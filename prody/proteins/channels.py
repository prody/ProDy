# -*- coding: utf-8 -*-

"""This module is called CaviTracer and defines functions for calculating 
channels, tunnels, and surface cavities within protein structure.
"""

__author__ = 'Karolina Mikulska-Ruminska', 'Jan Brezovsky', 'Eryk Trzcinski'
__credits__ = ['Karolina Mikulska-Ruminska', 'Jan Brezovsky', 'Eryk Trzcinski']
__email__ = ['karolamik@fizyka.umk.pl']

import logging
from contextlib import contextmanager

import numpy as np
from prody import LOGGER, PY3K
from prody.atomic import Atomic
from prody.utilities import checkCoords, getCoords, isListLike
from prody.proteins import writePDB, parsePDB, parsePQR
from prody.ensemble import Ensemble
from prody.measure import calcCenter, calcTransformation, calcDistance, calcRMSD, superpose


__all__ = ['getVmdModel', 'calcChannels', 'calcChannelsMultipleFrames', 
           'getChannelParameters', 'getChannelAtoms', 'showChannels', 
           'showCavities', 'showSurfaceCavities', 'selectChannelBySelection', 
           'getChannelResidueNames',
           'calcChannelSurfaceOverlaps', 'calcSurfaceCavities', 
           'calcSurfaceCavitiesMultipleFrames', 'getSurfaceCavityParameters',
           'getSurfaceCavityResidueNames', 'selectSurfaceCavityBySelection',
           'calcSurfaceCavityOverlaps',
           'getSurfaceCavityResidueNamesMultipleFrames',
           'getSurfaceCavityParametersMultipleFrames', 
           'getChannelParametersMultipleFrames', '_reportAtomsInputComposition',
           'getChannelResidueNamesMultipleFrames', 'calcPoresFromChannels',
           'showPores']

# Sampling of the enclosure test used to strip the moat (see
# ChannelCalculator.calcEnclosure). These are constants, not knobs: the enclosure
# of a point depends on how many directions are sampled and how far they are
# followed, so min_enclosure is only meaningful against a fixed sampling. Adding
# rays lowers every enclosure, since more directions find more of the thin ways
# out of a channel, and so invalidates the threshold rather than refining it.
ENCLOSURE_RAYS = 32
ENCLOSURE_RANGE = 25.0
ENCLOSURE_STEP = 0.75
# One radius for every atom, on the scale of a heavy-atom vdW radius. Enclosure is
# a burial heuristic, so resolving 1.52 A from 1.7 A would only shift every value
# by a little and be absorbed by min_enclosure; a single radius means a single
# tree and a plain nearest-neighbour test.
ENCLOSURE_RADIUS = 1.7

_OVERLAP_OFFSET_CACHE = {}

@contextmanager
def _warningsDelivered():
    """Let WARNING records through for the duration of the block.

    Importing ProDy installs a logging filter that drops every WARNING record from
    the package logger (prody.dynamics.adaptive2, at module scope), so every
    ``LOGGER.warn`` in ProDy is silently discarded. Whether that filter should exist
    at all is a question for the package as a whole; until it is settled, this
    module at least delivers its own warnings.

    Offending filters are found by asking them, rather than by importing the class
    and matching on type: any filter that rejects a synthetic WARNING record is
    detached for the block and reinstated afterwards. That keeps this working if the
    filter is renamed or moved, and it leaves the filter in force everywhere else,
    so nothing outside this module changes behaviour."""
    logger = LOGGER._logger
    probe = logging.LogRecord(logger.name, logging.WARNING, __file__, 0,
                              '', (), None)
    muting = [f for f in list(logger.filters)
              if not (f.filter(probe) if hasattr(f, 'filter') else f(probe))]
    for f in muting:
        logger.removeFilter(f)
    try:
        yield
    finally:
        for f in muting:
            logger.addFilter(f)


def _warn(message):
    """``LOGGER.warn``, but actually emitted. See :func:`_warningsDelivered`."""
    with _warningsDelivered():
        LOGGER.warn(message)


def checkAndImport(package_name):
    """Check for package and import it if possible and return **True**.
    Otherwise, return **False
        
    :arg package_name: name of package
    :type package_name: str

    :arg import_command: optional command to import submodules or with an alias
        default **None** means use "import {0}".format(package_name)
    :type import_command: None, str """
    
    if not isinstance(package_name, str):
        raise TypeError('package_name should be a string')

    if PY3K:
        import importlib.util
        if importlib.util.find_spec(package_name) is None:
            _warn("Package " + str(package_name) + " is not installed. "
            "Please install it to use this function.")
            return False
    else:
        try:
            __import__(package_name)
        except ImportError:
            _warn("Package " + str(package_name) + " is not installed. "
            "Please install it to use this function.")
            return False
    
    return True


def _getOverlapSphereOffsets(radius, resolution):
    """Return integer voxel offsets inside a sphere of a given radius."""

    key = (round(float(radius), 3), float(resolution))

    if key in _OVERLAP_OFFSET_CACHE:
        return _OVERLAP_OFFSET_CACHE[key]

    n = int(np.ceil(radius / resolution))
    grid = np.arange(-n, n + 1, dtype=int)

    dx, dy, dz = np.meshgrid(grid, grid, grid, indexing='ij')
    offsets = np.vstack((dx.ravel(), dy.ravel(), dz.ravel())).T
    xyz = offsets.astype(float) * resolution
    mask = np.sum(xyz * xyz, axis=1) <= radius * radius
    offsets = offsets[mask]
    _OVERLAP_OFFSET_CACHE[key] = offsets

    return offsets


def _surfaceFromPqrWorker(args):
    """Create voxelized FIL surface for one PQR file."""

    pqr_file, resolution = args
    atoms = parsePQR(pqr_file)
    fil = atoms.select('resname FIL')

    if fil is None:
        return set()

    coords = fil.getCoords()
    radii = fil.getRadii()
    surface = set()

    for center, radius in zip(coords, radii):
        center_idx = np.rint(center / resolution).astype(int)
        offsets = _getOverlapSphereOffsets(float(radius), resolution)
        voxels = offsets + center_idx
        surface.update(map(tuple, voxels))

    return surface


def _reportAtomsInputComposition(atoms):
    """Report the composition of atoms supplied for channel analysis.

    This function checks whether the input atomic structure contains only
    protein atoms or also includes water, non-water HETATM records, or other
    non-protein components. If non-protein atoms are present, a warning is
    issued indicating that all supplied atoms will be included in the channel
    calculation.

    The function does not modify or filter the input structure. To analyze only
    the protein, the user should provide an appropriate ProDy selection, for
    example ``atoms.select('protein')``. """
    
    if not isinstance(atoms, Atomic):
        raise TypeError(
            "atoms must be a ProDy Atomic object, such as an AtomGroup "
            "or Selection")

    protein = atoms.select('protein')
    water = atoms.select('water')
    hetero = atoms.select('hetero and not water')
    other = atoms.select('not protein and not hetero')
    nonprotein = atoms.select('not protein')

    if nonprotein is None:
        LOGGER.info("The atoms supplied to calcChannels contain protein atoms only.")
        return

    components = []

    if water is not None:
        components.append(
            "water: {0} atoms in {1} residues".format(
                water.numAtoms(),
                len(np.unique(water.getResindices()))))

    if hetero is not None:
        components.append(
            "non-water hetero components: {0} atoms "
            "(resnames: {1})".format(
                hetero.numAtoms(),
                ", ".join(sorted(np.unique(hetero.getResnames())))))

    if other is not None:
        components.append(
            "other non-protein components: {0} atoms "
            "(resnames: {1})".format(
                other.numAtoms(),
                ", ".join(sorted(np.unique(other.getResnames())))))

    _warn("The atoms supplied to calcChannels() contain non-protein components: "
        "{0}. All supplied atoms except waters will be used for channel analysis. "
        "To analyze only the protein structure, provide an appropriate "
        "selection, for example atoms.select('protein').".format(
            "; ".join(components)))


def getVmdModel(vmd_path, atoms, representation='NewCartoon'):
    """Generates a 3D model of molecular structures using VMD and returns 
    it as an Open3D TriangleMesh.

    This function creates a temporary PDB file from the provided atomic data
    and uses VMD (Visual Molecular Dynamics) to render this data into an STL
    file, which is then loaded into Open3D as a TriangleMesh. The function
    handles the creation and cleanup of temporary files and manages the
    subprocess call to VMD.
    
    To install Open3D use: 
    conda install open3d (for Anaconda users; version open3d-0.19.0 was used 
    during the development) or pip install open3d

    If problem with `ipykernel.comm.Comm` class appeared while using getVmdModel
    please update dash: python -m pip install -U dash

    :arg vmd_path: Path to the VMD executable. This is required to run VMD and 
        execute the TCL script.
    :type vmd_path: str

    :arg atoms: Atomic data to be written to a PDB file. This should be an 
        object or data structure that is compatible with the `writePDB` function.
    :type atoms: object

    :raises ImportError: If required libraries ('subprocess', 'pathlib', 
        'tempfile', 'open3d') are not installed, an ImportError is raised, 
        specifying which libraries are missing.

    :raises ValueError: If the STL file is not created or is empty, or if the 
        STL file cannot be read as a TriangleMesh,
        a ValueError is raised.

    :returns: An Open3D TriangleMesh object representing the 3D model generated
        from the PDB data.
    :rtype: open3d.geometry.TriangleMesh

    Example usage:
    model = getVmdModel('/path/to/vmd', atoms) """

    required = ['subprocess', 'pathlib', 'tempfile', 'open3d']
    missing = []
    errorMsg = None
    for name in required:
        if not checkAndImport(name):
            missing.append(name)
            if errorMsg is None:
                errorMsg = 'To run getVmdModel, ' \
                'please install {0}'.format(missing[0])
            else:
                errorMsg += ', ' + name

    if len(missing) > 0:
        if len(missing) > 1:
            errorMsg = ', '.join(errorMsg.split(', ')[:-1]) + ' and ' + errorMsg.split(', ')[-1]
        raise ImportError(errorMsg)

    import subprocess
    import tempfile
    import open3d as o3d
    import os

    representation_map = {
    'newcartoon': 'NewCartoon',
    'cartoon': 'NewCartoon',
    'vdw': 'VDW 1.0 20.0',
    'surf': 'Surf',
    'quicksurf': 'QuickSurf 1.0 0.5 0.25 2.0',
    'cpk': 'CPK 1.0 0.3 20.0 20.0'}

    rep_key = representation.lower()
    if rep_key not in representation_map:
        raise ValueError(
            "representation must be one of: 'NewCartoon', 'VDW', 'Surf', " \
            "'QuickSurf', or 'CPK'")
    representation_style = representation_map[rep_key]

    if PY3K:
        from pathlib import Path
    else:
        Path = lambda x: x 

    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as temp_pdb:
        temp_pdb_path = Path(temp_pdb.name)
        writePDB(temp_pdb.name, atoms)

    with tempfile.NamedTemporaryFile(suffix=".tcl", delete=False) as temp_script:
        temp_script_path = Path(temp_script.name)

        if PY3K:
            output_path = temp_script_path.parent / "output.stl"
        else:
            output_path = os.path.join(os.path.dirname(temp_script.name), 
                                       "output.stl")

        vmd_script = """
        set file_path [lindex $argv 0]
        set output_path [lindex $argv 1]

        mol new $file_path
        mol modstyle 0 0 %s

        set id_matrix {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}
        molinfo top set center_matrix [list $id_matrix]
        molinfo top set rotate_matrix [list $id_matrix]
        molinfo top set scale_matrix [list $id_matrix]

        rendering_method stl
        render STL $output_path

        exit
        """ % representation_style
        
        temp_script.write(vmd_script.encode('utf-8'))

    command = [vmd_path, '-e', str(temp_script_path), '-args', 
               str(temp_pdb_path), str(output_path)]

    try:
        if PY3K:
            subprocess.run(command, check=True)
        else:
            returncode = subprocess.call(command)
            if returncode != 0:
                LOGGER.info("VMD exited with status " + str(returncode) + ".")
    except Exception as e:
        _warn("An unexpected error occurred: " + str(e))
    finally:
        if os.path.exists(temp_script_path):
            os.unlink(temp_script_path)
        if os.path.exists(temp_pdb_path):
            os.unlink(temp_pdb_path)

        if not os.path.exists(output_path) or os.stat(output_path).st_size == 0:
            raise ValueError("STL file was not created or is empty.")

        stl_mesh = o3d.io.read_triangle_mesh(str(output_path))

        if stl_mesh.is_empty():
            raise ValueError("Failed to read the STL file as a TriangleMesh.")

        if os.path.exists(output_path):
            os.unlink(output_path)

        LOGGER.info("Model created successfully.")
        return stl_mesh


def showChannels(channels, model=None, surface=None):
    """Visualizes the channels, and optionally, the molecular model and 
    surface, using Open3D.
    
    This function renders a 3D visualization of molecular channels based on 
    their spline representations. It can also display a molecular model (e.g., 
    the protein structure) and a surface (e.g., cavity surface) in the same 
    visualization. The function utilizes the Open3D library to create and 
    render the 3D meshes.

    To install Open3D use: 
    conda install open3d (for Anaconda users; version open3d-0.19.0 was used 
    during the development) or pip install open3d
    
    :arg channels: A list of channel objects or a single channel object. Each 
        channel should have a `getSplines()` method that returns two 
        CubicSpline objects: one for the centerline and one for the radii.
    :type channels: list or single channel object
    
    :arg model: An optional Open3D TriangleMesh object representing the 
        molecular model, such as a protein. If provided, this model will be 
        rendered in the visualization.
        Model can be generated using getVmdModel() function.
    :type model: open3d.geometry.TriangleMesh, optional
    
    :arg surface: An optional list containing the surface data. The list should
         have two elements:
        - `points`: The coordinates of the vertices on the surface.
        - `simp`: The simplices that define the surface (e.g., triangles or 
           tetrahedra).
        If provided, the surface will be rendered as a wireframe overlay in the
        visualization.
    :type surface: list (with two numpy arrays), optional
    
    :raises ImportError: If the Open3D library is not installed, an ImportError
        is raised, prompting the user to install Open3D.
    
    :returns: None. This function only renders the visualization.
    
    Example usage:
    showChannels(channels, model=protein_mesh, surface=surface_data) """
    
    if not checkAndImport('open3d'):
        errorMsg = 'To run showChannels, please install open3d (version 0.19.0).'
        raise ImportError(errorMsg)
            
    import open3d as o3d
    
    def create_mesh_from_spline(centerline_spline, radius_spline, n=5):
        N = n * len(centerline_spline.x)
        t = np.linspace(centerline_spline.x[0], centerline_spline.x[-1], N)
        centers = centerline_spline(t)
        radii = radius_spline(t)

        spheres = [
            o3d.geometry.TriangleMesh.create_sphere(radius=r,
                                                    resolution=20).translate(c)
            for r, c in zip(radii, centers)
        ]
        mesh = spheres[0]
        for sphere in spheres[1:]:
            mesh += sphere

        return mesh
    
    if not isinstance(channels, list):
        channels = [channels]
    
    channel_meshes = [create_mesh_from_spline(*channel.getSplines()) for channel in channels]
    meshes_to_visualize = [o3d.geometry.TriangleMesh.create_coordinate_frame(size=0.1, origin=[0, 0, 0])]

    if model is not None:
        model.compute_vertex_normals()
        model.paint_uniform_color([0.1, 0.7, 0.3])
        meshes_to_visualize.append(model)
            
    if channel_meshes is not None:
        if not isinstance(channel_meshes, list):
            channel_meshes = [channel_meshes]
        for channel_mesh in channel_meshes:
            channel_mesh.compute_vertex_normals()
            channel_mesh.paint_uniform_color([0.5, 0.0, 0.5])
        meshes_to_visualize.extend(channel_meshes)
        
    if surface is not None:
        points = surface[0]
        simp = surface[1]
        
        triangles = []
        for tetra in simp:
            triangles.extend([sorted([tetra[0], tetra[1], tetra[2]]),
                                  sorted([tetra[0], tetra[1], tetra[3]]),
                                  sorted([tetra[0], tetra[2], tetra[3]]),
                                  sorted([tetra[1], tetra[2], tetra[3]])])
            
        triangles = np.array(triangles)
        triangles.sort(axis=1)
            
        triangles_tuple = [tuple(tri) for tri in triangles]
        unique_triangles, counts = np.unique(triangles_tuple, 
                                             return_counts=True, axis=0)
            
        surface_triangles = unique_triangles[counts == 1]
            
        lines = []
        for simplex in surface_triangles:
            for i in range(3):
                for j in range(i + 1, 3):
                    lines.append([simplex[i], simplex[j]])
            
        line_set = o3d.geometry.LineSet()
        line_set.points = o3d.utility.Vector3dVector(points)
        line_set.lines = o3d.utility.Vector2iVector(lines)
        
        meshes_to_visualize.append(line_set)
        
    if len(meshes_to_visualize) > 1:
        o3d.visualization.draw_geometries(meshes_to_visualize)
    else:
        LOGGER.info("Nothing to visualize.")


def showPores(pores, model=None, show_surface=False, surface=None, **kwargs):
    """Visualize pores calculated with :func:`calcPoresFromChannels`.

    :arg pores: Pore or sequence of Pore objects to visualize.
    :type pores: Pore or list """

    return showChannels(pores, model=model, surface=surface, **kwargs)


def showCavities(surface, show_surface=False):
    """Visualizes the cavities within a molecular surface using Open3D.

    This function displays a 3D visualization of cavities detected in a 
    molecular structure.
    It uses the Open3D library to render the cavities as a triangle mesh. 
    Optionally, it can also display the molecular surface as a wireframe 
    overlay.

    To install Open3D use: 
    conda install open3d (for Anaconda users; version open3d-0.19.0 was used 
    during the development) or pip install open3d

    :arg surface: A list containing three elements:
        - `points`: The coordinates of the vertices (atoms) in the molecular 
        structure.
        - `surf_simp`: The simplices that define the molecular surface.
        - `simp_cavities`: The simplices corresponding to the detected cavities.
    :type surface: list (with three numpy arrays)

    :arg show_surface: A boolean flag indicating whether to display the 
        molecular surface
        as a wireframe overlay in the visualization. If True, the surface will 
        be displayed in addition to the cavities. Default is False.
    :type show_surface: bool

    :raises ImportError: If the Open3D library is not installed, an ImportError
        is raised, prompting the user to install Open3D.

    :returns: None

    Example usage:
    showCavities(surface_data, show_surface=True) """
    
    if not checkAndImport('open3d'):
        errorMsg = 'To run showChannels, please install open3d.'
        raise ImportError(errorMsg)
            
    import open3d as o3d
    
    points = surface[0]
    surf_simp = surface[1]
    simp_cavities = surface[2]
    
    triangles = []
    for tetra in simp_cavities:
        triangles.extend([sorted([tetra[0], tetra[1], tetra[2]]),
                          sorted([tetra[0], tetra[1], tetra[3]]),
                          sorted([tetra[0], tetra[2], tetra[3]]),
                          sorted([tetra[1], tetra[2], tetra[3]])])
        
    surface_triangles = np.unique(np.array(triangles), axis=0, 
                                  return_counts=True)[0]
        
    mesh = o3d.geometry.TriangleMesh()
    mesh.vertices = o3d.utility.Vector3dVector(points)
    mesh.triangles = o3d.utility.Vector3iVector(surface_triangles)
        
    mesh.compute_vertex_normals()
    mesh.paint_uniform_color([0.1, 0.7, 0.3])
        
    vis = o3d.visualization.Visualizer()
    vis.create_window()
    vis.add_geometry(mesh)
        
    if show_surface == True:
        triangles = []
        for tetra in surf_simp:
            triangles.extend([sorted([tetra[0], tetra[1], tetra[2]]),
                                  sorted([tetra[0], tetra[1], tetra[3]]),
                                  sorted([tetra[0], tetra[2], tetra[3]]),
                                  sorted([tetra[1], tetra[2], tetra[3]])])
            
        triangles = np.array(triangles)
        triangles.sort(axis=1)
            
        triangles_tuple = [tuple(tri) for tri in triangles]
        unique_triangles, counts = np.unique(triangles_tuple, 
                                             return_counts=True, axis=0)
            
        surface_triangles = unique_triangles[counts == 1]
            
        lines = []
        for simplex in surface_triangles:
            for i in range(3):
                for j in range(i + 1, 3):
                    lines.append([simplex[i], simplex[j]])
            
        line_set = o3d.geometry.LineSet()
        line_set.points = o3d.utility.Vector3dVector(points)
        line_set.lines = o3d.utility.Vector2iVector(lines)
        
        vis.add_geometry(line_set)
            
    vis.get_render_option().mesh_show_back_face = True
    vis.get_render_option().background_color = np.array([1, 1, 1])
    vis.update_renderer()
    vis.run()
    vis.destroy_window()


def showSurfaceCavities(surface, cavities=None, model=None, show_surface=False, 
    cavity_atoms=None, mode='tetra', alpha=4.0, smoothing=0):
    """Visualize surface cavities together with an optional protein model.
    
    This function displays surface cavities calculated by :func:`calcSurfaceCavities`.
    Cavities can be visualized either directly from the tetrahedral/Voronoi
    representation returned by the calculation, or from pseudoatoms loaded from
    a PDB/PQR file through `cavity_atoms`.

    - `mode='tetra'` displays cavities as a tetrahedron-derived mesh.
      This mode follows the original geometric representation most closely,
      but the resulting surface may appear faceted.
    - `mode='smooth'` builds a smoother surface from Voronoi vertices
      assigned to each cavity using an alpha-shape reconstruction. This mode
      requires `cavities` to be provided.

    If `cavity_atoms` is provided, cavities are visualized from the
    coordinates of the supplied pseudoatoms instead of from `surface[2]` or
    `cavities`. 
    
    :arg surface: Surface data returned by :func:`calcSurfaceCavities`.
        Required for `mode='tetra'`, for `mode='smooth'` when
        `cavity_atoms` is not provided, and for displaying the molecular
        surface when `show_surface=True`. The expected list contains:

        - `surface[0]`: atomic coordinates used for the calculation,
        - `surface[1]`: simplices defining the molecular surface,
        - `surface[2]`: merged cavity simplices,
        - `surface[3]`: simplices after surface-layer removal,
        - `surface[4]`: Voronoi vertices used to represent surface cavities.

    :type surface: list or None

    :arg cavities: List of :class:`Cavity` objects returned by
        :func:`calcSurfaceCavities`. Required when `mode='smooth'` and
        `cavity_atoms` is not provided, because the function uses
        `cavity.tetrahedra` to select the corresponding Voronoi vertices.
    :type cavities: list or None

    :arg model: Optional Open3D `TriangleMesh` representing the protein or
        another molecular model. The model can be generated with
        :func:`getVmdModel`.
    :type model: open3d.geometry.TriangleMesh, or None

    :arg show_surface: If `True`, display the molecular surface wireframe
        derived from `surface[1]` in addition to the cavity representation.
        This requires `surface` to be provided. Default is `False`.
    :type show_surface: bool

    :arg mode: Visualization mode used when `cavity_atoms` is not provided.
        Accepted values are `'tetra'` and `'smooth'`. Default is `'tetra'`.
    :type mode: str

    :arg alpha: Alpha value used for alpha-shape surface reconstruction in
        `mode='smooth'` and when visualizing `cavity_atoms`. Smaller values
        produce tighter surfaces, while larger values may connect more distant
        points and generate broader surfaces. Default is 4.0.
    :type alpha: float

    :arg smoothing: Number of Taubin smoothing iterations applied to the
        reconstructed cavity mesh. If `0` or `None`, no smoothing is
        applied. Default is 0.
    :type smoothing: int or None

    :arg cavity_atoms: Optional pseudoatom representation of surface
        cavities. This can be either a path to a PDB/PQR file or a parsed ProDy
        `AtomGroup`, or an Open3D `TriangleMesh` generated, for example, with 
        :func:`getVmdModel`.
    :type cavity_atoms: str, :class:`.AtomGroup`, open3d.geometry.TriangleMesh,
         or None
    
    Examples:
    p = parsePDB('1tqn')
    protein = p.select('protein')
    cavities, surface = calcSurfaceCavities(protein, output_path='cavities.pqr')
    model_protein = getVmdModel(vmd_path, protein)
    cav_model = getVmdModel(vmd_path, parsePQR('cavities.pqr'), representation='QuickSurf')
    
    showSurfaceCavities(surface, model=model_protein, show_surface=True)
    or
    showSurfaceCavities(surface, model=model_protein, cavity_atoms='cavities.pqr', show_surface=True)
    or
    showSurfaceCavities(surface, model=model_protein, cavity_atoms=cav_model, show_surface=True)
    """

    if not checkAndImport('open3d'):
        raise ImportError('To run showSurfaceCavities, please install open3d.')

    import open3d as o3d
    import os

    points = surface[0]
    surf_simp = surface[1]
    simp_cavities = surface[2]
    
    if mode not in ['tetra', 'smooth']:
        raise ValueError("mode must be 'tetra' or 'smooth'")

    meshes_to_visualize = []

    if model is not None:
        model.compute_vertex_normals()
        model.paint_uniform_color([0.8, 0.8, 0.8])
        meshes_to_visualize.append(model)
    
    if cavity_atoms is not None:
        cavity_atoms_given = True
    
        if isinstance(cavity_atoms, o3d.geometry.TriangleMesh):
            cavity_atoms.compute_vertex_normals()
            cavity_atoms.paint_uniform_color([0.1, 0.7, 0.3])
            meshes_to_visualize.append(cavity_atoms)
        
        else:  
            if isinstance(cavity_atoms, str):
                ext = os.path.splitext(cavity_atoms)[1].lower()
                if ext == '.pqr':
                    cavity_atoms = parsePQR(cavity_atoms)
                else:
                    cavity_atoms = parsePDB(cavity_atoms)

            if not hasattr(cavity_atoms, 'getCoords'):
                raise TypeError("cavity_atoms must be a PDB/PQR filename, a " 
                                "ProDy AtomGroup,or an Open3D TriangleMesh.")

            resnums = np.unique(cavity_atoms.getResnums())

            for resnum in resnums:
                sele = cavity_atoms.select('resnum {0}'.format(resnum))
                if sele is None:
                    continue

                pts = sele.getCoords()
                if pts is None or len(pts) < 4:
                    continue

                pcd = o3d.geometry.PointCloud()
                pcd.points = o3d.utility.Vector3dVector(pts)
                cavity_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(pcd, alpha)

                if smoothing is not None and smoothing > 0:
                    cavity_mesh = cavity_mesh.filter_smooth_taubin(number_of_iterations=smoothing)

                cavity_mesh.compute_vertex_normals()
                cavity_mesh.paint_uniform_color([0.1, 0.7, 0.3])
                meshes_to_visualize.append(cavity_mesh)

    else:
        if surface is None:
            raise ValueError("surface must be provided when cavity_atoms is not given")

        points = surface[0]
        surf_simp = surface[1]
        simp_cavities = surface[2]
    
        if mode == 'tetra':
            triangles = []
            for tetra in simp_cavities:
                triangles.extend([
                    sorted([tetra[0], tetra[1], tetra[2]]),
                    sorted([tetra[0], tetra[1], tetra[3]]),
                    sorted([tetra[0], tetra[2], tetra[3]]),
                    sorted([tetra[1], tetra[2], tetra[3]])])

            surface_triangles = np.unique(np.array(triangles), axis=0, 
                                          return_counts=True)[0]
            cavity_mesh = o3d.geometry.TriangleMesh()
            cavity_mesh.vertices = o3d.utility.Vector3dVector(points)
            cavity_mesh.triangles = o3d.utility.Vector3iVector(surface_triangles)
            cavity_mesh.compute_vertex_normals()
            
            if smoothing is not None and smoothing > 0:
                cavity_mesh = cavity_mesh.filter_smooth_taubin(number_of_iterations=smoothing)
                cavity_mesh.compute_vertex_normals()
            
            cavity_mesh.paint_uniform_color([0.1, 0.7, 0.3])
            meshes_to_visualize.append(cavity_mesh)
            
        elif mode == 'smooth':
            if cavities is None:
                raise ValueError("cavities must be provided for mode='smooth'")
            
            verti = surface[4]
            for cavity in cavities:
                if cavity.tetrahedra is None or len(cavity.tetrahedra) < 4:
                    continue

                pts = verti[cavity.tetrahedra]
                pcd = o3d.geometry.PointCloud()
                pcd.points = o3d.utility.Vector3dVector(pts)
                cavity_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(pcd, alpha)

                if smoothing is not None and smoothing > 0:
                    cavity_mesh = cavity_mesh.filter_smooth_taubin(number_of_iterations=smoothing)

                cavity_mesh.compute_vertex_normals()
                cavity_mesh.paint_uniform_color([0.1, 0.7, 0.3])
                meshes_to_visualize.append(cavity_mesh)

    if show_surface == True:
        triangles = []
        for tetra in surf_simp:
            triangles.extend([sorted([tetra[0], tetra[1], tetra[2]]),
                                  sorted([tetra[0], tetra[1], tetra[3]]),
                                  sorted([tetra[0], tetra[2], tetra[3]]),
                                  sorted([tetra[1], tetra[2], tetra[3]])])
            
        triangles = np.array(triangles)
        triangles.sort(axis=1)
            
        triangles_tuple = [tuple(tri) for tri in triangles]
        unique_triangles, counts = np.unique(triangles_tuple, 
                                             return_counts=True, axis=0)
        surface_triangles = unique_triangles[counts == 1]
            
        lines = []
        for simplex in surface_triangles:
            for i in range(3):
                for j in range(i + 1, 3):
                    lines.append([simplex[i], simplex[j]])
            
        line_set = o3d.geometry.LineSet()
        line_set.points = o3d.utility.Vector3dVector(points)
        line_set.lines = o3d.utility.Vector2iVector(lines)
        meshes_to_visualize.append(line_set)

    o3d.visualization.draw_geometries(meshes_to_visualize)

def calcChannels(atoms, output_path=None, separate=False, start_point=None,
    restrict_channels_to_start_point=True, start_point_search=3.0,
    r1=3, r2=0.9, min_depth=5,
    min_volume=None, max_volume=None, max_depth=None, bottleneck=0.0,
    sparsity=1, min_tetrahedra=None, max_tetrahedra=None, cavities_only=False,
    diagram="homogenized", max_deviation=0.1, truncate_at_surface=True,
    similarity=0.8, route_tolerance=1.0, min_enclosure=0.70, max_peel_depth=None,
    weighted_cache=True, weighted_mouth_depth=2.5, edge_cost=None,
    return_details=False):
    """Computes and identifies channels within a molecular structure using 
    Voronoi and Delaunay tessellations.

    This function analyzes the provided atomic structure to detect channels, 
    which are voids or pathways within the molecular structure. It employs 
    Voronoi and Delaunay tessellations to identify these regions, then filters
    and refines the detected channels based on various parameters such as the 
    minimum depth and bottleneck size. The results can be saved to a PQR file 
    (PDB is optional) if an output path is provided. The `separate` parameter 
    controls whether each detected channel is saved to a separate file or if all 
    channels are saved in a single file.

    The implementation is inspired by the methods described in the following 
    publications:
    "MOLE 2.0: advanced approach for analysis of biomacromolecular channels" by
     D. Sehnal, et al., published in J Chemoinform, 5 (39) 2013.
     
     "CAVER: Algorithms for Analyzing Dynamics of Tunnels in Macromolecules". by
     A. Pavelka, et al., published in IEEE ACM T COMPUT BI, (13) 2016.

    "Software Tools for Identification, Visualization and Analysis of Protein 
    Tunnels and Channels". by J. Brezovsky, et al., Biotechnol Adv (31) 2013.


    :arg atoms: An object representing the molecular structure, typically 
        containing atomic coordinates and element types.
    :type atoms: `Atoms` object

    :arg output_path: Optional path to save the resulting channels and 
        associated data in PQR (or PDB) format. If None, results are not saved. 
        Default is None.
    :type output_path: str or None

    :arg separate: If True, each detected channel is saved to a separate PDB 
        file. If False, all channels are saved in a single PDB file. Default is 
        False.
    :type separate: bool

    :arg start_point: Optional starting point for channel search. This can be
        either a 3D coordinate point or an atomic selection/AtomGroup. If the
        3D coordinate point will be provided, the algorithm seeds the channel
        search near this point (overriding the default automatic seed selection
        based on the deepest tetrahedron); see ``start_point_search`` for how
        the seed tetrahedron itself is picked. Coordinates must be given in Å.
        If an atomic selection is provided, its geometric center is used as the
         starting point.
    :type start_point: list, tuple, or ndarray (length 3), :class:`.Atomic`, or None

    :arg restrict_channels_to_start_point: Only used when ``start_point`` is
        provided. If True (default), the channel search is restricted to the
        single cavity whose closest tetrahedron is globally nearest to
        ``start_point``, so  channels are computed only for the region around
        that point instead of one channel bundle per detected cavity. If False,
        ``start_point`` merely overrides the seed (starting) tetrahedron of
        every cavity and channels are still computed for all cavities.
    :type restrict_channels_to_start_point: bool

    :arg start_point_search: Only used when ``start_point`` is provided. Radius,
        in Angstrom, of the neighbourhood of ``start_point`` searched for the seed
        tetrahedron. The tetrahedron nearest ``start_point`` is often a tight one,
        and since every channel of the cavity starts there, its inscribed radius
        caps all of their bottlenecks and appears as one shared bottleneck at the
        joint beginning of the bundle. Seeded instead is the widest tetrahedron within
        ``start_point_search`` of ``start_point`` that belongs to the same cavity, is
        no shallower than the nearest one (so the seed cannot drift out towards the
        mouth) and is reachable from it through that neighbourhood (so it stays in the
        void the point sits in rather than crossing a wall). Default is 3.0; use 0 to
        seed the nearest tetrahedron as-is.
    :type start_point_search: float

    :arg r1: The first radius threshold used during the deletion of simplices, 
        which is used to define the outer surface of the channels. Default is 3
    :type r1: float

    :arg r2: The second radius threshold used to define the inner surface of
        the channels. Default is 0.9.

        Below about 1.2 Angstrom the probe is smaller than a water molecule, and
        then the structure must carry explicit hydrogens. Without them, every
        carbon keeps its full vdW radius while the space its hydrogens occupied is
        left empty, and a sub-water probe is small enough to thread those
        interstices: the interior percolates into a sponge and the channel count
        can rise several-fold. At 1.2 Angstrom and above, protonated and 
        unprotonated structures give the same channels, and an X-ray file may 
        be used as it comes. A warning is issued for the unsafe combination.
    :type r2: float

    :arg min_depth: The minimum depth, in Angstrom, a cavity must reach to be
        considered as a channel. Depth is the geodesic distance from the cavity's
        surface opening to its farthest point along the Voronoi network, so it is a
        physical length independent of the tessellation density. Default is 5.
    :type min_depth: float

    :arg max_depth: Maximum cavity depth, in Angstrom. Portions of a cavity deeper
        than this value are trimmed away. Default is None (no trimming).
    :type max_depth: float

    :arg bottleneck: Acts as secondary filter following channel identification.
        The minimum allowed bottleneck size (narrowest point) for the channels.
        It it critical when diagram=simple, as it partially corrects for wrong 
        diagram topology. Default is 0.0, no filtering applied. 
    :type bottleneck: float

    :arg min_volume: Minimum volume required for a channel/cavity to be 
        retained. Default is None.
    :type min_volume: float

    :arg max_volume: Maximum volume allowed for a channel/cavity to be 
        retained. Default is None.
    :type max_volume: float

    :arg sparsity: Size of a channel surface opening (mouth), in Angstrom: how far
        apart two exits must lie to count as separate openings. It is one quantity,
        reached by whichever branch of the search is running, and the two branches
        are mutually exclusive. With ``truncate_at_surface`` True (the default) it
        is a floor on the radius of a reported opening, so two channels leaving
        closer than ``sparsity`` are treated as sharing that opening and are merged
        if they also share a corridor (see ``similarity``); being applied *after*
        the search, it can only merge channels there, never hide one, and is a
        reporting preference rather than part of the geometry. With
        ``truncate_at_surface`` False it is instead the spacing at which exit
        tetrahedra are sampled as channel termini, and it does then decide which
        channels are found at all. Either way a higher value reports fewer channels.
        It has no effect on the cavities, which are found from the exit tetrahedra
        before any thinning. Default is 1.
    :type sparsity: float

    :arg diagram: 
        "homogenized" (default) - every atom is substituted by a set of 
        homogeneous balls whose common radius equals the smallest van der Waals
        radius present in the structure, before building the Voronoi and 
        Delaunay tessellations. This yields an accurate estimate of the 
        additively weighted Voronoi diagram from an ordinary one, as done in 
        MolAxis and CAVER 3
        "simple" - the original atoms are used with their individual van der 
        Waals radii directly. This is very inaccurate and should be avoided in 
        almost all cases.
        "weighted" - the true additively-weighted (Apollonius) Voronoi diagram is
        built directly from the atoms with their individual van der Waals radii,
        using the third-party ``vorpy`` package (with a compiled kernel when
        ``numba`` is available). This is the exact diagram the "homogenized" mode
        approximates, at a higher cost (~ 100x slower, for 4000 atoms). To use this
        approach install ``vorpy`` library using ``pip install vorpy3``. 
    :type diagram: str

    :arg max_deviation: Maximum tolerated deviation, in Angstrom, between the 
        union surface of the substitute balls and the original van der Waals 
        surface when ``diagram = homogenized`` . It controls the trade-off
        between surface accuracy and the number of balls generated: an atom 
        whose radius exceeds the smallest radius (``rho``) by more than 
        ``max_deviation`` is filled with several balls, otherwise it is kept as
         a single ``rho`` ball. Default is 0.1. Guideline values:

        * ``0.1`` fine accurate surface with minimal errors, but on 13-15x 
            more balls than original
        * ``0.15`` in heavy-atom-only structures it starts filling carbon, which
             is otherwise left as a single ball with a uniform ~0.18 A inset).
        * ``0.2`` speed optimized ; e.g. carbon fills to ~15 balls when 
            hydrogens are present (``rho``=1.2), resulting roughly to ~8 times 
            more balls. Without hydrogens, carbons are single balls.

        Only used when ``diagram = homogenized``.
    :type max_deviation: float

    :arg truncate_at_surface: If True (default), surface (exit) tetrahedra are
        made *absorbing*: a channel may end at one, but no channel may pass
        through one. A mouth is a surface tetrahedron a probe of the traversal
        radius ``r2`` can leave through. This forbids the cheapest path from
        surfacing at one mouth, running along the outside and re-entering at
        another - a surface hop, not a tunnel - which the width-rewarding cost
        would otherwise prefer, since surface grooves are the widest space
        available. Enforcing it during the search (rather than cutting the
        winning path afterwards) is what keeps genuine narrow interior corridors
        in the output: cut afterwards, such a corridor loses the cheapest-path
        race to the groove leading to the same mouth and is never enumerated at
        all. If False, paths run freely to their end tetrahedra, surface hops
        included.
    :type truncate_at_surface: bool

    :arg similarity: Only used when ``truncate_at_surface`` is True. Fraction
        (0-1) of the **longer** of two channels, measured in Angstrom along its
        centerline, that must run within ``route_tolerance`` of the other one for
        the two to count as the same corridor. Two channels are merged (cheapest
        kept) only when they take the same corridor **and** leave through the
        same opening (see ``sparsity``); a corridor that forks near the surface
        and exits twice is one tunnel, but two different corridors to one opening,
        or one corridor reaching two openings, are two tunnels. The comparison is
        geometric rather than a shared prefix of tetrahedra, so it is unaffected
        by *where* two routes diverge (variants that split and rejoin still count
        as one), by containment (a long route is not deleted as the "duplicate" of
        a short one it happens to start with), and by ``max_deviation`` (a
        tetrahedron count is not mesh-invariant; Angstrom are). ``1.0`` merges
        only routes that coincide along their whole length; ``0.0`` merges every
        channel that shares an opening. Default is 0.8.
    :type similarity: float

    :arg route_tolerance: Only used when ``truncate_at_surface`` is True. How far
        apart, in Angstrom, two centerlines may drift and still count as the same
        corridor when computing ``similarity``. Larger values merge more
        aggressively (nearby parallel routes read as one tunnel); smaller values
        report finer route variants separately. Default is 1.0.
    :type route_tolerance: float

    :arg min_enclosure: Fraction of directions that must be blocked by protein for
        a tetrahedron to count as interior, in ``[0, 1]``. Default is 0.70.

        Once the r1 surface is built it is eroded inward with the r2 probe, to
        strip the shell of true exterior that an r1 probe bridges over rather than
        enters (the "moat"); that shell would otherwise join the cavity and offer
        wide, low-cost routes along the outside of the protein. Erosion continues
        while the tetrahedra at the front are *open*, meaning that fewer than
        ``min_enclosure`` of the directions leaving them run into protein within
        :data:`ENCLOSURE_RANGE` Angstrom, and halts at the first buried layer.

        Bounding the erosion by size instead does not work. A count of tetrahedron
        layers is not mesh-invariant, as a layer is one tetrahedron thick and
        tetrahedra shrink as ``max_deviation`` is lowered. A depth in Angstrom is
        not ``r1``-invariant, as the moat is as deep as ``r1 - r2`` in a concavity
        but vanishes on a flat face, so a depth that clears it where it is thick
        also marches down the channel mouths and erodes the channels themselves.
        Testing burial locally instead leaves ``r1`` to decide only where the
        erosion starts, not where it stops, so results are independent of it, and
        ``r1`` is left doing the one job it should: capping the mouths.

        It is bounded from both sides, and the window is narrow.

        Too low and the moat is not fully stripped. Too high and the erosion 
        never meets a layer buried enough to stop it, so it percolates down the
          channels and eats the cavity. The default sits at the floor, which is
        the safe end: over-peeling deletes real channels, whereas
        under-peeling shows up as surface-riding routes that can be recognised.
    :type min_enclosure: float

    :arg max_peel_depth: Optional hard cap, **in Angstrom**, on how far the peel
        above may advance from the r1 surface. ``None`` (default) is uncapped, and
        the enclosure test alone decides where erosion stops. Set it only as a
        backstop on a structure where the peel misbehaves; it is deliberately not
        tied to ``r1``, since a cap that scales with ``r1`` reintroduces exactly
        the ``r1`` dependence that ``min_enclosure`` exists to remove.
    :type max_peel_depth: float or None

    :arg weighted_cache: Cache the raw additively-weighted Voronoi diagram to disk so
        that re-running ``diagram="weighted"`` on the same structure skips the
        expensive vorpy tessellation (which dominates the ~10 min run time). The
        diagram only depends on the atoms and ``r1``, so re-runs that change only
        ``r2``, ``bottleneck``, ``sparsity``, ``start_point`` etc. reuse it. ``True``
        (default) caches next to ``output_path`` (or, absent one, under the structure
        title in the current directory); pass a path to place it explicitly, or
        ``False`` to disable. The cache is keyed by content, so editing the structure
        or ``r1`` transparently forces a recompute. Only used for ``diagram="weighted"``.
    :type weighted_cache: bool or str

    :arg weighted_mouth_depth: Only used for ``diagram="weighted"``. The additively-
        weighted (Apollonius) tessellation is not a clean simplicial complex, so it
        leaves false interior boundary faces that the pipeline would misread as
        surface openings, truncating channels to stubs. To repair this a *homogenized*
        diagram of the same atoms is built as an interior/exterior oracle, and only
        exit tetrahedra whose Voronoi vertex lies within ``weighted_mouth_depth``
        Angstrom (geodesic distance below the molecular surface) are treated as mouths.
        Default 2.5 (the value at which the recovered channels match the
        ``"homogenized"`` result); ``None`` disables the relabeling.
    :type weighted_mouth_depth: float or None

    :arg edge_cost: How each Voronoi edge is priced in the Dijkstra tunnel search.
        ``"integral"`` prices each edge by the integral of its clearance profile
        along the edge, which is mesh-invariant. ``"bottleneck"`` uses the legacy
        ``length / (gate**2 + eps)``, charging the whole edge at its single
        narrowest point, whose value (and the routing it produces) drifts as
        ``max_deviation`` coarsens. ``None`` (default) selects ``"integral"`` for
        ``diagram="homogenized"``/``"simple"`` (straight edges, where the integral
        is exact) and ``"bottleneck"`` for ``diagram="weighted"``. ``"integral"``
        is rejected for ``diagram="weighted"`` (Apollonius edges are arcs the
        straight-chord integral cannot price). The reported bottleneck radius is
        unaffected by this choice.
    :type edge_cost: str or None

    :arg return_details: If True return an additional dictionary containing 
        internal calculation data, including the channel calculator, simplices,
        neighboring tetrahedra, Voronoi vertices, atomic coordinates, and van der
        Waals radii. Default is False.
    :type return_details: bool

    :returns: A tuple containing two elements:
        - `channels`: A list of detected channels, where each channel is an 
          object containing information about its path and geometry.
        - `surface`: A list containing additional information for further 
          visualization, including the atomic coordinates, simplices defining 
          the surface, and merged cavities.
    :rtype: tuple (list, list)

    This function performs the following steps:
    1. **Selection and Filtering:** Selects non-hetero atoms from the protein 
        and calculates van der Waals radii. When ``homogenize`` is True, each 
        atom is replaced by homogeneous balls of the smallest radius present 
        so that an ordinary tessellation approximates the additively weighted 
        Voronoi diagram. It then performs 3D Delaunay triangulation and Voronoi
         tessellation on the resulting coordinates.
    2. **State Management:** Creates and updates different stages of channel 
        detection of the protein structure to filter out simplices based on the
         given radii.
    3. **Surface Layer Calculation:** Determines the surface and second-layer 
        simplices from the filtered results.
    4. **Cavity and Channel Detection:** Finds and filters cavities based on 
        their depth and calculates channels using Dijkstra's algorithm.
    5. **Visualization and Saving:** Generates meshes for  detected channels, 
        filters them by bottleneck size, and either saves the results to a PDB 
        file or visualizes them based on the specified parameters.
       
    Example usage:
    channels, surface = calcChannels(atoms, output_path="channels", separate=True)
    
    channels, surface = calcChannels(atoms, output_path="all_channels.pdb", 
                                     start_point=[-22.312, -20.065, -11.144])
    
    start_sel = protein.select('resid 212 309 483')
    channels, surface = calcChannels(atoms, output_path="all_channels.pdb", 
                                     start_point=start_sel)
    
    To save the results as PDB file:
    channels, surface = calcChannels(atoms, output_path="channels.pdb",
                                     separate=False, r1=3, r2=0.9, min_depth=5,
                                     bottleneck=1, sparsity=3) """
    
    required = ['heapq', 'collections', 'scipy', 'pathlib', 'warnings']
    missing = []
    errorMsg = None
    for name in required:
        if not checkAndImport(name):
            missing.append(name)
            if errorMsg is None:
                errorMsg = 'To run calcChannels, please install {0}'.format(missing[0])
            else:
                errorMsg += ', ' + name

    if len(missing) > 0:
        if len(missing) > 1:
            errorMsg = ', '.join(errorMsg.split(', ')[:-1]) + ' and ' + errorMsg.split(', ')[-1]
        raise ImportError(errorMsg)

    if PY3K:
        from pathlib import Path
    else:
        from pathlib2 import Path
    
    if start_point is not None:
    
        if hasattr(start_point, 'getCoords'):
            if start_point.numAtoms() == 0:
                raise ValueError("start_point selection contains no atoms")
            start_point = calcCenter(start_point)

        elif not isListLike(start_point):
            raise TypeError("start_point must be a selection/AtomGroup or a list/tuple/ndarray "
                "with three numeric values")

        if len(start_point) != 3:
            raise ValueError(
                "start_point must be a selection/AtomGroup or a list of three numbers, e.g. "
                "start_point=[-12.312, 5.065, -1.144]")
        
        start_point = np.array(start_point, dtype=float)        
        
        LOGGER.info("Using user-provided start_point for channel seed: [{:.3f}, {:.3f}, {:.3f}] Å"
            .format(start_point[0], start_point[1], start_point[2]))

    LOGGER.timeit('_prody_calcChannels')

    # Edge-cost mode for the Dijkstra routing (buildSparseGraph). Default: the
    # mesh-invariant profile integral for the straight-edged homogenized/simple
    # diagrams, the legacy l/(d^2+b) for weighted (Apollonius edges are arcs,
    # where the straight-chord integral is only approximate; weighted is
    # experimental). An explicit value overrides the per-diagram default.
    if edge_cost is None:
        edge_cost = 'bottleneck' if diagram == 'weighted' else 'integral'
    elif edge_cost not in ('integral', 'bottleneck'):
        raise ValueError("edge_cost must be 'integral', 'bottleneck' or None, "
                         "got {0!r}".format(edge_cost))
    elif edge_cost == 'integral' and diagram == 'weighted':
        raise ValueError("edge_cost='integral' is only valid for the straight-edge "
                         "diagrams ('homogenized'/'simple'); the weighted "
                         "(Apollonius) diagram has arc edges the straight-chord "
                         "integral cannot price. Use edge_cost='bottleneck' (the "
                         "default for diagram='weighted') or None.")
    
    _reportAtomsInputComposition(atoms)
    atoms = atoms.select('not water') # water is excluded from the selection
    calculator = ChannelCalculator(atoms, r2=r2, sparsity=sparsity,
                                   route_tolerance=route_tolerance,
                                   edge_cost=edge_cost)

    elements = np.char.upper(np.asarray(atoms.getElements(), dtype=str))
    has_hydrogens = bool(np.any(elements == 'H'))

    # An experimental structure generally carries no hydrogens -- X-ray and cryo-EM
    # alike, since neither resolves them except at the very highest resolutions -- and
    # its carbons keep their full vdW radius, so the ~0.6 A the missing H occupied is
    # left as void, around every heavy atom at once, including buried contacts that
    # never come apart. That is usually harmless, and is often defended as standing in
    # for thermal motion: a probe of water size cannot enter those interstices anyway,
    # and protonated and unprotonated runs agree from about 1.2 A upwards. Below that
    # the probe is small enough to thread them and the interior percolates into a
    # sponge rather than merely widening. Those routes are fictitious, not the real
    # ones made wider. So a sub-water probe needs real hydrogens; above it, take the
    # file as it comes.
    if not has_hydrogens and r2 < 1.2:
        _warn("structure has no hydrogens and r2={0:.2f} is below 1.2 A: the space "
              "left by the missing H is then wide enough for the probe to pass, and "
              "channels will be found through interstices that do not exist in the "
              "real protein (their number can rise several-fold). Either add "
              "hydrogens, or raise r2 to 1.2 A or more, where protonated and "
              "unprotonated structures give the same channels.".format(r2))

    if diagram == "simple":
        # 'simple' builds an *unweighted* Delaunay of the atom centres, i.e. it
        # ignores the differences between atomic radii. That approximation is worst
        # when the radius spread is largest -- which is exactly when hydrogens (small
        # vdW) are present -- so warn there and steer the user to a radius-aware mode.
        # With H absent the heavy-atom radii are much closer, so 'simple' is more
        # defensible and matches the heavy-atom-only input most tools accept (at the
        # cost of over-large empty space where the missing H would sit).
        if has_hydrogens:
            _warn("diagram='simple' with hydrogens present: the unweighted "
                "Voronoi diagram ignores radius differences, which are largest when H "
                "are present, so its topology and clearances are significantly less "
                "accurate. Consider diagram='homogenized' (or 'weighted'), which "
                "account for per-atom radii.")

    coords = atoms.getCoords()
    vdw_radii = calculator.getVdwRadii(atoms.getElements())
    # Burial is a property of the protein, not of the tessellation, so the enclosure
    # test that strips the moat runs against the real atoms rather than the balls the
    # diagram happens to be built on. Homogenization would otherwise make it a
    # function of max_deviation, which min_enclosure must not be.
    atom_coords = coords
    # For diagram="weighted" only: a homogenized-surface depth oracle used to relabel
    # the additively-weighted diagram's leaky surface mouths (see getSurfaceCavities).
    mouth_oracle = None

    if diagram == "homogenized":
        LOGGER.timeit('_prody_channels_homogenize')
        coords, vdw_radii = calculator.homogenizeAtoms(coords, vdw_radii, max_deviation)
        LOGGER.report("Substituted {0} atoms with {1} homogeneous balls of radius {2:.2f} A in %.2fs.".format(
            atoms.numAtoms(), len(coords), float(vdw_radii[0])), '_prody_channels_homogenize')

    LOGGER.timeit('_prody_channels_tessellation')
    if diagram == "weighted":
        # True additively-weighted (Apollonius) Voronoi network via the third-party
        # vorpy package: van der Waals radii are baked into the diagram exactly,
        # instead of being approximated by homogenising atoms into uniform balls.
        # buildAwTessellation returns the same (simplices, neighbors, vertices) triple
        # a scipy Delaunay would, so the downstream erosion/cavity pipeline is
        # untouched. Because every AW vertex is equidistant (additively) to its 4
        # tangent atoms, the sum-based clearance test in deleteSimplices3d reduces
        # exactly to the per-atom clearance, so the returned clearances are not needed.
        if not checkAndImport('vorpy'):
            raise ImportError('diagram="weighted" requires the vorpy package for the '
                'additively-weighted (Apollonius) Voronoi diagram. Install vorpy, or '
                'use diagram="homogenized"/"simple".')
        # The compiled calc_vert kernel needs numba; without it the weighted path
        # still works but is ~5x slower, so fall back with a warning rather than fail.
        accelerate = checkAndImport('numba')
        if not accelerate:
            _warn('numba is not installed; the additively-weighted tessellation '
                'will run without the compiled kernel and may be very slow.')
        from ._vorpy_aw import buildAwTessellation, resolveCachePath
        try:
            title = atoms.getTitle()
        except Exception:
            title = None
        cache_path = resolveCachePath(weighted_cache, output_path, title)
        simplices, neighbors, verts, _ = buildAwTessellation(
            coords, vdw_radii, max_vert=max(2.0 * r1, 8), accelerate=accelerate,
            cache=cache_path)
        LOGGER.report('Additively-weighted (Apollonius) tessellation of {0} atoms '
            'constructed in %.2fs.'.format(len(coords)),
            '_prody_channels_tessellation')
        # The AW->simplicial mapping leaves false interior boundary faces that the
        # pipeline would misread as surface mouths (collapsing channels to stubs).
        # Build a homogenized diagram of the same atoms as an interior/exterior depth
        # oracle; getSurfaceCavities then keeps only exit tetrahedra within
        # weighted_mouth_depth Angstrom (geodesic) of the true molecular surface.
        if weighted_mouth_depth is not None:
            LOGGER.timeit('_prody_channels_mouth_oracle')
            mouth_oracle = calculator.buildSurfaceDepthOracle(
                coords, vdw_radii, r1, max_deviation, weighted_mouth_depth)
            LOGGER.report('Homogenized surface oracle (weighted mouth relabeling) '
                'built in %.2fs.', '_prody_channels_mouth_oracle')
    else:
        from scipy.spatial import Delaunay
        # We deliberately do NOT joggle/jitter the input (no QJ), unlike CAVER,
        # which perturbs by ~0.001 A to dodge the cospherical "T5" degeneracy it
        # reports in nearly every structure. scipy's default Qhull options
        # (Qbb Qc Qz) merge cospherical facets instead of joggling, so the
        # circumcenters stay finite even at degeneracies (verified: 0 NaN/inf
        # across millions of tetrahedra, even under heavy homogenized refinement).
        # Not joggling keeps the pipeline exactly reproducible (homogenizeAtoms is
        # a fixed Fibonacci lattice, and nothing here uses an RNG). The only
        # residual is coincident circumcenters at true degeneracies, handled where
        # it matters by the twin-tetrahedron guard in _edgeBottleneck.
        dela = Delaunay(coords)
        # circumcenters straight from the Delaunay paraboloid lifting, so we
        # skip the redundant second Qhull pass (scipy Voronoi). Numerically identical
        # to voro.vertices for points in general position.
        simplices = dela.simplices
        neighbors = dela.neighbors
        verts = calculator.calcCircumcenters(dela)
        LOGGER.report('Delaunay tessellation of {0} points constructed in %.2fs.'.format(
            len(coords)), '_prody_channels_tessellation')

    LOGGER.timeit('_prody_channels_surface')
    s_prt = State(simplices, neighbors, verts)
    
    if PY3K:
        s_tmp = State(*s_prt.getState())
        s_prv = State(None, None, None)
    else:
        s_tmp = apply(State, s_prt.getState())
        s_prv = State(None, None, None) 
        
    while True:
        s_prv.setState(*s_tmp.getState())
        
        if PY3K:
            s_tmp.setState(*calculator.deleteSimplices3d(coords, *(s_tmp.getState() + tuple([vdw_radii, r1, True]))))
        else:
            tmp_state = calculator.deleteSimplices3d(coords, *(s_tmp.getState() + [vdw_radii, r1, True]))
            s_tmp.setState(*tmp_state)

        if s_tmp == s_prv:
            break
        
    s_srf = State(*s_tmp.getState())

    # Moat removal: erode the r1 surface inward with the r2 probe, stripping the shell
    # of true exterior that a large r1 probe bridges over instead of entering (it would
    # otherwise join the cavity and offer wide, low-cost routes along the outside).
    # Erosion stops where the tetrahedra stop being open to the solvent, which is a
    # local criterion, so neither the mesh nor r1 sets how deep the peel goes.
    s_srf = State(*calculator.peelSurfaceByEnclosure(
        coords, *(s_srf.getState() + tuple([vdw_radii, r2, atom_coords,
                                            min_enclosure, max_peel_depth]))))

    s_inr = State(*calculator.deleteSimplices3d(coords, *(s_srf.getState() + tuple([vdw_radii, r2, False]))))

    l_first_layer_simp, l_second_layer_simp = calculator.surfaceLayer(s_srf.simp, s_inr.simp, s_srf.neigh)
    s_clr = State(*calculator.deleteSection(l_first_layer_simp, *s_inr.getState()))
    LOGGER.report('Surface and inner simplices filtered in %.2fs.', '_prody_channels_surface')

    LOGGER.timeit('_prody_channels_cavities')
    c_cavities = calculator.findGroups(s_clr.neigh)
    c_surface_cavities = calculator.getSurfaceCavities(c_cavities, s_clr.simp,
                                                       l_second_layer_simp,
                                                       s_clr, coords,
                                                       vdw_radii, sparsity,
                                                       mouth_oracle)

    calculator.findDeepestTetrahedra(c_surface_cavities, s_clr.neigh, s_clr.verti,
                                     coords, s_clr.simp)
    if start_point is not None:
        c_surface_cavities = calculator.setStartingTetrahedraFromPoint(
            c_surface_cavities, s_clr.verti, start_point, coords, vdw_radii,
            s_clr.simp, s_clr.neigh, restrict_channels_to_start_point,
            start_point_search)

    c_filtered_cavities = calculator.filterCavities(c_surface_cavities, min_depth)
    LOGGER.report('{0} surface cavities detected and filtered in %.2fs.'.format(
        len(c_filtered_cavities)), '_prody_channels_cavities')
    
    if cavities_only:
        if max_depth is not None:
            calculator.trimCavitiesByDepth(c_filtered_cavities, max_depth)

        if min_tetrahedra is not None or max_tetrahedra is not None:
            c_filtered_cavities = calculator.filterCavitiesByTetrahedra(
                c_filtered_cavities, min_tetrahedra, max_tetrahedra)

        calculator.calculate_cavity_volumes(c_filtered_cavities, s_clr.simp, coords)

        if min_volume is not None or max_volume is not None:
            c_filtered_cavities = calculator.filterCavitiesByVolume(
                c_filtered_cavities, min_volume, max_volume)
    
    merged_cavities = calculator.mergeCavities(c_filtered_cavities, s_clr.simp)

    # Early-return for the calcSurfaceCavities function:
    if cavities_only:
        LOGGER.info("Returning surface cavities")
        
        if output_path:
            output_path = Path(output_path)
            
            if output_path.is_dir():
                output_path = output_path / "surface_cavities.pqr"
            elif not (output_path.suffix == ".pdb" or output_path.suffix == ".pqr"):
                output_path = output_path.with_suffix(".pqr")

            if not separate:
                LOGGER.info("Saving surface cavities to " + str(output_path) + ".")
            else:
                LOGGER.info("Saving multiple surface cavities to directory " + str(output_path.parent) + ".")

            calculator.saveCavitiesToPdb(c_filtered_cavities, s_clr.verti, 
                                         output_path, separate)

        LOGGER.report('Surface cavity calculation completed in %.2fs.', '_prody_calcChannels')
        return c_filtered_cavities, [coords, s_srf.simp, merged_cavities, s_clr.simp, s_clr.verti]

    LOGGER.timeit('_prody_channels_pathfinding')
    # build the weighted adjacency matrix once for the whole cleared
    # state, then run a single multi-target Dijkstra per cavity (scipy csgraph),
    # instead of one heap Dijkstra per (seed, exit) pair.
    simplices, neighbors, vertices = s_clr.getState()
    graph = calculator.buildSparseGraph(simplices, neighbors, vertices, coords,
                                         vdw_radii)
    for cavity in c_filtered_cavities:
        calculator.dijkstra(cavity, graph, simplices, neighbors, vertices, 
                            coords, vdw_radii,
                            truncate_at_surface, similarity)
    LOGGER.report('Channel pathfinding (graph Dijkstra) over {0} cavities completed in %.2fs.'.format(
        len(c_filtered_cavities)), '_prody_channels_pathfinding')

    calculator.filterChannelsByBottleneck(c_filtered_cavities, bottleneck)
    
    if min_volume is not None or max_volume is not None:
        calculator.filterChannelsByVolume(c_filtered_cavities, min_volume, 
                                          max_volume)
    
    channels = [channel for cavity in c_filtered_cavities for channel in cavity.channels]
    # Order channels by ascending Dijkstra cost so that channel 0 is the best
    # tunnel (a short path through wide tetrahedra). This ordering drives both
    # the returned list and the channel numbering in the saved PQR/PDB files.
    channels.sort(key=lambda ch: ch.cost if ch.cost is not None else float('inf'))

    no_of_channels = len(channels)
    LOGGER.info("Detected " + str(no_of_channels) + " channels.")
        
    if output_path:
        output_path = Path(output_path)
        
        if output_path.is_dir():
            output_path = output_path / "output.pqr"
        
        elif not (output_path.suffix == ".pdb" or output_path.suffix == ".pqr"):
            output_path = output_path.with_suffix(".pqr")
    
        if not separate:
            LOGGER.info("Saving results to " + str(output_path) + ".")
        else:
            LOGGER.info("Saving multiple results to directory " + str(output_path.parent) + ".")
        calculator.saveChannelsToPdb(channels, output_path, separate)
    else:
        LOGGER.info("No output path given.")

    LOGGER.report('Channel calculation completed in %.2fs.', '_prody_calcChannels')

    # Additional information can be obtained
    if return_details:
        details = {'calculator': calculator, 
                    'simplices': s_clr.simp, 
                    'neighbors': s_clr.neigh, 
                    'vertices': s_clr.verti, 
                    'coords': coords, 
                    'vdw_radii': vdw_radii}
        
        return channels, [coords, s_srf.simp, merged_cavities, s_clr.simp], details

    return channels, [coords, s_srf.simp, merged_cavities, s_clr.simp]


def calcPoresFromChannels(channels, details):
    """Construct potential pores from previously identified channels using 
    :func:`calcChannels`. This function performs a post-processing analysis of 
    channels and requires ``return_details`` set to ``True`` in :func:`calcChannels`.
    
    The pore-construction procedure consists of the following steps:

    1. Group channels according to their starting tetrahedron.
    2. Generate all unique pairs of channels within each group.
    3. Identify the common initial segment and the last tetrahedron shared by
       each pair of channel paths.
    4. Join the non-overlapping parts of the two channels at their branching
       tetrahedron to obtain a surface-to-surface path.
    5. Reject paths containing loops or discontinuities between neighboring
       tetrahedra.
    6. Remove identical paths and paths differing only in direction.
    7. Recalculate the centerline spline, radius profile, length, bottleneck,
       and volume of each resulting pore using approach implemented for channels
       identification and visualization.
    
    :arg channels: A list of channel objects or a single channel object. Each 
        channel should have a `getSplines()` method that returns two 
        CubicSpline objects: one for the centerline and one for the radii.
    :type channels: list or single channel object
    
    :arg details: Additional calculation data returned by
        :func:`calcChannels` with ``return_details=True``. The dictionary must
        contain ``calculator``, ``simplices``, ``neighbors``, ``vertices``,
        ``coords``, and ``vdw_radii``.
    :type details: dict

    :returns: Potential pores constructed from compatible channel pairs.
    :rtype: list of Channel 
    
    Usage:
    channels, surface, details = calcChannels(protein, return_details=True)
    pores = calcPoresFromChannels(channels, details)    """
    
    calculator = details['calculator']
    simplices = details['simplices']
    vertices = details['vertices']
    coords = details['coords']
    vdw_radii = details['vdw_radii']
    neighbors = details['neighbors']
    
    pores = []
    pore_paths = []
    seen_paths = set()
    channel_groups = {}
    
    # Group channels using their starting tetrahedron
    for channel_index, channel in enumerate(channels):
        path = np.asarray(channel.tetrahedra, dtype=np.intp)
        # If channel is smaller than two tetrahedra (probably very rare)
        # Those channels should be excluded because they can not be connected with others
        if len(path) < 2:
            continue

        start_tetrahedron = int(path[0])
        channel_groups.setdefault(start_tetrahedron, []).append((channel_index, channel, path))
        
    from itertools import combinations
    # Generate all channel pairs within each group
    for start_tetrahedron, group in channel_groups.items():
        if len(group) < 2:
            continue

        for (channel1_index, channel1, path1), (channel2_index, channel2, path2) in combinations(group, 2):
            common_length = 0

            for tetrahedron1, tetrahedron2 in zip(path1, path2):
                if tetrahedron1 != tetrahedron2:
                    break
                common_length += 1

            if common_length == 0:
                continue

            # If we have for example: path1: start → A → B → C → mouth 1 and path2: start → A → B → D → mouth 2
            # it will create mouth 1 → C → B → D → mouth 2
            branch_index = common_length - 1
            pore_path = np.concatenate((path1[branch_index:][::-1], path2[branch_index + 1:]))
            
            # Reject paths containing loops
            if len(np.unique(pore_path)) != len(pore_path):
                continue

            # Continulity check of the pores (neighbours)
            is_continuous = True
            for tetrahedron1, tetrahedron2 in zip(pore_path[:-1], pore_path[1:]):
                if tetrahedron2 not in neighbors[tetrahedron1]:
                    is_continuous = False
                    break
            if not is_continuous:
                continue
                
            # Remove identical paths            
            path_key = tuple(int(tetrahedron) for tetrahedron in pore_path)
            canonical_key = min(path_key, path_key[::-1])

            if canonical_key in seen_paths:
                continue

            seen_paths.add(canonical_key)
            pore_paths.append(pore_path)
    
    # Pores reconstruction
    for pore_path in pore_paths:
        centerline_spline, radius_spline, length, bottleneck, volume = calculator.processChannel(
                                                pore_path, vertices, coords, vdw_radii, simplices)
        pore = Channel(pore_path, centerline_spline, radius_spline, length, bottleneck, volume, 0.0)
        pores.append(pore)

    return pores
            
                
def calcChannelsMultipleFrames(atoms, trajectory=None, output_path=None, 
    separate=False, start_point=None, **kwargs):
    """Compute channels for each frame in a given trajectory or multi-model 
    PDB file.

    This function calculates the channels for each frame in a trajectory or for
     each model in a multi-model PDB file. The `kwargs` can include parameters 
     necessary for channel calculation. If the `separate` parameter is set to 
     True, each detected channel will be saved in a separate PDB file.

    :arg atoms: Atomic data or object containing atomic coordinates and methods 
        for accessing them.
    :type atoms: object

    :arg trajectory: Trajectory object containing multiple frames or a 
        multi-model PDB file.
    :type trajectory: Atomic or Ensemble object

    :arg output_path: Optional path to save the resulting channels and 
        associated data in PDB format. If a directory is specified, each 
        frame/model will have its results saved in separate files. If None, 
        results are not saved. Default is None.
    :type output_path: str or None

    :arg separate: If True, each detected channel is saved to a separate PDB 
        file for each frame/model.
        If False, all channels for each frame/model are saved in a single file. 
        Default is False.
    :type separate: bool

    :arg start_point: Optional starting point for channel search. If provided, 
        the algorithm will use the tetrahedron whose Voronoi vertex is closest 
        to this point as the starting tetrahedron (overriding the default automatic 
        seed selection based on the deepest tetrahedron). Coordinates must be given in Å.
    :type start_point: list, tuple, or ndarray (length 3), or None 

    :arg kwargs: Additional parameters required for channel calculation. This can 
        include parameters such as radius values (r1, r2), minimum depth (min_depth), 
        bottleneck values, etc. 
        See the available parameters in calcChannels().
    :type kwargs: dict

    :returns: List of channels and surfaces computed for each frame or model. 
        Each entry in the list corresponds to a specific frame or model.
    :rtype: list of lists

    Example usage:
    channels_all, surfaces_all = calcChannelsMultipleFrames(atoms, trajectory=traj, 
                                    output_path="channels.pdb", separate=False, r1=3, 
                                    r2=0.9, min_depth=5, bottleneck=1, sparsity=3)
                                  
    channels_all, surfaces_all = calcChannelsMultipleFrames(atoms, trajectory=traj, 
                                    output_path="channels.pdb", separate=False, 
                                    start_point=[-10.353, -0.133, 5.608]) """

    
    if PY3K:
        if not checkAndImport('pathlib'):
            errorMsg = 'To run showChannels, please install open3d.'
            raise ImportError(errorMsg)
                
        from pathlib import Path
    else:
        if not checkAndImport('pathlib2'):
            errorMsg = 'To run showChannels, please install pathlib2 for Python 2.7.'
            raise ImportError(errorMsg)
        
        from pathlib2 import Path
        
    try:
        coords = getCoords(atoms)
    except AttributeError:
        try:
            checkCoords(coords)
        except TypeError:
            raise TypeError('coords must be an object with `getCoords` method')    

    channels_all = []
    surfaces_all = []
    start_frame = kwargs.pop('start_frame', 0)
    stop_frame = kwargs.pop('stop_frame', -1)
    
    if output_path:
        output_path = Path(output_path)
        if output_path.suffix == ".pqr":
            output_path = output_path.with_suffix('')

    if trajectory is not None:
        if isinstance(trajectory, Atomic):
            trajectory = Ensemble(trajectory)

        nfi = trajectory._nfi
        trajectory.reset()

        if stop_frame == -1:
            traj = trajectory[start_frame:]
        else:
            traj = trajectory[start_frame:stop_frame+1]

        atoms_copy = atoms.copy()
        for j0, frame0 in enumerate(traj, start=start_frame):
            LOGGER.info("Frame: {0}".format(j0))
            atoms_copy.setCoords(frame0.getCoords())
            if output_path:
                channels, surfaces = calcChannels(atoms_copy, str(output_path) + "{0}.pqr".format(j0), separate, start_point=start_point, **kwargs)
            else:
                channels, surfaces = calcChannels(atoms_copy, start_point=start_point, **kwargs)
            channels_all.append(channels)
            surfaces_all.append(surfaces)
        trajectory._nfi = nfi

    else:
        if atoms.numCoordsets() > 1:
            for i in range(len(atoms.getCoordsets()[start_frame:stop_frame])):
                LOGGER.info("Model: {0}".format(i+start_frame))
                atoms.setACSIndex(i+start_frame)
                if output_path:
                    channels, surfaces = calcChannels(atoms, str(output_path) + "{0}.pqr".format(i+start_frame), separate, start_point=start_point, **kwargs)
                else:
                    channels, surfaces = calcChannels(atoms, start_point=start_point, **kwargs)
                channels_all.append(channels)
                surfaces_all.append(surfaces)
        else:
            LOGGER.info("Include trajectory or use multi-model PDB file.")

    return channels_all, surfaces_all


def calcSurfaceCavitiesMultipleFrames(atoms, trajectory=None, output_path=None, 
    separate=False, **kwargs):
    """Compute surface cavities for each frame in a trajectory or multi-model PDB.

    This function calculates surface cavities for each frame of a trajectory or
    for each model of a multi-model PDB structure. For every frame/model, it
    calls :func:`calcSurfaceCavities` and stores the detected cavities together
    with the corresponding surface representation. The `kwargs` argument is
    passed directly to :func:`calcSurfaceCavities` and can include parameters
    controlling cavity detection, filtering, and output generation.

    :arg atoms: Atomic object containing the molecular structure. For trajectory 
        analysis, this object provides the reference topology and is updated with 
        coordinates from each frame. For multi-model PDB files, the individual 
        coordinate sets are analyzed one by one.
    :type atoms: :class:`.Atomic`

    :arg trajectory: Optional trajectory or ensemble object containing multiple
        coordinate frames. If provided, surface cavities are calculated for each
        selected trajectory frame. If not provided, the function attempts to use
        multiple coordinate sets stored in `atoms`.
    :type trajectory: :class:`.Atomic`, :class:`.Ensemble`, or trajectory-like object

    :arg output_path: Optional filename used to save detected surface cavities.
        If provided, one output file is generated for each frame/model by
        appending the frame/model index to the file name. If `None`, results are
        returned but not written in the folder. Default is `None`.
    :type output_path: str or None

    :arg separate: If `True`, each detected surface cavity is saved as a separate 
        PQR/PDB file for each frame/model. If `False`, all cavities detected 
        in a given frame/model are saved in a single file. Default is `False`.
    :type separate: bool

    :arg kwargs: Additional parameters passed to :func:`calcSurfaceCavities`.
        These can include `r1`, `r2`, `min_depth`, `max_depth`,
        `min_tetrahedra`, `max_tetrahedra`, `min_volume`, `max_volume`,
        `start_frame`, and `stop_frame`.
    :type kwargs: dict

    :returns: Two lists:
        - `cavities_all`: a list containing detected surface cavities for each
          analyzed frame/model,
        - `surfaces_all`: a list containing the corresponding surface representations 
          for each analyzed frame/model.
    :rtype: tuple (list, list)

    Example usage:
    protein = parsePDB('1tqn').select('protein')
    cavities_all, surfaces_all = calcSurfaceCavitiesMultipleFrames(protein, 
                                trajectory=traj, output_path="surface_cavities",
                                r1=4.5, r2=2.0, min_depth=1.5, max_depth=2.5, min_volume=50)

    cavities_all, surfaces_all = calcSurfaceCavitiesMultipleFrames(protein, start_frame=0, 
                                stop_frame=10, r1=4.5, r2=2.0) """

    if PY3K:
        if not checkAndImport('pathlib'):
            raise ImportError('To run calcSurfaceCavitiesMultipleFrames, please install pathlib.')
        from pathlib import Path
    else:
        if not checkAndImport('pathlib2'):
            raise ImportError('To run calcSurfaceCavitiesMultipleFrames, please install pathlib2 for Python 2.7.')
        from pathlib2 import Path

    try:
        coords = getCoords(atoms)
    except AttributeError:
        try:
            checkCoords(coords)
        except TypeError:
            raise TypeError('coords must be an object with `getCoords` method')

    cavities_all = []
    surfaces_all = []

    start_frame = kwargs.pop('start_frame', 0)
    stop_frame = kwargs.pop('stop_frame', -1)

    if output_path:
        output_path = Path(output_path)
        if output_path.suffix == ".pqr":
            output_path = output_path.with_suffix('')

    if trajectory is not None:
        if isinstance(trajectory, Atomic):
            trajectory = Ensemble(trajectory)

        nfi = trajectory._nfi
        trajectory.reset()
        if stop_frame == -1:
            traj = trajectory[start_frame:]
        else:
            traj = trajectory[start_frame:stop_frame + 1]

        atoms_copy = atoms.copy()
        for j0, frame0 in enumerate(traj, start=start_frame):
            LOGGER.info("Frame: {0}".format(j0))
            atoms_copy.setCoords(frame0.getCoords())

            if output_path:
                cavities, surface = calcSurfaceCavities(atoms_copy,
                    output_path=str(output_path) + "{0}.pqr".format(j0), 
                    separate=separate, **kwargs)
            else:
                cavities, surface = calcSurfaceCavities(atoms_copy, 
                                                        separate=separate, 
                                                        **kwargs)

            cavities_all.append(cavities)
            surfaces_all.append(surface)

        trajectory._nfi = nfi

    else:
        if atoms.numCoordsets() > 1:
            coordsets = atoms.getCoordsets()

            if stop_frame == -1:
                model_indices = range(start_frame, len(coordsets))
            else:
                model_indices = range(start_frame, stop_frame + 1)

            for i in model_indices:
                LOGGER.info("Model: {0}".format(i))
                atoms.setACSIndex(i)

                if output_path:
                    cavities, surface = calcSurfaceCavities(
                        atoms,
                        output_path=str(output_path) + "{0}.pqr".format(i),
                        separate=separate,
                        **kwargs)
                else:
                    cavities, surface = calcSurfaceCavities(
                        atoms,
                        separate=separate,
                        **kwargs)

                cavities_all.append(cavities)
                surfaces_all.append(surface)
        else:
            LOGGER.info("Include trajectory or use multi-model PDB file.")

    return cavities_all, surfaces_all


def parseParameters(channels, **kwargs):
    """Extracts and returns the lengths, bottlenecks, and volumes of each 
    channel in a given list of channels. """
    
    lengths = []
    bottlenecks = []
    volumes = []
    param_file_name = kwargs.pop('param_file_name', None)
    
    for nr_ch, channel in enumerate(channels):
        lengths.append(channel.length)
        bottlenecks.append(channel.bottleneck)
        volumes.append(channel.volume)
        
        if param_file_name is not None:
            with open(param_file_name+'_Parameters_All_channels.txt', "a") as f_par:
                f_par.write(("{0}_channel{1}: {2} {3} {4}\n".format(param_file_name, nr_ch, channel.length, channel.bottleneck, channel.volume)))
            
    return lengths, bottlenecks, volumes


def getChannelParameters(channels, **kwargs):
    """Extracts and returns the lengths, bottlenecks, and volumes of each 
    channel in a given list of channels.

    This function iterates through a list of channel objects, extracting the
    length, bottleneck, and volume of each channel. These values are collected
    into separate lists, which are returned as a tuple for further use.

    :arg channels: A list of channel objects, where each channel has attributes
      `length`, `bottleneck`,and `volume`. These attributes represent the 
      length of the channel, the minimum radius (bottleneck) along its path, 
      and the total volume of the channel, respectively.
    :type channels: list

    :arg param_file_name: The files with parameters will be saved in a text 
        file with the provided name. Use one word which will be added to
        '_Parameters_All_channels.txt' suffix. If further analysis will be
        performed with selectChannelBySelection() function, the preferable 
        param_file_name is PDB+chain for example: '1bbhA'.
    :type param_file_name: str 

    :returns: Three lists containing the lengths, bottlenecks, and volumes of 
        the channels.
    :rtype: tuple (list, list, list)

    Example usage:
    lengths, bottlenecks, volumes = getChannelParameters(channels) """
    
    multi_model_param = []
    param_file_name = kwargs.get('param_file_name', None)

    try:
        results_L_B_V = parseParameters(channels, **kwargs)
        lengths, bottlenecks, volumes = results_L_B_V
        LOGGER.info("Channel {0}: \t{1} \t{2} \t{3}".format('ID', 'Volume [Å³]',
                                                             'Length [Å]', 
                                                             'Bottleneck [Å]'))
        for i in range(len(lengths)):
            LOGGER.info("channel {0}: \t{1} \t\t{2} \t\t{3}".format(i, np.round(volumes[i],2), np.round(lengths[i], 2), np.round(bottlenecks[i], 2)))
        return results_L_B_V

    except:
        for nr_i,i in enumerate(channels):
            safe_param_file_name = param_file_name if param_file_name is not None else ""
            results = parseParameters(channels[nr_i], param_file_name=safe_param_file_name + str(nr_i))
            multi_model_param.append(results) 
            
        LOGGER.info("Channel {0}: \t{1} \t{2} \t{3}".format('ID', 'Volume [Å³]', 
                                                            'Length [Å]', 
                                                            'Bottleneck [Å]'))
        for frame_nr, frame in enumerate(multi_model_param):
            lengths, bottlenecks, volumes = frame
            LOGGER.info("Frame {0}".format(frame_nr))
            for i in range(len(lengths)):
                LOGGER.info("channel {0}: \t{1} \t\t{2} \t\t{3}".format(i, np.round(volumes[i],2), np.round(lengths[i], 2), np.round(bottlenecks[i], 2)))
        return multi_model_param


def getChannelParametersMultipleFrames(channels_all, **kwargs):
    """Extract channel parameters for multiple frames or models.

    This function is a multi-frame wrapper for :func:`getChannelParameters`.
    It extracts channel parameters for each model or trajectory frame separately.
    Each element of ``channels_all`` is treated as the list of channels calculated
    for one frame/model.

    This function should be used with channels returned by
    :func:`calcChannelsMultipleFrames`.

    :arg channels_all: list of channel lists returned by
        :func:`calcChannelsMultipleFrames`. Each element corresponds to one
        model or trajectory frame.
    :type channels_all: list

    :arg param_file_name: base name for the output parameter files. If provided,
        one file will be written for each model/frame with the frame/model index
        added to the file name.
    :type param_file_name: str

    :returns: A list of parameter tuples for each model/frame. Each tuple contains
        channel lengths, bottlenecks, and volumes.
    :rtype: list  """

    parameters_all = []
    for frame_nr, channels in enumerate(channels_all):
        LOGGER.info("Frame/model: {0}".format(frame_nr))
        params = getChannelParameters(channels, **kwargs)
        parameters_all.append(params)

    return parameters_all


def parseSurfaceCavityParameters(cavities, **kwargs):
    """Extract depths, volumes, and tetrahedra counts for surface cavities."""

    depths = []
    volumes = []
    tetrahedra_counts = []
    param_file_name = kwargs.pop('param_file_name', None)
    lines = []
    
    if param_file_name is not None:
        lines.append("# Cavity_id Volume [Å³] Depth [Å] Tetrahedra_count\n")
    
    for nr_cav, cavity in enumerate(cavities):
        depth = cavity.depth
        volume = cavity.volume
        tetrahedra_count = len(cavity.tetrahedra)
        depths.append(depth)
        volumes.append(volume)
        tetrahedra_counts.append(tetrahedra_count)

        if param_file_name is not None:
            lines.append("{0}_cavity{1}: {2:.3f} {3:.2f} {4}\n".format(
                param_file_name, nr_cav, volume, depth, tetrahedra_count))

    if param_file_name is not None:
        with open(param_file_name + '_Parameters_All_surface_cavities.txt', "w") as f_par:
            f_par.writelines(lines)
    
    return volumes, depths, tetrahedra_counts


def getSurfaceCavityParameters(cavities, **kwargs):
    """Extract volumes, depths, and tetrahedra counts of surface cavities.

    This function iterates through a list of surface cavity objects and extracts
    the volume, depth, and number of tetrahedra assigned to each cavity. These
    values are returned as separate lists and can optionally be saved to a text file.

    :arg cavities: A list of surface cavity objects returned by :func:`calcSurfaceCavities`.
    :type cavities: list

    :arg param_file_name: Optional name used to save cavity parameters to a text file. 
        The suffix '_Parameters_All_surface_cavities.txt' will be added.
    :type param_file_name: str

    :returns: Three lists containing volumes, depths, and tetrahedra counts.
    :rtype: tuple (list, list, list)

    Example usage:
    volumes, depths, tetrahedra_counts = getSurfaceCavityParameters(cavities)
    """

    multi_model_param = []
    param_file_name = kwargs.get('param_file_name', None)

    try:
        results_V_D_T = parseSurfaceCavityParameters(cavities, **kwargs)
        volumes, depths, tetrahedra_counts = results_V_D_T

        LOGGER.info("Cavity {0}: \t{1} \t{2} \t{3}".format('ID', 'Volume [Å³]', 
                                                           'Depth [Å]', 
                                                           'Tetrahedra count'))

        for i in range(len(volumes)):
            LOGGER.info("cavity {0}: \t{1} \t\t{2} \t\t{3}".format(i, np.round(volumes[i], 2), np.round(depths[i], 2),
                tetrahedra_counts[i]))

        return results_V_D_T

    except:
        for nr_i, i in enumerate(cavities):
            safe_param_file_name = param_file_name if param_file_name is not None else ""
            results = parseSurfaceCavityParameters(cavities[nr_i],
                param_file_name=safe_param_file_name + str(nr_i))
            multi_model_param.append(results)

        LOGGER.info("Cavity {0}: \t{1} \t{2} \t{3}".format('ID', 'Volume [Å³]', 
                                                           'Depth [Å]', 
                                                           'Tetrahedra count'))

        for frame_nr, frame in enumerate(multi_model_param):
            volumes, depths, tetrahedra_counts = frame
            LOGGER.info("Frame {0}".format(frame_nr))

            for i in range(len(volumes)):
                LOGGER.info("cavity {0}: \t{1} \t\t{2} \t\t{3}".format(i, np.round(volumes[i], 2), np.round(depths[i], 2),
                    tetrahedra_counts[i]))

        return multi_model_param


def getSurfaceCavityParametersMultipleFrames(cavities_all, **kwargs):
    """Provides surface cavity parameters for multiple frames or models.

    It analyzes surface cavities calculated for multi-model PDB files or
    trajectories and returns cavity parameters for each model/frame.

    :arg cavities_all: list of surface cavity lists returned by
        :func:`calcSurfaceCavitiesMultipleFrames`.
    :type cavities_all: list
    
    :arg param_file_name: base name for the output parameter files. If provided,
        one file will be written for each model/frame with the frame/model index
        added to the file name.
    :type param_file_name: str

    :returns: A list with surface cavity parameters for each frame/model.
    :rtype: list """

    parameters_all = []

    for i, cavities in enumerate(cavities_all):
        LOGGER.info("Model/frame: {0}".format(i))
        params = getSurfaceCavityParameters(cavities, **kwargs)
        parameters_all.append(params)

    return parameters_all


def getChannelAtoms(channels, protein=None, num_samples=5):
    """Generates an AtomGroup object representing the atoms along the paths of 
    the given channels and optionally combines them with an existing protein 
    structure.

    This function takes a list of channel objects and generates atomic 
    representations of the channels based on their centerline splines and 
    radius splines. The function samples points along each channel's centerline
     and assigns atom positions at these points with corresponding radii, 
    creating a list of PDB-formatted lines. These lines are then converted 
    into an AtomGroup object using the ProDy library. If a protein structure is
     provided, it is combined with the generated channel atoms by merging their
     respective PDB streams.

    :arg channels: A list of channel objects. Each channel has a method 
        `getSplines()` that
        returns the centerline spline and radius spline of the channel.
    :type channels: list

    :arg protein: An optional AtomGroup object representing a protein structure.
        If provided, it will be combined with the generated channel atoms.
    :type protein: prody.atomic.AtomGroup or None

    :arg num_samples: The number of atom samples to generate along each segment
         of the channel. More samples result in a finer representation of the
         channel. Default is 5.
    :type num_samples: int

    :returns: An AtomGroup object representing the combined atoms of the 
        channels and the protein, if a protein is provided.
    :rtype: prody.atomic.AtomGroup

    Example usage:
    atomic_structure = getChannelAtoms(channels, protein) """
    
    if PY3K:
        import io
    else:
        import StringIO as io
    
    from prody import parsePDBStream, writePDBStream

    def convert_lines_to_atomic(atom_lines):
        pdb_text = "\n".join(atom_lines)
        pdb_stream = io.StringIO(pdb_text)
        structure = parsePDBStream(pdb_stream)
        return structure

    atom_index = 1
    pdb_lines = []

    if not isinstance(channels, list):
        channels = [channels]
    
    for channel in channels:
        centerline_spline, radius_spline = channel.getSplines()
        samples = len(channel.tetrahedra) * num_samples
        t = np.linspace(centerline_spline.x[0], centerline_spline.x[-1], samples)
        centers = centerline_spline(t)
        radii = radius_spline(t)

        for i, (x, y, z, radius) in enumerate(zip(centers[:, 0], centers[:, 1], 
                                                  centers[:, 2], radii), 
                                                  start=atom_index):
            pdb_lines.append("ATOM  %5d  H   FIL T   1    %8.3f%8.3f%8.3f%6.2f%6.2f\n" % (i, x, y, z, 1.00, radius))
    
    if protein is not None:
        protein_stream = io.StringIO()
        writePDBStream(protein_stream, protein, csets=[protein.getACSIndex()])
        
        protein_stream.seek(0)

        protein_lines = protein_stream.readlines()
        if protein_lines[-1].strip() == 'END':
            protein_lines = protein_lines[:-1]
        
        combined_pdb_text = "".join(protein_lines) + "\n".join(pdb_lines) + "\nEND\n"
        combined_stream = io.StringIO(combined_pdb_text)
        combined_structure = parsePDBStream(combined_stream)

        return combined_structure

    channels_atomic = convert_lines_to_atomic(pdb_lines)
    return channels_atomic


def getChannelResidueNames(atoms, channels, **kwargs):
    '''Provides the resnames and resid of residues that are forming the channel(s). 
    Residues are extracted based on distA which is the distance between FIL atoms 
    (channel atoms) and protein residues.
    Results could be save as txt file by providing the `residues_file_name` parameter.
    
    :arg atoms: an Atomic object from which residues are selected 
    :type atoms: :class:`.Atomic`, :class:`.LigandInteractionsTrajectory`

    :arg channels: A list of channel objects. Each channel has a method 
        `getSplines()` that returns the centerline spline and radius spline of 
        the channel.
    :type channels: list

    :arg distA: Residues will be provided based on this value.
        default is 4 [Ang]
    :type distA: int, float 
    
    :arg residues_file_name: The file with residues will be saved in a text 
        file with the provided name. Use one word which will be added to 
        '_Residues_All_channels.txt' sufix. If further analysis will be 
        performed with selectChannelBySelection() function, the preferable 
        residues_file_name is PDB+chain for example: '1bbhA'.
    :type residues_file_name: str  
    
    :arg one_letter_aa: Whether to apply 1-latter code to residue name
        by defult is False
    :type one_letter_aa: bool  '''

    try:
        coords = (atoms._getCoords() if hasattr(atoms, '_getCoords') else
                    atoms.getCoords())
    except AttributeError:
        try:
            checkCoords(coords)
        except TypeError:
            raise TypeError('coords must be an object '
                            'with `getCoords` method')

    distA = kwargs.pop('distA', 4)
    residues_file_name = kwargs.pop('residues_file_name', None) 
    
    one_letter_aa = kwargs.pop('one_letter_aa', False)
    if one_letter_aa == True:
        from prody.atomic.atomic import AAMAP    
    
    if isinstance(channels, list):
        # Multiple channels
        selected_residues_ch = []
    
        for i, channel in enumerate(channels):
            atoms_protein = getChannelAtoms(channel, atoms)
            residues = atoms_protein.select('same residue as exwithin '+str(distA)+' of resname FIL')
    
            if residues is not None:
                resnames = residues.select('name CA').getResnames()
                if one_letter_aa == True:
                    resnames_1letter = [AAMAP["HIS"] if aa in ("HSD", "HSP") else AAMAP[aa] for aa in resnames]
                    resnames = resnames_1letter
                                    
                resnums = residues.select('name CA').getResnums()
                residues_info = ["{}{}".format(resname, resnum) for resname, resnum in zip(resnames, resnums)]
                residues_list = ", ".join(residues_info)
                residues_list = 'channel'+str(i)+': '+residues_list
                selected_residues_ch.append(residues_list)
            else:
                residues_list = "None"
            
    else:
        # Single channel analysis in case someone provide channels[0]
        atoms_protein = getChannelAtoms(channels, atoms)
        residues = atoms_protein.select('same residue as exwithin '+str(distA)+' of resname FIL')
        selected_residues_ch = []
        
        if residues is not None:
            resnames = residues.select('name CA').getResnames()
            if one_letter_aa == True:
                resnames_1letter = [AAMAP["HIS"] if aa in ("HSD", "HSP") else AAMAP[aa] for aa in resnames]
                resnames = resnames_1letter

            resnums = residues.select('name CA').getResnums()
            residues_info = ["{}{}".format(resname, resnum) for resname, resnum in zip(resnames, resnums)]
            residues_list = ", ".join(residues_info)
            selected_residues_ch.append(residues_list)
        else:
            residues_list = "None"

    if residues_file_name is not None:
        output_file = residues_file_name + '_Residues_All_channels.txt'
        with open(output_file, "a") as f_res:
            for k in selected_residues_ch:
                f_res.write(("{0}_{1}\n".format(residues_file_name, k)))
        
        LOGGER.info("Channel residues were saved to: {0}".format(output_file))
                
    return selected_residues_ch


def getChannelResidueNamesMultipleFrames(atoms, channels_all, trajectory=None, **kwargs):
    """Provides residue names for channels calculated for multiple frames/models.

    This function is a multi-frame wrapper for :func:`getChannelResidueNames`.
    For each model/frame, the atomic coordinates are matched with the
    corresponding channel prediction.

    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`

    :arg channels_all: list of channel lists returned by :func:`calcChannelsMultipleFrames`.
    :type channels_all: list

    :arg trajectory: optional trajectory object. If provided, coordinates are
        taken from trajectory frames. If None, a multi-model PDB is assumed and
        models are selected using ``setACSIndex``.
    :type trajectory: :class:`.Trajectory` or None

    :arg start_frame: first frame/model index. Default is 0.
    :type start_frame: int

    :arg stop_frame: last frame/model index. Default is -1, meaning all available
        frames/models in ``channels_all``.
    :type stop_frame: int

    :arg residues_file_name: base name for output residue files. If provided,
        one file will be written for each frame/model.
    :type residues_file_name: str

    :arg distA: maximal distance between channel FIL atoms and protein residues.
        Default is 4 Å.
    :type distA: int, float

    :arg one_letter_aa: whether to apply one-letter code to residue names.
        Default is False.
    :type one_letter_aa: bool

    :returns: A list of residue-name lists for each frame/model.
    :rtype: list  """

    start_frame = kwargs.pop('start_frame', 0)
    stop_frame = kwargs.pop('stop_frame', -1)
    residues_file_name = kwargs.pop('residues_file_name', None)
    selected_residues_all = []

    if trajectory is None:
        # multi-model PDB
        for frame_pos, channels in enumerate(channels_all):
            model_index = start_frame + frame_pos

            if stop_frame != -1 and model_index > stop_frame:
                break

            LOGGER.info("Model: {0}".format(model_index))
            atoms.setACSIndex(model_index)

            if residues_file_name is not None:
                frame_residues_file_name = residues_file_name + "_model{}".format(model_index)
            else:
                frame_residues_file_name = None

            residues = getChannelResidueNames(atoms, channels,
                residues_file_name=frame_residues_file_name, **kwargs)

            selected_residues_all.append(residues)

    else:
        # trajectory / DCD
        nfi = getattr(trajectory, '_nfi', None)

        if hasattr(trajectory, 'reset'):
            trajectory.reset()

        if stop_frame == -1:
            traj = trajectory[start_frame:]
        else:
            traj = trajectory[start_frame:stop_frame + 1]

        atoms_copy = atoms.copy()
        for frame_pos, frame in enumerate(traj):
            frame_index = start_frame + frame_pos

            if frame_pos >= len(channels_all):
                break

            LOGGER.info("Frame: {0}".format(frame_index))
            atoms_copy.setCoords(frame.getCoords())

            if residues_file_name is not None:
                frame_residues_file_name = residues_file_name + "_frame{}".format(frame_index)
            else:
                frame_residues_file_name = None

            residues = getChannelResidueNames(atoms_copy, channels_all[frame_pos],
                residues_file_name=frame_residues_file_name, **kwargs)

            selected_residues_all.append(residues)

        if nfi is not None:
            trajectory._nfi = nfi

    return selected_residues_all


def getSurfaceCavityResidueNames(atoms, cavities, surface, **kwargs):
    '''Provides the resnames and resid of residues that form surface cavities.

    Residues are extracted based on distA, which is the distance between surface
    cavity points and protein residues. Surface cavity points are taken from
    Voronoi vertices assigned to each cavity. Results can be saved as a txt file 
    by providing the `residues_file_name` parameter.

    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`, :class:`.LigandInteractionsTrajectory`

    :arg cavities: A list of surface cavity objects returned by :func:`calcSurfaceCavities`.
    :type cavities: list

    :arg surface: Surface data returned by :func:`calcSurfaceCavities`.
        The function uses `surface[4]`, which contains Voronoi vertices assigned
        to surface cavities.
    :type surface: list

    :arg distA: Residues will be provided based on this value.
        Default is 4 [Ang].
    :type distA: int, float

    :arg residues_file_name: The file with residues will be saved in a text file
        with the provided name. The suffix '_Residues_All_surface_cavities.txt' 
        will be added.
    :type residues_file_name: str

    :arg one_letter_aa: Whether to apply one-letter code to residue names.
        Default is False.
    :type one_letter_aa: bool

    :returns: A list of residue names and residue numbers for each surface cavity.
    :rtype: list
    '''

    try:
        coords = (atoms._getCoords() if hasattr(atoms, '_getCoords') else
                  atoms.getCoords())
    except AttributeError:
        try:
            checkCoords(coords)
        except TypeError:
            raise TypeError('coords must be an object '
                            'with `getCoords` method')

    distA = kwargs.pop('distA', 4)
    residues_file_name = kwargs.pop('residues_file_name', None)

    one_letter_aa = kwargs.pop('one_letter_aa', False)
    if one_letter_aa == True:
        from prody.atomic.atomic import AAMAP

    if surface is None or len(surface) < 5:
        raise ValueError('surface must contain Voronoi vertices in surface[4]')

    vertices = surface[4]
    if not isinstance(cavities, list):
        cavities = [cavities]

    selected_residues_cav = []

    for i, cavity in enumerate(cavities):
        if cavity.tetrahedra is None or len(cavity.tetrahedra) == 0:
            selected_residues_cav.append('cavity' + str(i) + ': None')
            continue

        points = vertices[cavity.tetrahedra]
        residues = atoms.select('same residue as exwithin ' + str(distA) + ' of center', center=points)

        if residues is not None:
            ca_residues = residues.select('name CA')

            if ca_residues is not None:
                resnames = ca_residues.getResnames()

                if one_letter_aa == True:
                    resnames_1letter = [AAMAP["HIS"] if aa in ("HSD", "HSE", "HSP", "HID", "HIE", "HIP")
                        else AAMAP.get(aa, aa) for aa in resnames]
                    resnames = resnames_1letter

                resnums = ca_residues.getResnums()
                residues_info = ["{}{}".format(resname, resnum)
                    for resname, resnum in zip(resnames, resnums)]

                residues_list = 'cavity' + str(i) + ': ' + ", ".join(residues_info)
                selected_residues_cav.append(residues_list)

            else:
                selected_residues_cav.append('cavity' + str(i) + ': None')
        else:
            selected_residues_cav.append('cavity' + str(i) + ': None')

    if residues_file_name is not None:
        output_file = residues_file_name + '_Residues_All_surface_cavities.txt'
        with open(output_file, "w") as f_res:
            f_res.write("# cavity_id residues_within_" + str(distA) + "_A\n")
            for k in selected_residues_cav:
                f_res.write("{0}_{1}\n".format(residues_file_name, k))
                
        LOGGER.info("Surface cavity residues were saved to: {0}".format(output_file))

    return selected_residues_cav


def getSurfaceCavityResidueNamesMultipleFrames(atoms, cavities_all, 
                                               surfaces_all, 
                                               trajectory=None, **kwargs):
    """Provides residue names for surface cavities calculated for multiple 
    frames/models.

    This function is a multi-frame wrapper for :func:`getSurfaceCavityResidueNames`. 
    For each model or trajectory frame, the atomic coordinates are matched with
      the corresponding surface cavity prediction. Thus, cavities calculated 
    for frame/model ``i`` are analyzed against the protein coordinates from 
    frame/model ``i``.

    This function should be used with results returned by 
    :func:`calcSurfaceCavitiesMultipleFrames`.

    :arg atoms: an Atomic object from which residues are selected.
    :type atoms: :class:`.Atomic`

    :arg cavities_all: list of surface cavity lists returned by
        :func:`calcSurfaceCavitiesMultipleFrames`. Each element corresponds
        to one model or trajectory frame.
    :type cavities_all: list

    :arg surfaces_all: list of surface data objects returned by
        :func:`calcSurfaceCavitiesMultipleFrames`. Each element corresponds
        to one model or trajectory frame and must contain Voronoi vertices in
        ``surface[4]``.
    :type surfaces_all: list

    :arg trajectory: optional trajectory object. If provided, coordinates are
        taken from trajectory frames. If None, a multi-model PDB is assumed.
    :type trajectory: :class:`.Trajectory` or None

    :arg start_frame: first frame/model index to analyze. Default is 0.
    :type start_frame: int

    :arg stop_frame: last frame/model index to analyze. Default is -1, meaning
        all available frames/models in ``cavities_all`` and ``surfaces_all``.
    :type stop_frame: int

    :arg residues_file_name: base name for output residue files. If provided,
        one file will be written for each model/frame with ``_modelX`` or
        ``_frameX`` added to the file name.
    :type residues_file_name: str

    :arg distA: maximal distance between surface cavity points and protein
        residues. Default is 4 Å.
    :type distA: int, float

    :arg one_letter_aa: whether to apply one-letter code to residue names.
        Default is False.
    :type one_letter_aa: bool  """

    start_frame = kwargs.pop('start_frame', 0)
    residues_file_name = kwargs.pop('residues_file_name', None)

    selected_residues_all = []
    
    if trajectory is None:
        # multi-model PDB
        for frame_pos, (cavities, surface) in enumerate(zip(cavities_all, 
                                                            surfaces_all)):
            model_index = start_frame + frame_pos
            atoms.setACSIndex(model_index)

            if residues_file_name is not None:
                frame_residues_file_name = residues_file_name + "_model{}".format(model_index)
            else:
                frame_residues_file_name = None

            residues = getSurfaceCavityResidueNames(atoms, cavities, surface,
                residues_file_name=frame_residues_file_name, **kwargs)

            selected_residues_all.append(residues)

    else:
        # trajectory / DCD
        if hasattr(trajectory, 'reset'):
            trajectory.reset()

        atoms_copy = atoms.copy()
        for frame_pos, frame in enumerate(trajectory):
            frame_index = start_frame + frame_pos
            atoms_copy.setCoords(frame.getCoords())

            if residues_file_name is not None:
                frame_residues_file_name = residues_file_name + "_frame{}".format(frame_index)
            else:
                frame_residues_file_name = None

            residues = getSurfaceCavityResidueNames(atoms_copy, cavities_all[frame_pos], surfaces_all[frame_pos],
                residues_file_name=frame_residues_file_name, **kwargs)

            selected_residues_all.append(residues)

    return selected_residues_all


def selectChannelBySelection(atoms, residue_sele, **kwargs):
    """Select PQR files with channels that are having FIL residues within 
    certain distance (distA) from selected residue (temporarily one residue).
    If not all files should be included use pqr_files to provide the new list. 
    For example:
    pqr_files = [file for file in os.listdir('.') if file.startswith('7lafA_') and file.endswith('.pqr')]
    pqr_files = [file for file in os.listdir('.') if '5kbd' in file and file.endswith('.pqr')]

    :arg atoms: an Atomic object from which residues are selected 
    :type atoms: :class:`.Atomic`, :class:`.LigandInteractionsTrajectory`

    :arg residue_sele: selection string
                        for example: 'resid 377 and chain A', 'resid 10 to 20'
    :type residue_sele: str
   
    :arg pqr_files: list of PQR files to analyze
                    default is False, which means that all .pqr files from the 
                    current directory will be analyzed.
    :type pqr_files: bool or list
    
    :arg folder_name: The name of the folder to which PDBs will be extracted
    :type folder_name: str

    :arg distA: non-zero value, maximal distance from selected region to 
        channel (FIL atoms). Default is 5.
    :type distA: int, float 
        
    :arg residues_file: File with residues forming the channel created by 
        getChannelResidues(), default is False 
    :type residues_file: bool

    :arg param_file: File with residues forming the channel created by 
        getChannelParameters(). Default is False.
    :type param_file: bool  """

    try:
        coords = (atoms._getCoords() if hasattr(atoms, '_getCoords') else
                    atoms.getCoords())
    except AttributeError:
        try:
            checkCoords(coords)
        except TypeError:
            raise TypeError('coords must be an object '
                            'with `getCoords` method')
    
    import os, shutil
    
    pqr_files = kwargs.pop('pqr_files', False)
    distA = kwargs.pop('distA', 5)
    folder_name = kwargs.pop('folder_name', 'selected_files')
    residues_file = kwargs.pop('residues_file', False)
    param_file = kwargs.pop('param_file', False)

    object_name = kwargs.pop('object_name', 'channel')
    residues_suffix = kwargs.pop('residues_suffix', '_Residues_All_channels.txt')
    parameters_suffix = kwargs.pop('parameters_suffix', '_Parameters_All_channels.txt')
    selected_residues_output = kwargs.pop('selected_residues_output', 
                                          'Selected_channel_residues.txt')
    selected_parameters_output = kwargs.pop('selected_parameters_output', 
                                            'Selected_channel_parameters.txt')

    copied_files_list = []
    
    if pqr_files == False:
        # take all PDBs from the current dir
        pqr_files = [file for file in os.listdir('.') if file.endswith('.pqr')]

    residue_sele = atoms.select(residue_sele)
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    for i in pqr_files:
        channel = parsePQR(i)
        if 'FIL' in np.unique(channel.getResnames()):
            sele_FIL = channel.select('same residue as exwithin '+str(distA)+' of center', center=residue_sele.getCoords())
        
            if sele_FIL is not None:
                shutil.copy(i, folder_name)
                LOGGER.info("Filtered files are now in: {0}".format(folder_name))
                copied_files_list.append(i)
            else:
                pass 

    # Extract parameters and/or residues with channel selection
    if residues_file == True:
        selected_residues = []
        for file in copied_files_list:
            try:
                PDB_id = file[:-4].split('_' + object_name)[0]
                selected_name = file[:-4].split('_')[-1]
                f = open(PDB_id + residues_suffix, 'r').readlines()

                for line in f:
                    if line.startswith(PDB_id + '_' + selected_name + ':'):
                        selected_residues.append(line)
            except:
                LOGGER.info("File {0} was not analyzed due to the lack of file or multiple channel file.".format(file))
                pass

        with open('Selected_channel_residues.txt', 'w') as f_out:
            f_out.writelines(selected_residues)

    if param_file == True:
        selected_param = []
        for file in copied_files_list:
            try:
                PDB_id = file[:-4].split('_channel')[0] 
                channel_name = file[:-4].split('_')[-1]
                f = open(PDB_id+'_Parameters_All_channels.txt', 'r').readlines()
                for line in f:
                    if line.startswith(PDB_id+'_'+channel_name+':'):
                        selected_param.append(line)
            except:
                LOGGER.info("File {0} was not analyzed due to the lack of file or multiple channel file.".format(file))
                pass

        with open('Selected_channel_parameters.txt', 'w') as f_out:
            f_out.writelines(selected_param)

    LOGGER.info("Selected files: ")
    LOGGER.info(' '.join(copied_files_list))
    LOGGER.info("If newly created files are empty please check whether the parameter names are: PDB_id+_Parameters_All_channels.txt")


def selectSurfaceCavityBySelection(atoms, residue_sele, **kwargs):
    """Select PQR files with surface cavities located close to a selected region.

    This function is a surface-cavity wrapper for :func:`selectChannelBySelection`.
    It selects PQR files containing surface cavities represented by FIL pseudoatoms
    that are located within a user-defined distance from a selected protein region.

    Surface cavity files should be generated by :func:`calcSurfaceCavities` or
    :func:`calcSurfaceCavitiesMultipleFrames`, preferably with ``separate=True``,
    so that each cavity is saved in an individual PQR file.

    :arg atoms: an Atomic object from which the reference region is selected
    :type atoms: :class:`.Atomic`, :class:`.LigandInteractionsTrajectory`

    :arg residue_sele: selection string defining the reference region.
        For example: ``'resid 173'``, ``'resid 173 and chain A'``, ``'resid 170 to 180'``.
    :type residue_sele: str

    :arg pqr_files: list of PQR files to analyze. If not provided, all PQR files
        from the current directory will be analyzed.
    :type pqr_files: bool or list

    :arg folder_name: name of the folder to which selected PQR files will be copied.
        Default is ``'selected_surface_cavities'``.
    :type folder_name: str

    :arg distA: maximal distance between the selected region and surface cavity
        FIL atoms. Default is 5 Å.
    :type distA: int, float

    :arg residues_file: if True, residue information for selected surface cavities
        will be extracted from ``*_Residues_All_surface_cavities.txt`` and saved
        to ``Selected_surface_cavity_residues.txt``.
        Default is False.
    :type residues_file: bool

    :arg param_file: if True, parameter information for selected surface cavities
        will be extracted from ``*_Parameters_All_surface_cavities.txt`` and saved
        to ``Selected_surface_cavity_parameters.txt``.
        Default is False.
    :type param_file: bool """

    kwargs.setdefault('object_name', 'cavity')
    kwargs.setdefault('residues_suffix', '_Residues_All_surface_cavities.txt')
    kwargs.setdefault('parameters_suffix', '_Parameters_All_surface_cavities.txt')
    kwargs.setdefault('selected_residues_output', 
                      'Selected_surface_cavity_residues.txt')
    kwargs.setdefault('selected_parameters_output', 
                      'Selected_surface_cavity_parameters.txt')

    return selectChannelBySelection(atoms, residue_sele, **kwargs)


def calcChannelSurfaceOverlaps(**kwargs):
    """Calculate overlapping parts of the predicted channels, tunnels, and 
    pores denote as 'FIL' atoms. Results are normalized within [0,1].

    :arg resolution: Surface sampling resolution.
        default is 0.5
    :type resolution: float
    
    :arg max_proc: Maximum number of parallel processes used to voxelize individual
        PQR files. If 1, files are processed serially. If None, all available CPU
        cores are used. Default is 2.
    :type max_proc: int or None

    :arg output_file_name: The name of the PDB file with overlapping surfaces.
    :type output_file_name: str

    :arg pqr_files: File with residues forming the channel created by
        getChannelResidues(). Default is False, then all the files from the
        current directory will be analyzed. When providing a list, only the
        PDBs from the list will be analyzed. When providing str, it will be
        treated as a folder path.
    :type pqr_files: bool, list or str
    
    Example usage:
    calcChannelSurfaceOverlaps() - all the files in the current directory will 
    be analyzed
    
    from pathlib import Path
    pqr_files = [str(f) for f in Path(".").glob("channels_*.pqr")]
    calcChannelSurfaceOverlaps(pqr_files=pqr_files, 
                output_file_name='results.pdb', max_proc=4)
    - files with the "channels_" prefix will be selected from the current folder 
    and analyzed using four parallel processes.
    
    calcChannelSurfaceOverlaps(pqr_files='./DATA', output_file_name='results.pdb') 
    - only files from the DATA folder will be analyzed and results will be saved 
    as results.pdb
    
    list_of_files = ['file1.pqr', 'file2.pqr', 'file3.pqr', ..]
    calcChannelSurfaceOverlaps(pqr_files=list_of_files, output_file_name='results.pdb') 
    - files from the list will be analyzed and results will be saved as results.pdb
    """
    
    import os
    import multiprocessing
    from collections import Counter

    resolution = kwargs.pop('resolution', 0.5)
    max_proc = kwargs.pop('max_proc', 2)
     
    pqr_files = kwargs.pop('pqr_files', False)
    if pqr_files == False or pqr_files is None:
        # take all PQRs from the current dir
        pqr_files = [file for file in os.listdir('.') if file.endswith('.pqr')]
    elif isinstance(pqr_files, str):
        # folder path
        if os.path.isdir(pqr_files):
            pqr_files = [os.path.join(pqr_files, file) for file in os.listdir(pqr_files) if file.endswith('.pqr')]
    elif isinstance(pqr_files, list):
        # list of PQRs
        pqr_files = [file for file in pqr_files if file.endswith('.pqr')]
    else:
        raise ValueError('Please provide list with PQR files, folder path, or nothing to analyze PQRs in the current folder')

    output_file_name = kwargs.pop('output_file_name','overlap_regions.pdb')
    
    if len(pqr_files) == 0:
        LOGGER.info("No PQR files found.")
        return None

    if os.path.exists(output_file_name):
        os.rename(output_file_name, output_file_name + '-old')

    if max_proc is None:
        max_proc = multiprocessing.cpu_count()

    max_proc = int(max_proc)
    if max_proc < 1:
        max_proc = 1

    max_proc = min(max_proc, len(pqr_files))

    LOGGER.info("Number of PQR files: {0}".format(len(pqr_files)))
    LOGGER.info("Resolution: {0}".format(resolution))
    LOGGER.info("max_proc: {0}".format(max_proc))

    merged_surface = Counter()
    tasks = [(pqr_file, resolution) for pqr_file in pqr_files]

    if max_proc > 1:
        LOGGER.info("Calculating overlaps using {0} processes.".format(max_proc))
        chunksize = max(1, len(tasks) // (max_proc * 4))
        with multiprocessing.Pool(processes=max_proc) as pool:
            for surface in pool.imap_unordered(_surfaceFromPqrWorker, tasks,
                                               chunksize=chunksize):
                merged_surface.update(surface)

    else:
        for pqr_file in pqr_files:
            LOGGER.info("Processing file: {0}".format(pqr_file))
            surface = _surfaceFromPqrWorker((pqr_file, resolution))
            merged_surface.update(surface)

    with open(output_file_name, 'w') as out:
        atom_id = 1

        for (ix, iy, iz), count in merged_surface.items():
            x = ix * resolution
            y = iy * resolution
            z = iz * resolution

            norm_count = float(count) / float(len(pqr_files))

            out.write("ATOM  {:5d}  H   FIL T   1    {:8.3f}{:8.3f}{:8.3f}{:6.2f}  1.00\n"
                .format(atom_id, x, y, z, norm_count))

            atom_id += 1

    LOGGER.info("Overlap written to: {0}".format(output_file_name))
    LOGGER.info("Number of occupied overlap voxels: {0}".format(len(merged_surface)))

    return output_file_name
    

def calcSurfaceCavityOverlaps(**kwargs):
    """Calculate overlapping regions of surface cavities represented as FIL atoms.

    It calculates spatial overlap between surface cavities saved as PQR files with
    FIL pseudoatoms, as generated by :func:`calcSurfaceCavities` or
    :func:`calcSurfaceCavitiesMultipleFrames`.

    Results are normalized within [0, 1], where the value corresponds to the
    fraction of analyzed PQR files contributing to a given spatial region.

    :arg resolution: surface sampling resolution. Default is 0.5.
    :type resolution: float
    
    :arg max_proc: Maximum number of parallel processes used to voxelize individual
        PQR files. If 1, files are processed serially. If None, all available CPU
        cores are used. Default is 2.
    :type max_proc: int or None

    :arg output_file_name: name of the output PDB file with overlapping cavity
        regions. Default is ``'surface_cavity_overlap_regions.pdb'``.
    :type output_file_name: str

    :arg pqr_files: PQR files with surface cavities represented as FIL atoms.
        If not provided, all PQR files from the current directory will be analyzed.
        A list of PQR files can also be provided.
    :type pqr_files: bool or list """

    kwargs.setdefault('output_file_name', 'surface_cavity_overlap_regions.pdb')

    return calcChannelSurfaceOverlaps(**kwargs)


def calcSurfaceCavities(atoms, output_path=None, r1=4.5, r2=2.0, min_depth=1.5,
                        max_depth=2.5, min_tetrahedra=None, max_tetrahedra=None,
                        min_volume=50, max_volume=None, sparsity=None,
                        separate=False):
    """Calculate surface cavities (pockets) on protein surface using CaviTracer 
    approach.

    :arg atoms: An object representing the molecular structure, typically 
        containing atomic coordinates and element types.
    :type atoms: `Atoms` object

    :arg output_path: Optional path to save the resulting cavities and 
        associated data in PQR (or PDB) format. If None, results are not saved.
         Default is None.
    :type output_path: str or None

    :arg separate: If True, each detected cavity is saved to a separate PQR 
        file. If False, all cavities are saved in a single PQR file. Default is 
        False.
    :type separate: bool

    :arg r1: The first radius threshold used during the deletion of simplices, 
        which is used to define the outer surface of the cavities. Default is 4.5.
    :type r1: float

    :arg r2: The second radius threshold used to define the inner surface of 
        the cavities. Default is 2.
    :type r2: float

    :arg min_depth: The minimum depth, in Angstrom, a cavity must reach to be
        considered. Depth is the geodesic distance from the surface opening along 
        the Voronoi network, a physical length independent of tessellation density.
        Default is 1.5.
    :type min_depth: float

    :arg max_depth: Maximum cavity depth, in Angstrom. Portions of a cavity deeper
        than this value are trimmed away, keeping the shallow surface shell that
        defines a pocket. Default is 2.5.
    :type max_depth: float

    :arg sparsity: Deprecated and ignored; accepted only so that existing calls
        keep working. It never affected surface cavities: it thinned the sampling
        of the mouth (exit) tetrahedra used as termini by the *channel* search,
        and no cavity property reads that thinned set. Cavity extent, depth,
        volume and filtering are all derived from the unthinned exit tetrahedra,
        so passing 1 or 15 returns the same cavities.
    :type sparsity: int

    :arg min_tetrahedra: Minimum number of tetrahedra required for a cavity to
        be retained. Smaller cavities are discarded. Default is None.
    :type min_tetrahedra: int

    :arg max_tetrahedra: Maximum number of tetrahedra allowed for a cavity to 
        be retained. Larger cavities are discarded. Default is None.
    :type max_tetrahedra: int

    :arg min_volume: Minimum volume required for a channel/cavity to be 
        retained. Default is 50.
    :type min_volume: float

    :arg max_volume: Maximum volume allowed for a channel/cavity to be 
        retained. Default is None.
    :type max_volume: float

    :returns: A tuple containing two elements:
        - `cavities`: A list of detected cavities, where each channel is an 
            object containing information about its path and geometry.
        - `surface`: A list containing additional information for further 
            visualization, including the atomic coordinates, simplices defining
            the surface, and merged cavities.
    :rtype: tuple (list, list)

    This function performs the following steps:
    1. **Selection and Filtering:** Selects non-hetero atoms from the protein, 
        calculates van der Waals radii, and performs 3D Delaunay triangulation 
        and Voronoi tessellation on the coordinates.
    2. **Surface and Interior Filtering:** Iteratively removes simplices based 
        on the user-defined radii (`r1` and `r2`) to distinguish the molecular 
        surface from the internal void space.
    3. **Surface Cavity Identification:** Detects connected void regions and 
        identifies those that remain connected to the protein surface, 
        corresponding to surface-accessible cavities and pockets.
    4. **Depth Calculation and Filtering:** Estimates cavity depth using a 
        graph-based traversal from the cavity openings, identifies the deepest
         tetrahedra, and filters cavities according to the specified depth criteria.
    5. **Output Generation:** Optionally trims cavities exceeding the specified
    	maximum depth, saves detected cavities to PDB/PQR files, and returns 
        cavity objects together with the surface representation for further 
        analysis and visualization.
       
    Example usage:
    p = parsePDB('1tqn')
    protein = p.select('protein')
    cavities, surface = calcSurfaceCavities(protein, output_path='test_surf_cav.pqr')   """

    if sparsity is not None:
        _warn("sparsity is deprecated in calcSurfaceCavities and is "
              "ignored. It thinned the mouth tetrahedra sampled as termini "
              "by the channel search; cavities are built from the unthinned "
              "ones, so it never changed them.")

    # No peel (min_enclosure=0). The enclosure peel strips the shell of true
    # exterior that a large r1 probe bridges over instead of entering, because it
    # offers a channel wide, low-cost routes along the outside of the protein. A
    # surface cavity *is* that shell: a pocket is shallow and open by definition,
    # so the peel deletes these cavities 
    cavities, surface = calcChannels(
            atoms,
            output_path=output_path,
            separate=separate,
            r1=r1, r2=r2,
            min_depth=min_depth, max_depth=max_depth,
            min_volume=min_volume, max_volume=max_volume,
            min_tetrahedra=min_tetrahedra, max_tetrahedra=max_tetrahedra,
            min_enclosure=0.0, cavities_only=True)
    
    return cavities, surface


class Channel:
    def __init__(self, tetrahedra, centerline_spline, radius_spline, length, 
                 bottleneck, volume, cost=None):
        self.tetrahedra = tetrahedra
        self.centerline_spline = centerline_spline
        self.radius_spline = radius_spline
        self.length = length
        self.bottleneck = bottleneck
        self.volume = volume
        # cost: accumulated Dijkstra path weight (sum of l / (d**2 + b) edge
        # costs) from the seed to the exit. Lower is better - a short path
        # through wide tetrahedra. Set by dijkstra(); None when not computed.
        self.cost = cost
        # curvature: path length / straight-line end-to-end distance
        # (dimensionless, >= 1; 1.0 == perfectly straight).
        self.curvature = self._computeCurvature()

    def _computeCurvature(self):
        """Path length divided by straight-line end-to-end distance."""
        x = self.centerline_spline.x
        start = np.asarray(self.centerline_spline(x[0]))
        end = np.asarray(self.centerline_spline(x[-1]))
        straight = float(np.linalg.norm(end - start))
        if straight <= 1e-9:
            return float('nan')
        return float(self.length / straight)

    def getSplines(self):
        return self.centerline_spline, self.radius_spline
    
    
class State:
    def __init__(self, simplices, neighbors, vertices):
        self.simp = simplices
        self.neigh = neighbors
        self.verti = vertices
        
    def __eq__(self, other):
        if not isinstance(other, State):
            return False
        return (np.array_equal(self.simp, other.simp) and
                np.array_equal(self.neigh, other.neigh) and
                np.array_equal(self.verti, other.verti))
        
    def setState(self, simplices, neighbors, vertices):
        self.simp = simplices
        self.neigh = neighbors
        self.verti = vertices
        
    def getState(self):
        return self.simp, self.neigh, self.verti

class Cavity:
    def __init__(self, tetrahedra, is_connected_to_surface):
        self.tetrahedra = tetrahedra
        self.is_connected_to_surface = is_connected_to_surface
        self.starting_tetrahedron = None
        self.channels = []
        self.depth = 0
        self.tetrahedra_depths = {}
        self.volume = 0.0
        
    def makeSurface(self):
        self.is_connected_to_surface = True
        
    def setExitTetrahedra(self, exit_tetrahedra, end_tetrahedra):
        self.exit_tetrahedra = exit_tetrahedra
        self.end_tetrahedra = end_tetrahedra
        
    def setStartingTetrahedron(self, tetrahedron):
        self.starting_tetrahedron = tetrahedron
        
    def setDepth(self, depth):
        self.depth = depth
        
    def addChannel(self, channel):
        self.channels.append(channel)
    
        
def _rowsIsin(a, b):
    """Boolean mask marking which rows of 2D integer array ``a`` occur as a row
    in 2D array ``b`` (exact, order-sensitive match).

    Uses a void-dtype view so each row is treated as a single hashable scalar,
    turning an O(len(a) x len(b)) row-by-row scan into an O(len(a) + len(b))
    hashed membership test.
    """
    a = np.ascontiguousarray(a)
    b = np.ascontiguousarray(b)
    if a.size == 0 or b.size == 0:
        return np.zeros(a.shape[0], dtype=bool)
    if a.dtype != b.dtype:
        b = b.astype(a.dtype)
    va = a.view(np.dtype((np.void, a.dtype.itemsize * a.shape[1]))).ravel()
    vb = b.view(np.dtype((np.void, b.dtype.itemsize * b.shape[1]))).ravel()
    return np.isin(va, vb)


class ChannelCalculator:
    def __init__(self, atoms, r2=0.9, sparsity=1, route_tolerance=1.0,
                 edge_cost='integral'):
        # Only the parameters the class actually consults are held here. r1,
        # min_depth and bottleneck are stages of the pipeline, applied to the
        # tessellation and to the finished channels by calcChannels; keeping copies
        # of them on the calculator suggested it filtered by them, which it does
        # not.
        self.atoms = atoms
        self.r2 = r2
        self.sparsity = sparsity
        self.route_tolerance = route_tolerance
        # 'integral' (clearance-profile integral) or 'bottleneck' (l/(d^2+b));
        # the Dijkstra edge weight in buildSparseGraph. Resolved per diagram by
        # calcChannels (weighted defaults to 'bottleneck').
        self.edge_cost = edge_cost
        # Filled once by buildSparseGraph and read by the channel geometry:
        # the per-simplex Voronoi-vertex clearance (the spline knots) and the
        # per-edge gate clearance on each shared Delaunay face (the reported
        # bottleneck and, later, the Dijkstra cost). Cached so the two consumers
        # share one definition of width instead of recomputing it apart.
        self._vertex_clearance = None
        self._edge_bottleneck = None

    def sphereFit(self, points, simplices, vertices, vdw_radii, r, rows=None):
        """Sum-based clearance test: for each tetrahedron, decide whether a probe
        of radius ``r`` fits at its Voronoi vertex.

        Compares the sum of the distances from the Voronoi vertex to the four
        atom centres against the sum of ``r + vdw_radius`` over the same four
        atoms. Summing over the four atoms instead of testing each one is the
        tangent-sphere test exactly when the vertex is equidistant from all four,
        and a mild relaxation of it otherwise; the erosion defaults are
        calibrated against that behaviour.

        :arg points: coordinates of all atoms, shape ``(n_atoms, 3)``.
        :type points: :class:`~numpy.ndarray`

        :arg simplices: atom indices of each tetrahedron, shape ``(n, 4)``.
        :type simplices: :class:`~numpy.ndarray`

        :arg vertices: Voronoi vertex of each tetrahedron, shape ``(n, 3)``.
        :type vertices: :class:`~numpy.ndarray`

        :arg vdw_radii: van der Waals radius of each atom.
        :type vdw_radii: :class:`~numpy.ndarray`

        :arg r: probe radius in Angstrom.
        :type r: float

        :arg rows: optional boolean mask over the tetrahedra. Rows outside it are
            reported as ``False`` without paying for the distance computation.
            :meth:`deleteSimplices3d` uses it to restrict the surface pass to the
            boundary shell.
        :type rows: :class:`~numpy.ndarray`, optional

        :returns: boolean array of length ``len(simplices)``, ``True`` where the
            probe fits.
        :rtype: :class:`~numpy.ndarray`
        """
        fits = np.zeros(len(simplices), dtype=bool)
        if rows is None:
            rows = slice(None)
        elif not rows.any():
            return fits

        atom_coords = points[simplices[rows]]                          # (m, 4, 3)
        d_sum = np.linalg.norm(
            atom_coords - vertices[rows][:, None, :], axis=2).sum(axis=1)
        r_sum = (r + vdw_radii[simplices[rows]]).sum(axis=1)
        fits[rows] = d_sum >= r_sum

        return fits

    def deleteSimplices3d(self, points, simplices, neighbors, vertices,
                          vdw_radii, r, surface):
        """Delete the tetrahedra that fail the :meth:`sphereFit` probe test and
        return the compacted tessellation.

        Which side of the test is deleted depends on the pass. The ``surface``
        pass erodes from the outside in: a boundary tetrahedron the probe fits
        into is open to the solvent, hence exterior, and goes. Only tetrahedra on
        the boundary (those with a ``-1`` neighbour) are reachable from outside,
        so the test is restricted to that shell (~n^(2/3) rows) instead of being
        evaluated over every tetrahedron on each erosion iteration, and the
        caller iterates to convergence. The inner pass instead drops every
        tetrahedron too tight for the probe, leaving the space it can actually
        occupy.

        :arg points: coordinates of all atoms, shape ``(n_atoms, 3)``.
        :type points: :class:`~numpy.ndarray`

        :arg simplices: atom indices of each tetrahedron, shape ``(n, 4)``.
        :type simplices: :class:`~numpy.ndarray`

        :arg neighbors: index of the tetrahedron opposite each vertex, ``-1`` on
            the boundary, shape ``(n, 4)``.
        :type neighbors: :class:`~numpy.ndarray`

        :arg vertices: Voronoi vertex of each tetrahedron, shape ``(n, 3)``.
        :type vertices: :class:`~numpy.ndarray`

        :arg vdw_radii: van der Waals radius of each atom.
        :type vdw_radii: :class:`~numpy.ndarray`

        :arg r: probe radius in Angstrom.
        :type r: float

        :arg surface: ``True`` for one erosion step of the surface pass, ``False``
            for the inner pass.
        :type surface: bool

        :returns: the surviving ``(simplices, neighbors, vertices)``, with
            neighbour indices remapped to the new numbering and deleted
            neighbours set to ``-1``.
        :rtype: tuple
        """
        simplices = np.asarray(simplices)
        neighbors = np.asarray(neighbors)
        vertices = np.asarray(vertices)

        n = len(simplices)
        if n == 0:
            return simplices, neighbors, vertices

        if surface:
            # Erode: a boundary tetrahedron the probe fits into is open to the
            # solvent, hence exterior. Only the boundary shell is reachable from
            # outside, so only it is tested.
            boundary = (neighbors == -1).any(axis=1)
            should_delete = self.sphereFit(points, simplices, vertices,
                                           vdw_radii, r, rows=boundary)
        else:
            # Carve: drop every tetrahedron too tight for the probe, anywhere,
            # leaving the space it can actually occupy.
            should_delete = ~self.sphereFit(points, simplices, vertices,
                                            vdw_radii, r)

        keep = ~should_delete
        simp = simplices[keep]
        neigh = neighbors[keep].copy()
        verti = vertices[keep]

        # Remap neighbour indices from the old numbering to the compacted one in a
        # single pass (deleted neighbours -> -1), replacing the previous
        # O(len(deleted) x len(neigh)) decrement loop.
        new_index = np.full(n, -1, dtype=neigh.dtype)
        new_index[keep] = np.arange(keep.sum(), dtype=neigh.dtype)
        neigh = np.where(neigh == -1, -1, new_index[neigh])

        return simp, neigh, verti

    def calcEnclosure(self, query, centers, tree=None):
        """Fraction of the directions seen from each point of ``query`` that are
        blocked by an atom within :data:`ENCLOSURE_RANGE` Angstrom.

        A local, probe-independent measure of burial: a point in the open solvent
        sees sky in most directions, a point inside a channel is surrounded
        whatever the channel's width.

        Rays are marched outwards and a ray is dropped as soon as it is blocked,
        which is what keeps this affordable: in a buried region most directions
        hit protein within the first few Angstrom, and only the few that escape
        are followed the whole way out. Marching costs
        ``ENCLOSURE_RAYS x steps`` tree queries per point and so is all but
        insensitive to how many atoms there are, whereas testing every atom in
        range against every ray costs a multiple of the atom count, and with rays
        this sparse nearly all of that work is wasted on atoms that lie near no
        ray at all.

        Pass the real atoms here, not the balls of a homogenized diagram: burial
        is a property of the protein, not of the tessellation. Atoms are all given
        the same :data:`ENCLOSURE_RADIUS`, so one tree and one plain
        nearest-neighbour test suffice. This is a burial heuristic and not a
        surface calculation, and the alternative -- a per-atom radius, which no
        nearest-neighbour query can express -- buys nothing that
        ``min_enclosure`` cannot absorb.

        :data:`ENCLOSURE_RAYS` is fixed rather than exposed, because it is part of
        the definition of the quantity and not an accuracy knob. Adding rays is
        not a free refinement: more directions discover more of the thin escape
        routes out of a channel, so the enclosure of every point drifts downwards
        and a threshold calibrated at one ray count does not carry over to
        another.

        :arg query: points to evaluate, ``(n, 3)``.
        :arg centers: atom centres.
        :arg tree: optional prebuilt :class:`~scipy.spatial.cKDTree` over
            ``centers``, to avoid rebuilding it on every call.
        :returns: ``n`` fractions in ``[0, 1]``."""
        from scipy.spatial import cKDTree

        query = np.asarray(query, dtype=float)
        if len(query) == 0:
            return np.empty(0)
        if tree is None:
            tree = cKDTree(np.asarray(centers, dtype=float))

        i = np.arange(ENCLOSURE_RAYS) + 0.5
        phi = np.arccos(1 - 2 * i / ENCLOSURE_RAYS)
        theta = np.pi * (1 + 5 ** 0.5) * i              # Fibonacci sphere
        directions = np.stack([np.sin(phi) * np.cos(theta),
                               np.sin(phi) * np.sin(theta),
                               np.cos(phi)], axis=1)

        blocked = np.zeros((len(query), ENCLOSURE_RAYS), dtype=bool)
        live = np.ones_like(blocked)
        for step in np.arange(ENCLOSURE_STEP,
                              ENCLOSURE_RANGE + ENCLOSURE_STEP, ENCLOSURE_STEP):
            point, ray = np.nonzero(live)
            if not len(point):
                break
            samples = query[point] + directions[ray] * step
            hit = tree.query(samples, distance_upper_bound=ENCLOSURE_RADIUS,
                             workers=-1)[0] <= ENCLOSURE_RADIUS
            blocked[point[hit], ray[hit]] = True
            live[point[hit], ray[hit]] = False

        return blocked.mean(axis=1)

    def peelSurfaceByEnclosure(self, points, simplices, neighbors, vertices,
                               vdw_radii, r, atom_coords, min_enclosure,
                               max_depth=None):
        """Erode the surface inward with a probe of radius ``r``, stopping where
        the tetrahedra stop being open to the solvent.

        This removes the "moat": the shell of true exterior that lies inside the
        ``r1`` surface, because an ``r1`` probe cannot enter the concavities it
        bridges over. Left in place the moat joins the cavity and offers wide,
        cheap routes along the outside of the protein.

        Neither obvious way of bounding the erosion works. A count of tetrahedron
        layers is not mesh-invariant, since a layer is one tetrahedron thick and
        tetrahedra shrink as the tessellation is refined. A depth in Angstrom is
        not ``r1``-invariant, since the moat has no constant thickness: it is as
        deep as ``r1 - r`` inside a concavity and vanishes on a flat face, so a
        depth large enough to clear it where it is thick also marches down the
        channel mouths and erodes the channels themselves. At ``r1 = 10``, a
        reasonable setting for a porin or a ribosome, that leaves almost nothing.

        The rule used here is local instead. A boundary tetrahedron is stripped
        only while it is *open*, that is while its enclosure is below
        ``min_enclosure`` (see :meth:`calcEnclosure`). The moat is open by
        construction and goes; erosion halts by itself at the first buried layer.
        ``r1`` then decides only where the erosion starts, not where it stops, so
        the result no longer depends on it, and ``r1`` is left doing the one job
        it should: capping the mouths.

        Since enclosure is a static field, the peel is really "delete the
        outside-connected component of ``{enclosure < min_enclosure}``". A
        threshold above the enclosure of a channel interior (empirically about
        0.93, as a channel is itself an escape direction) therefore percolates
        along the channels and erodes the cavity away entirely. ``max_depth`` is
        available as a hard backstop, but the default threshold leaves a wide
        margin and the failure mode is loud - no channels at all - rather than a
        plausible-looking result with the real channels missing.

        Note that the probe test and the enclosure test deliberately run against
        different spheres. ``points`` and ``vdw_radii`` are the balls the diagram
        is built on, and the probe test has to use them or its geometry stops
        agreeing with :meth:`deleteSimplices3d`. ``atom_coords`` are the real
        atoms, and the enclosure test has to use those, or burial would depend on
        the tessellation. Under ``diagram="homogenized"`` the two are not the same
        set, as some 4700 atoms become some 33000 equal balls, which would make
        the enclosure test both slower and a function of ``max_deviation``. Under
        ``"simple"`` and ``"weighted"`` they coincide.

        :arg atom_coords: the real atoms, for the enclosure test.
        :arg min_enclosure: fraction of directions that must be blocked for a
            tetrahedron to count as interior and stop the erosion. ``<= 0``
            returns the state unchanged.
        :arg max_depth: optional cap, in Angstrom, on how far the front may
            advance from the initial surface. ``None`` (default) is uncapped.
        :returns: ``(simplices, neighbors, vertices)``, compacted."""
        from scipy.spatial import cKDTree

        simplices = np.asarray(simplices)
        neighbors = np.asarray(neighbors)
        vertices = np.asarray(vertices)

        if min_enclosure <= 0 or len(simplices) == 0:
            return simplices, neighbors, vertices

        boundary = (neighbors == -1).any(axis=1)
        if not boundary.any():
            return simplices, neighbors, vertices

        # Fixed for the whole peel, so the cap bounds the total advance of the
        # front rather than its advance per pass.
        surface = cKDTree(vertices[boundary]) if max_depth is not None else None
        atoms = cKDTree(atom_coords)
        # Enclosure is a property of a point, not of the shrinking mesh, so a
        # tetrahedron re-examined on a later pass is never re-traced. Tetrahedra
        # are renumbered by the compaction below, but the four balls they are
        # built on are not, so those index the cache.
        traced = {}

        while True:
            n = len(simplices)
            if n == 0:
                break
            boundary = np.nonzero((neighbors == -1).any(axis=1))[0]
            if not len(boundary):
                break

            # The cheap tests first: does the probe fit (the same sum-based test
            # as deleteSimplices3d), and are we still inside the optional cap?
            ball_coords = points[simplices[boundary]]
            d_sum = np.linalg.norm(
                ball_coords - vertices[boundary][:, None, :], axis=2).sum(axis=1)
            r_sum = (r + vdw_radii[simplices[boundary]]).sum(axis=1)
            candidate = d_sum >= r_sum
            if max_depth is not None:
                candidate &= surface.query(vertices[boundary])[0] <= max_depth
            candidates = boundary[candidate]
            if not len(candidates):
                break

            # Ray tracing runs only on what survived those, and only once each.
            keys = [tuple(key) for key in simplices[candidates]]
            fresh = [i for i, key in enumerate(keys) if key not in traced]
            if fresh:
                values = self.calcEnclosure(vertices[candidates[fresh]],
                                            atom_coords, tree=atoms)
                for i, value in zip(fresh, values):
                    traced[keys[i]] = value
            enclosure = np.array([traced[key] for key in keys])

            should_delete = np.zeros(n, dtype=bool)
            should_delete[candidates[enclosure < min_enclosure]] = True
            if not should_delete.any():
                break

            keep = ~should_delete
            simplices = simplices[keep]
            neigh = neighbors[keep].copy()
            vertices = vertices[keep]

            new_index = np.full(n, -1, dtype=neigh.dtype)
            new_index[keep] = np.arange(keep.sum(), dtype=neigh.dtype)
            neighbors = np.where(neigh == -1, -1, new_index[neigh])

        return simplices, neighbors, vertices

    def deleteSection(self, simplices_subset, simplices, neighbors, vertices,
                      reverse=False):
        simplices = np.asarray(simplices)
        neighbors = np.asarray(neighbors)
        vertices = np.asarray(vertices)

        n = len(simplices)
        if n == 0:
            return simplices, neighbors, vertices

        # Which rows of `simplices` also appear in `simplices_subset` (exact,
        # order-sensitive row match via hashed membership instead of the former
        # O(n x len(subset)) per-row scan).
        matches = _rowsIsin(simplices, np.asarray(simplices_subset))
        keep = matches if reverse else ~matches

        simp = simplices[keep]
        neigh = neighbors[keep].copy()
        verti = vertices[keep]

        new_index = np.full(n, -1, dtype=neigh.dtype)
        new_index[keep] = np.arange(keep.sum(), dtype=neigh.dtype)
        neigh = np.where(neigh == -1, -1, new_index[neigh])

        return simp, neigh, verti

    def getVdwRadii(self, atoms):
        vdw_radii_dict = {
            'H': 1.20, 'HE': 1.40, 'LI': 1.82, 'BE': 1.53, 'B': 1.92, 'C': 1.70,
            'N': 1.55, 'O': 1.52, 'F': 1.47, 'NE': 1.54, 'NA': 2.27, 'MG': 1.73,
            'AL': 1.84, 'SI': 2.10, 'P': 1.80, 'S': 1.80, 'CL': 1.75, 'AR': 1.88,
            'K': 2.75, 'CA': 2.31, 'SC': 2.11, 'NI': 1.63, 'CU': 1.40, 'ZN': 1.39,
            'GA': 1.87, 'GE': 2.11, 'AS': 1.85, 'SE': 1.90, 'BR': 1.85, 'KR': 2.02,
            'RB': 3.03, 'SR': 2.49, 'PD': 1.63, 'AG': 1.72, 'CD': 1.58, 'IN': 1.93,
            'SN': 2.17, 'SB': 2.06, 'TE': 2.06, 'I': 1.98, 'XE': 2.16, 'CS': 3.43,
            'BA': 2.68, 'PT': 1.75, 'AU': 1.66, 'HG': 1.55, 'TL': 1.96, 'PB': 2.02,
            'BI': 2.07, 'PO': 1.97, 'AT': 2.02, 'RN': 2.20, 'FR': 3.48, 'RA': 2.83,
            'U': 1.86, 'FE': 2.44
        }
        
        return np.array([vdw_radii_dict[atom] for atom in atoms])

    def _fibonacciSphere(self, n):
        """Return ``n`` roughly evenly distributed unit vectors on a sphere using
        the Fibonacci (golden spiral) lattice."""
        n = int(np.maximum(1, n))
        indices = np.arange(n) + 0.5
        phi = np.arccos(1.0 - 2.0 * indices / n)
        theta = np.pi * (1.0 + 5.0 ** 0.5) * indices
        x = np.sin(phi) * np.cos(theta)
        y = np.sin(phi) * np.sin(theta)
        z = np.cos(phi)
        return np.stack([x, y, z], axis=1)

    def _shellPointCount(self, rad, rho, max_deviation):
        """Number of equal balls of radius ``rho`` to place on a shell of radius
        ``rad`` so that the outer envelope of their union stays within
        ``max_deviation`` of ``rad + rho``.

        A ball centered at radius ``rad`` only touches the target sphere of radius
        ``rad + rho`` at a single point, so the shell must be sampled densely
        enough that the "valleys" between neighbouring balls do not dip more than
        ``max_deviation``. Each ball covers a spherical cap of half-angle
        ``alpha`` on the ``rad + rho - max_deviation`` sphere; the count is the
        number of such caps needed to tile the sphere (with an overlap factor).
        """
        r = rad + rho - max_deviation
        cos_alpha = (rad * rad + r * r - rho * rho) / (2.0 * rad * r)
        cos_alpha = float(np.clip(cos_alpha, -1.0, 1.0))
        if cos_alpha >= 1.0:
            return 1
        # The exact number of caps to tile the sphere is 2 / (1 - cos_alpha); we use
        # a factor of 4 (a ~2x overlap margin). This is NOT slack to be trimmed: the
        # Fibonacci lattice is not an optimal packing and its coverage efficiency
        # degrades as the shell (and count) grows, so the factor needed to actually
        # hold the max_deviation bound increases with atom size. A single constant
        # must therefore be sized for the largest atoms (e.g. metals); 4 keeps the
        # measured dip below max_deviation across the whole range, whereas 3 already
        # fails for anything larger than the thinnest shell. Lowering it silently
        # breaks large-atom accuracy - tune max_deviation instead to change cost.
        return int(np.ceil(4.0 / (1.0 - cos_alpha)))

    def homogenizeAtoms(self, coords, vdw_radii, max_deviation=0.2):
        """Substitute every atom by a set of homogeneous balls whose common radius
        equals the smallest van der Waals radius present in the structure.

        Each atom of radius ``R`` is replaced by a collection of overlapping balls
        of radius ``rho = min(vdw_radii)`` arranged on concentric shells (plus a
        central ball) so that their union approximates the original atomic sphere
        to within ``max_deviation``. Because all resulting balls share the same
        radius, an ordinary Voronoi / Delaunay tessellation of their centers yields
        an accurate estimate of the additively weighted (power) Voronoi diagram of
        the original atoms. This is the approach used by MolAxis and CAVER 3 and
        avoids simply discarding the smaller (e.g. hydrogen) atoms.

        Atoms whose radius is within ``max_deviation`` of ``rho`` are kept as a
        single ball, so a structure of similarly sized atoms is left essentially
        unchanged while a structure containing hydrogens (small ``rho``) fills its
        larger atoms with several balls.

        :arg coords: atomic coordinates, shape ``(N, 3)``
        :arg vdw_radii: per-atom van der Waals radii, shape ``(N,)``
        :arg max_deviation: maximum tolerated deviation (in Angstrom) between the
            union surface of the substitute balls and the original atomic surface.
            Smaller values are more accurate but generate more balls. Default 0.2.
        :returns: a tuple ``(new_coords, new_radii)`` where every entry of
            ``new_radii`` equals ``rho``.
        """
        coords = np.asarray(coords, dtype=float)
        vdw_radii = np.asarray(vdw_radii, dtype=float)

        rho = float(np.min(vdw_radii))
        tol = 1e-6
        new_points = []

        for center, R in zip(coords, vdw_radii):
            # Atoms within max_deviation of the smallest radius stay a single ball.
            if R - rho <= max_deviation + tol:
                new_points.append(center)
                continue

            # Central ball plus concentric shells stepped by rho, with the
            # outermost shell at (R - rho) so the union surface reaches R.
            new_points.append(center)
            shell_radii = list(np.arange(rho, R - rho, rho))
            if not shell_radii or (R - rho) - shell_radii[-1] > tol:
                shell_radii.append(R - rho)

            for rad in shell_radii:
                if rad <= tol:
                    continue
                n = self._shellPointCount(rad, rho, max_deviation)
                new_points.extend(center + rad * self._fibonacciSphere(n))

        new_points = np.array(new_points)
        new_radii = np.full(len(new_points), rho)

        return new_points, new_radii

    def buildSurfaceDepthOracle(self, coords, vdw_radii, r1, max_deviation, max_depth):
        """Interior/exterior depth oracle for relabeling the additively-weighted
        diagram's surface mouths (``diagram="weighted"``).

        The AW tessellation is not a clean simplicial complex: many interior 3-ball
        faces are left unpaired and masquerade as surface boundaries, so channels
        truncate to stubs. This builds a *homogenized* Voronoi diagram of the same
        atoms (a clean simplicial complex), erodes it with an ``r1`` probe to separate
        solvent (exterior) from protein (interior), and labels every tetrahedron by its
        geodesic distance (A) below the molecular surface (0 = exterior/solvent).
        :meth:`getSurfaceCavities` then keeps only AW exit tetrahedra whose Voronoi
        vertex maps (via ``find_simplex``) to depth ``<= max_depth`` Angstrom.

        :returns: ``(delaunay, depth, max_depth)`` -- the homogenized
            :class:`~scipy.spatial.Delaunay`, its per-tetrahedron geodesic depth in
            Angstrom (points outside the hull are treated as depth 0), and the
            passed-through threshold.
        """
        from scipy.spatial import Delaunay

        hp, hrho = self.homogenizeAtoms(coords, vdw_radii, max_deviation)
        delaunay = Delaunay(hp)
        centers = self.calcCircumcenters(delaunay)
        clearance = (np.linalg.norm(hp[delaunay.simplices] - centers[:, None, :], axis=2)
                     - hrho[delaunay.simplices]).min(axis=1)
        neighbors = delaunay.neighbors
        n = len(delaunay.simplices)

        # r1 surface erosion: peel boundary tetrahedra wide enough for the probe, from
        # the hull inward, until nothing more can be removed. Survivors == interior.
        alive = np.ones(n, dtype=bool)
        while True:
            dead = np.zeros(n, dtype=bool)
            for k in range(4):
                col = neighbors[:, k]
                dead |= (col == -1) | ((col >= 0) & ~alive[col])
            peel = alive & dead & (clearance >= r1)
            if not peel.any():
                break
            alive[peel] = False

        # Geodesic depth (A) below the surface: shortest path from the exterior
        # (solvent) tetrahedra inward along Voronoi edges. A physical distance, not a
        # tetrahedron-layer count, so weighted_mouth_depth is a real Angstrom threshold
        # that does not drift with the homogenization density.
        degenerate = self._degenerateTetrahedra(delaunay.simplices, centers, hp)
        scratch = np.full(n, -1, dtype=np.intp)
        depth = self._geodesicDepth(np.arange(n), np.nonzero(~alive)[0], neighbors,
                                    centers, degenerate, scratch)
        # Enclosed pockets never reached from the exterior are deep (never a mouth).
        reached = depth[np.isfinite(depth)]
        depth[~np.isfinite(depth)] = (float(reached.max()) + 5.0) if reached.size \
            else float(max_depth + 1)

        return delaunay, depth, max_depth
    
    def surfaceLayer(self, shape_simplices, filtered_simplices, shape_neighbors):
        shape_simplices = np.asarray(shape_simplices)
        shape_neighbors = np.asarray(shape_neighbors)
        filtered_simplices = np.asarray(filtered_simplices)

        # Split simplices into those touching the boundary (a -1 neighbour) and
        # the interior ones, preserving order.
        boundary = (shape_neighbors == -1).any(axis=1)
        surface_simplices = shape_simplices[boundary]
        surface_neighbors = shape_neighbors[boundary]
        interior_simplices = shape_simplices[~boundary]

        # Row-membership tests replace the former (N, M, 4) broadcast temporaries.
        surf_keep = _rowsIsin(surface_simplices, filtered_simplices)
        filtered_surface_simplices = surface_simplices[surf_keep]
        filtered_surface_neighbors = surface_neighbors[surf_keep]

        filtered_surface_neighbors = np.unique(filtered_surface_neighbors)
        filtered_surface_neighbors = filtered_surface_neighbors[filtered_surface_neighbors != 0]

        filtered_interior_simplices = interior_simplices[
            _rowsIsin(interior_simplices, filtered_simplices)]

        surface_layer_neighbor_simplices = shape_simplices[filtered_surface_neighbors]

        second_layer = filtered_interior_simplices[
            _rowsIsin(filtered_interior_simplices, surface_layer_neighbor_simplices)]

        return filtered_surface_simplices, second_layer

            
    def findGroups(self, neigh, is_cavity=True):
        x = neigh.shape[0]
        visited = np.zeros(x, dtype=bool)
        groups = []

        def dfs(tetra_index):
            stack = [tetra_index]
            current_group = []
            while stack:
                index = stack.pop()
                if not visited[index]:
                    visited[index] = True
                    current_group.append(index)
                    stack.extend(neighbor for neighbor in neigh[index] if neighbor != -1 and not visited[neighbor])
            return np.array(current_group)

        for i in range(x):
            if not visited[i]:
                current_group = dfs(i)
                if is_cavity:
                    groups.append(Cavity(current_group, False))
                else:
                    groups.append(current_group)
            
        return groups

    def getSurfaceCavities(self, cavities, interior_simplices, second_layer, 
                           state, points, vdw_radii, sparsity, mouth_oracle=None):
        surface_cavities = []
        
        for cavity in cavities:
            tetrahedra = cavity.tetrahedra
            second_layer_mask = np.isin(interior_simplices[tetrahedra], second_layer).all(axis=1)
            
            if np.any(second_layer_mask):
                exit_tetrahedra = tetrahedra[second_layer_mask]
                if mouth_oracle is not None:
                    # diagram="weighted": drop the false (buried) mouths the leaky AW
                    # diagram produces. Keep an exit tetrahedron only if its Voronoi
                    # vertex lies within max_depth Angstrom (geodesic) of the true
                    # molecular surface, per a homogenized interior/exterior oracle.
                    delaunay, depth, max_depth = mouth_oracle
                    located = delaunay.find_simplex(state.verti[exit_tetrahedra])
                    surface_depth = np.where(located >= 0, depth[located.clip(0)], 0.0)
                    exit_tetrahedra = exit_tetrahedra[surface_depth <= max_depth]
                    if len(exit_tetrahedra) == 0:
                        continue
                cavity.makeSurface()
                end_tetrahedra = self.getEndTetrahedra(exit_tetrahedra, state.verti, points, vdw_radii, state.simp, sparsity)
                cavity.setExitTetrahedra(exit_tetrahedra, end_tetrahedra)
                surface_cavities.append(cavity)
                
        return surface_cavities


    def mergeCavities(self, cavities, simplices):
        if not cavities:
            # No cavities survived filtering (e.g. restrict_channels_to_start_point
            # selected a single cavity shallower than min_depth). Return an empty
            # (0, 4) slice so the pipeline yields zero channels instead of crashing
            # in np.concatenate on an empty list.
            return simplices[np.empty(0, dtype=np.intp)]
        merged_tetrahedra = np.concatenate([cavity.tetrahedra for cavity in cavities])
        return simplices[merged_tetrahedra]

    def _degenerateTetrahedra(self, simplices, vertices, points):
        # Scale-invariant flatness flag for the geodesic depth graph. A near-flat
        # tetrahedron has a runaway circumcenter, so its incident Voronoi edges can be
        # astronomically long (measured up to ~1e14 A) and would corrupt a shortest-path
        # depth. Flatness is the ratio of circumradius R to the tetrahedron's own atom
        # span L: well-shaped cells sit at R/L < ~4 whatever their absolute size - the
        # large fat cells that span a wide pore under a big probe included - while flat
        # slivers diverge to R/L -> infinity. Using the dimensionless R/L, never an
        # absolute length, keeps this correct for porins, ribosome tunnels and
        # large-probe (e.g. r1=20) runs alike. Edges touching a flagged tetrahedron are
        # dropped from the graph in _geodesicDepth.
        apex = points[simplices]                                  # (n, 4, 3)
        R = np.linalg.norm(apex[:, 0] - vertices, axis=1)         # circumradius
        L = np.zeros(len(simplices))
        for a in range(4):
            for b in range(a + 1, 4):
                L = np.maximum(L, np.linalg.norm(apex[:, a] - apex[:, b], axis=1))
        # Well-shaped cells measure R/L up to ~3.7 (p99.9); flat slivers reach ~1e14.
        # Flag above 5 - the ~14-order gap makes the exact cut irrelevant across [5, 100].
        return R / np.maximum(L, 1e-9) > 5.0

    def _geodesicDepth(self, tetrahedra, sources, neighbors, vertices, degenerate, scratch):
        # Shortest-path distance (A) from the source set to every tetrahedron in
        # `tetrahedra`, along Voronoi edges weighted by circumcenter-to-circumcenter
        # distance (a multi-source Dijkstra over the induced subgraph). Edges touching a
        # degenerate (runaway-circumcenter) tetrahedron are dropped. `scratch` is a
        # reusable global->local index buffer (-1 outside `tetrahedra`), reset before
        # return so the caller can pass it again. Returns the distance aligned to
        # `tetrahedra`, np.inf where a tetrahedron is unreachable from every source.
        from scipy.sparse import csr_matrix
        from scipy.sparse.csgraph import dijkstra

        T = np.asarray(tetrahedra, dtype=np.intp)
        m = len(T)
        if m == 0:
            return np.empty(0)
        scratch[T] = np.arange(m)
        nb = neighbors[T]                                         # (m, deg) global
        safe_nb = np.where(nb < 0, 0, nb)
        local = np.where(nb < 0, -1, scratch[safe_nb])
        keep = local >= 0
        keep &= ~degenerate[T][:, None]                          # drop from a flat tetra
        keep &= ~np.where(nb < 0, True, degenerate[safe_nb])     # drop into a flat tetra
        keep = keep.ravel()
        row = np.repeat(np.arange(m), nb.shape[1])[keep]
        col = local.ravel()[keep]
        w = np.linalg.norm(vertices[T[row]] - vertices[T[col]], axis=1)
        graph = csr_matrix((w, (row, col)), shape=(m, m))
        src = scratch[np.asarray(sources, dtype=np.intp)]
        src = src[src >= 0]
        scratch[T] = -1                                          # reset for reuse
        if len(src) == 0:
            return np.full(m, np.inf)
        return dijkstra(graph, indices=src, min_only=True)

    def findDeepestTetrahedra(self, cavities, neighbors, vertices, points, simplices):
        # Cavity depth = the geodesic distance (A) from the cavity's openings (its exit
        # tetrahedra) to its farthest point, as a shortest path along Voronoi edges. This
        # is a physical length independent of tetrahedron size, so min_depth is mesh
        # invariant; the old +1-per-tetrahedron layer count grew as the mesh refined.
        degenerate = self._degenerateTetrahedra(simplices, vertices, points)
        scratch = np.full(neighbors.shape[0], -1, dtype=np.intp)
        for cavity in cavities:
            tetra = np.asarray(cavity.tetrahedra, dtype=np.intp)
            dist = self._geodesicDepth(tetra, cavity.exit_tetrahedra, neighbors,
                                       vertices, degenerate, scratch)
            finite = np.isfinite(dist)
            if not finite.any():
                # No exit reached any tetrahedron (degenerate cavity); keep it minimal.
                cavity.setStartingTetrahedron(np.array([int(tetra[0])]))
                cavity.setDepth(0.0)
                cavity.tetrahedra_depths = {int(tetra[0]): 0.0}
                continue
            deepest = int(np.argmax(np.where(finite, dist, -np.inf)))
            cavity.setStartingTetrahedron(np.array([int(tetra[deepest])]))
            cavity.setDepth(float(dist[deepest]))
            cavity.tetrahedra_depths = {int(tetra[k]): float(dist[k])
                                        for k in np.nonzero(finite)[0]}
            
    def calcCircumcenters(self, dela):
        # per-simplex circumcenters recovered analytically from the Delaunay 
        # paraboloid lifting, avoiding a second Qhull pass. Identical to scipy 
        # Voronoi vertices in general position.
        eq = dela.equations
        scale = dela.paraboloid_scale
        centers = -eq[:, :-2] / (2 * scale * eq[:, -2][:, None])
        return centers

    def _edgeBottleneck(self, ci, cj, shared_atoms, points, vdw_radii):
        # Minimum clearance along the Voronoi edge - the segment between the two
        # circumcenters ci, cj, dual to the Delaunay face the two tetrahedra
        # share - measured against that face's atoms. This is the edge bottleneck
        # radius: the circumcenters are local clearance maxima, so the tightest
        # point of the segment (the gate) generally lies between them and is
        # narrower than either endpoint.
        #
        # min over t in [0, 1] and over the shared atoms of |p(t) - a| - vdw(a),
        # with p(t) = ci + t (cj - ci). Because the min over the (t, atom)
        # product equals the min of the per-atom minima, each atom reduces to an
        # independent clamped point-to-segment distance - a closed form, no
        # sampling. The foot clamps to an endpoint when it falls outside the
        # segment, so the gate value is always <= both vertex clearances. Exact
        # for the straight edges of the homogenized/simple diagrams; for the
        # weighted (Apollonius) diagram the true edge is a slight arc and the
        # chord is a local approximation.
        a = points[shared_atoms]
        r = vdw_radii[shared_atoms]
        u = cj - ci
        uu = float(u @ u)
        if uu <= 1e-12:
            # Twin tetrahedra: the circumcenters coincide, the edge is a point,
            # so return the shared vertex clearance directly. This guard is
            # load-bearing, not defensive boilerplate - do not remove it. Two
            # near-cospherical tetrahedra (the "T5" degeneracy CAVER's paper
            # jitters against) produce coincident circumcenters; measured on real
            # structures these are absent at the default max_deviation but do
            # appear - a handful per structure - as the diagram is refined
            # (max_deviation -> 0.02). scipy/Qhull's cospherical facet merging
            # keeps the geometry finite (no NaN/inf circumcenters) so we do not
            # need CAVER's jitter/random-rotation precautions, but the coincident
            # circumcenters would still divide by ~0 here. Just above the
            # threshold the clip below handles it: for a tiny but non-zero edge
            # the foot clamps to an endpoint and the gate degrades continuously
            # to the endpoint clearance.
            return float(np.min(np.linalg.norm(a - ci, axis=1) - r))
        t = np.clip((a - ci) @ u / uu, 0.0, 1.0)
        p = ci + t[:, None] * u
        return float(np.min(np.linalg.norm(p - a, axis=1) - r))

    def _edgeBottleneckBatch(self, ci, cj, shared_atoms, points, vdw_radii):
        # Vectorized _edgeBottleneck over F edges at once. ``ci``, ``cj`` are
        # (F, 3) circumcenters (``ci`` the lower-index endpoint, so the reverse
        # edge yields a bitwise-identical gate) and ``shared_atoms`` is (F, 3)
        # atom indices of the shared Delaunay face. Same closed-form clamped
        # point-to-segment as the scalar version, with the twin-tetrahedron
        # (coincident circumcenter) guard applied row-wise. See _edgeBottleneck.
        # Returns ``(gate, tstar)``: the min clearance and the parameter t in
        # [0, 1] where it is attained (the binding atom's clamped foot), the
        # latter used by _edgeCostIntegralBatch to force a quadrature node on
        # the pinch.
        a = points[shared_atoms]                          # (F, 3, 3)
        r = vdw_radii[shared_atoms]                       # (F, 3)
        u = cj - ci                                       # (F, 3)
        uu = np.einsum('fj,fj->f', u, u)                  # (F,)
        twin = uu <= 1e-12
        uu_safe = np.where(twin, 1.0, uu)
        diff = a - ci[:, None, :]                         # (F, 3, 3)
        t = np.clip(np.einsum('faj,fj->fa', diff, u) / uu_safe[:, None],
                    0.0, 1.0)                             # (F, 3)
        p = ci[:, None, :] + t[:, :, None] * u[:, None, :]
        clr = np.linalg.norm(p - a, axis=2) - r          # (F, 3)
        amin = clr.argmin(axis=1)
        rows = np.arange(len(u))
        gate = clr[rows, amin]
        tstar = t[rows, amin]
        if twin.any():
            # coincident circumcenters: the edge is a point, so the gate is the
            # shared vertex clearance measured at ci (== cj) and tstar is moot.
            twin_gate = (np.linalg.norm(diff, axis=2) - r).min(axis=1)
            gate = np.where(twin, twin_gate, gate)
            tstar = np.where(twin, 0.0, tstar)
        return gate, tstar

    def _edgeCostIntegralBatch(self, ci, cj, tstar, shared_atoms, points,
                               vdw_radii, z=2.0, delta=0.3, r_floor=1e-2):
        # Price the edge by the integral of its clearance profile r(t)^-z, not by
        # the whole length charged at its single narrowest point (the l/gate^2 MOLE
        # cost). This profile-integral idea follows CAVER (TCBB'15, Eq. 1), but this
        # is NOT a CAVER reimplementation - the quadrature deliberately differs (see
        # below). r(t) = min over the 3 shared-face balls of |ci + t (cj-ci) - a| -
        # vdw is the exact clearance profile of the straight (homogenized/simple)
        # edge. The integral is additive under subdivision, so unlike l/gate^2 the
        # cost is mesh-invariant; the vertex-only formula's error grows with edge
        # length (measured ~19% -> ~1500% p95 across length bins) and drifts the
        # routing as max_deviation coarsens.
        #
        # The quadrature differs from CAVER's: where CAVER samples a plain uniform
        # grid (and takes its grid-minimum as the bottleneck), we take the uniform
        # grid linspace(0, 1, K) AND force a node at the exact analytic gate t*, so
        # the narrowest point is never missed - and the reported edge bottleneck is
        # that exact gate (see _edgeBottleneckBatch), not a grid sample. K =
        # ceil(L/delta) is a fixed ARCLENGTH step in Angstrom, so the sample count
        # scales with physical length (what makes it mesh-invariant). A short edge
        # (L <= delta) collapses to {0, t*, 1}, the three clearances already in
        # hand. r(t) is floored at r_floor so the integrand cannot diverge on a
        # sub-r2 edge that dips through an atom (the integral's analog of the
        # l/(d^2+b) regularizer); traversable edges have r(t) >= r2 >> r_floor and
        # are untouched. Exact for straight edges; a chord approximation for the
        # weighted (Apollonius) diagram, whose edges are arcs.
        chunk = 20000
        u = cj - ci
        L = np.linalg.norm(u, axis=1)
        cost = np.zeros(len(ci))
        alive = np.nonzero(L > 1e-6)[0]                   # twin/zero-length -> 0
        K = np.maximum(2, np.ceil(L / delta).astype(int) + 1)
        for k in np.unique(K[alive]):
            bucket = alive[K[alive] == k]
            grid = np.linspace(0.0, 1.0, k)
            for s in range(0, len(bucket), chunk):
                ii = bucket[s:s + chunk]
                nodes = np.sort(np.concatenate(
                    [np.broadcast_to(grid, (len(ii), k)), tstar[ii, None]],
                    axis=1), axis=1)                      # (n, k+1)
                p = ci[ii][:, None, :] + nodes[:, :, None] * u[ii][:, None, :]
                a = points[shared_atoms[ii]]              # (n, 3, 3)
                rr = np.linalg.norm(p[:, :, None, :] - a[:, None, :, :], axis=3) \
                    - vdw_radii[shared_atoms[ii]][:, None, :]
                rr = np.maximum(rr.min(axis=2), r_floor)  # (n, k+1) profile
                cost[ii] = np.trapz(rr ** (-z), nodes, axis=1) * L[ii]
        return cost

    def buildSparseGraph(self, simplices, neighbors, vertices, points, vdw_radii):
        # One weighted CSR adjacency matrix for the whole cleared state, built
        # from array ops over the (N, deg) neighbour table - no Python loop.
        # Edge (tetra -> neigh) weight is l / (d**2 + b), where l is the
        # vertex-to-vertex distance and d is the gate clearance on the shared
        # Delaunay face (min clearance along the connecting Voronoi edge). The
        # gate is the width the cost should see: the clearance between the two
        # circumcenters, not the entered node's own vertex clearance, which is a
        # local maximum and lets the search prefer a route that is actually
        # narrower at a face it never measures. Because the gate is symmetric,
        # the cost no longer depends on traversal direction the way the
        # entered-node vertex clearance did.
        #
        # Two edge-cost modes (self.edge_cost): 'bottleneck' is the l/(d**2 + b)
        # above (d = gate); 'integral' replaces it on face edges with a clearance-
        # profile integral (_edgeCostIntegralBatch), which is mesh-invariant.
        # Either way the gate is still cached as the reported edge bottleneck.
        from scipy.sparse import csr_matrix

        simplices = np.asarray(simplices)
        neighbors = np.asarray(neighbors)
        N, deg = neighbors.shape

        tetra_points = points[simplices]
        distances = np.linalg.norm(tetra_points - vertices[:, None, :], axis=2)
        bottleneck = np.min(distances - vdw_radii[simplices], axis=1)
        # Cache the per-simplex vertex clearance: the channel geometry otherwise
        # recomputes this identical min over each path's tetrahedra.
        self._vertex_clearance = bottleneck

        # Directed edge list straight off the neighbour table.
        rows = np.repeat(np.arange(N), deg)
        cols = neighbors.ravel()
        keep = cols != -1
        rows = rows[keep]
        cols = cols[keep]

        # Drive the gate geometry from the lower-index endpoint so an edge's two
        # directed copies see identical inputs and get a bitwise identical,
        # direction-symmetric gate - i.e. compute each undirected edge once.
        lo = np.minimum(rows, cols)
        hi = np.maximum(rows, cols)

        # Shared 3 atoms of each Delaunay face, convention-agnostic (set
        # intersection, not scipy's opposite-vertex rule, so it holds for the
        # weighted diagram too): ``present`` marks which of the lower simplex's
        # four atoms also appear in the higher one; a face-adjacent edge has 3.
        slo = simplices[lo]
        shi = simplices[hi]
        present = (slo[:, :, None] == shi[:, None, :]).any(axis=2)   # (M, 4)
        face = present.sum(axis=1) == 3

        # Rare non-face links fall back to the entered node's vertex clearance;
        # face links (essentially all of them) overwrite it with the gate.
        d = bottleneck[cols].astype(float, copy=True)
        fi = tstar = shared = None
        if face.any():
            fi = np.nonzero(face)[0]
            shared = slo[fi][present[fi]].reshape(-1, 3)
            d[fi], tstar = self._edgeBottleneckBatch(vertices[lo[fi]],
                                                     vertices[hi[fi]], shared,
                                                     points, vdw_radii)

        l = np.linalg.norm(vertices[rows] - vertices[cols], axis=1)
        b = 1e-3
        weight = l / (d * d + b)          # 'bottleneck' cost; non-face fallback
        if self.edge_cost == 'integral' and fi is not None:
            # 'integral' cost: replace the l/(d^2+b) of every face edge (the
            # bottleneck-only MOLE cost) with the clearance-profile integral. The
            # l/(d^2+b) fallback stays on the rare non-face links.
            #
            # No R/L flatness guard is needed here (unlike buildSurfaceDepthOracle,
            # which runs on the full diagram): this graph is the *cleared interior*
            # state, and a runaway circumcenter means a huge clearance, so the r1
            # erosion has already stripped every flat boundary tetra - measured 0
            # degenerate tetra and max edge ~4 A in the cleared graph. That matters
            # because l/(d^2+b) *over*-prices a runaway edge (huge l) so the search
            # avoids it for free, whereas the integral would *under*-price it (r(t)
            # is huge over most of a runaway edge, so r(t)^-z ~ 0 there).
            weight[fi] = self._edgeCostIntegralBatch(
                vertices[lo[fi]], vertices[hi[fi]], tstar, shared, points, vdw_radii)

        # Per-edge gate cache (unordered key) read by _pathGates: each face edge
        # stored once, keyed (lo, hi). The reported bottleneck reads the same map.
        undirected = face & (rows < cols)
        self._edge_bottleneck = {
            (int(i), int(j)): float(v)
            for i, j, v in zip(rows[undirected], cols[undirected],
                               d[undirected])
        }

        return csr_matrix((weight, (rows, cols)), shape=(N, N))

    def dijkstra(self, cavity, graph, simplices, neighbors, vertices, points,
                 vdw_radii, truncate_at_surface=True, similarity=0.8):
        # a single multi-target Dijkstra from the seed over the cavity subgraph,
        # then every exit path reconstructed from the predecessor tree - 
        # instead of one heap search per (seed, exit) pair.
        # Channel geometry still goes through the current 
        # process_channel/Channel (Simpson-based volume).
        from scipy.sparse.csgraph import dijkstra
        from collections import defaultdict

        cavity_tetra = np.asarray(cavity.tetrahedra)
        if len(cavity_tetra) == 0:
            return
        global_to_local = {tetra: i for i, tetra in enumerate(cavity_tetra)}
        cavity_graph = graph[np.ix_(cavity_tetra, cavity_tetra)]

        # A tunnel ends at the surface, but the Dijkstra cost has no such term:
        # it rewards width, and the widest places are the surface grooves. Left
        # free, the cheapest path to a far exit leaves the pocket at one mouth,
        # runs along the outside and re-enters at another - which is not a
        # tunnel. So a mouth must not *conduct*, only *absorb*: we drop the
        # outgoing edges of every mouth before the search, making "a channel
        # ends at its first surface contact" a hard constraint of the search
        # rather than a cut applied afterwards to the winning path.
        # The ordering is the whole point. Truncating after selection cuts a
        # path that was itself chosen *because* it ran along the surface, while
        # a genuine narrow interior corridor to the same mouth loses the argmin
        # to that groove and is never enumerated at all - it vanishes from the
        # output even though it is open. Mesh refinement makes the groove
        # cheaper, so interior tunnels drop out one by one as max_deviation
        # shrinks; forbidding transit removes that dependence entirely.
        # A mouth is a surface (exit) tetrahedron a probe of the traversal
        # radius r2 can leave through. Note the gate is r2, not bottleneck:
        # bottleneck is a reporting filter, and letting it decide which mouths
        # absorb would let it silently re-open narrow mouths as transit nodes,
        # i.e. change the routes rather than filter them. For the homogenized
        # and weighted diagrams every surviving tetrahedron already clears r2 by
        # construction (equal radii + equidistant circumcenter collapse the
        # sum-based test in deleteSimplices3d to the per-atom clearance), so the
        # test is a no-op there; it earns its keep for diagram="simple", where
        # unequal radii break that identity.
        # Local indices of the tetrahedra a channel is allowed to end at.
        terminals_local = [global_to_local[int(t)]
                           for t in np.asarray(cavity.end_tetrahedra)
                           if int(t) in global_to_local]
        if truncate_at_surface:
            exit_tetra = np.asarray(getattr(cavity, 'exit_tetrahedra',
                                            np.empty(0, dtype=np.intp)))
            if len(exit_tetra):
                verts = vertices[exit_tetra]
                atom_pos = points[simplices[exit_tetra]]
                atom_rad = vdw_radii[simplices[exit_tetra]]
                clearance = (np.linalg.norm(atom_pos - verts[:, None, :],
                                            axis=2) - atom_rad).min(axis=1)
                seeds = set(int(s) for s in cavity.starting_tetrahedron)

                # Only the mouths themselves absorb. A tetrahedron that merely lies
                # inside a mouth's inscribed ball must NOT be absorbed: the ball's
                # radius is the clearance (up to ~2 A) and it reaches inward as
                # well as outward, so absorbing on it eats the corridors that
                # approach the surface and truncates real tunnels before they
                # arrive - measured to delete both known side tunnels at
                # max_deviation=0.1 while keeping them at 0.02, i.e. exactly the
                # silent, mesh-dependent tunnel loss this design exists to prevent.
                # A path can consequently still slip *past* a mouth through a twin
                # tetrahedron - a neighbour sharing almost the same circumcenter,
                # not itself in the second layer and so still conducting - and
                # surface again somewhere else. That leak is real but narrow (the
                # twins sit 0.1-0.7 A from a mouth, in the surface shell at depth
                # 1-4), and such a path always passes through the exit sphere of a
                # channel that is already reported. It is therefore handled in
                # _addDedupedChannels, which cuts a path at the first reported exit
                # sphere it enters - the point where it truly leaves the protein -
                # rather than walling the graph off against every mouth.
                absorbing = [global_to_local[int(t)]
                             for t, c in zip(exit_tetra, clearance)
                             if c >= self.r2 and int(t) in global_to_local
                             and int(t) not in seeds]
                if absorbing:
                    # Zero the mouths' rows: edges *into* a mouth survive (a
                    # channel may end there), edges *out of* it are gone.
                    cavity_graph = cavity_graph.tolil()
                    for i in absorbing:
                        cavity_graph.rows[i] = []
                        cavity_graph.data[i] = []
                    cavity_graph = cavity_graph.tocsr()
                # Every mouth is a terminus; the dedup decides which of them are
                # one opening. See the comment at the target loop below.
                terminals_local = absorbing

        candidates = []

        for start_global in cavity.starting_tetrahedron:
            if start_global not in global_to_local:
                continue
            start_local = global_to_local[start_global]
            # directed=True: edge (u -> v) keeps weight l / (d_v**2 + b), i.e.
            # clearance of the node being *entered* - exactly the current heap
            # Dijkstra's cost model. (directed=False would symmetrize each edge 
            # to l / (max(d_u, d_v)**2 + b) and pick slightly different paths)
            distances, predecessors = dijkstra(
                cavity_graph, directed=True, indices=start_local,
                return_predecessors=True)
            parent_to_children = defaultdict(list)

            for node, parent in enumerate(predecessors):
                if parent >= 0:
                    parent_to_children[parent].append(node)

            paths = {}
            stack = [(start_local, [start_local])]
            while stack:
                node, path = stack.pop()
                paths[node] = path
                for child in parent_to_children.get(node, []):
                    stack.append((child, path + [child]))

            # A channel ends where it first touches the surface, i.e. at whichever
            # mouth absorbed it - so emit a candidate for every *reachable* mouth,
            # not for a pre-sampled subset of them. Sampling the targets before the
            # search (the old `end_tetrahedra`, thinned by `sparsity`) can pick a
            # target that sits behind another mouth: the path is absorbed at that
            # nearer mouth and can go no further, the sampled target is never
            # reached, and because only sampled targets emit channels the tunnel is
            # reported nowhere at all. Which mouth happens to shadow which target is
            # a tessellation accident, so real tunnels vanished at some meshes and
            # not others. Every mouth is a legitimate terminus, so let every
            # reachable one produce a candidate and leave exit identity to the
            # dedup, where `sparsity` merges the mouths that share one opening.
            for exit_local in terminals_local:
                if exit_local == start_local:
                    continue
                if np.isinf(distances[exit_local]):
                    continue

                path_local = paths.get(exit_local)
                if path_local is None:
                    continue

                path_global = cavity_tetra[path_local]
                channel = Channel(path_global, *self.processChannel(
                    path_global, vertices, points, vdw_radii, simplices),
                    cost=float(distances[path_local[-1]]))
                # the Dijkstra cost at every node, so that a path cut short at a
                # reported exit can be re-costed at the node it was cut at
                node_costs = np.asarray(distances)[np.asarray(path_local)]
                candidates.append((channel, node_costs))

        if truncate_at_surface:
            self._addDedupedChannels(cavity, candidates, similarity, vertices,
                                     points, vdw_radii, simplices)
        else:
            for channel, _costs in candidates:
                cavity.addChannel(channel)

    def _addDedupedChannels(self, cavity, candidates, similarity, vertices,
                            points, vdw_radii, simplices):
        # Cheapest first, and each candidate is judged only against the channels
        # already kept - so the kept channel is always the cheapest of its group
        # and the result does not depend on the order candidates arrive in.
        #
        # Step 1, CUT. A reported channel's exit sphere (centred on its exit
        # vertex, radius the clearance there) is the volume of that opening. If a
        # candidate's route passes through it, the candidate has left the protein
        # at that opening: whatever it does afterwards is a hop across the outside,
        # not part of a tunnel. So cut it there. This is what stops a path from
        # slipping past a mouth through a twin tetrahedron and claiming some far
        # exit on the other side. Note the cut is made only against the handful of
        # *already reported* exits, never against all mouths - truncating against
        # every mouth is what used to demolish real tunnels.
        #
        # Step 2, COMPARE. Two channels are the same tunnel when they leave by the
        # same opening AND take the same corridor to get there. Both halves are
        # needed: a corridor that forks near the surface and exits twice through
        # one opening is one tunnel counted twice (merge), but two genuinely
        # different corridors that happen to surface at the same opening are two
        # tunnels (keep), and one corridor reaching two separate openings is also
        # two tunnels (keep). A candidate that was cut in step 1 is compared on its
        # cut route, which is the only part of it that is really a tunnel.
        #
        # Route identity is measured GEOMETRICALLY - how much of one centerline
        # runs alongside the other - not as a shared prefix of tetrahedra. A prefix
        # is the wrong instrument twice over: it is blind to rejoining (two routes
        # that split near the seed and then run together to the same exit share
        # almost no prefix, yet are plainly one tunnel) and it is fooled by
        # containment (a short route that is a prefix of a long one scores ~1.0 and
        # deletes the long one, which is the channel carrying the distinctive
        # route). Comparing the curves is immune to both, and to the tessellation:
        # a tetrahedron count is not mesh-invariant, so the same physical fork
        # scores differently at different max_deviation.
        prepared = sorted(candidates, key=lambda c: c[0].cost)
        if not prepared:
            return

        kept = []  # (channel, pts, exit_xyz, opening_radius)
        for channel, node_costs in prepared:
            tetra = np.asarray(channel.tetrahedra)
            pts = vertices[tetra]

            # step 1: cut at the first reported opening this route enters
            cut, cutter = None, None
            for i in range(1, len(pts)):
                for kxyz, kr in ((k[2], k[3]) for k in kept):
                    if np.linalg.norm(pts[i] - kxyz) < kr:
                        cut = i
                        break
                if cut is not None:
                    cutter = next(k for k in kept
                                  if np.linalg.norm(pts[cut] - k[2]) < k[3])
                    break
            if cut is not None:
                tetra = tetra[:cut + 1]
                pts = pts[:cut + 1]
                channel = Channel(tetra, *self.processChannel(
                    tetra, vertices, points, vdw_radii, simplices),
                    cost=float(node_costs[cut]))

            # step 2: same opening AND same corridor -> the same tunnel.
            # The corridor is compared OUTSIDE the shared opening. Inside it the
            # routes are already through the mouth and merely fanning out across
            # it, and that fan is not evidence of a different corridor: the
            # Voronoi network splays where a tunnel widens into its opening, so
            # sibling paths peel off in the last few Angstrom and end on
            # neighbouring exit tetrahedra. Counting that splay as divergence is
            # what used to report one tunnel as a bundle of near-copies.
            duplicate = False
            for _kc, kpts, kxyz, kr in kept:
                if np.linalg.norm(pts[-1] - kxyz) >= kr:
                    continue                        # a different opening
                if self._routeCoverage(pts, kpts, center=kxyz,
                                       radius=kr) >= similarity:
                    duplicate = True
                    break
            if not duplicate:
                if cut is not None:
                    # A cut channel stops inside an opening that is already
                    # reported, so it INHERITS that opening rather than declaring
                    # its own. Its last tetrahedron is an interior one that merely
                    # happens to lie in the exit volume, and its inscribed sphere
                    # is not a mouth - promoting it to a cutting surface would let
                    # an interior sphere start truncating other candidates. (It
                    # survives to here only when it reached that opening by a
                    # genuinely different corridor, which is a distinct tunnel and
                    # must be kept. Note its cost, taken at the cut node, is
                    # necessarily below that of the channel that cut it, since the
                    # cut lies upstream of that channel's mouth - so cost orders
                    # the output but does not mean the cut channel is "better".)
                    kept.append((channel, pts, cutter[2], cutter[3]))
                else:
                    # One radius stands for this opening everywhere: it cuts routes
                    # that pass through it, it decides which channels share it, and
                    # it is the region discounted when their corridors are compared.
                    # The clearance at the exit vertex measures the mouth, but on a
                    # coarse tessellation it is erratic and can collapse to almost
                    # nothing, fragmenting one physical mouth into several; the
                    # sparsity floor keeps it mesh-independent.
                    kept.append((channel, pts, pts[-1],
                                 max(self.calculateMaxRadius(
                                     pts[-1], points, vdw_radii,
                                     simplices[tetra[-1]]), self.sparsity)))
        for channel, _pts, _xyz, _r in kept:
            cavity.addChannel(channel)

    def _routeCoverage(self, a, b, tol=None, center=None, radius=0.0):
        """Fraction of the SHORTER centerline's arc length that runs within ``tol``
        Angstrom of the longer one.

        Answers "does the longer channel follow the shorter one's corridor?".
        ``1.0`` means the shorter route lies wholly inside the longer one's
        corridor, so they took the same way out - the longer one simply carried on
        past the point where the shorter one surfaced. That continuation is *not*
        counted as a difference, which is the point: two channels leaving through
        one opening are one tunnel even if one of them runs on and exits a few
        Angstrom further along. It is safe to ignore the continuation only because
        this is gated on the two channels sharing an opening; without that gate,
        scoring against the shorter route would delete long channels that head off
        to a quite different exit.

        ``center`` and ``radius`` describe that shared opening, and the part of
        either route lying inside it is discarded before the comparison. A tunnel
        splays as it widens into its mouth, so sibling paths peel apart over the
        last few Angstrom and land on neighbouring exit tetrahedra; that fan says
        nothing about which corridor they took, and counting it makes one tunnel
        look like several.

        Note this deliberately says nothing about *where* the routes differ, or how
        sharply the uncovered part turns away - only how much of the shorter route
        is shared. Where two corridors genuinely part company, they do so for a
        large fraction of the route, and the score falls."""
        from scipy.spatial import cKDTree

        if tol is None:
            tol = self.route_tolerance
        if center is not None and radius > 0:
            a = a[np.linalg.norm(a - center, axis=1) > radius]
            b = b[np.linalg.norm(b - center, axis=1) > radius]
            if len(a) < 2 or len(b) < 2:
                # Nothing survives outside the opening, so all either route ever
                # did was cross the mouth: there is no corridor to tell apart.
                return 1.0
        if len(a) < 2 or len(b) < 2:
            return 0.0

        def arclen(p):
            return float(np.linalg.norm(np.diff(p, axis=0), axis=1).sum())

        long_p, short_p = (a, b) if arclen(a) >= arclen(b) else (b, a)
        steps = np.linalg.norm(np.diff(short_p, axis=0), axis=1)
        total = steps.sum()
        if total <= 0:
            return 0.0

        # each node carries half of each adjacent segment, so its weight is the
        # arc length it stands for
        weight = np.zeros(len(short_p))
        weight[:-1] += steps / 2.0
        weight[1:] += steps / 2.0

        near = cKDTree(long_p).query(short_p)[0] <= tol
        return float(weight[near].sum() / total)

    def calculateMaxRadius(self, vertice, points, vdw_radii, simp):
        atom_positions = points[simp]
        radii = vdw_radii[simp]
        distances = np.linalg.norm(atom_positions - vertice, axis=1) - radii
        return np.min(distances)

    def _pathGates(self, tetrahedra, voronoi_vertices, points, vdw_radii,
                   simp, vertex_radii):
        # Per-edge gate clearance along the path: the minimum clearance on each
        # shared Delaunay face between consecutive circumcenters (the edge
        # bottleneck radius), read from the cache and recomputed for any edge the
        # map lacks. Each gate is <= the clearance at both its endpoints, so the
        # path minimum is the reported bottleneck and the gates are where the
        # tube pinches for the volume. Length len(tetrahedra) - 1; empty for a
        # single-tetrahedron path.
        n = len(tetrahedra)
        if n < 2:
            return np.empty(0)
        eb = self._edge_bottleneck
        gates = np.empty(n - 1)
        for k in range(n - 1):
            i, j = int(tetrahedra[k]), int(tetrahedra[k + 1])
            key = (i, j) if i < j else (j, i)
            g = eb.get(key) if eb is not None else None
            if g is None:
                shared = np.intersect1d(simp[i], simp[j], assume_unique=True)
                if len(shared) != 3:
                    # not face-adjacent (should not happen on a graph path);
                    # fall back to the tighter of the two endpoints
                    g = float(min(vertex_radii[k], vertex_radii[k + 1]))
                else:
                    g = self._edgeBottleneck(voronoi_vertices[i],
                                             voronoi_vertices[j], shared,
                                             points, vdw_radii)
            gates[k] = g
        return gates

    def calculateRadiusSpline(self, tetrahedra, voronoi_vertices, points,
                              vdw_radii, simp):
        tetrahedra = np.asarray(tetrahedra)
        # The per-vertex clearance is the same min buildSparseGraph already took
        # over every simplex; read it back instead of recomputing it per path.
        if self._vertex_clearance is not None:
            radii = self._vertex_clearance[tetrahedra]
        else:
            vertices = voronoi_vertices[tetrahedra]
            radii = np.array([self.calculateMaxRadius(v, points, vdw_radii, s)
                              for v, s in zip(vertices, simp[tetrahedra])])
        gates = self._pathGates(tetrahedra, voronoi_vertices, points,
                                vdw_radii, simp, radii)
        return radii, gates

    def processChannel(self, tetrahedra, voronoi_vertices, points, vdw_radii, 
                       simp):
        from scipy.interpolate import CubicSpline
        
        centers = voronoi_vertices[tetrahedra]
        radii, gates = self.calculateRadiusSpline(tetrahedra,
                                                  voronoi_vertices,
                                                  points, vdw_radii, simp)
        bottleneck = float(np.min(gates)) if len(gates) else float(np.min(radii))

        t = np.arange(len(centers))
        centerline_spline = CubicSpline(t, centers, bc_type='natural')
        # The tube pinches at the gates, not at the wide circumcenters, so give
        # the radius profile a knot at each gate (midway between its two
        # vertices) carrying the gate clearance. The centerline keeps only the
        # vertex knots; both splines share the same t domain, so the volume
        # integral samples them consistently and the endpoints (hence the cap
        # radii) are unchanged.
        if len(gates):
            knot_t = np.empty(2 * len(centers) - 1)
            knot_t[0::2] = t
            knot_t[1::2] = t[:-1] + 0.5
            knot_r = np.empty_like(knot_t)
            knot_r[0::2] = radii
            knot_r[1::2] = gates
            radius_spline = CubicSpline(knot_t, knot_r, bc_type='natural')
        else:
            radius_spline = CubicSpline(t, radii, bc_type='natural')

        length = self.calculateChannelLength(centerline_spline)
        volume = self.calculateChannelVolume(centerline_spline, radius_spline)
        
        return centerline_spline, radius_spline, length, bottleneck, volume

    def findBiggestTetrahedron(self, tetrahedra, voronoi_vertices, points, 
                               vdw_radii, simp):
        radii = np.array([self.calculateMaxRadius(voronoi_vertices[tetra], points, vdw_radii, simp[tetra]) for tetra in tetrahedra])
        max_radius_index = np.argmax(radii)
        return tetrahedra[max_radius_index]

    def getEndTetrahedra(self, tetrahedra, voronoi_vertices, points, vdw_radii, 
                         simp, sparsity):
        # Greedy sparse sampling of the mouth (exit) tetrahedra: seed with the
        # widest tetrahedron, then repeatedly add the widest tetrahedron that is
        # still at least `sparsity` away from every already-selected one, until
        # none qualify. Vectorized rewrite of the former O(N_exit x M^2) double
        # loop (which called np.linalg.norm once per candidate/selected pair and
        # re-scanned radii via findBiggestTetrahedron every pass):
        #   * the "far enough from all selected" test is exactly "running min
        #     distance to the selected set >= sparsity", so keep one min_dist
        #     array and fold in each new pick with a single vectorized norm;
        #   * inscribed radii are geometry-only, so precompute them once instead
        #     of recomputing find_biggest over the shrinking candidate set.
        # Selection order and argmax first-tie-break match the original, so the
        # returned end tetrahedra are identical.
        tetrahedra = np.asarray(tetrahedra)
        n = len(tetrahedra)
        if n == 0:
            return tetrahedra

        verts = voronoi_vertices[tetrahedra]            # (n, 3) circumcenters
        radii = np.array([
            self.calculateMaxRadius(voronoi_vertices[tetra], points, 
                                    vdw_radii, simp[tetra])
            for tetra in tetrahedra])

        min_dist = np.full(n, np.inf)
        selected = np.zeros(n, dtype=bool)
        order = []

        current = int(np.argmax(radii))             # widest tetrahedron (seed)
        while True:
            order.append(current)
            selected[current] = True
            min_dist = np.minimum(min_dist, np.linalg.norm(verts - verts[current], axis=1))

            feasible = (min_dist >= sparsity) & ~selected   # >= sparsity from every pick
            if not feasible.any():
                break
            # widest feasible tetrahedron; np.argmax breaks ties toward the
            # lowest index, matching the original input-order scan.
            current = int(np.argmax(np.where(feasible, radii, -np.inf)))

        return tetrahedra[order]

    def filterCavities(self, cavities, min_depth):
        return [cavity for cavity in cavities if cavity.depth >= min_depth]

    def filterChannelsByBottleneck(self, cavities, bottleneck):
        for cavity in cavities:
            cavity.channels = [channel for channel in cavity.channels if channel.bottleneck >= bottleneck]
    
    def filterChannelsByVolume(self, cavities, min_volume=None, max_volume=None):
        """Filter channels by volume."""
    
        for cavity in cavities:
            filtered_channels = []
            for channel in cavity.channels:
                if min_volume is not None and channel.volume < min_volume:
                    continue
                if max_volume is not None and channel.volume > max_volume:
                    continue
                filtered_channels.append(channel)
            cavity.channels = filtered_channels

    def filterCavitiesByTetrahedra(self, cavities, min_tetrahedra=None, 
                                   max_tetrahedra=None):
        """Filter cavities by cavity volume."""
    
        filtered = []
        for cavity in cavities:
            n = len(cavity.tetrahedra)
            if min_tetrahedra is not None and n < min_tetrahedra:
                continue
            if max_tetrahedra is not None and n > max_tetrahedra:
                continue
            filtered.append(cavity)
        return filtered

    def calculateTetrahedronVolume(self, a, b, c, d):
        return abs(np.dot(a - d, np.cross(b - d, c - d))) / 6.0

    def calculate_cavity_volumes(self, cavities, simplices, coords):
        """Calculate approximate cavity volumes from Delaunay tetrahedra."""

        for cavity in cavities:
            volume = 0.0
            for tetra in cavity.tetrahedra:
                atom_ids = simplices[tetra]
                a, b, c, d = coords[atom_ids]
                volume += self.calculateTetrahedronVolume(a, b, c, d)
            cavity.volume = volume

    def filterCavitiesByVolume(self, cavities, min_volume=None, max_volume=None):
        """Filter cavities by approximate volume."""

        filtered_cavities = []
        for cavity in cavities:
            if min_volume is not None and cavity.volume < min_volume:
                continue
            if max_volume is not None and cavity.volume > max_volume:
                continue
            filtered_cavities.append(cavity)
        return filtered_cavities

    def saveChannelsToPdb(self, channels, filename, separate=False, num_samples=5):
        # ``channels`` is a flat list, already ordered by cost - that order is
        # the order they are written and numbered here. Each channel is preceded
        # by a REMARK reporting its length, bottleneck radius, curvature and cost.
        filename = str(filename)

        # All channels in a single file, one after another in list (cost) order.
        with open(filename, 'w') as pqr_file:
            atom_index = 1
            for channel_index, channel in enumerate(channels):
                centerline_spline, radius_spline = channel.getSplines()
                samples = len(channel.tetrahedra) * num_samples
                t = np.linspace(centerline_spline.x[0], centerline_spline.x[-1], samples)
                centers = centerline_spline(t)
                radii = radius_spline(t)

                pqr_file.write(self._channelRemark(channel_index, channel))
                pdb_lines = []
                # Each channel gets its own residue number so the channels stay
                # separable at the record level, matching saveCavitiesToPdb.
                for i, (x, y, z, radius) in enumerate(zip(centers[:, 0], centers[:, 1], centers[:, 2], radii), start=atom_index):
                    pdb_lines.append("ATOM  %5d  H   FIL T%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n" % (i, channel_index + 1, x, y, z, 1.00, radius))

                # Bond consecutive samples of THIS channel only, using the global
                # atom serial numbers (start=atom_index). No CONECT spans two
                # channels, so each one is a separate strand in the viewer.
                for i in range(atom_index, atom_index + samples - 1):
                    pdb_lines.append("CONECT%5d%5d\n" % (i, i + 1))

                pqr_file.writelines(pdb_lines)
                pqr_file.write("\n")
                atom_index += samples

        # When separate is set to True also separate PDB/PQR files will be
        # created, one per channel, numbered by the same cost order.
        if separate:
            for channel_index, channel in enumerate(channels):
                # TODO channel suffix is rather long, making it hard to read in pymol, shorten or user defined?
                channel_filename = filename.replace('.pqr', '_chl{0}.pqr'.format(channel_index))
                channel_filename = channel_filename.replace('.pdb', '_chl{0}.pdb'.format(channel_index))

                with open(channel_filename, 'w') as pqr_file:
                    atom_index = 1
                    centerline_spline, radius_spline = channel.getSplines()
                    samples = len(channel.tetrahedra) * num_samples
                    t = np.linspace(centerline_spline.x[0], centerline_spline.x[-1], samples)
                    centers = centerline_spline(t)
                    radii = radius_spline(t)

                    pqr_file.write(self._channelRemark(channel_index, channel))
                    pdb_lines = []
                    for i, (x, y, z, radius) in enumerate(zip(centers[:, 0], centers[:, 1], centers[:, 2], radii), start=atom_index):
                        pdb_lines.append("ATOM  %5d  H   FIL T%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n" % (i, channel_index + 1, x, y, z, 1.00, radius))

                    for i in range(atom_index, atom_index + samples - 1):
                        pdb_lines.append("CONECT%5d%5d\n" % (i, i + 1))

                    pqr_file.writelines(pdb_lines)

    @staticmethod
    def _channelRemark(channel_index, channel):
        """One-line PQR/PDB REMARK with a channel's basic geometry:
        length, bottleneck radius, curvature and Dijkstra cost."""
        curv = 'n/a' if np.isnan(channel.curvature) else "%.3f" % channel.curvature
        cost = 'n/a' if channel.cost is None else "%.4g" % channel.cost
        return ("REMARK   Channel %d  length=%.3f A  bottleneck=%.3f A  "
                "curvature=%s  cost=%s\n" % (
                    channel_index, channel.length, channel.bottleneck, curv, cost))


    def saveCavitiesToPdb(self, cavities, vertices, filename, separate=False):
        """Save surface cavities to a PDB/PQR file as dummy atoms."""

        filename = str(filename)

        with open(filename, 'w') as pqr_file:
            atom_index = 1
            cavity_index = 0

            for cavity in cavities:
                tetrahedra = cavity.tetrahedra
                if tetrahedra is None or len(tetrahedra) == 0:
                    continue

                points = vertices[tetrahedra]
                for x, y, z in points:
                    pqr_file.write("ATOM  %5d  H   FIL T%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"
                        % (atom_index, cavity_index + 1, x, y, z, 1.00, 1.00))
                    atom_index += 1
                cavity_index += 1

        if separate:
            cavity_index = 0

            for cavity in cavities:
                tetrahedra = cavity.tetrahedra
                if tetrahedra is None or len(tetrahedra) == 0:
                    continue

                points = vertices[tetrahedra]
                cavity_filename = filename.replace('.pqr', '_cavity{0}.pqr'.format(cavity_index))
                cavity_filename = cavity_filename.replace('.pdb', '_cavity{0}.pdb'.format(cavity_index))

                with open(cavity_filename, 'w') as pqr_file:
                    atom_index = 1

                    for x, y, z in points:
                        pqr_file.write("ATOM  %5d  H   FIL T%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"
                            % (atom_index, cavity_index + 1, x, y, z, 1.00, 1.00))
                        atom_index += 1

                cavity_index += 1
            
    def calculateChannelLength(self, centerline_spline):
        t_values = np.linspace(centerline_spline.x[0], centerline_spline.x[-1], len(centerline_spline.x) * 10)
        points = centerline_spline(t_values)
        diffs = np.diff(points, axis=0)
        lengths = np.linalg.norm(diffs, axis=1)
        return np.sum(lengths)
    
    def calculateChannelVolume(self, centerline_spline, radius_spline):
        # Tube volume  V = \int pi r(t)^2 |x'(t)| dt  evaluated by vectorized
        # composite Simpson instead of an adaptive scipy.quad that called the
        # integrand thousands of times per channel. The centerline/radius are
        # piecewise cubic, so a uniform grid scaled to the number of spline
        # segments (~128 points/segment) reaches ~1e-8 relative accuracy - better
        # than quad's default tolerance - at a fraction of the cost.
        from scipy.integrate import simpson

        t_min = centerline_spline.x[0]
        t_max = centerline_spline.x[-1]
        n_segments = max(1, len(centerline_spline.x) - 1)
        n = n_segments * 128 + 1

        t = np.linspace(t_min, t_max, n)
        r = radius_spline(t)
        speed = np.linalg.norm(centerline_spline(t, 1), axis=1)
        volume = simpson(np.pi * r ** 2 * speed, x=t)

        r_start = radius_spline(t_min)
        r_end = radius_spline(t_max)

        hemisphere_volume_start = (2/3) * np.pi * r_start**3
        hemisphere_volume_end = (2/3) * np.pi * r_end**3

        total_volume = volume + hemisphere_volume_start + hemisphere_volume_end

        return total_volume
            
    def selectSeedTetrahedron(self, cavity, vertices, points, vdw_radii, simp,
                              neighbors, sp, search_radius):
        '''Map `sp` to the seed tetrahedron of one cavity.

        The tetrahedron nearest `sp` (the anchor) is frequently a tight one, and every
        channel of the cavity leaves through it: its inscribed radius then caps all of
        their bottlenecks, and the shared first links show up as one common bottleneck
        at the joint beginning of the bundle. So the anchor only says where to look.
        The seed is the widest (largest inscribed radius) tetrahedron of this cavity
        that lies within `search_radius` of `sp`, is no shallower than the anchor, and
        is reachable from the anchor through the tetrahedra within `search_radius`.
        Reachability is over the adjacency of the cleared tetrahedra, which is free
        space, so the seed can only move through the void the start point sits in and
        never hops across a wall into a lobe that merely passes nearby; the depth floor
        keeps it from sliding outward towards the mouth, where tetrahedra are wide but
        no longer inside the site. Note that the floor filters the seed, not the walk:
        a marginally shallower cell in between must not wall off the wider region
        behind it.

        `search_radius` <= 0 restores the plain nearest-vertex seed.

        :returns: dict of the seed and anchor properties (`seed`, `anchor`, and their
            `_vertex`, `_distance` from `sp`, inscribed `_radius` and `_depth`), plus
            the number of tetrahedra `searched` and how many of them were `eligible`'''

        from collections import deque

        depths = cavity.tetrahedra_depths
        tet = np.asarray(cavity.tetrahedra)
        d2 = np.sum((vertices[tet] - sp) ** 2, axis=1)
        anchor = int(tet[int(np.argmin(d2))])

        def properties(tetra):
            return dict(
                vertex=vertices[tetra],
                distance=float(np.linalg.norm(vertices[tetra] - sp)),
                radius=float(self.calculateMaxRadius(
                    vertices[tetra], points, vdw_radii, simp[tetra])),
                depth=float(depths.get(tetra, 0.0)))

        def report(seed, searched, eligible):
            info = {'seed': seed, 'anchor': anchor,
                    'searched': searched, 'eligible': eligible}
            for name, tetra in (('seed', seed), ('anchor', anchor)):
                for key, value in properties(tetra).items():
                    info['{0}_{1}'.format(name, key)] = value
            return info

        if not search_radius or search_radius <= 0:
            return report(anchor, 1, 1)

        near = set(int(t) for t, close in zip(tet, d2 <= search_radius ** 2) if close)

        # BFS from the anchor, staying inside the sphere and inside the cavity.
        reachable = [anchor]
        seen = {anchor}
        queue = deque([anchor])
        while queue:
            current = queue.popleft()
            for neighbor in neighbors[current]:
                neighbor = int(neighbor)
                if neighbor in near and neighbor not in seen:
                    seen.add(neighbor)
                    reachable.append(neighbor)
                    queue.append(neighbor)

        anchor_depth = depths.get(anchor, 0)
        eligible = [t for t in reachable if depths.get(t, 0) >= anchor_depth]

        reach = np.array(eligible, dtype=np.intp)
        atoms = simp[reach]
        clearance = (np.linalg.norm(points[atoms] - vertices[reach][:, None, :], axis=2)
                     - vdw_radii[atoms])
        radii = clearance.min(axis=1)
        # The anchor is eligible and comes first (BFS order), so argmax ties to it.
        best = int(np.argmax(radii))

        return report(int(reach[best]), len(reachable), len(eligible))

    def setStartingTetrahedraFromPoint(self, cavities, vertices, start_point,
                                       points, vdw_radii, simp, neighbors,
                                       restrict=False, search_radius=5.0):
        '''Set starting tetrahedra using a user-defined 3D point.
        The starting tetrahedron of a cavity is the widest one `selectSeedTetrahedron`
        finds in the neighbourhood of `start_point`; with ``search_radius=0`` it is
        simply the one whose Voronoi vertex is closest to `start_point`.

        :arg cavities: list of cavity objects
        :arg vertices: Voronoi vertices (array of shape (n, 3))
        :arg start_point: point [x, y, z] in Å (list/tuple/ndarray of length 3)
        :arg points: atom coordinates (array of shape (n_atoms, 3)), used to compute
            the inscribed radius of candidate tetrahedra
        :arg vdw_radii: per-atom van der Waals radii (array of shape (n_atoms,))
        :arg simp: simplices (tetrahedron -> its 4 atom indices)
        :arg neighbors: tetrahedron adjacency (tetrahedron -> its 4 neighbours, -1 none)
        :arg restrict: if True, only the single cavity whose closest tetrahedron is
            globally nearest to `start_point` is seeded and returned, so channels are
            computed exclusively for the region around `start_point`. If False (default),
            every cavity is seeded with its own seed tetrahedron and all cavities are
            returned unchanged.
        :type restrict: bool
        :arg search_radius: radius, in Angstrom, of the neighbourhood of `start_point`
            searched for a wider seed. 0 disables the search.
        :type search_radius: float
        :returns: list of cavities to search: all cavities when `restrict` is False, the
            single selected cavity when `restrict` is True, or an empty list if no cavity
            has any tetrahedra'''

        sp = np.asarray(start_point, dtype=float).reshape(3,)

        best_cavity = None
        best_info = None

        for i, cavity in enumerate(cavities):
            tet = cavity.tetrahedra
            if tet is None or len(tet) == 0:
                continue

            info = self.selectSeedTetrahedron(
                cavity, vertices, points, vdw_radii, simp, neighbors, sp, search_radius)

            if not restrict:
                cavity.setStartingTetrahedron(np.array([info['seed']]))
                self.reportSeedTetrahedron(info, search_radius, cavity_index=i)

            # The cavity is still chosen by proximity to start_point: widening moves the
            # seed inside a cavity, it must never decide between cavities.
            if best_info is None or info['anchor_distance'] < best_info['anchor_distance']:
                best_info = info
                best_cavity = cavity

        if not restrict:
            return cavities

        if best_cavity is None:
            _warn("start_point was provided but no cavity contains any "
                "tetrahedron; no channels will be computed.")
            return []

        best_cavity.setStartingTetrahedron(np.array([best_info['seed']]))
        self.reportSeedTetrahedron(best_info, search_radius)
        LOGGER.info("    restricting the channel search to the cavity that contains it "
            "({0} tetrahedra, depth {1:.1f} A).".format(len(best_cavity.tetrahedra),
                                                        float(best_cavity.depth)))

        return [best_cavity]

    def reportSeedTetrahedron(self, info, search_radius, cavity_index=None):
        '''Log the seed tetrahedron `selectSeedTetrahedron` picked, and, when it is not
        the one nearest the start point, the anchor it replaced -- the two radii are what
        tell the user whether the seed was capping the bottlenecks of the cavity.'''

        where = '' if cavity_index is None else ' of cavity {0}'.format(cavity_index)
        LOGGER.info("start_point seeded at tetrahedron {0}{1} (Voronoi vertex at "
            "[{2:.3f}, {3:.3f}, {4:.3f}], {5:.3f} A from start_point, inscribed radius "
            "{6:.3f} A, depth {7:.1f} A)."
            .format(info['seed'], where, info['seed_vertex'][0], info['seed_vertex'][1],
                    info['seed_vertex'][2], info['seed_distance'], info['seed_radius'],
                    info['seed_depth']))

        if info['seed'] != info['anchor']:
            LOGGER.info("    widened from the nearest tetrahedron {0} ({1:.3f} A away, "
                "inscribed radius {2:.3f} A, depth {3:.1f} A), the widest of the {4} tetrahedra "
                "no shallower than it among the {5} reachable within {6:.1f} A; seeding "
                "the narrow one would have capped every channel here at its radius."
                .format(info['anchor'], info['anchor_distance'], info['anchor_radius'],
                        info['anchor_depth'], info['eligible'], info['searched'],
                        float(search_radius)))
        elif search_radius and search_radius > 0:
            LOGGER.info("    already the widest of the {0} tetrahedra no shallower than "
                "it among the {1} reachable within {2:.1f} A."
                .format(info['eligible'], info['searched'], float(search_radius)))


    def trimCavitiesByDepth(self, cavities, max_depth):
        """Filtering cavities by max_depth."""
    
        for cavity in cavities:
            cavity.tetrahedra = np.array([
                tetra for tetra in cavity.tetrahedra
                if cavity.tetrahedra_depths.get(tetra, np.inf) <= max_depth])

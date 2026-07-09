# -*- coding: utf-8 -*-

"""This module is called CaviFinder and defines functions for calculating 
channels, tunnels and pores within protein structure.
"""

__author__ = 'Karolina Mikulska-Ruminska', 'Eryk Trzcinski'
__credits__ = ['Karolina Mikulska-Ruminska', 'Eryk Trzcinski']
__email__ = ['karolamik@fizyka.umk.pl']

import numpy as np
from prody import LOGGER, PY3K
from prody.atomic import Atomic
from prody.utilities import checkCoords, getCoords, isListLike
from prody.proteins import writePDB, parsePDB, parsePQR
from prody.ensemble import Ensemble
from prody.measure import calcCenter

from .fixer import *
from .compare import *
from prody.measure import calcTransformation, calcDistance, calcRMSD, superpose


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
           'getChannelParametersMultipleFrames',
           'getChannelResidueNamesMultipleFrames']


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
            LOGGER.warn("Package " + str(package_name) + " is not installed. "
            "Please install it to use this function.")
            return False
    else:
        try:
            __import__(package_name)
        except ImportError:
            LOGGER.warn("Package " + str(package_name) + " is not installed. "
            "Please install it to use this function.")
            return False
    
    return True


def getVmdModel(vmd_path, atoms, representation='NewCartoon'):
    """Generates a 3D model of molecular structures using VMD and returns it as
      an Open3D TriangleMesh.

    This function creates a temporary PDB file from the provided atomic data a
    nd uses VMD (Visual Molecular Dynamics) to render this data into an STL 
    file, which is then loaded into Open3D as a TriangleMesh. The function 
    handles the creation and cleanup of temporary files and manages the 
    subprocess call to VMD.
    
    To install Open3D use: 
    conda install open3d (for Anaconda users; version open3d-0.19.0 was used 
    during the developement) or pip install open3d

    :arg vmd_path: Path to the VMD executable. This is required to run VMD and 
    execute the TCL script.
    :type vmd_path: str

    :arg atoms: Atomic data to be written to a PDB file. This should be an 
    object or data structure
        that is compatible with the `writePDB` function.
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
        LOGGER.warn("An unexpected error occurred: " + str(e))
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
    during the developement) or pip install open3d
    
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

        spheres = [o3d.geometry.TriangleMesh.create_sphere(radius=r, 
                                                           resolution=20).translate(c) for r, c in zip(radii, centers)]
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


def showCavities(surface, show_surface=False):
    """Visualizes the cavities within a molecular surface using Open3D.

    This function displays a 3D visualization of cavities detected in a 
    molecular structure.
    It uses the Open3D library to render the cavities as a triangle mesh. 
    Optionally, it can also display the molecular surface as a wireframe 
    overlay.

    To install Open3D use: 
    conda install open3d (for Anaconda users; version open3d-0.19.0 was used 
    during the developement) or pip install open3d

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
      is raised,
        prompting the user to install Open3D.

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
                                          eturn_counts=True)[0]
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
    restrict_channels_to_start_point=False, r1=3, r2=0.9, min_depth=10, 
    min_volume=None, max_volume=None, max_depth=None, bottleneck=0.9, 
    sparsity=1, min_tetrahedra=None, max_tetrahedra=None, cavities_only=False, 
    diagram="homogenized", max_deviation=0.1, truncate_at_surface=True, 
    similarity=0.8, max_peel_depth=None):
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

    The implementation is inspired by the methods described in the publication:
    "MOLE 2.0: advanced approach for analysis of biomacromolecular channels" by
     D. Sehnal, et al., published in J Chemoinform, 5 (39) 2013.

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
        3D coordinate point will be provided, the algorithm will use the 
        tetrahedron whose Voronoi vertex is closest to this point as the 
        starting tetrahedron (overriding the default automatic seed selection 
        based on the deepest tetrahedron). Coordinates must be given in Å.
        If an atomic selection is provided, its geometric center is used as the
         starting point.
    :type start_point: list, tuple, or ndarray (length 3), :class:`.Atomic`, or None

    :arg restrict_channels_to_start_point: Only used when ``start_point`` is 
        provided. If True, the channel search is restricted to the single cavity
         whose closest tetrahedron is globally nearest to ``start_point``, so 
        channels are computed only for the region around that point instead of 
        one channel bundle per detected cavity. If False (default), 
        ``start_point`` merely overrides the seed (starting) tetrahedron of 
        every cavity and channels are still computed for all cavities.
    :type restrict_channels_to_start_point: bool 

    :arg r1: The first radius threshold used during the deletion of simplices, 
        which is used to define the outer surface of the channels. Default is 3
    :type r1: float

    :arg r2: The second radius threshold used to define the inner surface of 
        the channels. Default is 0.9.
    :type r2: float

    :arg min_depth: The minimum depth a cavity must have to be considered as a 
        channel. Default is 10.
    :type min_depth: int

    :arg max_depth: Maximum cavity depth. Cavities deeper than this value are 
        trimmed to the specified depth. Default is None.
    :type max_depth: int

    :arg bottleneck: The minimum allowed bottleneck size (narrowest point) for 
        the channels. Default is 0.9.
    :type bottleneck: float

    :arg min_volume: Minimum volume required for a channel/cavity to be 
        retained. Default is None.
    :type min_volume: float

    :arg max_volume: Maximum volume allowed for a channel/cavity to be 
        retained. Default is None.
    :type max_volume: float

    :arg sparsity: The sparsity parameter controls the sampling density when 
        analyzing the molecular surface. A higher value results in fewer 
        sampling points. Default is 1, which enables detection of most relevant 
         channel branches.
    :type sparsity: int
    
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
        "weighted" - TODO

    :type diagram: str

    :arg max_deviation: Maximum tolerated deviation, in Angstrom, between the 
        union surface of the substitute balls and the original van der Waals 
        surface when ``diagram = homogenized`` . It controls the trade-off
        between surface accuracy and the number of balls generated: an atom 
        whose radius exceeds the smallest radius (``rho``) by more than 
        ``max_deviation`` is filled with several balls, otherwise it is kept as
         a single ``rho`` ball. Default is 0.1. Guideline values:

        * ``0.1`` fine accurate surface with minimal errors, but on average 15 
            times more balls than original
        * ``0.15`` in heavy-atom-only structures it startsfilling carbon, which
             is otherwise left as a single ball with a uniform ~0.18 A inset).
        * ``0.2`` speed optimized ; e.g. carbon fills to ~15 balls when 
            hydrogens are present (``rho``=1.2), resulting roughly to ~8 times 
            more balls. Without hydrogens, carbons are single balls.

        Only used when ``diagram = homogenized``.
    :type max_deviation: float

    :arg truncate_at_surface: If True (default), each channel is terminated at 
        the first surface (exit)tetrahedron it reaches whose inscribed radius 
        is at least ``bottleneck``, instead of running all the way to its 
        assigned end tetrahedron. This prevents a cheapest path from surfacing 
        at one mouth and continuing on to another, and de-duplicates the 
        channels that collapse onto a shared mouth (keeping the cheapest per 
        terminal). If False, the original behaviour is kept (paths run to the 
        end tetrahedra).
    :type truncate_at_surface: bool

    :arg similarity: Only used when ``truncate_at_surface`` is True. Fraction 
        (0-1) of the shorter path that two channels must share, as a common 
        prefix from the seed, to be treated as the same tunnel when they leave
         through the same surface opening. Two channels are merged (cheapest 
         kept) only if their exit points coincide (the mouth spheres they leave
         through overlap) AND their shared-prefix fraction is at least 
         ``similarity``; channels that exit at distinct mouths, or reach one 
        exit by genuinely different corridors (diverging early, sharing a bit),
         are kept as separate tunnels. ``1.0`` merges only paths that share an 
         exit and are otherwise identical; ``0.0`` keeps one channel per 
         distinct exit. Default is 0.8.
    :type similarity: float

    :arg max_peel_depth: Safety cap on the bounded surface peel. After the r1 
        surface is built, it is eroded inward with the r2 probe by 
        ``round(r1 - r2)`` layers to strip the wide former-exterior shell (the
        "moat") that a large r1 probe bridges over; that shell would otherwise 
        act as a low-cost path on which channels truncate and collapse. The 
        peel is near-inert at the default ``r1``/``r2`` and grows with the
        gap ``r1 - r2``. ``max_peel_depth`` limits the number of eroded layers,
         as a guard against over-peeling into the interior at large ``r1``; 
         ``None`` (default) leaves the peel uncapped. Erosion also stops early
        on its own once no boundary tetrahedron wider than ``r2`` remains.
    :type max_peel_depth: int or None

    :returns: A tuple containing two elements:
        - `channels`: A list of detected channels, where each channel is an 
          object containing informationabout its path and geometry.
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
                                     separate=False, r1=3, r2=1.25, min_depth=10, 
                                     bottleneck=1, sparsity=15) """
    
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

    from scipy.spatial import Delaunay
    
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

    calculator = ChannelCalculator(atoms, r1, r2, min_depth, bottleneck, sparsity)

    if diagram not in ["homogenized", "weighted"]:
        atoms = atoms.select('not hetero and noh') # Excluding hydrogens
        # TODO in fact we should perhaps do the filtering outside, as you might
        #  want heteroatoms too, e.g., HEM in CYPs 

    coords = atoms.getCoords()
    vdw_radii = calculator.getVdwRadii(atoms.getElements())

    if diagram == "homogenized":
        LOGGER.timeit('_prody_channels_homogenize')
        coords, vdw_radii = calculator.homogenizeAtoms(coords, vdw_radii, max_deviation)
        LOGGER.report("Substituted {0} atoms with {1} homogeneous balls of radius {2:.2f} A in %.2fs.".format(
            atoms.numAtoms(), len(coords), float(vdw_radii[0])), '_prody_channels_homogenize')

    if diagram == "weighted":
        #TODO using vorpy3 package?
        pass

    LOGGER.timeit('_prody_channels_tessellation')
    dela = Delaunay(coords)
    # circumcenters straight from the Delaunay paraboloid lifting, so we
    # skip the redundant second Qhull pass (scipy Voronoi). Numerically identical
    # to voro.vertices for points in general position.
    verts = calculator.calcCircumcenters(dela)
    LOGGER.report('Delaunay tessellation of {0} points constructed in %.2fs.'.format(
        len(coords)), '_prody_channels_tessellation')

    LOGGER.timeit('_prody_channels_surface')
    s_prt = State(dela.simplices, dela.neighbors, verts)
    
    if PY3K:
        s_tmp = State(*s_prt.getState())
        s_prv = State(None, None, None)
    else:
        s_tmp = apply(State, s_prt.getState())
        s_prv = State(None, None, None) 
        
    while True:
        s_prv.setState(*s_tmp.getState())
        
        if PY3K:
            #s_tmp.setState(*calculator.deleteSimplices3d(coords, *(s_tmp.getState() + [vdw_radii, r1, True])))
            s_tmp.setState(*calculator.deleteSimplices3d(coords, *(s_tmp.getState() + tuple([vdw_radii, r1, True]))))
        else:
            tmp_state = calculator.deleteSimplices3d(coords, *(s_tmp.getState() + [vdw_radii, r1, True]))
            s_tmp.setState(*tmp_state)

        if s_tmp == s_prv:
            break
        
    s_srf = State(*s_tmp.getState())

    # Bounded r2 peel (moat removal): erode the r1 surface inward with the r2 probe
    # by round(r1 - r2) layers, stripping the wide former-exterior "moat" shell that
    # a large r1 bridges over (it would otherwise act as a low-cost path that truncates
    # channels). Stops early once erosion converges. max_peel_depth caps it (None = uncapped).
    peel_depth = int(round(r1 - r2))
    if max_peel_depth is not None:
        peel_depth = min(peel_depth, max_peel_depth)
    for _ in range(peel_depth):
        s_next = State(*calculator.deleteSimplices3d(coords, *(s_srf.getState() + tuple([vdw_radii, r2, True]))))
        if s_next == s_srf:
            break
        s_srf = s_next

    #s_inr = State(*calculator.deleteSimplices3d(coords, *(s_srf.getState() + [vdw_radii, r2, False])))
    s_inr = State(*calculator.deleteSimplices3d(coords, *(s_srf.getState() + tuple([vdw_radii, r2, False]))))

    l_first_layer_simp, l_second_layer_simp = calculator.surfaceLayer(s_srf.simp, s_inr.simp, s_srf.neigh)
    s_clr = State(*calculator.deleteSection(l_first_layer_simp, *s_inr.getState()))
    LOGGER.report('Surface and inner simplices filtered in %.2fs.', '_prody_channels_surface')

    LOGGER.timeit('_prody_channels_cavities')
    c_cavities = calculator.findGroups(s_clr.neigh)
    c_surface_cavities = calculator.getSurfaceCavities(c_cavities, s_clr.simp, 
                                                       l_second_layer_simp, 
                                                       s_clr, coords, 
                                                       vdw_radii, sparsity)

    calculator.findDeepestTetrahedra(c_surface_cavities, s_clr.neigh)
    if start_point is not None:
        c_surface_cavities = calculator.setStartingTetrahedraFromPoint(
            c_surface_cavities, s_clr.verti, start_point, restrict_channels_to_start_point)

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
    return channels, [coords, s_srf.simp, merged_cavities, s_clr.simp]

            
def calcChannelsMultipleFrames(atoms, trajectory=None, output_path=None, separate=False, start_point=None, **kwargs):
    """Compute channels for each frame in a given trajectory or multi-model PDB
      file.

    This function calculates the channels for each frame in a trajectory or for
     each model in a multi-model PDB file. The `kwargs` can include parameters 
     necessary for channel calculation. If the `separate` parameter is set to 
     True, each detected channel will be saved in a separate PDB file.

    :arg atoms: Atomic data or object containing atomic coordinates and methods for accessing them.
    :type atoms: object

    :arg trajectory: Trajectory object containing multiple frames or a 
        multi-model PDB file.
    :type trajectory: Atomic or Ensemble object

    :arg output_path: Optional path to save the resulting channels and 
        associated data in PDB format. If a directory is specified, each 
        frame/model will have its results saved in separate files. If None, 
        results are not saved. Default is None.
    :type output_path: str or None

    :arg separate: If True, each detected channel is saved to a separate PDB file for each frame/model.
        If False, all channels for each frame/model are saved in a single file. Default is False.
    :type separate: bool

    :arg start_point: Optional starting point for channel search. If provided, the algorithm will use 
        the tetrahedron whose Voronoi vertex is closest to this point as the starting tetrahedron (overriding 
        the default automatic seed selection based on the deepest tetrahedron). Coordinates must be given in Å.
    :type start_point: list, tuple, or ndarray (length 3), or None 

    :arg kwargs: Additional parameters required for channel calculation. This can include parameters such as
        radius values (r1, r2), minimum depth (min_depth), bottleneck values, etc. 
        See the available parameters in calcChannels().
    :type kwargs: dict

    :returns: List of channels and surfaces computed for each frame or model. Each entry in the list corresponds
        to a specific frame or model.
    :rtype: list of lists

    Example usage:
    channels_all, surfaces_all = calcChannelsMultipleFrames(atoms, trajectory=traj, output_path="channels.pdb", 
                                   separate=False, r1=3, r2=1.25, min_depth=10, bottleneck=1, sparsity=15) 
                                  
    channels_all, surfaces_all = calcChannelsMultipleFrames(atoms, trajectory=traj, output_path="channels.pdb", 
                                   separate=False, start_point=[-10.353, -0.133, 5.608]) """

    
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


def calcSurfaceCavitiesMultipleFrames(atoms, trajectory=None, output_path=None, separate=False, **kwargs):
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
        `sparsity`, `start_frame`, and `stop_frame`.
    :type kwargs: dict

    :returns: Two lists:
        - `cavities_all`: a list containing detected surface cavities for each
          analyzed frame/model,
        - `surfaces_all`: a list containing the corresponding surface representations 
          for each analyzed frame/model.
    :rtype: tuple (list, list)

    Example usage:
    protein = parsePDB('1tqn').select('protein')
    cavities_all, surfaces_all = calcSurfaceCavitiesMultipleFrames(protein, trajectory=traj, output_path="surface_cavities",
        r1=4.5, r2=2.0, min_depth=2, max_depth=3, min_volume=50)

    cavities_all, surfaces_all = calcSurfaceCavitiesMultipleFrames(protein, start_frame=0, stop_frame=10, r1=4.5, r2=2.0) """

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

    This functaaion iterates through a list of channel objects, extracting the 
    length, bottleneck,and volume of each channel. These values are collected 
    into separate lists, which are returned as a tuple for further use.

    :arg channels: A list of channel objects, where each channel has attributes
      `length`, `bottleneck`,and `volume`. These attributes represent the 
      length of the channel, the minimum radius (bottleneck) along its path, 
      and the total volume of the channel, respectively.
    :type channels: list

    :arg param_file_name: The files with parameters will be saved in a text 
        file with the provided name.Use one word which will be added to 
        '_Parameters_All_channels.txt' sufix. If further analysis will be 
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
            lines.append("{0}_cavity{1}: {2:.3f} {3} {4}\n".format(
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
         of the channel.More samples result in a finer representation of the 
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

    This function should be used with results returned by :func:`calcSurfaceCavitiesMultipleFrames`.

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
    certain distance (distA) from selected residue (temporarly one residue).
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
        channel (FIL atoms)default is 5
    :type distA: int, float 
        
    :arg residues_file: File with residues forming the channel created by 
        getChannelResidues(), default is False 
    :type residues_file: bool

    :arg param_file: File with residues forming the channel created by getChannelParameters()
                     default is False
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

    # Extract paramaters and/or residues with channel selection
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

    :arg output_file_name: The name of the PDB file with overlapping surfaces.
    :type output_file_name: str

    :arg pqr_files: File with residues forming the channel created by 
        getChannelResidues() default is False (then all the files from the 
        current directory will be analyzed)when providing a list, only the PDBs
         from list will be analyzedwhen providing str, it will be treated as a 
        folder path  
    :type pqr_files: bool, list or str
    
    Example usage:
    calcChannelSurfaceOverlaps() - all the files in the current directory will 
    be analyzed
    
    calcChannelSurfaceOverlaps(pqr_files='./DATA', output_file_name='results.pdb') - only files from 
    the DATA folder will be analyzed and results will be saved as results.pdb
    
    list_of_files = ['file1.pqr', 'file2.pqr', 'file3.pqr', ..]
    calcChannelSurfaceOverlaps(pqr_files=list_of_files, output_file_name='results.pdb') - files from 
    the list will be analyzed and results will be saved as results.pdb
    """
    
    import os

    resolution = kwargs.pop('resolution', 0.5)
     
    pqr_files = kwargs.pop('pqr_files', False)
    if pqr_files == False or pqr_files is None:
        # take all PQRs from the current dir
        pqr_files = [file for file in os.listdir('.') if file.endswith('.pqr')]
    elif isinstance(pqr_files, str):
        # folder path
        pqr_files = [file for file in os.listdir(pqr_files) if file.endswith('.pqr')]
    elif isinstance(pqr_files, list):
        # list of PQRs
        pqr_files = [file for file in pqr_files if file.endswith('.pqr')]
    else:
        raise ValueError('Please provide list with PQR files, folder path, or nothing to analyze PQRs in the current folder')

    output_file_name = kwargs.pop('output_file_name','overlap_regions.pdb')
    if os.path.exists(output_file_name):
        os.rename(output_file_name, output_file_name+'-old')

    def loadPDBdata(filepath):
        """Parse a PQR file and return a list of atom dictionaries for lines containing 'FIL'."""
        atoms_set = []
        FILatoms = parsePQR(filepath).select('resname FIL')
        
        if FILatoms == None:
            pass
        else:
            for nr_i, i in enumerate(FILatoms):
                FILatoms_coords = FILatoms.getCoords()[nr_i]
                FILBetas_value = FILatoms.getRadii()[nr_i]
                atoms_set.append({
                    'x': float(FILatoms_coords[0]),
                    'y': float(FILatoms_coords[1]),
                    'z': float(FILatoms_coords[2]),
                    'radius': float(FILBetas_value)
                })
        return atoms_set

    def create_surface(atoms, resolution=resolution):
        """Create a 3D grid representing the surface occupied by the atoms."""
        surface = {}
        Zr = 0
        for atom in atoms:
            x, y, z, radius = atom['x'], atom['y'], atom['z'], atom['radius']
            for i in np.arange(x - radius, x + radius, resolution):
                for j in np.arange(y - radius, y + radius, resolution):
                    for k in np.arange(z - radius, z + radius, resolution):
                        if (i - x) ** 2 + (j - y) ** 2 + (k - z) ** 2 <= radius ** 2:
                            key = (round(i, Zr), round(j, Zr), round(k, Zr))
                            surface[key] = surface.get(key, 0) + 1
        return surface

    def merge_surfaces(surfaces):
        """Merge multiple surfaces and calculate overlap counts."""
        merged_surface = {}
        for surface in surfaces:
            for key in surface:
                merged_surface[key] = merged_surface.get(key, 0) + 1
        return merged_surface

    def write_merge_surf_pdb(merged_surface, filename, nr_pdbs):
        """Write the merged surface into a PDB file."""
        with open(filename, 'w') as file:
            atom_id = 1
            for (x, y, z), count in merged_surface.items():
                norm_count = count/nr_pdbs
                file.write("ATOM  {:5d}  H   FIL T   1    {:8.3f}{:8.3f}{:8.3f}{:6.2f}  1.00\n".format(atom_id, x, y, z, norm_count))
                atom_id += 1

    surfaces = []
    for nr_pdbs,pqr_file in enumerate(pqr_files):
        LOGGER.info("Processing file: {0}".format(pqr_file))
        atoms = loadPDBdata(pqr_file)
        if atoms:
            surface = create_surface(atoms, resolution=resolution)
            surfaces.append(surface)
    
    nr_pdbs = nr_pdbs+1
    merged_surface = merge_surfaces(surfaces)
    write_merge_surf_pdb(merged_surface, output_file_name, nr_pdbs)


def calcSurfaceCavityOverlaps(**kwargs):
    """Calculate overlapping regions of surface cavities represented as FIL atoms.

    It calculates spatial overlap between surface cavities saved as PQR files with
    FIL pseudoatoms, as generated by :func:`calcSurfaceCavities` or
    :func:`calcSurfaceCavitiesMultipleFrames`.

    Results are normalized within [0, 1], where the value corresponds to the
    fraction of analyzed PQR files contributing to a given spatial region.

    :arg resolution: surface sampling resolution. Default is 0.5.
    :type resolution: float

    :arg output_file_name: name of the output PDB file with overlapping cavity
        regions. Default is ``'surface_cavity_overlap_regions.pdb'``.
    :type output_file_name: str

    :arg pqr_files: PQR files with surface cavities represented as FIL atoms.
        If not provided, all PQR files from the current directory will be analyzed.
        A list of PQR files can also be provided.
    :type pqr_files: bool or list """

    kwargs.setdefault('output_file_name', 'surface_cavity_overlap_regions.pdb')

    return calcChannelSurfaceOverlaps(**kwargs)


def calcSurfaceCavities(atoms, output_path=None, r1=4.5, r2=2.0, min_depth=2, 
                        max_depth=3, min_tetrahedra=None, max_tetrahedra=None, 
                        min_volume=50, max_volume=None, sparsity=15, 
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

    :arg min_depth: The minimum depth a cavity must have to be considered as a 
        cavity. Default is 2.
    :type min_depth: int

    :arg max_depth: Maximum cavity depth. Cavities deeper than this value are 
        trimmed to the specified depth. Default is 3.
    :type max_depth: int

    :arg sparsity: The sparsity parameter controls the sampling density when 
        analyzing the molecular surface.A higher value results in fewer 
        sampling points. Default is 15.
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
            visualization, includingthe atomic coordinates, simplices defining 
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


    cavities, surface = calcChannels(
            atoms, 
            output_path=output_path,
            separate=separate,
            r1=r1, r2=r2, 
            min_depth=min_depth, max_depth=max_depth,
            min_volume=min_volume, max_volume=max_volume,
            min_tetrahedra=min_tetrahedra, max_tetrahedra=max_tetrahedra,
            sparsity=sparsity, cavities_only=True)
    
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
    def __init__(self, atoms, r1=3, r2=1.25, min_depth=10, bottleneck=1, 
                 sparsity=15):
        self.atoms = atoms
        self.r1 = r1
        self.r2 = r2
        self.min_depth = min_depth
        self.bottleneck = bottleneck
        self.sparsity = sparsity
        
    # def sphereFit(self, vertices, tetrahedron, vertice, vdw_radii, r):
    #     center = vertice
    #     d_sum = sum(np.linalg.norm(center - vertices[atom]) for atom in tetrahedron)
    #     r_sum = sum(r + vdw_radii[atom] for atom in tetrahedron)
        
    #     return d_sum >= r_sum

    def deleteSimplices3d(self, points, simplices, neighbors, vertices, 
                          vdw_radii, r, surface):
        simplices = np.asarray(simplices)
        neighbors = np.asarray(neighbors)
        vertices = np.asarray(vertices)

        n = len(simplices)
        if n == 0:
            return simplices, neighbors, vertices

        # Vectorized sphereFit: for each tetrahedron compare the sum of distances
        # from its Voronoi vertex to its 4 atoms against the sum of (r + vdw_radius)
        # over those atoms. In the surface pass only boundary tetrahedra (those with
        # a -1 neighbour) can ever be deleted, so restrict the expensive norm to that
        # shell (~n^(2/3) rows) instead of evaluating it over every tetrahedron on
        # each erosion iteration.
        if surface:
            boundary = (neighbors == -1).any(axis=1)
            should_delete = np.zeros(n, dtype=bool)
            if boundary.any():
                atom_coords = points[simplices[boundary]]               # (m, 4, 3)
                d_sum = np.linalg.norm(
                    atom_coords - vertices[boundary][:, None, :], axis=2).sum(axis=1)
                r_sum = (r + vdw_radii[simplices[boundary]]).sum(axis=1)
                should_delete[boundary] = d_sum >= r_sum
        else:
            atom_coords = points[simplices]                             # (n, 4, 3)
            d_sum = np.linalg.norm(atom_coords - vertices[:, None, :], axis=2).sum(axis=1)
            r_sum = (r + vdw_radii[simplices]).sum(axis=1)
            should_delete = d_sum < r_sum

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
                           state, points, vdw_radii, sparsity):
        surface_cavities = []
        
        for cavity in cavities:
            tetrahedra = cavity.tetrahedra
            second_layer_mask = np.isin(interior_simplices[tetrahedra], second_layer).all(axis=1)
            
            if np.any(second_layer_mask):
                cavity.makeSurface()
                exit_tetrahedra = tetrahedra[second_layer_mask]
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

    def findDeepestTetrahedra(self, cavities, neighbors):
        from collections import deque
        
        for cavity in cavities:
            exit_tetrahedra = cavity.exit_tetrahedra
            # O(1) membership instead of scanning the tetrahedra array per edge.
            cavity_tetra_set = set(cavity.tetrahedra.tolist())
            visited = np.zeros(neighbors.shape[0], dtype=bool)
            visited[exit_tetrahedra] = True
            queue = deque([(tetra, 0) for tetra in exit_tetrahedra])
            max_depth = -1
            deepest_tetrahedron = None
            tetrahedra_depths = {}

            while queue:
                current, depth = queue.popleft()
                tetrahedra_depths[current] = depth

                if depth > max_depth:
                    max_depth = depth
                    deepest_tetrahedron = current

                for neighbor in neighbors[current]:
                    if neighbor != -1 and not visited[neighbor] and neighbor in cavity_tetra_set:
                        visited[neighbor] = True
                        queue.append((neighbor, depth + 1))

            cavity.setStartingTetrahedron(np.array([deepest_tetrahedron]))
            cavity.setDepth(max_depth)
            cavity.tetrahedra_depths = tetrahedra_depths
            
    def calcCircumcenters(self, dela):
        # per-simplex circumcenters recovered analytically from the Delaunay 
        # paraboloid lifting, avoiding a second Qhull pass. Identical to scipy 
        # Voronoi vertices in general position.
        eq = dela.equations
        scale = dela.paraboloid_scale
        centers = -eq[:, :-2] / (2 * scale * eq[:, -2][:, None])
        return centers

    def buildSparseGraph(self, simplices, neighbors, vertices, points, vdw_radii):
        # one weighted CSR adjacency matrix for the whole cleared state.
        # Edge (tetra -> neigh) weight is l / (d**2 + b) where
        # l is the vertex-to-vertex distance and d is the neighbour's clearance
        # (min over its 4 atoms of |vertex - atom| - vdw_radius)
        from scipy.sparse import csr_matrix

        tetra_points = points[simplices]
        distances = np.linalg.norm(tetra_points - vertices[:, None, :], axis=2)
        bottleneck = np.min(distances - vdw_radii[simplices], axis=1)

        rows = []
        cols = []
        data = []

        b = 1e-3
        for tetra, neighs in enumerate(neighbors):
            for neigh in neighs:
                if neigh == -1:
                    continue
                l = np.linalg.norm(vertices[tetra] - vertices[neigh])
                d = bottleneck[neigh]
                weight = l / (d * d + b)
                rows.append(tetra)
                cols.append(neigh)
                data.append(weight)
        graph = csr_matrix((data, (rows, cols)), shape=(len(simplices), 
                                                        len(simplices)))
        return graph

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

        # A tunnel physically ends at the surface, but the Dijkstra cost has 
        # no such term (it rewards width, and mouths are wide), so a cheapest 
        # path to a far exit can run through/past a nearer mouth. When 
        # truncate_at_surface is set we cut each reconstructed path at the 
        # first qualified mouth it reaches. A mouth is a surface (exit) 
        # tetrahedron whose inscribed clearance is >= bottleneck - one a probe
        #  of that radius can leave through. We test the Voronoi vertices 
        # geometrically, not tetra identity: near the surface many distinct 
        # exit tetra share almost the same circumcenter, so a path can be 
        # inside a mouth while its node is a neighbour,which a tetra-identity 
        # test would miss. Two truncated paths are then treated as the same 
        # channel only if they leave through overlapping mouths AND share most 
        # of their route (see _add_deduped_channels); distinct exits are kept.
        mouth_xyz = np.empty((0, 3))
        mouth_r = np.empty(0)
        if truncate_at_surface:
            exit_tetra = np.asarray(getattr(cavity, 'exit_tetrahedra', 
                                            np.empty(0, dtype=np.intp)))
            if len(exit_tetra):
                verts = vertices[exit_tetra]
                atom_pos = points[simplices[exit_tetra]]
                atom_rad = vdw_radii[simplices[exit_tetra]]
                clearance = (np.linalg.norm(atom_pos - verts[:, None, :], 
                                            axis=2) - atom_rad).min(axis=1)
                q = clearance >= self.bottleneck
                mouth_xyz = verts[q]
                mouth_r = clearance[q]

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

            for exit_global in cavity.end_tetrahedra:
                if exit_global == start_global:
                    continue
                if exit_global not in global_to_local:
                    continue
                exit_local = global_to_local[exit_global]
                if np.isinf(distances[exit_local]):
                    continue

                path_local = paths.get(exit_local)
                if path_local is None:
                    continue

                term_xyz = None
                term_r = 0.0
                if len(mouth_xyz):
                    # walk seed->exit; stop at the first tetra whose Voronoi 
                    # vertex lies inside some qualified mouth's sphere (skip 
                    # the seed). Record that entry point and the radius of the 
                    # mouth entered.. tunnel physically leaves  protein there.
                    for j in range(1, len(path_local)):
                        cc = vertices[cavity_tetra[path_local[j]]]
                        d = np.linalg.norm(mouth_xyz - cc, axis=1)
                        hit = np.nonzero(d < mouth_r)[0]
                        if len(hit):
                            path_local = path_local[:j + 1]
                            term_xyz = cc.copy()
                            term_r = float(mouth_r[hit[np.argmin(d[hit])]])
                            break

                path_global = cavity_tetra[path_local]
                channel = Channel(path_global, *self.processChannel(
                    path_global, vertices, points, vdw_radii, simplices),
                    cost=float(distances[path_local[-1]]))
                candidates.append((channel, term_xyz, term_r, list(path_local)))

        if truncate_at_surface:
            self._addDedupedChannels(cavity, candidates, similarity)
        else:
            for channel, _t, _r, _path in candidates:
                cavity.addChannel(channel)

    def _addDedupedChannels(self, cavity, candidates, similarity):
        # Keep one channel per (surface exit, distinct route). Two truncated 
        # channels are the same tunnel only if they leave through overlapping 
        # mouths (their exit spheres intersect, |Ti - Tj| < ri + rj) AND share 
        # most of their route (diverge late). Exits farther apart than their 
        # mouth radii are distinct openings and kept, even when the paths share
        #  a long trunk and split only near the surface; different corridors to
        #  one exit diverge early (low shared prefix) and are also kept. 
        # omparing the two actual exit points avoids the single-linkage chaining
        # of a mouth-cluster label, which can span many A and merge distinct exits.
        # Cost-sorted greedy, so the kept representative is always the cheapest
        # and the outcome is order-independent.
        kept = []  # (channel, term_xyz, term_r, path)
        for channel, term_xyz, term_r, path in sorted(candidates, key=lambda c: c[0].cost):
            duplicate = False
            if term_xyz is not None:
                for _kc, kxyz, kr, kpath in kept:
                    if kxyz is not None and \
                            np.linalg.norm(term_xyz - kxyz) < term_r + kr and \
                            self._sharedPrefixFraction(path, kpath) >= similarity:
                        duplicate = True
                        break
            if not duplicate:
                kept.append((channel, term_xyz, term_r, path))
        for channel, _t, _r, _path in kept:
            cavity.addChannel(channel)

    @staticmethod
    def _sharedPrefixFraction(a, b):
        # Fraction of the shorter path shared as a common prefix from the seed.
        # Robust to trunk-sharing: distinct tunnels share only the early trunk 
        # (small), a redundant wiggle shares almost everything (~1.0).
        n = 0
        for x, y in zip(a, b):
            if x == y:
                n += 1
            else:
                break
        m = min(len(a), len(b))
        return n / m if m else 0.0

    def calculateMaxRadius(self, vertice, points, vdw_radii, simp):
        atom_positions = points[simp]
        radii = vdw_radii[simp]
        distances = np.linalg.norm(atom_positions - vertice, axis=1) - radii
        return np.min(distances)

    def calculateRadiusSpline(self, tetrahedra, voronoi_vertices, points, 
                              vdw_radii, simp):
        vertices = voronoi_vertices[tetrahedra]
        radii = np.array([self.calculateMaxRadius(v, points, vdw_radii, s) for v, s in zip(vertices, simp[tetrahedra])])
        return radii, np.min(radii)

    def processChannel(self, tetrahedra, voronoi_vertices, points, vdw_radii, 
                       simp):
        from scipy.interpolate import CubicSpline
        
        centers = voronoi_vertices[tetrahedra]
        radii, bottleneck = self.calculateRadiusSpline(tetrahedra, 
                                                       voronoi_vertices, 
                                                       points, vdw_radii, simp)
        
        t = np.arange(len(centers))
        centerline_spline = CubicSpline(t, centers, bc_type='natural')
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
                for i, (x, y, z, radius) in enumerate(zip(centers[:, 0], centers[:, 1], centers[:, 2], radii), start=atom_index):
                    pdb_lines.append("ATOM  %5d  H   FIL T   1    %8.3f%8.3f%8.3f%6.2f%6.2f\n" % (i, x, y, z, 1.00, radius))

                for i in range(1, samples):
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
                        pdb_lines.append("ATOM  %5d  H   FIL T   1    %8.3f%8.3f%8.3f%6.2f%6.2f\n" % (i, x, y, z, 1.00, radius))

                    for i in range(1, samples):
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
            
    def setStartingTetrahedraFromPoint(self, cavities, vertices, start_point, 
                                       restrict=False):
        '''Set starting tetrahedra using a user-defined 3D point.
        The starting tetrahedron of a cavity is the one whose Voronoi vertex is closest
        to `start_point` (Euclidean distance).
        
        :arg cavities: list of cavity objects
        :arg vertices: Voronoi vertices (array of shape (n, 3))
        :arg start_point: point [x, y, z] in Å (list/tuple/ndarray of length 3)
        :arg restrict: if True, only the single cavity whose closest tetrahedron is
            globally nearest to `start_point` is seeded and returned, so channels are
            computed exclusively for the region around `start_point`. If False (default),
            every cavity is seeded with its own closest tetrahedron and all cavities are
            returned unchanged.
        :type restrict: bool
        :returns: list of cavities to search: all cavities when `restrict` is False, the
            single selected cavity when `restrict` is True, or an empty list if no cavity
            has any tetrahedra'''
        
        sp = np.asarray(start_point, dtype=float).reshape(3,)

        best_cavity = None
        best_tetra = None
        best_d2 = np.inf

        for cavity in cavities:
            tet = cavity.tetrahedra
            if tet is None or len(tet) == 0:
                continue

            # Voronoi vertex per tetrahedron: vertices[tetra_id] -> (x,y,z)
            v = vertices[tet]
            d2 = np.sum((v - sp) ** 2, axis=1)
            idx = int(np.argmin(d2))

            if not restrict:
                cavity.setStartingTetrahedron(np.array([tet[idx]]))

            if d2[idx] < best_d2:
                best_d2 = d2[idx]
                best_tetra = tet[idx]
                best_cavity = cavity

        if not restrict:
            return cavities

        if best_cavity is None:
            LOGGER.warn("start_point was provided but no cavity contains any "
                "tetrahedron; no channels will be computed.")
            return []

        best_cavity.setStartingTetrahedron(np.array([best_tetra]))
        LOGGER.info("start_point mapped to tetrahedron {0} (Voronoi vertex {1:.3f} A "
            "away); restricting channel search to the cavity that contains it."
            .format(int(best_tetra), float(np.sqrt(best_d2))))

        return [best_cavity]


    def trimCavitiesByDepth(self, cavities, max_depth):
        """Filtering cavities by max_depth."""
    
        for cavity in cavities:
            cavity.tetrahedra = np.array([
                tetra for tetra in cavity.tetrahedra
                if cavity.tetrahedra_depths.get(tetra, np.inf) <= max_depth])

# -*- coding: utf-8 -*-

"""This module called CaviFinder and defines functions for calculating channels, tunnels and pores
within protein structure.
"""

__author__ = 'Karolina Mikulska-Ruminska', 'Eryk Trzcinski'
__credits__ = ['Karolina Mikulska-Ruminska', 'Eryk Trzcinski']
__email__ = ['karolamik@fizyka.umk.pl']

import numpy as np
from numpy import *
from prody import LOGGER, SETTINGS, PY3K
from prody.atomic import AtomGroup, Atom, Atomic, Selection, Select
from prody.atomic import flags, sliceAtomicData
from prody.utilities import importLA, checkCoords, showFigure, getCoords, isListLike
from prody.measure import calcDistance, calcAngle, calcCenter
from prody.measure.contacts import findNeighbors
from prody.proteins import writePDB, parsePDB, parsePQR
from collections import Counter

from prody.trajectory import TrajBase, Trajectory, Frame
from prody.ensemble import Ensemble

import multiprocessing
from .fixer import *
from .compare import *
from prody.measure import calcTransformation, calcDistance, calcRMSD, superpose


__all__ = ['getVmdModel', 'calcChannels', 'calcChannelsMultipleFrames', 
           'getChannelParameters', 'getChannelAtoms', 'showChannels', 'showCavities',
           'showSurfaceCavities', 'selectChannelBySelection', 'getChannelResidueNames',
           'calcChannelSurfaceOverlaps', 'calcSurfaceCavities']



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
            LOGGER.warn("Package " + str(package_name) + " is not installed. Please install it to use this function.")
            return False
    else:
        try:
            __import__(package_name)
        except ImportError:
            LOGGER.warn("Package " + str(package_name) + " is not installed. Please install it to use this function.")
            return False
    
    return True


def getVmdModel(vmd_path, atoms, representation='NewCartoon'):
    """Generates a 3D model of molecular structures using VMD and returns it as an Open3D TriangleMesh.

    This function creates a temporary PDB file from the provided atomic data and uses VMD (Visual Molecular Dynamics)
    to render this data into an STL file, which is then loaded into Open3D as a TriangleMesh. The function handles
    the creation and cleanup of temporary files and manages the subprocess call to VMD.

    :param vmd_path: Path to the VMD executable. This is required to run VMD and execute the TCL script.
    :type vmd_path: str

    :param atoms: Atomic data to be written to a PDB file. This should be an object or data structure
        that is compatible with the `writePDB` function.
    :type atoms: object

    :raises ImportError: If required libraries ('subprocess', 'pathlib', 'tempfile', 'open3d') are not installed,
        an ImportError is raised, specifying which libraries are missing.

    :raises ValueError: If the STL file is not created or is empty, or if the STL file cannot be read as a TriangleMesh,
        a ValueError is raised.

    :returns: An Open3D TriangleMesh object representing the 3D model generated from the PDB data.
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
                errorMsg = 'To run getVmdModel, please install {0}'.format(missing[0])
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
            "representation must be one of: 'NewCartoon', 'VDW', 'Surf', 'QuickSurf', or 'CPK'")
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
            output_path = os.path.join(os.path.dirname(temp_script.name), "output.stl")

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

    command = [vmd_path, '-e', str(temp_script_path), '-args', str(temp_pdb_path), str(output_path)]

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
    """Visualizes the channels, and optionally, the molecular model and surface, using Open3D.
    
    This function renders a 3D visualization of molecular channels based on their spline representations.
    It can also display a molecular model (e.g., the protein structure) and a surface (e.g., cavity surface)
    in the same visualization. The function utilizes the Open3D library to create and render the 3D meshes.
    
    :arg channels: A list of channel objects or a single channel object. Each channel should have a 
        `get_splines()` method that returns two CubicSpline objects: one for the centerline and one for the radii.
    :type channels: list or single channel object
    
    :arg model: An optional Open3D TriangleMesh object representing the molecular model, such as a protein.
        If provided, this model will be rendered in the visualization.
        Model can be generated using getVmdModel() function.
    :type model: open3d.geometry.TriangleMesh, optional
    
    :arg surface: An optional list containing the surface data. The list should have two elements:
        - `points`: The coordinates of the vertices on the surface.
        - `simp`: The simplices that define the surface (e.g., triangles or tetrahedra).
        If provided, the surface will be rendered as a wireframe overlay in the visualization.
    :type surface: list (with two numpy arrays), optional
    
    :raises ImportError: If the Open3D library is not installed, an ImportError is raised,
        prompting the user to install Open3D.
    
    :returns: None. This function only renders the visualization.
    
    Example usage:
    showChannels(channels, model=protein_mesh, surface=surface_data) """
    
    if not checkAndImport('open3d'):
        errorMsg = 'To run showChannels, please install open3d.'
        raise ImportError(errorMsg)
            
    import open3d as o3d
    
    def create_mesh_from_spline(centerline_spline, radius_spline, n=5):
        N = n * len(centerline_spline.x)
        t = np.linspace(centerline_spline.x[0], centerline_spline.x[-1], N)
        centers = centerline_spline(t)
        radii = radius_spline(t)

        spheres = [o3d.geometry.TriangleMesh.create_sphere(radius=r, resolution=20).translate(c) for r, c in zip(radii, centers)]
        mesh = spheres[0]
        for sphere in spheres[1:]:
            mesh += sphere

        return mesh
    
    if not isinstance(channels, list):
        channels = [channels]
    
    channel_meshes = [create_mesh_from_spline(*channel.get_splines()) for channel in channels]
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
        unique_triangles, counts = np.unique(triangles_tuple, return_counts=True, axis=0)
            
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

    This function displays a 3D visualization of cavities detected in a molecular structure.
    It uses the Open3D library to render the cavities as a triangle mesh. Optionally, it can also
    display the molecular surface as a wireframe overlay.

    :arg surface: A list containing three elements:
        - `points`: The coordinates of the vertices (atoms) in the molecular structure.
        - `surf_simp`: The simplices that define the molecular surface.
        - `simp_cavities`: The simplices corresponding to the detected cavities.
    :type surface: list (with three numpy arrays)

    :arg show_surface: A boolean flag indicating whether to display the molecular surface
        as a wireframe overlay in the visualization. If True, the surface will be displayed
        in addition to the cavities. Default is False.
    :type show_surface: bool

    :raises ImportError: If the Open3D library is not installed, an ImportError is raised,
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
        
    surface_triangles = np.unique(np.array(triangles), axis=0, return_counts=True)[0]
        
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
        unique_triangles, counts = np.unique(triangles_tuple, return_counts=True, axis=0)
            
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
    
    :param surface: Surface data returned by :func:`calcSurfaceCavities`.
        Required for `mode='tetra'`, for `mode='smooth'` when
        `cavity_atoms` is not provided, and for displaying the molecular
        surface when `show_surface=True`. The expected list contains:

        - `surface[0]`: atomic coordinates used for the calculation,
        - `surface[1]`: simplices defining the molecular surface,
        - `surface[2]`: merged cavity simplices,
        - `surface[3]`: simplices after surface-layer removal,
        - `surface[4]`: Voronoi vertices used to represent surface cavities.

    :type surface: list or None

    :param cavities: List of :class:`Cavity` objects returned by
        :func:`calcSurfaceCavities`. Required when `mode='smooth'` and
        `cavity_atoms` is not provided, because the function uses
        `cavity.tetrahedra` to select the corresponding Voronoi vertices.
    :type cavities: list or None

    :param model: Optional Open3D `TriangleMesh` representing the protein or
        another molecular model. The model can be generated with
        :func:`getVmdModel`.
    :type model: open3d.geometry.TriangleMesh, or None

    :param show_surface: If `True`, display the molecular surface wireframe
        derived from `surface[1]` in addition to the cavity representation.
        This requires `surface` to be provided. Default is `False`.
    :type show_surface: bool

    :param mode: Visualization mode used when `cavity_atoms` is not provided.
        Accepted values are `'tetra'` and `'smooth'`. Default is `'tetra'`.
    :type mode: str

    :param alpha: Alpha value used for alpha-shape surface reconstruction in
        `mode='smooth'` and when visualizing `cavity_atoms`. Smaller values
        produce tighter surfaces, while larger values may connect more distant
        points and generate broader surfaces. Default is 4.0.
    :type alpha: float

    :param smoothing: Number of Taubin smoothing iterations applied to the
        reconstructed cavity mesh. If `0` or `None`, no smoothing is
        applied. Default is 0.
    :type smoothing: int or None

    :param cavity_atoms: Optional pseudoatom representation of surface
        cavities. This can be either a path to a PDB/PQR file or a parsed ProDy
        `AtomGroup`, or an Open3D `TriangleMesh` generated, for example, with 
        :func:`getVmdModel`.
    :type cavity_atoms: str, :class:`.AtomGroup`, open3d.geometry.TriangleMesh, or None
    
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
                raise TypeError("cavity_atoms must be a PDB/PQR filename, a ProDy AtomGroup, "
                                "or an Open3D TriangleMesh.")

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

            surface_triangles = np.unique(np.array(triangles), axis=0, return_counts=True)[0]
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
        unique_triangles, counts = np.unique(triangles_tuple, return_counts=True, axis=0)
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


def calcChannels(atoms, output_path=None, separate=False, start_point=None, r1=3, r2=1.25, min_depth=10, 
    min_volume=None, max_volume=None, max_depth=None, bottleneck=1, sparsity=15, 
    min_tetrahedra=None, max_tetrahedra=None, cavities_only=False):
    """Computes and identifies channels within a molecular structure using Voronoi and Delaunay tessellations.

    This function analyzes the provided atomic structure to detect channels, which are voids or pathways
    within the molecular structure. It employs Voronoi and Delaunay tessellations to identify these regions,
    then filters and refines the detected channels based on various parameters such as the minimum depth
    and bottleneck size. The results can be saved to a PQR file (PDB is optional) if an output path is provided. 
    The `separate` parameter controls whether each detected channel is saved to a separate file or if all 
    channels are saved in a single file.

    The implementation is inspired by the methods described in the publication:
    "MOLE 2.0: advanced approach for analysis of biomacromolecular channels" by D. Sehnal, et al., published in 
    J Chemoinform, 5 (39) 2013.

    :param atoms: An object representing the molecular structure, typically containing atomic coordinates
        and element types.
    :type atoms: `Atoms` object

    :param output_path: Optional path to save the resulting channels and associated data in PQR (or PDB) format.
        If None, results are not saved. Default is None.
    :type output_path: str or None

    :param separate: If True, each detected channel is saved to a separate PDB file. If False, all channels
        are saved in a single PDB file. Default is False.
    :type separate: bool

    :param start_point: Optional starting point for channel search. If provided, the algorithm will use 
        the tetrahedron whose Voronoi vertex is closest to this point as the starting tetrahedron (overriding 
        the default automatic seed selection based on the deepest tetrahedron). Coordinates must be given in Å.
    :type start_point: list, tuple, or ndarray (length 3), or None 

    :param r1: The first radius threshold used during the deletion of simplices, which is used to define 
        the outer surface of the channels. Default is 3.
    :type r1: float

    :param r2: The second radius threshold used to define the inner surface of the channels. Default is 1.25.
    :type r2: float

    :param min_depth: The minimum depth a cavity must have to be considered as a channel. Default is 10.
    :type min_depth: int

    :param max_depth: Maximum cavity depth. Cavities deeper than this value are trimmed to the specified depth. 
        Default is None.
    :type max_depth: int

    :param bottleneck: The minimum allowed bottleneck size (narrowest point) for the channels. Default is 1.
    :type bottleneck: float

    :param min_volume: Minimum volume required for a channel/cavity to be retained. Default is None.
    :type min_volume: float

    :param max_volume: Maximum volume allowed for a channel/cavity to be retained. Default is None.
    :type max_volume: float

    :param sparsity: The sparsity parameter controls the sampling density when analyzing the molecular surface.
        A higher value results in fewer sampling points. Default is 15.
    :type sparsity: int

    :returns: A tuple containing two elements:
        - `channels`: A list of detected channels, where each channel is an object containing information
          about its path and geometry.
        - `surface`: A list containing additional information for further visualization, including
          the atomic coordinates, simplices defining the surface, and merged cavities.
    :rtype: tuple (list, list)

    This function performs the following steps:
    1. **Selection and Filtering:** Selects non-hetero atoms from the protein, calculates van der Waals radii, 
        and performs 3D Delaunay triangulation and Voronoi tessellation on the coordinates.
    2. **State Management:** Creates and updates different stages of channel detection of the protein structure 
        to filter out simplices based on the given radii.
    3. **Surface Layer Calculation:** Determines the surface and second-layer simplices from the filtered results.
    4. **Cavity and Channel Detection:** Finds and filters cavities based on their depth and calculates channels 
       using Dijkstra's algorithm.
    5. **Visualization and Saving:** Generates meshes for the detected channels, filters them by bottleneck size, 
       and either saves the results to a PDB file or visualizes them based on the specified parameters.
       
    Example usage:
    channels, surface = calcChannels(atoms, output_path="channels", separate=True)
    
    channels, surface = calcChannels(atoms, output_path="all_channels.pdb", start_point=[-22.312, -20.065, -11.144])
    
    To save the results as PDB file:
    channels, surface = calcChannels(atoms, output_path="channels.pdb", separate=False, r1=3, r2=1.25, min_depth=10, 
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

    from scipy.spatial import Voronoi, Delaunay
    
    if PY3K:
        from pathlib import Path
    else:
        from pathlib2 import Path
    
    if start_point is not None:
        if not isListLike(start_point):
            raise TypeError("start_point must be a list/tuple/ndarray with three numeric values")
        if len(start_point) != 3:
            raise ValueError(
                "start_point must be a list of three numbers, e.g. "
                "start_point=[-12.312, 5.065, -1.144]")
        
        start_point = np.array(start_point, dtype=float)        
        
        LOGGER.info("Using user-provided start_point for channel seed: [{:.3f}, {:.3f}, {:.3f}] Å"
            .format(start_point[0], start_point[1], start_point[2]))


    calculator = ChannelCalculator(atoms, r1, r2, min_depth, bottleneck, sparsity)
    
    atoms = atoms.select('not hetero and noh') # Excluding hydrogens
    coords = atoms.getCoords()
    
    vdw_radii = calculator.get_vdw_radii(atoms.getElements())
        
    dela = Delaunay(coords)
    voro = Voronoi(coords)
        
    s_prt = State(dela.simplices, dela.neighbors, voro.vertices)
    
    if PY3K:
        s_tmp = State(*s_prt.get_state())
        s_prv = State(None, None, None)
    else:
        s_tmp = apply(State, s_prt.get_state())
        s_prv = State(None, None, None) 
        
    while True:
        s_prv.set_state(*s_tmp.get_state())
        
        if PY3K:
            #s_tmp.set_state(*calculator.delete_simplices3d(coords, *(s_tmp.get_state() + [vdw_radii, r1, True])))
            s_tmp.set_state(*calculator.delete_simplices3d(coords, *(s_tmp.get_state() + tuple([vdw_radii, r1, True]))))
        else:
            tmp_state = calculator.delete_simplices3d(coords, *(s_tmp.get_state() + [vdw_radii, r1, True]))
            s_tmp.set_state(*tmp_state)

        if s_tmp == s_prv:
            break
        
    s_srf = State(*s_tmp.get_state())
    #s_inr = State(*calculator.delete_simplices3d(coords, *(s_srf.get_state() + [vdw_radii, r2, False])))
    s_inr = State(*calculator.delete_simplices3d(coords, *(s_srf.get_state() + tuple([vdw_radii, r2, False]))))

    l_first_layer_simp, l_second_layer_simp = calculator.surface_layer(s_srf.simp, s_inr.simp, s_srf.neigh)
    s_clr = State(*calculator.delete_section(l_first_layer_simp, *s_inr.get_state()))
        
    c_cavities = calculator.find_groups(s_clr.neigh)
    c_surface_cavities = calculator.get_surface_cavities(c_cavities, s_clr.simp, l_second_layer_simp, s_clr, coords, vdw_radii, sparsity)
        
    calculator.find_deepest_tetrahedra(c_surface_cavities, s_clr.neigh)
    if start_point is not None:
        calculator.set_starting_tetrahedra_from_point(c_surface_cavities, s_clr.verti, start_point)
    
    c_filtered_cavities = calculator.filter_cavities(c_surface_cavities, min_depth)
    
    if cavities_only and max_depth is not None:
        calculator.trim_cavities_by_depth(c_filtered_cavities, max_depth)
    
    if cavities_only and (min_tetrahedra is not None or max_tetrahedra is not None):
        c_filtered_cavities = calculator.filter_cavities_by_tetrahedra(c_filtered_cavities, min_tetrahedra, max_tetrahedra)
    
    if cavities_only:
        calculator.calculate_cavity_volumes(c_filtered_cavities, s_clr.simp, coords)

    if cavities_only and (min_volume is not None or max_volume is not None):
        c_filtered_cavities = calculator.filter_cavities_by_volume(c_filtered_cavities, min_volume, max_volume)
    
    merged_cavities = calculator.merge_cavities(c_filtered_cavities, s_clr.simp)
    
    # Eraly-return for the calcSurfaceCavities function:
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

            calculator.save_cavities_to_pdb(c_filtered_cavities, s_clr.verti, output_path, separate)
        
        return c_filtered_cavities, [coords, s_srf.simp, merged_cavities, s_clr.simp, s_clr.verti]
        
    for cavity in c_filtered_cavities:
        #calculator.dijkstra(cavity, *(s_clr.get_state() + [coords, vdw_radii]))
        calculator.dijkstra(cavity, *(s_clr.get_state() + tuple([coords, vdw_radii])))
        
    calculator.filter_channels_by_bottleneck(c_filtered_cavities, bottleneck)
    
    if min_volume is not None or max_volume is not None:
        calculator.filter_channels_by_volume(c_filtered_cavities, min_volume, max_volume)
    
    channels = [channel for cavity in c_filtered_cavities for channel in cavity.channels]
        
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
        calculator.save_channels_to_pdb(c_filtered_cavities, output_path, separate)
    else:
        LOGGER.info("No output path given.")
                
    return channels, [coords, s_srf.simp, merged_cavities, s_clr.simp]

            
def calcChannelsMultipleFrames(atoms, trajectory=None, output_path=None, separate=False, start_point=None, **kwargs):
    """Compute channels for each frame in a given trajectory or multi-model PDB file.

    This function calculates the channels for each frame in a trajectory or for each model
    in a multi-model PDB file. The `kwargs` can include parameters necessary for channel calculation.
    If the `separate` parameter is set to True, each detected channel will be saved in a separate PDB file.

    :param atoms: Atomic data or object containing atomic coordinates and methods for accessing them.
    :type atoms: object

    :param trajectory: Trajectory object containing multiple frames or a multi-model PDB file.
    :type trajectory: Atomic or Ensemble object

    :param output_path: Optional path to save the resulting channels and associated data in PDB format.
        If a directory is specified, each frame/model will have its results saved in separate files. 
        If None, results are not saved. Default is None.
    :type output_path: str or None

    :param separate: If True, each detected channel is saved to a separate PDB file for each frame/model.
        If False, all channels for each frame/model are saved in a single file. Default is False.
    :type separate: bool

    :param start_point: Optional starting point for channel search. If provided, the algorithm will use 
        the tetrahedron whose Voronoi vertex is closest to this point as the starting tetrahedron (overriding 
        the default automatic seed selection based on the deepest tetrahedron). Coordinates must be given in Å.
    :type start_point: list, tuple, or ndarray (length 3), or None 

    :param kwargs: Additional parameters required for channel calculation. This can include parameters such as
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


def parseParameters(channels, **kwargs):
    """Extracts and returns the lengths, bottlenecks, and volumes of each channel in a given list of channels. """
    
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
    """Extracts and returns the lengths, bottlenecks, and volumes of each channel in a given list of channels.

    This functaaion iterates through a list of channel objects, extracting the length, bottleneck, 
    and volume of each channel. These values are collected into separate lists, which are returned 
    as a tuple for further use.

    :arg channels: A list of channel objects, where each channel has attributes `length`, `bottleneck`, 
        and `volume`. These attributes represent the length of the channel, the minimum radius 
        (bottleneck) along its path, and the total volume of the channel, respectively.
    :type channels: list

    :arg param_file_name: The files with parameters will be saved in a text file with the provided name.
        Use one word which will be added to '_Parameters_All_channels.txt' sufix.
        If further analysis will be performed with selectChannelBySelection() function, the preferable 
        param_file_name is PDB+chain for example: '1bbhA'.
    :type param_file_name: str 

    :returns: Three lists containing the lengths, bottlenecks, and volumes of the channels.
    :rtype: tuple (list, list, list)

    Example usage:
    lengths, bottlenecks, volumes = getChannelParameters(channels) """
    
    multi_model_param = []
    param_file_name = kwargs.get('param_file_name', None)

    try:
        results_L_B_V = parseParameters(channels, **kwargs)
        lengths, bottlenecks, volumes = results_L_B_V
        LOGGER.info("Channel {0}: \t{1} \t{2} \t{3}".format('ID', 'Volume [Å³]', 'Length [Å]', 'Bottleneck [Å]'))
        for i in range(len(lengths)):
            LOGGER.info("channel {0}: \t{1} \t\t{2} \t\t{3}".format(i, np.round(volumes[i],2), np.round(lengths[i], 2), np.round(bottlenecks[i], 2)))
        return results_L_B_V

    except:
        for nr_i,i in enumerate(channels):
            safe_param_file_name = param_file_name if param_file_name is not None else ""
            results = parseParameters(channels[nr_i], param_file_name=safe_param_file_name + str(nr_i))
            multi_model_param.append(results) 
            
        LOGGER.info("Channel {0}: \t{1} \t{2} \t{3}".format('ID', 'Volume [Å³]', 'Length [Å]', 'Bottleneck [Å]'))
        for frame_nr, frame in enumerate(multi_model_param):
            lengths, bottlenecks, volumes = frame
            LOGGER.info("Frame {0}".format(frame_nr))
            for i in range(len(lengths)):
                LOGGER.info("channel {0}: \t{1} \t\t{2} \t\t{3}".format(i, np.round(volumes[i],2), np.round(lengths[i], 2), np.round(bottlenecks[i], 2)))
        return multi_model_param


def getChannelAtoms(channels, protein=None, num_samples=5):
    """Generates an AtomGroup object representing the atoms along the paths of the given channels
    and optionally combines them with an existing protein structure.

    This function takes a list of channel objects and generates atomic representations of the
    channels based on their centerline splines and radius splines. The function samples points
    along each channel's centerline and assigns atom positions at these points with corresponding
    radii, creating a list of PDB-formatted lines. These lines are then converted into an AtomGroup
    object using the ProDy library. If a protein structure is provided, it is combined with the
    generated channel atoms by merging their respective PDB streams.

    :param channels: A list of channel objects. Each channel has a method `get_splines()` that
        returns the centerline spline and radius spline of the channel.
    :type channels: list

    :param protein: An optional AtomGroup object representing a protein structure. If provided, 
        it will be combined with the generated channel atoms.
    :type protein: prody.atomic.AtomGroup or None

    :param num_samples: The number of atom samples to generate along each segment of the channel.
        More samples result in a finer representation of the channel. Default is 5.
    :type num_samples: int

    :returns: An AtomGroup object representing the combined atoms of the channels and the protein,
        if a protein is provided.
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
        centerline_spline, radius_spline = channel.get_splines()
        samples = len(channel.tetrahedra) * num_samples
        t = np.linspace(centerline_spline.x[0], centerline_spline.x[-1], samples)
        centers = centerline_spline(t)
        radii = radius_spline(t)

        for i, (x, y, z, radius) in enumerate(zip(centers[:, 0], centers[:, 1], centers[:, 2], radii), start=atom_index):
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

    :param channels: A list of channel objects. Each channel has a method `get_splines()` that
        returns the centerline spline and radius spline of the channel.
    :type channels: list

    :arg distA: Residues will be provided based on this value.
        default is 4 [Ang]
    :type distA: int, float 
    
    :arg residues_file_name: The file with residues will be saved in a text file with the provided name.
        Use one word which will be added to '_Residues_All_channels.txt' sufix.
        If further analysis will be performed with selectChannelBySelection() function, the preferable 
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
        with open(residues_file_name+'_Residues_All_channels.txt', "a") as f_res:
            for k in selected_residues_ch:
                f_res.write(("{0}_{1}\n".format(residues_file_name, k)))
                
    return selected_residues_ch


def selectChannelBySelection(atoms, residue_sele, **kwargs):
    """Select PQR files with channels that are having FIL residues within certain distance (distA) from 
    selected residue (temporarly one residue).
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
                    default is False, which means that all .pqr files from the current directory will be analyzed.
    :type pqr_files: bool or list
    
    :arg folder_name: The name of the folder to which PDBs will be extracted
    :type folder_name: str

    :arg distA: non-zero value, maximal distance from selected region to channel (FIL atoms)
                default is 5
    :type distA: int, float 
        
    :arg residues_file: File with residues forming the channel created by getChannelResidues()
                        default is False 
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
    import numpy as np
    
    pqr_files = kwargs.pop('pqr_files', False)
    distA = kwargs.pop('distA', 5)
    folder_name = kwargs.pop('folder_name', 'selected_files')
    residues_file = kwargs.pop('residues_file', False)
    param_file = kwargs.pop('param_file', False)
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
                PDB_id = file[:-4].split('_channel')[0] 
                channel_name = file[:-4].split('_')[-1]
                f = open(PDB_id+'_Residues_All_channels.txt', 'r').readlines()
                for line in f:
                    if line.startswith(PDB_id+'_'+channel_name+':'):
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


def calcChannelSurfaceOverlaps(**kwargs):
    """Calculate overlapping parts of the predicted channels, tunnels, and pores denote as 'FIL' atoms.
    Results are normalized within [0,1].

    :arg resolution: Surface sampling resolution.
        default is 0.5
    :type resolution: float

    :arg output_file_name: The name of the PDB file with overlapping surfaces.
    :type output_file_name: str

    :arg pqr_files: File with residues forming the channel created by getChannelResidues()
        default is False (then all the files from the current directory will be analyzed)
        when providing a list, only the PDBs from list will be analyzed
        when providing str, it will be treated as a folder path  
    :type pqr_files: bool, list or str
    
    Example usage:
    calcChannelSurfaceOverlaps() - all the files in the current directory will be analyzed
    
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


def calcSurfaceCavities(atoms, output_path=None, r1=4.5, r2=2.0, min_depth=2, max_depth=3, 
    min_tetrahedra=None, max_tetrahedra=None, min_volume=50, max_volume=None,
    sparsity=15, separate=False):
    """Calculate surface cavities (pockets) on protein surface using CaviTracer approach.

    :param atoms: An object representing the molecular structure, typically containing atomic coordinates
        and element types.
    :type atoms: `Atoms` object

    :param output_path: Optional path to save the resulting cavities and associated data in PQR (or PDB) format.
        If None, results are not saved. Default is None.
    :type output_path: str or None

    :param separate: If True, each detected cavity is saved to a separate PQR file. If False, all cavities
        are saved in a single PQR file. Default is False.
    :type separate: bool

    :param r1: The first radius threshold used during the deletion of simplices, which is used to define 
        the outer surface of the cavities. Default is 4.5.
    :type r1: float

    :param r2: The second radius threshold used to define the inner surface of the cavities. Default is 2.
    :type r2: float

    :param min_depth: The minimum depth a cavity must have to be considered as a cavity. Default is 2.
    :type min_depth: int

    :param max_depth: Maximum cavity depth. Cavities deeper than this value are trimmed to the specified depth. 
        Default is 3.
    :type max_depth: int

    :param sparsity: The sparsity parameter controls the sampling density when analyzing the molecular surface.
        A higher value results in fewer sampling points. Default is 15.
    :type sparsity: int

    :param min_tetrahedra: Minimum number of tetrahedra required for a cavity to be retained. 
        Smaller cavities are discarded. Default is None.
    :type min_tetrahedra: int

    :param max_tetrahedra: Maximum number of tetrahedra allowed for a cavity to be retained.
        Larger cavities are discarded. Default is None.
    :type max_tetrahedra: int

    :param min_volume: Minimum volume required for a channel/cavity to be retained. Default is 50.
    :type min_volume: float

    :param max_volume: Maximum volume allowed for a channel/cavity to be retained. Default is None.
    :type max_volume: float

    :returns: A tuple containing two elements:
        - `cavities`: A list of detected cavities, where each channel is an object containing information
          about its path and geometry.
        - `surface`: A list containing additional information for further visualization, including
          the atomic coordinates, simplices defining the surface, and merged cavities.
    :rtype: tuple (list, list)

    This function performs the following steps:
    1. **Selection and Filtering:** Selects non-hetero atoms from the protein, calculates van der Waals radii, 
        and performs 3D Delaunay triangulation and Voronoi tessellation on the coordinates.
    2. **Surface and Interior Filtering:** Iteratively removes simplices based on the user-defined radii 
        (`r1` and `r2`) to distinguish the molecular surface from the internal void space.
    3. **Surface Cavity Identification:** Detects connected void regions and identifies those that remain 
        connected to the protein surface, corresponding to surface-accessible cavities and pockets.
    4. **Depth Calculation and Filtering:** Estimates cavity depth using a graph-based traversal from the cavity 
        openings, identifies the deepest tetrahedra, and filters cavities according to the specified depth criteria.
    5. **Output Generation:** Optionally trims cavities exceeding the specified	maximum depth, saves detected 
        cavities to PDB/PQR files, and returns cavity objects together with the surface representation for further 
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
    def __init__(self, tetrahedra, centerline_spline, radius_spline, length, bottleneck, volume):
        self.tetrahedra = tetrahedra
        self.centerline_spline = centerline_spline
        self.radius_spline = radius_spline
        self.length = length
        self.bottleneck = bottleneck
        self.volume = volume
            
    def get_splines(self):
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
        
    def set_state(self, simplices, neighbors, vertices):
        self.simp = simplices
        self.neigh = neighbors
        self.verti = vertices
        
    def get_state(self):
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
        
    def make_surface(self):
        self.is_connected_to_surface = True
        
    def set_exit_tetrahedra(self, exit_tetrahedra, end_tetrahedra):
        self.exit_tetrahedra = exit_tetrahedra
        self.end_tetrahedra = end_tetrahedra
        
    def set_starting_tetrahedron(self, tetrahedron):
        self.starting_tetrahedron = tetrahedron
        
    def set_depth(self, depth):
        self.depth = depth
        
    def add_channel(self, channel):
        self.channels.append(channel)
    
        
class ChannelCalculator:
    def __init__(self, atoms, r1=3, r2=1.25, min_depth=10, bottleneck=1, sparsity=15):
        self.atoms = atoms
        self.r1 = r1
        self.r2 = r2
        self.min_depth = min_depth
        self.bottleneck = bottleneck
        self.sparsity = sparsity
        
    def sphere_fit(self, vertices, tetrahedron, vertice, vdw_radii, r):
        center = vertice
        d_sum = sum(np.linalg.norm(center - vertices[atom]) for atom in tetrahedron)
        r_sum = sum(r + vdw_radii[atom] for atom in tetrahedron)
        
        return d_sum >= r_sum

    def delete_simplices3d(self, points, simplices, neighbors, vertices, vdw_radii, r, surface):
        simp, neigh, verti, deleted = [], [], [], []
        
        for i, tetrahedron in enumerate(simplices):
            should_delete = (-1 in neighbors[i] and self.sphere_fit(points, tetrahedron, vertices[i], vdw_radii, r)) if surface else not self.sphere_fit(points, tetrahedron, vertices[i], vdw_radii, r)
            
            if should_delete:
                deleted.append(i)
            else:
                simp.append(simplices[i])
                neigh.append(neighbors[i])
                verti.append(vertices[i])
        
        simp = np.array(simp)
        neigh = np.array(neigh)
        verti = np.array(verti)
        deleted = np.array(deleted)
        
        mask = np.isin(neigh, deleted)
        neigh[mask] = -1
        
        for i in reversed(deleted):
            mask = (neigh > i) & (neigh != -1)
            neigh[mask] -= 1
        
        return simp, neigh, verti

    def delete_section(self, simplices_subset, simplices, neighbors, vertices, reverse=False):
        simp, neigh, verti, deleted = [], [], [], []
      
        for i, tetrahedron in enumerate(simplices): 
            match = any((simplices_subset == tetrahedron).all(axis=1))
            if reverse:
                if match:
                    simp.append(tetrahedron)
                    neigh.append(neighbors[i])
                    verti.append(vertices[i])
                else:
                    deleted.append(i)
            else:
                if match:
                    deleted.append(i)
                else:
                    simp.append(tetrahedron)
                    neigh.append(neighbors[i])
                    verti.append(vertices[i])

        simp, neigh, verti = map(np.array, [simp, neigh, verti])
        deleted = np.array(deleted)
        
        mask = np.isin(neigh, deleted)
        neigh[mask] = -1
        
        for i in reversed(deleted):
            neigh = np.where((neigh > i) & (neigh != -1), neigh - 1, neigh)
        
        return simp, neigh, verti

    def get_vdw_radii(self, atoms):
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

    def surface_layer(self, shape_simplices, filtered_simplices, shape_neighbors):
        surface_simplices, surface_neighbors = [], []
        interior_simplices = []
        
        for i in range(len(shape_simplices)): 
            if -1 in shape_neighbors[i]:
                surface_simplices.append(shape_simplices[i])
                surface_neighbors.append(shape_neighbors[i])
            else:
                interior_simplices.append(shape_simplices[i])
        
        surface_simplices = np.array(surface_simplices)
        surface_neighbors = np.array(surface_neighbors)
        interior_simplices = np.array(interior_simplices)
        
        filtered_surface_simplices = surface_simplices[
            np.any(np.all(surface_simplices[:, None] == filtered_simplices, axis=2), axis=1)
        ]
        filtered_surface_neighbors = surface_neighbors[
            np.any(np.all(surface_simplices[:, None] == filtered_simplices, axis=2), axis=1)
        ]
        
        filtered_surface_neighbors = np.unique(filtered_surface_neighbors)
        filtered_surface_neighbors = filtered_surface_neighbors[filtered_surface_neighbors != 0]
        
        filtered_interior_simplices = interior_simplices[
            np.any(np.all(interior_simplices[:, None] == filtered_simplices, axis=2), axis=1)
        ]

        surface_layer_neighbor_simplices = shape_simplices[filtered_surface_neighbors]
        
        second_layer = filtered_interior_simplices[
            np.any(np.all(filtered_interior_simplices[:, None] == surface_layer_neighbor_simplices, axis=2), axis=1)
        ]

        return filtered_surface_simplices, second_layer

            
    def find_groups(self, neigh, is_cavity=True):
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

    def get_surface_cavities(self, cavities, interior_simplices, second_layer, state, points, vdw_radii, sparsity):
        surface_cavities = []
        
        for cavity in cavities:
            tetrahedra = cavity.tetrahedra
            second_layer_mask = np.isin(interior_simplices[tetrahedra], second_layer).all(axis=1)
            
            if np.any(second_layer_mask):
                cavity.make_surface()
                exit_tetrahedra = tetrahedra[second_layer_mask]
                end_tetrahedra = self.get_end_tetrahedra(exit_tetrahedra, state.verti, points, vdw_radii, state.simp, sparsity)
                cavity.set_exit_tetrahedra(exit_tetrahedra, end_tetrahedra)
                surface_cavities.append(cavity)
                
        return surface_cavities


    def merge_cavities(self, cavities, simplices):
        merged_tetrahedra = np.concatenate([cavity.tetrahedra for cavity in cavities])
        return simplices[merged_tetrahedra]

    def find_deepest_tetrahedra(self, cavities, neighbors):
        from collections import deque
        
        for cavity in cavities:
            exit_tetrahedra = cavity.exit_tetrahedra
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
                    if neighbor != -1 and not visited[neighbor] and neighbor in cavity.tetrahedra:
                        visited[neighbor] = True
                        queue.append((neighbor, depth + 1))

            cavity.set_starting_tetrahedron(np.array([deepest_tetrahedron]))
            cavity.set_depth(max_depth)
            cavity.tetrahedra_depths = tetrahedra_depths
            
    def dijkstra(self, cavity, simplices, neighbors, vertices, points, vdw_radii):
        import heapq
        
        def calculate_weight(current_tetra, neighbor_tetra):
            current_vertex = vertices[current_tetra]
            neighbor_vertex = vertices[neighbor_tetra]
            l = np.linalg.norm(current_vertex - neighbor_vertex)
            
            d = np.inf
            for atom, radius in zip(points[simplices[neighbor_tetra]], vdw_radii[simplices[neighbor_tetra]]):
                dist = np.linalg.norm(neighbor_vertex - atom) - radius
                if dist < d:
                    d = dist
            
            b = 1e-3
            return l / (d**2 + b)
        
        def dijkstra_algorithm(start, goal, tetrahedra_set):
            pq = [(0, start)]
            distances = {start: 0}
            previous = {start: None}

            while pq:
                current_distance, current_tetra = heapq.heappop(pq)
                
                if current_tetra == goal:
                    path = []
                    while current_tetra is not None:
                        path.append(current_tetra)
                        current_tetra = previous[current_tetra]
                    return path[::-1]
                
                if current_distance > distances[current_tetra]:
                    continue
                
                for neighbor in neighbors[current_tetra]:
                    if neighbor in tetrahedra_set:
                        weight = calculate_weight(current_tetra, neighbor)
                        distance = current_distance + weight
                        if distance < distances.get(neighbor, float('inf')):
                            distances[neighbor] = distance
                            previous[neighbor] = current_tetra
                            heapq.heappush(pq, (distance, neighbor))
            
            return None

        tetrahedra_set = set(cavity.tetrahedra)
        for exit_tetrahedron in cavity.end_tetrahedra:
            for starting_tetrahedron in cavity.starting_tetrahedron:
                if exit_tetrahedron != starting_tetrahedron:
                    path = dijkstra_algorithm(starting_tetrahedron, exit_tetrahedron, tetrahedra_set)
                    if path:
                        path_tetrahedra = np.array(path)
                        channel = Channel(path_tetrahedra, *self.process_channel(path_tetrahedra, vertices, points, vdw_radii, simplices))
                        cavity.add_channel(channel)

    def calculate_max_radius(self, vertice, points, vdw_radii, simp):
        atom_positions = points[simp]
        radii = vdw_radii[simp]
        distances = np.linalg.norm(atom_positions - vertice, axis=1) - radii
        return np.min(distances)

    def calculate_radius_spline(self, tetrahedra, voronoi_vertices, points, vdw_radii, simp):
        vertices = voronoi_vertices[tetrahedra]
        radii = np.array([self.calculate_max_radius(v, points, vdw_radii, s) for v, s in zip(vertices, simp[tetrahedra])])
        return radii, np.min(radii)

    def process_channel(self, tetrahedra, voronoi_vertices, points, vdw_radii, simp):
        from scipy.interpolate import CubicSpline
        
        centers = voronoi_vertices[tetrahedra]
        radii, bottleneck = self.calculate_radius_spline(tetrahedra, voronoi_vertices, points, vdw_radii, simp)
        
        t = np.arange(len(centers))
        centerline_spline = CubicSpline(t, centers, bc_type='natural')
        radius_spline = CubicSpline(t, radii, bc_type='natural')
        
        length = self.calculate_channel_length(centerline_spline)
        volume = self.calculate_channel_volume(centerline_spline, radius_spline)
        
        return centerline_spline, radius_spline, length, bottleneck, volume

    def find_biggest_tetrahedron(self, tetrahedra, voronoi_vertices, points, vdw_radii, simp):
        radii = np.array([self.calculate_max_radius(voronoi_vertices[tetra], points, vdw_radii, simp[tetra]) for tetra in tetrahedra])
        max_radius_index = np.argmax(radii)
        return tetrahedra[max_radius_index]

    def get_end_tetrahedra(self, tetrahedra, voronoi_vertices, points, vdw_radii, simp, sparsity):
        end_tetrahedra = []
        current_tetrahedron = self.find_biggest_tetrahedron(tetrahedra, voronoi_vertices, points, vdw_radii, simp)
        end_tetrahedra.append(current_tetrahedron)
        end_tetrahedra_set = {current_tetrahedron}
        
        while True:
            found_tetrahedra = []
            for tetra in tetrahedra:
                if tetra in end_tetrahedra_set:
                    continue

                all_far_enough = True  
                for selected_tetra in end_tetrahedra:
                    distance = np.linalg.norm(voronoi_vertices[selected_tetra] - voronoi_vertices[tetra])
                    if distance < sparsity:
                        all_far_enough = False
                        break
                
                if all_far_enough:
                    found_tetrahedra.append(tetra)

            if not found_tetrahedra:
                break
            
            biggest_tetrahedron = self.find_biggest_tetrahedron(found_tetrahedra, voronoi_vertices, points, vdw_radii, simp)
            end_tetrahedra.append(biggest_tetrahedron)
            end_tetrahedra_set.add(biggest_tetrahedron)

        return np.array(end_tetrahedra)

    def filter_cavities(self, cavities, min_depth):
        return [cavity for cavity in cavities if cavity.depth >= min_depth]

    def filter_channels_by_bottleneck(self, cavities, bottleneck):
        for cavity in cavities:
            cavity.channels = [channel for channel in cavity.channels if channel.bottleneck >= bottleneck]
    
    def filter_channels_by_volume(self, cavities, min_volume=None, max_volume=None):
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

    def filter_cavities_by_tetrahedra(self, cavities, min_tetrahedra=None, max_tetrahedra=None):
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

    def calculate_tetrahedron_volume(self, a, b, c, d):
        return abs(np.dot(a - d, np.cross(b - d, c - d))) / 6.0

    def calculate_cavity_volumes(self, cavities, simplices, coords):
        """Calculate approximate cavity volumes from Delaunay tetrahedra."""

        for cavity in cavities:
            volume = 0.0
            for tetra in cavity.tetrahedra:
                atom_ids = simplices[tetra]
                a, b, c, d = coords[atom_ids]
                volume += self.calculate_tetrahedron_volume(a, b, c, d)
            cavity.volume = volume

    def filter_cavities_by_volume(self, cavities, min_volume=None, max_volume=None):
        """Filter cavities by approximate volume."""

        filtered_cavities = []
        for cavity in cavities:
            if min_volume is not None and cavity.volume < min_volume:
                continue
            if max_volume is not None and cavity.volume > max_volume:
                continue
            filtered_cavities.append(cavity)
        return filtered_cavities

    def save_channels_to_pdb(self, cavities, filename, separate=False, num_samples=5):
        filename = str(filename)
        
        # All channels will be provided always when PDB/PQR will be created
        with open(filename, 'w') as pqr_file:
            atom_index = 1
            for cavity in cavities:
                for channel in cavity.channels:
                    centerline_spline, radius_spline = channel.get_splines()
                    samples = len(channel.tetrahedra) * num_samples
                    t = np.linspace(centerline_spline.x[0], centerline_spline.x[-1], samples)
                    centers = centerline_spline(t)
                    radii = radius_spline(t)

                    pdb_lines = []
                    for i, (x, y, z, radius) in enumerate(zip(centers[:, 0], centers[:, 1], centers[:, 2], radii), start=atom_index):
                        pdb_lines.append("ATOM  %5d  H   FIL T   1    %8.3f%8.3f%8.3f%6.2f%6.2f\n" % (i, x, y, z, 1.00, radius))

                    for i in range(1, samples):
                        pdb_lines.append("CONECT%5d%5d\n" % (i, i + 1))
                        
                    pqr_file.writelines(pdb_lines)
                    pqr_file.write("\n")
                    atom_index += samples
        
        # When separate is set to True also separate PDB/PQR files will be created
        if separate:
            channel_index = 0
            for cavity in cavities:
                for channel in cavity.channels:
                    channel_filename = filename.replace('.pqr', '_channel{0}.pqr'.format(channel_index))
                    
                    with open(channel_filename, 'w') as pqr_file:
                        atom_index = 1
                        centerline_spline, radius_spline = channel.get_splines()
                        samples = len(channel.tetrahedra) * num_samples
                        t = np.linspace(centerline_spline.x[0], centerline_spline.x[-1], samples)
                        centers = centerline_spline(t)
                        radii = radius_spline(t)
    
                        pdb_lines = []
                        for i, (x, y, z, radius) in enumerate(zip(centers[:, 0], centers[:, 1], centers[:, 2], radii), start=atom_index):
                            pdb_lines.append("ATOM  %5d  H   FIL T   1    %8.3f%8.3f%8.3f%6.2f%6.2f\n" % (i, x, y, z, 1.00, radius))
    
                        for i in range(1, samples):
                            pdb_lines.append("CONECT%5d%5d\n" % (i, i + 1))
                            
                        pqr_file.writelines(pdb_lines)
                        
                    channel_index += 1


    def save_cavities_to_pdb(self, cavities, vertices, filename, separate=False):
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
            
    def calculate_channel_length(self, centerline_spline):
        t_values = np.linspace(centerline_spline.x[0], centerline_spline.x[-1], len(centerline_spline.x) * 10)
        points = centerline_spline(t_values)
        diffs = np.diff(points, axis=0)
        lengths = np.linalg.norm(diffs, axis=1)
        return np.sum(lengths)
    
    def calculate_channel_volume(self, centerline_spline, radius_spline):
        import warnings
        from scipy.integrate import quad, IntegrationWarning  
        
        warnings.filterwarnings("ignore", category=IntegrationWarning)
        
        t_min = centerline_spline.x[0]
        t_max = centerline_spline.x[-1]
    
        def differential_volume(t):
            r = radius_spline(t)
            area = np.pi * r**2
            dx_dt = centerline_spline(t, 1)
            centerline_derivative = np.linalg.norm(dx_dt)
            return area * centerline_derivative
        
        volume, error = quad(differential_volume, t_min, t_max)
        
        r_start = radius_spline(t_min)
        r_end = radius_spline(t_max)
        
        hemisphere_volume_start = (2/3) * np.pi * r_start**3
        hemisphere_volume_end = (2/3) * np.pi * r_end**3
        
        total_volume = volume + hemisphere_volume_start + hemisphere_volume_end
        
        return total_volume
            
    def set_starting_tetrahedra_from_point(self, cavities, vertices, start_point):
        '''Set starting tetrahedra using a user-defined 3D point.
        The starting tetrahedron is selected as the one whose Voronoi vertex is closest
        to `start_point` (Euclidean distance).
        
        :arg cavities: list of cavity objects
        :arg vertices: Voronoi vertices (array of shape (n, 3))
        :arg start_point: point [x, y, z] in Å (list/tuple/ndarray of length 3)'''
        
        sp = np.asarray(start_point, dtype=float).reshape(3,)

        for cavity in cavities:
            tet = cavity.tetrahedra
            if tet is None or len(tet) == 0:
                continue

            # Voronoi vertex per tetrahedron: vertices[tetra_id] -> (x,y,z)
            v = vertices[tet]
            d2 = np.sum((v - sp) ** 2, axis=1)
            chosen = tet[int(np.argmin(d2))]

            cavity.set_starting_tetrahedron(np.array([chosen]))


    def trim_cavities_by_depth(self, cavities, max_depth):
        """Filtering cavities by max_depth."""
    
        for cavity in cavities:
            cavity.tetrahedra = np.array([
                tetra for tetra in cavity.tetrahedra
                if cavity.tetrahedra_depths.get(tetra, np.inf) <= max_depth])

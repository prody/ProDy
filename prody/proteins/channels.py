"""This module detects channels in molecules."""

import subprocess
import heapq
import numpy as np
import open3d as o3d
from collections import deque
from scipy.interpolate import CubicSpline
from scipy.spatial import Voronoi, Delaunay
from pathlib import Path

__all__ = ['run_vmd_script', 'detect_channels']

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
        
class Channel:
    def __init__(self, tetrahedra, centerline_spline, radius_spline, length, bottleneck):
        self.tetrahedra = tetrahedra
        self.centerline_spline = centerline_spline
        self.radius_spline = radius_spline
        self.length = length
        self.bottleneck = bottleneck
        
    def get_splines(self):
        return self.centerline_spline, self.radius_spline

def visualize_external_grid(points, simp, stl_file=None, other_mesh=None, channel_mesh=None, ret_lines=None):
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
    
    if ret_lines:
        return line_set
    
    geometries = [line_set, o3d.geometry.TriangleMesh.create_coordinate_frame(size=0.1, origin=[0, 0, 0])]

    if stl_file is not None:
        stl_mesh = o3d.io.read_triangle_mesh(stl_file)
        stl_mesh.compute_vertex_normals()
        stl_mesh.paint_uniform_color([0.1, 0.7, 0.3])
        geometries.append(stl_mesh)

    if channel_mesh is not None:
        if not isinstance(channel_mesh, list):
            channel_mesh = [channel_mesh]
        for mesh in channel_mesh:
            mesh.compute_vertex_normals()
            mesh.paint_uniform_color([0.5, 0.0, 0.5])
        geometries.extend(channel_mesh)

    o3d.visualization.draw_geometries(geometries)

def visualize_external_mesh(prot_coords, prot_simp, lines=None, alpha=0.5):
    triangles = []
    for tetra in prot_simp:
        triangles.extend([sorted([tetra[0], tetra[1], tetra[2]]),
                          sorted([tetra[0], tetra[1], tetra[3]]),
                          sorted([tetra[0], tetra[2], tetra[3]]),
                          sorted([tetra[1], tetra[2], tetra[3]])])
    
    surface_triangles = np.unique(np.array(triangles), axis=0, return_counts=True)[0]
    
    mesh = o3d.geometry.TriangleMesh()
    mesh.vertices = o3d.utility.Vector3dVector(prot_coords)
    mesh.triangles = o3d.utility.Vector3iVector(surface_triangles)
    
    mesh.compute_vertex_normals()
    mesh.paint_uniform_color([0.1, 0.7, 0.3])
    
    vis = o3d.visualization.Visualizer()
    vis.create_window()
    vis.add_geometry(mesh)
    
    if lines:
        vis.add_geometry(lines)
        
    vis.get_render_option().mesh_show_back_face = True
    vis.get_render_option().background_color = np.array([1, 1, 1])
    vis.update_renderer()
    vis.run()
    vis.destroy_window()

def visualize_channel(channel_meshes, stl_file=None):
    meshes_to_visualize = [o3d.geometry.TriangleMesh.create_coordinate_frame(size=0.1, origin=[0, 0, 0])]

    if stl_file is not None:
        mesh = o3d.io.read_triangle_mesh(stl_file)
        mesh.compute_vertex_normals()
        mesh.paint_uniform_color([0.1, 0.7, 0.3])
        meshes_to_visualize.append(mesh)
        
    if channel_meshes is not None:
        if not isinstance(channel_meshes, list):
            channel_meshes = [channel_meshes]
        for channel_mesh in channel_meshes:
            channel_mesh.compute_vertex_normals()
            channel_mesh.paint_uniform_color([0.5, 0.0, 0.5])
        meshes_to_visualize.extend(channel_meshes)
    
    if len(meshes_to_visualize) > 1:
        o3d.visualization.draw_geometries(meshes_to_visualize)
    else:
        print("No mesh to visualize.")

def sphere_fit(vertices, tetrahedron, vertice, vdw_radii, r):
    center = vertice
    d_sum = sum(np.linalg.norm(center - vertices[atom]) for atom in tetrahedron)
    r_sum = sum(r + vdw_radii[atom] for atom in tetrahedron)
    
    return d_sum >= r_sum

def delete_simplices3d(points, simplices, neighbors, vertices, vdw_radii, r, surface):
    simp, neigh, verti, deleted = [], [], [], []
    
    for i, tetrahedron in enumerate(simplices):
        should_delete = (-1 in neighbors[i] and sphere_fit(points, tetrahedron, vertices[i], vdw_radii, r)) if surface else not sphere_fit(points, tetrahedron, vertices[i], vdw_radii, r)
        
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

def delete_section(simplices_subset, simplices, neighbors, vertices, reverse=False):
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

def get_vdw_radii(atoms):
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

def surface_layer(shape_simplices, filtered_simplices, shape_neighbors):
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

        
def find_groups(neigh, is_cavity=True):
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

def get_surface_cavities(cavities, interior_simplices, second_layer, state, points, vdw_radii, sparsity):
    surface_cavities = []
    
    for cavity in cavities:
        tetrahedra = cavity.tetrahedra
        second_layer_mask = np.isin(interior_simplices[tetrahedra], second_layer).all(axis=1)
        
        if np.any(second_layer_mask):
            cavity.make_surface()
            exit_tetrahedra = tetrahedra[second_layer_mask]
            end_tetrahedra = get_end_tetrahedra(exit_tetrahedra, state.verti, points, vdw_radii, state.simp, sparsity)
            cavity.set_exit_tetrahedra(exit_tetrahedra, end_tetrahedra)
            surface_cavities.append(cavity)
            
    return surface_cavities


def merge_cavities(cavities, simplices):
    merged_tetrahedra = np.concatenate([cavity.tetrahedra for cavity in cavities])
    return simplices[merged_tetrahedra]

def find_deepest_tetrahedra(cavities, neighbors):
    for cavity in cavities:
        exit_tetrahedra = cavity.exit_tetrahedra
        visited = np.zeros(neighbors.shape[0], dtype=bool)
        visited[exit_tetrahedra] = True
        queue = deque([(tetra, 0) for tetra in exit_tetrahedra])
        max_depth = -1
        deepest_tetrahedron = None

        while queue:
            current, depth = queue.popleft()
            if depth > max_depth:
                max_depth = depth
                deepest_tetrahedron = current

            for neighbor in neighbors[current]:
                if neighbor != -1 and not visited[neighbor] and neighbor in cavity.tetrahedra:
                    visited[neighbor] = True
                    queue.append((neighbor, depth + 1))

        cavity.set_starting_tetrahedron(np.array([deepest_tetrahedron]))
        cavity.set_depth(max_depth)
        
def dijkstra(cavity, simplices, neighbors, vertices, points, vdw_radii):
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
                    channel = Channel(path_tetrahedra, *process_channel(path_tetrahedra, vertices, points, vdw_radii, simplices))
                    cavity.add_channel(channel)

def calculate_max_radius(vertice, points, vdw_radii, simp):
    atom_positions = points[simp]
    radii = vdw_radii[simp]
    distances = np.linalg.norm(atom_positions - vertice, axis=1) - radii
    return np.min(distances)

def calculate_radius_spline(tetrahedra, voronoi_vertices, points, vdw_radii, simp):
    vertices = voronoi_vertices[tetrahedra]
    radii = np.array([calculate_max_radius(v, points, vdw_radii, s) for v, s in zip(vertices, simp[tetrahedra])])
    return radii, np.min(radii)

def process_channel(tetrahedra, voronoi_vertices, points, vdw_radii, simp):
    centers = voronoi_vertices[tetrahedra]
    radii, bottleneck = calculate_radius_spline(tetrahedra, voronoi_vertices, points, vdw_radii, simp)
    
    t = np.arange(len(centers))
    centerline_spline = CubicSpline(t, centers, bc_type='natural')
    radius_spline = CubicSpline(t, radii, bc_type='natural')
    
    length = calculate_channel_length(centerline_spline)
    
    return centerline_spline, radius_spline, length, bottleneck

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

def find_biggest_tetrahedron(tetrahedra, voronoi_vertices, points, vdw_radii, simp):
    radii = np.array([calculate_max_radius(voronoi_vertices[tetra], points, vdw_radii, simp[tetra]) for tetra in tetrahedra])
    max_radius_index = np.argmax(radii)
    return tetrahedra[max_radius_index]

def get_end_tetrahedra(tetrahedra, voronoi_vertices, points, vdw_radii, simp, sparsity):
    end_tetrahedra = []

    current_tetrahedron = find_biggest_tetrahedron(tetrahedra, voronoi_vertices, points, vdw_radii, simp)
    end_tetrahedra.append(current_tetrahedron)
    
    end_tetrahedra_set = {current_tetrahedron}
    
    while True:
        found_tetrahedra = []

        for tetra in tetrahedra:
            if tetra in end_tetrahedra_set:
                continue

            if all(np.linalg.norm(voronoi_vertices[selected_tetra] - voronoi_vertices[tetra]) >= sparsity for selected_tetra in end_tetrahedra):
                found_tetrahedra.append(tetra)

        if not found_tetrahedra:
            break
        
        biggest_tetrahedron = find_biggest_tetrahedron(found_tetrahedra, voronoi_vertices, points, vdw_radii, simp)
        end_tetrahedra.append(biggest_tetrahedron)
        end_tetrahedra_set.add(biggest_tetrahedron)

    return np.array(end_tetrahedra)


def filter_cavities(cavities, min_depth):
    return [cavity for cavity in cavities if cavity.depth >= min_depth]

def filter_channels_by_bottleneck(cavities, bottleneck):
    for cavity in cavities:
        cavity.channels = [channel for channel in cavity.channels if channel.bottleneck >= bottleneck]

def save_channels_to_pdb(cavities, filename, num_samples=5):
    with open(filename, 'w') as pdb_file:
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
                    pdb_lines.append(f"ATOM  {i:5d}  H   FIL T   1    {x:8.3f}{y:8.3f}{z:8.3f}      {radius:6.2f}\n")
                
                for i in range(1, samples):
                    pdb_lines.append(f"CONECT{i:5d}{i + 1:5d}\n")

                pdb_file.writelines(pdb_lines)
                pdb_file.write("\n")
                atom_index += samples

def calculate_channel_length(centerline_spline):
    t_values = np.linspace(centerline_spline.x[0], centerline_spline.x[-1], len(centerline_spline.x) * 10)
    points = centerline_spline(t_values)
    diffs = np.diff(points, axis=0)
    lengths = np.linalg.norm(diffs, axis=1)
    return np.sum(lengths)

def run_vmd_script(vmd_path, file_path, script_path=None, output_path=None):
    """Executes a VMD script to create a mesh representation of a protein and save it as a .stl file.

    This function runs a VMD (Visual Molecular Dynamics) script using the specified VMD executable and input file.
    It also manages the paths for the script and output, ensuring they are correctly set up before execution.

    :arg vmd_path: Path to the VMD executable. This is required to run the VMD script.
    :type vmd_path: str

    :arg file_path: Path to the input file that will be processed by the VMD script.
    :type file_path: str

    :arg script_path: Path to the VMD script that will be executed. If **None**, defaults to 'script.tcl' in the
        current working directory. The script must be a valid Tcl script for VMD.
    :type script_path: str or None

    :arg output_path: Path where the output file will be saved. If **None**, defaults to 'output/protein.stl' in
        the current working directory. The output file will be created or overwritten at this location.
    :type output_path: str or None

    :returns: None

    This function performs the following steps:
    1. **Path Handling:** Resolves the paths for the VMD script and output file. If not provided, default paths are
       used. Creates the output directory if it does not exist.
    2. **Command Execution:** Constructs the command to run the VMD script with the specified arguments. Executes the
       command using `subprocess.run()` and checks for any errors during execution. The sript creates mesh representation 
       of the protein to be further visualized.
    3. **Error Handling:** Catches and prints errors if the VMD script fails or if any unexpected exceptions occur.

    Note: Ensure that VMD is correctly installed and accessible via the provided `vmd_path`, and that the script and
    file paths are valid for successful execution.
    """
    script_path = Path(script_path or 'script.tcl').resolve()
    if not script_path.is_file():
        raise FileNotFoundError("Script does not exist.")
    
    output_path = Path(output_path or 'output/protein.stl').resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    command = [vmd_path, '-e', str(script_path), '-args', str(file_path), str(output_path)]
    
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"VMD caused an error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

def detect_channels(protein, r1=3, r2=1.25, min_depth=10, bottleneck=1, sparsity=15, output_path=None, visualizer=None, stl_path=None):
    """Detects channels in a protein structure and visualizes or saves the results based on the given criteria.

    This function processes the provided protein structure to identify channels, cavities, and surface, and
    saves the results according to the specified parameters. The detection involves several steps, including
    filtering and visualization.

    :arg protein: The protein structure from which to detect channels. The protein object should support methods for
        selecting atoms and retrieving coordinates.
    :type protein: Protein object

    :arg r1: Radius for determining protein surface. Default is **3**. This parameter defines the cutoff for removing
        too big simplices based on distance from the center of the cavities.
    :type r1: float

    :arg r2: Radius for detecting cavities in protein. Default is **1.25**. This parameter defines the cutoff for
        removing interior small simplices.
    :type r2: float

    :arg min_depth: Minimum depth for filtering cavities. Default is **10**. Only cavities with a depth greater than
        or equal to this value will be retained.
    :type min_depth: int

    :arg bottleneck: Minimum bottleneck size for filtering channels. Default is **1**. Channels with a bottleneck size
        smaller than this value will be filtered out.
    :type bottleneck: float

    :arg sparsity: Controls the quantity of channels. The higher the sparsity, the less channels is detected. 
        Default is **15**. This parameter defines the minimum distance between end points of channels.
    :type sparsity: float

    :arg output_path: Path to save the results as a PDB file. If **None**, the results will not be saved. Default is
        **None**.
    :type output_path: str or None

    :arg visualizer: Type of visualization to perform. Options are **'surface'**, **'channels'**, or **'cavities'**.
        Default is **None**. Determines how the results are visualized.
    :type visualizer: str or None

    :arg stl_path: Path to an STL file for visualizing external meshes. If **None**, default visualization methods
        will be used. Otherwise, the results will be visualized with the protein structure on top. Default is **None**.
    :type stl_path: str or None

    :returns: None

    This function performs the following steps:
    1. **Selection and Filtering:** Selects non-hetero atoms from the protein, calculates van der Waals radii, and performs
       3D Delaunay triangulation and Voronoi tessellation on the coordinates.
    2. **State Management:** Creates and updates different states of the protein structure to filter out simplices based on
       the given radii.
    3. **Surface Layer Calculation:** Determines the surface and second-layer simplices from the filtered results.
    4. **Cavity and Channel Detection:** Finds and filters cavities based on their depth and calculates channels using
       Dijkstra's algorithm.
    5. **Visualization and Saving:** Generates meshes for the detected channels, filters them by bottleneck size, and either
       saves the results to a PDB file or visualizes them based on the specified parameters.

    Note: Ensure that the necessary external libraries and methods are properly imported and available for this function to
    execute correctly.
    """
    protein = protein.select('not hetero')
    coords = protein.getCoords()
    vdw_radii = get_vdw_radii(protein.getElements())
    
    dela = Delaunay(coords)
    voro = Voronoi(coords)
    
    s_prt = State(dela.simplices, dela.neighbors, voro.vertices)
    s_tmp = State(*s_prt.get_state())
    s_prv = State(None, None, None) 
    
    while True:
        s_prv.set_state(*s_tmp.get_state())
        s_tmp.set_state(*delete_simplices3d(coords, *s_tmp.get_state(), vdw_radii, r1, True))
        if s_tmp == s_prv:
            break
    
    s_srf = State(*s_tmp.get_state())
    s_inr = State(*delete_simplices3d(coords, *s_srf.get_state(), vdw_radii, r2, False))
    
    l_first_layer_simp, l_second_layer_simp = surface_layer(s_srf.simp, s_inr.simp, s_srf.neigh)
    s_clr = State(*delete_section(l_first_layer_simp, *s_inr.get_state()))
    
    c_cavities = find_groups(s_clr.neigh)
    c_surface_cavities = get_surface_cavities(c_cavities, s_clr.simp, l_second_layer_simp, s_clr, coords, vdw_radii, sparsity)
    
    find_deepest_tetrahedra(c_surface_cavities, s_clr.neigh)
    c_filtered_cavities = filter_cavities(c_surface_cavities, min_depth)
    merged_cavities = merge_cavities(c_filtered_cavities, s_clr.simp)
    
    for cavity in c_filtered_cavities:
        dijkstra(cavity, *s_clr.get_state(), coords, vdw_radii)
    
    filter_channels_by_bottleneck(c_filtered_cavities, bottleneck)
    channels = [create_mesh_from_spline(*channel.get_splines()) for cavity in c_filtered_cavities for channel in cavity.channels]
    
    no_of_channels = len(channels)
    print(f"Detected {no_of_channels} channels.")
    
    if output_path:
        print(f"Saving results to {output_path}")
        save_channels_to_pdb(c_filtered_cavities, Path(output_path), num_samples=5)
    else:
        print("No output path given.")
    
    if visualizer == 'surface':
        if stl_path:
            visualize_external_grid(coords, s_srf.simp, stl_path)
        else:
            visualize_external_mesh(coords, s_srf.simp)
            
    elif visualizer == 'channels':
        if stl_path:
            visualize_channel(channels, stl_path)
        else:
            visualize_external_grid(coords, s_srf.simp, channel_mesh=channels)
        
    elif visualizer == 'cavities':
        if stl_path:
            visualize_external_grid(coords, merged_cavities, stl_path)
        else:
            visualize_external_mesh(coords, merged_cavities, lines=visualize_external_grid(coords, s_srf.simp, ret_lines=True))




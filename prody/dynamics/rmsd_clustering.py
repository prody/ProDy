# -*- coding: utf-8 -*-

"""
This module defines functions for clustering. It provides tools for statistical analysis of the results
"""

__author__ = 'Konstantinos Theodosiadis'
__credits__ = ['James Krieger', 'Karolina Mikulska-Ruminska', 'Konstantinos Theodosiadis']
__email__ = ['karolamik@fizyka.umk.pl', 'jamesmkrieger@gmail.com', 'konst.theodosiadis@gmail.com']


import numpy as np
from prody.atomic import AtomGroup, Atom, Atomic, Selection, Select
from prody.utilities import getCoords
from prody.proteins import writePDB
from prody.trajectory import writeDCD
from prody.measure import calcTransformation

__all__ = ['alignTrajectory', 'calcPairwiseRMSD', 'calcClusterPopulations', 'getCluster', 'getClusterMedoid',
           'calcClusterStatistics', 'calcAllClusterStatistics', 'showClusterStatisticsTable', 'writeClusters',
           'clusterHierarchical', 'showDendrogram', 'clusterKMedoids', 'clusterDBSCAN', 'showReachabilityPlot']   


def alignTrajectory(atoms, trajectory, align='protein and backbone', select='all'):
    """
    Aligns each trajectory frame to the reference structure and returns a tuple of the 
    reference coordinates and the aligned coordinates of the selected atoms.
    
    The trajectory frames are aligned to the reference structure using the atoms specified by 
    ``align``. After alignment, the coordinates of the atoms specified by ``select`` are 
    extracted for every frame.
    
    
    :arg atoms: reference structure used for the alignment.
    :type atoms: :class:`prody.Atomic`
    
    :arg trajectory: trajectory containing the coordinate sets to align
    :type trajectory: :class:`prody.Trajectory`
    
    :arg align: atom selection used to calculate the alignment transformation.
                Must be a valid ProDy selection string.
                Default is ``"protein and backbone"``
    :type align: str
    
    :arg select: atom selection whose coordinates are returned.
                 Must be a valid ProDy selection string.
                 Default is ``"all"``
    :type select: str
                
    :returns: a tuple containing:
        * ref_coords (numpy.ndarray): coordinates of the selected atoms in the reference structure.
        * aligned_coords (numpy.ndarray): aligned coordinates of the selected atoms for every
        trajectory frame with shape ``(n_frames, n_atoms, 3)``.
    :rtype: tuple(:class:`numpy.ndarray`, :class:`numpy.ndarray`) 
    
    Example usage:
    >>> pdb = parsePDB("structure.pdb")
    >>> dcd = Trajectory("trajectory.dcd")
    >>> ref_coords, aligned_coords = alignTrajectory(pdb, dcd, select='resname IOA')
    """
    
    if trajectory.numAtoms() != atoms.numAtoms():
        raise ValueError("Trajectory atoms count does not match structure atoms count.")
    
    # Save original coordinates to restore in the end
    orig_coords = atoms.getCoords().copy()
    
    # Atoms we want to align
    atom_align = atoms.select(align)
    if atom_align is None:
        raise ValueError(f"No atoms match '{align}' in the structure.")
    ref_align = atom_align.copy()
    
    # Atoms we want to return
    atom_select = atoms.select(select)
    if atom_select is None:
        raise ValueError(f"No atoms match '{select}' in the structure.")
    ref_coords = atom_select.getCoords().copy()
    
    trajectory.link(atoms) # linking trajectory to update coordinates frame-by-frame
    
    n_frames = trajectory.numFrames()
    n_atoms = atom_select.numAtoms()
    trajectory.reset()
    aligned_coords = np.zeros((n_frames, n_atoms, 3))
    
    try:
        for i, frame in enumerate(trajectory):
            trans = calcTransformation(atom_align, ref_align)
            trans.apply(atom_select)
            aligned_coords[i] = atom_select.getCoords()
    finally:
        atoms.setCoords(orig_coords)
        trajectory.reset()
    
    return ref_coords, aligned_coords


def calcPairwiseRMSD(aligned_coords):
    """
    Calculates the frame-to-frame pairwise RMSD matrix using aligned structures
    and returns a symmetric distance matrix.
    Uses a vectorized approach for better efficiency.
    
    
    :arg aligned_coords: aligned coordinates with shape ``(n_frames, n_atoms, 3)``.
                         Recommended to generate them using :func:`alignTrajectory`.
    :type aligned_coords: :class:`numpy.ndarray`
    
    :returns: symmetric matrix containing the pairwise RMSD between all frames with shape 
              ``(n_frames, n_frames)``.
    :rtype: :class:`numpy.ndarray`
    
    Example usage:
    >>> ref_coords, aligned_coords = alignTrajectory(pdb, dcd, select='resname IOA')
    >>> distance_matrix = calcPairwiseRMSD(aligned_coords)
    """
    
    from scipy.spatial.distance import cdist
    
    
    aligned_coords = np.asarray(aligned_coords)
    
    if aligned_coords.ndim != 3 or aligned_coords.shape[2] != 3:
        raise ValueError(f"aligned_coords must have shape (n_frames, n_atoms, 3), but got {aligned_coords.shape}.")
        
    if aligned_coords.shape[0] == 0:
        raise ValueError("aligned_coords contains no frames.")

    if aligned_coords.shape[1] == 0:
        raise ValueError("aligned_coords contains no atoms.")
    
    
    n_frames, n_atoms, _ = aligned_coords.shape
    
    # Flatten frame coordinates (n_frames, n_atoms, 3) -> (n_frames, n_atoms * 3)
    coords_flat = aligned_coords.reshape(n_frames, 3 * n_atoms)
    
    euclidean_dists = cdist(coords_flat, coords_flat, metric = 'euclidean')
    distance_matrix = euclidean_dists / np.sqrt(n_atoms)
    
    return distance_matrix


def calcClusterPopulations(cluster_ids):
    """
    Uses the cluster IDs array, assigning each frame to a cluster, to calculate the population 
    of each cluster and its corresponding percentage.
    
    
    :arg cluster_ids: a one-dimensional array matching each frame to a cluster.
    :type cluster_ids: :class:`numpy.ndarray`
    
    :returns: a dictionary mapping each cluster ID to its population statistics.
              Each value is a dictionary with the keys ``"count"`` and ``"pct"``.
    :rtype: dict
    
    Example usage:
    >>> cluster_ids, _ = clusterHierarchical(distance_matrix)
    >>> populations = calcClusterPopulations(cluster_ids)
    >>> print(f"In cluster 1: {populations[1]['pct']:.2f}% of the objects")
    """
    
    cluster_ids = np.asarray(cluster_ids)
    if cluster_ids.ndim != 1:
        raise ValueError(f"cluster_ids must be a 1D array, but got shape {cluster_ids.shape}")
    if cluster_ids.size == 0:
        raise ValueError("cluster_ids is empty.")
    
    total_frames = len(cluster_ids)
    clusters, frequencies = np.unique(cluster_ids, return_counts=True)
    
    population_data = {}
    for c, n in zip(clusters, frequencies):
        population_data[int(c)] = {'count': int(n), 'pct': (n / total_frames) * 100}
        
    return population_data


def getCluster(cluster_ids, cluster_number):
    """
    Returns the indices of the frames grouped in a specific cluster.
    
    
    :arg cluster_ids: a one-dimensional array assigning each frame to a cluster.
    :type cluster_ids: :class:`numpy.ndarray`
    
    :arg cluster_number: the ID of the cluster whose members we want to return.
    :type cluster_number: int
    
    :returns: an array of the trajectory indices assigned to the specified cluster.
    :rtype: :class:`numpy.ndarray`
    
    Example usage:
    >>> cluster_ids, _ = clusterHierarchical(distance_matrix)
    >>> cluster1_indices = getCluster(cluster_ids, 1)
    """
    
    cluster_ids = np.asarray(cluster_ids)
    if cluster_ids.ndim != 1:
        raise ValueError(f"cluster_ids must be a 1D array, but got shape {cluster_ids.shape}.")
    if cluster_ids.size == 0:
        raise ValueError("cluster_ids is empty.")    
    
    if isinstance(cluster_number, bool) or not isinstance(cluster_number, (int, np.integer)):
        raise TypeError(f"cluster_number must be an integer, but got {type(cluster_number).__name__}.")
    
    cluster_indices = np.where(cluster_ids == cluster_number)[0]
    
    # Fail-safe
    if cluster_indices.size == 0:
        raise ValueError(f'No frames belong to Cluster {cluster_number}')
        
    return cluster_indices


def getClusterMedoid(distance_matrix, cluster_indices=None, cluster_ids=None, cluster_number=None):
    """
    Determines the medoid (representative) of the cluster from the minimum total pairwise 
    distance to all other cluster members.
    Returns a dictionary containing the medoid index within the cluster ("local") and in 
    the full trajectory ("global").
    
    The cluster can be specified either by:
        1. cluster_indices (recommended and can be generated by :func:`getCluster`) or
        2. cluster_ids and cluster_number
     
        
    :arg distance_matrix: a two-dimensional square, symmetric pairwise distance matrix.
    :type distance_matrix: :class:`numpy.ndarray`
    
    :arg cluster_indices: a one-dimensional array of the indices of the cluster members
                          Default is ``None``
    :type cluster_indices: :class:`numpy.ndarray`
    
    :arg cluster_ids: a one-dimensional array of the cluster IDs per element
                      Default is ``None``
    :type cluster_ids: :class:`numpy.ndarray`
    
    :arg cluster_number: the ID (number) of the cluster we want to investigate
                         Default is ``None``
    :type cluster_number: int
    
    :returns: a dictionary of the medoid index, with keys "local" for the index within the
              cluster and "global" for the index within the full trajectory
    :rtype: dict
    
    Example usage:
    >>> cluster_ids, _ = clusterHierarchical(distance_matrix)
    >>> cluster1_indices = getCluster(cluster_ids, 1)
    >>> cluster1_medoid = getClusterMedoid(distance_matrix, cluster1_indices)
    >>> print(f"Medoid of cluster 1 is # {cluster1_medoid['local']} in the cluster.")
    """
    
    cluster_indices = _resolveClusterIndices(cluster_indices=cluster_indices, 
                                             cluster_ids=cluster_ids, 
                                             cluster_number=cluster_number)
    
    distance_matrix = _validateDistanceMatrix(distance_matrix)

    n_frames = distance_matrix.shape[0]
    if np.any((cluster_indices < 0) | (cluster_indices >= n_frames)):
        raise ValueError("cluster_indices contain indices outside the bounds of the distance matrix.")
    
    cluster_distance_matrix = distance_matrix[np.ix_(cluster_indices, cluster_indices)]
    
    sum_dist = cluster_distance_matrix.sum(axis=1)
    medoid_local = int(np.argmin(sum_dist))
    medoid_global = int(cluster_indices[medoid_local])
    
    return {"local": medoid_local, "global": medoid_global}


def calcClusterStatistics(distance_matrix, cluster_indices=None, cluster_medoid=None,
                          cluster_ids=None, cluster_number=None):
    """
    Calculates descriptive statistics for a single cluster.
    Returns a dictionary with the statistics quantity as a key.
    
    The cluster can be specified either by:
        1. cluster_indices (recommended and can be generated by :func:`getCluster`) or
        2. cluster_ids and cluster_number
        
    
    :arg distance_matrix: a two-dimensional square, symmetric pairwise distance matrix.
    :type distance_matrix: :class:`numpy.ndarray`
    
    :arg cluster_indices: a one-dimensional array of the indices of the cluster members
                          Default is ``None``
    :type cluster_indices: :class:`numpy.ndarray`
    
    :arg cluster_medoid: the dictionary of the cluster medoid.
                         Expects the 'cluster_medoid' dictionary generated by the function
                         :func:`getClusterMedoid`.
                         Default is ``None``
    :type cluster_medoid: dict
    
    :arg cluster_ids: a one-dimensional array of the cluster IDs per element
                      Default is ``None``
    :type cluster_ids: :class:`numpy.ndarray`
    
    :arg cluster_number: the ID (number) of the cluster we want to investigate
                         Default is ``None``
    :type cluster_number: int
    
    :returns: a dictionary of statistics quantities with their names as keys
    :rtype: dict
    
    Example usage:
    >>> cluster_ids, _ = clusterHierarchical(distance_matrix)
    >>> cluster1_indices = getCluster(cluster_ids, 1)
    >>> cluster_stats1 = calcClusterStatistics(distance_matrix, cluster1_indices)
    >>> print(f"Mean distance from cluster Medoid: {cluster_stats1['mean']:.2f} ± {cluster_stats1['std']:.2f} [Å]")
    """   

    cluster_indices = _resolveClusterIndices(cluster_indices, cluster_ids, cluster_number)
    population = len(cluster_indices)
    
    distance_matrix = _validateDistanceMatrix(distance_matrix)
    
    if cluster_medoid is None:
        cluster_medoid = getClusterMedoid(distance_matrix, cluster_indices=cluster_indices)
    else:
        if not isinstance(cluster_medoid, dict):
            raise TypeError("cluster_medoid must be a dictionary with keys 'global' and 'local'.")
    
        required_keys = {"global", "local"}
        if not required_keys.issubset(cluster_medoid.keys()):
            missing = required_keys - cluster_medoid.keys()
            raise ValueError(f"Expected keys are 'global' and 'local'.\ncluster_medoid is missing required key(s): {missing}.")   
        
        if cluster_medoid["global"] not in cluster_indices:
            raise ValueError("Index not found.\ncluster_medoid does not belong to the specified cluster.")
            
        if not (0 <= cluster_medoid["local"] < population):
            raise ValueError("Local index is out of bounds.")
            
        if cluster_indices[cluster_medoid["local"]] != cluster_medoid["global"]:
            raise ValueError("cluster_medoid['local'] and cluster_medoid['global'] are inconsistent.")
        
    distance_to_medoid = distance_matrix[cluster_indices, cluster_medoid["global"]]
    total_frames = distance_matrix.shape[0]
    
    stats = {
        "population"    : population,
        "pct"           : (population / total_frames) * 100,
        "medoid_global" : cluster_medoid["global"],
        "medoid_local"  : cluster_medoid["local"],
        "distances"     : distance_to_medoid,
        "mean"          : np.mean(distance_to_medoid),
        "std"           : np.std(distance_to_medoid),
        "median"        : np.median(distance_to_medoid),
        "iqr"           : np.percentile(distance_to_medoid, 75) - np.percentile(distance_to_medoid, 25),
        "p95"           : np.percentile(distance_to_medoid, 95),
        "max"           : np.max(distance_to_medoid)
    }
    
    if cluster_number is not None:
        stats["cluster"] = cluster_number
        
    return stats


def calcAllClusterStatistics(distance_matrix, cluster_ids):
    """
    Calculates descriptive statistics for all clusters.
    Returns a list of dictionaries, each corresponding to a cluster with the statistic quantities
    as keys.
    
    
    :arg distance_matrix: a two-dimensional square, symmetric pairwise distance matrix.
    :type distance_matrix: :class:`numpy.ndarray`
    
    :arg cluster_ids: a one-dimensional array of the cluster IDs per element
    :type cluster_ids: :class:`numpy.ndarray`
    
    :returns: a list of dictionaries, one dictionary for each cluster with the statistic 
              quantities as keys.
    :rtype: list of dict
    
    Example usage:
    >>> cluster_ids, _= clusterHierarchical(distance_matrix, 4)
    >>> all_stats = calcAllClusterStatistics(distance_matrix, cluster_ids)
    """
    
    cluster_ids = np.asarray(cluster_ids)
    if cluster_ids.ndim != 1:
        raise ValueError(f"cluster_ids must be a 1D array, but got shape {cluster_ids.shape}.")
    
    distance_matrix = _validateDistanceMatrix(distance_matrix)
    
    if len(cluster_ids) != distance_matrix.shape[0]:
        raise ValueError(f"Dimension mismatch: cluster_ids length ({len(cluster_ids)}) "
                         f"does not match distance_matrix frames ({distance_matrix.shape[0]}).")
    
    clusters = np.unique(cluster_ids)
    all_stats = []
    
    for c in clusters:
        if c <= 0:
            continue
        
        cluster_stats = calcClusterStatistics(distance_matrix, cluster_ids = cluster_ids, cluster_number = int(c))
        all_stats.append(cluster_stats) 
        
    return all_stats


def _resolveClusterIndices(cluster_indices, cluster_ids, cluster_number):
    """
    Checks whether the input is cluster_indices or cluster_ids and cluster_number, and for the latter returns
    the cluster indices.
    """
   
    has_indices = cluster_indices is not None
    has_ids_info = (cluster_ids is not None) and (cluster_number is not None)
    
    if has_indices == has_ids_info:
        raise ValueError("Provide exactly one of 'cluster_indices' or both 'cluster_ids' and 'cluster_number'.")
    
    if has_indices:
        cluster_indices = np.asarray(cluster_indices)
        
        if cluster_indices.ndim != 1:
            raise ValueError(f"cluster_indices must be a 1D array, but shape {cluster_indices.shape} was given.")
        if cluster_indices.size == 0:
            raise ValueError("cluster_indices must not be empty.")
    else:
        cluster_ids = np.asarray(cluster_ids)
        cluster_indices = getCluster(cluster_ids, cluster_number)  
        
    return cluster_indices


def _validateDistanceMatrix(distance_matrix):
    """
    Validates that a distance matrix is a non-empty square 2D NumPy array.
    """
    
    distance_matrix = np.asarray(distance_matrix)
    if distance_matrix.ndim != 2:
        raise ValueError(f"The distance matrix must be a 2D array, but got shape {distance_matrix.shape}.")
    if distance_matrix.shape[0] != distance_matrix.shape[1]:
        raise ValueError(f"distance_matrix must be square, but got shape {distance_matrix.shape}.")
    if distance_matrix.size == 0:
        raise ValueError("distance_matrix cannot be empty.")
        
    return distance_matrix


def showClusterStatisticsTable(all_stats, dissimilarity="RMSD", units="Å", show=True, **kwargs):
    """
    Prints the cluster statistics in a table, where each column corresponds to a cluster.
    Expects either the 'stats' dictionary generated by :func:`calcClusterStatistics` or
    the 'all_stats' list of dictionaries generated by :func:`calcAllClusterStatistics`.
    
    
    :arg all_stats: a dictionary, or list of dictionaries, with the cluster descriptive statistics
                    Recommended to generate the dictionaries from the functions :func:`calcClusterStatistics`
                    of a single cluster or :func:`calcAllClusterStatistics`.
    :type all_stats: dict, or list of dict
                    
    :arg dissimilarity: the dissimilarity measure
                        Default is `"RMSD"`
    :type dissimilarity: str
    
    :arg units: the units of the dissimilarity measure
                Default is `"Å"`
    :type units: str
    
    :arg show: whether to print the table
               Default is `True`
    :type show: bool
    
    :arg **kwargs: keyword arguments passed directly to the ``tabulate`` function 
    :type **kwargs: dict
    
    :returns: a dictionary containing the raw numeric table data, the row labels, 
              and the headers, with their names as keys
    :rtype: dict
    
    Example usage:
    >>> cluster_ids, _ = clusterHierarchical(distance_matrix)
    >>> all_stats = calcAllClusterStatistics(distance_matrix, cluster_ids)
    >>> showClusterStatisticsTable(all_stats)
    """
    
    if isinstance(all_stats, dict):
        all_stats = [all_stats]
    elif not isinstance(all_stats, list):
        raise TypeError(f"all_stats must be a dict or list of dicts, but got {type(all_stats).__name__}.")
        
    if not all_stats:
        raise ValueError("all_stats cannot be empty.")
    
    required = {"population", "pct", "medoid_global", "medoid_local", "mean", "std",
                "median", "iqr", "p95", "max"}
    
    headers = []
    for i, stat in enumerate(all_stats):
        if not isinstance(stat, dict):
            raise TypeError(f"Expected a dictionary in all_stats, but got {type(stat).__name__} at index {i}.")
        
        cluster_id = stat.get('cluster', i + 1)
        headers.append(f"Cluster {int(cluster_id)}")
        
        missing = required - stat.keys()
        if missing:
            raise ValueError(f"Statistics dictionary is missing required keys: {sorted(missing)}.")
    
    
    floatfmt = kwargs.pop('floatfmt', '.2f')
    
    # The third element marks whether the metric is a floating point number (True) or an exact integer (False)
    metrics = [
        ("population",    "Total Frames",                  False),
        ("pct",           "Population Percentage (%)",     True),
        ("medoid_global", "Medoid Frame (Global)",         False),
        ("medoid_local",  "Medoid Frame (Within Cluster)", False),
        ("mean",          f"Mean {dissimilarity} [{units}]", True),
        ("std",           f"Std {dissimilarity} [{units}]",  True),
        ("median",        f"Median {dissimilarity} [{units}]", True),
        ("iqr",           f"IQR [{units}]",                  True),
        ("p95",           f"95th Percentile [{units}]",      True),
        ("max",           f"Max {dissimilarity} [{units}]",  True)
    ]
    
    table_raw = []
    table_display = []
    row_labels = []

    for key, label, is_float in metrics:
        row_labels.append(label)
        raw_row = [stat[key] for stat in all_stats]
        table_raw.append(raw_row)
        
        fmt = floatfmt if is_float else ".0f"
        
        formatted_row = []
        for val in raw_row:
            if val is None:
                formatted_row.append("N/A")
            else:
                try:
                    formatted_row.append(format(val, fmt))
                except ValueError:
                    formatted_row.append(str(val))
                    
        table_display.append(formatted_row)
        
    if not isinstance(show, bool):
        raise TypeError("show must be a bool.")
    
    if show:
        try:
            from tabulate import tabulate
        except ImportError:
            raise ImportError("The 'tabulate' package is required to display the table. "
                              "Please install it using 'pip install tabulate' or set show=False.")
                              
        kwargs.setdefault('tablefmt', 'fancy_grid')
        kwargs.setdefault('stralign', 'center')
        kwargs.setdefault('numalign', 'center')
        
        print(tabulate(table_display, headers=headers, showindex=row_labels, **kwargs))
    
    return {"table": table_raw, "row_labels": row_labels, "headers": headers}


def writeClusters(atoms, trajectory, distance_matrix, cluster_ids, write_dcd=True,
                  align ="protein and backbone", system="system", tag ="cluster"):
    """
    Aligns the trajectory once, then exports representative medoid structures as PDB files and 
    cluster-specific DCD trajectories.
    Returns the name of all exported .pdb and .dcd files.
    
    Note: This function loads all aligned coordinates into memory. It is highly optimized 
    for speed, provided the trajectory fits within available system RAM.
    
    :arg atoms: reference structure used for the alignment.
    :type atoms: :class:`prody.Atomic`
    
    :arg trajectory: trajectory containing the coordinate sets to align
    :type trajectory: :class:`prody.Trajectory`
    
    :arg distance_matrix: square, symmetric matrix of pairwise distances between all objects
    :type distance_matrix: :class:`numpy.ndarray`
    
    :arg cluster_ids: a one-dimensional array of the cluster IDs per element
    :type cluster_ids: :class:`numpy.ndarray`   
    
    :arg write_dcd: determines whether to save a .dcd file of the cluster frames
                    Default is ``True``
    :type write_dcd: bool
    
    :arg align: atom selection used to calculate the alignment transformation.
                Must be a valid ProDy selection string.
                Default is ``"protein and backbone"``
    :type align: str
    
    :arg system: the name of the system under investigation
                 Default is ``"system"``
    :type system: str
    
    :arg tag: the name of the clustering method
              Default is ``"cluster"``
    :type tag: str
    
    :returns: list of exported filenames
    :rtype: list[str]
    
    Example usage:
    >>> pdb = parsePDB("structure.pdb")
    >>> dcd = Trajectory("trajectory.dcd")
    >>> _, aligned_coords = alignTrajectory(pdb, dcd, select='resname IOA')
    >>> distance_matrix = calcPairwiseRMSD(aligned_coords)
    >>> cluster_ids, _ = clusterHierarchical(distance_matrix)
    >>> writeClusters(pdb, dcd, distance_matrix, cluster_ids,
                      system="type1_RUN23", tag="hier")
    """
    
    distance_matrix = _validateDistanceMatrix(distance_matrix)
    
    cluster_ids = np.asarray(cluster_ids)
    if cluster_ids.ndim != 1:
        raise ValueError(f"cluster_ids must be a 1D array, but got shape {cluster_ids.shape}.")
    if cluster_ids.size == 0:
        raise ValueError("cluster_ids is empty.")
    if len(cluster_ids) != trajectory.numFrames():
        raise ValueError("cluster_ids must have one entry per trajectory frame.")
        
    if not isinstance(write_dcd, bool):
        raise TypeError(f"write_dcd must be a bool, but got {type(write_dcd).__name__}")

    clusters = np.unique(cluster_ids)
    clusters = clusters[clusters > 0]
    num_clusters = clusters.size
    if num_clusters == 0:
        raise ValueError("No clusters were found. All frames are labeled as noise")
    
    exported_files = []
    
    # NOTE: Loads the entire aligned trajectory into memory
    _, aligned_coords = alignTrajectory(atoms, trajectory, align=align, select="all")
    
    for cluster in clusters:
        cluster_indices = getCluster(cluster_ids, cluster)
        medoid = getClusterMedoid(distance_matrix, cluster_indices = cluster_indices)

        cluster_coords = aligned_coords[cluster_indices]

        cluster_atoms = atoms.copy()
        cluster_atoms.setCoords(cluster_coords[0])

        if len(cluster_coords) > 1:
            cluster_atoms.addCoordset(cluster_coords[1:])

        if write_dcd:
            dcd_filename = f"{system}_{tag}_n{num_clusters}_cluster{cluster}.dcd"
            writeDCD(dcd_filename, cluster_atoms)
            exported_files.append(dcd_filename)

        medoid_atoms = atoms.copy()
        medoid_atoms.setCoords(aligned_coords[medoid["global"]])

        pdb_filename = f"{system}_{tag}_n{num_clusters}_cluster{cluster}_medoid.pdb"
        writePDB(pdb_filename, medoid_atoms)
        exported_files.append(pdb_filename)

    return exported_files
    

def clusterHierarchical(distance_matrix, method='average', cutoff=None):
    """
    Performs bottom-up hierarchical clustering from a pairwise distance matrix.
    
    Note that the resulting cluster IDs are 1-based (starting at 1, not 0).
    
    If cutoff is ``None`` or ``'auto'``, it is automatically chosen as the midpoint 
    of the largest gap between consecutive linkage distances.
    
    Recommendation: Plot the dendrogram using showDendrogram() and choose the cutoff 
                    manually whenever possible.                
        
        
    :arg distance_matrix: either a one-dimensional condensed distance matrix or 
                          a two-dimensional pairwise distance matrix.
    :type distance_matrix: :class:`numpy.ndarray`
    
    :arg method: linkage criterion used when constructing the linkage matrix.
                 Default is ``'average'``
    :type method: str
    
    :arg cutoff: the cutoff distance that determines the number of clusters.
                 If ``None`` or ``'auto'``, the midpoint of the largest linkage gap
                 is calculated automatically.
                 Default is ``None``
    :type cutoff: float, str, None
    
    :returns:
        * a one-dimensional array containing the cluster ID for each object.
          IDs are 1-indexed (1 to number of clusters for the cutoff)
        * the hierarchical clustering linkage matrix
    :rtype: tuple(:class:`numpy.ndarray`, :class:`numpy.ndarray`)
    
    
    .. seealso::
        :func:`clusterMatrix`
            Hierarchical clustering of a distance matrix using SciPy.

    
    Example usage:
    >>> distance_matrix = calcPairwiseRMSD(aligned_coords)
    >>> cluster_ids, linkage_matrix = clusterHierarchical(distance_matrix, cutoff='auto')
    """

    from scipy.cluster.hierarchy import linkage, fcluster


    condensed_distance_matrix = _condenseDistanceMatrix(distance_matrix)    
    linkage_matrix = linkage(condensed_distance_matrix, method=method)
    
    if cutoff is None or cutoff == "auto":
        cutoff = _calcAutoCutoff(linkage_matrix)
    elif isinstance(cutoff, str):
        raise ValueError(f"Invalid string for cutoff: '{cutoff}'. Use 'auto', None, or a numeric value.")
        
    cluster_ids = fcluster(linkage_matrix, t=cutoff, criterion='distance')
    
    return cluster_ids, linkage_matrix


def showDendrogram(distance_matrix=None, linkage_matrix=None, *args, **kwargs):
    """
    Plots a hierarchical clustering dendrogram using either a pairwise distance matrix 
    or a pre-computed linkage matrix.

    By default, the complete dendrogram is shown. To display a truncated dendrogram, 
    pass ``truncate_mode`` and ``p`` via kwargs.

    Exactly one of ``distance_matrix`` or ``linkage_matrix`` must be provided.
    
    
    :arg distance_matrix: one-dimensional condensed or two-dimensional square pairwise 
                          distance matrix. 
                          Default is ``None``
    :type distance_matrix: :class:`numpy.ndarray`

    :arg linkage_matrix: pre-computed hierarchical clustering linkage matrix.
                         Default is ``None``
    :type linkage_matrix: :class:`numpy.ndarray`
    
    :arg *args: positional arguments passed directly to SciPy's ``dendrogram`` function.
    :type *args: tuple

    :arg method: linkage criterion used when constructing the linkage matrix from 
                 ``distance_matrix``. (Passed via kwargs).
                 Default is ``"average"``
    :type method: str

    :arg cutoff: dendrogram color threshold. If ``"auto"``, the cutoff is placed at the midpoint
                 of the largest linkage gap. (Passed via kwargs).
                 Default is ``None``
    :type cutoff: float, str
    
    :arg title: title of the generated plot. (Passed via kwargs).
                Default is ``"Hierarchical Clustering Dendrogram"``
    :type title: str
    
    :arg xlabel: label for the x-axis. (Passed via kwargs).
                 Default is ``"Index"``
    :type xlabel: str
    
    :arg ylabel: label for the y-axis. (Passed via kwargs).
                 Default is ``"Distance"``
    :type ylabel: str

    :arg ax: axes on which to draw the plot. (Passed via kwargs).
             Default is ``None`` and the current axes are used.
    :type ax: :class:`matplotlib.axes.Axes`
    
    :arg **kwargs: additional keyword arguments passed to SciPy's ``dendrogram`` 
                   (e.g., ``truncate_mode``, ``p``, ``color_threshold``, ``no_labels``).
    :type **kwargs: dict

    :returns: the matplotlib axes and the linkage matrix.
    :rtype: tuple(:class:`matplotlib.axes.Axes`, :class:`numpy.ndarray`)

    Example usage:
    >>> import matplotlib.pyplot as plt
    >>> distance_matrix = calcPairwiseRMSD(aligned_coords)
    >>> plt.figure(figsize=(12, 8))
    >>> ax, linkage = showDendrogram(distance_matrix, cutoff='auto', truncate_mode='lastp', p=30)
    >>> plt.show()
    """
    
    import matplotlib.pyplot as plt
    from scipy.cluster.hierarchy import dendrogram, linkage


    if (distance_matrix is None) == (linkage_matrix is None):
        raise ValueError("Provide exactly one of 'distance_matrix' or 'linkage_matrix'.")

    method = kwargs.pop('method', 'average')
    cutoff = kwargs.pop('cutoff', None)

    if linkage_matrix is None:
        condensed_distance_matrix = _condenseDistanceMatrix(distance_matrix)
        linkage_matrix = linkage(condensed_distance_matrix, method=method)
    else:
        linkage_matrix = np.asarray(linkage_matrix)
        if linkage_matrix.ndim != 2 or linkage_matrix.shape[1] != 4 or linkage_matrix.shape[0] == 0:
            raise ValueError("linkage_matrix must be a non-empty array with shape (n_samples-1, 4).")

    if cutoff == "auto":
        cutoff = _calcAutoCutoff(linkage_matrix)
    elif isinstance(cutoff, str):
        raise ValueError(f"Invalid string for cutoff: '{cutoff}'. Use 'auto' or a numeric value.")

    ax = kwargs.pop('ax', None)
    if ax is None:
        ax = plt.gca()

    title = kwargs.pop('title', "Hierarchical Clustering Dendrogram")
    xlabel = kwargs.pop('xlabel', "Index")
    ylabel = kwargs.pop('ylabel', "Distance")

    # SciPy Dendogram Defaults
    kwargs.setdefault('truncate_mode', None)
    kwargs.setdefault('p', 30)
    kwargs.setdefault('no_labels', True)
    
    if cutoff is not None:
        kwargs.setdefault('color_threshold', cutoff)

    with plt.rc_context({"lines.linewidth": 0.6}):
        dendrogram(linkage_matrix, *args, ax=ax, **kwargs)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    if cutoff is not None:
        ax.axhline(cutoff, color="black", linestyle="--", linewidth=1.2, label=f"Cutoff ({cutoff:.2f})")
        ax.legend()

    return ax, linkage_matrix


def _calcAutoCutoff(linkage_matrix):
    """
    Calculates the optimal distance cutoff based on the midpoint of the largest gap
    in the linkage matrix.
    """
    
    linkage_matrix = np.asarray(linkage_matrix)
    merge_dist = linkage_matrix[:, 2] # locations of the merges
    
    if len(merge_dist) == 0:
        raise ValueError("At least two frames are required for clustering")
    elif len(merge_dist) == 1:
        return merge_dist[0] + 1.0

    gaps = np.diff(merge_dist)
    largest_gap_idx = np.argmax(gaps)
    
    return (merge_dist[largest_gap_idx] + merge_dist[largest_gap_idx + 1]) / 2.0


def _condenseDistanceMatrix(distance_matrix):
    """
    Converts a pairwise distance matrix into condensed form.

    Accepts either:
        * a condensed one-dimensional distance matrix, or
        * a square two-dimensional distance matrix.
    
    Returns the condensed distance matrix suitable for scipy.cluster.hierarchy.linkage.
    """
    from scipy.spatial.distance import squareform
    
    
    distance_matrix = np.asarray(distance_matrix)
    
    if distance_matrix.size == 0:
        raise ValueError("distance_matrix cannot be empty.")
    
    if distance_matrix.ndim == 2:
        if distance_matrix.shape[0] != distance_matrix.shape[1]:
            raise ValueError(f"distance_matrix must be square, but got shape {distance_matrix.shape}.")
        condensed_distance_matrix = squareform(distance_matrix)
    elif distance_matrix.ndim == 1:
        condensed_distance_matrix = distance_matrix
    else:
        raise ValueError("distance_matrix must be either a condensed 1D array or a square 2D array.")
        
    return condensed_distance_matrix


def clusterKMedoids(distance_matrix, k, method='alternate', method_sklearn='pam', 
                    initial_medoids=None, seed=None, max_iter=100, n_init=10):
    """
    Performs K-Medoids clustering using various algorithms.

    This function acts as a facade, routing the clustering task to the specified backend 
    ('alternate', 'pam', or 'sklearn'). 'alternate' is generally faster, while 'pam' typically
    finds solutions with lower cost. 'sklearn' needs the sklearn_extra.cluster module to be installed.
    
    Note that the resulting cluster IDs are 1-based (starting at 1, not 0).


    :arg distance_matrix: square, symmetric matrix of pairwise distances between all objects
    :type distance_matrix: :class:`numpy.ndarray`
    
    :arg k: prespecified number of clusters to form
    :type k: int
    
    :arg method: the clustering algorithm to use.
                 Options are 'alternate' (custom Alternating Medoids), 'pam' (custom PAM), 'sklearn' (sklearn_extra)
                 Default is 'alternate'
    :type method: str
    
    :arg method_sklearn: the specific method to pass to sklearn_extra if method = 'sklearn' is chosen.
                         Options are 'pam', 'alternate'
                         Default is 'pam'
    :type method_sklearn: str
    
    :arg initial_medoids: one-dimensional array of indices to use as starting medoids.
                          Default is ``None`` and the starting medoids are picked randomly.
                          It is not supported for the 'sklearn' method.
    :type initial_medoids: :class:`numpy.ndarray`
    
    :arg seed: random seed for reproducibility
               Default is ``None``
    :type seed: int
    
    :arg max_iter: maximum number of iterations per a single run
                   Default is ``100``
    :type max_iter: int
    
    :arg n_init: number of times the algorithm will be run with different initial medoids
                 and the best result with the lowest cost is returned
    
    :returns: a tuple of:
        * an one-dimensional array containing the cluster ID for each object.
          IDs are 1-indexed (1 to k)
        * a one-dimensional array of shape (k,) containing the indices of the final cluster
          medoids
        * the final sum of distances from each point to its nearest medoid
    :rtype: tuple(:class:`numpy.ndarray`, :class:`numpy.ndarray`, float)
    
    
    .. seealso::
        :func:`prody.utilities.catchall.calcKmedoidClusters`
            Direct K-Medoids clustering on coordinate sets using ``scikit-learn-extra``.
    
    
    Example usage:
    >>> distance_matrix = calcPairwiseRMSD(aligned_coords)
    >>> cluster_ids, medoids, _ = clusterKMedoids(distance_matrix, 4, method = 'alternate', 
                                                        seed = 42, n_init = 30)
    """
    
    distance_matrix = _validateDistanceMatrix(distance_matrix)
    n = distance_matrix.shape[0]
    
    if not isinstance(k, (int, np.integer)):
        raise TypeError(f"k must be an integer, got {type(k).__name__}")
    if k <= 0 or k > n:
        raise ValueError(f"k must be between 1 and {n} (the number of points).")
        
    if not isinstance(max_iter, (int, np.integer)) or max_iter <= 0:
        raise ValueError("max_iter must be a positive integer.")

    if not isinstance(n_init, (int, np.integer)) or n_init <= 0:
        raise ValueError("n_init must be a positive integer.")
        
    if seed is not None and not isinstance(seed, (int, np.integer)):
        raise TypeError("seed must be an integer or None.")
        
    if initial_medoids is not None:
        if method == 'sklearn':
            raise ValueError("The 'sklearn' method does not support custom 'initial_medoids'. Use 'alternate' or 'pam' instead.")
                             
        initial_medoids = np.asarray(initial_medoids, dtype=int)
        if initial_medoids.ndim != 1:
            raise ValueError("initial_medoids must be a 1D array.")
        if len(initial_medoids) != k:
            raise ValueError(f"Number of initial medoids ({len(initial_medoids)}) must equal k ({k}).")
        if len(np.unique(initial_medoids)) != k:
            raise ValueError("initial_medoids must be unique.")
        if np.any((initial_medoids < 0) | (initial_medoids >= n)):
            raise ValueError(f"initial_medoids contain indices out of bounds (must be 0 to {n-1}).")
            
        n_init = 1

    if method == 'alternate':
        return _clusterKMedoidsAlternating(distance_matrix, k, initial_medoids, seed, max_iter, n_init)
    elif method == 'pam':
        return _clusterKMedoidsPAM(distance_matrix, k, initial_medoids, seed, max_iter, n_init)
    elif method == 'sklearn':
        if method_sklearn not in ['alternate', 'pam']:
            raise ValueError(f"method_sklearn '{method_sklearn}' is not valid. Options: 'alternate', 'pam'")
        return _clusterKMedoidsSklearn(distance_matrix, k, method_sklearn, seed, max_iter, n_init)
    else:
        raise ValueError(f"Method '{method}' is not valid. Options: 'alternate', 'pam', 'sklearn'")


def _cost(distance_matrix, medoids):
    """Calculates the sum of distances from all points to their nearest medoid."""
    distances = distance_matrix[:, medoids] 
    closest_distance = np.min(distances, axis=1) 
    return np.sum(closest_distance)


def _clusterKMedoidsAlternating(distance_matrix, k, initial_medoids, seed, max_iter, n_init):
    """K-Medoids clustering with an Alternating Medoids algorithm"""
    
    n = distance_matrix.shape[0]
    rng = np.random.default_rng(seed)

    best_cost = np.inf
    best_medoids = None
    
    for run in range(n_init): 
        if initial_medoids is None:
            medoids = np.asarray(rng.choice(n, size = k, replace = False))
        else:
            medoids = initial_medoids.copy()

        for iteration in range(max_iter):
            distances_to_medoids = distance_matrix[:, medoids] 
            cluster_ids = np.argmin(distances_to_medoids, axis = 1) 

            new_medoids = np.zeros_like(medoids)
            
            for i in range(len(medoids)):
                cluster_members = np.where(cluster_ids == i)[0]
                
                # Handles empty clusters by retaining the old medoid
                if len(cluster_members) == 0:
                    new_medoids[i] = medoids[i]
                    continue
                
                intra_cluster_distances = distance_matrix[np.ix_(cluster_members, cluster_members)]
                sum_distances = intra_cluster_distances.sum(axis=1)
                best_medoid_idx_in_cluster = np.argmin(sum_distances)
                new_medoids[i] = cluster_members[best_medoid_idx_in_cluster]
            
            current_cost = _cost(distance_matrix, new_medoids)

            if np.array_equal(medoids, new_medoids):
                break
                
            medoids = new_medoids
    
        if current_cost < best_cost:
            best_cost = current_cost
            best_medoids = medoids.copy()

    distances = distance_matrix[:, best_medoids]
    cluster_ids = np.argmin(distances, axis=1) + 1 

    return cluster_ids, best_medoids, best_cost


def _clusterKMedoidsPAM(distance_matrix, k, initial_medoids, seed, max_iter, n_init):
    """K-Medoids clustering with PAM (Partitioning Around Medoids)"""
    
    n = distance_matrix.shape[0]
    rng = np.random.default_rng(seed)

    best_cost = np.inf
    best_medoids = None

    for run in range(n_init): 
        if initial_medoids is None:
            medoids = np.asarray(rng.choice(n, size=k, replace=False))
        else:
            medoids = initial_medoids.copy()
            
        # Using a set for removal/addition operations
        non_medoids = set(i for i in range(n) if i not in medoids)
        current_cost = _cost(distance_matrix, medoids)

        for iteration in range(max_iter):
            best_cost_swap = current_cost
            best_swap = None
            
            for medoid_idx, old_medoid in enumerate(medoids):
                for new_medoid in non_medoids:
                    candidate = medoids.copy()
                    candidate[medoid_idx] = new_medoid
                    candidate_cost = _cost(distance_matrix, candidate)

                    if candidate_cost < best_cost_swap:
                        best_cost_swap = candidate_cost
                        best_swap = (medoid_idx, old_medoid, new_medoid)

            if best_swap is None:
                break

            idx, old_medoid, new_medoid = best_swap
            medoids[idx] = new_medoid
            non_medoids.remove(new_medoid)
            non_medoids.add(old_medoid)
            current_cost = best_cost_swap

        if current_cost < best_cost:
            best_cost = current_cost
            best_medoids = medoids.copy()

    distances = distance_matrix[:, best_medoids]
    cluster_ids = np.argmin(distances, axis=1) + 1 

    return cluster_ids, best_medoids, best_cost


def _clusterKMedoidsSklearn(distance_matrix, k, method_sklearn, seed, max_iter, n_init):
    """K-Medoids clustering with sklearn_extra"""
    try:
        from sklearn_extra.cluster import KMedoids
    except ImportError:
        raise ImportError("The 'sklearn_extra' package is required for this K-Medoids approach. "
                          "Please install it using 'pip install scikit-learn-extra'.")
    
    
    best_cost = np.inf
    best_labels = None
    best_medoids = None

    rng = np.random.default_rng(seed)

    for run in range(n_init): 
        kmedoids = KMedoids(n_clusters=k, metric='precomputed', method=method_sklearn, init='random', 
                            max_iter=max_iter, random_state=int(rng.integers(0, 1000000)))

        kmedoids.fit(distance_matrix)
        medoids = kmedoids.medoid_indices_

        cost = np.sum(np.min(distance_matrix[:, medoids], axis=1))

        if cost < best_cost:
            best_cost = cost
            best_labels = kmedoids.labels_ + 1 
            best_medoids = medoids.copy()

    return best_labels, best_medoids, best_cost


def clusterDBSCAN(distance_matrix, eps=None, minPts=None, method='custom'):
    """
    Performs DBSCAN clustering using various algorithms.

    This function acts as a facade, routing the clustering task to the specified backend 
    ('custom' or 'sklearn'). 'custom' uses a built-in implementation without additional 
    dependencies beyond NumPy. 'sklearn' needs the sklearn.cluster module to be installed.
    
    Note that the resulting cluster IDs are 1-based (starting at 1, not 0).
    Noise points are labeled -1.


    :arg distance_matrix: square, symmetric matrix of pairwise distances between all objects
    :type distance_matrix: :class:`numpy.ndarray`
    
    :arg eps: the "radius" of the neighborhood within which we count neighbors
              Default is ``None`` and automatically the median of the distances 
              in the distance matrix is used. 
              Ideally use :func:`showReachabilityPlot` to determine
              manually the most suitable eps.
    :type eps: float
    
    :arg minPts: the minimum number of neighbors required for a point to be considered 
                 a core point
                 Default is ``None`` and automatically the 5% of the total objects,
                 or for less than 20 2 is used.
                 Ideally choose manually the most suitable minPts
    :type minPts: int
    
    :arg method: the clustering algorithm to use
                 Options are 'custom' and 'sklearn'
                 Default is 'custom' because it needs no module installation
    :type method: str
    
    :returns: a tuple of:
        * a one-dimensional array containing the cluster ID for each object.
          IDs are 1-indexed and noise corresponds to -1.
        * a one-dimensional array of the frame indices corresponding to noise
    :rtype: tuple(:class:`numpy.ndarray`, :class:`numpy.ndarray`)
    
    Example usage:
    >>> distance_matrix = calcPairwiseRMSD(aligned_coords)
    >>> cluster_ids, _ = clusterDBSCAN(distance_matrix, eps = 1.8, minPts = 30)
    """
    
    distance_matrix = _validateDistanceMatrix(distance_matrix)
    
    if eps is None:
       eps = _calcAutoEps(distance_matrix)
    _validateEps(eps, distance_matrix)
        
    if minPts is None:
        minPts = _calcAutoMinPts(distance_matrix)
    _validateMinPts(minPts, distance_matrix)
    
    if method == 'custom':
        return _clusterDBSCANCustom(distance_matrix, eps, minPts)
    elif method == 'sklearn':
        return _clusterDBSCANSklearn(distance_matrix, eps, minPts)
    else:
        raise ValueError(f"method can be either 'custom' or 'sklearn', but got {method}.")


def _clusterDBSCANCustom(distance_matrix, eps, minPts):
    """DBSCAN clustering with custom algorithm"""
    total_points = distance_matrix.shape[0]
    
    labels = np.zeros(total_points, dtype=int)
    cluster_id = 0
    
    for p in range(total_points):

        if labels[p] != 0:
            continue
            
        neighbors, = np.where(distance_matrix[p] <= eps)
        
        if len(neighbors) < minPts:
            labels[p] = -1
        else:
            cluster_id += 1
            labels[p] = cluster_id
            
            candidate_set = [n for n in neighbors if n != p]
            while candidate_set:
                q = candidate_set.pop() 
                
                if labels[q] == -1:
                    labels[q] = cluster_id
                    
                if labels[q] != 0:
                    continue
                    
                labels[q] = cluster_id
                
                q_neighbors, = np.where(distance_matrix[q] <= eps)
                if len(q_neighbors) >= minPts:
                    for n in q_neighbors:
                        if labels[n] == 0:
                            candidate_set.append(n)
                        elif labels[n] == -1:
                            labels[n] = cluster_id
                            
    noise_frames, = np.where(labels == -1)
    
    return labels, noise_frames        


def _clusterDBSCANSklearn(distance_matrix, eps, minPts):
    """DBSCAN clustering with sklearn"""    
    
    try:
        from sklearn.cluster import DBSCAN
    except ImportError:
        raise ImportError("The 'sklearn' package is required for this DBSCAN approach. "
                          "Please install it using 'pip install scikit-learn'.")
    
    
    dbscan = DBSCAN(eps = eps, min_samples = minPts, metric = "precomputed")
    labels = dbscan.fit_predict(distance_matrix)
    
    cluster_ids = np.copy(labels)
    cluster_ids[cluster_ids >= 0] += 1

    noise_frames, = np.where(labels == -1)
    
    return cluster_ids, noise_frames


def _calcAutoEps(distance_matrix):
    """Automatically determine the DBSCAN eps parameter"""
    
    import warnings
    
    
    pairwise_distances = distance_matrix[np.triu_indices_from(distance_matrix, k = 1)]
    eps = float(np.median(pairwise_distances))
    warnings.warn(f"No eps provided. Automatically chosen at {eps:.3f}. Ideally provide your own value.")
    return eps


def _validateEps(eps, distance_matrix):
    """Validate the DBSCAN eps parameter"""
    
    import warnings
    
    
    maxDistance = np.max(distance_matrix)    
        
    if not isinstance(eps, (float, int, np.floating, np.integer)):
        raise TypeError(f"eps must be a numeric value, but got {type(eps).__name__}")
    
    if eps <= 0:
        raise ValueError("eps must be positive")
    elif eps > maxDistance:
        warnings.warn(f"eps ({eps}) is greater than the maximum pairwise distance ({maxDistance:.3f}).\n"
                      "All frames will be clustered together with no noise.")


def _calcAutoMinPts(distance_matrix):
    """ Automatically determine the minPts parameter"""
    
    import warnings
    
    
    total_points = distance_matrix.shape[0]
    minPts = max(2, int(total_points // 20))
    warnings.warn(f"No minPts provided. Automatically chosen at {minPts}. Ideally provide your own value.")
    return minPts
  
  
def _validateMinPts(minPts, distance_matrix):
    """Validate minPts parameter"""
    total_points = distance_matrix.shape[0]    
        
    if not isinstance(minPts, (int, np.integer)):
        raise TypeError(f"minPts must be a positive integer, but got {type(minPts).__name__}")
        
    if minPts <= 0 or minPts > total_points:
        raise ValueError(f"minPts must be between 1 and {total_points}.")
        

def showReachabilityPlot(distance_matrix, *args, minPts=None, method='custom', eps=None, 
                         fill=True, colors=None, **kwargs):
    """
    Plots the reachability plot using a simplified OPTICS algorithm.
    The reachability plot should be used to determine the most suitable eps 
    parameter for DBSCAN.
    
    This function acts as a facade, routing the ordering task to the specified backend 
    ('custom' or 'sklearn'). 'custom' uses a built-in implementation without additional 
    dependencies beyond NumPy. 'sklearn' needs the sklearn.cluster module to be installed.
    
    
    :arg distance_matrix: a two-dimensional square, symmetric pairwise distance matrix.
    :type distance_matrix: :class:`numpy.ndarray`
    
    :arg *args: positional arguments passed directly to Matplotlib's ``plot`` function.
    :type *args: tuple
    
    :arg minPts: the minimum number of neighbors required for a point to be considered 
                 a core point.
                 Default is ``None`` and automatically the 5% of the total objects,
                 or for less than 20, 2 is used.
    :type minPts: int
    
    :arg method: the OPTICS algorithm to use ('custom' or 'sklearn').
                 Default is 'custom'.
    :type method: str
    
    :arg eps: the "radius" of the neighborhood within which we count neighbors.
              If `'auto'`, the median of the pairwise distances is used.
              Default is ``None``.
    :type eps: float, str
    
    :arg fill: whether to fill the area beneath the curve and the valleys below eps.
               Default is ``True``.
    :type fill: bool
    
    :arg colors: custom color mapping for the valleys (list, dict, or colormap name).
                 Default is ``None`` and uses Matplotlib's tab10 palette.
    :type colors: list, dict, str, or None
    
    :arg **kwargs: additional keyword arguments passed to Matplotlib's ``plot`` function 
    :type **kwargs: dict
    
    :returns: the Matplotlib axes containing the plot.
    :rtype: :class:`matplotlib.axes.Axes`
    
    Example usage:
    >>> import matplotlib.pyplot as plt
    >>> plt.figure()
    >>> showReachabilityPlot(distance_matrix, minPts = 20, eps = 2.8)
    >>> plt.show()
    """
    
    import matplotlib.pyplot as plt

    
    distance_matrix = _validateDistanceMatrix(distance_matrix)
    
    if minPts is None:
        minPts = _calcAutoMinPts(distance_matrix) 
    _validateMinPts(minPts, distance_matrix)
    
    if not isinstance(fill, bool):
        raise TypeError(f"fill must be a bool, but got {type(fill).__name__}")
    
    if method == 'custom':
        reachability, ordering = _orderOPTICSCustom(distance_matrix, minPts)
    elif method == 'sklearn':
        reachability, ordering = _orderOPTICSSklearn(distance_matrix, minPts)
    else:
        raise ValueError(f"method must be either 'custom' or 'sklearn', but got {method}")
    
    ax = kwargs.pop('ax', None)
    if ax is None:
        ax = plt.gca()
        
    title = kwargs.pop('title', "OPTICS Reachability Plot")
    xlabel = kwargs.pop('xlabel', "Frames (Sorted by OPTICS)")
    ylabel = kwargs.pop('ylabel', "Reachability Distance")
    
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    y = reachability[ordering]
    x = np.arange(len(y))
    
    kwargs.setdefault('color', "#36454F")
    kwargs.setdefault('lw', 1.5)
    
    ax.plot(x, y, *args, **kwargs)
    
    if fill:
        ax.fill_between(x, 0, y, color='black', alpha=0.4)
    
    if eps == 'auto':
        eps = _calcAutoEps(distance_matrix)
    
    if eps is not None:
        _validateEps(eps, distance_matrix)
        
        line_width = kwargs.get('lw', kwargs.get('linewidth', 1.5))
        ax.axhline(y=eps, color='black', linestyle='--', linewidth=line_width)
        
        if fill:
            below_eps_mask = (y <= eps)
            segments = []
            start = None

            for i, below in enumerate(below_eps_mask):
                if below and start is None:
                    start = i
                elif not below and start is not None:
                    segments.append((start, i))
                    start = None
            
            if start is not None:
                segments.append((start, len(below_eps_mask)))
            
            if colors is None:
                cmap = plt.get_cmap('tab10')
            elif isinstance(colors, str):
                cmap = plt.get_cmap(colors)
            else:
                cmap = colors

            for i, (start, end) in enumerate(segments):
                if isinstance(cmap, dict):
                    valley_color = cmap.get(i, 'gray')
                elif isinstance(cmap, (list, tuple)):
                    valley_color = cmap[i % len(cmap)]
                else:
                    num_colors = getattr(cmap, 'N', 10) 
                    valley_color = cmap(i % num_colors)
                
                ax.fill_between(x[start:end], y[start:end], eps, color=valley_color)
    
    if 'label' in kwargs:
        ax.legend()
        
    ax.grid(axis='y', linestyle='--', linewidth=0.8, alpha=0.3)

    return ax    

    
def _orderOPTICSCustom(distance_matrix, minPts):
    """OPTICS algorithm with built-in modules and NumPy"""
    
    import heapq
    
    
    total_points = distance_matrix.shape[0]
    
    sorted_distances = np.sort(distance_matrix, axis=1)
    core_distances = sorted_distances[:, minPts - 1] 
    
    reachability = np.full(total_points, np.inf)
    processed = np.zeros(total_points, dtype=bool) 
    ordering = []
    
    def _updateSeeds(idx):
        new_reaches = np.maximum(core_distances[idx], distance_matrix[idx, :]) # Reachability distance definition
        update_mask = (~processed) & (new_reaches < reachability)
        points_to_update, = np.where(update_mask)
        reachability[update_mask] = new_reaches[update_mask]
        
        for j in points_to_update:
            heapq.heappush(seeds, (reachability[j], j))
    
    for i in range(total_points):
        if processed[i]:
            continue
        
        processed[i] = True
        ordering.append(i)
        seeds = []
        
        _updateSeeds(i)

        while seeds:
            current_reach, q = heapq.heappop(seeds)
            
            if processed[q]:
                continue
    
            processed[q] = True
            ordering.append(q)
            
            _updateSeeds(q)
                
    return reachability, np.array(ordering)


def _orderOPTICSSklearn(distance_matrix, minPts):
    """OPTICS algorithm with sklearn"""
    
    try:
        from sklearn.cluster import OPTICS
    except ImportError:
        raise ImportError("The 'sklearn' package is required for this OPTICS approach. "
                          "Please install it using 'pip install scikit-learn'.")
        
    optics = OPTICS(min_samples = minPts, metric = 'precomputed')
    optics.fit(distance_matrix)
    reachability = optics.reachability_
    ordering = optics.ordering_
    
    return reachability, ordering
    
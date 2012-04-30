"""
The KD tree data structure can be used for all kinds of searches that
involve N-dimensional vectors. For example, neighbor searches (find all points
within a radius of a given point) or finding all point pairs in a set
that are within a certain radius of each other. See "Computational Geometry: 
Algorithms and Applications" (Mark de Berg, Marc van Kreveld, Mark Overmars, 
Otfried Schwarzkopf).
"""


def getKDTree(coords, bucket_size=1):
    """Internal function to get KDTree for coordinates without any checks."""

    from KDTree import KDTree
    return KDTree(coords, bucket_size)
    

"""This module defines utility functions for drawing."""
# Note for developers: this module should be imported only a specific function is needed. 
# Please do not import this module when initializing ProDy for robustness reasons.

import numpy as np

from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

__all__ = ['Arrow3D', 'drawArrow3D']

class Arrow3D(FancyArrowPatch):
    """This function is implemented by tacaswell on stackoverflow: 
    https://stackoverflow.com/a/29188796."""
    def __init__(self, xs, ys, zs, *args, **kwargs):
        super(Arrow3D, self).__init__((0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

def drawArrow3D(*xyz, **kwargs):
    from matplotlib.pyplot import gca

    allargs = xyz
    xs = []; ys = []; zs = []; args = []

    for arg in allargs:
        if not np.isscalar(arg):
            if len(arg) == 3:
                xs.append(arg[0])
                ys.append(arg[1])
                zs.append(arg[2])
                continue
        args.append(arg)

    ax = gca()
    if 'arrowstyle' not in kwargs:
        kwargs['arrowstyle'] = '-|>'
    if 'mutation_scale' not in kwargs:
        kwargs['mutation_scale'] = 20
    if 'color' not in kwargs:
        kwargs['color'] = 'k'

    arrow = Arrow3D(xs, ys, zs, *args, **kwargs)
    ax.add_artist(arrow)
    return arrow

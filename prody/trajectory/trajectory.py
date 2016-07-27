# -*- coding: utf-8 -*-
"""This module defines a class for handling multiple trajectories."""

import os.path

import numpy as np

from .trajbase import TrajBase
from .frame import Frame

from prody.trajectory import openTrajFile

__all__ = ['Trajectory']

class Trajectory(TrajBase):

    """A class for handling trajectories in multiple files."""

    def __init__(self, name, **kwargs):
        """Trajectory can be instantiated with a *name* or a filename. When
        name is a valid path to a trajectory file it will be opened for
        reading."""

        TrajBase.__init__(self, name)
        self._trajectory = None
        self._trajectories = []
        self._filenames = set()
        self._n_files = 0
        self._cfi = 0 # current file index
        assert 'mode' not in kwargs, 'mode is an invalid keyword argument'
        self._kwargs = kwargs
        if os.path.isfile(name):
            self.addFile(name)

    def __repr__(self):

        if self._closed:
            return '<Trajectory: {0} (closed)>'.format(self._title)

        link = ''
        if self._ag is not None:
            link = 'linked to ' + str(self._ag) + '; '

        files = ''
        if self._n_files > 1:
            files = '{0} files; '.format(self._n_files)

        next = 'next {0} of {1} frames; '.format(self._nfi, self._n_csets)

        if self._indices is None:
            atoms = '{0} atoms'.format(self._n_atoms)
        else:
            atoms = 'selected {0} of {1} atoms'.format(
                                            self.numSelected(), self._n_atoms)

        return '<Trajectory: {0} ({1}{2}{3}{4})>'.format(
                                        self._title, link, files, next, atoms)

    def _nextFile(self):

        self._cfi += 1
        if self._cfi < self._n_files:
            self._trajectory = self._trajectories[self._cfi]
            if self._trajectory.nextIndex() > 0:
                self._trajectory.reset()

    def _gotoFile(self, i):

        if i < self._n_files:
            self._cfi = i
            self._trajectory = self._trajectories[i]
            if self._trajectory.nextIndex() > 0:
                self._trajectory.reset()

    def setAtoms(self, atoms):

        for traj in self._trajectories:
            traj.setAtoms(atoms)
        TrajBase.setAtoms(self, atoms)

    setAtoms.__doc__ = TrajBase.setAtoms.__doc__

    def link(self, *ag):

        if ag:
            TrajBase.link(self, *ag)
            for traj in self._trajectories:
                traj.link(*ag)
        else:
            return self._ag

    link.__doc__ = TrajBase.link.__doc__

    def addFile(self, filename, **kwargs):
        """Add a file to the trajectory instance. Currently only DCD files
        are supported."""

        if not isinstance(filename, str):
            raise ValueError('filename must be a string')
        if os.path.abspath(filename) in self._filenames:
            raise IOError('{0} is already added to the trajectory'
                          .format(filename))
        assert 'mode' not in kwargs, 'mode is an invalid keyword argument'
        traj = openTrajFile(filename, **(kwargs or self._kwargs))
        n_atoms = self._n_atoms
        if n_atoms != 0 and n_atoms != traj.numAtoms():
            raise IOError('{0} must have same number of atoms as '
                            'previously loaded files'.format(traj.getTitle()))

        if self._n_files == 0:
            self._trajectory = traj
            self._title = traj.getTitle()
        if n_atoms == 0:
            self._n_atoms = traj.numAtoms()
            self._coords = traj._coords
        self._trajectories.append(traj)
        self._n_csets += traj.numFrames()
        self._n_files += 1
        if self._ag is not None:
            traj.setAtoms(self._ag)

    def numFiles(self):
        """Returns number of open trajectory files."""

        return self._n_files

    def getFilenames(self, absolute=False):
        """Returns list of filenames opened for reading."""

        return [traj.getFilename(absolute) for traj in self._trajectories]

    def getFrame(self, index):

        if self._closed:
            raise ValueError('I/O operation on closed file')
        self.goto(index)
        return next(self)

    def getCoordsets(self, indices=None):

        if self._closed:
            raise ValueError('I/O operation on closed file')
        if indices is None:
            indices = np.arange(self._n_csets)
        elif isinstance(indices, int):
            indices = np.array([indices])
        elif isinstance(indices, slice):
            indices = np.arange(*indices.indices(self._n_csets))
            indices.sort()
        elif isinstance(indices, (list, np.ndarray)):
            indices = np.unique(indices)
        else:
            raise TypeError('indices must be an integer or a list of '
                            'integers')

        nfi = self._nfi
        self.reset()
        coords = np.zeros((len(indices), self.numSelected(), 3),
                          self._trajectories[0]._dtype)
        prev = -1
        next = self.nextCoordset
        for i, index in enumerate(indices):
            diff = index - prev
            if diff > 1:
                self.skip(diff - 1)
            coords[i] = next()
            prev = index
        self.goto(nfi)
        return coords

    getCoordsets.__doc__ = TrajBase.getCoordsets.__doc__

    def __next__(self):

        if self._closed:
            raise ValueError('I/O operation on closed file')
        nfi = self._nfi
        if nfi < self._n_csets:
            traj = self._trajectory
            while traj._nfi == traj._n_csets:
                self._nextFile()
                traj = self._trajectory
            unitcell = traj._nextUnitcell()
            coords = traj._nextCoordset()

            if self._ag is None:
                frame = Frame(self, nfi, coords, unitcell)
            else:
                frame = self._frame
                Frame.__init__(frame, self, nfi, None, unitcell)
                self._ag.setACSLabel(self._title + ' frame ' + str(self._nfi))
            self._nfi += 1
            return frame

    __next__.__doc__ = TrajBase.__next__.__doc__
    next = __next__

    def nextCoordset(self):

        if self._closed:
            raise ValueError('I/O operation on closed file')
        if self._nfi < self._n_csets:
            traj = self._trajectory
            while traj._nfi == traj._n_csets:
                self._nextFile()
                traj = self._trajectory
            if self._ag is not None:
                self._ag.setACSLabel(self._title + ' frame ' + str(self._nfi))
                traj = self._trajectory
            self._nfi += 1
            #if self._indices is None:
            return traj.nextCoordset()
            #else:
            #    return traj.nextCoordset()[self._indices]

    nextCoordset.__doc__ = TrajBase.nextCoordset.__doc__

    def goto(self, n):

        if self._closed:
            raise ValueError('I/O operation on closed file')
        if not isinstance(n, int):
            raise ValueError('n must be an integer')
        n_csets = self._n_csets
        if n == 0:
            self.reset()
        else:
            if n < 0:
                n = n_csets + n
            if n < 0:
                n = 0
            elif n > n_csets:
                n = n_csets
            nfi = n
            for which, traj in enumerate(self._trajectories):
                if traj._n_csets >= nfi:
                    break
                else:
                    nfi -= traj._n_csets
            self._gotoFile(which)
            self._trajectory.goto(nfi)
            self._nfi = n

    goto.__doc__ = TrajBase.goto.__doc__

    def skip(self, n):

        if self._closed:
            raise ValueError('I/O operation on closed file')
        if not isinstance(n, int):
            raise ValueError('n must be an integer')
        left = self._n_csets - self._nfi
        if n > left:
            n = left
        while self._nfi < self._n_csets and n > 0:
            traj = self._trajectory
            skip = min(n, traj.numFrames() - traj.nextIndex())
            traj.skip(skip)
            if n > skip:
                self._nextFile()
            self._nfi += skip
            n -= skip

    skip.__doc__ = TrajBase.skip.__doc__

    def reset(self):

        if self._closed:
            raise ValueError('I/O operation on closed file')
        if self._trajectories:
            for traj in self._trajectories:
                traj.reset()
            self._trajectory = self._trajectories[0]
            self._cfi = 0
            self._nfi = 0

    reset.__doc__ = TrajBase.reset.__doc__

    def close(self):

        for traj in self._trajectories:
            traj.close()
        self._closed = True

    close.__doc__ = TrajBase.close.__doc__

    def hasUnitcell(self):

        return np.all([traj.hasUnitcell() for traj in self._trajectories])

    hasUnitcell.__doc__ = TrajBase.hasUnitcell.__doc__

    def getTimestep(self):
        """Returns list of timestep sizes, one number from each file."""

        return [traj.getTimestep() for traj in self._trajectories]

    def getFirstTimestep(self):
        """Returns list of first timestep values, one number from each file."""

        return [traj.getFirstTimestep() for traj in self._trajectories]

    def getFrameFreq(self):
        """Returns list of timesteps between frames, one number from each file.
        """

        return [traj.getFrameFreq() for traj in self._trajectories]

    def numFixed(self):
        """Returns a list of fixed atom numbers, one from each file."""

        return [traj.numFixed() for traj in self._trajectories]

# -*- coding: utf-8 -*-
# ProDy: A Python Package for Protein Dynamics Analysis
# 
# Copyright (C) 2010-2012 Ahmet Bakan
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>

"""This module defines classes for handling trajectory files in DCD format."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import os
from time import time
from struct import calcsize, unpack, pack

import numpy as np

from prody.atomic import Atomic
from prody.ensemble import Ensemble
from prody.tools import getsize, checkCoords, now

from frame import Frame
from trajbase import TrajBase
from trajfile import TrajFile
from trajectory import Trajectory

__all__ = ['DCDFile', 'parseDCD', 'writeDCD']

pkg = __import__(__package__)
LOGGER = pkg.LOGGER

DEBUG = False

PISQUARE = np.pi ** 2
RECSCALE32BIT = 1
RECSCALE64BIT = 2

class DCDFile(TrajFile):
    
    """A class for reading DCD files. DCD header and first frame is parsed at 
    instantiation.  Coordinates from the first frame is set as the reference 
    coordinates."""
    
    def __init__(self, filename, mode='r'):
        
        TrajFile.__init__(self, filename, mode)
        if self._mode != 'w':
            self._parseHeader()
            
    __init__.__doc__ = TrajFile.__init__.__doc__
        
    def _parseHeader(self):
        """Read the header information from a dcd file.
        Input: fd - a file struct opened for binary reading.
        Output: 0 on success, negative error code on failure.
        Side effects: *natoms set to number of atoms per frame
                      *nsets set to number of frames in dcd file
                      *istart set to starting timestep of dcd file
                      *nsavc set to timesteps between dcd saves
                      *delta set to value of trajectory timestep
                      *nfixed set to number of fixed atoms 
                      *freeind may be set to heap-allocated space
                      *reverse set to one if reverse-endian, zero if not.
                      *charmm set to internal code for handling charmm data.
        """
        
        dcd = self._file
        endian = '' #'=' # native endian
        rec_scale = RECSCALE32BIT
        charmm = None
        dcdcordmagic = unpack(endian+'i', 'CORD')[0]
        # Check magic number in file header and determine byte order
        bits = dcd.read(calcsize('ii'))
        temp = unpack(endian+'ii', bits)
        if DEBUG: print len(temp), temp
        if temp[0] + temp[1] == 84:
            LOGGER.info('Detected CHARMM -i8 64-bit DCD file of native '
                        'endianness.')
            rec_scale = RECSCALE64BIT
        elif temp[0] == 84 and temp[1] == dcdcordmagic:
            LOGGER.info('Detected standard 32-bit DCD file of native '
                        'endianness.')
        else:
            if unpack('>ii', bits) == temp:
                endian = '>'
            else:
                endian = '<'
            temp = unpack(endian+'ii', bits)
            if temp[0] + temp[1] == 84:
                rec_scale = RECSCALE64BIT
                LOGGER.info('Detected CHARMM -i8 64-bit DCD file of opposite '
                            'endianness.')
            else:
                endian = ''
                temp = unpack(endian+'ii', bits)
                if temp[0] == 84 and temp[1] == dcdcordmagic:
                    LOGGER.info('Detected standard 32-bit DCD file of '
                                'opposite endianness.')
                else:
                    raise IOError('Unrecognized DCD header or unsupported '
                                  'DCD format.')
                    
        
        # check for magic string, in case of long record markers
        if rec_scale == RECSCALE64BIT:
            raise IOError('CHARMM 64-bit DCD files are not yet supported.');
            temp = unpack('I', dcd.read(calcsize('I')))
            if temp[0] != dcdcordmagic: 
                raise IOError('Failed to find CORD magic in CHARMM -i8 64-bit '
                              'DCD file.');
        
        # Buffer the entire header for random access
        bits = dcd.read(80)
        # CHARMm-genereate DCD files set the last integer in the
        # header, which is unused by X-PLOR, to its version number.
        # Checking if this is nonzero tells us this is a CHARMm file
        # and to look for other CHARMm flags.
        temp = unpack(endian + 'i'*20 , bits)
        if DEBUG: print len(temp), temp
        if temp[-1] != 0:
            charmm = True

        if charmm:
            LOGGER.info('CHARMM format DCD file (also NAMD 2.1 and later).')
            temp = unpack(endian + 'i'*9 + 'f' + 'i'*10 , bits)
        else:
            LOGGER.info('X-PLOR format DCD file (also NAMD 2.0 and earlier) '
                        'is not supported.')
            return None
        
        # Store the number of sets of coordinates (NSET)
        self._n_csets = temp[0]
        # Store ISTART, the starting timestep
        self._first_ts = temp[1]
        # Store NSAVC, the number of timesteps between dcd saves
        self._framefreq = temp[2]
        # Store NAMNF, the number of fixed atoms
        self._n_fixed = temp[8]
        
        if self._n_fixed > 0:
            raise IOError('DCD files with fixed atoms is not yet supported.')
        
        # Read in the timestep, DELTA
        # Note: DELTA is stored as double with X-PLOR but as float with CHARMm
        self._timestep = temp[9]
        self._unitcell = temp[10] == 1
        
        # Get the end size of the first block
        if unpack(endian+'i', dcd.read(rec_scale * calcsize('i')))[0] != 84:
            raise IOError('Unrecognized DCD format.')
        
        # Read in the size of the next block
        temp = unpack(endian+'i', dcd.read(rec_scale * calcsize('i')))
        if DEBUG: print len(temp), temp
        if temp[0] != 164:
            raise IOError('Unrecognized DCD format.')

        # Read NTITLE, the number of 80 character title strings there are
        temp = unpack(endian+'i', dcd.read(rec_scale * calcsize('i')))
        if DEBUG: print len(temp), temp
        self._dcdtitle = dcd.read(80)
        if DEBUG: print self._dcdtitle
        self._remarks = dcd.read(80)
        if DEBUG: print self._remarks
        # Get the ending size for this block
        temp = unpack(endian+'i', dcd.read(rec_scale * calcsize('i')))
        if DEBUG: print len(temp), temp
        if temp[0] != 164:
            raise IOError('Unrecognized DCD format.')


        # Read in an integer '4'
        if unpack(endian+'i', dcd.read(rec_scale * calcsize('i')))[0] != 4:
            raise IOError('Unrecognized DCD format.')

        # Read in the number of atoms
        self._n_atoms = unpack(endian+'i',dcd.read(rec_scale*calcsize('i')))[0]
        # Read in an integer '4'
        if unpack(endian+'i', dcd.read(rec_scale * calcsize('i')))[0] != 4:
            raise IOError('Bad DCD format.')

        self._is64bit = rec_scale == RECSCALE64BIT
        self._endian = endian
        self._n_floats = (self._n_atoms + 2) * 3
        
        if self._is64bit:
            if self._unitcell:
                self._bytes_per_frame = 56 + self._n_floats * 8
            else:
                self._bytes_per_frame = self._n_floats * 8
            LOGGER.warning('Reading of 64 bit DCD files has not been tested. '
                           'Please report any problems that you may find.')
            self._dtype = np.float64
        else: 
            if self._unitcell:
                self._bytes_per_frame = 56 + self._n_floats * 4
            else:
                self._bytes_per_frame = self._n_floats * 4
            self._dtype = np.float32
        
        self._first_byte = self._file.tell()
        n_csets = (getsize(self._filename) - self._first_byte
                                                    ) / self._bytes_per_frame
        if n_csets != self._n_csets: 
            LOGGER.warning('DCD header claims {0:d} frames, file size '
                           'indicates there are actually {1:d} frames.'
                           .format(self._n_csets, n_csets))
            self._n_csets = n_csets
        self._coords = self.nextCoordset()
        self._file.seek(self._first_byte)
        self._nfi = 0
   
    def hasUnitcell(self):
        
        return self._unitcell
    
    hasUnitcell.__doc__ = TrajBase.hasUnitcell.__doc__ 
   
    
    def getRemarks(self):
        """Return remarks parsed from DCD file."""
        
        return self._remarks
        
    def next(self):
        
        if self._closed: 
            raise ValueError('I/O operation on closed file')
        nfi = self._nfi
        if nfi < self._n_csets:
            unitcell = self._nextUnitcell()
            coords = self._nextCoordset()
            frame = Frame(self, nfi, coords, unitcell)
            return frame
    
    next.__doc__ = TrajBase.next.__doc__  
        
    def nextCoordset(self):
        """Return next coordinate set."""
        
        if self._closed: 
            raise ValueError('I/O operation on closed file')
        if self._nfi < self._n_csets:
            #Skip extended system coordinates (unit cell data)
            if self._unitcell:
                self._file.seek(56, 1)
            if self._sel is None:
                return self._nextCoordset()
            else:            
                return self._nextCoordset()[self._indices]

    def _nextCoordset(self):
    
        n_floats = self._n_floats
        n_atoms = self._n_atoms
        xyz = np.fromfile(self._file, dtype=self._dtype, count=n_floats)
        if len(xyz) != n_floats:
            return None
        xyz = xyz.reshape((3, n_atoms+2)).T[1:-1,:]
        xyz = xyz.reshape((n_atoms, 3))
        if self._ag is not None:
            self._ag._setCoords(xyz, self._title + ' frame ' + str(self._nfi))
        self._nfi += 1
        return xyz

    nextCoordset.__doc__ = TrajBase.nextCoordset.__doc__  

    def _nextUnitcell(self):
        
        if self._unitcell:
            self._file.read(4)
            unitcell = np.fromfile(self._file, dtype=np.float64, count=6)
            unitcell = unitcell[[0,2,5,1,3,4]]
            if np.all(abs(unitcell[3:]) <= 1):
                # This file was generated by CHARMM, or by NAMD > 2.5, with the angle */
                # cosines of the periodic cell angles written to the DCD file.        */ 
                # This formulation improves rounding behavior for orthogonal cells    */
                # so that the angles end up at precisely 90 degrees, unlike acos().   */
                unitcell[3:] = 90. - np.arcsin(unitcell[3:]) * 90 / PISQUARE  
            self._file.read(4)
            return unitcell

    def getCoordsets(self, indices=None):
        """Returns coordinate sets at given *indices*. *indices* may be an 
        integer, a list of integers or ``None``. ``None`` returns all 
        coordinate sets."""
                
        if self._closed: 
            raise ValueError('I/O operation on closed file')
        if self._indices is None and \
            (indices is None or indices == slice(None)):
            nfi = self._nfi
            self.reset()
            n_floats = self._n_floats + self._unitcell * 14
            n_atoms = self._n_atoms
            n_csets = self._n_csets
            data = np.fromfile(self._file, self._dtype, 
                               n_floats * n_csets)
            if len(data) > n_floats * n_csets:
                n_csets = len(data)/n_floats
                data = data[:n_csets]
                LOGGER.warning('DCD is corrupt, {0:d} out of {1:d} frames '
                               'were parsed.'.format(n_csets, self._n_csets))
            data = data.reshape((n_csets, n_floats))
            if self._unitcell:
                data = data[:, 14:]
            data = data.reshape((n_csets, 3, n_atoms+2))
            data = data[:, :, 1:-1]
            data = data.transpose(0, 2, 1)
            self.goto(nfi)
            return data
        else:            
            return TrajFile.getCoordsets(self, indices)
    
    getCoordsets.__doc__ = TrajBase.getCoordsets.__doc__

    def write(self, coords, unitcell=None, **kwargs):
        """Write *coords* to the file.  Number of atoms will be determined 
        based on the size of the first coordinate set.  If *unitcell* is 
        provided for the first coordinate set, it will be expected for the
        following coordinate sets as well.
        
        The following keywords are used when writing the first coordinate set 
        for files open at 'w' mode:        
            
        :arg timestep: timestep used for integration, default is 1
        :arg firsttimestep: number of the first timestep, default is 0
        :arg framefreq: number of timesteps between frames, default is 1"""
        
        if self._closed:
            raise ValueError('I/O operation on closed file')
        if self._mode == 'r':
            raise IOError('File not open for writing')
        coords = checkCoords(coords, 'coords', True, dtype=np.float32)
        if coords.ndim == 2:
            n_atoms = coords.shape[0]
            coords = [coords]
        else:
            n_atoms = coords.shape[1]
        if self._n_atoms == 0:
            self._n_atoms = n_atoms
        else:
            if self._n_atoms != n_atoms:
                raise ValueError('coords to not have correct number '
                                 'of atoms')
        dcd = self._file
        pack_i_4N = pack('i', self._n_atoms * 4)
        if self._n_csets == 0:
            if unitcell is None:
                self._unitcell = False
            else:
                self._unitcell = True
            timestep = float(kwargs.get('timestep', 1.0))
            first_ts = int(kwargs.get('firsttimestep', 0))
            framefreq = int(kwargs.get('framefreq', 1))
            n_fixed = 0

            pack_i_0 = pack('i', 0)
            pack_ix4_0x4 = pack('i'*4, 0, 0, 0, 0)
            pack_i_1 = pack('i', 1)
            pack_i_2 = pack('i', 2)
            pack_i_4 = pack('i', 4)
            pack_i_84 = pack('i', 84)
            pack_i_164 = pack('i', 164)

            dcd.write(pack_i_84)
            dcd.write('CORD')
            dcd.write(pack_i_0) # 0 Number of frames in file, none written yet
            dcd.write(pack('i', first_ts)) # 1 Starting timestep
            dcd.write(pack('i', framefreq)) # 2 Timesteps between frames
            dcd.write(pack_i_0) # 3 Number of timesteps in simulation
            dcd.write(pack_i_0) # 4 NAMD writes NSTEP or ISTART - NSAVC here?
            dcd.write(pack_ix4_0x4) # 5, 6, 7, 8
            dcd.write(pack('f', timestep)) # 9 timestep
            dcd.write(pack('i', int(self._unitcell))) # 10 with unitcell
            dcd.write(pack_ix4_0x4) # 11, 12, 13, 14
            dcd.write(pack_ix4_0x4) # 15, 16, 17, 18
            dcd.write(pack('i', 24)) # 19 Pretend to be CHARMM version 24
            dcd.write(pack_i_84)
            dcd.write(pack_i_164)
            dcd.write(pack_i_2)
            dcd.write('{0:80s}'.format('Created by ProDy'))
            dcd.write('{0:80s}'.format('REMARKS Created ' + 
                                       now().strftime('%d %B, %Y at %H:%M')))
            dcd.write(pack_i_164)
            
            dcd.write(pack_i_4)
            dcd.write(pack('i', n_atoms))
            dcd.write(pack_i_4)
            self._first_byte = dcd.tell()
        if self._unitcell: 
            if unitcell is None:
                raise TypeError('unitcell data is expected')
            else:
                uc = unitcell
                uc[3:] = np.sin((PISQUARE/90) * (90-uc[3:]))
                uc = uc[[0,3,1,4,5,2]]
                pack_i_48 = pack('i', 48)
        dcd.seek(0, 2)
        for xyz in coords:
            if self._unitcell:
                dcd.write(pack_i_48)
                uc.tofile(dcd)
                dcd.write(pack_i_48)
            xyz = xyz.T
            dcd.write(pack_i_4N)
            xyz[0].tofile(dcd)
            dcd.write(pack_i_4N)
            dcd.write(pack_i_4N)
            xyz[1].tofile(dcd)
            dcd.write(pack_i_4N)
            dcd.write(pack_i_4N)
            xyz[2].tofile(dcd)
            dcd.write(pack_i_4N)
            self._n_csets += 1
            dcd.seek(8, 0)
            dcd.write(pack('i', self._n_csets))
            dcd.seek(0, 2)
        self._nfi = self._n_csets

    def flush(self):
        """Flush the internal output buffer."""
        
        if self._mode != 'r':
            self._file.flush()
            os.fsync(self._file.fileno())
            
def parseDCD(filename, start=None, stop=None, step=None):
    """Parse CHARMM format DCD files (also NAMD 2.1 and later).  Returns an 
    :class:`Ensemble` instance. Conformations in the ensemble will be ordered 
    as they appear in the trajectory file.  Use :class:`DCDFile` class for 
    parsing  coordinates of a subset of atoms.
    
    :arg filename: DCD filename
    :type filename: str
    
    :arg start: index of first frame to read
    :type start: int
        
    :arg stop: index of the frame that stops reading
    :type stop: int
        
    :arg step: steps between reading frames, default is 1 meaning every frame
    :type step: int"""
    
    dcd = DCDFile(filename)
    time_ = time()
    n_frames = dcd.numFrames()
    LOGGER.info('DCD file contains {0:d} coordinate sets for {1:d} atoms.'
                .format(n_frames, dcd.numAtoms()))
    ensemble = dcd[slice(start,stop,step)]    
    dcd.close()
    time_ = time() - time_ or 0.01
    dcd_size = 1.0 * dcd.numFrames() * dcd._bytes_per_frame / (1024*1024)
    LOGGER.info('DCD file was parsed in {0:.2f} seconds.'.format(time_))
    LOGGER.info('{0:.2f} MB parsed at input rate {1:.2f} MB/s.'
                .format(dcd_size, dcd_size/time_))
    LOGGER.info('{0:d} coordinate sets parsed at input rate {1:d} frame/s.'
                .format(n_frames, int(n_frames/time_)))
    return ensemble



def writeDCD(filename, trajectory, start=None, stop=None, step=None, 
             align=False):
    """Write 32-bit CHARMM format DCD file (also NAMD 2.1 and later).
    *trajectory can be an :class:`Trajectory`, :class:`DCDFile`, or 
    :class:`Ensemble` instance. *filename* is returned upon successful
    output of file."""
    
    if not isinstance(trajectory, (TrajBase, Ensemble, Atomic)):
        raise TypeError('{0:s} is not a valid type for trajectory'
                        .format(type(trajectory)))
    
    irange = range(*slice(start, stop, 
                          step).indices(trajectory.numCoordsets()))
    n_csets = len(irange)
    if n_csets == 0:
        raise ValueError('trajectory does not have any coordinate sets, or '
                         'no coordinate sets are selected')
    
    if isinstance(trajectory, Atomic):
        isEnsemble = False
        isAtomic = True
        n_atoms = trajectory.numAtoms()
    else:
        isEnsemble = True
        isAtomic = False
        n_atoms = trajectory.numSelected()
    if n_atoms == 0:
        raise ValueError('no atoms are selected in the trajectory')
    if isinstance(trajectory, TrajBase):
        isTrajectory = True
        unitcell = trajectory.hasUnitcell()
        nfi = trajectory.getNextIndex() 
        trajectory.reset()
        pack_i_48 = pack('i', 48)
        if isinstance(trajectory, Trajectory):
            timestep = trajectory.getTimestep()[0]
            first_ts = trajectory.getFirstTimestep()[0]
            framefreq = trajectory.getFrameFreq()[0]
            n_fixed = trajectory.numFixed()[0]
        else:
            timestep = trajectory.getTimestep()
            first_ts = trajectory.getFirstTimestep()
            framefreq = trajectory.getFrameFreq()
            n_fixed = trajectory.numFixed()
    else:
        isTrajectory = False
        unitcell = False
        if isinstance(trajectory, Ensemble):
            frame = trajectory[0]
        else:
            frame = trajectory
            acsi = trajectory.getACSIndex()
        timestep = 1
        first_ts = 0
        framefreq = 1
        n_fixed = 0
        
    dcd = DCDFile(filename, mode='w')
    LOGGER.progress('Writing DCD', len(irange))
    prev = -1
    uc = None
    time_ = time()
    for j, i in enumerate(irange):
        diff = i - prev
        if diff > 1:
            trajectory.skip(diff-1)
        prev = i
        if isTrajectory:
            frame = trajectory.next()
            if frame is None:
                break
            if unitcell:
                uc = frame._getUnitcell()
                uc[3:] = np.sin((PISQUARE/90) * (90-uc[3:]))
                uc = uc[[0,3,1,4,5,2]]
        elif isEnsemble:
            frame._index = i
        else:
            frame.setACSIndex(i) 
        if align:
            frame.superpose()
        if j == 0:
            dcd.write(frame._getCoords(), uc, timestep=timestep, 
                      firsttimestep=first_ts, framefreq=framefreq)
        else:
            dcd.write(frame._getCoords(), uc)
        LOGGER.update(i)
    if isAtomic:
        trajectory.setACSIndex(acsi)
    j += 1
    LOGGER.clear()
    dcd.close()
    time_ = time() - time_ or 0.01
    dcd_size = 1.0 * (56 + (n_atoms * 3 + 6) * 4 ) * n_csets / (1024*1024)
    LOGGER.info('DCD file was written in {0:.2f} seconds.'.format(time_))
    LOGGER.info('{0:.2f} MB written at input rate {1:.2f} MB/s.'
                .format(dcd_size, dcd_size/time_))
    LOGGER.info('{0:d} coordinate sets written at output rate {1:d} frame/s.'
                .format(n_csets, int(n_csets/time_)))
    if j != n_csets:
        LOGGER.warning('Warning: {0:d} frames expected, {1:d} written.'
                       .format(n_csets, j))
    if isTrajectory:
        trajectory.goto(nfi)
    return filename

# Druggability: Python Package and VMD GUI for Druggability Index Analysis
# Copyright (C) 2010  Ahmet Bakan <ahb12@pitt.edu>
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
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""This module defines a class to analyze probe grid file.

This class is used by :class:`druggability.dia.DruggabilityIndexAnalysis`, 
which will frequently be referred to as DIA.

Classes:
    
* :class:`Probe`

Functions:

* :func:`get_expected_occupancy`

For each probe type, a probe card is defined as a dictionary. Probe cards 
contain information on the chemical identity and physical properties of the 
probe molecules:
    
* name: full chemical name
* radius: average distance of the central atom other (moleule) heavy atoms
* atomname: name of the central atom, used when writing PDB files
* n_atoms: number of heavy atoms
* charge: charge of the probe

For example, isopropanol card is defined as follows::
    
    IPRO = {'name': 'isopropanol', 
        'radius': 3.99,
        'atomname': 'C2',
        'n_atoms': 4,
        'charge': 0}

And then this card is registered as::    

    PROBE_CARDS['IPRO'] = IPRO
    
"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010  Ahmet Bakan'
__version__ = '0.5.2'


import time
import os.path
from os.path import isfile, join

import numpy as np

from druggability import functions
from druggability.grid import XPLOR, OpenDX
from druggability.exceptions import ProbeError

# Probe Cards are defined
IPRO = {'name': 'isopropanol', 
        'radius': 2.564, #V=70.602 3.24,
        'atomname': 'C2',
        'n_atoms': 4,
        'charge': 0}
IBUT = {'name': 'isobutane', 
        'radius': 2.664, #V=79.146 3.40,
        'atomname': 'C2',
        'n_atoms': 4,
        'charge': 0}
IPAM = {'name': 'isopropylamine', 
        'radius': 2.603, #V=73.873 3.18,
        'atomname': 'C2',
        'n_atoms': 4,
        'charge': +1}
ACET = {'name': 'acetate', 
        'radius': 2.376, #V=56.197 2.94,
        'atomname': 'C2',
        'n_atoms': 4,
        'charge': -1}
ACAM = {'name': 'acetamide',
        'radius': 2.421, #V=59.468 3.14,
        'atomname': 'C2',
        'n_atoms': 4,
        'charge': 0}
ALL = {'name': 'all-probes',
        'radius': 2.564*.6 + 2.664*.1 + 2.603*.1 + 2.376*.1 + 2.421*.1,
        'atomname': 'C2',
        'n_atoms': None,
        'charge': 0}

# JYL 20141009 different btw druggability and drugui
PROBE_CARDS = {'IPRO': IPRO, 'IBUT': IBUT, 'IPAM': IPAM, 'ACET': ACET, 
               'ACAM': ACAM, 'IBTN': IBUT, 'ACTT': ACET, 'ALL': ALL,
               'PRO2': IPRO} # included for backwards compatibility

import os
path = os.path.dirname(os.path.abspath(__file__))
    
with open(join(path, 'probe.dat')) as dat:
    for line in dat:
        items = line.split()
        if len(items) > 2 and items[1] == 'name' and items[0] not in PROBE_CARDS:
            PROBE_CARDS[items[0]] = {'name': ' '.join(items[2:]),
                                     'radius': 2.376,
                                     'atomname': 'C', 'charge': 0,
                                     'n_atoms': 4}
             
resi = None
with open(join(path, 'top_all36_cgenff.rtf')) as top: 
    for line in top:
        if line.startswith('RESI'):
            resi = None
            items = line.split()
            if items[1] in PROBE_CARDS:
                resi = items[1]
                PROBE_CARDS[resi]['charge'] = float(items[2])
                PROBE_CARDS[resi]['n_atoms'] = 0
        if resi and line.startswith('ATOM'):
            aname = line.split()[1]
            if not aname.startswith('H'):
                PROBE_CARDS[resi]['n_atoms'] += 1

def get_expected_occupancy(grid_spacing):
    """Return number of expected probes in a rectangular grid element.
    
    Calculation is based on a reference simulation of water and probe mixture
    (see :ref:`mixture`).
    
    """
    gs3 = grid_spacing * grid_spacing * grid_spacing
    n_probes = 343                # of isopropanols
    #nW = 6860                    # of waters
    ave_vol = 240934.57           # from a 10 ns reference simulation 
    #return (nW*gs3/ave_vol, nI*gs3/ave_vol)
    return n_probes * gs3 / ave_vol

__all__ = ['Probe']

class Probe(object):
    
    """A class to manipulate grid based probe count data.
    
    Probe count data typically comes from molecular dynamic simulations, and 
    may be present in different data file formats. This class recognizes grid 
    file formats for which a grid class is defined in the 
    :mod:`druggability.grid` module.

    .. attribute:: dia

       DIA instance that the probe grid data belongs to. 

    .. attribute:: type

       Type of the probe molecule, which is the residue name used in the force 
       field definition of the probe.

    .. attribute:: grid

       A :class:`druggability.grid.Grid` instance that holds probe grid data.

    .. attribute:: name 

       Full chemical name of the probe molecule.

    .. attribute:: radius

       An effective radius for the probe molecule.  

    .. attribute:: atomname

       Name of the atom that represents the probe. This is usually the name
       of the central carbon atom.

    .. attribute:: charge

       Charge of the probe molecule.
    
    """
    
    
    def __init__(self, dia, probe_type, grid):
        """Instantiate with a DIA instance, probe type and grid.
        
        At instantiation, passing a DIA instance is required. Many operations 
        over a probe object depends on the attributes of the DIA instance.

        :arg dia: DIA instance
        :type dia: :class:`druggability.dia.DruggabilityIndexAnalysis`
        
        :arg probe_type: probe type, see the following table
        :type probe_type: str
        
        :arg grid: grid filename or Grid instance
        :type grid: str or :class:`druggability.grid.Grid`
        
        ========== =============
        Probe Type Chemical name
        ========== =============
        IPRO, PRO2 isopropanol
        IBUT       isobutane
        IPAM       isopropylamine
        ACET       acetate
        ACAM       acetamide
        ========== =============
            
        Grid file types handled based on their extensions are:
            
        * .xplor (ptraj output grid file) 
        * .dx (VMD output grid file)
        
        """
        self.dia = dia
        self.type = probe_type
        try:
            card = PROBE_CARDS[self.type]
        except:
            self.dia.logger.info('An unknown probe is encountered: ' 
                                 + probe_type)
            card = {'name': probe_type.lower(),
                    'radius': 2.376,
                    'atomname': 'C',
                    'n_atoms': 4,
                    'charge': 0}
        self.name = card['name']
        self.radius = card['radius']
        self.atomname = card['atomname']
        self.charge = card['charge']
        
        if isinstance(grid, str):
            start = time.time()
            if grid.endswith('.xplor'):
                self.dia.logger.info('Parsing Xplor file {0:s}.'
                                     .format(os.path.relpath(grid)))
                self.grid = XPLOR(grid)
            elif grid.endswith('.dx'):
                self.dia.logger.info('Parsing OpenDX file {0:s}.'
                                     .format(os.path.relpath(grid)))
                self.grid = OpenDX(grid)
            else:
                raise ProbeError('Grid file {0:s} is not in a recognized '
                    'format (.xplor or .dx)'.format(grid))
            self.dia.logger.info('{0:s} was parsed in {1:.2f}s'
                                 .format(self.grid.name, time.time()-start))                
        else:
            self.grid = grid
        
        self.expected_population = get_expected_occupancy(self.grid.spacing[0]
                                    ) #* self.fraction  
        
    def __repr__(self):
        return 'Probe {0:s} ({1:s}) (grid: {2:s})'.format(self.name, self.type,
                                                          self.grid.filename)
        
    def write_grid(self, what=None, filename=None, smooth=True):
        """Write original, smoothed, free-energy, or enrichment grid.
        
        Output file content is determined by *what* argument. If *None*, is 
        provided grid original data is outputed. 
        
        By default, grid data is smoothed by calling the 
        :meth:`druggability.grid.Grid.smooth` method of grid instance.
        
        :arg what: original (None), enrichment, free-energy
        :type what: str or None
        
        :arg filename: output filename
        :type filename: str
         
        :arg smooth: smooth grid data before writing grid file
        :type smooth: bool, default is True
        
        """
        if what is not None and not what in ['enrichment', 'free-energy']:
            raise ProbeError('{0:s} is not recognized as a valid grid type'
                             .format(what))
        if smooth and self.grid.state == 'original':
            self.grid.smooth()
        if what is None:
            what = self.grid.state
        
        start = time.time()
        self.dia.logger.info('Writing {0:s} {1:s} grid in {2:s}.'
                         .format(self.grid.name, what, self.grid.format))

        if what in ['origional', 'smooth']:
            filename = self.grid.write(filename)
        else:
            array = self.grid.array
            state = self.grid.state
            self.grid.array = self.grid.array / \
                              self.dia.parameters['n_frames'] / \
                              self.expected_population
            if what == 'enrichment':
                self.grid.state = 'enrichment'
                filename = self.grid.write(filename)
            else:
                self.grid.state = 'free-energy'
                which = self.grid.array <= 1
                self.grid.array[which] = 0
                which = self.grid.array > 1
                self.grid.array[which] = -functions.convert_Kd_to_dG(
                                            self.grid.array[which])
                filename = self.grid.write(filename)
            self.grid.array = array
            self.grid.state = state
        self.dia.logger.info('{0:s} was written in {1:.2f}s'
                         .format(os.path.relpath(filename), time.time()-start))

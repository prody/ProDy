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

"""Class to analyze druggability simulation trajectories.

Classes:
    
* :class:`DruggabilityIndexAnalysis`

Algorithmic details and theoretical basis of DIA is explained in 
:ref:`methodology`. 


DIA Parameters
-------------------------------------------------------------------------------

A table of parameters required for druggability index analysis is provided 
below. User may set these parameters in a number of ways, including passing
them in a file, or passing them using a method in an interactive Python 
session. See the documentation on :class:`DruggabilityIndexAnalysis` 
for available options to set these parameters.

============ ======= ======== =================================================
Parameter    Default Unit     Description
============ ======= ======== =================================================
n_frames                      number of simulation frames used when
                              generating probe grid data 
                              
                              When VMD volmap plug-in is used, this parameter
                              needs to be set to 1, since volmap 
                              calculates a frame averaged value. When
                              AmberTools ptraj module is used, the
                              actual number of frames is required, 
                              since ptraj does not average over frames.

temperature          K        simulation temperature

delta_g      -1.239  kcal/mol binding free energy limit for hotspots

                              Hotspots with binding free energy below the value
                              set by :attr:`delta_g` will be used for assessing
                              druggability. A way to think of this parameter 
                              may be in terms of ligand efficiency.
                              Currently available probes have 4 heavy atoms:
                              :math:`-delta_g / 4 = 0.31 LE`.
                                
n_probes     6                number of probes to merge to make a drug-size 
                              molecule mimic
                              
                              Merging 6 hotspots puts together 24 heavy atoms.
                              Drug-like molecules have have 32-40 heavy atoms.
                              The difference is assumed to be atoms serving as
                              scaffolds without contributing to affinity.

min_n_probes 5                minimum number of probes in an acceptable 
                              solution

merge_radius 6.0     A        hotspot merge distance

                              This radius is used when clustering hotspots
                              and merging them to make up a drug-size molecule.

low_affinity 10      uM       lower binding affinity limit 

                              Potential binding sites with predicted affinity 
                              better than :attr:`low_affinity` will be 
                              reported.
                              
n_solutions  3                number of drug-size solutions to report for
                              each potential binding site
                              
                              This parameter may be set indefinetely large 
                              in order to see all solutions in a potential 
                              binding site.
                              
max_charge   2       e        maximum absolute charge on the drug-size molecule

n_charged    3                number of charged hotspots in a solutions 
============ ======= ======== =================================================

Note that some of the parameters do not have default values. User needs 
to set at least their (:attr:`n_frames` and :attr:`temperature`) values 
using :meth:`DruggabilityIndexAnalysis.set_parameters` method. 
Where a unit is specified, user provided parameteer value will be assumed 
to be in the specified unit. 

"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010  Ahmet Bakan'
__version__ = '0.5.2'


import os.path
import time
from glob import glob
from code import interact

import numpy as np

from druggability.abase import ABase
from druggability import functions
from druggability.probe import Probe, PROBE_CARDS
from druggability.exceptions import DIAError

__all__ = ['DruggabilityIndexAnalysis']

def _grow_molecules(alist, seed, dist, cutoff, size, dG):
    """Grow molecules starting from a seed hotspot.
    
    :arg alist: list of grown molecules
    :arg seed: list of indices
    :arg dist: distance matrix
    :arg cutoff: merge distance
    :arg size: maximum size
    
    """
    #mkp print alist, seed
    #mkp raw_input()
    if len(seed) >= size:
        alist.append(seed)
        return alist
    else:
        which = np.zeros(dist.shape[0])
        for i in seed:
            which += dist[i, :] <= cutoff
        which[seed] = False
        which = which.nonzero()[0]
        if len(which) == 0:
            return alist
        #i = which[(dG[which] == np.min(dG[which])).nonzero()[0]][0]
        #_grow_molecules(alist, seed + [i], dist, cutoff, size, dG)
        for i in which[:3]:
            _grow_molecules(alist, seed + [i], dist, cutoff, size, dG)
    return alist

class DruggabilityIndexAnalysis(ABase):
    
    """Analyze druggability of proteins based on MD simulation data.

    .. attribute:: parameters

       A dictionary that keeps parameter values.

    .. attribute:: hotspots

       A Numpy array instance that keeps hotspot data. Each row corresponds to 
       a hotspot. Columns correspond to:

       0. binding free energy (kcal/mol)
       1. X coordinate (A)
       2. Y coordinate (A)
       3. Z coordinate (A)
       4. probe identifier
       5. fractional occupancy
       6. charge (averaged by occupancy)
       7. effective radius
       8. cluster identifier

       First 8 columns are obtained by calling :meth:`identify_hotspots` 
       method. Last column is obtained by calling :meth:`assess_druggability` 
       method.
        

    """
    
    
    def __init__(self, name, **kwargs):
        """Initialize with a name and keyword arguments.
        
        All keyword arguments are passed to :class:`druggability.abase.ABase`.
        
        """
        ABase.__init__(self, name, **kwargs)
       
        self._probes = []
        self._dict = {}
        self._all = None
        self.logger.info('Druggability Analysis {0:s} is initialized.'
                         .format(name))
        
        self.parameters = {}
        
        self.parameters['temperature'] = None
        self.parameters['n_frames'] = None
        
        self.parameters['delta_g'] = -1.0
        
        self.parameters['n_probes'] = 6
        self.parameters['min_probes'] = 5
        self.parameters['merge_radius'] = 6.5
        self.parameters['n_solutions'] = 3
        self.parameters['low_affinity'] = -6.8589778920184257
        self.parameters['max_charge'] = 2.0
        self.parameters['n_charged'] = 3
        
        
        self.hotspots = None
            
    def __repr__(self):
        return 'Druggability analysis {0:s} containing {1:d} probes.'.format(
                self.name, len(self._probes))
    
    def set_parameters(self, **kwargs):
        """Set druggability analysis parameters.
        
        :mod:`druggability.dia` module documentation describes parameters 
        and their default values.       
            
        """
        if 'temperature' in kwargs:
            self.parameters['temperature'] = float(kwargs['temperature'])
            kwargs.pop('temperature')
            self.logger.info('Parameter: temperature {0:.2f} K'.format(
                                            self.parameters['temperature']))

        if 'n_frames' in kwargs:
            self.parameters['n_frames'] = int(kwargs['n_frames'])
            kwargs.pop('n_frames')
            self.logger.info('Parameter: n_frames {0:d}'.format(
                                                self.parameters['n_frames']))

        if 'delta_g' in kwargs: 
            self.parameters['delta_g'] = float(kwargs['delta_g'])
            kwargs.pop('delta_g')
            self.logger.info('Parameter: delta_g {0:.3f} kcal/mol'
                             .format(self.parameters['delta_g']))
        
        if 'n_probes' in kwargs:
            self.parameters['n_probes'] = int(kwargs['n_probes'])
            kwargs.pop('n_probes')
            self.logger.info('Parameter: n_probes {0:d}'.format(
                                                self.parameters['n_probes']))

        if 'min_n_probes' in kwargs:
            self.parameters['min_n_probes'] = int(kwargs['min_n_probes'])
            kwargs.pop('min_n_probes')
            self.logger.info('Parameter: min_n_probes {0:d}'.format(
                                            self.parameters['min_n_probes']))
        if 'n_solutions' in kwargs:
            self.parameters['n_solutions'] = int(kwargs['n_solutions'])
            kwargs.pop('n_solutions')
            self.logger.info('Parameter: n_solutions {0:d}'.format(
                                            self.parameters['n_solutions']))

        
        if 'merge_radius' in kwargs:
            self.parameters['merge_radius'] = float(kwargs['merge_radius'])
            kwargs.pop('merge_radius')
            self.logger.info('Parameter: merge_radius {0:.1f} A'.format(
                                            self.parameters['merge_radius']))
            
        if 'low_affinity' in kwargs:
            self.parameters['low_affinity'] = functions.convert_Kd_to_dG( 
                                float(kwargs['low_affinity']) * 1e-6)
            kwargs.pop('low_affinity')
            self.logger.info('Parameter: low_affinity {0[0]:.2f} '
            '{0[1]:s}'.format(functions.format_Kd(functions.convert_dG_to_Kd(
                                            self.parameters['low_affinity']))))
        if 'max_charge' in kwargs:
            self.parameters['max_charge'] = float(kwargs['max_charge'])
            kwargs.pop('max_charge')
            self.logger.info('Parameter: max_charge {0:.1f} e'.format(
                                                self.parameters['max_charge']))

        if 'n_charged' in kwargs:
            self.parameters['n_charged'] = int(kwargs['n_charged'])
            kwargs.pop('n_charged')
            self.logger.info('Parameter: n_charged {0:d}'.format(
                                                self.parameters['n_charged']))
            
                    
        while kwargs:
            key, value = kwargs.popitem()
            self.logger.warning('Parameter: {0:s} ({1:s}) is not a valid '
                                'parameter.'.format(key, str(value)))

    
    def parse_parameter_file(self, filename):
        """Parse parameters and probe grid types and filenames from a file.

        Example file content::
            
            n_frames 8000
            # Some comment
            temperature 300 K
            probe IPRO grid_IPRO.xplor 0.5
            probe ACAM grid_ACAM.xplor 0.5
            delta_g -1.0 kcal/mol
            n_probes 6
            # Another comment
        
        Note that all possible parameters are not included in this example.
        See documentation for "set_parameters" method.
        
        """
        kwargs = {}
        inp = open(filename)
        probes = []
        for line in inp:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            items = line.split()
            if items[0] == 'probe':
                probes.append((items[1], items[2]))
            else:
                if kwargs.has_key(items[0]):
                    raise DIAError('{0:s} parameter is repeated in {1:s}'
                                   .format(items[0], filename))
                kwargs[items[0]] = items[1]
        self.logger.info('Setting parameters from file {0:s}.'
                         .format(filename))
        self.set_parameters(**kwargs)
        self.logger.info('{0:d} probe files will be parsed.'
                         .format(len(probes)))
        for probe in probes:
            self.add_probe(probe[0], probe[1])

    
    def add_probe(self, probe_type, filename):
        """Add a probe grid file to the analysis.
        
        :arg probe_type: probe type is the residue name of the probe used by  
                         the force field
        :type probe_type: str
        :arg filename: Path to the file that contains grid data for the probe. 
                       Acceptable file formats are XPLOR and OpenDX.        
        :type filename: str
        
        """
        if probe_type not in PROBE_CARDS:
            raise DIAError('Probe type {0:s} is not recognized.'
                           .format(str(probe_type)))
        for probe in self._probes:
            if probe.grid.filename == filename:
                raise DIAError('Probe grid file {0:s} has already been added '
                               'to the analysis.'.format(filename))
        if probe_type in self._dict:
            raise DIAError('Probe type {0:s} has already been added '
                           'to the analysis.'.format(probe_type))
            
        
        self._probes.append(Probe(self, probe_type, filename))
        self._dict[probe_type] = self._probes[-1]
        self._dict[self._probes[-1].name] = self._probes[-1]

        
    def get_probe(self, probe=None):
        """Return probe instance for a given probe name or type.
        
        Probe instance for given name (e.g. *isopropanol*) or probe type (e.g.
        *IPRO*) will be return. If *None* is provided, probe instance for
        combined probe data will be returned if it is available.
        
        """
        
        if probe is None:
            if self._all is not None:
                return self._all
            else:
                DIAError('Combined probe data is not available yet.')
        else:
            if self._dict.has_key(probe):
                return self._dict[probe]   
            else:
                DIAError('{0:s} is not recognized as a probe name or type'
                         .format(str(probe)))
                
                
                                  
    def identify_hotspots(self, **kwargs):
        """Identify probe binding hotspots.
        
        Keyword arguments are passed to :meth:`set_parameters` method.
        
        A probe binding hotspot is a grid element with binding free energy 
        lower than the value of :attr:`delta_g` parameter and than the 
        binding free energies of the surrounding grid elements.

        .. seealso::

            For algorithmic details see :ref:`identifyhotspots`.
                    
        """
        if not self._probes:
            raise DIAError('You need to add probes before you '
                           'can identify binding hotspots.')
        self.set_parameters(**kwargs)
        n_frames = self.parameters['n_frames'] 
        if n_frames is None:
            raise DIAError('n_frames is a required parameter, and is not set'
                           'to identify binding hotspots.')
        delta_g = self.parameters['delta_g']
        if delta_g is None: 
            raise DIAError('delta_g parameter is not set, and is required '
                           'to identify binding hotspots.')
        temperature = self.parameters['temperature']
        # Calculate enrichment for given delta_g
        class ProbeError(Exception):
            pass
        if delta_g >= 0: 
            raise ProbeError('{0:.2f} kcal/mol is not a valid binding '
                                 'free energy limit'.format(delta_g))
        enrichment = functions.convert_dG_to_Kd(-delta_g, temperature)
        self.logger.info('Searching probe binding hotspots with deltaG less ' 
                         'than {0:.2f} kcal/mol (~{1:.0f} folds enrichment).'
                         .format(delta_g, enrichment))

        start = time.time()
        for probe in self._probes:
            if probe.grid.state == 'original':
                probe.grid.smooth()
            elif probe.grid.state != 'original': 
                raise DIAError('{0:s} grid status {1:s} is not recognized'
                               .format(probe.type, probe.grid.state))
        # Get all hotspots from the combined grid (get_hotspots smooths grid)
        if len(self._probes) == 1:
            self._all = self._probes[0]
        else:
            # Merge individual grids
            self._all = Probe(self, 'ALL', np.sum(
                              [probe.grid for probe in self._probes]))
            self._all.grid.filename = self.name + '_ALL' + os.path.splitext(
                                           self._probes[0].grid.filename)[1] 
        probes = self._all

        # Calculate population for given delta_g
        population = enrichment * n_frames * probes.expected_population 
        #print 'population', population
        #population = population / 5
        which = (probes.grid.array >= population).nonzero()
        #print 'n', len(which[0])
        if len(which[0]) == 0:
            raise ProbeError('No hotspots found with binding free energy '
                             'below {0:.2f} kcal/mol (~{1:.0f} folds '
                             'enrichment)'.format(delta_g, enrichment))
        spots = sorted(zip(-probes.grid.array[which], which[0], which[1], which[2]))
        #spots.sort()
        spots = np.array(spots)
        """
        # write all spots for MDM2
        spots[:, 0] = -spots[:, 0] / n_frames / probes.expected_population
        spots[:, 0] = -functions.convert_Kd_to_dG(spots[:, 0], temperature)
        self.hotspots = spots
        self._atomnames = np.array(['CA'] * len(spots))
        self._resnames = np.array(['IPRO'] * len(spots))
        self.hotspots[:, 1:4] = (self.hotspots[:, 1:4] * self._all.grid.spacing
                                 ) + self._all.grid.offset
        self.hotspots = np.concatenate([self.hotspots, 
                                        np.ones((len(spots), 1))], 1)
        self.hotspots = np.concatenate([self.hotspots, 
                                        np.ones((len(spots), 1))], 1)
        self.hotspots = np.concatenate([self.hotspots, 
                                        np.ones((len(spots), 1))], 1)
        self.hotspots = np.concatenate([self.hotspots, 
                                        np.ones((len(spots), 1))], 1)
        self.write_hotspots()
        self.hotspots = None
        stop
        """
        radii = np.array([probe.radius for probe in self._probes])
        charges = np.array([probe.charge for probe in self._probes])
        i = 0
        type_and_occupancy = []
        full_occupancy = []
        charge_list = []
        radii_list = []
        while i < len(spots):
            ijk = tuple(spots[i, 1:].astype('i'))
            counts = [(probe.grid.array[ijk], k) 
                      for k, probe in enumerate(self._probes)]
            weights = np.array(counts)[:, 0].T
            weightsum = weights.sum()
            counts.sort(reverse=True)
            counts = np.array(counts)
            counts[:, 0] = counts[:, 0] / weightsum

            full_occupancy.append(counts)
            type_and_occupancy.append((counts[0, 1], counts[0, 0]))
            
            # eliminate points within weighted radius average
            radius = (radii * weights).sum() / weightsum
            charge_list.append([(charges * weights).sum() / weightsum])
            radii_list.append([radius])
            # an average distance corresponding C-X bond is added
            # to the radius to account for some buffer space
            # when connecting hotspots
            #cutoff = ((radius + 1.4)  / probes.grid.spacing[0]) ** 2
            vdw_compensation = 0
            cutoff = ((radius + vdw_compensation) / probes.grid.spacing[0]) ** 2
            keep = ((spots[:, 1:]-spots[i, 1:]) ** 2).sum(1) > cutoff
            keep[:i+1] = True
            spots = spots[keep]
            i += 1
        spots[:, 0] = -spots[:, 0] / n_frames / probes.expected_population 
        spots[:, 0] = -functions.convert_Kd_to_dG(spots[:, 0], temperature)
        
        self.hotspots = spots
        
        
        self.logger.info('{0:d} {1:s} binding spots were identified in '
                         '{2:.2f}s.'.format(len(spots), probes.name, 
                         time.time()-start))
        self.logger.info('Minimum binding free energy is {0:.2f} '
                         'kcal/mol.'.format(spots[0,0]))

        # print details of occupancy in debugging mode
        for j, counts in enumerate(full_occupancy):         
            self.logger.debug('Hotspot {0:3d} {1:5.2f} kcal/mol {2:s}'.format(
                    j+1, self.hotspots[j, 0],
                    ''.join(
                    [('{0:5.1f}% {1:s} '.format(count[0]*100, 
                                            self._probes[int(count[1])].type)) 
                    for count in counts[(counts[:,0] > 0).nonzero()[0]]]
                    )))

             
        self.hotspots = np.concatenate([self.hotspots, 
                                        np.array(type_and_occupancy)], 1)
        self.hotspots = np.concatenate([self.hotspots, 
                                        np.array(charge_list)], 1)
        self.hotspots = np.concatenate([self.hotspots, 
                                        np.array(radii_list)], 1)
        # Convert indices to coordinates
        # Example
        #       nX = 10     first = -4       last = 5   spacing = 0.5
        #       0    1    2    3    4    5    6    7    8    9 Python_Index
        #   ----|----|----|----|----|----|----|----|----|----|
        #      -4   -3   -2   -1    0    1    2    3    4    5 Xplor_Index
        # -2.0 -1.5 -1.0 -0.5   0   0.5  1.0  1.5  2.0  2.5 X_coordinate
        # X_coordinate = (Python_Index + first) * spacing
        # Python_Index = (X_coordinate / spacing) - first -> is shifting
        self.hotspots[:, 1:4] = (self.hotspots[:, 1:4] * self._all.grid.spacing
                                 ) + self._all.grid.offset
                        
        # Write a short report for each probe type
        if len(self._probes) > 1:
            for i, probe in enumerate(self._probes):
                which = (self.hotspots[:, 4] == i).nonzero()[0]
                self.logger.info('{0:s}: {1:d} {2:s} binding hotspots '
                                 'were identified.'
                                 .format(probe.type, len(which), probe.name))
                if len(which) > 0:
                    self.logger.info('{0:s}: lowest binding free energy is '
                                     '{1:.2f} kcal/mol.'
                                     .format(probe.type, 
                                             self.hotspots[which[0],0]))
        self._resnames = np.array([self._probes[i].type 
                                     for i in self.hotspots[:, 4].astype('i')])
        self._atomnames = np.array([self._probes[i].atomname 
                                     for i in self.hotspots[:, 4].astype('i')])


    def assess_druggability(self, **kwargs):
        """Identify potential binding sites and estimate achievable affinity.

        Keyword arguments are passed to :meth:`set_parameters` method.
        
        """
        if self.hotspots is None or not isinstance(self.hotspots, np.ndarray):
            raise DIAError('Binding hotspots have not been identified.')
        self.set_parameters(**kwargs)
        
        merge_radius = self.parameters.get('merge_radius', None)
        if not isinstance(merge_radius, float):
            raise DIAError('merge_radius of type float must be provided')
        low_affinity = self.parameters.get('low_affinity', None)
        if not isinstance(low_affinity, float):
            raise DIAError('low_affinity of type float must be provided')
        n_probes = self.parameters.get('n_probes', None)
        if not isinstance(n_probes, int):
            raise DIAError('n_probes of type int must be provided')
        min_n_probes = self.parameters.get('min_n_probes', None)
        if not isinstance(min_n_probes, int):
            raise DIAError('min_n_probes of type int must be provided')
        if min_n_probes > n_probes:
            raise DIAError('min_n_probes cannot be greater than n_probes')
        n_solutions = self.parameters.get('n_solutions', None)
        if not isinstance(n_solutions, int):
            raise DIAError('n_solutions of type int must be provided')
        max_charge = self.parameters.get('max_charge', None)
        if not isinstance(max_charge, float):
            raise DIAError('max_charge of type float must be provided')
        n_charged = self.parameters.get('n_charged', None)
        if not isinstance(n_charged, int):
            raise DIAError('n_charged of type int must be provided')

        # Calculate distances        
        length = self.hotspots.shape[0]
        dist2 = np.zeros((length, length))
        cutoff2 = merge_radius**2
        coords = self.hotspots[:, 1:4]
        for i in range(length):
            dist2[i+1:, i] = ((coords[i+1:, :] - coords[i])**2).sum(1)
        
        # Perform single-linkage hierarchical clustering
        start = time.time()
        self.logger.info('Clustering probe binding hotspots.')
        clusters = functions.single_linkage(dist2, cutoff2)
        self.hotspots = np.concatenate([self.hotspots, 
                                        clusters.reshape((length, 1))], 1) 
        self.logger.info('Clustering completed in {0:.2f}ms.'
                         .format((time.time()-start)*1000))
                         

        # Determine potential sites
        sites = []
        for i in set(clusters):
            which = (self.hotspots[:, -1] == i).nonzero()[0]
            if len(which) >= min_n_probes:
                delta_g = self.hotspots[which, 0].sum()
                #if delta_g <= low_affinity + 5:
                if delta_g < low_affinity:
                    if len(which) <= n_probes:
                        if abs(self.hotspots[which, 6].sum()) > max_charge:
                            continue
                    sites.append((delta_g, len(which), i))
        sites.sort()
        self.logger.info('{0:d} potential sites are identified.'
                         .format(len(sites)))
        #mkp print self.hotspots[0]
        #mkp print sites
                                 
        self.logger.info('Calculating achievable affinity ranges.')
        sdict = {}
        for k, (delta_g, size, i) in enumerate(sites):
            if delta_g < low_affinity:
                which = (self.hotspots[:, -1] == i).nonzero()[0]
                hotspots = self.hotspots[which, :]
                mlist = []
                mdict = {}
                n_which = len(which)
                if n_which > n_probes:
                    #mkp print 'len(which) > n_probes'
                    dist2sub = dist2[which, :][:, which]
                    dist2sub = dist2sub + dist2sub.T
                    ###print len(which), 'probes'
                    tempc = 0
                    for i in range(n_which - n_probes):
                        # grow molecules for all hotspots in site except those n_probes with highest dG
                        #mkp print 'for i in range(n_which - n_probes):'
                        grownups = _grow_molecules([], [i], dist2sub, 
                                                   cutoff2, n_probes,
                                                   hotspots[:,0])
                        #mkp stop
                        #grownups = _grow_molecules([], [n_which-i-1], 
                        #                           dist2sub[:n_which-i, :][:, :n_which-i], 
                        #                           cutoff2, n_probes)
                        tempc += len(grownups)
                        #mkp print 'tempc += len(grownups)'
                        for molecule in grownups:
                            mtuple = list(molecule)
                            mtuple.sort()
                            mtuple = tuple(mtuple)
                            if mtuple not in mdict: 
                                mdict[mtuple] = molecule
                                mol_delta_g = hotspots[molecule, 0].sum()
                                mol_charge = hotspots[molecule, 6].sum()
                                mol_n_charged = np.round(np.abs(
                                                  hotspots[molecule, 6])).sum()
                                #if mol_charge != 0:
                                #interact(local=locals())
                                if mol_delta_g <= low_affinity and \
                                   abs(mol_charge) <= max_charge and \
                                   mol_n_charged <= n_charged:
                                    mlist.append((mol_delta_g, list(mtuple), 
                                                  mol_charge))
                    ###print tempc, 'grownups'
                    if not mlist:
                        continue
                    mlist.sort()
                    drugindices = []
                    for i in mlist:
                        drugindices.append(i[0])
                    drugindices = np.array(drugindices)
                    affinities = functions.convert_dG_to_Kd(drugindices)
                else:
                    mol_charge = hotspots[:, 6].sum()
                    if abs(mol_charge) <= max_charge:
                        mlist.append((delta_g, range(len(which)), mol_charge))
                    drugindices = np.array([delta_g])
                    affinities = functions.convert_dG_to_Kd(drugindices)
                sdict[k] = {}
                sdict[k]['hotspots'] = hotspots
                sdict[k]['which'] = which
                sdict[k]['size'] = size
                sdict[k]['min_delta_g'] = hotspots[0, 0]
                sdict[k]['mean_delta_g'] = delta_g / size      
                sdict[k]['mlist'] = mlist
                sdict[k]['drugindices'] = drugindices
                sdict[k]['affinities'] = affinities
                    
        slist = []
        for key, value in sdict.items():
            slist.append((value['drugindices'][0], key))
        slist.sort()


        for filename in glob(os.path.join(self.workdir, 
                                          self.name + '_site_*.pdb')):
            os.remove(filename)
        
        for k, (temp, key) in enumerate(slist):
            hotspots = sdict[key]['hotspots']
            which = sdict[key]['which']
            atomnames = self._atomnames[which]
            resnames = self._resnames[which]
            # write all hotspots in this site
            filename = os.path.join(self.workdir, self.name + 
                                    '_site_{0:d}.pdb'.format(k+1))
            out = open(filename, 'w')                                               
            pdb = functions.get_pdb_lines(hotspots[:, 1:4], 
                atomnames=atomnames, resnames=resnames, 
                occupancy=hotspots[:, 5], bfactor=hotspots[:, 0])
            out.write(pdb)
            out.close()

            
            drugindices = sdict[key]['drugindices']
            affinities = sdict[key]['affinities']
            mlist = sdict[key]['mlist']
            self.logger.info('Site {0:d}: {1:d} probe binding hotspots'
                     .format(k+1, sdict[key]['size']))
            self.logger.info('Site {0:d}: Lowest probe binding free energy'
                ' {1:.2f} kcal/mol'.format(k+1, sdict[key]['min_delta_g']))                         
            self.logger.info('Site {0:d}: Average probe binding free energy'
                '{1:.2f} kcal/mol'.format(k+1, sdict[key]['mean_delta_g']))

            if len(drugindices) > 1:
                self.logger.info('Site {0:d}: Total of {1:d} '
                    'solutions.'.format(k+1, len(drugindices)))

                self.logger.debug('\n' + 
                    functions.get_histr(-np.log10(affinities),
                                        label='-log10(affinity)',
                                        title=('Achievable '
                                        'affinities for site {0:d}')
                                        .format(k+1)))
            self.logger.info('Site {0:d}: Lowest drug-like binding '
                             'free energy {1:.2f} kcal/mol'
                             .format(k+1, drugindices[0]))
            self.logger.info('Site {0:d}: Highest drug-like affinity '
                             '{1[0]:.3f} {1[1]:s}'
                     .format(k+1, functions.format_Kd(affinities[0])))  


            for m in range(min(n_solutions, len(mlist))):
                which = mlist[m][1]
                mol_charge = mlist[m][2]
                filename = os.path.join(self.workdir, self.name + 
                        '_site_{0:d}_soln_{1:d}.pdb'.format(k+1, m+1))
                out = open(filename, 'w')                    
                pdb = functions.get_pdb_lines(hotspots[which, 1:4], 
                                        atomnames=atomnames[which], 
                                        resnames=resnames[which], 
                                        occupancy=hotspots[which, 5], 
                                        bfactor=hotspots[which, 0])
                out.write(pdb)
                out.close()
                self.logger.debug('Site {0:d}: Solution {2:d} binding '
                                 'free energy {1:.2f} kcal/mol'
                                 .format(k+1, drugindices[m], m+1))
                self.logger.debug('Site {0:d}: Solution {2:d} affinity'
                                 ' {1[0]:.3f} {1[1]:s}'.format(k+1, 
                                 functions.format_Kd(affinities[m]), m+1))                                         
                self.logger.debug('Site {0:d}: Solution {2:d} total '
                                 'charge {1:.2f} e'
                                 .format(k+1, mol_charge, m+1))
                self.logger.debug('Site {0:d}: Solution {2:d} number of '
                                 'hotspots {1:d}'
                                 .format(k+1, len(which), m+1))
                
                volume = 0
                volcor = 0
                for a, w in enumerate(which): 
                    h1 = hotspots[w]
                    R = h1[7]
                    volume += R ** 3
                    # apply pairwise overlap correction
                    for b in range(a+1, len(which)):
                        h2 = hotspots[which[b]]
                        r = h2[7]
                        d = ((h1[1:4] - h2[1:4]) ** 2).sum() ** 0.5
                        if (R + r) > d:
                            volcor += (np.pi * (d**2 + 2*d*r - 3*r**2 + 
                                                2*d*R - 3*R**2 + 6*r*R) * 
                                       (R + r - d)**2) / 12 / d
                volume = volume * 4 / 3 * np.pi - volcor
                self.logger.debug('Site {0:d}: Solution {2:d} approximate '
                                 'volume {1:.2f} A^3'
                                 .format(k+1, volume, m+1))
                


                    
    
    def write_hotspots(self, probe_type=None):
        """Write probe binding hotspots in PDB format.
        
        :arg probe_type: probe type to output, if None all probes
        :type probe_type: str or None, default is None
        
        If :meth:`assess_druggability` method is called, cluster ids of 
        hotspots will be printed as the residue number of the hotspot.
        
        """
        if self.hotspots is None:
            raise DIAError('Binding hotspots have not been identified.')
        which = None
        if probe_type is None:
            which = range(self.hotspots.shape[0])
            filename = os.path.join(self.workdir, 
                                    self.name + '_all_hotspots.pdb')
        else:
            for i, probe in enumerate(self._probes):
                if probe.type == probe_type:
                    which = self.hotspots[:, 4] == i
                    break
        if which is None:
            raise DIAError('{0:s} is not recognized as a valid probe type'
                           .format(probe_type))
            filename = os.path.join(self.workdir, self.name + '_' + 
                                                  probe_type + '_hotspots.pdb')
        out = open(filename, 'w')
        if probe_type is None:
            out.write('REMARK Hotspots are sorted according to their binding '
            'free energies in ascending order.\n'
            'REMARK Bindiing free energies are written in b-factor column.\n'
            'REMARK Each hotspot is represented by the probe type that is '
            'most frequently observed at its location.\n'
            'REMARK Residue names correspond to so selected probe type.\n'
            'REMARK Occupancy values are the fraction of frames that the '
            'was occupied by the most frequent probe type at that location.\n'
            'REMARK Cluster numbers are written to residue numbers column.\n'
            'REMARK Single-linkage clustering cutoff was {0:.1f}A.\n'
            .format(self.parameters['merge_radius']))
        if self.hotspots.shape[1] > 6:
            resids = self.hotspots[which, -1].flatten().astype('int64') + 1
        else:
            resids = None
        pdb = functions.get_pdb_lines(self.hotspots[:, 1:4], 
                                      atomnames=self._atomnames[which], 
                                      resnames=self._resnames[which], 
                                      resids=resids, 
                                      occupancy=self.hotspots[:, 5], 
                                      bfactor=self.hotspots[:, 0])
        out.write(pdb)
        out.close()
        self.logger.info('Hotspots are written into file {0:s}.'
                         .format(os.path.relpath(filename)))


    def perform_analysis(self):
        """Perform all analysis steps at once.
        
        All parameters without a default value must have been set by the user.
        
        This method runs the following steps:
            
        * :meth:`identify_hotspots`
        * :meth:`assess_druggability`
        * :meth:`write_hotspots`
          
        """
        self.identify_hotspots()
        self.assess_druggability()
        self.write_hotspots()


    def evaluate_ligand(self, filename, **kwargs):
        """Predict affinity of a site bound to a given ligand.
        
        This method selects hotspots within given radius of the ligand atoms
        and calculates affinity based on selected hotspots.

        """
        """
        :arg filename: ligand pdb filename, should contain only ligand atoms
        :type filename: str
        
        :arg radius: the distance (A) at which hotspots will be considered as 
          overlaping with the ligand 
        :type radius: float, default is 1.5
        
        """
        if not self._probes:
            raise DIAError('You need to add probes before you '
                           'can identify binding hotspots.')
        n_frames = self.parameters['n_frames'] 
        if n_frames is None:
            raise DIAError('n_frames is a required parameter, and is not set'
                           'to identify binding hotspots.')
        delta_g = kwargs.get('delta_g', -0.01)
        if delta_g is None: 
            raise DIAError('delta_g parameter is not set, and is required '
                           'to identify binding hotspots.')
        temperature = self.parameters['temperature']
        # Calculate enrichment for given delta_g
        if delta_g >= 0: 
            raise ProbeError('{0:.2f} kcal/mol is not a valid binding '
                                 'free energy limit'.format(delta_g))
        enrichment = functions.convert_dG_to_Kd(-delta_g, temperature)
        self.logger.info('Searching probe binding hotspots with deltaG less ' 
                         'than {0:.2f} kcal/mol (~{1:.0f} folds enrichment).'
                         .format(delta_g, enrichment))

        probes = self._all

        # Calculate population for given delta_g
        population = enrichment * n_frames * probes.expected_population 
        #which = (probes.grid.array > population).nonzero()
        

        self.logger.info('Evaluating binding site of ligand in {0:s}.'.format(
                                                    os.path.relpath(filename)))
        pdb_file = open(filename)
        coords = []
        for line in pdb_file:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                coords.append((float(line[30:38]), float(line[38:46]),
                                    float(line[46:54])))
        pdb_file.close()
        coords = np.array(coords, 'd')
        indices = np.array((coords - self._all.grid.offset) 
                            / self._all.grid.spacing, 'i')
        ranges = np.array((kwargs.get('radius', 1.5) / self._all.grid.spacing).round(), 'i')
        ilists = [[], [], []]
        for index in indices:
            for i in range(ranges[0]):
                for j in range(ranges[1]):
                    for k in range(ranges[0]):
                        ilists[0].append(index[0]-i)
                        ilists[0].append(index[0]+i)
                        ilists[1].append(index[1]-j)
                        ilists[1].append(index[1]+j)
                        ilists[2].append(index[2]-k)
                        ilists[2].append(index[2]+k)

        spots = zip(-self._all.grid.array[ilists], ilists[0], ilists[1], ilists[2])
        spots.sort()
        spots = np.array(spots)
        spots = spots[spots[:, 0] < -population, :]
        
        radii = np.array([probe.radius for probe in self._probes])
        charges = np.array([probe.charge for probe in self._probes])
        i = 0
        type_and_occupancy = []
        full_occupancy = []
        charge_list = []
        radii_list = []
        while i < len(spots):
            ijk = tuple(spots[i, 1:].astype('i'))
            counts = [(probe.grid.array[ijk], k) 
                      for k, probe in enumerate(self._probes)]
            weights = np.array(counts)[:, 0].T
            weightsum = weights.sum()
            counts.sort(reverse=True)
            counts = np.array(counts)
            counts[:, 0] = counts[:, 0] / weightsum

            full_occupancy.append(counts)
            type_and_occupancy.append((counts[0, 1], counts[0, 0]))
            
            # eliminate points within weighted radius average
            radius = (radii * weights).sum() / weightsum
            charge_list.append([(charges * weights).sum() / weightsum])
            radii_list.append([radius])
            # an average distance corresponding C-X bond is added
            # to the radius to account for some buffer space
            # when connecting hotspots
            #cutoff = ((radius + 1.4)  / probes.grid.spacing[0]) ** 2
            cutoff = (radius / probes.grid.spacing[0]) ** 2
            keep = ((spots[:, 1:]-spots[i, 1:]) ** 2).sum(1) > cutoff
            keep[:i+1] = True
            spots = spots[keep]
            i += 1
        spots[:, 0] = -spots[:, 0] / n_frames / probes.expected_population 
        spots[:, 0] = -functions.convert_Kd_to_dG(spots[:, 0], temperature)
        hotspots = spots

        # print details of occupancy in debugging mode
        for j, counts in enumerate(full_occupancy):         
            self.logger.debug('Hotspot {0:3d} {1:5.2f} kcal/mol {2:s}'.format(
                    j+1, hotspots[j, 0],
                    ''.join(
                    [('{0:5.1f}% {1:s} '.format(count[0]*100, 
                                            self._probes[int(count[1])].type)) 
                    for count in counts[(counts[:,0] > 0).nonzero()[0]]]
                    )))

             
        hotspots = np.concatenate([hotspots, np.array(type_and_occupancy)], 1)
        hotspots = np.concatenate([hotspots, np.array(charge_list)], 1)
        hotspots = np.concatenate([hotspots, np.array(radii_list)], 1)
        hotspots[:, 1:4] = (hotspots[:, 1:4] * self._all.grid.spacing
                                 ) + self._all.grid.offset
        n_hotspots = len(hotspots)
        self.logger.info('{0:d} hotspots were identified.'
                         .format(n_hotspots))
        resnames = [self._probes[i].type for i in hotspots[:, 4].astype('i')]
        atomnames = [self._probes[i].atomname 
                     for i in hotspots[:, 4].astype('i')]
        filename = os.path.join(self.workdir, self.name + '_' +
                                        os.path.split(filename)[-1])
        delta_g = hotspots[:, 0].sum()
        affinity = functions.format_Kd(functions.convert_dG_to_Kd(delta_g))                                           

        self.logger.info('Estimated binding free energy is {0:.2f} kcal/mol.'
                         .format(delta_g))
        self.logger.info('Estimated affinity is {0[0]:.3f} {0[1]:s}.'
                         .format(affinity))                                         

        volume = 0
        volcor = 0
        for a, h1 in enumerate(hotspots): 
            R = h1[7]
            volume += R ** 3
            # apply pairwise overlap correction
            for b in range(a+1, n_hotspots):
                h2 = hotspots[b]
                r = h2[7]
                d = ((h1[1:4] - h2[1:4]) ** 2).sum() ** 0.5
                if (R + r) > d:
                    volcor += (np.pi * (d**2 + 2*d*r - 3*r**2 + 
                                        2*d*R - 3*R**2 + 6*r*R) * 
                               (R + r - d)**2) / 12 / d
        volume = volume * 4 / 3 * np.pi - volcor
        self.logger.info('Solution approximate volume {0:.2f} A^3'
                          .format(volume))

        pdb_file = open(filename, 'w')
        pdb_file.write('REMARK Estimated binding free energy is {0:.2f} '
                       'kcal/mol\n'.format(delta_g))
        pdb_file.write('REMARK Estimated affinity is {0[0]:.3f} {0[1]:s}\n'
                         .format(affinity))           
        pdb = functions.get_pdb_lines(hotspots[:, 1:4], atomnames=atomnames, 
                                      resnames=resnames, 
                                      occupancy=hotspots[:, 5], 
                                      bfactor=hotspots[:, 0])
        pdb_file.write(pdb)
        pdb_file.close()
        self.logger.info('Selected hotspots are written into {0:s}.'.format(
                                                    os.path.relpath(filename)))

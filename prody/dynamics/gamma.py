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

"""This module defines customized gamma functions for elastic network model
analysis."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import numpy as np

from prody.atomic import Atomic 

__all__ = ['Gamma', 'GammaStructureBased', 'GammaVariableCutoff',]

class Gamma(object):
    
    """Base class for facilitating use of atom type, residue type, or residue
    property dependent force constants (γ).
    
    Derived classes:
        
    * :class:`GammaStructureBased`
    * :class:`GammaVariableCutoff`"""
    
    
    def __init__(self):
        pass
    
    def gamma(self, dist2, i, j):
        """Return force constant.
        
        For efficiency purposes square of the distance between interacting
        atom/residue (node) pairs is passed to this function. In addition, 
        node indices are passed."""
        
        pass
    
    
class GammaStructureBased(Gamma):
    
    """Facilitate setting the spring constant based on the secondary structure 
    and connectivity of the residues.
    
    A recent systematic study [LT10]_ of a large set of NMR-structures analyzed 
    using a method based on entropy maximization showed that taking into 
    consideration properties such as sequential separation between 
    contacting residues and the secondary structure types of the interacting 
    residues provides refinement in the ENM description of proteins.
    
    This class determines pairs of connected residues or pairs of proximal 
    residues in a helix or a sheet, and assigns them a larger user defined 
    spring constant value.
    
     DSSP single letter abbreviations are recognized: 
       * **H**: α-helix
       * **G**: 3-10-helix
       * **I**: π-helix
       * **E**: extended part of a sheet
    
    *helix*: 
        Applies to residue (or Cα atom) pairs that are in the same helical 
        segment, at most 7 Å apart, and separated by at most 
        3 (3-10-helix), 4 (α-helix), or 5 (π-helix) residues.
        
    *sheet*:  
        Applies to Cα atom pairs that are in different β-strands and at most 
        6 Å apart.
        
    *connected*:
        Applies to Cα atoms that are at most 4 Å apart.
        
    Note that this class does not take into account insertion codes.        
    
    **Example**:

    Let's parse coordinates and header data from a PDB file, and then
    assign secondary structure to the atoms. 
        
    >>> from prody import *
    >>> ubi, header = parsePDB('1aar', chain='A', subset='calpha', header=True)
    >>> assignSecstr(header, ubi)
    <AtomGroup: 1aar_A_ca (76 atoms)>

    In the above we parsed only the atoms needed for this calculation, i.e.
    Cα atoms from chain A. 
    
    We build the Hessian matrix using structure based force constants as 
    follows;
    
    >>> gamma = GammaStructureBased(ubi)
    >>> anm = ANM('')
    >>> anm.buildHessian(ubi, gamma=gamma)
    
    We can obtain the force constants assigned to residue pairs from the 
    Kirchhoff matrix as follows: 
    
    >>> k = anm.getKirchhoff()
    >>> k[0,1] # a pair of connected residues
    -10.0
    >>> k[0,16] # a pair of residues from a sheet
    -6.0"""
    
    def __init__(self, atoms, gamma=1.0, helix=6.0, sheet=6.0, connected=10.0):
        """Setup the parameters.
        
        :arg atoms: A set of atoms with chain identifiers, residue numbers,
            and secondary structure assignments are set.
        :type atoms: :class:`~prody.atomic.Atomic`

        :arg gamma: Force constant in arbitrary units. Default is 1.0.
        :type gamma: float
            
        :arg helix: Force constant factor for residues hydrogen bonded in 
            α-helices, 3,10-helices, and π-helices. Default is 6.0, i.e.
            ``6.0`*gamma``.
        :type helix: float

        :arg sheet: Force constant factor for residue pairs forming a hydrogen 
            bond in a β-sheet. Default is 6.0, i.e. ``6.0`*gamma``.
        :type sheet: float
            
        :arg connected: Force constant factor for residue pairs that are
            connected. Default is 10.0, i.e. ``10.0`*gamma``.
        :type connected: float"""
        
        if not isinstance(atoms, Atomic):
            raise TypeError('atoms must be an Atomic instance')
        n_atoms = atoms.numAtoms()
        sstr = atoms.getSecstrs()
        assert sstr is not None, 'secondary structure assignments must be set'
        chid = atoms.getChids()
        assert chid is not None, 'chain identifiers must be set'
        rnum = atoms.getResnums()
        assert rnum is not None, 'residue numbers must be set'
        gamma = float(gamma)
        assert gamma > 0, 'gamma must be greater than 0'
        helix = float(helix)
        assert helix > 0, 'helix must be greater than 0'
        sheet = float(sheet)
        assert sheet > 0, 'sheet must be greater than 0'
        connected = float(connected)
        assert connected > 0, 'connected must be greater than 0'
        
        ssid = np.zeros(n_atoms)
        for i in range(1, n_atoms):
            if (sstr[i-1] == sstr[i] and chid[i-1] == chid[i] and
               rnum[i]-rnum[i-1] == 1): 
                ssid[i] = ssid[i-1]
            else:
                ssid[i] = ssid[i-1] + 1
        self._sstr = sstr
        self._chid = chid
        self._rnum = rnum
        self._ssid = ssid
        self._gamma = gamma
        self._helix = gamma * helix
        self._sheet = gamma * sheet
        self._connected = gamma * connected
    
    def getSecstrs(self):
        """Return a copy of secondary structure assignments."""
        
        return self._sstr.copy()    
    
    def getChids(self):
        """Return a copy of chain identifiers."""
        
        return self._chid.socopypy()    

    def getResnums(self):
        """Return a copy of residue numbers."""
        
        return self._rnum.copy()    


    def gamma(self, dist2, i, j):
        """Return force constant."""
        
        if dist2 <= 16:
            return self._connected
        sstr = self._sstr
        ssid = self._ssid
        rnum = self._rnum
        if ssid[i] == ssid[j]:
            i_j = abs(rnum[j] - rnum[i])
            if ((i_j <= 4 and sstr[i] == 'H') or 
                (i_j <= 3 and sstr[i] == 'G') or 
                (i_j <= 5 and sstr[i] == 'I')) and dist2 <= 49: 
                return self._helix
        elif sstr[i] == sstr[j] == 'E' and dist2 <= 36:
            return self._sheet
        
        return self._gamma
    
    
class GammaVariableCutoff(Gamma):
    
    """Facilitate setting the cutoff distance based on user defined 
    atom/residue (node) radii.
    
    Half of the cutoff distance can be thought of as the radius of a node. 
    This class enables setting different radii for different node types.
    
    **Example**:
    
    Let's think of a protein-DNA complex for which we want to use different
    radius for different residue types. Let's say, for protein Cα atoms we
    want to set the radius to 7.5 Å, and for nucleic acid phosphate atoms to 
    10 Å. We use the HhaI-DNA complex structure :file:`1mht`.

    >>> hhai = parsePDB('1mht')
    >>> ca_p = hhai.select('(protein and name CA) or (nucleic and name P)')
    >>> print( ca_p.getNames() ) # doctest: +ELLIPSIS
    ['P' 'P' 'P' 'P' 'P' 'P' 'P' 'P' 'P' 'P' 'P' 'P' 'P' 'P' 'P' 'P' 'P' 'P'
     'P' 'P' 'P' 'P' 'CA' 'CA' 'CA' 'CA' 'CA' 'CA' 'CA' 'CA' 'CA' 'CA' 'CA'
     ...
     'CA' 'CA' 'CA' 'CA' 'CA' 'CA' 'CA' 'CA' 'CA' 'CA' 'CA' 'CA' 'CA' 'CA' 'CA'
     'CA']
         
    We set the radii of atoms: 
     
    >>> varcutoff = GammaVariableCutoff(ca_p.getNames(), gamma=1,  
    ... default_radius=7.5, debug=True, P=10)
    >>> print( varcutoff.getRadii() ) # doctest: +ELLIPSIS
    [ 10.   10.   10.   10.   10.   10.   10.   10.   10.   10.   10.   10.
      10.   10.   10.   10.   10.   10.   10.   10.   10.   10.    7.5   7.5
      ...
       7.5   7.5   7.5   7.5   7.5   7.5   7.5   7.5   7.5   7.5   7.5   7.5
       7.5]
    
    The above shows that for phosphate atoms radii is set to 10 Å, because
    we passed the ``P=10`` argument.  As for Cα atoms, the default 7.5 Å
    is set as the radius (``default_radius=7.5``).  Note we also passed
    ``debug=True`` argument for demonstration purposes. This argument 
    allows printing debugging information on the screen.
    
    We build :class:`ANM` Hessian matrix as follows:  
        
    >>> anm = ANM('HhaI-DNA')
    >>> anm.buildHessian(ca_p, gamma=varcutoff, cutoff=20) # doctest: +ELLIPSIS
    P_0 -- P_1 effective cutoff: 20.0 distance: 7.00984971308 gamma: 1.0
    P_0 -- P_2 effective cutoff: 20.0 distance: 13.0828573714 gamma: 1.0
    P_0 -- P_3 effective cutoff: 20.0 distance: 17.7058394322 gamma: 1.0
    P_0 -- P_16 effective cutoff: 20.0 distance: 16.8254374386 gamma: 1.0
    P_0 -- P_17 effective cutoff: 20.0 distance: 16.4588128977 gamma: 1.0
    P_0 -- P_18 effective cutoff: 20.0 distance: 17.3123524109 gamma: 1.0
    ...
    CA_345 -- CA_347 effective cutoff: 15.0 distance: 6.0211884209 gamma: 1.0
    CA_345 -- CA_348 effective cutoff: 15.0 distance: 9.71260639581 gamma: 1.0
    CA_346 -- CA_347 effective cutoff: 15.0 distance: 3.80736523071 gamma: 1.0
    CA_346 -- CA_348 effective cutoff: 15.0 distance: 6.73513808322 gamma: 1.0
    CA_347 -- CA_348 effective cutoff: 15.0 distance: 3.80721748788 gamma: 1.0
        
    Note that we passed ``cutoff=20.0`` to the :meth:`ANM.buildHessian` 
    method.  This is equal to the largest possible cutoff distance (between 
    two phosphate atoms) for this system, and ensures that all of the 
    potential interactions are evaluated. 
    
    For pairs of atoms for which the actual distance is larger than the 
    effective cutoff, the :meth:`GammaVariableCutoff.gamma` method returns 
    ``0``.  This annuls the interaction between those atom pairs."""
    
    def __init__(self, identifiers, gamma=1., default_radius=7.5, debug=False, 
                 **kwargs):
        """Set the radii of atoms.
        
        :arg identifiers: List of atom names or types, or residue names.
        :type identifiers: list or :class:`numpy.ndarray`
        
        :arg gamma: Uniform force constant value. Default is 1.0.
        :type gamma: float
        
        :arg default_radius: Default radius for atoms whose radii is not set
            as a keyword argument. Default is 7.5
        :value default_radius: float
        
        :arg debug: Print debugging information. Default is ``False``.
        :type debug: bool
        
        Keywords in keyword arguments must match those in *atom_identifiers*.
        Values of keyword arguments must be :class:`float`."""
        
        self._identifiers = identifiers
        radii = np.ones(len(identifiers)) * default_radius
        
        for i, identifier in enumerate(identifiers): 
            radii[i] = kwargs.get(identifier, default_radius)
        self._radii = radii
        self._gamma = float(gamma)
        self._debug = bool(debug)
    
    def getRadii(self):
        """Return a copy of radii array."""
        
        return self._radii.copy()

    def getGamma(self):
        """Return the uniform force constant value."""
        
        return self._gamma

    def gamma(self, dist2, i, j):
        """Return force constant."""
        
        cutoff = (self._radii[i] + self._radii[j])
        cutoff2 = cutoff ** 2
        
        if dist2 < cutoff2: 
            gamma = self._gamma
        else:
            gamma = 0
        if self._debug:
            print self._identifiers[i]+'_'+str(i), '--', \
                  self._identifiers[j]+'_'+str(j), \
                  'effective cutoff:', cutoff, 'distance:', dist2**0.5, \
                  'gamma:', gamma    
        return gamma

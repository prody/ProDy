# -*- coding: utf-8 -*-
"""This module defines functions for parsing and writing `PDB files`_.

.. _EMD maps: http://emdatabank.org/mapformat.html"""

from collections import defaultdict
import os.path

from prody.atomic import AtomGroup
from prody.atomic import flags
from prody.utilities import openFile
from prody import LOGGER, SETTINGS

import struct as st
import numpy as np

__all__ = ['parseEMDStream', 'parseEMD']

class EMDParseError(Exception):
	pass

""" For  documentation"""

def parseEMD(emd, **kwargs):
	"""Returns an :class:`.AtomGroup` containing the information parsed from EMD file. 

	This function extends :func:`.parseEMDStream`.

	See :ref:`parseEMD` for a detailed usage example. 

	:arg emd: an EMD identifier or a file name, EMD files should be locally available. 
	"""

	title = kwargs.get('title', None)
	cutoff = float(kwargs.get('cutoff', 1.20))
	if not os.path.isfile(emd):
		raise IOError('EMD file {0} is not available in the directory {1}'
						.format(emd),os.getcwd())
	if title is None:
		title, ext = os.path.splitext(os.path.split(emd)[1])
		kwargs['title'] = title
	emd = openFile(emd, 'rt')
	result = parseEMDStream(emd, **kwargs)
	emd.close()
	return result

def _parseEMDLines(atomgroup, stream, format='EMD'):
	""" Returns an AtomGroup. see also :func:`.parseEMDStream()`.

	:arg stream: stream from parser.
	"""

	format = format.upper()
	if format == 'EMD':
		isEMD = True
	else:
		isEMD = False

	# Number of columns, rows, and sections (3 words, 12 bytes, 1-12)
	NC = st.unpack('<L', stream.read(4))[0]
	NR = st.unpack('<L', stream.read(4))[0]
	NS = st.unpack('<L', stream.read(4))[0]
	Ntot = NC * NR * NS

	# Mode (1 word, 4 bytes, 13-16)
	mode = st.unpack('<L', stream.read(4))[0]

	# Number of first column, row, section (3 words, 12 bytes, 17-28)
	ncstart = st.unpack('<L', stream.read(4))[0]
	nrstart = st.unpack('<L', stream.read(4))[0]
	nsstart = st.unpack('<L', stream.read(4))[0]

	# Number of intervals along x, y, z (3 words, 12 bytes, 29-40)
	Nx = st.unpack('<L', stream.read(4))[0]
	Ny = st.unpack('<L', stream.read(4))[0]
	Nz = st.unpack('<L', stream.read(4))[0]

	# Cell dimensions (Angstroms) (3 words, 12 bytes, 41-52)
	Lx = st.unpack('<f', f.read(4))[0]
    Ly = st.unpack('<f', f.read(4))[0]
    Lz = st.unpack('<f', f.read(4))[0]
	
	# Cell angles (Degrees) (3 words, 12 bytes, 53-64)
    a = st.unpack('<f', f.read(4))[0]
    b = st.unpack('<f', f.read(4))[0]
    c = st.unpack('<f', f.read(4))[0]

    # Which axis corresponds to column, row, and sections (1, 2, 3 for x, y ,z)
    # (3 words, 12 bytes, 65-76)
    mapc = st.unpack('<L', f.read(4))[0]
    mapr = st.unpack('<L', f.read(4))[0]
    maps = st.unpack('<L', f.read(4))[0]

    # Density values (min, max, mean) (3 words, 12 bytes, 77-88)
    dmin = st.unpack('<f', f.read(4))[0]
    dmax = st.unpack('<f', f.read(4))[0]
    dmean = st.unpack('<f', f.read(4))[0]

    # Not interested (1 word, 4 bytes, 89-92)
    stream.read(1*4)

	# Not interested (1 word, 4 bytes, 93-96)
    nsym = st.unpack('<f', f.read(4))[0]

    # Not interested (25 word, 4 bytes, 97-196)
    stream.read(25*4)




def parseEMDStream(stream, **kwargs):
	""" Returns an :class:`.AtomGroup` containing EMD data parsed from a stream of EMD file.

	:arg stream: Anything that implements the method ``readlines``
        (e.g. :class:`file`, buffer, stdin)"""

    ag = None
    if 'ag' in kwargs:
    	ag = kwargs['ag']
    	if not isinstance(ag, AtomGroup):
            raise TypeError('ag must be an AtomGroup instance')
        n_csets = ag.numCoordsets()

    biomol = kwargs.get('biomol', False)
	hd = None
	LOGGER.timeit()
	try:



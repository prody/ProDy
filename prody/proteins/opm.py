# -*- coding: utf-8 -*-
"""This module defines functions for accessing wwPDB servers."""

from os import getcwd
from os.path import join

from prody import LOGGER
from prody.utilities import openURL, sympath

from .localpdb import pathPDBFolder

__all__ = ['fetchPDBfromOPM']

def fetchPDBfromOPM(pdb, filename=None):
    url = 'https://opm-assets.storage.googleapis.com/pdb/%s.pdb' % pdb.lower()

    try:
        handle = openURL(url)
    except Exception as err:
        LOGGER.warn('{0} download failed ({1}).'.format(pdb, str(err)))
        filename = ''
    else:
        data = handle.read()
        if len(data):

            if filename is None:
                output_folder = pathPDBFolder()
                if output_folder is None:
                    output_folder = getcwd()

                filename = join(output_folder, pdb + '-opm.pdb')

            with open(filename, 'w+b') as pdbfile:
                pdbfile.write(data)

            LOGGER.debug('{0} downloaded ({1})'
                            .format(pdb, sympath(filename)))
        else:
            LOGGER.warn('{0} download failed, reason unknown.'
                        .format(pdb))
            filename = ''

    return filename

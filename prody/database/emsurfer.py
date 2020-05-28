# -*- coding: utf-8 -*-
"""This module defines functions for Emsurfer searching Protein Data Bank."""

import re
import numpy as np
from prody.proteins.emdfile import EMDMAP
from prody.utilities import createStringIO
from prody import LOGGER, PY3K
from prody import parseEMD, writeEMD
# if PY3K:
    # import urllib.parse as urllib
    # import urllib.request as urllib2
# else:
    # import urllib
    # import urllib2
import os

__all__ = ['EmsurferRecord', 'searchEmsurfer']

def searchEmsurfer(emd, **kwargs):
    """Search with the EM-Surfer server with input of EMD ID (or local EMD file).
    EM-Surfer server: http://kiharalab.org/em-surfer/
    
    :arg emd: EMD code or local EMD map file for the query protein
    """
    
    import requests
    from requests.models import Request
    
    LOGGER.timeit('_emsurfer')
    # timeout = 120
    timeout = kwargs.pop('timeout', 120)
    
    emsurferURL = "http://kiharalab.org/em-surfer/cgi-bin/listResults.cgi"
    
    volumeFilter = kwargs.get('volumeFilter', 'on')
    representation = kwargs.get('representation','recommend')
    minResolution = kwargs.get('minResolution', 0.5)
    maxResolution = kwargs.get('maxResolution', 30.)

    if isinstance(emd, EMDMAP):
        emdmap = emd
        stream = createStringIO()
        writeEMD(stream, emdmap)
        data = stream.getvalue()
        stream.close()
        files = {"file1" : data}

        emdId = emdmap.getTitle()
        emdId = ''
        emsurfer_title = 'Title_'+emdId
    elif isinstance(emd, str):
        if os.path.isfile(emd):
            emdmap = parseEMD(emd)
            filename = os.path.basename(emd)
            filename, ext = os.path.splitext(filename)
            if ext.lower() == '.gz':
                filename2, ext2 = os.path.splitext(filename)
                if ext2.lower() == '.emd':
                    filename = filename2
            emdId = filename
            files = {"file1" : data}
            emdId = ''
            emsurfer_title = 'Title_' + emdId
        else:
            emdId = emd
            emsurfer_title = 'Title_' + emdId
            files = ''

    method='post'
    url=emsurferURL
    params = { 'emdbid' : emdId, 'volumefilter' : volumeFilter, 'representation' : representation, 
               'minresolution': minResolution, 'maxresolution': maxResolution }

    # Generate request url deep inside 
    data=None; headers=None; cookies=None; files=None
    auth=None; timeout=None; allow_redirects=True; proxies=None
    hooks=None; stream=None; verify=None; cert=None; json=None
    req = Request(
        method=method.upper(),
        url=url,
        headers=headers,
        files=files,
        data=data or {},
        json=json,
        params=params or {},
        auth=auth,
        cookies=cookies,
        hooks=hooks,
    )
    session = requests.sessions.Session()
    prep = session.prepare_request(req)
    resp = session.send(prep)
    url = resp.url

    LOGGER.debug('Submitted Emsurfer search for EMD "{0}".'.format(emdId))
    LOGGER.info(url)
    LOGGER.clear()
    obj = EmsurferRecord(url, emdId, timeout=timeout, **kwargs)
        
    return obj


class EmsurferRecord(object):

    """A class to store results from Emsurfer EMD search."""

    def __init__(self, url, emdId, localFile=False, **kwargs):
        """Instantiate a emsurferEMD object instance.

        :arg url: url of Emsurfer results page or local emsurfer results file
        :arg emdId: EMD code for searched protein
        :arg localFile: provided url is a path for local emsurfer results file
        """

        self._url = url
        self._emdId = emdId
        timeout = kwargs.pop('timeout', 120)

        self._title = emdId
        self.isSuccess = self.fetch(self._url, localFile=localFile, timeout=timeout, **kwargs)

    def fetch(self, url=None, localFile=False, **kwargs):
        if localFile:
            emsurfer_file = open(url, 'r')
            data = emsurfer_file.read()
            emsurfer_file.close()
        else:
            import requests
            
            if url == None:
                url = self._url

            html = requests.get(url).content

            if PY3K:
                html = html.decode()

            LOGGER.clear()
            LOGGER.report('Emsurfer results were fetched in %.1fs.', '_emsurfer')
            data = html.strip().split('\n')
        
        data_list = []
        for line in data[3:-2]:
            data_list.append(tuple(line.split('\t')))

        # Rank	EMDB_ID	EUC_D	RESOLUTION
        emsurferInfo = np.array(data_list, dtype=[('Rank', '<i4'), ('EMDB_ID', '<U70'),
                                                  ('EUC_D', '<f4'), ('RESOLUTION', '<f4')])
        emdListAll = []
        self._emsurferInfo = emsurferInfo
        emsurfer_temp_dict = dict()
        for temp in self._emsurferInfo:
            temp_dict = dict()
            temp_dict['Rank'] = temp[0]
            temp_dict['EMDB_ID'] = emdbId = temp[1]
            temp_dict['EUC_D'] = temp[2]
            temp_dict['RESOLUTION'] = temp[3]
            emsurfer_temp_dict[emdbId] = temp_dict
            emdListAll.append(emdbId)
        self._emdListAll = tuple(emdListAll)
        self._emdList = self._emdListAll
        self._alignEMD = emsurfer_temp_dict
        LOGGER.info('Obtained ' + str(len(emdListAll)) + ' EMD matches from Emsurfer for '+self._emdId+'.')
        return True
        
    def getEMDs(self, filtered=True):
        """Returns EMD list (filters may be applied)"""
        if filtered:
            return self._emdList
        return self._emdListAll
        
    def getHits(self):
        return self._alignEMD
        
    def getFilterList(self):
        filterDict = self._filterDict
        temp_str = ', '.join([str(len(filterDict['len'])), str(len(filterDict['rmsd'])), str(len(filterDict['Z'])), str(len(filterDict['identity']))])
        LOGGER.info('Filter out [' + temp_str + '] for [length, RMSD, Z, identity]')
        return self._filterList
    
    def getMapping(self, key):
        try:
            info = self._alignEMD[key]
            mapping = [info['map_ref'], info['map_sel']]
        except KeyError:
            return None
        return mapping

    def getMappings(self):
        map_dict = {}
        for key in self._alignEMD:
            info = self._alignEMD[key]
            mapping = [info['map_ref'], info['map_sel']]
            map_dict[key] = mapping
        return map_dict

    mappings = property(getMappings)

    def filter(self, cutoff_len=None, cutoff_rmsd=None, cutoff_Z=None, cutoff_identity=None):
        """Filters out EMDs from the EMDList and returns the EMD list.
        EMDs satisfy any of following criterion will be filtered out.
        (1) Length of aligned residues < cutoff_len (must be an integer or a float between 0 and 1);
        (2) RMSD < cutoff_rmsd (must be a positive number);
        (3) Z score < cutoff_Z (must be a positive number);
        (4) Identity > cutoff_identity (must be an integer or a float between 0 and 1).
        """
        if cutoff_len == None:
            # cutoff_len = int(0.8*self._max_index)
            cutoff_len = 0
        elif not isinstance(cutoff_len, (float, int)):
            raise TypeError('cutoff_len must be a float or an integer')
        elif cutoff_len <= 1 and cutoff_len > 0:
            cutoff_len = int(cutoff_len*self._max_index)
        elif cutoff_len <= self._max_index and cutoff_len > 0:
            cutoff_len = int(cutoff_len)
        else:
            raise ValueError('cutoff_len must be a float between 0 and 1, or an int not greater than the max length')
            
        if cutoff_rmsd == None:
            cutoff_rmsd = 0
        elif not isinstance(cutoff_rmsd, (float, int)):
            raise TypeError('cutoff_rmsd must be a float or an integer')
        elif cutoff_rmsd >= 0:
            cutoff_rmsd = float(cutoff_rmsd)
        else:
            raise ValueError('cutoff_rmsd must be a number not less than 0')
            
        if cutoff_Z == None:
            cutoff_Z = 0
        elif not isinstance(cutoff_Z, (float, int)):
            raise TypeError('cutoff_Z must be a float or an integer')
        elif cutoff_Z >= 0:
            cutoff_Z = float(cutoff_Z)
        else:
            raise ValueError('cutoff_Z must be a number not less than 0')
            
        if cutoff_identity == None or cutoff_identity == 0:
            cutoff_identity = 100
        elif not isinstance(cutoff_identity, (float, int)):
            raise TypeError('cutoff_identity must be a float or an integer')
        elif cutoff_identity <= 1 and cutoff_identity > 0:
            cutoff_identity = float(cutoff_identity*100)
        elif cutoff_identity <= 100 and cutoff_identity > 0:
            cutoff_identity = float(cutoff_identity)
        else:
            raise ValueError('cutoff_identity must be a float between 0 and 1, or a number between 0 and 100')
            
        # debug:
        # print('cutoff_len: ' + str(cutoff_len) + ', ' + 'cutoff_rmsd: ' + str(cutoff_rmsd) + ', ' + 'cutoff_Z: ' + str(cutoff_Z) + ', ' + 'cutoff_identity: ' + str(cutoff_identity))
        
        emsurferInfo = self._alignEMD
        emdListAll = self._emdListAll
        missing_ind_dict = dict()
        ref_indices_set = set(range(self._max_index))
        filterListLen = []
        filterListRMSD = []
        filterListZ = []
        filterListIdentiry = []
        
        # keep the first EMD (query EMD)
        for emdId in emdListAll[1:]:
            temp_dict = emsurferInfo[emdId]
            # filter: len_align, identity, rmsd, Z
            if temp_dict['len_align'] < cutoff_len:
                # print('Filter out ' + emdId + ', len_align: ' + str(temp_dict['len_align']))
                filterListLen.append(emdId)
                continue
            if temp_dict['rmsd'] < cutoff_rmsd:
                # print('Filter out ' + emdId + ', rmsd: ' + str(temp_dict['rmsd']))
                filterListRMSD.append(emdId)
                continue
            if temp_dict['Z'] < cutoff_Z:
                # print('Filter out ' + emdId + ', Z: ' + str(temp_dict['Z']))
                filterListZ.append(emdId)
                continue
            if temp_dict['identity'] > cutoff_identity:
                # print('Filter out ' + emdId + ', identity: ' + str(temp_dict['identity']))
                filterListIdentiry.append(emdId)
                continue
            temp_diff = list(ref_indices_set - set(temp_dict['map_ref']))
            for diff_i in temp_diff:
                if not diff_i in missing_ind_dict:
                    missing_ind_dict[diff_i] = 1
                else:
                    missing_ind_dict[diff_i] += 1
        self._missing_ind_dict = missing_ind_dict
        filterList = filterListLen + filterListRMSD + filterListZ + filterListIdentiry
        filterDict = {'len': filterListLen, 'rmsd': filterListRMSD, 'Z': filterListZ, 'identity': filterListIdentiry}
        self._filterList = filterList
        self._filterDict = filterDict
        self._emdList = [self._emdListAll[0]] + list(set(list(self._emdListAll[1:])) - set(filterList))
        LOGGER.info(str(len(filterList)) + ' EMDs have been filtered out from '+str(len(emdListAll))+' Emsurfer hits (remaining: '+str(len(emdListAll)-len(filterList))+').')
        return self._emdList
    
    def getTitle(self):
        """Return the title of the record"""

        return self._title


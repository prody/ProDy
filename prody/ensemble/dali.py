# -*- coding: utf-8 -*-
"""This module defines functions for Dali searching Protein Data Bank."""

import re
from prody import LOGGER, PY3K
if PY3K:
    import urllib.parse as urllib
    import urllib.request as urllib2
else:
    import urllib
    import urllib2

import numpy as np
from prody import parsePDB, Ensemble, PDBEnsemble, AtomMap

__all__ = ['daliRecord', 'daliSearchPDB']
    
# from prody import *

def daliSearchPDB(pdbId, chainId, daliURL=None, subset='fullPDB', **kwargs):
    """Search Dali server with input of PDB ID and chain ID.
    Dali server: http://ekhidna2.biocenter.helsinki.fi/dali/
    
    :arg subset: fullPDB, PDB25, PDB50, PDB90
    :type sequence: str"""
    LOGGER.timeit('_dali')
    timeout = 120
    
    if daliURL is None:
        daliURL = "http://ekhidna2.biocenter.helsinki.fi/cgi-bin/sans/dump.cgi"
    pdbId = pdbId.lower()
    pdb_chain = pdbId + chainId
    parameters = { 'cd1' : pdb_chain, 'method': 'search', 'title': 'Title_'+pdb_chain, 'address': '' }
    enc_params = urllib.urlencode(parameters).encode('utf-8')
    request = urllib2.Request(daliURL, enc_params)
    try_error = 3
    while try_error >= 0:
        try:
            url = urllib2.urlopen(request).url
            break
        except:
            try_error -= 1
            if try_error >= 0:
                LOGGER.sleep(2, '. Connection error happened. Trying to reconnect...')
                continue
            else:
                url = urllib2.urlopen(request).url
                break
    if url.split('.')[-1].lower() in ['html', 'php']:
        url = url.strip(url.split('/')[-1])
    LOGGER.debug('Submitted Dali search for PDB and chain "{0} and {1}".'.format(pdbId, chainId))
    LOGGER.info(url)
    LOGGER.clear()
    obj = daliRecord(url, pdbId, chainId, subset=subset)
    if obj.isSuccess:
        return obj
    else:
        return None


class daliRecord(object):

    """A class to store results from Dali PDB search."""

    def __init__(self, url, pdbId=None, chainId=None, subset='fullPDB'):
        """Instantiate a daliPDB object instance.

        :arg subset: fullPDB, PDB25, PDB50, PDB90
        :type sequence: str"""

        self._url = url
        self._pdbId = pdbId
        # self._chainId = chainId
        self._chainId = chainId.upper()
        subset = subset.upper()
        if subset == "FULLPDB" or subset not in ["PDB25", "PDB50", "PDB90"]:
            self._subset = ""
        else:
            self._subset = "-"+subset[3:]
        
        self.isSuccess = self.getRecord(self._url)

    def getRecord(self, url):
        sleep = 2
        timeout = 120
        LOGGER.timeit('_dali')
        log_message = ''
        try_error = 3
        while True:
            LOGGER.sleep(int(sleep), 'to reconnect Dali '+log_message)
            LOGGER.clear()
            LOGGER.write('Connecting Dali for search results...')
            LOGGER.clear()
            try:
                html = urllib2.urlopen(url).read()
            except:
                try_error -= 1
                if try_error >= 0:
                    LOGGER.sleep(2, '. Connection error happened. Trying to reconnect...')
                    continue
                else:
                    html = urllib2.urlopen(url).read()
            if html.find('Status: Queued') > -1:
                log_message = '(Dali searching is queued)...'
            elif html.find('Status: Running') > -1:
                log_message = '(Dali searching is running)...'
            elif html.find('Your job') == -1 and html.find('.txt') > -1:
                break
            elif html.find('ERROR:') > -1:
                LOGGER.warn(': Dali search reported an ERROR!')
                return None
                break
            sleep = 20 if int(sleep * 1.5) >= 20 else int(sleep * 1.5)
            if LOGGER.timing('_dali') > timeout:
                LOGGER.warn(': Dali search is time out. \nThe results can be obtained using getRecord() function later.')
                return None
                break
            LOGGER.clear()
        LOGGER.clear()
        LOGGER.report('Dali results completed in %.1fs.', '_dali')
        lines = html.strip().split('\n')
        with open("temp.txt", "w") as file_temp: file_temp.write(html + '\n')
        file_name = re.search('=.+-90\.txt', html).group()[1:]
        file_name = file_name[:-7]
        # LOGGER.info(url+file_name+self._subset+'.txt')
        data = urllib2.urlopen(url+file_name+self._subset+'.txt').read()

        with open("temp.txt", "a+") as file_temp: file_temp.write(url+file_name + '\n' + data)
        data_list = data.strip().split('# ')
        # No:  Chain   Z    rmsd lali nres  %id PDB  Description -> data_list[3]
        # Structural equivalences -> data_list[4]
        # Translation-rotation matrices -> data_list[5]
        map_temp_dict = dict()
        mapping = []
        lines = data_list[4].strip().split('\n')
        self._lines_4 = lines
        mapping_temp = np.genfromtxt(lines[1:], delimiter = (4,1,14,6,2,4,4,5,2,4,4,3,5,4,3,5,6,3,5,4,3,5,28), usecols = [0,3,5,7,9,12,15,15,18,21], dtype='|i4')
        # [0,3,5,7,9,12,15,15,18,21] -> [index, residue_a, residue_b, residue_i_a, residue_i_b, resid_a, resid_b, resid_i_a, resid_i_b]
        for map_i in mapping_temp:
            if not map_i[0] in map_temp_dict:
                map_temp_dict[map_i[0]] = [[map_i[1], map_i[2], map_i[3], map_i[4]]]
            else:
                map_temp_dict[map_i[0]].append([map_i[1], map_i[2], map_i[3], map_i[4]])
        self._max_index = max(mapping_temp[:,2])
        self._mapping = map_temp_dict
        self._data = data_list[3]
        lines = data_list[3].strip().split('\n')
        daliInfo = np.genfromtxt(lines[1:], delimiter = (4,3,6,5,5,5,6,5,57), usecols = [0,2,3,4,5,6,7], dtype=[('id', '<i4'), ('pdb_chain', '|S6'), ('Z', '<f4'), ('rmsd', '<f4'), ('len_align', '<i4'), ('res_num', '<i4'), ('identity', '<i4')])
        if daliInfo.ndim == 0:
            daliInfo = np.array([daliInfo])
        pdbListAll = []
        self._daliInfo = daliInfo
        dali_temp_dict = dict()
        for temp in self._daliInfo:
            temp_dict = dict()
            pdb_chain = temp[1].strip()[0:6]
            temp_dict['pdbId'] = pdb_chain[0:4]
            temp_dict['chainId'] = pdb_chain[5:6]
            temp_dict['pdb_chain'] = pdb_chain
            temp_dict['Z'] = temp[2]
            temp_dict['rmsd'] = temp[3]
            temp_dict['len_align'] = temp[4]
            temp_dict['res_num'] = temp[5]
            temp_dict['identity'] = temp[6]
            temp_dict['mapping'] = (np.array(map_temp_dict[temp[0]])-1).tolist()
            temp_dict['map_ref'] = [x for map_i in (np.array(map_temp_dict[temp[0]])-1).tolist() for x in range(map_i[0], map_i[1]+1)]
            temp_dict['map_sel'] = [x for map_i in (np.array(map_temp_dict[temp[0]])-1).tolist() for x in range(map_i[2], map_i[3]+1)]
            dali_temp_dict[temp_dict['pdb_chain']] = temp_dict
            pdbListAll.append(pdb_chain)
        self._pdbListAll = tuple(pdbListAll)
        self._pdbList = self._pdbListAll
        self._alignPDB = dali_temp_dict

        LOGGER.info(str(len(pdbListAll)) + ' Dali results have been searched.')
        return True
        
    def getPDBList(self):
        """Returns PDB list (filters may applied)"""
        return self._pdbList
        
    def getPDBListAll(self):
        """Returns all PDB list"""
        return self._pdbListAll
        
    def getHits(self):
        return self._alignPDB
        
    def filter(self, cutoff_len=None, cutoff_rmsd=0.5, cutoff_Z=20, cutoff_identity=20):
        """Filters out PDBs from the PDBList and returns the PDB list.
        PDBs satisfy any of following criterion will be filtered out.
        (1) Length of aligned residues < cutoff_len;
        (2) RMSD < cutoff_rmsd;
        (3) Z score < cutoff_Z;
        (4) Identity < cutoff_identity.
        By default, cutoff_len is None and a cutoff of 0.8*Length will be applied.
        """
        daliInfo = self._alignPDB
        pdbListAll = self._pdbListAll
        missing_ind_dict = dict()
        ref_indices_set = set(range(self._max_index))
        filterList = []
        filterMissing = []
        if cutoff_len == None:
            cutoff_len = int(0.8*self._max_index)
        for pdb_chain in pdbListAll:
            temp_dict = daliInfo[pdb_chain]
            # filter: len_align, identity, rmsd, Z
            if temp_dict['len_align'] < cutoff_len:
                # print('Filter out ' + pdb_chain + ', len_align: ' + str(temp_dict['len_align']))
                filterList.append(pdb_chain)
                continue
            if temp_dict['rmsd'] < cutoff_rmsd:
                # print('Filter out ' + pdb_chain + ', rmsd: ' + str(temp_dict['rmsd']))
                filterList.append(pdb_chain)
                continue
            if temp_dict['Z'] < cutoff_Z:
                # print('Filter out ' + pdb_chain + ', Z: ' + str(temp_dict['Z']))
                filterList.append(pdb_chain)
                continue
            if temp_dict['identity'] < cutoff_identity:
                # print('Filter out ' + pdb_chain + ', identity: ' + str(temp_dict['identity']))
                filterList.append(pdb_chain)
                continue
            temp_diff = list(ref_indices_set - set(temp_dict['map_ref']))
            for diff_i in temp_diff:
                if not diff_i in missing_ind_dict:
                    missing_ind_dict[diff_i] = 1
                else:
                    missing_ind_dict[diff_i] += 1
        self._filterList = filterList
        self._missing_ind_dict = missing_ind_dict
        self._pdbList = list(set(list(self._pdbListAll)) - set(filterList))
        LOGGER.info(str(len(filterList)) + ' PDBs have been filtered out from '+str(len(pdbListAll))+' Dali hits.')
        return self._pdbList
        
    def buildDaliEnsemble(self):
        daliInfo = self._alignPDB
        pdbList = self._pdbList
        ref_pdb = parsePDB(self._pdbId).select('chain '+self._chainId).copy()
        ref_pdb_ca = ref_pdb.select("protein and name CA").copy()
        ref_chain = ref_pdb_ca.getHierView().getChain(self._chainId)
        ref_indices_set = set(range(len(ref_chain)))

        ensemble = PDBEnsemble('3h5v-A')
        ensemble.setAtoms(ref_chain)
        ensemble.setCoords(ref_chain)

        for pdb_chain in pdbList:
            print pdb_chain
            temp_dict = daliInfo[pdb_chain]
            sel_pdb = parsePDB(pdb_chain[0:4]).select('chain '+pdb_chain[5:6]).copy()
            sel_pdb_ca = sel_pdb.select("protein and name CA").copy()
            map_ref = []
            map_sel = []
            for i in range(len(temp_dict['map_ref'])):
                map_ref.append(temp_dict['map_ref'][i])
                map_sel.append(temp_dict['map_sel'][i])
            dum_sel = list(ref_indices_set - set(map_ref))
            atommap = AtomMap(sel_pdb_ca, indices=map_sel, mapping=map_ref, dummies=dum_sel)
            ensemble.addCoordset(atommap, weights=atommap.getFlags('mapped'))
            # ensemble.addCoordset(atommap, weights=atommap.getFlags('mapped'), degeneracy=True)
        ensemble.iterpose()
        RMSDs = ensemble.getRMSDs()
        return ensemble
    

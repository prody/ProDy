# -*- coding: utf-8 -*-
"""This module defines classes and functions for browsing QuartataWeb.

Based on code written by the CHARMM-GUI team (http://charmm-gui.org) and modified by James Krieger

This suite uses the following softwares:
- python Splinter package (https://splinter.readthedocs.org/en/latest/)
- a web browser, such as Google Chrome or Mozilla Firefox
- the corresponding driver such as chromedriver (https://sites.google.com/a/chromium.org/chromedriver/downloads)
   for Chrome or geckodriver (https://github.com/mozilla/geckodriver/releases) for Firefox
"""

from prody import PY3K, LOGGER
from prody.utilities import openFile

import numpy as np
import os

__all__ = ['QuartataWebBrowser', 'QuartataChemicalRecord', 'searchQuartataWeb',
           'initializeBrowser']


class QuartataWebBrowser(object):
    """Class to browse the QuartataWeb website.

    :arg data_source: source database for QuartataWeb analysis
        options are ``"DrugBank"`` or ``"STITCH"``. Default is ``"DrugBank"``
    :type data_source: str

    :arg drug_group: group of drugs if using DrugBank
        options are ``"Approved"`` or ``"All"``. Default is ``"All"``
    :type drug_group: str

    :arg input_type: number corresponding to the input type, options are 
        ``1`` (Chemical and/or target) or 
        ``2`` (A list of chemicals, targets or chemical combinations). 
        Default is ``1``
    :type input_type: int

    :arg query_type: number corresponding to the query type. Options are 
        dependent on input_type. 
        
        With input_type 1, they are:
        * ``1`` (chemical-target interaction)
        * ``2`` (chemical-chemical similarity)
        * ``3`` (target-target similarity)

        With input_type 2, they are:
        * ``1`` (chemicals)
        * ``2`` (targets)
        * ``3`` (chemical combinations)

        Default is ``1``
    :type query_type: int

    :arg data: data to enter into the box or boxes. This varies depending on input type 
        and query type, but will always be a list of strings.
        
        For input_type 1, a list with two items is expected. These will be one of the 
        following depending on query_type:
        * With query_type 1, the first would be a chemical and the second a target. 
            One of these can also be left blank.
        * With query_type 2, the first would be a chemical and the second a chemical.
        * With query_type 3, the first would be a target and the second a target.

        For input_type 2, a list with any length is expected. These will be one of the 
        following depending on query_type:
        * With query_type 1, these would be chemicals. 
        * With query_type 2, these would be targets.
        * With query_type 3, these would be pairs of chemicals, separated by semicolons.
    :type data: list

    :arg num_predictions: number of predictions to show or consider in addition to 
        known interactions. Default is ``0``. 
        With DrugBank and input_type 1, a second number can be provided in a list 
        for secondary interactions.
    :type num_predictions: int, list

    :arg browser_type: browser type for navigation
        Default is ``"Chrome"``
    :type browser_type: str

    :arg job_id: job ID for accessing previous jobs
        Default is **None**
    :type job_id: int        

    :arg tsv: a filename for a file that contains the results 
        or a file to save the results in tsv format
    :type tsv: str
    """

    def __init__(self, data_source=None, drug_group=None, input_type=None, query_type=None, 
                 data=None, num_predictions=None, browser_type=None, job_id=None, 
                 tsv=None, chem_type='known'):

        self.browser_type = None
        self.browser = None

        self.data_source = None
        self.drug_group = None
        self.input_type = None
        self.query_type = None
        self.data = None
        self.num_predictions = None

        self.chemical_data = {}
        self.fields = {}
        self.num_fields = {}
        self.num_rows = {}

        self.job_id = job_id

        self.filename = None
        self.no_data = {'known': True, 'predicted': True}
        if tsv is not None:
            try:
                self.parseChemicals(tsv, chem_type)
            except:
                raise ValueError('please provide a valid filename')

        self.setBrowserType(browser_type)
        self.setDataSource(data_source)
        self.setDrugGroup(drug_group)
        self.setInputType(input_type)
        self.setQueryType(query_type)
        self.setData(data)
        self.setNumPredictions(num_predictions)

    def updateHomePage(self):
        """Update the home page with data from setting variables"""
        url = "http://quartata.csb.pitt.edu"

        if self.data_source == 'DrugBank':
            url += '/index'
        else:
            url += '/index_stitch'

        if self.input_type == 2:
            url += '_2'

        url += '.php'
        self.browser.visit(url)

        if self.data_source == 'DrugBank':
            if self.drug_group == 'Approved':
                self.browser.find_by_name('db_type')[0].click()
            else:
                self.browser.find_by_name('db_type')[1].click()

        if self.query_type is not None:
            self.browser.find_by_name('pattern')[self.query_type - 1].click()

        if self.data is not None:
            if self.input_type == 1:
                if self.query_type == 1:
                    self.browser.find_by_name('q_drug_1')[0].fill(self.data[0])
                    self.browser.find_by_name('q_target_1')[
                        0].fill(self.data[1])
                elif self.query_type == 2:
                    self.browser.find_by_name('q_drug_1')[0].fill(self.data[0])
                    self.browser.find_by_name('q_drug_2')[0].fill(self.data[1])
                else:
                    self.browser.find_by_name('q_target_1')[0].fill(self.data[0])
                    self.browser.find_by_name('q_target_2')[0].fill(self.data[1])
            else:
                if self.query_type == 1:
                    self.browser.find_by_name('q_drugs')[0].fill('\n'.join(self.data))
                if self.query_type == 2:
                    self.browser.find_by_name('q_targets')[0].fill('\n'.join(self.data))
                else:
                    self.browser.find_by_name('q_drug_pairs')[0].fill('\n'.join(self.data))

        if self.num_predictions is not None:
            self.browser.find_by_name('pred_n')[0].fill(self.num_predictions[0])
            if self.data_source == 'DrugBank' and self.input_type == 1:
                self.browser.find_by_name('pred_n_2nd')[0].fill(self.num_predictions[1])

    def setDataSource(self, data_source):
        """Set data_source and update home page
        
        :arg data_source: source database for QuartataWeb analysis
            options are ``"DrugBank"`` or ``"STITCH"``. Default is ``"DrugBank"``
        :type data_source: str
        """
        if data_source is None:
            data_source = 'DrugBank'
        elif not isinstance(data_source, str):
            raise TypeError('data_source should be a string or None')
        elif data_source.lower() == 'drugbank':
            data_source = 'DrugBank'
        elif data_source.lower() == 'stitch':
            data_source = 'STITCH'
        else:
            raise ValueError('data_source should be DrugBank, STITCH or None')

        self.data_source = data_source
        if self.no_data:
            self.updateHomePage()

    def setDrugGroup(self, group):
        """Set drug_group and update home page
        
        :arg group: group of drugs if using DrugBank
            options are ``"Approved"`` or ``"All"``. Default is ``"All"``
        :type group: str
        """
        if self.data_source == 'DrugBank':
            if group is None:
                group = 'All'
            elif not isinstance(group, str):
                raise TypeError('group must be string or None')
            elif group.lower() == 'all':
                group = 'All'
            elif group.lower() == 'approved':
                group = 'Approved'
            else:
                raise ValueError('group should be approved, all or None')

            self.drug_group = group
            if self.no_data:
                self.updateHomePage()

        elif group is not None:
            LOGGER.warn('there are no groups when using STITCH')

    def setInputType(self, input_type):
        """Set input_type and update home page
        
        :arg input_type: number corresponding to the input type, options are 
            ``1`` (Chemical and/or target) or 
            ``2`` (A list of chemicals, targets or chemical combinations). 
            Default is ``1``
        :type input_type: int
        """        
        if input_type is None:
            input_type = 1
        elif not isinstance(input_type, int):
            raise TypeError('input_type should be an integer (1 or 2) or None')
        elif not input_type in [1, 2]:
            raise ValueError('input_type should be 1, 2 or None')

        self.input_type = input_type
        if self.no_data:
            self.updateHomePage()

    def setQueryType(self, query_type):
        """Set query_type and update home page
        
        :arg query_type: number corresponding to the query type. Options are 
            dependent on input_type. 
            
            With input_type 1, they are:
            * ``1`` (chemical-target interaction)
            * ``2`` (chemical-chemical similarity)
            * ``3`` (target-target similarity)

            With input_type 2, they are:
            * ``1`` (chemicals)
            * ``2`` (targets)
            * ``3`` (chemical combinations)

            Default is ``1``
        :type query_type: int
        """
        if query_type is None:
            query_type = 1
        elif not isinstance(query_type, int):
            raise TypeError(
                'query_type should be an integer (1, 2 or 3) or None')
        elif not query_type in [1, 2, 3]:
            raise ValueError('query_type should be 1, 2, 3 or None')

        self.query_type = query_type
        if self.no_data:
            self.updateHomePage()

    def setData(self, data):
        """Set data and update home page
        
        :arg data: data to enter into the box or boxes. This varies depending on input type 
            and query type, but will always be a list of strings.
            
            For input_type 1, a list with two items is expected. These will be one of the 
            following depending on query_type:
            * With query_type 1, the first would be a chemical and the second a target. 
                One of these can also be left blank.
            * With query_type 2, the first would be a chemical and the second a chemical.
            * With query_type 3, the first would be a target and the second a target.

            For input_type 2, a list with any length is expected. These will be one of the 
            following depending on query_type:
            * With query_type 1, these would be chemicals. 
            * With query_type 2, these would be targets.
            * With query_type 3, these would be pairs of chemicals, separated by semicolons.
        :type data: list
        """
        if data is None:
            LOGGER.warn('data is not set')
        elif not isinstance(data, list):
            raise TypeError('data should be a list')
        else:
            for item in data:
                if not isinstance(item, str):
                    raise TypeError('data should be a list of strings')

            if self.input_type == 1:
                if len(data) > 2:
                    raise ValueError(
                        'data can only have two values with input_type 1')

                if len(data) == 1:
                    data.append('')

            if self.input_type == 3:
                for item in data:
                    if item.find(';') == -1:
                        raise ValueError(
                            'each item in data must be a pair with ; as delimiter')

        self.data = data
        if self.no_data:
            self.updateHomePage()

    def setNumPredictions(self, num_predictions):
        """Set num_predictions and update home page
        
        :arg num_predictions: number of predictions to show or consider in addition to 
            known interactions. Default is ``0``. 
            With DrugBank and input_type 1, a second number can be provided in a list 
            for secondary interactions.
        :type num_predictions: int, list
        """
        if num_predictions is None:
            num_predictions = 0

        if not isinstance(num_predictions, (int, list)):
            raise TypeError(
                'num_predictions should be an integer, a list or None')

        if isinstance(num_predictions, int):
            num_predictions = [num_predictions, 0]

        if num_predictions[0] > 100:
            raise ValueError('1st num_predictions must be <= 100')

        if num_predictions[1] > 20:
            raise ValueError('2nd num_predictions must be <= 20')

        self.num_predictions = num_predictions
        if self.no_data:
            self.updateHomePage()

    def setBrowserType(self, browser_type):
        """Set browser_type and update home page
        
        :arg browser_type: browser type for navigation
            Default is ``"Chrome"``
        :type browser_type: str
        """
        if self.no_data:
            self.browser_type, self.browser = initializeBrowser(browser_type)
            self.updateHomePage()


    def setJObID(self, job_id):
        """Set job_id and view results
        
        :arg job_id: job ID for accessing previous jobs
            Default is **None**
        :type job_id: int
        """
        self.job_id = job_id
        if self.no_data:
            self.viewResults()


    def viewResults(self):
        """View results by clicking submit or using a job_id"""
        if self.job_id is None or self.browser.url.find('index') != -1:
            self.browser.find_by_name('submit')[0].click()
            self.job_id = self.browser.url.split('_')[-1].split('=')[-1]

        else:
            if self.data_source == 'DrugBank':
                url = '_'.join(['http://quartata.csb.pitt.edu/quartata_result.php?job',
                                'id={0}'.format(self.job_id)])
            else:
                url = '_'.join(['http://quartata.csb.pitt.edu/quartata_result',
                                'stitch.php?job', 'id={0}'.format(self.job_id)])

            self.browser.visit(url)


    def goToDownloads(self):
        """Go to downloads page"""
        if self.job_id is None:
            self.viewResults()

        if self.data_source == 'DrugBank':
            url = '_'.join(['http://quartata.csb.pitt.edu/quartata_download.php?job',
                            'id={0}'.format(self.job_id)])
        else:
            url = '_'.join(['http://quartata.csb.pitt.edu/quartata_download',
                            'stitch.php?job', 'id={0}'.format(self.job_id)])

        self.browser.visit(url)


    def goToWorkDir(self):
        """Go to working directory"""
        if self.job_id is None:
            self.viewResults()

        url = 'http://quartata.csb.pitt.edu/work/{0}'.format(self.job_id)
        self.browser.visit(url)


    def parseChemicals(self, filename=None, chem_type='known'):
        """Go to working directory and parse chemicals for query protein.
        Updates self.chemical_data"""
        
        if filename is None:
            filename = self.filename

        try:
            if filename is not None:
                if not self.no_data[chem_type]:
                    return True

                if not isinstance(filename, str):
                    raise TypeError('filename should be a string')

                if os.path.isfile(filename):
                    # read the contents
                    LOGGER.info('reading chemicals from {0}'.format(filename))
                    stream = openFile(filename, 'rt')
                    lines = stream.readlines()
                    stream.close()
                    self.no_data[chem_type] = False
                else:
                    # filename contains a filename for writing
                    self.no_data[chem_type] = True

                self.filename = filename

            if self.no_data[chem_type]:
                self.goToWorkDir()
                
                if self.data_source == 'DrugBank':
                    data_filename = '%s_drugs_for_query_protein.txt' % chem_type
                else:
                    data_filename = '%s_chemicals_for_query_protein.txt' % chem_type

                self.browser.find_by_text(data_filename)[0].click()
                
                import requests
                html = requests.get(self.browser.url).content
                if PY3K:
                    html = html.decode()

                if filename is not None:
                    LOGGER.info('writing chemicals to {0}'.format(filename))
                    out = open(filename, 'w')
                    out.write(html)
                    out.close()

                lines = html.split('\n')

            self.fields[chem_type] = lines[0].split('\t')
            self.num_fields[chem_type] = len(self.fields[chem_type])

            self.num_rows[chem_type] = len(lines[1:])
            if lines[-1].strip() == '':
                self.num_rows[chem_type] -= 1

            dtypes = []
            for i, item in enumerate(lines[1].split('\t')):
                if item.isnumeric():
                    dtypes.append((self.fields[chem_type][i], int))
                elif item.find('.') != -1 and item.replace('.','0').isnumeric():
                    dtypes.append((self.fields[chem_type][i], float))
                else:
                    dtypes.append((self.fields[chem_type][i], object))

            self.chemical_data[chem_type] = np.empty(self.num_rows[chem_type], dtype=dtypes)

            for i, line in enumerate(lines[1:self.num_rows[chem_type]+1]):
                items = line.strip().split('\t')
                if len(items) != self.num_fields[chem_type]:
                    raise ValueError('line {0} has the wrong number of fields'.format(i+1))

                for j, item in enumerate(items):
                    self.chemical_data[chem_type][i][j] = item
        except:
            self.no_data[chem_type] = True
        else:
            self.no_data[chem_type] = False

        return not self.no_data[chem_type]


    def quit(self):
        if self.browser is not None:
            self.browser.quit()


class QuartataChemicalRecord(object):
    """Class for handling chemical data from QuartataWebBrowser"""

    def __init__(self, data_source=None, drug_group=None, input_type=None, query_type=None, 
                 data=None, num_predictions=None, browser_type=None, job_id=None, 
                 filename=None):
        """Instantiate a QuartataChemicalRecord object instance.
        Inputs are the same as QuartataWebBrowser.
        """
        self._chemData = None
        self._filterDict = None
        self.data_source = data_source 
        self.drug_group = drug_group
        self.input_type = input_type
        self.query_type = query_type
        self.data = data
        self.num_predictions = num_predictions
        self.job_id = job_id
        self.filename = filename

        self.isSuccess = self.fetch(data_source, drug_group, input_type, query_type,
                                    data, num_predictions, browser_type, job_id, filename)


    def fetch(self, data_source=None, drug_group=None, input_type=None, query_type=None, 
              data=None, num_predictions=None, browser_type=None, job_id=None, filename=None):
        """Fetch data"""
        if data_source is None:
            data_source = self.data_source
        if drug_group is None:
            drug_group = self.drug_group
        if input_type is None:
            input_type = self.input_type
        if query_type is None:
            query_type = self.query_type
        if data is None:
            data = self.data

        if data is None:
            raise ValueError('data cannot be None')

        if num_predictions is None:
            num_predictions = self.num_predictions
        if job_id is None:
            job_id = self.job_id
        if filename is None:
            filename = self.filename

        self.qwb = QuartataWebBrowser(data_source, drug_group, input_type, query_type,
                                      data, num_predictions, browser_type, job_id, filename)
        
        isSuccess = self.qwb.parseChemicals()
        if self.qwb.num_predictions[0] > 0:
            isSuccess = self.qwb.parseChemicals(chem_type='predicted')

        self.qwb.quit()

        self._chemData = self.qwb.chemical_data
        if self._chemData is None:
            raise ValueError('')
        chem_temp_dict = dict()
        listAll = []
        for key in self._chemData:
            for temp in self._chemData[key]:
                temp_dict = dict()
                chem_name = temp[1]

                temp_dict['DB_ID'] = temp[0]
                temp_dict['chemical_name'] = chem_name
                temp_dict['mol_weight'] = temp[2]
                temp_dict['SMILES'] = temp[3]
                temp_dict['conf_score'] = temp[4]

                chem_temp_dict[chem_name] = temp_dict
                listAll.append(chem_name)

            self._listAll = tuple(listAll)
            self._list = self._listAll
            self._chemDict = chem_temp_dict
        
        return isSuccess


    def getChemicalList(self, filtered=True):
        """Returns chemical list (filters may be applied)"""
        if not self.isSuccess:
            LOGGER.warn('Quartata Chemical Record does not have any data yet. '
                        'Please run fetch again, possibly with different parameters.')
        
        if filtered:
            return self._list
        return self._listAll
        

    def getSMILESList(self, filtered=True):
        """Returns SMILES list (filters may be applied)"""
        if not self.isSuccess:
            LOGGER.warn('Quartata Chemical Record does not have any data yet.'
                        'Please run fetch again, possibly with different parameters.')
        
        if filtered:
            return [self._chemDict[key]['SMILES'] for key in self._list]
        return self._chemData['SMILES']
        

    def getParticularSMILES(self, key):
        """Returns SMILES for a particular chemical"""
        if not self.isSuccess:
            LOGGER.warn('Quartata Chemical Record does not have any data yet.'
                        'Please run fetch again, possibly with different parameters.')

        return self._chemDict[key]['SMILES']


    def getFilterList(self):
        """Returns a list of chemicals for the entries that were filtered out"""
        
        filterDict = self._filterDict
        if filterDict is None:
            raise ValueError('You cannot obtain the list of filtered out entries before doing any filtering.')

        temp_str = ', '.join([str(len(filterDict['lower_MW'])), str(len(filterDict['upper_MW'])), 
                              str(len(filterDict['conf_score']))])
        LOGGER.info('Filtered out [' + temp_str + '] for [lower weight, upper weight, confidence score]')
        return self._filterList


    def filter(self, lower_weight=None, upper_weight=None, cutoff_score=None):
        """Filters out chemicals from the list and returns the updated list.
        Chemicals that satisfy any of the following criterion will be filtered out.
        (1) Molecular weight < lower_weight (must be a positive number);
        (2) Molecular weight > upper_weight (must be a positive number);
        (3) Confidence score < cutoff_score (must be a positive number);

        Please note that every time this function is run, this overrides any previous runs.
        Therefore, please provide all filters at once.
        """
        if not self.isSuccess:
            LOGGER.warn('Quartata Chemical Record does not have any data yet.'
                        'Please run fetch again, possibly with different parameters.')
            return None

        if lower_weight == None:
            lower_weight = 0
        elif not isinstance(lower_weight, (float, int)):
            raise TypeError('lower_weight must be a float or an integer')
        if lower_weight >= 0:
            lower_weight = float(lower_weight)
        else:
            raise ValueError('lower_weight must be a number not less than 0')
            
        if upper_weight == None:
            upper_weight = 0
        elif not isinstance(upper_weight, (float, int)):
            raise TypeError('upper_weight must be a float or an integer')
        if upper_weight >= 0:
            upper_weight = float(upper_weight)
        else:
            raise ValueError('upper_weight must be a number not less than 0')
            
        if cutoff_score == None:
            cutoff_score = 0
        elif not isinstance(cutoff_score, (float, int)):
            raise TypeError('cutoff_score must be a float or an integer')
        elif cutoff_score >= 0:
            cutoff_score = float(cutoff_score)
        else:
            raise ValueError('cutoff_score must be a number not less than 0')

        quartataInfo = self._chemDict
        if quartataInfo is None:
            raise ValueError("Quartata Chemical Record does not have any data yet. Please run fetch.")

        listAll = self._listAll
        ref_indices_set = set(range(self.qwb.num_rows))
        filterListLowerMW = []
        filterListUpperMW = []
        filterListConf = []
        
        for chem in listAll:
            temp_dict = quartataInfo[chem]

            if temp_dict['mol_weight'] < lower_weight:
                filterListLowerMW.append(chem)
                continue

            if upper_weight > 0 and temp_dict['mol_weight'] > upper_weight:
                filterListUpperMW.append(chem)
                continue

            if temp_dict['conf_score'] < cutoff_score:
                filterListConf.append(chem)
                continue

        filterList = filterListLowerMW + filterListUpperMW + filterListConf
        filterDict = {'lower_MW': filterListLowerMW, 'upper_MW': filterListUpperMW, 'conf_score': filterListConf}
        self._filterList = filterList
        self._filterDict = filterDict
        self._list = [item for item in self._listAll if not item in filterList]
        LOGGER.info(str(len(self._listAll)-len(self._list)) + ' chemicals have been filtered out from '+str(len(self._listAll))+' QuartataWeb hits (remaining: '+str(len(self._list))+').')
        return self._list
    

def searchQuartataWeb(data_source=None, drug_group=None, input_type=None, query_type=None, 
                      data=None, num_predictions=None, browser_type=None, job_id=None, 
                      filename=None, result_type='Chemical'):
    """Wrapper function for searching QuartataWeb.

    :arg result_type: type of results to get from QuartataWeb.
        So far only ``'Chemical'`` is supported.
    :type result_type: str

    All other arguments are the same as :class:`.QuartataWebBrowser`.
    """
    if result_type == 'Chemical':
        return QuartataChemicalRecord(data_source, drug_group, input_type, query_type,
                                      data, num_predictions, browser_type, job_id,
                                      filename)
    else:
        LOGGER.warn('No other result types are supported yet')
        return None


def initializeBrowser(browser_type, url):
    try:
        from splinter import Browser
    except ImportError:
        raise ImportError('Browser module could not be imported. '
                            'install splinter package to solve the problem.')
    else:
        from selenium.webdriver.common.service import WebDriverException
    
    if url is None:
        url = "http://quartata.csb.pitt.edu"

    if browser_type is None:
        try:
            browser = Browser('chrome')
            browser.visit(url)
        except WebDriverException:
            try:
                browser = Browser('firefox')
                browser.visit(url)
            except WebDriverException:
                raise ValueError('No web driver found for Chrome or Firefox. '
                                    'Please specify a different browser type or download an appropriate driver.')
            else:
                browser_type = 'firefox'
        else:
            browser_type = 'chrome'

    elif not isinstance(browser_type, str):
        raise TypeError('browser_type should be a string or None')
    else:
        try:
            browser = Browser(browser_type)
            browser.visit(url)
        except WebDriverException:
            raise ValueError('No web driver found for browser_type. '
                                'Please specify a different browser type or download an appropriate driver.')
        else:
            browser_type = browser_type

    return browser_type, browser
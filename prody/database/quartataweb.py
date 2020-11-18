# -*- coding: utf-8 -*-
"""This module defines classes and functions for browsing QuartataWeb.

----------------------------------------------------------------------------------------
Based on code written by the CHARMM-GUI team (http://charmm-gui.org) and modified by James Krieger

This suite uses the following softwares:
a) python Splinter package (https://splinter.readthedocs.org/en/latest/)
b) a web browser, such as Google Chrome or Mozilla Firefox
c) the corresponding driver such as chromedriver (https://sites.google.com/a/chromium.org/chromedriver/downloads)
   for Chrome or geckodriver (https://github.com/mozilla/geckodriver/releases) for Firefox
----------------------------------------------------------------------------------------
"""

from prody import PY3K, LOGGER
import numpy as np

try:
    from splinter import Browser
except ImportError:
    raise ImportError('Browser module could not be imported. '
                      'install splinter package to solve the problem.')
else:
    from selenium.webdriver.common.service import WebDriverException

import requests

__all__ = ['QuartataWebBrowser']


class QuartataWebBrowser(object):
    """A class to browse the QuartataWeb website."""

    def __init__(self, data_source=None, drug_group=None, input_type=None, query_type=None, 
                 data=None, num_predictions=None, browser_type=None, job_id=None):
        """Instantiate a QuartataWebBrowser object instance.

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
            Default is ``None``
        :type browser_type: int        
        """

        self.browser_type = None
        self.browser = None

        self.data_source = None
        self.drug_group = None
        self.input_type = None
        self.query_type = None
        self.data = None
        self.num_predictions = None

        self.job_id = job_id

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
        else:
            LOGGER.warn('there are no groups when using STITCH')

        self.drug_group = group
        self.updateHomePage()

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
        self.updateHomePage()

    def setBrowserType(self, browser_type):
        """Set browser_type and update home page
        
        :arg browser_type: browser type for navigation
            Default is ``"Chrome"``
        :type browser_type: str
        """
        if browser_type is None:
            try:
                browser = Browser('chrome')
                url = "http://quartata.csb.pitt.edu"
                browser.visit(url)
            except WebDriverException:
                try:
                    browser = Browser('firefox')
                    url = "http://quartata.csb.pitt.edu"
                    browser.visit(url)
                except WebDriverException:
                    raise ValueError(
                        'No web driver found for Chrome or Firefox. Please specify a browser type or download an appropriate driver.')
                else:
                    self.browser_type = 'firefox'
            else:
                self.browser_type = 'chrome'

        elif not isinstance(browser_type, str):
            raise TypeError('browser_type should be a string or None')
        else:
            try:
                browser = Browser(browser_type)
                url = "http://quartata.csb.pitt.edu"
                browser.visit(url)
            except WebDriverException:
                raise ValueError(
                    'No web driver found for browser_type. Please specify a different browser type or download an appropriate driver.')
            else:
                self.browser_type = browser_type

        self.browser = browser
        self.updateHomePage()


    def setJObID(self, job_id):
        """Set job_id and view results
        
        :arg job_id: job ID for accessing previous jobs
            Default is ``None``
        :type browser_type: int
        """
        self.job_id = job_id
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


    def parseChemicals(self):
        """Go to working directory and parse chemicals for query protein.
        Updates self.chemical_data"""
        self.goToWorkDir()
        
        if self.data_source == 'DrugBank':
            filename = 'known_drugs_for_query_protein.txt'
        else:
            filename = 'known_chemicals_for_query_protein.txt'

        self.browser.find_by_text(filename)[0].click()
        
        html = requests.get(self.browser.url).content
        if PY3K:
            html = html.decode()

        lines = html.split('\n')

        self.fields = lines[0].split('\t')
        self.num_fields = len(self.fields)

        self.num_rows = len(lines[1:])
        if lines[-1].strip() == '':
            self.num_rows -= 1

        dtypes = []
        for i, item in enumerate(lines[1].split('\t')):
            if item.isnumeric():
                dtypes.append((self.fields[i], int))
            elif item.find('.') != -1 and item.replace('.','0').isnumeric():
                dtypes.append((self.fields[i], float))
            else:
                dtypes.append((self.fields[i], object))

        self.chemical_data = np.empty(self.num_rows, dtype=dtypes)

        for i, line in enumerate(lines[1:self.num_rows]):
            items = line.split('\t')
            if len(items) != self.num_fields:
                raise ValueError('line {0} has the wrong number of fields'.format(i+1))

            for j, item in enumerate(items):
                self.chemical_data[i][j] = item


    def quit(self):
        self.browser.quit()


class QuartataWebRecord(object):
    pass
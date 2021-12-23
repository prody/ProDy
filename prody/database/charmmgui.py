# -*- coding: utf-8 -*-
"""This module defines classes and functions for using CHARMM-GUI.

Written by CHARMM-GUI team (www.charmm-gui.org) and modified by James Krieger

This suite uses the following softwares:
a) python splinter (https://splinter.readthedocs.org/en/latest/) and selenium (https://selenium-python.readthedocs.io/) packages
b) a web browser, such as Google Chrome or Mozilla Firefox
c) the corresponding driver such as chromedriver (https://sites.google.com/a/chromium.org/chromedriver/downloads)
   for Chrome or geckodriver (https://github.com/mozilla/geckodriver/releases) for Firefox
"""

from prody.proteins.pdbfile import parsePDB, fetchPDB
from prody import LOGGER
from .quartataweb import initializeBrowser

import os
import numpy as np
from numbers import Integral
import time


__all__ = ['CharmmGUIBrowser']


class CharmmGUIBrowser(object):
    """Class for running CHARMM-GUI in a browser

    :arg fname: filename for starting PDB file. A PDB ID can be provided and ProDy will fetch the PDB file first.
    :type fname: str

    :arg cwd: current working directory for CHARMM-GUI. Default is the current working directory where Python is running.
    :type cwd: str

    :arg ndir: relative path for next directory for accessing and saving files. Default is an empty string meaning use cwd.
    :type ndir: str

    :arg segids: segment IDs for selection in CHARMM-GUI. Default is to select all protein chains.
    :type segids: list

    :arg job_type: type of job for CHARMM-GUI. Currently `'solution'` and `'membrane.bilayer'` are supported.
        Default is `'solution'`.
    :type job_type: str

    :arg job_id: job ID for retrieving previous runs. Default is empty string.
    :type job_id: str

    :arg timeout: number of seconds before timing out when waiting for each step.
        Default is 12000, corresponding to 3 hours and 20 minutes, which seems high enough to not time out.
    :type timeout: int

    :arg saveas: tgz folder name for saving the CHARMM-GUI outputs
    :type saveas: str

    :arg browser_type: type of browser for running CHARMM-GUI. 
        Default is to find one that has compatible drivers.
    :type browser_type: str

    :arg max_warn: maximum number of warnings before raising an error.
        Default is 0.
    :type max_warn: int 
    """
    def __init__(self, fname=None, **kwargs):
        if fname is None:
            raise ValueError('Please enter a starting PDB filename.')
        if not isinstance(fname, str):
            raise TypeError('fname should be a string')

        if not fname.endswith('.pdb'):
            fname += '.pdb'

        self.fname = fname

        self.cwd = kwargs.get('cwd', os.getcwd())
        if not isinstance(self.cwd, str):
            raise TypeError('cwd should be a string')

        self.ndir = kwargs.get('ndir', '')
        if not isinstance(self.ndir, str):
            raise TypeError('ndir should be a string')

        if not os.path.isfile("%s/%s/%s" % (self.cwd, self.ndir, fname)):
            fpath = fetchPDB(fname[:-4], compressed=False)
            if not os.path.isdir("%s/%s" % (self.cwd, self.ndir)):
                os.mkdir("%s/%s" % (self.cwd, self.ndir))
            os.rename(fpath, "%s/%s/%s" % (self.cwd, self.ndir, fname))

        self.ag = parsePDB("%s/%s/%s" % (self.cwd, self.ndir, fname),
                           subset='ca')

        self.segids = kwargs.get(
            'segids', ['PRO'+ch.getChid() for ch in self.ag.iterChains()])
        if not isinstance(self.segids, list):
            raise TypeError('segids should be a list')

        self.job_id = kwargs.get('job_id', '')

        self.job_type = kwargs.get('job_type', 'solution')
        if not isinstance(self.job_type, str):
            raise TypeError('job_type should be a string')
        if not self.job_type in ['solution', 'membrane.bilayer']:
            raise ValueError(
                "job_type should be `'solution'` or `'membrane.bilayer'`. No other types are supported yet.")

        self.timeout = kwargs.get('timeout', 12000)
        if not isinstance(self.timeout, Integral):
            raise TypeError('timeout should be an integer')

        self.browser = None
        self.browser_type = kwargs.get('browser_type', None)

        self.link = None

        self.saveas = kwargs.get('saveas', 'charmm-gui.tgz')
        if not isinstance(self.saveas, str):
            raise TypeError('saveas should be a string')

        self.isSuccess = self.run(**kwargs)

    def download(self, browser=None, link=None, saveas=None):
        """Download CHARMM-GUI files from *link* to *saveas* using *browser*. 
        If not set, these taken from the CharmmGUIBrowser object."""

        if browser is None:
            browser = self.browser

        if link is None:
            link = self.link

        if saveas is None:
            saveas = self.saveas

        saveas = self.cwd + '/' + self.ndir + '/' + saveas

        import requests
        LOGGER.info("downloading %s to %s" % (link, saveas))
        response = requests.get(link, stream=True)
        if response.status_code == 200:
            with open(saveas, "wb") as code:
                code.write(response.raw.read(decode_content=True))

        fsize = float(os.stat(saveas).st_size)/(1024.0*1024.0)
        LOGGER.info("download complete, file size is %5.2f MB" % fsize)

    def get_download_link(self, browser=None):
        """Get download link for CHARMM-GUI files using *browser*. 
        If not set, *browser* is taken from the CharmmGUIBrowser object."""
        if browser is None:
            browser = self.browser

        sleep = 2
        while True:
            try:
                downlink = browser.links.find_by_partial_href("download")[-1]
                d = downlink['href']
                if d == None:
                    continue
                else:
                    return d
            except:
                LOGGER.info(
                    "Waiting for download link, sleep for %s s" % sleep)
                LOGGER.sleep(int(sleep), 'to check for download link.')
                sleep = int(sleep * 1.5)

    def check_error(self, browser=None, **kwargs):
        """Check errors on current *browser* page. If not set, this is taken from the CharmmGUIBrowser object.

        :arg max_warn: number of warnings to allow
            Default is 0
        :type max_warn: int
        """
        if browser is None:
            browser = self.browser

        max_warn = kwargs.get('max_warn', 0)
        if not isinstance(max_warn, Integral):
            raise TypeError('max_warn should be an integer')

        sleep = 2
        while True:
            try:
                res = not (browser.is_text_present("CHARMM was terminated abnormally")
                           or len(browser.find_by_id("error_msg")[:]) > max_warn)
                if res:
                    return res, None
                else:
                    return res, browser.find_by_id("error_msg")[0].text
            except:
                LOGGER.info("Waiting for error check, sleep for %s s" % sleep)
                LOGGER.sleep(int(sleep), 'to check for error.')
                sleep = int(sleep * 1.5)

    def wait_for_text(self, browser=None, text=None, timeout=None, **kwargs):
        """Wait for text (using the content of a button) when going to the next page."""
        if browser is None:
            browser = self.browser

        if text is None:
            text = self.text

        if timeout is None:
            timeout = self.timeout

        LOGGER.info("    Next step: %s" % text)
        sleep_time = 0
        sleep = 2
        while sleep_time < timeout:
            try:
                res = browser.is_text_present(text)
                if res == True:
                    status, err = self.check_error(browser, **kwargs)
                    if status == True or err != 'The CHARMM-GUI process is still running. Please reload this page again a few minutes later.':
                        break
                else:
                    LOGGER.sleep(int(sleep), 'for %s.' % text)
                    sleep = int(sleep * 1.5)
                    sleep_time += sleep
            except:
                LOGGER.sleep(int(sleep), 'for %s.' % text)
                sleep = int(sleep * 1.5)
                sleep_time += sleep

            if sleep_time > 10 and self.job_id is not None:
                if browser.url.find(self.job_id) == -1:
                    browser.visit('http://charmm-gui.org/?doc=input/retriever')
                    browser.find_by_name('jobid')[0].fill(self.job_id)
                    browser.find_by_value("Submit")[0].click()
                    browser.links.find_by_text("Go")[-1].click()
                else:
                    browser.reload()
                    status, _ = self.check_error(browser, **kwargs)
                    if status == True:
                        browser.execute_script("proceed()")

            if sleep > 180:
                sleep = 180

        if sleep_time >= timeout:
            LOGGER.warn('wait for %s step timed out' % text)
            return False
        
        return True

    def next_step(self, browser=None, text=None, timeout=None, **kwargs):
        """Go to next step including waiting for text and checking for errors."""
        if browser is None:
            browser = self.browser

        if text is None:
            text = self.text

        browser.execute_script("proceed()")
        status = True

        LOGGER.info("    Going to next step")
        status = self.wait_for_text(browser, text, timeout)
        if status == False:
            return False

        status, err = self.check_error(browser, **kwargs)
        if status == False:
            LOGGER.warn(err)
            return False

        return True

    def run(self, cwd=None, ndir=None, fname=None, segids=None, job_type=None, saveas=None, timeout=None, browser_type=None, **kwargs):
        """Main running function that runs all the other steps. Arguments are as in initialization and can overwrite them."""
        if cwd is None:
            cwd = self.cwd

        if ndir is None:
            ndir = self.ndir

        if fname is None:
            fname = self.fname

        if segids is None:
            segids = self.segids

        if job_type is None:
            job_type = self.job_type

        if not job_type in ['solution', 'membrane.bilayer']:
            raise ValueError('job_type must be solution or membrane.bilayer')

        LOGGER.info("Running CHARMM-GUI in %s directory"
                    % ndir)

        if saveas is None:
            saveas = self.saveas

        if timeout is None:
            timeout = self.timeout
        else:
            self.timeout = timeout

        if browser_type is None:
            browser_type = self.browser_type

        url = "http://www.charmm-gui.org/input/%s" % job_type
        self.browser_type, self.browser = initializeBrowser(browser_type, url=url)
        browser = self.browser
        browser.visit(url)
        
        browser.find_by_id('email').first.fill('kriegerj@pitt.edu') # we will have to change this
        browser.find_by_id('password').first.fill('ProDy123')
        browser.find_by_value("Submit").click()
        
        browser.visit(url)
        browser.attach_file('file', "%s/%s/%s" % (cwd, ndir, fname))
        browser.find_by_value("PDB").click()

        status = self.next_step(browser, "Manipulate PDB", **kwargs)
        if status == False:
            return False

        for item in browser.find_by_value("1"):
            item.uncheck()

        for segid in segids:
            browser.find_by_name("chains[%s][checked]" % segid.upper()).check()

        status = self.next_step(browser, "Generate PDB", **kwargs)
        if status == False:
            return False

        if job_type == 'membrane.bilayer':
            status = self.next_step(browser, "Calculate Cross-Sectional Area", **kwargs)
            if status == False:
                return False

            browser.find_by_name("align_option").first.click()
            status = self.next_step(browser, "Determine the System Size", **kwargs)
            if status == False:
                return False

            maxes = np.max(self.ag.getCoords(), axis=0)
            minis = np.min(self.ag.getCoords(), axis=0)
            len_x = maxes[0] - minis[0]
            len_y = maxes[1] - minis[1]
            max_xy = np.ceil(np.max([len_x, len_y]))
            browser.find_by_name("hetero_lx").first.fill(str(max_xy+10))

            browser.find_by_id("lipid_ratio[upper][chl1]").first.fill("1")
            browser.find_by_id("lipid_ratio[lower][chl1]").first.fill("1")

            [y for y in browser.find_by_css("img") if y.outer_html.find(
                "liptype_ratio_pc") != -1][0].click()
            time.sleep(5)
            browser.find_by_id("lipid_ratio[upper][popc]").first.fill("3")
            browser.find_by_id("lipid_ratio[lower][popc]").first.fill("3")

            browser.find_by_value("Show the system info").first.click()

            while len(browser.find_by_id("error_msg")[:]) > 0:
                if browser.find_by_id("error_msg")[0].text.find("more lipids") != -1:
                    browser.find_by_name("hetero_xy_option")[1].click()
                    time.sleep(5)
                    [y for y in browser.find_by_css("img")
                     if y.outer_html.find("liptype_number_pc") != -1][0].click()
                    time.sleep(5)

                if browser.find_by_id("error_msg")[0].text.find("lower") != -1:
                    lnum_popc = browser.find_by_id("lipid_number[lower][popc]").first
                    lnum_popc.fill(str(int(lnum_popc.value)+1))
                elif browser.find_by_id("error_msg")[0].text.find("upper") != -1:
                    unum_popc = browser.find_by_id("lipid_number[upper][popc]").first
                    unum_popc.fill(str(int(unum_popc.value)+1))
                else:
                    LOGGER.warn(browser.find_by_id("error_msg")[0].text)
                    return False

                browser.find_by_value("Show the system info")[1].click()

            status = self.next_step(browser, "Build Components", max_warn=1)
            if status == False:
                return False
        else:
            status = self.next_step(browser, "Solvate Molecule", **kwargs)
            if status == False:
                return False

        browser.find_option_by_text('NaCl').first.click()
        browser.find_option_by_text('Monte-Carlo').first.click()

        self.job_id = browser.find_by_name("jobid")[0].value
        if job_type == 'membrane.bilayer':
            status = self.next_step(browser, "Assemble Components", **kwargs)
            if status == False:
                return False
        else:
            status = self.next_step(browser, "Setup Periodic Boundary Condition", **kwargs)
            if status == False:
                return False

        status = self.next_step(browser, "Generate Equilibration and Dynamics", **kwargs)
        if status == False:
            return False

        browser.find_by_name("namd_checked").check()

        status = self.next_step(browser, "production.inp", **kwargs)
        if status == False:
            return False

        link = self.get_download_link(browser)
        LOGGER.info("Build success")
        status = True

        self.browser = browser
        self.link = link

        self.download(saveas=saveas)
        self.browser.quit()
        return status


if __name__ == "__main__":

    cgb = CharmmGUIBrowser(fname='1ake',
                           saveas='1ake-charmm-gui.tgz',
                           ndir='test')

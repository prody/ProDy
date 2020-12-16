#!/usr/bin/env python
# ----------------------------------------------------------------------------------------
# Written by CHARMM-GUI team (www.charmm-gui.org) and modified by James Krieger
# This suite uses the following softwares:
# a) google chrome browser
# b) python Splinter package (https://splinter.readthedocs.org/en/latest/)
# c) chromedriver (https://sites.google.com/a/chromium.org/chromedriver/downloads).
#    chromedriver is required by Splinter. Make sure it is installed in your PATH
# ----------------------------------------------------------------------------------------

from prody.proteins.pdbfile import parsePDB
import time
import os
import sys
from prody import PY3K, LOGGER

if PY3K:
    import urllib.parse as urllib
    import urllib.request as urllib2
else:
    import urllib
    import urllib2

import subprocess


__all__ = ['CharmmGUIBrowser']


class CharmmGUIBrowser(object):
    def __init__(self, cwd=None, ndir=None, fname=None, segids=None):
        if cwd is None:
            cwd = os.getcwd()
        self.cwd = cwd

        if ndir is None:
            ndir = ''
        self.ndir = ndir

        if not isinstance(fname, str):
            raise TypeError('please provide a filename')
        self.fname = fname

        ag = parsePDB("%s/%s/%s" % (cwd, ndir, fname),
                      subset='ca')

        if segids is None:
            segids = ['PRO'+ch.getChid() for ch in ag.iterChains()]
        self.segids = segids

        self.browser = None
        self.link = None
        self.status = False

        self.saveas = 'charmm-gui'

        # self.run()

    def download(self, browser=None, link=None, saveas=None):
        if browser is None:
            browser = self.browser

        if link is None:
            link = self.link

        if saveas is None:
            saveas = self.saveas

        LOGGER.info("downloading %s to %s" % (link, saveas))
        url = "http://www.charmm-gui.org/?doc=input/download"
        user_agent = 'Mozilla/4.0 (compatible; MSIE 5.5; Windows NT)'
        headers = {'User-Agent': user_agent}
        link1 = link.split("&")
        tag = link1[1].split("=")[1]
        step = link1[2].split("=")[1]
        values = {'time': tag,
                  'step': step,
                  'archive': 'tgz'}
        data = urllib.urlencode(values)
        req = urllib2.Request(url, data, headers)
        response = urllib2.urlopen(req)
        the_page = response.read()
        with open(saveas, "wb") as code:
            code.write(the_page)
        fsize = float(os.stat(saveas).st_size)/(1024.0*1024.0)
        LOGGER.info("download complete, file size is %5.2f MB" % fsize)

    def get_download_link(self, browser=None):
        if browser is None:
            browser = self.browser

        while True:
            try:
                downlink = browser.find_link_by_partial_href("archive=tgz")[-1]
                d = downlink['href']
                if d == None:
                    time.sleep(1)
                    continue
                else:
                    return d
            except:
                LOGGER.info("Waiting for download link, sleep for 2 s")
                time.sleep(2)

    def check_error(self, browser=None):
        if browser is None:
            browser = self.browser

        while True:
            try:
                res = browser.is_text_present(
                    "CHARMM was terminated abnormally")
                if res == True:
                    LOGGER.info("Build Failed")
                    return False
                if res == False:
                    return True
            except:
                LOGGER.info("check_error sleep for 1 second")
                time.sleep(1)

    def wait_for_text(self, browser=None, text=None):
        if browser is None:
            browser = self.browser

        if text is None:
            text = self.text

        LOGGER.info("    Waiting for %s button" % text)
        while True:
            try:
                res = browser.is_text_present(text)
                if res == True:
                    break
                #err = browser.
            except:
                LOGGER.info("    sleep for 2 seconds")
                time.sleep(2)

    def next_step(self, browser=None, text=None):
        if browser is None:
            browser = self.browser

        if text is None:
            text = self.text

        browser.execute_script("proceed()")
        if len(browser.find_by_id("error_msg")[:]) > 0:
            raise Exception(browser.find_by_id("error_msg")[0].text)

        LOGGER.info("    Going to next step")
        self.wait_for_text(browser, text)

    def run(self, cwd=None, ndir=None, fname=None, segids=None):
        if cwd is None:
            cwd = self.cwd

        if ndir is None:
            ndir = self.ndir

        if fname is None:
            fname = self.fname

        if segids is None:
            segids = self.segids

        LOGGER.info(
            "Running CHARMM-GUI Solution Builder in %s directory" % ndir)

        try:
            from splinter import Browser
        except ImportError:
            raise ImportError('Browser module could not be imported. '
                              'install splinter package to solve the problem.')

        self.browser = browser = Browser('chrome')
        url = "http://www.charmm-gui.org/input/solution"
        browser.visit(url)

        browser.attach_file('file', "%s/%s/%s" % (cwd, ndir, fname))
        browser.find_by_value("PDB").click()

        self.next_step(browser, "Manipulate PDB")

        for item in browser.find_by_value("1"):
            item.uncheck()

        for segid in segids:
            browser.find_by_name("chains[%s][checked]" % segid.upper()).check()

        self.next_step(browser, "Generate PDB")

        self.next_step(browser, "Solvate Molecule")
        if self.check_error(browser) == False:
            return browser, self.get_download_link(browser)  # , status
        browser.find_option_by_text('NaCl').first.click()
        browser.find_option_by_text('Monte-Carlo').first.click()

        self.next_step(browser, "Setup Periodic Boundary Condition")
        if self.check_error(browser) == False:
            return browser, self.get_download_link(browser)  # , status

        self.next_step(browser, "Generate Equilibration and Dynamics Inputs")
        if self.check_error(browser) == False:
            return browser, self.get_download_link(browser)  # , status

        self.next_step(browser, "step5_production.inp")

        self.check_error(browser)
        link = self.get_download_link(browser)
        LOGGER.info("Build success")
        status = "Success"

        self.browser = browser
        self.link = link
        self.status = status
        return


if __name__ == "__main__":

    from prody import fetchPDB
    fname = fetchPDB('1ake')
    cgb = CharmmGUIBrowser(fname=fname)
    browser, link, status = cgb.run()

    if status == "Failed":
        LOGGER.info("Failed")
    else:
        LOGGER.info("Success")
        cgb.download

    browser.quit()

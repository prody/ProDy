#!/usr/bin/env python
# ----------------------------------------------------------------------------------------
# Written by CHARMM-GUI team (www.charmm-gui.org) and modified by James Krieger
# This suite uses the following softwares:
# a) google chrome browser
# b) python Splinter package (https://splinter.readthedocs.org/en/latest/)
# c) chromedriver (https://sites.google.com/a/chromium.org/chromedriver/downloads).
#    chromedriver is required by Splinter. Make sure it is installed in your PATH
# ----------------------------------------------------------------------------------------

from prody.proteins.pdbfile import parsePDB, fetchPDB
import time
import os
import numpy as np
from prody import LOGGER



__all__ = ['CharmmGUIBrowser']


class CharmmGUIBrowser(object):
    def __init__(self, cwd=None, ndir=None, fname=None, segids=None, job_type=None, timeout=None, saveas=None):
        if cwd is None:
            cwd = os.getcwd()
        self.cwd = cwd

        if ndir is None:
            ndir = ''
        self.ndir = ndir

        self.job_id = ''

        if job_type is None:
            job_type = 'solution'
        self.job_type = job_type

        if timeout is None:
            timeout = 15
        self.timeout = timeout

        if not isinstance(fname, str):
            raise TypeError('please provide a filename')

        if not fname.endswith('.pdb'):
            fname += '.pdb'

        self.fname = fname

        if not os.path.isfile("%s/%s/%s" % (cwd, ndir, fname)):
            fpath = fetchPDB(fname[:-4], compressed=False)
            if not os.path.isdir("%s/%s" % (cwd, ndir)):
                os.mkdir("%s/%s" % (cwd, ndir))
            os.rename(fpath, "%s/%s/%s" % (cwd, ndir, fname))

        self.ag = parsePDB("%s/%s/%s" % (cwd, ndir, fname),
                        subset='ca')

        if segids is None:
            segids = ['PRO'+ch.getChid() for ch in self.ag.iterChains()]
        self.segids = segids

        self.browser = None
        self.link = None
        self.status = False

        if saveas is None:
            saveas = 'charmm-gui'
        self.saveas = saveas

        self.run()

    def download(self, browser=None, link=None, saveas=None):
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
                LOGGER.info("Waiting for download link, sleep for %s s" % sleep)
                LOGGER.sleep(int(sleep), 'to check for download link.')
                sleep = int(sleep * 1.5)


    def check_error(self, browser=None, max_warn=0):
        if browser is None:
            browser = self.browser

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

    def wait_for_text(self, browser=None, text=None, timeout=None, max_warn=0):
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
                    status, err = self.check_error(browser, max_warn)
                    if status == True:
                        browser.execute_script("proceed()")                    

            if sleep > 180:
                sleep = 180

        if sleep_time >= timeout:
            raise TimeoutError('wait for %s step timed out' % text)

    def next_step(self, browser=None, text=None, timeout=None, max_warn=0):
        if browser is None:
            browser = self.browser

        if text is None:
            text = self.text

        browser.execute_script("proceed()")
        status = True

        LOGGER.info("    Going to next step")
        self.wait_for_text(browser, text, timeout)

        status, err = self.check_error(browser, max_warn)
        if status == False:
            raise Exception(err)

    def run(self, cwd=None, ndir=None, fname=None, segids=None, job_type=None, saveas=None, timeout=None):
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

        if saveas is None:
            saveas = self.saveas

        if timeout is None:
            timeout = self.timeout
        else:
            self.timeout = timeout

        if not job_type in ['solution', 'membrane.bilayer']:
            raise ValueError('job_type must be solution or membrane.bilayer')

        LOGGER.info("Running CHARMM-GUI Solution Builder in %s directory"
                    % ndir)

        try:
            from splinter import Browser
        except ImportError:
            raise ImportError('Browser module could not be imported. '
                              'install splinter package to solve the problem.')

        self.browser = browser = Browser('chrome')
        url = "http://www.charmm-gui.org/input/%s" % job_type
        browser.visit(url)

        browser.attach_file('file', "%s/%s/%s" % (cwd, ndir, fname))
        browser.find_by_value("PDB").click()

        self.next_step(browser, "Manipulate PDB")

        for item in browser.find_by_value("1"):
            item.uncheck()

        for segid in segids:
            browser.find_by_name("chains[%s][checked]" % segid.upper()).check()

        self.next_step(browser, "Generate PDB")

        if job_type == 'membrane.bilayer':
            self.next_step(browser, "Calculate Cross-Sectional Area")

            browser.find_by_name("align_option").first.click()
            self.next_step(browser, "Determine the System Size")

            maxes = np.max(self.ag.getCoords(), axis=0)
            minis = np.min(self.ag.getCoords(), axis=0)
            len_x = maxes[0] - minis[0]
            len_y = maxes[1] - minis[1]
            max_xy = np.ceil(np.max([len_x, len_y]))
            browser.find_by_name("hetero_lx").first.fill(str(max_xy+10))

            browser.find_by_id("lipid_ratio[upper][chl1]").first.fill("1")
            browser.find_by_id("lipid_ratio[lower][chl1]").first.fill("1")

            [y for y in browser.find_by_css("img") if y.outer_html.find("liptype_ratio_pc") != -1][0].click()
            time.sleep(5)
            browser.find_by_id("lipid_ratio[upper][popc]").first.fill("3")
            browser.find_by_id("lipid_ratio[lower][popc]").first.fill("3")

            browser.find_by_value("Show the system info").first.click()

            while len(browser.find_by_id("error_msg")[:]) > 0:
                if browser.find_by_id("error_msg")[0].text.find("more lipids") != -1:
                    browser.find_by_name("hetero_xy_option")[1].click()
                    time.sleep(5)
                    [y for y in browser.find_by_css("img") if y.outer_html.find("liptype_number_pc") != -1][0].click()
                    time.sleep(5)

                if browser.find_by_id("error_msg")[0].text.find("lower") != -1:
                    lower_chl1_num = browser.find_by_id("lipid_number[lower][popc]").first
                    lower_chl1_num.fill(str(int(lower_chl1_num.value)+1))
                elif browser.find_by_id("error_msg")[0].text.find("upper") != -1:
                    lower_chl1_num = browser.find_by_id("lipid_number[upper][popc]").first
                    lower_chl1_num.fill(str(int(lower_chl1_num.value)+1))
                else:
                    raise Exception(browser.find_by_id("error_msg")[0].text)

                browser.find_by_value("Show the system info")[1].click()

            self.next_step(browser, "Build Components", max_warn=1)
        else:
            self.next_step(browser, "Solvate Molecule")

        browser.find_option_by_text('NaCl').first.click()
        browser.find_option_by_text('Monte-Carlo').first.click()

        if job_type == 'membrane.bilayer':
            self.job_id = browser.find_by_name("jobid")[0].value
            self.next_step(browser, "Assemble Components")
        else:
            self.next_step(browser, "Setup Periodic Boundary Condition")

        self.next_step(browser, "Generate Equilibration and Dynamics", timeout=12000)

        browser.find_by_name("namd_checked").check()

        self.next_step(browser, "production.inp")

        link = self.get_download_link(browser)
        LOGGER.info("Build success")
        status = "Success"

        self.browser = browser
        self.link = link
        self.status = status
        
        self.download(saveas=saveas)
        self.browser.quit()
        return


if __name__ == "__main__":

    cgb = CharmmGUIBrowser(fname='1ake', saveas='1ake-charmm-gui.tgz', ndir='test')

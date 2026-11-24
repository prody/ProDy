# -*- coding: utf-8 -*-

"""RefSeq database operations."""

import os
import re
from concurrent.futures import ThreadPoolExecutor

import requests
import requests_cache
from datetime import timedelta
from prody import LOGGER
from prody.utilities.motifhelpers import downloadFile

__all__ = ["RefSeq"]

requests_cache.install_cache(expire_after=timedelta(weeks=1))


class RefSeq:
    """RefSeq database."""

    RELEASE_FILE = "release.txt"
    CHECK_SUMS = "files_installed.txt"

    @classmethod
    def getCurrentRelease(cls) -> str:
        """Get current RefSeq db release version.

        Raises:
            Exception: No release version found.

        Returns:
            str: RefSeq db release version.
        """
        rs_release = ""
        url = "https://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER"
        try:
            with requests_cache.disabled():
                response = requests.get(url)
        except requests.exceptions.RequestException as exception:
            LOGGER.error(str(exception))
        else:
            rs_release = response.text.strip()
        if not re.match(r"\d+", rs_release):
            LOGGER.error("Could't determine release version.")
        LOGGER.debug("RefSeq current release: {}".format(rs_release))
        return rs_release

    @classmethod
    def getInstalledFiles(cls) -> str:
        """Reads installed files with corresponding check sums."""
        LOGGER.debug("Downloading installed files with check sums.")
        current_release = cls.getCurrentRelease()
        url = "https://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/release{}.files.installed".format(current_release)
        try:
            with requests_cache.disabled():
                response = requests.get(url)
        except requests.exceptions.RequestException as exception:
            LOGGER.error(str(exception))
        else:
            result = ""
            lines = response.text.split("\n")
            for line in lines:
                if re.search(r"complete\.(\d+\.){1,2}protein\.faa\.gz", line):
                    result += line + "\n"
            return result

    @classmethod
    def saveInstalledFiles(cls, data_dir: str) -> None:
        """Saves installed files list with check sums locally."""
        installed_files = cls.getInstalledFiles()
        path = os.path.join(data_dir, cls.__name__)
        os.makedirs(path, exist_ok=True)
        with open("{}/{}".format(path, cls.CHECK_SUMS), "w", encoding="utf-8") as file:
            file.write(installed_files)
        LOGGER.debug("{} file saved locally.".format(cls.CHECK_SUMS))

    @classmethod
    def getLocalFiles(cls, data_dir: str) -> dict:
        """Lists local RefSeq FASTA protein files and corresponding check sums.

        Returns:
            dict: file names and corresponding check sums.
        """
        LOGGER.debug("Getting local RefSeq protein FASTA files")
        path = os.path.join(data_dir, cls.__name__)
        os.makedirs(path, exist_ok=True)
        filename = os.path.join(path, cls.CHECK_SUMS)
        if os.path.exists(filename):
            with open(filename, "r", encoding="utf-8") as file:
                text = file.read()
                results = re.findall(r"^(\d+)\s+(\S+)$", text, re.M)
                return {result[1]: result[0] for result in results}
        else:
            return {}

    @classmethod
    def getFiles(cls) -> dict:
        """Lists all FASTA protein files on RefSeq ftp server.

        Returns:
            dict: FASTA protein files with check sums.
        """
        LOGGER.debug("Getting protein FASTA file list from RefSeq.")
        installed_files = cls.getInstalledFiles()
        file_matcher = re.compile(r"^(\d+)\s+(\S+)$", re.M)
        files = re.findall(file_matcher, installed_files)
        return {file[1]: file[0] for file in files} if files else {}

    @classmethod
    def pepareDownloadFileList(cls, data_dir: str) -> list:
        """Prepare file list to be downloaded"""
        LOGGER.debug("Preparing file list to be downloaded.")
        remote_files = cls.getFiles()
        local_files = cls.getLocalFiles(data_dir)
        download_list = []
        for filename in remote_files.keys():
            if filename in local_files.keys():
                if remote_files[filename] != local_files[filename]:
                    download_list.append(filename)
            else:
                download_list.append(filename)
        return download_list

    @classmethod
    def downloadRelease(cls, data_dir: str) -> None:
        """Download latest RefSeq database release."""

        def report_done(future):
            LOGGER.info("Downloaded file {}".format(future.result()))

        files = cls.pepareDownloadFileList(data_dir)
        url = "https://ftp.ncbi.nlm.nih.gov/refseq/release/complete/"
        LOGGER.timeit("downloadRelease")
        with ThreadPoolExecutor(max_workers=100) as executor:
            for file in files:
                directory = os.path.join(data_dir, cls.__name__)
                future = executor.submit(downloadFile, url, directory, file)
                future.add_done_callback(report_done)
        LOGGER.report(msg="downloadRelease completed in %.2fs.", label="downloadRelease")

    @classmethod
    def saveRelease(cls, data_dir: str) -> None:
        """Write current release version to disk."""
        current_release = cls.getCurrentRelease()
        path = os.path.join(data_dir, cls.__name__)
        os.makedirs(path, exist_ok=True)
        with open("{}/{}".format(path, cls.RELEASE_FILE), "w", encoding="utf-8") as file:
            file.write(current_release)
        LOGGER.debug("RefSeq release {} saved.".format(current_release))

    @classmethod
    def updateRelease(cls, data_dir) -> None:
        """Update release to the most recent one."""
        LOGGER.debug("Updating RefSeq release version.")
        cls.downloadRelease(data_dir)
        cls.saveRelease(data_dir)

    @classmethod
    def getLocalRelease(cls, data_dir) -> str:
        """Get release version from local disk."""
        path = os.path.join(data_dir, cls.__name__)
        os.makedirs(path, exist_ok=True)
        version = ""
        release = os.path.join(path, cls.RELEASE_FILE)
        if os.path.exists(release):
            with open(release, "r", encoding="utf-8") as file:
                version = file.readline()
                LOGGER.debug("RefSeq local release: {}".format(version))
        else:
            LOGGER.debug("RefSeq local release not found.")
        return version

    @classmethod
    def checkForUpdates(cls, data_dir) -> None:
        """Check if local version is the recent one."""
        LOGGER.debug("RefSeq checking for updates.")
        if cls.getCurrentRelease() != cls.getLocalRelease(data_dir):
            cls.updateRelease(data_dir)

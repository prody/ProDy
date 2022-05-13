# -*- coding: utf-8 -*-

"""RefSeq database operations."""

import os
import re
from concurrent.futures import ThreadPoolExecutor

import requests
from prody import LOGGER
from prody.utilities.helpers import downloadFile
from prody.utilities.pathtools import PRODY_DATA

__all__ = ["RefSeq"]


class RefSeq:
    """RefSeq database."""

    RELEASE = "release.txt"
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
            response = requests.get(url)
        except requests.exceptions.RequestException as exception:
            LOGGER.error(str(exception))
        else:
            rs_release = response.text
        if not re.match(r"\d+", rs_release):
            LOGGER.error("Could't determine release version.")
        LOGGER.debug("RefSeq current release: {}".format(rs_release))
        return rs_release

    @classmethod
    def getInstalledFiles(cls) -> str:
        """Reads installed files with corresponding check sums."""
        LOGGER.debug("Downloading installed files with check sums.")
        current_release = cls.getCurrentRelease()
        url = f"https://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/release{current_release}.files.installed"
        try:
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
    def saveInstalledFiles(cls) -> None:
        """Saves installed files list with check sums locally."""
        installed_files = cls.getInstalledFiles()
        path = os.path.join(PRODY_DATA, cls.__name__)
        os.makedirs(path, exist_ok=True)
        with open(f"{PRODY_DATA}/{cls.__name__}/{cls.CHECK_SUMS}", "w", encoding="utf-8") as file:
            file.write(installed_files)
        LOGGER.debug(f"{cls.CHECK_SUMS} file saved locally.")

    @classmethod
    def getLocalFiles(cls) -> dict:
        """Lists local RefSeq FASTA protein files and corresponding check sums.

        Returns:
            dict: file names and corresponding check sums.
        """
        LOGGER.debug("Getting local RefSeq protein FASTA files")
        path = os.path.join(PRODY_DATA, cls.__name__)
        os.makedirs(path, exist_ok=True)
        with open(f"{PRODY_DATA}/{cls.__name__}/{cls.CHECK_SUMS}", "r", encoding="utf-8") as file:
            text = file.read()
            results = re.findall(r"^(\d+)\s+(\S+)$", text, re.M)
            return {result[1]: result[0] for result in results}

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
    def pepareDownloadFileList(cls) -> list:
        """Prepare file list to be downloaded"""
        LOGGER.debug("Preparing file list to be downloaded.")
        remote_files = cls.getFiles()
        local_files = cls.getLocalFiles()
        download_list = []
        for filename in remote_files.keys():
            if filename in local_files.keys():
                if remote_files[filename] != local_files[filename]:
                    download_list.append(filename)
            else:
                download_list.append(filename)
        return download_list

    @classmethod
    def downloadRelease(cls) -> None:
        """Download latest RefSeq database release."""
        files = cls.pepareDownloadFileList()
        url = ""
        LOGGER.timeit()
        with ThreadPoolExecutor(max_workers=3) as executor:
            for file in files:
                future = executor.submit(downloadFile, url, cls.__name__, file)
                LOGGER.info(future.result())
        LOGGER.report()

    @classmethod
    def saveRelease(cls) -> None:
        """Write current release version to disk."""
        current_release = cls.getCurrentRelease()
        path = os.path.join(PRODY_DATA, cls.__name__)
        os.makedirs(path, exist_ok=True)
        with open(f"{PRODY_DATA}/{cls.__name__}/{cls.RELEASE}", "w", encoding="utf-8") as file:
            file.write(current_release)
        LOGGER.debug("RefSeq release {} saved.".format(current_release))

    @classmethod
    def updateRelease(cls) -> None:
        """Update release to the most recent one."""
        LOGGER.debug("Updating RefSeq release version.")
        cls.downloadRelease()
        cls.saveRelease()

    @classmethod
    def getLocalRelease(cls) -> str:
        """Get release version from local disk."""
        LOGGER.debug("Getting RefSeq local release version.")
        path = os.path.join(PRODY_DATA, cls.__name__)
        os.makedirs(path, exist_ok=True)
        with open(f"{PRODY_DATA}/{cls.__name__}/{cls.RELEASE}", "r", encoding="utf-8") as file:
            return file.readline()

    @classmethod
    def checkForUpdates(cls) -> None:
        """Check if local version is the recent one."""
        LOGGER.debug("RefSeq checking for updates.")
        if cls.getCurrentRelease() != cls.getLocalRelease():
            cls.updateRelease()

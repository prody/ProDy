# -*- coding: utf-8 -*-

"""Swiss-Prot database operations."""

import os
import re
import requests
import requests_cache
from datetime import timedelta
from concurrent.futures import ThreadPoolExecutor

from prody import LOGGER
from prody.utilities.motifhelpers import downloadFile

__all__ = ["SwissProt"]

requests_cache.install_cache(expire_after=timedelta(weeks=1))


class SwissProt:
    """Swiss-Prot database."""

    RELEASE_FILE = "release.txt"
    HEADERS = {"User-Agent": "Python pattern search agent", "Contact": "mkonstanty@gmail.com"}

    @classmethod
    def getCurrentRelease(cls) -> str:
        """Get current swiss-prot db release version.

        Raises:
            Exception: No release version found.

        Returns:
            str: Swiss-prot db release version.
        """
        sp_release = ""
        url = "https://ftp.expasy.org/databases/swiss-prot/release/reldate.txt"
        try:
            with requests_cache.disabled():
                response = requests.get(url, headers=cls.HEADERS)
        except requests.exceptions.RequestException as exception:
            LOGGER.error(str(exception))
        else:
            sp_release = re.search(r"Swiss-Prot Release (\d{4}_\d{2})", response.text)
        if not sp_release:
            LOGGER.error("Could't determine release version.")
        LOGGER.debug("Swiss-Prot current release: {}".format(sp_release[1]))
        return sp_release[1]

    @classmethod
    def downloadRelease(cls, data_dir: str, types=None) -> None:
        """Download latest swiss-prot database release.

        Args:
            types (list, optional): Database file types. Defaults to None.
        """

        def report_done(future):
            LOGGER.info("Downloaded file {}".format(future.result()))

        types = types if types else ["xml", "dat", "fasta"]
        files = ["uniprot_sprot.{}.gz".format(type) for type in types]
        url = "https://ftp.expasy.org/databases/swiss-prot/release/"
        LOGGER.timeit("downloadRelease")
        with ThreadPoolExecutor(max_workers=3) as executor:
            for file in files:
                directory = os.path.join(data_dir, cls.__name__)
                future = executor.submit(downloadFile, url, directory, file, headers=cls.HEADERS)
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
        LOGGER.debug("Swiss-Prot release {} saved.".format(current_release))

    @classmethod
    def updateRelease(cls, data_dir) -> None:
        """Update release to the most recent one."""
        LOGGER.debug("Updating Swiss-Prot release version.")
        cls.downloadRelease(data_dir, types=["fasta"])
        cls.saveRelease(data_dir)

    @classmethod
    def getLocalRelease(cls, data_dir: str) -> str:
        """Get release version from local disk."""
        LOGGER.debug("Getting Swiss-Prot local release version.")
        path = os.path.join(data_dir, cls.__name__)
        os.makedirs(path, exist_ok=True)
        version = ""
        release = os.path.join(path, cls.RELEASE_FILE)
        if os.path.exists(release):
            with open(release, "r", encoding="utf-8") as file:
                version = file.readline()
                LOGGER.debug("Swiss-Prot local release: {}".format(version))
        else:
            LOGGER.debug("Swiss-Prot local release not found.")
        return version

    @classmethod
    def checkForUpdates(cls, data_dir) -> None:
        """Check if local version is the recent one."""
        LOGGER.debug("Swiss-Prot checking for updates.")
        if cls.getCurrentRelease() != cls.getLocalRelease(data_dir):
            cls.updateRelease(data_dir)

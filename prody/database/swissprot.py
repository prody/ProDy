# -*- coding: utf-8 -*-

"""Swiss-Prot database operations."""

import os
import re
from concurrent.futures import ThreadPoolExecutor

import requests
from prody import LOGGER
from prody.utilities.helpers import downloadFile
from prody.utilities.pathtools import PRODY_DATA

__all__ = ["SwissProt"]


class SwissProt:
    """Swiss-Prot database."""

    RELEASE = "release.txt"
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
    def downloadRelease(cls, types=None) -> None:
        """Download latest swiss-prot database release.

        Args:
            types (list, optional): Database file types. Defaults to None.
        """
        types = types if types else ["xml", "dat", "fasta"]
        files = [f"uniprot_sprot.{type}.gz" for type in types]
        url = "https://ftp.expasy.org/databases/swiss-prot/release/"
        LOGGER.timeit()
        with ThreadPoolExecutor(max_workers=3) as executor:
            for file in files:
                future = executor.submit(downloadFile, url, cls.__name__, file, headers=cls.HEADERS)
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
        LOGGER.debug("Swiss-Prot release {} saved.".format(current_release))

    @classmethod
    def updateRelease(cls) -> None:
        """Update release to the most recent one."""
        LOGGER.debug("Updating Swiss-Prot release version.")
        cls.downloadRelease()
        cls.saveRelease()

    @classmethod
    def getLocalRelease(cls) -> str:
        """Get release version from local disk."""
        LOGGER.debug("Getting Swiss-Prot local release version.")
        path = os.path.join(PRODY_DATA, cls.__name__)
        os.makedirs(path, exist_ok=True)
        with open(f"{PRODY_DATA}/{cls.__name__}/{cls.RELEASE}", "r", encoding="utf-8") as file:
            return file.readline()

    @classmethod
    def checkForUpdates(cls) -> None:
        """Check if local version is the recent one."""
        LOGGER.debug("Swiss-Prot checking for updates.")
        if cls.getCurrentRelease() != cls.getLocalRelease():
            cls.updateRelease()

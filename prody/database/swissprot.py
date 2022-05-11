# -*- coding: utf-8 -*-

"""Swiss-Prot database operations."""

import re

import requests
from prody import LOGGER

__all__ = ["SwissProt"]


class SwissProt:
    """Swiss-Prot database."""

    headers = {"User-Agent": "Python pattern search agent", "Contact": "mkonstanty@gmail.com"}

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
            response = requests.get(url, headers=cls.headers)
        except requests.exceptions.RequestException as e:
            LOGGER.error(str(e))
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
        for file in files:
            LOGGER.info("Downloading file {}.".format(file))
            try:
                with requests.get(url + file, headers=cls.headers, stream=True) as request:
                    with open(f"./data/{file}", "wb") as output:
                        for chunk in request.iter_content(chunk_size=16 * 1024):
                            output.write(chunk)
            except Exception as exception:
                LOGGER.error("Failed to download file {}: {}".format(file, exception))
        LOGGER.report()

    @classmethod
    def saveRelease(cls) -> None:
        """Write current release version to disk."""
        current_release = cls.getCurrentRelease()
        with open("./data/sp_release.txt", "w", encoding="utf-8") as file:
            file.write(current_release)
        LOGGER.debug("Swiss-Prot release {} saved.".format(current_release))

    @classmethod
    def updateRelease(cls) -> None:
        """Update release to the most recent one."""
        LOGGER.debug("Updating Swiss-Prot release version.")
        cls.downloadRelease()
        cls.saveRelease()

    @staticmethod
    def getLocalRelease() -> str:
        """Get release version from local disk."""
        LOGGER.debug("Getting Swiss-Prot local release version.")
        with open("./data/sp_release.txt", "r", encoding="utf-8") as file:
            return file.readline()

    @classmethod
    def checkForUpdates(cls) -> None:
        """Check if local version is the recent one."""
        LOGGER.debug("Swiss-Prot checking for updates.")
        if cls.getCurrentRelease() != cls.getLocalRelease():
            cls.updateRelease()

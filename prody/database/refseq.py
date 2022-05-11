# -*- coding: utf-8 -*-

"""RefSeq database operations."""

import re

import requests
from prody import LOGGER

__all__ = ["RefSeq"]


class RefSeq:
    """RefSeq database."""

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
        except requests.exceptions.RequestException as e:
            LOGGER.error(str(e))
        else:
            rs_release = response.text
        if not re.match(r"\d+", rs_release):
            LOGGER.error("Could't determine release version.")
        LOGGER.debug("RefSeq current release: {}".format(rs_release))
        return rs_release

    @classmethod
    def downloadRelease(cls) -> None:
        """Download latest RefSeq database release."""
        pass

    @classmethod
    def saveRelease(cls) -> None:
        """Write current release version to disk."""
        current_release = cls.getCurrentRelease()
        with open("./data/rs_release.txt", "w", encoding="utf-8") as file:
            file.write(current_release)
        LOGGER.debug("RefSeq release {} saved.".format(current_release))

    @classmethod
    def updateRelease(cls) -> None:
        """Update release to the most recent one."""
        LOGGER.debug("Updating RefSeq release version.")
        cls.downloadRelease()
        cls.saveRelease()

    @staticmethod
    def getLocalRelease() -> str:
        """Get release version from local disk."""
        LOGGER.debug("Getting RefSeq local release version.")
        with open("./data/rs_release.txt", "r", encoding="utf-8") as file:
            return file.readline()

    @classmethod
    def checkForUpdates(cls) -> None:
        """Check if local version is the recent one."""
        LOGGER.debug("RefSeq checking for updates.")
        if cls.getCurrentRelease() != cls.getLocalRelease():
            cls.updateRelease()

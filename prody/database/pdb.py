# -*- coding: utf-8 -*-

"""Protein Data Bank database operations."""

import os

from prody import LOGGER
from prody.utilities.motifhelpers import downloadFile

__all__ = ["ProteinDataBank"]


class ProteinDataBank:
    """Protein Data Bank"""

    @classmethod
    def downloadRelease(cls, data_dir: str) -> None:
        """Download latest pdb database.

        Args:
            data_dir (str): Directory to download to.
        """
        url = "https://ftp.wwpdb.org/pub/pdb/derived_data/"
        file = "pdb_seqres.txt.gz"
        LOGGER.timeit("downloadRelease")
        directory = os.path.join(data_dir, cls.__name__)
        downloadFile(url, directory, file)
        LOGGER.report(msg="downloadRelease completed in %.2fs.", label="downloadRelease")

# -*- coding: utf-8 -*-

"""This module defines functions and constants related to MOTIF."""

from typing import Any, Optional
import re
import os
import requests
import requests_cache
from datetime import timedelta
from .logger import LOGGER

__all__ = [
    "prositeToRegEx",
    "validateMotif",
    "downloadFile",
    "openFile",
    "getLocalDBs",
    "getGenericDBs",
]

CHUNK_SIZE = 16 * 1024

requests_cache.install_cache(expire_after=timedelta(weeks=1))


def prositeToRegEx(motif: str) -> str:
    """Change PROSITE Motif to python regular expression.

    Args:
        pattern (str): PROSITE Motif

    Returns:
        str: regular expression
    """
    pattern = (
        motif.replace("{", "[^")
        .replace("}", "]")
        .replace("(", "{")
        .replace(")", "}")
        .replace("-", "")
        .replace("x", ".")
        .replace(">", "$")
        .replace("<", "*")
    )
    return pattern


def validateMotif(prosite: str) -> bool:
    """_summary_

    Args:
        prosite (str): _description_

    Returns:
        bool: _description_
    """
    motif = prosite.replace("-", "")
    motif_pattern = r"^([ARDNCEQGHILKMFPSTWYVx]*(\[[ARDNCEQGHILKMFPSTWYV]+\])*({{[{{ARDNCEQGHILKMFPSTWYV}}]+}})*([ARDNCEQGHILKMFPSTWYVx]\(\d(,\d)?\))*)+$"  # noqa
    motif_matcher = re.compile(motif_pattern)
    matcher = re.match(motif_matcher, motif)
    return True if matcher else False


def downloadFile(url: str, dir: str, file: str, **kwargs: Any) -> str:
    """Downloads a file via http protocol.

    Args:
        url (str): file url for download
        dir (str): directory to download to
        file (str): filename to download
        kwargs (Any): keyword arguments for GET request
    """
    kwargs.update({"stream": True})
    LOGGER.info("Downloading file {}.".format(file))
    try:
        with requests_cache.disabled():
            with requests.get(url + file, **kwargs) as request:
                with open(os.path.join(dir, file), "wb") as output:
                    for chunk in request.iter_content(chunk_size=CHUNK_SIZE):
                        output.write(chunk)
    except requests.exceptions.RequestException as exception:
        LOGGER.error(str(exception))
    return file


def _getCompressionType(filename: str) -> Optional[str]:
    """Attempts to guess the compression (if any) on a file using the first few bytes.

    Args:
        filename (str): file to analize compression type

    Returns:
        Optional[str]: compression type
    """
    magic_dict = {
        "gz": (b"\x1f", b"\x8b", b"\x08"),
        "bz2": (b"\x42", b"\x5a", b"\x68"),
        "zip": (b"\x50", b"\x4b", b"\x03", b"\x04"),
    }
    max_len = max(len(x) for x in magic_dict)
    with open(filename, "rb") as unknown_file:
        file_start = unknown_file.read(max_len)
    compression_type = None
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    return compression_type


def openFile(filename: str) -> Any:
    """Open file and return handler.

    Args:
        filename (str): file to open

    Raises:
        NotImplementedError: compression not supported

    Returns:
        Any: file handler
    """
    compression = _getCompressionType(filename)
    if compression == "gz":
        import gzip

        return gzip.open(filename, "rt")
    elif not compression:
        return open(filename, "r")
    else:
        raise NotImplementedError


def getLocalDBs(upload_dir: str) -> list:
    """Get all local custom databases.

    Args:
        upload_dir (str): Directory to search

    Returns:
        list: databases
    """
    return os.listdir(upload_dir)


def getGenericDBs(data_dir: str) -> dict:
    """Get all buildin databases.

    Args:
        data_dir (str): Directory to search

    Returns:
        dict: databases
    """
    result = {}
    sp = f"{data_dir}/SwissProt/release.txt"
    if os.path.isfile(sp):
        result.update({"sp-local": "Swiss-Prot (local, buildin)"})
    rs = f"{data_dir}/RefSeq/release.txt"
    if os.path.isfile(rs):
        result.update({"rs-local": "RefSeq (local, buildin)"})
    pdb = f"{data_dir}/PDB/release.txt"
    if os.path.isfile(pdb):
        result.update({"pdb-local": "Protein Data Bank (local, buildin)"})
    return result

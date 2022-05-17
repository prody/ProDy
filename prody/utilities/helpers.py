# -*- coding: utf-8 -*-
"""This module defines functions and constants related to MOTIF."""

from typing import Any

import requests
from .logger import LOGGER
from prody.utilities.pathtools import PRODY_DATA

__all__ = ["prositeToRegEx", "downloadFile"]

CHUNK_SIZE = 16 * 1024


def prositeToRegEx(pattern: str) -> str:
    """Change PROSITE Motif to python regular expression.

    Args:
        pattern (str): PROSITE Motif

    Returns:
        str: regular expression
    """
    pattern = (
        pattern.replace("{", "[^")
        .replace("}", "]")
        .replace("(", "{")
        .replace(")", "}")
        .replace("-", "")
        .replace("x", ".")
        .replace(">", "$")
        .replace("<", "*")
    )
    return pattern


def downloadFile(url: str, dir: str, file: str, **kwargs: Any) -> None:
    """Downloads a file via http protocol.

    Args:
        url (str): website address
        file (str): filename to download
        kwargs (Any): keyword arguments for GET request
    """
    kwargs.update({"stream": True})
    LOGGER.info("Downloading file {}.".format(file))
    try:
        with requests.get(url + file, **kwargs) as request:
            with open(f"{PRODY_DATA}/{dir}/{file}", "wb") as output:
                for chunk in request.iter_content(chunk_size=CHUNK_SIZE):
                    output.write(chunk)
    except requests.exceptions.RequestException as exception:
        LOGGER.error(str(exception))

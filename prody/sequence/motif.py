#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""This module gets PDB code from MOTIF."""

import argparse
import re

import requests
from prody import LOGGER

__author__ = "Mariusz Konstanty"

__all__ = ["getPdbFromMotif", "expasySearchMotif", "_expasyTextToStructure"]

DATABASES = {
    "sp": "Swiss-Prot",
    # "rs": "RefSeq",
    "pdb": "PDB",
    "tr": "TrEMBL",
    "local": "local database",
}
MOTIF_PATTERN = (
    r"^([ARDNCEQGHILKMFPSTWYVx]*"
    r"(\[[ARDNCEQGHILKMFPSTWYV]+\])*"
    r"({{[{{ARDNCEQGHILKMFPSTWYV}}]+}})*"
    r"([ARDNCEQGHILKMFPSTWYVx]\(\d(,\d)?\))*)+$"
)
MOTIF_MATCHER = re.compile(MOTIF_PATTERN)


def getPdbFromMotif(motif: str, database: str) -> str:
    if database not in DATABASES:
        raise ValueError(f"Database must be one of: {DATABASES.keys()}.")
    if not re.match(MOTIF_MATCHER, motif.replace("-", "")):
        raise ValueError(f"{motif} is not valid PROSITE motif.")
    result = expasySearchMotif(motif, database)
    pdb_code = ""
    return pdb_code


def _expasyTextToStructure(response: str) -> list:
    RESULT_PATTERN = (
        r">(sp|tr|pdb)\|"
        r"(\w+)\|(\w+).*\n"
        r"(.*)\n"
        r"([ARDNCEQGHILKMFPSTWYVx\n]+)\W+"
        r"(\d{,3}) - (\d{,3}):\W+(\w+)"
    )
    motifs = re.compile(RESULT_PATTERN, re.M)
    results = re.findall(motifs, response)
    structure = []
    names = ("database", "accession", "id", "description", "sequence", "start", "end", "match")
    for result in results:
        result = map(lambda x: x.replace("\n", ""), result)
        structure.append(dict(zip(names, result)))
    return structure


def expasySearchMotif(motif, database):
    """Search motif in remote Swiss-Prot database.

    https://prosite.expasy.org/scanprosite/scanprosite_doc.html#rest_intro

    Args:
        motif (str): Motif in PROSITE format.
        database (str): 'sp' (UniProtKB/Swiss-Prot) or 'tr' (UniProtKB/TrEMBL) or 'pdb' (PDB).

    Returns:
        str: Expasy response from motif search.
    """
    headers = {"User-Agent": "Python pattern search agent", "Contact": "mkonstanty@gmail.com"}
    url = "https://prosite.expasy.org/cgi-bin/prosite/PSScan.cgi"
    payload = {"sig": motif, "db": database}
    try:
        result = requests.get(url, params=payload, headers=headers)
    except Exception as exception:
        print(f"Remote search for motif {motif} in Swiss-Prot database failed: {exception}")
        return []
    else:
        return _expasyTextToStructure(result.text)


def _argMotif(motif: str) -> str:
    """Create motif type for argparse.ArgumentParser.

    Args:
        value (str): motif to be validated

    Returns:
        str: valid Motif
    """
    new_motif = str(motif).replace("-", "")
    if not re.match(MOTIF_MATCHER, new_motif):
        msg = f"{motif} is not a valid PROSITE MOTIF."
        raise argparse.ArgumentTypeError(msg)
    return new_motif


def _parseArgs():
    """Parse arguments from the command line.

    Returns:
        argparse.Namespace: parsed script arguments
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-m",
        "--motif",
        type=_argMotif,
        required=True,
        help="motif to search in a databases",
    )
    parser.add_argument(
        "-d",
        "--databases",
        choices=DATABASES,
        action="append",
        required=True,
        help="choose the databases to search from",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = _parseArgs()

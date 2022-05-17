# -*- coding: utf-8 -*-
"""This module gets PDB code from MOTIF."""

import os
import re
import xml.dom.minidom

import requests
from prody import LOGGER
from prody.utilities.pathtools import PRODY_DATA

__all__ = [
    "getPdbCodesFromMotif",
    "getUniprotId",
    "pdbSearch",
    "saveSearchResults",
    "expasySearchMotif",
    "_expasyMotifToStructure",
]

DATABASES = {
    "sp": "Swiss-Prot",
    "rs": "RefSeq",
    "pdb": "PDB",
    "local": "local database",
}
MOTIF_PATTERN = (
    r"^([ARDNCEQGHILKMFPSTWYVx]*"
    r"(\[[ARDNCEQGHILKMFPSTWYV]+\])*"
    r"({{[{{ARDNCEQGHILKMFPSTWYV}}]+}})*"
    r"([ARDNCEQGHILKMFPSTWYVx]\(\d(,\d)?\))*)+$"
)
MOTIF_MATCHER = re.compile(MOTIF_PATTERN)


def getPdbCodesFromMotif(motif: str) -> list:
    """Get PDB code from MOTIF.

    Args:
        motif (str): PROSITE motif

    Returns:
        list: List of results
    """
    LOGGER.debug("Getting PDB codes from MOTIF.")
    url = "http://www.rcsb.org/pdb/rest/search"
    headers = {"Content-Type": "application/x-www-form-urlencoded"}
    data = f"""
    <xml version="1.0" encoding="UTF-8">
        <orgPdbQuery version="1.0">
            <queryType>org.pdb.query.simple.MotifQuery</queryType>
            <description>Motif Query For: {motif}</description>
            <motif>{motif}</motif>
        </orgPdbQuery>
    </xml>
    """
    try:
        response = requests.post(url, headers=headers, data=data)
    except requests.exceptions.RequestException as exception:
        LOGGER.debug(str(exception))
    else:
        response.encoding = "UTF-8"
        return response.text.split("\n")


def getUniprotId(pdbId, pdbEntityNr) -> str:
    """Get UniProt ID from PDB code.

    Args:
        pdbId (str): PDB code
        pdbEntityNr (str): PDB entity number

    Returns:
        str: comma-separated accession ids
    """
    LOGGER.debug("Getting UniProt ID from PDB code.")
    url = f"http://www.rcsb.org/pdb/rest/describeMol?structureId={pdbId}"
    accession_ids = ""
    try:
        response = requests.get(url)
    except requests.exceptions.RequestException as exception:
        LOGGER.debug(str(exception))
    else:
        response.encoding = "UTF-8"
        mol_description = xml.dom.minidom.parseString(response.text)
        polymers = mol_description.getElementsByTagName("polymer")
        for polymer in polymers:
            if polymer.getAttribute("entityNr") == pdbEntityNr:
                accessions = polymer.getElementsByTagName("accession")
                for counter, accession in enumerate(accessions):
                    accession_ids += accession.getAttribute("id")
                    if len(accessions) > counter + 1:
                        accession_ids += ","
    return accession_ids


def pdbSearch() -> list:
    """Search for a motif in PDB database.

    Returns:
        list: list of proteins
    """
    proteins = []
    codes = getPdbCodesFromMotif()
    for code in codes:
        if code:
            code_and_nr = code.split(":")
            pdbEntityNr = code_and_nr.pop()
            pdbId = code_and_nr.pop()
            uniprotIds = getUniprotId(pdbId, pdbEntityNr)
            protein = {
                "pdbEntityNr": pdbEntityNr,
                "pdbId": pdbId,
                "uniprotId": uniprotIds,
            }
            proteins.append(protein)
    return proteins


def saveSearchResults(proteins: list) -> None:
    """Save protein information to the file.

    Args:
        proteins (list): List of proteins to be saved.
    """
    output_dir = "output"
    path = os.path.join(PRODY_DATA, "output")
    os.makedirs(path, exist_ok=True)
    with open(f"{PRODY_DATA}/{output_dir}/proteins.csv", "w", encoding="utf-8") as file:
        file.write(";".join(proteins[0].keys()) + "\n")
        for protein in proteins:
            file.write(";".join(protein.values()) + "\n")


def _expasyMotifToStructure(response: str) -> list:
    """Parse HTML response from Expasy.

    Args:
        response (str): HTML response from Expasy

    Returns:
        list: List of result dicts
    """
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
    keys = ("database", "accession", "id", "description", "sequence", "start", "end", "match")
    for result in results:
        result = map(lambda x: x.replace("\n", ""), result)
        structure.append(dict(zip(keys, result)))
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
        return _expasyMotifToStructure(result.text)

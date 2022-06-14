# -*- coding: utf-8 -*-

"""This module prase protein sequences using PROSITE motifs."""

import os
import re
import requests
import requests_cache
from datetime import timedelta
import json
from typing import Optional
from multiprocessing import Pool
from prody import LOGGER
from prody.utilities.motifhelpers import prositeToRegEx, openFile


__all__ = ["motifSearch", "getPdbCodesFromMotif", "saveMotifResults", "getUniprot"]

MOTIF_PATTERN = (
    r"^([ARDNCEQGHILKMFPSTWYVx]*"
    r"(\[[ARDNCEQGHILKMFPSTWYV]+\])*"
    r"({{[{{ARDNCEQGHILKMFPSTWYV}}]+}})*"
    r"([ARDNCEQGHILKMFPSTWYVx]\(\d(,\d)?\))*)+$"
)
MOTIF_MATCHER = re.compile(MOTIF_PATTERN)

requests_cache.install_cache(expire_after=timedelta(weeks=1))


def getUniprot(uniprot: str) -> str:
    headers = {"User-Agent": "Python pattern search agent", "Contact": "mkonstanty@gmail.com"}
    url = "https://www.uniprot.org/"
    payload = {"file": uniprot, "format": "txt"}
    response = requests.post(url, headers=headers, data=json.dumps(payload))
    return response.text


def getPdbCodesFromMotif(motif: str, results: Optional[int] = None, page: Optional[int] = None) -> list:
    """Get PDB code from motif.

    Args:
        motif (str): FASTA motif
        results (Optional[int], optional): How many results to return. Defaults to None.
        page (Optional[int], optional): Which result page to return. Defaults to None.

    Returns:
        list: PDB ids
    """
    LOGGER.debug("Get PDB code from motif {}".format(motif))
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    headers = {"Content-Type": "application/json;charset=utf-8"}
    if results and page:
        paginate = {"paginate": {"start": results * (page - 1), "rows": results * page - 1}}
    else:
        paginate = {"return_all_hits": True}
    data = {
        "return_type": "entry",
        "query": {
            "type": "terminal",
            "service": "seqmotif",
            "parameters": {"value": motif, "pattern_type": "prosite", "target": "pdb_protein_sequence"},
        },
        "request_options": paginate,
    }
    try:
        LOGGER.timeit("getPdbCodesFromMotif")
        response = requests.post(url, headers=headers, data=json.dumps(data))
        LOGGER.debug(f"Result taken from cache: {response.from_cache}")
        LOGGER.report(msg="getPdbCodesFromMotif completed in %.2fs.", label="getPdbCodesFromMotif")
    except requests.exceptions.RequestException as exception:
        LOGGER.error(str(exception))
    else:
        response.encoding = "UTF-8"
        result = json.loads(response.text)
        result["result_set"]
        return [x["identifier"] for x in result["result_set"]]


def saveMotifResults(motif: str, data: list[dict], directory: str) -> str:
    """Save search data to a file.

    Args:
        motif (str): FASTA motif
        data (list[dict]): List of proteins to be saved.
        directory (str): directory to save results

    Returns:
        str: saved file name
    """
    os.makedirs(directory, exist_ok=True)
    filename = f"{motif.replace('-', '')}.txt"
    with open(f"{directory}/{filename}", "w", encoding="utf-8") as file:
        for protein in data:
            for key, value in protein.items():
                if key == "sequence":
                    value = value.replace("\n", "")
                file.write(f"{key}: {value}\n")
            file.write("\n")
    LOGGER.info(f"Search written to a file {filename}.")
    return filename


def _expasyMotifToStructure(response: str) -> list:
    """Parse HTML response from Expasy.

    Args:
        response (str): HTML response from Expasy

    Returns:
        list: List of result dicts
    """
    LOGGER.timeit("_expasyMotifToStructure")
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
        structure.append(dict(zip(keys, result)))
    LOGGER.report(msg="_expasyMotifToStructure completed in %.2fs.", label="_expasyMotifToStructure")
    return structure


def expasySearchMotif(motif: str, database: str) -> list:
    """Search motif in remote Swiss-Prot database.

    https://prosite.expasy.org/scanprosite/scanprosite_doc.html#rest_intro

    Args:
        motif (str): Motif in PROSITE format.
        database (str): 'sp' (UniProtKB/Swiss-Prot) or 'tr' (UniProtKB/TrEMBL) or 'pdb' (PDB).

    Returns:
        list: List of result dicts
    """
    LOGGER.debug("Searching motif {} in {} database via Expasy".format(motif, database))
    headers = {"User-Agent": "Python pattern search agent", "Contact": "mkonstanty@gmail.com"}
    url = "https://prosite.expasy.org/cgi-bin/prosite/PSScan.cgi"
    payload = {"sig": motif, "db": database}
    try:
        LOGGER.timeit("expasySearchMotif")
        result = requests.get(url, params=payload, headers=headers)
        LOGGER.debug(f"Result taken from cache: {result.from_cache}")
        LOGGER.report(msg="expasySearchMotif completed in %.2fs.", label="expasySearchMotif")
    except Exception as exception:
        LOGGER.error(f"Remote search for motif {motif} in Swiss-Prot database failed: {exception}")
    else:
        return _expasyMotifToStructure(result.text)


def spOnlineMotifSearch(motif: str, *args) -> list:
    """Search motif in remote Swiss-Prot database.

    Args:
        motif (str): Motif in PROSITE format.

    Returns:
        list: List of result dicts
    """
    LOGGER.debug("Search for motif in Swiss-Prot online.")
    return expasySearchMotif(motif, "sp")


def rsOnlineMotifSearch(motif: str, *args) -> list:
    """ "Search motif in remote RefSeq database.

    Args:
        motif (str): Motif in PROSITE format.

    Raises:
        NotImplementedError: not implemented yet

    Returns:
        list: List of result dicts
    """
    LOGGER.debug("Search for motif in RefSeq online.")
    # return remoteSearchMotif(motif, "refseq")
    raise NotImplementedError


def pdbOnlineMotifSearch(motif: str, *args) -> list:
    """Search motif in remote PDB database.


    Args:
        motif (str): Motif in PROSITE format.

    Returns:
        list: list: List of result dicts
    """
    LOGGER.debug("Search for motif in Protein Data Bank online.")
    return expasySearchMotif(motif, "pdb")


def _parseProtein(pattern: str, protein: list) -> list:
    """Parse protein sequence.

    Args:
        pattern (str): motif regex
        protein (list): FASTA protein

    Returns:
        Optional[list]: Parsed protein.
    """
    line = protein[0].replace("\n", "").replace(">", "")
    sequence = "".join(protein[1:]).replace("\n", "")
    tmp = line.split(" ")
    database, accession, id = tmp[0].split("|")
    description = " ".join(tmp[1:])
    motif_matcher = re.compile(pattern)
    results = re.finditer(motif_matcher, sequence)
    elements = []
    for result in results:
        elements.append(
            {
                "database": database,
                "accession": accession,
                "id": id,
                "description": description,
                "sequence": sequence,
                "start": result.start(),
                "end": result.end(),
                "match": result.group(0),
            }
        )
    return elements


def localMotifSearch(filename: str, motif: str) -> list:
    """Search for motif in local database file.

    Args:
        filename (str): path to filename
        motif (str): PROSITE motif

    Returns:
        list: Parsed proteins
    """
    data = []
    LOGGER.timeit("localMotifSearch")
    pattern = prositeToRegEx(motif)
    protein = []
    try:
        file = openFile(filename)
        line = file.readline()
        if line.startswith(">"):
            protein.append(line)
        else:
            LOGGER.error("Only FASTA files supported.")
        # pool = Pool(processes=os.cpu_count())
        # results = []
        for line in file:
            if line.startswith(">"):
                # results.append(pool.apply_async(_parseProtein, (pattern, protein,)))
                data.extend(_parseProtein(pattern, protein))
                protein.clear()
                protein.append(line)
            else:
                protein.append(line)
        # data = [result.get() for result in results]
        # pool.close()
        # pool.join()
    except Exception as e:
        LOGGER.error(str(e))
    finally:
        file.close()
    LOGGER.report(msg="localMotifSearch completed in %.2fs.", label="localMotifSearch")
    return data


def spLocalMotifSearch(motif: str, data_dir: str, *args) -> list:
    """Search for motif in local Swiss-Prot database file.

    Args:
        motif (str): PROSITE motif
        data_dir (str): local Swiss-Prot database directory
    Returns:
        list: Parsed proteins
    """
    LOGGER.debug("Running local Swiss-Prot motif search...")
    filename = os.path.join(data_dir, "SwissProt/uniprot_sprot.fasta.gz")
    return localMotifSearch(filename, motif)


def rsLocalMotifSearch(motif: str, data_dir: str, upload_dir: str) -> list:
    """Search for motif in local RefSeq database file.

    Args:
        motif (str): PROSITE motif
        data_dir (str): local RefSeq database directory
        upload_dir (str): directory for local custom databases

    Raises:
        NotImplementedError: not implemented yet

    Returns:
        list: list: Parsed proteins
    """
    raise NotImplementedError


def pdbLocalMotifSearch(motif: str, data_dir: str, upload_dir: str) -> list:
    """Search for motif in local PDB database file.

    Args:
        motif (str): PROSITE motif
        data_dir (str): local RefSeq database directory
        upload_dir (str): directory for local custom databases

    Raises:
        NotImplementedError: not implemented yet

    Returns:
        list: Parsed proteins
    """
    raise NotImplementedError


def motifSearch(database: str, motif: str, upload_dir: str, data_dir: str) -> Optional[list]:
    """Search for motif in protein database.

    Args:
        database (str): selected database
        motif (str): motif to be searched
        data_dir (str): local builin database directory
        upload_dir (str): local custom database directory

    Returns:
        Optional[list]: Parsed proteins
    """

    def custom(motif: str, data_dir: str, upload_dir: str) -> Optional[list]:
        LOGGER.debug("Running local motif search...")
        filename = os.path.join(upload_dir, database)
        return localMotifSearch(filename, motif)

    LOGGER.debug(f"Motif search: {database=}\t{motif=}")
    return {
        "sp-online": spOnlineMotifSearch,
        "rs-online": rsOnlineMotifSearch,
        "pdb-online": pdbOnlineMotifSearch,
        "sp-local": spLocalMotifSearch,
        "rs-local": rsLocalMotifSearch,
        "pdb-local": pdbLocalMotifSearch,
    }.get(database, custom)(motif, data_dir, upload_dir)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Residue-level evolutionary conservation from multiple sequence alignments.

This module provides a convenience function to compute wild-type and
type-based conservation metrics for each residue of a target protein
sequence using a multiple sequence alignment (MSA) generated from BLAST
XML results or provided directly as a pre-aligned MSA file.

The workflow:

1. Reads the target (wild-type) amino-acid sequence from a text file
   (single line, one-letter codes).
2. Optionally parses a single-chain PDB file whose sequence matches
   the target and extracts per-residue identifiers (chain, residue
   number, insertion code) to enable mapping conservation scores into
   PDB B-factors.
3. Obtains an MSA either by:
   - Parsing BLAST XML output and writing a padded `.seqmsa` file
     with sequences filtered by coverage and e-value, or
   - Loading a pre-computed `.seqmsa` file (one aligned sequence
     per line) whose length matches the target sequence.
4. For each residue position, computes:
   - Wild-type–based scores such as:
     * Identity fraction
     * Shannon entropy and normalized entropy
     * Consensus residue and its frequency
     * Gap fraction
     * BLOSUM62-based scores (WT vs others and pairwise)
   - Type-based scores using a standard 4-group classification
     (hydrophobic, polar, negative, positive).
5. Writes a tab-delimited text file containing one row per residue
   with sequence index, optional PDB identifiers, consensus residue
   and type, and all selected wild-type and type metrics.
6. Optionally:
   - Generates one or more PNG heatmaps visualizing selected metrics
     in a ConSurf-style layout.
   - Writes out one or two PDB files in which per-residue conservation
     scores are stored in the B-factor column (occupancy set to 1.00).

The main entry point is:

* :func:`computeConservationFromMSA` – run the full workflow for a
  given sequence and MSA (from BLAST XML or pre-computed file) and
  write text, PNG, and PDB outputs as requested.

The same module can also be used as a standalone script from the
command line; see the CLI options in :func:`main`.
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from collections import Counter
from typing import Iterable, List, Optional, Sequence, Tuple, Dict

from Bio.Blast import NCBIXML
from Bio.PDB import PDBParser
from Bio.Align import substitution_matrices

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np

__all__ = ["computeConservationFromMSA"]

BLOSUM62 = substitution_matrices.load("BLOSUM62")

# ------------------------- Amino-acid helpers -------------------------

THREE_TO_ONE = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLU": "E",
    "GLN": "Q",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "SEC": "U",
    "PYL": "O",
    "UNK": "X",
}

# Standard 4-group classification for type conservation
# H: Hydrophobic, P: Polar Neutral, N: Negative/Acidic, B: Positive/Basic, X: Special/Unknown/Gap
AA_TYPES = {
    "A": "H",
    "I": "H",
    "L": "H",
    "M": "H",
    "F": "H",
    "P": "H",
    "W": "H",
    "V": "H",
    "G": "H",  # Hydrophobic (H)
    "S": "P",
    "T": "P",
    "C": "P",
    "N": "P",
    "Q": "P",
    "Y": "P",  # Polar Neutral (P)
    "D": "N",
    "E": "N",  # Negative/Acidic (N)
    "K": "B",
    "R": "B",
    "H": "B",  # Positive/Basic (B)
    "U": "X",
    "O": "X",
    "X": "X",
    "-": "X",  # Special/Unknown/Gap (X)
}

WT_HEADER_MAP = {
    "identity": "WT_frac_id",
    "entropy": "WT_entropy_bits",
    "norm_entropy": "WT_conservation",
    "consensus_freq": "WT_frac_consensus",
    "gap_fraction": "WT_frac_gaps",
    "blosum_wt": "WT_blosum_to_WT",
    "blosum_pairwise": "WT_blosum_pairwise",
}

TYPE_HEADER_MAP = {
    "type_freq": "TYPE_frac_class",
    "type_entropy": "TYPE_class_entropy",
    "type_consensus_freq": "TYPE_frac_consensus_cls",
}

WT_TITLE_MAP = {
    "identity": "WT identity fraction",
    "entropy": "WT entropy (bits)",
    "norm_entropy": "WT conservation (1 − norm. entropy)",
    "consensus_freq": "WT consensus residue frequency",
    "gap_fraction": "WT gap fraction",
    "blosum_wt": "WT BLOSUM62 to WT",
    "blosum_pairwise": "WT pairwise BLOSUM62",
}

TYPE_TITLE_MAP = {
    "type_freq": "Class conservation (WT class frequency)",
    "type_entropy": "Class entropy (H/P/N/B)",
    "type_consensus_freq": "Consensus class frequency",
}


# ------------------------- Core helpers -------------------------

def read_target_sequence(seq_file: str) -> str:
    """
    Read the target (wild-type) sequence from a text file.

    The file is expected to contain a single line of one-letter amino-acid
    codes; leading/trailing whitespace is stripped.

    Parameters
    ----------
    seq_file : str
        Path to the sequence file.

    Returns
    -------
    str
        The target sequence as a single string.

    Raises
    ------
    FileNotFoundError
        If the sequence file cannot be found.
    ValueError
        If the file is empty.
    """
    with open(seq_file, "r") as f:
        sequence = f.read().strip()
    if not sequence:
        raise ValueError("Sequence file is empty.")
    return sequence


def read_pdb_residue_ids(
    pdb_file: str, target_sequence: str
) -> List[Tuple[str, int, str]]:
    """
    Parse a PDB file, verify the sequence against the target, and extract
    per-residue identifiers (chain, residue number, insertion code).

    Parameters
    ----------
    pdb_file : str
        Path to a single-chain PDB file.
    target_sequence : str
        Target sequence to validate against the PDB.

    Returns
    -------
    list of (str, int, str)
        List of (chain_id, resseq, insertion_code) for each residue,
        in the same order as the target_sequence.

    Raises
    ------
    ValueError
        If the PDB has no chains, the sequence does not match, or the
        lengths differ.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    model = structure[0]

    try:
        chain = next(model.get_chains())
    except StopIteration:
        raise ValueError("PDB file does not contain any chains.")

    pdb_residue_ids: List[Tuple[str, int, str]] = []
    seq_index = 0
    for residue in chain.get_residues():
        res_name = residue.get_resname()
        if res_name not in THREE_TO_ONE:
            continue
        one_letter = THREE_TO_ONE[res_name]
        if seq_index >= len(target_sequence) or one_letter != target_sequence[seq_index]:
            res_num_str = f"{residue.get_id()[1]}{residue.get_id()[2].strip()}"
            raise ValueError(
                f"Sequence mismatch at position {seq_index + 1}: "
                f"PDB residue {chain.id}:{res_num_str} is '{one_letter}', "
                f"but target residue is '{target_sequence[seq_index]}'."
            )
        res_id_tuple = (chain.id, residue.get_id()[1], residue.get_id()[2].strip())
        pdb_residue_ids.append(res_id_tuple)
        seq_index += 1

    if seq_index != len(target_sequence):
        raise ValueError(
            f"PDB sequence length ({seq_index}) does not match target "
            f"sequence length ({len(target_sequence)})."
        )

    return pdb_residue_ids


def generate_msa_from_blast_xml(
    blast_file: str, seq: str, msa_out_path: str
) -> int:
    """
    Generate an MSA from BLAST XML results.

    Sequences are filtered by e-value (<= 0.001) and query coverage
    (>= 30 %). Subject sequences are padded with '-' characters to
    match the exact length of the query sequence.

    Parameters
    ----------
    blast_file : str
        Path to BLAST XML file.
    seq : str
        Target (query) sequence.
    msa_out_path : str
        Path to the output `.seqmsa` file.

    Returns
    -------
    int
        Number of sequences written to the MSA.

    Raises
    ------
    FileNotFoundError
        If the BLAST XML file is missing.
    ValueError
        If the target sequence is empty.
    RuntimeError
        For general parsing errors.
    """
    seq_length = len(seq)
    if seq_length == 0:
        raise ValueError("Target query sequence is empty; cannot parse BLAST results.")

    hits_processed = 0
    try:
        with open(blast_file, "r") as blast_handle:
            blast_records = NCBIXML.parse(blast_handle)
            with open(msa_out_path, "w") as fp:
                for blast_record in blast_records:
                    for ali in blast_record.alignments:
                        for hsp in ali.hsps:
                            query_len_aligned = hsp.query_end - hsp.query_start + 1
                            query_cover = query_len_aligned / seq_length * 100.0
                            if query_cover < 30.0 or hsp.expect > 0.001:
                                continue
                            pre_pad = "-" * (hsp.query_start - 1)
                            aligned_region = hsp.sbjct
                            post_pad = "-" * (seq_length - hsp.query_end)
                            aligned_seq = (pre_pad + aligned_region + post_pad)
                            aligned_seq = aligned_seq[:seq_length].ljust(seq_length, "-")
                            fp.write(aligned_seq + "\n")
                            hits_processed += 1
    except FileNotFoundError:
        raise
    except Exception as e:
        raise RuntimeError(f"Error during BLAST XML parsing: {e}") from e

    return hits_processed


def load_msa(msa_file_path: str) -> List[str]:
    """
    Load sequences from an MSA file.

    Parameters
    ----------
    msa_file_path : str
        Path to an MSA file (one sequence per line).

    Returns
    -------
    list of str
        List of aligned sequences.
    """
    with open(msa_file_path, "r") as f:
        return [line.strip() for line in f if line.strip()]


def calculate_conservation_scores(
    msa_sequences: Sequence[str],
    target_sequence: str,
    total_hits: int,
    wt_scores: Iterable[str],
    type_scores: Iterable[str],
) -> Tuple[List[Dict[str, float]], List[Dict[str, float]], str, str]:
    """
    Calculate multiple wild-type–based and type-based conservation scores
    per residue.

    Parameters
    ----------
    msa_sequences : sequence of str
        Aligned sequences (all same length as target).
    target_sequence : str
        Target (wild-type) sequence.
    total_hits : int
        Number of sequences in the MSA.
    wt_scores : iterable of str
        Names of wild-type metrics to compute.
    type_scores : iterable of str
        Names of type-based metrics to compute.

    Returns
    -------
    wt_metrics : list of dict
        One dict per residue, containing wild-type metrics.
    type_metrics : list of dict
        One dict per residue, containing type-based metrics.
    consensus_seq_res : str
        Consensus residue sequence.
    consensus_seq_type : str
        Consensus type sequence.
    """
    if not msa_sequences:
        return [], [], "", ""

    wt_scores = set(wt_scores)
    type_scores = set(type_scores)

    seq_len = len(target_sequence)
    wt_metrics: List[Dict[str, float]] = []
    type_metrics: List[Dict[str, float]] = []
    consensus_residues: List[str] = []
    consensus_types: List[str] = []

    for i in range(seq_len):
        wt_residue = target_sequence[i]
        column = [seq[i] for seq in msa_sequences]
        counts = Counter(column)

        wt_dict: Dict[str, float] = {}
        type_dict: Dict[str, float] = {}

        # Consensus residue
        if counts:
            cons_residue, cons_count = max(counts.items(), key=lambda x: x[1])
        else:
            cons_residue, cons_count = "-", 0
        cons_res_freq = cons_count / total_hits if total_hits else 0.0
        wt_dict["cons_residue"] = cons_residue
        if "consensus_freq" in wt_scores:
            wt_dict["consensus_freq"] = cons_res_freq
        consensus_residues.append(cons_residue)

        # WT-based scores
        if "identity" in wt_scores:
            wt_dict["identity"] = counts.get(wt_residue, 0) / total_hits if total_hits else 0.0

        if "entropy" in wt_scores or "norm_entropy" in wt_scores:
            probs = []
            for aa, c in counts.items():
                if aa == "-":
                    continue
                p = c / total_hits
                if p > 0:
                    probs.append(p)
            H = -sum(p * math.log(p, 2) for p in probs) if probs else 0.0
            if "entropy" in wt_scores:
                wt_dict["entropy"] = H
            if "norm_entropy" in wt_scores:
                H_max = math.log(len(probs), 2) if probs else 1.0
                wt_dict["norm_entropy"] = 1.0 - (H / H_max if H_max > 0 else 0.0)

        if "gap_fraction" in wt_scores:
            wt_dict["gap_fraction"] = counts.get("-", 0) / total_hits if total_hits else 0.0

        if "blosum_wt" in wt_scores:
            score_sum = 0.0
            for aa, c in counts.items():
                if aa == "-" or aa == wt_residue:
                    continue
                try:
                    s = BLOSUM62[wt_residue, aa]
                except (KeyError, IndexError):
                    try:
                        s = BLOSUM62[aa, wt_residue]
                    except (KeyError, IndexError):
                        s = 0
                score_sum += s * c
            wt_dict["blosum_wt"] = score_sum / total_hits if total_hits else 0.0

        if "blosum_pairwise" in wt_scores:
            aas = [(aa, c) for aa, c in counts.items() if aa != "-"]
            tot_pairs = 0
            score_sum = 0.0
            for idx_a, (aa1, c1) in enumerate(aas):
                for aa2, c2 in aas[idx_a:]:
                    if aa1 == aa2:
                        num_pairs = c1 * (c1 - 1) // 2 if idx_a == 0 else c1 * c2
                    else:
                        num_pairs = c1 * c2
                    try:
                        s = BLOSUM62[aa1, aa2]
                    except (KeyError, IndexError):
                        try:
                            s = BLOSUM62[aa2, aa1]
                        except (KeyError, IndexError):
                            s = 0
                    score_sum += s * num_pairs
                    tot_pairs += num_pairs
            wt_dict["blosum_pairwise"] = score_sum / tot_pairs if tot_pairs else 0.0

        wt_metrics.append(wt_dict)

        # Type-based scores
        wt_type = AA_TYPES.get(wt_residue.upper(), "X")

        if "type_freq" in type_scores:
            type_count = 0
            for aa, c in counts.items():
                aa_type = AA_TYPES.get(aa.upper(), "X")
                if aa_type == wt_type:
                    type_count += c
            type_dict["type_freq"] = type_count / total_hits if total_hits else 0.0

        if (
            "type_entropy" in type_scores
            or "type_consensus_freq" in type_scores
        ):
            class_counts: Dict[str, int] = {}
            for aa, c in counts.items():
                aa_type = AA_TYPES.get(aa.upper(), "X")
                class_counts[aa_type] = class_counts.get(aa_type, 0) + c

            if class_counts:
                cons_type, cons_type_count = max(class_counts.items(), key=lambda x: x[1])
            else:
                cons_type, cons_type_count = "X", 0
            cons_type_freq = cons_type_count / total_hits if total_hits else 0.0
            type_dict["cons_type"] = cons_type
            if "type_consensus_freq" in type_scores:
                type_dict["type_consensus_freq"] = cons_type_freq

            if "type_entropy" in type_scores:
                probs_t = []
                for _, c in class_counts.items():
                    p = c / total_hits
                    if p > 0:
                        probs_t.append(p)
                Ht = -sum(p * math.log(p, 2) for p in probs_t) if probs_t else 0.0
                type_dict["type_entropy"] = Ht

            consensus_types.append(cons_type)
        else:
            consensus_types.append("X")

        type_metrics.append(type_dict)

    return wt_metrics, type_metrics, "".join(consensus_residues), "".join(consensus_types)


def save_conservation_scores(
    wt_metrics: Sequence[Dict[str, float]],
    type_metrics: Sequence[Dict[str, float]],
    seq: str,
    pdb_ids: Optional[Sequence[Tuple[str, int, str]]],
    wt_scores: Iterable[str],
    type_scores: Iterable[str],
    output_file: str,
) -> None:
    """
    Save wild-type and type-based metrics per residue to a tab-delimited file.

    Parameters
    ----------
    wt_metrics, type_metrics : sequence of dict
        Per-residue metric dictionaries.
    seq : str
        Target sequence.
    pdb_ids : sequence of (chain_id, resseq, inscode), optional
        If provided, PDB identifiers are written for each residue.
    wt_scores, type_scores : iterable of str
        Names of metrics to write.
    output_file : str
        Output text file path.
    """
    wt_scores = list(wt_scores)
    type_scores = list(type_scores)

    with open(output_file, "w") as f:
        header_cols = ["seq_index", "res_name", "AA_consensus", "Class_consensus"]
        if pdb_ids is not None:
            header_cols.insert(1, "PDB_chain")
            header_cols.insert(2, "PDB_resnum")

        for s in wt_scores:
            header_cols.append(WT_HEADER_MAP.get(s, f"WT_{s}"))

        for s in type_scores:
            header_cols.append(TYPE_HEADER_MAP.get(s, f"TYPE_{s}"))


        f.write("\t".join(header_cols) + "\n")

        if pdb_ids is not None:
            iterable = zip(wt_metrics, type_metrics, seq, pdb_ids)
        else:
            iterable = zip(wt_metrics, type_metrics, seq)

        for i, items in enumerate(iterable, 1):
            if pdb_ids is not None:
                wt_dict, type_dict, residue, pdb_id = items
                chain_id, res_num, _ = pdb_id
                cons_res = wt_dict.get("cons_residue", "X")
                cons_type = type_dict.get("cons_type", "X")
                row = [str(i), chain_id, str(res_num), residue, cons_res, cons_type]
            else:
                wt_dict, type_dict, residue = items
                cons_res = wt_dict.get("cons_residue", "X")
                cons_type = type_dict.get("cons_type", "X")
                row = [str(i), residue, cons_res, cons_type]

            for s in wt_scores:
                row.append(f"{wt_dict.get(s, 0.0):.4f}")
            for s in type_scores:
                row.append(f"{type_dict.get(s, 0.0):.4f}")

            f.write("\t".join(row) + "\n")


def write_conservation_pdb(
    original_pdb_file: str,
    conservation_scores_wt: Sequence[float],
    pdb_ids: Sequence[Tuple[str, int, str]],
    output_file: str,
) -> None:
    """
    Write a new PDB file with conservation scores in the B-factor column.

    ATOM/HETATM lines are updated so that:
    - Occupancy (cols 55–60) is set to 1.00
    - B-factor (cols 61–66) is set to the corresponding conservation score

    Parameters
    ----------
    original_pdb_file : str
        Input PDB file.
    conservation_scores_wt : sequence of float
        Per-residue conservation scores for the wild-type (or selected metric).
    pdb_ids : sequence of (chain_id, resseq, inscode)
        Residue identifiers corresponding to scores.
    output_file : str
        Output PDB path.
    """
    score_map: Dict[Tuple[str, int, str], float] = {}
    for score, pdb_id in zip(conservation_scores_wt, pdb_ids):
        score_map[pdb_id] = score

    with open(original_pdb_file, "r") as f_in, open(output_file, "w") as f_out:
        for line in f_in:
            if line.startswith(("ATOM", "HETATM")):
                if len(line) < 66:
                    f_out.write(line)
                    continue
                chain_id = line[21:22].strip()
                try:
                    res_num = int(line[22:26].strip())
                except ValueError:
                    f_out.write(line)
                    continue
                ins_code = line[26:27].strip()
                pdb_key = (chain_id, res_num, ins_code)
                score = score_map.get(pdb_key)

                if score is not None:
                    new_b_factor = f"{score:.2f}".rjust(6)
                    new_occupancy = f"{1.00:.2f}".rjust(6)
                    modified_line = (
                        line[:54] + new_occupancy + new_b_factor + line[66:]
                    )
                    f_out.write(modified_line)
                else:
                    f_out.write(line)
            else:
                f_out.write(line)


def extract_metric_list(metric_dicts: Sequence[Dict[str, float]], name: str) -> List[float]:
    """Extract a list of values for a named metric from a list of dicts."""
    return [d.get(name, 0.0) for d in metric_dicts]


def should_plot(metric_name: str, selection: str) -> bool:
    """
    Decide whether to plot a metric based on a selection string.

    selection can be 'all', 'none', or a comma-separated list of metric names.
    """
    if selection == "all":
        return True
    if selection == "none":
        return False
    return metric_name in selection.split(",")


def generate_matplotlib_heatmap_png(
    conservation_scores_1d: Sequence[float],
    seq: str,
    output_file: str,
    title_suffix: str,
    base_name: str,
    metric_name: str,
    residues_per_row: int = 50,
    cmap_name="Blues",
) -> None:
    """
    Generate a ConSurf-style conservation heatmap PNG using Matplotlib.

    Parameters
    ----------
    conservation_scores_1d : sequence of float
        1D list/array of scores (WT or type).
    seq : str
        Sequence used for residue labels.
    output_file : str
        Output PNG path.
    title_suffix : str
        Text to append to the heatmap title.
    base_name : str
        Base name of the dataset (used in the title).
    metric_name : str
        Metric name (affects color scaling).
    residues_per_row : int, optional
        Number of residues per row in the heatmap grid.
    """
    scores = np.array(conservation_scores_1d, dtype=float)
    seq_len = len(seq)

    # Color scaling
    n_grades = 10
    vmin, vmax = 0.0, 1.0
    if metric_name == "entropy":
        vmin, vmax = 0.0, math.log(20, 2)
    elif metric_name == "type_entropy":
        vmin, vmax = 0.0, math.log(5, 2)
    elif metric_name in ("blosum_wt", "blosum_pairwise"):
        vmin, vmax = -4.0, 11.0

    bounds = np.linspace(vmin, vmax, n_grades + 1)
    cmap = plt.cm.get_cmap(cmap_name, n_grades)
    norm = mcolors.BoundaryNorm(bounds, cmap.N)

    rows = (seq_len + residues_per_row - 1) // residues_per_row
    padded_len = rows * residues_per_row
    padded_scores = np.pad(scores, (0, padded_len - seq_len), constant_values=np.nan)
    heatmap_data = padded_scores.reshape(rows, residues_per_row)

    fig_width = 10
    fig_height = 1 + rows * 0.3
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=300)

    cax = ax.imshow(heatmap_data, cmap=cmap, norm=norm, aspect="equal")

    # Residue labels
    for r in range(rows):
        for c in range(residues_per_row):
            linear_index = r * residues_per_row + c
            if linear_index < seq_len:
                residue = seq[linear_index]
                score = heatmap_data[r, c]
                if np.isnan(score):
                    text_color = "black"
                elif score < 0.3 or score > 0.7:
                    text_color = "white"
                else:
                    text_color = "black"
                ax.text(
                    c,
                    r,
                    residue,
                    ha="center",
                    va="center",
                    color=text_color,
                    fontsize=7,
                    fontweight="bold",
                )

    col_labels = [str(i + 1) if (i + 1) % 10 == 0 else "" for i in range(residues_per_row)]
    ax.set_xticks(np.arange(residues_per_row))
    ax.set_xticklabels(col_labels, fontsize=10)
    ax.tick_params(axis="x", length=0)
    ax.xaxis.set_label_position("top")
    ax.xaxis.tick_top()

    row_labels = [str(r * residues_per_row + 1) for r in range(rows)]
    ax.set_yticks(np.arange(rows))
    ax.set_yticklabels(row_labels, fontsize=10)
    ax.tick_params(axis="y", length=0)
    for side in ("top", "right", "bottom", "left"):
        ax.spines[side].set_visible(False)

    cbar = fig.colorbar(cax, orientation="horizontal", fraction=0.065, pad=0.1)
    cbar.set_ticks(bounds)
    if metric_name in {
        "identity",
        "type_freq",
        "consensus_freq",
        "type_consensus_freq",
        "gap_fraction",
        "norm_entropy",
    }:
        tick_labels = [f"{int(b * 100)}%" for b in bounds]
    else:
        tick_labels = [f"{b:.2f}" for b in bounds]
    cbar.set_ticklabels(tick_labels)
    cbar.set_label(f"{title_suffix}", fontsize=10)
    cbar.ax.tick_params(labelsize=10)

    ax.set_title(
        f"{title_suffix} Heatmap for {base_name}",
        fontsize=10,
        pad=20,
        fontweight="bold",
    )

    plt.tight_layout(pad=0.5)
    plt.savefig(output_file)
    plt.close(fig)


# ------------------------- ProDy-style callable -------------------------

def computeConservationFromMSA(
    sequence_file: str,
    pdb_file: Optional[str] = None,
    xml_file: Optional[str] = None,
    msa_file: Optional[str] = None,
    wt_scores: str = "identity,type_freq",
    type_scores: str = "type_freq",
    wt_bfactor: str = "identity",
    type_bfactor: str = "type_freq",
    wt_heatmaps: str = "identity",
    type_heatmaps: str = "type_freq",
    residues_per_row: int = 50,
    out_prefix: Optional[str] = None,
    return_data: bool = False,
    wt_cmap: str = "Blues",
    type_cmap: str = "Blues",    
) -> Union[int, Dict[str, Any]]:
    """
    Compute residue-level evolutionary conservation from an MSA.

    This high-level function loads the target sequence (and optional PDB),
    obtains an MSA from BLAST XML or a precomputed MSA file, computes
    wild-type and type-based conservation metrics per residue, writes
    results to disk (TXT, PNG, PDB), and optionally returns all data in
    a structured dictionary for programmatic use.

    Parameters
    ----------
    sequence_file : str
        Path to the file containing the target (wild-type) sequence,
        one-letter amino-acid codes on a single line.
    pdb_file : str, optional
        Path to a single-chain PDB file to extract residue IDs.
        If provided, conservation scores can be mapped into B-factors.
    xml_file : str, optional
        Path to a BLAST XML results file. If given and ``msa_file`` is
        not provided, an MSA is generated from the BLAST hits.
    msa_file : str, optional
        Path to a pre-computed MSA file (one sequence per line,
        aligned to the target). If provided, BLAST XML parsing is
        skipped.
    wt_scores : str, optional
        Comma-separated wild-type metrics to compute, e.g.:
        ``"identity,entropy,norm_entropy,consensus_freq,gap_fraction,"
        "blosum_wt,blosum_pairwise"``.
    type_scores : str, optional
        Comma-separated type-based metrics to compute, e.g.:
        ``"type_freq,type_entropy,type_consensus_freq"``.
    wt_bfactor : str, optional
        Wild-type metric name to write into B-factors, or ``"none"``.
    type_bfactor : str, optional
        Type metric name to write into B-factors, or ``"none"``.
    wt_heatmaps : str, optional
        Comma-separated wild-type metrics for heatmaps, or ``"all"`` /
        ``"none"``.
    type_heatmaps : str, optional
        Comma-separated type metrics for heatmaps, or ``"all"`` /
        ``"none"``.
    residues_per_row : int, optional
        Number of residues per row in heatmap layout.
    out_prefix : str, optional
        Optional prefix for output files. If not provided, the prefix
        is derived from ``xml_file`` or ``msa_file``.
    return_data : bool, optional
        If ``False`` (default), return only the number of MSA sequences
        used. If ``True``, also return a dictionary containing all
        per-residue metrics, consensus sequences, PDB mapping, and
        output file paths.

    Returns
    -------
    int or dict
        If ``return_data`` is ``False``, returns the number of sequences
        used in the conservation calculation.

        If ``return_data`` is ``True``, returns a dictionary:

        ``{
            "n_sequences": int,
            "target_sequence": str,
            "msa_sequences": list[str],
            "wt_scores": list[str],
            "type_scores": list[str],
            "wt_metrics": list[dict],
            "type_metrics": list[dict],
            "consensus": {
                "residues": str,
                "types": str,
            },
            "pdb": {
                "pdb_file": str or None,
                "pdb_ids": list[(chain, resnum, inscode)] or None,
            },
            "output_files": {
                "scores_txt": str,
                "wt_pdb": list[str],
                "type_pdb": list[str],
                "wt_heatmaps": list[str],
                "type_heatmaps": list[str],
            },
        }``
    """
    if xml_file is None and msa_file is None:
        raise ValueError("Either 'xml_file' or 'msa_file' must be provided.")

    wt_score_set = set(wt_scores.split(",")) if wt_scores else set()
    type_score_set = set(type_scores.split(",")) if type_scores else set()

    # 1. Load target sequence
    target_sequence = read_target_sequence(sequence_file)

    # 2. Optional PDB mapping
    pdb_ids: Optional[List[Tuple[str, int, str]]] = None
    if pdb_file is not None:
        pdb_ids = read_pdb_residue_ids(pdb_file, target_sequence)

    # 3. Obtain MSA (from BLAST XML or precomputed file)
    if xml_file is not None and msa_file is None:
        base_name = os.path.splitext(os.path.basename(xml_file))[0]
        msa_output_file = f"{base_name}.seqmsa"
        total_hits = generate_msa_from_blast_xml(xml_file, target_sequence, msa_output_file)
        msa_sequences = load_msa(msa_output_file)
    else:
        if msa_file is None:
            raise ValueError("msa_file must be provided when xml_file is None.")
        base_name = os.path.splitext(os.path.basename(msa_file))[0]
        msa_output_file = msa_file
        msa_sequences = load_msa(msa_file)
        total_hits = len(msa_sequences)

    if total_hits == 0:
        return 0 if not return_data else {
            "n_sequences": 0,
            "target_sequence": target_sequence,
            "msa_sequences": [],
            "wt_scores": sorted(wt_score_set),
            "type_scores": sorted(type_score_set),
            "wt_metrics": [],
            "type_metrics": [],
            "consensus": {"residues": "", "types": ""},
            "pdb": {"pdb_file": pdb_file, "pdb_ids": pdb_ids},
            "output_files": {
                "scores_txt": None,
                "wt_pdb": [],
                "type_pdb": [],
                "wt_heatmaps": [],
                "type_heatmaps": [],
            },
        }

    if len(msa_sequences[0]) != len(target_sequence):
        raise ValueError(
            f"Length mismatch: MSA sequence length ({len(msa_sequences[0])}) "
            f"must match target length ({len(target_sequence)})."
        )

    prefix = out_prefix or base_name
    scores_output_file = f"{prefix}_conservation_scores.txt"

    # 4. Per-residue metrics
    wt_metrics, type_metrics, consensus_seq_res, consensus_seq_type = (
        calculate_conservation_scores(
            msa_sequences,
            target_sequence,
            total_hits,
            wt_score_set,
            type_score_set,
        )
    )

    # 5. Save scores table
    save_conservation_scores(
        wt_metrics,
        type_metrics,
        target_sequence,
        pdb_ids,
        wt_score_set,
        type_score_set,
        scores_output_file,
    )

    # Track created files for optional return
    wt_pdb_paths: List[str] = []
    type_pdb_paths: List[str] = []
    wt_png_paths: List[str] = []
    type_png_paths: List[str] = []

    # 6. Heatmaps (WT)
    for s in wt_score_set:
        if should_plot(s, wt_heatmaps):
            vals = extract_metric_list(wt_metrics, s)
            label_seq = consensus_seq_res if s == "consensus_freq" else target_sequence
            from_header = WT_HEADER_MAP.get(s, f"WT_{s}")
            out_png = f"{prefix}_{from_header}_heatmap.png"
            title_suffix = WT_TITLE_MAP.get(s, f"WT {s}")
            generate_matplotlib_heatmap_png(
                vals,
                label_seq,
                out_png,
                title_suffix,
                prefix,
                metric_name=s,
                residues_per_row=residues_per_row,
                cmap_name=wt_cmap,
            )
            wt_png_paths.append(out_png)

    # 7. Heatmaps (TYPE)
    for s in type_score_set:
        if should_plot(s, type_heatmaps):
            vals = extract_metric_list(type_metrics, s)
            label_seq = consensus_seq_res if s == "type_consensus_freq" else target_sequence
            out_png = f"{prefix}_TYPE_{s}_heatmap.png"
            title_suffix = TYPE_TITLE_MAP.get(s, f"TYPE {s}")
            generate_matplotlib_heatmap_png(
                vals,
                label_seq,
                out_png,
                title_suffix,
                prefix,
                metric_name=s,
                residues_per_row=residues_per_row,
                cmap_name=type_cmap,
            )
            type_png_paths.append(out_png)

    # 8. PDB outputs
    if pdb_file is not None and pdb_ids is not None:
        if wt_bfactor and wt_bfactor != "none":
            wt_list = extract_metric_list(wt_metrics, wt_bfactor)
            pdb_output_file_wt = f"{prefix}_WT_{wt_bfactor}_bfactor.pdb"
            write_conservation_pdb(pdb_file, wt_list, pdb_ids, pdb_output_file_wt)
            wt_pdb_paths.append(pdb_output_file_wt)

        if type_bfactor and type_bfactor != "none":
            type_list = extract_metric_list(type_metrics, type_bfactor)
            pdb_output_file_type = f"{prefix}_TYPE_{type_bfactor}_bfactor.pdb"
            write_conservation_pdb(pdb_file, type_list, pdb_ids, pdb_output_file_type)
            type_pdb_paths.append(pdb_output_file_type)

    if not return_data:
        return total_hits

    # 9. Rich return object
    return {
        "n_sequences": total_hits,
        "target_sequence": target_sequence,
        "msa_sequences": msa_sequences,
        "wt_scores": sorted(wt_score_set),
        "type_scores": sorted(type_score_set),
        "wt_metrics": wt_metrics,
        "type_metrics": type_metrics,
        "consensus": {
            "residues": consensus_seq_res,
            "types": consensus_seq_type,
        },
        "pdb": {
            "pdb_file": pdb_file,
            "pdb_ids": pdb_ids,
        },
        "output_files": {
            "scores_txt": scores_output_file,
            "wt_pdb": wt_pdb_paths,
            "type_pdb": type_pdb_paths,
            "wt_heatmaps": wt_png_paths,
            "type_heatmaps": type_png_paths,
        },
    }



# ------------------------- CLI wrapper -------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Compute residue-level conservation metrics from an MSA generated "
            "from BLAST XML or from a pre-computed MSA file, and optionally "
            "write PDB B-factors and PNG heatmaps."
        )
    )

    parser.add_argument(
        "-s",
        "--sequence-file",
        required=True,
        help="Path to the file containing the target sequence (one-line, one-letter codes).",
    )
    parser.add_argument(
        "-p",
        "--pdb-file",
        required=False,
        help=(
            "Path to a single-chain PDB file to extract residue IDs. "
            "Optional if only sequence-based scores are needed."
        ),
    )
    parser.add_argument(
        "--wt-scores",
        type=str,
        default="identity,type_freq",
        help=(
            "Comma-separated wild-type metrics to compute: "
            "identity,entropy,norm_entropy,consensus_freq,"
            "gap_fraction,blosum_wt,blosum_pairwise"
        ),
    )
    parser.add_argument(
        "--type-scores",
        type=str,
        default="type_freq",
        help=(
            "Comma-separated type-based metrics to compute: "
            "type_freq,type_entropy,type_consensus_freq"
        ),
    )
    parser.add_argument(
        "--wt-bfactor",
        type=str,
        default="identity",
        help="Wild-type metric to write into B-factors (or 'none').",
    )
    parser.add_argument(
        "--type-bfactor",
        type=str,
        default="type_freq",
        help="Type metric to write into B-factors (or 'none').",
    )
    parser.add_argument(
        "--wt-heatmaps",
        type=str,
        default="identity",
        help="Comma-separated WT metrics for heatmaps, or 'all'/'none'.",
    )
    parser.add_argument(
        "--type-heatmaps",
        type=str,
        default="type_freq",
        help="Comma-separated type metrics for heatmaps, or 'all'/'none'.",
    )
    parser.add_argument(
        "--residues-per-row",
        type=int,
        default=50,
        help="Number of residues per row in heatmap layout.",
    )
    parser.add_argument(
        "--wt-cmap",
        type=str,
        default="Blues",
        help="Matplotlib colormap name for WT heatmaps (default: Blues).",
    )

    parser.add_argument(
        "--type-cmap",
        type=str,
        default="Blues",
        help="Matplotlib colormap name for TYPE heatmaps (default: Blues).",
    )
    parser.add_argument(
        "--out-prefix",
        type=str,
        default=None,
        help="Optional prefix for all output files.",
    )

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "-x",
        "--xml-file",
        help="Path to BLAST XML results file; an MSA will be generated.",
    )
    input_group.add_argument(
        "-m",
        "--msa-file",
        help=(
            "Path to a pre-computed MSA file (.seqmsa format, one sequence per line), "
            "aligned to the target sequence."
        ),
    )

    args = parser.parse_args()

    try:
        total = computeConservationFromMSA(
            sequence_file=args.sequence_file,
            pdb_file=args.pdb_file,
            xml_file=args.xml_file,
            msa_file=args.msa_file,
            wt_scores=args.wt_scores,
            type_scores=args.type_scores,
            wt_bfactor=args.wt_bfactor,
            type_bfactor=args.type_bfactor,
            wt_heatmaps=args.wt_heatmaps,
            type_heatmaps=args.type_heatmaps,
            residues_per_row=args.residues_per_row,
            out_prefix=args.out_prefix,
            wt_cmap=args.wt_cmap,
            type_cmap=args.type_cmap,
        )
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Processed {total} MSA sequences.")


if __name__ == "__main__":
    main()


#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Ligand efficiency and ADME scoring tools for small molecules.

This module provides a convenience function to evaluate ligand efficiency,
basic ADME/physicochemical descriptors, and a composite drug-likeness
score from a multi-record SDF file using RDKit. It is intended for
quick triage of virtual screening results or docking hits.

For each molecule, the workflow:

1. Reads structures from an SDF file without initial sanitization.
2. Optionally standardizes each molecule using RDKit MolStandardize
   (cleanup, metal disconnection, reionization, normalization, tautomer
   canonicalization, and final sanitization).
3. Retrieves per-ligand binding energy and minimized RMSD from SD
   properties if present.
4. Computes a set of ADME-related descriptors, including:
   - Molecular weight (MW)
   - Crippen cLogP
   - Topological polar surface area (TPSA)
   - H-bond donors and acceptors (HBD, HBA)
   - Rotatable bond count (RB)
   - Fraction of sp3 carbons (fSP3)
   - Formal charge and heavy atom count
   - Aromatic ring, ring count, heteroatom count, and aromatic atom
     fraction
   - Rule-of-five style filters (Lipinski, Veber, Egan, lead-like)
   - Simple HIA / BBB likelihood classifiers (TPSA-based)
5. Runs several structural alert filters using RDKit filter catalogs:
   PAINS (A/B/C), Brenk, NIH, and ZINC.
6. Optionally computes:
   - QED score (Quantitative Estimate of Drug-likeness)
   - Synthetic accessibility (SA) score (if RDKit SA_Score is available)
   - ESOL-like aqueous solubility estimate (logS)
   - Simple aggregator-risk flag using a couple of SMARTS patterns.
7. Computes a per-ligand ligand efficiency (LE) defined as:
   LE = |energy_kcal_per_mol| / heavy_atoms
   when both energy and heavy atom count are available.
8. Maps the descriptors and alert flag into a normalized 0–10 "drug_score"
   using weighted heuristic scoring functions for:
   - LE
   - cLogP
   - TPSA
   - MW
   - HBA/HBD
   - RB
   - Presence or absence of structural alerts
9. Derives a qualitative summary label from the score and key filters:
   "Good (drug-like)", "Medium (optimize)", or "Poor (low priority)".

The main entry point is:

* :func:`calcLigandEfficiencyFromSDF` – run the full workflow on an SDF
  file and write a CSV table with one row per molecule.

The same module can also be used as a standalone script from the command
line; in that case, see the docstring of :func:`calcLigandEfficiencyFromSDF`
for details of the generated CSV columns and how to control the input
and output paths.
"""
__author__ = 'Anupam Banerjee'
__email__ = ['anupam.banerjee@stonybrook.edu']

from __future__ import annotations

import argparse
import csv
import math
import sys
from typing import List, Tuple, Optional

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors as rdmd, Lipinski
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem import SDWriter
from rdkit.Chem import QED
from rdkit.Chem.MolStandardize import rdMolStandardize as std

__all__ = ["calcLigandEfficiencyFromSDF"]

# Optional SA score
try:
    from rdkit.Contrib.SA_Score import sascorer

    HAVE_SA = True
except Exception:
    HAVE_SA = False

# ------------------------- Controls -------------------------

WEIGHTS_FULL = {
    "LE": 0.25,
    "cLogP": 0.15,
    "TPSA": 0.15,
    "MW": 0.15,
    "HBAHBD": 0.10,
    "RB": 0.10,
    "Alerts": 0.10,
}


def renormalize_weights(le_present: bool):
    """If LE is missing, drop LE and renormalize remaining weights."""
    if le_present:
        return WEIGHTS_FULL
    w = WEIGHTS_FULL.copy()
    del w["LE"]
    s = sum(w.values())
    for k in w:
        w[k] = w[k] / s
    return w


# ------------------------- Helpers -------------------------


def maybe_quiet_rdkit(enable: bool):
    if enable:
        RDLogger.DisableLog("rdApp.warning")


def pick_id(mol: Chem.Mol, use_full_title: bool = False) -> str:
    for key in [
        "zinc_id",
        "ZINC_ID",
        "ZINC",
        "ZINCID",
        "id",
        "ID",
        "name",
        "Name",
        "MOL_NAME",
        "MOL_ID",
        "MolPort",
        "MolPortID",
    ]:
        if mol.HasProp(key):
            v = mol.GetProp(key).strip()
            if v:
                return v
    title = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
    if not title:
        return ""
    return title if use_full_title else title.split()[0]


def get_float_prop(mol: Chem.Mol, prop: str) -> Optional[float]:
    if not mol.HasProp(prop):
        return None
    try:
        return float(mol.GetProp(prop).strip())
    except Exception:
        return None


def build_catalog(cat_type) -> FilterCatalog:
    params = FilterCatalogParams()
    params.AddCatalog(cat_type)
    return FilterCatalog(params)


# Filter catalogs
FC_PAINS_A = build_catalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
FC_PAINS_B = build_catalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
FC_PAINS_C = build_catalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
FC_BRENK = build_catalog(FilterCatalogParams.FilterCatalogs.BRENK)
FC_NIH = build_catalog(FilterCatalogParams.FilterCatalogs.NIH)
FC_ZINC = build_catalog(FilterCatalogParams.FilterCatalogs.ZINC)


def list_alerts(
    mol: Chem.Mol,
) -> Tuple[List[str], List[str], List[str], List[str]]:
    pains: List[str] = []
    for fc in (FC_PAINS_A, FC_PAINS_B, FC_PAINS_C):
        pains.extend([e.GetDescription() for e in fc.GetMatches(mol)])
    brenk = [e.GetDescription() for e in FC_BRENK.GetMatches(mol)]
    nih = [e.GetDescription() for e in FC_NIH.GetMatches(mol)]
    zinc = [e.GetDescription() for e in FC_ZINC.GetMatches(mol)]
    return pains, brenk, nih, zinc


def standardize(m: Chem.Mol) -> Chem.Mol:
    if m is None:
        return None
    try:
        clean = std.Cleanup(m)
        clean = std.MetalDisconnector().Disconnect(clean)
        clean = std.Reionize(clean)
        clean = std.Normalize(clean)
        can = std.TautomerCanonicalizer()
        clean = can.Canonicalize(clean)
        Chem.SanitizeMol(clean)
        return clean
    except Exception:
        return m


def esol_logs(m: Chem.Mol) -> Optional[float]:
    try:
        mw = Descriptors.MolWt(m)
        logp = Crippen.MolLogP(m)
        rb = Lipinski.NumRotatableBonds(m)
        # robust aromatic atom count
        try:
            arom = rdmd.CalcNumAromaticAtoms(m)  # may not exist in some RDKit builds
        except Exception:
            arom = sum(1 for a in m.GetAtoms() if a.GetIsAromatic())
        heavy = m.GetNumHeavyAtoms() or 1
        arom_prop = arom / heavy
        logs = 0.16 - 0.63 * logp - 0.0062 * mw + 0.066 * rb - 0.74 * arom_prop
        return round(float(logs), 2)
    except Exception:
        return None


def aggregator_risk(m: Chem.Mol) -> bool:
    try:
        long_chain = Chem.MolFromSmarts("CCCCCCCC")
        catechol = Chem.MolFromSmarts("c1ccc(c(c1)O)O")
        if m.HasSubstructMatch(long_chain):
            return True
        if m.HasSubstructMatch(catechol) and Crippen.MolLogP(m) > 2.5:
            return True
        return False
    except Exception:
        return False


def calc_qed(m: Chem.Mol) -> Optional[float]:
    try:
        return round(float(QED.qed(m)), 3)
    except Exception:
        return None


def calc_sa(m: Chem.Mol) -> Optional[float]:
    if not HAVE_SA:
        return None
    try:
        return round(float(sascorer.calculateScore(m)), 2)
    except Exception:
        return None


def compute_adme(mol: Chem.Mol):
    Chem.SanitizeMol(mol)
    mw = Descriptors.MolWt(mol)
    clogp = Crippen.MolLogP(mol)
    tpsa = rdmd.CalcTPSA(mol)
    hba = Lipinski.NumHAcceptors(mol)
    hbd = Lipinski.NumHDonors(mol)
    rb = Lipinski.NumRotatableBonds(mol)
    fsp3 = rdmd.CalcFractionCSP3(mol)
    charge = Chem.GetFormalCharge(mol)
    heavy = mol.GetNumHeavyAtoms()
    arom_rings = rdmd.CalcNumAromaticRings(mol)

    # rings via ring info (robust across RDKit versions)
    ri = mol.GetRingInfo()
    rings = ri.NumRings()

    try:
        arom_atoms = rdmd.CalcNumAromaticAtoms(mol)
    except Exception:
        arom_atoms = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())

    arom_prop = (arom_atoms / heavy) if heavy else 0.0

    try:
        hetero = rdmd.CalcNumHeteroatoms(mol)
    except Exception:
        hetero = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() not in (1, 6))

    lipinski = (hbd <= 5) and (hba <= 10) and (mw <= 500) and (clogp <= 5)
    veber = (tpsa <= 140) and (rb <= 10)
    egan = (tpsa <= 131) and (clogp <= 5.88)
    leadlike = (
        (mw < 350)
        and (clogp <= 3)
        and (hbd <= 3)
        and (hba <= 6)
        and (rb <= 7)
    )
    hia_likely = (tpsa <= 130)
    bbb_likely = (tpsa <= 90)

    pains, brenk, nih, zinc = list_alerts(mol)
    qed_val = calc_qed(mol)
    sa_val = calc_sa(mol)
    logs = esol_logs(mol)
    agg = aggregator_risk(mol)
    has_alert = any((pains, brenk, nih, zinc))

    return {
        "MW": round(mw, 3),
        "cLogP": round(clogp, 3),
        "TPSA": round(tpsa, 2),
        "HBA": int(hba),
        "HBD": int(hbd),
        "RB": int(rb),
        "fSP3": round(fsp3, 3),
        "formal_charge": int(charge),
        "heavy_atoms": int(heavy),
        "aromatic_rings": int(arom_rings),
        "Rings": int(rings),
        "HeteroAtoms": int(hetero),
        "AromaticAtoms": int(arom_atoms),
        "AromaticProp": round(arom_prop, 3),
        "Lipinski_pass": lipinski,
        "Veber_pass": veber,
        "Egan_pass": egan,
        "Leadlike_pass": leadlike,
        "HIA_likely": hia_likely,
        "BBB_likely": bbb_likely,
        "PAINS_alerts": "; ".join(pains) if pains else "",
        "Brenk_alerts": "; ".join(brenk) if brenk else "",
        "NIH_alerts": "; ".join(nih) if nih else "",
        "ZINC_alerts": "; ".join(zinc) if zinc else "",
        "QED": qed_val,
        "SA_score": sa_val,
        "ESOL_logS": logs,
        "AggregatorRisk": bool(agg),
        "_has_alert": has_alert,
    }


def compute_le(energy_kcal, n_heavy):
    if energy_kcal is None or n_heavy is None or n_heavy <= 0:
        return None
    return abs(energy_kcal) / n_heavy


# ------------------------- Scoring -------------------------


def clamp01(x: float) -> float:
    return max(0.0, min(1.0, x))


def score_le(le: Optional[float]) -> float:
    if le is None:
        return 0.0
    return clamp01(le / 0.4)


def score_clogp(clogp: Optional[float]) -> float:
    if clogp is None:
        return 0.0
    return clamp01(1.0 - abs(clogp - 2.0) / 3.0)


def score_tpsa(tpsa: Optional[float]) -> float:
    if tpsa is None:
        return 0.0
    return clamp01(math.exp(-((tpsa - 80.0) ** 2) / (2 * 30.0**2)))


def score_mw(mw: Optional[float]) -> float:
    if mw is None:
        return 0.0
    return clamp01(1.0 - abs(mw - 325.0) / 225.0)


def score_hba_hbd(hba: Optional[int], hbd: Optional[int]) -> float:
    if hba is None or hbd is None:
        return 0.0
    penalty = 0.0
    if hba > 10:
        penalty += 0.1 * (hba - 10)
    if hbd > 5:
        penalty += 0.15 * (hbd - 5)
    return clamp01(1.0 - penalty)


def score_rb(rb: Optional[int]) -> float:
    if rb is None:
        return 0.0
    return clamp01(1.0 - max(0, rb - 8) / 8.0)


def score_alerts(has_alert: bool) -> float:
    return 0.0 if has_alert else 1.0


def compute_summary_score(
    le,
    clogp,
    tpsa,
    mw,
    hba,
    hbd,
    rb,
    has_alert,
) -> float:
    le_present = le is not None
    weights = renormalize_weights(le_present)
    s = 0.0
    if le_present:
        s += weights["LE"] * score_le(le)
    s += weights["cLogP"] * score_clogp(clogp)
    s += weights["TPSA"] * score_tpsa(tpsa)
    s += weights["MW"] * score_mw(mw)
    s += weights["HBAHBD"] * score_hba_hbd(hba, hbd)
    s += weights["RB"] * score_rb(rb)
    s += weights["Alerts"] * score_alerts(has_alert)
    return round(10.0 * s, 2)


def summarize(
    drug_score: float,
    le: Optional[float],
    lipinski_pass: bool,
    has_alert: bool,
) -> str:
    if drug_score >= 8.0 and lipinski_pass and not has_alert:
        if le is None or le >= 0.30:
            return "Good (drug-like)"
    if drug_score >= 6.0:
        return "Medium (optimize)"
    return "Poor (low priority)"


# ------------------------- ProDy-style callable -------------------------


def calcLigandEfficiencyFromSDF(
    sdf: str,
    out: Optional[str] = None,
    energy_field: str = "minimizedAffinity",
    rmsd_field: str = "minimizedRMSD",
    id_full: bool = False,
    limit: Optional[int] = None,
    extract_sdf: Optional[str] = None,
    quiet: bool = False,
) -> int:
    """
    Calculate ligand efficiency, ADME descriptors, and a composite
    drug-likeness score from an SDF file.

    This function reads a multi-record SDF file, optionally standardizes
    each molecule using RDKit MolStandardize, and computes for each
    record a panel of physicochemical descriptors, structural alerts,
    a ligand-efficiency value, and a heuristic 0–10 "drug_score"
    summarizing drug-likeness. Results are written to a CSV file with
    one row per molecule.

    Parameters
    ----------
    sdf : str
        Path to the input SDF file. The file may contain multiple
        molecular records. Structures are read without initial
        sanitization and are subsequently standardized inside this
        function.
    out : str, optional
        Path to the output CSV file. If not provided, the CSV name is
        derived from the SDF path by appending ``.le_adme.csv`` before
        the extension.
    energy_field : str, optional
        Name of the SD property holding the per-ligand binding energy
        in kcal/mol. The default is ``"minimizedAffinity"``. If the
        property is missing or cannot be parsed for a given molecule,
        its energy and LE are reported as ``None`` for that row.
    rmsd_field : str, optional
        Name of the SD property holding the per-ligand minimized RMSD,
        typically from a docking or refinement step. The default is
        ``"minimizedRMSD"``. If missing or invalid, the RMSD entry is
        reported as ``None``.
    id_full : bool, optional
        If ``False`` (default), the output ID field is derived from a
        subset of common SD properties (e.g. ``ZINC_ID``, ``id``,
        ``name``) or the first token of the internal title line
        (``_Name``). If ``True``, the full title line is used when
        falling back to ``_Name``.
    limit : int, optional
        If provided and positive, limit processing to the first
        ``limit`` records in the SDF. When used together with
        ``extract_sdf``, those same records can also be written into a
        smaller SDF file for inspection.
    extract_sdf : str, optional
        Path to a new SDF file to which the first ``limit`` standardized
        molecules are written. If ``limit`` is given but
        ``extract_sdf`` is omitted, the output SDF name is formed by
        appending ``.headN.sdf`` (where N = ``limit``) to the input
        SDF basename.
    quiet : bool, optional
        If ``True``, suppress RDKit warning messages (e.g. sanitization
        warnings) using ``RDLogger.DisableLog("rdApp.warning")``.
        Default is ``False``.

    Returns
    -------
    int
        The number of molecules successfully processed and written to
        the CSV (and optional extracted SDF).

    Output
    ------
    The CSV file contains (at minimum) the following columns, one row
    per molecule:

    - ``index``: 1-based index of the record in the original SDF.
    - ``id``: chosen identifier (ZINC ID, other SD property, or title).
    - ``energy_kcal_per_mol``: parsed binding energy from ``energy_field``.
    - ``heavy_atoms``: heavy-atom count used in LE calculation.
    - ``LE``: ligand efficiency, defined as
      ``|energy_kcal_per_mol| / heavy_atoms`` when both are available.
    - ``minimizedRMSD``: parsed value from ``rmsd_field``.
    - ``MW``, ``cLogP``, ``TPSA``, ``HBA``, ``HBD``, ``RB``, ``fSP3``,
      ``formal_charge``, ``aromatic_rings``: standard RDKit descriptors.
    - ``Lipinski_pass``, ``Veber_pass``, ``Egan_pass``, ``Leadlike_pass``:
      simple rule-of-five / lead-likeness flags.
    - ``HIA_likely``, ``BBB_likely``: approximate oral absorption and
      blood–brain-barrier likelihood flags based on TPSA thresholds.
    - ``PAINS_alerts``, ``Brenk_alerts``, ``NIH_alerts``, ``ZINC_alerts``:
      concatenated descriptions of matching alert filters (empty string
      if no alerts of that type are found).
    - ``drug_score``: composite 0–10 score derived from per-descriptor
      sub-scores and an alert penalty. Higher values indicate more
      drug-like profiles according to this heuristic.
    - ``summary``: a qualitative label summarizing the score and key
      filters: ``"Good (drug-like)"``, ``"Medium (optimize)"``, or
      ``"Poor (low priority)"``.
    - ``Rings``, ``HeteroAtoms``, ``AromaticAtoms``, ``AromaticProp``:
      expanded topological descriptors (ring count, heteroatom count,
      aromatic atom count, and fraction of aromatic atoms).
    - ``QED``: RDKit QED score, if available.
    - ``SA_score``: synthetic accessibility score, if the RDKit
      SA_Score contribution is installed; otherwise ``None``.
    - ``ESOL_logS``: approximate aqueous solubility prediction.
    - ``AggregatorRisk``: boolean flag for simple aggregator-risk
      motifs (long hydrophobic chains or catechol-like fragments with
      high cLogP).

    Notes
    -----
    - Molecule standardization is performed on each record before
      descriptor calculation. If standardization fails, the original
      molecule is used as a fallback.
    - In the presence of missing or invalid numeric properties, the
      corresponding descriptors and LE are set to ``None`` and
      contribute neutrally (score 0) to the drug_score.
    - This function does not modify the input SDF on disk; any
      standardized structures are written only to the optional
      ``extract_sdf`` output.
    """
    ...

    maybe_quiet_rdkit(quiet)

    out_csv = out or (sdf.rsplit(".", 1)[0] + ".le_adme.csv")
    try:
        csv_f = open(out_csv, "w", newline="", encoding="utf-8")
    except Exception as e:
        print(f"ERROR: cannot open CSV for writing '{out_csv}': {e}", file=sys.stderr)
        raise

    fieldnames = [
        "index",
        "id",
        "energy_kcal_per_mol",
        "heavy_atoms",
        "LE",
        "minimizedRMSD",
        "MW",
        "cLogP",
        "TPSA",
        "HBA",
        "HBD",
        "RB",
        "fSP3",
        "formal_charge",
        "aromatic_rings",
        "Lipinski_pass",
        "Veber_pass",
        "Egan_pass",
        "Leadlike_pass",
        "HIA_likely",
        "BBB_likely",
        "PAINS_alerts",
        "Brenk_alerts",
        "NIH_alerts",
        "ZINC_alerts",
        "drug_score",
        "summary",
        "Rings",
        "HeteroAtoms",
        "AromaticAtoms",
        "AromaticProp",
        "QED",
        "SA_score",
        "ESOL_logS",
        "AggregatorRisk",
    ]

    writer = csv.DictWriter(csv_f, fieldnames=fieldnames)
    writer.writeheader()

    suppl = Chem.SDMolSupplier(sdf, removeHs=False, sanitize=False)
    if suppl is None:
        print(f"ERROR: cannot open SDF: {sdf}", file=sys.stderr)
        csv_f.close()
        raise IOError(f"Cannot open SDF: {sdf}")

    sdw = None
    n_to_process = limit if (limit and limit > 0) else None

    if n_to_process:
        extract_path = extract_sdf or (
            sdf.rsplit(".", 1)[0] + f".head{n_to_process}.sdf"
        )
        try:
            sdw = SDWriter(extract_path)
        except Exception as e:
            print(f"ERROR: cannot open SDF writer '{extract_path}': {e}", file=sys.stderr)
            csv_f.close()
            raise

    processed = 0
    try:
        for idx, mol in enumerate(suppl, start=1):
            if n_to_process and processed >= n_to_process:
                break
            if mol is None:
                continue

            mol_std = standardize(mol)
            energy_kcal = get_float_prop(mol_std, energy_field)
            rmsd_val = get_float_prop(mol_std, rmsd_field)

            adme = compute_adme(mol_std)
            le_val = compute_le(energy_kcal, adme["heavy_atoms"])
            drug_score = compute_summary_score(
                le=le_val,
                clogp=adme["cLogP"],
                tpsa=adme["TPSA"],
                mw=adme["MW"],
                hba=adme["HBA"],
                hbd=adme["HBD"],
                rb=adme["RB"],
                has_alert=adme["_has_alert"],
            )
            label = summarize(
                drug_score, le_val, adme["Lipinski_pass"], adme["_has_alert"]
            )

            row = {
                "index": idx,
                "id": pick_id(mol_std, id_full),
                "energy_kcal_per_mol": energy_kcal,
                "heavy_atoms": adme["heavy_atoms"],
                "LE": le_val,
                "minimizedRMSD": rmsd_val,
                **{
                    k: v
                    for k, v in adme.items()
                    if k not in ("heavy_atoms", "_has_alert")
                },
                "drug_score": drug_score,
                "summary": label,
            }

            writer.writerow(row)

            if sdw is not None:
                title = row["id"]
                if title:
                    mol_std.SetProp("_Name", title)
                sdw.write(mol_std)

            processed += 1
    finally:
        csv_f.close()
        if sdw is not None:
            sdw.close()

    return processed


# ------------------------- CLI wrapper -------------------------


def main():
    ap = argparse.ArgumentParser(
        description="Ligand Efficiency + ADME + Summary scoring from SDF "
        "(expanded descriptors/alerts)."
    )
    ap.add_argument("sdf", help="input multi-record SDF file")
    ap.add_argument(
        "-o", "--out", default=None, help="output CSV (default: .le_adme.csv)"
    )
    ap.add_argument(
        "--energy-field",
        default="minimizedAffinity",
        help="SD property with energy in kcal/mol (default: minimizedAffinity)",
    )
    ap.add_argument(
        "--rmsd-field",
        default="minimizedRMSD",
        help="SD property with minimized RMSD (default: minimizedRMSD)",
    )
    ap.add_argument("--id-full", action="store_true", help="use full title line as ID")
    ap.add_argument(
        "--limit",
        type=int,
        default=None,
        help=(
            "if set, process ONLY the first N records and also extract those N to a new SDF"
        ),
    )
    ap.add_argument(
        "--extract-sdf",
        default=None,
        help="custom path for extracted SDF when --limit is used (default: .headN.sdf)",
    )
    ap.add_argument(
        "--quiet", action="store_true", help="suppress RDKit warning logs"
    )

    args = ap.parse_args()

    processed = calcLigandEfficiencyFromSDF(
        sdf=args.sdf,
        out=args.out,
        energy_field=args.energy_field,
        rmsd_field=args.rmsd_field,
        id_full=args.id_full,
        limit=args.limit,
        extract_sdf=args.extract_sdf,
        quiet=args.quiet,
    )

    if args.limit:
        print(f"Processed and extracted first {processed} records.")
    else:
        print(f"Processed {processed} records.")


if __name__ == "__main__":
    main()


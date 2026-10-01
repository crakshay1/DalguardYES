#!/usr/bin/env python3
"""
RiboGuard GA Engine
===================

Expected input:
    A list of candidate dictionaries, one dictionary per candidate.

Minimum candidate format:
    {
        "rbs": "UACAAG",
        "spacer": "AAUAAA"
    }

Recommended candidate format:
    {
        "five_prime_flank": "AAUAAU",
        "rbs": "UACAAG",
        "spacer": "AAUAAA",
        "cds_start": "AUGGCUACUAAAGAAAACGCU"
    }

The GA evaluates candidates using a Salis-inspired objective:
    T-score = orthogonal_TIR - wt_penalty_constant * WT_TIR

It outputs dashboard-ready files:
    - rbs_dataset.json     React-friendly combined dataset
    - candidates.csv       ranked candidates
    - fitness.csv          generation history
    - landscape.csv        scatter plot data
    - binding_sites.csv    all binding-site results for top candidates

Notes:
    - Folding uses LinearFold if available, otherwise ViennaRNA.
    - Delta-G duplex and standby calculations require ViennaRNA.
    - No pseudo-folding or approximate duplex fallback is used in this version.
    - The code uses RNA internally. DNA T is normalized to RNA U.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import multiprocessing as _mp
import os
import random
import subprocess
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Any

# Force fork so worker processes inherit already-imported ViennaRNA state
# instead of re-importing under spawn.
_FORK_CTX = _mp.get_context("fork")

# -------------------------
# Constants / defaults
# -------------------------

BASES = ["A", "U", "G", "C"]

DEFAULT_ORTH_ASD = "ACTTGTATA"
DEFAULT_WT_ASD = "ACCTCCTTA"
DEFAULT_FLANK = "AAUAAU"  # legacy alias for older callers
DEFAULT_FIVE_PRIME_UTR = DEFAULT_FLANK
DEFAULT_STANDBY = ""
DEFAULT_CDS_START = "AUGGCUACUAAAGAAAACGCUACUGCU"

# Physical constants requested by team
R_KCAL = 1.987e-3       # kcal/mol/K
TEMP_K = 310.15         # 37 C
RT_PHYSICAL = R_KCAL * TEMP_K

# Salis-style prefactor
R0 = 2500.0

# Objective:
# T-score = orthogonal_TIR - WT_PENALTY_CONSTANT * WT_TIR
DEFAULT_WT_PENALTY_CONSTANT = 0.35

# Spacing penalty constants
OPTIMAL_SPACING = 5.0

# Start codon placeholder/lookup
USE_START_LOOKUP = True
START_CODON_DG = {
    "AUG": -1.194,
    "GUG": -0.0748,
    "UUG": -0.0435,
    "CUG": -0.03406,
}

# Long-range filter
LONG_RANGE_THRESHOLD_NT = 35


# -------------------------
# Data classes
# -------------------------

@dataclass
class BindingSiteResult:
    candidate_id: str
    anti_sd_type: str
    site_rank: int
    mrna_start_0based: int
    mrna_end_0based_exclusive: int
    mrna_start_1based: int
    mrna_end_1based: int
    asd_start_5p_0based: int
    asd_start_5p_1based: int
    asd_end_5p_0based_exclusive: int
    asd_end_5p_1based: int
    overlap: int
    pairing_score: float
    aligned_spacing: float
    d: float
    dG_spacing: float
    dG_mrna_unfolding: float
    dG_duplex_candidate: float
    dG_start: float
    dG_standby: float
    dG_total: float
    tir: float
    is_best_site: bool


@dataclass
class CandidateEval:
    candidate_id: str

    # Canonical sequence schema
    five_prime_utr: str
    standby: str
    rbs_left: str
    rbs_core: str
    rbs_right: str
    rbs: str
    spacer: str
    cds_start: str

    # Assembled sequences / folding
    full_seq: str
    utr_seq: str
    structure: str
    mfe: float
    backend: str
    rbs_start: int
    rbs_end: int
    aug_start: int
    aug_end: int
    rbs_access: float
    aug_access: float
    long_range_flag: bool
    long_range_pairs: str

    # Candidate-level thermodynamics
    dG_duplex_orth: float
    dG_duplex_wt: float
    dG_start: float
    dG_standby_orth: float
    dG_standby_wt: float

    # Best orthogonal binding-site terms
    best_dG_spacing_orth: float
    best_dG_mrna_unfolding_orth: float
    dG_total_orth: float

    # Best WT binding-site terms
    best_dG_spacing_wt: float
    best_dG_mrna_unfolding_wt: float
    dG_total_wt: float

    # Expression / ranking
    orth_tir: float
    wt_tir: float
    t_score: float
    fitness_for_selection: float

    # Legacy/dashboard aliases retained internally for compatibility.
    # They are intentionally omitted from the clean workbook export.
    dG_standby: float
    best_dG_spacing: float
    best_dG_mrna_unfolding: float
    dG_total: float
    orthScore: float
    wtLeakage: float
    fitness: float


# -------------------------
# Basic utilities
# -------------------------

def normalize_rna(seq: str) -> str:
    """Uppercase, DNA-to-RNA, keep only A/U/G/C."""
    seq = (seq or "").upper().replace("T", "U")
    return "".join(ch for ch in seq if ch in "AUGC")


def reverse_complement(seq: str) -> str:
    seq = normalize_rna(seq)
    comp = {"A": "U", "U": "A", "G": "C", "C": "G"}
    return "".join(comp[b] for b in reversed(seq))


def clamp01(x: float) -> float:
    try:
        x = float(x)
    except Exception:
        return 0.0
    return max(0.0, min(1.0, x))


def safe_log10(x: float) -> float:
    if x <= 0 or math.isnan(x):
        return 0.0
    if math.isinf(x):
        return 308.0
    return math.log10(x)


def signed_log10_score(x: float) -> float:
    """Stable scalar for GA selection if raw T-score is huge or negative."""
    if math.isinf(x):
        return 308.0 if x > 0 else -308.0
    if x >= 0:
        return math.log10(x + 1.0)
    return -math.log10(abs(x) + 1.0)


def gc_fraction(seq: str) -> float:
    seq = normalize_rna(seq)
    return 0.0 if not seq else (seq.count("G") + seq.count("C")) / len(seq)


# -------------------------
# Folding backend
# -------------------------

def fold_sequence(sequence: str) -> Tuple[str, float, str]:
    """
    Fold a single RNA sequence.

    Backend priority:
      1) LinearFold if ./LinearFold/bin/linearfold_c exists
      2) ViennaRNA Python package

    No pseudo-folding fallback is used. If neither backend is available,
    this function raises RuntimeError so the result cannot silently rely on
    fake folding.
    """
    seq = normalize_rna(sequence)
    if not seq:
        return "", 0.0, "empty"

    # 1) LinearFold, if available.
    linear_path = "./LinearFold/bin/linearfold_c"
    if os.path.exists(linear_path) and os.access(linear_path, os.X_OK):
        try:
            cmd = [linear_path, "100", "0", "0", "0", "0", "0", "5.0", "", "0", "2"]
            process = subprocess.Popen(
                cmd,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            out, err = process.communicate(input=f"{seq}\n", timeout=20)
            toks = out.split()
            if len(toks) >= 3:
                dot = toks[1]
                mfe = float(toks[2].strip("()"))
                if len(dot) == len(seq):
                    return dot, mfe, "LinearFold"
            # If LinearFold exists but output was invalid, continue to ViennaRNA.
        except Exception:
            # Continue to ViennaRNA if LinearFold fails.
            pass

    # 2) ViennaRNA. Required if LinearFold is not available/usable.
    try:
        import RNA  # type: ignore
    except Exception as exc:
        raise RuntimeError(
            "No RNA folding backend available. Install ViennaRNA or provide "
            "./LinearFold/bin/linearfold_c."
        ) from exc

    dot, mfe = RNA.fold(seq)
    if len(dot) != len(seq):
        raise RuntimeError("ViennaRNA returned a structure with unexpected length.")
    return dot, float(mfe), "ViennaRNA"


# -------------------------
# Secondary structure helpers
# -------------------------

def pair_table_from_dotbracket(dotbracket: str) -> List[int]:
    pair_table = [-1] * len(dotbracket)
    stack: List[int] = []
    for i, ch in enumerate(dotbracket):
        if ch == "(":
            stack.append(i)
        elif ch == ")" and stack:
            j = stack.pop()
            pair_table[i] = j
            pair_table[j] = i
    return pair_table


def accessibility_score(dotbracket: str, start: int, end: int) -> float:
    if not dotbracket or start >= end:
        return 0.0
    start = max(0, start)
    end = min(len(dotbracket), end)
    window = dotbracket[start:end]
    if not window:
        return 0.0
    return window.count(".") / len(window)


def rbs_long_range_interactions(
    dotbracket: str,
    full_seq: str,
    rbs_start: int,
    rbs_end: int,
    threshold: int = LONG_RANGE_THRESHOLD_NT,
) -> Tuple[bool, str]:
    """
    Return True if an RBS base is paired with a far-away mRNA region.
    """
    pair_table = pair_table_from_dotbracket(dotbracket)
    details = []
    for i in range(rbs_start, rbs_end):
        j = pair_table[i] if 0 <= i < len(pair_table) else -1
        if j != -1 and abs(j - i) > threshold:
            details.append(f"{i+1}{full_seq[i]}-{j+1}{full_seq[j]} dist={abs(j-i)}")
    return (len(details) > 0), "; ".join(details)


# -------------------------
# Energetics
# -------------------------

def pair_energy(rbs_base: str, asd_base: str) -> float:
    """
    Lightweight energy-like base-pair scoring.
    Negative = favorable.
    """
    a = normalize_rna(rbs_base)
    b = normalize_rna(asd_base)
    if not a or not b:
        return 1.5
    pair = (a[0], b[0])
    if pair in [("G", "C"), ("C", "G")]:
        return -3.0
    if pair in [("A", "U"), ("U", "A")]:
        return -2.0
    if pair in [("G", "U"), ("U", "G")]:
        return -0.8
    return 1.2


def dG_start_codon(codon: str) -> float:
    if not USE_START_LOOKUP:
        return 0.0
    codon = normalize_rna(codon)[:3]
    return float(START_CODON_DG.get(codon, 0.0))


# -------------------------
# dG standby term
# -------------------------
# This is your friend's ViennaRNA-based standby calculation, wrapped so the
# GA still runs if ViennaRNA is not installed. It returns <= 0 kcal/mol.
# If no SD:aSD pairing is detected, standby is irrelevant and returns 0.

STANDBY_LEN = 4  # nt of standby site immediately 5' of SD binding region


def dG_standby(mRNA_upstream: str, anti_sd: str) -> float:
    """
    Standby contribution in kcal/mol.

    Inputs:
        mRNA_upstream: 5' UTR sequence up to but not including AUG
        anti_sd: anti-Shine-Dalgarno sequence, 5'->3'

    Returns:
        dG_standby <= 0.0. If no SD:aSD interaction is found, returns 0.0.
        Requires ViennaRNA; no approximate fallback is used.
    """
    try:
        import RNA  # type: ignore
    except Exception as exc:
        raise RuntimeError("dG_standby requires ViennaRNA. Install with: pip install ViennaRNA") from exc

    mRNA = normalize_rna(mRNA_upstream)
    rRNA = normalize_rna(anti_sd)
    if not mRNA or not rRNA:
        return 0.0

    md = RNA.md()
    m_len = len(mRNA)

    def _fold_alone(seq: str) -> Tuple[float, str]:
        if not seq:
            return 0.0, ""
        fc = RNA.fold_compound(seq, md)
        ss, dG = fc.mfe()
        return float(dG), ss

    def _parse_pairs(ss: str) -> List[Tuple[int, int]]:
        stack: List[int] = []
        pairs: List[Tuple[int, int]] = []
        for i, c in enumerate(ss):
            if c == "(":
                stack.append(i)
            elif c == ")" and stack:
                j = stack.pop()
                pairs.append((j, i))
        return pairs

    def _pairs_to_dotbracket(n: int, bp_x: List[int], bp_y: List[int]) -> str:
        ss = ["."] * n
        for x, y in zip(bp_x, bp_y):
            if 0 <= x < n and 0 <= y < n:
                ss[x] = "("
                ss[y] = ")"
        return "".join(ss)

    def _cofold(seq1: str, seq2: str) -> Tuple[str, float]:
        fc = RNA.fold_compound(seq1 + "&" + seq2, md)
        ss, dG = fc.mfe()
        # ViennaRNA may include &, but the original friend's code expects no separator.
        ss = ss.replace("&", "")
        return ss, float(dG)

    def _eval_structure(mRNA_seq: str, rRNA_seq: str, bp_x: List[int], bp_y: List[int]) -> float:
        n = len(mRNA_seq) + len(rRNA_seq)
        ss = _pairs_to_dotbracket(n, bp_x, bp_y)
        fc = RNA.fold_compound(mRNA_seq + "&" + rRNA_seq, md)
        return float(fc.eval_structure(ss))

    # Step 1: cofold mRNA upstream with anti-SD.
    ss_co, energy_before = _cofold(mRNA, rRNA)

    all_pairs = _parse_pairs(ss_co)
    bp_x = [p[0] for p in all_pairs]
    bp_y = [p[1] for p in all_pairs]

    # Step 2: find mRNA:rRNA intermolecular pairs.
    mRNA_rRNA_pairs: List[Tuple[int, int]] = []
    for x, y in all_pairs:
        if x < m_len <= y:
            mRNA_rRNA_pairs.append((x, y))
        elif y < m_len <= x:
            mRNA_rRNA_pairs.append((y, x))

    if not mRNA_rRNA_pairs:
        return 0.0

    # Left edge of SD:aSD interaction on mRNA.
    most_5p_mRNA = min(p[0] for p in mRNA_rRNA_pairs)

    # Step 3: define pre-standby region.
    sb_start = max(0, most_5p_mRNA - STANDBY_LEN - 1)
    pre_standby_seq = mRNA[0:sb_start]

    # Step 4: keep all SD+rRNA pairs downstream of binding start.
    bp_x_3p = [x for x, y in zip(bp_x, bp_y) if x >= most_5p_mRNA]
    bp_y_3p = [y for x, y in zip(bp_x, bp_y) if x >= most_5p_mRNA]

    # Step 5: fold the pre-standby region alone.
    _, ss_5p = _fold_alone(pre_standby_seq)
    pairs_5p = _parse_pairs(ss_5p)
    bp_x_5p = [p[0] for p in pairs_5p]
    bp_y_5p = [p[1] for p in pairs_5p]

    # Step 6: evaluate structure where standby bases are left unpaired.
    bp_x_after = bp_x_5p + bp_x_3p
    bp_y_after = bp_y_5p + bp_y_3p

    energy_after = _eval_structure(mRNA, rRNA, bp_x_after, bp_y_after)

    # Free standby site is not a bonus beyond zero.
    return float(min(0.0, energy_before - energy_after))


def dG_spacing_penalty(aligned_spacing: float, optimal_spacing: float = OPTIMAL_SPACING) -> Tuple[float, float]:
    """
    Salis-style spacing penalty.

    aligned_spacing = (AUG index - mRNA binding start) - aSD binding start from 5'
    d = aligned_spacing - optimal_spacing
    """
    d = float(aligned_spacing) - float(optimal_spacing)
    if d >= 0:
        penalty = 0.048 * d * d + 0.24 * d
    else:
        penalty = 12.2 / ((1.0 + math.exp(2.5 * (d + 2.0))) ** 3.0)
    return d, float(penalty)


def calculate_tir(dG_total: float, r0: float = R0, rt: float = RT_PHYSICAL) -> float:
    """
    TIR = r0 * exp(-dG_total / RT)
    Uses physical RT by default as requested.
    """
    if math.isnan(dG_total):
        return 0.0
    arg = -float(dG_total) / float(rt)
    # avoid overflow
    arg = max(min(arg, 700), -700)
    return float(r0 * math.exp(arg))


def dG_mrna_unfolding_site(full_seq: str, site_start: int, site_end: int, window: int = 35) -> float:
    """
    Positive local unfolding penalty around a binding site.
    Approximation: penalty = -MFE(local region), clipped at >= 0.
    """
    full_seq = normalize_rna(full_seq)
    left = max(0, int(site_start) - window)
    right = min(len(full_seq), int(site_end) + window)
    local_region = full_seq[left:right]
    _dot, mfe, _backend = fold_sequence(local_region)
    return max(0.0, -float(mfe))


def cofold_deltaG(utr: str, anti_sd: str) -> float:
    """
    Candidate-level duplex term:
      dG_binding = dG_complex - dG_utr - dG_asd

    Requires ViennaRNA. No approximate fallback is used.
    """
    try:
        import RNA  # type: ignore
    except Exception as exc:
        raise RuntimeError("cofold_deltaG requires ViennaRNA. Install with: pip install ViennaRNA") from exc

    utr = normalize_rna(utr)
    anti_sd = normalize_rna(anti_sd)
    if not utr or not anti_sd:
        return 0.0

    struct_complex, mfe_complex = RNA.cofold(f"{utr}&{anti_sd}")
    struct_utr, mfe_utr = RNA.fold(utr)
    struct_asd, mfe_asd = RNA.fold(anti_sd)
    return float(mfe_complex - mfe_utr - mfe_asd)


def delta_g_duplex_whole_utr(utr: str, anti_sd: str) -> Tuple[float, str]:
    """
    Calculate dG duplex for the whole 5' UTR with anti-SD.
    Uses ViennaRNA cofold only. No approximate scan fallback.
    """
    return cofold_deltaG(utr, anti_sd), "ViennaRNA cofold"


def normalized_affinity_from_dG(dg: float, scale: float = 12.0) -> float:
    """Dashboard-friendly 0-1 affinity; more negative dG = higher score."""
    return clamp01((-dg) / scale)


# -------------------------
# Candidate construction and binding-site scanning
# -------------------------

def canonicalize_candidate(
    candidate: Dict[str, Any],
    default_five_prime_utr: str = DEFAULT_FIVE_PRIME_UTR,
    default_standby: str = DEFAULT_STANDBY,
    default_cds_start: str = DEFAULT_CDS_START,
) -> Dict[str, Any]:
    """Return a candidate in the canonical DalGuardYES sequence schema.

    Canonical sequence order:
        five_prime_utr + standby + rbs_left + rbs_core + rbs_right + spacer + cds_start

    Older ``five_prime_flank`` / ``core`` inputs are accepted only as a
    compatibility bridge.  New GA populations always retain the explicit
    canonical fields.
    """
    c = dict(candidate or {})

    # New schema first; fall back to the legacy flank only when necessary.
    if c.get("five_prime_utr") is not None:
        five_prime_utr = normalize_rna(c.get("five_prime_utr", ""))
    elif c.get("five_prime_flank") is not None:
        five_prime_utr = normalize_rna(c.get("five_prime_flank", ""))
    else:
        five_prime_utr = normalize_rna(default_five_prime_utr)

    standby = normalize_rna(c.get("standby", default_standby))
    spacer = normalize_rna(c.get("spacer", ""))
    cds_start = normalize_rna(c.get("cds_start", default_cds_start))
    if not cds_start.startswith("AUG"):
        cds_start = "AUG" + cds_start

    rbs = normalize_rna(c.get("rbs", ""))
    rbs_core = normalize_rna(c.get("rbs_core", c.get("core", "")))
    rbs_left = normalize_rna(c.get("rbs_left", ""))
    rbs_right = normalize_rna(c.get("rbs_right", ""))

    # Explicit left/right fields mean the canonical components are the source
    # of truth.  A legacy input may provide only ``core`` plus a full ``rbs``;
    # in that case we must split the full RBS rather than collapse it to core.
    has_explicit_sides = ("rbs_left" in c) or ("rbs_right" in c)
    if has_explicit_sides:
        component_rbs = rbs_left + rbs_core + rbs_right
        if component_rbs:
            rbs = component_rbs
    elif rbs and rbs_core and rbs_core in rbs:
        idx = rbs.find(rbs_core)
        rbs_left = rbs[:idx]
        rbs_right = rbs[idx + len(rbs_core):]
    elif not rbs and rbs_core:
        rbs = rbs_core
    # Last-resort compatibility: protect the full supplied RBS.
    elif rbs and not rbs_core:
        rbs_core = rbs

    out: Dict[str, Any] = dict(c)
    out.update({
        "five_prime_utr": five_prime_utr,
        "standby": standby,
        "rbs_left": rbs_left,
        "rbs_core": rbs_core,
        "rbs_right": rbs_right,
        "rbs": rbs,
        "spacer": spacer,
        "cds_start": cds_start,
        "mutable_regions": c.get(
            "mutable_regions",
            ["standby", "rbs_left", "rbs_right", "spacer"],
        ),
        "protected_regions": c.get(
            "protected_regions",
            ["five_prime_utr", "rbs_core", "cds_start"],
        ),
    })
    return out


def build_full_sequence(
    candidate: Dict[str, Any],
    default_five_prime_utr: str = DEFAULT_FIVE_PRIME_UTR,
    default_standby: str = DEFAULT_STANDBY,
    default_cds_start: str = DEFAULT_CDS_START,
) -> Tuple[str, str, str, str, int, int, int, int]:
    c = canonicalize_candidate(
        candidate,
        default_five_prime_utr=default_five_prime_utr,
        default_standby=default_standby,
        default_cds_start=default_cds_start,
    )
    prefix = c["five_prime_utr"] + c["standby"]
    rbs = c["rbs"]
    spacer = c["spacer"]
    cds = c["cds_start"]

    full_seq = prefix + rbs + spacer + cds
    rbs_start = len(prefix)
    rbs_end = rbs_start + len(rbs)
    aug_start = rbs_end + len(spacer)
    aug_end = aug_start + 3
    utr_seq = full_seq[:aug_start]

    return full_seq, utr_seq, rbs, spacer, rbs_start, rbs_end, aug_start, aug_end


def scan_binding_sites(
    utr: str,
    anti_sd: str,
    aug_start: int,
    min_overlap: int = 5,
    max_sites: int = 12,
) -> List[Dict[str, Any]]:
    """
    Scan possible anti-SD binding sites in the UTR.
    We use this for binding-site coordinates, spacing penalty, and unfolding windows.

    Note:
    - dG_duplex is treated as candidate-level in this project.
    - This scanner is primarily for likely binding coordinates.
    """
    utr = normalize_rna(utr)
    anti_sd = normalize_rna(anti_sd)
    anti_rev = anti_sd[::-1]
    raw: List[Dict[str, Any]] = []

    for mrna_start in range(len(utr)):
        for asd_start_rev in range(len(anti_rev)):
            energy = 0.0
            overlap = 0
            i, j = mrna_start, asd_start_rev

            while i < len(utr) and j < len(anti_rev):
                e = pair_energy(utr[i], anti_rev[j])
                energy += e
                overlap += 1
                i += 1
                j += 1

            if overlap < min_overlap:
                continue

            # Favorable pairing score: more favorable = higher positive score.
            pairing_score = max(0.0, -energy)
            if pairing_score <= 0:
                continue

            # Convert reversed anti-SD scan coordinates to original 5' coordinates.
            # Reversed index block: [asd_start_rev, asd_start_rev+overlap)
            # Original 5' index block:
            #   start = len(asd) - (asd_start_rev + overlap)
            #   end   = len(asd) - asd_start_rev
            asd_start_5p = len(anti_sd) - (asd_start_rev + overlap)
            asd_end_5p_exclusive = len(anti_sd) - asd_start_rev
            asd_start_5p = max(0, asd_start_5p)
            asd_end_5p_exclusive = min(len(anti_sd), asd_end_5p_exclusive)

            aligned_spacing = (aug_start - mrna_start) - asd_start_5p
            d, spacing = dG_spacing_penalty(aligned_spacing)

            raw.append({
                "mrna_start": mrna_start,
                "mrna_end": i,
                "asd_start_5p": asd_start_5p,
                "asd_end_5p_exclusive": asd_end_5p_exclusive,
                "overlap": overlap,
                "pairing_score": pairing_score,
                "aligned_spacing": aligned_spacing,
                "d": d,
                "dG_spacing": spacing,
            })

    # Sort by pairing quality, then spacing penalty; remove near-duplicates.
    raw.sort(key=lambda x: (x["pairing_score"], -x["dG_spacing"]), reverse=True)

    filtered: List[Dict[str, Any]] = []
    used_starts: List[int] = []
    for site in raw:
        if all(abs(site["mrna_start"] - s) > 2 for s in used_starts):
            filtered.append(site)
            used_starts.append(site["mrna_start"])
        if len(filtered) >= max_sites:
            break

    return filtered


# -------------------------
# Evaluation
# -------------------------


def score_sites_for_asd(
    candidate_id: str,
    anti_sd_type: str,
    full_seq: str,
    utr_seq: str,
    anti_sd: str,
    aug_start: int,
    dG_duplex_candidate: float,
    dG_start: float,
    dG_standby: float,
    collect_sites: bool = True,
) -> Tuple[float, float, float, float, List[BindingSiteResult]]:
    """
    Score all possible binding sites for one anti-SD.

    Candidate-level constants:
      - dG_duplex_candidate
      - dG_start
      - dG_standby

    Site-specific terms:
      - dG_spacing
      - dG_mrna_unfolding

    Returns:
      best_tir, best_dG_total, best_dG_spacing, best_dG_unfolding, site_results
    """
    sites = scan_binding_sites(utr_seq, anti_sd, aug_start, min_overlap=5, max_sites=12)
    if not sites:
        return 0.0, 999.0, float("nan"), float("nan"), []

    rows: List[Dict[str, Any]] = []
    for site in sites:
        dG_unfold = dG_mrna_unfolding_site(
            full_seq=full_seq,
            site_start=site["mrna_start"],
            site_end=site["mrna_end"],
            window=35,
        )

        dG_total_site = (
            dG_duplex_candidate
            + dG_start
            + dG_standby
            + site["dG_spacing"]
            + dG_unfold
        )
        tir_site = calculate_tir(dG_total_site)

        rows.append({
            **site,
            "dG_mrna_unfolding": dG_unfold,
            "dG_total": dG_total_site,
            "tir": tir_site,
        })

    rows.sort(key=lambda x: x["tir"], reverse=True)
    best = rows[0]

    site_results: List[BindingSiteResult] = []
    if collect_sites:
        for rank, site in enumerate(rows, 1):
            site_results.append(BindingSiteResult(
                candidate_id=candidate_id,
                anti_sd_type=anti_sd_type,
                site_rank=rank,
                mrna_start_0based=site["mrna_start"],
                mrna_end_0based_exclusive=site["mrna_end"],
                mrna_start_1based=site["mrna_start"] + 1,
                mrna_end_1based=site["mrna_end"],
                asd_start_5p_0based=site["asd_start_5p"],
                asd_start_5p_1based=site["asd_start_5p"] + 1,
                asd_end_5p_0based_exclusive=site["asd_end_5p_exclusive"],
                asd_end_5p_1based=site["asd_end_5p_exclusive"],
                overlap=site["overlap"],
                pairing_score=round(float(site["pairing_score"]), 4),
                aligned_spacing=round(float(site["aligned_spacing"]), 4),
                d=round(float(site["d"]), 4),
                dG_spacing=round(float(site["dG_spacing"]), 4),
                dG_mrna_unfolding=round(float(site["dG_mrna_unfolding"]), 4),
                dG_duplex_candidate=round(float(dG_duplex_candidate), 4),
                dG_start=round(float(dG_start), 4),
                dG_standby=round(float(dG_standby), 4),
                dG_total=round(float(site["dG_total"]), 4),
                tir=float(site["tir"]),
                is_best_site=(rank == 1),
            ))

    return (
        float(best["tir"]),
        float(best["dG_total"]),
        float(best["dG_spacing"]),
        float(best["dG_mrna_unfolding"]),
        site_results,
    )


def evaluate_candidate(
    candidate: Dict[str, Any],
    orth_anti_sd: str,
    wt_anti_sd: str,
    default_five_prime_utr: str = DEFAULT_FIVE_PRIME_UTR,
    default_standby: str = DEFAULT_STANDBY,
    default_cds_start: str = DEFAULT_CDS_START,
    wt_penalty_constant: float = DEFAULT_WT_PENALTY_CONSTANT,
    candidate_id: Optional[str] = None,
    default_flank: Optional[str] = None,
) -> Tuple[CandidateEval, List[BindingSiteResult]]:
    """Evaluate one canonical candidate without changing the thermodynamic model."""
    if default_flank is not None:
        default_five_prime_utr = default_flank

    cid = str(candidate_id or candidate.get("id", f"cand_{random.randrange(10**9)}"))
    c = canonicalize_candidate(
        candidate,
        default_five_prime_utr=default_five_prime_utr,
        default_standby=default_standby,
        default_cds_start=default_cds_start,
    )

    full_seq, utr_seq, rbs, spacer, rbs_start, rbs_end, aug_start, aug_end = build_full_sequence(
        c,
        default_five_prime_utr=default_five_prime_utr,
        default_standby=default_standby,
        default_cds_start=default_cds_start,
    )

    structure, mfe, backend = fold_sequence(full_seq)
    rbs_access = accessibility_score(structure, rbs_start, rbs_end)
    aug_access = accessibility_score(structure, max(0, aug_start - 3), min(len(full_seq), aug_end + 6))

    long_range_flag, long_range_pairs = rbs_long_range_interactions(
        structure, full_seq, rbs_start, rbs_end, LONG_RANGE_THRESHOLD_NT
    )

    dG_duplex_orth, duplex_backend_orth = delta_g_duplex_whole_utr(utr_seq, orth_anti_sd)
    dG_duplex_wt, _duplex_backend_wt = delta_g_duplex_whole_utr(utr_seq, wt_anti_sd)
    dG_start = dG_start_codon(full_seq[aug_start:aug_end])

    # Keep the existing standby-energy model unchanged; retain both orthogonal
    # and WT values instead of discarding the WT value.
    dG_standby_orth = dG_standby(utr_seq, orth_anti_sd)
    dG_standby_wt = dG_standby(utr_seq, wt_anti_sd)

    orth_tir, dG_total_orth, orth_spacing, orth_unfold, orth_sites = score_sites_for_asd(
        candidate_id=cid,
        anti_sd_type="orthogonal",
        full_seq=full_seq,
        utr_seq=utr_seq,
        anti_sd=orth_anti_sd,
        aug_start=aug_start,
        dG_duplex_candidate=dG_duplex_orth,
        dG_start=dG_start,
        dG_standby=dG_standby_orth,
        collect_sites=True,
    )

    wt_tir, dG_total_wt, wt_spacing, wt_unfold, wt_sites = score_sites_for_asd(
        candidate_id=cid,
        anti_sd_type="wt",
        full_seq=full_seq,
        utr_seq=utr_seq,
        anti_sd=wt_anti_sd,
        aug_start=aug_start,
        dG_duplex_candidate=dG_duplex_wt,
        dG_start=dG_start,
        dG_standby=dG_standby_wt,
        collect_sites=True,
    )

    site_results = orth_sites + wt_sites
    t_score = safe_log10(orth_tir) - wt_penalty_constant * safe_log10(wt_tir)
    fitness_for_selection = signed_log10_score(t_score)

    # Existing hard filter retained exactly for selection/ranking.
    if long_range_flag or not orth_sites:
        fitness_for_selection = -1e9
        t_score = -1e9

    orthScore = normalized_affinity_from_dG(dG_duplex_orth)
    wtLeakage = normalized_affinity_from_dG(dG_duplex_wt)

    def _round_or_nan(x: float) -> float:
        return round(float(x), 4) if not math.isnan(float(x)) else float("nan")

    eval_obj = CandidateEval(
        candidate_id=cid,
        five_prime_utr=c["five_prime_utr"],
        standby=c["standby"],
        rbs_left=c["rbs_left"],
        rbs_core=c["rbs_core"],
        rbs_right=c["rbs_right"],
        rbs=rbs,
        spacer=spacer,
        cds_start=c["cds_start"],
        full_seq=full_seq,
        utr_seq=utr_seq,
        structure=structure,
        mfe=round(float(mfe), 4),
        backend=f"{backend}; duplex={duplex_backend_orth}",
        rbs_start=rbs_start,
        rbs_end=rbs_end,
        aug_start=aug_start,
        aug_end=aug_end,
        rbs_access=round(float(rbs_access), 4),
        aug_access=round(float(aug_access), 4),
        long_range_flag=bool(long_range_flag),
        long_range_pairs=long_range_pairs,
        dG_duplex_orth=round(float(dG_duplex_orth), 4),
        dG_duplex_wt=round(float(dG_duplex_wt), 4),
        dG_start=round(float(dG_start), 4),
        dG_standby_orth=round(float(dG_standby_orth), 4),
        dG_standby_wt=round(float(dG_standby_wt), 4),
        best_dG_spacing_orth=_round_or_nan(orth_spacing),
        best_dG_mrna_unfolding_orth=_round_or_nan(orth_unfold),
        dG_total_orth=round(float(dG_total_orth), 4),
        best_dG_spacing_wt=_round_or_nan(wt_spacing),
        best_dG_mrna_unfolding_wt=_round_or_nan(wt_unfold),
        dG_total_wt=round(float(dG_total_wt), 4),
        orth_tir=float(orth_tir),
        wt_tir=float(wt_tir),
        t_score=float(t_score),
        fitness_for_selection=round(float(fitness_for_selection), 6),
        dG_standby=round(float(dG_standby_orth), 4),
        best_dG_spacing=_round_or_nan(orth_spacing),
        best_dG_mrna_unfolding=_round_or_nan(orth_unfold),
        dG_total=round(float(dG_total_orth), 4),
        orthScore=round(float(orthScore), 4),
        wtLeakage=round(float(wtLeakage), 4),
        fitness=0.0,
    )
    return eval_obj, site_results


# -------------------------
# GA operators
# -------------------------

def _mutate_one_base(seq: str) -> str:
    if not seq:
        return seq
    chars = list(seq)
    i = random.randrange(len(chars))
    chars[i] = random.choice([b for b in BASES if b != chars[i]])
    return "".join(chars)


def mutate_candidate(candidate: Dict[str, Any], mutate_flank: bool = True) -> Dict[str, Any]:
    """Mutate only standby, non-core RBS segments, and spacer.

    ``mutate_flank`` is retained as a legacy argument but now controls standby
    mutation only.  The fixed 5′ UTR, protected RBS core, and CDS are never
    mutated.
    """
    new = canonicalize_candidate(candidate)
    mutable = set(new.get("mutable_regions") or ["standby", "rbs_left", "rbs_right", "spacer"])

    ops: List[str] = []
    if mutate_flank and "standby" in mutable and new["standby"]:
        ops.extend(["standby_sub", "standby_sub"])
    if ("rbs_left" in mutable or "rbs" in mutable) and new["rbs_left"]:
        ops.extend(["rbs_left_sub", "rbs_left_sub"])
    if ("rbs_right" in mutable or "rbs" in mutable) and new["rbs_right"]:
        ops.extend(["rbs_right_sub", "rbs_right_sub"])
    if "spacer" in mutable:
        ops.extend(["spacer_sub", "spacer_sub", "spacer_insert", "spacer_delete"])

    if not ops:
        return new

    op = random.choice(ops)
    spacer = new["spacer"]
    spacer_min = int(new.get("spacer_len_min", 3) or 3)
    spacer_max = int(new.get("spacer_len_max", 14) or 14)

    if op == "standby_sub":
        new["standby"] = _mutate_one_base(new["standby"])
    elif op == "rbs_left_sub":
        new["rbs_left"] = _mutate_one_base(new["rbs_left"])
    elif op == "rbs_right_sub":
        new["rbs_right"] = _mutate_one_base(new["rbs_right"])
    elif op == "spacer_sub" and spacer:
        new["spacer"] = _mutate_one_base(spacer)
    elif op == "spacer_insert" and len(spacer) < spacer_max:
        i = random.randrange(len(spacer) + 1)
        new["spacer"] = spacer[:i] + random.choice(BASES) + spacer[i:]
    elif op == "spacer_delete" and len(spacer) > spacer_min:
        i = random.randrange(len(spacer))
        new["spacer"] = spacer[:i] + spacer[i + 1:]

    new["rbs"] = new["rbs_left"] + new["rbs_core"] + new["rbs_right"]
    return new


def crossover(a: Dict[str, Any], b: Dict[str, Any]) -> Dict[str, Any]:
    """Modular crossover that preserves all protected sequence fields."""
    ca = canonicalize_candidate(a)
    cb = canonicalize_candidate(b)
    child: Dict[str, Any] = {
        "five_prime_utr": ca["five_prime_utr"],
        "standby": random.choice([ca["standby"], cb["standby"]]),
        "rbs_left": random.choice([ca["rbs_left"], cb["rbs_left"]]),
        "rbs_core": ca["rbs_core"],
        "rbs_right": random.choice([ca["rbs_right"], cb["rbs_right"]]),
        "spacer": random.choice([ca["spacer"], cb["spacer"]]),
        "cds_start": ca["cds_start"],
        "mutable_regions": ca.get("mutable_regions", ["standby", "rbs_left", "rbs_right", "spacer"]),
        "protected_regions": ca.get("protected_regions", ["five_prime_utr", "rbs_core", "cds_start"]),
        "spacer_len_min": ca.get("spacer_len_min", cb.get("spacer_len_min", 3)),
        "spacer_len_max": ca.get("spacer_len_max", cb.get("spacer_len_max", 14)),
        "source": "ga_crossover",
    }
    child["rbs"] = child["rbs_left"] + child["rbs_core"] + child["rbs_right"]
    return child


def generate_guided_seed_candidates(
    orth_anti_sd: str,
    n: int,
    default_five_prime_utr: str = DEFAULT_FIVE_PRIME_UTR,
    default_standby: str = DEFAULT_STANDBY,
    default_cds_start: str = DEFAULT_CDS_START,
    default_flank: Optional[str] = None,
) -> List[Dict[str, Any]]:
    """Fallback guided seeds when the supplied initial population is too small."""
    if default_flank is not None:
        default_five_prime_utr = default_flank
    core = reverse_complement(orth_anti_sd) or "UACAAG"

    def random_spacer(length: int) -> str:
        return "".join(random.choice(BASES) for _ in range(length))

    out: List[Dict[str, Any]] = []
    while len(out) < n:
        spacer = random_spacer(random.randint(4, 12))
        out.append({
            "five_prime_utr": normalize_rna(default_five_prime_utr),
            "standby": normalize_rna(default_standby),
            "rbs_left": "",
            "rbs_core": core,
            "rbs_right": "",
            "rbs": core,
            "spacer": spacer,
            "cds_start": normalize_rna(default_cds_start),
            "mutable_regions": ["standby", "rbs_left", "rbs_right", "spacer"],
            "protected_regions": ["five_prime_utr", "rbs_core", "cds_start"],
            "spacer_len_min": 4,
            "spacer_len_max": 12,
            "source": "guided_seed",
        })
    return out


# -------------------------
# GA main loop
# -------------------------

def _eval_worker(args: Tuple) -> Tuple[CandidateEval, List[BindingSiteResult]]:
    """Module-level so ProcessPoolExecutor can pickle it."""
    (
        candidate, orth_anti_sd, wt_anti_sd, default_five_prime_utr,
        default_standby, default_cds_start, wt_penalty_constant, cid,
    ) = args
    return evaluate_candidate(
        candidate,
        orth_anti_sd=orth_anti_sd,
        wt_anti_sd=wt_anti_sd,
        default_five_prime_utr=default_five_prime_utr,
        default_standby=default_standby,
        default_cds_start=default_cds_start,
        wt_penalty_constant=wt_penalty_constant,
        candidate_id=cid,
    )


def run_ga(
    initial_candidates: List[Dict[str, Any]],
    orth_anti_sd: str = DEFAULT_ORTH_ASD,
    wt_anti_sd: str = DEFAULT_WT_ASD,
    default_five_prime_utr: str = DEFAULT_FIVE_PRIME_UTR,
    default_standby: str = DEFAULT_STANDBY,
    default_cds_start: str = DEFAULT_CDS_START,
    generations: int = 30,
    population_size: int = 80,
    elite_fraction: float = 0.20,
    wt_penalty_constant: float = DEFAULT_WT_PENALTY_CONSTANT,
    seed: int = 7,
    parallel: bool = False,
    max_workers: Optional[int] = None,
    default_flank: Optional[str] = None,
) -> Tuple[List[CandidateEval], List[Dict[str, Any]], List[BindingSiteResult], List[Dict[str, Any]]]:
    """Run the GA and retain scientifically meaningful per-generation statistics."""
    if default_flank is not None:
        default_five_prime_utr = default_flank

    random.seed(seed)
    orth_anti_sd = normalize_rna(orth_anti_sd)
    wt_anti_sd = normalize_rna(wt_anti_sd)

    population: List[Dict[str, Any]] = []
    for i, c in enumerate(initial_candidates or []):
        cc = canonicalize_candidate(
            c,
            default_five_prime_utr=default_five_prime_utr,
            default_standby=default_standby,
            default_cds_start=default_cds_start,
        )
        cc["id"] = str(c.get("id", f"seed_{i+1:04d}"))
        cc["source"] = c.get("source", "initial_seed")
        population.append(cc)

    if len(population) < population_size:
        population.extend(generate_guided_seed_candidates(
            orth_anti_sd,
            population_size - len(population),
            default_five_prime_utr=default_five_prime_utr,
            default_standby=default_standby,
            default_cds_start=default_cds_start,
        ))

    random.shuffle(population)
    population = population[:population_size]

    all_evals_by_key: Dict[Tuple[str, str, str, str], CandidateEval] = {}
    all_binding_sites: List[BindingSiteResult] = []
    history: List[Dict[str, Any]] = []
    n_workers = max_workers or max(1, (os.cpu_count() or 2) - 1)

    for gen in range(generations):
        scored: List[Tuple[CandidateEval, Dict[str, Any]]] = []

        if parallel and len(population) > 1:
            tasks = [
                (
                    ind, orth_anti_sd, wt_anti_sd, default_five_prime_utr,
                    default_standby, default_cds_start, wt_penalty_constant,
                    ind.get("id", f"g{gen}_i{idx}"),
                )
                for idx, ind in enumerate(population)
            ]
            with ProcessPoolExecutor(max_workers=n_workers, mp_context=_FORK_CTX) as executor:
                results = list(executor.map(_eval_worker, tasks))
            pairs = zip(population, results)
            for ind, (ev, sites) in pairs:
                scored.append((ev, ind))
                all_binding_sites.extend(sites)
                key = (ev.five_prime_utr, ev.standby, ev.rbs, ev.spacer)
                if key not in all_evals_by_key or ev.fitness_for_selection > all_evals_by_key[key].fitness_for_selection:
                    all_evals_by_key[key] = ev
        else:
            for idx, ind in enumerate(population):
                cid = ind.get("id", f"g{gen}_i{idx}")
                ev, sites = evaluate_candidate(
                    ind,
                    orth_anti_sd=orth_anti_sd,
                    wt_anti_sd=wt_anti_sd,
                    default_five_prime_utr=default_five_prime_utr,
                    default_standby=default_standby,
                    default_cds_start=default_cds_start,
                    wt_penalty_constant=wt_penalty_constant,
                    candidate_id=cid,
                )
                scored.append((ev, ind))
                all_binding_sites.extend(sites)
                key = (ev.five_prime_utr, ev.standby, ev.rbs, ev.spacer)
                if key not in all_evals_by_key or ev.fitness_for_selection > all_evals_by_key[key].fitness_for_selection:
                    all_evals_by_key[key] = ev

        scored.sort(key=lambda x: x[0].fitness_for_selection, reverse=True)
        evals = [x[0] for x in scored]
        best = evals[0]
        valid = [
            e for e in evals
            if e.fitness_for_selection > -1e8
            and not math.isnan(e.fitness_for_selection)
            and not math.isinf(e.fitness_for_selection)
        ]

        if valid:
            best_valid = max(valid, key=lambda e: e.fitness_for_selection)
            mean_valid_t_score = sum(e.t_score for e in valid) / len(valid)
            mean_valid_selection = sum(e.fitness_for_selection for e in valid) / len(valid)
            mean_valid_orth_tir = sum(e.orth_tir for e in valid) / len(valid)
            mean_valid_wt_tir = sum(e.wt_tir for e in valid) / len(valid)
            best_t_score = best_valid.t_score
            best_selection = best_valid.fitness_for_selection
            best_candidate_orth_tir = best_valid.orth_tir
            best_candidate_wt_tir = best_valid.wt_tir
            best_rbs_access = best_valid.rbs_access
        else:
            best_t_score = float("nan")
            mean_valid_t_score = float("nan")
            best_selection = best.fitness_for_selection
            mean_valid_selection = float("nan")
            best_candidate_orth_tir = float("nan")
            mean_valid_orth_tir = float("nan")
            best_candidate_wt_tir = float("nan")
            mean_valid_wt_tir = float("nan")
            best_rbs_access = float("nan")

        history.append({
            "generation": gen,
            "best_t_score": best_t_score,
            "mean_valid_t_score": mean_valid_t_score,
            "best_selection_score": best_selection,
            "mean_valid_selection_score": mean_valid_selection,
            "best_candidate_orth_tir": best_candidate_orth_tir,
            "mean_valid_orth_tir": mean_valid_orth_tir,
            "best_candidate_wt_tir": best_candidate_wt_tir,
            "mean_valid_wt_tir": mean_valid_wt_tir,
            "best_rbs_access": best_rbs_access,
            "valid_candidates": len(valid),
            "filtered_candidates": len(evals) - len(valid),
        })

        elite_count = max(1, min(len(scored), max(4, int(population_size * elite_fraction))))
        elites = [ind for _ev, ind in scored[:elite_count]]

        next_pop: List[Dict[str, Any]] = []
        for i, e in enumerate(elites):
            kept = dict(e)
            kept["id"] = f"g{gen+1}_elite_{i}"
            next_pop.append(kept)

        while len(next_pop) < population_size:
            if random.random() < 0.35 and len(elites) >= 2:
                p1, p2 = random.sample(elites, 2)
                child = crossover(p1, p2)
            else:
                child = dict(random.choice(elites))
            child = mutate_candidate(child, mutate_flank=True)
            child["id"] = f"g{gen+1}_child_{len(next_pop)}"
            next_pop.append(child)
        population = next_pop

    final_evals = list(all_evals_by_key.values())
    final_evals.sort(key=lambda e: e.fitness_for_selection, reverse=True)

    # Legacy dashboard-only normalized candidate fitness.  Valid scientific
    # outputs use T-score / selection score and are exported separately.
    finite_scores = [
        e.fitness_for_selection for e in final_evals
        if e.fitness_for_selection > -1e8
        and not math.isinf(e.fitness_for_selection)
        and not math.isnan(e.fitness_for_selection)
    ]
    lo, hi = (min(finite_scores), max(finite_scores)) if finite_scores else (0.0, 1.0)
    for e in final_evals:
        if e.fitness_for_selection <= -1e8:
            e.fitness = 0.0
        elif hi > lo:
            e.fitness = round((e.fitness_for_selection - lo) / (hi - lo), 4)
        else:
            e.fitness = 1.0

    # Legacy history aliases retained without the old -1e9-contaminated average.
    valid_hist_scores = [
        v for h in history
        for v in (h.get("best_selection_score"), h.get("mean_valid_selection_score"))
        if isinstance(v, (int, float)) and not math.isnan(v) and not math.isinf(v)
    ]
    hmin, hmax = (min(valid_hist_scores), max(valid_hist_scores)) if valid_hist_scores else (0.0, 1.0)
    for h in history:
        h["bestRaw"] = h.get("best_selection_score")
        h["avgRaw"] = h.get("mean_valid_selection_score")
        h["bestTScore"] = h.get("best_t_score")
        h["avgValidTScore"] = h.get("mean_valid_t_score")
        h["bestTIR"] = h.get("best_candidate_orth_tir")
        h["avgTIR"] = h.get("mean_valid_orth_tir")
        for source, target in (("best_selection_score", "best"), ("mean_valid_selection_score", "avg")):
            val = h.get(source)
            if not isinstance(val, (int, float)) or math.isnan(val) or math.isinf(val):
                h[target] = float("nan")
            elif hmax > hmin:
                h[target] = round((val - hmin) / (hmax - hmin), 4)
            else:
                h[target] = 1.0

    return final_evals, history, all_binding_sites, population


# -------------------------
# Input/output adapters
# -------------------------

def load_initial_candidates(path: Optional[str]) -> List[Dict[str, str]]:
    """
    Load friend's candidates from JSON or CSV.
    JSON can be:
      - list of dicts
      - {"candidates": [...]}
    CSV should have columns:
      rbs, spacer, optional five_prime_flank, cds_start, id, source
    """
    if not path:
        return []

    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Initial candidates file not found: {path}")

    if p.suffix.lower() == ".json":
        data = json.loads(p.read_text())
        if isinstance(data, list):
            return data
        if isinstance(data, dict):
            return data.get("candidates", [])
        raise ValueError("JSON input must be a list or an object with a 'candidates' key.")

    if p.suffix.lower() == ".csv":
        out: List[Dict[str, str]] = []
        with p.open(newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                out.append({k: v for k, v in row.items() if v is not None})
        return out

    raise ValueError("Unsupported input format. Use .json or .csv")


def candidate_to_dashboard_dict(e: CandidateEval) -> Dict[str, Any]:
    """Compact dashboard representation; the workbook uses full CandidateEval rows."""
    return {
        "candidateId": e.candidate_id,
        "fivePrimeUTR": e.five_prime_utr,
        "standby": e.standby,
        "rbsLeft": e.rbs_left,
        "rbsCore": e.rbs_core,
        "rbsRight": e.rbs_right,
        "rbs": e.rbs,
        "spacer": e.spacer,
        "cdsStart": e.cds_start,
        "fullSeq": e.full_seq,
        "utrSeq": e.utr_seq,
        "augStart": e.aug_start,
        "structure": e.structure,
        "rbsAccess": f"{e.rbs_access:.3f}",
        "orthScore": f"{e.orthScore:.3f}",
        "wtLeakage": f"{e.wtLeakage:.3f}",
        "fitness": f"{e.fitness:.3f}",
        "dGDuplexOrth": e.dG_duplex_orth,
        "dGDuplexWT": e.dG_duplex_wt,
        "dGStart": e.dG_start,
        "dGStandbyOrth": e.dG_standby_orth,
        "dGStandbyWT": e.dG_standby_wt,
        "dGSpacingOrth": e.best_dG_spacing_orth,
        "dGmRNAUnfoldingOrth": e.best_dG_mrna_unfolding_orth,
        "dGTotalOrth": e.dG_total_orth,
        "dGSpacingWT": e.best_dG_spacing_wt,
        "dGmRNAUnfoldingWT": e.best_dG_mrna_unfolding_wt,
        "dGTotalWT": e.dG_total_wt,
        # Legacy orthogonal aliases used by existing dashboard cards.
        "dGTotal": e.dG_total,
        "dGStandby": e.dG_standby,
        "dGSpacing": e.best_dG_spacing,
        "dGmRNAUnfolding": e.best_dG_mrna_unfolding,
        "orthTIR": e.orth_tir,
        "wtTIR": e.wt_tir,
        "tScore": e.t_score,
        "backend": e.backend,
    }


def build_dashboard_dataset(
    evals: List[CandidateEval],
    history: List[Dict[str, Any]],
    binding_sites: List[BindingSiteResult],
    inputs: Dict[str, Any],
    top_n: int = 20,
) -> Dict[str, Any]:
    top = evals[:top_n]
    best_ids = {e.candidate_id for e in top[:5]}
    site_payload = [asdict(s) for s in binding_sites if s.candidate_id in best_ids]

    return {
        "inputs": {
            "orthogonalAntiSD": inputs.get("orthogonalAntiSD", DEFAULT_ORTH_ASD),
            "wtAntiSD": inputs.get("wtAntiSD", DEFAULT_WT_ASD),
            "fivePrimeUTR": inputs.get("fivePrimeUTR", DEFAULT_FIVE_PRIME_UTR),
            "standby": inputs.get("standby", DEFAULT_STANDBY),
            "cdsStart": inputs.get("cdsStart", DEFAULT_CDS_START),
            "targetExpression": inputs.get("targetExpression", "High"),
            "wtPenaltyConstant": inputs.get("wtPenaltyConstant", DEFAULT_WT_PENALTY_CONSTANT),
        },
        "candidates": [candidate_to_dashboard_dict(e) for e in top],
        "fitnessData": [dict(h) for h in history],
        "scatterPoints": [
            {
                "wtLeakage": max(1e-8, e.wt_tir),
                "binding": max(1e-8, e.orth_tir),
                "access": e.rbs_access,
                "rbs": e.rbs,
                "spacer": e.spacer,
                "fitness": e.fitness,
            }
            for e in evals
        ],
        "bindingSites": site_payload,
    }


def write_outputs(
    out_dir: str,
    evals: List[CandidateEval],
    history: List[Dict[str, Any]],
    binding_sites: List[BindingSiteResult],
    dataset: Dict[str, Any],
) -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    (out / "rbs_dataset.json").write_text(json.dumps(dataset, indent=2))

    # Candidates CSV
    cand_rows = [asdict(e) for e in evals]
    if cand_rows:
        with (out / "candidates.csv").open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(cand_rows[0].keys()))
            writer.writeheader()
            writer.writerows(cand_rows)

    # Fitness CSV
    if history:
        with (out / "fitness.csv").open("w", newline="") as f:
            keys = list(history[0].keys())
            writer = csv.DictWriter(f, fieldnames=keys)
            writer.writeheader()
            writer.writerows(history)

    # Landscape CSV, React-style
    landscape_rows = [
        {
            "wtLeakage": max(1e-8, e.wt_tir),
            "binding": e.orth_tir,
            "access": e.rbs_access,
            "rbs": e.rbs,
            "spacer": e.spacer,
            "fitness": e.fitness,
        }
        for e in evals
    ]
    if landscape_rows:
        with (out / "landscape.csv").open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(landscape_rows[0].keys()))
            writer.writeheader()
            writer.writerows(landscape_rows)

    # Binding sites CSV
    site_rows = [asdict(s) for s in binding_sites]
    if site_rows:
        with (out / "binding_sites.csv").open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(site_rows[0].keys()))
            writer.writeheader()
            writer.writerows(site_rows)


# -------------------------
# CLI
# -------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="RiboGuard GA engine: initial candidates -> optimized RBS/spacer outputs")
    parser.add_argument("--seeds", type=str, default=None, help="Friend's initial candidates JSON/CSV")
    parser.add_argument("--out", type=str, default="riboguard_outputs", help="Output directory")
    parser.add_argument("--orth-asd", type=str, default=DEFAULT_ORTH_ASD)
    parser.add_argument("--wt-asd", type=str, default=DEFAULT_WT_ASD)
    parser.add_argument("--flank", type=str, default=DEFAULT_FLANK)
    parser.add_argument("--cds-start", type=str, default=DEFAULT_CDS_START)
    parser.add_argument("--generations", type=int, default=25)
    parser.add_argument("--population", type=int, default=80)
    parser.add_argument("--wt-penalty", type=float, default=DEFAULT_WT_PENALTY_CONSTANT)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--parallel", action="store_true", help="Evaluate candidates across multiple CPU cores")
    parser.add_argument("--workers", type=int, default=None, help="Worker processes when --parallel is set (default: CPU count - 1)")
    args = parser.parse_args()

    initial_candidates = load_initial_candidates(args.seeds)

    # If friend's candidate file has fewer candidates than population size,
    # the GA auto-fills using guided seeds from reverse complement.
    evals, history, binding_sites, _final_pop = run_ga(
        initial_candidates=initial_candidates,
        orth_anti_sd=args.orth_asd,
        wt_anti_sd=args.wt_asd,
        default_five_prime_utr=args.flank,
        default_standby=DEFAULT_STANDBY,
        default_cds_start=args.cds_start,
        generations=args.generations,
        population_size=args.population,
        wt_penalty_constant=args.wt_penalty,
        seed=args.seed,
        parallel=args.parallel,
        max_workers=args.workers,
    )

    dataset = build_dashboard_dataset(
        evals=evals,
        history=history,
        binding_sites=binding_sites,
        inputs={
            "orthogonalAntiSD": args.orth_asd,
            "wtAntiSD": args.wt_asd,
            "cdsStart": args.cds_start,
            "targetExpression": "High",
            "wtPenaltyConstant": args.wt_penalty,
        },
        top_n=20,
    )
    write_outputs(args.out, evals, history, binding_sites, dataset)

    print(f"Done. Outputs written to: {args.out}")
    if evals:
        best = evals[0]
        print("Best candidate:")
        print(f"  RBS: {best.rbs}")
        print(f"  Spacer: {best.spacer}")
        print(f"  Orth TIR: {best.orth_tir:.4g}")
        print(f"  WT TIR: {best.wt_tir:.4g}")
        print(f"  T-score: {best.t_score:.4g}")


if __name__ == "__main__":
    main()

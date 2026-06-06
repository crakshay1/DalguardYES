"""
spacing_analyser.py
───────────────────
Given an mRNA and an aSD sequence, computes and displays:

  - ASCII diagram of the SD:aSD duplex with base-pair connectors
  - All positional indices (last_mRNA_nt, farthest_3p_rRNA, distance_to_start)
  - Full step-by-step aligned spacing calculation
  - dG_spacing with equation expansion

Requires: ViennaRNA Python bindings
    pip install ViennaRNA

Usage:
    from spacing_analyser import analyse_spacing

    analyse_spacing(
        mRNA_input = "AAAGGAGAAAAATGAAACGT",
        aSD_input  = "ACCUCCUUA",           # E. coli default
    )
"""

import RNA
import math

# ── Spacing penalty parameters (Salis 2009) ───────────────────────────────────
OPTIMAL_SPACING = 5
PUSH = [12.2, 2.5, 2.0, 3.0]  # sigmoidal under-spacing: [A, B, C, D]
PULL = [0.048, 0.24, 0.0]  # quadratic over-spacing:  [a, b, c]

ECOLI_ASD = "ACCUCCUUA"


def to_rna(s: str) -> str:
    return s.upper().replace("T", "U")


def parse_pairs(structure: str) -> list[tuple[int, int]]:
    """Extract all (i, j) base pairs from dot-bracket. Both 0-indexed."""
    stack, pairs = [], []
    for i, c in enumerate(structure):
        if c == "(":
            stack.append(i)
        elif c == ")":
            j = stack.pop()
            pairs.append((j, i))
    return pairs


def calc_dG_spacing(spacing: float) -> float:
    if spacing >= 1e10:
        return 1e10
    ds = spacing - OPTIMAL_SPACING
    if spacing < OPTIMAL_SPACING:
        A, B, C, D = PUSH
        return A / (1.0 + math.exp(B * (ds + C))) ** D
    else:
        a, b, c = PULL
        return a * ds**2 + b * ds + c


def analyse_spacing(mRNA_input: str, aSD_input: str = ECOLI_ASD) -> dict:
    """
    Analyse SD:aSD spacing for a given mRNA and aSD sequence.

    Parameters
    ----------
    mRNA_input : mRNA sequence (DNA or RNA). Must contain AUG.
    aSD_input  : 3' tail of 16S rRNA, 5'->3' (DNA or RNA).
                 Defaults to E. coli sequence.

    Returns
    -------
    dict with all computed values, and prints a full annotated report.
    """

    mRNA = to_rna(mRNA_input)
    aSD = to_rna(aSD_input)
    r_len = len(aSD)

    # ── Find AUG ──────────────────────────────────────────────────────────────
    aug_pos = mRNA.find("AUG")
    if aug_pos == -1:
        raise ValueError("No AUG start codon found in mRNA sequence.")

    mRNA_upstream = mRNA[:aug_pos]
    m_len = len(mRNA_upstream)

    # ── Cofold mRNA upstream with aSD ─────────────────────────────────────────
    # ViennaRNA cofold: concatenate with '&', output is m_len + r_len chars (no separator)
    fc = RNA.fold_compound(mRNA_upstream + "&" + aSD, RNA.md())
    ss, dG_cofold = fc.mfe()

    mRNA_struct = ss[:m_len]
    rRNA_struct = ss[m_len:]

    # ── Extract intermolecular pairs ──────────────────────────────────────────
    # In the concatenated structure, rRNA positions have global index >= m_len
    all_pairs = parse_pairs(ss)

    mRNA_rRNA_pairs = []
    for x, y in all_pairs:
        if x < m_len <= y:
            mRNA_rRNA_pairs.append((x, y))
        elif y < m_len <= x:
            mRNA_rRNA_pairs.append((y, x))

    if not mRNA_rRNA_pairs:
        print("No SD:aSD intermolecular pairing detected.")
        return {}

    # ── Convert rRNA global indices to 5'->3' 1-indexed positions ─────────────
    # global_idx = m_len + (0-based offset into aSD)
    # 5'->3' position (1-indexed) = global_idx - m_len + 1
    def rRNA_5p_pos(global_idx: int) -> int:
        return global_idx - m_len + 1

    pair_map = {}  # mRNA 0-idx -> aSD 5'->3' 1-idx
    for mx, ry in mRNA_rRNA_pairs:
        pair_map[mx] = rRNA_5p_pos(ry)

    paired_rRNA_positions = list(pair_map.values())
    paired_mRNA_positions = list(pair_map.keys())

    # ── Salis indices ──────────────────────────────────────────────────────────
    # farthest_3p_rRNA: highest 5'->3' index among paired rRNA nts
    # (3' end of aSD = highest index = farthest from 5' anchor in 30S body)
    farthest_3p_rRNA = max(paired_rRNA_positions)

    # The mRNA nt paired with farthest_3p_rRNA
    last_mRNA_nt = paired_mRNA_positions[paired_rRNA_positions.index(farthest_3p_rRNA)]

    # distance_to_start: mRNA nt count from last_mRNA_nt to AUG (exclusive)
    distance_to_start = aug_pos - last_mRNA_nt

    # aligned_spacing: geometry-corrected distance
    aligned_spacing = distance_to_start - farthest_3p_rRNA

    dG_spacing = calc_dG_spacing(aligned_spacing)
    ds = aligned_spacing - OPTIMAL_SPACING

    # ── Print report ──────────────────────────────────────────────────────────
    W = 70
    div = "─" * W
    hdiv = "═" * W

    print()
    print(hdiv)
    print("  COFOLD ANALYSIS")
    print(hdiv)
    print()

    # mRNA line
    print("mRNA 5'- ", end="")
    for nt in mRNA_upstream:
        print(f"{nt} ", end="")
    print(f"· {mRNA[aug_pos:aug_pos+3]} -3'")

    # Connector line
    print("         ", end="")
    for i in range(m_len):
        print("| " if i in pair_map else "  ", end="")
    print()

    # aSD line — printed antiparallel (3'->5' left to right under mRNA)
    # For each mRNA position, show the paired aSD nucleotide if any
    print("         ", end="")
    for i in range(m_len):
        if i in pair_map:
            rpos = pair_map[i]  # 5'->3' 1-indexed
            print(f"{aSD[rpos-1]} ", end="")
        else:
            print(". ", end="")

    min_rpos = min(paired_rRNA_positions)
    max_rpos = max(paired_rRNA_positions)
    print(f"  aSD 3'←5' (paired: pos {max_rpos}→{min_rpos})")

    print()
    print(div)
    print("  INDICES")
    print(div)
    print()

    # mRNA position ruler
    print("  mRNA positions (0-indexed):")
    print("  pos: ", end="")
    for i in range(m_len):
        print(f"{i:<2}", end="")
    print(f"  AUG={aug_pos}")

    print("  seq: ", end="")
    for nt in mRNA_upstream:
        print(f"{nt:<2}", end="")
    print()
    print()

    print(
        f"  last_mRNA_nt     = {last_mRNA_nt}"
        f"  ← 0-idx, nt='{mRNA_upstream[last_mRNA_nt]}'"
        f"  (paired with farthest_3p_rRNA)"
    )
    print(f"  aug_pos          = {aug_pos}" f"  ← 0-idx, A of AUG")
    print()

    # aSD position ruler
    print("  aSD positions (5'->3', 1-indexed):")
    print("  pos: ", end="")
    for i in range(1, r_len + 1):
        print(f"{i:<2}", end="")
    print()
    print("  seq: ", end="")
    for nt in aSD:
        print(f"{nt:<2}", end="")
    print()
    print()

    print(
        f"  farthest_3p_rRNA = {farthest_3p_rRNA}"
        f"  ← 5'->3' pos, nt='{aSD[farthest_3p_rRNA-1]}'"
        f"  (most 3' paired rRNA nt)"
    )
    print(f"  all paired aSD positions: {sorted(paired_rRNA_positions)}")
    print()

    print(div)
    print("  SPACING CALCULATION")
    print(div)
    print()

    print(f"  distance_to_start  =  aug_pos - last_mRNA_nt")
    print(f"                     =  {aug_pos} - {last_mRNA_nt}")
    print(f"                     =  {distance_to_start}")
    print()
    print(f"  aligned_spacing    =  distance_to_start - farthest_3p_rRNA")
    print(f"                     =  {distance_to_start} - {farthest_3p_rRNA}")
    print(f"                     =  {aligned_spacing}")
    print()
    print(f"  ds                 =  aligned_spacing - optimal_spacing")
    print(f"                     =  {aligned_spacing} - {OPTIMAL_SPACING}")
    print(f"                     =  {ds}")
    print()

    if aligned_spacing < OPTIMAL_SPACING:
        A, B, C, D = PUSH
        inner = math.exp(B * (ds + C))
        bracket = 1.0 + inner
        penalty = bracket**D
        print(f"  under-spacing → sigmoidal wall")
        print(f"  dG_spacing  =  A / (1 + exp(B*(ds+C)))^D")
        print(f"              =  {A} / (1 + exp({B}*({ds}+{C})))^{D}")
        print(f"              =  {A} / (1 + {inner:.6f})^{D}")
        print(f"              =  {A} / ({bracket:.6f})^{D}")
        print(f"              =  {A} / {penalty:.6f}")
        print(f"              =  {dG_spacing:.4f} kcal/mol")
    else:
        a, b, c = PULL
        term1 = a * ds**2
        term2 = b * ds
        print(f"  over-spacing → quadratic")
        print(f"  dG_spacing  =  a*ds^2 + b*ds + c")
        print(f"              =  {a}*({ds})^2 + {b}*({ds}) + {c}")
        print(f"              =  {term1:.4f} + {term2:.4f} + {c}")
        print(f"              =  {dG_spacing:.4f} kcal/mol")

    print()
    print(f"  dG_cofold          =  {dG_cofold:.4f} kcal/mol")
    print(f"  dG_spacing         =  {dG_spacing:.4f} kcal/mol")
    print(hdiv)

    return {
        "mRNA": mRNA,
        "aSD": aSD,
        "aug_pos": aug_pos,
        "last_mRNA_nt": last_mRNA_nt,
        "farthest_3p_rRNA": farthest_3p_rRNA,
        "distance_to_start": distance_to_start,
        "aligned_spacing": aligned_spacing,
        "dG_cofold": dG_cofold,
        "dG_spacing": dG_spacing,
        "pair_map": pair_map,
    }


# ── Demo ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":

    print("\n>>> CASE 1: Strong SD, suboptimal spacing")
    analyse_spacing(
        mRNA_input="AAAGGAGAAAATGAAACGT",
        aSD_input=ECOLI_ASD,
    )

    print("\n>>> CASE 2: Strong SD, optimal spacing")
    analyse_spacing(
        mRNA_input="AAAGGAGAAAAAAAATGAAACGT",
        aSD_input=ECOLI_ASD,
    )

    print("\n>>> CASE 3: Custom o-aSD, full pairing (explodes)")
    analyse_spacing(
        mRNA_input="AAATAAGGAGGTGATAAATGAAACGT",
        aSD_input="ATCACCTCCTTA",
    )

    print("\n>>> CASE 4: Custom o-aSD, core only, correct spacing")
    analyse_spacing(
        mRNA_input="AAATAAGGAGAAAAAAAATGAAACGT",
        aSD_input="ATCACCTCCTTA",
    )

"""
standby_analyser.py
───────────────────
Computes and displays the standby site penalty (dG_standby) following
the Salis (2009) algorithm exactly.

The standby site is the STANDBY_LEN nucleotides immediately 5' of the
SD binding region. The penalty measures the energetic cost of keeping
those nucleotides single-stranded so the 30S can dock initially.

Algorithm:
  1. Cofold mRNA upstream + aSD  →  energy_before (unconstrained MFE)
  2. Identify most_5p_mRNA  =  left edge of SD:aSD binding on mRNA
  3. Define standby site    =  STANDBY_LEN nt immediately 5' of most_5p_mRNA
  4. Construct constrained structure:
       - pre-standby region folds freely
       - standby site forced single-stranded (absent from pair list)
       - SD+rRNA pairs preserved from cofold exactly
  5. Evaluate energy of constrained structure  →  energy_after
  6. dG_standby = energy_before - energy_after, capped at 0

Requires: ViennaRNA Python bindings
    pip install ViennaRNA

Usage:
    from standby_analyser import calc_standby, print_standby
    print_standby("CCGGCCGGAAGGAGAAAAAA", "ACCUCCUUA")
"""

import RNA

md = RNA.md()
STANDBY_LEN = 4
ECOLI_ASD = "ACCUCCUUA"


# ── Utilities ─────────────────────────────────────────────────────────────────


def to_rna(s: str) -> str:
    return s.upper().replace("T", "U")


def fold_alone(seq: str) -> tuple[float, str]:
    if not seq:
        return 0.0, ""
    fc = RNA.fold_compound(seq, md)
    ss, dG = fc.mfe()
    return dG, ss


def parse_pairs(ss: str) -> list[tuple[int, int]]:
    """Extract (i, j) base pairs from dot-bracket. Both 0-indexed."""
    stack, pairs = [], []
    for i, c in enumerate(ss):
        if c == "(":
            stack.append(i)
        elif c == ")":
            j = stack.pop()
            pairs.append((j, i))
    return pairs


def pairs_to_dotbracket(n: int, bp_x: list, bp_y: list) -> str:
    """Reconstruct dot-bracket string from explicit pair lists."""
    ss = ["."] * n
    for x, y in zip(bp_x, bp_y):
        ss[x] = "("
        ss[y] = ")"
    return "".join(ss)


def cofold(seq1: str, seq2: str) -> tuple[str, float]:
    """
    Cofold two RNA strands. Output structure is len(seq1)+len(seq2) chars
    with no separator — ViennaRNA strips the '&' from the output.
    """
    fc = RNA.fold_compound(seq1 + "&" + seq2, md)
    ss, dG = fc.mfe()
    return ss, dG


def eval_constrained_energy(mRNA: str, rRNA: str, bp_x: list, bp_y: list) -> float:
    """
    Evaluate the energy of a specific base pair assignment without
    re-folding. Reconstructs dot-bracket and uses ViennaRNA eval_structure.
    """
    n = len(mRNA) + len(rRNA)
    ss = pairs_to_dotbracket(n, bp_x, bp_y)
    fc = RNA.fold_compound(mRNA + "&" + rRNA, md)
    return fc.eval_structure(ss)


# ── Core calculation ──────────────────────────────────────────────────────────


def calc_standby(mRNA_upstream: str, aSD: str = ECOLI_ASD) -> tuple[float, dict]:
    """
    Compute dG_standby for a given mRNA upstream region and aSD sequence.

    Parameters
    ----------
    mRNA_upstream : sequence from 5' end to (not including) AUG
    aSD           : 3' tail of 16S rRNA, 5'->3'

    Returns
    -------
    (dG_standby, result_dict)
    dG_standby is 0 or negative (penalty is a cost, sign convention matches
    Salis: it is subtracted in the TIR equation via the standby term).
    """
    mRNA = to_rna(mRNA_upstream)
    rRNA = to_rna(aSD)
    m_len = len(mRNA)

    # ── Step 1: cofold ─────────────────────────────────────────────────────────
    ss_co, energy_before = cofold(mRNA, rRNA)
    all_pairs = parse_pairs(ss_co)
    bp_x = [p[0] for p in all_pairs]
    bp_y = [p[1] for p in all_pairs]

    # Identify intermolecular mRNA:rRNA pairs
    mRNA_rRNA_pairs = []
    for x, y in all_pairs:
        if x < m_len <= y:
            mRNA_rRNA_pairs.append((x, y))
        elif y < m_len <= x:
            mRNA_rRNA_pairs.append((y, x))

    if not mRNA_rRNA_pairs:
        return 0.0, {}

    # ── Step 2: most_5p_mRNA ───────────────────────────────────────────────────
    # Left edge of SD:aSD binding region on the mRNA
    most_5p_mRNA = min(p[0] for p in mRNA_rRNA_pairs)

    # ── Step 3: define regions ─────────────────────────────────────────────────
    standby_start = max(0, most_5p_mRNA - STANDBY_LEN - 1)
    standby_end = most_5p_mRNA
    mRNA_subsequence = mRNA[0:standby_start]  # pre-standby: folds freely
    standby_site = mRNA[standby_start:standby_end]  # forced single-stranded

    # ── Step 4: extract SD+rRNA pairs (preserved) ──────────────────────────────
    # All pairs at or downstream of most_5p_mRNA
    bp_x_3p = [x for x, y in zip(bp_x, bp_y) if x >= most_5p_mRNA]
    bp_y_3p = [y for x, y in zip(bp_x, bp_y) if x >= most_5p_mRNA]

    # ── Step 5: fold pre-standby alone ─────────────────────────────────────────
    dG_5p, ss_5p = fold_alone(mRNA_subsequence)
    pairs_5p = parse_pairs(ss_5p)
    bp_x_5p = [p[0] for p in pairs_5p]
    bp_y_5p = [p[1] for p in pairs_5p]

    # ── Step 6: assemble constrained pair list ──────────────────────────────────
    # pre-standby pairs (free)  +  SD/rRNA pairs (preserved)
    # standby site nucleotides are absent  =  forced single-stranded
    bp_x_after = bp_x_5p + bp_x_3p
    bp_y_after = bp_y_5p + bp_y_3p

    # ── Step 7: evaluate constrained structure energy ───────────────────────────
    energy_after = eval_constrained_energy(mRNA, rRNA, bp_x_after, bp_y_after)

    # ── Step 8: penalty ─────────────────────────────────────────────────────────
    # energy_before is the unconstrained MFE — always <= energy_after
    # so raw dG_standby = energy_before - energy_after <= 0
    # cap at 0: a free standby site costs nothing
    dG_standby_raw = energy_before - energy_after
    dG_standby = min(0.0, dG_standby_raw)

    constrained_struct = pairs_to_dotbracket(
        len(mRNA) + len(rRNA), bp_x_after, bp_y_after
    )

    return dG_standby, {
        "mRNA": mRNA,
        "rRNA": rRNA,
        "cofold_struct": ss_co,
        "energy_before": energy_before,
        "most_5p_mRNA": most_5p_mRNA,
        "standby_start": standby_start,
        "standby_end": standby_end,
        "mRNA_subsequence": mRNA_subsequence,
        "standby_site": standby_site,
        "ss_5p": ss_5p,
        "dG_5p": dG_5p,
        "constrained_struct": constrained_struct,
        "energy_after": energy_after,
        "dG_standby_raw": dG_standby_raw,
        "dG_standby": dG_standby,
    }


# ── Display ───────────────────────────────────────────────────────────────────


def print_standby(mRNA_upstream: str, aSD: str = ECOLI_ASD, label: str = "") -> None:
    """Compute and print annotated standby site analysis."""

    dG_standby, r = calc_standby(mRNA_upstream, aSD)
    if not r:
        print("No SD:aSD pairing detected — dG_standby = 0")
        return

    mRNA = r["mRNA"]
    rRNA = r["rRNA"]
    m_len = len(mRNA)
    ss = r["cofold_struct"]

    # Build mRNA->rRNA pair map for diagram
    all_pairs = parse_pairs(ss)
    mRNA_rRNA_map = {}
    for x, y in all_pairs:
        if x < m_len <= y:
            mRNA_rRNA_map[x] = y - m_len  # rRNA 0-idx
        elif y < m_len <= x:
            mRNA_rRNA_map[y] = x - m_len

    sb = r["standby_start"]
    sbe = r["standby_end"]
    m5p = r["most_5p_mRNA"]

    W = 68
    print()
    print("═" * W)
    if label:
        print(f"  {label}")
    print("═" * W)
    print()

    # mRNA line
    print("  mRNA 5'- ", end="")
    for nt in mRNA:
        print(f"{nt} ", end="")
    print("-3'")

    # Connector line
    print("           ", end="")
    for i in range(m_len):
        print("| " if i in mRNA_rRNA_map else "  ", end="")
    print()

    # aSD line (antiparallel — 3'->5' left to right)
    print("  aSD  3'- ", end="")
    for i in range(m_len):
        if i in mRNA_rRNA_map:
            rpos = mRNA_rRNA_map[i]
            print(f"{rRNA[rpos]} ", end="")
        else:
            print(". ", end="")
    print("-5'")

    print()

    # Region annotation line
    print("           ", end="")
    for i in range(m_len):
        if i < sb:
            print("· ", end="")
        elif sb <= i < sbe:
            print("S ", end="")
        elif i in mRNA_rRNA_map:
            print("D ", end="")
        else:
            print("  ", end="")
    print("  · = pre-standby  S = standby  D = SD")

    print()
    print("─" * W)
    print("  REGIONS")
    print("─" * W)
    pre = r["mRNA_subsequence"] or "(empty)"
    print(f"  pre-standby:      '{pre}'")
    print(f"                    pos 0–{sb-1 if sb > 0 else 'n/a'}  folds freely")
    print(f"  standby site:     '{r['standby_site']}'")
    print(f"                    pos {sb}–{sbe-1}  forced single-stranded")
    print(f"  SD left edge:     pos {m5p}  (most_5p_mRNA)")

    print()
    print("─" * W)
    print("  ENERGY CALCULATION")
    print("─" * W)
    print()
    print(f"  Step 1  cofold (unconstrained MFE)")
    print(f"          structure:    {r['cofold_struct']}")
    print(f"          energy_before = {r['energy_before']:.4f} kcal/mol")
    print()
    print(f"  Step 2  fold pre-standby alone")
    print(f"          sequence:  '{r['mRNA_subsequence'] or '(empty)'}'")
    print(f"          structure: '{r['ss_5p'] or '(empty)'}'")
    print(f"          dG_5p      = {r['dG_5p']:.4f} kcal/mol")
    print()
    print(f"  Step 3  assemble constrained structure")
    print(f"          pre-standby pairs  +  SD/rRNA pairs (from cofold)")
    print(f"          standby site nucleotides absent  →  forced unpaired")
    print(f"          constrained: {r['constrained_struct']}")
    print()
    print(f"  Step 4  evaluate constrained energy")
    print(f"          energy_after  = {r['energy_after']:.4f} kcal/mol")
    print()
    print(f"  Step 5  penalty")
    print(f"          dG_standby  =  energy_before - energy_after")
    print(
        f"                      =  {r['energy_before']:.4f} - ({r['energy_after']:.4f})"
    )
    print(f"                      =  {r['dG_standby_raw']:.4f}")
    print(f"                         capped at min(value, 0)")
    print(f"                      =  {dG_standby:.4f} kcal/mol")
    print("═" * W)


# ── Demo ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":

    print_standby(
        "AAAGGAGAAAAAA",
        ECOLI_ASD,
        label="Case 1: clean standby — poly-A context, no structure",
    )

    print_standby(
        "CCGGCCGGAAGGAGAAAAAA",
        ECOLI_ASD,
        label="Case 2: structured pre-standby — GC stem upstream of standby",
    )

    print_standby(
        "CCCCAAGGAGGGGGAAAAAA",
        ECOLI_ASD,
        label="Case 3: structured standby — GGGG standby pairs with upstream CCCC",
    )

    print_standby(
        "AAATAAGGAGGTGATAAAAAAA",
        "ATCACCTCCTTA",
        label="Case 4: custom o-aSD — full pairing, poly-A standby",
    )

    print_standby(
        "UUUACCUCCGACAGUAAAAAAGGAGAAAAAAAAUG",
        "ACCUCCUUA",
        label="Case 5: custom o-aSD",
    )

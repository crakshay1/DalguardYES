"""

Typical usage
-------------
    from dg_standby import dG_standby

    mRNA_upstream = "CCGGCCGGAAGGAGAAAAAA"   # everything upstream of AUG
    aSD           = "ACCUCCUUA"              # E. coli, or your o-aSD

    penalty = dG_standby(mRNA_upstream, aSD)
"""

import RNA

# ── Module-level constants ─────────────────────────────────────────────────────

STANDBY_LEN = 4  # nt of standby site immediately 5' of SD binding region
ECOLI_ASD = "ACCUCCUUA"

_md = RNA.md()  # ViennaRNA model details — instantiated once


# ── Private utilities ──────────────────────────────────────────────────────────


def _to_rna(s: str) -> str:
    return s.upper().replace("T", "U")


def _fold_alone(seq: str) -> tuple[float, str]:
    """MFE fold a single RNA strand. Returns (dG, dot-bracket)."""
    if not seq:
        return 0.0, ""
    fc = RNA.fold_compound(seq, _md)
    ss, dG = fc.mfe()
    return dG, ss


def _parse_pairs(ss: str) -> list[tuple[int, int]]:
    """
    Extract all (i, j) base pairs from a dot-bracket string.
    Both indices are 0-based.
    """
    stack, pairs = [], []
    for i, c in enumerate(ss):
        if c == "(":
            stack.append(i)
        elif c == ")":
            j = stack.pop()
            pairs.append((j, i))
    return pairs


def _pairs_to_dotbracket(n: int, bp_x: list[int], bp_y: list[int]) -> str:
    """Reconstruct dot-bracket from explicit pair lists."""
    ss = ["."] * n
    for x, y in zip(bp_x, bp_y):
        ss[x] = "("
        ss[y] = ")"
    return "".join(ss)


def _cofold(seq1: str, seq2: str) -> tuple[str, float]:
    """
    Cofold two RNA strands via ViennaRNA.
    Output structure is len(seq1)+len(seq2) characters — no separator.
    ViennaRNA strips the '&' from the structure output.
    """
    fc = RNA.fold_compound(seq1 + "&" + seq2, _md)
    ss, dG = fc.mfe()
    return ss, dG


def _eval_structure(mRNA: str, rRNA: str, bp_x: list[int], bp_y: list[int]) -> float:
    """
    Evaluate the free energy of a specific base pair assignment.
    Uses ViennaRNA eval_structure — a single traversal, no DP.
    The dot-bracket must not contain '&'; pass the full concatenated sequence.
    """
    n = len(mRNA) + len(rRNA)
    ss = _pairs_to_dotbracket(n, bp_x, bp_y)
    fc = RNA.fold_compound(mRNA + "&" + rRNA, _md)
    return fc.eval_structure(ss)


# ── Public function ────────────────────────────────────────────────────────────


def dG_standby(
    mRNA_upstream: str,
    aSD: str,
    *,
    cofold_struct: str | None = None,
    cofold_energy: float | None = None,
) -> float:
    """
    mRNA_upstream : str
        mRNA sequence from the 5' end up to (not including) the A of AUG.
        AKA +1 -> (A)UG

    aSD : str
        3' tail of 16S rRNA (anti-SD), 5'->3', DNA or RNA.

    cofold_struct : str, optional
        Pre-computed cofold dot-bracket structure (len(mRNA)+len(aSD) chars,
        no '&' separator). If provided together with cofold_energy, skips the
        cofold step — useful when the TIR calculator has already run cofold
        and you want to avoid a redundant fold call.

    cofold_energy : float, optional
        Pre-computed cofold MFE energy (kcal/mol). Must be provided together
        with cofold_struct to take effect.

    Returns
    -------
        dG_standby in kcal/mol. Always <= 0.
        Returns 0.0 if no SD:aSD pairing is detected.

    """

    mRNA = _to_rna(mRNA_upstream)
    rRNA = _to_rna(aSD)
    m_len = len(mRNA)

    # Step 1: cofold (or use pre-computed dG_cofold result)
    if cofold_struct is not None and cofold_energy is not None:
        ss_co = cofold_struct
        energy_before = cofold_energy
    else:
        ss_co, energy_before = _cofold(mRNA, rRNA)

    all_pairs = _parse_pairs(ss_co)
    bp_x = [p[0] for p in all_pairs]
    bp_y = [p[1] for p in all_pairs]

    # Step 2: find most_5p_mRNA
    # Left edge of SD:aSD binding region on the mRNA.
    mRNA_rRNA_pairs = []
    for x, y in all_pairs:
        if x < m_len <= y:
            mRNA_rRNA_pairs.append((x, y))
        elif y < m_len <= x:
            mRNA_rRNA_pairs.append((y, x))  # normalise: mRNA index first

    if not mRNA_rRNA_pairs:
        return 0.0  # no SD:aSD pairing — standby is irrelevant

    most_5p_mRNA = min(p[0] for p in mRNA_rRNA_pairs)

    # Step 3: define regions
    #
    #   5'─[pre_standby]─[standby STANDBY_LEN nt]─[SD...]─...─3'
    #       0            sb_start              sb_end=most_5p_mRNA
    #
    sb_start = max(0, most_5p_mRNA - STANDBY_LEN - 1)
    pre_standby_seq = mRNA[0:sb_start]

    # Step 4: SD+rRNA pairs — preserved from cofold
    # All pairs whose mRNA index is at or downstream of most_5p_mRNA.
    bp_x_3p = [x for x, y in zip(bp_x, bp_y) if x >= most_5p_mRNA]
    bp_y_3p = [y for x, y in zip(bp_x, bp_y) if x >= most_5p_mRNA]

    # Step 5: fold pre-standby alone
    _, ss_5p = _fold_alone(pre_standby_seq)
    pairs_5p = _parse_pairs(ss_5p)
    bp_x_5p = [p[0] for p in pairs_5p]
    bp_y_5p = [p[1] for p in pairs_5p]

    # Step 6: assemble constrained pair list
    # pre-standby pairs (free to fold)  +  SD/rRNA pairs (fixed from cofold)
    # Standby site nucleotides are absent from both lists
    bp_x_after = bp_x_5p + bp_x_3p
    bp_y_after = bp_y_5p + bp_y_3p

    # Step 7: evaluate constrained structure
    energy_after = _eval_structure(mRNA, rRNA, bp_x_after, bp_y_after)

    # Step 8: calculate penalty
    # a free standby site contributes nothing (not a bonus).
    return min(0.0, energy_before - energy_after)


# ── Smoke test ────────────────────────────────────────────────────────────────

if __name__ == "__main__":

    cases = [
        ("clean standby (poly-A)", "AAAGGAGAAAAAA", ECOLI_ASD),
        ("structured pre-standby (GC stem)", "CCGGCCGGAAGGAGAAAAAA", ECOLI_ASD),
        ("custom o-aSD", "AAATAAGGAGGTGATAAAAAAA", "ATCACCTCCTTA"),
    ]

    print(f"\n{'Case':<45} {'dG_standby':>12}")
    print("─" * 59)
    for label, mRNA_up, asd in cases:
        penalty = dG_standby(mRNA_up, asd)
        print(f"{label:<45} {penalty:>12.4f} kcal/mol")
    print()

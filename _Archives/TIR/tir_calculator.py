"""
tir_calculator.py
─────────────────
Translation Initiation Rate (TIR) calculator based on the Salis (2009)
thermodynamic model. Supports custom anti-SD sequences for orthogonal
ribosome design.

Requires: ViennaRNA Python bindings
    pip install ViennaRNA

Usage:
    from tir_calculator import calc_tir, print_result

    result = calc_tir("AAGTTAAGAGGCAAGCTTATGAAAGGTTTTGTT")
    print_result(result)

    # Orthogonal ribosome — supply custom aSD
    result = calc_tir(mRNA, rrna="ATCACCTCCTTA")

Reference:
    Salis HM, Mirsky EA, Voigt CA (2009) Automated design of synthetic
    ribosome binding sites to control protein expression. Nature Biotechnology.
"""

import math
import RNA  # ViennaRNA

# ── Model parameters (Salis RBS Calculator v1) ────────────────────────────────

OPTIMAL_SPACING = 5          # nt, optimal aligned SD-to-AUG distance

# Under-spacing penalty: sigmoidal wall
# dG = A / (1 + exp(B*(ds + C)))^D   where ds = spacing - optimal
PUSH = [12.2, 2.5, 2.0, 3.0]        # [A, B, C, D]

# Over-spacing penalty: quadratic
# dG = a*ds^2 + b*ds + c
PULL = [0.048, 0.24, 0.0]            # [a, b, c]

CUTOFF = 35                  # nt window around AUG used for folding
STANDBY_LEN = 4              # nt upstream of SD that must be single-stranded

# Start codon hybridization energies (kcal/mol) — empirical, Salis 2009
START_ENERGIES = {
    'AUG': -1.194,
    'GUG': -0.0748,
    'UUG': -0.0435,
    'CUG': -0.03406,
}

# E. coli 16S rRNA 3' tail (anti-SD), RNA, 5'->3'
ECOLI_RRNA = 'ACCUCCUUA'

RT = 1.987e-3 * 310.15      # kcal/mol at 37°C
INFINITY = 1e12

_md = RNA.md()               # ViennaRNA model details (default parameters)


# ── Utility ───────────────────────────────────────────────────────────────────

def to_rna(seq: str) -> str:
    return seq.upper().replace('T', 'U')


def fold_alone(seq: str) -> tuple[float, str]:
    """MFE fold a single RNA strand. Returns (dG, dot-bracket structure)."""
    if not seq:
        return 0.0, ''
    fc = RNA.fold_compound(seq, _md)
    ss, dG = fc.mfe()
    return dG, ss


def cofold(seq1: str, seq2: str) -> tuple[str, float]:
    """
    Cofold two RNA strands using ViennaRNA.
    Returns (structure, dG).
    Structure is len(seq1)+len(seq2) characters — no separator in output.
    """
    fc = RNA.fold_compound(seq1 + '&' + seq2, _md)
    ss, dG = fc.mfe()
    return ss, dG


def parse_pairs(structure: str, offset: int = 0) -> list[tuple[int, int]]:
    """
    Extract all (i, j) base pairs from a dot-bracket string.
    Positions are 0-indexed, with optional offset applied to both indices.
    """
    stack, pairs = [], []
    for i, c in enumerate(structure):
        if c == '(':
            stack.append(i + offset)
        elif c == ')':
            j = stack.pop()
            pairs.append((j, i + offset))
    return pairs


# ── Spacing ───────────────────────────────────────────────────────────────────

def calc_dG_spacing(spacing: float) -> float:
    """
    Compute the spacing penalty given aligned SD-to-AUG spacing in nt.

    Under optimal (< OPTIMAL_SPACING): sigmoidal wall — steep penalty
    for compressing mRNA against ribosome geometry.

    Over optimal (>= OPTIMAL_SPACING): quadratic — gentler penalty
    for stretching the spacer.
    """
    if spacing >= INFINITY:
        return INFINITY

    ds = spacing - OPTIMAL_SPACING

    if spacing < OPTIMAL_SPACING:
        A, B, C, D = PUSH
        return A / (1.0 + math.exp(B * (ds + C))) ** D
    else:
        a, b, c = PULL
        return a * ds**2 + b * ds + c


def find_aligned_spacing(
    mRNA_struct: str,
    rRNA_struct: str,
    start_pos_in_upstream: int,
) -> tuple[float, int, int]:
    """
    Compute aligned spacing from cofold dot-bracket structures.

    The aligned spacing is the distance from the farthest 3' mRNA nucleotide
    involved in SD:aSD base pairing to the AUG, minus the position of the
    farthest 3' rRNA nucleotide involved in that pairing.

    This accounts for where in the aSD the pairing terminates — not simply
    the sequence distance from SD motif start to AUG.

    Returns:
        (aligned_spacing, last_mRNA_nt_0idx, last_rRNA_nt_0idx)
        aligned_spacing = INFINITY if no SD:aSD pairing detected.
    """
    m_len = len(mRNA_struct)

    # All base pairs across the concatenated cofold structure
    # rRNA positions are offset by m_len in the combined structure
    all_pairs = parse_pairs(mRNA_struct + rRNA_struct)

    # Isolate mRNA:rRNA intermolecular pairs
    # One partner index < m_len (mRNA), other >= m_len (rRNA)
    mRNA_rRNA_pairs = []
    for x, y in all_pairs:
        if x < m_len <= y:
            mRNA_rRNA_pairs.append((x, y))
        elif y < m_len <= x:
            mRNA_rRNA_pairs.append((y, x))

    if not mRNA_rRNA_pairs:
        return INFINITY, -1, -1

    # Farthest 3' mRNA nucleotide in the SD:aSD duplex
    last_mRNA_nt = max(p[0] for p in mRNA_rRNA_pairs)

    # Farthest 3' rRNA nucleotide (0-indexed within rRNA)
    last_rRNA_nt = max(p[1] - m_len for p in mRNA_rRNA_pairs)

    distance_to_start = start_pos_in_upstream - last_mRNA_nt
    farthest_3p_rRNA  = last_rRNA_nt + 1   # convert to 1-indexed

    aligned_spacing = distance_to_start - farthest_3p_rRNA
    return aligned_spacing, last_mRNA_nt, last_rRNA_nt


# ── Main calculation ──────────────────────────────────────────────────────────

def calc_tir(mRNA_input: str, rrna: str = ECOLI_RRNA) -> dict:
    """
    Compute all thermodynamic components of the Salis TIR model for a
    given mRNA sequence and ribosomal anti-SD sequence.

    Args:
        mRNA_input: mRNA sequence (DNA or RNA, any case). Must contain AUG.
        rrna:       3' tail of 16S rRNA (anti-SD), 5'->3', RNA or DNA.
                    Defaults to E. coli sequence. Supply your o-aSD here
                    for orthogonal ribosome calculations.

    Returns:
        dict with keys:
            mRNA, start_pos, start_codon, aligned_spacing,
            dG_mRNA, dG_SD_aSD, dG_spacing, dG_standby, dG_start,
            dG_total, rel_rate, ss_mRNA, cofold_struct

    dG_total = dG_SD_aSD + dG_start + dG_spacing + dG_standby - dG_mRNA
    rel_rate = exp(-dG_total / RT)   [arbitrary units]
    """
    mRNA = to_rna(mRNA_input)
    rRNA = to_rna(rrna)

    # Locate first AUG
    start_pos = mRNA.find('AUG')
    if start_pos == -1:
        raise ValueError('No AUG start codon found in mRNA sequence.')
    start_codon = mRNA[start_pos:start_pos + 3]

    # Folding window
    win_start     = max(0, start_pos - CUTOFF)
    win_end       = min(len(mRNA), start_pos + CUTOFF)
    mRNA_window   = mRNA[win_start:win_end]
    mRNA_upstream = mRNA[win_start:start_pos]

    # ── 1. dG_mRNA ─────────────────────────────────────────────────────────────
    # MFE of mRNA window folding on its own — represents stored structural
    # energy that the ribosome must pay to access the SD and AUG.
    dG_mRNA, ss_mRNA = fold_alone(mRNA_window)

    # ── 2. dG_SD:aSD ───────────────────────────────────────────────────────────
    # Net hybridization energy between mRNA upstream region and rRNA.
    # = dG(cofold) - dG(mRNA alone) - dG(rRNA alone)
    dG_mRNA_up_alone, _ = fold_alone(mRNA_upstream)
    dG_rRNA_alone,    _ = fold_alone(rRNA)
    cofold_struct, dG_cofold = cofold(mRNA_upstream, rRNA)
    dG_SD_aSD = dG_cofold - dG_mRNA_up_alone - dG_rRNA_alone

    # ── 3. dG_spacing ──────────────────────────────────────────────────────────
    # Penalty for SD-to-AUG distance deviating from ribosome geometry optimum.
    m_len = len(mRNA_upstream)
    mRNA_up_struct = cofold_struct[:m_len]
    rRNA_struct    = cofold_struct[m_len:]

    start_pos_in_upstream = start_pos - win_start
    aligned_spacing, last_mRNA_nt, last_rRNA_nt = find_aligned_spacing(
        mRNA_up_struct, rRNA_struct, start_pos_in_upstream
    )
    dG_spacing = calc_dG_spacing(aligned_spacing)

    # ── 4. dG_standby ──────────────────────────────────────────────────────────
    # Cost of keeping the standby site (STANDBY_LEN nt upstream of SD) free
    # for initial 30S docking. Defined as <= 0 (never destabilizing).
    if last_mRNA_nt > STANDBY_LEN:
        pre_standby    = mRNA_upstream[:last_mRNA_nt - STANDBY_LEN]
        dG_pre, _      = fold_alone(pre_standby)
        dG_standby     = dG_cofold - (dG_pre + dG_rRNA_alone)
        if dG_standby > 0:
            dG_standby = 0.0
    else:
        dG_standby = 0.0

    # ── 5. dG_start ────────────────────────────────────────────────────────────
    # Empirical energy of initiator tRNA accommodation at the start codon.
    # Fixed lookup — not sequence-computed.
    dG_start = START_ENERGIES.get(start_codon, 0.0)

    # ── 6. Assemble dG_total and translation rate ───────────────────────────────
    dG_total = dG_SD_aSD + dG_start + dG_spacing + dG_standby - dG_mRNA
    rel_rate = math.exp(-dG_total / RT)

    return {
        'mRNA':            mRNA,
        'start_pos':       start_pos,
        'start_codon':     start_codon,
        'aligned_spacing': aligned_spacing,
        'dG_mRNA':         dG_mRNA,
        'dG_SD_aSD':       dG_SD_aSD,
        'dG_spacing':      dG_spacing,
        'dG_standby':      dG_standby,
        'dG_start':        dG_start,
        'dG_total':        dG_total,
        'rel_rate':        rel_rate,
        'ss_mRNA':         ss_mRNA,
        'cofold_struct':   cofold_struct,
    }


def print_result(r: dict, label: str = '') -> None:
    """Pretty-print all dG components and the relative translation rate."""
    if label:
        print(f'\n{"═"*44}')
        print(f'  {label}')
        print(f'{"═"*44}')
    print(f"mRNA:            {r['mRNA']}")
    print(f"AUG position:    {r['start_pos']}")
    print(f"Start codon:     {r['start_codon']}")
    print(f"Aligned spacing: {r['aligned_spacing']} nt")
    print(f"mRNA structure:  {r['ss_mRNA']}")
    print(f"\n{'Term':<22} {'kcal/mol':>10}")
    print('─' * 34)
    print(f"{'dG_SD:aSD':<22} {r['dG_SD_aSD']:>10.4f}")
    print(f"{'dG_start':<22} {r['dG_start']:>10.4f}")
    print(f"{'dG_spacing':<22} {r['dG_spacing']:>10.4f}")
    print(f"{'dG_standby':<22} {r['dG_standby']:>10.4f}")
    print(f"{'-dG_mRNA':<22} {-r['dG_mRNA']:>10.4f}")
    print('─' * 34)
    print(f"{'dG_total':<22} {r['dG_total']:>10.4f}")
    print(f"\nRelative rate:   {r['rel_rate']:.2f}  (arbitrary units)")


# ── Demo ──────────────────────────────────────────────────────────────────────

if __name__ == '__main__':

    cases = [
        ('Strong RBS (consensus SD)',
         'AAGTTAAGAGGCAAGCTTATGAAAGGTTTTGTTTTAGTTTT',
         ECOLI_RRNA),

        ('Medium RBS (partial SD: GAGG)',
         'AAGTTAAAGAGACAAGCTTATGAAAGGTTTTGTTTTAGTTTT',
         ECOLI_RRNA),

        ('Weak RBS (no SD)',
         'AAGTTAAAGTTTCAAGCTTATGAAAGGTTTTGTTTTAGTTTT',
         ECOLI_RRNA),

        ('Orthogonal RBS (custom o-aSD)',
         'AAGTTAATAAGGAGGTGATCAAGCTTATGAAAGGTTTTGTT',
         'ATCACCTCCTTA'),   # your o-aSD
    ]

    for label, mRNA, rrna in cases:
        r = calc_tir(mRNA, rrna=rrna)
        print_result(r, label=label)

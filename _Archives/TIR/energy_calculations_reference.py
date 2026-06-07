import math

# Salis RBS Calculator Implementation

# ── Constants ──────────────────────────────────────────────────────────────

infinity = 1e12

# E. coli 16S rRNA 3' tail (the anti-SD sequence)
rRNA = "acctcctta"

# Start codon hybridization energies (kcal/mol)
start_codon_energies = {
    "ATG": -1.194,
    "AUG": -1.194,
    "GTG": -0.0748,
    "GUG": -0.0748,
    "TTG": -0.0435,
    "UUG": -0.0435,
    "CTG": -0.03406,
    "CUG": -0.03406,
}

# Spacing penalty parameters
optimal_spacing = 5  # nt, optimum SD-to-AUG distance

# Under-spacing: sigmoidal wall
# dG = A / (1 + exp(B * (ds + C))) ** D
dG_spacing_constant_push = [12.2, 2.5, 2.0, 3.0]  # [A, B, C, D]

# Over-spacing: quadratic
# dG = a*ds^2 + b*ds + c
dG_spacing_constant_pull = [0.048, 0.24, 0.0]  # [a, b, c]

# Folding window: nt upstream and downstream of start codon to consider
cutoff = 35

# Standby site: nt upstream of SD that must be single-stranded
standby_site_length = 4

# ── Spacing penalty ────────────────────────────────────────────────────────


def calc_dG_spacing(self, aligned_spacing):

    ds = aligned_spacing - self.optimal_spacing

    if aligned_spacing < self.optimal_spacing:
        # Sigmoidal — steep wall for under-spacing
        A, B, C, D = self.dG_spacing_constant_push
        return A / (1.0 + math.exp(B * (ds + C))) ** D

    else:
        # Quadratic — gentler rise for over-spacing
        a, b, c = self.dG_spacing_constant_pull
        return a * ds**2 + b * ds + c


# ── Aligned spacing ────────────────────────────────────────────────────────


def calc_aligned_spacing(self, mRNA, start_pos, bp_x, bp_y):
    """
    Finds the farthest 3' rRNA nucleotide that base-pairs with
    the mRNA upstream of the start codon, then computes the
    distance from that mRNA partner to the AUG.

    bp_x, bp_y: paired lists of mRNA and rRNA positions
                as returned by the folding engine.
    """
    seq_len = len(mRNA) + self.rRNA_len
    Ok = False

    for rRNA_nt in range(seq_len, seq_len - self.rRNA_len, -1):
        if rRNA_nt in bp_y:
            rRNA_pos = bp_y.index(rRNA_nt)
            if bp_x[rRNA_pos] < start_pos:
                Ok = True
                farthest_3_prime_rRNA = rRNA_nt - len(mRNA)
                mRNA_nt = bp_x[rRNA_pos]
                distance_to_start = start_pos - mRNA_nt + 1
                break
            else:
                break

    if Ok:
        return distance_to_start - farthest_3_prime_rRNA
    else:
        return self.infinity

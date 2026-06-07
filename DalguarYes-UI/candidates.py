from __future__ import annotations
"""
    Generates RBS candidates from a mRNA (json format).
    @crakshay1
"""


"""
Strategy
─────────────────────────────────────────────────────────────────────────────
Stop-codon filter 
─────────────────────────────────────────────────────────────────────────────
Only FORWARD reading frames are checked, mRNA is single-stranded and
only translated in the 5′→3′ direction. Three forward frames are checked
in both the mutated spacer and in the junction window
(RBS tail + spacer + CDS head). The CDS itself is excluded from checking
because it naturally contains stop codons in off-frames.

Mutation strategy
─────────────────────────────────────────────────────────────────────────────
Three mutable regions per candidate:

- five_prime_flank  (partial)
    Positions 0 .. standby_start-1  → FROZEN
    Positions standby_start .. end  → independently randomised per position
    (each position: 50% chance of synonymous change, 50% kept)

- Core motif inside the full RBS is ALWAYS intact, but the other parts of the RBS are mutated.

- Spacer is freshly generated per candidate: random RNA of length drawn uniformly
in [spacer_len_min, spacer_len_max]. Filtered so no triplet in any of
the 3 forward reading frames is a stop codon.

Parameters
─────────────────────────────────────────────────────────────────────────────
  --input  / -i    Input JSON
  --output / -o    Output JSON path
  --n              Number of candidates to generate per seed  (default: 100)
  --standby-start  0-based index into five_prime_flank; positions BEFORE this
                   are frozen  (default: 0 → entire flank is mutable)
  --spacer-len     MIN MAX  tuple for random spacer length  (default: 4 - 7)
  --max-tries      Internal retry cap per candidate (default: 10000)
  --seed           RNG seed for reproducibility

Input JSON format
─────────────────────────────────────────────────────────────────────────────
{
  "name": "seed_002",
  "five_prime_flank": "UUUAAA",
  "rbs": "AAGGUACAAGUCU",
  "core": "UACAAG",
  "cds_start": "AUGGCUACUAAAGAAAACGCU",
  "mutable_regions": ["five_prime_flank", "rbs", "spacer"]
}

Output JSON format
─────────────────────────────────────────────────────────────────────────────
[
  {
    "name": "seed_002_c0001",
    "five_prime_flank": "UUUCGA",
    "rbs": "AAGGUACAAGUCU",
    "spacer": "AAUAAA",
    "cds_start": "AUGGCUACUAAAGAAAACGCU",
    "mutable_regions": ["five_prime_flank", "rbs", "spacer"]
  },
  ...
]
"""
import argparse
import json
import random
import sys
from pathlib import Path
from typing import Any



# Constants
# ─────────────────────────────────────────────────────────────────────────────

NUCLEOTIDES: tuple[str, ...] = ("A", "U", "G", "C")
STOP_CODONS: frozenset[str]  = frozenset({"UAA", "UAG", "UGA"})

# Sequence utilities
# ─────────────────────────────────────────────────────────────────────────────

def normalise(seq: str) -> str:
    """Uppercase + T→U."""
    return seq.upper().replace("T", "U")


def has_stop_fwd(seq: str) -> bool:
    """
        True if any triplet in any of the 3 FORWARD reading frames is a stop codon.
        Reverse complement is intentionally NOT checked, mRNA is single-stranded
        and only translated 5'→3'.
    """
    for frame in range(3):
        for i in range(frame, len(seq) - 2, 3):
            if seq[i:i + 3] in STOP_CODONS:
                return True
    return False


def junction_has_stop(flank: str, rbs: str, spacer: str, cds_start: str, window: int = 15) -> bool:
    """
        Check the junction zone (end of RBS + spacer) for stop codons in all 3
        forward frames. CDS is excluded, it naturally contains off-frame stops.
    """
    region = rbs[-window:] + spacer  
    return has_stop_fwd(region)


# Mutators
# ─────────────────────────────────────────────────────────────────────────────

def _randomise_positions(seq: str) -> str:
    """
        Per-position independent randomisation:
        each nucleotide has a 50% chance of being replaced by a different nt.
    """
    result = []
    for nt in seq:
        if random.random() < 0.5:
            result.append(random.choice([n for n in NUCLEOTIDES if n != nt]))
        else:
            result.append(nt)
    return "".join(result)


def mutate_flank(flank: str, standby_start: int) -> str:
    """
        Positions 0 .. standby_start-1  => FROZEN.
        Positions standby_start .. end  => randomised.
    """
    return flank[:standby_start] + _randomise_positions(flank[standby_start:])


def mutate_rbs_non_core(full_rbs: str, core: str) -> str | None:
    """
        Locate the core inside the full RBS and mutate everything around it.
        Returns the mutated full RBS (core preserved), or None if core not found.
    """
    core_idx = full_rbs.find(core)
    if core_idx == -1:
        return None
    prefix = _randomise_positions(full_rbs[:core_idx])
    suffix = _randomise_positions(full_rbs[core_idx + len(core):])
    return prefix + core + suffix


def random_spacer(length: int, max_tries: int = 1000) -> str | None:
    """
        Generate a random RNA of `length` with no stop codon in any of the
        3 forward reading frames. Returns None if max_tries is exhausted.
    """
    for _ in range(max_tries):
        seq = "".join(random.choice(NUCLEOTIDES) for _ in range(length))
        if not has_stop_fwd(seq):
            return seq
    return None



# Per-seed generation
# ─────────────────────────────────────────────────────────────────────────────

def generate_candidates(seed: dict[str, Any], n: int, standby_start: int, spacer_len_min: int, spacer_len_max: int, max_tries: int) -> list[dict[str, Any]]:
    """
        Generates a list of candidates with the methods created above.
    """
    name = seed["name"]
    flank = normalise(seed["five_prime_flank"])
    rbs = normalise(seed["rbs"])
    core = normalise(seed["core"])
    cds_start = normalise(seed["cds_start"])
    mutable = seed.get("mutable_regions", [])

    # Validation 
    if core not in rbs:
        raise ValueError(
            f"[{name}] Core '{core}' not found inside RBS '{rbs}'."
        )
    # Catch an unmutatable core
    for frame in range(3):
        for i in range(frame, len(core) - 2, 3):
            if core[i:i+3] in STOP_CODONS:
                raise ValueError(f"[{name}] Core '{core}' contains stop codon '{core[i:i+3]}': all candidates would be rejected.")
    if standby_start > len(flank):
        raise ValueError(f"[{name}] standby_start={standby_start} exceeds five_prime_flank length ({len(flank)}).")

    candidates: list[dict[str, Any]] = []
    attempts   = 0

    while len(candidates) < n and attempts < max_tries:
        attempts += 1
        # five_prime_flank
        new_flank = mutate_flank(flank, standby_start) if "five_prime_flank" in mutable else flank

        # RBS non-core region
        if "rbs" in mutable:
            new_rbs = mutate_rbs_non_core(rbs, core)
            if new_rbs is None:
                continue
        else:
            new_rbs = rbs

        # Spacer
        if "spacer" in mutable:
            slen = random.randint(spacer_len_min, spacer_len_max)
            new_spacer = random_spacer(slen)
            if new_spacer is None:
                continue
        else:
            new_spacer = normalise(seed.get("spacer", ""))

        # Junction stop-codon check
        if junction_has_stop(new_flank, new_rbs, new_spacer, cds_start):
            continue

        # Build record
        candidates.append({
            "name": f"{name}_c{len(candidates) + 1:04d}",
            "five_prime_flank": new_flank,
            "rbs": new_rbs,
            "spacer": new_spacer,
            "cds_start": cds_start,
            "mutable_regions":  ["five_prime_flank", "rbs", "spacer"]
        })

    return candidates



# Main
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="RBS candidate generator, DalguardYES", 
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--input",  "-i", required=True,
                        help="Input JSON (single seed or list).")
    parser.add_argument("--output", "-o", required=True,
                        help="Output JSON path.")
    parser.add_argument("--n", type=int, default=100,
                        help="Candidates per seed (default: 100).")
    parser.add_argument("--standby-start", type=int, default=0, metavar="N",
                        help="0-based index: positions before N in five_prime_flank "
                             "are frozen (default: 0).")
    parser.add_argument("--spacer-len", type=int, nargs=2, default=[4, 7],
                        metavar=("MIN", "MAX"),
                        help="Spacer length range (default: 4 7).")
    parser.add_argument("--max-tries", type=int, default=10000,
                        help="Max attempts per seed before giving up (default: 10000).")
    parser.add_argument("--seed", type=int, default=None,
                        help="RNG seed for reproducibility.")

    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    sp_min, sp_max = args.spacer_len
    if sp_min > sp_max:
        parser.error(f"--spacer-len MIN ({sp_min}) must be ≤ MAX ({sp_max}).")

    with open(args.input) as fh:
        raw = json.load(fh)
    seeds: list[dict] = raw if isinstance(raw, list) else [raw]

    print(f"[INFO] {len(seeds)} seed(s) loaded.")

    all_candidates: list[dict] = []
    for seed in seeds:
        print(f"\n[INFO] Seed '{seed['name']}' → generating {args.n} candidates …")
        cands = generate_candidates(
            seed          = seed,
            n             = args.n,
            standby_start = args.standby_start,
            spacer_len_min= sp_min,
            spacer_len_max= sp_max,
            max_tries     = args.max_tries,
        )
        all_candidates.extend(cands)
        print(f"       └─ {len(cands)} candidate(s) OK.")

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as fh:
        json.dump(all_candidates, fh, indent=2, ensure_ascii=False)

    print(f"\n Done. {len(all_candidates)} total candidate(s) → '{out}'")


if __name__ == "__main__":
    main()
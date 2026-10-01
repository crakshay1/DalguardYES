from __future__ import annotations
"""
Generate RBS candidates using one explicit candidate schema across the pipeline.

Canonical candidate sequence order
─────────────────────────────────────────────────────────────────────────────
    five_prime_utr + standby
    + rbs_left + rbs_core + rbs_right
    + spacer + cds_start

Mutation policy
─────────────────────────────────────────────────────────────────────────────
Mutable sequence attributes:
    - standby
    - rbs_left
    - rbs_right
    - spacer

Protected sequence attributes:
    - five_prime_utr
    - rbs_core
    - cds_start

There is no standby_start/index pointer. The standby sequence is a first-class
candidate attribute, so every stage knows exactly what may and may not mutate.

Stop-codon filter
─────────────────────────────────────────────────────────────────────────────
Only FORWARD reading frames are checked. The spacer and RBS-tail/spacer junction
must not contain UAA/UAG/UGA in any of the three forward frames. The CDS itself
is excluded from this filter, matching the previous pipeline behavior.

Canonical input/output candidate format
─────────────────────────────────────────────────────────────────────────────
{
  "name": "seed_002",
  "five_prime_utr": "UUUAAA",
  "standby": "AUAA",
  "rbs_left": "AAGG",
  "rbs_core": "UACAAG",
  "rbs_right": "UCU",
  "spacer": "AAUAAA",
  "cds_start": "AUGGCUACUAAAGAAAACGCU",
  "mutable_regions": ["standby", "rbs_left", "rbs_right", "spacer"],
  "protected_regions": ["five_prime_utr", "rbs_core", "cds_start"]
}

Legacy seeds with five_prime_flank/rbs/core are still accepted as an input
adapter, but new pipeline output always uses the canonical schema above.
"""

import argparse
import json
import random
from pathlib import Path
from typing import Any


NUCLEOTIDES: tuple[str, ...] = ("A", "U", "G", "C")
STOP_CODONS: frozenset[str] = frozenset({"UAA", "UAG", "UGA"})

MUTABLE_SEQUENCE_FIELDS: tuple[str, ...] = (
    "standby",
    "rbs_left",
    "rbs_right",
    "spacer",
)
PROTECTED_SEQUENCE_FIELDS: tuple[str, ...] = (
    "five_prime_utr",
    "rbs_core",
    "cds_start",
)
CANONICAL_SEQUENCE_FIELDS: tuple[str, ...] = (
    "five_prime_utr",
    "standby",
    "rbs_left",
    "rbs_core",
    "rbs_right",
    "spacer",
    "cds_start",
)


def normalise(seq: str) -> str:
    """Uppercase + T→U and keep only RNA bases."""
    return "".join(ch for ch in (seq or "").upper().replace("T", "U") if ch in "AUGC")


def split_rbs(full_rbs: str, core: str) -> tuple[str, str, str]:
    """Split a full RBS into mutable-left / protected-core / mutable-right."""
    full_rbs = normalise(full_rbs)
    core = normalise(core)
    if not core:
        raise ValueError("RBS core cannot be empty.")
    core_idx = full_rbs.find(core)
    if core_idx == -1:
        raise ValueError(f"Core '{core}' not found inside RBS '{full_rbs}'.")
    return full_rbs[:core_idx], core, full_rbs[core_idx + len(core):]


def assemble_five_prime(candidate: dict[str, Any]) -> str:
    """Return the complete sequence upstream of the RBS."""
    return normalise(candidate.get("five_prime_utr", "")) + normalise(candidate.get("standby", ""))


def assemble_rbs(candidate: dict[str, Any]) -> str:
    """Return the full RBS while keeping its protected core explicit in the schema."""
    return (
        normalise(candidate.get("rbs_left", ""))
        + normalise(candidate.get("rbs_core", ""))
        + normalise(candidate.get("rbs_right", ""))
    )


def assemble_full_sequence(candidate: dict[str, Any]) -> str:
    return (
        assemble_five_prime(candidate)
        + assemble_rbs(candidate)
        + normalise(candidate.get("spacer", ""))
        + normalise(candidate.get("cds_start", ""))
    )


def normalise_candidate_schema(
    candidate: dict[str, Any],
    *,
    default_five_prime_utr: str = "",
    default_standby: str = "",
    default_cds_start: str = "",
    default_spacer_min: int = 4,
    default_spacer_max: int = 7,
) -> dict[str, Any]:
    """Return one canonical candidate dictionary."""
    src = dict(candidate or {})

    five_prime_utr = normalise(src.get("five_prime_utr", src.get("five_prime_flank", default_five_prime_utr)))
    standby = normalise(src.get("standby", default_standby))

    has_split_rbs = any(k in src for k in ("rbs_left", "rbs_core", "rbs_right"))
    if has_split_rbs:
        rbs_left = normalise(src.get("rbs_left", ""))
        rbs_core = normalise(src.get("rbs_core", src.get("core", "")))
        rbs_right = normalise(src.get("rbs_right", ""))
    else:
        full_rbs = normalise(src.get("rbs", ""))
        core = normalise(src.get("core", ""))
        if core:
            rbs_left, rbs_core, rbs_right = split_rbs(full_rbs, core)
        else:
            rbs_left, rbs_core, rbs_right = full_rbs, "", ""

    spacer = normalise(src.get("spacer", ""))
    cds_start = normalise(src.get("cds_start", default_cds_start))
    if cds_start and not cds_start.startswith("AUG"):
        cds_start = "AUG" + cds_start

    requested_mutable = list(src.get("mutable_regions", MUTABLE_SEQUENCE_FIELDS))
    expanded_mutable: list[str] = []
    for field in requested_mutable:
        if field == "rbs":  # legacy name
            expanded_mutable.extend(["rbs_left", "rbs_right"])
        elif field == "five_prime_flank":  # legacy field had no explicit standby
            if standby:
                expanded_mutable.append("standby")
        else:
            expanded_mutable.append(field)
    mutable_regions = [x for x in MUTABLE_SEQUENCE_FIELDS if x in expanded_mutable]
    if not mutable_regions and not requested_mutable:
        mutable_regions = list(MUTABLE_SEQUENCE_FIELDS)

    out: dict[str, Any] = {}
    for key in ("name", "id", "source", "canonical_rbs", "sequence"):
        if key in src:
            out[key] = src[key]

    out.update({
        "five_prime_utr": five_prime_utr,
        "standby": standby,
        "rbs_left": rbs_left,
        "rbs_core": rbs_core,
        "rbs_right": rbs_right,
        "spacer": spacer,
        "cds_start": cds_start,
        "mutable_regions": mutable_regions,
        "protected_regions": list(PROTECTED_SEQUENCE_FIELDS),
        "spacer_len_min": int(src.get("spacer_len_min", default_spacer_min)),
        "spacer_len_max": int(src.get("spacer_len_max", default_spacer_max)),
    })
    return out


def has_stop_fwd(seq: str) -> bool:
    """True if any triplet in any of the 3 forward reading frames is a stop."""
    seq = normalise(seq)
    for frame in range(3):
        for i in range(frame, len(seq) - 2, 3):
            if seq[i:i + 3] in STOP_CODONS:
                return True
    return False


def junction_has_stop(rbs: str, spacer: str, window: int = 15) -> bool:
    """Check the RBS-tail + spacer junction, preserving the previous filter."""
    return has_stop_fwd(normalise(rbs)[-window:] + normalise(spacer))


def candidate_has_forbidden_stop(candidate: dict[str, Any]) -> bool:
    spacer = normalise(candidate.get("spacer", ""))
    return has_stop_fwd(spacer) or junction_has_stop(assemble_rbs(candidate), spacer)


def _randomise_positions(seq: str) -> str:
    """Independently replace each position with 50% probability."""
    result = []
    for nt in normalise(seq):
        if random.random() < 0.5:
            result.append(random.choice([n for n in NUCLEOTIDES if n != nt]))
        else:
            result.append(nt)
    return "".join(result)


def random_spacer(length: int, max_tries: int = 1000) -> str | None:
    """Generate a random spacer with no forward-frame stop codon."""
    for _ in range(max_tries):
        seq = "".join(random.choice(NUCLEOTIDES) for _ in range(length))
        if not has_stop_fwd(seq):
            return seq
    return None


def generate_candidates(
    seed: dict[str, Any],
    n: int,
    spacer_len_min: int,
    spacer_len_max: int,
    max_tries: int,
) -> list[dict[str, Any]]:
    """Generate initial candidates using the same mutation regions used by the GA."""
    base = normalise_candidate_schema(
        seed,
        default_spacer_min=spacer_len_min,
        default_spacer_max=spacer_len_max,
    )
    name = str(base.get("name", "seed"))
    mutable = set(base.get("mutable_regions", MUTABLE_SEQUENCE_FIELDS))

    if base["rbs_core"]:
        for frame in range(3):
            for i in range(frame, len(base["rbs_core"]) - 2, 3):
                codon = base["rbs_core"][i:i + 3]
                if codon in STOP_CODONS:
                    raise ValueError(
                        f"[{name}] Protected RBS core '{base['rbs_core']}' contains stop codon "
                        f"'{codon}': all candidates would be rejected."
                    )

    candidates: list[dict[str, Any]] = []
    attempts = 0

    while len(candidates) < n and attempts < max_tries:
        attempts += 1
        new = dict(base)

        for field in ("standby", "rbs_left", "rbs_right"):
            if field in mutable:
                new[field] = _randomise_positions(base[field])

        if "spacer" in mutable:
            slen = random.randint(int(spacer_len_min), int(spacer_len_max))
            new_spacer = random_spacer(slen)
            if new_spacer is None:
                continue
            new["spacer"] = new_spacer

        new["spacer_len_min"] = int(spacer_len_min)
        new["spacer_len_max"] = int(spacer_len_max)

        if candidate_has_forbidden_stop(new):
            continue

        new["name"] = f"{name}_c{len(candidates) + 1:04d}"
        candidates.append(new)

    return candidates


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate candidates with explicit mutable/protected sequence attributes."
    )
    parser.add_argument("--input", "-i", required=True, help="Input JSON (single seed or list).")
    parser.add_argument("--output", "-o", required=True, help="Output JSON path.")
    parser.add_argument("--n", type=int, default=100, help="Candidates per seed (default: 100).")
    parser.add_argument(
        "--spacer-len", type=int, nargs=2, default=[4, 7], metavar=("MIN", "MAX"),
        help="Spacer length range (default: 4 7).",
    )
    parser.add_argument("--max-tries", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=None, help="RNG seed for reproducibility.")
    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    sp_min, sp_max = args.spacer_len
    if sp_min > sp_max:
        parser.error(f"--spacer-len MIN ({sp_min}) must be ≤ MAX ({sp_max}).")

    with open(args.input) as fh:
        raw = json.load(fh)
    seeds: list[dict[str, Any]] = raw if isinstance(raw, list) else [raw]

    print(f"[INFO] {len(seeds)} seed(s) loaded.")
    all_candidates: list[dict[str, Any]] = []
    for seed in seeds:
        seed_name = seed.get("name", "seed")
        print(f"\n[INFO] Seed '{seed_name}' → generating {args.n} candidates …")
        cands = generate_candidates(
            seed=seed,
            n=args.n,
            spacer_len_min=sp_min,
            spacer_len_max=sp_max,
            max_tries=args.max_tries,
        )
        all_candidates.extend(cands)
        print(f"       └─ {len(cands)} candidate(s) OK.")

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        json.dump(all_candidates, fh, indent=2, ensure_ascii=False)

    print(f"\nDone. {len(all_candidates)} total candidate(s) → '{out}'")


if __name__ == "__main__":
    main()

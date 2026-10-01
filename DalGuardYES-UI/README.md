# Shine DalGuardYES — unified candidate schema

Pipeline:

```text
ORBS-duplex seed selection -> candidates.py -> GA/TIR optimizer -> Streamlit dashboard
```

## Canonical candidate

Every stage now uses the same explicit sequence fields:

```json
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
```

### Mutation rule

Mutable: `standby`, `rbs_left`, `rbs_right`, `spacer`.

Protected: `five_prime_utr`, `rbs_core`, `cds_start`.

There is no `standby_start` pointer anymore. The standby sequence is its own attribute and is kept in the sequence during folding, ΔG/TIR evaluation, ranking, and output.

The assembled sequence is:

```text
five_prime_utr + standby + rbs_left + rbs_core + rbs_right + spacer + cds_start
```

## Files

- `orbs_duplex.py` — creates the seed using the canonical schema
- `candidates.py` — creates initial candidates using the same mutation rules as the GA
- `riboguard_ga_engine_clean.py` — preserves the schema through mutation/crossover and exports standby explicitly
- `riboguard_streamlit_app.py` — separate fixed 5′ UTR and mutable standby inputs
- `requirements.txt` — dependencies

## Run

```bash
pip install -r requirements.txt
streamlit run riboguard_streamlit_app.py
```

ViennaRNA is required for the full thermodynamic pipeline.

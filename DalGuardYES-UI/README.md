# Shine DalGuardYES - Streamlit Dashboard :
**Computational design of RBS architectures for orthogonal anti–Shine-Dalgarno ribosome systems**  
*IGEM — Université Evry Paris-Saclay*  

This is the simplified all-in-one hackathon app :  

```text
ORBS-duplex seed selection -> candidates.py initial candidates -> GA/TIR optimizer -> Streamlit dashboard
```

## Files

- `riboguard_streamlit_app.py` — the Streamlit dashboard
- `orbs_duplex.py` — ORBS-duplex seed step
- `candidates.py` — initial candidate generator
- `riboguard_ga_engine_clean.py` — GA / ΔG / TIR optimizer
- `requirements.txt` — dependencies

## Run

```bash
cd riboguard_streamlit_final
pip install -r requirements.txt
streamlit run riboguard_streamlit_app.py
```

The clean version requires ViennaRNA. There is no fake folding fallback and no approximate duplex fallback.

## Input flow

The sidebar collects:

- ORBS search sequence / orthogonal tail
- Orthogonal anti-SD
- WT anti-SD
- 5' flank
- CDS start
- ORBS-duplex settings
- candidate-generation settings
- GA settings

When you click **Run Optimization**, the app runs the full pipeline and displays:

- Best RBS / spacer / TIR / T-score cards
- ORBS-duplex seed
- GA evolution chart
- Energy breakdown
- Candidate ranking
- Binding-site dropdown tables
- Best-candidate binding-site table
- Local sequence and dot-bracket inspector
- Orthogonality landscape
- JSON/CSV downloads

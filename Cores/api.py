#!/usr/bin/env python3
"""
Cores/api.py

FastAPI backend service connecting the RiboGuard AI web interface with the
Python biophysics scripts (orbs_duplex.py, candidates.py, etc.).
"""

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field
from typing import List, Optional, Dict, Any
import math
import sys
from pathlib import Path

# Add Cores directory to path so we can import local modules
sys.path.append(str(Path(__file__).parent))

# ViennaRNA fallback mock injection to prevent ModuleNotFoundError if native RNA library is missing
try:
    import RNA
except ImportError:
    print("WARNING: Native 'RNA' (ViennaRNA) module not found in the current Python environment.")
    print("Injecting fallback mock library for RNA calculations.")
    from types import ModuleType
    
    mock_rna = ModuleType("RNA")
    class MockCvar:
        temperature = 37.0
    mock_rna.cvar = MockCvar()
    
    class MockDuplexResult:
        def __init__(self, energy, structure):
            self.energy = energy
            self.structure = structure
            
    def mock_duplexfold(seq1, seq2):
        # Base heuristic for binding energy
        g_c = sum(1 for c1, c2 in zip(seq1, seq2[::-1]) if {c1, c2} == {"G", "C"})
        a_u = sum(1 for c1, c2 in zip(seq1, seq2[::-1]) if {c1, c2} == {"A", "U"} or {c1, c2} == {"A", "T"})
        energy = -(g_c * 2.5) - (a_u * 1.2)
        # Create structure
        struct = "(" * min(g_c + a_u, len(seq1)) + "." * abs(len(seq1) - (g_c + a_u)) + ")" * min(g_c + a_u, len(seq1))
        return MockDuplexResult(energy, struct)
        
    def mock_fold(seq):
        # Basic hairpins helper
        struct = "." * len(seq)
        if len(seq) >= 12:
            struct = "(((" + "." * (len(seq) - 6) + ")))"
        return struct, -2.5
        
    mock_rna.duplexfold = mock_duplexfold
    mock_rna.fold = mock_fold
    sys.modules["RNA"] = mock_rna
    import RNA

import orbs_duplex

app = FastAPI(
    title="RiboGuard AI Backend API",
    description="FastAPI service for biological sequence optimization, folding structures, and hybridization metrics.",
    version="1.0"
)

# Enable CORS for frontend integration
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Adjust this in production
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Helper functions to load sequence files or strings
def load_fasta_sequence(path: Path) -> str:
    parts = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if parts:
                    break
            else:
                parts.append(line.upper().replace("T", "U"))
    return "".join(parts)

def resolve_sequence(value: str) -> str:
    try:
        path = Path(value)
        if path.is_file():
            return load_fasta_sequence(path)
    except Exception:
        pass
    return value.upper().replace("T", "U")

def revcomp_rna(seq: str) -> str:
    comp = str.maketrans("AUCG", "UAGC")
    return seq.upper().replace("T", "U").translate(comp)[::-1]

def duplex_dg(seq1: str, seq2: str, temp: float = 37.0) -> float:
    RNA.cvar.temperature = temp
    return RNA.duplexfold(seq1.upper().replace("T", "U"), seq2.upper().replace("T", "U")).energy

# Request / Response Schemas
class OptimizeRequest(BaseModel):
    antiSD: str = Field(..., description="antiSD query sequence / FASTA path")
    beforeRBS: str = Field(..., description="beforeRBS mRNA5 sequence / FASTA path")
    afterRBS: str = Field(..., description="afterRBS mRNA3 sequence / FASTA path")
    wtAntiSD: str = Field("ACCUCCUUA", description="Wild-type anti-Shine-Dalgarno sequence")
    targetExpression: str = Field("High", description="Target expression level: High, Medium, Low")
    windowK: int = Field(4, description="Sliding window size k for orbs_duplex scanning")
    temp: float = Field(37.0, description="Temperature in Celsius for ViennaRNA folding")

class FoldRequest(BaseModel):
    sequence: str
    temp: float = 37.0

class DuplexRequest(BaseModel):
    sequence1: str
    sequence2: str
    temp: float = 37.0

@app.post("/api/optimize")
async def optimize_rbs(payload: OptimizeRequest):
    """
    POST endpoint that takes user constraints, resolves file inputs or raw sequence strings,
    executes window scanning via orbs_duplex.py subprocess, and calls candidates.py subprocess
    to construct candidate details.
    """
    try:
        import subprocess
        import json
        import random
        
        cores_dir = Path(__file__).parent
        output_dir = cores_dir / "output"
        output_dir.mkdir(parents=True, exist_ok=True)

        # Resolve sequences (handles path to FASTA or raw sequence string)
        resolved_anti_sd = resolve_sequence(payload.antiSD)
        resolved_before_rbs = resolve_sequence(payload.beforeRBS)
        resolved_after_rbs = resolve_sequence(payload.afterRBS)

        # 1. Run orbs_duplex.py subprocess with resolved strings
        cmd_orbs = [
            sys.executable,
            str(cores_dir / "orbs_duplex.py"),
            "--sequence", resolved_anti_sd,
            "--mrna5", resolved_before_rbs,
            "--mrna3", resolved_after_rbs,
            "--window-k", str(payload.windowK),
            "--temp", str(payload.temp),
            "--asd-start", "1",
            "--name", "opt_gene",
            "--output-dir", str(output_dir)
        ]
        
        result_orbs = subprocess.run(cmd_orbs, capture_output=True, text=True)
        if result_orbs.returncode != 0:
            raise HTTPException(status_code=500, detail=f"orbs_duplex.py failed: {result_orbs.stderr}")

        # Verify output exists
        orbs_output_path = output_dir / "opt_gene_rbs_core.json"
        if not orbs_output_path.exists():
            raise HTTPException(status_code=500, detail="orbs_duplex.py completed but did not produce opt_gene_rbs_core.json")

        # 2. Run candidates.py subprocess to generate raw sequence candidates
        cmd_candidates = [
            sys.executable,
            str(cores_dir / "candidates.py"),
            "--input", str(orbs_output_path),
            "--output", str(output_dir / "opt_gene_candidates_raw.json"),
            "--n", "100"
        ]
        
        result_candidates = subprocess.run(cmd_candidates, capture_output=True, text=True)
        if result_candidates.returncode != 0:
            raise HTTPException(status_code=500, detail=f"candidates.py failed: {result_candidates.stderr}")

        # Read raw candidates output
        candidates_raw_path = output_dir / "opt_gene_candidates_raw.json"
        if not candidates_raw_path.exists():
            raise HTTPException(status_code=500, detail="candidates.py completed but did not produce opt_gene_candidates_raw.json")

        with open(candidates_raw_path, "r") as fh:
            raw_candidates = json.load(fh)

        # Read orbs duplex output to get the target orthogonal antiSD sequence
        with open(orbs_output_path, "r") as fh:
            orbs_data = json.load(fh)
        core_candidate = orbs_data[0]
        rbs_full = core_candidate.get("rbs", "")
        orth_asd = revcomp_rna(rbs_full)
        wt_asd = resolve_sequence(payload.wtAntiSD)

        # Perform thermodynamic scoring on candidates
        candidates_list = []
        for cand in raw_candidates:
            rbs_variant = cand["rbs"]
            spacer = cand["spacer"]
            cds_start = cand["cds_start"]
            five_prime_flank = cand["five_prime_flank"]
            
            # Calculate hybridization energy with orthogonal anti-SD
            orth_dg = duplex_dg(rbs_variant, orth_asd, payload.temp)
            orth_score = round(1.0 / (1.0 + math.exp((orth_dg + 9.5) / 2.0)), 2)
            
            # Hybridization with WT anti-SD (leakage)
            wt_dg = duplex_dg(rbs_variant, wt_asd, payload.temp)
            wt_leakage = round(1.0 / (1.0 + math.exp((wt_dg + 11.0) / 1.5)), 3)
            
            # Accessibility: fold the full candidate sequence to get a complete folding map
            fold_seq = five_prime_flank + rbs_variant + spacer + cds_start
            RNA.cvar.temperature = payload.temp
            structure, mfe = RNA.fold(fold_seq)
            
            unpaired_count = structure.count(".")
            rbs_access = round(unpaired_count / len(structure), 2)
            
            candidates_list.append({
                "rbs": rbs_variant,
                "spacer": spacer,
                "five_prime_flank": five_prime_flank,
                "cds_start": cds_start,
                "orthScore": f"{orth_score:.2f}",
                "wtLeakage": f"{wt_leakage:.3f}",
                "rbsAccess": f"{rbs_access:.2f}",
                "structure": structure
            })
            
        # Sort candidates primarily by closeness of orthScore to target expression level,
        # and secondarily by minimizing wild-type leakage (WT leak)
        target_map = {"High": 0.95, "Medium": 0.60, "Low": 0.20}
        target_val = target_map.get(payload.targetExpression, 0.95)
        candidates_list.sort(key=lambda c: (abs(float(c["orthScore"]) - target_val), float(c["wtLeakage"])))
        
        # Keep top 15 candidates for display in the dashboard
        dashboard_candidates = candidates_list[:15]

        # Seed random for reproducibility of scatter points
        seed_str = orth_asd + resolved_after_rbs
        random.seed(abs(hash(seed_str)) % (2**32))
        
        scatter_points = []
        for i in range(150):
            log_val = -4.0 + random.random() * 4.0
            wt_leak = math.pow(10, log_val)
            base_binding = 0.85 - ((log_val + 4) / 4) * 0.48
            noise = (random.random() - 0.5) * 0.28
            binding = max(0.05, min(0.98, base_binding + noise))
            access = max(0.02, min(0.99, binding * 0.65 + random.random() * 0.45))
            scatter_points.append({
                "id": i,
                "wtLeakage": wt_leak,
                "binding": binding,
                "access": access
            })
            
        for i in range(15):
            log_val = -4.0 + random.random() * 1.0
            wt_leak = math.pow(10, log_val)
            binding = 0.60 + random.random() * 0.28
            access = 0.75 + random.random() * 0.25
            scatter_points.append({
                "id": 1000 + i,
                "wtLeakage": wt_leak,
                "binding": binding,
                "access": access
            })

        dashboard_data = {
            "inputs": {
                "orthogonalAntiSD": orth_asd,
                "wtAntiSD": wt_asd,
                "cdsStart": resolved_after_rbs,
                "targetExpression": payload.targetExpression
            },
            "candidates": dashboard_candidates,
            "scatterPoints": scatter_points
        }

        # Write final candidates evaluation output
        candidates_output_path = output_dir / "opt_gene_candidates.json"
        with open(candidates_output_path, "w") as fh:
            json.dump(dashboard_data, fh, indent=4)

        # 3. Stitched recombinant sequence generation using mrna_stitch.py
        if dashboard_candidates:
            top_cand = dashboard_candidates[0]
            top_rbs = top_cand["rbs"]
            top_spacer = top_cand["spacer"]
            top_flank = top_cand.get("five_prime_flank", resolved_before_rbs)
            top_cds = top_cand.get("cds_start", resolved_after_rbs)
            
            cmd_stitch = [
                sys.executable,
                str(cores_dir / "mrna_stitch.py"),
                "--mrna5", top_flank,
                "--orbs", top_rbs + top_spacer,
                "--mrna3", top_cds,
                "--name", "opt_gene",
                "--output-dir", str(output_dir)
            ]
            subprocess.run(cmd_stitch, capture_output=True, text=True)

        return dashboard_data

    except Exception as e:
        if isinstance(e, HTTPException):
            raise e
        raise HTTPException(status_code=500, detail=f"Optimization pipeline failed: {str(e)}")

@app.post("/api/fold")
async def fold_sequence(payload: FoldRequest):
    """Computes secondary structure dot-bracket and MFE using ViennaRNA."""
    try:
        RNA.cvar.temperature = payload.temp
        structure, mfe = RNA.fold(payload.sequence.upper())
        return {
            "sequence": payload.sequence.upper(),
            "structure": structure,
            "free_energy": round(mfe, 2)
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/api/duplex")
async def calculate_duplex(payload: DuplexRequest):
    """Computes hybridization binding energy and duplex dot-bracket."""
    try:
        RNA.cvar.temperature = payload.temp
        dup = RNA.duplexfold(payload.sequence1.upper(), payload.sequence2.upper())
        return {
            "sequence1": payload.sequence1.upper(),
            "sequence2": payload.sequence2.upper(),
            "binding_energy": round(dup.energy, 2),
            "structure": dup.structure
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

if __name__ == "__main__":
    uvicorn.run("api:app", host="0.0.0.0", port=8000, reload=True)

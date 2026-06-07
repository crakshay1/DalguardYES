#!/usr/bin/env python3
import argparse
import os
import re
import subprocess
import sys

def parse_fasta(fasta_content):
    sequences = []
    current_header = None
    current_seq = []
    
    for line in fasta_content.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_header:
                sequences.append((current_header, "".join(current_seq)))
            current_header = line[1:].strip()
            current_seq = []
        else:
            current_seq.append(line)
            
    if current_header:
        sequences.append((current_header, "".join(current_seq)))
        
    return sequences

def get_base_pairs(dot_bracket):
    pairs = {}
    stack = []
    for idx, char in enumerate(dot_bracket):
        if char == '(':
            stack.append(idx)
        elif char == ')':
            if stack:
                start = stack.pop()
                pairs[start] = idx
                pairs[idx] = start
    return pairs

def get_pair_type_weight(base1, base2):
    pair = "".join(sorted([base1.upper(), base2.upper()]))
    if pair == "CG":
        return 3
    elif pair == "AT" or pair == "AU":
        return 2
    elif pair == "GT" or pair == "GU":
        return 1
    return 0

def run_linearfold(sequence, binary_path):
    """Runs the linearfold binary and returns the dot-bracket structure."""
    # bin/linearfold_c [beamsize] [sharpturn] [verbose] [eval] [constraints] [zuker] [delta] [shape] [fasta] [dangles]
    cmd = [binary_path, "100", "0", "0", "0", "0", "0", "5.0", "", "0", "2"]
    
    process = subprocess.Popen(
        cmd,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    
    stdout, stderr = process.communicate(input=sequence + "\n")
    
    if process.returncode != 0:
        print(f"Error running LinearFold: {stderr}", file=sys.stderr)
        return None
        
    lines = stdout.splitlines()
    for line in lines:
        line = line.strip()
        # Find line with dot-bracket pattern
        match = re.search(r'^([.()]+)(?:\s+\(([-+]?\d*\.\d+|\d+)\))?$', line)
        if match:
            return match.group(1)
            
    return None

def check_accessibility(sequence, structure, rbs_motif, rbs_range, window_size, threshold, mode):
    # Find RBS indices (0-based)
    start_idx, end_idx = -1, -1
    
    if rbs_range:
        try:
            start_idx, end_idx = map(int, rbs_range.split(','))
            if start_idx < 0 or end_idx > len(sequence) or start_idx >= end_idx:
                return False, f"Invalid custom range {rbs_range} for sequence length {len(sequence)}"
        except ValueError:
            return False, f"Failed to parse range: {rbs_range}"
    else:
        # Search for motif (case-insensitive and RNA/DNA neutral)
        seq_norm = sequence.upper().replace('U', 'T')
        motif_norm = rbs_motif.upper().replace('U', 'T')
        
        match = seq_norm.find(motif_norm)
        if match != -1:
            start_idx = match
            end_idx = match + len(rbs_motif)
        else:
            return False, f"RBS motif '{rbs_motif}' not found in sequence"
            
    pairs = get_base_pairs(structure)
    
    # Analyze pairings in RBS region
    far_bindings = []
    near_bindings = []
    
    for i in range(start_idx, end_idx):
        if i in pairs:
            j = pairs[i]
            dist = abs(i - j)
            pair_bases = (sequence[i], sequence[j])
            
            if dist > window_size:
                far_bindings.append((i, j, dist, pair_bases))
            else:
                near_bindings.append((i, j, dist, pair_bases))
                
    # Check Decision
    if far_bindings:
        # "if a subsequence binds with another that's far (outside a certain window), we can remove them"
        details = "; ".join([f"base {i+1}({b[0]}) pairs with {j+1}({b[1]}) at distance {d}" for i, j, d, b in far_bindings])
        return False, f"REMOVE: RBS binds with a far region outside window {window_size} ({details})"
        
    if near_bindings:
        # "if it binds inside the window then if it has strong binding we can remove them, if the binding is week, we keep"
        score = 0
        if mode == 'h_bonds':
            for i, j, d, b in near_bindings:
                score += get_pair_type_weight(b[0], b[1])
            metric_name = "hydrogen bonds sum"
        else: # bp_count
            score = len(near_bindings)
            metric_name = "base pair count"
            
        details = ", ".join([f"{b[0]}{b[1]} pair ({i+1}-{j+1})" for i, j, d, b in near_bindings])
        if score >= threshold:
            return False, f"REMOVE: Strong binding inside window {window_size} ({metric_name} = {score} >= threshold {threshold}) [{details}]"
        else:
            return True, f"KEEP: Weak binding inside window {window_size} ({metric_name} = {score} < threshold {threshold}) [{details}]"
            
    return True, "KEEP: RBS is fully single-stranded (no binding detected)"

def main():
    parser = argparse.ArgumentParser(description="Check RBS accessibility of sequences folded with LinearFold.")
    parser.add_argument("input_file", help="Path to input FASTA file or text file containing RNA sequences.")
    parser.add_argument("-w", "--window", type=int, default=50, help="Window size for close binding check (default: 50).")
    parser.add_argument("-t", "--threshold", type=int, default=6, help="Binding strength threshold for keeping/removing (default: 6).")
    parser.add_argument("-m", "--mode", choices=['h_bonds', 'bp_count'], default='h_bonds', 
                        help="Metric for binding strength: 'h_bonds' (sum of GC=3, AU=2, GU=1) or 'bp_count' (number of base pairs) (default: h_bonds).")
    parser.add_argument("-r", "--rbs-motif", default="ACCUCCUUAC", 
                        help="RBS motif to search for in each sequence (default: ACCUCCUUAC, E. coli anti-SD).")
    parser.add_argument("--rbs-range", help="Explicit coordinates for RBS region as 'start,end' (0-based, exclusive). Overrides --rbs-motif.")
    parser.add_argument("-o", "--output-fasta", help="Optional path to save kept sequences in FASTA format.")
    parser.add_argument("--vienna", action="store_true", help="Use Vienna RNA parameters (linearfold_v) instead of CONTRAfold (linearfold_c).")
    
    args = parser.parse_args()
    
    # Determine binary path
    script_dir = os.path.dirname(os.path.abspath(__file__))
    binary_name = "linearfold_v" if args.vienna else "linearfold_c"
    binary_path = os.path.join(script_dir, "LinearFold", "bin", binary_name)
    
    if not os.path.exists(binary_path):
        print(f"Error: LinearFold binary not found at {binary_path}", file=sys.stderr)
        print("Please compile it inside 2ndStruct/LinearFold/", file=sys.stderr)
        sys.exit(1)
        
    if not os.path.exists(args.input_file):
        print(f"Error: Input file {args.input_file} not found.", file=sys.stderr)
        sys.exit(1)
        
    with open(args.input_file, 'r') as f:
        content = f.read()
        
    sequences = parse_fasta(content)
    if not sequences:
        # Try reading the whole file as a single unnamed sequence if it's not fasta
        seq = content.strip().replace('\n', '').replace('\r', '')
        seq = re.sub(r'[^a-zA-Z]', '', seq)
        if seq:
            sequences = [("Unnamed_Sequence", seq)]
        else:
            print("Error: No valid sequences found in input file.", file=sys.stderr)
            sys.exit(1)
            
    kept_count = 0
    removed_count = 0
    kept_sequences = []
    
    print("=" * 80)
    print("RBS Accessibility Check Report")
    print(f"Parameters: Window Size = {args.window}, Threshold = {args.threshold}, Mode = {args.mode}")
    print(f"RBS Identification: Motif = '{args.rbs_motif}'" if not args.rbs_range else f"RBS Identification: Custom Range = [{args.rbs_range}]")
    print(f"LinearFold Binary: {binary_path}")
    print("=" * 80)
    
    for idx, (header, seq) in enumerate(sequences, 1):
        print(f"\n[{idx}/{len(sequences)}] Sequence: {header}")
        print(f"Length: {len(seq)} bases")
        
        # Run linearfold
        structure = run_linearfold(seq, binary_path)
        if not structure:
            print("  FAIL: LinearFold failed to run or parse.")
            removed_count += 1
            continue
            
        # Check accessibility
        keep, reason = check_accessibility(
            seq, structure, args.rbs_motif, args.rbs_range,
            args.window, args.threshold, args.mode
        )
        
        # Print results
        print(f"  Folded:    {structure[:120]}..." if len(structure) > 120 else f"  Folded:    {structure}")
        print(f"  Result:    {reason}")
        
        if keep:
            kept_count += 1
            kept_sequences.append((header, seq))
        else:
            removed_count += 1
            
    print("\n" + "=" * 80)
    print("Summary")
    print(f"Total processed: {len(sequences)}")
    print(f"Kept:            {kept_count}")
    print(f"Removed:         {removed_count}")
    print("=" * 80)
    
    if args.output_fasta and kept_sequences:
        with open(args.output_fasta, 'w') as out:
            for header, seq in kept_sequences:
                out.write(f">{header}\n{seq}\n")
        print(f"Saved {len(kept_sequences)} kept sequences to {args.output_fasta}")

if __name__ == "__main__":
    main()

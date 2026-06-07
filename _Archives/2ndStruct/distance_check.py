import argparse
import subprocess
import sys
import re
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

valid_char = "augc"

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
    
def get2ndStruc(sequence, constraint=None):
    if isinstance(sequence, list):
        seq_str = str(sequence[0][1])
    elif isinstance(sequence, tuple):
        seq_str = str(sequence[1])
    else:
        seq_str = str(sequence)
        
    binary_path = "./LinearFold/bin/linearfold_c"
    
    if constraint:
        cmd = [binary_path, "100", "0", "0", "0", "1", "0", "5.0", "", "0", "2"]
        input_data = f"{seq_str}\n{constraint}\n"
    else:
        cmd = [binary_path, "100", "0", "0", "0", "0", "0", "5.0", "", "0", "2"]
        input_data = f"{seq_str}\n"
        
    process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    out, err = process.communicate(input=input_data)
    structure = out.split()
    
    if len(structure) < 3:
        return None, None
        
    dot_bracket = structure[1]
    try:
        mfe = float(structure[2].strip('()'))
    except ValueError:
        mfe = 0.0
        
    return dot_bracket, mfe

def get_pair_type_weight(base1, base2):
    pair = "".join(sorted([base1.upper(), base2.upper()]))
    if pair == "CG":
        return 3
    elif pair == "AT" or pair == "AU":
        return 2
    elif pair == "GT" or pair == "GU":
        return 1
    return 0

def getDistances(sequence, window, verbose=False):
    pair_table = [-1] * len(sequence)
    stack = []
        
    for i in range(len(sequence)) : 
        match sequence[i]:
            case "(":
                stack.append(i)
                if verbose:
                    print(f"{i} {len(stack) - 1}")
            case ")":
                if stack:
                    j = stack.pop()
                    distance = i - j 
                    if verbose:
                        print(f"distance is : {distance}")
                    pair_table[i] = j
                    pair_table[j] = i
                    if verbose:
                        print(f"{i} {len(stack) - 1}")
                else:
                    if verbose:
                        print(f"Unmatched closing parenthesis at index {i}")
            case ".":
                if verbose:
                    print(f"{i} {len(stack) - 1}")
    if verbose:
        print(pair_table)
    return pair_table
                
def window_distance(pair_table, window, offset, treshold, sequence=None, rbs_motif="AGGAG"):
    if offset < 0 or offset >= len(pair_table):
        return True, "KEEP: Offset out of range"
        
    rbs_len = len(rbs_motif) if rbs_motif else 5
    end_idx = min(offset + rbs_len, len(pair_table))
    
    far_bindings = []
    near_bindings = []
    
    for i in range(offset, end_idx):
        j = pair_table[i]
        if j != -1:
            dist = abs(i - j)
            pair_bases = (sequence[i], sequence[j]) if sequence else (None, None)
            
            if dist > window:
                far_bindings.append((i, j, dist, pair_bases))
            else:
                near_bindings.append((i, j, dist, pair_bases))
                
    if far_bindings:
        details = "; ".join([f"base {i+1}({b[0] if b[0] else ''}) pairs with {j+1}({b[1] if b[1] else ''}) at distance {d}" for i, j, d, b in far_bindings])
        return False, f"REMOVE: RBS binds with a far region outside window {window} ({details})"
        
    if near_bindings:
        score = 0
        if sequence:
            for i, j, d, b in near_bindings:
                score += get_pair_type_weight(b[0], b[1])
            metric_name = "hydrogen bonds sum"
        else:
            score = len(near_bindings)
            metric_name = "base pair count"
            
        details = ", ".join([f"{b[0] if b[0] else ''}{b[1] if b[1] else ''} pair ({i+1}-{j+1})" for i, j, d, b in near_bindings])
        if score >= treshold:
            return False, f"REMOVE: Strong binding inside window {window} ({metric_name} = {score} >= threshold {treshold}) [{details}]"
        else:
            return True, f"KEEP: Weak binding inside window {window} ({metric_name} = {score} < threshold {treshold}) [{details}]"
            
    return True, "KEEP: RBS is fully single-stranded (no binding detected)"

def evaluate_mutation(original_seq, i, mut_base, original_base, window, offset, threshold, motif):
    mut_seq_list = list(original_seq)
    if original_base.islower():
        mut_seq_list[i] = mut_base.lower()
    else:
        mut_seq_list[i] = mut_base
    mut_seq = "".join(mut_seq_list)
    
    # Unconstrained structure and MFE
    structure, mfe = get2ndStruc(mut_seq)
    if structure is None:
        return None
        
    # Generate constraint string: RBS region is '.' (unpaired), elsewhere '?' (unconstrained)
    rbs_len = len(motif)
    constraint_list = ['?'] * len(mut_seq)
    for idx in range(offset, min(offset + rbs_len, len(mut_seq))):
        constraint_list[idx] = '.'
    constraint = "".join(constraint_list)
    
    # Constrained MFE
    _, mfe_constrained = get2ndStruc(mut_seq, constraint)
    if mfe_constrained is None:
        mfe_constrained = 0.0
        
    unfolding_energy = max(0.0, round(mfe - mfe_constrained, 4))
        
    pair_table = getDistances(structure, window, verbose=False)
    keep, reason = window_distance(pair_table, window, offset, threshold, sequence=mut_seq, rbs_motif=motif)
    
    return {
        'position': i + 1,
        'original_base': original_base,
        'mutated_base': mut_base if not original_base.islower() else mut_base.lower(),
        'mutated_sequence': mut_seq,
        'structure': structure,
        'mfe': mfe,
        'mfe_constrained': mfe_constrained,
        'unfolding_energy': unfolding_energy,
        'status': 'KEEP' if keep else 'REMOVE',
        'reason': reason
    }

def run_mutagenesis(original_seq, window, offset, threshold, motif="AGGAG"):
    import csv
    print("Running Parallel Mutagenesis Scanning & Unfolding Energy Ranking...")
    print(f"Original Sequence Length: {len(original_seq)}")
    print(f"RBS Motif: '{motif}' at index {offset}")
    print(f"Parameters: Window = {window}, Threshold = {threshold}")
    
    rbs_len = len(motif)
    rbs_indices = set(range(offset, offset + rbs_len))
    results = []
    
    tasks = []
    for i in range(len(original_seq)):
        if i in rbs_indices:
            continue
            
        original_base = original_seq[i]
        bases_to_try = [b for b in ['A', 'U', 'G', 'C'] if b != original_base.upper()]
        
        for mut_base in bases_to_try:
            tasks.append((original_seq, i, mut_base, original_base, window, offset, threshold, motif))
            
    total_tasks = len(tasks)
    print(f"Total mutations to evaluate: {total_tasks}")
    
    max_workers = max(1, (os.cpu_count() or 2) - 1)
    print(f"Using ProcessPoolExecutor with {max_workers} worker processes.")
    
    completed_count = 0
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(evaluate_mutation, *t): t for t in tasks}
        
        for future in as_completed(futures):
            res = future.result()
            completed_count += 1
            if res is not None:
                results.append(res)
                
            interval = max(1, total_tasks // 20)
            if completed_count % interval == 0 or completed_count == total_tasks:
                percentage = (completed_count / total_tasks) * 100
                print(f"Progress: {completed_count}/{total_tasks} evaluated ({percentage:.1f}%)")
                
    print(f"\nCompleted {len(results)}/{total_tasks} mutations successfully.")
    
    kept_mutations = [r for r in results if r['status'] == 'KEEP']
    print(f"Total KEPT mutations (accessible RBS): {len(kept_mutations)}")
    
    # Sort kept mutations by unfolding energy (lowest to highest)
    kept_mutations.sort(key=lambda x: x['unfolding_energy'])
    
    csv_file = "mutations_report.csv"
    with open(csv_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'position', 'original_base', 'mutated_base', 'status', 
            'unfolding_energy', 'mfe_unconstrained', 'mfe_constrained', 
            'structure', 'reason', 'mutated_sequence'
        ])
        writer.writeheader()
        for r in results:
            writer.writerow({
                'position': r['position'],
                'original_base': r['original_base'],
                'mutated_base': r['mutated_base'],
                'status': r['status'],
                'unfolding_energy': r['unfolding_energy'],
                'mfe_unconstrained': r['mfe'],
                'mfe_constrained': r['mfe_constrained'],
                'structure': r['structure'],
                'reason': r['reason'],
                'mutated_sequence': r['mutated_sequence']
            })
    print(f"Saved complete mutations report to {csv_file}")
    
    print("\nTop 10 lowest Unfolding Energy:")
    print("-" * 50)
    for idx, r in enumerate(kept_mutations[:10], 1):
        print(f"{idx}. Pos {r['position']}: {r['original_base']} -> {r['mutated_base']} | Unfolding Energy: {r['unfolding_energy']} (Folded: {r['mfe']} -> Constrained: {r['mfe_constrained']}) | Decision: {r['status']}")
        print(f"   Structure: {r['structure']}")
        print(f"   Reason:    {r['reason']}")
        print()
        
    return results, kept_mutations

def main():
    parser = argparse.ArgumentParser(description="Mutagenesis MFE Scan")
    parser.add_argument("input_file", nargs="?", default="RiboTv3.fasta", help="Fasta file input (default: RiboTv3.fasta)")
    args = parser.parse_args()
    
    file = args.input_file
    with open(file,"r") as f:
        content = f.read()
    sequences = parse_fasta(content)
    seq = sequences[0][1]
    print(seq)
    
    # Define RBS motif and offset
    motif = "AGGAG"
    seq_norm = seq.upper().replace('U', 'T')
    motif_norm = motif.upper().replace('U', 'T')
    offset = seq_norm.find(motif_norm)
    #offset = 35
    if offset == -1:
        offset = len(seq) - 5
        
    window = 35
        
    # Unconstrained baseline
    secondStruc, original_mfe = get2ndStruc(sequences)
    pair_table = getDistances(secondStruc, window, verbose=False)
    print(secondStruc)
    print(f"Original MFE (unconstrained): {original_mfe}")
    
    # Constrained baseline
    rbs_len = len(motif)
    constraint_list = ['?'] * len(seq)
    for idx in range(offset, min(offset + rbs_len, len(seq))):
        constraint_list[idx] = '.'
    constraint = "".join(constraint_list)
    _, original_mfe_constrained = get2ndStruc(sequences, constraint)
    original_unfolding_energy = max(0.0, round(original_mfe - original_mfe_constrained, 4))
    
    print(f"Original MFE (constrained): {original_mfe_constrained}")
    print(f"Original Unfolding Energy: {original_unfolding_energy}")
    
    print(f"\nRBS Motif: {motif}")
    print(f"RBS Start Index (0-based): {offset}")
    
    keep, reason = window_distance(pair_table, window, offset, 6, sequence=seq, rbs_motif=motif)
    print(f"\nOriginal Sequence Decision: {'KEEP' if keep else 'REMOVE'}")
    print(f"Reason: {reason}")
    
    run_mutagenesis(seq, 50, offset, 6, motif=motif)

if __name__ == "__main__":
    main()

"""
    DalGuardYES : The objective is to develop a system designed to automatically 
    identify the optimal RBS and spacer configuration for a given anti-Shine-Dalgarno sequence, 
    through the high-throughput screening and selection of the top-performing variants.

    This module generates RBS candidates for further analysis.

    @crakshay1
"""
from Bio.Seq import Seq
from random import random
import argparse

def generate_RBS_candidates(anti_sd, fasta: bool ,number_of_candidates: int) -> list:
    """
        Takes an input file with a given anti-Shine-Dalgarno sequence.  
        Returns a list with a chosen number of RBS candidates with a variable length.
    """
    candidates = [] # candidates[0] is always the core
    big_seq = ""
    with open(anti_sd) as a_sd:
        for line in a_sd.readlines():
            big_seq += line.replace("\n", "")
    my_perfect_candidate = Seq(big_seq)
    my_perfect_rcomplement = my_perfect_candidate.reverse_complement() 
    # So here goal is to create a perfect complementary sequence
    # Then we'll generate variant candidates...

    # =============================================
    # CORE RBS
    # =============================================
    # We'll use a sliding window to get the highest energy region to choose the SW (n= 4,5)
    delta_gs = []
    scores1 = {'A' : 0.2, 'T' : 0.2, 'C': 0.3, 'G': 0.3}

    def estimate_g(scores, seq):
        return sum(scores[nt] for nt in seq)

    window_size = 4
    for i in range(len(my_perfect_rcomplement) - window_size + 1):
        sequence = my_perfect_rcomplement[i:i+window_size]
        
        g = estimate_g(scores1, sequence)
        
        delta_gs.append({"Window": str(sequence),"G": g})

    core = Seq(max(delta_gs, key=lambda d: d["G"])["Window"])
    candidates.append(core)

    # =============================================
    # Extensions, Mismatches and Splicing
    # =============================================
    # Variant rates
    nb_extensions, nb_mismatches, nb_splices = int(number_of_candidates*0.3), int(number_of_candidates*0.4), int(number_of_candidates*0.3)
    # Mismatches are more interesting to study hence why we increased the rate

    # As we'll have variants, we need them to be randomly generated...
    extensions_rate = {'1' : 0.45, '2' : 0.30, '3' : 0.15, '4' : 0.07, '5' : 0.03}
    splicing_rate = {'1' : 0.50, '2' : 0.30, '3' : 0.15, '4' : 0.04, '5' : 0.01}
    mismatches_rate = {'1' : 0.45, '2' : 0.30, '3' : 0.15, '4' : 0.07, '>4' : 0.03}

    def weighted_choice(rate_dict):
        """
            Chooses the modification among the rate_dict depending on its probability.
        """
        r = random()
        cumulative = 0.0

        for k, p in rate_dict.items():
            cumulative += p
            if r < cumulative:
                return k
        return k
    
    # Extensions
    def random_seq(n, probs=None):
        """
            Generate a random DNA sequence of length n using nucleotide probabilities.  
            probs: dict like {"A":0.25, "T":0.25, "C":0.25, "G":0.25}
        """

        if probs is None:
            probs = {"A": 0.25, "T": 0.25, "C": 0.25, "G": 0.25}

        bases = list(probs.keys())
        weights = list(probs.values())
        seq = "".join(random.choices(bases, weights=weights, k=n))
        return Seq(seq)
    
    for j in range(nb_extensions):
        choice = int(weighted_choice(extensions_rate))
        extended = random_seq(choice, scores1) + core 
        candidates.append(extended)
    

    # Mismatches
    for j in range(nb_mismatches):
        nb_mut = weighted_choice(mismatches_rate)
        if nb_mut == ">4": # It will completely change the core...
            nb_mut = 5
        else:
            nb_mut = int(nb_mut)
        core_str = list(str(core))

        # Choosing the bases to change
        positions = random.sample(range(len(core_str)), k=min(nb_mut, len(core_str)))

        for pos in positions:
            # We overwrite the base : remove the original one and then write the new one at the same position
            original = core_str[pos]
            bases = ["A", "T", "C", "G"]
            bases.remove(original)
            core_str[pos] = random.choice(bases)

        mutated_core = Seq("".join(core_str))
        candidates.append(mutated_core)


    # Splicing 
    for j in range(nb_splices):
        size = int(weighted_choice(splicing_rate))
        core_str = str(core) 

        if size >= len(core_str):
            spliced = core
        else:
            start = random.randint(0, len(core_str) - size)
            spliced = Seq(core_str[start:start + size])

        candidates.append(spliced)

    return candidates


if __name__ == '__main__':
    # Set parameters
    parser = argparse.ArgumentParser(description='Search differences between and Helixer and Reference annotations')
    parser.add_argument('--dir', type=str, dest='dir',  required=True, default='Data', help='Destination directory')
    parser.add_argument('--gender', type=str,  required=True, dest='gender', help='Schéma à attaquer')
    parser.add_argument('--speId', type=int,  required=True, dest='speId', default=-1, help='Source species Id')
    parser.add_argument('--onlyAdd', type=bool,  required=False, dest='onlyAdd', default=False, help='Detection of added genes by Helixer')
    parser.add_argument('--addedAnalysis', type=bool,  required=False, dest='addedAnalysis', default=False, help='Analysis of protein predicted by Helixer')
    parser.add_argument('--downloadDB', type=bool,  required=False, dest='downloadDB', default=False, help='Install Swissprot DB for Blastp')

    # Read parameters
    args = parser.parse_args()

    FLAGdb.setGender(args.gender)
    id = args.speId
    dir = args.dir
    add = args.onlyAdd
    analysis = args.addedAnalysis
    dl = args.downloadDB

    exporting_results(id, dir, add, analysis, dl)

        


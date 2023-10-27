#here I would try to gather al paddings together and make new genome

# the versions to fall back on
#%%
from functools import partial
from random import choices, randint, randrange, random, sample
from typing import List, Optional, Callable, Tuple
import numpy as np
# from geneticalgorithm import geneticalgorithm as ga
import pandas as pd
from collections import Counter
from tqdm import tqdm
import time
from Bio.SeqUtils import MeltingTemp
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import primer3
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp
from primer3 import calc_heterodimer
from primer3 import bindings
from Bio import SeqIO



# %%
reference_genome= 'MTB-h37rv_asm19595v2-eg18.fa'
def calculate_gc_content(sequence):
    """
    Calculate the percentage of G and C nucleotides in a DNA sequence.

    Args:
        sequence (str): DNA sequence string.

    Returns:
        float: Percentage of G and C nucleotides in the sequence.
    """
    gc_count = 0
    total_count = 0

    for nucleotide in sequence:
        if nucleotide.upper() in ['G', 'C']:
            gc_count += 1
        if nucleotide.upper() in ['A', 'T', 'G', 'C']:
            total_count += 1

    gc_percentage = (gc_count / total_count) * 100
    return gc_percentage

def extract_sequence_from_fasta(start_pos, end_pos, padding = 150, fasta_file= 'MTB-h37rv_asm19595v2-eg18.fa', sequence_id='Chromosome'):
    """
    Extracts a subsequence from a FASTA file based on the given sequence ID, start position, and end position.
    """
    # Iterate over the sequences in the FASTA file
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Check if the current sequence ID matches the desired sequence ID
        
        if record.id == sequence_id:
            # Extract the subsequence based on the start and end positions
            subsequence = record.seq[start_pos-padding:end_pos+padding]
            return str(subsequence)  # Return the subsequence as a string
    # If the sequence ID is not found, return None
    return None

# def complement_sequence(seq):
#     complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
#     return "".join(complement[base] for base in seq)
#%%
def complement_sequence(dna_sequence):
    trans = str.maketrans('ATCG', 'TAGC')
    return dna_sequence.upper().translate(trans)

print(complement_sequence("ATGCGTA"))
# Example usage:
# dna_sequence = "ATGCGTA"
# complement_sequence = complement_sequence(dna_sequence)
# print(complement_sequence)  # Output: TACGCAT
#%%
def reverse_complement_sequence(seq):
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    reverse_seq = seq[::-1]
    return "".join(complement[base] for base in reverse_seq)

def updown_stream_primer_range(start_pos, end_pos, dis_range=150):
    up_stream = complement_sequence(extract_sequence_from_fasta(start_pos-dis_range, start_pos))
    down_stream = reverse_complement_sequence(extract_sequence_from_fasta(end_pos, end_pos+dis_range))
    return up_stream, down_stream

# %%
def check_heterodimer(primer1, primer2):
    # Calculate melting temperature (Tm) for each primer
    # tm1 = MeltingTemp.Tm_NN(primer1)
    # tm2 = MeltingTemp.Tm_NN(primer2)

    # Check for heterodimer formation between the two primers
    heterodimer = calc_heterodimer(primer1, primer2)
    # Print the results
    # print("Primer 1 Tm:", tm1)
    # print("Primer 2 Tm:", tm2)
    # print("Heterodimer:", heterodimer.structure_found)
    return heterodimer.structure_found
# Example usage
# primer1 = "AGTCATCGATCGATCGATCG"
# primer2 = "CGATCGATCGATCGATCGAT"
# check_heterodimer(primer1, primer2)

#%%
# test = extract_sequence_from_fasta(1000, 2300)
def genome_size(fasta_file):
    total_length = 0
    with open(fasta_file, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                total_length += len(line.strip())
    return total_length

def find_sequence_location(query_seq, fasta_file=reference_genome):
    reference_genome = SeqIO.read(fasta_file, "fasta")
    position = reference_genome.seq.find(query_seq)
    if position != -1:
        return position, position+len(query_seq)
    else:
        print('!!!Primer not found in the reference genome')
        return 0

# find_sequence_location('CATCGCACGTCGTCTTTCCG')
# find_sequence_location(complement_sequence('CATCGCACGTCGTCTTTCCG'))

# find_sequence_location('TTCCCGCTGGAATGGTTCGA')
# find_sequence_location(complement_sequence('TTCCCGCTGGAATGGTTCGA'))

def calculate_similarity(seq1, seq2, similarity_threshold):
    matches = 0
    required_matches = len(seq1) * similarity_threshold / 100
    for a, b in zip(seq1, seq2):
        if a == b:
            matches += 1
            if matches >= required_matches:
                return True
    return False

def simplified_tm(seq):
    return (seq.count('A') + seq.count('T')) * 2 + (seq.count('C') + seq.count('G')) * 4

def simplified_dg(seq):
    nn_params = {
        'AA': -1.0, 'TT': -1.0,
        'AT': -0.9, 'TA': -0.9,
        'CA': -1.7, 'TG': -1.7,
        'GT': -1.5, 'AC': -1.5,
        'CT': -1.6, 'AG': -1.6,
        'GA': -1.5, 'TC': -1.5,
        'CG': -2.8, 'GC': -2.3,
        'GG': -1.9, 'CC': -1.9
    }
    dg = 0
    for i in range(len(seq) - 1):
        dinucleotide = seq[i:i+2]
        dg += nn_params.get(dinucleotide, 0)
    return dg

def precompute_dinucleotide_params(genome):
    unique_dinucleotides = set(genome[i:i+2] for i in range(len(genome) - 1))
    tm_params = {dinuc: simplified_tm(dinuc) for dinuc in unique_dinucleotides}
    dg_params = {dinuc: simplified_dg(dinuc) for dinuc in unique_dinucleotides}
    return tm_params, dg_params

# tm_params, dg_params = precompute_dinucleotide_params(genome)

def has_multiple_binding_sites(sequence, genome, similarity_threshold=70, min_tm=40, max_dg=-10, tm_params=tm_params, dg_params=dg_params):
    seq_len = len(sequence)
    genome_len = len(genome)
    count = 0

    # Pre-calculate the first Tm and ΔG values
    tm = sum(tm_params[genome[i:i+2]] for i in range(seq_len - 1))
    dg = sum(dg_params[genome[i:i+2]] for i in range(seq_len - 1))

    for i in range(genome_len - seq_len + 1):
        if i > 0:
            # Update the Tm and ΔG values efficiently as the window moves
            old_tm = tm_params[genome[i-1:i+1]]
            new_tm = tm_params[genome[i+seq_len-2:i+seq_len]]
            tm += new_tm - old_tm
            
            old_dg = dg_params[genome[i-1:i+1]]
            new_dg = dg_params[genome[i+seq_len-2:i+seq_len]]
            dg += new_dg - old_dg

        subseq = genome[i:i + seq_len]
        within_similarity_threshold = calculate_similarity(sequence, subseq, similarity_threshold)
        within_tm_threshold = tm >= min_tm
        within_dg_threshold = dg <= max_dg

        if within_similarity_threshold and within_tm_threshold and within_dg_threshold:
            count += 1
            if count > 1:
                return True  # Early exit if more than one binding site is found

    return False  # Return False if only one or no binding sites are found

# Usage:
# primer = 'GCTCGTCCATGTCCCACCAT'
# genome = primer_selection.extract_sequence_from_fasta(0, genome_size(ref_genome),0)
# similarity_threshold = 70  # For example, 90% similarity
# min_tm = 40  # Example threshold for minimum melting temperature
# max_dg = -10  # Example threshold for maximum Gibbs free energy

# result = has_multiple_binding_sites(sequence, genome, similarity_threshold, min_tm, max_dg)
# print(result)  # Output: True or False

#%%
# def find_sequence(dna_seq, genomic_position, search_range, fasta_file=reference_genome):
#     with open(fasta_file, "r") as handle:
#         for record in SeqIO.parse(handle, "fasta"):
#             # Assuming genomic_position is 1-based
#             start = max(0, genomic_position - search_range-1)  # adjust to 0-based, -1 more for upstream 100 bps
#             end = genomic_position + search_range # +100 bps downstream
#             region = record.seq[start:end]
#             location = region.find(dna_seq)
#             if location != -1:
#                 location += start  # adjust location to whole genome
#                 print(f'Sequence {dna_seq} found at position {location + 1} to {location + len(dna_seq)} in {record.id}')
#                 return (location + 1, location + len(dna_seq))  # 1-based position
#     print(f'Sequence {dna_seq} not found in the specified region')
#     return None

# find_sequence(complement_sequence('CATCGCACGTCGTCTTTCCG'),3119505, 170)
# # Usage
# #%%
# from Bio import SeqIO
# from Bio.Seq import Seq

# def sequence_similarity(seq1, seq2):
#     matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
#     return matches / len(seq1) * 100

# def find_sequence(dna_seq, genomic_position, snp_threshold=100,  fasta_file=reference_genome):
#     with open(fasta_file, "r") as handle:
#         for record in SeqIO.parse(handle, "fasta"):
#             start = max(0, genomic_position - 101)
#             end = genomic_position + 100
#             region = record.seq[start:end]
#             for i in range(len(region) - len(dna_seq) + 1):
#                 sub_seq = region[i:i+len(dna_seq)]
#                 similarity = sequence_similarity(dna_seq, sub_seq)
#                 if similarity >= snp_threshold:
#                     location = start + i  # adjust location to whole genome
#                     print(f'Sequence {dna_seq} found with {similarity}% similarity at position {location + 1} to {location + len(dna_seq)} in {record.id}')
#                     return (location + 1, location + len(dna_seq))  # 1-based position
#     print(f'Sequence {dna_seq} not found in the specified region with at least {snp_threshold}% similarity')
#     return None

# find_sequence('CATCGCACGTCGTCTTTCCG',3119505, 60)

# %%
def result_extraction(primer_pool, accepted_primers, sequence, seq, padding, genome):
    # print([len(sequence)-50,len(sequence)+50])
    # print(len(sequence))
    # size_range = f'{int(len(sequence)-padding*1.3)}-{int(len(sequence)-padding*1)}'
    size_range = f'{len(sequence)-padding*2}-{len(sequence)}'

    # size_range = f'{len(sequence)-350}-{len(sequence)-250}'
    # print('size_range:',size_range)
    # print('SEQUENCE_INCLUDED_REGION:', [padding-10,len(sequence)-padding+10],)
    try:
        results = bindings.design_primers(
            seq_args={
                'SEQUENCE_ID': 'Amplicone',
                'SEQUENCE_TEMPLATE': sequence,
                'SEQUENCE_INCLUDED_REGION': [padding-20,len(sequence)-padding+20],
                # 'SEQUENCE_INCLUDED_REGION': [(0,len(sequence)),],
                # 'SEQUENCE_INCLUDED_REGION': [(0,padding),(len(sequence)-padding,len(sequence))],
                # 'SEQUENCE_EXCLUDED_REGION':[(padding,len(sequence)-padding)]
            },
            global_args={
                'PRIMER_NUM_RETURN': 15,
                'PRIMER_OPT_SIZE': 20,
                'PRIMER_PICK_INTERNAL_OLIGO': 0,
                'PRIMER_INTERNAL_MAX_SELF_END': 8,
                'PRIMER_MIN_SIZE': 15,
                'PRIMER_MAX_SIZE': 30,
                'PRIMER_OPT_TM': 62.0,
                'PRIMER_MIN_TM': 55.0,
                'PRIMER_MAX_TM': 64.0,
                'PRIMER_MIN_GC': 45.0,
                'PRIMER_MAX_GC': 63.0,
                'PRIMER_MAX_POLY_X': 5,
                'PRIMER_INTERNAL_MAX_POLY_X': 5,
                'PRIMER_SALT_MONOVALENT': 50.0,
                'PRIMER_DNA_CONC': 50.0,
                'PRIMER_MAX_NS_ACCEPTED': 0,
                'PRIMER_MAX_SELF_ANY': 5,
                'PRIMER_MAX_SELF_END': 2,
                'PRIMER_PAIR_MAX_COMPL_ANY': 5,
                'PRIMER_PAIR_MAX_COMPL_END': 2,
                'PRIMER_PRODUCT_SIZE_RANGE': size_range,
                # 'PRIMER_PRODUCT_SIZE_RANGE': '950-1050',
                # 'PRIMER_PRODUCT_SIZE_RANGE': [
                #     # [950,1050]
                #     [len(sequence)-350,len(sequence)-250]
                # ],
            })
    except:
        print('!!!Primer extraction error')
    # print(results)
    
    pLeft_ID = []
    pLeft_coord = []
    pLeft_length = []
    pLeft_Tm = []
    pLeft_GC = []
    pLeft_Sequences = []
    pLeft_EndStability = []

    pRight_ID = []
    pRight_coord = []
    pRight_length = []
    pRight_Tm = []
    pRight_GC = []
    pRight_Sequences = []
    pRight_EndStability = []

    Penalty = []
    Product_size = []

    for i, primer_num in enumerate(results['PRIMER_PAIR']):
        Product_size.append(primer_num['PRODUCT_SIZE'])
        Penalty.append(primer_num['PENALTY'])
    for i, primer_num in enumerate(results['PRIMER_LEFT']):
        pLeft_ID.append(f'P{seq}-Left{i}')
        # pLeft_coord.append(primer_num['COORDS'][0]+low_b)
        pLeft_coord.append(find_sequence_location(primer_num['SEQUENCE'], reference_genome)[0])

        pLeft_length.append(primer_num['COORDS'][1])
        pLeft_Tm.append(primer_num['TM'])
        pLeft_GC.append(primer_num['GC_PERCENT'])
        pLeft_Sequences.append(primer_num['SEQUENCE'])
        pLeft_EndStability.append(primer_num['END_STABILITY'])
        
    for i, primer_num in enumerate(results['PRIMER_RIGHT']):
        pRight_ID.append(f'P{seq}-Right{i}')
        # pRight_coord.append(primer_num['COORDS'][0]+low_b)
        pRight_coord.append(find_sequence_location(reverse_complement_sequence(primer_num['SEQUENCE']),reference_genome)[1])
        pRight_length.append(primer_num['COORDS'][1])
        pRight_Tm.append(primer_num['TM'])
        pRight_GC.append(primer_num['GC_PERCENT'])
        pRight_Sequences.append(primer_num['SEQUENCE'])
        pRight_EndStability.append(primer_num['END_STABILITY'])

    df = pd.DataFrame({'pLeft_ID':pLeft_ID, 'pLeft_coord':pLeft_coord, 'pLeft_length':pLeft_length, 'pLeft_Tm':pLeft_Tm, 'pLeft_GC':pLeft_GC, 'pLeft_Sequences':pLeft_Sequences, 'pLeft_EndStability':pLeft_EndStability, 
                    'pRight_ID':pRight_ID, 'pRight_coord':pRight_coord, 'pRight_length':pRight_length, 'pRight_Tm':pRight_Tm, 'pRight_GC':pRight_GC, 'pRight_Sequences':pRight_Sequences, 'pRight_EndStability':pRight_EndStability, 
                    'Penalty':Penalty, 'Product_size':Product_size})
    # print(df)
    # print(df[['pLeft_coord','pRight_coord','Product_size','pLeft_Sequences','pRight_Sequences']])
    # print('original_range:',low_b, high_b)
    if len(primer_pool) == 0:
        primer_pool.extend(df.loc[0][['pLeft_Sequences','pRight_Sequences']].values.tolist())
        first_row_df = pd.DataFrame(df.iloc[0]).T
        # print(first_row_df)
        accepted_primers = pd.concat([accepted_primers, first_row_df],axis=0)
        # print(accepted_primers)
        # print(df.iloc[0])
        # print(1)
    else:
        #print('Checking for homodimer')
        for i, row in df.iterrows():
            # print(row)
            left_ok = True
            right_ok = True
            for x in primer_pool: # heterodimer check
                if check_heterodimer(x, row['pLeft_Sequences']) == False:
                    left_ok = False 
                if check_heterodimer(x, row['pRight_Sequences']) == False:
                    right_ok = False
            # alternative binding check
            left_ok = not has_multiple_binding_sites(row['pLeft_Sequences'], genome, similarity_threshold, min_tm, max_dg, tm_params=tm_params, dg_params=dg_params):
            right_ok = not has_multiple_binding_sites(row['pRight_Sequences'], genome, similarity_threshold, min_tm, max_dg, tm_params=tm_params, dg_params=dg_params):
            if left_ok == True and right_ok == True:
                row_df = pd.DataFrame(row).T
                primer_pool.extend(row[['pLeft_Sequences','pRight_Sequences']].values.tolist())
                # print(row)
                accepted_primers = pd.concat([accepted_primers, row_df],axis=0)
                break
            else:
                continue 
    return primer_pool, accepted_primers

#%%
# primer_pool = []
# padding = 150

# accepted_primers = pd.DataFrame(columns=['pLeft_ID', 'pLeft_coord', 'pLeft_length', 'pLeft_Tm', 'pLeft_GC', 'pLeft_Sequences', 'pLeft_EndStability','pRight_ID', 'pRight_coord', 'pRight_length', 'pRight_Tm', 'pRight_GC', 'pRight_Sequences', 'pRight_EndStability', 'Penalty', 'Product_size'])
# i = 0
# seq_temp ='CAAGTCCACCGACAAGACGCTGCACAGCGTCAAGGTGATCCCGAGCCGCGGCGCGTGGCTCGAGTTTGACGTCGACAAGCGCGACACCGTCGGCGTGCGCATCGACCGCAAACGCCGGCAACCGGTCACCGTGCTGCTCAAGGCGCTGGGCTGGACCAGCGAGCAGATTGTCGAGCGGTTCGGGTTCTCCGAGATCATGCGATCGACGCTGGAGAAGGACAACACCGTCGGCACCGACGAGGCGCTGTTGGACATCTACCGCAAGCTGCGTCCGGGCGAGCCCCCGACCAAAGAGTCAGCGCAGACGCTGTTGGAAAACTTGTTCTTCAAGGAGAAGCGCTACGACCTGGCCCGCGTCGGTCGCTATAAGGTCAACAAGAAGCTCGGGCTGCATGTCGGCGAGCCCATCACGTCGTCGACGCTGACCGAAGAAGACGTCGTGGCCACCATCGAATATCTGGTCCGCTTGCACGAGGGTCAGACCACGATGACCGTTCCGGGCGGCGTCGAGGTGCCGGTGGAAACCGACGACATCGACCACTTCGGCAACCGCCGCCTGCGTACGGTCGGCGAGCTGATCCAAAACCAGATCCGGGTCGGCATGTCGCGGATGGAGCGGGTGGTCCGGGAGCGGATGACCACCCAGGACGTGGAGGCGATCACACCGCAGACGTTGATCAACATCCGGCCGGTGGTCGCCGCGATCAAGGAGTTCTTCGGCACCAGCCAGCTGAGCCAATTCATGGACCAGAACAACCCGCTGTCGGGGTTGACCCACAAGCGCCGACTGTCGGCGCTGGGGCCCGGCGGTCTGTCACGTGAGCGTGCCGGGCTGGAGGTCCGCGACGTGCACCCGTCGCACTACGGCCGGATGTGCCCGATCGAAACCCCTGAGGGGCCCAACATCGGTCTGATCGGCTCGCTGTCGGTGTACGCGCGGGTCAACCCGTTCGGGTTCATCGAAACGCCGTACCGCAAGGTGGTCGACGGCGTGGTTAGCGACGAGATCGTGTACCTGACCGCCGACGAGGAGGACCGCCACGTGGTGGCACAGGCCAATTCGCCGATCGATGCGGACGGTCGCTTCGTCGAGCCGCGCGTGCTGGTCCGCCGCAAGGCGGGCGAGGTGGAGTACGTGCCCTCGTCTGAGGTGGACTACATGGACGTCTCGCCCCGCCAGATGGTGTCGGTGGCCACCGCGATGATTCCCTTCCTGGAGCACGACGACGCCAACCGTGCCCTCATGGGGGCAAACATGCAGCGCCAGGCGGTGCCGCTGGTCCGTAGCGAGGCCCCGCTGGTGGGCACCGGGATGGAGCTGCGCGCGGCGATCGACGCCGGCGACGTCGTCGTCGCCGAAGAAAGCGGCGTCATCGAGGAGGTGTCGGCCGACTACATCACTGTGATGCACGACAACGGCACCCGGCGTACCTACCGGATGCGCAAGTTT'
# primer_pool, accepted_primers = result_extraction(primer_pool, accepted_primers, seq_temp, 1, padding)

# %% trial with the whole genome at once
# df.loc[0][['pLeft_Sequences','pRight_Sequences']].values.tolist()

# %%

def extract_sequence_from_string(start_pos, end_pos, padding=150, sequence_string='', sequence_id='Chromosome'):
    """
    Extracts a subsequence from a given sequence string based on the start position and end position.
    """
    # Convert the sequence string to a Bio.Seq object
    sequence = Seq(sequence_string)
    
    # Extract the subsequence based on the start and end positions
    subsequence = sequence[start_pos - padding:end_pos + padding]
    
    return str(subsequence)  # Return the subsequence as a string

# selected_seq = extract_sequence_from_string(249, 270, padding = 0, sequence_string=sequence[:250]+'N'*100+sequence[-250:])
# print(selected_seq)
# %%
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
ref_genome= 'MTB-h37rv_asm19595v2-eg18.fa'
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

def genome_size(fasta_file):
    total_length = 0
    with open(fasta_file, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                total_length += len(line.strip())
    return total_length
#%%
def complement_sequence(dna_sequence):
    trans = str.maketrans('ATCG', 'TAGC')
    return dna_sequence.upper().translate(trans)

# print(complement_sequence("ATGCGTA"))
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

def find_sequence_location(query_seq, fasta_file=ref_genome):
    ref_genome = SeqIO.read(fasta_file, "fasta")
    position = ref_genome.seq.find(query_seq)
    if position != -1:
        return [position, position+len(query_seq)]
    else:
        print('!!!Primer not found in the reference genome')
        return 0

# pLeft_coord = []
# pLeft_coord.append(find_sequence_location('CATCGCACGTCGTCTTTCCG', ref_genome)[0])

# find_sequence_location('CATCGCACGTCGTCTTTCCG')[0]
# find_sequence_location(complement_sequence('CATCGCACGTCGTCTTTCCG'))
#%%
def simplified_tm(seq):
    # Simplified melting temperature calculation based on base pair count
    return (seq.count('A') + seq.count('T')) * 2 + (seq.count('C') + seq.count('G')) * 4

def simplified_dg(seq):
    # Simplified Gibbs free energy calculation based on the nearest-neighbor model (not accurate)
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

def calculate_similarity(seq1, seq2):
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    return matches / len(seq1) * 100

def has_multiple_binding_sites(sequence, genome, similarity_threshold=81, min_tm=60, max_dg=-10):
    seq_len = len(sequence)
    genome_len = len(genome)
    count = 0

    for i in range(genome_len - seq_len + 1):
        subseq = genome[i:i + seq_len]
        if calculate_similarity(subseq, sequence)>similarity_threshold:
            tm = simplified_tm(subseq)
            dg = simplified_dg(subseq)
            
            within_tm_threshold = tm >= min_tm
            within_dg_threshold = dg <= max_dg

            if within_tm_threshold and within_dg_threshold:
                # print(genome[i:i + seq_len])
                count += 1
                if count > 1:
                    return True  # Early exit if more than one binding site is found

    return False  # Return False if only one or no binding sites are found

# # Usage:
# primer = 'CGAACTCGAGGCTGCCTACT'
# # primer = 'GCTCGTCCATGTCCCACCAT'
# genome = primer_selection.extract_sequence_from_fasta(0, genome_size(ref_genome),0)
# result = has_multiple_binding_sites(primer, genome, 81)
# print(result)  # Output: True or False

# %%
def result_extraction(primer_pool, accepted_primers, sequence, seq, padding, ref_genome, high_b, low_b, read_size):
    # print([len(sequence)-50,len(sequence)+50])
    # print(len(sequence))
    # size_range = f'{int(len(sequence)-padding*1.3)}-{int(len(sequence)-padding*1)}'
    size_range = f'{len(sequence)-padding*2}-{len(sequence)}'
    genome = extract_sequence_from_fasta(0, genome_size(ref_genome),0)
    no_primer = '-'
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
        pLeft_ID.append(f'P{seq}-L{i}')
        # pLeft_coord.append(primer_num['COORDS'][0]+low_b)
        # print(primer_num['SEQUENCE'])
        pLeft_coord.append(find_sequence_location(primer_num['SEQUENCE'], ref_genome)[0])
        
        pLeft_length.append(primer_num['COORDS'][1])
        pLeft_Tm.append(primer_num['TM'])
        pLeft_GC.append(primer_num['GC_PERCENT'])
        pLeft_Sequences.append(complement_sequence(primer_num['SEQUENCE']))
        pLeft_EndStability.append(primer_num['END_STABILITY'])
        
    for i, primer_num in enumerate(results['PRIMER_RIGHT']):
        pRight_ID.append(f'P{seq}-R{i}')
        # pRight_coord.append(primer_num['COORDS'][0]+low_b)
        pRight_coord.append(find_sequence_location(reverse_complement_sequence(primer_num['SEQUENCE']),ref_genome)[1])
        pRight_length.append(primer_num['COORDS'][1])
        pRight_Tm.append(primer_num['TM'])
        pRight_GC.append(primer_num['GC_PERCENT'])
        pRight_Sequences.append(reverse_complement_sequence(primer_num['SEQUENCE']))
        pRight_EndStability.append(primer_num['END_STABILITY'])

    df = pd.DataFrame({'pLeft_ID':pLeft_ID, 'pLeft_coord':pLeft_coord, 'pLeft_length':pLeft_length, 'pLeft_Tm':pLeft_Tm, 'pLeft_GC':pLeft_GC, 'pLeft_Sequences':pLeft_Sequences, 'pLeft_EndStability':pLeft_EndStability, 
                    'pRight_ID':pRight_ID, 'pRight_coord':pRight_coord, 'pRight_length':pRight_length, 'pRight_Tm':pRight_Tm, 'pRight_GC':pRight_GC, 'pRight_Sequences':pRight_Sequences, 'pRight_EndStability':pRight_EndStability, 
                    'Penalty':Penalty, 'Product_size':Product_size})
    # print(df)
    # print(df[['pLeft_coord','pRight_coord','Product_size','pLeft_Sequences','pRight_Sequences']])
    # print('original_range:',low_b, high_b)
    # tm_params, dg_params = precompute_dinucleotide_params(genome)
    
    if len(primer_pool) == 0:
        print(f'{df.shape[0]} primers designed')
        for i, row in df.iterrows():
            # print(i, df.shape[0]-1)
            left_ok = True
            right_ok = True
            too_far = False
            #runs first to avoid meaningless running of has_multiple_binding_sites which take a long time
            if abs(low_b - row['pLeft_coord']) > read_size/2:
                left_ok = False
                too_far = True
            if abs(high_b - row['pRight_coord']) > read_size/2:
                right_ok = False
                too_far = True
            if (not left_ok or not right_ok) and i != df.shape[0]-1:
                print(f'Primer pair #{i} has alternative binding site')
                continue
            else:
                pass
                # print(f'Primer pair #{i} has alternative binding site')
                
            if too_far == False:
                left_ok = not has_multiple_binding_sites(row['pLeft_Sequences'], genome)
                right_ok = not has_multiple_binding_sites(reverse_complement_sequence(row['pRight_Sequences']), genome)
                if (not left_ok or not right_ok) and i != df.shape[0]-1:
                    print(f'Primer pair #{i} has alternative binding site')
                    continue  
                else:
                    pass
                    # print(f'Primer pair #{i} has alternative binding site')
            else:
                pass
            
            # print(too_far, left_ok, right_ok)

            #if pass all filtering
            if left_ok == True and right_ok == True:
                primer_pool.append(complement_sequence(row['pLeft_Sequences']))
                primer_pool.append(row['pRight_Sequences'])
                row_df = pd.DataFrame(row).T
                accepted_primers = pd.concat([accepted_primers, row_df],axis=0)
                print(f'Primer pair #{i} accepted')
                print(low_b, row['pLeft_coord'])
                print(abs(low_b - row['pLeft_coord']) > read_size/2)
                no_primer = '-'
                print('***')
                break
            else:
                if i == df.shape[0]-1:
                    print(f'!!!No suitable primer found:, please manually inspect the sequence')
                    if abs(low_b - row['pLeft_coord']) > read_size/2:
                        left_ok = False
                        print('!Problem with left(forward) primer')
                        print('How should I moved the left range? (e.g. -50 = moving start of covered range 50bp upstream)')
                        change = input('Where to move (+/-bps):')
                        low_b = low_b+int(change)
                        seq_template = extract_sequence_from_fasta(low_b, high_b, padding=0)
                        print(f'Redesigning primers for the new range ({change}bps): {low_b, high_b} = {low_b-int(change), high_b}')
                        primer_pool, accepted_primers, no_primer = result_extraction(primer_pool, accepted_primers, seq_template, i+1, padding, ref_genome = ref_genome, high_b = high_b, low_b = low_b, read_size = read_size)
                    elif abs(high_b - row['pRight_coord']) > read_size/2:
                        right_ok = False
                        print('!Problem with Right(backward) primer')
                        print('How should I moved the right range? (e.g. -50 = moving start of covered range 50bp upstream)')
                        change = input('Where to move (+/-bps):')
                        high_b = high_b+int(change)
                        seq_template = extract_sequence_from_fasta(low_b, high_b+int(change), padding=0)
                        print(f'Redesigning primers for the new range ({change}bps): {low_b, high_b} instead of {low_b, high_b-int(change)}')
                        primer_pool, accepted_primers, no_primer = result_extraction(primer_pool, accepted_primers, seq_template, i+1, padding, ref_genome = ref_genome, high_b = high_b, low_b = low_b, read_size = read_size)
                    no_primer='Redesigned'
                else:
                    print(f'Primer pair #{i} has alternative binding site')
                    continue
        # print(primer_pool, accepted_primers)
        # print(accepted_primers)
        # print(df.iloc[0])
        # print(1)
    else:
        #print('Checking for homodimer')
        print(f'{df.shape[0]} primers designed')
        for i, row in df.iterrows():
            # print(row)
            left_ok = True
            right_ok = True
            too_far = False
            hetero = False
            # checking if the primer is too far away from the original range
            if abs(low_b - row['pLeft_coord']) > read_size/2:
                left_ok = False
                too_far = True
            if abs(high_b - row['pRight_coord']) > read_size/2:
                right_ok = False
                too_far = True
            if (not left_ok or not right_ok) and i != df.shape[0]-1:
                print(f'Primer pair #{i} has alternative binding site')
                continue
            else:
                pass
                # print(f'Primer pair #{i} has alternative binding site')

            if too_far == False:
                for x in primer_pool: # heterodimer check
                    if check_heterodimer(x, complement_sequence(row['pLeft_Sequences'])) == False:
                        left_ok = False
                        hetero = True
                    if check_heterodimer(x, row['pRight_Sequences']) == False:
                        right_ok = False
                        hetero = True

                if (not left_ok or not right_ok) and i != df.shape[0]-1:
                    print(f'Primer pair #{i} has homodimer')
                    continue  
                else:
                    pass
                    # print(f'Primer pair #{i} has alternative binding site')
            else:
                pass
            
            if too_far or hetero:
                pass
            else:  
                left_ok = not has_multiple_binding_sites(row['pLeft_Sequences'], genome)
                right_ok = not has_multiple_binding_sites(reverse_complement_sequence(row['pRight_Sequences']), genome)
                if (not left_ok or not right_ok) and i != df.shape[0]-1:
                    print(f'Primer pair #{i} has alternative binding site')
                    continue
                else:
                    pass
                    # print(f'Primer pair #{i} has alternative binding site')
            # print(too_far, hetero, left_ok, right_ok)
            if left_ok == True and right_ok == True:
                row_df = pd.DataFrame(row).T
                # primer_pool.extend(row[['pLeft_Sequences','pRight_Sequences']].values.tolist())
                primer_pool.append(complement_sequence(row['pLeft_Sequences']))
                primer_pool.append(row['pRight_Sequences'])
                # print(row)
                accepted_primers = pd.concat([accepted_primers, row_df],axis=0)
                print(f'Primer pair #{i} accepted')
                # print(low_b, row['pLeft_coord'])
                # print(abs(low_b - row['pLeft_coord']) > read_size/2)
                no_primer='-'
                print('***')
                break
            else:
                # print(i, df.shape[0]-1)
                if i == df.shape[0]-1:
                    print(f'!!!No suitable primer found: please manually inspect the sequence')
                    # no_primer.append(seq)
                    if abs(low_b - row['pLeft_coord']) > read_size/2:
                        left_ok = False
                        print('!Problem with left(forward) primer')
                        print('How should I moved the left range? (e.g. -50 = moving start of covered range 50bp upstream)')
                        change = input('Where to move (+/-bps):')
                        low_b = low_b+int(change)
                        seq_template = extract_sequence_from_fasta(low_b, high_b, padding=0)
                        print(f'Redesigning primers for the new range ({change}bps): {low_b, high_b} instead of {low_b-int(change), high_b}')
                        primer_pool, accepted_primers, no_primer = result_extraction(primer_pool, accepted_primers, seq_template, i+1, padding, ref_genome = ref_genome, high_b = high_b, low_b = low_b, read_size = read_size)

                    if abs(high_b - row['pRight_coord']) > read_size/2:
                        right_ok = False
                        print('!Problem with Right(backward) primer')
                        print('How should I moved the right range? (e.g. -50 = moving start of covered range 50bp stream)')
                        change = input('Where to move (+/-bps):')
                        high_b = high_b+int(change)
                        seq_template = extract_sequence_from_fasta(low_b, high_b, padding=0)
                        print(f'Redesigning primers for the new range ({change}bps): {low_b, high_b} instead of {low_b, high_b-int(change)}')
                        primer_pool, accepted_primers, no_primer = result_extraction(primer_pool, accepted_primers, seq_template, i+1, padding, ref_genome = ref_genome, high_b = high_b, low_b = low_b, read_size = read_size)
                    no_primer='Redesigned'
                else:
                    print(f'Primer pair #{i} has alternative binding site')
                    continue
            #Alternative binding check
    return primer_pool, accepted_primers, no_primer

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

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


# %%
reference_genome= '/mnt/storage10/lwang/Projects/Amplicone_design_tool/MTB-h37rv_asm19595v2-eg18.fa'
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

def extract_sequence_from_fasta(start_pos, end_pos, padding = 150, fasta_file= '/mnt/storage10/lwang/Projects/Amplicone_design_tool/MTB-h37rv_asm19595v2-eg18.fa', sequence_id='Chromosome'):
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

def complement_sequence(seq):
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(complement[base] for base in seq)

def reverse_complement_sequence(seq):
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    reverse_seq = seq[::-1]
    return "".join(complement[base] for base in reverse_seq)

def updown_stream_primer_range(start_pos, end_pos, dis_range=150):
    up_stream = complement_sequence(extract_sequence_from_fasta(start_pos-dis_range, start_pos))
    down_stream = reverse_complement_sequence(extract_sequence_from_fasta(end_pos, end_pos+dis_range))
    return up_stream, down_stream

# %%
# updown_stream_primer_range(1000, 2000)

# %%
# from Bio.Blast.Applications import NcbiblastnCommandline

# def blast_primers(primer_sequences, reference_genome_file):
#     # Create a temporary FASTA file for the primer sequences
#     primer_file = "primers.fasta"
#     with open(primer_file, "w") as file:
#         for i, seq in enumerate(primer_sequences):
#             file.write(f">Primer_{i}\n{seq}\n")

#     # Perform BLAST search using primer sequences against the reference genome
#     blastn_cline = NcbiblastnCommandline(query=primer_file, subject=reference_genome_file, out="m_cold.xml", outfmt=5)#, out="blast_results.txt")
#     stdout, stderr = blastn_cline()

#     # Process the BLAST results as desired
#     # Here, we simply print the results
#     with open("blast_results.txt") as file:
#         blast_results = file.read()
#         print(blast_results)

#     # Remove the temporary primer file and BLAST result file
#     # os.remove(primer_file)
#     # os.remove("blast_results.txt")
    
# #test
# blast_primers(['AGAAGGTGTGGAAGTTGTGGAACGTGTTACGGTTGTTCGTTTAGC', 'GTCGAGGGGTTTGCTCTGTT','CGATCGATCGATCGATCGAT'], reference_genome)
#%%
# primer_file = "primers.fasta"
# reference_genome_file= '/mnt/storage10/lwang/Projects/Amplicone_design_tool/MTB-h37rv_asm19595v2-eg18.fa'

# blastn_cline = NcbiblastnCommandline(query=primer_file, subject=reference_genome_file, outfmt=6, out="blast_results.txt")
# # query="m_cold.fasta", db="nt", strand="plus",evalue=0.001, out="m_cold.xml", outfmt=5
# stdout, stderr = blastn_cline()
# %%
def check_heterodimer(primer1, primer2):
    # Calculate melting temperature (Tm) for each primer
    tm1 = MeltingTemp.Tm_NN(primer1)
    tm2 = MeltingTemp.Tm_NN(primer2)

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


from Bio import SeqIO

def find_sequence_location(query_sequence, position, fasta_file=reference_genome):
    closest_start = float('inf')
    closest_end = float('inf')
    closest_record = None

    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        index = sequence.find(query_sequence)
        
        while index != -1:
            start_index = index
            end_index = index + len(query_sequence)
            
            if abs(position - start_index) < abs(position - closest_start):
                closest_start = start_index
                closest_end = end_index
                closest_record = record.id

            index = sequence.find(query_sequence, index + 1)

    if closest_record is not None:
        return closest_start, closest_end
    else:
        print("Primer Sequence not found in FASTA file.")
        return 0,0

# %%
def result_extraction(primer_pool, accepted_primers, sequence, seq, padding, low_b, high_b):
    # print([len(sequence)-50,len(sequence)+50])
    # print(len(sequence))
    size_range = f'{int(len(sequence)-padding*1.3)}-{int(len(sequence)-padding*1)}'

    # size_range = f'{len(sequence)-350}-{len(sequence)-250}'
    print('size_range:',size_range)
    print('SEQUENCE_INCLUDED_REGION:', [padding-10,len(sequence)-padding+10],)
    try:
        results = bindings.design_primers(
            seq_args={
                'SEQUENCE_ID': 'MH1000',
                'SEQUENCE_TEMPLATE': sequence,
                # 'SEQUENCE_INCLUDED_REGION': [150,len(sequence)-150]
                # 'SEQUENCE_INCLUDED_REGION': [(0,len(sequence)),],
                'SEQUENCE_INCLUDED_REGION': [(0,padding),(len(sequence)-padding,len(sequence))],
                'SEQUENCE_EXCLUDED_REGION':[(padding,len(sequence)-padding)]
            },
            global_args={
                'PRIMER_NUM_RETURN': 15,
                'PRIMER_OPT_SIZE': 20,
                'PRIMER_MIN_3_PRIME_DISTANCE': 20,
                'PRIMER_MIN_5_PRIME_DISTANCE': 20,
                'PRIMER_PICK_INTERNAL_OLIGO': 0,
                'PRIMER_INTERNAL_MAX_SELF_END': 8,
                'PRIMER_MIN_SIZE': 15,
                'PRIMER_MAX_SIZE': 35,
                'PRIMER_OPT_TM': 62.0,
                'PRIMER_MIN_TM': 55.0,
                'PRIMER_MAX_TM': 68.0,
                'PRIMER_MIN_GC': 42.0,
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
        pLeft_ID.append(f'Seq{seq}-Left{i}')
        # pLeft_coord.append(primer_num['COORDS'][0]+low_b)
        pLeft_coord.append(find_sequence_location(primer_num['SEQUENCE'], low_b,reference_genome)[0])

        pLeft_length.append(primer_num['COORDS'][1])
        pLeft_Tm.append(primer_num['TM'])
        pLeft_GC.append(primer_num['GC_PERCENT'])
        pLeft_Sequences.append(primer_num['SEQUENCE'])
        pLeft_EndStability.append(primer_num['END_STABILITY'])
        
    for i, primer_num in enumerate(results['PRIMER_RIGHT']):
        pRight_ID.append(f'Seq{seq}-Right{i}')
        # pRight_coord.append(primer_num['COORDS'][0]+low_b)
        pRight_coord.append(find_sequence_location(reverse_complement_sequence(primer_num['SEQUENCE']), high_b,reference_genome)[1])
        pRight_length.append(primer_num['COORDS'][1])
        pRight_Tm.append(primer_num['TM'])
        pRight_GC.append(primer_num['GC_PERCENT'])
        pRight_Sequences.append(primer_num['SEQUENCE'])
        pRight_EndStability.append(primer_num['END_STABILITY'])

    df = pd.DataFrame({'pLeft_ID':pLeft_ID, 'pLeft_coord':pLeft_coord, 'pLeft_length':pLeft_length, 'pLeft_Tm':pLeft_Tm, 'pLeft_GC':pLeft_GC, 'pLeft_Sequences':pLeft_Sequences, 'pLeft_EndStability':pLeft_EndStability, 
                    'pRight_ID':pRight_ID, 'pRight_coord':pRight_coord, 'pRight_length':pRight_length, 'pRight_Tm':pRight_Tm, 'pRight_GC':pRight_GC, 'pRight_Sequences':pRight_Sequences, 'pRight_EndStability':pRight_EndStability, 
                    'Penalty':Penalty, 'Product_size':Product_size})
    # print(df)
    print(df[['pLeft_coord','pRight_coord','Product_size','pLeft_Sequences','pRight_Sequences']])
    print('original_range:',low_b, high_b)
    
    if len(primer_pool) == 0:
        primer_pool.extend(df.loc[0][['pLeft_Sequences','pRight_Sequences']].values.tolist())
        first_row_df = pd.DataFrame(df.iloc[0]).T
        accepted_primers = pd.concat([accepted_primers, first_row_df],axis=0)
        # print(accepted_primers)
        # print(df.iloc[0])
    else:
        for i, row in df.iterrows():
            # print(row)
            left_ok = True
            right_ok = True
            for x in primer_pool:
                if check_heterodimer(x, row['pLeft_Sequences']) == False:
                    left_ok = False 
                if check_heterodimer(x, row['pRight_Sequences']) == False:
                    right_ok = False
            if left_ok == True and right_ok == True:
                row_df = pd.DataFrame(row).T
                primer_pool.extend(row[['pLeft_Sequences','pRight_Sequences']].values.tolist())
                # print(row)
                accepted_primers = pd.concat([accepted_primers, row_df],axis=0)
            break
    return primer_pool, accepted_primers

def result_extraction_sides(primer_pool, accepted_primers, sequence, seq, padding, low_b, high_b):
    # print([len(sequence)-50,len(sequence)+50])
    # print(len(sequence))
    primer_pool = primer_pool[:-2]
    accepted_primers = accepted_primers[:-1]
    left_sequence = sequence[:padding]
    right_sequence = sequence[-padding:]
    size_range = f'0-{padding}'

    # size_range = f'{len(sequence)-350}-{len(sequence)-250}'
    print('size_range:',size_range)
    print('SEQUENCE_INCLUDED_REGION:', [0,len(left_sequence)])
    try:
        left_results = bindings.design_primers(
            seq_args={
                'SEQUENCE_ID': 'leftprimer',
                'SEQUENCE_TEMPLATE': left_sequence,
                # 'SEQUENCE_INCLUDED_REGION': [150,len(sequence)-150]
                # 'SEQUENCE_INCLUDED_REGION': [0,len(left_sequence)],
            },
            global_args={
                'PRIMER_NUM_RETURN': 20,
                'PRIMER_PICK_RIGHT_PRIMER':0,
                'PRIMER_OPT_SIZE': 20,
                'PRIMER_PICK_INTERNAL_OLIGO': 0,
                'PRIMER_INTERNAL_MAX_SELF_END': 8,
                'PRIMER_MIN_THREE_PRIME_DISTANCE':10,
                'PRIMER_MIN_FIVE_PRIME_DISTANCE':10,
                'PRIMER_MIN_SIZE': 18,
                'PRIMER_MAX_SIZE': 30,
                'PRIMER_OPT_TM': 62.0,
                'PRIMER_MIN_TM': 60.0,
                'PRIMER_MAX_TM': 64.0,
                'PRIMER_MIN_GC': 45.0,
                'PRIMER_MAX_GC': 60.0,
                'PRIMER_MAX_POLY_X': 5,
                'PRIMER_INTERNAL_MAX_POLY_X': 5,
                'PRIMER_SALT_MONOVALENT': 50.0,
                'PRIMER_DNA_CONC': 50.0,
                'PRIMER_MAX_NS_ACCEPTED': 0,
                'PRIMER_MAX_SELF_ANY': 5,
                'PRIMER_MAX_SELF_END': 2,
                'PRIMER_PAIR_MAX_COMPL_ANY': 5,
                'PRIMER_PAIR_MAX_COMPL_END': 2,
                # 'PRIMER_PRODUCT_SIZE_RANGE': size_range,
                # 'PRIMER_PRODUCT_SIZE_RANGE': '950-1050',
                # 'PRIMER_PRODUCT_SIZE_RANGE': [
                #     # [950,1050]
                #     [len(sequence)-350,len(sequence)-250]
                # ],
            })
    except:
        print('!!!Left Primer extraction error')
        
    try:
        right_results = bindings.design_primers(
        seq_args={
            'SEQUENCE_ID': 'rightprimer',
            'SEQUENCE_TEMPLATE': right_sequence,
            # 'SEQUENCE_INCLUDED_REGION': [150,len(sequence)-150]
            # 'SEQUENCE_INCLUDED_REGION': [0, len(right_sequence)],
        },
        global_args={
            'PRIMER_NUM_RETURN': 20,
            'PRIMER_PICK_LEFT_PRIMER':0,            
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_PICK_INTERNAL_OLIGO': 0,
            'PRIMER_INTERNAL_MAX_SELF_END': 8,
            'PRIMER_MIN_THREE_PRIME_DISTANCE':10,
            'PRIMER_MIN_FIVE_PRIME_DISTANCE':10,
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_MAX_SIZE': 30,
            'PRIMER_OPT_TM': 62.0,
            'PRIMER_MIN_TM': 60.0,
            'PRIMER_MAX_TM': 64.0,
            'PRIMER_MIN_GC': 45.0,
            'PRIMER_MAX_GC': 60.0,
            'PRIMER_MAX_POLY_X': 5,
            'PRIMER_INTERNAL_MAX_POLY_X': 5,
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_DNA_CONC': 50.0,
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_MAX_SELF_ANY': 5,
            'PRIMER_MAX_SELF_END': 2,
            'PRIMER_PAIR_MAX_COMPL_ANY': 5,
            'PRIMER_PAIR_MAX_COMPL_END': 2,
            # 'PRIMER_PRODUCT_SIZE_RANGE': size_range,
            # 'PRIMER_PRODUCT_SIZE_RANGE': '950-1050',
            # 'PRIMER_PRODUCT_SIZE_RANGE': [
            #     # [950,1050]
            #     [len(sequence)-350,len(sequence)-250]
            # ],
        })
    except:
        print('!!!Right Primer extraction error')
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

    # left_Penalty = []
    Product_size = []
    # right_Penalty = []
    # for i, primer_num in enumerate(left_results['PRIMER_PAIR']):
    #     Product_size.append(primer_num['PRODUCT_SIZE'])
    #     Penalty.append(primer_num['PENALTY'])
    for i, primer_num in enumerate(left_results['PRIMER_LEFT']):
        pLeft_ID.append(f'Seq{seq}-Left{i}')
        # pLeft_coord.append(primer_num['COORDS'][0]+low_b)
        pLeft_coord.append(find_sequence_location(primer_num['SEQUENCE'], low_b,reference_genome)[0])

        pLeft_length.append(primer_num['COORDS'][1])
        pLeft_Tm.append(primer_num['TM'])
        pLeft_GC.append(primer_num['GC_PERCENT'])
        pLeft_Sequences.append(primer_num['SEQUENCE'])
        pLeft_EndStability.append(primer_num['END_STABILITY'])
        
    for i, primer_num in enumerate(right_results['PRIMER_RIGHT']):
        pRight_ID.append(f'Seq{seq}-Right{i}')
        # pRight_coord.append(primer_num['COORDS'][0]+low_b)
        pRight_coord.append(find_sequence_location(reverse_complement_sequence(primer_num['SEQUENCE']), high_b,reference_genome)[1])
        pRight_length.append(primer_num['COORDS'][1])
        pRight_Tm.append(primer_num['TM'])
        pRight_GC.append(primer_num['GC_PERCENT'])
        pRight_Sequences.append(primer_num['SEQUENCE'])
        pRight_EndStability.append(primer_num['END_STABILITY'])

    df_left = pd.DataFrame({'pLeft_ID':pLeft_ID, 'pLeft_coord':pLeft_coord, 'pLeft_length':pLeft_length, 'pLeft_Tm':pLeft_Tm, 'pLeft_GC':pLeft_GC, 'pLeft_Sequences':pLeft_Sequences, 'pLeft_EndStability':pLeft_EndStability, })
    df_right = pd.DataFrame({'pRight_ID':pRight_ID, 'pRight_coord':pRight_coord, 'pRight_length':pRight_length, 'pRight_Tm':pRight_Tm, 'pRight_GC':pRight_GC, 'pRight_Sequences':pRight_Sequences, 'pRight_EndStability':pRight_EndStability, })
    print(df_left.shape)
    print(df_right.shape)
    # print(df)
    # print(df[['pLeft_coord','pRight_coord','Product_size','pLeft_Sequences','pRight_Sequences']])
    # print('original_range:',low_b, high_b)
    
    if len(primer_pool) == 0:
        found = False
        for x in range(df_left.shape[0]):
            if found:
                break
            for y in range(df_right.shape[0]):
                if check_heterodimer(df_left.loc[x]['pLeft_Sequences'], df_right.loc[y]['pRight_Sequences']) == False:
                    primer_pool.extend(df_left.loc[x]['pLeft_Sequences'].values.tolist())
                    primer_pool.extend(df_left.loc[x]['pLeft_Sequences'].values.tolist())
                    data = {'Penalty': ['re-assembled'], 'Product_size': [df_right.loc[x]['pRight_coord'] - df_left.loc[0]['pLeft_coord']]}
                    df_insert = pd.DataFrame(data)
                    first_row_df = pd.concat([df_left.loc[x], df_right.loc[x], df_insert], axis=1)
                    accepted_primers = pd.concat([accepted_primers, first_row_df],axis=0)
                    print('found1',x,y)
                    break
                else:
                    print('homodimer1',x,y)
        # first_row_df = pd.DataFrame(df_left.iloc[0]).T
        # print(accepted_primers)
        # print(df.iloc[0])
    else:
        for i in range(df_left.shape[0]):
            # print(row)
            left_ok = True
            right_ok = True
            for x in primer_pool:
                if check_heterodimer(x, df_left.loc[i]['pLeft_Sequences']) != False:
                    left_ok = False 
                if check_heterodimer(x, df_right.loc[i]['pRight_Sequences']) != False:
                    right_ok = False
            if left_ok == True and right_ok == True:
                if check_heterodimer([df_left.loc[x]['pLeft_Sequences'], df_right.loc[x]['pRight_Sequences']]) == False:
                    primer_pool.extend(df_left.loc[x]['pLeft_Sequences'].values.tolist())
                    primer_pool.extend(df_left.loc[x]['pLeft_Sequences'].values.tolist())
                    data = {'Penalty': ['re-assembled'], 'Product_size': [df_right.loc[x]['pRight_coord'] - df_left.loc[0]['pLeft_coord']]}
                    df_insert = pd.DataFrame(data)
                    first_row_df = pd.concat([df_left.loc[x], df_right.loc[x], df_insert], axis=1)
                    accepted_primers = pd.concat([accepted_primers, first_row_df],axis=0)
                    print('found2',x)
                    break

                else:
                    print('homodimer2',x)
                    continue         
    return primer_pool, accepted_primers
# primer_pool, accepted_primers = result_extraction(primer_pool, accepted_primers, test, seq)

# %% trial with the whole genome at once
df.loc[0][['pLeft_Sequences','pRight_Sequences']].values.tolist()

#%%
tb_genome = extract_sequence_from_fasta(0, 4411532, padding = 0, fasta_file= '/mnt/storage10/lwang/Projects/Amplicone_design_tool/MTB-h37rv_asm19595v2-eg18.fa', sequence_id='Chromosome')

# covered_ranges = [(4326008, 4326858), (4043869, 4044719), (4247642, 4248492), (1472644, 1473494), (2284366, 2285216), (2287239, 2288089), (4407543, 4408393), (1917939, 1918789), (760314, 761164), (2518310, 2519160), (2154492, 2155342), (2715339, 2716189), (3071476, 3072326), (4246734, 4247584), (765572, 766422), (4327029, 4327879), (1475956, 1476806), (1416212, 1417062), (764363, 765213), (762089, 762939), (3087136, 3087986), (4325194, 4326044), (4240592, 4241442), (3072537, 3073387), (2282776, 2283626), (3073768, 3074618), (4243203, 4244053), (2285354, 2286204), (6575, 7425), (2288412, 2289262)]
# covered_ranges = [(4326008, 4327008)]
covered_ranges = [(2226008, 2226858)]
formatted_ranges = [f"{start}-{end}" for start, end in covered_ranges]

results = bindings.design_primers(
            seq_args={
                'SEQUENCE_ID': 'MH1000',
                'SEQUENCE_TEMPLATE':tb_genome,
                # 'SEQUENCE_INCLUDED_REGION': [150,len(sequence)-150]
                # 'SEQUENCE_INCLUDED_REGION': [(0,len(sequence)),],
                'SEQUENCE_INCLUDED_REGION': covered_ranges,
                # 'SEQUENCE_EXCLUDED_REGION':[(padding,len(sequence)-padding)]
            },
            global_args={
                'PRIMER_NUM_RETURN': 5,
                'PRIMER_OPT_SIZE': 20,
                'PRIMER_MIN_3_PRIME_DISTANCE': 20,
                'PRIMER_MIN_5_PRIME_DISTANCE': 20,
                'PRIMER_PICK_INTERNAL_OLIGO': 0,
                'PRIMER_INTERNAL_MAX_SELF_END': 8,
                'PRIMER_MIN_SIZE': 15,
                'PRIMER_MAX_SIZE': 35,
                'PRIMER_OPT_TM': 62.0,
                'PRIMER_MIN_TM': 55.0,
                'PRIMER_MAX_TM': 68.0,
                'PRIMER_MIN_GC': 42.0,
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
                # 'PRIMER_PRODUCT_SIZE_RANGE': [500,1100],
                'PRIMER_PRODUCT_SIZE_RANGE': '850-1250',
                # 'PRIMER_PRODUCT_SIZE_RANGE': [
                #     # [950,1050]
                #     [len(sequence)-350,len(sequence)-250]
                # ],
            })
#%%
len(extract_sequence_from_fasta(0, 4411532, padding = 0, fasta_file= '/mnt/storage10/lwang/Projects/Amplicone_design_tool/MTB-h37rv_asm19595v2-eg18.fa', sequence_id='Chromosome'))
# %%
sequence = extract_sequence_from_fasta(48586, 48606, padding = 0, fasta_file= '/mnt/storage10/lwang/Projects/Amplicone_design_tool/MTB-h37rv_asm19595v2-eg18.fa', sequence_id='Chromosome')
# sequence[:250]+'N'*100+sequence[-250:]
print(sequence)
#%%
left_results = bindings.design_primers(
            seq_args={
                'SEQUENCE_ID': 'leftprimer',
                # 'SEQUENCE_TEMPLATE': sequence[:250]+sequence[-250:],
                'SEQUENCE_TEMPLATE': sequence[:250]+sequence[-250:],
                'SEQUENCE_INCLUDED_REGION': [100,400]
                # 'SEQUENCE_INCLUDED_REGION': sequence[:251],
            },
            global_args={
                'PRIMER_NUM_RETURN': 15,
                'PRIMER_OPT_SIZE': 20,
                'PRIMER_NUM_RETURN': 15,
                'PRIMER_PICK_INTERNAL_OLIGO': 0,
                'PRIMER_INTERNAL_MAX_SELF_END': 8,
                'PRIMER_MIN_THREE_PRIME_DISTANCE':10,
                'PRIMER_MIN_FIVE_PRIME_DISTANCE':10,
                'PRIMER_MIN_SIZE': 18,
                'PRIMER_MAX_SIZE': 30,
                'PRIMER_OPT_TM': 62.0,
                'PRIMER_MIN_TM': 60.0,
                'PRIMER_MAX_TM': 64.0,
                'PRIMER_MIN_GC': 45.0,
                'PRIMER_MAX_GC': 60.0,
                'PRIMER_MAX_POLY_X': 5,
                'PRIMER_INTERNAL_MAX_POLY_X': 5,
                'PRIMER_SALT_MONOVALENT': 50.0,
                'PRIMER_DNA_CONC': 50.0,
                'PRIMER_MAX_NS_ACCEPTED': 0,
                'PRIMER_MAX_SELF_ANY': 5,
                'PRIMER_MAX_SELF_END': 2,
                'PRIMER_PAIR_MAX_COMPL_ANY': 5,
                'PRIMER_PAIR_MAX_COMPL_END': 2,
                # 'PRIMER_PRODUCT_SIZE_RANGE': [490,510],
                'PRIMER_PRODUCT_SIZE_RANGE': '300-400',
                # 'PRIMER_PRODUCT_SIZE_RANGE': [
                #     # [950,1050]
                #     [len(sequence)-350,len(sequence)-250]
                # ],
            })
# %%
low_b = 4326008
seq = 1
pLeft_ID = []
pLeft_coord = []
pLeft_length = []
pLeft_Tm = []
pLeft_GC = []
pLeft_Sequences = []
pLeft_EndStability = []
Penalty = []
Product_size = []
for i, primer_num in enumerate(left_results['PRIMER_LEFT']):
    pLeft_ID.append(f'Seq{seq}-Left{i}')
    # pLeft_coord.append(primer_num['COORDS'][0]+low_b)
    pLeft_coord.append(find_sequence_location(primer_num['SEQUENCE'], low_b,reference_genome)[0])

    pLeft_length.append(primer_num['COORDS'][1])
    pLeft_Tm.append(primer_num['TM'])
    pLeft_GC.append(primer_num['GC_PERCENT'])
    pLeft_Sequences.append(primer_num['SEQUENCE'])
    pLeft_EndStability.append(primer_num['END_STABILITY'])
    
# for i, primer_num in enumerate(left_results['PRIMER_PAIR']):
#     Product_size.append(primer_num['PRODUCT_SIZE'])
#     Penalty.append(primer_num['PENALTY'])

df_left = pd.DataFrame({'pLeft_ID':pLeft_ID, 'pLeft_coord':pLeft_coord, 'pLeft_length':pLeft_length, 'pLeft_Tm':pLeft_Tm, 'pLeft_GC':pLeft_GC, 'pLeft_Sequences':pLeft_Sequences, 'pLeft_EndStability':pLeft_EndStability})


# %%
data = {'Penalty': ['re-assembled'], 'Product_size': ['Product_size']}  # Replace value1, value2, value3 with your actual values
df_insert = pd.DataFrame(data)
df_left_1 = pd.concat([df_left, df_left, df_insert], axis=1)

# %%
seq = 1
primer_pool = []
accepted_primers = []
sequence = extract_sequence_from_fasta(4326008, 4326858, padding = 250, fasta_file= '/mnt/storage10/lwang/Projects/Amplicone_design_tool/MTB-h37rv_asm19595v2-eg18.fa', sequence_id='Chromosome')
primer_pool1, accepted_primers1 = result_extraction_sides(primer_pool, accepted_primers, sequence, seq, 250, low_b=4326008, high_b=4326858)

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

selected_seq = extract_sequence_from_string(249, 270, padding = 0, sequence_string=sequence[:250]+'N'*100+sequence[-250:])
print(selected_seq)
# %%
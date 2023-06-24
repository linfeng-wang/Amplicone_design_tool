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
#%%
def print_full(x):
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', 2000)
    pd.set_option('display.float_format', '{:20,.2f}'.format)
    pd.set_option('display.max_colwidth', None)
    print(x)
    pd.reset_option('display.max_rows')
    pd.reset_option('display.max_columns')
    pd.reset_option('display.width')
    pd.reset_option('display.float_format')
    pd.reset_option('display.max_colwidth')

#%%
# raw_data = pd.read_csv("https://raw.githubusercontent.com/GaryNapier/tb-lineages/main/fst_results_clean_fst_1_for_paper.csv")
# raw_data['weight'] = sample(range(1, 10000), raw_data.shape[0])
# raw_data['weight'] = raw_data['weight'] / 100
# raw_data = raw_data.sort_values(by=['Pos'])
# raw_data = raw_data.reset_index(drop=True)
full_data = pd.read_csv('/mnt/storage10/jody/projects/variant_dump/variants.csv')
full_data = full_data[~full_data['drugs'].isna()]
full_data = full_data.sort_values(by=['genome_pos'])
full_data = full_data.reset_index(drop=True)
full_data['weight'] = full_data['freq']
#%%
# snp_weight = full_data['sample_id'].value_counts()/full_data['sample_id'].value_counts().max()
# snp_weight = snp_weight.to_frame()
#%%
# full_data['weight'] = 0
# full_data = full_data.merge(snp_weight, left_on='sample_id', right_index=True, how='left')
# merged_data = merged_data.sort_values(by=['genome_pos'])
# merged_data = merged_data.reset_index(drop=True)

# merged_data['weight'].fillna(value=0, inplace=True)

#%%
def rolling_sum(df, weight, window_size, genomic_pos):
    """
    Calculates the rolling sum of a list with a given window size.

    Parameters:
        lst (list): The list to calculate the rolling sum for.
        window_size (int): The size of the rolling window.

    Returns:
        list: A list containing the rolling sum values.
    """
    # Calculate the rolling sum using a list comprehension
    rolling_sum = []
    pos = np.unique(genomic_pos).tolist()
    # print(pos)
    for x in pos:
        start = x
        in_range = [i for i in pos if i <= start+window_size]
        end = min(in_range, key=lambda x:abs(x-(start+window_size)))
        freq_sum = df[(df['genome_pos']>=start) & (df['genome_pos']<=end)][f'{weight}'].sum()
        rolling_sum.append(freq_sum)
    return rolling_sum

# test = rolling_sum(full_data, 400, full_data['genome_pos'].tolist())

#%% trial
# window_size = 400
# pos = full_data['genome_pos'].unique()
# rolling_sum = []
# for x in pos:
#     start = x
#     in_range = [i for i in pos if i <= start+window_size]
#     end = min(in_range, key=lambda x:abs(x-(start+window_size)))
#     freq_sum = full_data[(full_data['genome_pos']>=start) & (full_data['genome_pos']<=end)]['freq'].sum()
#     rolling_sum.append(freq_sum)
full_data['weight'].min()
#%%
#Full data trial
read_number = 30
read_size = 1000
window_size = read_size
# priorities = []

weight_window_sum = rolling_sum(full_data, 'weight', window_size, full_data['genome_pos'].tolist())
pos = full_data['genome_pos'].unique()
covered_positions = {}
covered_ranges = []
for run in tqdm(range(0,read_number)):
    start = pos[np.argmax(weight_window_sum)] # find the index of the max value in the rolling sum
    in_range = [i for i in pos if i <= start+window_size]
    end = min(in_range, key=lambda x:abs(x-(start+window_size)))

    # if len(covered_ranges) != 0:
    #     pass
    # break
    # elif len(covered_ranges) > 0:
    #     for i in range(len(covered_ranges)):
    #         if start > covered_ranges[i][0] and start < covered_ranges[i][1]:
    #             start_index = pos.index(covered_ranges[i][1])+1
    #             start = pos[start_index]
    #             if covered_ranges[i][1]+1 + read_size > full_data.shape[0]:
    #                 end = full_data.shape[0]
    #             end = covered_ranges[i][1]+1 + read_size
    #         if end > covered_ranges[i][0] and end < covered_ranges[i][1]:
    #             end = covered_ranges[i][0]-1
    #             if covered_ranges[i][0]-1 - read_size < 0:
    #                 start = 0
    #             start = covered_ranges[i][0]-1 - read_size
    # else:
    #     print('error')
    # covered_positions[f'Amplicon_{run+1}'] = {'Range':{'Start': start, 'End': end}, 'Markers':full_data[['genome_pos','gene','sublin','drtype','drugs','weight']][start:end].sort_values(by=['weight']).to_dict('records')} # verbose version of output
    covered_positions[f'Amplicon_{run+1}'] = {'Range':{'Start': start, 'End': end}}  # concise version of output
    covered_ranges.append([start, end])
    # 
    # full_data.loc[(full_data['genome_pos']>=start) & (full_data['genome_pos']<=end), 'weight'] = 0 # set the weight of the covered positions to 0
    full_data.loc[(full_data['genome_pos']>=start) & (full_data['genome_pos']<=end), 'weight'] = full_data.loc[(full_data['genome_pos']>=start) & (full_data['genome_pos']<=end), 'weight']/10 # set the weight of the covered positions to 0
    weight_window_sum = rolling_sum(full_data, 'weight', window_size, full_data['genome_pos'].tolist())

#%%
# output covered position into a bed file
bed_file_path = "intervals1kbps.bed"

with open(bed_file_path, "w") as bed_file:
    for interval in covered_ranges:
        start = interval[0]
        end = interval[1]
        bed_line = f"chr1\t{start}\t{end}\n"  # Modify "chr1" with the appropriate chromosome name
        bed_file.write(bed_line)
        
#%% checking
full_data[(full_data['genome_pos']>=761004) & (full_data['genome_pos']<761284)]
full_data[full_data['genome_pos']<window_size]
full_data[(full_data['genome_pos']>=183247) & (full_data['genome_pos']<=761284)]








#%%
#Full data trial with coverage percentage
percentage_cover = 0.95
read_size = 1000
window_size = read_size
# priorities = []
covered_percentage = 0

weight_window_sum = rolling_sum(full_data, 'weight', window_size, full_data['genome_pos'].tolist())
pos = full_data['genome_pos'].unique()
covered_positions = {}
covered_ranges = []
covered_snps = []
# for run in tqdm(range(0,read_number)):
run = 0
while covered_percentage < percentage_cover:
    print(f'{np.round(covered_percentage,3)}-->{percentage_cover}')
    start = pos[np.argmax(weight_window_sum)] # find the index of the max value in the rolling sum
    in_range = [i for i in pos if i <= start+window_size]
    end = min(in_range, key=lambda x:abs(x-(start+window_size)))

    covered_positions[f'Amplicon_{run+1}'] = {'Range':{'Start': start, 'End': end}}  # concise version of output
    covered_ranges.append([start, end])
    run +=1
    # 
    # full_data.loc[(full_data['genome_pos']>=start) & (full_data['genome_pos']<=end), 'weight'] = 0 # set the weight of the covered positions to 0
    full_data.loc[(full_data['genome_pos']>=start) & (full_data['genome_pos']<=end), 'weight'] = full_data.loc[(full_data['genome_pos']>=start) & (full_data['genome_pos']<=end), 'weight']/10 # set the weight of the covered positions to 0
    weight_window_sum = rolling_sum(full_data, 'weight', window_size, full_data['genome_pos'].tolist())
    covered_snps.extend(list(full_data.loc[(full_data['genome_pos']>=start) & (full_data['genome_pos']<=end), 'sample_id'].values))

    covered_snps = list(set(covered_snps))
    covered_percentage = len(covered_snps)/full_data['sample_id'].nunique()
print(f'{np.round(covered_percentage,3)}-->{percentage_cover}')
print(f'{len(covered_positions)} amplicons needed to cover {percentage_cover} of the knowns SNPs')

#%%
# output covered position into a bed file
bed_file_path = "intervals1kbps95.bed"

with open(bed_file_path, "w") as bed_file:
    for interval in covered_ranges:
        start = interval[0]
        end = interval[1]
        bed_line = f"chr1\t{start}\t{end}\n"  # Modify "chr1" with the appropriate chromosome name
        bed_file.write(bed_line)
            
# %%
list(full_data.loc[(full_data['genome_pos']>=1) & (full_data['genome_pos']<=10000), 'sample_id'].values)



# %%
# Extracting the sequence from a ref genome


def extract_sequence_from_fasta(fasta_file, sequence_id, start_pos, end_pos):
    """
    Extracts a subsequence from a FASTA file based on the given sequence ID, start position, and end position.
    """
    # Iterate over the sequences in the FASTA file
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Check if the current sequence ID matches the desired sequence ID
        
        if record.id == sequence_id:
            print(record.id)
            # Extract the subsequence based on the start and end positions
            subsequence = record.seq[start_pos:end_pos]
            return str(subsequence)  # Return the subsequence as a string

    # If the sequence ID is not found, return None
    return None


# Example usage
fasta_file = '/mnt/storage10/lwang/Projects/Amplicone_design_tool/MTB-h37rv_asm19595v2-eg18.fa'
sequence_id = 'Chromosome'
start_pos = 0
end_pos = 20

# Extract the sequence
sequence = extract_sequence_from_fasta(fasta_file, sequence_id, start_pos, end_pos)
print(MeltingTemp.Tm_NN(sequence))

# %%
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



# %%

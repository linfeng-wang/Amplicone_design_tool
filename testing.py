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



# full_data = pd.read_csv('/mnt/storage10/jody/projects/variant_dump/variants.csv')
# full_data = full_data[~full_data['drugs'].isna()]
# full_data = full_data.sort_values(by=['genome_pos'])
# full_data = full_data.reset_index(drop=True)
# full_data['weight'] = full_data['freq']
#%%
# snp_weight = full_data['sample_id'].value_counts()/full_data['sample_id'].value_counts().max()
# snp_weight = snp_weight.to_frame()
#%%
# Test
# full_data['weight'] = 0
# full_data = full_data.merge(snp_weight, left_on='sample_id', right_index=True, how='left')
# merged_data = merged_data.sort_values(by=['genome_pos'])
# merged_data = merged_data.reset_index(drop=True)
# merged_data['weight'].fillna(value=0, inplace=True)
def test_coverage(bed_file_path, full_data,tb_genes):
    intervals = []
    with open(bed_file_path, "r") as bed_file:
        for line in bed_file:
            line = line.strip().split("\t")
            start = int(line[1])
            end = int(line[2])
            intervals.append([start, end])
    # Print the intervals
    destination_df = pd.DataFrame(columns=full_data.columns)
    for interval in intervals:
        segment = full_data[(full_data['genome_pos']>= interval[0])&(full_data['genome_pos']<= interval[1])]
        destination_df = pd.concat([destination_df, segment])
    destination_df = destination_df.drop_duplicates(keep='first')
    
    for i, x in full_data['gene'].value_counts().items():
        if i in destination_df['gene'].value_counts().keys():
            if i in tb_genes:
                print(i, np.round(destination_df['gene'].value_counts()[i]/full_data['gene'].value_counts()[i],2), '<--')
            else:
                print(i, np.round(destination_df['gene'].value_counts()[i]/full_data['gene'].value_counts()[i],2))
        else:
            if i in tb_genes:
                print(i, 0, '<--')
            else:
                print(i, 0)

# test_coverage('/mnt/storage10/lwang/Projects/Amplicone_design_tool/intervals-30-1000.bed', full_data)
# %% debug
# for x in intervals:
#     print(full_data[(full_data['genome_pos']>= x[0])&(full_data['genome_pos']<= x[1])])
#     break
# %%
def calculate_overlap(ranges):
    total_overlap = 0

    for i in range(len(ranges)):
        start_i, end_i = ranges[i]

        for j in range(i + 1, len(ranges)):
            start_j, end_j = ranges[j]

            overlap = max(0, min(end_i, end_j) - max(start_i, start_j))
            total_overlap += overlap

    return total_overlap

# ranges = [[760314, 761314], [2154853, 2155853], [2288257, 2289257], [4247186, 4248186],
#           [781585, 782585], [1472358, 1473358], [1673373, 1674373], [6620, 7620],
#           [4407543, 4408543], [2715339, 2716339], [761089, 762089], [2154169, 2155169],
#           [4243203, 4244203], [2288257, 2289257], [4247186, 4248186], [781585, 782585],
#           [1472272, 1473272], [1674048, 1675048], [764363, 765363], [2726112, 2727112],
#           [6620, 7620], [1673373, 1674373], [4248747, 4249747], [1917939, 1918939],
#           [4407543, 4408543], [2155125, 2156125], [2715339, 2716339], [4326081, 4327081],
#           [760314, 761314], [6575, 7575]]

# overlap = calculate_overlap(ranges)
# print("Total overlap:", overlap)

# %%primer found

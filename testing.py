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
bed_file_path = "intervals.bed"
intervals = []

with open(bed_file_path, "r") as bed_file:
    for line in bed_file:
        line = line.strip().split("\t")
        start = int(line[1])
        end = int(line[2])
        intervals.append([start, end])
#%%
# Print the intervals
destination_df = pd.DataFrame()

for interval in intervals:
    segment = full_data[(full_data['genome_pos']>= interval[0])&(full_data['genome_pos']<= interval[1])]
    destination_df = pd.concat([destination_df, segment])
destination_df = destination_df.drop_duplicates(keep='first')

# %%
destination_df['gene'].value_counts()

# %%
full_data['gene'].value_counts()

# %%
for i, x in full_data['gene'].value_counts().items():
    if i in destination_df['gene'].value_counts().keys():
        print(i, np.round(destination_df['gene'].value_counts()[i]/full_data['gene'].value_counts()[i],2))
    else:
        print(i, 0)
        

#%%
for i, x in destination_df['gene'].value_counts().items():
    print(i,x)
# %%

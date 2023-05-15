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

#%%
#Full data trial
read_number = 40
read_size = 400
window_size = 400
# priorities = []

weight_window_sum = rolling_sum(full_data, 'weight', 400, full_data['genome_pos'].tolist())
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
    full_data.loc[(full_data['genome_pos']>=start) & (full_data['genome_pos']<=end), 'weight'] = 0 # set the weight of the covered positions to 0
    weight_window_sum = rolling_sum(full_data, 'weight', 400, full_data['genome_pos'].tolist())

#%% checking
full_data[(full_data['genome_pos']>=761004) & (full_data['genome_pos']<761284)]
full_data[full_data['genome_pos']<400]
full_data[(full_data['genome_pos']>=183247) & (full_data['genome_pos']<=761284)]
#version to fall back on
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
from plotly import graph_objects as go
import json
from imp import reload
import primer_selection
reload(primer_selection)
import testing
reload(testing)
#%%
def value_counts_list(lst):
    """
    Computes the frequency count of unique elements in a list and returns a dictionary, sorted by frequency count in
    descending order.

    Args:
    - lst (list): List of elements

    Returns:
    - dict: Dictionary with unique elements as keys and their frequency count as values, sorted by frequency count
    in descending order
    """
    value_counts = {}
    for item in lst:
        if item in value_counts:
            value_counts[item] += 1
        else:
            value_counts[item] = 1
    sorted_value_counts = dict(sorted(value_counts.items(), key=lambda x: x[1], reverse=True))
    return sorted_value_counts

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
tb_drug_resistance_genes = {
    'gyrB': ['Levofloxacin'],
    'gyrA': ['Levofloxacin'],
    'mshA': ['Isoniazid'],
    'rpoB': ['Rifampicin'],
    'rpoC': ['Rifampicin'],
    'rpsL': ['Streptomycin'],
    'embR': ['Ethambutol'],
    'rrs': ['Kanamycin', 'Capreomycin', 'Amikacin', 'Streptomycin'],
    'fabG1': ['Isoniazid'],
    'inhA': ['Isoniazid'],
    'rpsA': ['Pyrazinamide'],
    'tlyA': ['Capreomycin'],
    'ndh': ['Isoniazid'],
    'katG': ['Isoniazid'],
    'pncA': ['Pyrazinamide'],
    'kasA': ['Isoniazid'],
    'eis': ['Kanamycin', 'Amikacin'],
    'ahpC': ['Isoniazid'],
    'rpoA': ['Rifampicin'],
    'panD': ['Pyrazinamide'],
    'embC': ['Ethambutol'],
    'embA': ['Ethambutol'],
    'embB': ['Ethambutol'],
    'ubiA': ['Ethambutol'],
    'gid': ['Streptomycin']
}

full_data.loc[full_data['gene'].isin(tb_drug_resistance_genes.keys()), 'weight'] += 0.5
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
# full_data['weight'].min()
#%%
#Full data trial
read_number = 30
read_size = 1000
window_size = read_size-150
# priorities = []
graph_output = False # set to True to see the graph output
weight_window_sum = rolling_sum(full_data, 'weight', window_size, full_data['genome_pos'].tolist())
pos = full_data['genome_pos'].unique()
covered_positions = {}
covered_ranges = []
for run in tqdm(range(0,read_number)):
    while graph_output:
        trace = go.Scatter(
        x=list(range(1, len(weight_window_sum) + 1)),
        y=weight_window_sum,
        mode='lines',
        line=dict(color='blue'),
        fill='tozeroy',
        fillcolor='rgba(0, 0, 255, 0.3)')
        # Create the layout
        layout = go.Layout(
            title='Line Density Plot',
            xaxis=dict(title='Genomic Position'),
            yaxis=dict(title='Density'),
            shapes=[
            # Add a vertical line at x=8
            dict(
                type='line',
                x0=np.argmax(weight_window_sum),
                x1=np.argmax(weight_window_sum),
                y0=0,
                y1=max(weight_window_sum),
                line=dict(color='red', width=2)
            )
        ]
        )
        # Create the figure
        fig = go.Figure(data=[trace], layout=layout)
        # Display the plot
        fig.show()
    
    start = pos[np.argmax(weight_window_sum)] # find the index of the max value in the rolling sum
    # in_range = [i for i in pos if i <= start+window_size] # find all values in the window
    # end = min(in_range, key=lambda x:abs(x-(start+window_size))) # find the closest value to the end of the window
    end = start+window_size
    if end > 4485058:
        end = 4485058-200
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
    # covered_positions[f'Amplicon_{run+1}'] = {'Range':{'Start': start, 'End': end}, 'Markers':full_data[['genome_pos','gene','sublin','drtype','drugs','weight']][start:end].sort_values(by=['weight']).to_dict('records')} 
    covered_positions[f'Amplicon_{run+1}'] = {'Range':{'Start': start, 'End': end}, 'Markers':full_data[(full_data['genome_pos']>= start) & (full_data['genome_pos']<=end)][['genome_pos','gene','sublin','drtype','drugs','weight']].sort_values(by=['weight']).to_dict('records')}# verbose version of output
    # covered_positions[f'Amplicon_{run+1}'] = {'Range':{'Start': start, 'End': end}}  # concise version of output
    covered_ranges.append([start, end])
    # 
    # full_data.loc[(full_data['genome_pos']>=start) & (full_data['genome_pos']<=end), 'weight'] = 0 # set the weight of the covered positions to 0
    full_data.loc[(full_data['genome_pos']>=start) & (full_data['genome_pos']<=end), 'weight'] = full_data.loc[(full_data['genome_pos']>=start) & (full_data['genome_pos']<=end), 'weight']/10 # set the weight of the covered positions to 0
    weight_window_sum = rolling_sum(full_data, 'weight', window_size, full_data['genome_pos'].tolist())
#! full data needs to be reloaded after each run
#%%
# output
seq = 1
padding = 250
accepted_primers = pd.DataFrame(columns=['pLeft_ID', 'pLeft_coord', 'pLeft_length', 'pLeft_Tm', 'pLeft_GC', 'pLeft_Sequences', 'pLeft_EndStability','pRight_ID', 'pRight_coord', 'pRight_length', 'pRight_Tm', 'pRight_GC', 'pRight_Sequences', 'pRight_EndStability', 'Penalty', 'Product_size'])
primer_pool = []
for i, x in enumerate(covered_ranges):
    low_b = x[0]
    high_b = x[1]
    if low_b <= 150:
        low_b = 151
    elif high_b >= 4411532:
        high_b = 4411532-151
    # if high_b - low_b < 300:
    #     high_b+= 150
    #     low_b-= 150
    if (high_b - low_b) < 350:
        # print('======')
        high_b+= 450

    seq_template = primer_selection.extract_sequence_from_fasta(low_b, high_b, padding=padding)
    
    # print(low_b, high_b)
    # print(seq_template)
    # print(x)
    primer_pool, accepted_primers = primer_selection.result_extraction(primer_pool, accepted_primers, seq_template, i+1, padding, low_b, high_b)
    # print(accepted_primers)
    primer_pos = accepted_primers.iloc[-1][['pLeft_coord','pRight_coord']].values
    amplified_segment = full_data[(full_data['genome_pos']>= primer_pos[0])&(full_data['genome_pos']<= primer_pos[1])]
    covered_segment = full_data[(full_data['genome_pos']>= low_b)&(full_data['genome_pos']<= high_b)]
    # print('======')
    # print(amplified_segment.shape[0]/(covered_segment.shape[0]+0.00001))
    # ratio = amplified_segment.shape[0]/(covered_segment.shape[0]+0.00001)
    # # if ratio < 0.85:
        
    # print('======')
    
    # if amplified_segment.shape[0]/covered_segment.shape[0] < 0.5:
    #     print()
    # break
accepted_primers.to_csv(f'output/accepted_primers-{read_number}-{read_size}.csv')

# output covered position info into a json file
class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)

# output covered coordinates into a bed file
bed_file_path = f"output/intervals-{read_number}-{read_size}.bed"

with open(bed_file_path, "w") as bed_file:
    for start, end in zip(accepted_primers['pLeft_coord'],accepted_primers['pRight_coord']+accepted_primers['pRight_length']):
        bed_line = f"chr1\t{start}\t{end}\n"  # Modify "chr1" with the appropriate chromosome name
        bed_file.write(bed_line)

run = 0
covered_positions = {}
for start, end in zip(accepted_primers['pLeft_coord'],accepted_primers['pRight_coord']+accepted_primers['pRight_length']):
    # print(start, end)
    covered_positions[f'Amplicon_{run+1}'] = {'Range':{'Start': start, 'End': end}, 'Markers':full_data[(full_data['genome_pos']>= start) & (full_data['genome_pos']<=end)][['genome_pos','gene','change','drugs','weight']].sort_values(by=['weight']).to_dict('records')}# verbose version of output
    run+=1
out_file = open(f"output/covered_positions-{read_number}-{read_size}.json", "w")

json.dump(covered_positions, out_file, cls=NpEncoder, indent = 6)

out_file.close()
print(f'Output_file:{read_number}-{read_size}')

testing.test_coverage(bed_file_path, full_data, tb_drug_resistance_genes.keys())



#%%
tb_genes = tb_drug_resistance_genes.keys()
destination_df = pd.DataFrame(columns=full_data.columns)
for interval in covered_ranges:
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












#%%
run = 0
covered_positions = {}
for start, end in zip(accepted_primers['pLeft_coord'],accepted_primers['pRight_coord']+accepted_primers['pRight_length']):
    print(start, end)
    covered_positions[f'Amplicon_{run+1}'] = {'Range':{'Start': start, 'End': end}, 'Markers':full_data[(full_data['genome_pos']>= start) & (full_data['genome_pos']<=end)][['genome_pos','gene','sublin','drtype','drugs','weight']].sort_values(by=['weight']).to_dict('records')}# verbose version of output
    run+=1
#%%
covered_positions


#%% checking
full_data[(full_data['genome_pos']>=761004) & (full_data['genome_pos']<761284)]
full_data[full_data['genome_pos']<window_size]
full_data[(full_data['genome_pos']>=183247) & (full_data['genome_pos']<=761284)]

#%% checking
# for x in covered_ranges:
#     print(x)
#     print(x[1]-x[0])

#%%
for i,e in covered_positions.items():
    print(i, e['Range']['Start'], e['Range']['End'])
    print(len(e['Markers']))
    # print(e['Markers'])
    weight = 0
    gene = []
    for x in e['Markers']:
        weight += x['weight']
        gene.append(x['gene'])
    print(weight)
    print(value_counts_list(gene))
    print('------------------')
#%%
# full_data[(full_data['genome_pos']>= 4247186) & (full_data['genome_pos']<=4248186)]['genome_pos','gene','sublin','drtype','drugs','weight'].sort_values(by=['weight'])
full_data[(full_data['genome_pos']>= 4247186) & (full_data['genome_pos']<=4248186)][['genome_pos','gene','sublin','drtype','drugs','weight']].sort_values(by=['weight']).to_dict('records')
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

# # Example usage
# fasta_file = '/mnt/storage10/lwang/Projects/Amplicone_design_tool/MTB-h37rv_asm19595v2-eg18.fa'
# sequence_id = 'Chromosome'
# start_pos = 0
# end_pos = 20



# %%
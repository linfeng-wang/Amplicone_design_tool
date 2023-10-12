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
import Amplicone_no
reload(Amplicone_no)
import argparse
from functools import reduce
import os

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
full_data = pd.read_csv('/mnt/storage10/jody/projects/variant_dump/variants.csv')
full_data = full_data[~full_data['drugs'].isna()]
full_data = full_data.sort_values(by=['genome_pos'])
full_data = full_data.reset_index(drop=True)
full_data['weight'] = full_data['freq']
ref_genome = 'MTB-h37rv_asm19595v2-eg18.fa'

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

def genome_size(fasta_file):
    total_length = 0
    with open(fasta_file, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                total_length += len(line.strip())
    return total_length


#Full data trial
def place_amplicone(full_data, read_number, read_size, graphic_output=False, ref_size = genome_size(ref_genome)):
    read_number = read_number
    read_size = read_size
    window_size = read_size
    run = 0
    full_data_cp = full_data.copy()
    # priorities = []
    graph_output = False # set to True to see the graph output
    weight_window_sum = rolling_sum(full_data_cp, 'weight', window_size, full_data_cp['genome_pos'].tolist())
    pos = full_data_cp['genome_pos'].unique()
    covered_positions = {}
    covered_ranges = []
    reduce_amplicone = 0
    print('Placing Amplicones...')
    # for run in tqdm(range(0,read_number)):
    while run < read_number:
        print(f'Amplicone #{run+1}')
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
        if end >  ref_size:
            end =  ref_size-200

        full_data_cp.loc[(full_data_cp['genome_pos']>=start) & (full_data_cp['genome_pos']<=end), 'weight'] = full_data_cp.loc[(full_data_cp['genome_pos']>=start) & (full_data_cp['genome_pos']<=end), 'weight']/10 # set the weight of the covered positions to 0
        if [start, end] in covered_ranges:
            # print(full_data_cp)
            # print(full_data_cp.loc[(full_data_cp['genome_pos']>=start) & (full_data_cp['genome_pos']<=end)]['weight'].values)
            print('Already covered, consider reducing amplicone number.. Finding alternative sequence...')
            reduce_amplicone += 1
        else:
            run += 1
            covered_ranges.append([start, end])
            covered_positions[f'Amplicon_{run+1}'] = {'Range':{'Start': start, 'End': end}, 'Markers':full_data_cp[(full_data_cp['genome_pos']>= start) & (full_data_cp['genome_pos']<=end)][['genome_pos','gene','sublin','drtype','drugs','weight']].sort_values(by=['weight']).to_dict('records')}# verbose version of output
            # print(full_data_cp.loc[(full_data_cp['genome_pos']>=start) & (full_data_cp['genome_pos']<=end)]['weight'].values)
            # print('==============')
        weight_window_sum = rolling_sum(full_data_cp, 'weight', window_size, full_data_cp['genome_pos'].tolist())
    print('====================')
    print(f'Consider reducing number of amplicones by: {reduce_amplicone}')
    print('====================')

    return covered_positions, covered_ranges

#%%
read_size = 1000
specific_gene = ['katG']
specific_gene_amplicone = 5
non_specific_amplicone = 25
specific_gene_data = full_data[full_data['gene'].isin(specific_gene)]
non_specific_gene_data = full_data[~full_data['gene'].isin(specific_gene)]

covered_positions, covered_ranges = [], []
if len(specific_gene)>0:
    covered_positions_sp, covered_ranges_sp = place_amplicone(specific_gene_data, specific_gene_amplicone, read_size, genome_size(ref_genome))
    covered_positions_nosp, covered_ranges_nosp = place_amplicone(non_specific_gene_data, non_specific_amplicone, read_size, genome_size(ref_genome))
    covered_positions = {**covered_positions_sp, **covered_positions_nosp}
    covered_ranges = covered_ranges_sp + covered_ranges_nosp

else:
    covered_positions_nosp, covered_ranges_nosp = place_amplicone(non_specific_gene_data, non_specific_amplicone, read_size, genome_size(ref_genome))
    covered_positions = covered_positions_nosp
    covered_ranges = covered_ranges_nosp

#%%
spoligotype = True
if spoligotype:
    spacers = pd.read_csv('spacers.bed', sep='\t', header=None)
    spacers = np.array(spacers)
    spacers = spacers[:, 1:3]
    spacers = spacers.tolist()
    flattened_data = [item for sublist in spacers for item in sublist]
    spacer_max = max(flattened_data)
    spacer_min = min(flattened_data)
    spol_list = np.arange(spacer_min-400,spacer_max+400,1)
    weight = [0.01]*len(spol_list) 
    spol_data = pd.DataFrame({'genome_pos':spol_list,'weight':weight})
    # Create a list of boolean masks, one for each range
    masks = [(spol_data['genome_pos'] >= start) & (spol_data['genome_pos'] <= end) for start, end in spacers]
    # Use reduce and the | operator to combine the masks into a single mask
    combined_mask = reduce(lambda x, y: x | y, masks)
    # Use .loc and the combined mask to update the weight column
    spol_data.loc[combined_mask, 'weight'] = 1
    covered_ranges_spol  = Amplicone_no.place_amplicone_spol(spol_data, 1, read_size, graphic_output=False, ref_size = genome_size(ref_genome))
    covered_ranges.extend(covered_ranges_spol)
    read_number = specific_gene_amplicone + non_specific_amplicone + len(covered_ranges_spol)

#%%
# output
# seq = 1
output_path = '.'
op = f'{output_path}/Amplicone_design_output'
os.makedirs(op, exist_ok=True) #output path
ref_size = genome_size(ref_genome)
padding = 150
accepted_primers = pd.DataFrame(columns=['pLeft_ID', 'pLeft_coord', 'pLeft_length', 'pLeft_Tm', 'pLeft_GC', 'pLeft_Sequences', 'pLeft_EndStability','pRight_ID', 'pRight_coord', 'pRight_length', 'pRight_Tm', 'pRight_GC', 'pRight_Sequences', 'pRight_EndStability', 'Penalty', 'Product_size'])
primer_pool = []
no_primer = []
print('Designing primers...')
for i, x in tqdm(enumerate(covered_ranges)):
    low_b = x[0]
    high_b = x[1]
    if low_b <= padding:
        low_b = 151
    elif high_b >= ref_size-padding:
        high_b = ref_size-padding
    # if high_b - low_b < 300:
    #     high_b+= 150
    #     low_b-= 150
    if (high_b - low_b) < padding*2+50:
        # print('======')
        high_b+= 450
    else:
        high_b+= padding
        low_b-= padding

    # seq_template = primer_selection.extract_sequence_from_fasta(low_b, high_b, padding=padding)
    seq_template = primer_selection.extract_sequence_from_fasta(low_b, high_b, padding=0) #set to 0 to avoid padding again
    # print(low_b, high_b)
    # print(seq_template)
    # print(x)
    # print(seq_template)
    primer_pool, accepted_primers, no_primer_ = primer_selection.result_extraction(primer_pool, accepted_primers, seq_template, i+1, padding, ref_genome = ref_genome)

    no_primer.extend(no_primer_)
    
primer_label = ['Gene_specific']*specific_gene_amplicone + ['Non_specific']*non_specific_amplicone + ['Spoligotype']*len(covered_ranges_spol)

for x in no_primer:
    {f'Designed covered_ranges[{x}] cannot be found'}
    my_list.pop(x)

accepted_primers['Amplicone_type'] = primer_label
accepted_primers.to_csv(f'{op}/Primer_design-accepted_primers-{read_number}-{read_size}.csv')

primer_pos = accepted_primers[['pLeft_coord','pRight_coord']].values
columns = ['pLeft_ID', 'pRight_ID', 'pLeft_coord', 'pRight_coord', 'SNP_inclusion']

# Create an empty DataFrame with the specified column headings
primer_inclusion =pd.DataFrame(columns=columns)
for i, row in accepted_primers.iterrows():
    data = full_data[(full_data['genome_pos']>= row['pLeft_coord']) & (full_data['genome_pos']<= row['pRight_coord'])]    
    info = row[['pLeft_ID', 'pRight_ID', 'pLeft_coord', 'pRight_coord']]
    SNP = data['gene'].str.cat(data['change'], sep='-').unique()
    info['SNP_inclusion'] = ','.join(SNP)
    primer_inclusion.loc[len(primer_inclusion)] = info.tolist()

primer_inclusion.to_csv(f'{op}/SNP_inclusion.csv')
#%%
# Evaluation
    # print(accepted_primers)
columns = ['sample_id', 'genome_pos', 'gene', 'change', 'freq', 'type', 'sublin', 'drtype', 'drugs', 'weight']

# Create an empty DataFrame with the specified column headingsref_genome
designed =pd.DataFrame(columns=columns)

# Create DataFrame
covered = pd.DataFrame(columns=columns)
primer_pos = accepted_primers[['pLeft_coord','pRight_coord']].values
for x in primer_pos:
    covered = pd.concat([covered, full_data[(full_data['genome_pos']>= x[0])&(full_data['genome_pos']<= x[1])]])
for x in covered_ranges:
    designed = pd.concat([covered, full_data[(full_data['genome_pos']>= x[0])&(full_data['genome_pos']<= x[1])]])

print(designed['gene'].value_counts())
print(covered['gene'].value_counts())
print(covered['gene'].value_counts()/designed['gene'].value_counts())

#%%
columns = ['sample_id', 'genome_pos', 'gene', 'change', 'freq', 'type', 'sublin', 'drtype', 'drugs', 'weight']
# Create an empty DataFrame with the specified column headings
tested =pd.DataFrame(columns=columns)

variants_dr = [] # all the drugs that have variants
tested_dr = [] # all the drugs that have been covered by the primers
variants = pd.read_csv('/mnt/storage10/lwang/Projects/Amplicone_design_tool/model/variants.txt') # put in the reference genome to test for
variants = variants[~variants['drugs'].isna()]

for x in variants['drugs'].values:
    variants_dr.extend(x.split(','))

inclusion  = []
for x in primer_pos:
    drop_ind = variants[(variants['genome_pos']>= x[0])&(variants['genome_pos']<= x[1])].index
    inclusion.append(len(drop_ind))
    tested = pd.concat([tested, variants[(variants['genome_pos']>= x[0])&(variants['genome_pos']<= x[1])]])
    variants.drop(drop_ind, inplace=True)
    # print(variants.shape)

for x in tested['drugs'].values:
    tested_dr.extend(x.split(','))
    
variants_dr_ = value_counts_list(variants_dr)
tested_dr_ = value_counts_list(tested_dr)

dr = []
ratio = []
percent = []
for x,y in variants_dr_.items():
    dr.append(x)
    if x in tested_dr_.keys():
        percent.append(round(tested_dr_[x]/variants_dr_[x]*100,1))
        # print(tested_dr_[x]/variants_dr_[x])
        ratio.append(f'{tested_dr_[x]}/{variants_dr_[x]}')
    else:
        percent.append(0/variants_dr_[x])
        ratio.append(f'0/{variants_dr_[x]}')

evaluation_df = pd.DataFrame({'Drug':dr, 'Ratio detected':ratio, 'Percentage detected':percent})
evaluation_df.to_csv(f'{op}/Primer_design-Percentage_gene_covered.csv')

#%%
def main():
    print(args.argument_name)

# %%
if __name__ == "_main__":
    parser = argparse.ArgumentParser(description='Amplicone design',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # in
    parser.add_argument('--SNP_priority', type = str, help = 'SNP_priority CSV files', default='variants.csv')
    parser.add_argument('--Amplicon_size', type = int, help = 'Amplicon size', default=1000)
    parser.add_argument('--Reference_genome', type = str, help = 'reference file', default='MTB-h37rv_asm19595v2-eg18.fa')
    parser.add_argument('--specific_amplicone_no', type = int, help = 'number of amplicone dedicated to amplifying specific genes', default=0)
    parser.add_argument('--specific_amplicone_gene', type = str, help = 'give a list of fene names separated by Lineage ', default='')
    parser.add_argument('--Non_specific_amplicone_no', type = int, help = 'number of amplicone dedicated to amplifying all SNPs in all genes according the importants list', default=30)
    
    parser.add_argument('--Spoligo_sequencing', type = bool, help = 'Whether to amplify Spoligotype', default=False)
    # out
    parser.add_argument('--output_folder_path', default = '', type = str, help = 'output_folder_path (accepted_primers, SNP_inclusion, gene_covered)')
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)
    main(args)
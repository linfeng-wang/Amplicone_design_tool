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
import time
import plotting
reload(plotting)
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


full_data = pd.read_csv('/mnt/storage10/jody/projects/variant_dump/variants.csv')
full_data = full_data[~full_data['drugs'].isna()]
full_gene = full_data[~full_data['type'].isin(['synonymous_variant','non_coding_transcript_exon_variant'])]
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
full_data.to_csv('snp_priority.csv', index=False)
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
    print(len(pos))
    for x in tqdm(pos):
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

def extraction_prep(x, ref_size= genome_size(ref_genome), padding=150):
    padding = int(padding/2)
    low_b = x[0]
    high_b = x[1]
    if low_b <= padding:
        low_b = padding+1
    elif high_b >= ref_size-padding:
        high_b = ref_size-padding
    # if high_b - low_b < 300:
    #     high_b+= 150
    #     low_b-= 150
    elif (high_b - low_b) < padding*2+50:
        # print('======')
        high_b+= padding*3
    else:
        high_b+= padding
        low_b-= padding
    # seq_template = primer_selection.extract_sequence_from_fasta(low_b, high_b, padding=padding)
    seq_template = primer_selection.extract_sequence_from_fasta(low_b, high_b, padding=padding/2) #set to 0 to avoid padding again
    # print(low_b, high_b)
    # print(seq_template)
    return seq_template, low_b, high_b

#%%
ideal_range = []
#place_amplicone function
def place_amplicone(full_data, read_number, read_size, primer_pool, accepted_primers, no_primer_, graphic_output=False, ref_size = genome_size(ref_genome),padding=150, output_path = '.'):
    start = time.time()

    read_number = read_number
    read_size = read_size
    window_size = read_size
    run = 0
    full_data_cp = full_data.copy()
    # priorities = []
    weight_window_sum = rolling_sum(full_data_cp, 'weight', window_size, full_data_cp['genome_pos'].tolist())
    pos = full_data_cp['genome_pos'].unique()
    covered_positions = {}
    covered_ranges = []
    reduce_amplicone = 0
    
    # primer design storage
    accepted_primers = pd.DataFrame(columns=['pLeft_ID', 'pLeft_coord', 'pLeft_length', 'pLeft_Tm', 'pLeft_GC', 'pLeft_Sequences', 'pLeft_EndStability','pRight_ID', 'pRight_coord', 'pRight_length', 'pRight_Tm', 'pRight_GC', 'pRight_Sequences', 'pRight_EndStability', 'Penalty', 'Product_size'])
    primer_pool = []
    no_primer_ = no_primer_
    print('Placing Amplicones...')
    # for run in tqdm(range(0,read_number)):
    while run < read_number:
        print(f'Amplicone #{run+1}')
        print(f'Designing primers...for {read_size}bps Amplicons...with {padding}bps padding')

        if graphic_output == True:
            output_path = '.'
            op = f'{output_path}/Running_graphs'
            os.makedirs(op, exist_ok=True) #output path
            # print(np.argmax(weight_window_sum))

            trace = go.Scatter(
            x=[item*window_size for item in list(range(1, len(weight_window_sum) + 1))],
            y=weight_window_sum,
            mode='lines',
            line=dict(color='blue'),
            fill='tozeroy',
            fillcolor='rgba(0, 0, 255, 0.3)')
            # Create the layout§
            layout = go.Layout(
                title=f'Weight sum calculation by {window_size}bps sliding windows - #{run+1}',
                xaxis=dict(title='Sliding Window Genomic Position(bps)'),
                yaxis=dict(title='Window Weight Sum'),
                shapes=[
                # Add a vertical line at x=8
                dict(
                    type='line',
                    x0=np.argmax(weight_window_sum)*window_size,
                    x1=np.argmax(weight_window_sum)*window_size+window_size,
                    y0=0,
                    y1=max(weight_window_sum),
                    line=dict(color='rgba(255, 0, 0, 0.5)', width=20)
                )
            ]
            )
            # Create the figure
            fig = go.Figure(data=[trace], layout=layout)
            # Display the plot
            fig.show()
            fig.write_image(f"{op}/{window_size}bps sliding windows-#{run+1}.png")

        
        start_r = pos[np.argmax(weight_window_sum)] # find the index of the max value in the rolling sum
        # in_range = [i for i in pos if i <= start+window_size] # find all values in the window
        # end = min(in_range, key=lambda x:abs(x-(start+window_size))) # find the closest value to the end of the window
        # print(start_r)
        end_r = start_r+window_size
        if end_r > ref_size:
            end_r = ref_size-200
        
        ideal_range.append([start_r, end_r])

        seq_template, low_b, high_b = extraction_prep([start_r, end_r], ref_size=genome_size(ref_genome), padding=padding)
        primer_pool, accepted_primers, no_primer = primer_selection.result_extraction(primer_pool, accepted_primers, seq_template, run+1, padding, ref_genome = ref_genome, high_b = high_b, low_b = low_b, read_size = read_size)
        if accepted_primers.shape[0] != 0:
            start_p, end_p = accepted_primers.iloc[accepted_primers.shape[0]-1][['pLeft_coord','pRight_coord']].values
        else:
            run= max(0,run-1)
            print('No suitable primers found')
            break
        # print('---WindowMin', min(weight_window_sum))
        # print('---WindowMean', np.mean(weight_window_sum))
        # print('---Min', full_data_cp['weight'].min())
        # print('---Mean', full_data_cp['weight'].mean())
        # c = full_data_cp[(full_data_cp['genome_pos']>=start_p) & (full_data_cp['genome_pos']<=end_p)].shape[0]
        c = full_data_cp.shape[0]
        # full_data_cp.loc[(full_data_cp['genome_pos']>=start_p) & (full_data_cp['genome_pos']<=end_p), 'weight'] = full_data_cp['weight'].min()/10/c  # set the weight of the covered positions smaller
        full_data_cp.loc[(full_data_cp['genome_pos']>=start_p) & (full_data_cp['genome_pos']<=end_p), 'weight'] = 0 # set the weight of the covered positions smaller
        # print(start_p, end_p)
        # full_data_cp.loc[(full_data_cp['genome_pos']>=start_p) & (full_data_cp['genome_pos']<=end_p), 'weight'] = full_data_cp.loc[(full_data_cp['genome_pos']>=start_p) & (full_data_cp['genome_pos']<=end_p), 'weight']/100 # set the weight of the covered positions smaller
        # print(full_data_cp[(full_data_cp['genome_pos']>=start_p) & (full_data_cp['genome_pos']<=end_p)]['weight'])
        # print('min',full_data_cp[(full_data_cp['genome_pos']>=start_p) & (full_data_cp['genome_pos']<=end_p)]['genome_pos'].min())
        # # print(full_data_cp.loc[(full_data_cp['genome_pos']>=start_p) & (full_data_cp['genome_pos']<=end_p), 'weight'].sum())
        # print(f'Cover ranges {run+1}/{read_number}: {[start_p, end_p]}==========')
# full_data_cp[(full_data_cp['genome_pos']>=1673440) & (full_data_cp['genome_pos']<=1674487)]
# full_data_cp[(full_data_cp['genome_pos']>=1673373) & (full_data_cp['genome_pos']<=1674373)]
        if [start_p, end_p] in covered_ranges:
            if pos[np.argmax(weight_window_sum)] == start_r:
                # c = full_data_cp[(full_data_cp['genome_pos']>=start_r) & (full_data_cp['genome_pos']<=end_r)].shape[0]
                c = full_data_cp.shape[0]
                
                # full_data_cp.loc[(full_data_cp['genome_pos']>=start_r) & (full_data_cp['genome_pos']<=end_r), 'weight'] = full_data_cp['weight'].min()/10/c # set the weight of the covered positions smaller
                full_data_cp.loc[(full_data_cp['genome_pos']>=start_r) & (full_data_cp['genome_pos']<=end_r), 'weight'] = 0 # set the weight of the covered positions smaller
            # this is when problem comes, there is a difference in range coverage according to the design by weighted sum, however the actual range obtained from designed primers are dont cover the same range, hence the sae primers are repeatedly designed 
            print('Already covered, consider reducing amplicone number...Finding alternative sequence...')
            reduce_amplicone += 1
            accepted_primers = accepted_primers.iloc[:-1]
            print('***')
            
        else:
            run += 1
            covered_ranges.append([start_p, end_p])
            no_primer_.append(no_primer)

            covered_positions[f'Amplicon_{run+1}'] = {'Range':{'Start': start_p, 'End': end_p}, 'Markers':full_data_cp[(full_data_cp['genome_pos']>=start_p) & (full_data_cp['genome_pos']<=end_p)][['genome_pos','gene','sublin','drtype','drugs','weight']].sort_values(by=['weight']).to_dict('records')}# verbose version of output
            print('***')

            # print(full_data_cp.loc[(full_data_cp['genome_pos']>=start_p) & (full_data_cp['genome_pos']<=end_p)]['weight'].values)
            # print('==============')
        # print(weight_window_sum)
        weight_window_sum = rolling_sum(full_data_cp, 'weight', window_size, full_data_cp['genome_pos'].tolist())
    if reduce_amplicone > 0:
        print('====================')
        print(f'Consider reducing number of amplicones by: {reduce_amplicone}')
        print('====================')
        end = time.time()
        print(f'Programme ran for {round((end - start)/60,1)} min')

    return covered_positions, covered_ranges, full_data_cp, primer_pool, accepted_primers, no_primer_


#%%
# read_size = 1000
read_size = 400
# specific_gene = ['katG']
specific_gene = []
# padding = int(read_size/10)
# padding = 152
padding = 100
# read_size = read_size - padding*2
non_specific_gene = []
specific_gene_amplicone = 0
non_specific_amplicone = 35
# this is way we separate the specific and non-specific amplicones
# specific_gene_data = full_data[full_data['gene'].isin(specific_gene)]
# non_specific_gene_data = full_data[~full_data['gene'].isin(specific_gene)]
# this way we still have non specific amplicones incorporate the specific genes
specific_gene_data = full_data[full_data['gene'].isin(specific_gene)]
non_specific_gene_data = full_data
ref_size = genome_size(ref_genome)
specific_gene_data_count = 0
output_path = '.'


# Calculating number of amplicone needed
# target_coverage = 1
# gene_coverage = Amplicone_no.place_amplicone_search(full_data, target_coverage, read_size, genome_size(ref_genome))

covered_positions, covered_ranges = [], []
primer_pool, no_primer_ = [], []
accepted_primers = pd.DataFrame(columns = ['pLeft_ID', 'pLeft_coord', 'pLeft_length', 'pLeft_Tm', 'pLeft_GC',
    'pLeft_Sequences', 'pLeft_EndStability', 'pRight_ID', 'pRight_coord',
    'pRight_length', 'pRight_Tm', 'pRight_GC', 'pRight_Sequences',
    'pRight_EndStability', 'Penalty', 'Product_size'])

if len(specific_gene)>0:
    print('=====Specific amplicon=====')
    covered_positions_sp, covered_ranges_sp, specific_gene_data, primer_pool, accepted_primers, no_primer_ = place_amplicone(specific_gene_data, specific_gene_amplicone, read_size, primer_pool, accepted_primers, no_primer_, genome_size(ref_genome),padding=padding, output_path =output_path)
    specific_gene_data_count = accepted_primers.shape[0]                                                  
    print('=====Non-specific amplicon=====')
    non_specific_gene_data.update(specific_gene_data) # add the specific gene data to the non-specific gene data
    covered_positions_nosp, covered_ranges_nosp, full_data_cp, primer_pool, accepted_primers, no_primer_ = place_amplicone(non_specific_gene_data, non_specific_amplicone, read_size, primer_pool, accepted_primers, no_primer_, genome_size(ref_genome),padding=padding, output_path =output_path)
    covered_positions = {**covered_positions_sp, **covered_positions_nosp}
    covered_ranges = covered_ranges_sp + covered_ranges_nosp
    non_specific_gene_data_count = accepted_primers.shape[0] - specific_gene_data_count
    
else:
    covered_positions_nosp, covered_ranges_nosp, full_data_cp, primer_pool, accepted_primers, no_primer_ = place_amplicone(non_specific_gene_data, non_specific_amplicone, read_size, primer_pool, accepted_primers, no_primer_,genome_size(ref_genome),padding=padding, output_path =output_path)
    covered_positions = covered_positions_nosp
    covered_ranges = covered_ranges_nosp
    non_specific_gene_data_count = accepted_primers.shape[0]
spoligotype = False
covered_ranges_spol = []
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
    covered_ranges_spol = Amplicone_no.place_amplicone_spol(spol_data, 1, read_size, graphic_output=False, ref_size = genome_size(ref_genome))
    covered_ranges.extend(covered_ranges_spol)

read_number = specific_gene_data_count + non_specific_gene_data_count + len(covered_ranges_spol)

# output
primer_label = ['Gene_specific']*specific_gene_data_count + ['Non_specific']*non_specific_gene_data_count + ['Spoligotype']*len(covered_ranges_spol)

accepted_primers['Amplicone_type'] = primer_label
accepted_primers['Redesign'] = no_primer_
accepted_primers['Designed_ranges'] = covered_ranges

op = f'{output_path}/Amplicon_design_output'
os.makedirs(op, exist_ok=True) #output path
accepted_primers.to_csv(f'{op}/Primer_design-accepted_primers-{read_number}-{read_size}.csv',index=False)

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

primer_inclusion.to_csv(f'{op}/SNP_inclusion-{read_number}-{read_size}.csv',index=False)

print('Primer design output files:')
print(f'{op}/SNP_inclusion-{read_number}-{read_size}.csv')
print(f'{op}/Primer_design-accepted_primers-{read_number}-{read_size}.csv')

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

columns = ['sample_id', 'genome_pos', 'gene', 'change', 'freq', 'type', 'sublin', 'drtype', 'drugs', 'weight']
# Create an empty DataFrame with the specified column headings
tested =pd.DataFrame(columns=columns)

variants_dr = [] # all the drugs that have variants
tested_dr = [] # all the drugs that have been covered by the primers
variants = pd.read_csv('/mnt/storage10/lwang/Projects/Amplicone_design_tool/model2in1/variants.txt') # put in the reference samples to test for
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


op = f'{output_path}/Eval'
os.makedirs(op, exist_ok=True) #output path

evaluation_df = pd.DataFrame({'Drug':dr, 'Ratio detected':ratio, 'Percentage detected':percent})
evaluation_df.to_csv(f'{op}/Primer_design-Percentage_gene_covered-{read_number}-{read_size}.csv',index=False)
print('Evalutation output files:')
print(f'{op}/Primer_design-Percentage_gene_covered-{read_number}-{read_size}.csv')

#%%
#Plotting coverage
plotting.gene_coverage_plot(read_size, full_data, gene_range='regions.bed', tb_drug_resistance_genes= tb_drug_resistance_genes, output_path = output_path)
    

#%%
def main(args):
    # Rationale:
    # This program is designed to perform primer design for PCR amplification,
    # taking into consideration various parameters such as amplicon size, SNP priority,
    # reference genome, and more. It also allows users to specify particular amplicons
    # based on gene names and offers graphical output options. Users can also choose to
    # include spoligotyping sequencing information. The program outputs results to a 
    # specified folder.

    # Displaying the User Settings for Verification:
    print("========== User Settings ==========")
    print(f"Amplicon Size: {args.amplicone_size}")
    print(f"SNP Priority: {args.snp_priority}")
    print(f"Reference Genome: {args.reference_genome}")
    if args.padding_size == None:
        print(f"Padding_size: {args.amplione_size/4}")
    else:
        print(f"Padding_size: {args.padding_size}")
    print(f"Specific Amplicon Number: {args.specific_amplicone_no}")
    if args.specific_amplicone_no > 0:
        print(f"Specific Amplicon Gene: {args.specific_amplicone_gene}")
    print(f"Non-specific Amplicon Number: {args.non_specific_amplicone_no}")
    print(f"Graphic Option: {args.graphic_option}")
    print(f"Spoligo Sequencing: {args.spoligo_sequencing}")
    if args.spoligo_sequencing:
        if args.spoligo_sequencing_file == None:
            print(f"Spoligo Sequencing File: default MTB spoligotype spacer range file")
        else:
            print(f"Spoligo Sequencing File: {args.spoligo_sequencing_file}")
    print(f"Output Folder Path: {args.output_folder_path}")
    
    print("==================================")


    # read_size
    read_size = args.amplicone_size
    # full data - priority file with weights&frequency for each snp
    if args.snp_priority:
        full_data= pd.read_csv(args.snp_priority)
    else:
        full_data = pd.read_csv('snp_priority.csv')
    # paddding size
    if args.padding_size == None:
        padding = int(read_size/4)
    else:
        padding = args.padding_size

    # Reference Genome
    ref_genome = args.reference_genome
    
    #specific_genes
    if args.pecific_amplicone_gene:
            specific_gene = args.pecific_amplicone_gene.split(',')
    else:
        specific_gene = []
        
    non_specific_gene = []
    # specific_gene_amplicone
    specific_gene_amplicone = args.specific_amplicone_no
    # non_specific_amplicone
    non_specific_amplicone = args.non_specific_amplicone_no
    # this is way we separate the specific and non-specific amplicones
    # specific_gene_data = full_data[full_data['gene'].isin(specific_gene)]
    # non_specific_gene_data = full_data[~full_data['gene'].isin(specific_gene)]
    # this way we still have non specific amplicones incorporate the specific genes
    specific_gene_data = full_data[full_data['gene'].isin(specific_gene)]
    non_specific_gene_data = full_data
    ref_size = genome_size(ref_genome)
    specific_gene_data_count = 0
    #main output folder path
    output_path = args.output_folder_path


    # Calculating number of amplicone needed
    # target_coverage = 1
    # gene_coverage = Amplicone_no.place_amplicone_search(full_data, target_coverage, read_size, genome_size(ref_genome))

    covered_positions, covered_ranges = [], []
    primer_pool, no_primer_ = [], []
    accepted_primers = pd.DataFrame(columns = ['pLeft_ID', 'pLeft_coord', 'pLeft_length', 'pLeft_Tm', 'pLeft_GC',
        'pLeft_Sequences', 'pLeft_EndStability', 'pRight_ID', 'pRight_coord',
        'pRight_length', 'pRight_Tm', 'pRight_GC', 'pRight_Sequences',
        'pRight_EndStability', 'Penalty', 'Product_size'])

    if len(specific_gene)>0:
        print('=====Specific amplicon=====')
        covered_positions_sp, covered_ranges_sp, specific_gene_data, primer_pool, accepted_primers, no_primer_ = place_amplicone(specific_gene_data, specific_gene_amplicone, read_size, primer_pool, accepted_primers, no_primer_, genome_size(ref_genome),padding=padding, output_path =output_path)
        specific_gene_data_count = accepted_primers.shape[0]                                                  
        print('=====Non-specific amplicon=====')
        non_specific_gene_data.update(specific_gene_data) # add the specific gene data to the non-specific gene data
        covered_positions_nosp, covered_ranges_nosp, full_data_cp, primer_pool, accepted_primers, no_primer_ = place_amplicone(non_specific_gene_data, non_specific_amplicone, read_size, primer_pool, accepted_primers, no_primer_, genome_size(ref_genome),padding=padding, output_path =output_path)
        covered_positions = {**covered_positions_sp, **covered_positions_nosp}
        covered_ranges = covered_ranges_sp + covered_ranges_nosp
        non_specific_gene_data_count = accepted_primers.shape[0] - specific_gene_data_count
        
    else:
        covered_positions_nosp, covered_ranges_nosp, full_data_cp, primer_pool, accepted_primers, no_primer_ = place_amplicone(non_specific_gene_data, non_specific_amplicone, read_size, primer_pool, accepted_primers, no_primer_,genome_size(ref_genome),padding=padding, output_path =output_path)
        covered_positions = covered_positions_nosp
        covered_ranges = covered_ranges_nosp
        non_specific_gene_data_count = accepted_primers.shape[0]
    
    # whether or not you wanna include spoligotyping sequencing
    spoligotype = args.spoligo_sequencing
    covered_ranges_spol = []
    if spoligotype:
        
        if args.spoligo_sequencing_file == None:
            spacers = pd.read_csv('spacers.bed', sep='\t', header=None)
        else:
            spacers = pd.read_csv(args.spoligo_sequencing_file, sep='\t', header=None)
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
        covered_ranges_spol = Amplicone_no.place_amplicone_spol(spol_data, 1, read_size, graphic_output=False, ref_size = genome_size(ref_genome))
        covered_ranges.extend(covered_ranges_spol)

    read_number = specific_gene_data_count + non_specific_gene_data_count + len(covered_ranges_spol)

    # output
    primer_label = ['Gene_specific']*specific_gene_data_count + ['Non_specific']*non_specific_gene_data_count + ['Spoligotype']*len(covered_ranges_spol)

    accepted_primers['Amplicone_type'] = primer_label
    accepted_primers['Redesign'] = no_primer_
    accepted_primers['Designed_ranges'] = covered_ranges

    op = f'{output_path}/Amplicon_design_output'
    os.makedirs(op, exist_ok=True) #output path
    accepted_primers.to_csv(f'{op}/Primer_design-accepted_primers-{read_number}-{read_size}.csv',index=False)

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

    primer_inclusion.to_csv(f'{op}/SNP_inclusion-{read_number}-{read_size}.csv',index=False)

    print('Primer design output files:')
    print(f'{op}/SNP_inclusion-{read_number}-{read_size}.csv')
    print(f'{op}/Primer_design-accepted_primers-{read_number}-{read_size}.csv')
# %%
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Amplicone design',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # in
    parser.add_argument('-s','--snp_priority', type = str, help = 'SNP priority CSV files (default: collated global 50k clinical TB samples)', default='variants.csv', Default=None)
    parser.add_argument('-a','--amplicon_size', type = int, help = 'Amplicon size', default=400)
    parser.add_argument('-p','--padding_size', type = int, help = 'Size of padding on each side of the target sequence during primer design', default=None)
    parser.add_argument('-ref','--reference_genome', type = str, help = 'reference fasta file (default: MTB-h37rv genome)', default='MTB-h37rv_asm19595v2-eg18.fa')
    parser.add_argument('-sn','--specific_amplicone_no', type = int, help = 'number of amplicone dedicated to amplifying specific genes', default=0)
    parser.add_argument('-sg','--specific_amplicone_gene', type = str, help = 'give a list of gene names separated by Lineage ', default='')
    parser.add_argument('-nn','--non-specific_amplicone_no', type = int, help = 'number of amplicone dedicated to amplifying all SNPs in all genes according the importants list', default=30)
    parser.add_argument('-g','--graphic_option', action='store_true', default = False, help = 'output graphic on amplicone coverage')
    
    parser.add_argument('-sp','--spoligo_sequencing', action='store_true', help = 'Whether to amplify Spoligotype')
    parser.add_argument('-sp_f','--spoligo_sequencing_file', type = str, help = 'Custom spoligotype range files (default: TB spligotype space ranges)', default = None)
    # out
    parser.add_argument('-op','--output_folder_path', default = '', type = str, help = 'output_folder_path (accepted_primers, SNP_inclusion, gene_covered)')
    parser.set_defaults(func=main)
    args = parser.parse_args()
    main(args)
    
    
# running
# python your_script.py -a 400 -ref "MTB-h37rv_asm19595v2-eg18.fa" -sn 0 -nn 30 -g -sp -op "/mnt/storage10/lwang/Projects/Amplicone_design_tool/test_runs"

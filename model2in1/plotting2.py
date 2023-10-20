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
import matplotlib.pyplot as plt
import pandas as pd


#%%
# Sample data
full_data = pd.read_csv('/mnt/storage10/jody/projects/variant_dump/variants.csv')
full_data = full_data.sort_values(by=['genome_pos'])
full_data = full_data.reset_index(drop=True)
full_data['weight'] = full_data['freq']
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

susana_design = pd.read_excel('/mnt/storage10/lwang/Projects/Amplicone_design_tool/model2in1/amplicon_TB_2023_V2.xlsx', header = 0)
susana_design = susana_design.iloc[:-3]
# susana_design = susana_design['primer name'].contains('')
#%%

gene_range = pd.read_csv('regions.bed', sep = '\t', header = 0)
read_size = 400

for g in tb_drug_resistance_genes:
    print(g, '===============')
    gene = g
    g_range = gene_range[gene_range['gene']==gene]
    # full_gene = full_data[full_data['gene']==gene]
    full_gene = full_data
    full_gene = full_gene[~full_gene['drugs'].isna()]
    full_gene = full_gene[~full_gene['type'].isin(['synonymous_variant','non_coding_transcript_exon_variant'])]
    full_gene = full_gene[(full_gene['genome_pos']>=g_range['start'].values[0]) & (full_gene['genome_pos']<=g_range['end'].values[0])]
    # full_gene = full_gene[(full_gene['genome_pos']>=6400) & (full_gene['genome_pos']<=7100)]
    
    full_gene = full_gene[['genome_pos','gene','change','freq']]
    snp_names = full_gene['change'].tolist()
    snp_names = [item for item in snp_names if not item.startswith('c.')]

    full_gene = full_gene[full_gene['change'].isin(snp_names)]

    full_gene['change'] = full_gene['change'].str.extract(r'([a-zA-Z]+\d+)')[0]
    full_gene = full_gene.groupby(['genome_pos', 'gene', 'change']).agg({'freq': 'sum'}).reset_index()

    snp_names = full_gene['change'].tolist()
    snp_frequency = full_gene['freq'].tolist()
    genomic_positions = full_gene['genome_pos'].tolist()
    primers = pd.read_csv('./Amplicon_design_output/Primer_design-accepted_primers-10-1000.csv', index_col=False)
    primers = pd.read_csv('./Amplicon_design_output/Primer_design-accepted_primers-35-400.csv', index_col=False)
    primers.drop(columns=['Unnamed: 0'], inplace=True)

    df = pd.DataFrame({
        'Gene': snp_names,
        'SNP_Frequency': snp_frequency,
        'Genomic_Position': genomic_positions
    })
    if len(df) == 0:
        continue
        

    amplicon_df = pd.DataFrame(columns=primers.columns.tolist())
    for i, row in primers.iterrows():
        in_range = False
        start1 = row['pLeft_coord']
        end1 = row['pRight_coord']
        start2 = df['Genomic_Position'].min()
        end2 = df['Genomic_Position'].max()
        if (start1 <= end2) and (end1 >= start2):
            in_range = True
            # print(row.to_frame())
            amplicon_df = pd.concat([amplicon_df, row.to_frame().T],axis = 0)
    amplicon_positions = []
    amplicon_names = []

    for i, row in amplicon_df.iterrows():
        amplicon_positions.append((row['pLeft_coord'], row['pRight_coord']))
        amplicon_names.append(amplicon_df['pLeft_ID'])


    susana_amplicon_df = pd.DataFrame(columns=susana_design.columns.tolist())
    for i, row in susana_design.iterrows():
        in_range = False
        start1 = row['Start']
        end1 = row['End']
        start2 = df['Genomic_Position'].min()
        end2 = df['Genomic_Position'].max()
        if (start1 <= end2) and (end1 >= start2):
            in_range = True    
            # print(row.to_frame())
            susana_amplicon_df = pd.concat([susana_amplicon_df, row.to_frame().T],axis = 0)
    susana_amplicon_positions = []
    susana_amplicon_names = []

    for i, row in susana_amplicon_df.iterrows():
        susana_amplicon_positions.append((row['Start'], row['End']))
        susana_amplicon_names.append(row['primer name'])

    #getting horizontal bars

    # Create the figure and axis
    fig, ax1 = plt.subplots(figsize=(20, 16), frameon=False, dpi=50)

    # Create the bar chart
    bars = ax1.bar(df['Genomic_Position'], df['SNP_Frequency'], color='b', alpha=0.6)

    # Set y-axis label
    ax1.set_ylabel('SNP Frequency')

    # Set x-axis limits
    ax1.set_xlim([min(df['Genomic_Position'])-10, max(df['Genomic_Position'])+10])
    # ax1.set_xlim([2153898, 2156121])
    # print(min(df['Genomic_Position'])-10, max(df['Genomic_Position'])+10)
    # Set y-axis limits
    ax1.set_ylim([-len(amplicon_positions+susana_amplicon_positions)-10, 100])

    # ax1.set_xticks(df['Genomic_Position'])
    # ax1.set_xticklabels(df['Genomic_Position'])
    ax1.set_xticks(ax1.get_xticks()[::2])  # Only keep every 2nd tick
    ax1.tick_params(axis='x', labelsize=26)
    ax1.tick_params(axis='y', labelsize=26)

    ax1.set_xlabel('Genomic Positions(bps)',fontsize=30)
    ax1.set_ylabel('SNP Frequency',fontsize=30)

    # Add a horizontal line at y=0
    ax1.axhline(0, color='black', linewidth=5, alpha=0.8)

    # Add ticks at y=0 for each bar
    for x in df['Genomic_Position']:
        ax1.plot([x, x], [0.3, -0.3], color='black')

    # # Rotate x-axis tick labels by 90 degrees
    for label in ax1.get_xticklabels():
        label.set_rotation(30)

    # Set title
    ax1.set_title(f'{gene} - Amplicon Coverage {read_size}bps', fontsize=40)

    for i, ((start, end), primer_name) in enumerate(zip(amplicon_positions, amplicon_names)):
        ax1.barh(-i-2, end-start, left=start, color='r', alpha=0.6, label='Designer Amplicons' if i == 0 else "")
        # ax1.axhline(0, color='black', linewidth=5, alpha=0.8)
        # ax1.text((start + end)/2, -5, primer_name, va='bottom', ha='center', color='black')
    d = i
    for i, ((start, end), primer_name) in enumerate(zip(susana_amplicon_positions, susana_amplicon_names)):
        ax1.barh(-i-2-3-d, end-start, left=start, color='g', alpha=0.6, label='Reference Amplicons' if i == 0 else "")

    ax1.legend(fontsize=20, loc='upper right')

    output_path = '.'
    op = f'{output_path}/amplicon_coverage/'
    os.makedirs(op, exist_ok=True) #output path
    plt.savefig(f'{op}/{read_size}bps amplicon_coverage-{gene}.png', bbox_inches='tight')
    plt.show()




# %%

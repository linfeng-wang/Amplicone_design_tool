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

#%%
gene = 'gyrA'
full_gene = full_data[full_data['gene']==gene]
full_gene = full_gene[~full_gene['drugs'].isna()]
full_gene = full_gene[~full_gene['type'].isin(['synonymous_variant','non_coding_transcript_exon_variant'])]

full_gene = full_gene[['genome_pos','gene','change','freq']]
snp_names = full_gene['change'].tolist()
snp_names = [item for item in snp_names if not item.startswith('c.')]

full_gene = full_gene[full_gene['change'].isin(snp_names)]

full_gene['change'] = full_gene['change'].str.extract(r'([a-zA-Z]+\d+)')[0]
full_gene = full_gene.groupby(['genome_pos', 'gene', 'change']).agg({'freq': 'sum'}).reset_index()

snp_names = full_gene['change'].tolist()
snp_frequency = full_gene['freq'].tolist()
genomic_positions = full_gene['genome_pos'].tolist()
primers = pd.read_csv('/mnt/storage10/lwang/Projects/Amplicone_design_tool/model/Amplicone_design_output/Primer_design-accepted_primers-30-1000.csv', index_col=False)
primers.drop(columns=['Unnamed: 0'], inplace=True)

df = pd.DataFrame({
    'Gene': snp_names,
    'SNP_Frequency': snp_frequency,
    'Genomic_Position': genomic_positions
})

amplicon_df = pd.DataFrame(columns=primers.columns.tolist())
for i, row in primers.iterrows():
    in_range = False
    if row['pLeft_coord'] >= df['Genomic_Position'].min() and row['pLeft_coord'] <= df['Genomic_Position'].max():
        in_range = True
    elif row['pRight_coord'] >= df['Genomic_Position'].min() and row['pRight_coord'] <= df['Genomic_Position'].max():
        in_range = True
    if in_range == True:
        # print(row.to_frame())
        amplicon_df = pd.concat([amplicon_df, row.to_frame().T],axis = 0)
amplicon_positions = []
amplicon_names = []

for i, row in amplicon_df.iterrows():
    amplicon_positions.append((row['pLeft_coord'], row['pRight_coord']))
    amplicon_names.append(amplicon_df['pLeft_ID'])

df = pd.read_csv('filenamedata.csv')
df = df[df['SNP_Frequency']> 0]

#getting horizontal bars

# Create the figure and axis
fig, ax1 = plt.subplots(figsize=(20, 16), frameon=False, dpi=50)

# Create the bar chart
bars = ax1.bar(df['Genomic_Position'], df['SNP_Frequency'], color='b', alpha=0.6, width=10)

# Set y-axis label
ax1.set_ylabel('SNP Frequency')

# Set x-axis limits
ax1.set_xlim([min(df['Genomic_Position'])-10, max(df['Genomic_Position'])+10])
# print(min(df['Genomic_Position'])-10, max(df['Genomic_Position'])+10)
# Set y-axis limits
ax1.set_ylim([-10, 100])

# ax1.set_xticks(df['Genomic_Position'])
# ax1.set_xticklabels(df['Genomic_Position'])
ax1.set_xticks(ax1.get_xticks()[::2])  # Only keep every 2nd tick
ax1.tick_params(axis='x', labelsize=16)
ax1.tick_params(axis='y', labelsize=16)

ax1.set_xlabel('Genomic Positions(bps)',fontsize=20)
ax1.set_ylabel('SNP Frequency',fontsize=20)

# Add a horizontal line at y=0
ax1.axhline(0, color='black', linewidth=5, alpha=0.8)

# Add ticks at y=0 for each bar
for x in df['Genomic_Position']:
    ax1.plot([x, x], [0.3, -0.3], color='black')

# # Rotate x-axis tick labels by 90 degrees
for label in ax1.get_xticklabels():
    label.set_rotation(30)

# Set title
ax1.set_title(f'{gene} - Amplicon Cover', fontsize=30)

for i, ((start, end), primer_name) in enumerate(zip(amplicon_positions, amplicon_names)):
    ax1.barh(-i-2, end-start, left=start, color='r', alpha=0.8)
    # ax1.axhline(0, color='black', linewidth=5, alpha=0.8)

    # ax1.text((start + end)/2, -5, primer_name, va='bottom', ha='center', color='black')

plt.show()



# %%

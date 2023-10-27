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
import matplotlib.pyplot as plt
import pandas as pd

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
gene = 'katG'
full_gene = full_data[full_data['gene']==gene]
# full_gene['change'] = full_gene['change'].str.extract(r'(\d+)')[0]
# full_gene = full_gene[full_gene['weight']>1.3]
full_gene = full_gene[~full_gene['drugs'].isna()]
full_gene = full_gene[~full_gene['type'].isin(['synonymous_variant','non_coding_transcript_exon_variant'])]
# full_gene = full_gene.drop_duplicates(subset='change',)
# full_gene = full_gene.drop_duplicates(subset='genome_pos')
full_gene = full_gene[['genome_pos','gene','change','freq']]
snp_names = full_gene['change'].tolist()
snp_names = [item for item in snp_names if not item.startswith('c.')]

full_gene = full_gene[full_gene['change'].isin(snp_names)]

full_gene['change'] = full_gene['change'].str.extract(r'([a-zA-Z]+\d+)')[0]
full_gene = full_gene.groupby(['genome_pos', 'gene', 'change']).agg({'freq': 'sum'}).reset_index()
#%%
snp_names = full_gene['change'].tolist()
snp_frequency = full_gene['freq'].tolist()
genomic_positions = full_gene['genome_pos'].tolist()
primers = pd.read_csv('/mnt/storage10/lwang/Projects/Amplicone_design_tool/model/Amplicone_design_output/Primer_design-accepted_primers-30-1000.csv', index_col=False)
primers.drop(columns=['Unnamed: 0'], inplace=True)

amplicon_df = pd.DataFrame(columns=primers.columns.tolist())
for i, row in primers.iterrows():
    in_range = False
    if row['pLeft_coord'] > full_gene['genome_pos'].min() and row['pLeft_coord'] < full_gene['genome_pos'].max():
        in_range = True
    if row['pRight_coord'] < full_gene['genome_pos'].min() and row['pRight_coord'] < full_gene['genome_pos'].max():
        in_range = True
    if in_range == True:
        # print(row.to_frame())
        amplicon_df = pd.concat([amplicon_df, row.to_frame().T],axis = 0)
amplicon_positions = []
amplicon_names = []

for i, row in amplicon_df.iterrows():
    amplicon_positions.append((row['pLeft_coord'], row['pRight_coord']))
    amplicon_names.append(amplicon_df['pLeft_ID'])
print(len(amplicon_positions))
print(len(amplicon_names))
print(len(snp_names))

#%%
# Create DataFrame
df = pd.DataFrame({
    'Gene': snp_names,
    'SNP_Frequency': snp_frequency,
    'Genomic_Position': genomic_positions
})

#%%
# Create the figure and axis
fig, ax1 = plt.subplots(figsize=(100, 60))

# Create the bar chart
bars = ax1.bar(df['Genomic_Position'], df['SNP_Frequency'], color='b', alpha=0.6)

# Add gene names on top of bars
for bar, gene_name in zip(bars, df['Gene']):
    yval = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2.0, yval, gene_name, va='bottom', ha='center')

# Set x-axis to represent genomic positions
ax1.set_xticks(df['Genomic_Position'])
ax1.set_xticklabels(df['Genomic_Position'])
ax1.set_xlabel('Genomic Positions')
ax1.set_ylabel('SNP Frequency')
ax1.set_xlim([min(df['Genomic_Position'])-10, max(df['Genomic_Position'])+10])
# Set y-axis limits
ax1.set_ylim([10, 100])
plt.title('Gene SNP Frequency and Primer Positions')

# Add a horizontal line at y=0
ax1.axhline(0, color='black', linewidth=0.8)

# Overlay horizontal bars for primer positions
ax2 = ax1.twiny()  # Create new x-axis, sharing y-axis
ax2.set_xlim(ax1.get_xlim())  # Same x-axis limits
ax2.set_xticks([])  # No tick labels for the new x-axis

# # Add primers as horizontal bars and label them
# for (start, end), primer_name in zip(amplicon_positions, amplicon_names):
#     ax2.barh(-1, end-start, left=start, height=0.5, color='r', alpha=0.5)
#     ax2.text((start + end)/2, -1, primer_name, va='bottom', ha='center', color='black')

plt.show()
#%%
df.to_csv('filenamedata.csv', index=False)

#%%
df = pd.read_csv('filenamedata.csv')
df = df[df['SNP_Frequency']> 10]
# Create the figure and axis
fig, ax1 = plt.subplots(figsize=(15, 6), frameon=False, dpi=100)

# Create the bar chart
bars = ax1.bar(df['Genomic_Position'], df['SNP_Frequency'], color='b', alpha=0.6, width=5)

# Set y-axis label
ax1.set_ylabel('SNP Frequency')

# Set x-axis limits
ax1.set_xlim([min(df['Genomic_Position'])-10, max(df['Genomic_Position'])+10])

# Set y-axis limits
ax1.set_ylim([10, 100])

ax1.set_xticks(df['Genomic_Position'])
ax1.set_xticklabels(df['Genomic_Position'])
ax1.set_xlabel('Genomic Positions')
ax1.set_ylabel('SNP Frequency')
plt.title('Gene SNP Frequency and Primer Positions')

# Add a horizontal line at y=0
ax1.axhline(0, color='black', linewidth=0.8)

# Overlay horizontal bars for primer positions
ax2 = ax1.twiny()  # Create new x-axis, sharing y-axis
ax2.set_xlim(ax1.get_xlim())  # Same x-axis limits
ax2.set_xticks([])  # No tick labels for the new x-axis

# Add primers as horizontal bars and label them
for (start, end), primer_name in zip(amplicon_positions, amplicon_names):
    ax2.barh(-1, end-start, left=start, height=0.5, color='r', alpha=0.5)
    ax2.text((start + end)/2, -1, primer_name, va='bottom', ha='center', color='black')

plt.show()
# # Rotate x-axis tick labels by 90 degrees
# for label in ax1.get_xticklabels():
#     label.set_rotation(90)

# Set title
plt.title('Gene SNP Frequency and Primer Positions')

# Show the plot
plt.show()










# %%
import matplotlib.pyplot as plt
import pandas as pd

# Sample data
# snp_names = ['gene1', 'gene2', 'gene3', 'gene4']
# snp_frequency = [5, 10, 7, 3]
# genomic_positions = [10, 20, 30, 40]  # Replace with actual genomic positions
# amplicon_positions = [(8, 25), (18, 22), (28, 32)]
# amplicon_names = ['primer1', 'primer2', 'primer3']

# Create DataFrame
df = pd.DataFrame({
    'Gene': snp_names,
    'SNP_Frequency': snp_frequency,
    'Genomic_Position': genomic_positions
})

# Create the figure and axis
fig, ax1 = plt.subplots(figsize=(10, 6),frameon=False, dpi=100)

# Create the bar chart
bars = ax1.bar(df['Genomic_Position'], df['SNP_Frequency'], color='b', alpha=0.6)

# Add gene names on top of bars
for bar, gene_name in zip(bars, df['Gene']):
    yval = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2.0, yval, gene_name, va='bottom', ha='center')

# Set x-axis to represent genomic positions
ax1.set_xticks(df['Genomic_Position'])
ax1.set_xticklabels(df['Genomic_Position'])
# ax1.set_xscale('log')
# ax1.set_xlabel('Genomic Positions (log scale)')
ax1.set_ylabel('SNP Frequency')
plt.title('Gene SNP Frequency and Primer Positions')

# Add a horizontal line at y=0
ax1.axhline(0, color='black', linewidth=2)

# Add ticks at y=0 for each bar
for x in df['Genomic_Position']:
    ax1.plot([x, x], [0.3, -0.3], color='black')

# Overlay horizontal bars for primer positions
ax2 = ax1.twiny()  # Create new x-axis, sharing y-axis
ax2.set_xlim(ax1.get_xlim())  # Same x-axis limits
ax2.set_xticks([])  # No tick labels for the new x-axis

# Add primers as horizontal bars and label them
for (start, end), primer_name in zip(amplicon_positions, amplicon_names):
    ax2.barh(-1, end-start, left=start, height=0.5, color='r', alpha=0.5)
    ax2.text((start + end)/2, -1, primer_name, va='bottom', ha='center', color='black')

plt.savefig('filename.png')

# %%
# Create DataFrame
df = pd.DataFrame({
    'Gene': gene_names,
    'SNP_Frequency': snp_frequency,
    'Genomic_Position': genomic_positions
})

# Create the figure and axis
fig, ax1 = plt.subplots(figsize=(10, 6))

# Create the bar chart
bars = ax1.bar(df['Genomic_Position'], df['SNP_Frequency'], color='b', alpha=0.6)

# Add gene names on top of bars
for bar, gene_name in zip(bars, df['Gene']):
    yval = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2.0, yval, gene_name, va='bottom', ha='center')

# Set x-axis to represent genomic positions with exponential scale
ax1.set_xscale('log')
ax1.set_xlabel('Genomic Positions (log scale)')
ax1.set_ylabel('SNP Frequency')
plt.title('Gene SNP Frequency and Primer Positions')

# Add a horizontal line at y=0
ax1.axhline(0, color='black', linewidth=0.8)

# Add ticks at y=0 for each bar
for x in df['Genomic_Position']:
    ax1.plot([x, x], [0, -0.5], color='black')

# Overlay horizontal bars for primer positions
ax2 = ax1.twiny()  # Create new x-axis, sharing y-axis
ax2.set_xlim(ax1.get_xlim())  # Same x-axis limits
ax2.set_xticks([])  # No tick labels for the new x-axis

# Add primers as horizontal bars and label them
for (start, end), primer_name in zip(primer_positions, primer_names):
    ax2.barh(-1, end-start, left=start, height=0.5, color='r', alpha=0.5)
    ax2.text((start + end)/2, -1, primer_name, va='bottom', ha='center', color='black')

plt.savefig('filename.png')
# %%
#%%
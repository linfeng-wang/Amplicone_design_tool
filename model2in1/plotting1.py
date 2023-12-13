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
import Amplicon_no
reload(Amplicon_no)
import argparse
from functools import reduce
import os
import matplotlib.pyplot as plt
import pandas as pd


#%%
# Sample data
def plotting(priority, read_size, accepted_primers, gff, reference_design, output_dir):
    full_data = pd.read_csv(priority)
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
    if reference_design is not None:
        reference_design = pd.read_csv(reference_design, header = 0)
    # reference_design = reference_design.iloc[:-3]
    # reference_design = reference_design['primer name'].contains('')

    #%%
    # ideal_ranges = [[761004, 761404],
    # [2154961, 2155361],
    # [4247393, 4247793],
    # [2288820, 2289220],
    # [781585, 781985],
    # [4247186, 4247586],
    # [1673373, 1673773],
    # [7247, 7647],
    # [2288528, 2288928],
    # [1472358, 1472758],
    # [4407756, 4408156],
    # [1473003, 1473403],
    # [2715339, 2715739],
    # [2715339, 2715739],
    # [1472035, 1472435],
    # [4243203, 4243603],
    # [4243203, 4243603],
    # [4247730, 4248130],
    # [4326074, 4326474],
    # [1674048, 1674448],
    # [4327084, 4327484],
    # [4326602, 4327002],
    # [1674048, 1674448],
    # [2726112, 2726512],
    # [764724, 765124],
    # [6575, 6975],
    # [4249389, 4249789],
    # [1674481, 1674881],
    # [2155462, 2155862],
    # [2747141, 2747541],
    # [4407543, 4407943],
    # [4327484, 4327884],
    # [3067958, 3068358],
    # [3067958, 3068358],
    # [1918292, 1918692],
    # [6575, 6975],
    # [2726112, 2726512],
    # [4325761, 4326161],
    # [3841083, 3841483],
    # [762089, 762489],
    # [762089, 762489],
    # [3841083, 3841483],
    # [4326265, 4326665]]

    #%%
    genes = pd.read_csv(gff, sep = '\t', header = None, skiprows=[0,1,2,3,4,5,6])
    genes = genes[genes[2]=='gene']
    genes['Name'] = genes[8].str.extract(r'Name=([^;]*)')
    genes = genes[[3,4,'Name']]
    genes.columns = ['start', 'end', 'gene']
    # print(genes)
    #%%
    gene_range = pd.read_csv('regions.bed', sep = '\t', header = 0)
        
    # primers = pd.read_csv('./Amplicon_design_output/Primer_design-accepted_primers-10-1000.csv', index_col=False)
    primers = pd.read_csv(accepted_primers, index_col=False)
    # primers.drop(columns=['Unnamed: 0'], inplace=True)

    for g in tb_drug_resistance_genes:
        # print(g, '===============')
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

        full_gene['freq'] = full_gene['freq'].apply(lambda x: max(x, 1))

        # print(full_gene)
        snp_names = full_gene['change'].tolist()
        snp_frequency = full_gene['freq'].tolist()
        genomic_positions = full_gene['genome_pos'].tolist()


        df = pd.DataFrame({
            'Gene': snp_names,
            'SNP_Frequency': np.log10(snp_frequency),
            'Genomic_Position': genomic_positions
        })

        if len(df) == 0:
            print(f'---> No SNPs found in: {gene}')
            continue
        else:
            # print(f'---> SNPs found in: {gene}')
            pass
            
    # getting gene ranges   
        genes_df = pd.DataFrame(columns=genes.columns.tolist())
        for i, row in genes.iterrows():
            in_range = False
            start1 = row['start']
            end1 = row['end']
            start2 = df['Genomic_Position'].min()
            end2 = df['Genomic_Position'].max()
            if (start1 <= end2) and (end1 >= start2):
                in_range = True    
                # print(row.to_frame())
                genes_df = pd.concat([genes_df, row.to_frame().T],axis = 0)
        genes_positions = []
        genes_names = []
        # print(genes_positions)
        for i, row in genes_df.iterrows():
            genes_positions.append((row['start'], row['end']))
            genes_names.append(row['gene'])

    # getting designed amplicon ranges

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
            amplicon_names.append(row['pLeft_ID'])
        
    # getting reference amplicon ranges
        if reference_design is not None:
            reference_amplicon_df = pd.DataFrame(columns=reference_design.columns.tolist())
            for i, row in reference_design.iterrows():
                in_range = False
                start1 = row['Start']
                end1 = row['End']
                start2 = df['Genomic_Position'].min()
                end2 = df['Genomic_Position'].max()
                if (start1 <= end2) and (end1 >= start2):
                    in_range = True    
                    # print(row.to_frame())
                    reference_amplicon_df = pd.concat([reference_amplicon_df, row.to_frame().T],axis = 0)
            reference_amplicon_positions = []
            reference_amplicon_names = []

            for i, row in reference_amplicon_df.iterrows():
                reference_amplicon_positions.append((row['Start'], row['End']))
                reference_amplicon_names.append(row['primer name'])

    # getting designed amplicon ranges

        # ideal_amplicons = []
        # for row in ideal_ranges:
        #     in_range = False
        #     start1 = row[0]
        #     end1 = row[1]
        #     start2 = df['Genomic_Position'].min()
        #     end2 = df['Genomic_Position'].max()
        #     if (start1 <= end2) and (end1 >= start2):
        #         in_range = True    
        #         # print(row.to_frame())
        #         ideal_amplicons.append(row)

        #getting horizontal bars

        # Create the figure and axis
        fig, ax1 = plt.subplots(figsize=(20, 16), frameon=True, dpi=50)

        # Create the bar chart
        bars = ax1.bar(df['Genomic_Position'], df['SNP_Frequency'], color='b', alpha=0.6)

        # Set y-axis label
        ax1.set_ylabel('SNP Frequency')

        # Set x-axis limits
        ax1.set_xlim([min(df['Genomic_Position'])-10, max(df['Genomic_Position'])+10])
        # ax1.set_xlim([2153898, 2156121])
        # print(min(df['Genomic_Position'])-10, max(df['Genomic_Position'])+10)
        # Set y-axis limits
        if reference_design is not None:
            # ax1.set_ylim([-len(amplicon_positions+reference_amplicon_positions)-10, None])
            ax1.set_ylim(ymin=-len(amplicon_positions+reference_amplicon_positions)-5, ymax=len(genes_names)+3)
            # ax1.set_ylim(ymin=-len(amplicon_positions+reference_amplicon_positions+ideal_amplicons)-8, ymax=None)


        ax1.set_xticks(ax1.get_xticks()[::2])  # Only keep every 2nd tick
        ax1.tick_params(axis='x', labelsize=26)
        ax1.tick_params(axis='y', labelsize=26)

        ax1.set_xlabel('Genomic Positions(bps)',fontsize=30)
        ax1.set_ylabel('log10(SNP Frequency)',fontsize=30)

        # Add a horizontal line at y=0
        ax1.axhline(0, color='black', linewidth=5, alpha=0.8)

        # Set title
        ax1.set_title(f'{gene} - Amplicon Coverage {read_size}bps', fontsize=40)
        
        for i, ((start, end), primer_name) in enumerate(zip(genes_positions, genes_names)):
            ax1.barh(i+0.7, end-start, left=start, color='grey', alpha=0.3, label='Gene ranges' if i == 0 else "")
            # ax1.axhline(0, color='black', linewidth=5, alpha=0.8)
            if min(df['Genomic_Position'])-10 >= (start + end)/2:
                ax1.text(min(df['Genomic_Position'])+100, i+0.7, primer_name, va='bottom', ha='center', color='black', fontsize=25)
            elif (start + end)/2 >= max(df['Genomic_Position'])+10:
                ax1.text(max(df['Genomic_Position'])-100, i+0.7, primer_name, va='bottom', ha='center', color='black', fontsize=25)
            else:
                ax1.text((start + end)/2, i+0.7, primer_name, va='bottom', ha='center', color='black', fontsize=25)
        d = i

        for i, ((start, end), primer_name) in enumerate(zip(amplicon_positions, amplicon_names)):
            ax1.barh(-i-2, end-start, left=start, color='r', alpha=0.4, label='Designer Amplicons' if i == 0 else "")
            # ax1.axhline(0, color='black', linewidth=5, alpha=0.8)
            primer_name = primer_name.split('-')
            primer_name = primer_name[:-1]
            primer_name = '-'.join(primer_name)
            ax1.text((start + end)/2, -i-2.2, primer_name, va='bottom', ha='center', color='black', fontsize=25)
        if reference_design is not None:
            d = i
            for i, ((start, end), primer_name) in enumerate(zip(reference_amplicon_positions, reference_amplicon_names)):
                ax1.barh(-i-2-3-d, end-start, left=start, color='g', alpha=0.6, label='Reference Amplicons' if i == 0 else "")
        # print(reference_amplicon_positions)
        # e = i
        # for i, (start, end) in enumerate(ideal_amplicons):
        #     ax1.barh(-i-2-6-d-e, end-start, left=start, color='b', alpha=0.6, label='Ideal Amplicons' if i == 0 else "")
        ax1.legend(fontsize=20, loc='lower left')

        output_path = f'{output_dir}'
        op = f'{output_path}/Amplicon_coverage_plots/'
        os.makedirs(op, exist_ok=True) #output path
        fig.savefig(f'{op}/{read_size}bps-amplicon_coverage-{gene}.png', bbox_inches='tight')
    print(f'**Graphic output saved to: {op}')
        # plt.show()

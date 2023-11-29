from functools import partial
from random import choices, randint, randrange, random, sample
from typing import List, Optional, Callable, Tuple
import numpy as np
# from geneticalgorithm import geneticalgorithm as ga
from rich_argparse import ArgumentDefaultsRichHelpFormatter
import pandas as pd
from collections import Counter
from tqdm import tqdm
import time
from Bio.SeqUtils import MeltingTemp
from Bio import SeqIO
from plotly import graph_objects as go
import json
from importlib import reload
import primer_selection
reload(primer_selection)
import Amplicone_no
reload(Amplicone_no)
import plotting
reload(plotting)
import working_algo_gene_2in1 as w
reload(w)
import argparse
from functools import reduce
import os
import matplotlib.pyplot as plt
import pandas as pd

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

    gene_names = [
        "rpoB",
        "katG",
        "embB",
        "pncA",
        "rpsL",
        "rrs",
        "ethA",
        "fabG1",
        "gyrA",
        "gid",
        "inhA",
        "ethR",
        "rpoC",
        "ahpC",
        "gyrB",
        "folC",
        "tlyA",
        "alr",
        "embA",
        "thyA",
        "eis"
    ]
    # full_data = w.load_data(args.snp_priority, gene_names)
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
    ref_size = w.genome_size(ref_genome)
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
        covered_positions_sp, covered_ranges_sp, specific_gene_data, primer_pool, accepted_primers, no_primer_ = w.place_amplicone(specific_gene_data, specific_gene_amplicone, read_size, primer_pool, accepted_primers, no_primer_, w.genome_size(ref_genome),padding=padding, output_path =output_path)
        specific_gene_data_count = accepted_primers.shape[0]                                                  
        print('=====Non-specific amplicon=====')
        non_specific_gene_data.update(specific_gene_data) # add the specific gene data to the non-specific gene data
        covered_positions_nosp, covered_ranges_nosp, full_data_cp, primer_pool, accepted_primers, no_primer_ = w.place_amplicone(non_specific_gene_data, non_specific_amplicone, read_size, primer_pool, accepted_primers, no_primer_, w.genome_size(ref_genome),padding=padding, output_path =output_path)
        covered_positions = {**covered_positions_sp, **covered_positions_nosp}
        covered_ranges = covered_ranges_sp + covered_ranges_nosp
        non_specific_gene_data_count = accepted_primers.shape[0] - specific_gene_data_count
        
    else:
        covered_positions_nosp, covered_ranges_nosp, full_data_cp, primer_pool, accepted_primers, no_primer_ = w.place_amplicone(non_specific_gene_data, non_specific_amplicone, read_size, primer_pool, accepted_primers, no_primer_,w.genome_size(ref_genome),padding=padding, output_path =output_path)
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
        covered_ranges_spol = Amplicone_no.place_amplicone_spol(spol_data, 1, read_size, graphic_output=False, ref_size = w.genome_size(ref_genome))
        covered_ranges.extend(covered_ranges_spol)

    read_number = specific_gene_data_count + non_specific_gene_data_count + len(covered_ranges_spol)

    # output
    primer_label = ['Gene_specific']*specific_gene_data_count + ['Non_specific']*non_specific_gene_data_count + ['Spoligotype']*len(covered_ranges_spol)

    accepted_primers['Amplicone_type'] = primer_label
    accepted_primers['Redesign'] = no_primer_
    accepted_primers['Designed_ranges'] = covered_ranges

    #accepted_primers - change primer to iupac
    threshold = 0.001
    all_snps = pd.read_csv(args.all_snps, sep = '\t', header = None)
    all_snps.drop_duplicates(inplace=True)
    for i, row in accepted_primers.iterrows():
        #left primer
        primer_seq = ''
        for x,y in zip(range(row['pLeft_coord']+len(row['pLeft_Sequences'])), row['pLeft_Sequences']):
            if (all_snps[all_snps[0] == x][[1,2]].shape[0] > 0) and (all_snps[all_snps[0] == x][[3]].values.astype('float').item()>= threshold):
                alleles = ''.join(all_snps[all_snps[0] == x][[1,2]].values[0])
                primer_seq+w.nucleotide_to_iupac(alleles)
            else:
                primer_seq+y
        if row['pLeft_Sequences'] != primer_seq:
            print('SNP in Left primer')
            accepted_primers.loc[i, 'pLeft_Sequences'] = primer_seq
        #right primer
        primer_seq = ''
        for x,y in zip(range(row['pRight_coord']+len(row['pRight_Sequences'])), primer_selection.reverse_complement_sequence(row['pRight_Sequences'])):
            if all_snps[all_snps[0] == x][[1,2]].shape[0] > 0 and all_snps[all_snps[0] == x][[3]].values.astype('float').item()>= threshold:
                alleles = ''.join(all_snps[all_snps[0] == x][[1,2]].values[0])
                primer_seq+w.nucleotide_to_iupac(alleles)
            else:
                primer_seq+y
        if row['pRight_Sequences'] != primer_seq:
            print('SNP in Left primer')
            accepted_primers.loc[i, 'pRight_Sequences'] = primer_seq


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
    
    
def main_amplicon_no(args):
    full_data = pd.read_csv(args.snp_priority)
    full_data = full_data[~full_data['drugs'].isna()]
    full_data = full_data.sort_values(by=['genome_pos'])
    full_data = full_data.reset_index(drop=True)
    full_data['weight'] = full_data['freq']
    ref_genome = args.reference_genome
    target_coverage = args.target_coverage
    read_size = args.amplicon_size
    gene_coverage = Amplicone_no.place_amplicone_search(full_data, target_coverage, read_size, w.genome_size(ref_genome))
    print(gene_coverage)
    
    return 0
    
def main_plotting(args):
    priority = pd.read_csv(args.snp_priority)
    accepted_primers = args.accepted_primers
    gff = args.gff
    reference_design = args.reference_design
    output_dir = args.output_dir
    plotting(priority, accepted_primers, gff, reference_design, output_dir)
    return 0
    
# %%
if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='Amplicon designer for TB', 
                                    description='Amplicone design, list of specific genes that can be priotised: rpoB,katG,embB,pncA,rpsL,rrs,ethA,fabG1,gyrA,gid,inhA,ethR,rpoC,ahpC,gyrB,folC,tlyA,alr,embA,thyA,eis',
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    subparsers = parser.add_subparsers(dest="command", help="Task to perform")


    ###### Design pipeline
    parser_sub = subparsers.add_parser('design', help='Run whole design pipeline', formatter_class=ArgumentDefaultsRichHelpFormatter)
    input=parser_sub.add_argument_group("Input options")
    # parser.add_argument('-c','--country_file', type = str, help = 'SNP priority CSV files (default: collated global 50k clinical TB samples)', default='variants.csv', Default=None)
    # in
    input.add_argument('-s','--snp_priority', type = str, help = 'SNP priority CSV files (default: collated global 50k clinical TB samples)', default='variants.csv', Default='/mnt/storage10/lwang/Projects/Amplicone_design_tool/model2in1/variants.csv')
    input.add_argument('-ref','--reference_genome', type = str, help = 'reference fasta file (default: MTB-h37rv genome)', default='MTB-h37rv_asm19595v2-eg18.fa')
    input.add_argument('-sp_f','--spoligo_sequencing_file', type = str, help = 'Custom spoligotype range files (default: TB spligotype space ranges)', default = None)
    input.add_argument('-as','--all_snps', type = str, help = 'All SNPs in the reference genome', default = 'all_snps.txt')
    
    setting=parser_sub.add_argument_group("Design options")
    setting.add_argument('-a','--amplicon_size', type = int, help = 'Amplicon size', default=400)
    setting.add_argument('-p','--padding_size', type = int, help = 'Size of padding on each side of the target sequence during primer design', default=200)
    setting.add_argument('-sn','--specific_amplicone_no', type = int, help = 'number of amplicone dedicated to amplifying specific genes', default=0)
    setting.add_argument('-sg','--specific_amplicone_gene', type = str, help = 'give a list of gene names separated by Lineage ', default='')
    setting.add_argument('-nn','--non-specific_amplicone_no', type = int, help = 'number of amplicone dedicated to amplifying all SNPs in all genes according the importants list', default=30)
    setting.add_argument('-g','--graphic_option', action='store_true', default = False, help = 'output graphic on amplicone coverage')
    setting.add_argument('-sp','--spoligo_sequencing', action='store_true', help = 'Whether to amplify Spoligotype', default = False)

    # out
    output=parser_sub.add_argument_group("Output options")
    output.add_argument('-op','--output_folder_path', default = '', type = str, help = 'output_folder_path (accepted_primers, SNP_inclusion, gene_covered)', required=True)
    parser_sub.set_defaults(func=main)
    
    ###### Amplicon number estimates
    parser_sub = subparsers.add_parser('amplicon_no', help='Amplicon number estimates', formatter_class=ArgumentDefaultsRichHelpFormatter)
    input=parser_sub.add_argument_group("input options")
    input.add_argument('-s','--snp_priority', type = str, help = 'SNP priority CSV files (default: collated global 50k clinical TB samples)', default='variants.csv')
    input.add_argument('-op','--spliogo range', type = str, help = 'Genomic ranges for spligotype region', default='spacers.bed')
    input.add_argument('-ref','--reference_genome', type = str, help = 'reference fasta file (default: MTB-h37rv genome)', default='MTB-h37rv_asm19595v2-eg18.fa')
    
    setting=parser_sub.add_argument_group("Amplicon options")
    setting.add_argument('-a','--amplicon_size', type = int, help = 'Amplicon size', default=400)
    setting.add_argument('-c','--target_coverage', type = int, help = 'target coverage of SNPs default: full coverage', default=1)
    setting.add_argument('-s','--spoligotype', action='store_true', help = 'Whether to do a separate run on chekcing the number of amplicone needed to cover the spligotypes', default= False)
    parser_sub.set_defaults(func=main_amplicon_no)

    # output=parser_sub.add_argument_group("Output options")
    # output.add_argument('-op','--output_folder_path', default = './', type = str, help = 'output folder path for covered ranges')
    # output.add_argument('-op','--output_folder_path', default = '', type = str, help = 'output folder path for covered ranges')

    ###### Visualisation of designed amplicones
    parser_sub = subparsers.add_parser('plotting', help='Amplicon number estimates', formatter_class=ArgumentDefaultsRichHelpFormatter)
    input=parser_sub.add_argument_group("input options")
    input.add_argument('-s','--snp_priority', type = str, help = 'SNP priority CSV files (default: collated global 50k clinical TB samples)', default='variants.csv')
    input.add_argument('-gff','--gff_features', type = str, help = 'genomic feature file .gff for the corresponding genome', default='MTB-h37rv_asm19595v2-eg18.gff', required=True)
    input.add_argument('-ap','--accepted_primers', type = str, help = 'primer design output file from desgin function', required=True)
    input.add_argument('-rp','--reference_design', type = str, help = '(reference) design that can be plotted against the designed amplicons for comparision', default=None)
    output=parser_sub.add_argument_group("output options")
    output.add_argument('-op','--output_folder_path', default = '', type = str, help = 'output_folder_path (accepted_primers, SNP_inclusion, gene_covered)', required=True)
    parser_sub.set_defaults(func=main_plotting)

    args = parser.parse_args()
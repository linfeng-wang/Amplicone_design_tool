#%%
from functools import partial
from random import choices, randint, randrange, random, sample
from typing import List, Optional, Callable, Tuple
import numpy as np
# from geneticalgorithm import geneticalgorithm as ga
import pandas as pd
from collections import Counter
from tqdm import tqdm
#%%
# Genome = List[int] # 0 or 1
# Population = List[Genome] # list of genomes
# PopulateFunc = Callable[[], Population] # function that returns a population
# FitnessFunc = Callable[[Genome], int] # function that returns a fitness score for a genome
# SelectionFunc = Callable[[Population, FitnessFunc], Tuple[Genome, Genome]] # function that returns a pair of genomes
# CrossoverFunc = Callable[[Genome, Genome], Tuple[Genome, Genome]] # function that returns a pair of genomes
# MutationFunc = Callable[[Genome], Genome] # function that returns a genome
# PrinterFunc = Callable[[Population, int, FitnessFunc], None] # function that prints stats

#%%
# raw_data = pd.read_csv("https://raw.githubusercontent.com/GaryNapier/tb-lineages/main/fst_results_clean_fst_1_for_paper.csv")
# raw_data['weight'] = sample(range(1, 10000), raw_data.shape[0])
# raw_data['weight'] = raw_data['weight'] / 100
# raw_data = raw_data.sort_values(by=['Pos'])
# raw_data = raw_data.reset_index(drop=True)
full_data = pd.read_csv('/mnt/storage10/jody/projects/variant_dump/variants.csv')
full_data = full_data[~full_data['drugs'].isna()]

#%%
snp_weight = full_data['sample_id'].value_counts()/full_data['sample_id'].value_counts().max()
snp_weight = snp_weight.to_frame()
#%%
full_data['weight'] = 0
merged_data = full_data.merge(snp_weight, left_on='sample_id', right_index=True, how='left')
merged_data = merged_data.sort_values(by=['genome_pos'])
merged_data = merged_data.reset_index(drop=True)

# merged_data['weight'].fillna(value=0, inplace=True)
#%%
def rolling_sum(lst, window_size):
    """
    Calculates the rolling sum of a list with a given window size.

    Parameters:
        lst (list): The list to calculate the rolling sum for.
        window_size (int): The size of the rolling window.

    Returns:
        list: A list containing the rolling sum values.
    """
    # Calculate the rolling sum using a list comprehension
    rolling_sum = [sum(lst[i:i+window_size]) for i in range(len(lst)-window_size+1)]
    return rolling_sum
#%%
#Full data trial
read_number = 40
read_size = 400
# priorities = []

weight_window_sum = rolling_sum(merged_data['weight_y'].tolist(), read_size)

covered_positions = {}
covered_ranges = []
for x in tqdm(range(read_number)):
    if len(covered_ranges) == 0:
        start = np.argmax(weight_window_sum)
        end = start + read_size
    if len(covered_ranges) > 0:
        start = np.argmax(weight_window_sum)
        end = start + read_size
        for i in range(len(covered_ranges)):
            if start > covered_ranges[i][0] and start < covered_ranges[i][1]:
                start = covered_ranges[i][1]+1
                if covered_ranges[i][1]+1 + read_size > merged_data.shape[0]:
                    end = merged_data.shape[0]
                end = covered_ranges[i][1]+1 + read_size
            if end > covered_ranges[i][0] and end < covered_ranges[i][1]:
                end = covered_ranges[i][0]-1
                if covered_ranges[i][0]-1 - read_size < 0:
                    start = 0
                start = covered_ranges[i][0]-1 - read_size

    # covered_positions[f'Amplicon_{x+1}'] = {'Range':{'Start': start, 'End': end}, 'Markers':merged_data[['genome_pos','gene','sublin','drtype','drugs','weight']][start:end].sort_values(by=['weight']).to_dict('records')} # verbose version of output
    covered_positions[f'Amplicon_{x+1}'] = {'Range':{'Start': start, 'End': end}}  # concise version of output
    covered_ranges.append([start, end])
    # 
    merged_data.loc[(merged_data['genome_pos']>=merged_data.loc[start,:]['genome_pos']) & (merged_data['genome_pos']<=merged_data.loc[end,:]['genome_pos']), 'weight_y'] = 0 # set the weight of the covered positions to 0
    weight_window_sum = rolling_sum(merged_data['weight_y'].tolist(), read_size)


#%%
merged_data.loc[start,:]

#%%
covered_positions
#%%
covered_ranges

#%%
# Raw data trial
weight_window_sum = rolling_sum(raw_data['weight'].tolist(), read_size)

covered_positions = {}
covered_ranges = []
for x in range(read_number):
    if len(covered_ranges) == 0:
        start = np.argmax(weight_window_sum)
        end = start + read_size
    if len(covered_ranges) > 0:
        start = np.argmax(weight_window_sum)
        end = start + read_size
        for i in range(len(covered_ranges)):
            if start > covered_ranges[i][0] and start < covered_ranges[i][1]:
                start = covered_ranges[i][1]+1
                if covered_ranges[i][1]+1 + read_size > raw_data.shape[0]:
                    end = raw_data.shape[0]
                end = covered_ranges[i][1]+1 + read_size
            if end > covered_ranges[i][0] and end < covered_ranges[i][1]:
                end = covered_ranges[i][0]-1
                if covered_ranges[i][0]-1 - read_size < 0:
                    start = 0
                start = covered_ranges[i][0]-1 - read_size

    # covered_positions[f'Amplicon_{x+1}'] = {'Range':{'Start': start, 'End': end}, 'Markers':raw_data[['Pos','Gene','Lin','change_type','biotype','weight']][start:end].sort_values(by=['weight']).to_dict('records')} # verbose version of output
    covered_positions[f'Amplicon_{x+1}'] = {'Range':{'Start': start, 'End': end}}  # concise version of output
    covered_ranges.append([start, end])
    # 
    raw_data.loc[(raw_data['Pos']>=raw_data.loc[start,:]['Pos']) & (raw_data['Pos']<=raw_data.loc[end,:]['Pos']), 'weight'] = 0
    weight_window_sum = rolling_sum(raw_data['weight'].tolist(), read_size)



#%%
# def generate_genome(length: int) -> Genome:
#     return choices([0, 1], k=length)


# def generate_population(size: int, genome_length: int) -> Population:
#     return [generate_genome(genome_length) for _ in range(size)]


# def single_point_crossover(a: Genome, b: Genome) -> Tuple[Genome, Genome]:
#     if len(a) != len(b):
#         raise ValueError("Genomes a and b must be of same length")

#     length = len(a)
#     if length < 2:
#         return a, b

#     p = randint(1, length - 1)
#     return a[0:p] + b[p:], b[0:p] + a[p:]


# def mutation(genome: Genome, num: int = 1, probability: float = 0.5) -> Genome:
#     for _ in range(num):
#         index = randrange(len(genome))
#         genome[index] = genome[index] if random() > probability else abs(genome[index] - 1)
#     return genome


# def population_fitness(population: Population, fitness_func: FitnessFunc) -> int:
#     return sum([fitness_func(genome) for genome in population])


# def selection_pair(population: Population, fitness_func: FitnessFunc) -> Population:
#     return choices(
#         population=population,
#         weights=[fitness_func(gene) for gene in population],
#         k=2
#     )


def sort_population(population: Population, fitness_func: FitnessFunc) -> Population:
    return sorted(population, key=fitness_func, reverse=True)


def genome_to_string(genome: Genome) -> str:
    return "".join(map(str, genome))


def print_stats(population: Population, generation_id: int, fitness_func: FitnessFunc):
    print("GENERATION %02d" % generation_id)
    print("=============")
    print("Population: [%s]" % ", ".join([genome_to_string(gene) for gene in population]))
    print("Avg. Fitness: %f" % (population_fitness(population, fitness_func) / len(population)))
    sorted_population = sort_population(population, fitness_func)
    print(
        "Best: %s (%f)" % (genome_to_string(sorted_population[0]), fitness_func(sorted_population[0])))
    print("Worst: %s (%f)" % (genome_to_string(sorted_population[-1]),
                              fitness_func(sorted_population[-1])))
    print("")

    return sorted_population[0]

def run_evolution(
        populate_func: PopulateFunc,
        fitness_func: FitnessFunc,
        fitness_limit: int,
        selection_func: SelectionFunc = selection_pair,
        crossover_func: CrossoverFunc = single_point_crossover,
        mutation_func: MutationFunc = mutation,
        generation_limit: int = 100,
        printer: Optional[PrinterFunc] = None) \
        -> Tuple[Population, int]:
    population = populate_func()
    
    for i in range(generation_limit): 
        population = sorted(population, key=lambda genome: fitness_func(genome), reverse=True)
        
        if printer is not None:
            printer(population, i, fitness_func)

        if fitness_func(population[0]) >= fitness_limit:
            break

        next_generation = population[0:2]

        for j in range(int(len(population) / 2) - 1):
            parents = selection_func(population, fitness_func)
            offspring_a, offspring_b = crossover_func(parents[0], parents[1])
            offspring_a = mutation_func(offspring_a)
            offspring_b = mutation_func(offspring_b)
            next_generation += [offspring_a, offspring_b]

        population = next_generation

    return population, i
# %%
population, generations = run_evolution(
    PopulateFunc= partial(generate_population, size=10, genome_length=len(things)),
    fitness_func= partial(fitness, things=things, weight_limit=40),
    fitness_limit=40,
    generation_limit=100)
print(f"number of generations: {generations}")
# %%

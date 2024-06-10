#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 16:14:23 2023

@author: nathan
"""

import os
import tskit
import msprime
import csv
import numpy as np
from pathlib import Path


# hardcoding some parameters from the original fwdpy11 simulation
num_demes = 100
N0 = 10000
Nf = 500

# getting parameters fed to the original simulations
params = []
with open("/home/nathan/Documents/GitHub/path_integral/simulations/params.txt", 'r') as file:
    csvreader = csv.reader(file)
    for row in csvreader:
        params.append(row)
        
seeds = [int(params[i][0]) for i in range(len(params))]
mut_rate = [float(params[i][1]) for i in range(len(params))]
effect_sizes = [float(params[i][2]) for i in range(len(params))]

del(csvreader, file, params, row)

# working directory is where i exported all the trees from the original simulations
os.chdir("/media/nathan/T7/path_integral/simulations/out/trees(neut)")

# path = os.listdir()[1]
count = 0
for path in os.listdir():
    count += 1
    print(count,":",path)
    # match current tree with its parameters
    cur_seed = int(path.split("_")[1])
    cur_index = seeds.index(cur_seed)
    cur_mut_rate = mut_rate[cur_index]
    cur_effect_size = effect_sizes[cur_index]
    
    # loading tree
    ts = tskit.load(path + "/" + os.listdir(path)[0])
    # ts.num_mutations
    
    # replacing mutations
    ts = msprime.sim_mutations(tree_sequence=ts, 
                                  rate=cur_mut_rate * 1e-9,
                                  random_seed=cur_seed,
                                  model=msprime.BinaryMutationModel(),
                                  discrete_genome=False,
                                  keep=False)
    print(ts.num_mutations)
    
    # calculating allele frequencies
    times = set()
    for s in ts.samples():
        times.add(ts.node(s).time)
    
    allele_frequencies = {i: np.zeros((ts.num_sites, len(times))) for i in range(num_demes)}
    for j, t in enumerate(times):
        print(str(j) + " : " + str(t))
        samples = [s for s in ts.samples() if ts.node(s).time == t]
        ts_slice = ts.simplify(samples, filter_sites=False)
        G = ts_slice.genotype_matrix()
        if t == max(times):
            # initial generation before shift in lines
            afs = G.sum(axis=1) / 2 / N0
            for i in range(num_demes):
                allele_frequencies[i][:, j] = afs
        else:
            for i in range(num_demes):
                G_sub = G[:, i * 2 * Nf : (i + 1) * 2 * Nf]
                afs = G_sub.sum(axis=1) / 2 / Nf
                allele_frequencies[i][:, j] = afs
    
    out_dir = Path(path)
    
    out_files = {
        "start_freqs": out_dir / "start_freqs.csv",
        "end_freqs": out_dir / "end_freqs.csv",
    }
    
    with open(out_files["start_freqs"], "w") as outfile:
        writer = csv.writer(outfile)
        deme_list = list(allele_frequencies.keys())
        
        writer.writerow(deme_list)
        writer.writerow("start")
        # iterate each column and assign the
        # corresponding values to the column
        for i in range(ts.num_mutations):
                writer.writerow([allele_frequencies[x][i][1] for x in deme_list])
                
    with open(out_files["end_freqs"], "w") as outfile:
        writer = csv.writer(outfile)
        deme_list = list(allele_frequencies.keys())
        
        writer.writerow(deme_list)
        writer.writerow("end")
        # iterate each column and assign the
        # corresponding values to the column
        for i in range(ts.num_mutations):
                writer.writerow([allele_frequencies[x][i][0] for x in deme_list])
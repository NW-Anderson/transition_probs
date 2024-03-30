import fwdpy11
import numpy as np
from dataclasses import dataclass
from typing import List

import time
import demes

import os
import sys
import tarfile

## Set up parameters
# seed = 1
seed = int(sys.argv[1])
## TODO uncomment for cluster

# 1 Morgan chrosome
r = 1e-8
L = 1e8
assert r * L == 1

# VS always set to 1
VS = 1.0

# Change U in order to adjust VG
# U = 2.5e-3
U = float(sys.argv[2])
## TODO uncomment for cluster

expectedVG = 4 * U * VS
optimum = 0.0
shift = 1.0

# Set up selected mutations
# a = 0.01
a = float(sys.argv[3])
## TODO uncomment for cluster

print(f"Using the parameters: seed: {seed}, U: {U}, a: {a}")

sregions = [fwdpy11.ConstantS(0, L, 1, a), fwdpy11.ConstantS(0, L, 1, -a)]

# Load demographic model using demes
# os.chdir("/home/nathan/Documents/GitHub/path_integral")
g = demes.load("/opt/demo.yaml")

# import demesdraw
# demesdraw.tubes(g, colours=("blue"))

burnin = 10
model = fwdpy11.discrete_demography.from_demes(g, burnin=burnin)
simlen = model.metadata["total_simulation_length"]

# Setting up the population
initial_sizes = [
    model.metadata["initial_sizes"][i]
    for i in sorted(model.metadata["initial_sizes"].keys())
]
N0 = initial_sizes[0]
assert len(initial_sizes) == 1
Nf = g.demes[-1].epochs[0].start_size

pop = fwdpy11.DiploidPopulation(initial_sizes, L)

# after the burnin, the population is split into
# the replicate lines and the optimum shift occurs
GSSmo = fwdpy11.GSSmo(
    [
        fwdpy11.Optimum(when=0, optimum=optimum, VS=VS),
        fwdpy11.Optimum(when=model.metadata["burnin_time"], optimum=shift, VS=VS),
    ]
)

pdict = {
    "nregions": [],
    "sregions": sregions,
    "recregions": [fwdpy11.BinomialInterval(0, L, 1)],
    "rates": (0.0, U, None),
    "gvalue": fwdpy11.Additive(scaling=2, gvalue_to_fitness=GSSmo),
    "simlen": simlen,
    "demography": model,
    "prune_selected": False,
}
params = fwdpy11.ModelParams(**pdict)

## set up recorders
@dataclass
class SimData:
    generation: int
    demes_ids: List[int]
    mean_phenotype: List[float]
    mean_fitness: List[float]
    var_phenotype: List[float]


@dataclass
class Recorder:
    data: list

    def __call__(self, pop, sampler):
        md = np.array(pop.diploid_metadata)
        # general properties of the population
        # store lists of mean phenotypes and fitness, and var(pheno)
        deme_ids = sorted(list(set(md["deme"])))
        mean_pheno = [md[md["deme"] == i]["g"].mean() for i in deme_ids]
        mean_fitness = [md[md["deme"] == i]["w"].mean() for i in deme_ids]
        var_pheno = [md[md["deme"] == i]["g"].var() for i in deme_ids]
        self.data.append(
            SimData(pop.generation, deme_ids, mean_pheno, mean_fitness, var_pheno)
        )
        # record last generation of full population and of the replicate lines
        if pop.generation == pop.N * burnin or pop.generation == simlen:
            sampler.assign(np.arange(0, pop.N))


recorder = Recorder(data=[])

## Initialize and evolve full population
pop = fwdpy11.DiploidPopulation(N0, L)
rng = fwdpy11.GSLrng(seed)


time1 = time.time()
fwdpy11.evolvets(rng, pop, params, 100, recorder=recorder, suppress_table_indexing=True)
time2 = time.time()
print(f"simulation took {(time2 - time1)/60:0.2f} minutes")
assert pop.generation == simlen

ts = pop.dump_tables_to_tskit()
line_time = g.demes[0].end_time

## find mean VG, var VG and mean pheno for the final N0 generations of the burnin
gen = []
mean_pheno = []
mean_fit = []
var_pheno = []
# for i in range(N0 * (burnin - 1)+1, N0 * burnin+1): 
for i in range(simlen):
    gen.append(float(recorder.data[i].generation))
    mean_pheno.append(recorder.data[i].mean_phenotype[0])
    mean_fit.append(recorder.data[i].mean_fitness[0])
    var_pheno.append(recorder.data[i].var_phenotype[0])


# initial_VG = np.mean(var_pheno[N0 * (burnin - 1):N0 * burnin])
# print(initial_VG)

## The tree sequence should record the last generation
## prior to the split into replicate lines, and the final generation
print(f"The tree sequence now has {ts.num_mutations} mutations, at "
      f"{ts.num_sites} distinct sites.")
mut_times = []
for m in ts.mutations():
    mut_times.append(m.time)

# infinite sites, so this should be true
assert ts.num_mutations == ts.num_sites

# get times to compute frequencies at
times = set()
for s in ts.samples():
    times.add(ts.node(s).time)

times = sorted(list(times))[::-1]

# get number of replicate populations, minus 1 cause we dont want ancestral
num_demes = len(model.metadata["deme_labels"]) - 1

# calculating allele frequencies in each deme at the two generations of interest
allele_frequencies = {i: np.zeros((ts.num_sites, len(times))) for i in range(num_demes)}
for j, t in enumerate(times):
    print(j)
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

# get effect sizes of each mutation
effectsizes = []
for i in range(ts.num_mutations):
    effectsizes.append(ts.tables.mutations[i].metadata['s'])
    
# import matplotlib.pyplot as plt 
# f, ax = plt.subplots()
# ax.plot(gen[19900:20100], var_pheno[19900:20100], label="Genetic Variance")
# ax.set_xlabel("Generation")
# ax.set_ylabel("Value")
# plt.legend()
# plt.show()


# exporting results
import csv
from pathlib import Path

out_dir = Path(f"result_{seed}")
out_dir.mkdir()

out_files = {
    "start_freqs": out_dir / "start_freqs.csv",
    "end_freqs": out_dir / "end_freqs.csv",
    "effect_sizes": out_dir / "effect_sizes.csv",
    "popstats": out_dir / "popStats.csv",
    "trees": out_dir / "simulation.trees",
}

with open(out_files["start_freqs"], "w") as outfile:
    writer = csv.writer(outfile)
    deme_list = list(allele_frequencies.keys())
    
    writer.writerow(deme_list)
    writer.writerow("start")
    # iterate each column and assign the
    # corresponding values to the column
    for i in range(ts.num_mutations):
            writer.writerow([allele_frequencies[x][i][0] for x in deme_list])
            
with open(out_files["end_freqs"], "w") as outfile:
    writer = csv.writer(outfile)
    deme_list = list(allele_frequencies.keys())
    
    writer.writerow(deme_list)
    writer.writerow("end")
    # iterate each column and assign the
    # corresponding values to the column
    for i in range(ts.num_mutations):
            writer.writerow([allele_frequencies[x][i][1] for x in deme_list])

with open(out_files["effect_sizes"],"w") as outfile:
    writer = csv.writer(outfile)
    writer.writerow(effectsizes)
    
with open(out_files["popstats"], "w") as outfile:
    writer = csv.writer(outfile)
    writer.writerow(["gen","mean_pheno",  "var_pheno","mean_fit"])
    for i in range(len(gen)):
        writer.writerow([gen[i], mean_pheno[i], var_pheno[i], mean_fit[i]])
        
ts.dump(str(out_files["trees"]))

with tarfile.open(f"result_{seed}.tar.gz", mode="w:gz") as tar:
    for out_name in out_files.values():
        tar.add(out_name)

for out_name in out_files.values():
    out_name.unlink()

out_dir.rmdir()

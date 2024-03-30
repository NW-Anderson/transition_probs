import numpy as np
import os


def get_fprime(freqs, es, optimum, VS):
    # get expected frequencies of each mutation in the next generation
    two_f_a = 2 * freqs * es
    VG_contributions = 2 * freqs * (1 - freqs) * es * es
    EG = two_f_a.sum()
    VG = VG_contributions.sum()
    VGprime = VG - VG_contributions
    EGprime = EG - two_f_a

    delta_p = (
        freqs
        * (1 - freqs)
        * np.sqrt((VG + VS) / (VGprime + VS))
        * np.exp(((EG - optimum) ** 2) / (2 * (VG + VS)))  # nearly 1
        * (  # also nearly 1
            freqs * np.exp(-((EGprime + 2 * es - optimum) ** 2) / (2 * (VGprime + VS)))
            - (1 - freqs) * np.exp(-((EGprime - optimum) ** 2) / (2 * (VGprime + VS)))
            + (1 - 2 * freqs)
            * np.exp(-((EGprime + es - optimum) ** 2) / (2 * (VGprime + VS)))
        )
    )
    return freqs + delta_p


def sample_generation(Ne, fprime):
    # drift through binomial sampling, given expected frequencies in the next generation
    return np.random.binomial(2 * Ne, p=fprime) / 2 / Ne


def new_mutations(freqs, es, ss, mu, Ne, a):
    num_muts = np.random.poisson(2 * Ne * mu)
    extension = num_muts - (ss == 0).sum()
    if extension > 0:
        freqs = np.concatenate((freqs, np.zeros(extension)))
        es = np.concatenate((es, np.zeros(extension)))
        ss = np.concatenate((ss, np.zeros(extension, dtype=int)))
    free_sites = np.where(ss == 0)[0]
    for i in range(num_muts):
        site_idx = free_sites[i]
        # add to G
        freqs[site_idx] = 1 / 2 / Ne
        # add to seg sites
        ss[site_idx] = 1
        # add effect size
        es[site_idx] = (-1) ** np.random.randint(2) * a
    return freqs, es, ss


def evolve(freqs, es, ss, a, optimum, VS, mu, Ne):
    # get marginal fitnesses of alleles at all segregating sites
    fprime = get_fprime(freqs, es, optimum, VS)

    # binomial sampling of alleles to create offspring
    freqs = sample_generation(Ne, fprime)

    # update to remove lost mutations
    ss[freqs == 0] = 0

    # introduce new mutations
    freqs, es, ss = new_mutations(freqs, es, ss, mu, Ne, a)

    # clean up if many nonseg sites
    if np.sum(ss) < 0.9 * len(ss):
        freqs, es, ss = cleanup(freqs, es, ss)

    return freqs, es, ss


def evolve_line(freqs, es, ss, a, optline, VS, mu, Nline):
    # get marginal fitnesses of alleles at all segregating sites
    fprime = get_fprime(freqs, es, optline, VS)

    # binomial sampling of alleles to create offspring
    freqs = sample_generation(Nline, fprime)

    # introduce new mutations
    freqs, es, ss = new_mutations(freqs, es, ss, mu, Nline, a)

    return freqs, es, ss


def cleanup(freqs, es, ss):
    es = es.compress(ss == 1)
    if freqs.ndim == 1:
        freqs = freqs.compress(ss == 1)
    else:
        freqs = freqs.compress(ss == 1, axis=1)
    ss = ss.compress(ss == 1)
    return freqs, es, ss


def report(freqs, es, ss, VGs=None, report=True):
    EP = 2 * np.sum(freqs * es)
    VP = 2 * np.sum(freqs * (1 - freqs) * es * es)
    if report:
        print(f" num sites: {ss.sum()}; VG={VP:0.5f}; mean phenotype={EP:0.5f}")
    if VGs is not None:
        VGs.append(VP)


def write_data(j, meanVG, es, all_freqs, fname):
    if not os.path.exists(fname):
        fout = open(fname, "w+")
        # write header line (rep number, source VG, effect size, starting freq, all line freqs)
        fout.write(
            "rep,meanVG,a,x0,"
            + ",".join(["xf" + str(i) for i in range(len(all_freqs[0]) - 1)])
            + "\n"
        )
    else:
        fout = open(fname, "a+")
    for e, fs in zip(es, all_freqs):
        line_out = f"{j},{meanVG},{e}," + ",".join([str(f) for f in fs]) + "\n"
        fout.write(line_out)
    fout.close()


if __name__ == "__main__":
    # draw mutations from a normal
    optimum = 0
    VS = 1
    optline = 1

    # per-generation per-individual mutation rate
    mu = 0.025
    Ne = 10000
    Nline = 500

    a = 0.01
    numlines = 100
    gens = 50

    # filename to write data to
    fname = f"e_and_r.mu_{mu}.a_{a}.txt"

    # SHOC, just to compare
    expectedVG = 4 * mu * VS / (1 + VS / Ne / a ** 2)
    print("Expected VG from SHOC:", expectedVG)

    # initialize population
    freqs = np.array([])  # frequencies
    es = np.array([])  # effect sizes
    ss = np.array([])  # indicator of segregating sites

    # track VGs in full population
    VGs = []

    report(freqs, es, ss)
    for i in range(20 * Ne):
        # burn in from initial state
        freqs, es, ss = evolve(freqs, es, ss, a, optimum, VS, mu, Ne)
        report(freqs, es, ss, VGs=VGs, report=False)
        if i % Ne == 0:
            report(freqs, es, ss)

    report(freqs, es, ss)

    for j in range(100):
        print("Ave VG over past Ne generations:", np.mean(VGs[-Ne:]))
        freqs, es, ss = cleanup(freqs, es, ss)
        # run E&R lines, write the data, and then
        # advance the base population some number of generations

        # stores all frequencies of seg sites in initial population and the freq
        # of each mutation in each line
        all_freqs = np.empty((len(freqs), numlines + 1))
        all_freqs[:, 0] = freqs

        # evolve {gens} generations from the existing
        for rep in range(numlines):
            # copy to create a line
            freqs_rep = freqs.copy()
            es_rep = es.copy()
            ss_rep = ss.copy()
            # evolve the line
            for gen in range(gens):
                freqs_rep, es_rep, ss_rep = evolve_line(
                    freqs_rep, es_rep, ss_rep, a, optline, VS, mu, Nline
                )
            # store those final frequencies
            all_freqs[:, rep + 1] = freqs_rep[: len(freqs)]

        write_data(j, np.mean(VGs[-1000]),es, all_freqs, fname)

        # advance Ne generations
        for i in range(Ne):
            freqs, es, ss = evolve(freqs, es, ss, a, optimum, VS, mu, Ne)
            report(freqs, es, ss, VGs=VGs, report=False)
        report(freqs, es, ss)

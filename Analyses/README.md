Here's all the scripts that involve evaluating an analytic solution to the diffusion. The scripts are not elegant but they work. I output each datapoint as its own file that are merged later because my poor laptop has a habit of crashing mid analysis, and this is an easy way to save your work as you go. You'll need to put all of the output files from each script into their own directory

## ComparisonToExisting
Evaluates the polygenic adaptation, genic and neutral models under the same parameter regime and outputs the density of ending frequencies at grid points. This data is used to create Figure 2 in the V1 manuscript

## ComparisonToSims
Evaluate the different analytic solutions under the same parameter regimes as our simulations and outputs the density of ending frequencies. This data is used in Figures 3, S5-S7 of the V1 manuscript

## Convergence
Evaluate the polygenic selection analytic solution for different values of kmax and mmax, as well as the solution from NDSolve. Outputs the density of ending freqs. This data is used in figures S1-S4 of the V1 manuscript

## ErrorCalculations
Computes statistical distance between numerical and path integral solution under a specific error regime. Data used in Figure S10 of the V1 manuscript

## ProbDetection
Finds a neutral 99% afc threshold, then finds the probability an allele exceeds this threshold both numerically and using the analytic solution. Also computes statistical distance between the numerical and analytic solutions for each parameter regime. Data used in Figures 4, 5, S9, S11, S12 in the V1 manuscript

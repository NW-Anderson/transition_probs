# transition_probs
Code used in "A Path Integral Approach for Allele Frequency Dynamics Under Polygenic Selection"

Scripts are organized based on their purpose (e.g. simulations, or finding probability of extreme frequency changes), some results require running scripts located in multiple files (such as running a simulation in Simulations, then prepping the data to read into R using DataCleaning, then plotting the results in Figures).

## Analyses
Contains Mathematica scripts used to evaluate numerical or analytic solutions to the diffusion. 

## DataCleaning
Contains some miscellaneous scripts to clean the output of the simulation or mathematica scripts to use in plotting

## Figures 
Contains ggplot scripts to make figures, as well as pdfs of the figures used in the manuscript

## PathIntegrals
Contains scripts to compute the path integral solution to the diffusion (genic, underdominance and polygenic adaptation). These expressions are saved and then evaluated in "Analyses"

## Simulations
Contains scripts to run all of the simulations.

TODO: 
* Add docs for path_integrals
* Add docs for simulations
* Do something about the sim and Sim files merging or not depending on who looks at them
* 

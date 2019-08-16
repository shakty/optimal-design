# Parameter-Sampled GPUCB-PE

This code samples model parameters and simulates likely datasets that you would observe if experiment participants were parameterized as such. After sampling datasets, each outcome is assigned a likelihood based on its frequency in the sample. Using these likelihoods, the information gain is calculated by taking the Kullback-Leibler Divergence between the likelihoods of the datasets generated under several competing models.  


## Scripts

Start with the ps_gpucbpe_testbed.R file. It initializes the process that was used in "Fast Model-Selection through Adaptive Design of Experiments Maximizing Information Gain", recreating Figure 3c. 

ps_gpucbpe_testbed.R - 
ps_gpucbpe_simulateDatasets.R - this script houses the Parameter-Sampled GPUCB-PE code
ps_gpucbpe_models.R - the four models used in our model comparison, including three from El-Gamal & Palfrey (1995)
ps_gpucbpe_process.R - Gaussian Process script, using functions from GPfit, adapted for optimizing experimental design 
ps_gpucbpe_histories.R - enumerates all possible game histories 
ps_gpucbpe_helper.R - plotting, statistics, and file naming
ps_gpucbpe_calc_likelihoods.R - used for calculating the likelihoods of each dataset--assuming they have not been sampled
matches.R - enumerates possible pairings in the experiment

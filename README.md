# Parameter-sampled GPUCB-PE

R code for adaptively constructing and searching an "information landscape"
in order to find the particular experimental design that will maximize the
information gained as a result of running the experiment. This code accompanies
the recent paper: 

**Optimal design of experiments to identify latent behavioral types**\
Stefano Balietti, Brennan Klein, and Christoph Riedl, 2020.\
[arXiv:1807.07024](https://arxiv.org/abs/1807.07024)

This code samples model parameters and simulates likely datasets that you would 
observe if experiment participants were parameterized as such. After sampling 
datasets, each outcome is assigned a likelihood based on its frequency in the 
sample. Using these likelihoods, the information gain is calculated by taking 
the Kullback-Leibler Divergence between the likelihoods of the datasets 
generated under several competing models.  

## Scripts

Start with the ps_gpucbpe_testbed.R file. It initializes the process that was 
used in "Fast Model-Selection through Adaxptive Design of Experiments Maximizing 
Information Gain", recreating Figure 3c. 

1. [ps_gpucbpe_testbed.R](https://github.com/jkbren/ps-gpucbpe/blob/master/ps-gpucbpe/PS_GPUCBPE/R_CODE/ps_gpucbpe_testbed.R) - start here!
2. [ps_gpucbpe_simulateDatasets.R](https://github.com/jkbren/ps-gpucbpe/blob/master/ps-gpucbpe/PS_GPUCBPE/R_CODE/ps_gpucbpe_simulateDatasets.R) - this script houses the Parameter-Sampled G
PUCB-PE code
3. [ps_gpucbpe_models.R](https://github.com/jkbren/ps-gpucbpe/blob/master/ps-gpucbpe/PS_GPUCBPE/R_CODE/ps_gpucbpe_models.R) - the four models used in our model comparison, 
including three from El-Gamal & Palfrey (1995)
4. [ps_gpucbpe_process.R](https://github.com/jkbren/ps-gpucbpe/blob/master/ps-gpucbpe/PS_GPUCBPE/R_CODE/ps_gpucbpe_process.R) - Gaussian Process script, using functions from GPfit, 
adapted for optimizing experimental design 
5. [ps_gpucbpe_histories.R](https://github.com/jkbren/ps-gpucbpe/blob/master/ps-gpucbpe/PS_GPUCBPE/R_CODE/ps_gpucbpe_histories.R) - enumerates all possible game histories 
6. [ps_gpucbpe_helper.R](https://github.com/jkbren/ps-gpucbpe/blob/master/ps-gpucbpe/PS_GPUCBPE/R_CODE/ps_gpucbpe_helper.R) - plotting, statistics, and file naming
7. [ps_gpucbpe_calc_likelihoods.R](https://github.com/jkbren/ps-gpucbpe/blob/master/ps-gpucbpe/PS_GPUCBPE/R_CODE/ps_gpucbpe_calc_likelihoods.R) - used for calculating the likelihoods of each
dataset--assuming they have not been sampled
8. [matches.R](https://github.com/jkbren/ps-gpucbpe/blob/master/ps-gpucbpe/PS_GPUCBPE/R_CODE/matches.R) - enumerates possible pairings in the experiment

## Citation   <a name="citation"/>

If you use these methods and this code in your own research, 
please cite our paper:

Balietti, S., Klein, B. & Riedl, C. (2020). **Optimal design of experiments to identify latent behavioral types**. 
[arXiv:1807.07024](https://arxiv.org/abs/1807.07024)

Bibtex: 
```text
@article{balietti_optimaldesign_2020,
  title={Optimal Design of Experiments to Identify Latent Behavioral Types},
  author={Balietti, S.; Klein, B. and Riedl, C.},
  journal={arXiv preprint arXiv:1807.07024},
  year={2020}
}
```

## See also:

* El-Gamal, M. A., & Palfrey, T. R. (1996). **Economical experiments: Bayesian 
efficient experimental design**. *International Journal of Game Theory*, 25(4), 
495-517. doi: [10.1007/BF01803953](https://link.springer.com/article/10.1007/BF01803953).
    + motivating work from which many of these ideas derive

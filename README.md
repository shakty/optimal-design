# Optimal design of experiments to identify latent behavioral types

This repository contains R code for adaptively constructing and searching an "information landscape" in order to find the particular experimental design that will maximize the information gained as a result of running the experiment. This code accompanies the recent paper:

**Optimal design of experiments to identify latent behavioral types**\
Stefano Balietti, Brennan Klein, and Christoph Riedl, \
[Experimental Economics](https://link.springer.com/article/10.1007%2Fs10683-020-09680-w) 24, pages 772–799 (2021)\
[arXiv:1807.07024](https://arxiv.org/abs/1807.07024)

## Parameter-sampled GPUCB-PE

This code samples model parameters and simulates likely datasets that you would
observe if experiment participants were parameterized as such. After sampling
datasets, each outcome is assigned a likelihood based on its frequency in the
sample. Using these likelihoods, the information gain is calculated by taking
the Kullback-Leibler Divergence between the likelihoods of the datasets
generated under several competing models.  

### Scripts

Start with the ps_gpucbpe_testbed.R file. It initializes the process that was
used in "Fast Model-Selection through Adaxptive Design of Experiments Maximizing
Information Gain", recreating Figure 3c.

1. [main.R](https://github.com/shakty/optimal-design/blob/master/R_CODE/main.R) - start here!
2. [simulateDatasets.R](https://github.com/shakty/optimal-design/blob/master/R_CODE/ps_gpucbpatasets.R) - this script houses the Parameter-Sampled G
PUCB-PE code
3. [models.R](https://github.com/shakty/optimal-design/blob/master/R_CODE/models.R) - the four models used in our model comparison,
including three from El-Gamal & Palfrey (1995)
4. [process.R](https://github.com/shakty/optimal-design/blob/master/R_CODE/process.R) - Gaussian Process script, using functions from GPfit,
adapted for optimizing experimental design
5. [histories.R](https://github.com/shakty/optimal-design/blob/master/R_CODE/psstories.R) - enumerates all possible game histories
6. [helper.R](https://github.com/shakty/optimal-design/blob/master/R_CODE_helper.R) - plotting, statistics, and file naming
7. [calc_likelihoods.R](https://github.com/shakty/optimal-design/blob/master/R_CODE/ps_gpucbplihoods.R) - used for calculating the likelihoods of each
dataset--assuming they have not been sampled
8. [matches.R](https://github.com/shakty/optimal-design/blob/master/ps-gpucmatches.R) - enumerates possible pairings in the experiment


## Citation   <a name="citation"/>

If you use these methods and this code in your own research,
please cite our paper:

Balietti, S., Klein, B. & Riedl, C. (2021). **Optimal design of experiments to identify latent behavioral types**, _[Experimental Economics](https://link.springer.com/article/10.1007%2Fs10683-020-09680-w)_ 24, pages 772–799


Bibtex:
```text
@article{balietti_optimaldesign_2021,
  title={Optimal Design of Experiments to Identify Latent Behavioral Types},
  author={Balietti, S.; Klein, B. and Riedl, C.},
  journal={Experimental Economics},
  issue={24},
  year={2021},
  pages={772--799}
}
```

## See also:

* El-Gamal, M. A., & Palfrey, T. R. (1996). **Economical experiments: Bayesian
efficient experimental design**. *International Journal of Game Theory*, 25(4),
495-517. doi: [10.1007/BF01803953](https://link.springer.com/article/10.1007/BF01803953).
    + motivating work from which many of these ideas derive

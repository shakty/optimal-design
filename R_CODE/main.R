##############################################
###                                        ###
###   Adaptive Design of Experiments for   ###
###       Maximizing Information Gain      ###
###        code by: Brennan Klein          ###
###                                        ###
##############################################

## Associated Manuscript:
## Balietti, S., Klein, B. & Riedl, C. (2020).
## Optimal design of experiments to identify latent behavioral types.
## https://arxiv.org/abs/1807.07024


# libraries
library(ggplot2)
library(stringr)
library(gridExtra)
library(GPfit)
library(RColorBrewer)
library(lhs)
library(latex2exp)
library(scales)
library(lhs)
suppressMessages(library(data.table))
suppressMessages(library(randtoolbox))
library(Matrix, quietly=T, warn.conflicts=F)

# Set working directory as needed.
## wd_path <- "~/Desktop/code_ParamSampledGPUCBPE/"
## setwd(wd_path)

# Source other files.
source("helper.R")
source("simulateDatasets.R")
source("process.R")
source("histories.R")
source("models.R")
source("calc_likelihoods.R")

# Initialize your experiment.
#############################

# How precise do you want your information landscape to be?
## (creates a 41x41 grid of possible experiments).
n_dim <- 41

expParam1_range <- c(2.0, 6.0)    # "A" in the manuscript.
expParam2_range <- c(0.2, 0.8)    # "pi" in the manuscript.

# Space of all possible experiments.
seqA <-  seq(expParam1_range[1], expParam1_range[2], length.out=n_dim)
seqPI <- seq(expParam2_range[1], expParam2_range[2], length.out=n_dim)
experiment_grid <- expand.grid(A = seqA, PI = seqPI)

sampling_method <- "Uniform"        # How to sample model parameters.
sampling <- "SamplePost_SampleHist" # Only value accepted for now.

# GPUCB-PE parameters.
######################

L <- 4            # Number of previous searches back to look for stopping rule
failsafe <- 250   # Number of iterations, in case the stopping rule doesnt work
rho0 <- 0.001     # Stopping Threshold for GPUCB_PE
k <- 1            # Number of "explore" searches per run of the algorithm

# Information calculation.
##########################

## When calculating the KL divergence, is it symmetric or not?
## Note: KL distance from A to B is generally different from distance
## from B to A. The symmetric version takes the average of the two distances.
asymmetric <- T


# Initializing parameters for search process.
#############################################
n_init      <- 20           # How many seed points to sample?
sobol       <- "Sobol"      # Distribution initial seeds: "Random" or "Sobol"
search      <- "GPUCBPE"    # Algorithm to search with: Random, Grid, GPUCBPE
n_samples   <- 1000         # Num. histories to sample per search: 1000, 5000, 10000
seed        <- 666          # Random seed
set.seed(seed)              # Set the random seed

## Initialize game details.
model_nums <- c(1,2,3, "all") # (Unsure if needed).
models_used <- c(1,2,3)       # (Unsure if needed).

## If asymmetric, this model is taken as reference.
current_top_model <- model_nums[1]

n_players <- 10
if (n_players%%2 !=0) {
    stop("HEY YOU NEED n_players TO BE EVEN")
}
n_pairs <- n_players/2
n_rounds <- 3

posteriors <- read.csv("../DATA_INPUT/posteriors_elgml.csv")
jet.colors <- rev(colorRampPalette(brewer.pal(11, "Spectral"))(100))

#################################
# Now let's practice on a grid. #
#################################

history <- initialize_process(posteriors, n_init, expParam1_range,
                              expParam2_range, n_dim, sampling_method,
                              n_players, n_rounds, n_samples, sobol, sampling,
                              asymmetric, L, models_used, current_top_model)


plot_landscape(history, n_samples, n_init, sobol, n_dim, expParam1_range,
               expParam2_range)

history <- run_algorithm_once(history, experiment_grid, posteriors, rho0,
                              failsafe, L, n_samples, k, asymmetric, sampling,
                              search, n_init, sampling_method, sobol, seed,
                              models_used, current_top_model)

plot_landscape(history, n_samples, n_init, sobol, n_dim, expParam1_range,
               expParam2_range)

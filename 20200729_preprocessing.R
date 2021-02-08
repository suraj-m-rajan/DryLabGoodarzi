################################################################################
# 20200729_preprocessing.R
# Author: Mitchel Cole
# Using bed files generated in 20200713_m6a_feature_exploration.Rmd, this script
# obtains 201 bp sequence and conservation scores associated with each entry in
# the bed file (HEK, CD8, and A549 miCLIP data from Linders 2015 and Ke 2015)
# The output will be used for modeling in a subsequent ipynb
################################################################################

library(purrr)
library(readr)

source("R/20200730_m6a_project_functions.R")

samples <- c("a549","cd8t","hek293")

isBackground_full <- rep(c(TRUE, FALSE), times = length(samples))
samples_full <- rep(samples, each = 2)

map2_dfr(samples_full, isBackground_full, createInputData) %>%
  write_csv("processed_data/201_bp_and_phylo_ml_data.csv.gz")
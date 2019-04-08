# This script peforms stepwise QTL mapping with R/QTL for one of seven traits.
# Which trait to map is determined by changing the index for the variable ix.
# This is based on input R/QTL files created with the script /R/rqtl_files.R

library("snow")
library("qtl")
# parameters to input
set.seed(249) # random seed.

source('R/perform_mapping.R')
perform_mapping("ffit_reduced", nperms = 1000, nclusters = 7)
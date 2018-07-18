# This script peforms MQM mapping with R/QTL for one of seven traits.
# Which trait to map is determined by changing the index for the variable ix.
# This is based on input R/QTL files created with the script /R/rqtl_files.R,
# Also on genome-wide significance thresholds determined by permutation.
# Permutations are performed by a corresponding script in the same folder.

# This is designed to run on the Uppsala University computing cluster Rackham.
# This does not allow using R studio projects, so setting the working directory one level up keeps paths the same.
setwd('..')

# There are six phenotypes:
traits <- c(
  "mass", # 1. seed mass (mass)
  "frut", # 2. fruits/plant ()
  "seed", # 3. seeds/fruit
  "tofu", # 4. Total fecundity (estimated seeds/plant = frut*seed)
  "ffit", # 5. Fruits per planted seedling; fitness measure used by Agren etal 2013
  "tfit", # 6. Seeds per planted seedling; fitness incorporating seed number.
  "surv") # 7. Survival
# parameters to input
cl <- 16
pm <- 10000

cols <- c("it2010", "it2011", "sw2010", "sw2011")

source("R/perform_mqm.R")
library("qtl",  lib.loc='~/R/x86_64-redhat-linux-gnu-library/3.5')
library("snow", lib.loc='~/R/x86_64-redhat-linux-gnu-library/3.5')

perform_mqm(traits[4], phenocol = cols, nclusters = cl, nperms = pm)


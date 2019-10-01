# This script peforms stepwise QTL mapping with R/QTL for one of seven traits.
# Which trait to map is determined by changing the index for the variable ix.
# This is based on input R/QTL files created with the script /R/rqtl_files.R

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
set.seed(249) # random seed.

source('R/perform_mapping.R')
perform_mapping(traits[4], nperms = 10000, nclusters = 16)

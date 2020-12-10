#' Tom Ellis
#' 
#' Script to create two sets of CSV files in a consistent format:
#' 1. A set of files giving RIL means for each trait, primarily to calculate 
#' genetic correlations.
#' 2. A set of files that can be read by R/QTL for linkage mapping.

library(qtl)
library(arghqtl)

# genotype data
geno <- read.csv('data_raw/RIL genotypes.csv')
# seed trait data
seed <- mass <- read.csv("data_raw/RIL_means_seedmass_seednumber.csv")

# subset columns.
seed <- seed[, 1:5]
mass <- mass[, c(1,10:13)]

# Fitness component data
surv <- frut <- read.csv("data_raw/Agren2013_fitness_components.csv") # fruit and survival
frut <- frut[, c("id",
                 "F10_IT_NFrFrPr_noedge_mean","F11_IT_NFrFrPr_noedge_mean",
                 "F10_SW_NFrFrPr_noedge_mean","F11_SW_NFrFrPr_noedge_mean")]
surv <- surv[, c("id",
                 "F10_IT_Surv_noedge_mean","F11_IT_Surv_noedge_mean",
                 "F10_SW_Surv_noedge_mean", "F11_SW_Surv_noedge_mean")]
# fruits per seedling
ffit <- read.csv('data_raw/Agren2013_LSMfitness.csv')
ffit <- ffit[,c("id","F10_NFruitIt_noedge_lsm", "F11_NFruitIt_noedge_lsm", "F10_NFruitSw_noedge_lsm", "F11_NFruitSw_noedge_lsm")]

# rename columns
colnms <- c("id", "it2010", "it2011", "sw2010", "sw2011")
colnames(seed) <- colnames(mass) <- colnames(frut) <- colnames(ffit) <- colnames(surv) <- colnms

# reorder rows to a common order.
rilnames <- sort(unique(c(frut$id, seed$id, ffit$id)))
mass <- mass[match(rilnames, mass$id),]
seed <- seed[match(rilnames, seed$id),]
frut <- frut[match(rilnames, frut$id),]
ffit <- ffit[match(rilnames, ffit$id),]
surv <- surv[match(rilnames, surv$id),]

# seeds per adult
tofu <- seed[,2:5] * frut[,2:5]
tofu <- cbind(genotype = frut$id, tofu)
# seeds per seedling
tfit <- seed[,2:5] * ffit[,2:5]
tfit <- cbind(genotype = ffit$id, tfit)

# write phenotypes to disk
write.csv(surv, file = "data_derived/RIL_survival.csv", row.names = F)
write.csv(seed, file = "data_derived/RIL_seeds_per_fruit.csv", row.names = F)
write.csv(mass, file = "data_derived/RIL_seed_mass.csv", row.names = F)
write.csv(frut, file = "data_derived/RIL_fruits_per_adult.csv", row.names = F)
write.csv(tofu, file = "data_derived/RIL_seeds_per_adult.csv", row.names = F)
write.csv(ffit, file = "data_derived/RIL_fruits_per_seedling.csv", row.names = F)
write.csv(tfit, file = "data_derived/RIL_seeds_per_seedling.csv", row.names = F)

# For QTL mapping of components of fitness, ensure that all traits are based
# only on lines for which seed numbers are available.
nax <- is.na(mass) + is.na(seed) + is.na(frut) + is.na(tofu) + is.na(surv)
nax <- nax > 0

mass[nax] <- NA
seed[nax] <- NA
frut[nax] <- NA
tofu[nax] <- NA
surv[nax] <- NA
#ffit[nax] <- NA

# For fruits per seedling, remove lines that have no data on seeds per seedling
ffit[is.na(tfit)] <- NA

# Files for R/QTL for new traits, and total seed number traits.
# This will throw a warning about NA labels.
# This is because not all lines are included across datasets.
write.csv(qtl_dataframe(mass, geno), file = "data_derived/rqtl_mass.csv", row.names = F)
write.csv(qtl_dataframe(seed, geno), file = "data_derived/rqtl_seed.csv", row.names = F)
write.csv(qtl_dataframe(frut, geno), file = "data_derived/rqtl_frut.csv", row.names = F)
write.csv(qtl_dataframe(tofu, geno), file = "data_derived/rqtl_tofu.csv", row.names = F)
write.csv(qtl_dataframe(ffit, geno), file = "data_derived/rqtl_ffit.csv", row.names = F)
write.csv(qtl_dataframe(tfit, geno), file = "data_derived/rqtl_tfit.csv", row.names = F)
write.csv(qtl_dataframe(surv, geno), file = "data_derived/rqtl_surv.csv", row.names = F)
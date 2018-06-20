library(qtl)
library(qtltools)

# import large dataframe with all the data.
bigtable <- read.csv("data_files/RILFieldExptItalySweden_SeedSizeNumber_corr2_SixYrs_398RIL By Genotype.csv")
surv <- frut <- read.csv("data_files/MeansFitnessComp_noedge_20092014.csv") # fruit and survival
geno <- read.csv('data_files/RIL genotypes.csv') # genotype data
# subset columns.
numb <- bigtable[, c(1, 3,4,8,9,13,14,18,19)]
mass <- bigtable[, c(1,23,24,28,29,33,34,38,39)]
frut <- frut[, c(1,16,18,17,19)]
surv <- surv[, c(1,4,6,5,7)]
# set values for lines with <3 replicates to NA
for(c in 2:5) numb[numb[,c+4] < 3 ,c] <- NA
for(c in 2:5) mass[mass[,c+4] < 3 ,c] <- NA
# rename columns
colnms <- c("genotype", "mean_it2010", "mean_it2011", "mean_sw2010", "mean_sw2011")
colnames(numb)[1:5] <- colnames(mass)[1:5] <- colnames(frut) <- colnames(frut) <- colnames(surv) <- colnms
# multiply seed numb per fruit by fruit number
tofu <- numb[match(frut$genotype, numb$genotype),2:5] * frut[,2:5]
tofu <- cbind(genotype = frut$genotype, tofu)

# create copies of the fruit and tofu arrays for fitness; i.e. coerce NAs to zero
ffit <- frut; tfit <- tofu
ffit[is.na(ffit)] <- 0
tfit[is.na(tfit)] <- 0

# write phenotypes to disk
write.csv(surv, file = "data_files/RIL_survival.csv", row.names = F)
write.csv(numb, file = "data_files/RIL_seeds_per_fruit.csv", row.names = F)
write.csv(mass, file = "data_files/RIL_seed_mass.csv", row.names = F)
write.csv(frut, file = "data_files/RIL_fruits_per_adult.csv", row.names = F)
write.csv(tofu, file = "data_files/RIL_seeds_per_adult.csv", row.names = F)
write.csv(ffit, file = "data_files/RIL_fruits_per_seedling.csv", row.names = F)
write.csv(tfit, file = "data_files/RIL_seeds_per_seedling.csv", row.names = F)

write.csv(qtl_dataframe(surv, geno), file = "data_files/rqtl_surv.csv", row.names = F)
write.csv(qtl_dataframe(mass, geno), file = "data_files/rqtl_mass.csv", row.names = F)
write.csv(qtl_dataframe(numb, geno), file = "data_files/rqtl_numb.csv", row.names = F)
write.csv(qtl_dataframe(frut, geno), file = "data_files/rqtl_frut.csv", row.names = F)
write.csv(qtl_dataframe(tofu, geno), file = "data_files/rqtl_tofu.csv", row.names = F)
write.csv(qtl_dataframe(ffit, geno), file = "data_files/rqtl_ffit.csv", row.names = F)
write.csv(qtl_dataframe(tfit, geno), file = "data_files/rqtl_tfit.csv", row.names = F)
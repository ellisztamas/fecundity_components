# This script peforms stepwise multiple QTL mapping with R/QTL for one of six traits.
# Which trait to map is determined by changing the index for the variable ix.
# This is based on input R/QTL files created with the script /R/rqtl_files.R
# This is designed to run on the Uppsala University computing cluster Milou.

# There are six phenotypes:
traits <- c(
  "mass",# 1. seed mass (mass)
  "frut", # 2. fruits/plant ()
  "seed", # 3. seeds/fruit
  "tofu", # 4. Total fecundity (estimated seeds/plant = frut*seed)
  "ffit", # 5. Fruits per planted seedling; fitness measure used by Agren etal 2013
  "tfit") # 6. Seeds per planted seedling; fitness incorporating seed number.
# parameters to input
ix <- traits[5] # pull out the name of the trait we are mapping here
set.seed(13) # random seed.
nperms <- 10000

path <- paste("data_files/rqtl_", ix, ".csv", sep="")
cat("Stepwise QTL analysis for R/QTL file", path,"\n")
cat("Started", format(Sys.time(),usetz = TRUE), "\n\n")

t0 <- proc.time() # save the system time.
###############
# Import data #
###############
library("qtl",  lib.loc='~/R/x86_64-unknown-linux-gnu-library/3.0/')
library("snow", lib.loc='~/R/x86_64-unknown-linux-gnu-library/3.0/')
# create r/qtl object from the file created above.
cross <- read.cross("csv", "", path, na.strings=NA, genotypes = c("a", "b"), alleles=c("It","Sw"))
cross <- convert2riself(cross) # state that this is RIL data from a selfer

# calculate the genotype probabilities
cross <- calc.genoprob(cross, step=2, error.prob=.0001, map.function="kosambi",stepwidth="max", off.end=0)

# perform quantile normal transformations on each variable.
cat("Performing quantile normal transformations.\n")
source('R/quantnorm.R')
for(c in 2:5) cross$pheno[,c] <- quantnorm(cross$pheno[,c])

################
# Permutations #
################
# Perumtations for seed fruter.
cat("Starting permutations.\n")
# Empty list to store permutations for each year.
permutations <- list(it2010=NULL, it2011=NULL, sw2010=NULL, sw2011=NULL)
# loop across each site-year combination and perform permutations for each.
for(y in 1:length(permutations)){
  permutations[[y]] <- scantwo(cross, pheno.col = y+1, method = "hk", n.perm = nperms, n.cluster = 8, verbose = F)
  cat("Permutations completed for", names(permutations)[y], "after", round(((proc.time()-t0)[3] /3600),2),"hours.\n")
}
# store permutations for later
save(permutations, file= paste("output/permutations_", ix, ".csv", sep=""))

#######################
# Calculate penalties #
#######################
cat("Calculating penalties.\n")
# Empty list to store penalties for each site-year combination.
penalties <- list(it2010=NULL, it2011=NULL, sw2010=NULL, sw2011=NULL)
# get penalties from the permutations
for(y in 1:length(penalties)) penalties[[y]] <- calc.penalties(permutations[[y]], alpha = 0.05)

########################
# Stepwise QTL fitting #
########################
cat("Fitting stepwise models.\n")
# Empty list to store penalties for each site-year combination.
stepwise <- list(it2010=NULL, it2011=NULL, sw2010=NULL, sw2011=NULL)
# loop across each site-year combination and perform QTL search for each.
for(y in 1:length(stepwise)){
  stepwise[[y]] <- stepwiseqtl(cross, pheno.col=y+1, penalties=penalties[[y]], method="hk", model="normal",
                                    max.qtl=10, covar=NULL, scan.pairs=T, additive.only=F, keeplodprofile=T, keeptrace=T,
                                    refine.locations=T, verbose=F)
  cat("QTL search completed for", names(stepwise)[y], "after",  round(((proc.time()-t0)[3] /3600),2),"hours.\n")
}

# save the stepwise QTL models.
save(stepwise, file=paste("output/stepwise_", ix, ".Rdata", sep=""))

cat("Analyses completed.")

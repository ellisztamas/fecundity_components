# This script peforms stepwise multiple QTL mapping with R/QTL for one of seven traits.
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
ix <- traits[1] # pull out the name of the trait we are mapping here

t0 <- proc.time() # save the system time.
cat("Stepwise QTL analysis for R/QTL file", path,"\n")
cat("Started", format(Sys.time(),usetz = TRUE), "\n\n")

###############
# Import data #
###############
library("qtl",  lib.loc='~/R/x86_64-redhat-linux-gnu-library/3.5')
library("snow", lib.loc='~/R/x86_64-redhat-linux-gnu-library/3.5')
# import cross object
path <- paste("data_files/rqtl_", ix, ".csv", sep="")
# create r/qtl object from the file created above.
cross <- read.cross("csv", "", path, na.strings=NA, genotypes = c("a", "b"), alleles=c("It","Sw"))
cross <- convert2riself(cross) # state that this is RIL data from a selfer
# calculate the genotype probabilities
cross <- calc.genoprob(cross, step=2, error.prob=.0001, map.function="kosambi",stepwidth="max", off.end=0)

# perform quantile normal transformations on each variable.
cat("Performing quantile normal transformations.\n")
source('R/quantnorm.R')
for(c in 2:5) cross$pheno[,c] <- quantnorm(cross$pheno[,c])

penalties <- readRDS(file = "output/penalties_surv.rds")


###################
# Stepwise models #
###################
cat("Fitting stepwise models.\n")
# Empty list to store penalties for each site-year combination.
stepwise <- list(it2010=NULL, it2011=NULL, sw2010=NULL, sw2011=NULL)
# loop across each site-year combination and perform QTL search for each.
for(y in 1:length(stepwise)){
  stepwise[[y]] <- stepwiseqtl(cross, pheno.col=y+1, penalties=penalties[[y]], method="hk", model="normal",
                               max.qtl=20, covar=NULL, scan.pairs=T, additive.only=F, keeplodprofile=T, keeptrace=T,
                               refine.locations=T, verbose=F)
  cat("QTL search completed for", names(stepwise)[y], "after",  round(((proc.time()-t0)[3] /3600),2),"hours.\n")
}

# save the stepwise QTL models.
saveRDS(stepwise, file=paste("output/stepwise_", ix, ".rds", sep=""))

cat("Analyses completed.")

#' Tom Ellis
#' 
#' Script to estimate confidence intervals around selection coefficients between
#' parental ecotypes in each experiment by non-parametric bootstrapping. This 
#' defines function `sboot` to take a single bootstrap subsample of the 
#' observed data on individual plants and calculate the new selection 
#' coefficient. This is then repeated many times.

source("R/001.parents.R")

siteyear <- sort(unique(seed_par$siteyear))
frut_par$block <- paste(frut_par$year, frut_par$site, frut_par$tray)
seed_par$block <- paste(seed_par$year, seed_par$site, seed_par$tray)

# function to take a single bootstrap resample of the observed data
sboot <- function(group, nreps){
  dfrut <- frut_par[frut_par$siteyear == group,]
  dseed <- seed_par[seed_par$siteyear == group,]
  dt <- matrix(0, nreps, 5)
  for(i in 1:nreps){
    ix <- tapply(1:nrow(dfrut), dfrut$block, function(x) sample(x, size=length(x), replace = T))
    ix <- do.call('c', ix)
    this_exp1 <- dfrut[ix,]
    this_exp2 <- dseed[ix,]
    pmeans <- data.frame(surv = tapply(this_exp1$survival,        this_exp1$genotype, mean, na.rm=T),
                         frut = tapply(this_exp1$fruit_per_plant, this_exp1$genotype, mean, na.rm=T),
                         seed = tapply(this_exp2$nseeds,          this_exp2$genotype, mean, na.rm=T),
                         ffit = tapply(this_exp1$fruit_per_seed,  this_exp1$genotype, mean, na.rm=T),
                         tfit = tapply(this_exp1$tfit,            this_exp1$genotype, mean, na.rm=T))
    pmeans <- t(pmeans) / apply(pmeans,2,max)
    if("italy" %in% strsplit(group, " ")[[1]])  pmeans <- (pmeans[,1]-pmeans[,2])
    if("sweden" %in% strsplit(group, " ")[[1]]) pmeans <- (pmeans[,2]-pmeans[,1])
    dt[i, ] <- as.numeric(pmeans)
  }
  return(dt)
}

# perform bootsraps
nreps <- 1000
sel_boots <- lapply(siteyear, function(x) sboot(x, nreps))
names(sel_boots) <- c("it2010", "it2011", "sw2010", "sw2011")

saveRDS(sel_boots, "output/bootstrap_selection_coefficients.rds")

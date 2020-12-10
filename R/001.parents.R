#' Tom Ellis
#' 
#' Script to summarise means, standard errors and sample sizes of parental
#' ecotypes for each trait. Selection coefficients are then calculated.
#' This script contains a lot of very ugly base-R calls to `tapply` because I
#' wrote it before I became aware of tidyverse tools. If you have downloaded 
#' this code at some point in the future, please don't judge me. 

# Data on individual plants of the parental ecotypes.
seed_par <- read.csv(file = 'data_raw/individual_parents_massnumber.csv')
frut_par <- read.csv('data_raw/individual_parents_nfruit.csv')
# site-year identifers
frut_par$siteyear <- paste(frut_par$site, frut_par$year)
seed_par$siteyear <- paste(seed_par$site, seed_par$year)
# survival
frut_par$survival <- as.integer(!is.na(frut_par$fruit_per_plant))

parents <- data.frame(
  variable = rep(c("seed_mass", "seeds_per_fruit", "fruit_per_RP","fruit_per_seedling","survival"), each=4),
  site     = rep(c("italy", "sweden"), 5, each=2),
  year     = rep(rep(c(2010:2011), 2), 5),
  it_parent= c(
    tapply(seed_par$mean_mass_ug   [seed_par$genotype=="italy"], seed_par$siteyear[seed_par$genotype=="italy"], mean, na.rm=T),
    tapply(seed_par$nseeds         [seed_par$genotype=="italy"], seed_par$siteyear[seed_par$genotype=="italy"], mean, na.rm=T),
    tapply(frut_par$fruit_per_plant[frut_par$genotype=="italy"], frut_par$siteyear[frut_par$genotype=="italy"], mean, na.rm=T),
    tapply(frut_par$fruit_per_seed [frut_par$genotype=="italy"], frut_par$siteyear[frut_par$genotype=="italy"], mean, na.rm=T),
    tapply(frut_par$survival       [frut_par$genotype=="italy"], frut_par$siteyear[frut_par$genotype=="italy"], mean, na.rm=T)),
  it_se = c(
    tapply(seed_par$mean_mass_ug   [seed_par$genotype=="italy"], seed_par$siteyear[seed_par$genotype=="italy"], sd, na.rm=T),
    tapply(seed_par$nseeds         [seed_par$genotype=="italy"], seed_par$siteyear[seed_par$genotype=="italy"], sd, na.rm=T),
    tapply(frut_par$fruit_per_plant[frut_par$genotype=="italy"], frut_par$siteyear[frut_par$genotype=="italy"], sd, na.rm=T),
    tapply(frut_par$fruit_per_seed [frut_par$genotype=="italy"], frut_par$siteyear[frut_par$genotype=="italy"], sd, na.rm=T),
    tapply(frut_par$survival       [frut_par$genotype=="italy"], frut_par$siteyear[frut_par$genotype=="italy"], sd, na.rm=T)),
  it_n = c(
    tapply(!is.na(seed_par$mean_mass_ug   [seed_par$genotype=="italy"]), seed_par$siteyear[seed_par$genotype=="italy"], sum),
    tapply(!is.na(seed_par$nseeds         [seed_par$genotype=="italy"]), seed_par$siteyear[seed_par$genotype=="italy"], sum),
    tapply(!is.na(frut_par$fruit_per_plant[frut_par$genotype=="italy"]), frut_par$siteyear[frut_par$genotype=="italy"], sum),
    tapply(!is.na(frut_par$fruit_per_seed [frut_par$genotype=="italy"]), frut_par$siteyear[frut_par$genotype=="italy"], sum),
    tapply(!is.na(frut_par$survival       [frut_par$genotype=="italy"]), frut_par$siteyear[frut_par$genotype=="italy"], sum)),
  sw_parent= c(
    tapply(seed_par$mean_mass_ug   [seed_par$genotype=="sweden"], seed_par$siteyear[seed_par$genotype=="sweden"], mean, na.rm=T),
    tapply(seed_par$nseeds         [seed_par$genotype=="sweden"], seed_par$siteyear[seed_par$genotype=="sweden"], mean, na.rm=T),
    tapply(frut_par$fruit_per_plant[frut_par$genotype=="sweden"], frut_par$siteyear[frut_par$genotype=="sweden"], mean, na.rm=T),
    tapply(frut_par$fruit_per_seed [frut_par$genotype=="sweden"], frut_par$siteyear[frut_par$genotype=="sweden"], mean, na.rm=T),
    tapply(frut_par$survival       [frut_par$genotype=="sweden"], frut_par$siteyear[frut_par$genotype=="sweden"], mean, na.rm=T)),
  sw_se= c(
    tapply(seed_par$mean_mass_ug   [seed_par$genotype=="sweden"], seed_par$siteyear[seed_par$genotype=="sweden"], sd, na.rm=T),
    tapply(seed_par$nseeds         [seed_par$genotype=="sweden"], seed_par$siteyear[seed_par$genotype=="sweden"], sd, na.rm=T),
    tapply(frut_par$fruit_per_plant[frut_par$genotype=="sweden"], frut_par$siteyear[frut_par$genotype=="sweden"], sd, na.rm=T),
    tapply(frut_par$fruit_per_seed [frut_par$genotype=="sweden"], frut_par$siteyear[frut_par$genotype=="sweden"], sd, na.rm=T),
    tapply(frut_par$survival       [frut_par$genotype=="sweden"], frut_par$siteyear[frut_par$genotype=="sweden"], sd, na.rm=T)),
  sw_n= c(
    tapply(!is.na(seed_par$mean_mass_ug   [seed_par$genotype=="sweden"]), seed_par$siteyear[seed_par$genotype=="sweden"], sum),
    tapply(!is.na(seed_par$nseeds         [seed_par$genotype=="sweden"]), seed_par$siteyear[seed_par$genotype=="sweden"], sum),
    tapply(!is.na(frut_par$fruit_per_plant[frut_par$genotype=="sweden"]), frut_par$siteyear[frut_par$genotype=="sweden"], sum),
    tapply(!is.na(frut_par$fruit_per_seed [frut_par$genotype=="sweden"]), frut_par$siteyear[frut_par$genotype=="sweden"], sum),
    tapply(!is.na(frut_par$survival       [frut_par$genotype=="sweden"]), frut_par$siteyear[frut_par$genotype=="sweden"], sum))
)
# correct SD to SE
parents$it_se <- parents$it_se / sqrt(parents$it_n -1)
parents$sw_se <- parents$sw_se / sqrt(parents$sw_n -1)

# insert mean seed numbers for each line in each site-year combination.
ix <- which(frut_par$genotype=="italy")
frut_par$seed_numb[ix] <- parents$it_parent[match(frut_par$siteyear, paste(parents$site, parents$year))+4][ix]
sx <- which(frut_par$genotype=="sweden")
frut_par$seed_numb[sx] <- parents$sw_parent[match(frut_par$siteyear, paste(parents$site, parents$year))+4][sx]
# Estimates of total fecundity and total fitness.
frut_par$tofu <- frut_par$fruit_per_plant * frut_par$seed_numb
frut_par$tfit <- frut_par$fruit_per_seed  * frut_par$seed_numb

# put together a data frame of data on seeds per reproductive plant.
tofu <- data.frame(
  variable = rep("total_fecundity", 4),
  site     = rep(c("italy", "sweden"), each=2),
  year     = rep(2010:2011, 2)
)
tofu$it_parent <- tapply(frut_par$tofu[frut_par$genotype=="italy"],  frut_par$siteyear[frut_par$genotype=="italy"],  mean, na.rm=T)
tofu$it_se     <- tapply(frut_par$tofu[frut_par$genotype=="italy"], frut_par$siteyear[frut_par$genotype=="italy"], sd, na.rm=T)/
  sqrt(tapply(!is.na(frut_par$tofu[frut_par$genotype=="italy"]),  frut_par$siteyear[frut_par$genotype=="italy"],  sum)-1)
tofu$it_n      <- tapply(!is.na(frut_par$tofu[frut_par$genotype=="italy"]),  frut_par$siteyear[frut_par$genotype=="italy"],  sum)
tofu$sw_parent <- tapply(frut_par$tofu[frut_par$genotype=="sweden"], frut_par$siteyear[frut_par$genotype=="sweden"], mean, na.rm=T)
tofu$sw_se     <- tapply(frut_par$tofu[frut_par$genotype=="sweden"], frut_par$siteyear[frut_par$genotype=="sweden"], sd, na.rm=T)/
  sqrt(tapply(!is.na(frut_par$tofu[frut_par$genotype=="sweden"]),  frut_par$siteyear[frut_par$genotype=="sweden"],  sum)-1)
tofu$sw_n      <- tapply(!is.na(frut_par$tofu[frut_par$genotype=="sweden"]),  frut_par$siteyear[frut_par$genotype=="sweden"],  sum)

# data frame on seeds per planted seed.
tfit <- data.frame(
  variable = rep("total_fitness", 4),
  site     = rep(c("italy", "sweden"), each=2),
  year     = rep(2010:2011, 2)
)
tfit$it_parent <- tapply(frut_par$tfit[frut_par$genotype=="italy"],  frut_par$siteyear[frut_par$genotype=="italy"],  mean, na.rm=T)
tfit$it_se     <- tapply(frut_par$tfit[frut_par$genotype=="italy"], frut_par$siteyear[frut_par$genotype=="italy"], sd, na.rm=T)/
  sqrt(tapply(!is.na(frut_par$tfit[frut_par$genotype=="italy"]),  frut_par$siteyear[frut_par$genotype=="italy"],  sum)-1)
tfit$it_n      <- tapply(!is.na(frut_par$tfit[frut_par$genotype=="italy"]),  frut_par$siteyear[frut_par$genotype=="italy"],  sum)
tfit$sw_parent <- tapply(frut_par$tfit[frut_par$genotype=="sweden"], frut_par$siteyear[frut_par$genotype=="sweden"], mean, na.rm=T)
tfit$sw_se     <- tapply(frut_par$tfit[frut_par$genotype=="sweden"], frut_par$siteyear[frut_par$genotype=="sweden"], sd, na.rm=T)/
  sqrt(tapply(!is.na(frut_par$tfit[frut_par$genotype=="sweden"]),  frut_par$siteyear[frut_par$genotype=="sweden"],  sum)-1)
tfit$sw_n      <- tapply(!is.na(frut_par$tfit[frut_par$genotype=="sweden"]),  frut_par$siteyear[frut_par$genotype=="sweden"],  sum)

# bind table together
parents <- rbind(parents, tofu, tfit)
# difference between parental ecotypes
parents$delta <- abs(parents$it_parent - parents$sw_parent)

# FITNESS DIFFERENCES AMONG PARENTS
maxvals <- apply(parents[,c(4,7)], 1, max) # get the fittest phenotype in each site-year combination.
# relative fitnesses of the Italian and Swedish parents.
parents$it_relw <- parents$it_parent / maxvals
parents$sw_relw <- parents$sw_parent / maxvals
# selective advantage of the local morph
parents$s[parents$site == 'italy']  <- parents$it_relw[parents$site == 'italy']  - parents$sw_relw[parents$site == 'italy']
parents$s[parents$site == 'sweden'] <- parents$sw_relw[parents$site == 'sweden'] - parents$it_relw[parents$site == 'sweden']

# get the differences between parental means in each expeirment.
parent_diffs <- data.frame(
  it2010 = parents$delta[parents$site == "italy"  & parents$year == "2010"],
  it2011 = parents$delta[parents$site == "italy"  & parents$year == "2011"],
  sw2010 = parents$delta[parents$site == "sweden" & parents$year == "2010"],
  sw2011 = parents$delta[parents$site == "sweden" & parents$year == "2011"])
parent_diffs <- parent_diffs[c(1,3,2,6,4,7,5),]
parent_diffs$trait <- c("mass", "seed", "frut", "ffit", "surv", "tofu", "tfit")[c(1,3,2,6,4,7,5)]

rm(ix,sx)
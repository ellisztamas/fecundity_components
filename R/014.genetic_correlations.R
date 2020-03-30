# RIL phenotypes
ril_numb <- read.csv('data_derived/RIL_seeds_per_fruit.csv')
ril_mass <- read.csv('data_derived/RIL_seed_mass.csv')
ril_frut <- read.csv("data_derived/RIL_fruits_per_adult.csv")
ril_tofu <- read.csv('data_derived/RIL_seeds_per_adult.csv')
ril_surv <- read.csv('data_derived/RIL_survival.csv')

cortests <- list(
  # tofu vs survival
  cor.test(ril_tofu$it2010, ril_surv$it2010, method = 'p'),
  cor.test(ril_tofu$it2011, ril_surv$it2011, method = 'p'),
  cor.test(ril_tofu$sw2010, ril_surv$sw2010, method = 'p'),
  cor.test(ril_tofu$sw2011, ril_surv$sw2011, method = 'p'),
  # frut vs seeds\fruit
  cor.test(ril_frut$it2010, ril_numb$it2010, method = 'p'),
  cor.test(ril_frut$it2011, ril_numb$it2011, method = 'p'),
  cor.test(ril_frut$sw2010, ril_numb$sw2010, method = 'p'),
  cor.test(ril_frut$sw2011, ril_numb$sw2011, method = 'p'),
  # tofu vs seed mass
  cor.test(ril_tofu$it2010, ril_mass$it2010, method = 'p'),
  cor.test(ril_tofu$it2011, ril_mass$it2011, method = 'p'),
  cor.test(ril_tofu$sw2010, ril_mass$sw2010, method = 'p'),
  cor.test(ril_tofu$sw2011, ril_mass$sw2011, method = 'p'))

cor_obs <- sapply(cortests, function(x) x$estimate)


nreps <- 1000
cor_bootstraps <- matrix(NA, nrow=nreps, ncol=12)

for(r in 1:nreps){
  ix <- sample(1:nrow(ril_frut), replace = T)
  # Correlations on resampled data
  cortests_boot <- list(
    # tofu vs survival
    cor.test(ril_tofu$it2010[ix], ril_surv$it2010[ix], method = 'p'),
    cor.test(ril_tofu$it2011[ix], ril_surv$it2011[ix], method = 'p'),
    cor.test(ril_tofu$sw2010[ix], ril_surv$sw2010[ix], method = 'p'),
    cor.test(ril_tofu$sw2011[ix], ril_surv$sw2011[ix], method = 'p'),
    # frut vs seed number
    cor.test(ril_frut$it2010[ix], ril_numb$it2010[ix], method = 'p'),
    cor.test(ril_frut$it2011[ix], ril_numb$it2011[ix], method = 'p'),
    cor.test(ril_frut$sw2010[ix], ril_numb$sw2010[ix], method = 'p'),
    cor.test(ril_frut$sw2011[ix], ril_numb$sw2011[ix], method = 'p'),
    # todu vs seed mass
    cor.test(ril_tofu$it2010[ix], ril_mass$it2010[ix], method = 'p'),
    cor.test(ril_tofu$it2011[ix], ril_mass$it2011[ix], method = 'p'),
    cor.test(ril_tofu$sw2010[ix], ril_mass$sw2010[ix], method = 'p'),
    cor.test(ril_tofu$sw2011[ix], ril_mass$sw2011[ix], method = 'p'))
  cor_bootstraps[r,] <- sapply(cortests_boot, function(x) x$estimate)
}

rhos <- data.frame(
  traits= rep(c("tofu_vs_surv", "frut_vs_numb", "tofu_vs_mass"), each=4),
  expt  = rep(colnames(ril_mass[-1]), 3),
  n     = sapply(cortests, function(x) x$parameter) + 2,
  obs   = sapply(cortests, function(x) x$estimate),
  p     = sapply(cortests, function(x) x$p.value),
  lower = apply(cor_bootstraps, 2, quantile, 0.025),
  upper = apply(cor_bootstraps, 2, quantile, 0.975)
)

write.csv(rhos, 'output/genetic_correlations.csv', row.names = F)

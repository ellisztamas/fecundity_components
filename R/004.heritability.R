#' Tom Ellis
#' 
#' Script to estimate broad-sense heritabilities for survival, fruits/RP,
#' seeds/fruit, seed mass and fruits/seedling.

library(lme4)

# Import and subset data
ind_data <- read.csv("data_raw/RILFieldExptItalySweden_SeedSizeNumber_SixYrs_398RIL Excl Edge.csv")
ind_data <- ind_data[ind_data$Year %in% 2010:2011,]
ind_data$Tray <- as.factor(ind_data$Tray)
# Add a column for survival
ind_data$survival <- ifelse(ind_data$NFruits > 0, 1, 0)
# Vector of traits of interest
traits <- c("survival", "Fecundity","X..Seeds..good.", "NFruits", "MSeedMass.ug.")
# Vector of site_year identifiers.
ind_data$site_year <- paste(ind_data$Site, ind_data$Year)
site_year <- unique(ind_data$site_year)

# Function to get broad-sense heritability from a mixed model.
herit <- function(model){
  vc <- as.data.frame(VarCorr(model))$vcov
  vc[1] / sum(vc)
}

# Observed heritabilities
# Empty df to store values
H2 <- list(
  obs = data.frame(matrix(NA, nrow = length(traits), ncol = length(site_year))),
  upper = data.frame(matrix(NA, nrow = length(traits), ncol = length(site_year))),
  lower = data.frame(matrix(NA, nrow = length(traits), ncol = length(site_year)))
)
# Label columns and rows
colnames(H2$obs)  <- colnames(H2$upper)  <- colnames(H2$lower)  <- site_year
row.names(H2$obs) <- row.names(H2$upper) <- row.names(H2$lower) <- traits

nreps <- 200
t0 <- proc.time()
for(t in traits){
  cat("Running",nreps,"bootstraps for the heritability of", t,"in..." )
  
  model <- as.formula(paste(t, "~ Tray + (1| Genotype2)"))
  for(sy in site_year){
    cat(sy, ", ", sep="")
    # Fit the model
    m <- lmer(model,
              data   = ind_data,
              subset = site_year == sy)
    # Pull out variance components and get H2.
    vc <- as.data.frame(VarCorr(m))$vcov
    H2$obs[t, sy] <- vc[1] / sum(vc)
    # Parametric bootstrapping from the model
    para_boot <- bootMer(m, herit, nsim = nreps)
    H2$lower[t, sy] <- quantile(para_boot$t, 0.025)
    H2$upper[t, sy] <- quantile(para_boot$t, 0.975)
  }
  cat("Done.\n")
}
paste("Parametric bootstraps completed after", round((proc.time()-t0)[3] / 60, 2), "minutes.")

expt <- split(ind_data, ind_data$site_year)
H2$n <- sapply(expt, function(x) colSums(!is.na(x[, traits])))
H2$n <- as.data.frame(H2$n)

saveRDS(H2, file = "output/heritabilities.rds")

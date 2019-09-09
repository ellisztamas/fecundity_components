# Data formatting
source("R/parents.R")
source("R/rqtl_files.R")
# Analyses of phenotypes
source("R/bootstrap_selection_coefficients.R")
source("R/heritability.R")
source("R/mass_wilcox_tests.R")
# QTL mapping
source("R/mapping_mass.R")# 1. seed mass
source("R/mapping_frut.R")# 2. fruits/plant
source("R/mapping_seed.R")# 3. seeds/fruit
source("R/mapping_tofu.R")# 4. Total fecundity (estimated seeds/plant = frut*seed)
source("R/mapping_ffit.R")# 5. Fruits per planted seedling
source("R/mapping_tfit.R")# 6. Seeds per planted seedling; fitness incorporating seed number.
source("R/mapping_surv.R")# 7. Survival
source("R/format_QTL_models.R")

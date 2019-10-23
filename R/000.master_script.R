# Data formatting
source("R/001.parents.R")
source("R/002.rqtl_files.R")
# Analyses of phenotypes
source("R/003.bootstrap_selection_coefficients.R")
source("R/004.heritability.R")
source("R/005.mass_wilcox_tests.R")
# QTL mapping
source("R/006.mapping_mass.R")# 1. seed mass
source("R/007.mapping_frut.R")# 2. fruits/plant
source("R/008.mapping_seed.R")# 3. seeds/fruit
source("R/009.mapping_tofu.R")# 4. Total fecundity (estimated seeds/plant = frut*seed)
source("R/010.mapping_ffit.R")# 5. Fruits per planted seedling
source("R/011.mapping_tfit.R")# 6. Seeds per planted seedling; fitness incorporating seed number.
source("R/012.mapping_surv.R")# 7. Survival to reproduction
# Stitch QTL results together and cluster them.
source("R/013.format_QTL_models.R")

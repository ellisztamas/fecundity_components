setwd("..")

library(qtl)



sweden_italy_cross <- read.cross(
  "csvs",
  dir = "",
  genotypes = c("a", "b"),
  genfile = "data_files/geno_SurvFec.csv",
  phefile = "data_files/pheno_SurvFec_wQN.csv",
  na.strings = c("NA", "-"),
  crosstype = "riself"
)

sweden_italy.calc <- calc.genoprob(
  sweden_italy_cross,
  step = 2,
  error.prob = .0001,
  map.function = "kosambi",
  stepwidth = "max",
  off.end = 0
)

set.seed(9274)

cnms <-
  c(
    "QN_F10_IT_FecFrt_noedge_mean",
    "QN_F11_IT_FecFrt_noedge_mean",
    "QN_F10_SW_FecFrt_noedge_mean",
    "QN_F11_SW_FecFrt_noedge_mean"
  )
qtlmodels <- lapply(1:4, function(x)
  NULL)
names(qtlmodels) <- c("it2010", "it2011", "sw2010", "sw2011")

for (c in 1:4) {
  EQPC_perm <- scantwo(
    sweden_italy.calc,
    pheno.col = cnms[c],
    model = "normal",
    method = "hk",
    n.perm = 10000,
    n.cluster = 16,
    addcovar = NULL,
    intcovar = NULL,
    weights = NULL,
    clean.output = FALSE,
    clean.nmar = FALSE,
    clean.distance = FALSE,
    incl.markers = TRUE,
    assumeCondIndep = FALSE,
    verbose = FALSE
  )
  
  EQPC_pen1.05 <- calc.penalties(EQPC_perm, alpha = .05)
  
  EPI_HEAVY_stepwise_EQPC <-
    stepwiseqtl(
      sweden_italy.calc,
      method = "hk",
      model = "normal",
      pheno.col = cnms[c],
      penalties = EQPC_pen1.05[1:2],
      max.qtl = 25,
      covar = NULL,
      scan.pairs = FALSE,
      additive.only = F,
      keeplodprofile = TRUE,
      keeptrace = TRUE,
      refine.locations = TRUE,
      verbose = TRUE
    )
}
saveRDS(qtlmodels, file = 'output/oakley_stepwise.rds')
perform_mapping <- function(ix, nclusters, nperms){
  path <- paste("data_files/rqtl_", ix, ".csv", sep="")
  cat("Stepwise QTL matching for R/QTL file", path,".\n")
  cat("Started", format(Sys.time(),usetz = TRUE), "\n\n")
  
  t0 <- proc.time() # save the system time.
  ###############
  # Import data #
  ###############
    library("qtl",  lib.loc='~/R/x86_64-redhat-linux-gnu-library/3.5')
    library("snow", lib.loc='~/R/x86_64-redhat-linux-gnu-library/3.5')
  # create r/qtl object from the file created above.
  cross <- read.cross(
    "csv",
    "",
    path,
    na.strings = NA,
    genotypes = c("a", "b"),
    alleles = c("It", "Sw"),
    crosstype="riself")
  # calculate the genotype probabilities
  cross <- calc.genoprob(
    cross,
    step = 2,
    error.prob = .0001,
    map.function = "kosambi",
    stepwidth = "max",
    off.end = 0
  )
  
  # perform quantile normal transformations on each variable.
  # cat("Performing quantile normal transformations.\n")
  # source('R/quantnorm.R')
  # for(c in 2:5) cross$pheno[,c] <- quantnorm(cross$pheno[,c])
  
  ################
  # Permutations #
  ################
  # Perumtations for seed fruter.
  cat("Starting permutations. ", nperms, "permutations will be performed.\n")
  # Empty list to store permutations for each year.
  permutations <- list(it2010=NULL, it2011=NULL, sw2010=NULL, sw2011=NULL)
  # loop across each site-year combination and perform permutations for each.
  for(y in 1:length(permutations)){
    t1 <- proc.time()
    permutations[[y]] <- 
      scantwo(
        cross,
        pheno.col = y + 1,
        model = "normal",
        method = "hk",
        n.perm = nperms,
        n.cluster = nclusters,
        verbose = F,
        clean.nmar = FALSE,
        clean.distance = FALSE,
        incl.markers = TRUE
      )
    cat("Permutations completed for", names(permutations)[y], "in", round(((proc.time()-t1)[3] /60),2),"minutes.\n")
  }
  # store permutations for later
  saveRDS(permutations, file= paste("output/permutations_", ix, ".rds", sep=""))
  
  #######################
  # Calculate penalties #
  #######################
  cat("Calculating penalties.\n")
  # Empty list to store penalties for each site-year combination.
  penalties <- list(it2010=NULL, it2011=NULL, sw2010=NULL, sw2011=NULL)
  # get penalties from the permutations
  for(y in 1:length(penalties)) penalties[[y]] <- calc.penalties(permutations[[y]], alpha = 0.05)
  saveRDS(penalties, file= paste("output/penalties_", ix, ".rds", sep=""))
  
  ###################
  # Stepwise models #
  ###################
  # perform quantile normal transformations on each variable.
  cat("Performing quantile normal transformations.\n")
  source('R/quantnorm.R')
  for(c in 2:5) cross$pheno[,c] <- quantnorm(cross$pheno[,c])
  
  cat("Fitting stepwise models.\n")
  # Empty list to store models for each site-year combination.
  stepwise <- list(it2010=NULL, it2011=NULL, sw2010=NULL, sw2011=NULL)
  # loop across each site-year combination and perform QTL search for each.
  for(y in 1:length(stepwise)){
    t2 <- proc.time()
    stepwise[[y]] <-
      stepwiseqtl(
        cross,
        pheno.col = y + 1,
        penalties = penalties[[y]],
        method = "hk",
        model = "normal",
        max.qtl = 15,
        keeptrace = TRUE,
        verbose = FALSE
      )
    cat("QTL search completed for", names(stepwise)[y], "in",  round(((proc.time()-t2)[3] /60),2),"minutes\n")
  }
  
  # save the stepwise QTL models.
  saveRDS(stepwise, file=paste("output/stepwise_", ix, ".rds", sep=""))
  
  cat("Analyses completed in", round(((proc.time()-t2)[3] /3600),2),"hours.\n")
  
}
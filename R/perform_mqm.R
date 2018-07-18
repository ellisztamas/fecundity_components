perform_mqm <- function(ix, phenocol, nperms=10000, nclusters= 16){
  path <- paste("data_files/rqtl_", ix, ".csv", sep="")
  
  t0 <- proc.time() # save the system time.
  
  cat("Multiple-QTL-Mapping analysis for R/QTL file", path,"\n")
  cat("Started", format(Sys.time(),usetz = TRUE), "\n\n")
  
  ###############
  # Import data #
  ###############
  # create r/qtl object from the file created above.
  cross <-
    read.cross(
      "csv",
      "",
      path,
      na.strings = NA,
      genotypes = c("a", "b"),
      alleles = c("It", "Sw"),
      crosstype="riself"
    )
  
  # perform quantile normal transformations on each variable.
  cat("Performing quantile normal transformations.\n")
  source('R/quantnorm.R')
  for(c in 2:5) cross$pheno[,c] <- quantnorm(cross$pheno[,c])
  
  
  # Ideally one would use mqmaugment, but this cannot handle this cross object.
  # Instead use full.geno; this is jutsified by this comment by Karl Broman: https://groups.google.com/forum/#!topic/Rqtl-disc/qtlh7DoM-sk
  # cross <- fill.geno(cross, error.prob=.0001, map.function="kosambi")
  
  # Data augmentation
  augcross <- mqmaugment(cross)
  saveRDS(augcross, file=paste("output/augment_", ix, ".rds", sep=""))

  
  cat("Beginning multiple QTL mapping.\n")
  setcofactors <- mqmsetcofactors(augcross, each = 5)
  mqm <- lapply(1:length(cols), function(x) NULL)
  names(mqm) <- cols
  
  for(i in cols) mqm[[i]] <- mqmscan(augcross, setcofactors, pheno.col = i, n.clusters = nclusters)
  saveRDS(mqm, file=paste("output/mqm_", ix, ".rds", sep=""))
  
  cat("Beginning permutations.", nperms, "permutations will be performed on",nclusters,".\n")
  permutations <- lapply(1:length(cols), function(x) NULL)
  names(permutations) <- cols
  for(i in cols){
    permutations[[i]] <-
      mqmpermutation(
        augcross,
        scanfunction = mqmscan,
        cofactors = setcofactors,
        pheno.col = i,
        n.cluster = nclusters,
        n.perm = nperms,
        batchsize = nclusters
      )
    cat("Permutations completed for", i, "after", round(((proc.time()-t0)[3] /3600),2),"hours.\n")
  }
  # store permutations for later
  saveRDS(permutations, file= paste("output/permutations_", ix, ".rds", sep=""))
  
  cat("Exporting MQM results as qtl objects.\n")
  qtlmodels <- lapply(mqm, mqmgetmodel)
  # save the QTL models.
  saveRDS(qtlmodels, file=paste("output/qtlmodels_", ix, ".rds", sep=""))
  
  cat("Analyses completed at", format(Sys.time(),usetz = TRUE), "\n\n")
  cat("======================================================================\n\n\n")
}
  
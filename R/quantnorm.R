quantnorm<-function(x){
  # apply quantile normalisation to a vector of floats.
  # Note that this will suppress warnings about converting to NaN
  n=sum(!is.na(x),na.rm=T)
  x=rank(x)/(n+1)
  x=suppressWarnings(qnorm(x))
  x[is.infinite(x)]=NA
  #x[is.nan(x)]=NA
  return(x)
}



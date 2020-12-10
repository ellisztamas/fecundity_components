#' Tom Ellis
#' 
#' Perform Wilcox-signed-rank tests on differences in seed mass between 
#' parental ecotypes in ecah experiment.

parent_mass <- read.csv('data_raw/individual_parents_massnumber.csv')
# unique identifier for each fite year combination
parent_mass$siteyear <- paste(parent_mass$site, parent_mass$year,sep="_")
siteyear <- sort(unique(parent_mass$siteyear))

# empty lists to store tests
mass_wilcox <- vector('list', length(siteyear))
names(mass_wilcox) <- siteyear
# run tests
for(y in siteyear){ # seed number
  this_exp <- parent_mass[parent_mass$siteyear == y,]
  mass_wilcox[[y]] <- wilcox.test(mean_mass_ug ~ genotype, data=this_exp)
}

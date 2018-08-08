# Positive correlations among components of fitness boost estimates of local adaptation in *Arabidopsis thaliana*

This repository documents analysis associated with the manuscript "Positive correlations among components of fitness boost estimates of local adaptation in *Arabidopsis thaliana*" by Tom Ellis, Froukje Postma, Christopher Oakley and Jon Ågren.

## Table of contents

1. [Experimental set up](#experimental-set-up)
2. [Analysis workflow](#analysis-workflow)
3. [Dependencies](#Dependencies)
4. [Data](#data)
    1. [Parental lines](#parental-lines)
    2. [RIL data](#ril-data)
    3. [Derived data](#derived-data)
    4. [R/QTL objects](#rqtl-objects)
5. [To-do list](#to-do)

## Experimental set up
We compared fitness components in parental lines and a recombinant inbred line (RIL) population derived from those parents.
Plants were transplanted as seedlings in reciprocal transplant experiments at the native sites of the parental lines in Italy and Sweden in 2010 and 2011.

We consider a total of seven traits. These have four letter abbreviations to keep the code straightforward:

1. **mass**: Mean seed mass.
2. **seed**: Seed number per fruit, a component of fecundity
3. **frut**: Fruit number per reproductive plant, a component of fecundity.
4. **tofu**: Seed number per reproductive plant, proxy for overall fecundity (estimated from 2 and 3).
5. **ffit**: Fruit number per planted seedling, a proxy for fitness
6. **tfit**: Fruit number per planted seedling, a proxy for fitness that includes seed number (estimated from 2 and 5).
7. **surv**: Survival to reproduction.

Traits 3, 5 and 7 are take directly from [Ågren *et al.*, 2013](http://www.pnas.org/content/110/52/21077/).
Traits 4 and 6 are not directly observed, but estimated as the product of other variables.


## Analysis workflow

* Parental lines are formatted in `R/parents.R`
* Wilcox-tests on parental seed mass data are run in `R/mass_wilcox_tests.R`
* QTL mapping:
	1. `R/rqtl_files.R` formats RIL data into R/QTL files.
	2. `R/perform_mapping.R` is a generic function to import an R/QTL input file, run permutations, and fit `stepwiseqtl()` models to the data. Output is saved as as RDS file in `output` folder with a generic names like `stepwise_trait.rds`.
	3. For each trait there is a script `R/mapping_trait.R` that called `perform_mapping.R` on that trait.
	4. `R/format_qtl_models.R` runs `fitqtl()` on each `stepwiseqtl()` model, and clusters colocalising QTL into groups. It groups cross objects, QTL models, model fits and clusters for a single trait into one list of objects, which is saved as an RDS file.
* The manuscript text and code to create figures is found in `manuscript/fecundity_components.Rmd` (see also `manuscript/supporting_information.Rmd`).

## Dependencies
All analyses were done in R using RStudio. The following additional packages are required, which can be installed from CRAN using `install.packages()`:

* `magittr`
* `knitr`
* `vioplot`
* `kableExtra`
* `qtl`

In addition, plotting and clustering QTL requires the custom package `qtltools` written by Tom Ellis for this manuscript, hosted on [GitHub](https://github.com/ellisztamas/qtltools). The easiest way to install it is to use `devtools`:

```
library('devtools')
install_github('ellisztamas/qtltools')

```

## Data
### Parental lines

* **individual_parents_massnumber.csv**: Values for seed mass and seed number for individual parental plants.
* **individual_parents_nfruit.csv**: Values for fruit number per plant and per seedling for individual parental plants.

Column headers used:
	* site: Italian or Swedish site
	* year: Year the experiment was transplanted.
	* tray: Block ID.
	* row: Row within block.
	* column: column within block
	* pos: Well position within block.
	* genotype: Italian or Swedish parent
	* subline: Italian and Swedish parents are divided into four sublines each; this is ignored in analysis
	* fruit_per_plant: Number of fruits per reproductive plant
	* fruit_per_seed: Number of fruits per planted seedling
	* nseeds: Number of seeds per fruit
	* mean_mass_ug: Mean seed mass
	* total_mass_mg: Total mass of all seeds in the fruit

### RIL data

* **Agren2013_fitness_components.csv**: RIL mean survival (Surv) and fruits per reproductive plant (NFrFrPr) from [Ågren *et al.*, 2013](http://www.pnas.org/content/110/52/21077/)
* **Agren2013_LSMfitness.csv**: RIL least-square mean fruits per planted seedling (NFruit) from [Ågren *et al.*, 2013](http://www.pnas.org/content/110/52/21077/)
* **ril_genotypes.csv**: RIL genotypes from [Ågren *et al.*, 2013](http://www.pnas.org/content/110/52/21077/)
* **genetic_map_agren_etal_2013.csv**: Genetic map from [Ågren *et al.*, 2013](http://www.pnas.org/content/110/52/21077/), used for plotting QTL positions in the main text.

Phenotype data column headers show RIL name (**id**), and then a code that shows abbreviations used by Ågren *et al.*:

* year (F09, F10, F11)
* site (IT, SW)
* trait (NFrFrPr, Surv, NFruit)
* whether line means are based on arithmetic means (mean) or least-square means (lsm).
* noedge means edge plants from the last three rows of each tray were removed to avoid edge effects

Genotype data are in R/QTL format, showing RIL names, following by each locus. Lines 2 and 3 indicate chromosome number an linkage map position.

### Derived data

The script `R/rqtl_files.R` formats RIL data. This outputs .csv for each trait with RIL names in a common order for performing genetic correlations with five columns:

* **id**: RIL name.
* **it2010**: RIL means for Italy in 2010.
* **it2011**: RIL means for Italy in 2011.
* **sw2010**: RIL means for Sweden in 2010.
* **sw2011**: RIL means for Sweden in 2011.

It also outputs R/QTL files for mapping with the same columns, but with genotype information included.

### RQTL objects
R/QTL output is saved to `./output`. See [the section on workflow](#analysis-workflow) for more on how analyses were run. Because there are so many trait-year-site levels, I have bundled objects together where possible, and loop over them. This avoids copy-paste errors, but requires some explanation. 

There is one object of class `qtl` derived from the function `stepwiseqtl()` for each of the seven traits with obvious names like `stepwise_surv.rds` etc that follow the abbreviations [above](#experimental-set-up). Models for `frut`, `ffit` and `surv` are taken straight from [Ågren *et al.*, 2013](http://www.pnas.org/content/110/52/21077/). Each of these is a list of four `qtl` objects:

1. **it2010**: QTL detected in Italy in 2010
2. **it2011**: QTL detected in Italy in 2011
3. **sw2010**: QTL detected in Sweden in 2010
4. **sw2011**: QTL detected in Sweden in 2011

There are three other RDS objects that are nested lists. The first level has seven elements (one for each trait) and the second level has four elements (one for each site-year combination):

1. **qtl_stepwise_models.rds**: objects from `stepwiseqtl`.
2. **qtl_model_fits**: `fitqtl` objects derived from each stepwise model.
3. **qtl_clusters.rds**: clusters of colocalising QTL.

Note that no QTL were detected for survival to reproduction in Sweden in 2010.

## To do
Finalise bibtex file.
Reformat qtl figure.
Document workflow and variables.

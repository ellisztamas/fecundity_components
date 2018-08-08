# Positive correlations among components of fitness boost estimates of local adaptation in *Arabidopsis thaliana*

This repository documents analysis associated with the manuscript "Positive correlations among components of fitness boost estimates of local adaptation in *Arabidopsis thaliana*" by Tom Ellis, Froukje Postma, Christopher Oakley and Jon Ågren.



## Table of contents

## Traits measured
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

## Data
### Raw data
#### Parental lines

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

#### RIL data

* **Agren2013_fitness_components.csv**: RIL mean survival (Surv) and fruits per reproductive plant (NFrFrPr) from [Ågren *et al.*, 2013](http://www.pnas.org/content/110/52/21077/)
* **Agren2013_LSMfitness.csv**: RIL least-square mean fruits per planted seedling (NFruit) from [Ågren *et al.*, 2013](http://www.pnas.org/content/110/52/21077/)** [Ågren *et al.*, 2013](http://www.pnas.org/content/110/52/21077/)

Column headers show RIL name (**id**), and then a code that shows abbreviations used by Ågren *et al.*:

* year (F09, F10, F11)
* site (IT, SW)
* trait (NFrFrPr, Surv, NFruit)
* whether line means are based on arithmetic means (mean) or least-square means (lsm).
* noedge means edge plants from the last three rows of each tray were removed to avoid edge effects

### Derived data

### R/QTL objects

## Analysis workflow

* Parental lines are formatted in `R/parents.R`
* Wilcox-tests on parental seed mass data are run in `R/mass_wilcox_tests.R`



`R/rqtl_files.R`

## Dependencies

## To do
Finalise bibtex file.
Reformat qtl figure.
Document workflow and variables.

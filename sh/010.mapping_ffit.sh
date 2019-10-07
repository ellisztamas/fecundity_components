#!/usr/bin/env bash

# SLURM
# This job requests from SLURM to allocate 1 node.
#SBATCH --nodes=16
# On that node, it will run 4 tasks, each with 1 core and 1 GB of memory.
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1G
#SBATCH --time=0-24:00:00
#SBATCH --output=./output/rqtl_logs/010.mapping_ffit.log

# ENVIRONMENT #
module load r/3.5.1-foss-2018b

Rscript R/010.mapping_ffit.R

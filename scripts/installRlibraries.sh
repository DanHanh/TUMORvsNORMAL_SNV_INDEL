#!/bin/bash

#SBATCH --mem=16000
#SBATCH --output=R-%x_%j.stdout
#SBATCH --error=R-%x_%j.stderr
#SBATCH --cpus-per-task=8

module load R/4.1.0-foss-2021a

libPath='Rlibrary'

mkdir Rlibrary

Rscript ./scripts/InstallLibraries.R $libPath


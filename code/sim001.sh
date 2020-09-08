#!/bin/bash
#SBATCH -o slurm.out
#SBATCH -e slurm.err 
#SBATCH --mem=32G
#SBATCH -p herringlab
module load R/x86_64-pc-linux-gnu-library/3.6
R CMD BATCH sim001.R
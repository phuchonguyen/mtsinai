#!/bin/bash
#SBATCH -o slurm.out
#SBATCH -e slurm.err 
#SBATCH --account=herringlab
#SBATCH --mem=8G
#SBATCH -p herringlab,statdept-low,volfovskylab-low
export R_LIBS_USER=~/R/x86_64-pc-linux-gnu-library/3.6
module load R
R CMD BATCH sim001.R
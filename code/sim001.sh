#!/bin/bash
#SBATCH -o slurm.out
#SBATCH -e slurm.err 
#SBATCH --mem=8G
#SBATCH -p herringlab,statdept-low,volfovskylab-low
module load R/3.6
R CMD BATCH sim001.R
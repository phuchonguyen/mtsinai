#!/bin/bash
#SBATCH -o slurm.out
#SBATCH -e slurm.err 
#SBATCH -–mem=32G
#SBATCH –p herringlab
module load R/3.6.0
R CMD BATCH sim_long.R
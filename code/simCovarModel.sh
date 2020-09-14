#!/bin/bash
#SBATCH -o slurm_covar.out
#SBATCH -e slurm_covar.err 
#SBATCH --account=herringlab
#SBATCH --mem=8G
#SBATCH -p herringlab,statdept-low,volfovskylab-low
export R_LIBS_USER=~/R/x86_64-pc-linux-gnu-library/4.0
module load R
Rscript code/simCovarModel.R
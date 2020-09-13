#!/bin/bash
#SBATCH -o slurm.out
#SBATCH -e slurm.err 
#SBATCH --account=herringlab
#SBATCH --mem=8G
#SBATCH -p herringlab,statdept-low,volfovskylab-low
export R_LIBS_USER=~/R/x86_64-pc-linux-gnu-library/4.0
module load R
Rscript /work/phn5/mtsinai/mtsinai/code/sim001.R

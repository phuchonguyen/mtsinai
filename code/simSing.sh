#!/bin/bash
#SBATCH -o slurm.out
#SBATCH -e slurm.err 
#SBATCH --account=herringlab
#SBATCH -c32
#SBATCH -p herringlab,statdept-low,volfovskylab-low

singularity exec --bind samples mtsinai.sif ls -l /samples
singularity exec mtsinai.sif Rscript simSing.R
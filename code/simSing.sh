#!/bin/bash
#SBATCH -o slurm.out
#SBATCH -e slurm.err 
#SBATCH --account=herringlab
#SBATCH --mem=8G
#SBATCH -p herringlab,statdept-low,volfovskylab-low

singularity exec --bind code/samples code/mtsinai.sif ls -l code/samples
singularity run code/mtsinai.sif
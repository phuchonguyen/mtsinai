#!/bin/bash
#SBATCH -o slurm.out
#SBATCH -e slurm.err 
#SBATCH --account=herringlab
#SBATCH -c32
#SBATCH -p herringlab,statdept-low,volfovskylab-low

singularity exec --bind /work/phn5/mtsinai/mtsinai/code/samples mtsinai.sif ls samples
singularity exec mtsinai.sif Rscript simSing.R
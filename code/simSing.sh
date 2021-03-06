#!/bin/bash
#SBATCH -o slurm.out
#SBATCH -e slurm.err 
#SBATCH --account=herringlab
#SBATCH -c32
#SBATCH -p herringlab,statdept-low,volfovskylab-low
#SBATCH --mail-user=phn5@duke.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE

singularity exec --bind /hpc/home/phn5/mtsinai/code/samples mtsinai.sif ls samples
singularity run mtsinai.sif simSing.R
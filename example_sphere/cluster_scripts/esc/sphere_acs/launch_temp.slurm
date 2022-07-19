#!/bin/bash
#SBATCH --mail-user=caw97@cam.ac.uk
#SBATCH --job-name=tdda
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4G
##SBATCH --mail-type=START,FAIL,END
#SBATCH --output=slurm-%j.out  

module load anaconda3

/home/caw97/source_code/g-dda/source_code/ddscat
wait
python sphere_qabs.py
wait
/home/caw97/source_code/t-dda/source_code/Lattice_Diffusion /home/caw97/source_code/t-dda/myGreen.num_300 var.par tdda_input temp.out




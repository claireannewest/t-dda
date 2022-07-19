#!/bin/bash
#SBATCH --mail-user=caw97@cam.ac.uk
#SBATCH --job-name=tdda
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4G
##SBATCH --mail-type=START,FAIL,END
#SBATCH --output=slurm-%j.out

/home/caw97/source_code/g-dda/source_code/ddscat
wait
mv tdda_input_w000_ddscat.par tdda_input
wait
/home/caw97/source_code/t-dda/source_code/Lattice_Diffusion /home/caw97/source_code/t-dda/myGreen.num_300 var.par tdda_input temp.out




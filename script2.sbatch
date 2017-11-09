#!/bin/bash
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=sing5957@vandals.uidaho.edu

. /usr/modules/init/bash
module load gcc

cd $SLURM_SUBMIT_DIR

make

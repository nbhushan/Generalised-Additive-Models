#!/bin/bash
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --job-name=gam_final_p275030
#SBATCH --mem=120GB
#SBATCH --mail-user=n.bhushan@rug.nl
#SBATCH --mail-type=BEGIN,END,FAIL

pwd
module load R
cd $HOME/GAMM/
Rscript gamFinalscript.r

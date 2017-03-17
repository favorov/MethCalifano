#! /bin/bash
#SBATCH
#SBATCH --time=45:0:0
#SBATCH --partition=shared
#skel for marcc
#SBATCH --mem=40G
module load R/3.3.1
Rscript prepare.gw.noodles.C.R

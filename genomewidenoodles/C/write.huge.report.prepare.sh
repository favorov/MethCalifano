#! /bin/bash
#SBATCH
#SBATCH --time=45:0:0
#SBATCH --partition=shared
#skel for marcc
#SBATCH --mem=15G
module load R/3.3.1
Rscript write.huge.report.prepare.R

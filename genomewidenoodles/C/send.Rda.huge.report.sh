#! /bin/bash
#SBATCH
#SBATCH --time=45:0:0
#SBATCH --partition=shared
#skel for marcc
#SBATCH --mem=10G
scp huge.report.frame.Rda favorov@favorov.bioinfolab.net:/mnt/storage/favorov/reports

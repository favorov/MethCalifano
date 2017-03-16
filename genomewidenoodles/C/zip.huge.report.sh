#! /bin/bash
#SBATCH
#SBATCH --time=45:0:0
#SBATCH --partition=shared
#skel for marcc
#SBATCH --mem=10G
zip -9 report.HNSCC.C.noodles.complete.zip noodles.C.complete.annotaion.tsv

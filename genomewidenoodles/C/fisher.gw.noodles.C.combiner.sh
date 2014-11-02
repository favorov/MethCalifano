#/usr/bin/bash
workers=50
Rscript fisher.gw.noodles.C.R combiner $workers || exit 1

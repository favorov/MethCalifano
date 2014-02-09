Which cytobands are differentially methylated in the 2013-11-25_RO1_Batch1_ClinicalData_DeID data?
========================================================

## Alexander Favorov 
### 26 Dec 2013



We load the cytoband list (karyotype) from UCSC DAS server and we load our MACS from our /Clinical/2013-11-25_RO1_Batch1_ClinicalData_DeID.xls.csv file. 



Then, for each cytoband we calculate the net length of methylated regions (we read them from to MACS bedfiles) that ovelap with the cytoband in each of the samples.




Thus, for each interval (cytodand) we have a methylated length (coverage) for each sample).
Then, we calculate the Wicoxon p-value and anova test Pr(>F) for each cytoband's tumor and normal DM coverages.




And, finally, we make list of all the cytobands with one of the p-values that is <=0.05 and output its coords, p-values and the samplewise methylated regions coverage. 


```
## RangedData with 1 row and 3 value columns across 24 spaces
##      space               ranges |       id wilcoxon.result anova.result
##   <factor>            <IRanges> | <factor>       <numeric>    <numeric>
## 1     chr3 [44100001, 44200000] |   p21.32         0.03311      0.03323
```

```
## RangedData with 1 row and 17 value columns across 24 spaces
##      space               ranges |       id    X37941    X35667    X25441    X35675     X4931    X18732    X31124    X35663    X27281    X27279    X30553    X30596    X30577    X30580    X30018    X32913
##   <factor>            <IRanges> | <factor> <integer> <integer> <integer> <integer> <integer> <integer> <integer> <integer> <integer> <integer> <integer> <integer> <integer> <integer> <integer> <integer>
## 1     chr3 [44100001, 44200000] |   p21.32         0       241         0      1611       406         0       623       534       558         0       426       589      1055      1311      1610       893
```


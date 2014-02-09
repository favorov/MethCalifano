Trying to test whether chrY methylation is represented in female samples less than in male
========================================================

### Alexander Favorov 
### 19 Dec 2013



We will test three measures: the coverage of methylation MACS peaks, their score (summed over chromosome) and their product.



Here is coverage of Y-chromome of all the samples by the methylated peaks, the sum of all the scores and the product of length and scores for all the samples.


```
##    code DNAid tumor  male Peak length sum Peak score sum Product sum
## 1     1 37941  TRUE  TRUE          413742         101598   207614610
## 2     2 35667  TRUE FALSE          376294         102885   150700149
## 3     3 25441  TRUE  TRUE          313265          52486   143475969
## 4     4 35675  TRUE FALSE          413638         104126   147113371
## 5     5  4931  TRUE  TRUE          651050         162215   220521307
## 6     6 18732  TRUE  TRUE          564458         145679   230736378
## 7     7 31124  TRUE  TRUE          403700          91285   121400575
## 8     8 35663  TRUE  TRUE          413478         103125   132697472
## 9     9 27281  TRUE  TRUE          467948         128386   192180775
## 10   10 27279  TRUE  TRUE          430185          97185   107257870
## 11  101 30553 FALSE  TRUE          304730          62312    98240747
## 12  102 30596 FALSE  TRUE          763616         207376   366795930
## 13  103 30577 FALSE  TRUE          462897         114264   164465789
## 14  104 30580 FALSE  TRUE          518039         127544   187512987
## 15  105 30018 FALSE FALSE          239088          35748    32169528
## 16  106 32913 FALSE FALSE          392266          93407   155932321
```


Wilcoxon test of the methylation peaks coverage the Y-chromosome for males/females 


```
## 
## 	Wilcoxon rank sum test
## 
## data:  Y.chrom.stats$"Peak length sum"[Y.chrom.stats$male] and Y.chrom.stats$"Peak length sum"[!Y.chrom.stats$male]
## W = 40, p-value = 0.05824
## alternative hypothesis: true location shift is not equal to 0
```

Wilcoxon test of the scores sum of Y-chromosome methylation peaks for males/females 


```
## 
## 	Wilcoxon rank sum test
## 
## data:  Y.chrom.stats$"Peak score sum"[Y.chrom.stats$male] and Y.chrom.stats$"Peak score sum"[!Y.chrom.stats$male]
## W = 34, p-value = 0.2615
## alternative hypothesis: true location shift is not equal to 0
```


Wilcoxon test of the scores*length product sum Y-chromosome for males/females 

```
## 
## 	Wilcoxon rank sum test
## 
## data:  Y.chrom.stats$"Product sum"[Y.chrom.stats$male] and Y.chrom.stats$"Product sum"[!Y.chrom.stats$male]
## W = 33, p-value = 0.3165
## alternative hypothesis: true location shift is not equal to 0
```


### The best is just to calculate lenghts.

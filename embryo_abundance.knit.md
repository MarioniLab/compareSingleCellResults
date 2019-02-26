---
title: Differential abundance analyses with the embryo chimera data set
author:
- name: Aaron T. L. Lun
  affiliation: &CRUK Cancer Research UK Cambridge Institute, Li Ka Shing Centre, Robinson Way, Cambridge CB2 0RE, United Kingdom
- name: Jonathan A. Griffiths
  affiliation: *CRUK
date: "2019-02-26"
vignette: >
  %\VignetteIndexEntry{04. Embryo differential abundance}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
output:
  BiocStyle::html_document:
    toc_float: true
    titlecaps: false
bibliography: ref.bib
---

<!--
AL: asterisks below remove weird highlighting on my editor.
****
-->





# Introduction

We previously merged the samples together onto the same coordinate system to generate common clusters.
Each cluster represents a molecular phenotype (in terms of gene expression profiles), which in turn is a proxy for cell states or types.
Now, our aim is to identify differences between conditions for cells in each cluster.
This allows us to quantify the effects of losing _Tal1_ during early embryonic development.


```r
library(SingleCellExperiment)
sce <- readRDS("embryo_merged.rds")
```

The most obvious differences between conditions are those of changes in the per-cluster cell abundance.
This will reveal which cell types are depleted or enriched upon loss of _Tal1_ function.
Specifically, our aim is to detect significant changes in cell abundance between genotypes.
This is done using the table of cell counts per cluster per sample, as shown below.


```r
cluster.counts <- table(sce$cluster, sce$sample)
cluster.counts    
```

```
##     
##        1   2   3   4
##   1  388 338 389 423
##   2    0   1 833 873
##   3  511 472 301 308
##   4    0   0   0   0
##   5    0   1 157 199
##   6  138 114  47  47
##   7  180 151 116 115
##   8    0   0   0   0
##   9   85 109 152 152
##   10 393 338 162 174
##   11 104  82  76  63
##   12 259 284 157 179
##   13 203 178  61  81
##   14 297 253  96 107
##   15 152 138 136 128
##   16  89  79  34  41
##   17 131 118  94 113
##   18 155 153  52  43
##   19 182 152 136 130
##   20  78  61  42  27
##   21   0   1  28  30
##   22 138 155  43  39
##   23  96 102  18  24
##   24  49  49  20  25
##   25  31  31  22  22
##   26  12  13  14   8
##   27  10  10  15   4
##   28   6   5   9  15
```

# Setting up for statistical modelling 

Our differential abundance analysis will be performed using methods from the *[edgeR](https://bioconductor.org/packages/3.9/edgeR)* package [@robinson2010edgeR].
This uses a negative binomial generalized linear model (NB GLM) to handle overdispersed count data in experiments with limited replication.
In our case, we have biological variation with only two replicates per condition, so *[edgeR](https://bioconductor.org/packages/3.9/edgeR)* (or its contemporaries) is a natural choice for the analysis.
The same strategy is used to analyze cell count data in mass cytometry experiments [@lun2017testing].


```r
library(edgeR)
y.ab <- DGEList(cluster.counts)
y.ab
```

```
## An object of class "DGEList"
## $counts
##    
##       1   2   3   4
##   1 388 338 389 423
##   2   0   1 833 873
##   3 511 472 301 308
##   4   0   0   0   0
##   5   0   1 157 199
## 23 more rows ...
## 
## $samples
##   group lib.size norm.factors
## 1     1     3687            1
## 2     1     3388            1
## 3     1     3210            1
## 4     1     3370            1
```

Typical applications of *[edgeR](https://bioconductor.org/packages/3.9/edgeR)* (for differential _expression_ analyses) perform a normalization step with `calcNormFactors()` to remove composition biases.
This requires the assumption that most of the input features (i.e., genes) are not differentially expressed between conditions.
While this assumption is reasonable for most gene expression data sets, we cannot say the same thing for cell abundance -
it is generally too strong to assume that most cell types do not change in abundance between conditions.
Thus, we will _not_ run `calcNormFactors()` here - at least, not at first (see below).
Any changes we detect between conditions will represent differences in the proportion of cells in each cluster.

Another typical step in *[edgeR](https://bioconductor.org/packages/3.9/edgeR)* is to filter out low-abundance features.
This aims to reduce computational work, improve the accuracy of the trend fitting and reduce the severity of the multiple testing correction.
However, our features (i.e., clusters) are data-defined and - depending on the clustering algorithm - will generally not be of low-abundance^[Otherwise there would not have been enough evidence to define it in the first place!].
Thus, this step can also be skipped unless the clustering algorithm tends to output many small clusters.

# Modelling the biological variation

We will use the quasi-likelihood (QL) framework [@chen2016reads], which involves estimating both the NB dispersion and the QL dispersion for each cluster. 
We set up the design matrix to block on the batch differences between replicates.


```r
genotype <- c("ko", "ko", "wt", "wt")
batch <- factor(c(1,2,1,2))
design <- model.matrix(~0 + genotype + batch)
design
```

```
##   genotypeko genotypewt batch2
## 1          1          0      0
## 2          1          0      1
## 3          0          1      0
## 4          0          1      1
## attr(,"assign")
## [1] 1 1 2
## attr(,"contrasts")
## attr(,"contrasts")$genotype
## [1] "contr.treatment"
## 
## attr(,"contrasts")$batch
## [1] "contr.treatment"
```

We use the `estimateDisp()` function to estimate the NB dipersion for each cluster (Figure \@ref(fig:abplotbcv)).
The role of the NB dispersion is to model the mean-variance trend, which is not easily accommodated by QL dispersions alone^[Due to the quadratic nature of the NB mean-variance trend.].


```r
y.ab <- estimateDisp(y.ab, design)
summary(y.ab$trended.dispersion)
```

```
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
## 9.766e-05 9.766e-05 9.766e-05 9.766e-05 9.766e-05 9.766e-05
```

```r
plotBCV(y.ab, cex=1)
```

<div class="figure">
<img src="embryo_abundance_files/figure-html/abplotbcv-1.png" alt="Biological coefficient of variation (BCV) for each cluster with respect to its average abundance. BCVs are defined as the square root of the NB dispersion. Trended and common dispersion estimates are shown in blue and red, respectively." width="100%" />
<p class="caption">(\#fig:abplotbcv)Biological coefficient of variation (BCV) for each cluster with respect to its average abundance. BCVs are defined as the square root of the NB dispersion. Trended and common dispersion estimates are shown in blue and red, respectively.</p>
</div>

The QL dispersion models the uncertainty and variability of the per-cluster variance (Figure \@ref(fig:.
This is, conversely, not well handled by the NB dispersions, so the two dispersion types complement each other in the final analysis.
We use `glmQLFit()` to fit a GLM to the counts for each cluster and estimate the QL dispersion from the GLM deviance. 
We set `robust=TRUE` to avoid distortions from highly variable clusters [@phipson2016robust].
We also turn off the abundance trend as there are not enough features for a stable trend fit, nor is there evidence for a strong trend (Figure \@ref(fig:abplotql)).


```r
fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
summary(fit.ab$var.prior)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   1.463   1.463   1.463   1.463   1.463   1.463
```

```r
summary(fit.ab$df.prior)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##     Inf     Inf     Inf     Inf     Inf     Inf
```

```r
plotQLDisp(fit.ab, cex=1)
```

<div class="figure">
<img src="embryo_abundance_files/figure-html/abplotql-1.png" alt="QL dispersion estimates for each cluster with respect to its average abundance. Quarter-root values of the raw estimates are shown in black while the shrunken estimates are shown in red. Shrinkage is performed towards the common dispersion in blue." width="100%" />
<p class="caption">(\#fig:abplotql)QL dispersion estimates for each cluster with respect to its average abundance. Quarter-root values of the raw estimates are shown in black while the shrunken estimates are shown in red. Shrinkage is performed towards the common dispersion in blue.</p>
</div>

# Testing for significant differences

We test for differences in abundance between genotypes using `glmQLFTest()`.
Many clusters are differentially abundant at a false discovery rate (FDR) of 5%.


```r
con <- makeContrasts(genotypeko - genotypewt, levels=design)
res <- glmQLFTest(fit.ab, contrast=con)
summary(decideTests(res))
```

```
##        1*genotypeko -1*genotypewt
## Down                            5
## NotSig                         10
## Up                             13
```

The top few clusters are strongly depleted in the KO condition.


```r
tab <- topTags(res, n=nrow(res))
head(tab$table, 10)
```

```
##         logFC   logCPM          F       PValue          FDR
## 2  -10.508916 16.93894 1640.30343 2.119701e-26 5.935164e-25
## 5   -8.239736 14.70494  342.55738 3.066234e-17 4.292727e-16
## 14   1.330383 15.76405   94.92404 1.705040e-10 1.591371e-09
## 10   1.014395 16.26237   81.26162 9.027742e-10 6.319419e-09
## 22   1.735478 14.77373   75.76610 1.875017e-09 1.050010e-08
## 18   1.589765 14.87550   70.05746 4.181802e-09 1.951508e-08
## 23   2.133777 14.14676   67.84085 5.783914e-09 2.313566e-08
## 13   1.319240 15.24500   65.37685 8.370438e-09 2.929653e-08
## 21  -5.633614 12.25890   51.66209 8.010719e-08 2.492224e-07
## 6    1.311483 14.66014   42.96199 4.171633e-07 1.168057e-06
```

Further examination indicates that they are derived from the erythroid lineage, with high expression of hemoglobins and related genes (Figure \@ref(abplothemo)).
This is consistent with the expected function of _Tal1_ in promoting erythroid differentiation.


```r
library(scater)
plotExpression(sce, x="cluster", "Hba-x", colour_by="tomato")
```

<div class="figure">
<img src="embryo_abundance_files/figure-html/abplothemo-1.png" alt="Distribution of _Hba-x_ expression in each cluster. Each point represents a cell and is coloured by the genotype." width="100%"  class="widefigure" />
<p class="caption">(\#fig:abplothemo)Distribution of _Hba-x_ expression in each cluster. Each point represents a cell and is coloured by the genotype.</p>
</div>



A more subtle example is that of cluster 13, which increases in abundance in the KO condition. 
Based on the expression of _Plac1_ and _Lgals1_ (Figure \@ref(subtleheat)), this cluster probably contains placental cells or their precursors.
From this, it is tempting to speculate that the loss of _Tal1_ diverts cells into lineages that they would not normally differentiate towards.


```r
chosen <- "13"
tab$table[chosen,]
```

```
##      logFC logCPM        F       PValue          FDR
## 13 1.31924 15.245 65.37685 8.370438e-09 2.929653e-08
```

```r
markers <- readRDS("embryo_markers.rds")
marker.set <- markers[[chosen]]
head(rownames(marker.set), 20)
```

```
##  [1] "H19"     "Krt8"    "Lgals1"  "Acp5"    "Plac1"   "Tmsb10"  "Prrx2"  
##  [8] "Pmp22"   "Spin2c"  "Gypa"    "Phlda2"  "Cdkn1c"  "Tmsb4x"  "Pitx1"  
## [15] "Hbb-bh1" "Ptn"     "Hoxaas3" "Tuba1a"  "Tspan32" "Ppic"
```

```r
logFCs <- as.matrix(marker.set[1:50,-(1:3)])
colnames(logFCs) <- sub("logFC.", "", colnames(logFCs))

library(pheatmap)
max.lfc <- max(abs(range(logFCs)))
pheatmap(logFCs, breaks=seq(-5, 5, length.out=101))
```

<div class="figure">
<img src="embryo_abundance_files/figure-html/subtleheat-1.png" alt="Heatmap of the log-fold changes for the top 50 genes expressed in cluster 13 compared to all other clusters. Each column represents another cluster." width="100%" />
<p class="caption">(\#fig:subtleheat)Heatmap of the log-fold changes for the top 50 genes expressed in cluster 13 compared to all other clusters. Each column represents another cluster.</p>
</div>



# Handling composition effects

## Assuming most clusters do not change

This data set contains a large change in abundance as the erythroid lineage is completely lost in the KO condition.
An interesting question, then, is whether the _increases_ in abundance in the KO condition are (i) caused by concomitant compositional effects from the loss of the erythroid lineage,
or (ii) they are a separate result of another biological process involving _Tal1_.
This is not a question that is generally answerable without additional assumptions, one of which is that most clusters are not differentially abundant.
Under this assumption, we use `calcNormFactors()` to compute normalization factors^[Not the same as size factors!] for each sample.


```r
y.ab2 <- calcNormFactors(y.ab)
y.ab2$samples
```

```
##   group lib.size norm.factors
## 1     1     3687    1.2534398
## 2     1     3388    1.1487510
## 3     1     3210    0.8328802
## 4     1     3370    0.8338503
```

We then proceed with the remainder of the *[edgeR](https://bioconductor.org/packages/3.9/edgeR)* analysis, shown below in condensed format.
The top hits do not change but many of the positive log-fold changes are shrunk towards zero.
This suggests that, while still significant, the increases in abundance observed in our original analysis were magnified by composition effects 
(conditional on the correctness of our above assumption).


```r
y.ab2 <- estimateDisp(y.ab2, design)
fit.ab2 <- glmQLFit(y.ab2, design, robust=TRUE, abundance.trend=FALSE)
res2 <- glmQLFTest(fit.ab2, contrast=con)
topTags(res2, n=10)
```

```
## Coefficient:  1*genotypeko -1*genotypewt 
##          logFC   logCPM          F       PValue        FDR
## 1   -0.7931553 16.75459 4064.89135 0.0009077037 0.01686157
## 2  -10.9882202 16.91665 2938.94263 0.0016682395 0.01686157
## 23   1.6165238 14.11056 1749.64395 0.0018065967 0.01686157
## 10   0.4882554 16.22599  974.95264 0.0029115669 0.02038097
## 6    0.7846628 14.62405  550.12642 0.0051036946 0.02171596
## 5   -8.7127229 14.67216  644.82662 0.0053394741 0.02171596
## 14   0.8046856 15.72744  630.98719 0.0054289905 0.02171596
## 21  -6.1097159 12.22372   94.03414 0.0232307322 0.07906466
## 16   0.5342232 14.12833   83.53067 0.0254136407 0.07906466
## 17  -0.3598262 15.01553   63.98595 0.0310894559 0.08705048
```

## Removing the offending clusters

Another approach to avoiding composition effects is to repeat the analysis after removing differentially abundant clusters with many cells.
This provides a clearer picture of the changes in abundance among the remaining clusters.
In this case, we would like to remove the blood-related clusters with strong hemoglobin expression.


```r
offenders <- c("2", "5", "4", "8") # TODO: replace with a more automated mechanism.
y.ab3 <- DGEList(cluster.counts[setdiff(rownames(cluster.counts), offenders),])
y.ab3
```

```
## An object of class "DGEList"
## $counts
##    
##       1   2   3   4
##   1 388 338 389 423
##   3 511 472 301 308
##   6 138 114  47  47
##   7 180 151 116 115
##   9  85 109 152 152
## 19 more rows ...
## 
## $samples
##   group lib.size norm.factors
## 1     1     3687            1
## 2     1     3386            1
## 3     1     2220            1
## 4     1     2298            1
```

Note how the "library sizes" (i.e., the total number of cells in each sample) are much lower for the WT samples 3 and 4.
This reflects the fact that the WT-only cluster of erythroid cells has been removed.
The differential analysis will subsequently test for proportional differences in the non-blood cells.


```r
y.ab3 <- estimateDisp(y.ab3, design)
fit.ab3 <- glmQLFit(y.ab3, design, robust=TRUE, abundance.trend=FALSE)
res3 <- glmQLFTest(fit.ab3, contrast=con)
topTags(res3, n=10)
```

```
## Coefficient:  1*genotypeko -1*genotypewt 
##         logFC   logCPM        F       PValue          FDR
## 1  -0.8078762 17.02635 86.73368 1.936518e-09 4.647643e-08
## 21 -6.1223245 12.49550 74.87789 7.652775e-09 8.197762e-08
## 9  -1.2887352 15.41361 72.52985 1.024720e-08 8.197762e-08
## 23  1.5930356 14.38228 39.28521 1.770323e-06 1.062194e-05
## 22  1.1940726 15.00913 37.58612 2.472628e-06 1.186861e-05
## 14  0.7892022 15.99919 35.31604 3.923221e-06 1.569288e-05
## 18  1.0488587 15.11092 32.04962 7.880944e-06 2.702038e-05
## 13  0.7779225 15.48045 24.03432 5.323241e-05 1.596972e-04
## 10  0.4730759 16.49772 18.81664 2.237308e-04 5.966154e-04
## 6   0.7710781 14.89581 15.69649 5.801561e-04 1.392375e-03
```

A similar strategy can be used to focus on proportional changes within a single subpopulation of a very heterogeneous data set.
For example, if we collected a whole blood data set, we could subset to T cells and test for changes in T cell subtypes (memory, killer, regulatory, etc.) 
using the total number of T cells in each sample as the library size.
This avoids detecting changes in T cell subsets that are driven by compositional effects from changes in abundance of, say, B cells in the same sample.

# Concluding remarks 

This workflow focuses on finding clusters that change in abundance, which are the most obvious manifestation of differences between conditions.
The complementary approach is to look for changes in expression between conditions _within_ each cluster.
As we will see, these two analyses represent two sides of the same coin when dealing with clusters defined by expression data.

All software packages used in this workflow are publicly available from the Comprehensive R Archive Network (https://cran.r-project.org) or the Bioconductor project (http://bioconductor.org). 
The specific version numbers of the packages used are shown below, along with the version of the R installation.


```r
sessionInfo()
```

```
## R Under development (unstable) (2019-02-19 r76128)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 16.04.5 LTS
## 
## Matrix products: default
## BLAS: /home/cri.camres.org/lun01/Software/R/trunk/lib/libRblas.so
## LAPACK: /home/cri.camres.org/lun01/Software/R/trunk/lib/libRlapack.so
## 
## locale:
##  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] pheatmap_1.0.12             scater_1.11.11             
##  [3] ggplot2_3.1.0               edgeR_3.25.3               
##  [5] limma_3.39.12               SingleCellExperiment_1.5.2 
##  [7] SummarizedExperiment_1.13.0 DelayedArray_0.9.8         
##  [9] BiocParallel_1.17.15        matrixStats_0.54.0         
## [11] Biobase_2.43.1              GenomicRanges_1.35.1       
## [13] GenomeInfoDb_1.19.2         IRanges_2.17.4             
## [15] S4Vectors_0.21.10           BiocGenerics_0.29.1        
## [17] BiocStyle_2.11.0           
## 
## loaded via a namespace (and not attached):
##  [1] viridis_0.5.1            BiocSingular_0.99.12    
##  [3] viridisLite_0.3.0        splines_3.6.0           
##  [5] DelayedMatrixStats_1.5.2 assertthat_0.2.0        
##  [7] statmod_1.4.30           BiocManager_1.30.4      
##  [9] highr_0.7                GenomeInfoDbData_1.2.0  
## [11] vipor_0.4.5              yaml_2.2.0              
## [13] pillar_1.3.1             lattice_0.20-38         
## [15] glue_1.3.0               digest_0.6.18           
## [17] RColorBrewer_1.1-2       XVector_0.23.0          
## [19] colorspace_1.4-0         cowplot_0.9.4           
## [21] htmltools_0.3.6          Matrix_1.2-16           
## [23] plyr_1.8.4               pkgconfig_2.0.2         
## [25] bookdown_0.9             zlibbioc_1.29.0         
## [27] purrr_0.3.0              scales_1.0.0            
## [29] processx_3.2.1           tibble_2.0.1            
## [31] withr_2.1.2              lazyeval_0.2.1          
## [33] magrittr_1.5             crayon_1.3.4            
## [35] evaluate_0.13            ps_1.3.0                
## [37] beeswarm_0.2.3           tools_3.6.0             
## [39] stringr_1.4.0            munsell_0.5.0           
## [41] locfit_1.5-9.1           irlba_2.3.3             
## [43] callr_3.1.1              compiler_3.6.0          
## [45] rsvd_1.0.0               rlang_0.3.1             
## [47] grid_3.6.0               RCurl_1.95-4.11         
## [49] BiocNeighbors_1.1.11     labeling_0.3            
## [51] bitops_1.0-6             rmarkdown_1.11          
## [53] gtable_0.2.0             codetools_0.2-16        
## [55] R6_2.4.0                 gridExtra_2.3           
## [57] knitr_1.21               dplyr_0.8.0.1           
## [59] stringi_1.3.1            ggbeeswarm_0.6.0        
## [61] compareSingleCell_0.98.0 Rcpp_1.0.0              
## [63] tidyselect_0.2.5         xfun_0.5
```

# References

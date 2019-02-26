---
title: Merging the embryo chimera data set
author:
- name: Aaron T. L. Lun
  affiliation: &CRUK Cancer Research UK Cambridge Institute, Li Ka Shing Centre, Robinson Way, Cambridge CB2 0RE, United Kingdom
- name: Jonathan A. Griffiths
  affiliation: *CRUK
date: "2019-02-26"
vignette: >
  %\VignetteIndexEntry{03. Embryo merging}
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

Once a processed `SingleCellExperiment` is available, the next step is to merge different samples to place them on the same coordinate system.
This allows us to perform downstream procedures like dimensionality reduction and clustering without needing to model sample-to-sample variation in expression.
We will merge using the mutual nearest neighbors (MNN) method [@haghverdi2018batch] as described [elsewhere](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/batch.html),
using the `SingleCellExperiment` object that we constructed in the [previous workflow](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/embryo_merge.html).


```r
library(SingleCellExperiment)
sce <- readRDS("embryo_processed.rds")
```

# Modelling the technical variability

We model the technical noise in each sample using the `multiBlockVar()` function as described [here](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/var.html#fitting-batch-specific-trends).
As this data set does not contain spike-ins, we set `make.tech.trend=TRUE` to estimate the technical component based on Poisson noise -  
see [here](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/tenx.html#modelling-the-mean-variance-trend) for details.


```r
library(scran)
dec.out <- multiBlockVar(sce, block=sce$sample, 
    make.tech.trend=TRUE, weighted=FALSE)
dec.out[,-7] # don't show per-block results.
```

```
## DataFrame with 29453 rows and 6 columns
##                                mean                total
##                           <numeric>            <numeric>
## Xkr4           5.29782915802558e-05 4.13932024379949e-05
## Gm1992                            0                    0
## Gm37381                           0                    0
## Rp1            0.000364210237139479 0.000428754187376142
## Sox17            0.0812588334691563    0.122990549478257
## ...                             ...                  ...
## AC149090.1        0.130897470873326    0.135827036918682
## DHRSX             0.262823762140174     0.24302992192422
## Vmn2r122       8.87345580461688e-05 0.000106737528207588
## CAAA01147332.1   0.0015295605947613  0.00163468918662797
## tomato-td         0.108275242720852     0.10373134041948
##                                  bio                 tech
##                            <numeric>            <numeric>
## Xkr4           -1.57099248267328e-05 5.71031272647278e-05
## Gm1992                             0                    0
## Gm37381                            0                    0
## Rp1             5.23629276668067e-05 0.000376391259709335
## Sox17             0.0403341523186133   0.0826563971596436
## ...                              ...                  ...
## AC149090.1       0.00603041734862433    0.129796619570058
## DHRSX          -0.000951030047235964    0.243980951971456
## Vmn2r122        1.50503096507976e-05 9.16872185567908e-05
## CAAA01147332.1  8.82386582789878e-09  0.00163468036276215
## tomato-td        0.00100594066801961     0.10272539975146
##                              p.value                   FDR
##                            <numeric>             <numeric>
## Xkr4                               1                     1
## Gm1992                             1                     1
## Gm37381                            1                     1
## Rp1             8.04439454819487e-27  6.38113527142428e-26
## Sox17          8.11961936219155e-299 3.86344344223954e-297
## ...                              ...                   ...
## AC149090.1      0.000326312835877702   0.00102156589658864
## DHRSX              0.736539614509996                     1
## Vmn2r122        1.60195121762685e-07  6.36051081326014e-07
## CAAA01147332.1      0.74681609334657                     1
## tomato-td       2.06028029843154e-10  9.46520599433851e-10
```

We also turn off weighting to ensure that each sample contributes equally to the results, regardless of the number of cells.
This is desirable here because different conditions may have different sets of highly variable genes (HVGs).
If one condition contains more cells, weighting by the number of cells would bias the combined statistics in favour of that condition's HVGs.
This differs from more typical applications of `multiBlockVar()` where all samples are replicates.
In such cases, the underlying set of HVGs should be similar and thus weighting could be used to improve precision without biasing the results.

The trend lies close to the lower bound of the per-gene variances for each sample (Figure \@ref(fig:trendplots)).
This indicates that the assumption of Poisson noise is reasonable.


```r
par(mfrow=c(2,2))
per.block.stats <- dec.out$per.block
for (i in colnames(per.block.stats)) {
    cur.stats <- per.block.stats[,i]
    plot(cur.stats$mean, cur.stats$total, xlab="Mean log-expression", 
        ylab="Variance", main=sprintf("Sample %s", i), pch=16)
    curve(metadata(cur.stats)$trend(x), add=TRUE, col="dodgerblue", lwd=2)
}
```

<div class="figure">
<img src="embryo_merge_files/figure-html/trendplots-1.png" alt="Mean-variance relationship in each sample, using statistics computed from log-normalized expression values. Each point represents a gene and the blue line represents the fitted Poisson trend in each sample." width="100%"  class="widefigure" />
<p class="caption">(\#fig:trendplots)Mean-variance relationship in each sample, using statistics computed from log-normalized expression values. Each point represents a gene and the blue line represents the fitted Poisson trend in each sample.</p>
</div>

**Comments from Aaron:**

- Turning off weighting in `multiBlockVar()` (and indirectly, `combineVar()`) assumes that there are no samples with very few cells.
Such samples will have imprecise variance estimates that, without weighting, would add extra uncertainty to the combined statistics.

## Identifying genes of interest

We define the features of interest as those with net biological components greater than zero.
This enriches for genes that contain some biological variation, reducing the effect of uninteresting Poisson noise on downstream analyses.


```r
to.use <- dec.out$bio > 0
summary(to.use)
```

```
##    Mode   FALSE    TRUE 
## logical   14108   15345
```

The injected KO cells were derived from a single male cell line while the host embryos were a mix of male or female mice.
Thus, there is a (largely uninteresting) sex effect between the WT and KO samples.
To mitigate this, we remove _Xist_ and genes on the Y chromosome.


```r
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
loc <- mapIds(TxDb.Mmusculus.UCSC.mm10.ensGene, key=rowData(sce)$ENSEMBL, 
    keytype="GENEID", column="CDSCHROM")
is.y <- loc=="chrY" & !is.na(loc)
to.use <- to.use & !is.y 
to.use[which(rowData(sce)$SYMBOL=="Xist")] <- FALSE
sum(to.use)
```

```
## [1] 15338
```

We also remove the tdTomato marker used to sort for KO cells.


```r
to.use[rowData(sce)$SYMBOL=="tomato-td"] <- FALSE
sum(to.use)
```

```
## [1] 15337
```

# Dimensionality reduction with PCA

We use principal components analysis (PCA) to perform dimensionality reduction prior to merging.
This reduces computational work and removes some high-dimensional noise, as [previously discussed](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/reads.html#denoising-expression-values-using-pca).
We perform PCA across all samples by using `multiBatchPCA()` on our selected subset of genes in `to.use`.
This ensures each sample contributes equally to the definition of the coordinate space, as described [here](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/batch.html#hierarchical-merging).


```r
library(batchelor)
library(BiocSingular)

set.seed(101) # for irlba.
pc.out <- batchelor::multiBatchPCA(sce, batch=sce$sample, 
    subset.row=to.use, get.variance=TRUE,
    BSPARAM=IrlbaParam(deferred=TRUE))
pc.out # one output matrix per level of 'sample'.
```

```
## List of length 4
## names(4): 1 2 3 4
```

By default, `multiBatchPCA()` will return 50 PCs for all samples.
We use `denoisePCANumber()` to choose the number of PCs to retain based on our previous estimates of technical noise in this data set.
This discards later PCs until the variance lost is equal to the total technical component.


```r
to.retain <- denoisePCANumber(
    metadata(pc.out)$var.explained, # variance explained per PC.
    sum(dec.out$tech[to.use]), # technical noise in subset of genes.
    metadata(pc.out)$var.total # total variance in the data
)
to.retain
```

```
## [1] 20
```

We then subset the matrices to only retain the first `to.retain` PCs in each sample.


```r
for (i in seq_along(pc.out)) {
    pc.out[[i]] <- pc.out[[i]][,seq_len(to.retain),drop=FALSE]
}
```

**Comments from Aaron:**

- The `BSPARAM=` argument instructs `multiBatchPCA()` to use methods from *[irlba](https://CRAN.R-project.org/package=irlba)* to speed up the PCA.
This is done by using an approximate algorithm, with deferred centered and scaling to preserve sparsity.
For large data sets, it is possible to achieve further speed gains with parallelization via the `BPPARAM=` argument; 
or by switching to randomized SVD in *[rsvd](https://CRAN.R-project.org/package=rsvd)* via `BiocSingular::RandomParam()`.

# Batch correction with MNN

We use the `fastMNN()` function in a hierarchical manner as described [elsewhere](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/batch.html#hierarchical-merging).
This involves merging the most similar samples before merging those that are more different, to weaken the assumption of shared populations required for correct MNN detection.
In this case, we first merge samples from the same genotype to remove the batch effect.
Note the use of `pc.input=TRUE` to specify that the input values are already in low-dimensional PC space.


```r
ko.out <- batchelor::fastMNN(pc.out[["1"]], pc.out[["2"]], 
    pc.input=TRUE, compute.variances=TRUE)
metadata(ko.out)$lost.var
```

```
## [1] 0.01787954 0.01787954
```

```r
wt.out <- batchelor::fastMNN(pc.out[["3"]], pc.out[["4"]], 
    pc.input=TRUE, compute.variances=TRUE)
metadata(wt.out)$lost.var
```

```
## [1] 0.02629509 0.02629509
```

The `lost.var` represents the proportion of variance in each batch that is removed by the batch correction^[Specifically, the orthogonalization step, as described [here](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/batch.html#with-diagnostics).].
If the assumptions underlying the MNN approach hold, this should be low and represent removal of noise along the batch vector.
A large proportion of lost variance (>10%) may be cause for concern as it suggests that biological structure within each batch is being discarded.

Our next step is to merge samples across genotypes to remove technical and uninteresting biological differences.
This includes sex and changes in expression induced by injection or cell sorting.


```r
overall <- batchelor::fastMNN(ko.out$corrected, wt.out$corrected,
    pc.input=TRUE, compute.variances=TRUE)
metadata(overall)$lost.var
```

```
## [1] 0.0419255 0.0419255
```

We store the result in the `reducedDims` slot of our `SingleCellExperiment` object, for use in downstream functions.


```r
# Cell order is the same between sce and corrected: see comments.
reducedDim(sce, "corrected") <- overall$corrected
```

Note that the merge to create `overall` will also eliminate changes in expression caused by loss of _Tal1_.
This is necessary for the purposes of mapping all cells onto a common coordinate system.
Without such correction, cells of the same type would be separated by the genotype difference across samples, precluding common annotation in downstream analyses.
Nonetheless, differential expression upon KO is of substantial biological interest and will be recovered in our downstream differential testing.

**Comments from Aaron:**

- The `reducedDim` assignment assumes that the order of cells in `overall$corrected` is the same as the order of cells in `sce`.
This is already the case here, as `ko.out$corrected` contains cells from samples 1 and 2 (in that order) while `wt.out$corrected` contains cells from samples 3 and 4 (in that order).
However, this may not be true for arbitrary merge orders!
In such cases, we suggest assigning unique column names to each cell so that one can `match()` the row names of `corrected` with the column names of `sce` prior to assignment.

# Visual inspection of merge results

We use $t$-stochastic neighbor embedding (t-SNE) plots [@van2008visualizing] to perform further dimensionality reduction for visualization.
In the uncorrected data, cells clearly separate by genotype with no visible batch effects between replicates (Figure \@ref(fig:beforeplot)).


```r
library(scater)
old <- sce
reducedDim(sce, "PCA") <- do.call(rbind, pc.out)
plotTSNE(sce, rerun=TRUE, run_args=list(use_dimred="PCA"), colour_by="sample")
```

<div class="figure">
<img src="embryo_merge_files/figure-html/beforeplot-1.png" alt="t-SNE plot of the embryo data before MNN correction. Each point represents a cell and is coloured by the sample of origin (KO: 1 and 2, WT: 3 and 4)." width="100%" />
<p class="caption">(\#fig:beforeplot)t-SNE plot of the embryo data before MNN correction. Each point represents a cell and is coloured by the sample of origin (KO: 1 and 2, WT: 3 and 4).</p>
</div>

After correction, cells from all samples are merged together in the majority of subpopulations (Figure \@ref(fig:afterplot)).
This is consistent with the removal of the inter-genotype differences and simplifies annotation and interpretation in downstream analyses.


```r
sce <- runTSNE(sce, use_dimred="corrected", colour_by="sample")
plotTSNE(sce, colour_by="sample")
```

<div class="figure">
<img src="embryo_merge_files/figure-html/afterplot-1.png" alt="t-SNE plot of the embryo data after MNN correction. Each point represents a cell and is coloured by the sample of origin." width="100%" />
<p class="caption">(\#fig:afterplot)t-SNE plot of the embryo data after MNN correction. Each point represents a cell and is coloured by the sample of origin.</p>
</div>

# Defining common annotation

## Clustering cells

We use the shared nearest-neighbour approach [@xu2015identification] to cluster cells in the corrected space, as described
[here](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/umis.html#clustering-cells-into-putative-subpopulations)
and [here](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/batch.html#using-the-corrected-values-in-downstream-analyses).


```r
snn.gr <- buildSNNGraph(sce, use.dimred="corrected")
clusters <- igraph::cluster_walktrap(snn.gr)
table(clusters$membership, sce$sample)
```

```
##     
##        1   2   3   4
##   1  388 338 389 423
##   2    0   1 833 873
##   3  511 472 301 308
##   4    0   1 129 181
##   5    0   1 157 199
##   6  138 114  47  47
##   7  180 151 116 115
##   8    0   0  41  23
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

We visually examine the clusters on a t-SNE plot to confirm that a sensible partitioning was generated (Figure \@ref(fig:clusterplot)).


```r
sce$cluster <- factor(clusters$membership)
plotTSNE(sce, colour_by="cluster", text_by="cluster", text_colour="red")
```

<div class="figure">
<img src="embryo_merge_files/figure-html/clusterplot-1.png" alt="t-SNE plot of the MNN-corrected embryo data, where each point represents a cell and is coloured by the identity of the cluster to which it is assigned. The cluster number is also shown at the median coordinates across all cells in that cluster." width="100%" />
<p class="caption">(\#fig:clusterplot)t-SNE plot of the MNN-corrected embryo data, where each point represents a cell and is coloured by the identity of the cluster to which it is assigned. The cluster number is also shown at the median coordinates across all cells in that cluster.</p>
</div>

**Comments from Aaron:**

- More sophisticated diagnostics for graph-based clustering are possible with `clusterModularity()`,
see [this section](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/umis.html#evaluating-graph-based-clusters) for more details.
For brevity's sake, we will skip that step here.
- We do not run `library(igraph)`, but instead use `igraph::` to extract methods from the *[igraph](https://CRAN.R-project.org/package=igraph)* package. 
This is because *[igraph](https://CRAN.R-project.org/package=igraph)* contains a normalize method that will override its counterpart from *[scater](https://bioconductor.org/packages/3.9/scater)*, resulting in some unusual bugs.

## Annotating cell types

We use `findMarkers()` to identify the genes that define each cluster.
This is done by testing for differential expression between each pair of clusters and consolidating the results into a single table per cluster -
see [here](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/reads.html#detecting-marker-genes-between-clusters) for details.
We also block on the sample of origin to avoid confounding effects from sample-to-sample variability or differential expression between genotypes.


```r
markers <- findMarkers(sce, sce$cluster, block=sce$sample)
```



Of particular interest is cluster 5.
This upregulates a range of hemoglobin genes (Figure \@ref(fig:bloodheat)) and probably represents cells in the erythroid lineage.


```r
blood.set <- markers[["5"]]
as.data.frame(blood.set[1:20,1:3])
```

```
##          Top       p.value           FDR
## Hbb-bh1    1  0.000000e+00  0.000000e+00
## Blvrb      1  0.000000e+00  0.000000e+00
## Gypa       1  0.000000e+00  0.000000e+00
## Mest       1  0.000000e+00  0.000000e+00
## Alas2      1 3.304393e-297 2.949221e-294
## Smim1      1 6.197642e-287 5.215404e-284
## Hbb-y      1 2.156206e-262 1.512065e-259
## Klf1       1 2.517351e-232 1.300764e-229
## Acp5       1 1.376621e-188 4.885015e-186
## Mrap       1 1.774735e-123 2.694395e-121
## Hba-x      2  0.000000e+00  0.000000e+00
## Prdx2      2  0.000000e+00  0.000000e+00
## Cited4     2  0.000000e+00  0.000000e+00
## Hebp1      2 2.561026e-254 1.676220e-251
## Marcksl1   2 2.219265e-202 8.715200e-200
## Trim10     2 1.212626e-118 1.725385e-116
## Hba-a1     3  0.000000e+00  0.000000e+00
## Tmem14c    3  0.000000e+00  0.000000e+00
## Mgst3      3  0.000000e+00  0.000000e+00
## Tmsb10     3 9.852223e-221 4.534024e-218
```

```r
logFCs <- as.matrix(blood.set[1:50,-(1:3)])
colnames(logFCs) <- sub("logFC.", "", colnames(logFCs))

library(pheatmap)
max.lfc <- max(abs(range(logFCs)))
pheatmap(logFCs, breaks=seq(-5, 5, length.out=101))
```

<div class="figure">
<img src="embryo_merge_files/figure-html/bloodheat-1.png" alt="Heatmap of the log-fold changes for the top 50 genes expressed in cluster 5 compared to all other clusters. Each column represents another cluster." width="100%" />
<p class="caption">(\#fig:bloodheat)Heatmap of the log-fold changes for the top 50 genes expressed in cluster 5 compared to all other clusters. Each column represents another cluster.</p>
</div>

One could repeat this procedure with all of the other clusters.
However, it is more efficient to use our differential analyses to prioritize clusters of interest, so we will delay further annotation until that point.

## Detecting doublets

One potentially problematic feature of this data set is its high doublet frequency.
Thus, we want a measure of how "doublet-like" each cell is, in order to assist downstream interpretation of our results.
This is achieved using the `doubletCells()` function on each sample, as described [here](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/doublets.html#doublet-detection-by-simulation).


```r
set.seed(1000)
doublets <- numeric(ncol(sce))
for (i in levels(sce$sample)) {
    keep <- sce$sample==i
    cur.sce <- sce[,keep]
    scores <- doubletCells(cur.sce, BSPARAM=BiocSingular::IrlbaParam(deferred=TRUE))
    doublets[keep] <- scores
}
summary(doublets)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 0.00000 0.01083 0.03907 0.16942 0.11573 7.06882
```

Some clusters are entirely composed of putative doublets, while others simply contain a small proportion of doublets (Figure \@ref(fig:doubletplot)).
The former are not of any interest and can be dismissed immediately.
The latter are salvageable but require some care in interpretation during downstream analyses.


```r
sce$doublet <- doublets
plotColData(sce, "doublet", "cluster", colour_by="sample")
```

<div class="figure">
<img src="embryo_merge_files/figure-html/doubletplot-1.png" alt="Distribution of per-cell doublet scores within each cluster. Each cell is a point coloured by sample of origin." width="100%" />
<p class="caption">(\#fig:doubletplot)Distribution of per-cell doublet scores within each cluster. Each cell is a point coloured by sample of origin.</p>
</div>

We perform a complementary analysis based on the clusters directly, as described [here](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/doublets.html#doublet-detection-with-clusters).
The top-ranking doublet-like clusters are defined by the absence of any uniquely expressed genes distinguishing them from a putative pair of source clusters.
These are consistent with the most obvious offenders in Figure \@ref(fig:doubletplot).


```r
by.clust <- doubletCluster(sce, sce$cluster, block=sce$sample)
by.clust[,1:8]
```

```
## DataFrame with 28 rows and 8 columns
##         source1     source2         N        best               p.value
##     <character> <character> <integer> <character>             <numeric>
## 4            10           2         0      Ndufa4      0.14761091364786
## 8            14           2         0        Rhoj    0.0609629735308206
## 16           25          12        20       Sfrp5  8.79778267485982e-09
## 23           24           8        34        Nnat  6.21436797217934e-17
## 18           24          12        35        Cdx2  1.06444750786317e-13
## ...         ...         ...       ...         ...                   ...
## 7            27           4       304     Hoxb5os 3.38336457413179e-111
## 10            8           1       309     Hoxaas3  1.24953668286433e-37
## 25           27           8       325        Myl7  3.59662188287309e-72
## 1            19          17       375        Lix1  4.18887342887083e-85
## 9            28          10       411      Tfap2b 4.75114616938885e-120
##             lib.size1         lib.size2                prop
##             <numeric>         <numeric>           <numeric>
## 4   0.853759467947534 0.530944023646776  0.0221667854597292
## 8   0.571674912241086 0.589986245971301 0.00456165359942979
## 16  0.817886018536778 0.997633602839677   0.017320028510335
## 23   1.14547051831947  2.51733760529172  0.0171062009978617
## 18  0.887634150248278 0.810347589300016  0.0287241625089095
## ...               ...               ...                 ...
## 7      0.655635341867  2.41425449355515  0.0400570206699929
## 10   1.05407452287186 0.522417449257801  0.0760513186029936
## 25  0.886075949367089  2.93628691983122  0.0075552387740556
## 1   0.845710972124425 0.968603736072568   0.109622238061297
## 9   0.767270879213959  1.56578011180756  0.0354953670705631
```

We remove the obvious doublet clusters prior to any downstream analysis.
We give the benefit of the doubt to clusters containing but not dominated by doublets, provided that their results are treated with caution.


```r
offenders <- c("4", "8")
retain <- !sce$cluster %in% offenders
sce <- sce[,retain]
summary(retain)
```

```
##    Mode   FALSE    TRUE 
## logical     375   13655
```



# Concluding remarks

We save the `SingleCellExperiment` object with the merged coordinates to file for downstream use in differential testing.
We also save the test results for later annotation.


```r
saveRDS(sce, "embryo_merged.rds")
saveRDS(markers, "embryo_markers.rds")
```

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
## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] pheatmap_1.0.12                       
##  [2] scater_1.11.11                        
##  [3] ggplot2_3.1.0                         
##  [4] BiocSingular_0.99.12                  
##  [5] batchelor_0.99.0                      
##  [6] TxDb.Mmusculus.UCSC.mm10.ensGene_3.4.0
##  [7] GenomicFeatures_1.35.7                
##  [8] AnnotationDbi_1.45.0                  
##  [9] scran_1.11.20                         
## [10] SingleCellExperiment_1.5.2            
## [11] SummarizedExperiment_1.13.0           
## [12] DelayedArray_0.9.8                    
## [13] BiocParallel_1.17.15                  
## [14] matrixStats_0.54.0                    
## [15] Biobase_2.43.1                        
## [16] GenomicRanges_1.35.1                  
## [17] GenomeInfoDb_1.19.2                   
## [18] IRanges_2.17.4                        
## [19] S4Vectors_0.21.10                     
## [20] BiocGenerics_0.29.1                   
## [21] BiocStyle_2.11.0                      
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-6             bit64_0.9-7             
##  [3] RColorBrewer_1.1-2       httr_1.4.0              
##  [5] progress_1.2.0           dynamicTreeCut_1.63-1   
##  [7] tools_3.6.0              R6_2.4.0                
##  [9] irlba_2.3.3              vipor_0.4.5             
## [11] DBI_1.0.0                lazyeval_0.2.1          
## [13] colorspace_1.4-0         withr_2.1.2             
## [15] tidyselect_0.2.5         gridExtra_2.3           
## [17] prettyunits_1.0.2        processx_3.2.1          
## [19] bit_1.1-14               compiler_3.6.0          
## [21] BiocNeighbors_1.1.11     labeling_0.3            
## [23] rtracklayer_1.43.1       bookdown_0.9            
## [25] scales_1.0.0             callr_3.1.1             
## [27] Rsamtools_1.99.2         stringr_1.4.0           
## [29] digest_0.6.18            rmarkdown_1.11          
## [31] XVector_0.23.0           pkgconfig_2.0.2         
## [33] htmltools_0.3.6          highr_0.7               
## [35] limma_3.39.12            rlang_0.3.1             
## [37] RSQLite_2.1.1            DelayedMatrixStats_1.5.2
## [39] dplyr_0.8.0.1            RCurl_1.95-4.11         
## [41] magrittr_1.5             GenomeInfoDbData_1.2.0  
## [43] Matrix_1.2-16            Rcpp_1.0.0              
## [45] ggbeeswarm_0.6.0         munsell_0.5.0           
## [47] viridis_0.5.1            stringi_1.3.1           
## [49] yaml_2.2.0               edgeR_3.25.3            
## [51] zlibbioc_1.29.0          plyr_1.8.4              
## [53] grid_3.6.0               blob_1.1.1              
## [55] crayon_1.3.4             lattice_0.20-38         
## [57] cowplot_0.9.4            Biostrings_2.51.2       
## [59] hms_0.4.2                locfit_1.5-9.1          
## [61] knitr_1.21               ps_1.3.0                
## [63] pillar_1.3.1             igraph_1.2.4            
## [65] codetools_0.2-16         biomaRt_2.39.2          
## [67] XML_3.98-1.17            glue_1.3.0              
## [69] evaluate_0.13            BiocManager_1.30.4      
## [71] gtable_0.2.0             purrr_0.3.0             
## [73] assertthat_0.2.0         xfun_0.5                
## [75] rsvd_1.0.0               viridisLite_0.3.0       
## [77] tibble_2.0.1             compareSingleCell_0.98.0
## [79] GenomicAlignments_1.19.1 beeswarm_0.2.3          
## [81] memoise_1.1.0            statmod_1.4.30
```

# References

---
title: Differential expression analyses with the embryo chimera data set
author:
- name: Aaron T. L. Lun
  affiliation: &CRUK Cancer Research UK Cambridge Institute, Li Ka Shing Centre, Robinson Way, Cambridge CB2 0RE, United Kingdom
- name: Jonathan A. Griffiths
  affiliation: *CRUK
date: "2019-02-26"
vignette: >
  %\VignetteIndexEntry{05. Differential expression}
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

Another obvious differential analysis is to look for changes in expression between conditions within each cluster.
This allows us to identify cell type-specific transcriptional differences caused by the loss of _Tal1_.
Differential expression (DE) analyses are very much complementary to differential abundance analyses;
these represent two sides of the same coin when dealing with clusters defined from expression data.

To illustrate this duality, consider a scRNA-seq experiment involving a heterogeneous population and two biological conditions.
Assume we have a cell type $X$ present in both conditions.
Within $X$, some genes are differentially expressed between conditions.
This leads to two possible outcomes:

- The DE between conditions causes $X$ to form two separate clusters in expression space (or the low-dimensional space derived from the gene expression values).
This manifests as differential abundance where one cluster is enriched in each condition.
- The DE between conditions is not sufficient to split $X$ into two separate clusters, e.g., because `fastMNN()` identifies them as corresponding cell types and merges them together.
This means that the differences between conditions manifests directly as DE within the single cluster corresponding to $X$.

Thus, it is often necessary to consider both differential expression and abundance to fully characterize the differences between conditions in scRNA-seq experiments.
In this workflow, we will be identifying differentially expressed genes (DEGs) upon loss of _Tal1_ function.
We will focus on DEGs in the cluster that we [previously](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/embryo_merge.html) annotated as placenta precursors.


```r
library(SingleCellExperiment)
sce <- readRDS("embryo_merged.rds")
placenta <- 13
```



# Setting up the count matrix

We create a new count matrix where we sum counts together for all cells in the same cluster and sample.
The summation creates "pseudo-bulk" count vectors that are more amenable to downstream analysis with standard tools such as *[edgeR](https://bioconductor.org/packages/3.9/edgeR)*.
More importantly, it reflects the fact that our biological replication occurs at the sample level [@lun2017overcoming].
Supplying the per-cell counts directly would indicate that each cell is a biological replicate^[Unless one is using a mixture model, though this has its own problems.], which is not true in an experimental setting.


```r
library(scater)
cluster.sample <- sprintf("Cluster%i.Sample%i", sce$cluster, sce$sample)
summed <- sumCountsAcrossCells(sce, cluster.sample)
dim(summed)
```

```
## [1] 29453   101
```

We will subset this matrix to our cluster of interest, i.e., cluster 13.
We do not use all clusters in the DE analysis as the strong DE between clusters makes it difficult to compute a sensible average abundance to model the mean-dispersion trend.
Batch effects may also differ between clusters, which would not be easily handled with a single additive term in the design matrix for the batch.
Finally, it is usually more convenient to fit a number of small generalized linear models (one per cluster) than to try to fit one large model involving all clusters.


```r
keep <- grep(paste0("^Cluster", placenta), colnames(summed))
summed <- summed[,keep]
dim(summed)
```

```
## [1] 29453     4
```

We extract the sample of origin for each column from the column names of `summed`.


```r
head(colnames(summed))
```

```
## [1] "Cluster13.Sample1" "Cluster13.Sample2" "Cluster13.Sample3"
## [4] "Cluster13.Sample4"
```

```r
sample.id <- as.integer(sub(".*([0-9]+)$", "\\1", colnames(summed)))
head(sample.id)
```

```
## [1] 1 2 3 4
```

Finally, we construct a `DGEList` object for use in *[edgeR](https://bioconductor.org/packages/3.9/edgeR)* [@robinson2010edgeR].


```r
library(edgeR)
y.exp <- DGEList(summed)
```

# Testing for differential expression

## Filtering and normalization

A typical step in bulk RNA-seq data analyses is to remove samples with very low library sizes corresponding to failed library preparation or sequencing.
In our situation, this is equivalent to removing cluster-sample combinations that have very few or lowly-sequenced cells.
The corresponding summed count vectors are likely to be highly variable and thus reduce power for DE detection.
The exact definition of "very low" will vary, but in this case, we define it to be log-library sizes that are more than 3 median absolute deviations from the median. 


```r
discarded <- isOutlier(y.exp$samples$lib.size, log=TRUE, type="lower")
y.exp <- y.exp[,!discarded]
```

Another typical step in bulk RNA-seq analyses is to remove genes that are lowly expressed.
This reduces computational work, improves the accuracy of mean-variance trend modelling and decreases the severity of the multiple testing correction.
Here, we remove genes with an average count below 1 across samples^[This choice of threshold is somewhat arbitrary, but the exact choice is not critical.].


```r
keep <- aveLogCPM(y.exp) > aveLogCPM(1, lib.size=mean(y.exp$samples$lib.size))
y.exp <- y.exp[keep,]
summary(keep)
```

```
##    Mode   FALSE    TRUE 
## logical   17738   11715
```

Finally, we correct for composition biases by computing normalization factors with the trimmed mean of M-values method [@robinson2010scaling].


```r
y.exp <- calcNormFactors(y.exp)
y.exp$samples
```

```
##                   group lib.size norm.factors
## Cluster13.Sample1     1  2358841    0.9961578
## Cluster13.Sample2     1  2133587    0.9895658
## Cluster13.Sample3     1   708039    0.9977105
## Cluster13.Sample4     1   850959    1.0167699
```

## Modelling biological variability 

We set up the design matrix with one term for each genotype/cluster combination and an additive term for the batch effect between replicates.
Modelling the batch effect is necessary as `summed` is derived from the original count matrix, i.e., before batch correction.


```r
batch <- factor(c(1,2,1,2))[sample.id]
genotype <- rep(c("KO", "WT"), each=2)[sample.id]
design <- model.matrix(~0 + genotype + batch)
design
```

```
##   genotypeKO genotypeWT batch2
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

We estimate the negative binomial (NB) dispersions with `estimateDisp()`.
As previously mentioned, this models the mean-variance trend in count data (Figure \@ref(fig:bcvplot)).


```r
y.exp <- estimateDisp(y.exp, design)
summary(y.exp$trended.dispersion)
```

```
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## 0.001261 0.002084 0.003664 0.013743 0.018763 0.061121
```

```r
plotBCV(y.exp)
```

<div class="figure">
<img src="embryo_expression_files/figure-html/bcvplot-1.png" alt="Biological coefficient of variation (BCV) for each gene as a function of the average abundance. The BCV is computed as the square root of the NB dispersion after empirical Bayes shrinkage towards the trend. Trended and common BCV estimates are shown in blue and red, respectively." width="100%" />
<p class="caption">(\#fig:bcvplot)Biological coefficient of variation (BCV) for each gene as a function of the average abundance. The BCV is computed as the square root of the NB dispersion after empirical Bayes shrinkage towards the trend. Trended and common BCV estimates are shown in blue and red, respectively.</p>
</div>

We also estimate the quasi-likelihood dispersions with `glmQLFit()` [@chen2018reads].
This accounts for the uncertainty and variability of gene-specific dispersions (Figure \@ref(fig:qlplot)). 


```r
fit.exp <- glmQLFit(y.exp, design, robust=TRUE)
summary(fit.exp$var.prior)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.8595  0.9812  0.9964  1.1795  1.1420  2.4917
```

```r
summary(fit.exp$df.prior)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   3.726     Inf     Inf     Inf     Inf     Inf
```

```r
plotQLDisp(fit.exp)    
```

<div class="figure">
<img src="embryo_expression_files/figure-html/qlplot-1.png" alt="QL dispersion estimates for each gene as a function of abundance. Raw estimates (black) are shrunk towards the trend (blue) to yield squeezed estimates (red)." width="100%" />
<p class="caption">(\#fig:qlplot)QL dispersion estimates for each gene as a function of abundance. Raw estimates (black) are shrunk towards the trend (blue) to yield squeezed estimates (red).</p>
</div>

## Hypothesis testing for DEGs

We use the `glmQLFTest()` function to identify DEGs in the KO condition compared to the WT.
DEGs are defined as those with non-zero log-fold changes at a false discovery rate of 5%.


```r
con <- makeContrasts(genotypeKO - genotypeWT, levels=design)
res.exp <- glmQLFTest(fit.exp, contrast=con)
summary(decideTests(res.exp))
```

```
##        1*genotypeKO -1*genotypeWT
## Down                          103
## NotSig                      11475
## Up                            137
```

```r
topTags(res.exp, n=10)
```

```
## Coefficient:  1*genotypeKO -1*genotypeWT 
##              logFC    logCPM         F        PValue           FDR
## Hbb-bh1  -8.518881  8.487174 3128.9527  0.000000e+00  0.000000e+00
## Cdkn1c   -7.857900 11.663836 8026.0055  0.000000e+00  0.000000e+00
## Hbb-y    -7.765848  7.469557 1824.8134  0.000000e+00  0.000000e+00
## Hba-x    -7.142118  7.734239 2047.3651  0.000000e+00  0.000000e+00
## Phlda2   -2.081694 11.937407 1417.8290 5.307829e-293 1.243624e-289
## Xist    -11.898570  6.553760 1142.0563 5.679260e-239 1.108876e-235
## Hba-a1   -6.992932  5.739749  649.3852 1.837117e-139 3.074547e-136
## Hba-a2   -8.055647  4.551025  308.9038  2.836557e-68  4.153783e-65
## Erdr1     2.228181  7.599153  228.5864  1.370458e-40  1.783879e-37
## Acta2    -1.004403  8.794472  171.1849  7.801128e-39  9.139021e-36
```

Amusingly, the top DEGs are hemoglobins that are downregulated in the KO condition.
One would expect that these genes should not have been expressed at all in non-blood WT cells^[Though this may not be entirely unreasonable, as we shall discuss below.].
This result may be caused by contamination of each droplet from the ambient pool of extracellular RNA [@lun2018distinguishing;@young2018soupx],
which is filled with hemoglobin mRNAs in the WT samples but not in the KO samples.

Ambient contamination is a phenomenon that is generally most pronounced in massively multiplexed scRNA-seq protocols.
It tends to be less of an issue for plate-based experiments where the ambient solution is greatly diluted to ensure separation of individual cells in the cell sorter.
Thus, for data from protocols such as SMART-seq2, it is generally satisfactory to stop at this point.
For droplet-based experiments, some further work is required to reassure ourselves that the DE is not caused by differences in ambient RNA.

# Avoiding problems with ambient RNA

## Defining the ambient profile 

To eliminate the effects of ambient RNA, we need to obtain an estimate of the ambient "expression" profile.
This cannot be obtained from the count matrix in `sce`, so instead we need to go back to the raw count matrix for each sample.
We download and cache these matrices using *[BiocFileCache](https://bioconductor.org/packages/3.9/BiocFileCache)* as previousy described.


```r
library(BiocFileCache)
bfc <- BiocFileCache("raw_data", ask=FALSE)
sample.paths <- character(length(sample.id))
for (i in seq_along(sample.id)) {
    fname <- sprintf("sample_%s_unswapped.mtx.gz", sample.id[i]) 
    sample.paths[i] <- bfcrpath(bfc, file.path("https://content.cruk.cam.ac.uk/",
       "jmlab/chimera_tal1_data/unfiltered", fname))
}
```

We follow the approach used in `emptyDrops()` [@lun2018distinguishing] and consider all barcodes with total counts below 100 to represent empty droplets.


```r
collected <- vector("list", length(sample.id))
for (i in seq_along(collected)) {
    curmat <- Matrix::readMM(sample.paths[i])
    is.empty <- Matrix::colSums(curmat) < 100
    collected[[i]] <- Matrix::rowSums(curmat[,is.empty])
}
collected <- do.call(cbind, collected)
colnames(collected) <- sample.id
dim(collected)
```

```
## [1] 29453     4
```

We supply it with the original row names for the entire data set.


```r
gene.path <- bfcrpath(bfc, file.path("https://content.cruk.cam.ac.uk/",
    "jmlab/chimera_tal1_data/genes.tsv.gz"))
gene.tab <- read.delim(gene.path, header=FALSE, stringsAsFactors=FALSE)
rownames(collected) <- gene.tab[,1]
```

This allows us to easily match it with the rows in `sce`.
In this case, we did not do any subsetting of `sce` so the correspondence should be exactly 1:1.


```r
m <- match(rowData(sce)$ENSEMBL, rownames(collected))
stopifnot(all(!is.na(m))) # sanity check!
collected <- collected[m,]    
```

## Repeating the DE analysis

We combine `summed` with `collected` before creating a `DGEList` object.


```r
overall <- cbind(summed, collected)
y.amb <- DGEList(overall)
y.amb$samples 
```

```
##                   group lib.size norm.factors
## Cluster13.Sample1     1  2358841            1
## Cluster13.Sample2     1  2133587            1
## Cluster13.Sample3     1   708039            1
## Cluster13.Sample4     1   850959            1
## 1                     1  6302283            1
## 2                     1  6104236            1
## 3                     1  5330601            1
## 4                     1  6102033            1
```

We filter out low-abundance genes and estimate normalization factors as previously described.


```r
keep <- aveLogCPM(y.amb) > aveLogCPM(1, lib.size=mean(y.amb$samples$lib.size))
y.amb <- y.amb[keep,]
summary(keep)
```

```
##    Mode   FALSE    TRUE 
## logical   15402   14051
```

```r
y.amb <- calcNormFactors(y.amb)
y.amb$samples
```

```
##                   group lib.size norm.factors
## Cluster13.Sample1     1  2358841    1.0372072
## Cluster13.Sample2     1  2133587    1.0207244
## Cluster13.Sample3     1   708039    1.0341107
## Cluster13.Sample4     1   850959    1.0423806
## 1                     1  6302283    1.0228940
## 2                     1  6104236    1.0177315
## 3                     1  5330601    0.9230598
## 4                     1  6102033    0.9118828
```

Here, we use a new design matrix that accounts for the relationship between cluster 13's expression and the ambient pool in the same sample.
The first four coefficients represent the log-expression of the ambient pool in each sample,
while the last two coefficients represent the log-fold change of the summed cell expression profiles over the ambient pool in the KO or WT conditions. 


```r
s <- factor(rep(sample.id, 2))
new.geno <- rep(genotype, 2)
ambient <- rep(c("N", "Y"), each=4)
design.amb <- model.matrix(~0 + s + new.geno:ambient)

# Get to full rank:
design.amb <- design.amb[,!grepl("ambientY", colnames(design.amb))] 

# Syntactically valid colnames:
colnames(design.amb) <- make.names(colnames(design.amb)) 
design.amb
```

```
##   s1 s2 s3 s4 new.genoKO.ambientN new.genoWT.ambientN
## 1  1  0  0  0                   1                   0
## 2  0  1  0  0                   1                   0
## 3  0  0  1  0                   0                   1
## 4  0  0  0  1                   0                   1
## 5  1  0  0  0                   0                   0
## 6  0  1  0  0                   0                   0
## 7  0  0  1  0                   0                   0
## 8  0  0  0  1                   0                   0
```

We use this new design matrix to estimate the NB and QL dispersions.


```r
y.amb <- estimateDisp(y.amb, design.amb)
summary(y.amb$trended.dispersion)
```

```
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
## 0.0003460 0.0005035 0.0008534 0.0040607 0.0046828 0.0214139
```

```r
fit.amb <- glmQLFit(y.amb, design.amb, robust=TRUE)    
summary(fit.amb$var.prior)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.9519  1.0252  1.0522  1.2477  1.3522  2.1796
```

```r
summary(fit.amb$df.prior)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    2327     Inf     Inf     Inf     Inf     Inf
```

Finally, we test for differences between WT and KO.
The key here is to identify genes that have different cell/ambient log-fold changes between conditions.
This corresponds to a non-zero two-way interaction effect between genotype and the cell/ambient factors. 
(Equivalently, the interaction term can be interpreted as the difference in the KO/WT log-fold change computed from the cell profile compared to that from the ambient profiles.)
By doing so, we can mitigate the effect of differences in the ambient pool between conditions.


```r
con <- makeContrasts(new.genoKO.ambientN - new.genoWT.ambientN, levels=design.amb)
res.amb <- glmQLFTest(fit.amb, contrast=con)
summary(decideTests(res.amb))
```

```
##        1*new.genoKO.ambientN -1*new.genoWT.ambientN
## Down                                            537
## NotSig                                        12981
## Up                                              533
```

```r
topTags(res.amb, n=10)
```

```
## Coefficient:  1*new.genoKO.ambientN -1*new.genoWT.ambientN 
##              logFC    logCPM        F        PValue           FDR
## Cdkn1c  -4.2247076  9.606981 631.7684 6.892576e-138 9.684758e-134
## Acta2   -1.9892717  7.778484 496.7825 4.209652e-109 2.957491e-105
## Fth1     1.1398851 10.985979 484.4775 1.801908e-106 8.439535e-103
## Tmem14c  1.9357006  8.911486 438.1727  1.477921e-96  5.191567e-93
## Gpx1     1.6373299  8.779213 315.8539  2.812054e-70  7.902435e-67
## Dlk1    -1.4329059  8.052810 285.3870  1.038621e-63  2.432278e-60
## Krt8    -1.0306325  9.321042 260.0062  3.124343e-58  6.271449e-55
## Smim1    3.3618487  6.596244 240.8483  4.309657e-54  7.569375e-51
## Hmbs     2.0508429  7.756492 233.2395  1.904176e-52  2.972842e-49
## Prdx2    0.8987387 10.249319 225.4598  9.175036e-51  1.289184e-47
```

## Combining with the standard DE analysis

Unfortunately, testing the interaction term is not sufficient to avoid problems due to differences in the ambient pool.
Consider a situation where the ambient pool contains a transcript $G$ in one condition but not the other.
Further assume that there are droplets containing cells that also express $G$ at the same level between conditions.
If the ambient pool contaminates the droplets in all samples, a non-zero log-fold change for $G$ will be introduced between conditions for the cell expression profiles.
However, due to the presence of endogenous $G$, the log-fold change in the cells will not be as large as the log-fold change in the ambient profiles between conditions.
This causes the interaction term to be significantly non-zero.

Thus, some extra work is required to ensure that we are not detecting spurious DEGs.
Specifically, we only consider DEGs where the direct KO/WT comparison is significant; the interaction term is significant; 
and the KO/WT log-fold change is of the same sign as the interaction effect.
This focuses on genes with changes in expression beyond that expected from ambient contamination.

<!--
As formulated above, the interaction term is equivalent to the KO/WT log-fold change in cells minus the KO/WT log-fold change in the ambient pool.
Consider all the cases:

```
cells < ambient < 0 # Interesting
cells < 0 < ambient # Interesting
0 < cells < ambient # Possible artifact
0 < ambient < cells # Interesting
ambient < 0 < cells # Interesting
ambient < cells < 0 # Possible artifact
```

So, as long as `cells` is of the same sign as `cells - ambient`, we're in business!
-->



```r
common <- intersect(rownames(res.amb), rownames(res.exp))
tab.exp <- res.exp$table[common,]
tab.amb <- res.amb$table[common,]
okay <- sign(tab.exp$logFC)==sign(tab.amb$logFC) 
summary(okay)
```

```
##    Mode   FALSE    TRUE 
## logical    2702    8979
```

We compute a single $p$-value by taking the larger of the two $p$-values from the direct and interaction contrasts.
This is equivalent to an intersection-union test [@berger1978bioequivalence] and ensures that we can properly control the FDR later.
For all genes that are not of interest, we set their $p$-values to unity.


```r
iut.p <- pmax(tab.exp$PValue, tab.amb$PValue)
iut.p[!okay] <- 1
```

We use these statistics to create a final table of results.
This conservatively defines a set of DEGs in the presence of differences in the ambient pool between conditions.


```r
final <- data.frame(row.names=common,
    logFC=tab.exp$logFC, interaction=tab.amb$logFC,
    PValue=iut.p, FDR=p.adjust(iut.p, method="BH"))
final <- final[order(final$PValue),]
sum(final$FDR <= 0.05)
```

```
## [1] 99
```

```r
head(final, 10)
```

```
##               logFC interaction        PValue           FDR
## Cdkn1c   -7.8578999  -4.2247076 6.892576e-138 8.051218e-134
## Phlda2   -2.0816944  -0.9132656  8.184493e-50  4.780153e-46
## Hbb-bh1  -8.5188806  -3.5182419  6.550576e-42  2.550576e-38
## Acta2    -1.0044031  -1.9892717  7.801128e-39  2.278124e-35
## Asb4     -0.9044369  -1.6572831  3.084862e-27  7.206854e-24
## Xist    -11.8985704  -8.8029510  3.230173e-26  6.288608e-23
## Actc1    -1.0680425  -1.7986714  1.856506e-23  3.097977e-20
## Hba-x    -7.1421180  -2.3830414  1.012139e-22  1.477850e-19
## Tagln    -1.3879801  -1.9732654  1.871630e-16  2.429168e-13
## Hbb-y    -7.7658479  -2.5733576  2.416959e-16  2.823250e-13
```

# Visualizing the results

The `final` table contains a number of hemoglobin genes that one would not expect to be expressed in placental precursors.
Further examination indicates that the hemoglobins are downregulated in the KO beyond what would be expected from their loss in the ambient pool (Figure \@ref(fig:hemoplot)).
This suggests that there may be some hemoglobin expression in the placental precursors in the WT condition, 
possibly representing erythropoietic potential [@zeigler2006allantois;@corbel2007hematopoietic] that is lost upon knocking out _Tal1_.


```r
hemo <- head(grep("Hb[ab]-", rownames(final)), 3)
hemo <- rownames(final)[hemo]

# Computing relative to WT baseline.
as.cpms <- cpm(y.amb[hemo,], log=TRUE, prior.count=3)
as.cpms[,1:4] <- as.cpms[,1:4] - rowMeans(as.cpms[,3:4])
as.cpms[,5:8] <- as.cpms[,5:8] - rowMeans(as.cpms[,7:8])

par(mfrow=c(3,1))
for (i in hemo) {
    plot(as.cpms[i,], main=i, xlab="", cex=2, cex.lab=1.5,
        cex.main=1.5, ylab="Adjusted log-expression",
        pch=ifelse(ambient=="Y", 16, 4),
        col=ifelse(new.geno=="WT", "black", "red")
    )
}
```

<div class="figure">
<img src="embryo_expression_files/figure-html/hemoplot-1.png" alt="Relative log-expression of hemoglobin genes in the cells (crosses) and ambient pool (closed circles) for each sample, for the WT (black) and KO genotypes (red)." width="100%" />
<p class="caption">(\#fig:hemoplot)Relative log-expression of hemoglobin genes in the cells (crosses) and ambient pool (closed circles) for each sample, for the WT (black) and KO genotypes (red).</p>
</div>

Differences between the cell and ambient log-fold changes are more obvious in other genes.
A few muscle-related genes are downregulated in the KO condition in spite of their upregulation in the KO ambient pool (Figure \@ref(fig:muscleplot)).
This may indicate some loss of differentiation potential for muscle or connective tissue in placental precursors when _Tal1_ is lost.


```r
muscle <- c("Tagln", "Acta2", "Actc1")

# Computing relative to WT baseline.
as.cpms <- cpm(y.amb[muscle,], log=TRUE, prior.count=3)
as.cpms[,1:4] <- as.cpms[,1:4] - rowMeans(as.cpms[,3:4])
as.cpms[,5:8] <- as.cpms[,5:8] - rowMeans(as.cpms[,7:8])

par(mfrow=c(3,1))
for (i in muscle) {
    plot(as.cpms[i,], main=i, xlab="", cex=2, cex.lab=1.5,
        cex.main=1.5, ylab="Adjusted log-expression",
        pch=ifelse(ambient=="Y", 16, 4),
        col=ifelse(new.geno=="WT", "black", "red")
    )
}
```

<div class="figure">
<img src="embryo_expression_files/figure-html/muscleplot-1.png" alt="Relative log-expression of muscle-related genes in the cells (crosses) and ambient pool (closed circles) for each sample, for the WT (black) and KO genotypes (red)." width="100%" />
<p class="caption">(\#fig:muscleplot)Relative log-expression of muscle-related genes in the cells (crosses) and ambient pool (closed circles) for each sample, for the WT (black) and KO genotypes (red).</p>
</div>

For particularly interesting changes, it is always worthwhile to return to the per-cell expression profiles to inspect the magnitude of the effect.
Figure \@ref(percellplot) suggests that the downregulation of muscle-related genes is genuine but small relative to the cell-to-cell heterogeneity.
This is consistent with the small fold changes (barely 2-fold) between conditions reported in `final`.


```r
subsce <- sce[,sce$cluster==placenta]
plotExpression(subsce, x="sample", colour_by="tomato", 
    features=muscle, show_median=TRUE, ncol=3)
```

<div class="figure">
<img src="embryo_expression_files/figure-html/percellplot-1.png" alt="Distribution of log-transformed normalized expression values for the muscle-related genes in the placenta cluster. Each point represents a cell and is coloured according to tdTomato expression (i.e., positive for KO cells). The black bard represents the median across all cells in each sample." width="100%" />
<p class="caption">(\#fig:percellplot)Distribution of log-transformed normalized expression values for the muscle-related genes in the placenta cluster. Each point represents a cell and is coloured according to tdTomato expression (i.e., positive for KO cells). The black bard represents the median across all cells in each sample.</p>
</div>

# Concluding remarks 

The DE analysis shown above can be repeated for each cluster, provided there are enough cells in that cluster from each sample to obtain a summed count matrix.
This complements the differential abundance analysis by providing another perspective into the differences between conditions.
Indeed, examination of the DEGs may suggest that a single cluster should actually be treated as two separate cell types that have been merged together by `fastMNN()`.
(Conversely, examination of markers for two differentially abundant clusters may suggest that they should actually be a single cell type with DE between conditions.)

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
##  [1] BiocFileCache_1.7.0         dbplyr_1.3.0               
##  [3] edgeR_3.25.3                limma_3.39.12              
##  [5] scater_1.11.11              ggplot2_3.1.0              
##  [7] SingleCellExperiment_1.5.2  SummarizedExperiment_1.13.0
##  [9] DelayedArray_0.9.8          BiocParallel_1.17.15       
## [11] matrixStats_0.54.0          Biobase_2.43.1             
## [13] GenomicRanges_1.35.1        GenomeInfoDb_1.19.2        
## [15] IRanges_2.17.4              S4Vectors_0.21.10          
## [17] BiocGenerics_0.29.1         BiocStyle_2.11.0           
## 
## loaded via a namespace (and not attached):
##  [1] httr_1.4.0               viridis_0.5.1           
##  [3] BiocSingular_0.99.12     bit64_0.9-7             
##  [5] viridisLite_0.3.0        splines_3.6.0           
##  [7] DelayedMatrixStats_1.5.2 assertthat_0.2.0        
##  [9] statmod_1.4.30           BiocManager_1.30.4      
## [11] highr_0.7                blob_1.1.1              
## [13] GenomeInfoDbData_1.2.0   vipor_0.4.5             
## [15] yaml_2.2.0               RSQLite_2.1.1           
## [17] pillar_1.3.1             lattice_0.20-38         
## [19] glue_1.3.0               digest_0.6.18           
## [21] XVector_0.23.0           colorspace_1.4-0        
## [23] cowplot_0.9.4            htmltools_0.3.6         
## [25] Matrix_1.2-16            plyr_1.8.4              
## [27] pkgconfig_2.0.2          bookdown_0.9            
## [29] zlibbioc_1.29.0          purrr_0.3.0             
## [31] scales_1.0.0             processx_3.2.1          
## [33] tibble_2.0.1             withr_2.1.2             
## [35] lazyeval_0.2.1           magrittr_1.5            
## [37] crayon_1.3.4             memoise_1.1.0           
## [39] evaluate_0.13            ps_1.3.0                
## [41] beeswarm_0.2.3           tools_3.6.0             
## [43] stringr_1.4.0            munsell_0.5.0           
## [45] locfit_1.5-9.1           irlba_2.3.3             
## [47] callr_3.1.1              compiler_3.6.0          
## [49] rsvd_1.0.0               rlang_0.3.1             
## [51] grid_3.6.0               RCurl_1.95-4.11         
## [53] BiocNeighbors_1.1.11     rappdirs_0.3.1          
## [55] labeling_0.3             bitops_1.0-6            
## [57] rmarkdown_1.11           gtable_0.2.0            
## [59] curl_3.3                 DBI_1.0.0               
## [61] R6_2.4.0                 gridExtra_2.3           
## [63] knitr_1.21               dplyr_0.8.0.1           
## [65] bit_1.1-14               stringi_1.3.1           
## [67] ggbeeswarm_0.6.0         compareSingleCell_0.98.0
## [69] Rcpp_1.0.0               tidyselect_0.2.5        
## [71] xfun_0.5
```

# References

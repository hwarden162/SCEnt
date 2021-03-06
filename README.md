
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SCEnt

<!-- badges: start -->
<!-- badges: end -->

SCEnt is a package for single cell entropy analysis. It can calculate
metrics for the heterogeneity or homogeneity of a gene within a cell
population. It can use these metrics to perform feature selection for
scRNA-seq data.

From the work of Michael J. Casey, Ruben J. Sanchez-Garcia and Ben D.
MacArthur (2020) &lt;doi:
<https://doi.org/10.1101/2020.10.01.322255>&gt;. Package written by Hugh
Warden.

## Installation

SCEnt is available to install from CRAN:

``` r
install.packages("SCEnt")
```

And the development version of SCEnt is available from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("hwarden162/SCEnt")
```

## Entropy Analysis

The expression of each gene can be considered as a probability
distribution. Where given the information that a gene has been
expressed, the expression counts can be used to infer the probability
that it was expressed in a given cell.

Each of these probability distributions for each gene contiain a certain
amount of information, known in mathematics as entropy. Entropy has been
historically described as the amount of ‘surprise’ encoded in the
system. If a given gene is only ever expressed in one cell then this
gene is said to have low entropy, as if that gene is expressed then the
cell expressing it is known and there is no ‘surprise’. However, if a
given gene is expressed equally in every cell, then this gene will have
high entropy. As expression of this gene does not narrow down the amount
of cells it could have come from, meaning the answer will be a
‘surprise’.

What this means, is that the more homogeneously a gene is expressed
within a cell population, the more entropy the expression distrbution
will have. Therefore, allowing the entropy to be a metric of homogeneity
within the cell population.

Every probability distribution has an entropy value and there are many
ways to compare entropy values between probability distributions. One of
these is called the Kullback-Liebler Divergence (or KL Divergence for
short). KL Divergence measures the amount of entropy that would be lost
if one distribution were used to approximate another.

The KL Divergence can then be used in this case to see how much
information is lost if a uniform distribution were used to represent the
expression distribution of a gene. Remembering that the uniform
distribution would be the expression of a completely homogeneous gene,
the KL Divergence gives us a measure of how much extra information we
are getting from the heterogeneity of the gene. Essentially giving a
metric for the heterogeneity of a gene within a cell population.

This is intended to be a very brief overview of the underpinning
mathematics of SCEnt. For a more detailed explanation please see the
paper [Measuring the Information Obtained from a Single-Cell Sequencing
Experiment](https://www.biorxiv.org/content/biorxiv/early/2020/10/01/2020.10.01.322255.full.pdf).

## Using SCEnt to Quantify Homogeneity and Heterogeneity

Here is some synthetic scRNA-seq data:

``` r
gene_counts
#>       gene1 gene2 gene3 gene4 gene5
#> cell1     0     5     2     3     0
#> cell2     0     5     0     3     0
#> cell3     0     3     2     3     0
#> cell4     0     2     1     3     0
#> cell5     1     0     3     3     5
#> cell6     2     0     0     3     0
#> cell7     3     0     1     3     0
```

The `gene_hom()` and `gene_het()` functions can be used to calculate the
homogeneity or heterogeneity of a gene, respectively. Each of these can
be passed a gene and it will return a value.

``` r
(gene1 <- gene_counts[,1])
#> cell1 cell2 cell3 cell4 cell5 cell6 cell7 
#>     0     0     0     0     1     2     3
gene_hom(gene1)
#> [1] 1.459148
gene_het(gene1)
#> [1] 1.348207

(gene2 <- gene_counts[,2])
#> cell1 cell2 cell3 cell4 cell5 cell6 cell7 
#>     5     5     3     2     0     0     0
gene_hom(gene2)
#> [1] 1.908613
gene_het(gene2)
#> [1] 0.8987422

(gene3 <- gene_counts[,3])
#> cell1 cell2 cell3 cell4 cell5 cell6 cell7 
#>     2     0     2     1     3     0     1
gene_hom(gene3)
#> [1] 2.19716
gene_het(gene3)
#> [1] 0.6101952

(gene4 <- gene_counts[,4])
#> cell1 cell2 cell3 cell4 cell5 cell6 cell7 
#>     3     3     3     3     3     3     3
gene_hom(gene4)
#> [1] 2.807355
gene_het(gene4)
#> [1] 0

(gene5 <- gene_counts[,5])
#> cell1 cell2 cell3 cell4 cell5 cell6 cell7 
#>     0     0     0     0     5     0     0
gene_hom(gene5)
#> [1] 0
gene_het(gene5)
#> [1] 2.807355
```

Note that gene4 is a maximally homogeneous gene and gene5 is a maximally
heterogeneous gene. Therefore, their entropy calculations always return
the extreme values of 0 or 2.807355, the non-zero extreme value will be
equal to log<sub>2</sub>\[*N*\] where *N* is the number of cells in the
sample.

Rather than submitting each gene individually, the whole matrix of gene
expressions can be passed. However, the data needs to be in the format
of having genes represented as rows and cells represented as columns.

``` r
gene_hom(t(gene_counts))
#>    gene1    gene2    gene3    gene4    gene5 
#> 1.459148 1.908613 2.197160 2.807355 0.000000
```

Rather than transposing the matrix as an input, there is also a
transpose parameter for convenience

``` r
gene_het(gene_counts, transpose = TRUE)
#>     gene1     gene2     gene3     gene4     gene5 
#> 1.3482070 0.8987422 0.6101952 0.0000000 2.8073549
```

The `gene_hom()` and `gene_het()` functions have other parameters not
covered in detail here. The `unit` parameter is passed to the entropy
calculations to change the units by which entropy is calclated, by
default this is set to `"log2"` such that all outputs are in bits. The
`normalise` parameter is there for if the gene expression is already in
probability form. This parameter should generally not be used except in
some cases to speed up calculations.

## Using SCEnt For Feature Selection

When wanting to use gene expression for prediction, homogeneously
expressed genes within the population are not going to be useful for any
sort of classification. Conversely, heterogeneously expressed genes are
much more likely to be useful to make predictions from.

The function `scent_select()` will carry out feature selection on
scRNA-seq data by finding gene heterogeneity and applying some user
defined threshold. The thresholds that can be selected are
`bit_threshold`, `count_threshold` and `perc_threshold`. `bit_threshold`
takes a numeric value and SCEnt will only select genes with a
heterogenity value greater than the given bit value. `count_threshold`
takes an integer value, SCEnt will return the data with only the top *n*
heteogeneously expressed genes, where *n* is the `count_threshold`.
`perc_threshold` takes a value between 0 and 1, SCEnt will only return
genes with heterogeneity greater than the percentile of population
heterogeneity as given by `perc_threshold`. Unlike the entropy
calculations above, the matrix of gene expressions should be in the
format of cells as rows and genes as columns. Again, there is a
`transpose` option built in, as well as the `unit` and `normalise`
parameters to be passed to the entropy calculations.

Here is the sample data again:

``` r
gene_counts
#>       gene1 gene2 gene3 gene4 gene5
#> cell1     0     5     2     3     0
#> cell2     0     5     0     3     0
#> cell3     0     3     2     3     0
#> cell4     0     2     1     3     0
#> cell5     1     0     3     3     5
#> cell6     2     0     0     3     0
#> cell7     3     0     1     3     0
```

`scent_select()` can carry out each of these three forms of feature
selection depending on which threshold is supplied.

``` r
scent_select(gene_counts, bit_threshold = 0.85)
#>       gene1 gene2 gene5
#> cell1     0     5     0
#> cell2     0     5     0
#> cell3     0     3     0
#> cell4     0     2     0
#> cell5     1     0     5
#> cell6     2     0     0
#> cell7     3     0     0
scent_select(gene_counts, count_threshold = 2)
#>       gene1 gene5
#> cell1     0     0
#> cell2     0     0
#> cell3     0     0
#> cell4     0     0
#> cell5     1     5
#> cell6     2     0
#> cell7     3     0
scent_select(gene_counts, perc_threshold = 0.25)
#>       gene1 gene2 gene3 gene5
#> cell1     0     5     2     0
#> cell2     0     5     0     0
#> cell3     0     3     2     0
#> cell4     0     2     1     0
#> cell5     1     0     3     5
#> cell6     2     0     0     0
#> cell7     3     0     1     0
```

Trying to threshold on multiple values will throw an error

``` r
scent_select(gene_counts, bit_threshold = 0.85, count_threshold = 2)
#> Error in scent_select(gene_counts, bit_threshold = 0.85, count_threshold = 2): 
#>  Only one threshold can be set at a time
```

This is to avoid problems with orders of operation, and instead
`scent_select()` should be piped back into itself so that the order in
which the thresholds are applied is explicit

``` r
gene_counts %>%
  scent_select(bit_threshold = 0.85) %>%
  scent_select(count_threshold = 2)
#>       gene1 gene5
#> cell1     0     0
#> cell2     0     0
#> cell3     0     0
#> cell4     0     0
#> cell5     1     5
#> cell6     2     0
#> cell7     3     0
```

## Tidy Implementation of Feature Selection

There is a tidy wrapper for `scent_select()` called
`scent_select_tidy()`. This is currently being used to develop a new
recipe step for use within the `tidymodels` framework.

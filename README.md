
<!-- README.md is generated from README.Rmd. Please edit that file -->

# beer <img src='man/figures/logo.png' align="right" height="138" />

<!-- badges: start -->
<!-- badges: end -->

Phage immuno-precipitation sequencing (PhIP-seq) is a high-throughput
approach for characterizing antibody responses to a variety of target
antigens. A typical component of PhIP-seq analyses involves identifying
which peptides elicit enriched antibody responses. `beer` provides two
approaches for identifying peptide enrichments.

The first approach is based on
[`edgeR`](https://bioconductor.org/packages/release/bioc/html/edgeR.html)’s
standard pipeline for identifying differential expression from read
count data[^1][^2][^3]. Though
[`edgeR`](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
is remarkably effective at quickly identifying enriched antibody
responses, it is less likely to pick up enriched peptides at the lower
fold-change range.

The second approach, Bayesian Estimation in R (BEER) was developed
specifically for the PhIP-seq setting and implements a Bayesian model to
identify peptide enrichments as described in Chen et. al[^4]. Though
BEER is more likely to identify enriched peptides at the lower
fold-change range, it tends to take much longer to run.

Below we give a brief overview of the two approaches. For more
information, see the package vignette using `browseVignettes("beer")`.
Both methods can be run in synchronously or asynchronously as supported
by
[`future`](https://cran.r-project.org/web/packages/future/index.html).

## Installation

### `rjags`

For Bayesian MCMC modeling, `beer` relies on
[`rjags`](https://cran.r-project.org/web/packages/rjags/index.html) to
interface [Just Another Gibbs Sampler
(JAGS)](https://mcmc-jags.sourceforge.io/). JAGS can be downloaded from
[this link](https://sourceforge.net/projects/mcmc-jags/files/).
[Homebrew](https://brew.sh/) users can install JAGS using,

    brew install jags

For M1 Mac users using Rosetta emulation of intel, Homebrew installation
of JAGS will likely work. However, we recommend installing JAGS from
source for all other M1 Mac users.

Once JAGS has been installed,
[`rjags`](https://cran.r-project.org/web/packages/rjags/index.html) can
be installed in `R` via `install.packages("rjags")`.

### `beer`

Once `rjags` has been installed, the stable release version of `beer` in
Bioconductor can be installed using `BiocManager`:

``` r
if (!require("BiocManager"))
    install.packages("BiocManager")
    
BiocManager::install("beer")
```

To load the package:

``` r
library(beer)
```

## edgeR

Differentially enriched peptides between a particular serum sample and
all beads-only samples indicate enriched antibody responses to those
peptides. Thus, to identify enriched peptides, we can run the standard
[`edgeR`](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
pipeline for differential expression.

The `runEdgeR()` function estimates peptide-specific dispersion
parameters then runs the exact test proposed by Robinson and Smyth[^5]
for the difference in mean between two groups of negative binomial
random variables. Since peptides are enriched only if average proportion
of reads pulled in the serum sample is higher than the average
proportion of reads pulled in a beads-only samples, two-sided p-values
are converted to one-sided p-values.

``` r
## Load data
data_path <- system.file("extdata/sim_data.RDS", package = "beer")
sim_data <- readRDS(data_path)
```

``` r
edgeR_out <- runEdgeR(sim_data, 
                   assay.names = c(logfc = "edgeR_logfc", 
                                   prob = "edgeR_logpval"))
```

Using BH correction to adjust for multiple testing, enriched peptides
are given by the matrix,

``` r
assay(edgeR_out, "edgeR_hits") <- apply(
  assay(edgeR_out, "edgeR_logpval"), 2, 
  function(sample){
    pval <- 10^(-sample)
    p.adjust(pval, method = "BH") < 0.05
  })

colSums(assay(edgeR_out, "edgeR_hits"))
#>  1  2  3  4  5  6  7  8  9 10 
#> NA NA NA NA  2  1  1  1  2  0
```

## BEER (Bayesian Estimation Enrichment in R)

BEER uses a Bayesian hierarchical model to derive posterior
probabilities of enrichment and estimated fold-changes. Briefly, each
sample is run individually in comparison to all beads-only samples as
follows:

1.  **Define prior parameters.** Though most prior parameters are
    supplemented by the user (or use the defaults), prior parameters for
    non-enriched peptides are first approximated using all beads-only
    samples.
2.  **Identify super enriched peptides.** Based on the prior parameters,
    super enriched peptides are first excluded as these peptides should
    always have posterior probabilities of enrichment of 1.
3.  **Re-estimate beads-only prior parameters.** Prior parameters are
    then re-estimated from the beads-only samples for the remaining
    peptides.
4.  **Initialize and run the MCMCs.** To reduce convergence time, MLE
    estimates are used to initialize the MCMC sampler, and samples are
    drawn from the posterior distributions of the unknown parameters.
5.  **Summarize and store results.** Posterior samples are summarized
    using the means of the posterior distribution and are stored in the
    PhIPData object.

BEER can be easily run with `brew()`:

``` r
## Named vector specifying where we want to store the summarized MCMC output
## NULL indicates that the output should not be stored.
assay_locations <- c(
  phi = "beer_fc_marg", 
  phi_Z = "beer_fc_cond", 
  Z = "beer_prob", 
  c = "sampleInfo", 
  pi = "sampleInfo"
)

beer_out <- brew(sim_data, assay.names = assay_locations)
```

Thus, supposing peptides with posterior probability above 0.5 are
enriched and noting that super enriched peptides were not run (and thus
are missing entries in the posterior probability matrix), the matrix of
enriched peptides is given by,

``` r
## Define matrix of peptides that were run in BEER
was_run <- matrix(rep(beer_out$group != "beads", each = nrow(beer_out)), 
                  nrow = nrow(beer_out))

## Identify super-enriched peptides
## These peptides were in samples that were run, but have missing posterior 
## probabilities
are_se <- was_run & is.na(assay(beer_out, "beer_prob"))

## Enriched peptides are peptides with:
## - posterior probability > 0.5, OR
## - super-enriched peptides
assay(beer_out, "beer_hits") <- assay(beer_out, "beer_prob") > 0.5 | are_se

colSums(assay(beer_out, "beer_hits"))
#>  1  2  3  4  5  6  7  8  9 10 
#> NA NA NA NA  5  7  3  3  3  2
```

## Beads-only round robin

To approximate the false positive rate, we often run each of the
beads-only samples against all other beads-only samples. This beads-only
round robin also provides a sense of how similar the beads-only samples
are to each other.

The beads-only round robin can be included in `brew()` and `runEdgeR()`
by specifying `beadsRR = TRUE`.

``` r
## edgeR with beadsRR
edgeR_beadsRR <- runEdgeR(sim_data, beadsRR = TRUE, 
                          assay.names = c(logfc = "edgeR_logfc", 
                                          prob = "edgeR_logpval"))
## Calculate hits
assay(edgeR_beadsRR, "edgeR_hits") <- apply(
  assay(edgeR_beadsRR, "edgeR_logpval"), 2, 
  function(sample){
    pval <- 10^(-sample)
    p.adjust(pval, method = "BH") < 0.05
  })

## Note samples 1-4 have 0 instead of NA now
colSums(assay(edgeR_beadsRR, "edgeR_hits"))
#>  1  2  3  4  5  6  7  8  9 10 
#>  0  0  0  0  2  1  1  1  2  0
```

``` r
## BEER with beadsRR added to edgeR output
beer_beadsRR <- brew(edgeR_beadsRR, beadsRR = TRUE, 
                     assay.names = assay_locations)

## Check BEER hits like before
was_run <- matrix(rep(beer_beadsRR$group != "beads", each = nrow(beer_beadsRR)), 
                  nrow = nrow(beer_beadsRR))
are_se <- was_run & is.na(assay(beer_beadsRR, "beer_prob"))
beer_hits <- assay(beer_beadsRR, "beer_prob") > 0.5 | are_se

## Note again that samples 1-4 are not NA 
colSums(beer_hits)
#>  1  2  3  4  5  6  7  8  9 10 
#>  2  3  0  0  5  7  3  3  3  2
```

Alternatively, one can run `beadsRR()` separately,

``` r
## edgeR with beadsRR
edgeR_beadsRR <- beadsRR(sim_data, method = "edgeR", 
                         assay.names = c(logfc = "edgeR_logfc", 
                                         prob = "edgeR_logpval"))
## Calculate hits
assay(edgeR_beadsRR, "edgeR_hits") <- apply(
  assay(edgeR_beadsRR, "edgeR_logpval"), 2, 
  function(sample){
    pval <- 10^(-sample)
    p.adjust(pval, method = "BH") < 0.05
  })

## Note samples 5-10 are NA now
colSums(assay(edgeR_beadsRR, "edgeR_hits"))
#>  1  2  3  4  5  6  7  8  9 10 
#>  0  0  0  0 NA NA NA NA NA NA
```

## References

[^1]: Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a
    Bioconductor package for differential expression analysis of digital
    gene expression data. Bioinformatics 26, 139-140

[^2]: McCarthy DJ, Chen Y and Smyth GK (2012). Differential expression
    analysis of multifactor RNA-Seq experiments with respect to
    biological variation. Nucleic Acids Research 40, 4288-4297

[^3]: Chen Y, Lun ATL, Smyth GK (2016). From reads to genes to pathways:
    differential expression analysis of RNA-Seq experiments using
    Rsubread and the edgeR quasi-likelihood pipeline. F1000Research 5,
    1438

[^4]: Chen A, Kammers K, Larman HB, Scharpf R, Ruczinski I. Detecting
    antibody reactivities in phage immunoprecipitation sequencing data
    (2022). *bioRxiv*.
    <https://www.biorxiv.org/content/10.1101/2022.01.19.476926v1>

[^5]: Robinson MD and Smyth GK. Small-sample estimation of negative
    binomial dispersion, with applications to SAGE data (2008).
    *Biostatistics*, 9, 321-332.
    <https://doi.org/10.1093/biostatistics/kxm030>

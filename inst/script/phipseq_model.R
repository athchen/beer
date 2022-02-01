#' Description of phipseq_model.bugs
#'
#' phipseq_model.bugs encodes the Bayesian model for identifying enriched
#' peptides in PhIP-Seq data as described in Chen et al. The model is written
#' in JAGS. Each sample is run individually (rather than on a per-plate basis).
#'
#' Briefly, the model describes the observed read counts using a Binomial
#' distribution. The probability of drawing a read depends on (1) the
#' attenuation constant, (2) the underlying fold-change in comparison to
#' negative controls or beads-only samples, and (3) the distribution of a
#' peptide in a beads-only sample.
#'
#' The underlying fold-change depends on whether the peptide of interest is
#' enriched in the particular sample. If the peptide is enriched, then the fold-
#' change is modeled using a shifted Gamma distribution. If the peptide is not
#' enriched, then the fold-change is 1. The attenuation constant and likelihood
#' of being enriched are modeled using Beta distributions.
#'
#' Default prior parameters can be found in the \code{\link{brew}} documentation.
#'
#' @seealso  \code{\link{brew}}
#' @seealso [Detecting Antibody Reactivities in Phage ImmunoPrecipitation Sequencing Data](https://www.biorxiv.org/content/10.1101/2022.01.19.476926v1)

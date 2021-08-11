#' Guess super-enriched peptides based on edgeR fold-change estimates
#'
#' @param object \code{\link[PhIPData]{PhIPData}} object.
#' @param threshold minimum estimated fc for a peptide to be considered
#' super-enriched.
#' @param fc.name assay name corresponding to the assay that stores the edgeR
#' estimated log2 fold-changes.
#'
#' @return logical matrix of the with the same dimensions as \code{object}
#' indicating which peptides are considered super-enriched.
.guessEnriched_edgeR <- function(object, threshold = 15, fc.name = "logfc"){

    ## Check that assay is present
    if(!fc.name %in% assayNames(object)){
        stop(paste0(fc.name, " assay does not exist in the PhIPData object."))
    }

    ## Check that assay is not empty
    if(all(is.na(assay(object, fc.name)))){
        stop(paste0(fc.name, " is empty. edgeR estimated fold-changes can be ",
                    "added to the object using the `edgeR()` function."))
    }

    assay(object, fc.name) > log2(threshold)
}

#' Guess enriched peptides based on MLE estimates of the true fold-change
#'
#' @param object \code{\link[PhIPData]{PhIPData}} object.
#' @param threshold minimum estimated fc for a peptide to be considered
#' super-enriched.
#' @param beads.prior data.frame of prior parameters for beads-only samples.
#'
#' @return logical matrix of the with the same dimensions as \code{object}
#' indicating which peptides are considered super-enriched.
.guessEnriched_MLE <- function(object, threshold = 15, beads.prior){

    n <- librarySize(object)

    ## Calculate expected proportion of reads pulled based on beads-only samples
    expected_prop <- beads.prior[["a_0"]]/
        (beads.prior[["a_0"]] + beads.prior[["b_0"]])
    expected_rc <- vapply(n, function(n_i) n_i*expected_prop,
                          numeric(length(expected_prop)))

    ## Calculate fold change based on an attenuation constant of 1
    fc_beads <- counts(object)/expected_rc

    ## Rough estimate of attenuation constant
    ## Require an observed FC of 5 (w/o accounting for the attenuation constant)
    guess_e <- apply(fc_beads, c(1, 2), function(x) ifelse(x > 5, 1, 0 ))
    c_est <- vapply(seq(ncol(guess_e)), function(col) {
        if(object$group[col] != getBeadsName()){
            ne_peps <- !guess_e[, col]
            coef(lm(counts(object)[ne_peps, col] ~
                        expected_rc[ne_peps, col] - 1))
        } else 1
    }, numeric(1))

    ## Guess enriched peptides based on estimated attenuation constant
    expected_rc_attn <- vapply(n*c_est, function(n_i) n_i*expected_prop,
                               numeric(length(expected_prop)))
    fc_attn <- counts(object)/expected_rc_attn

    fc_attn > threshold
}

#' @title Identifying clearly enriched peptides
#'
#' @description As clearly enriched peptides will always have a 100\%
#' posterior probability of enrichment, BEER removes these peptides a priori to
#' running the model. Clearly enriched peptides can be identified using edgeR
#' estimated fold-changes or maximum likelihood estimates based on the specified
#' prior parameters. Additional parameters for each method can be found in the
#' details below.
#'
#' @details \strong{edgeR}. Identification of clearly enriched peptides relies
#' on edgeR fold-change estimates, so \code{\link[beer]{edgeR}} must be run on
#' the \code{\link[PhIPData]{PhIPData}} object beforehand. Additional parameters
#' for identifying clearly enriched peptides based on edgeR estimated
#' fold-changes are listed below:
#'
#' \itemize{
#'      \item \code{object}: a \code{\link[PhIPData]{PhIPData}} object.
#'      \item \code{threshold}: minimum estimated fc for a peptide to be
#'      considered super-enriched. The default value is 15.
#'      \item \code{fc.name}: assay name corresponding to the assay that stores
#'      the edgeR estimated log2 fold-changes.
#' }
#'
#' \strong{MLE}. As the number of reads tends to be quite large, the estimates
#' for the proportion of reads pulled are generally accurate. Clearly enriched
#' peptides are identified by first comparing the observed read count to the
#' expected read count based on the beads-only prior parameters. Peptides with
#' observed read counts larger than 5 times the expected read counts are
#' temporarily labeled as enriched, and attenuation constants are estimated by
#' regressing the observed read counts on the expected read counts for all
#' non-enriched peptides. Using this attenuation constant, peptides with
#' fold-changes above some predefined threshold after adjusting for the
#' attenuation constant are considered enriched. Parameters for identifying
#' clearly enriched peptides using the MLE approach are listed below.
#'
#' \itemize{
#'      \item \code{object}: a \code{\link[PhIPData]{PhIPData}} object.
#'      \item \code{threshold}: minimum estimated fc for a peptide to be
#'      considered super-enriched.
#'      \item \code{beads.prior}: data.frame of prior parameters for beads-only
#'      samples.
#' }
#'
#' @param object a \code{\link[PhIPData]{PhIPData}} object
#' @param method one of "mle" or "edgeR", specifying which method to use to
#' identify clearly enriched peptides
#' @param ... additional parameters dependent on the method used. See details
#' for more information
#'
#' @return a logical matrix of the with the same dimensions as \code{object}
#' indicating which peptides are considered super-enriched.
#'
#' @examples
#' sim_data <- readRDS(system.file("extdata", "sim_data.rds", package = "beer"))
#' edgeR_out <- edgeR(sim_data)
#'
#' guessEnriched(edgeR_out, method = "edgeR", threshold = 15, fc.name = "logfc")
#' guessEnriched(edgeR_out, method = "mle", threshold = 15,
#'     beads.prior = getAB(edgeR_out, method = "edgeR"))
#'
#' @export
guessEnriched <- function(object, method = "mle", ...){

    if(method == "mle"){ .guessEnriched_MLE(object, ...)
    } else if (method == "edgeR") {
        .guessEnriched_edgeR(object, ...)
    } else {
        stop("Invalid method. Valid methods include 'mle' and 'edgeR'")
    }
}

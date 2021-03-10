#' @export
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

    assay(object, fc.name) > threshold
}

#' @export
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

#' @export
guessEnriched <- function(object, method = "mle", ...){

    if(method == "mle"){ .guessEnriched_MLE(object, ...)
    } else .guessEnriched_edgeR(object, ...)
}

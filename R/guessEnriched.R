guessEnriched_edgeR <- function(object, threshold = 15, fc.name = "logfc"){

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

guessEnriched_MLE <- function(object, threshold = 15, beads.prior){

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
    guess_e <- apply(fc_beads, c(1, 2), function(x) if(x > 5) 0 else 1)
    new_n <- colSums(guess_e*counts(object))
    guess_e <- which(fc_beads > 1.5)
    c_est <- new_n/n

    ## Guess enriched peptides based on estimated attenuation constant
    expected_rc_attn <- vapply(new_n, function(n_i) n_i*expected_prop,
                               numeric(length(expected_prop)))
    fc_attn <- counts(object)/expected_rc_attn

    fc_attn > threshold
}

guessEnriched <- function(object, ..., method = "mle"){

    if(method == "edgeR"){ guessEnriched_edgeR(object, ...)
    } else guessEnriched_MLE(object, ...)
}

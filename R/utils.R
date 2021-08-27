#' Function to check that the counts matrix only contains integers
#'
#' @param object PhIPData object
#'
#' @return nothing if all counts are integers, and error otherwise
#'
#' @import PhIPData
.checkCounts <- function(object) {
    is_int <- vapply(PhIPData::counts(object), is.integer, logical(1))
    if (!all(is_int)) {
        stop("`counts` entries must be integers.")
    }
}

#' Function to check whether an assay will be overwritten
#'
#' If the an assay is not specified (e.g. with \code{NA}), then
#' \code{.checkOverwrite()} will return \code{FALSE} (rather than NA).
#'
#' @param object PhIPData object
#' @param assay.names character vector of assay names
#'
#' @return logical vector indicating whether data in an assay will be
#' overwritten
#'
#' @import PhIPData SummarizedExperiment
.checkOverwrite <- function(object, assay.names) {
    ## check sampleInfo objects
    in_sample <- ifelse(assay.names == "sampleInfo",
        names(assay.names) %in% colnames(sampleInfo(object)),
        FALSE
    )

    in_assay <- vapply(assay.names, function(name) {
        if (name %in% assayNames(object)) {
            any(!is.na(assay(object, name)))
        } else {
            FALSE
        }
    }, logical(1))

    ifelse(is.na(in_sample | in_assay), FALSE, in_sample | in_assay)
}

#' Function to tidy inputs pertaining to parallelization
#'
#' This function ensures that the specified plan is in the list of valid
#' plans for \code{[future::plan]}. The function also ensures that
#' the maximum number of workers is one less that number of available number
#' of cpus.
#'
#' @param parallel list of parameters passed to the
#' \code{[future::plan]} function.
#'
#' @return tidied list of parameters passed to \code{[future::plan]}
#'
#' @importFrom future availableCores
.tidyParallel <- function(parallel) {
    valid_strat <- c(
        "sequential", "transparent", "multisession",
        "multicore", "cluster", "remote"
    )

    ## Convert to list if parallel is a character vector
    ## The first entry always defines the plan
    if (is.vector(parallel)) {
        strat <- ifelse("strategy" %in% names(parallel),
            parallel[["strategy"]],
            parallel[1]
        )
        workers <- if ("workers" %in% names(parallel)) {
            parallel[["workers"]]
        } else {
            NULL
        }
        parallel <- list(strategy = strat, workers = as.numeric(workers))
    }

    ## Check that strategy is valid
    if (!all(parallel[["strategy"]] %in% valid_strat)) {
        stop(paste0(
            "The specified parallel strategy is not supported. ",
            "Valid strategies are '",
            paste0(valid_strat, collapse = "', '"), "'."
        ))
    }

    ## Get number of workers, ensure that the number of workers is at least one
    ## and less than the number of available cores - 1
    if (parallel[["strategy"]] == "sequential") {
        parallel[["workers"]] <- NULL
    } else {
        parallel[["workers"]] <- max(1, min(
            future::availableCores() - 1,
            parallel[["workers"]]
        ))
    }

    parallel
}

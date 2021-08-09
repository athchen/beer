#' Function to check that the counts matrix only contains integers
#'
#' @param object PhIPData object
#' @return nothing if all counts are integers, and error otherwise
.checkCounts <- function(object){
    is_int <- vapply(PhIPData::counts(object), is.integer, logical(1))
    if(!all(is_int)){
        stop("`counts` entries must be integers.")
    }
}

#' Function to tidy inputs pertaining to parallelization
#'
#' This function ensures that the specified plan is in the list of valid
#' plans for \code{\link{future::plan}{future}}. The function also ensures that
#' the maximum number of workers is one less that number of available number
#' of cpus.
#'
#' @param parallel list of parameters passed to the
#' \code{\link{future::plan}{future::plan()}} function.
#'
#' @importFrom future availableCores
.tidyParallel <- function(parallel){
    valid_strat <- c("sequential", "transparent", "multisession",
                     "multicore", "cluster", "remote")

    if(is.null(parallel[["strategy"]])) parallel[["strategy"]] <- "sequential"

    ## Check that strategy is valid
    if (!all(parallel[["strategy"]] %in% valid_strat)) {
        stop(paste0("The specified parallel strategy is not supported. ",
                    "Valid strategies are '",
                    paste0(valid_strat, collapse = "', '"), "'."))
    }

    ## Get number of workers, ensure that the number of workers is at least one
    ## and less than the number of available cores - 1
    if(parallel[["strategy"]] == "sequential"){
        parallel[["workers"]] <- NULL
    } else {
        parallel[["workers"]] <- max(1, min(future::availableCores() - 1,
                                            parallel[["workers"]]))
    }

    parallel
}

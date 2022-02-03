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


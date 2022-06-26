#' Calculate expected read counts or proportion of reads
#'
#' @param object PhIPData object
#' @param type one of `rc` or `prop` indicating whether the function should
#' return the expected read counts or expected proportion of reads, respectively
#' @param assay.name name indicating where the results should be stored in the
#' PhIPData object
#'
#' @return PhIPData object with the results stored in the location specified
#' by \code{assay.name}.
#'
#' @examples
#' sim_data <- readRDS(system.file("extdata", "sim_data.rds", package = "beer"))
#'
#' ## Calculate expected read counts
#' getExpected(sim_data, "rc", "expected_rc")
#'
#' ## Calculate expected proportion of reads
#' getExpected(sim_data, "prop", "expected_prop")
#'
#' @importFrom cli cli_alert_warning
#' @export
getExpected <- function(object, type = "rc", assay.name = "expected_rc") {

    ## Check inputs
    if (!type %in% c("rc", "prop")) {
        stop(
            "Invalid `type` specified. `type` must be one of `rc` or `prop`."
        )
    }

    ## Check assay overwrite
    if (.checkOverwrite(object, assay.name)) {
        cli::cli_alert_warning("The specified assay exists and will be overwritten.")
    }

    beads_prop <- rowMeans(propReads(subsetBeads(object)))
    prop_mat <- matrix(
        rep(beads_prop, each = ncol(object)),
        nrow = length(beads_prop),
        byrow = TRUE
    )

    lib_mat <- matrix(
        rep(librarySize(object), each = ncol(object)),
        nrow = nrow(object)
    )

    assay_out <- if (type == "rc") {
        lib_mat * prop_mat
    } else {
        prop_mat
    }
    colnames(assay_out) <- colnames(object)
    rownames(assay_out) <- rownames(object)

    ## Store in phipdata_obj
    assay(object, assay.name) <- assay_out

    ## Return object
    object
}

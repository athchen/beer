#' Calculate expected read counts or proportion of reads
#'
#' @param object PhIPData object
#' @param type any of `rc` or `prop` indicating whether the function should
#' return the expected read counts or expected proportion of reads, respectively
#' @param assay.names name(s) indicating where the results should be stored in the
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
#' ## Calculate both
#' getExpected(sim_data)
#'
#' @importFrom cli cli_alert_warning
#' @export
getExpected <- function(object, type = c("rc", "prop"), 
    assay.names = c("expected_rc", "expected_prop")) {

    ## Check inputs
    if(all(!c("rc", "prop") %in% type)){
        stop("Invalid `type` specified. `type` must be any of `rc` or `prop`.")
    }

    ## Check assay overwrite
    overwrite <- vapply(assay.names, .checkOverwrite, logical(1), object = object)
    overwrite_assays <- names(overwrite)[overwrite]
    if (any(overwrite)) {
        cli::cli_alert_warning(
            paste0(
                "{.overwrite_assays {overwrite_assays}} assay{?s} {?exists/exist} ", 
                "and will be overwritten."
            )
        )
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
    expected_rc <- lib_mat * prop_mat
    
    ## Tidy names
    colnames(prop_mat) <- colnames(expected_rc) <- colnames(object)
    rownames(prop_mat) <- rownames(expected_rc) <- rownames(object)

    ## Store in phipdata_obj
    if("prop" %in% type){
        prop_index <- which(type == "prop")
        assay(object, assay.names[prop_index]) <- prop_mat
    }
    if("rc" %in% type){
        rc_index <- which(type == "rc")
        assay(object, assay.names[rc_index]) <- expected_rc
    }
    
    ## Return object
    object
}

#' Calculate Bayes Factors
#'
#' @param object PhIPData object
#' @param assay.postprob string indicating the assay where posterior
#' probabilities are stored. 
#' @param assay.name name indicating where the results should be stored in the
#' PhIPData object
#' @param prior.params prior parameters for the probability of enrichment 
#' (a_pi, b_pi)
#'
#' @return PhIPData object with the results stored in the location specified
#' by \code{assay.name}.
#'
#' @examples
#' sim_data <- readRDS(system.file("extdata", "sim_data.rds", package = "beer"))
#'
#' ## Calculate Bayes Factors
#' getBF(sim_data, "prob", "bayes_factor")
#'
#' @importFrom cli cli_alert_warning
#' @export
getBF <- function(object, assay.postprob = "prob", assay.name = "bayes_factors", 
    prior.params = list(a_pi = 2, b_pi = 300)) {
    
    ## Check inputs
    if (!assay.postprob %in% assayNames(object)) {
        stop(
            "Posterior probability assay not found."
        )
    }
    
    ## Check assay overwrite
    if (.checkOverwrite(object, assay.name)) {
        cli::cli_alert_warning("The specified assay exists and will be overwritten.")
    }
    
    ## Calculate Bayes Factors
    post_odds <- assay(object, assay.postprob)/(1 - assay(object, assay.postprob))
    prior_prob <- prior.params$a_pi/(prior.params$a_pi + prior.params$b_pi)
    prior_odds <- matrix(
        prior_prob/(1 - prior_prob),
        nrow = nrow(object), ncol = ncol(object)
    )

    assay(object, assay.name) <- post_odds/prior_odds
    
    object
}

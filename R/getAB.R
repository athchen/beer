#' @include edgeR.R utils.R

### edgeR estimates  ==============================================
#' @title Derive beta shape parameters using edgeR dispersion estimates
#'
#' Given a \code{\link[PhIPData]{PhIPData}} object, beads-only shape parameters
#' are estimated by first deriving the peptide-specific edgeR dispersion
#' estimate \eqn{\phi^{edgeR}}. \eqn{\phi^{edgeR}} corresponds to the squared
#' coefficient of variation for the proportion of reads pulled for a given
#' peptide. Using \eqn{\phi^{edgeR}} to derive an estimate of the variance for
#' the proportion of reads pulled by a single peptide, the mean and variance are
#' converted to shape parameters of a beta distribution.
#'
#' @param object a \code{\link[PhIPData]{PhIPData}} object.
#' @param threshold.cpm CPM threshold to be considered present in a sample.
#' @param threshold.prevalence proportion of beads-only samples that surpass
#' \code{threshold.cpm}.
#' @param lower minimum value of the beta shape parameters.
#' @param upper maximum value of the beta shape parameters.
#'
#' @return dataframe with rows corresponding to peptides and columns
#' corresponding to estimated shape parameters of the beta distribution.
#'
#' @seealso [.edgeRBeads()] for estimating \eqn{\phi^{edgeR}}
.getAB_edgeR <- function(object, threshold.cpm = 0, threshold.prevalence = 0,
    lower = 1, upper = Inf) {
    edgeR_beads <- .edgeRBeads(object, threshold.cpm, threshold.prevalence)

    ## Transform edgeR estimates to BB parameters
    mean_prop <- 2^edgeR_beads$AveLogCPM / 1e6
    var_prop <- edgeR_beads$tagwise.dispersion * mean_prop^2

    ## In these settings a << b, so we first adjust a, then calculate b.
    alpha_0 <- (1 - mean_prop) * mean_prop^2 / var_prop - mean_prop
    alpha_0 <- vapply(
        alpha_0, function(a_est) {
            if (a_est < lower) {
                lower
            } else if (a_est > upper) {
                upper
            } else {
                a_est
            }
        },
        numeric(1)
    )

    beta_0 <- alpha_0 * (1 / mean_prop - 1)

    data.frame(a_0 = alpha_0, b_0 = beta_0)
}

### Method of moments estimates  ==============================================
#' Helper function to derive MOM estimates of a, b from a vector of proportions
#'
#' @param prop vector of proportions.
#' @param offsets vector defining the offset to use when the mean and/or
#' variance are zero.
#' @param lower lowerbound for the shape parameters.
#' @param upper upper bound for the shape parameters.
#' @param ... parameters passed to \code{[base::mean]} and
#' \code{[stats::var]}.
#'
#' @importFrom stats var
#' @return a data frame with MOM estimates of a, b
.getAB_MOM_prop <- function(prop, offsets = c(mean = 1e-8, var = 1e-8),
    lower = 1, upper = Inf, ...) {
    mean_prop <- mean(prop, ...)
    var_prop <- var(prop, ...)

    ## If mean is negative, return an error
    ## if mean is 0, add small offset else use mean
    mean_prop <- if (mean_prop < 0) {
        stop("Mean cannot be negative.")
    } else if (mean_prop == 0) {
        mean_prop + offsets[["mean"]]
    } else {
        mean_prop
    }

    ## If variance is 0, add small offset
    var_prop <- if (var_prop < 0) {
        stop("Variance cannot be negative.")
    } else if (var_prop == 0) {
        var_prop + offsets[["var"]]
    } else {
        var_prop
    }

    ## estimate a and b
    ## In these settings a << b, so we first adjust a, then calculate b.
    a_est <- (1 - mean_prop) * mean_prop^2 / var_prop - mean_prop
    a_est <- if (a_est < lower) {
        lower
    } else if (a_est > upper) {
        upper
    } else {
        a_est
    }
    b_est <- a_est * (1 / mean_prop - 1)

    c(a = a_est, b = b_est)
}

#' Wrapper function to derive MOM estimates of a, b from beads-only samples
#'
#' @param object a \code{\link[PhIPData]{PhIPData}} object.
#' @param offsets vector defining the offset to use when the mean and/or
#' variance are zero.
#' @param lower lowerbound for the shape parameters.
#' @param upper upper bound for the shape parameters.
#' @param ... parameters passed to \code{[base::mean]} and
#' \code{[stats::var]}.
#'
#' @return a data frame with MOM estimates of a, b
.getAB_MOM <- function(object, offsets = c(mean = 1e-8, var = 1e-8),
    lower = 1, upper = Inf, ...) {
    prop_dat <- PhIPData::propReads(object)
    params <- apply(prop_dat, 1, .getAB_MOM_prop,
        offsets = offsets,
        lower = lower, upper = upper, ...
    )

    data.frame(a_0 = params[1, ], b_0 = params[2, ])
}

### MLE estimates  ==============================================
#' Helper function to derive MLE estimates of a, b from a vector of proportions
#'
#' @param prop vector of proportions.
#' @param prop.offset offset to use when the proportion of reads is 0.
#' @param optim.method optimization method passed to \code{[stats::optim]}.
#' @param lower lowerbound for the shape parameters.
#' @param upper upper bound for the shape parameters.
#'
#' @return a data frame of MLE estimates of a, b
#'
#' @seealso \code{[stats::optim]} for available optimization methods
#'
#' @import PhIPData
#' @importFrom stats dbeta optim
.getAB_MLE_prop <- function(prop, prop.offset = 1e-8, optim.method = "default",
    lower = 1, upper = Inf) {

    ## Add small offset when the proportion equals to 0
    prop <- prop + prop.offset * (prop == 0)

    ## Use MOM as initial values
    start <- .getAB_MOM_prop(prop)

    optim.method <- if (optim.method == "default") {
        optim.method <- "L-BFGS-B"
    } else {
        optim.method
    }

    nll <- function(prop, par) {
        log_lik <- -sum(dbeta(prop, par[["a"]], par[["b"]], log = TRUE))
        if (is.infinite(log_lik)) {
            sign(log_lik) * .Machine$double.xmax
        } else {
            log_lik
        }
    }

    opt <- stats::optim(
        par = start, fn = nll, prop = prop,
        method = optim.method,
        lower = lower, upper = upper
    )

    c(a_0 = opt$par[["a"]], b_0 = opt$par[["b"]])
}

#' Wrapper function to derive MLE estimates of a, b from beads-only samples
#'
#' @param object a \code{\link[PhIPData]{PhIPData}} object
#' @param prop.offset offset to use when the proportion of reads is 0.
#' @param optim.method optimization method passed to \code{[stats::optim]}.
#' @param lower lowerbound for the shape parameters.
#' @param upper upper bound for the shape parameters.
#'
#' @return a data frame of MLE estimates of a, b
#'
#' @seealso \code{[stats::optim]} for available optimization methods
.getAB_MLE <- function(object, prop.offset = 1e-8, optim.method = "default",
    lower = 1, upper = Inf) {
    prop_dat <- propReads(object)

    params <- apply(prop_dat, 1, .getAB_MLE_prop,
        prop.offset = prop.offset,
        optim.method = optim.method, lower = lower, upper = upper
    )

    data.frame(a_0 = params[1, ], b_0 = params[2, ])
}

### Estimating beads-only shape parameters  ==================================
#' Estimate beads-only shape parameters
#'
#' @description Beta shape parameters are estimated using the proportion of
#' reads-pulled per petide across the beads-only samples. Currently, only three
#' estimation methods are supported: edgeR, method of moments (MOM),
#' maximum likelihood (MLE). Note that edgeR can only be used on
#' \code{\link[PhIPData]{PhIPData}} objects while MOM and MLE methods can also
#' be applied to vectors of values between 0 and 1. Parameters that can be
#' passed to each method are listed in the details.
#'
#' @details \strong{edgeR} derived estimates rely on edgeR's peptide-specific
#' dispersion estimates, denoted \eqn{\phi^{edgeR}}. \eqn{\phi^{edgeR}}
#' corresponds to the squared coefficient of variation for the proportion of
#' reads pulled for a given peptide. Using \eqn{\phi^{edgeR}} to derive an
#' estimate of the variance for the proportion of reads pulled by a single
#' peptide, the mean and variance are transformed into shape parameters
#' satisfying the lower and upper bounds. When \code{method = "edgeR"}, the
#' following additional parameters can be specified.
#'
#' \itemize{
#'     \item \code{threshold.cpm}: CPM threshold to be considered present in
#'     a sample.
#'     \item \code{threshold.prevalence}: proportion of beads-only samples
#'     that surpass \code{threshold.cpm}.
#'     \item \code{lower}: minimum value of the beta shape parameters.
#'     \item \code{upper}: maximum value of the beta shape parameters.
#' }
#'
#' \strong{Method of Moments (MOM)} estimates are derived by transforming the
#' sample mean and variance to shape parameters of the beta distribution. For
#' \code{method = "mom"}, the following parameters can be adjusted:
#'
#' \itemize{
#'     \item \code{offsets}: vector defining the offset to use when the mean
#'     and/or variance are zero.
#'     \item \code{lower}: lowerbound for the shape parameters.
#'     \item \code{upper}: upper bound for the shape parameters.
#'     \item \code{...}: parameters passed to \code{[base::mean]} and
#'     \code{[stats::var]}.
#' }
#'
#' \strong{Maximum Likelihood (MLE)} estimates rely on
#' \code{[stats::optim]} to derive shape parameters that maximize the
#' likelihood of observed data. By default the L-BFGS-B optimization method is
#' used. Parameters for MLE estimates include:
#'
#' \itemize{
#'      \item \code{prop.offset}: offset to use when the proportion of reads
#'      is 0.
#'      \item \code{optim.method}: optimization method passed to
#'      \code{[stats::optim]}.
#'      \item \code{lower}: lowerbound for the shape parameters.
#'      \item \code{upper}: upper bound for the shape parameters.
#' }
#'
#' @param object a \code{\link[PhIPData]{PhIPData}} object or a vector
#' @param method one of \code{c("edgeR", "mle", "mom")} designating which method
#' to use to estimate beads-only prior parameters. MOM is the default method
#' used to estimate shape parameters.
#' @param ... parameters passed to specific estimating functions. See details
#' for more information
#'
#' @return a data frame of beta shape parameters where each row corresponds to
#' a peptide.
#'
#' @examples
#' ## PhIPData object
#' sim_data <- readRDS(system.file("extdata", "sim_data.rds", package = "beer"))
#'
#' getAB(sim_data, method = "edgeR")
#' getAB(sim_data, method = "mle")
#' getAB(sim_data, method = "mom")
#'
#' ## Vector of proportions
#' prop <- rbeta(100, 2, 8)
#' getAB(prop, method = "mle")
#' getAB(prop, method = "mom")
#' @importFrom methods is
#' @export
getAB <- function(object, method = "mom", ...) {

    ## Check that specified method is valid
    if (!method %in% c("edgeR", "mle", "mom")) {
        stop(
            "Invalid specified method for estimating a, b. ",
            "Valid methods are 'edgeR', 'mle', and 'mom'."
        )
    }

    if (is.vector(object)) {
        ## Check that valid inputs are between 0 and 1
        if (any(!is.numeric(object)) | max(object) > 1 | min(object) < 0) {
            stop(
                "Invalid inputs. Vectors can only contain numeric ",
                "values between 0 and 1"
            )
        }

        if (method == "mle") {
            .getAB_MLE_prop(object, ...)
        } else if (method == "mom") {
            .getAB_MOM_prop(object, ...)
        } else {
            stop("edgeR is not a valid method for vectors.")
        }
    } else if (is(object, "PhIPData")) {
        if (method == "mom") {
            .getAB_MOM(object, ...)
        } else if (method == "mle") {
            .getAB_MLE(object, ...)
        } else {
            .getAB_edgeR(object, ...)
        }
    } else {
        stop(
            "Invalid inputs. 'object' must be a vector of values ",
            "between 0 and 1 or a PhIPData object."
        )
    }
}

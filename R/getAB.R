#' @include edgeR.R utils.R

### edgeR estimates  ==============================================
#' @export
.getAB_edgeR <- function(object, threshold.cpm = 0, threshold.prevalence = 0,
                         lower = 1, upper = Inf){

    edgeR_beads <- .edgeRBeads(object, threshold.cpm, threshold.prevalence)

    ## Transform edgeR estimates to BB parameters
    mean_prop <- 2^edgeR_beads$AveLogCPM/1e6
    var_prop <- edgeR_beads$tagwise.dispersion*mean_prop^2

    ## In these settings a << b, so we first adjust a, then calculate b.
    alpha_0 <- (1-mean_prop)*mean_prop^2/var_prop - mean_prop
    alpha_0 <- vapply(alpha_0, function(a_est) {
        if(a_est < lower){ lower
            } else if (a_est > upper){ upper
            } else a_est},
        numeric(1))

    beta_0 <- alpha_0*(1/mean_prop - 1)

    data.frame(a_0 = alpha_0, b_0 = beta_0)
}

### Method of moments estimates  ==============================================
#' @export
.getAB_MOM_prop <- function(prop, offsets = c(mean = 1e-8, var = 1e-8),
                            lower = 1, upper = Inf, ...){

    mean_prop <- mean(prop, ...)
    var_prop <- var(prop, ...)

    ## If mean is negative, return an error
    ## if mean is 0, add small offset else use mean
    mean_prop <- if(mean_prop < 0){
        stop("Mean cannot be negative.")
    } else if (mean_prop == 0) {
        mean_prop + offsets[["mean"]]
    } else mean_prop

    ## If variance is 0, add small offset
    var_prop <- if(var_prop < 0){
        stop("Variance cannot be negative.")
    } else if(var_prop == 0) {
        var_prop + offsets[["var"]]
    } else var_prop

    ## estimate a and b
    ## In these settings a << b, so we first adjust a, then calculate b.
    a_est <- (1-mean_prop)*mean_prop^2/var_prop - mean_prop
    a_est <- if(a_est < lower){
        lower
    } else if (a_est > upper){
        upper
    } else a_est
    b_est <- a_est*(1/mean_prop - 1)

    c(a = a_est, b = b_est)
}

#' @export
.getAB_MOM <- function(object, offsets = c(mean = 1e-8, var = 1e-8),
                      lower = 1, upper = Inf, ...){

    prop_dat <- getPropReads(object)
    params <- apply(prop_dat, 1, .getAB_MOM_prop, offsets = offsets,
                    lower = lower, upper = upper, ...)

    data.frame(a_0 = params[1, ], b_0 = params[2, ])
}

### MLE estimates  ==============================================
#' @export
.getAB_MLE_prop <- function(prop, prop.offset = 1e-8, optim.method = "default",
                       lower = 1, upper = Inf){

    ## Add small offset when the proportion equals to 0
    prop <- prop + prop.offset*(prop == 0)

    ## Use MOM as initial values
    start <- .getAB_MOM_prop(prop)

    optim.method <- if(optim.method == "default"){
        optim.method <- "L-BFGS-B"
    } else optim.method

    nll <- function(prop, par){

        log_lik <- -sum(dbeta(prop, par[["a"]], par[["b"]], log = TRUE))
        if(is.infinite(log_lik)){
            sign(log_lik)*.Machine$double.xmax
        } else log_lik
    }

    opt <- stats::optim(par = start, fn = nll, prop = prop,
                        method = optim.method,
                        lower = lower, upper = upper)

    c(a_0 = opt$par[["a"]], b_0 = opt$par[["b"]])
}

#' @export
.getAB_MLE <- function(object, prop.offset = 1e-8, optim.method = "default",
                      lower = 1, upper = Inf){

    prop_dat <- getPropReads(object)

    params <- apply(prop_dat, 1, .getAB_MLE_prop, prop.offset = prop.offset,
                    optim.method = optim.method, lower = lower, upper = upper)

    data.frame(a_0 = params[1, ], b_0 = params[2, ])
}

### MLE estimates  ==============================================
#' @export
getAB <- function(object, method = "mom", ...){

    ## Check that specified method is valid
    if(!method %in% c("edgeR", "mle", "mom")){
        stop(paste0("Invalid specified method for estimating a, b. ",
                    "Valid methods are 'edgeR', 'mle', and 'mom'."))
    }

    if(is.vector(object)){
        if (method == "mle"){ .getAB_MLE_prop(object, ...)
        } else if (method == "mom") { .getAB_MOM_prop(object, ....)
        } else stop("edgeR is not a valid method for vectors.")
    } else if (is(object, "PhIPData")){
        if(method == "mom"){ .getAB_MOM(object, ...)
        } else if (method == "mle"){ .getAB_MLE(object, ...)
        } else .getAB_edgeR(object, ...)
    } else {
        stop("Invalid inputs. 'object' must be a vector or a PhIPData object.")
    }
}

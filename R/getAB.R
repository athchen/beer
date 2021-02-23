#' @include utils.R

### edgeR estimates  ==============================================
#' @import edgeR
getAB_edgeR <- function(object, threshold.cpm = 1, threshold.prevalence = 2){
    dge_obj <- as(object, "DGEList")

    # filter data to keep only peptides with cpm > 1
    # in at least 2 beads-only samples
    keep_ind <- rowSums(edgeR::cpm(dge_obj) > threshold.cpm) >=
        threshold.prevalence

    # calculate common dispersion in beads only data
    dge_obj <- dge_obj[keep_ind, , keep.lib.size = FALSE]
    dge_obj <- edgeR::calcNormFactors(dge_obj)
    dge_obj <- edgeR::estimateCommonDisp(dge_obj)
    dge_obj <- edgeR::estimateTagwiseDisp(dge_obj)

    # Transform edgeR estimates to BB parameters
    phat <- 2^dge_obj$AveLogCPM/1e6
    disp_BB <- ((1 + dge_obj$pseudo.lib.size*phat*dge_obj$tagwise.dispersion)/
                    (1-phat) - 1)/(dge_obj$pseudo.lib.size - 1)
    beta_0 <- (1 - phat) * (1 - disp_BB)/(disp_BB)
    alpha_0 <- phat*beta_0/(1 - phat)

    data.frame(a_0 = alpha_0, b_0 = beta_0)
}

### Method of moments estimates  ==============================================
getAB_MOM_ <- function(prop, offsets = c(mean = 1e-8, var = 1e-8),
                       lower = 0, upper = Inf, ...){

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

getAB_MOM <- function(object, offsets = c(mean = 1e-8, var = 1e-8),
                      lower = 0, upper = Inf, ...){

    prop_dat <- get_phat(object)
    mv_dat <- apply(prop_dat, 1, function(row){
        c(mean = mean(row), var = var(row))
    })

    params <- apply(prop_dat, 1, getAB_MOM_, offsets = offsets,
                    lower = lower, upper = upper, ...)

    data.frame(a_0 = params[1, ], b_0 = params[2, ])
}

### MLE estimates  ==============================================
getAB_MLE_ <- function(prop, prop.offset = 1e-8, optim.method = "default",
                       lower = 0, upper = Inf){

    ## Add small offset when the proportion equals to 0
    prop <- prop + prop.offset*(prop == 0)

    ## Use MOM as initial values
    start <- getAB_MOM_(prop)

    optim.method <- if(optim.method == "default"){
        optim.method <- "L-BFGS-B"
    } else optim.method

    nll <- function(prop, par){

        log_lik <- -sum(dbeta(prop, par[["a"]], par[["b"]], log = TRUE))
        if(is.infinite(log_lik)){
            sign(log_lik)*Machine$double.xmax
        } else log_lik
    }

    opt <- stats::optim(par = start, fn = nll, prop = prop,
                        method = optim.method,
                        lower = lower, upper = upper)

    c(a_0 = opt$par[["a"]], b_0 = opt$par[["b"]])
}

getAB_MLE <- function(object, prop.offset = 1e-8, optim.method = "default",
                      lower = 0, upper = Inf){
    prop_dat <- get_phat(object)

    params <- apply(prop_dat, 1, getAB_MLE_, prop.offset = prop.offset,
                    optim.method = optim.method, lower = lower, upper = upper)

    data.frame(a_0 = params[1, ], b_0 = params[2, ])
}

### MLE estimates  ==============================================
getAB <- function(object, method = "mom", ...){

    if(method == "edgeR"){ getAB_edgeR(object)
    } else if (method == "mle"){ getAB_MLE(object, ...)
    } else getAB_MOM(object, ...)

}

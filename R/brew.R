#' @include utils.R beadsRR.R getAB.R summarizeRun.R
NULL

#' Clean up inputs for prior estimation
#'
#' Tidy inputs related to `prior.parameters`. Supplies default values for
#' missing parameters and ensures that all required parameters are present.
#'
#' @param prior.params named list of prior parameters
#' @param object PhIPData object
#' @param beads.args parameters used to estimate a_0, b_0
#'
#' @return tidied list of prior parameters.
.tidyInputsPrior <- function(prior.params, object, beads.args){
    params <- c("method", "a_0", "b_0", "a_pi", "b_pi",
                "a_phi", "b_phi", "a_c", "b_c", "fc")

    ## Set missing parameters to defaults
    if(!"method" %in% names(prior.params)){
        if(all(c("a_0", "b_0") %in% names(prior.params))){
            prior.params$method <- "custom"
        } else prior.params$method <- "edgeR"
    }

    if(!all(c("a_0", "b_0") %in% names(prior.params))){
        beads_ab <- do.call(getAB, c(list(object = subsetBeads(object),
                                          method = prior.params$method),
                                     beads.args))
        prior.params$a_0 <- beads_ab[["a_0"]]
        prior.params$b_0 <- beads_ab[["b_0"]]
    }

    if(!"fc" %in% names(prior.params)) prior.params$fc <- 1

    ## Check that all necessary prior params are specified
    if(!all(params %in% names(prior.params))){
        missing_params <- params[!params %in% names(prior.params)]
        stop(paste0("The following elements are missing in prior.params: ",
                    paste0(missing_params, collapse = ", "), "."))
    }

    ## Check that method is valid
    if(!prior.params$method %in% c("edgeR", "mle", "mom", "custom")){
        stop(paste0("Invalid specified method. ",
                    "Valid methods are 'custom', 'edgeR', 'mle', and 'mom'."))
    }

    ## Return only needed parameters
    prior.params[params]
}

#' Clean up inputs for identifying super-enriched peptides
#'
#' Tidy inputs related to `se.params`. Supplies default values for
#' missing parameters and ensures that all required parameters are present.
#'
#' @param se.params named list of parameters for super-enriched estimation
#' @param beads.prior data.frame with beads-only parameters
#'
#' @return tidied list of parameters for identifying super-enriched peptides.
.tidyInputsSE <- function(se.params, beads.prior){

    ## Set method to default if it is missing
    if(!"method" %in% names(se.params)) se.params$method <- "mle"

    ## Check that method is valid
    if(!se.params$method %in% c("edgeR", "mle")){
        stop(paste0("Invalid specified method. ",
                    "Valid methods are 'edgeR' and 'mle'"))
    }

    ## Set missing parameters to defaults and return only needed parameters
    if(se.params$method == "edgeR"){
        if(!"threshold" %in% names(se.params)) se.params$threshold <- 15
        if(!"fc.name" %in% names(se.params)) se.params$fc.name <- "logfc"

        se.params[c("method", "threshold", "fc.name")]

    } else {
        if(!"threshold" %in% names(se.params)) se.params$threshold <- 15
        if(!"beads.prior" %in% names(se.params))
            se.params$beads.prior <- beads.prior

        se.params[c("method", "threshold", "beads.prior")]
    }
}

#' Clean inputs for JAGS parameters
#'
#' Tidy inputs related to `jags.params`. Supplies default values for
#' missing parameters and ensures that all required parameters are present.
#'
#' @param jags.params named list of JAGS parameters
#'
#' @return tidied list of JAGS parameters.
.tidyInputsJAGS <- function(jags.params){

    default <- list(n.chains = 1, n.adapt = 1e3,
                    n.iter = 1e4, thin = 1, na.rm = TRUE,
                    burn.in = 0, post.thin = 1,
                    seed = as.numeric(format(Sys.Date(), "%Y%m%d")))

    params <- c("n.chains", "n.adapt",
                "n.iter", "thin", "na.rm",
                "burn.in", "post.thin", "seed")

    ## Add missing inputs
    missing_params <- params[!params %in% names(jags.params)]
    jags.params[missing_params] <- default[missing_params]

    ## Check that all jags params are specified
    if(!all(params %in% names(jags.params))){
        missing_params <- params[!params %in% names(jags.params)]
        stop(paste0("The following elements are missing in jags.params: ",
                    paste0(missing_params, collapse = ", "), "."))
    }

    ## Return tidied jags needed parameters
    jags.params[params]
}

#' Clean-up specified assay names
#'
#' Tidy inputs related to `assay.names`. Supplies default values for
#' missing parameters and ensures that all required parameters are present.
#'
#' @param assay.names named list specifying where to store each assay.
#'
#' @return tidied list of assay.names
.tidyAssayNames <- function(assay.names){

    default <- c(phi = NA, phi_Z = "logfc", Z = "prob",
                 c = "sampleInfo", pi = "sampleInfo")

    assays <- c("phi", "phi_Z", "Z", "c", "pi")

    ## Set missing parameters to defaults
    missing_assays <- assays[!assays %in% names(assay.names)]
    assay.names[missing_assays] <- default[missing_assays]

    ## Check that c and pi assay options are valid
    valid <- c(is.na(assay.names["c"]) | assay.names["c"] == "sampleInfo",
               is.na(assay.names["pi"]) | assay.names["pi"] == "sampleInfo")
    if(sum(valid) != 2){
        stop(paste0("Invalid location specified. ",
                    paste0(c("c", "pi")[!valid], collapse = " and "),
                    " can only be stored in the sampleInfo of a ",
                    "PhIPData object (or not stored)."))
    }

    ## Check that phi, phi_Z, Z are unique assays
    if (length(unique((assay.names[c("phi", "phi_Z", "Z")]))) != 3) {
        stop("phi, phi_Z, and Z must have unique assay names.")
    }

    ## Return only needed parameters
    assay.names[assays]
}

#' Run BEER for one sample
#'
#' This function is not really for external use. It's exported for
#' parallelization purposes. For more detailed descriptions see
#' \code{\link{brew}}.
#'
#' @param object PhIPData object
#' @param sample sample name
#' @param prior.params vector of prior parameters
#' @param n.chains number of chains to run
#' @param n.adapt number of iterations to use as burn-in.
#' @param quiet logical value specifying whether warnings/informative messages
#' should be printed
#' @param n.iter number of iterations for the MCMC chain to run (after n.adapt)
#' @param thin thinning parameter
#' @param na.rm what to do with NA values (for JAGS)
#' @param ... extra params for JAGS
#' @param seed number/string for reproducibility purposes.
#'
#' @return nothing, saves the the results to an RDS in either a temp directory
#' or the specified directory.
#'
#' @examples
#' \dontrun{
#' sim_data <- readRDS(system.file("extdata", "sim_data.rds", package = "beer"))
#'
#' beads_prior <- getAB(subsetBeads(sim_data), "edgeR")
#' .brew_one(sim_data, "9", list(a_0 = beads_prior[["a_0"]],
#'           b_0 = beads_prior[["b_0"]],
#'           a_pi = 2, b_pi = 300,
#'           a_phi = 1.25, b_phi = 0.1,
#'           a_c = 80, b_c = 20,
#'           fc = 1))
#' }
#'
#' @export
.brew_one <- function(object, sample, prior.params,
                     n.chains = 1, n.adapt = 1e3,
                     n.iter = 1e4, thin = 1, na.rm = TRUE, ...,
                     seed = as.numeric(format(Sys.Date(), "%Y%m%d"))){
    ## Define data
    data_list <- list(N = 1,
                      P = nrow(object),
                      B = 0,
                      n = librarySize(object[, sample], withDimnames = FALSE),
                      Y = unname(counts(object)[, sample, drop = FALSE]))
    data_list <- c(data_list, prior.params)

    ## Define initial values
    inits_list <- guessInits(object[, sample],
                             list(a_0 = prior.params[["a_0"]],
                                  b_0 = prior.params[["b_0"]]))
    inits_list$`.RNG.name` <- "base::Wichmann-Hill"
    inits_list$`.RNG.seed` <- seed

    model_file <- system.file("extdata/phipseq_model.bugs", package = "beer")

    #Compile and run model
    capture.output({
        jags_model <- rjags::jags.model(file = model_file,
                                        data = data_list,
                                        inits = inits_list,
                                        n.chains = n.chains,
                                        n.adapt = n.adapt)
        mcmc <- coda.samples(jags_model,
                             variable.names = c("c", "pi", "Z", "phi"),
                             n.iter = n.iter, thin = thin,
                             na.rm = na.rm, ...)
    })

    mcmc
}

#' Run BEER for all samples
#'
#' Encapsulated function to run each sample against all beads-only samples.
#' The code is wrapped in this smaller function to (1) modularize the code and
#' (2) make sure the cli output colors don't change.
#'
#' @param object PhIPData object
#' @param sample_id vector of sample IDs to iterate over
#' @param beads_id vector of IDs corresponding to beads-only samples
#' @param se.matrix matrix indicating which peptides are clearly enriched
#' @param prior.params list of prior parameters
#' @param beads.prior data frame of beads-only prior parameters
#' @param jags.params list of JAGS parameters
#' @param tmp.dir directory to store JAGS samples
#'
#' @return vector of process id's for internal checking of whether functions
#' were parallized correctly.
.brew_samples <- function(object, sample_id, beads_id, se.matrix,
                          prior.params, beads.prior, beads.args, jags.params,
                          tmp.dir){
    cli::cli_text("Sample runs")
    progressr::handlers("txtprogressbar")
    p <- progressr::progressor(along = sample_id)

    jags_out <- future.apply::future_lapply(sample_id, function(sample){
        sample_counter <- paste0(which(sample_id == sample), " of ",
                                 length(sample_id))
        p(sample_counter, class = "sticky", amount = 1)

        ## Subset super-enriched and only single sample
        one_sample <- object[!se.matrix[, sample], c(beads_id, sample)]

        ## Calculate new beads-only priors
        new_beads <- if(prior.params$method == "custom"){
            lapply(beads.prior, function(x) x[!se.matrix[, sample]])
        } else {
            do.call(getAB, c(list(object = subsetBeads(one_sample),
                                  method = prior.params$method),
                             beads.args))
        }

        new_prior <- c(list(a_0 = new_beads[["a_0"]],
                            b_0 = new_beads[["b_0"]]),
                       prior.params[c("a_pi", "b_pi", "a_phi", "b_phi",
                                      "a_c", "b_c", "fc")])

        jags_run <- do.call(.brew_one, c(list(object = one_sample,
                                              sample = sample,
                                              prior.params = new_prior,
                                              quiet = TRUE),
                                         jags.params))

        saveRDS(jags_run, paste0(tmp.dir,"/", sample, ".rds"))

        Sys.getpid()
    }, future.conditions = "message")

    jags_out
}

#' Bayesian Enrichment Estimation in R (BEER)
#'
#' Run BEER to estimate posterior probabilities of enrichment, sample-specific
#' attenuation constants, relative fold-changes in comparison to beads-only
#' samples, and proportion of peptides enriched per sample as described in
#' Chen et. al. See \emph{Details} for more information on input parameters.
#'
#' @details \strong{\code{prior.params}}. List of prior parameters. Parameters
#' include,
#' \itemize{
#'     \item \code{method}: method used to estimate beads-only prior parameters
#'     a_0, b_0. Valid methods include 'custom' or any of the methods specified
#'     in \code{\link{getAB}}. If \code{method = 'custom'} is specified,
#'     \code{a_0} and \code{b_0} must be included in the list of prior
#'     parameters. \code{'edgeR'} is used as the default method for estimating
#'     a_0, b_0.
#'     \item \code{a_pi}, \code{b_pi}: prior shape parameters for the proportion
#'     of peptides enriched in a sample. Defaults to 2 and 300, respectively.
#'     \item \code{a_phi}, \code{b_phi}: prior shape parameters of the gamma
#'     distribution that describe the valid range of enriched-fold changes. The
#'     shift is specified by \code{fc}. The default values of \code{a_phi} and
#'     \code{b_phi} are 1.25 and 0.1, respectively.
#'     \item \code{a_c}, \code{b_c}: prior shape parameters for the attenuation
#'     constant. Default values for \code{a_c} and \code{b_c} are 80 and 20.
#'     \item{\code{fc}}: minimum fold change for an enriched-peptide. \code{fc}
#'     describes the shift in the gamma distribution.
#' }
#'
#' \strong{\code{beads.args}}. Named list of parameters supplied to
#' \code{\link[beer]{getAB}}. The estimation method used is specified in
#' \code{prior.params}, but other valid parameters include lower and upper
#' bounds for elicited parameters. As JAGS recommends that \eqn{a, b > 1} for
#' the beta distribution, \code{beads.args} defaults to \code{list(lower = 1)}.
#'
#' \strong{\code{se.params}}. Named list of parameters supplied to
#' \code{\link[beer]{guessEnriched}}. By default \code{list(method = 'mle')} is
#' used to identify clearly enriched peptides.
#'
#' \strong{\code{jags.params}}. Named list of parameters for MCMC sampling. By
#' default, BEER only runs one chain with 1,000 adaptation iteration and 10,000
#' sampling iterations. If unspecified, BEER uses the current date as the seed.
#'
#' \strong{\code{sample.dir}}. Path specifying where to store the intermediate
#' results. If \code{NULL}, then results are stored in the default temporary
#' directory. Otherwise, the MCMC samples for running BEER on each sample is
#' stored as a single \code{RDS} file in the specified directory.
#'
#' \strong{\code{assay.names}}. Named list specifying where to store the point
#' estimates. If \code{NULL}, estimates are not added to the PhIPData object.
#' Valid exported estimates include,
#' \itemize{
#'     \item \code{phi}: fold-change estimate after marginalizing over the
#'     posterior probability of enrichment. By default point estimates are not
#'     exported.
#'     \item \code{phi_Z}: fold-change estimate presuming the peptide is
#'     enriched. By default \code{phi_Z} estimates are stored in \code{'logfc'}
#'     assay.
#'     \item \code{Z}: posterior probability of enrichment. Estimates are stored
#'     in the \code{'prob'} assay by default.
#'     \item \code{c}: attenuation constant estimates. Stored in
#'     \code{'sampleInfo'} by default.
#'     \item \code{pi}: point estimates for the proportion of peptides enriched
#'     in a sample. Stored in \code{'sampleInfo'} by default.
#' }
#'
#' @param object PhIPData object
#' @param prior.params named list of prior parameters
#' @param beads.args named list of parameters supplied to estimating beads-only
#' prior parameters (a_0, b_0)
#' @param se.params named list of parameters specific to identifying clearly
#' enriched peptides
#' @param jags.params named list of parameters for running MCMC using JAGS
#' @param sample.dir path to temporarily store RDS files for each sample run
#' @param assay.names named vector indicating where MCMC results should be
#' stored in the PhIPData object
#' @param beadsRR logical value specifying whether each beads-only sample
#' should be compared to all other beads-only samples.
#' @param parallel character indicating which parallelization strategy to use.
#' Alternatively, a named list of parameters available in
#' \code{\link{future::plan}{future::plan()}}.
#'
#' @return A PhIPData object with BEER results stored in the locations specified
#' by \code{assay.names}.
#'
#' @seealso \code{\link[future]{plan}} for parallelization options,
#' \code{\link{beadsRR}} for running each beads-only sample against all
#' remaining samples, \code{\link{getAB}} for more information about valid parameters for estimating beads-only prior parameters,
#' \code{\link{guessEnriched}} for more information about how clearly
#' enriched peptides are identified, and \code{\link[rjags]{jags.model}} for
#' MCMC sampling parameters.
#'
#' @examples
#' sim_data <- readRDS(system.file("extdata", "sim_data.rds", package = "beer"))
#'
#' ## Sequential evaluation
#' brew(sim_data)
#'
#' ## Multisession evaluation
#' \dontrun{
#' beer(sim_data, parallel = "multisession")
#' }
#'
#' @importFrom future.apply future_lapply
#' @importFrom future plan
#' @importFrom cli cli_alert_warning
#'
#' @export
brew <- function(object,
                 prior.params = list(method = "edgeR",
                                     a_pi = 2, b_pi = 300,
                                     a_phi = 1.25, b_phi = 0.1,
                                     a_c = 80, b_c = 20,
                                     fc = 1),
                 beads.args = list(lower = 1),
                 se.params = list(method = "mle"),
                 jags.params = list(n.chains = 1, n.adapt = 1e3,
                                    n.iter = 1e4, thin = 1, na.rm = TRUE,
                                    burn.in = 0, post.thin = 1,
                                    seed = as.numeric(format(Sys.Date(),
                                                             "%Y%m%d"))),
                 sample.dir = NULL,
                 assay.names = c(phi = NULL, phi_Z = "logfc", Z = "prob",
                                 c = "sampleInfo", pi = "sampleInfo"),
                 beadsRR = FALSE,
                 parallel = "sequential"){

    .checkCounts(object)

    ## Check and tidy inputs
    prior.params <- .tidyInputsPrior(prior.params, object, beads.args)
    beads.prior <- prior.params[c("a_0", "b_0")]
    se.params <- if(!is.null(se.params)) { .tidyInputsSE(se.params, beads.prior)
    } else { NULL }
    jags.params <- .tidyInputsJAGS(jags.params)
    assay.names <- .tidyAssayNames(assay.names)
    parallel <- .tidyParallel(parallel)

    ## Define plan and return current plan to oplan
    oplan <- do.call(plan, parallel)
    ## Rest to original plan on exit
    on.exit(plan(oplan))

    ## Get sample names
    beads_id <- colnames(object[, object$group == getBeadsName()])
    sample_id <- colnames(object[, object$group != getBeadsName()])

    ## Guess which peptides are super-enriched
    se.matrix <- if(is.null(se.params)){
        matrix(FALSE, nrow(object), ncol(object),
               dimnames = dimnames(object))
    } else {
        do.call(guessEnriched, c(list(object = object), se.params))
    }

    ## Create temporary directory for output of JAGS models
    tmp.dir <- if(is.null(sample.dir)) {
        paste0(tempdir(), "/beer_run",
               as.numeric(format(Sys.Date(), "%Y%m%d")),
               as.numeric(Sys.time()))
    } else normalizePath(sample.dir, mustWork = FALSE)

    if(dir.exists(tmp.dir)){
        delete <- menu(c("Yes", "No"),
                       title = paste0("Specified directory for samples exists. ",
                                      "Delete the existing directory?"))
        if(delete == 1) unlink(tmp.dir, recursive = TRUE)
    }
    dir.create(tmp.dir, recursive = TRUE)

    ## Check whether assays will be overwritten
    beads_over <- .checkOverwrite(object[, object$group == getBeadsName()],
                                  assay.names)
    sample_over <- .checkOverwrite(object[ ,object$group != getBeadsName()],
                                   assay.names)
    msg <- if(beadsRR & any(beads_over | sample_over, na.rm = TRUE)){
        paste0("Values in the following assays will be overwritten: ",
               paste0(unique(assay.names[beads_over | sample_over]),
                      collapse = ", "))
    } else if (!beadsRR & any(sample_over, na.rm = TRUE)) {
        paste0("Values in the following assays will be overwritten: ",
               paste0(unique(assay.names[beads_over & sample_over]),
                      collapse = ", "))
    } else character(0)
    if(length(msg) > 0) cli::cli_alert_warning(msg)

    ## Run model one sample at a time
    cli::cli_h1("Running JAGS")
    if(beadsRR){
        cli::cli_text("Beads-only round robin")
        beadsRR(subsetBeads(object), method = "beer",
                prior.params, beads.args, se.matrix,
                jags.params, tmp.dir, parallel, parallel.params)
    }

    pids <- .brew_samples(object, sample_id, beads_id, se.matrix,
                          prior.params, beads.prior, beads.args, jags.params,
                          tmp.dir)

    ## Summarize output
    cli::cli_h1("Summarizing results")
    run_out <- summarizeRun(object, tmp.dir, se.matrix,
                            burn.in = jags.params$burn.in,
                            post.thin = jags.params$post.thin,
                            assay.names, quiet = FALSE)

    ## Clean-up if necessary
    if(is.null(sample.dir)){
        unlink(sample.dir, recursive = TRUE)
    }

    run_out
}

#' @include utils.R beadsRR.R getAB.R summarizeRun.R

.tidyInputsPrior <- function(prior.params, object, beads.args){
    params <- c("method", "a_0", "b_0", "a_pi", "b_pi",
                "a_phi", "b_phi", "a_c", "b_c", "fc")

    ## Set missing parameters to defaults
    if(!"method" %in% names(prior.params)){
        if(all(c("a_0", "b_0") %in% names(prior.params))){
            prior.params$method <- "custom"
        } else prior.params$method <- "mom"
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
    if(!all(names(prior.params) %in% params)){
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

.tidyInputsJAGS <- function(jags.params){

    default <- list(n.chains = 1, n.adapt = 1e3, quiet = FALSE,
                    n.iter = 1e4, thin = 1, na.rm = TRUE,
                    burn.in = 0, post.thin = 1,
                    seed = as.numeric(format(Sys.Date(), "%Y%m%d")))

    params <- c("n.chains", "n.adapt", "quiet",
                "n.iter", "thin", "na.rm",
                "burn.in", "post.thin", "seed")

    ## Add missing inputs
    missing_params <- params[!params %in% names(jags.params)]
    jags.params[missing_params] <- default[missing_params]

    ## Check that all jags params are specified
    if(!all(names(jags.params) %in% params)){
        missing_params <- params[!params %in% names(prior.params)]
        stop(paste0("The following elements are missing in jags.params: ",
                    paste0(missing_params, collapse = ", "), "."))
    }

    ## Return tidied jags needed parameters
    jags.params[params]
}

.tidyAssayNames <- function(assay.names){

    default <- c(phi = NA, phi_Z = "logfc", Z = "prob",
                 c = "metadata", pi = "metadata")

    assays <- c("phi", "phi_Z", "Z", "c", "pi")

    ## Set missing parameters to defaults
    missing_assays <- assays[!assays %in% names(assay.names)]
    assay.names[missing_assays] <- default[missing_assays]

    ## Check that c and pi assay options are valid
    valid <- c(is.na(assay.names["c"]) | assay.names["c"] == "metadata",
               is.na(assay.names["pi"]) | assay.names["pi"] == "metadata")
    if(sum(valid) != 2){
        stop(paste0("Invalid location specified. ",
                    paste0(c("c", "pi")[!valid], collapse = " and "),
                    " can only be stored in the metadata of a ",
                    "PhIPData object (or not stored)."))
    }

    ## Return only needed parameters
    assay.names[assays]
}

brew_one <- function(object, sample, prior.params,
                     n.chains = 1, n.adapt = 1e3, quiet = FALSE,
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

    # Compile and run model
    jags_model <- rjags::jags.model(file = model_file, data = data_list,
                                    inits = inits_list, n.chains = n.chains,
                                    n.adapt = n.adapt)

    rjags::coda.samples(jags_model,
                        variable.names = c("c", "pi", "Z", "phi"),
                        n.iter = n.iter,
                        thin = thin,
                        na.rm = na.rm, ...)
}

#' @export
brew <- function(object,
                 prior.params = list(method = "mom",
                                     a_pi = 2, b_pi = 300,
                                     a_phi = 1.25, b_phi = 0.1,
                                     a_c = 80, b_c = 20,
                                     fc = 1),
                 beads.args = list(lower = 1),
                 se.params = list(method = "mle"),
                 jags.params = list(n.chains = 1, n.adapt = 1e3, quiet = FALSE,
                                    n.iter = 1e4, thin = 1, na.rm = TRUE,
                                    burn.in = 0, post.thin = 1,
                                    seed = as.numeric(format(Sys.Date(), "%Y%m%d"))),
                 sample.dir = NULL,
                 assay.names = c(phi = NULL, phi_Z = "logfc", Z = "prob",
                                 c = "metadata", pi = "metadata"),
                 beadsRR = FALSE,
                 parallel = TRUE,
                 parallel.params = list()){

    .checkCounts(object)

    ## Check and tidy inputs
    prior.params <- .tidyInputsPrior(prior.params, object, beads.args)
    beads.prior <- prior.params[c("a_0", "b_0")]
    se.params <- if(!is.null(se.params)) { .tidyInputsSE(se.params, beads.prior)
    } else { NULL }
    jags.params <- .tidyInputsJAGS(jags.params)
    assay.names <- .tidyAssayNames(assay.names)

    ## Get sample names
    beads_id <- colnames(object[, object$group == getBeadsName()])
    sample_id <- colnames(object[, object$group != getBeadsName()])

    ## Guess which peptides are super-enriched
    se_peps <- if(is.null(se.params)){
        matrix(FALSE, nrow(object), ncol(object),
               dimnames = dimnames(object))
    } else {
        do.call(guessEnriched, c(list(object = object), se.params))
    }

    ## Create temporary directory for output of JAGS models
    tmp_dir <- if(is.null(sample.dir)) {
        paste0(tempdir(), "/beer_run", as.numeric(format(Sys.Date(), "%Y%m%d")))
    } else normalizePath(sample.dir)

    if(dir.exists(tmp_dir)){
        delete <- menu(c("Yes", "No"),
                       title = paste0("Specified directory for samples exists. ",
                                      "Delete the existing directory?"))
        if(delete == 1) unlink(tmp_dir, recursive = TRUE)
    }
    dir.create(tmp_dir, recursive = TRUE)

    ## Run model one sample at a time
    if(!jags.params$quiet) cli::cli_h1("Running JAGS")

    ## Run
    if(beadsRR){
       if(!jags.params$quiet) cli::cli_h2("Beads-only round robin")
        beadsRR(subsetBeads(object), method = "beer",
                prior.params, beads.args, se_peps,
                jags.params, tmp_dir, parallel, parallel.params)
    }

    if(!jags.params$quiet) cli::cli_h2("Sample runs")
    for(sample in sample_id){
        if(!jags.params$quiet) {
            print(paste0(which(sample_id == sample), " of ", length(sample_id)))
        }

        ## Subset super-enriched and only single sample
        one_sample <- object[!se_peps[, sample], c(beads_id, sample)]

        ## Calculate new beads-only priors
        new_beads <- if(prior.params$method == "custom"){
            lapply(beads.prior, function(x) x[!se_peps[, sample]])
        } else {
            do.call(getAB, c(list(object = subsetBeads(one_sample),
                                  method = prior.params$method),
                             beads.args))
        }

        new_prior <- c(list(a_0 = new_beads[["a_0"]],
                            b_0 = new_beads[["b_0"]]),
                       prior.params[c("a_pi", "b_pi", "a_phi", "b_phi",
                                      "a_c", "b_c", "fc")])
        jags_run <- do.call(brew_one, c(list(object = one_sample,
                                             sample = sample,
                                             prior.params = new_prior),
                                        jags.params))

        saveRDS(jags_run, paste0(tmp_dir,"/", sample, ".rds"))
    }

    ## Summarize output
    if(!jags.params$quiet) cli::cli_h1("Summarizing results")
    run_out <- summarizeRun(object, tmp_dir, se_peps,
                            burn.in = jags.params$burn.in,
                            post.thin = jags.params$post.thin,
                            assay.names, quiet = FALSE)

    ## Clean-up if necessary
    if(is.null(sample.dir)){
        unlink(sample.dir, recursive = TRUE)
    }

    run_out
}

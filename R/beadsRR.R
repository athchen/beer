#' Function to run the beads-only round robin using BEER
#'
#' Each sample is run in comparison to all other beads-only samples to
#' approximate the false positive rate of detecting enrichments.
#'
#' @param object PhIPData object
#' @param prior.params named list of prior parameters
#' @param beads.args named list of parameters supplied to estimating beads-only
#' prior parameters (a_0, b_0)
#' @param jags.params named list of parameters for running MCMC using JAGS
#' @param sample.dir path to temporarily store RDS files for each sample run,
#' if \code{NULL} then \code{[base::tempdir]} is used to temporarily store
#' MCMC output and cleaned afterwards.
#' @param assay.names named vector indicating where MCMC results should be
#' stored in the PhIPData object
#' @param summarize logical indicating whether to return a PhIPData object.
#' @param bp.param \code{[BiocParallel::BiocParallelParam]} passed to
#' BiocParallel functions.
#'
#' @return vector of process IDs or a PhIPData object
#'
#' @import PhIPData
#' @importFrom BiocParallel bplapply
#' @importFrom progressr handlers progressor
.beadsRRBeer <- function(object,
    prior.params = list(
        method = "edgeR",
        a_pi = 2, b_pi = 300,
        a_phi = 1.25, b_phi = 0.1,
        a_c = 80, b_c = 20,
        fc = 1
    ),
    beads.args = list(lower = 1),
    jags.params = list(
        n.chains = 1, n.adapt = 1e3,
        n.iter = 1e4, thin = 1, na.rm = TRUE,
        burn.in = 0, post.thin = 1,
        seed = as.numeric(format(
            Sys.Date(),
            "%Y%m%d"
        ))
    ),
    sample.dir = NULL,
    assay.names = c(
        phi = NULL, phi_Z = "logfc", Z = "prob",
        c = "sampleInfo", pi = "sampleInfo"
    ),
    summarize = TRUE,
    bp.param = bpparam()) {
    .checkCounts(object)

    ## Check and tidy inputs
    prior.params <- .tidyInputsPrior(prior.params, object, beads.args)
    beads.prior <- prior.params[c("a_0", "b_0")]
    jags.params <- .tidyInputsJAGS(jags.params)
    assay.names <- .tidyAssayNames(assay.names)

    ## Create temporary directory for output of JAGS models
    tmp.dir <- if (is.null(sample.dir)) {
        file.path(
            tempdir(),
            paste0("beer_run", as.numeric(format(Sys.Date(), "%Y%m%d")))
        )
    } else {
        normalizePath(sample.dir, mustWork = FALSE)
    }

    ## Check if any sample files exist
    beads_id <- colnames(object[, object$group == getBeadsName()])
    if (any(file.exists(file.path(tmp.dir, paste0(beads_id, ".rds"))))) {
        cli::cli_alert_warning(
            paste0(
                "Sample files already exist in the specified directory. ",
                "Some files may be overwritten."
            )
        )
    } else if (!dir.exists(tmp.dir)) {
        dir.create(tmp.dir, recursive = TRUE)
    }

    progressr::handlers("txtprogressbar")
    p <- progressr::progressor(along = colnames(object))

    jags_out <- bplapply(beads_id, function(sample) {
        sample_counter <- paste0(
            which(colnames(object) == sample), " of ",
            length(beads_id)
        )
        p(sample_counter, class = "sticky", amount = 1)

        # recode existing object
        one_beads <- object
        one_beads$group[colnames(one_beads) == sample] <- "sample"

        ## Calculate new beads-only priors
        new_beads <- if (prior.params$method == "custom") {
            beads.prior
        } else {
            do.call(getAB, c(
                list(
                    object = subsetBeads(one_beads),
                    method = prior.params$method
                ),
                beads.args
            ))
        }

        new_prior <- c(
            list(
                a_0 = new_beads[["a_0"]],
                b_0 = new_beads[["b_0"]]
            ),
            prior.params[c(
                "a_pi", "b_pi", "a_phi", "b_phi",
                "a_c", "b_c", "fc"
            )]
        )

        jags_run <- do.call(brewOne, c(
            list(
                object = one_beads,
                sample = sample,
                prior.params = new_prior
            ),
            jags.params
        ))

        saveRDS(jags_run, file.path(tmp.dir, paste0(sample, ".rds")))

        Sys.getpid()
    }, BPPARAM = bp.param)

    if (summarize) {
        out <- summarizeRun(object, file.path(tmp.dir, paste0(beads_id, ".rds")),
            matrix(FALSE, nrow(object), ncol(object),
                dimnames = dimnames(object)
            ),
            burn.in = jags.params$burn.in,
            post.thin = jags.params$post.thin,
            assay.names,
            bp.param
        )

        ## Clean-up after summarization
        if (is.null(sample.dir)) {
            unlink(tmp.dir, recursive = TRUE)
        }
    } else {
        out <- unlist(jags_out)
    }

    out
}

#' Function to run the beads-only round robin using edgeR
#'
#' Each sample is run in comparison to all other beads-only samples to
#' approximate the false positive rate of detecting enrichments.
#'
#' @param object A PhIPData object of only beads-only samples.
#' @param threshold.cpm CPM threshold to be considered present in a sample
#' @param threshold.prevalence proportion of beads-only samples that surpass
#' \code{threshold.cpm}.
#' @param assay.names named vector specifying the assay names for the
#' log2(fold-change) and exact test p-values. If the vector is not names,
#' the first and second entries are used as defaults
#' @param bp.param \code{[BiocParallel::BiocParallelParam]} passed to
#' BiocParallel functions.
#'
#' @return vector of process IDs
#'
#' @import PhIPData SummarizedExperiment
#' @importFrom BiocParallel bplapply
.beadsRREdgeR <- function(object, threshold.cpm = 0, threshold.prevalence = 0,
    assay.names = c(logfc = "logfc", prob = "prob"),
    bp.param = BiocParallel::bpparam()) {

    ## Set-up output matrices ----------
    ## Make empty matrix for the cases where fc and prob do not exist
    empty_mat <- matrix(nrow = nrow(object), ncol = ncol(object))
    colnames(empty_mat) <- colnames(object)
    rownames(empty_mat) <- rownames(object)

    edgeR_fc <- if (assay.names[["logfc"]] %in% assayNames(object)) {
        assay(object, assay.names[["logfc"]])
    } else {
        empty_mat
    }

    edgeR_pval <- if (assay.names[["prob"]] %in% assayNames(object)) {
        assay(object, assay.names[["prob"]])
    } else {
        empty_mat
    }

    beads_names <- colnames(subsetBeads(object))

    output <- bplapply(beads_names, function(sample) {
        ## Recode existing object
        one_beads <- object
        one_beads$group[colnames(one_beads) == sample] <- "sample"

        ## Derive dispersion estimates from beads-only samples
        edgeR_beads <- .edgeRBeads(object, threshold.cpm, threshold.prevalence)
        common_disp <- edgeR_beads$common.dispersion
        tagwise_disp <- edgeR_beads$tagwise.dispersion
        trended_disp <- edgeR_beads$trended.dispersion

        ## Run edgeR for the one beads-only sample ------------
        edgeROne(
            one_beads, sample, beads_names[beads_names != sample],
            common_disp, tagwise_disp, trended_disp
        )
    })

    ## Unnest output
    for (result in output) {
        edgeR_fc[, result$sample] <- result$logfc
        edgeR_pval[, result$sample] <- result$log10pval
    }

    ## append to object
    assay(object, assay.names[["logfc"]]) <- edgeR_fc
    assay(object, assay.names[["prob"]]) <- edgeR_pval

    object
}

#' Beads-only round robin
#'
#' To approximate the false positive rate of each approach, each beads-only
#' sample is run individually against all other samples. For BEER, this means
#' that the sample to be compared is encoded as an actual sample, and prior
#' parameters for beads-only samples are re-estimated. Thus, the beads-only
#' round robin also serves to assess how similar the beads-only samples are to
#' one another.
#'
#' @details If \strong{\code{method == 'beer'}}, then valid parameters include
#' \code{prior.params}, \code{beads.args}, \code{jags.params},
#' \code{sample.dir}, \code{assay.names}, and \code{summarize}. A description
#' of the first four parameters can be found in \code{\link{brew}}.
#' \code{summarize} is a logical value indicating whether a PhIPData object
#' with the summarized results should be returned. When running \code{beadsRR},
#' \code{summarize} typically does not need to be changed.
#'
#' When \strong{\code{method == 'edgeR'}}, \code{threshold.cpm},
#' \code{threshold.prevalence}, and \code{assay.names} are valid additional
#' parameters that can be supplied to \code{beadsRR}. See \code{\link{edgeR}}
#' for additional details on each of these parameters.
#'
#' @param object PhIPData object
#' @param method one of \code{'beer'} or \code{'edgeR'} specifying which method
#' to use.
#' @param bp.param \code{[BiocParallel::BiocParallelParam]} passed to
#' BiocParallel functions.
#' @param ... parameters passed to the method specific functions. See the
#' \emph{Details} section below for additional information.
#'
#' @return a PhIPData object
#'
#' @seealso \code{\link{brew}} for BEER parameters, \code{\link{edgeR}} for
#' edgeR parameters, and \code{[BiocParallel::BiocParallelParam]} for
#' parallelization.
#'
#' @examples
#' sim_data <- readRDS(system.file("extdata", "sim_data.rds", package = "beer"))
#'
#' beadsRR(sim_data, method = "beer")
#' beadsRR(sim_data, method = "edgeR")
#' @export
beadsRR <- function(object, method, bp.param = BiocParallel::bpparam(), ...) {
    if (!method %in% c("edgeR", "beer")) {
        stop("Invalid specified method for beads-only round robin.")
    } else if (method == "edgeR") {
        .beadsRREdgeR(object, bp.param = bp.param, ...)
    } else {
        .beadsRRBeer(object, bp.param = bp.param, ...)
    }
}

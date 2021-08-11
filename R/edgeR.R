#' Estimate edgeR dispersion parameters from the beads-only data
#'
#' Wrapper function to estimate edgeR dispersion parameters from beads-only
#' samples. Peptides can be pre-filtered based on a minimum read count per
#' million (cpm) and the proportion of beads-only samples that surpass the cpm
#' threshold.
#'
#' @param object \code{\link[PhIPData]{PhIPData}} object (can have actual serum
#' samples)
#' @param threshold.cpm CPM threshold to be considered present in a sample
#' @param threshold.prevalence proportion of beads-only samples that surpass
#' \code{threshold.cpm}.
#'
#' @return a DGEList object with common, trended, and tagwise dispersion
#' estimates
.edgeRBeads <- function(object, threshold.cpm = 0, threshold.prevalence = 0){
    edgeR_beads <- as(PhIPData::subsetBeads(object), "DGEList")

    ## For dispersion estimates, keep only peptides with cpm above the
    ## specified threshold in the given number of samples.
    keep_ind <- rowSums(edgeR::cpm(edgeR_beads) >= threshold.cpm) >=
        threshold.prevalence

    ## Estimate common, trended, and tagwise dispersion in the beads-only data
    edgeR_beads <- edgeR_beads[keep_ind, , keep.lib.size = FALSE]
    edgeR_beads <- edgeR::calcNormFactors(edgeR_beads)
    edgeR_beads <- suppressMessages(edgeR::estimateDisp(edgeR_beads))

    edgeR_beads
}

#' Run edgeR for one sample against all the beads-only samples.
#'
#' @param object \code{\link[PhIPData]{PhIPData}} object
#' @param sample sample name of the sample to compare against beads-only samples
#' @param beads sample names for beads-only samples
#' @param common.disp edgeR estimated common disperion parameter
#' @param tagwise.dsip edgeR estimated tagwise dispersion parameter
#' @param trended.disp edgeR estimated trended dispersion parameter
.edgeR_one <- function(object, sample, beads,
                       common.disp, tagwise.disp, trended.disp){
    ## Coerce into edgeR object
    ## Set common dispersion to disp estimated from beads-only samples
    data <- as(object[, c(beads, sample)], "DGEList")
    data$common.dispersion <- common.disp
    data$tagwise.dispersion <- tagwise.disp
    data$trended.disp <- trended.disp

    ## edgeR output
    output <- edgeR::exactTest(data)$table
    # Convert to one-sided p-values and take -log10
    log10pval <- ifelse(output$logFC > 0, -log10(output$PValue/2),
                        -log10(1 - output$PValue/2))

    list(sample = sample,
         logfc = output$logFC,
         log10pval = log10pval)
}

#' Run edgeR on PhIP-Seq data
#'
#'
#' @param object \code{\link[PhIPData]{PhIPData}} object
#' @param threshold.cpm CPM threshold to be considered present in a sample
#' @param threshold.prevalence proportion of beads-only samples that surpass
#' \code{threshold.cpm}.
#' @param assay.names named vector specifying the assay names for the log2(fold-change) and exact test p-values. The vector must have entries for `logfc` and
#' `prob`.
#' @param parallel list of parameters specifying the parallelization strategy
#' as described in \code{\link{future::plan}{future::plan()}}.
#' @param beadsRR logical value specifying whether each beads-only sample
#' should be compared to all other beads-only samples.
#'
#' @return PhIPData object with log2 estimated fold-changes and p-values for
#' enrichment stored in the assays specified by `assay.names`.
#'
#' @seealso \code{\link[future]{plan}}
#' @seealso \code{\link{beadsRR}}
#'
#' @examples
#' sim_data <- readRDS(system.file("extdata", "sim_data.rds", package = "beer"))
#'
#' ## Sequential evaluation
#' edgeR(sim_data)
#'
#' ## Multisession evaluation
#' \dontrun{
#' edgeR(sim_data, parallel = list(strategy = "multisession"))
#' }
#'
#' @importFrom future.apply future_lapply
#' @importFrom future plan
#' @export
edgeR <- function(object, threshold.cpm = 0, threshold.prevalence = 0,
                  assay.names = c(logfc = "logfc", prob = "prob"),
                  parallel = list(strategy = "sequential"),
                  beadsRR = FALSE){

    ## Check assay names
    if(is.null(assay.names[["logfc"]])) assay.names[["logfc"]] <- "logfc"
    if(is.null(assay.names[["prob"]])) assay.names[["prob"]] <- "prob"

    ## Tidy parallelization inputs
    parallel <- .tidyParallel(parallel)
    ## Define plan and return current plan to oplan
    oplan <- do.call(plan, parallel)
    ## Rest to original plan on exit
    on.exit(plan(oplan))

    edgeR_beads <- .edgeRBeads(object, threshold.cpm, threshold.prevalence)
    common_disp <- edgeR_beads$common.dispersion
    tagwise_disp <- edgeR_beads$tagwise.dispersion
    trended_disp <- edgeR_beads$trended.dispersion

    ## Do beadsRR if necessary
    if(beadsRR){
        object <- beadsRR(object, method = "edgeR",
                          threshold.cpm, threshold.prevalence, assay.names)
    }

    ## Set-up output matrices ----------
    ## Make empty matrix for the cases where fc and prob do not exist
    empty_mat <- matrix(nrow = nrow(object), ncol = ncol(object))
    colnames(empty_mat) <- colnames(object)
    rownames(empty_mat) <- rownames(object)

    edgeR_fc <- if(assay.names[["logfc"]] %in% assayNames(object)){
        assay(object, assay.names[["logfc"]])
    } else { empty_mat }

    edgeR_pval <- if(assay.names[["prob"]] %in% assayNames(object)){
        assay(object, assay.names[["prob"]])
    } else { empty_mat }

    ## Run edgeR one-sample at a time ------------
    sample_names <- colnames(object[, object$group != getBeadsName()])
    beads_names <- colnames(object[, object$group == getBeadsName()])

    output <- future.apply::future_lapply(sample_names, function(sample){
        .edgeR_one(object, sample, beads_names,
                   common_disp, tagwise_disp, trended_disp)
    })

    # Unnest output items
    for(result in output){
        edgeR_fc[, result$sample] <- result$logfc
        edgeR_pval[, result$sample] <- result$log10pval
    }

    ## append to object
    assay(object, assay.names[["logfc"]]) <- edgeR_fc
    assay(object, assay.names[["prob"]]) <- edgeR_pval

    object
}

.getCommonDisp <- function(object, threshold.cpm = 1, threshold.prevalence = 2){
    edgeR_beads <- as(PhIPData::subsetBeads(object), "DGEList")

    ## For dispersion estimates, keep only peptides with cpm above the
    ## specified threshold in the given number of samples.
    keep_ind <- rowSums(edgeR::cpm(edgeR_beads) >= threshold.cpm) >=
        threshold.prevalence

    ## Estimate common dispersion in the beads-only data
    edgeR_beads <- edgeR_beads[keep_ind, , keep.lib.size = FALSE]
    edgeR_beads <- edgeR::calcNormFactors(edgeR_beads)
    edgeR_beads <- edgeR::estimateCommonDisp(edgeR_beads)
    edgeR_beads <- edgeR::estimateTagwiseDisp(edgeR_beads)

    edgeR_beads$common.dispersion
}

edgeR_one <- function(object, sample, beads, common.disp){
    ## Coerce into edgeR object
    ## Set common dispersion to disp estimated from beads-only samples
    data <- as(object[, c(beads, sample)], "DGEList")
    data <- edgeR::calcNormFactors(data)
    data$common.dispersion <- common.disp

    ## edgeR output
    output <- edgeR::exactTest(data)$table

    log10pval <- -log10(output$PValue)
    ## Truncate large p-values for storage reasons
    log10pval <- ifelse(log10pval > 250.1, 250.1, log10pval)
    ## Multiply by the sign of logfc so that negatively enriched peptides
    ## do not appear as significant
    log10pval <- sign(2^output$logFC)*log10pval

    list(sample = sample,
         logfc = output$logFC,
         log10pval = log10pval)
}

edgeR <- function(object, threshold.cpm = 1, threshold.prevalence = 2,
                  assay.names = c(logfc = "logfc", prob = "prob")){

    ## Derive dispersion estimates from beads-only samples
    common_disp <- .getCommonDisp(object, threshold.cpm, threshold.prevalence)

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

    for(sample in sample_names){
        output <- edgeR_one(object, sample, beads_names, common_disp)

        edgeR_fc[, output$sample] <- output$logfc
        edgeR_pval[, output$sample] <- output$log10pval
    }

    ## append to object
    assay(object, assay.names[["logfc"]]) <- edgeR_fc
    assay(object, assay.names[["prob"]]) <- edgeR_pval

    object
}

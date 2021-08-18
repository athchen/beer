.beadsRR_beer <- function(object, prior.params, beads.args, se.matrix,
                          jags.params, sample.dir){
    for(sample in colnames(object)){
        if(!jags.params$quiet) {
            print(paste0(which(colnames(object) == sample), " of ", ncol(object)))
        }

        # recode existing object
        one_beads <- object
        one_beads$group[colnames(one_beads) == sample] <- "sample"

        ## Calculate new beads-only priors
        new_beads <- if(prior.params$method == "custom"){
            lapply(beads.prior, function(x) x[!se.matrix[, sample]])
        } else {
            do.call(getAB, c(list(object = subsetBeads(one_beads),
                                  method = prior.params$method),
                             beads.args))
        }

        new_prior <- c(list(a_0 = new_beads[["a_0"]],
                            b_0 = new_beads[["b_0"]]),
                       prior.params[c("a_pi", "b_pi", "a_phi", "b_phi",
                                      "a_c", "b_c", "fc")])
        jags_run <- do.call(brew_one, c(list(object = one_beads,
                                             sample = sample,
                                             prior.params = new_prior),
                                        jags.params))

        saveRDS(jags_run, paste0(sample.dir,"/", sample, ".rds"))
    }
}

.beadsRR_edgeR <- function(object, threshold.cpm = 0, threshold.prevalence = 0,
                           assay.names = c(logfc = "logfc", prob = "prob")){

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

    beads_names <- colnames(subsetBeads(object))
    for(sample in beads_names){

        ## Recode existing object
        one_beads <- object
        one_beads$group[colnames(one_beads) == sample] <- "sample"

        ## Derive dispersion estimates from beads-only samples
        edgeR_beads <- .edgeRBeads(object, threshold.cpm, threshold.prevalence)
        common_disp <- edgeR_beads$common.dispersion
        tagwise_disp <- edgeR_beads$tagwise.dispersion
        trended_disp <- edgeR_beads$trended.dispersion

        ## Run edgeR for the one beads-only sample ------------
        output <- edgeR_one(one_beads, sample, beads_names[beads_names != sample],
                            common_disp, tagwise_disp, trended_disp)

        edgeR_fc[, output$sample] <- output$logfc
        edgeR_pval[, output$sample] <- output$log10pval
    }

    ## append to object
    assay(object, assay.names[["logfc"]]) <- edgeR_fc
    assay(object, assay.names[["prob"]]) <- edgeR_pval

    object
}

beadsRR <- function(object, method, ...){
    if(!method %in% c("edgeR", "beer")){
        stop("Invalid specified method for beads-only round robin.")
    } else if (method == "edgeR"){ .beadsRR_edgeR(object, ...)
    } else {.beadsRR_beer(object, ...)}
}

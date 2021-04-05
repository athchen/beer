beadsRR <- function(object, ...){
    for(sample in colnames(object)){
        if(!jags.params$quiet) {
            print(paste0(which(colnames(object) == sample), " of ", ncol(object)))
        }

        # recode existing object
        one_beads <- object
        one_beads$group[colnames(one_beads) == sample] <- "sample"

        ## Calculate new beads-only priors
        new_beads <- if(prior.params$method == "custom"){
            lapply(beads.prior, function(x) x[!se_peps[, sample]])
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

        saveRDS(jags_run, paste0(tmp_dir,"/", sample, ".rds"))
    }
}

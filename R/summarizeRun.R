summarizeRun_one <- function(object, directory, sample, se.matrix,
                             burn.in = 0, post.thin = 1){
    file <- paste0(directory, "/", sample, ".rds")

    rds <- readRDS(file)
    mcmc_matrix <- as.matrix(rds)
    iter_ind <- seq(burn.in + 1, nrow(mcmc_matrix), by = post.thin)
    mcmc_matrix <- mcmc_matrix[iter_ind, ]

    ## translate peptide indices
    pep_ind <- rep(NA, nrow(object))
    pep_ind[which(!se.matrix[, sample])] <- seq(sum(!se.matrix[, sample]))
    names(pep_ind) <- rownames(object)

    # for convenience extract parameter specific samples
    samples_c <- mcmc_matrix[, grepl("c", colnames(mcmc_matrix))]
    samples_pi <- mcmc_matrix[, grepl("pi", colnames(mcmc_matrix))]
    samples_phi <- mcmc_matrix[, grepl("phi\\[", colnames(mcmc_matrix))]
    samples_Z <- mcmc_matrix[, grepl("Z\\[", colnames(mcmc_matrix))]

    # summarize info
    point_c <- data.frame(parameter = "c",
                          sample = sample,
                          est_value = mean(samples_c))
    point_pi <- data.frame(parameter = "pi",
                           sample = sample,
                           est_value = mean(samples_pi))
    point_phi <- data.frame(parameter = "phi",
                            sample = sample,
                            peptide = rownames(object),
                            est_value = unname(colMeans(samples_phi)[pep_ind]),
                            est_enriched = unname(colMeans(samples_phi*
                                                               samples_Z)
                                                  [pep_ind]))

    point_Z <- data.frame(parameter = "Z",
                          sample = sample,
                          peptide = rownames(object),
                          est_value = unname(colMeans(samples_Z)[pep_ind]))

    return(list(point_c = point_c,
                point_pi = point_pi,
                point_phi = point_phi,
                point_Z = point_Z))
}

summarizeRun <- function(object, directory, se.matrix,
                         burn.in = 0, post.thin = 1,
                         assay.names,
                         quiet = TRUE){

    jags_files <- list.files(directory, full.names = TRUE)

    ## Pre-allocate containers
    point_c <- if(assay.names["c"] %in% names(metadata(object))) {
        metadata(object)[[assay.names["c"]]]
    } else rep(NA, ncol(object))
    point_pi <- if(assay.names["pi"] %in% names(metadata(object))) {
        metadata(object)[[assay.names["pi"]]]
    } else rep(NA, ncol(object))
    point_phi <- if(!is.null(assay.names["phi"]) &
                    assay.names["phi"] %in% assayNames(object)){
        assay(object, assay.names["phi"])
    } else matrix(NA, nrow = nrow(object), ncol = ncol(object))
    point_phi <- if(!is.null(assay.names["phi"]) &
                    assay.names["phi"] %in% assayNames(object)){
        assay(object, assay.names["phi"])
    } else matrix(NA, nrow = nrow(object), ncol = ncol(object))
    point_phi_Z <- if(!is.null(assay.names["phi_Z"]) &
                    assay.names["phi_Z"] %in% assayNames(object)){
        assay(object, assay.names["phi_Z"])
    } else matrix(NA, nrow = nrow(object), ncol = ncol(object))
    point_Z <- if(!is.null(assay.names["Z"]) &
                    assay.names["Z"] %in% assayNames(object)){
        assay(object, assay.names["Z"])
    } else matrix(NA, nrow = nrow(object), ncol = ncol(object))

    names(point_c) <- names(point_pi) <- colnames(point_phi) <-
        colnames(point_phi_Z) <- colnames(point_Z) <-
        colnames(object)

    rownames(point_phi) <- rownames(point_phi_Z) <- rownames(point_Z) <-
        rownames(object)

    for(file in jags_files){
        if(!quiet){
            print(paste0(which(file == jags_files), " of ",
                         length(jags_files)))
        }

        pattern <- paste0(directory, "/(.*)\\.rds")
        sample <- regmatches(file, regexec(pattern, file))[[1]][2]

        out <- summarizeRun_one(object, directory, sample, se.matrix,
                                burn.in, post.thin)

        point_c[sample] <- out$point_c$est_value
        point_pi[sample] <- out$point_pi$est_value
        point_phi[, sample] <- out$point_phi$est_value
        point_phi_Z[, sample] <- out$point_phi$est_enriched
        point_Z[, sample] <- out$point_Z$est_value
    }

    ## Assign c and pi to metadata
    add <- !vapply(assay.names[c("c", "pi")], is.na, logical(1))
    metadata(object) <- c(metadata(object),
                          list(c = point_c, pi = point_pi)[add])

    ## Assign phi, phi_Z, and Z to assays
    assay <- c("phi", "phi_Z", "Z")[!is.na(assay.names[c("phi", "phi_Z", "Z")])]
    assays(object)[assay] <- list(phi = point_phi, phi_Z = point_phi_Z,
                                  Z = point_Z)[assay]
    object
}

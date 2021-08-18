#' Derive point estimates for c, pi, phi, and Z for a particular sample
#'
#' Posterior means are used as point estimates for \eqn{c}, \eqn{\pi},
#' \eqn{\phi}, and \eqn{Z}. As super-enriched peptides are tossed out before
#' MCMC sampling, super-enriched peptides return \code{NA} for the \eqn{\phi}
#' and \eqn{Z} point estimates. Indices corresponding to a particular peptide in
#' the MCMC sampler are mapped back to the original peptide names.
#'
#' @param object a \code{\link[PhIPData]{PhIPData}} object
#' @param file path to rds file
#' @param se.matrix logical matrix indicating which peptides were identified as
#' super-enriched peptides
#' @param burn.in number of iterations to be burned
#' @param post.thin thinning parameter
#'
#' @return list of point estimates for c, pi, phi and Z
summarizeRun_one <- function(object, file, se.matrix,
                             burn.in = 0, post.thin = 1){
    sample <- regmatches(file, regexec("/([^/]*)\\.rds", file))[[1]][2]

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
                            est_enriched = unname((colSums(samples_phi*samples_Z)/
                                                      colSums(samples_Z))
                                                  [pep_ind]))
    point_Z <- data.frame(parameter = "Z",
                          sample = sample,
                          peptide = rownames(object),
                          est_value = unname(colMeans(samples_Z)[pep_ind]))

    list(point_c = point_c,
         point_pi = point_pi,
         point_phi = point_phi,
         point_Z = point_Z)
}

#' Summarize MCMC chain and return point estimates for BEER parameters
#'
#' Posterior means are used as point estimates for \eqn{c}, \eqn{\pi},
#' \eqn{\phi}, and \eqn{Z}. As super-enriched peptides are tossed out before
#' MCMC sampling, super-enriched peptides return \code{NA} for the \eqn{\phi}
#' and \eqn{Z} point estimates. Indices corresponding to a particular peptide in
#' the MCMC sampler are mapped back to the original peptide names.
#'
#' @param object a \code{\link[PhIPData]{PhIPData}} object
#' @param directory path of the directory containing JAGS output
#' @param se.matrix logical matrix indicating which peptides were identified as
#' super-enriched peptides
#' @param burn.in number of iterations to be burned
#' @param post.thin thinning parameter
#' @param assay.names named vector of specifying where to store point estimates
#'
#' @return PhIPData object with point estimates stored in the assays specified
#' by `assay.names`.
summarizeRun <- function(object, directory, se.matrix,
                         burn.in = 0, post.thin = 1,
                         assay.names, quiet = FALSE){

    jags_files <- list.files(directory, full.names = TRUE)
    samples <- vapply(
        regmatches(jags_files, regexec("/([^/]*)\\.rds", jags_files)),
        function(x) x[[2]],
        character(1))
    names(jags_files) <- samples

    ## Pre-allocate containers
    point_c <- if(assay.names["c"] %in% colnames(sampleInfo(object))) {
        sampleInfo(object)[[assay.names["c"]]]
    } else rep(NA, ncol(object))
    point_pi <- if(assay.names["pi"] %in% colnames(sampleInfo(object))) {
        sampleInfo(object)[[assay.names["pi"]]]
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

    ## Summarize each file, use future_lapply here to enable reading/import
    ## to be parallelized
    progressr::handlers("txtprogressbar")
    p <- progressr::progressor(along = jags_files)
    files_out <- future_lapply(jags_files, function(file){
        file_counter <- paste0(which(file == jags_files), " of ",
                               length(jags_files))
        p(file_counter, class = "sticky", amount = 1)
        summarizeRun_one(object, file, se.matrix, burn.in, post.thin)
    })

    for(out in files_out){
        sample <- out$point_c$sample

        point_c[sample] <- out$point_c$est_value
        point_pi[sample] <- out$point_pi$est_value
        point_phi[, sample] <- out$point_phi$est_value
        point_phi_Z[, sample] <- out$point_phi$est_enriched
        point_Z[, sample] <- out$point_Z$est_value
     }

    ## Assign c and pi to sampleInfo
    if(!is.na(assay.names["c"])) object$c <- point_c
    if(!is.na(assay.names["pi"])) object$pi <- point_pi

    ## Assign phi, phi_Z, and Z to assays
    assay <- c("phi", "phi_Z", "Z")[!is.na(assay.names[c("phi", "phi_Z", "Z")])]
    assays(object)[assay.names[assay]] <- list(phi = point_phi, phi_Z = point_phi_Z,
                                  Z = point_Z)[assay]
    object
}

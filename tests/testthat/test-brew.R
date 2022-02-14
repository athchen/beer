sim_data <- readRDS(system.file("extdata", "sim_data.rds", package = "beer"))

test_that("Prior parameters are correctly tidied", {
    prior.params <- list(
        a_pi = 2, b_pi = 300, a_phi = 1.25, b_phi = 0.1,
        a_c = 80, b_c = 20, fc = 1
    )
    beads_args <- list(lower = 1)
    tidied_inputs <- .tidyInputsPrior(prior.params, sim_data, beads_args)

    ## Missing specified method defaults to edgeR
    expect_equal(tidied_inputs$method, "edgeR")

    ## a_0, b_0 are calculated accordingly
    expect_true(all(c("a_0", "b_0") %in% names(tidied_inputs)))

    ## fc defaults to 1
    prior.params$fc <- NULL
    expect_equal(.tidyInputsPrior(prior.params, sim_data, beads_args)$fc, 1)

    ## Missing priors returns an error
    prior.params$a_phi <- NULL
    expect_error(
        .tidyInputsPrior(prior.params, sim_data, beads_args),
        "elements are missing"
    )
    prior.params$a_phi <- 1.25

    ## Invalid specified method
    prior.params$method <- "test"
    expect_error(
        .tidyInputsPrior(prior.params, sim_data, beads_args),
        "Invalid specified method"
    )

    ## No extra parameters are returned
    prior.params$method <- "edgeR"
    prior.params$extra <- "extra"
    tidied_inputs <- .tidyInputsPrior(prior.params, sim_data, beads_args)
    expect_false("extra" %in% names(tidied_inputs))
})

test_that("SE parameters are correctly tidied", {
    se_params <- list()
    beads_prior <- getAB(subsetBeads(sim_data), method = "edgeR")

    ## Missing method defaults to mle
    expect_equal(.tidyInputsSE(se_params, beads_prior)$method, "mle")

    ## Error for invalid method
    se_params$method <- "test"
    expect_error(
        .tidyInputsSE(se_params, beads_prior),
        "Invalid specified method"
    )

    ## Add correct defaults
    se_params$method <- "edgeR"
    expect_equal(
        names(.tidyInputsSE(se_params, beads_prior)),
        c("method", "threshold", "fc.name")
    )
    se_params$method <- "mle"
    expect_equal(
        names(.tidyInputsSE(se_params, beads_prior)),
        c("method", "threshold", "beads.prior")
    )

    ## No extra parameters
    se_params <- list(
        method = "edgeR", threshold = 15, fc.name = "prob",
        extra = "extra"
    )
    expect_equal(
        names(.tidyInputsSE(se_params, beads_prior)),
        c("method", "threshold", "fc.name")
    )
    se_params <- list(
        method = "mle", threshold = 5,
        beads.prior = "beasd.prior", extra = "extra"
    )
    expect_equal(
        names(.tidyInputsSE(se_params, beads_prior)),
        c("method", "threshold", "beads.prior")
    )
})

test_that("JAGS parameters are correctly tidied", {
    jags.params <- list(
        n.chains = 1, n.adapt = 1e3,
        n.iter = 1e4, thin = 1, na.rm = TRUE,
        burn.in = 0, post.thin = 1,
        seed = as.numeric(format(Sys.Date(), "%Y%m%d"))
    )

    ## Missing params
    jags.params$post.thin <- NULL
    expect_equal(.tidyInputsJAGS(jags.params)$post.thin, 1)

    ## Extra params
    jags.params$extra <- "extra"
    expect_equal(
        names(.tidyInputsJAGS(jags.params)),
        c(
            "n.chains", "n.adapt", "n.iter", "thin", "na.rm",
            "burn.in", "post.thin", "seed"
        )
    )
})

test_that("Assay names are correctly tidied", {
    assay.names <- c(
        phi = NULL, phi_Z = "logfc", Z = "prob",
        c = "sampleInfo", pi = "sampleInfo"
    )
    ## Missing parameters
    expect_equal(.tidyAssayNames(assay.names)[["phi"]], NA_character_)

    ## Invalid c/pi options
    assay.names[["c"]] <- "test"
    expect_error(.tidyAssayNames(assay.names), "Invalid location specified")
    assay.names[["pi"]] <- "test"
    expect_error(.tidyAssayNames(assay.names), "Invalid location specified")
    assay.names[c("c", "pi")] <- rep("sampleInfo", 2)

    ## Duplicate assay names
    assay.names[c("phi", "phi_Z")] <- rep("logfc", 2)
    expect_error(.tidyAssayNames(assay.names), "unique assay names")

    ## Extra parameters
    assay.names[["extra"]] <- "extra"
    assay.names["phi"] <- NA
    expect_equal(
        names(.tidyAssayNames(assay.names)),
        c("phi", "phi_Z", "Z", "c", "pi")
    )
})

test_that("warns when overwriting sampleInfo", {
    sim_data$c <- "overwrite"
    expect_true(any(
        grepl(
            "Values in the following assays will be overwritten: sampleInfo",
            cli::cli_format_method(
                brew(sim_data, jags.params = list(seed = 123))
            )
        )
    ))

    ## Clean sampleInfo
    sampleInfo(sim_data)[, "c"] <- NULL
})

test_that("brew can run with different BiocParallelParam classes", {

    ## Serial
    ## Also check that files are saved in the right directory
    ex_dir <- paste0(system.file("extdata", package = "beer"), "/ex_dir")
    # if exists delete
    if (dir.exists(ex_dir)) unlink(ex_dir, recursive = TRUE)
    brew_seq <- brew(sim_data,
        sample.dir = ex_dir,
        jags.params = list(seed = 123),
        BPPARAM = BiocParallel::SerialParam()
    )
    expect_equal(list.files(ex_dir), paste0(c(10, 5:9), ".rds"))
    # Clean directory
    unlink(ex_dir, recursive = TRUE)

    ## Snow
    suppressWarnings(brew_snow <- brew(sim_data,
        jags.params = list(seed = 123),
        BPPARAM = BiocParallel::SnowParam()
    ))

    ## Check that there's nothing different
    expect_identical(brew_seq, brew_snow)
})

test_that("brew works with beadsRR", {
    brew_rr <- brew(sim_data, jags.params = list(seed = 123), beadsRR = TRUE)
    expect_s4_class(brew_rr, "PhIPData")
})

test_that("Files are saved to the correct directory", {
    ## Also check that files are saved in the right directory
    ex_dir <- paste0(system.file("extdata", package = "beer"), "/ex_dir")
    # if exists delete
    if (dir.exists(ex_dir)) unlink(ex_dir, recurive = TRUE)
    brew(sim_data,
        jags.params = list(seed = 123),
        sample.dir = ex_dir, beadsRR = TRUE
    )

    expect_equal(length(list.files(ex_dir)), ncol(sim_data))
    # Clean directory
    unlink(ex_dir, recursive = TRUE)
})

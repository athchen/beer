sim_data <- readRDS(system.file("extdata", "sim_data.rds", package = "beer"))
edgeR_out <- edgeR(sim_data)

test_that("Clearly enriched peptides are accurately identified with edgeR", {

    ## Test that it works
    expect_equal(guessEnriched(edgeR_out, method = "edgeR"),
                 assay(edgeR_out, "logfc") > log2(15))

    ## Missing assay
    expect_error(guessEnriched(edgeR_out, method = "edgeR", fc.name = "test"),
                 "assay does not exist")

    ## Empty assay
    assay(edgeR_out, "test") <- matrix(NA, nrow = nrow(edgeR_out),
                                       ncol = ncol(edgeR_out))
    expect_error(guessEnriched(edgeR_out, method = "edgeR", fc.name = "test"),
                 "is empty.")
})

test_that("MLE super-enriched peptides are all enriched", {
    mle_se <- guessEnriched(sim_data, method = "mle",
                            beads.prior = getAB(sim_data, method = "edgeR"))
    expect_true(all(which(mle_se) %in% which(assay(sim_data, "true_Z") == 1)))
})

test_that("Wrapper function catched method specification error", {
    expect_error(guessEnriched(sim_data, method = "test"), "Invalid method.")
})

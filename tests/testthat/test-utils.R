sim_data <- readRDS(system.file("extdata", "sim_data.rds", package = "beer"))

test_that("Counts matrices are properly checked for integer values", {
    expect_silent(.checkCounts(sim_data))

    counts(sim_data) <- matrix(runif(prod(dim(sim_data)), 0, 1000),
        nrow = nrow(sim_data)
    )
    expect_error(.checkCounts(sim_data), "`counts` entries must be integers.")
})

test_that("Non-empty assays are correctly identified", {
    expect_equal(
        unname(.checkOverwrite(sim_data, c("counts", "logfc"))),
        c(TRUE, FALSE)
    )
})

test_that("Overwritten assays are identified", {
    tmp_data <- sim_data
    assay.names <- c(
        phi = NA, phi_Z = "logfc", Z = "prob",
        c = "sampleInfo", pi = "sampleInfo"
    )

    ## Non-empty sampleInfo
    tmp_data$c <- tmp_data$pi <- rep("test", ncol(sim_data))
    expect_equal(
        unname(.checkOverwrite(tmp_data, assay.names)),
        c(FALSE, FALSE, FALSE, TRUE, TRUE)
    )

    # Non-empty assay
    logfc(tmp_data) <- matrix(1, nrow = nrow(tmp_data), ncol = ncol(tmp_data))
    expect_equal(
        unname(.checkOverwrite(tmp_data, assay.names)),
        c(FALSE, TRUE, FALSE, TRUE, TRUE)
    )

    # No problems
    tmp_data$c <- tmp_data$pi <- NULL
    logfc(tmp_data) <- matrix(NA, nrow = nrow(tmp_data), ncol = ncol(tmp_data))
    expect_equal(
        unname(.checkOverwrite(tmp_data, assay.names)),
        c(FALSE, FALSE, FALSE, FALSE, FALSE)
    )
})

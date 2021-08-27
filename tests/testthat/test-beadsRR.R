sim_data <- readRDS(system.file("extdata", "sim_data.rds", package = "beer"))

test_that("beadsRR only accepts two methods", {
    expect_error(beadsRR(sim_data, method = "test"), "Invalid specified method")
})

test_that("beadsRR works with BEER", {
    expect_snapshot(beadsRR(sim_data,
        jags.params = list(seed = 123),
        method = "beer"
    ))

    # Check that it works after changing plans
    beadsRR_seq <- beadsRR(sim_data,
        method = "beer",
        jags.params = list(seed = 123),
        summarize = FALSE
    )
    plan("multisession", workers = future::availableCores() - 1)
    beadsRR_par <- beadsRR(sim_data,
        method = "beer",
        jags.params = list(seed = 123),
        summarize = FALSE
    )
    expect_equal(length(unique(beadsRR_seq)), 1)
    expect_equal(
        unname(length(unique(beadsRR_par))),
        unname(future::availableCores() - 1)
    )

    ## reset plan
    plan("sequential")
})

test_that("beadsRR works with edgeR", {
    beadsRR_seq <- beadsRR(sim_data, method = "edgeR")
    plan("multisession")
    beadsRR_par <- beadsRR(sim_data, method = "edgeR")
    expect_identical(beadsRR_seq, beadsRR_par)

    ## reset plan
    plan("sequential")
})

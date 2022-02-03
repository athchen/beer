sim_data <- readRDS(system.file("extdata", "sim_data.rds", package = "beer"))

test_that("beadsRR only accepts two methods", {
    expect_error(beadsRR(sim_data, method = "test"), "Invalid specified method")
})

test_that("beadsRR works with BEER", {
    expect_snapshot(beadsRR(sim_data,
        jags.params = list(seed = 123),
        method = "beer"
    ))

    ## Test parallelization
    beer_ser <- beadsRR(sim_data,
                        method = "beer",
                        jags.params = list(seed = 123),
                        bp.param = BiocParallel::SerialParam())
    suppressWarnings(beer_snow <- beadsRR(sim_data, method = "beer",
                                          jags.params = list(seed = 123),
                                          bp.param = BiocParallel::SnowParam()))
    expect_identical(beer_ser, beer_snow)
})

test_that("beadsRR works with edgeR", {
    expect_snapshot(beadsRR(sim_data, method = "edgeR"))

    ## Test parallelization
    beadsRR_ser <- beadsRR(sim_data, method = "edgeR",
                           bp.param = BiocParallel::SerialParam())
    beadsRR_snow <- beadsRR(sim_data, method = "edgeR",
                            bp.param = BiocParallel::SnowParam())
    expect_identical(beadsRR_ser, beadsRR_snow)
})

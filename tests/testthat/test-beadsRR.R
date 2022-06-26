sim_data <- readRDS(system.file("extdata", "sim_data.rds", package = "beer"))

test_that("beadsRR only accepts two methods", {
    expect_error(beadsRR(sim_data, method = "test"), "Invalid specified method")
})

test_that("beadsRR works with BEER", {
    expect_s4_class(beadsRR(sim_data,
        jags.params = list(seed = 123),
        method = "beer"
    ), "PhIPData")

    ## Test parallelization
    beer_ser <- beadsRR(sim_data,
        method = "beer",
        jags.params = list(seed = 123),
        BPPARAM = BiocParallel::SerialParam()
    )
    suppressWarnings(beer_snow <- beadsRR(sim_data,
        method = "beer",
        jags.params = list(seed = 123),
        BPPARAM = BiocParallel::SnowParam()
    ))
    expect_identical(beer_ser, beer_snow)
})

test_that("beadsRR works with edgeR", {
    ## Test exactTest
    expect_s4_class(beadsRR(sim_data, method = "edgeR"), "PhIPData")

    ## Test parallelization
    beadsRR_ser <- beadsRR(sim_data,
        method = "edgeR",
        BPPARAM = BiocParallel::SerialParam()
    )
    beadsRR_snow <- beadsRR(sim_data,
        method = "edgeR",
        BPPARAM = BiocParallel::SnowParam()
    )
    expect_identical(beadsRR_ser, beadsRR_snow)

    ## Test glmQLF
    expect_s4_class(
        beadsRR(sim_data, method = "edgeR", de.method = "glmQLFTest"),
        "PhIPData"
    )

    ## Test parallelization
    beadsRR_ser <- beadsRR(sim_data,
        method = "edgeR", de.method = "glmQLFTest",
        BPPARAM = BiocParallel::SerialParam()
    )
    beadsRR_snow <- beadsRR(sim_data,
        method = "edgeR", de.method = "glmQLFTest",
        BPPARAM = BiocParallel::SnowParam()
    )
    expect_identical(beadsRR_ser, beadsRR_snow)
})

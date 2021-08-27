sim_data <- readRDS(system.file("extdata", "sim_data.rds", package = "beer"))

test_that("Estimated initial values are consistent", {
    expect_snapshot(guessInits(sim_data, getAB(sim_data, "edgeR")))
})

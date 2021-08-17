sim_data <- readRDS(system.file("extdata", "sim_data.rds", package = "beer"))

test_that("Counts matrices are properly checked for integer values", {
    expect_silent(.checkCounts(sim_data))

    counts(sim_data) <- matrix(runif(prod(dim(sim_data)), 0, 1000),
                               nrow = nrow(sim_data))
    expect_error(.checkCounts(sim_data), "`counts` entries must be integers.")
})

test_that("Non-empty assays are correctly identified", {
    expect_equal(unname(.checkOverwrite(sim_data, c("counts", "logfc"))),
                 c(TRUE, FALSE))
})

test_that("Parallelization paramers are correctly tidied", {

    ## Convert vector to list
    parallel <- "sequential"
    expect_equal(.tidyParallel(parallel)[["strategy"]], "sequential")

    ## Invalid plan
    parallel[["strategy"]] <- "test"
    expect_error(.tidyParallel(parallel),
                 "The specified parallel strategy is not supported.")
    parallel[["strategy"]] <- "sequential"

    ## Check workers default for sequential plan
    max_workers <- unname(future::availableCores() - 1)
    expect_equal(.tidyParallel(parallel)[["workers"]], NULL)

    ## Check workers default for multisession plan
    parallel[["strategy"]] <- "multisession"
    expect_equal(.tidyParallel(parallel)[["workers"]], max_workers)

    ## Poorly specified # of workers
    parallel[["workers"]] <- -1
    expect_equal(.tidyParallel(parallel)[["workers"]], 1)
    parallel[["workers"]] <- max_workers + 4
    expect_equal(.tidyParallel(parallel)[["workers"]], max_workers)
    parallel[["workers"]] <- 2
    expect_equal(.tidyParallel(parallel)[["workers"]], 2)

})

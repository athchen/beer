sim_data <- readRDS(file.path("..", "..", "inst/extdata/sim_data.rds"))

test_that("edgeR runs with multisession and sequential evaluation", {
    curr_plan <- plan()

    ## Sequential, check that current plan is properly reset
    edgeR_seq <- edgeR(sim_data, parallel = list(strategy = "sequential"))
    expect_identical(plan(), curr_plan)

    ## Multisession, check that current plan is properly reset
    edgeR_multi <- edgeR(sim_data, parallel = list(strategy = "multisession"))
    expect_identical(plan(), curr_plan)

    ## Check that there's nothing different
    expect_identical(edgeR_seq, edgeR_multi)
})

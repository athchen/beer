sim_data <- readRDS(system.file("extdata", "sim_data.rds", package = "beer"))

test_that("edgeR corrects for poorly specified assay.names", {

    ## Unnamed assay.names vector
    edgeR_out <- runEdgeR(sim_data, assay.names = c("logfc", "prob"))
    expect_equal(
        unname(.checkOverwrite(edgeR_out, c("logfc", "prob"))),
        c(TRUE, TRUE)
    )

    ## Shorter vector length
    edgeR_out <- runEdgeR(sim_data, assay.names = c("logfc"))
    expect_equal(
        unname(.checkOverwrite(edgeR_out, c("logfc", "prob"))),
        c(TRUE, TRUE)
    )
})

cli::test_that_cli("warns when overwriting matrices", {
    expect_snapshot(runEdgeR(sim_data, assay.names = c("logfc", "counts")))
})

test_that("edgeR runs with multisession and sequential evaluation", {
    curr_plan <- plan()

    ## Sequential, check that current plan is properly reset
    edgeR_seq <- runEdgeR(sim_data, parallel = "sequential")
    expect_identical(plan(), curr_plan)

    ## Multisession, check that current plan is properly reset
    edgeR_multi <- runEdgeR(sim_data, parallel = "multisession")
    expect_identical(plan(), curr_plan)

    ## Check that there's nothing different
    expect_identical(edgeR_seq, edgeR_multi)
})

test_that("edgeR works with beadsRR", {
    expect_snapshot(runEdgeR(sim_data, beadsRR = TRUE))
})

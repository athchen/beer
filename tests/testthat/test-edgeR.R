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

test_that("warns when overwriting matrices", {
    expect_true(
        grepl(
            "Values in the following assays will be overwritten: counts",
            cli::cli_format_method(
                runEdgeR(sim_data, assay.names = c("logfc", "counts"))
            )
        )
    )
})

test_that("edgeR runs with different BiocParallelParam classes", {

    ## Serial
    edgeR_ser <- runEdgeR(sim_data, bp.param = BiocParallel::SerialParam())
    ## Snow
    suppressWarnings(
        edgeR_snow <- runEdgeR(sim_data, bp.param = BiocParallel::SnowParam())
    )

    ## Check that there's nothing different
    expect_identical(edgeR_ser, edgeR_snow)
})

test_that("edgeR works with beadsRR", {
    expect_s4_class(runEdgeR(sim_data, beadsRR = TRUE), "PhIPData")
})

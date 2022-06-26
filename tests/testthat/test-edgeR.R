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
    edgeR_ser <- runEdgeR(sim_data, BPPARAM = BiocParallel::SerialParam())
    ## Snow
    suppressWarnings(
        edgeR_snow <- runEdgeR(sim_data, BPPARAM = BiocParallel::SnowParam())
    )

    ## Check that there's nothing different
    expect_identical(edgeR_ser, edgeR_snow)
})

test_that("edgeR works with beadsRR", {
    exact <- runEdgeR(sim_data, beadsRR = TRUE)
    glmQLF <- runEdgeR(sim_data, beadsRR = TRUE, de.method = "glmQLFTest")

    expect_snapshot(exact)
    expect_snapshot(glmQLF)

    expect_false(identical(exact, glmQLF))
})

test_that("runEdgeR captures invalid inputs", {

    # Error for invalid method
    expect_error(
        runEdgeR(sim_data, de.method = "invalid"),
        "Invalid edgeR method for identifying DE peptides."
    )
})

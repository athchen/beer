sim_data <- readRDS(system.file("extdata", "sim_data.rds", package = "beer"))

test_that("Expected read counts and proportion of reads can be calculated", {

    ## Test invalid option
    expect_error(
        getExpected(sim_data, "test"),
        "Invalid `type` specified. `type` must be any of `rc` or `prop`."
    )

    ## Test read counts
    rc_out <- getExpected(sim_data, "rc")
    expect_snapshot(rc_out)
    expect_true("expected_rc" %in% assayNames(rc_out) &
            !"expected_prop" %in% assayNames(rc_out))

    ## Test overwrite
    expect_true(any(
        grepl(
            "assay exists and will be overwritten",
            cli::cli_format_method(
                getExpected(rc_out, "rc")
            )
        )
    ))

    ## Test prop
    prop_out <- getExpected(sim_data, "prop", "expected_prop")
    expect_true(all(assay(prop_out, "expected_prop") <= 1 &
        assay(prop_out, "expected_prop") >= 0))
    
    ## Test both
    both_out <- getExpected(sim_data)
    expect_true(all(c("expected_rc", "expected_prop") %in% assayNames(both_out)))
    
    ## Test one invalid option
    one_invalid <- getExpected(sim_data, c("rc", "test"))
    expect_true(
        "expected_rc" %in% assayNames(one_invalid) &
            !"expected_prop" %in% assayNames(one_invalid)
    )
})

test_that("Bayes factors can be calculated", {
    
    ## Test invalid assay name
    expect_error(
        getBF(sim_data, assay.postprob = "test"),
        "Posterior probability assay not found."
    )
    
    ## Test Bayes factor
    brew_out <- brew(sim_data, beadsRR = TRUE)
    bf_out <- getBF(brew_out, assay.postprob = "prob", "bayes_factors")
    expect_snapshot(bf_out)
    
    ## Test different priors
    expect_false(identical(bf_out, getBF(brew_out, prior.params = list(a_pi = 10, b_pi = 100))))
    
    ## Test overwrite
    expect_true(any(
        grepl(
            "The specified assay exists and will be overwritten",
            cli::cli_format_method(
                getBF(bf_out, assay.name = "bayes_factors")
            )
        )
    ))
})

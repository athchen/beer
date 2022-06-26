sim_data <- readRDS(system.file("extdata", "sim_data.rds", package = "beer"))

test_that("Expected read counts and proportion of reads can be calculated", {

    ## Test invalid option
    expect_error(
        getExpected(sim_data, "test"),
        "Invalid `type` specified. `type` must be one of `rc` or `prop`."
    )

    ## Test read counts
    rc_out <- getExpected(sim_data, "rc")
    expect_snapshot(rc_out)

    ## Test overwrite
    expect_true(any(
        grepl(
            "The specified assay exists and will be overwritten",
            cli::cli_format_method(
                getExpected(rc_out, "rc")
            )
        )
    ))

    ## Test prop
    prop_out <- getExpected(sim_data, "prop", "expected_prop")
    expect_true(all(assay(prop_out, "expected_prop") <= 1 &
        assay(prop_out, "expected_prop") >= 0))
})

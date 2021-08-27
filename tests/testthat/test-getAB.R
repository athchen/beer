sim_data <- readRDS(system.file("extdata", "sim_data.rds", package = "beer"))

test_that("Valid beta parameter estimates are derived from PhIPData objects", {
    expect_silent(getAB(sim_data, "mom"))
    expect_silent(getAB(sim_data, "mle"))
    expect_silent(getAB(sim_data, "edgeR"))

    # Test invalid method
    expect_error(
        getAB(sim_data, "test"),
        "Invalid specified method for estimating a, b."
    )
})

test_that("Valid beta parameter estimates are derived from vectors", {
    prop <- rbeta(100, 2, 8)

    expect_type(getAB(prop, "mom"), "double")
    expect_type(getAB(prop, "mle"), "double")

    # Test invalid method
    expect_error(
        getAB(prop, "edgeR"),
        "edgeR is not a valid method for vectors."
    )
    expect_error(
        getAB(prop, "test"),
        "Invalid specified method for estimating a, b."
    )
})

test_that("getAB fails for misspecified vectors and incorrect object types", {
    ## Above 1
    expect_error(getAB(1:10), "Invalid inputs. Vectors ")

    ## Below 0
    expect_error(getAB(-10:-1), "Invalid inputs. Vectors ")

    ## Invalid vector type
    expect_error(getAB(c("a", "b", "c")), "Invalid inputs. Vectors ")

    ## Invalid object type
    expect_error(
        getAB(data.frame(x = 1:10, y = -10:-1)),
        "Invalid inputs. 'object'"
    )
})

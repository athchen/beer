#' @title Derive initial estimates of unknown model parameters
#'
#' @description To reduce converge time and to reduce the likelihood of the
#' slice sampler getting stuck, we use maximum likelihood to derive initial
#' estimates for unknown model parameters.
#'
#' @details Briefly initial values are defined as follows:
#' \enumerate{
#'     \item \code{theta_guess[i, j] = Y[i, j]/n[j]}, or the the MLE for theta.
#'     \item \code{Z_guess[i, j] = 1} if \eqn{j} is a serum sample, and the
#'     observed read count is >2x the expected read count assuming
#'     \code{c[j] = 1}.
#'     \item \code{pi_guess[j]} is the mean of column \eqn{j} in \code{Z_guess}.
#'     \item \code{c_guess[j]} is the estimated slope from regressing the
#'     observed read counts against the expected read counts (without adjusting
#'     for the attenuation constant) for non-enriched peptides only.
#'     \item \code{phi_guess[i,j]} is the ratio of the observed read counts to
#'     the expected read counts multiplied by the attenuation constant.
#' }
#'
#' @param object a \code{\link[PhIPData]{PhIPData}} object
#' @param beads.prior a data frame with two columns (named a_0, b_0) containing
#' estimated shape parameters from beads-only samples.
#'
#' @return a list of estimated initial values.
#'
#' @seealso \emph{Methods} in [Chen et. al 2022](https://www.biorxiv.org/content/10.1101/2022.01.19.476926v1)
#'
#' @import PhIPData
#' @importFrom stats coef lm
guessInits <- function(object, beads.prior) {
    ## Extract convenient parameters
    N <- ncol(object)
    P <- nrow(object)
    B <- ifelse(object$group == getBeadsName(), 1, 0)
    n <- librarySize(object)
    Y <- counts(object)

    a_0 <- beads.prior[["a_0"]]
    b_0 <- beads.prior[["b_0"]]

    # Guess proportions, add offset of 1e-8 when theta is 0
    theta_guess <- Y / matrix(rep(n, P), nrow = P, byrow = TRUE)
    theta_guess <- theta_guess + 1e-8 * (theta_guess == 0)

    ## Guess enriched peptides
    expected_prop <- a_0 / (a_0 + b_0)
    expected_rc <- vapply(
        n, function(n_i) n_i * expected_prop,
        numeric(length(expected_prop))
    )
    Z_guess <- (Y > expected_rc * 2) *
        (1 - matrix(rep(B, P), nrow = P, byrow = TRUE))

    ## Guess proportions of enriched peptides
    pi_guess <- colMeans(Z_guess)

    ## Guess attenuation constant
    c_guess <- vapply(seq(N), function(index) {
        if (B[index] != 1) {
            coef(lm(Y[, index] ~ expected_rc[, index] - 1))
        } else {
            1
        }
    }, numeric(1))

    ## Guess fold change
    expected_rc_attn <- vapply(
        c_guess * n, function(n_i) n_i * expected_prop,
        numeric(length(expected_prop))
    )
    phi_guess <- (1 - Z_guess) + Y / expected_rc_attn * Z_guess

    list(
        theta = theta_guess,
        c_sample = vapply(c_guess, min, numeric(1), 1 - 1e-6), # ensure c < 1
        pi_sample = vapply(pi_guess, max, numeric(1), 1e-6), # ensure pi > 0
        Z = Z_guess,
        phi_sample = phi_guess
    )
}

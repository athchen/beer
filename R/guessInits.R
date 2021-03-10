# ----- guess_inits.R --------------
# Function to give crude initial values for JAGS to avoid slice sampler getting stuck.
# Initial values are found as follows:
#     1. theta_guess[i, j] = Y[i, j]/n[i], or the is the MLE for theta
#     2. Z_guess[i, j] = 1 if the i is a sample, and the observed rc is >2x the expected rc assuming c[i] = 1
#     3. pi_guess[i] is the mean of row i
#     3. c_guess[i] is the estimated slope from regressing the observed rc to expected rc assuming c[i] = 1 for non-enriched peptides
#     4. phi_guess[i,j] = Y[i, j]/(n[i]*c_guess[i]*a_0[j]/(a_0[j] + b_0[j]))
# Note: we can imagine re-estimating Z_guess and c_guess given phi_guess to get closer estimates.
# ----- depends -----
# tidyverse
# stats
# ----- input -----
# data_list       list of simulation parameters fed to JAGS. This includes:
#                     N, P, B, n, Y, a_0, b_0, a_c,
#                     b_c, a_pi, b_pi, a_phi, b_phi, fc
# ----- output -----
# inits_list      list of initial values fed to JAGS. This includes:
#                     theta, c, pi, Z, phi
# -----------------------------
guessInits <- function(object, beads.prior){
    ## Extract convenient parameters
    N <- ncol(object)
    P <- nrow(object)
    B <- ifelse(object$group == getBeadsName(), 1, 0)
    n <- librarySize(object)
    Y <- counts(object)

    a_0 <- beads.prior[["a_0"]]
    b_0 <- beads.prior[["b_0"]]

    # Guess proportions, add offset of 1e-8 when theta is 0
    theta_guess <- Y/matrix(rep(n, P), nrow = P, byrow = TRUE)
    theta_guess <- theta_guess + 1e-8*(theta_guess == 0)

    ## Guess enriched peptides
    expected_prop <- a_0/(a_0 + b_0)
    expected_rc <- vapply(n, function(n_i) n_i*expected_prop,
                          numeric(length(expected_prop)))
    Z_guess <- (Y > expected_rc*2)*
        (1 - matrix(rep(B, P), nrow = P, byrow = TRUE))

    ## Guess proportions of enriched peptides
    pi_guess <- colMeans(Z_guess)

    ## Guess attenuation constant
    c_guess <- sapply(seq(N), function(index){
        if(B[index] != 1){
            coef(lm(Y[, index] ~ expected_rc[, index] - 1))
        } else 1
    })

    ## Guess fold change
    expected_rc_attn <- vapply(c_guess*n, function(n_i) n_i*expected_prop,
                               numeric(length(expected_prop)))
    phi_guess <- (1-Z_guess) + Y/expected_rc_attn*Z_guess

    list(theta = theta_guess,
         c_sample = vapply(c_guess, min, numeric(1), 1-1e-6), # ensure c < 1
         pi_sample = vapply(pi_guess, max, numeric(1), 1e-6), # ensure pi > 0
         Z = Z_guess,
         phi_sample = phi_guess)
}

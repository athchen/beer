#' Simulate a small data set for testing purposes.
#'
#' Beads-only prior parameters are estimated based on EC data provided below in
#' `hiv_url`.

hiv_url <- "https://github.com/athchen/beer_manuscript/blob/master/data_raw/hiv.rds"
hiv <- readRDS(url(hiv_url, "rb"))

# Simulate data -----------
set.seed(20210223)

# Constants
N <- 10
P <- 10
B <- c(rep(1, 4), rep(0, 6))
n <- round(abs(rnorm(N, 1.5e6, 1e5)))

# beads_only parameters, arrange peptides by proportion
beads_params <- beer::getAB(hiv, lower = 1) |>
    dplyr::as_tibble() |>
    dplyr::mutate(mean = a_0/(a_0 + b_0)) |>
    dplyr::arrange(mean)

# Randomly select peptides
pep_ind <- sort(sample(1:nrow(hiv), P))
a_0 <- beads_params$a_0[pep_ind]
b_0 <- beads_params$b_0[pep_ind]

# proportion of peptides enriched
Z <- cbind(matrix(0, nrow = P, ncol = sum(B)),
           sapply(1:(N - sum(B)), function(x) sample(c(1, rep(0, P-1)))))
pi <- colMeans(Z)

# fold-changes
phi <- apply(Z, 2, function(x) {
    phi_j <- rep(1, P)
    phi_cat <- runif(1, 2, 32)
    phi_j[which(x!= 0)] <- phi_cat
    phi_j
})

# Generate remaining data
a <- matrix(NA, nrow = P, ncol = N)
b <- matrix(NA, nrow = P, ncol = N)
theta <- matrix(NA, nrow = P, ncol = N)
Y <- matrix(NA, nrow = P, ncol = N)

for(j in 1:N){
    for(i in 1:P) {
        if (Z[i, j] == 1 & B[j] == 0) {
            mean_e <- min(phi[i, j]*a_0[i]/(a_0[i] + b_0[i]), 1)
            var_e <- a_0[i]*b_0[i]/((a_0[i] + b_0[i])^2*(a_0[i] + b_0[i] + 1))

            a[i, j] <- max((1-mean_e)*mean_e^2/var_e - mean_e, 1)
            b[i, j] <- a[i, j]*(1/mean_e - 1)

        } else if (Z[i, j] == 0 & B[j] == 1) {

            a[i, j] <- a_0[i]
            b[i, j] <- b_0[i]

        } else {
            mean_ne <- a_0[i]/(a_0[i] + b_0[i])
            var_ne <- (a_0[i]*b_0[i])/((a_0[i]+b_0[i])^2*(a_0[i] + b_0[i] + 1))

            a[i, j] <- max((1-mean_ne)*mean_ne^2/var_ne - var_ne, 1)
            b[i, j] <- a[i, j]*(1/mean_ne - 1)

        }

        theta[i, j] <- rbeta(1, a[i, j], b[i, j])
        if(is.na(theta[i, j])){
            print(paste0("Sample ", j, " and peptide ", i, " failed to create theta."))
        }

        Y[i, j] <- rbinom(1, n[j], theta[i, j])
    }
}

# Define new n
n_init <- n
n <- colSums(Y)

# Define c
c <- sapply(1:N, function(x){
    if(x %in% (sum(B) + 1):N) {
        ne_index <- which(Z[, x] == 0)
        data_sub <- data.frame(Y = Y[ne_index, x],
                               expected_rc = a_0[ne_index]/(a_0[ne_index] +
                                                                b_0[ne_index])*n[x])
        lm_fit <- lm(Y ~ expected_rc - 1, data = data_sub)
        coef(lm_fit)[[1]]
    } else 1
})

# Define phip_obj ------------
sample_params <- data.frame(group = c(rep("beads", sum(B)),
                                      rep("sample", N - sum(B))),
                            n_init = n_init,
                            n = n, true_c = c, true_pi = pi)
pep_params <- data.frame(a_0 = a_0, b_0 = b_0)
prior_params <- list(seed = 20210223,
                     a_pi = 1000, b_pi = 2.4e4,
                     a_phi = 1.25, b_phi = 0.1,
                     a_c = 33.02, b_c = 17.92,
                     fc = 1)
sim_data <- PhIPData(counts = Y,
                     peptideInfo = pep_params,
                     sampleInfo = sample_params,
                     metadata = prior_params)
assays(sim_data)[c("true_Z", "true_phi", "true_a", "true_b", "true_theta")] <-
    list(true_Z = Z, true_phi = phi, true_a = a, true_b = b, true_theta = theta)

saveRDS(sim_data, "inst/extdata/sim_data.rds")

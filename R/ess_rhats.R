# Effective sample size (ESS) and Gelman–Rubin R-hat

# Requires: sns, asbio
# Assumes an rstan fit object `fit` is in the workspace, with:
#   chains = 4, iter = 3000, warmup = 1000  (adjust below if different)

library(sns)
library(asbio)

iterations <- 3000   # total iterations per chain (including warmup)
warmup     <- 1000
chains     <- 4
post_start <- warmup + 1

c1 <- fit@sim[["samples"]][[1]]
c2 <- fit@sim[["samples"]][[2]]
c3 <- fit@sim[["samples"]][[3]]
c4 <- fit@sim[["samples"]][[4]]

# Parameter names as stored in the Stan samples list
par_list <- c(
  "beta[1]", "beta[2]", "beta[3]", "beta[4]",
  "beta[5]", "beta[6]", "beta[7]",
  "sigmasq", "phi"
)

# Friendly labels for printing
par_labels <- c(
  "beta[1]" = "beta_0 (ASC)",
  "beta[2]" = "beta_t (time)",
  "beta[3]" = "beta_c (cost)",
  "beta[4]" = "gamma_Vehicle",
  "beta[5]" = "gamma_Income",
  "beta[6]" = "gamma_Manhattan",
  "beta[7]" = "gamma_Male",
  "sigmasq" = "sigmasq",
  "phi"     = "phi"
)

ess_rhat <- data.frame(
  parameter = character(0),
  ESS = numeric(0),
  R_hat = numeric(0),
  stringsAsFactors = FALSE
)

for (par in par_list) {
  post1 <- c1[[par]][post_start:iterations]
  post2 <- c2[[par]][post_start:iterations]
  post3 <- c3[[par]][post_start:iterations]
  post4 <- c4[[par]][post_start:iterations]

  ess_total <- ess(post1, method = c("coda", "ise")) +
               ess(post2, method = c("coda", "ise")) +
               ess(post3, method = c("coda", "ise")) +
               ess(post4, method = c("coda", "ise"))

  rhat <- R.hat(cbind(post1, post2, post3, post4))

  label <- if (!is.null(par_labels[[par]])) par_labels[[par]] else par
  cat(label, "\n")
  cat("  ESS  :", ess_total, "\n")
  cat("  R_hat:", rhat, "\n\n")

  ess_rhat <- rbind(
    ess_rhat,
    data.frame(parameter = label, ESS = ess_total, R_hat = rhat)
  )
}

print(ess_rhat)

# Optional: write diagnostics to disk
# write.csv(ess_rhat, "results/ess_rhat.csv", row.names = FALSE)

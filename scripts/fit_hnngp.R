# End-to-end binary HNNGP estimation with rstan
# Mode choice: transit (1) vs driving (0), NYC travel-diary style data
#
# Usage (from the hnngp-tutorial/ directory):
#   source("scripts/fit_hnngp.R")

library(rstan)
library(HDInterval)

source("R/NNMatrix.R")

# ---------------------------------------------------------------------------
# 1. Data
# ---------------------------------------------------------------------------
my_data <- read.csv("data/mode_binary.csv", header = TRUE, sep = ",",
                    dec = ".", check.names = FALSE)
my_data <- my_data[complete.cases(my_data), ]

# Binary outcome: transit (PMODE_R3 == 2) vs driving
y    <- as.integer(my_data$PMODE_R3 == 2)
male <- as.integer(my_data$GENDER == "Male")
t    <- (my_data$T_Transit - my_data$T_Driving) / 60   # hours
c    <- my_data$Cost_Transit - my_data$Cost_Driving
veh  <- as.integer(my_data$HHVEH_R == 0)
hinc <- as.integer(my_data$INCOM_R == 4)
manh <- as.integer(my_data$GEO_GROUP1_D == 1)

# Spatial coordinates: trip direction (destination - origin)
s_o <- my_data[, c("INTPTLAT", "INTPTLONG")]
s_d <- my_data[, c("DLAT", "DLONG")]
s   <- as.matrix(s_d - s_o)

X <- cbind(1, t, c, veh, hinc, manh, male)
N <- nrow(X)

# Optional: subsample for a faster demo
#set.seed(42); idx <- sample(N, 500); y <- y[idx]; X <- X[idx, ]; s <- s[idx, ]; N <- nrow(X)

# Train / test mask (1 = training)
train_prop <- 0.8
set.seed(42)
mask <- sample(c(1L, 0L), N, replace = TRUE, prob = c(train_prop, 1 - train_prop))

# ---------------------------------------------------------------------------
# 2. Nearest-neighbor matrix (NNGP)
# ---------------------------------------------------------------------------
M <- 20
NN.matrix <- NNMatrix(
  coords = s,
  n.neighbors = M,
  n.omp.threads = 2,
  ord = order(s[, 1])
)

# ---------------------------------------------------------------------------
# 3. Priors
# ---------------------------------------------------------------------------
P  <- ncol(X) - 1
uB <- rep(0, P + 1)
uB[4] <- 3.42                 # informative prior mean for no-vehicle indicator
VB <- diag(P + 1)
ss <- 0.5                     # sigma ~ Normal(1, ss)  [half via <lower=0>]
st <- 0.001                   # unused when tausq is fixed in the Stan file
ap <- 3; bp <- 0.5            # phi ~ Gamma(ap, bp)

# ---------------------------------------------------------------------------
# 4. Stan data (observations ordered by NN ordering)
# ---------------------------------------------------------------------------
stan_data <- list(
  N = N, M = M, P = P,
  Y = y[NN.matrix$ord],
  X = X[NN.matrix$ord, , drop = FALSE],
  NN_ind = NN.matrix$NN_ind,
  NN_dist = NN.matrix$NN_dist,
  NN_distM = NN.matrix$NN_distM,
  mask = mask[NN.matrix$ord],
  uB = uB, VB = VB, ss = ss, st = st, ap = ap, bp = bp
)

# ---------------------------------------------------------------------------
# 5. MCMC
# ---------------------------------------------------------------------------
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

parameters <- c("beta", "sigmasq", "phi", "log_lik", "y_est")

fit <- stan(
  file = "stan/nngp_binary.stan",
  data = stan_data,
  pars = parameters,
  chains = 4,
  iter = 3000,
  warmup = 1000,
  refresh = 50,
  verbose = FALSE
)

dir.create("results", showWarnings = FALSE)
saveRDS(fit, "results/hnngp_binary_fit.rds")
# save.image("results/hnngp_binary.RData")

print(fit, pars = c("beta", "sigmasq", "phi"))

# ---------------------------------------------------------------------------
# 6. Posterior summaries of interest
# ---------------------------------------------------------------------------
post <- rstan::extract(fit, pars = c("beta", "sigmasq", "phi"))

# Value of time (USD/hr): (beta_t / beta_c) * 60
# beta[,2] = time, beta[,3] = cost
vot <- (post$beta[, 2] / post$beta[, 3]) * 60
cat("\nVOT posterior mean:", mean(vot),
    " median:", median(vot), "\n")
cat("VOT 95% HDI:", paste(round(hdi(vot, credMass = 0.95), 2), collapse = " to "), "\n")

# ---------------------------------------------------------------------------
# 7. ESS and R-hat (see R/ess_rhats.R)
# ---------------------------------------------------------------------------
source("R/ess_rhats.R")

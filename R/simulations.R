# simulations for sequestration potential paper
library(withr)
library(dplyr)


############ core DGP function ###########
dgp_soc <- function(
    N,
    tau = 0.2,
    scenario = c("none", "saturation", "complex"),
    beta = 0.1,
    beta1 = -0.1,
    beta2 = 0
) {
  scenario <- match.arg(scenario)
  
  # ---- Baseline %SOC ----
  baseline <- pmax(rnorm(N, mean = 3, sd = 1), 0)
  
  # Control potential outcome
  y0 <- baseline
  
  # Normalize baseline to [-1,1]
  B <- 2 * (baseline - min(baseline)) /
    (max(baseline) - min(baseline)) - 1
  
  # Additional moderators
  S <- sample(rep(c(0, 1), each = N/2))    # soil type
  C <- runif(N, -1, 1)          # NDVI proxy
  
  # ---- Treatment potential outcome ----
  if (scenario == "none") {
    y1 <- y0 + tau
    
  } else if (scenario == "saturation") {
    y1 <- y0 + tau + beta * B
    
  } else if (scenario == "complex") {
    nonlinear_term <- exp(B * C)
    nonlinear_term <- nonlinear_term / max(nonlinear_term)
    
    y1 <- y0 +
      tau +
      beta1 * B * S +
      beta2 * nonlinear_term
  }
  
  # Randomized assignment
  Z <- rep(c(0, 1), each = N/2)
  Z <- sample(Z)
  
  Y <- ifelse(Z == 1, y1, y0)
  
  X <- cbind(
    baseline = baseline,
    soil = S,
    ndvi = C
  )
  
  list(
    X = X,
    Z = Z,
    Y = Y,
    y0 = y0,
    y1 = y1,
    tau_true = y1 - y0
  )
}

########## simulation runner function #######

run_soc_simulation <- function(N, tau, scenario_params, num_sims, seed) {
  
  results <- vector("list", num_sims)
  
  for (m in seq_len(num_sims)) {
    
    # remove metadata field before calling dgp_soc
    dgp_args <- scenario_params[names(scenario_params) != "label"]
    
    dat <- with_seed(seed + m, {
      do.call(
        dgp_soc,
        c(
          list(N = N, tau = tau),
          dgp_args
        )
      )
    })
    
    y_pop <- cbind(dat$y1, dat$y0)
    
    dim_res <- estimate_optimal_policy(
      Y = dat$Y,
      Z = dat$Z,
      X = dat$X,
      learner = learner_dim,
      y_pop = y_pop
    )
    
    ols_res <- estimate_optimal_policy(
      Y = dat$Y,
      Z = dat$Z,
      X = dat$X,
      learner = learner_ols,
      y_pop = y_pop
    )
    
    rf_res <- estimate_optimal_policy(
      Y = dat$Y,
      Z = dat$Z,
      X = dat$X,
      learner = learner_rf,
      y_pop = y_pop
    )
    
    K <- floor(0.10 * N)
    
    ols_10 <- estimate_optimal_policy(
      Y = dat$Y,
      Z = dat$Z,
      X = dat$X,
      learner = learner_ols,
      y_pop = y_pop,
      budget = K
    )
    
    rf_10 <- estimate_optimal_policy(
      Y = dat$Y,
      Z = dat$Z,
      X = dat$X,
      learner = learner_rf,
      y_pop = y_pop,
      budget = K
    )
    
    results[[m]] <- tibble(
      sim = m,
      N = N,
      tau = tau,
      scenario = scenario_params$label,
      blanket = dim_res$policy_value,
      ols = ols_res$policy_value,
      rf = rf_res$policy_value,
      ols_10 = ols_10$policy_value,
      rf_10 = rf_10$policy_value,
      oracle = ols_res$oracle_value
    )
  }
  
  bind_rows(results)
}








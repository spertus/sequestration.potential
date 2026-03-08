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
  
  # ---- Normalized Baseline %SOC ----
  baseline <- pmax(rnorm(N, mean = 3, sd = 1), 0)
  
  # Control potential outcome
  y0 <- baseline
  
  # Normalize baseline to standard gaussian
  B <- (baseline - mean(baseline)) / sd(baseline)
  
  # Additional moderators
  S <- sample(rep(c(0, 1), each = N/2))    # soil type proxy
  C <- runif(N, -1, 1)          # NDVI proxy
  
  # ---- Treatment potential outcome ----
  if (scenario == "none") {
    y1 <- y0 + tau
    X <- cbind(
      baseline = B
    )
  } else if (scenario == "saturation") {
    y1 <- y0 + tau + beta * B
    X <- cbind(
      baseline = B
    )
  } else if (scenario == "complex") {
    
    y1 <- y0 +
      tau +
      beta1 * B * S +
      beta2 * B * C
    X <- cbind(
      baseline = B,
      soil = S,
      ndvi = C
    )
  }
  
  # Randomized assignment
  Z <- rep(c(0, 1), each = N/2)
  Z <- sample(Z)
  
  Y <- ifelse(Z == 1, y1, y0)
  
  
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

run_soc_simulation <- function(
    N,
    tau,
    scenario_params,
    num_sims,
    seed,
    N_eval = 500
) {
  
  results <- vector("list", num_sims)
  
  # remove metadata field before calling dgp_soc
  dgp_args <- scenario_params[names(scenario_params) != "label"]
  
  for (m in seq_len(num_sims)) {
    
    # ----------------------------
    # Training study draw
    # ----------------------------
    dat_train <- with_seed(seed + m, {
      do.call(
        dgp_soc,
        c(
          list(N = N, tau = tau),
          dgp_args
        )
      )
    })
    
    # ----------------------------
    # Held-out evaluation population
    # ----------------------------
    dat_eval <- with_seed(seed + 100000 + m, {
      do.call(
        dgp_soc,
        c(
          list(N = N_eval, tau = tau),
          dgp_args
        )
      )
    })
    
    x_pop <- dat_eval$X
    y_pop <- cbind(dat_eval$y1, dat_eval$y0)
    
    # Blanket policy evaluated on held-out population
    dim_res <- estimate_optimal_policy(
      Y = dat_train$Y,
      Z = dat_train$Z,
      X = dat_train$X,
      learner = learner_dim,
      x_pop = x_pop,
      y_pop = y_pop
    )
    
    # OLS optimal policy evaluated on held-out population
    ols_res <- estimate_optimal_policy(
      Y = dat_train$Y,
      Z = dat_train$Z,
      X = dat_train$X,
      learner = learner_ols,
      x_pop = x_pop,
      y_pop = y_pop
    )
    
    # RF optimal policy evaluated on held-out population
    rf_res <- estimate_optimal_policy(
      Y = dat_train$Y,
      Z = dat_train$Z,
      X = dat_train$X,
      learner = learner_rf,
      x_pop = x_pop,
      y_pop = y_pop
    )
    
    # Treat-best-10% on held-out population
    K_eval <- floor(0.10 * N_eval)
    
    ols_10 <- estimate_optimal_policy(
      Y = dat_train$Y,
      Z = dat_train$Z,
      X = dat_train$X,
      learner = learner_ols,
      x_pop = x_pop,
      y_pop = y_pop,
      budget = K_eval
    )
    
    rf_10 <- estimate_optimal_policy(
      Y = dat_train$Y,
      Z = dat_train$Z,
      X = dat_train$X,
      learner = learner_rf,
      x_pop = x_pop,
      y_pop = y_pop,
      budget = K_eval
    )
    
    results[[m]] <- tibble(
      sim = m,
      N = N,
      N_eval = N_eval,
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








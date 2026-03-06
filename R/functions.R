library(withr)
library(testthat)
library(tidyverse)
library(grf)




# helper functions
`%||%` <- function(x, y) if (is.null(x)) y else x


# main functions
estimate_cate <- function(
    Y,
    Z,
    X,
    x_pop = NULL,
    learner = learner_ols,
    ci_level = 0.95,
    ...
) {
  # takes data from a size n RCT with uniform assignment and returns an OLS estimate of the science table either for that RCT, or for a size N population with covariates x_pop
  # inputs:
    # Y = length-n vector; outcomes in the RCT 
    # Z = length-n binary vector; treatment assignment generically, 0 denotes control and 1 denotes treatment
    # X = n-by-p matrix; covariates in RCT
    # x_pop = N-by-p matrix; covariates in population; may be left out, in which case the science table for the experiment is returned instead. must be in exact same order of covariates as X 
    # learner = function; 
    # any additional arguments to the learner function (must be named)
  # outputs: 
  # a n-by-2 (default) science table for the experiment or N-by-2 (if x_pop) is given science table for the population
  
  # ensure correct dimensions
  stopifnot(
    nrow(X) == length(Y),
    length(Z) == length(Y)
  )
  
  if (is.null(x_pop)) x_pop <- X
  
  # run the CATE learner on the data
  res <- learner(
    Y = Y,
    Z = Z,
    X = X,
    x_pop = x_pop,
    ...
  )
  
  if (!is.list(res) || !("cate" %in% names(res))) {
    stop("learner must return a list with element 'cate'")
  }
  
  cate <- res$cate # retrieve CATE estimate
  se   <- res$se %||% rep(NA, length(cate)) # retrieve SE if produced, otherwise return NA
  
  # normal theory confidence intervals
  alpha <- 1 - ci_level
  z <- qnorm(1 - alpha / 2)
  
  # return dataframe with estimates, standard errors, CIs for each population unit
  data.frame(
    cate  = cate,
    se    = se,
    lower = cate - z * se,
    upper = cate + z * se
  )
}

learner_ols <- function(Y, Z, X, x_pop, ...) {
  # function to learn CATE using OLS. 
  # A detailed discussion of this method appears in https://www.tandfonline.com/doi/full/10.1080/01621459.2017.1407322
  
  # check that the rank conditions are satisfied
  p <- ncol(X)
  if (sum(Z == 0) <= p || sum(Z == 1) <= p) {
    stop("Not enough observations in one treatment arm to fit OLS models.")
  }
  
  # ensure matrices
  X <- as.matrix(X)
  x_pop <- if (is.null(x_pop)) X else as.matrix(x_pop)
  if(ncol(X) != ncol(x_pop)) stop("Number of covariates (columns) in x_pop does not match number in X")
  
  # add intercept explicitly
  X_ctl <- cbind(1, X[Z == 0, , drop = FALSE])
  X_trt <- cbind(1, X[Z == 1, , drop = FALSE])
  
  Y_ctl <- Y[Z == 0]
  Y_trt <- Y[Z == 1]
  
  # OLS fits
  beta_ctl <- solve(t(X_ctl) %*% X_ctl, t(X_ctl) %*% Y_ctl)
  beta_trt <- solve(t(X_trt) %*% X_trt, t(X_trt) %*% Y_trt)
  
  # prediction matrix
  Xp <- cbind(1, x_pop)
  
  mu_0 <- as.vector(Xp %*% beta_ctl)
  mu_1 <- as.vector(Xp %*% beta_trt)
  
  cate <- mu_1 - mu_0
  
  # conservative SE from squared residuals (correct under homoskedasticity)
  sigma2_ctl <- mean((Y_ctl - X_ctl %*% beta_ctl)^2)
  sigma2_trt <- mean((Y_trt - X_trt %*% beta_trt)^2)
  
  V_ctl <- sigma2_ctl * solve(t(X_ctl) %*% X_ctl)
  V_trt <- sigma2_trt * solve(t(X_trt) %*% X_trt)
  
  se <- sqrt(
    rowSums((Xp %*% V_trt) * Xp) +
      rowSums((Xp %*% V_ctl) * Xp)
  )
  
  list(cate = cate, se = se)
}



learner_rf <- function(Y, Z, X, x_pop, ...) {
  # function to learn CATE using 'causal forest' as described in https://grf-labs.github.io/grf/reference/causal_forest.html
  # see also https://projecteuclid.org/journals/annals-of-statistics/volume-47/issue-2/Generalized-random-forests/10.1214/18-AOS1709.full
  # and https://www.tandfonline.com/doi/full/10.1080/01621459.2017.1319839
  if (!requireNamespace("grf", quietly = TRUE)) {
    stop("Package 'grf' is required to run learner_rf")
  }
  # check that the rank conditions are satisfied
  # ensure matrices
  X <- as.matrix(X)
  x_pop <- if (is.null(x_pop)) X else as.matrix(x_pop)
  if(ncol(X) != ncol(x_pop)) stop("Number of covariates (columns) in x_pop does not match number in X")
  
  p <- ncol(X)
  if (sum(Z == 0) < p || sum(Z == 1) < p) {
    stop("Not enough observations in one treatment arm to fit OLS models.")
  }
  
  prop_score <- mean(Z) # propensity scores for RCT, for use with causal_forest
  # conversion of the two covariate matrices to dataframes with common structure
  X_df <- as.data.frame(X)
  colnames(X_df) <- paste0("X", seq_len(ncol(X_df)))
  if (is.null(x_pop)) {
    x_pop_df <- X_df
  } else {
    x_pop_df <- as.data.frame(x_pop)
    colnames(x_pop_df) <- colnames(X_df)
  }
  
  # causal forest is fit using default parameters from the 'grf' package
  cf <- grf::causal_forest(
    X = X,
    Y = Y,
    W = Z,
    W.hat = prop_score,
    ...
  )
  
  cate_preds <- predict(cf, newdata = x_pop, estimate.variance = TRUE)
  
  list(
    cate = as.numeric(cate_preds$predictions),
    se   = sqrt(cate_preds$variance.estimates)
  )
}

learner_dim <- function(Y, Z, X, x_pop, ...) {
  # function that estimates the CATE using the difference-in-means, that is, ignoring the covariates and using the standard ATE estimate
  
  # this just records the size of the population or the study
  if(is.null(x_pop)){
    N <- nrow(X)
  } else{
    N <- nrow(x_pop)
  }
  
  mu_0 <- mean(Y[Z == 0])
  mu_1 <- mean(Y[Z == 1])
  
  cate <- rep(mu_1 - mu_0, N)
  se   <- rep(sqrt(var(Y[Z == 1]) / sum(Z == 1) + var(Y[Z == 0]) / sum(Z == 0)), N)
  
  list(cate = cate, se = se)
}


solve_lp_policy <- function(
    cate,
    budget = NULL,
    costs = NULL
) {
  # function to solve optimal policy problem with (estimated) CATEs as a general linear program with budget constraints
  # inputs:
    # cate: length-N vector, conditional average treatment effects for the population
    # budget: scalar, the overall budget for the problem; defaults to NULL, i.e. an unconstrained portfolio; to treat the K-best, set budget = K and do not provide a budget 
    # costs: length-N vector, additive costs of treatment for each unit of the population
  # outputs:
    # an (estimated) optimal policy within the portfolio
  if (!requireNamespace("lpSolve", quietly = TRUE)) {
    stop("Package 'lpSolve' is required.")
  }
  
  n <- length(cate)
  
  # Unconstrained case: treat if cate > 0
  if (is.null(budget) && is.null(costs)) {
    return(as.integer(cate > 0))
  }
  
  # Defaults
  if (is.null(costs)) {
    costs <- rep(1, n)
  }
  
  stopifnot(length(costs) == n)
  
  # Objective
  obj <- cate
  
  # Budget constraint
  A <- matrix(costs, nrow = 1)
  b <- budget
  
  res <- lpSolve::lp(
    direction = "max",
    objective.in = obj,
    const.mat = A,
    const.dir = "<=",
    const.rhs = b,
    all.bin = TRUE
  )
  
  as.integer(res$solution > 0.5)
}


estimate_optimal_policy <- function(
    Y,
    Z,
    X,
    learner = learner_ols,
    x_pop = NULL,
    y_pop = NULL,
    budget = NULL,
    costs = NULL
) {
  # takes data from a size n RCT with uniform assignment, estimates CATES, and returns an estimate of a "first-best" (unconstrained) policy, either for the study subjects or for a population with covariates x_pop
  # the "first-best" policy over an unconstrained portfolio of binary treatments just assigns units with positive CATEs to treatment and negative CATEs to control
  # inputs:
    # Y = length-n vector; outcomes in the RCT 
    # Z = length-n binary vector; treatment assignment generically, 0 denotes control and 1 denotes treatment
    # X = n-by-p matrix; covariates in RCT
    # x_pop = N-by-p matrix; covariates in population; may be left out, in which case the science table for the experiment is returned instead
    # y_pop = N-by-2 matrix; optional true science table for the population as a matrix of potential outcomes on treatment (1st column) and control (2nd column); if ommitted, only the policy estimate is returned without the acutal return
    # budget = scalar; the budget for the portfolio; defaults to NULL, which is an unconstrained portfolio; for 'treat the K-best' set this to K and leave the costs as NULL
    # costs = length-N vector of positive reals; the cost of treating each unit in an additive cost model
  # outputs: 
    # list with elements: 
    #required 
      # policy = an estimate of the "first-best" (unconstrained) optimal policy  
      # cate = the conditional average treatment effect for every unit in the population (if x_pop is provided) or study   
    #optional (if y_pop, the science table for the population, is provided)
      # policy_value = the realized value of the estimated optimal policy
      # oracle_value = the realized value of the true optimal policy
      # regret = the realized regret of using the estimated optimal policy
  cate_res <- estimate_cate(
    Y = Y,
    Z = Z,
    X = X,
    learner = learner,
    x_pop = x_pop
  )
  
  # extract the CATE point estimates
  cate <- cate_res$cate
  
  policy <- solve_lp_policy(
    cate = cate,
    budget = budget,
    costs = costs
  )
  
  out <- list(
    policy = policy,
    cate = cate
  )
  
  if (!is.null(y_pop)) {
    stopifnot(
      is.matrix(y_pop),
      nrow(y_pop) == length(policy),
      ncol(y_pop) == 2
    )
    
    policy_value <- mean(
      ifelse(policy == 1, y_pop[,1], y_pop[,2])
    )
    
    oracle_value <- mean(pmax(y_pop[,1], y_pop[,2]))
    
    out$policy_value <- policy_value
    out$oracle_value <- oracle_value
    out$regret       <- oracle_value - policy_value
  }
  
  out
}

########### generic simulation functions ###########

run_policy_simulation <- function(
    dgp,
    n,
    learner,
    num_sims = 100,
    budget = NULL,
    costs = NULL,
    seed = NULL
) {
  # code to run simulate regret of optimal policy estimation
  # inputs:
    # dgp = function, the data generating process
    # n = integer, the size of the study
    # learner = function, the base learner for estimating response functions
    # num_sims = integer, the number of simulation iterations
    # budget = scalar, the budget for a constrained policy 
    # costs = length-n vector of positive reals, the additive costs for each unit in the DGP
    # seed = integer, an optional seed for the pseudo-random numbers
  # outputs: 
    # a num_sims-by-4 data.frame with a row for each simulation run and columns for estimates and oracle values 
  if (!is.null(seed)) set.seed(seed)
  
  results <- vector("list", num_sims)
  
  for (m in seq_len(num_sims)) {
    # this locks the seed for each DGP, which allows DGPs to be replicated when the seed is set, even if learners use auxiliary randomness that advance the seed
    if(!is.null(seed)){
      seed_each <- seed + m
      withr::with_seed(seed_each, {
        dat <- dgp(n)
      })
    } else{
      dat <- dgp(n)
    }
    
    
    res <- estimate_optimal_policy(
      Y = dat$Y,
      Z = dat$Z,
      X = dat$X,
      learner = learner,
      y_pop = cbind(dat$y_1, dat$y_0),
      budget = budget,
      costs = costs
    )
    
    results[[m]] <- data.frame(
      policy_value = res$policy_value,
      oracle_value = res$oracle_value,
      regret       = res$regret,
      treat_rate   = mean(res$policy)
    )
  }
  
  do.call(rbind, results)
}

dgp_zero_tau <- function(n) {
  # trivial dgp with no treatment effect whatsoever (Fisher's sharp null is satisfied)
  X <- matrix(rnorm(n), n, 1)
  
  #Z <- rbinom(n, 1, 0.5) # bernoulli experiment
  Z <- rep(0, n)
  Z[sample(1:n, ceiling(n/2))] <- 1 # completely randomized experiment
  
  y_0 <- rnorm(n)
  y_1 <- y_0
  
  Y <- y_0
  
  list(
    X = X,
    Z = Z,
    Y = Y,
    y_0 = y_0,
    y_1 = y_1,
    tau = rep(0, n)
  )
}


dgp_constant_tau <- function(n, tau = 1) {
  # an example of a data generating process with constant treatmnent effect that defaults to 1
  X <- matrix(rnorm(n), n, 1)
  
  #Z <- rbinom(n, 1, 0.5) # bernoulli experiment
  Z <- rep(0, n)
  Z[sample(1:n, ceiling(n/2))] <- 1 # completely randomized experiment
  
  y_0 <- rnorm(n)
  y_1 <- y_0 + tau
  
  Y <- ifelse(Z == 1, y_1, y_0)
  
  list(
    X = X,
    Z = Z,
    Y = Y,
    y_0 = y_0,
    y_1 = y_1,
    tau = rep(tau, n)
  )
}



dgp_saturation <- function(n, saturation_effect = 1){
  # data generating process function under the saturation hypothesis
  # input: 
    # n = the size of the study 
  # output:
    #  list with covariates, treatments, observed outcomes, potential outcomes, and individual treatment effects
  
  X <- matrix(rnorm(n * 2), n, 2) # proxy for baseline SOC
  
  tau <- 1 - saturation_effect * X[,1] # heterogeneous treatment effect
  mu_0  <- X[,2]
  
  y_0 <- mu_0 + rnorm(n, sd = sigma)
  y_1 <- y_0 + tau
  
  #Z <- rbinom(n, 1, 0.5) # bernoulli experiment
  Z <- rep(0, n)
  Z[sample(1:n, ceiling(n/2))] <- 1 # completely randomized experiment
  Y <- ifelse(Z == 1, y_1, y_0)
  
  list(
    X = X,
    Z = Z,
    Y = Y,
    y_0 = y_0,
    y_1 = y_1,
    tau = tau
  )
}



dgp_linear_heterogeneous <- function(n, p = 3, sigma = 1){
  # an example of a data generating process function
  # input: 
    # n = the size of the study 
    # p = the number of covariates, defaults to 3
    # sigma = the 'noise' in the potential outcomes
  # output:
    #  list with covariates, treatments, observed outcomes, potential outcomes, and individual treatment effects
  X <- matrix(rnorm(n * p), n, p)
  
  tau <- 1 + X[,1] # heterogeneous treatment effect
  mu_0  <- X[,2]
  
  y_0 <- mu_0 + rnorm(n, sd = sigma)
  y_1 <- y_0 + tau
  
  #Z <- rbinom(n, 1, 0.5) # bernoulli experiment
  Z <- rep(0, n)
  Z[sample(1:n, ceiling(n/2))] <- 1 # completely randomized experiment
  Y <- ifelse(Z == 1, y_1, y_0)
  
  list(
    X = X,
    Z = Z,
    Y = Y,
    y_0 = y_0,
    y_1 = y_1,
    tau = tau
  )
}


dgp_nonlinear_interactions <- function(
    n,
    p = 5,
    sigma = 1,
    interaction_scale = 1,
    nonlinear_scale = 1
) {
  # a complex data generating process with non-linear covariate effects and interactions between covariates 
  stopifnot(p >= 3)
  
  # Covariates
  X <- matrix(rnorm(n * p), n, p)
  
  # Nonlinear baseline
  mu_0 <- nonlinear_scale * (
    sin(X[,1]) +
      0.5 * X[,2]^2 -
      0.3 * exp(-X[,3])
  )
  
  # CATE with nonlinearity and interactions
  tau <- interaction_scale * (
    1 +
      X[,1] * X[,2] -
      0.5 * X[,3]^2 +
      0.3 * X[,4] * X[,5]
  )
  
  # Potential outcomes
  y_0 <- mu_0 + rnorm(n, sd = sigma)
  y_1 <- y_0 + tau
  
  # Randomized treatment
  Z <- rep(0, n)
  Z[sample(1:n, ceiling(n/2))] <- 1 # completely randomized experiment
  
  Y <- ifelse(Z == 1, y_1, y_0)
  
  list(
    X = X,
    Z = Z,
    Y = Y,
    y_0 = y_0,
    y_1 = y_1,
    tau = tau,
    mu_0 = mu_0
  )
}






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
  mu_0 <- res$mu_0 %||% rep(NA_real_, length(cate)) # per-arm potential outcome estimates
  mu_1 <- res$mu_1 %||% rep(NA_real_, length(cate)) # (NA when the learner does not provide them)

  # normal theory confidence intervals
  alpha <- 1 - ci_level
  z <- qnorm(1 - alpha / 2)

  # return dataframe with estimates, standard errors, CIs, and (when available)
  # per-arm potential outcome predictions for each population unit
  data.frame(
    cate  = cate,
    se    = se,
    lower = cate - z * se,
    upper = cate + z * se,
    mu_0  = mu_0,
    mu_1  = mu_1
  )
}

learner_ols <- function(Y, Z, X, x_pop, ...) {
  # Estimates the CATE via the interaction-regression approach of
  # Ding, Feller, and Miratrix (2019, JASA;
  # https://www.tandfonline.com/doi/full/10.1080/01621459.2017.1407322).
  # Fits a single OLS of Y on (1, Z, X_centered, Z * X_centered), then forms
  # CATEs and HC2 (Eicker-Huber-White) standard errors at each row of x_pop.
  #
  # The CATE point estimates are mathematically identical to those from fitting
  # a separate OLS in each arm (the previous implementation), but the HC2 SEs
  # are valid under arbitrary heteroskedasticity and treatment-effect
  # heterogeneity (homoskedastic SE is not). HC2 also
  # approximates the Neyman-style design-based variance for CATE estimators.

  # ensure matrices
  X     <- as.matrix(X)
  x_pop <- if (is.null(x_pop)) X else as.matrix(x_pop)
  if (ncol(X) != ncol(x_pop)) {
    stop("Number of covariates (columns) in x_pop does not match number in X")
  }

  n <- length(Y)
  p <- ncol(X)
  if (sum(Z == 0) <= p || sum(Z == 1) <= p) {
    stop("Not enough observations in one treatment arm to fit OLS models.")
  }

  # Center X at the study mean. Centering does not change CATE estimates or
  # SEs at any x_pop, but it makes the coefficient on Z directly interpretable
  # as the estimated ATE at the study mean.
  X_center <- colMeans(X)
  X_c      <- sweep(X,2,X_center,FUN = "-")
  ZX_c     <- Z * X_c

  # Design matrix columns: intercept, Z, X_centered (p cols), Z * X_centered (p cols)
  D <- cbind(1, Z, X_c, ZX_c)

  # OLS coefficients
  XtX_inv <- solve(crossprod(D))
  beta    <- as.vector(XtX_inv %*% crossprod(D, Y))

  # Residuals and leverages for HC2
  e <- as.vector(Y - D %*% beta)
  h <- rowSums((D %*% XtX_inv) * D)             # leverage h_ii
  w <- e^2 / pmax(1 - h, .Machine$double.eps)   # HC2 weights

  # Sandwich variance: V_HC2 = (D'D)^{-1} (sum_i w_i D_i D_i') (D'D)^{-1}
  meat  <- crossprod(D, D * w)
  V_hc2 <- XtX_inv %*% meat %*% XtX_inv

  # CATE at each population point: mu(1, x) - mu(0, x) =
  #   beta_Z + beta_{Z:X}^T (x - X_center).
  # In the design parameterization, this is the contrast c = (0, 1, 0_p, x_centered).
  X_pop_c <- sweep(x_pop, 2, X_center, FUN = "-")
  N       <- nrow(X_pop_c)
  C       <- cbind(0, 1, matrix(0, N, p), X_pop_c)

  cate <- as.vector(C %*% beta)
  se   <- sqrt(rowSums((C %*% V_hc2) * C))

  # Per-arm potential outcome predictions at each population point.
  # mu_0(x) corresponds to the design contrast c0 = (1, 0, x_centered, 0)
  #   => mu_0(x) = beta_intercept + beta_X^T (x - X_center)
  # mu_1(x) corresponds to c1 = (1, 1, x_centered, x_centered)
  #   => mu_1(x) = (beta_intercept + beta_Z) + (beta_X + beta_{Z:X})^T (x - X_center)
  # mu_1 - mu_0 reproduces `cate` exactly (and equals the per-arm OLS fit).
  C0   <- cbind(1, 0, X_pop_c, matrix(0, N, p))
  C1   <- cbind(1, 1, X_pop_c, X_pop_c)
  mu_0 <- as.vector(C0 %*% beta)
  mu_1 <- as.vector(C1 %*% beta)

  list(
    cate = cate,
    se   = se,
    mu_0 = mu_0,
    mu_1 = mu_1,
    inference_method = "HC2 (Ding-Feller-Miratrix interaction regression)"
  )
}



learner_rf <- function(Y, Z, X, x_pop, ...) {
  # Learn CATE and per-arm potential outcomes via grf::causal_forest, with a
  # companion grf::regression_forest used to estimate the marginal regression
  # E[Y | X]. Per-arm potential outcomes are backed out coherently from the
  # causal forest's CATE and the regression forest's E[Y | X], using the known
  # (constant) propensity for an RCT:
  #     mu_1(x) = E[Y | X = x] + (1 - p) * cate(x)
  #     mu_0(x) = E[Y | X = x] -      p  * cate(x)
  # See https://grf-labs.github.io/grf/reference/causal_forest.html and
  # Athey & Wager (2019, Annals of Statistics; arXiv:1610.01271).
  if (!requireNamespace("grf", quietly = TRUE)) {
    stop("Package 'grf' is required to run learner_rf")
  }

  # ensure matrices
  X     <- as.matrix(X)
  x_pop <- if (is.null(x_pop)) X else as.matrix(x_pop)
  if (ncol(X) != ncol(x_pop)) stop("Number of covariates (columns) in x_pop does not match number in X")

  p <- ncol(X)
  if (sum(Z == 0) < p || sum(Z == 1) < p) {
    stop("Not enough observations in one treatment arm to fit OLS models.")
  }

  prop_score <- mean(Z) # known propensity for an RCT

  # Marginal regression E[Y | X]. We use the OOB predictions on training X as
  # Y.hat for causal_forest (avoiding the in-bag overfit), and in-sample
  # predictions at x_pop to compute mu_0 and mu_1. The two prediction calls use
  # the same forest, so mu_1 - mu_0 is consistent with `cate` from the causal
  # forest at training points and is a smooth extrapolation elsewhere.
  rf_y        <- grf::regression_forest(X = X, Y = Y, ...)
  Y_hat_train <- as.numeric(predict(rf_y)$predictions)                # OOB
  Y_hat_pop   <- as.numeric(predict(rf_y, newdata = x_pop)$predictions)

  # Causal forest. Passing Y.hat explicitly ensures the cate is residualised
  # against the same regression forest used below to compute mu_0 and mu_1.
  cf <- grf::causal_forest(
    X     = X,
    Y     = Y,
    W     = Z,
    Y.hat = Y_hat_train,
    W.hat = prop_score,
    ...
  )

  cate_preds <- predict(cf, newdata = x_pop, estimate.variance = TRUE)
  cate       <- as.numeric(cate_preds$predictions)

  mu_1 <- Y_hat_pop + (1 - prop_score) * cate
  mu_0 <- Y_hat_pop -      prop_score  * cate

  list(
    cate = cate,
    se   = sqrt(cate_preds$variance.estimates),
    mu_0 = mu_0,
    mu_1 = mu_1
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

  mu_0_hat <- mean(Y[Z == 0])
  mu_1_hat <- mean(Y[Z == 1])

  cate <- rep(mu_1_hat - mu_0_hat, N)
  se   <- rep(sqrt(var(Y[Z == 1]) / sum(Z == 1) + var(Y[Z == 0]) / sum(Z == 0)), N)

  list(
    cate = cate,
    se   = se,
    mu_0 = rep(mu_0_hat, N),
    mu_1 = rep(mu_1_hat, N)
  )
}


akm_winner_ci <- function(estimates, cov_matrix, ci_level = 0.95) {
  # Lee-Sun-Sun-Taylor / Andrews-Kitagawa-McCloskey conditional confidence
  # interval for the value V_{k_hat} of the "winning" candidate among M
  # candidates, where k_hat = argmax(estimates) and the candidates are
  # jointly approximately normal: estimates ~ N(true_values, cov_matrix).
  #
  # The CI corrects the "winner's curse": E[max_k V_hat_k] is upward-biased
  # for max_k V_k, and a naive Wald CI around the winner under-covers V_{k_hat}.
  # The AKM/Lee CI is built by inverting tests for V_{k_hat} based on the
  # conditional distribution of V_hat_{k_hat} given the polyhedral selection
  # event {V_hat_{k_hat} >= V_hat_j for all j != k_hat}; this conditional
  # distribution is a truncated normal (Lee et al. 2016, Annals of Statistics
  # 44(3); Andrews, Kitagawa, McCloskey 2024, QJE 139(1)).
  #
  # Inputs:
  #   estimates  - length-M numeric vector of point estimates V_hat_1,...,V_hat_M
  #   cov_matrix - M x M (asymptotic / bootstrap) covariance matrix of estimates
  #   ci_level   - nominal coverage level (default 0.95)
  # Output: list with
  #   estimate     - V_hat_{k_hat}
  #   ci           - 2-vector (lower, upper) AKM conditional CI for V_{k_hat}
  #   sigma        - sd(V_hat_{k_hat}) from cov_matrix
  #   truncation   - 2-vector (V_minus, V_plus) implied by the polyhedral
  #                  selection event "k_hat won"
  #   winner_index - which.max(estimates)
  #   method       - string description

  M <- length(estimates)
  if (M < 2) stop("akm_winner_ci needs at least 2 candidates.")
  if (!is.matrix(cov_matrix) || nrow(cov_matrix) != M || ncol(cov_matrix) != M) {
    stop("cov_matrix must be M x M, where M = length(estimates)")
  }

  k_hat <- which.max(estimates)
  V_hat <- estimates[k_hat]
  sigma <- sqrt(cov_matrix[k_hat, k_hat])
  if (!is.finite(sigma) || sigma <= 0) {
    return(list(
      estimate = V_hat, ci = c(NA_real_, NA_real_),
      sigma = sigma, truncation = c(NA_real_, NA_real_),
      winner_index = k_hat,
      method = "AKM CI not available (winner has non-positive variance)"
    ))
  }

  # Build A so the selection event "k_hat won" is A V <= 0.
  others <- setdiff(seq_len(M), k_hat)
  A <- matrix(0, length(others), M)
  for (i in seq_along(others)) {
    A[i, others[i]] <-  1
    A[i, k_hat]     <- -1
  }

  # Lee polyhedral-lemma decomposition for test stat eta'V = V_hat_{k_hat}.
  eta <- numeric(M); eta[k_hat] <- 1
  c_vec <- as.vector(cov_matrix %*% eta) / sigma^2     # length M
  Ac <- as.vector(A %*% c_vec)                         # length M-1
  AV <- as.vector(A %*% estimates)                     # length M-1
  u  <- AV - Ac * V_hat                                # nuisance, fixed given V_hat

  upper_idx <- which(Ac > 0)
  lower_idx <- which(Ac < 0)
  V_plus  <- if (length(upper_idx) > 0) min(-u[upper_idx] / Ac[upper_idx]) else  Inf
  V_minus <- if (length(lower_idx) > 0) max(-u[lower_idx] / Ac[lower_idx]) else -Inf

  # Numerically stable CDF of N(mu, sigma^2) truncated to [a, b], at v.
  # Uses log-survival differencing so it stays accurate when mu << a (deeply
  # truncated case) where direct Phi differencing underflows.
  trunc_cdf <- function(v, mu, sigma, a, b) {
    if (sigma <= 0) return(NA_real_)
    if (v <= a) return(0)
    if (is.finite(b) && v >= b) return(1)
    log_Sa <- if (a == -Inf) 0       else pnorm((a - mu) / sigma, lower.tail = FALSE, log.p = TRUE)
    log_Sv <-                            pnorm((v - mu) / sigma, lower.tail = FALSE, log.p = TRUE)
    log_Sb <- if (b ==  Inf) -Inf    else pnorm((b - mu) / sigma, lower.tail = FALSE, log.p = TRUE)
    if (log_Sv > log_Sa) return(0)   # numerical: v at the boundary
    log_num   <- log_Sa + log1p(-exp(log_Sv - log_Sa))
    log_denom <- if (b == Inf) log_Sa else log_Sa + log1p(-exp(log_Sb - log_Sa))
    if (!is.finite(log_num)   || log_num   == -Inf) return(0)
    if (!is.finite(log_denom) || log_denom == -Inf) return(NA_real_)
    pmin(pmax(exp(log_num - log_denom), 0), 1)
  }

  alpha <- 1 - ci_level

  find_endpoint <- function(target_p) {
    f <- function(V) trunc_cdf(V_hat, V, sigma, V_minus, V_plus) - target_p
    lo <- V_hat - 50 * sigma
    hi <- V_hat + 50 * sigma
    f_lo <- f(lo); f_hi <- f(hi)
    if (is.na(f_lo) || is.na(f_hi)) return(NA_real_)
    if (f_lo * f_hi > 0) {
      return(if (abs(f_lo) < abs(f_hi)) lo else hi)
    }
    tryCatch(uniroot(f, c(lo, hi), tol = 1e-8)$root, error = function(e) NA_real_)
  }

  list(
    estimate     = V_hat,
    ci           = c(find_endpoint(1 - alpha / 2), find_endpoint(alpha / 2)),
    sigma        = sigma,
    truncation   = c(V_minus, V_plus),
    winner_index = k_hat,
    method       = "AKM/Lee polyhedral conditional CI"
  )
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
    costs = NULL,
    inference = c("none", "bootstrap"),
    n_boot = 1000,
    boot_seed = NULL,
    ci_level = 0.95
) {
  # takes data from a size n RCT with uniform assignment, estimates CATES, and returns an estimate of the optimal policy and return
  # inputs:
    # Y = length-n vector; outcomes in the RCT
    # Z = length-n binary vector; treatment assignment generically, 0 denotes control and 1 denotes treatment
    # X = n-by-p matrix; covariates in RCT
    # x_pop = N-by-p matrix; covariates in population; may be left out, in which case the science table for the experiment is returned instead
    # y_pop = N-by-2 matrix; optional true science table for the population as a matrix of potential outcomes on treatment (1st column) and control (2nd column); if ommitted, only the policy estimate is returned without the acutal return
    # budget = scalar; the budget for the portfolio; defaults to NULL, which is an unconstrained portfolio; for 'treat the K-best' set this to K and leave the costs as NULL
    # costs = length-N vector of positive reals; the cost of treating each unit in an additive cost model
    # inference = "none" (default) or "bootstrap"; if "bootstrap", a stratified
    #   bootstrap (with replacement, within each treatment arm) is used to
    #   construct percentile CIs for the value of the estimated optimal policy
    #   and the values of the two blanket policies. The whole pipeline (CATE
    #   re-estimation, LP re-solution, value re-computation) is bootstrapped
    #   to capture both estimation and policy-selection uncertainty.
    # n_boot = number of bootstrap replicates (default 1000)
    # boot_seed = optional integer seed for reproducible bootstrap draws
    # ci_level = nominal coverage level for bootstrap CIs (default 0.95)
  # outputs:
    # list with elements:
    # required
      # policy = an estimate of the "first-best" (unconstrained) optimal policy
      # cate = the conditional average treatment effect for every unit in the population (if x_pop is provided) or study
      # mu_0, mu_1 = predicted potential outcomes under control / treatment
      #   (NA-filled if the learner does not provide them, e.g. learner_rf)
      # predicted_value = mean(policy * mu_1 + (1 - policy) * mu_0); the
      #   plug-in average potential outcome under the estimated optimal policy
      # treat_all_value = mean(mu_1); value of the "treat every plot" blanket
      # treat_none_value = mean(mu_0); value of the "treat no plot" blanket
      # best_blanket_value = max(treat_all_value, treat_none_value)
    # optional (if y_pop, the science table for the population, is provided)
      # policy_value = the realized value of the estimated optimal policy
      # oracle_value = the realized value of the true optimal policy
      # regret = the realized regret of using the estimated optimal policy
    # optional (if inference = "bootstrap")
      # inference = list with bootstrap CI results (see end of function)
  inference <- match.arg(inference)

  cate_res <- estimate_cate(
    Y = Y,
    Z = Z,
    X = X,
    learner = learner,
    x_pop = x_pop
  )

  cate <- cate_res$cate
  mu_0 <- cate_res$mu_0
  mu_1 <- cate_res$mu_1

  policy <- solve_lp_policy(
    cate = cate,
    budget = budget,
    costs = costs
  )

  out <- list(
    policy = policy,
    cate   = cate,
    mu_0   = mu_0,
    mu_1   = mu_1
  )

  # Plug-in policy / blanket values (only if learner provides per-arm predictions)
  has_potentials <- !any(is.na(mu_0)) && !any(is.na(mu_1))
  if (has_potentials) {
    out$predicted_value    <- mean(policy * mu_1 + (1 - policy) * mu_0)
    out$treat_all_value    <- mean(mu_1)
    out$treat_none_value   <- mean(mu_0)
    out$best_blanket_value <- max(out$treat_all_value, out$treat_none_value)
  }

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

  if (inference == "bootstrap") {
    if (!is.null(boot_seed)) set.seed(boot_seed)
    n         <- length(Y)
    N_pop     <- length(cate)
    is_X_mat  <- is.matrix(X) || is.data.frame(X)

    # Pin x_pop to a concrete matrix so each bootstrap replicate evaluates the
    # CATE / mu_0 / mu_1 at the SAME population locations as the original fit.
    x_pop_eval <- if (is.null(x_pop)) X else x_pop

    idx_ctl <- which(Z == 0); n_ctl <- length(idx_ctl)
    idx_trt <- which(Z == 1); n_trt <- length(idx_trt)

    boot_cates    <- matrix(NA_real_, n_boot, N_pop)
    boot_mu_0     <- matrix(NA_real_, n_boot, N_pop)
    boot_mu_1     <- matrix(NA_real_, n_boot, N_pop)
    boot_policies <- matrix(NA_integer_, n_boot, N_pop)
    boot_predicted_value  <- rep(NA_real_, n_boot)
    boot_treat_all_value  <- rep(NA_real_, n_boot)
    boot_treat_none_value <- rep(NA_real_, n_boot)
    boot_best_blanket     <- rep(NA_real_, n_boot)
    boot_policy_value     <- if (!is.null(y_pop)) rep(NA_real_, n_boot) else NULL

    for (b in seq_len(n_boot)) {
      # Stratified resample with replacement: preserves treatment arm sizes,
      # mirroring the design's marginal assignment probabilities.
      idx_b <- c(
        sample(idx_ctl, n_ctl, replace = TRUE),
        sample(idx_trt, n_trt, replace = TRUE)
      )

      Y_b <- Y[idx_b]
      Z_b <- Z[idx_b]
      X_b <- if (is_X_mat) X[idx_b, , drop = FALSE] else X[idx_b]

      cate_b_res <- tryCatch(
        estimate_cate(
          Y = Y_b, Z = Z_b, X = X_b,
          learner = learner,
          x_pop = x_pop_eval
        ),
        error = function(e) NULL
      )
      if (is.null(cate_b_res)) next

      cate_b <- cate_b_res$cate
      mu_0_b <- cate_b_res$mu_0
      mu_1_b <- cate_b_res$mu_1
      if (any(!is.finite(cate_b))) next

      policy_b <- solve_lp_policy(
        cate = cate_b, budget = budget, costs = costs
      )

      boot_cates[b, ]    <- cate_b
      boot_mu_0[b, ]     <- mu_0_b
      boot_mu_1[b, ]     <- mu_1_b
      boot_policies[b, ] <- policy_b

      if (all(!is.na(mu_0_b)) && all(!is.na(mu_1_b))) {
        boot_predicted_value[b]  <- mean(policy_b * mu_1_b + (1 - policy_b) * mu_0_b)
        boot_treat_all_value[b]  <- mean(mu_1_b)
        boot_treat_none_value[b] <- mean(mu_0_b)
        boot_best_blanket[b]     <- max(boot_treat_all_value[b], boot_treat_none_value[b])
      }

      if (!is.null(y_pop)) {
        boot_policy_value[b] <- mean(
          ifelse(policy_b == 1, y_pop[,1], y_pop[,2])
        )
      }
    }

    alpha <- 1 - ci_level
    qs    <- c(alpha / 2, 1 - alpha / 2)
    # quantile function
    qfn   <- function(x) {
      if (all(is.na(x))) return(c(NA_real_, NA_real_))
      stats::quantile(x, probs = qs, names = FALSE, na.rm = TRUE)
    }

    inf <- list(
      method     = "stratified bootstrap (resample within each treatment arm)",
      n_boot     = n_boot,
      n_completed = sum(!is.na(boot_predicted_value) | !is.na(boot_policy_value %||% NA)),
      ci_level   = ci_level,
      # Per-unit CATE bootstrap CIs (alternative to the analytic SEs)
      cate_lower = apply(boot_cates, 2, function(x) qfn(x)[1]),
      cate_upper = apply(boot_cates, 2, function(x) qfn(x)[2]),
      # Treatment-stability score: proportion of bootstrap reps in which each
      # population unit is treated under the LP-solved optimal policy. Useful
      # for reporting which plots are robustly in / out of the optimal policy.
      treatment_stability = colMeans(boot_policies, na.rm = TRUE),
      # Bootstrap raw samples (rows = reps, cols = population units)
      boot_cates    = boot_cates,
      boot_mu_0     = boot_mu_0,
      boot_mu_1     = boot_mu_1,
      boot_policies = boot_policies
    )

    if (has_potentials) {
      inf$predicted_value_ci    <- qfn(boot_predicted_value)
      inf$treat_all_value_ci    <- qfn(boot_treat_all_value)
      inf$treat_none_value_ci   <- qfn(boot_treat_none_value)
      inf$best_blanket_value_ci <- qfn(boot_best_blanket)
      inf$boot_predicted_value  <- boot_predicted_value
      inf$boot_treat_all_value  <- boot_treat_all_value
      inf$boot_treat_none_value <- boot_treat_none_value
      inf$boot_best_blanket_value <- boot_best_blanket

      # Andrews-Kitagawa-McCloskey conditional CIs ("Inference on Winners")
      # using the bootstrap-estimated covariance of the policy-value vector.
      # Two natural winner problems arise here:
      #   1. "Best blanket" = argmax(V_TA, V_TN). The naive percentile/quantile
      #      CI for max(V_hat_TA, V_hat_TN) may be upward-biased; AKM corrects it.
      #   2. "Best policy among {treat-all, treat-none, LP-solved optimal}".
      #      For the unconstrained LP this is essentially trivial because the
      #      LP is mathematically guaranteed to dominate both blankets, so the
      #      AKM CI usually collapses to the naive Wald CI for V_opt. For the
      #      budget-constrained LP, V_opt may not strictly dominate the
      #      blankets and the AKM correction can matter.
      #
      # Note: AKM here treats the LP-output policy as fixed; it does NOT
      # correct for the within-LP selection of which plots to treat. For that,
      # a more elaborate selective-inference argument over the LP polyhedron
      # would be needed (not yet implemented).
      V_boot_mat <- cbind(boot_treat_all_value, boot_treat_none_value, boot_predicted_value)
      complete   <- stats::complete.cases(V_boot_mat)
      if (sum(complete) >= 30) {
        Sigma_full <- stats::cov(V_boot_mat[complete, , drop = FALSE])

        # AKM CI for "best blanket" using the (V_TA, V_TN) sub-covariance
        Sigma_bb <- Sigma_full[1:2, 1:2]
        akm_bb <- akm_winner_ci(
          c(out$treat_all_value, out$treat_none_value),
          Sigma_bb,
          ci_level = ci_level
        )
        inf$best_blanket_akm_ci    <- akm_bb$ci
        inf$best_blanket_akm_winner <- if (akm_bb$winner_index == 1L) "treat_all" else "treat_none"

        # AKM CI for predicted-value winner among {V_TA, V_TN, V_opt}
        akm_opt <- akm_winner_ci(
          c(out$treat_all_value, out$treat_none_value, out$predicted_value),
          Sigma_full,
          ci_level = ci_level
        )
        inf$predicted_value_akm_ci    <- akm_opt$ci
        inf$predicted_value_akm_winner <- c("treat_all", "treat_none", "lp_optimal")[akm_opt$winner_index]
      } else {
        inf$best_blanket_akm_ci   <- c(NA_real_, NA_real_)
        inf$predicted_value_akm_ci <- c(NA_real_, NA_real_)
      }
    }
    if (!is.null(y_pop)) {
      inf$policy_value_ci   <- qfn(boot_policy_value)
      inf$boot_policy_value <- boot_policy_value
    }

    out$inference <- inf
  }

  out
}


# NB: experimental, not currently used
policy_value_cv <- function(
    Y,
    Z,
    X,
    learner = learner_ols,
    x_pop = NULL,
    budget = NULL,
    costs = NULL,
    K = 2,
    n_repeats = 1,
    split_seed = NULL,
    ci_level = 0.95
) {
  # Cross-fit AIPW evaluation of the value of the (procedure that estimates the)
  # optimal policy. For each of K folds:
  #   1. Estimate the CATE / mu_0 / mu_1 (via `learner`) on the K-1 training folds.
  #   2. Solve the LP for the optimal policy at the held-out fold's covariates.
  #   3. Score each held-out unit with the AIPW estimator
  #        V_i = mu_pi(X_i) + (Z_i == pi(X_i)) * (Y_i - mu_pi(X_i)) / e(Z_i)
  #      where mu_pi(x) = mu_{pi(x)}(x), e(1) = mean(Z), e(0) = 1 - mean(Z).
  # The mean of V_i over all units is an unbiased estimator (in the RCT) of the
  # value of the policy-selection procedure when applied to a fresh sample of
  # size n_train, evaluated on a fresh sample. Because each unit's score uses
  # only training-fold information for mu_hat and pi, the standard error
  # var(V_i) / n is approximately valid (ignoring weak fold dependence).
  #
  # When n_repeats > 1, K-fold splits are drawn n_repeats times and each unit's
  # score is averaged across the splits before forming the mean / SE.
  #
  # Unlike the bootstrap inference in estimate_optimal_policy(), this CI does
  # NOT rely on the bootstrap being consistent for the value of the optimum
  # (which it generally is not, due to non-differentiability of the max);
  # validity follows directly from sample splitting and the known propensity
  # in an RCT. See e.g. Athey & Wager (2021, Econometrica) for the AIPW score.
  #
  # Inputs:
  #   Y, Z, X    - study data (as in estimate_optimal_policy)
  #   learner    - CATE learner; must return $mu_0 and $mu_1 (learner_ols,
  #                learner_dim, learner_rf all do)
  #   x_pop      - currently unused; the cross-fit estimator scores observed
  #                study units (so the population is implicitly the study).
  #   budget     - if non-NULL, scaled proportionally on each held-out fold:
  #                budget_test = ceiling(budget * n_test / n).
  #   costs      - if non-NULL, sliced to the held-out fold each iteration.
  #   K          - number of folds (default 2; use larger K for larger n).
  #   n_repeats  - number of independent K-fold splits to average over.
  #   split_seed - optional integer seed for the fold assignments.
  #   ci_level   - nominal coverage of the returned Wald CI.
  #
  # Output: list with
  #   estimate, se, ci             - cross-fit AIPW value and Wald CI
  #   method                       - string description
  #   K, n_repeats                 - settings
  #   full_policy, full_predicted_value
  #                                - the policy and plug-in value computed on
  #                                  the full data, for reference
  #   treatment_stability          - per-unit fraction of folds (across all
  #                                  repeats) in which the unit was treated
  #                                  by the held-out policy
  #   scores                       - per-unit averaged AIPW scores

  if (!is.null(split_seed)) set.seed(split_seed)
  if (!is.null(x_pop)) {
    warning("x_pop is currently ignored by policy_value_cv (the held-out fold ",
            "is implicitly the population). The policy is still re-fit on the ",
            "full data with x_pop = NULL for the `full_policy` field.")
  }

  n <- length(Y)
  prop_score <- mean(Z)
  is_X_mat   <- is.matrix(X) || is.data.frame(X)

  # Reference: full-data policy + plug-in value (no inference).
  full_res <- estimate_optimal_policy(
    Y = Y, Z = Z, X = X,
    learner = learner,
    budget  = budget,
    costs   = costs
  )
  if (any(is.na(full_res$mu_0)) || any(is.na(full_res$mu_1))) {
    stop("Learner must return mu_0 and mu_1 for policy_value_cv. ",
         "(learner_ols, learner_dim, learner_rf all support this.)")
  }

  score_mat  <- matrix(NA_real_,    n_repeats, n)
  policy_mat <- matrix(NA_integer_, n_repeats, n)

  for (rep_idx in seq_len(n_repeats)) {

    # Stratified K-fold assignment: balance treatment groups across folds
    fold_id <- integer(n)
    for (z in c(0, 1)) {
      ix          <- which(Z == z)
      fold_id[ix] <- sample(rep(seq_len(K), length.out = length(ix)))
    }

    for (k in seq_len(K)) {
      train <- which(fold_id != k)
      test  <- which(fold_id == k)

      Y_tr <- Y[train]; Z_tr <- Z[train]
      X_tr <- if (is_X_mat) X[train, , drop = FALSE] else X[train]
      Y_te <- Y[test];  Z_te <- Z[test]
      X_te <- if (is_X_mat) X[test,  , drop = FALSE] else X[test]

      res_tr <- learner(Y = Y_tr, Z = Z_tr, X = X_tr, x_pop = X_te)

      cate_te <- res_tr$cate
      mu_0_te <- res_tr$mu_0
      mu_1_te <- res_tr$mu_1
      if (is.null(mu_0_te) || is.null(mu_1_te) ||
          any(is.na(mu_0_te)) || any(is.na(mu_1_te))) {
        stop("Learner must return mu_0 and mu_1 for policy_value_cv.")
      }

      # Scale the budget proportionally on the held-out fold; this matches
      # the "treat top X% of population" interpretation, but does change the
      # absolute number of treated units relative to the full-data LP.
      budget_te <- if (!is.null(budget)) max(0L, ceiling(budget * length(test) / n)) else NULL
      costs_te  <- if (!is.null(costs))  costs[test] else NULL

      policy_te <- solve_lp_policy(
        cate = cate_te, budget = budget_te, costs = costs_te
      )

      mu_pi <- ifelse(policy_te == 1, mu_1_te, mu_0_te)
      e_obs <- ifelse(Z_te == 1, prop_score, 1 - prop_score)

      score_mat[rep_idx, test]  <- mu_pi + (Z_te == policy_te) * (Y_te - mu_pi) / e_obs
      policy_mat[rep_idx, test] <- policy_te
    }
  }

  # Per-unit averaged scores across repeats; mean and Wald CI.
  unit_scores         <- colMeans(score_mat,  na.rm = TRUE)
  treatment_stability <- colMeans(policy_mat, na.rm = TRUE)

  est   <- mean(unit_scores)
  se    <- sqrt(stats::var(unit_scores) / n)
  alpha <- 1 - ci_level
  zcrit <- qnorm(1 - alpha / 2)

  list(
    estimate              = est,
    se                    = se,
    ci                    = c(est - zcrit * se, est + zcrit * se),
    method                = sprintf(
      "%d-fold AIPW (cross-fit, averaged over %d split%s)",
      K, n_repeats, if (n_repeats == 1) "" else "s"
    ),
    K                     = K,
    n_repeats             = n_repeats,
    full_policy           = full_res$policy,
    full_predicted_value  = full_res$predicted_value,
    treatment_stability   = treatment_stability,
    scores                = unit_scores
  )
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






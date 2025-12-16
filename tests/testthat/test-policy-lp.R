# test linear program solver for optimal policies

test_that("unconstrained policy matches sign of CATE", {
  cate <- c(-2, -1, 0, 1, 3)
  
  policy <- solve_lp_policy(cate)
  
  expect_equal(policy, c(0, 0, 0, 1, 1))
})


test_that("budget constraint treats K best units", {
  cate <- c(5, 4, 3, 2, 1)
  K <- 2
  
  policy <- solve_lp_policy(cate, budget = K)
  
  expect_equal(sum(policy), K)
  
  treated_indices <- which(policy == 1)
  expect_equal(treated_indices, c(1, 2))
})


test_that("cost constraint respects budget", {
  cate  <- c(10, 8, 6, 4)
  costs <- c(5, 4, 3, 1)
  budget <- 6
  
  policy <- solve_lp_policy(
    cate = cate,
    costs = costs,
    budget = budget
  )
  
  total_cost <- sum(costs * policy)
  
  expect_true(total_cost <= budget)
  
  # Should pick the best affordable combo
  # Unit 1 (cost 5, cate 10) dominates unit 2 (cost 4, cate 8)
  expect_true(all(policy %in% c(0, 1))) # ensure policy is binary
  expect_true(sum(costs * policy) <= budget) # ensure budget is met
  expect_equal(policy, c(1, 0, 0, 1)) # check resulting policy
})


test_that("oracle value >= constrained policy value", {
  set.seed(1)
  n <- 100
  X <- matrix(rnorm(n), n, 1)
  Z <- rbinom(n, 1, 0.5)
  
  tau <- rnorm(n)
  y0 <- rnorm(n)
  y1 <- y0 + tau
  Y  <- ifelse(Z == 1, y1, y0)
  
  res_unconstrained <- estimate_optimal_policy(
    Y = Y,
    Z = Z,
    X = X,
    learner = learner_ols,
    y_pop = cbind(y1, y0)
  )
  
  res_constrained <- estimate_optimal_policy(
    Y = Y,
    Z = Z,
    X = X,
    learner = learner_ols,
    budget = 30,
    y_pop = cbind(y1, y0)
  )
  
  expect_true(res_unconstrained$oracle_value >= res_constrained$oracle_value)
})

# test simulation function 

test_that("run_policy_simulation returns correct structure", {
  res <- run_policy_simulation(
    dgp = dgp_constant_tau,
    n = 100,
    learner = learner_ols,
    num_sims = 10,
    seed = 1
  )
  
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 10)
  
  expect_true(all(c(
    "policy_value",
    "oracle_value",
    "regret",
    "treat_rate"
  ) %in% names(res)))
})

test_that("simulation is reproducible with fixed seed", {
  res1 <- run_policy_simulation(
    dgp = dgp_constant_tau,
    n = 80,
    num_sims = 5,
    learner = learner_ols,
    seed = 123
  )
  
  res2 <- run_policy_simulation(
    dgp = dgp_constant_tau,
    n = 80,
    num_sims = 5,
    learner = learner_ols,
    seed = 123
  )
  
  expect_equal(res1, res2)
})


test_that("zero-treatment-effect DGP yields zero regret", {
  res <- run_policy_simulation(
    dgp = dgp_zero_tau,
    n = 200,
    num_sims = 20,
    learner = learner_ols,
    seed = 1
  )
  
  expect_true(all(abs(res$regret) < 1e-6))
})


test_that("oracle value is always >= policy value", {
  res <- run_policy_simulation(
    dgp = dgp_constant_tau,
    n = 150,
    num_sims = 20,
    learner = learner_ols,
    seed = 2
  )
  
  expect_true(all(res$oracle_value >= res$policy_value))
  expect_true(all(res$regret >= 0))
})


test_that("budget constraint respected in simulation", {
  K <- 30
  
  res <- run_policy_simulation(
    dgp = dgp_constant_tau,
    n = 100,
    num_sims = 10,
    learner = learner_ols,
    policy_args = list(budget = K),
    seed = 3
  )
  
  expect_true(all(res$treat_rate <= K / 100 + 1e-6))
})

test_that("simulation runs for very small n", {
  res <- run_policy_simulation(
    dgp = dgp_constant_tau,
    n = 10,
    num_sims = 5,
    learner = learner_ols
  )
  
  expect_equal(nrow(res), 5)
})




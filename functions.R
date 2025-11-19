library(tidyverse)
library(data.table)
library(gstat)
#permuter can be downloaded from Github by following README: https://github.com/statlab/permuter
library(permuter)
#library(kolmim)




num_sims <- 200

run_simulation_simple <- function(param_list, ate, moderator_effect, idiosyncratic_teh, measurement_noise, n, num_sims = 200){
  # takes a variety of parameters as inputs, generates potential outcomes, runs a completely randomized experiments, and estimates those parameters
  # inputs:
    #param_list:
      # baseline_avg = the average SOC at baseline across plots
      # baseline_ap_var = the variability of baseline SOC across plots
      # mean_change_ctl = the average change between baseline and follow-up for control plots
      # ap_var_change_ctl = the across-plot variance among control plots
      # ate = average treatment effect
      # moderator_effect = the effect of a single continuous (gaussian) moderator on treatment effects, which could represent baseline SOC for instance
      # idiosyncratic_teh = the variability of 'random' treatment effect heterogeneity _not_ due to the moderator
      # measurement_noise = the amount of measurement noise, representing both assay error and within-plot heterogeneity
      # n = the size of the experiment (must be divisible by 2)
    # num_sims:
      # the number of simulation iterations
  # returns:
    # a data.frame of various results summarizing the performance of estimators over the num_sims simulations
    
  attach(param_list)
  
  #simulate at population level
  b <- rnorm(N, mean = baseline_avg, sd = baseline_ap_var) #baseline %TC
  teh_noise_i <- rnorm(N, mean = 0, sd = sqrt(idiosyncratic_teh))
  y_0 <- b + rnorm(N, mean = mean_change_ctl, sd = sqrt(ap_var_change_ctl))
  y_1 <- y_0 + ate + moderator_effect * stdize(b) + teh_noise_i
  #parameters to be estimated
  pate <- mean(y_1) - mean(y_0) #(finite) population average treatment effect
  fpme <- coef(lm((y_1 - y_0) ~ stdize(b)))['stdize(b)'] #finite-population moderator effect
  #for now, assume the policy is not budget constrained
  oracle_mean <- sum(pmax(y_0, y_1))/N #using actual ITEs
  restricted_oracle_mean <- max(sum(y_1), sum(y_0))/N #total for best restricted regime (treat or don't treat whole population)
  
  #storage for results
  #each result is a matrix with columns ['estimate','lower_ci','upper_ci']
  diff_means_results <- matrix(0, nrow = num_sims, ncol = 3) #difference-in-means estimate of PATE
  diff_diffs_results <- matrix(0, nrow = num_sims, ncol = 3) #difference-in-differences estimate of PATE
  ols_est_results <- matrix(0, nrow = num_sims, ncol = 3) #OLS with full baseline interaction estimate of PATE
  mod_est_results <- matrix(0, nrow = num_sims, ncol = 3) #OLS estimate of moderator effect
  naive_mod_est_results <- matrix(0, nrow = num_sims, ncol = 3) #"naive" estimate from regressing followup SOCs on differences
  optimal_policy_results <- rep(0, num_sims) #OLS estimate of optimal regime by ITE prediction
  restricted_optimal_policy_results <- rep(0, num_sims)
  for(i in 1:num_sims){
    #study level potential outcomes and covariates
    sample_index <- sample(1:N, size = n, replace = FALSE)
    #sample from population and add measurement noise from core sampling (currently commented out)
    baseline_noise <-  rnorm(n, mean = 0, sd = sqrt(measurement_noise))
    B <- b[sample_index] + baseline_noise
    B_std <- stdize(b)[sample_index] + baseline_noise
    Y_0 <- y_0[sample_index] + rnorm(n, mean = 0, sd = sqrt(measurement_noise))
    Y_1 <- y_1[sample_index] + rnorm(n, mean = 0, sd = sqrt(measurement_noise))
    #balanced complete randomization
    trt_index <- sample(1:n, size = n/2, replace = FALSE)
    
    #estimate average treatment effect
    #difference in means
    diff_means <- mean(Y_1[trt_index]) - mean(Y_0[-trt_index])
    diff_means_se <- sqrt(var(Y_1[trt_index])/(n/2) + var(Y_0[trt_index])/(n/2)) 
    diff_means_ci <- c(diff_means - qnorm(0.975) * diff_means_se, diff_means + qnorm(0.975) * diff_means_se)
    diff_means_results[i,] <- c(diff_means, diff_means_ci)
    #difference in differences
    diff_diffs <- mean(Y_1[trt_index] - B[trt_index]) - mean(Y_0[-trt_index] - B[-trt_index])
    diff_diffs_se <- sqrt(var(Y_1[trt_index] - B[trt_index])/(n/2) + var(Y_0[-trt_index] - B[-trt_index])/(n/2))
    diff_diffs_ci <- c(diff_diffs - qnorm(0.975) * diff_diffs_se, diff_diffs + qnorm(0.975) * diff_diffs_se)
    diff_diffs_results[i,] <- c(diff_diffs, diff_diffs_ci)
    #OLS estimate
    Z <- rep(0, n) #treatment dummy
    Z[trt_index] <- 1
    Y <- Y_0 #pooled observed data
    Y[trt_index] <- Y_1[trt_index]
    ols <- lm(Y ~ Z*B_std) #fully interacted OLS, ala Lin (2013) "Agnostic notes..."
    ols_est <- summary(ols)$coefficients["Z","Estimate"]
    ols_se <- summary(ols)$coefficients["Z","Std. Error"] #can use usual formula for SE (not sandwich SE) because experiment is balanced; see remark (ii), page 307 of Lin (2013)
    ols_ci <- c(ols_est - qnorm(0.975) * ols_se, ols_est + qnorm(0.975) * ols_se)
    ols_est_results[i,] <- c(ols_est, ols_ci)
    study_data <- data.frame(Y=Y,Z=Z,B_std=B_std,B=B)
    
    #estimate moderator effects
    #using the Ding et al package "hettx" 
    sys_model <- summary(estimate_systematic(Y ~ Z, interaction.formula = ~ B_std, method = "OLS", data = study_data))
    # ols_est <- sys_model$ATE
    # ols_se <- sys_model$SE.ATE
    # ols_ci <- c(ols_est - qnorm(0.975) * ols_se, ols_est + qnorm(0.975) * ols_se)
    # ols_est_results[i,] <- c(ols_est, ols_ci)
    #moderator effects
    mod_est <- sys_model$coefficients["B_std"]
    mod_se <- sqrt(sys_model$vcov[2,2])
    mod_ci <- c(mod_est - qnorm(0.975) * mod_se, mod_est + qnorm(0.975) * mod_se)
    mod_est_results[i,] <- c(mod_est, mod_ci)
    #naive estimates of "baseline effects" 
    D <- Y - B
    naive_model <- summary(lm(Y ~ D))
    naive_mod_est <- naive_model$coefficients["D","Estimate"]
    naive_mod_se <- naive_model$coefficients["D","Std. Error"]
    naive_mod_ci <- c(naive_mod_est - qnorm(0.975) * naive_mod_se, naive_mod_est + qnorm(0.975) * naive_mod_se)
    naive_mod_est_results[i,] <- c(naive_mod_est, naive_mod_ci)
    
    
    #estimate optimal policy
    #predictions of ITEs and estimates of oracle and restricted oracle policies
    hat_y_1 <- predict(ols, newdata = data.frame(Z = 1, B_std = stdize(b)), interval = "prediction")
    hat_y_0 <- predict(ols, newdata = data.frame(Z = 0, B_std = stdize(b)), interval = "prediction")
    est_po_matrix <- cbind(hat_y_1[,1], hat_y_0[,1]) #estimated potential outcome matrix
    estimated_policy <- apply(est_po_matrix, 1, which.max) #pick out the estimated best treatment for each plot
    estimated_optimal_mean <- sum(cbind(y_1, y_0)[cbind(1:N,estimated_policy)])/N #return for estimated policy
    optimal_policy_results[i] <- estimated_optimal_mean 
    
    estimated_restricted_optimal_mean <- c(mean(y_0), mean(y_1))[which.max(c(mean(Y_0), mean(Y_1)))] #return for estimated restricted optimal policy
    restricted_optimal_policy_results[i] <- estimated_restricted_optimal_mean
  }
  results <- data.frame(
    spate = ate, #super-population ATE
    fpate = pate, #finite-population ATE
    spme = moderator_effect,
    fpme = fpme,
    oracle_mean = oracle_mean,
    restricted_oracle_mean = restricted_oracle_mean,
    teh = idiosyncratic_teh, 
    measurementnoise = measurement_noise,
    n = n,
    diffmeans_bias = mean(diff_means_results[,1]) - pate,
    diffmeans_rmse = sqrt(mean((diff_means_results[,1] - pate)^2)),
    diffmeans_coverage = mean((diff_means_results[,2]) < pate & (diff_means_results[,3] > pate)),
    diffmeans_ciwidth = mean(diff_means_results[,3] - diff_means_results[,2]),
    diffmeans_power = mean(diff_means_results[,2] > 0),
    diffdiffs_bias = mean(diff_diffs_results[,1]) - pate,
    diffdiffs_rmse = sqrt(mean((diff_diffs_results[,1] - pate)^2)),
    diffdiffs_coverage = mean((diff_diffs_results[,2]) < pate & (diff_diffs_results[,3] > pate)),
    diffdiffs_ciwidth = mean(diff_diffs_results[,3] - diff_diffs_results[,2]),
    diffdiffs_power = mean(diff_diffs_results[,2] > 0),
    olsest_bias = mean(ols_est_results[,1]) - pate,
    olsest_rmse = sqrt(mean((ols_est_results[,1] - pate)^2)),
    olsest_coverage = mean((ols_est_results[,2]) < pate & (ols_est_results[,3] > pate)),
    olsest_ciwidth = mean(ols_est_results[,3] - ols_est_results[,2]),
    olsest_power = mean(ols_est_results[,2] > 0),
    modest_bias = mean(mod_est_results[,1]) - fpme,
    modest_rmse = sqrt(mean((mod_est_results[,1] - moderator_effect)^2)),
    modest_coverage = mean((mod_est_results[,2]) < moderator_effect & (mod_est_results[,3] > moderator_effect)),
    modest_ciwidth = mean(mod_est_results[,3] - mod_est_results[,2]),
    naive_modest_bias = mean(naive_mod_est_results[,1]) - fpme,
    naive_modest_rmse = sqrt(mean((naive_mod_est_results[,1] - moderator_effect)^2)),
    naive_modest_coverage = mean((naive_mod_est_results[,2]) < moderator_effect & (naive_mod_est_results[,3] > moderator_effect)),
    naive_modest_ciwidth = mean(naive_mod_est_results[,3] - naive_mod_est_results[,2]),
    avg_optimal_policy_estimate = mean(optimal_policy_results),
    avg_restricted_optimal_policy_estimate = mean(restricted_optimal_policy_results)
  )
  results
}


run_simulation_complex <- function(param_list, num_sims = 200){
  #generate the potential outcomes one time as a finite population
  b <- rnorm(N, mean = baseline_avg, sd = baseline_ap_var) # baseline %SOC
  c <- rbinom(N, size = 1, prob = 0.5) # bernoulli covariate
  d <- rnorm(N, mean = 0, sd = 1) # continuous covariate
  teh_noise_i <- rnorm(N, mean = 0, sd = sqrt(idiosyncratic_teh)) # idosyncratic treatment effect heterogneeity 
  y_0 <- b + rnorm(N, mean = mean_change_ctl, sd = sqrt(ap_var_change_ctl))
  y_1 <- y_0 + ate + moderator_effect * stdize(b) - c + c*d + teh_noise_i
  #parameters to be estimated
  pate <- mean(y_1) - mean(y_0) # (finite) population average treatment effect
  fpme <- coef(lm((y_1 - y_0) ~ stdize(b)))['stdize(b)'] #finite-population moderator effect
  #for now, assume the policy is not budget constrained
  oracle_mean <- sum(pmax(y_0, y_1))/N #using actual ITEs
  restricted_oracle_mean <- max(sum(y_1), sum(y_0))/N #total for best restricted regime (treat or don't treat whole population)
  
  #storage for results
  #each result is a matrix with columns ['estimate','lower_ci','upper_ci']
  diff_means_results <- matrix(0, nrow = num_sims, ncol = 3) #difference-in-means estimate of PATE
  diff_diffs_results <- matrix(0, nrow = num_sims, ncol = 3) #difference-in-differences estimate of PATE
  ols_est_results <- matrix(0, nrow = num_sims, ncol = 3) #OLS with full baseline interaction estimate of PATE
  mod_est_results <- matrix(0, nrow = num_sims, ncol = 3) #OLS estimate of moderator effect
  naive_mod_est_results <- matrix(0, nrow = num_sims, ncol = 3) #"naive" estimate from regressing followup SOCs on differences
  optimal_policy_results <- rep(0, num_sims) #OLS estimate of optimal regime by ITE prediction
  restricted_optimal_policy_results <- rep(0, num_sims)
  for(i in 1:num_sims){
    #study level potential outcomes and covariates
    sample_index <- sample(1:N, size = n, replace = FALSE)
    #sample from population and add measurement noise from core sampling (currently commented out)
    baseline_noise <-  rnorm(n, mean = 0, sd = sqrt(measurement_noise))
    B <- b[sample_index] + baseline_noise
    B_std <- stdize(b)[sample_index] + baseline_noise
    Y_0 <- y_0[sample_index] + rnorm(n, mean = 0, sd = sqrt(measurement_noise))
    Y_1 <- y_1[sample_index] + rnorm(n, mean = 0, sd = sqrt(measurement_noise))
    #balanced complete randomization
    trt_index <- sample(1:n, size = n/2, replace = FALSE)
    
    #estimate average treatment effect
    #difference in means
    diff_means <- mean(Y_1[trt_index]) - mean(Y_0[-trt_index])
    diff_means_se <- sqrt(var(Y_1[trt_index])/(n/2) + var(Y_0[trt_index])/(n/2)) 
    diff_means_ci <- c(diff_means - qnorm(0.975) * diff_means_se, diff_means + qnorm(0.975) * diff_means_se)
    diff_means_results[i,] <- c(diff_means, diff_means_ci)
    #difference in differences
    diff_diffs <- mean(Y_1[trt_index] - B[trt_index]) - mean(Y_0[-trt_index] - B[-trt_index])
    diff_diffs_se <- sqrt(var(Y_1[trt_index] - B[trt_index])/(n/2) + var(Y_0[-trt_index] - B[-trt_index])/(n/2))
    diff_diffs_ci <- c(diff_diffs - qnorm(0.975) * diff_diffs_se, diff_diffs + qnorm(0.975) * diff_diffs_se)
    diff_diffs_results[i,] <- c(diff_diffs, diff_diffs_ci)
    #OLS estimate
    Z <- rep(0, n) #treatment dummy
    Z[trt_index] <- 1
    Y <- Y_0 #pooled observed data
    Y[trt_index] <- Y_1[trt_index]
    ols <- lm(Y ~ Z*B_std) #fully interacted OLS, ala Lin (2013) "Agnostic notes..."
    ols_est <- summary(ols)$coefficients["Z","Estimate"]
    ols_se <- summary(ols)$coefficients["Z","Std. Error"] #can use usual formula for SE (not sandwich SE) because experiment is balanced; see remark (ii), page 307 of Lin (2013)
    ols_ci <- c(ols_est - qnorm(0.975) * ols_se, ols_est + qnorm(0.975) * ols_se)
    ols_est_results[i,] <- c(ols_est, ols_ci)
    study_data <- data.frame(Y=Y,Z=Z,B_std=B_std,B=B)
    
    #estimate moderator effects
    #using the Ding et al package "hettx" 
    sys_model <- summary(estimate_systematic(Y ~ Z, interaction.formula = ~ B_std, method = "OLS", data = study_data))
    # ols_est <- sys_model$ATE
    # ols_se <- sys_model$SE.ATE
    # ols_ci <- c(ols_est - qnorm(0.975) * ols_se, ols_est + qnorm(0.975) * ols_se)
    # ols_est_results[i,] <- c(ols_est, ols_ci)
    #moderator effects
    mod_est <- sys_model$coefficients["B_std"]
    mod_se <- sqrt(sys_model$vcov[2,2])
    mod_ci <- c(mod_est - qnorm(0.975) * mod_se, mod_est + qnorm(0.975) * mod_se)
    mod_est_results[i,] <- c(mod_est, mod_ci)
    #naive estimates of "baseline effects" 
    D <- Y - B
    naive_model <- summary(lm(Y ~ D))
    naive_mod_est <- naive_model$coefficients["D","Estimate"]
    naive_mod_se <- naive_model$coefficients["D","Std. Error"]
    naive_mod_ci <- c(naive_mod_est - qnorm(0.975) * naive_mod_se, naive_mod_est + qnorm(0.975) * naive_mod_se)
    naive_mod_est_results[i,] <- c(naive_mod_est, naive_mod_ci)
    
    
    #estimate optimal policy
    #predictions of ITEs and estimates of oracle and restricted oracle policies
    hat_y_1 <- predict(ols, newdata = data.frame(Z = 1, B_std = stdize(b)), interval = "prediction")
    hat_y_0 <- predict(ols, newdata = data.frame(Z = 0, B_std = stdize(b)), interval = "prediction")
    est_po_matrix <- cbind(hat_y_1[,1], hat_y_0[,1]) #estimated potential outcome matrix
    estimated_policy <- apply(est_po_matrix, 1, which.max) #pick out the estimated best treatment for each plot
    estimated_optimal_mean <- sum(cbind(y_1, y_0)[cbind(1:N,estimated_policy)])/N #return for estimated policy
    optimal_policy_results[i] <- estimated_optimal_mean 
    
    estimated_restricted_optimal_mean <- c(mean(y_0), mean(y_1))[which.max(c(mean(Y_0), mean(Y_1)))] #return for estimated restricted optimal policy
    restricted_optimal_policy_results[i] <- estimated_restricted_optimal_mean
  }
  results <- data.frame(
    spate = ate, #super-population ATE
    fpate = pate, #finite-population ATE
    spme = moderator_effect,
    fpme = fpme,
    oracle_mean = oracle_mean,
    restricted_oracle_mean = restricted_oracle_mean,
    teh = idiosyncratic_teh, 
    measurementnoise = measurement_noise,
    n = n,
    diffmeans_bias = mean(diff_means_results[,1]) - pate,
    diffmeans_rmse = sqrt(mean((diff_means_results[,1] - pate)^2)),
    diffmeans_coverage = mean((diff_means_results[,2]) < pate & (diff_means_results[,3] > pate)),
    diffmeans_ciwidth = mean(diff_means_results[,3] - diff_means_results[,2]),
    diffmeans_power = mean(diff_means_results[,2] > 0),
    diffdiffs_bias = mean(diff_diffs_results[,1]) - pate,
    diffdiffs_rmse = sqrt(mean((diff_diffs_results[,1] - pate)^2)),
    diffdiffs_coverage = mean((diff_diffs_results[,2]) < pate & (diff_diffs_results[,3] > pate)),
    diffdiffs_ciwidth = mean(diff_diffs_results[,3] - diff_diffs_results[,2]),
    diffdiffs_power = mean(diff_diffs_results[,2] > 0),
    olsest_bias = mean(ols_est_results[,1]) - pate,
    olsest_rmse = sqrt(mean((ols_est_results[,1] - pate)^2)),
    olsest_coverage = mean((ols_est_results[,2]) < pate & (ols_est_results[,3] > pate)),
    olsest_ciwidth = mean(ols_est_results[,3] - ols_est_results[,2]),
    olsest_power = mean(ols_est_results[,2] > 0),
    modest_bias = mean(mod_est_results[,1]) - fpme,
    modest_rmse = sqrt(mean((mod_est_results[,1] - moderator_effect)^2)),
    modest_coverage = mean((mod_est_results[,2]) < moderator_effect & (mod_est_results[,3] > moderator_effect)),
    modest_ciwidth = mean(mod_est_results[,3] - mod_est_results[,2]),
    naive_modest_bias = mean(naive_mod_est_results[,1]) - fpme,
    naive_modest_rmse = sqrt(mean((naive_mod_est_results[,1] - moderator_effect)^2)),
    naive_modest_coverage = mean((naive_mod_est_results[,2]) < moderator_effect & (naive_mod_est_results[,3] > moderator_effect)),
    naive_modest_ciwidth = mean(naive_mod_est_results[,3] - naive_mod_est_results[,2]),
    avg_optimal_policy_estimate = mean(optimal_policy_results),
    avg_restricted_optimal_policy_estimate = mean(restricted_optimal_policy_results)
  )
  results
}
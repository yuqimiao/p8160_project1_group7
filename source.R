library(MASS)
library(glmnet)
library(tidyverse)

# this function is to calculate lasso function by cv
lasso_function = function(data, response){
  lambda_best = glmnet::cv.glmnet(x = data, y = response, alpha = 1, nfolds = 5) #cv to find best lambda
  lambda = lambda_best$lambda.min #best lambda
  lasso_best_fit = glmnet::glmnet(x = data, y = response, alpha = 1, nlambda = 1, lambda = lambda) #fit the lasso model
  lasso_best_fit
}

# this function is to calculate correlation matrix
correlation_process = function(parameters, type) {
  
  # function of density type
  correlation_density = function(parameters) {
    corr_matrix = parameters$corr_matrix
    strong = as.integer(parameters$strong * nrow(corr_matrix))
    weak_corr = as.integer(parameters$weak_corr * nrow(corr_matrix))
    weak_no_corr = as.integer(parameters$weak_no_corr * nrow(corr_matrix))
    density_num = as.integer(min(parameters$density, 1) * weak_corr * strong)
    for (i in (strong + 1):(strong + weak_corr)) {
      changed = FALSE
      for (j in 1:strong) {
        changed = rbernoulli(1)
        if (changed) {
          corr_matrix[i, j] = parameters$corr
          corr_matrix[j, i] = parameters$corr
          density_num = density_num - 1
          break
        }
      }
      if (!changed) {
        corr_matrix[i, 1] = parameters$corr
        corr_matrix[1, i] = parameters$corr
        density_num = density_num - 1
      }
    }
    changed_rate = max(density_num / (strong * weak_corr - weak_corr), 0)
    for (i in (strong + 1):(strong + weak_corr)) {
      for (j in 1:strong) {
        if (corr_matrix[i, j] == 0.5) {
          next
        }
        changed = rbernoulli(1, changed_rate)
        if (changed) {
          corr_matrix[i, j] = parameters$corr
          corr_matrix[j, i] = parameters$corr
          density_num = density_num - 1
        }
      }
    }
    return(corr_matrix)
  }
  
  # function of unique correlation.
  correlation_unique = function(parameters) {
    corr_matrix = parameters$corr_matrix
    strong = as.integer(parameters$strong * nrow(corr_matrix))
    weak_corr = as.integer(parameters$weak_corr * nrow(corr_matrix))
    weak_no_corr = as.integer(parameters$weak_no_corr * nrow(corr_matrix))
    for (i in (strong + 1):(strong + weak_corr)) {
      changed = FALSE
      for (j in 1:strong) {
        changed = rbernoulli(1)
        if (changed) {
          corr_matrix[i, j] = parameters$corr
          corr_matrix[j, i] = parameters$corr
          break
        }
      }
      if (!changed) {
        corr_matrix[i, 1] = parameters$corr
        corr_matrix[1, i] = parameters$corr
      }
    }
    return(corr_matrix)
  }
  
  
  correlation_definite = function(parameters) {
    corr_matrix = parameters$corr_matrix
    strong = as.integer(parameters$strong * nrow(corr_matrix))
    weak_corr = as.integer(parameters$weak_corr * nrow(corr_matrix))
    weak_no_corr = as.integer(parameters$weak_no_corr * nrow(corr_matrix))
    if (strong < weak_corr) {
      stop("#strong must be greater than #weak_corr")
    }
    for (i in 1:weak_corr) {
      corr_matrix[i + strong, i] = parameters$corr
      corr_matrix[i, i + strong] = parameters$corr
    }
    return(corr_matrix)
  }
  
  if (type == "density") {
    return(correlation_density(parameters))
  }
  
  if (type == "unique") {
    return(correlation_unique(parameters))
  }
  
  if (type == "definite") {
    return(correlation_definite(parameters))
  }
}




# This is the main function of simulation
simulation <- function(time, n, p, percentage, c = 1, mu = 1, sigma = 0.5, beta = NULL, seed = 100, error_sigma = 3, correlation = 0.5, type = "density", weak_corr_density = 0.5) {
  set.seed(seed) 
  
  # check input
  if (length(percentage) != 3) {
    stop("please send the 3 samples percentage")
  }
  if (p > n) {
    warning("#Features are more than #samples")
  }
  if (length(beta) > 0 && length(beta) != p) {
    stop("#Features is not consistent with #coefficient")
  }
  if (sum(percentage < 0)) {
    stop("percentage cannot be negative")
  }
  if (weak_corr_density < 0 || weak_corr_density > 1) {
    stop("weak_corr_density must be in [0, 1]")
  }
  
  # scale percentage
  strong = percentage[1] / sum(percentage)
  weak_corr = percentage[2] / sum(percentage)
  weak_no_corr = percentage[3] / sum(percentage)
  corr_feature = as.integer((strong + weak_corr) * p)
  no_corr_feature = p - corr_feature
  
  num_corr_feature = as.integer(weak_corr_density * weak_corr * strong * p * p)
  
  sigma_matrix = matrix(rep(0, p * p), nrow = p, ncol = p)
  diag(sigma_matrix) = 1
  
  parameters = list(corr_matrix = sigma_matrix, strong = strong, weak_corr = weak_corr, weak_no_corr = weak_no_corr, density = weak_corr_density, corr = correlation)
  
  sigma_matrix = correlation_process(parameters, type)
  
  #for (i in 1:as.integer(strong * p)) {
  #  sigma_matrix[i][i + as.integer(strong * p)] = 0.5
  #  sigma_matrix[i + as.integer(strong * p)][i] = 0.5
  #}
  #print(sigma_matrix)
  
  # simulate correlated data from multivariate Guassian distribution
  #print(correlation_process(parameters, "density"))
  #return(correlation_process(parameters, "density"))
  #print(sigma_matrix * sigma)
  #return(sigma_matrix)
  
  
  
  # if beta is specified,  calculate c
  if (length(beta) != 0) {
    strong_beta = sort(beta, decreasing = TRUE)[as.integer(strong * p) + 1]
    c = strong_beta / sqrt(log(p) / n)
  }
  
  # if beta is not specified, simulate beta by c
  if (length(beta) == 0) {
    strong_beta = runif(as.integer(strong * 0.5* p), min = c * sqrt(log(p) / n), max = 5*c * sqrt(p / n))
    
    weak_corr_beta = runif(as.integer(weak_corr *0.5* p), min = 0, max = c * sqrt(log(p) / n))
    
    weak_beta = runif(0.5*p - as.integer(strong * 0.5*p) - as.integer(weak_corr * 0.5* p), min = 0, max = c * sqrt(log(p) / n))
    
    null_beta = rep(0, p - length(strong_beta)- length(weak_corr_beta) - length(weak_beta))
    
    beta = c(strong_beta, weak_corr_beta, weak_beta, null_beta)
    
  }
  
  # calculate response variable
  beta_list = list()
  measures = list()
  fitted = list()
  
  
  for (j in 1:time) {
    set.seed(time + j)
    print(j)
    data = mvrnorm(n = n, mu = rep(mu, p), Sigma = sigma_matrix * sigma)
    
    
    
    # set the variable names
    data_columns_names = paste("V", 1:p, sep = "")
    colnames(data) <- data_columns_names
    
    
    error_term = rnorm(n, mean = 0, sd = error_sigma)
    response = data %*% beta + error_term
    # regression
    lr = lm(response ~ data)
    step_lr = step(lm(response ~ 1, data = data.frame(data)), 
                   scope = formula(lm(response ~ ., data = data.frame(data))), 
                   k = 2,
                   trace = 0,
                   direction = "forward")
    lasso = lasso_function(data, response)
    
    # get fitted values
    step_fitted = fitted.values(step_lr)
    lasso_fitted = predict(lasso, data) %>% as.vector()
    fitted_values = tibble(idx = j,
                           true_y = as.vector(response), 
                           step_fitted = step_fitted,
                           lasso_fitted = lasso_fitted,
                           n = n, p = p,
                           correlation = sigma)
    
    # get coefficient estimates of the stepwise model 
    coef_step = coefficients(step_lr)[2:length(coefficients(step_lr))]
    coef_step_full = rep(0, p)
    names(coef_step_full) = paste0("V", 1:p)
    coef_step_full[match(names(coef_step), names(coef_step_full))] <- coef_step
    
    # collect beta as tibble
    beta_result = tibble(
      "idx" = j,
      "true_beta" = beta,
      "lr" = coefficients(lr)[2:length(coefficients(lr))],
      "step" = coef_step_full,
      "lasso" = as.vector(coefficients(lasso))[2:length(coefficients(lasso))],
      "type" = c(rep("strong", as.integer(strong * 0.5* p)), rep("weak_corr", as.integer(weak_corr *0.5* p)), rep("weak", 0.5*p - as.integer(strong * 0.5* p) - as.integer(weak_corr *0.5* p)), rep("null", length(null_beta)))
    ) %>% 
      replace(is.na(.), 0) 
    
    # tidy beta result for later calculation of recall, precision, f1 score
    model_measures = beta_result %>% 
      dplyr::select(- lr) %>% 
      pivot_longer(step:lasso, names_to  = "method", values_to = "coefficient") %>% 
      mutate(correct = ifelse((type != "null"& coefficient != 0) | 
                                (type == "null" & coefficient == 0), 1, 0 )) %>% 
      mutate(n = n, sigma = sigma, idx = j) 
    
    beta_list[[j]] = beta_result
    measures[[j]] = model_measures
    fitted[[j]] = fitted_values
  }
  
  
  # collect output
  #result = list(beta = beta_list, data = data, response = response, c = c, percentage = c(strong, weak_corr, weak_no_corr), lr_model = lr, step_model = step_lr, lasso_model = lasso, cov_matrix = sigma_matrix * sigma)
  result = list(beta = beta_list, measures = measures, fitted = fitted, c = c, percentage = c(strong, weak_corr, weak_no_corr), cov_matrix = sigma_matrix * sigma)
  return(result)
}


Simulation
================
Jungang Zou
2/4/2020

# correlation\_process function

Do not call this function out of simulation function

``` r
# this function is to calculate lasso function by cv
lasso_function = function(data, response){
  lambda_best = glmnet::cv.glmnet(x = data, y = response, alpha = 1, nfolds = 10) #cv to find best lambda
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
```

# Simulation Function:

input:

  - time(required) : \#simulation related to same parameters
  - n(required) : \#sample
  - p(required) : \#feature
  - percentage(required): vector of length 3, percentage of strong,
    weak\_correlated and weak data, must be positive, not need to be sum
    up to 1
  - c : to calculate the zone of beta, if beta is specified, this
    parameter has no use
  - mu : the mean of strong and weak\_correlated data
  - sigma : the covariance of strong and weak\_correlated data
  - beta : true beta. If beta is not specified, beta will be sampled by
    formulas of c
  - seed : set the random seed.
  - error\_sigma : the variance of error term
  - correlation : the correlation between strong variables and
    weak\_corr variables
  - type : the type of correlation, you can specify it in {“unique”,
    “density”, “definite”}
  - weak\_corr\_density : the density of correlation. The number of
    correlation will be density \* \#strong \* \# weak\_corr. Randomly
    select correlation. If type is not “density”, this parameter will be
    no use.

output:

  - beta : a tibble consist of true beta, lr beta, step beta, lasso
    beta, beta type
  - c : c is specified or calculated
  - percentage: scaled percentage
  - cov\_matrix: covariance matrix

The relationship between c and beta: If c is specified and beta is NULL,
then the function will calculate different zone for beta, and simulate
beta from the calculated zone by percentage. If beta is specified, c
will be no use. And new c will be calculate in function.

Correlation type: If type is “unique”, each weak\_corr variable will be
randomly correlate to only 1 strong variables. If type is “density”,
each weak\_corr variable will be randomly correlate to at least 1 strong
variables. All the correlation will be sum up to weak\_corr\_density \*
\#strong \* \#weak\_corr. If type is “definite”, each weak\_corr
variable will be assign to only 1 different strong variables. If
\#strong \< \#weak\_corr, some error will be occurred.

``` r
# This is the main function of simulation
simulation <- function(time, n, p, percentage, c = 1, mu = 1, sigma = 0.5, beta = NULL, seed = 100, error_sigma = 1, correlation = 0.5, type = "density", weak_corr_density = 0.5) {
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
    strong_beta = runif(as.integer(strong * p), min = c * sqrt(log(p) / n), max = c * sqrt(p / n))
    
    weak_corr_beta = runif(as.integer(weak_corr * p), min = 0, max = c * sqrt(log(p) / n))
    
    weak_beta = runif(p - as.integer(strong * p) - as.integer(weak_corr * p), min = 0, max = c * sqrt(log(p) / n))
    
    beta = c(strong_beta, weak_corr_beta, weak_beta)
    
    # randomly make some beta negative 
    beta_sign = runif(length(beta), -1, 1)
    beta_sign = beta_sign / abs(beta_sign)
    beta = beta * beta_sign
  }
  
  # calculate response variable
  beta_list = list()
  
  for (j in 1:time) {
    set.seed(time + j)
    
    data = mvrnorm(n = n, mu = rep(mu, p), Sigma = sigma_matrix * sigma, tol = 1)
  
  
  
    # set the variable names
    data_columns_names = paste("V", 1:p, sep = "")
    colnames(data) <- data_columns_names
    
    
    error_term = rnorm(n, mean = 0, sd = error_sigma)
    response = data %*% beta + error_term
    # regression
    lr = lm(response ~ data)
    step_lr = step(lm(response ~ data), direction = "forward")
    lasso = lasso_function(data, response)

  
    # collect beta as tibble
    beta_result = tibble(
        "true beta" = beta,
        "lr coefficient" = coefficients(lr)[2:length(coefficients(lr))],
        "step coefficient" = coefficients(step_lr)[2:length(coefficients(step_lr))],
        "lasso coefficient" = as.vector(coefficients(lasso))[2:length(coefficients(lasso))],
        "type" = c(rep("strong", as.integer(strong * p)), rep("weak_corr", as.integer(weak_corr * p)), rep("weak", p - as.integer(strong * p) - as.integer(weak_corr * p)))
        ) %>% 
      replace(is.na(.), 0)
    
    beta_list[[j]] = beta_result
  }
  
  
  
  
  
  # collect output
  #result = list(beta = beta_list, data = data, response = response, c = c, percentage = c(strong, weak_corr, weak_no_corr), lr_model = lr, step_model = step_lr, lasso_model = lasso, cov_matrix = sigma_matrix * sigma)
  result = list(beta = beta_list, c = c, percentage = c(strong, weak_corr, weak_no_corr), cov_matrix = sigma_matrix * sigma)
  return(result)
}





a = simulation(10, 100, 40, c(6, 3, 2))
```

    ## Start:  AIC=19.7
    ## response ~ data
    ## 
    ## Start:  AIC=46.63
    ## response ~ data
    ## 
    ## Start:  AIC=-4.09
    ## response ~ data
    ## 
    ## Start:  AIC=10.58
    ## response ~ data
    ## 
    ## Start:  AIC=47.27
    ## response ~ data
    ## 
    ## Start:  AIC=28.47
    ## response ~ data
    ## 
    ## Start:  AIC=26.03
    ## response ~ data
    ## 
    ## Start:  AIC=39.82
    ## response ~ data
    ## 
    ## Start:  AIC=27.65
    ## response ~ data
    ## 
    ## Start:  AIC=46.11
    ## response ~ data

``` r
a
```

    ## $beta
    ## $beta[[1]]
    ## # A tibble: 40 x 5
    ##    `true beta` `lr coefficient` `step coefficient` `lasso coefficien… type 
    ##          <dbl>            <dbl>              <dbl>              <dbl> <chr>
    ##  1      -0.457           0.0269             0.0269            -0.483  stro…
    ##  2       0.611           1.25               1.25               0.190  stro…
    ##  3       0.311           0.559              0.559              0      stro…
    ##  4       0.482           0.622              0.622              0.0358 stro…
    ##  5      -0.226          -0.710             -0.710             -0.272  stro…
    ##  6      -0.223           0.400              0.400              0      stro…
    ##  7       0.355           0.654              0.654              0      stro…
    ##  8       0.323           0.106              0.106              0      stro…
    ##  9       0.435           0.488              0.488              0.232  stro…
    ## 10       0.355           0.322              0.322              0.0476 stro…
    ## # … with 30 more rows
    ## 
    ## $beta[[2]]
    ## # A tibble: 40 x 5
    ##    `true beta` `lr coefficient` `step coefficient` `lasso coefficien… type 
    ##          <dbl>            <dbl>              <dbl>              <dbl> <chr>
    ##  1      -0.457          -0.980             -0.980            -0.721   stro…
    ##  2       0.611           0.278              0.278             0.262   stro…
    ##  3       0.311           0.0212             0.0212            0.0668  stro…
    ##  4       0.482           0.206              0.206             0.593   stro…
    ##  5      -0.226           0.109              0.109            -0.00308 stro…
    ##  6      -0.223          -0.426             -0.426            -0.0790  stro…
    ##  7       0.355          -0.0642            -0.0642            0.433   stro…
    ##  8       0.323           0.0763             0.0763            0.154   stro…
    ##  9       0.435           0.874              0.874             0.308   stro…
    ## 10       0.355           0.541              0.541             0       stro…
    ## # … with 30 more rows
    ## 
    ## $beta[[3]]
    ## # A tibble: 40 x 5
    ##    `true beta` `lr coefficient` `step coefficient` `lasso coefficien… type 
    ##          <dbl>            <dbl>              <dbl>              <dbl> <chr>
    ##  1      -0.457          -0.314             -0.314             -0.484  stro…
    ##  2       0.611           1.04               1.04               0.716  stro…
    ##  3       0.311           0.456              0.456              0.169  stro…
    ##  4       0.482          -0.0501            -0.0501             0      stro…
    ##  5      -0.226          -0.404             -0.404             -0.145  stro…
    ##  6      -0.223          -0.330             -0.330             -0.0330 stro…
    ##  7       0.355           0.192              0.192              0.333  stro…
    ##  8       0.323           0.518              0.518              0.360  stro…
    ##  9       0.435           0.697              0.697              0.206  stro…
    ## 10       0.355           0.330              0.330              0.226  stro…
    ## # … with 30 more rows
    ## 
    ## $beta[[4]]
    ## # A tibble: 40 x 5
    ##    `true beta` `lr coefficient` `step coefficient` `lasso coefficien… type 
    ##          <dbl>            <dbl>              <dbl>              <dbl> <chr>
    ##  1      -0.457          -0.492             -0.492              0      stro…
    ##  2       0.611           0.564              0.564              0.751  stro…
    ##  3       0.311           0.119              0.119              0.473  stro…
    ##  4       0.482          -0.0400            -0.0400             0      stro…
    ##  5      -0.226          -0.252             -0.252              0      stro…
    ##  6      -0.223          -0.450             -0.450             -0.0313 stro…
    ##  7       0.355           0.383              0.383              0.581  stro…
    ##  8       0.323           0.431              0.431              0.360  stro…
    ##  9       0.435           0.673              0.673              0.470  stro…
    ## 10       0.355           0.305              0.305              0      stro…
    ## # … with 30 more rows
    ## 
    ## $beta[[5]]
    ## # A tibble: 40 x 5
    ##    `true beta` `lr coefficient` `step coefficient` `lasso coefficien… type 
    ##          <dbl>            <dbl>              <dbl>              <dbl> <chr>
    ##  1      -0.457          -1.05              -1.05              -0.634  stro…
    ##  2       0.611           0.605              0.605              0.467  stro…
    ##  3       0.311          -0.286             -0.286              0      stro…
    ##  4       0.482           0.356              0.356              0.0202 stro…
    ##  5      -0.226           0.0152             0.0152             0      stro…
    ##  6      -0.223          -0.197             -0.197             -0.0578 stro…
    ##  7       0.355           0.626              0.626              0.238  stro…
    ##  8       0.323           0.164              0.164              0      stro…
    ##  9       0.435           0.0658             0.0658             0.0630 stro…
    ## 10       0.355           0.0867             0.0867             0.0541 stro…
    ## # … with 30 more rows
    ## 
    ## $beta[[6]]
    ## # A tibble: 40 x 5
    ##    `true beta` `lr coefficient` `step coefficient` `lasso coefficien… type 
    ##          <dbl>            <dbl>              <dbl>              <dbl> <chr>
    ##  1      -0.457          -0.0154            -0.0154            -0.0798 stro…
    ##  2       0.611           0.941              0.941              0.540  stro…
    ##  3       0.311           0.365              0.365              0.289  stro…
    ##  4       0.482           0.276              0.276              0      stro…
    ##  5      -0.226          -0.288             -0.288              0      stro…
    ##  6      -0.223          -0.122             -0.122             -0.129  stro…
    ##  7       0.355           0.573              0.573              0.209  stro…
    ##  8       0.323           0.103              0.103              0      stro…
    ##  9       0.435           0.598              0.598              0.454  stro…
    ## 10       0.355           0.186              0.186              0      stro…
    ## # … with 30 more rows
    ## 
    ## $beta[[7]]
    ## # A tibble: 40 x 5
    ##    `true beta` `lr coefficient` `step coefficient` `lasso coefficien… type 
    ##          <dbl>            <dbl>              <dbl>              <dbl> <chr>
    ##  1      -0.457           0.771              0.771              0      stro…
    ##  2       0.611           1.68               1.68               0.130  stro…
    ##  3       0.311           1.22               1.22               0.454  stro…
    ##  4       0.482           1.34               1.34               0.274  stro…
    ##  5      -0.226          -0.727             -0.727              0      stro…
    ##  6      -0.223           0.852              0.852              0      stro…
    ##  7       0.355           1.62               1.62               0.0214 stro…
    ##  8       0.323           0.635              0.635              0      stro…
    ##  9       0.435           0.0562             0.0562             0.364  stro…
    ## 10       0.355           0.0746             0.0746             0      stro…
    ## # … with 30 more rows
    ## 
    ## $beta[[8]]
    ## # A tibble: 40 x 5
    ##    `true beta` `lr coefficient` `step coefficient` `lasso coefficien… type 
    ##          <dbl>            <dbl>              <dbl>              <dbl> <chr>
    ##  1      -0.457          -0.988             -0.988             -0.575  stro…
    ##  2       0.611           0.532              0.532              0.645  stro…
    ##  3       0.311           0.335              0.335              0.464  stro…
    ##  4       0.482           0.249              0.249              0.308  stro…
    ##  5      -0.226          -0.229             -0.229             -0.0727 stro…
    ##  6      -0.223          -0.513             -0.513             -0.0955 stro…
    ##  7       0.355           0.0421             0.0421             0.421  stro…
    ##  8       0.323           0.498              0.498              0.468  stro…
    ##  9       0.435           0.691              0.691              0.235  stro…
    ## 10       0.355           0.601              0.601              0.0789 stro…
    ## # … with 30 more rows
    ## 
    ## $beta[[9]]
    ## # A tibble: 40 x 5
    ##    `true beta` `lr coefficient` `step coefficient` `lasso coefficien… type 
    ##          <dbl>            <dbl>              <dbl>              <dbl> <chr>
    ##  1      -0.457           0.0104             0.0104           -0.264   stro…
    ##  2       0.611           0.877              0.877             0.319   stro…
    ##  3       0.311           0.991              0.991             0.762   stro…
    ##  4       0.482           0.883              0.883             0.147   stro…
    ##  5      -0.226          -0.622             -0.622            -0.00861 stro…
    ##  6      -0.223          -0.208             -0.208            -0.415   stro…
    ##  7       0.355           0.874              0.874             0       stro…
    ##  8       0.323           0.310              0.310             0       stro…
    ##  9       0.435           0.217              0.217             0.321   stro…
    ## 10       0.355          -0.164             -0.164             0.0768  stro…
    ## # … with 30 more rows
    ## 
    ## $beta[[10]]
    ## # A tibble: 40 x 5
    ##    `true beta` `lr coefficient` `step coefficient` `lasso coefficien… type 
    ##          <dbl>            <dbl>              <dbl>              <dbl> <chr>
    ##  1      -0.457           -0.223             -0.223              0     stro…
    ##  2       0.611            0.493              0.493              0.511 stro…
    ##  3       0.311           -0.194             -0.194              0     stro…
    ##  4       0.482            0.334              0.334              0.148 stro…
    ##  5      -0.226           -0.266             -0.266              0     stro…
    ##  6      -0.223           -0.486             -0.486             -0.245 stro…
    ##  7       0.355            0.386              0.386              0.125 stro…
    ##  8       0.323            0.336              0.336              0     stro…
    ##  9       0.435            0.559              0.559              0.498 stro…
    ## 10       0.355            0.287              0.287              0.304 stro…
    ## # … with 30 more rows
    ## 
    ## 
    ## $c
    ## [1] 1
    ## 
    ## $percentage
    ## [1] 0.5454545 0.2727273 0.1818182
    ## 
    ## $cov_matrix
    ##       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
    ##  [1,] 0.50 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00  0.00  0.00  0.00  0.00
    ##  [2,] 0.00 0.50 0.00 0.00 0.00 0.00 0.00 0.00 0.00  0.00  0.00  0.00  0.00
    ##  [3,] 0.00 0.00 0.50 0.00 0.00 0.00 0.00 0.00 0.00  0.00  0.00  0.00  0.00
    ##  [4,] 0.00 0.00 0.00 0.50 0.00 0.00 0.00 0.00 0.00  0.00  0.00  0.00  0.00
    ##  [5,] 0.00 0.00 0.00 0.00 0.50 0.00 0.00 0.00 0.00  0.00  0.00  0.00  0.00
    ##  [6,] 0.00 0.00 0.00 0.00 0.00 0.50 0.00 0.00 0.00  0.00  0.00  0.00  0.00
    ##  [7,] 0.00 0.00 0.00 0.00 0.00 0.00 0.50 0.00 0.00  0.00  0.00  0.00  0.00
    ##  [8,] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.50 0.00  0.00  0.00  0.00  0.00
    ##  [9,] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.50  0.00  0.00  0.00  0.00
    ## [10,] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00  0.50  0.00  0.00  0.00
    ## [11,] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00  0.00  0.50  0.00  0.00
    ## [12,] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00  0.00  0.00  0.50  0.00
    ## [13,] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00  0.00  0.00  0.00  0.50
    ## [14,] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00  0.00  0.00  0.00  0.00
    ## [15,] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00  0.00  0.00  0.00  0.00
    ## [16,] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00  0.00  0.00  0.00  0.00
    ## [17,] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00  0.00  0.00  0.00  0.00
    ## [18,] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00  0.00  0.00  0.00  0.00
    ## [19,] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00  0.00  0.00  0.00  0.00
    ## [20,] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00  0.00  0.00  0.00  0.00
    ## [21,] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00  0.00  0.00  0.00  0.00
    ## [22,] 0.25 0.25 0.25 0.00 0.00 0.25 0.25 0.25 0.00  0.00  0.25  0.00  0.25
    ## [23,] 0.25 0.25 0.25 0.25 0.00 0.25 0.25 0.00 0.00  0.00  0.00  0.00  0.00
    ## [24,] 0.25 0.25 0.25 0.00 0.00 0.00 0.00 0.00 0.25  0.00  0.00  0.25  0.25
    ## [25,] 0.25 0.25 0.25 0.00 0.25 0.25 0.00 0.00 0.25  0.25  0.00  0.00  0.00
    ## [26,] 0.25 0.00 0.00 0.25 0.25 0.25 0.25 0.00 0.25  0.25  0.25  0.00  0.00
    ## [27,] 0.00 0.25 0.25 0.00 0.00 0.25 0.25 0.25 0.00  0.00  0.25  0.25  0.00
    ## [28,] 0.25 0.00 0.25 0.00 0.25 0.25 0.00 0.00 0.00  0.00  0.25  0.00  0.00
    ## [29,] 0.00 0.00 0.00 0.25 0.00 0.00 0.25 0.25 0.25  0.00  0.00  0.25  0.25
    ## [30,] 0.25 0.25 0.25 0.00 0.00 0.00 0.00 0.00 0.25  0.25  0.00  0.00  0.25
    ## [31,] 0.25 0.25 0.25 0.00 0.25 0.25 0.25 0.00 0.00  0.00  0.25  0.00  0.00
    ## [32,] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00  0.00  0.00  0.00  0.00
    ## [33,] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00  0.00  0.00  0.00  0.00
    ## [34,] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00  0.00  0.00  0.00  0.00
    ## [35,] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00  0.00  0.00  0.00  0.00
    ## [36,] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00  0.00  0.00  0.00  0.00
    ## [37,] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00  0.00  0.00  0.00  0.00
    ## [38,] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00  0.00  0.00  0.00  0.00
    ## [39,] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00  0.00  0.00  0.00  0.00
    ## [40,] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00  0.00  0.00  0.00  0.00
    ##       [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24]
    ##  [1,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.25  0.25  0.25
    ##  [2,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.25  0.25  0.25
    ##  [3,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.25  0.25  0.25
    ##  [4,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.25  0.00
    ##  [5,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
    ##  [6,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.25  0.25  0.00
    ##  [7,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.25  0.25  0.00
    ##  [8,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.25  0.00  0.00
    ##  [9,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.25
    ## [10,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
    ## [11,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.25  0.00  0.00
    ## [12,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.25
    ## [13,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.25  0.00  0.25
    ## [14,]  0.50  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.25  0.25  0.25
    ## [15,]  0.00  0.50  0.00  0.00  0.00  0.00  0.00  0.00  0.25  0.00  0.25
    ## [16,]  0.00  0.00  0.50  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.25
    ## [17,]  0.00  0.00  0.00  0.50  0.00  0.00  0.00  0.00  0.25  0.00  0.25
    ## [18,]  0.00  0.00  0.00  0.00  0.50  0.00  0.00  0.00  0.25  0.25  0.25
    ## [19,]  0.00  0.00  0.00  0.00  0.00  0.50  0.00  0.00  0.00  0.00  0.00
    ## [20,]  0.00  0.00  0.00  0.00  0.00  0.00  0.50  0.00  0.00  0.00  0.00
    ## [21,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.50  0.25  0.25  0.25
    ## [22,]  0.25  0.25  0.00  0.25  0.25  0.00  0.00  0.25  0.50  0.00  0.00
    ## [23,]  0.25  0.00  0.00  0.00  0.25  0.00  0.00  0.25  0.00  0.50  0.00
    ## [24,]  0.25  0.25  0.25  0.25  0.25  0.00  0.00  0.25  0.00  0.00  0.50
    ## [25,]  0.25  0.00  0.00  0.00  0.00  0.25  0.00  0.00  0.00  0.00  0.00
    ## [26,]  0.25  0.00  0.25  0.25  0.00  0.25  0.00  0.00  0.00  0.00  0.00
    ## [27,]  0.25  0.00  0.00  0.00  0.00  0.25  0.25  0.00  0.00  0.00  0.00
    ## [28,]  0.25  0.25  0.25  0.25  0.00  0.25  0.25  0.25  0.00  0.00  0.00
    ## [29,]  0.25  0.00  0.25  0.25  0.00  0.25  0.00  0.25  0.00  0.00  0.00
    ## [30,]  0.25  0.00  0.00  0.25  0.00  0.25  0.00  0.25  0.00  0.00  0.00
    ## [31,]  0.00  0.25  0.25  0.25  0.00  0.00  0.25  0.00  0.00  0.00  0.00
    ## [32,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
    ## [33,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
    ## [34,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
    ## [35,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
    ## [36,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
    ## [37,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
    ## [38,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
    ## [39,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
    ## [40,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
    ##       [,25] [,26] [,27] [,28] [,29] [,30] [,31] [,32] [,33] [,34] [,35]
    ##  [1,]  0.25  0.25  0.00  0.25  0.00  0.25  0.25   0.0   0.0   0.0   0.0
    ##  [2,]  0.25  0.00  0.25  0.00  0.00  0.25  0.25   0.0   0.0   0.0   0.0
    ##  [3,]  0.25  0.00  0.25  0.25  0.00  0.25  0.25   0.0   0.0   0.0   0.0
    ##  [4,]  0.00  0.25  0.00  0.00  0.25  0.00  0.00   0.0   0.0   0.0   0.0
    ##  [5,]  0.25  0.25  0.00  0.25  0.00  0.00  0.25   0.0   0.0   0.0   0.0
    ##  [6,]  0.25  0.25  0.25  0.25  0.00  0.00  0.25   0.0   0.0   0.0   0.0
    ##  [7,]  0.00  0.25  0.25  0.00  0.25  0.00  0.25   0.0   0.0   0.0   0.0
    ##  [8,]  0.00  0.00  0.25  0.00  0.25  0.00  0.00   0.0   0.0   0.0   0.0
    ##  [9,]  0.25  0.25  0.00  0.00  0.25  0.25  0.00   0.0   0.0   0.0   0.0
    ## [10,]  0.25  0.25  0.00  0.00  0.00  0.25  0.00   0.0   0.0   0.0   0.0
    ## [11,]  0.00  0.25  0.25  0.25  0.00  0.00  0.25   0.0   0.0   0.0   0.0
    ## [12,]  0.00  0.00  0.25  0.00  0.25  0.00  0.00   0.0   0.0   0.0   0.0
    ## [13,]  0.00  0.00  0.00  0.00  0.25  0.25  0.00   0.0   0.0   0.0   0.0
    ## [14,]  0.25  0.25  0.25  0.25  0.25  0.25  0.00   0.0   0.0   0.0   0.0
    ## [15,]  0.00  0.00  0.00  0.25  0.00  0.00  0.25   0.0   0.0   0.0   0.0
    ## [16,]  0.00  0.25  0.00  0.25  0.25  0.00  0.25   0.0   0.0   0.0   0.0
    ## [17,]  0.00  0.25  0.00  0.25  0.25  0.25  0.25   0.0   0.0   0.0   0.0
    ## [18,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00   0.0   0.0   0.0   0.0
    ## [19,]  0.25  0.25  0.25  0.25  0.25  0.25  0.00   0.0   0.0   0.0   0.0
    ## [20,]  0.00  0.00  0.25  0.25  0.00  0.00  0.25   0.0   0.0   0.0   0.0
    ## [21,]  0.00  0.00  0.00  0.25  0.25  0.25  0.00   0.0   0.0   0.0   0.0
    ## [22,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00   0.0   0.0   0.0   0.0
    ## [23,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00   0.0   0.0   0.0   0.0
    ## [24,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00   0.0   0.0   0.0   0.0
    ## [25,]  0.50  0.00  0.00  0.00  0.00  0.00  0.00   0.0   0.0   0.0   0.0
    ## [26,]  0.00  0.50  0.00  0.00  0.00  0.00  0.00   0.0   0.0   0.0   0.0
    ## [27,]  0.00  0.00  0.50  0.00  0.00  0.00  0.00   0.0   0.0   0.0   0.0
    ## [28,]  0.00  0.00  0.00  0.50  0.00  0.00  0.00   0.0   0.0   0.0   0.0
    ## [29,]  0.00  0.00  0.00  0.00  0.50  0.00  0.00   0.0   0.0   0.0   0.0
    ## [30,]  0.00  0.00  0.00  0.00  0.00  0.50  0.00   0.0   0.0   0.0   0.0
    ## [31,]  0.00  0.00  0.00  0.00  0.00  0.00  0.50   0.0   0.0   0.0   0.0
    ## [32,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00   0.5   0.0   0.0   0.0
    ## [33,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00   0.0   0.5   0.0   0.0
    ## [34,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00   0.0   0.0   0.5   0.0
    ## [35,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00   0.0   0.0   0.0   0.5
    ## [36,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00   0.0   0.0   0.0   0.0
    ## [37,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00   0.0   0.0   0.0   0.0
    ## [38,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00   0.0   0.0   0.0   0.0
    ## [39,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00   0.0   0.0   0.0   0.0
    ## [40,]  0.00  0.00  0.00  0.00  0.00  0.00  0.00   0.0   0.0   0.0   0.0
    ##       [,36] [,37] [,38] [,39] [,40]
    ##  [1,]   0.0   0.0   0.0   0.0   0.0
    ##  [2,]   0.0   0.0   0.0   0.0   0.0
    ##  [3,]   0.0   0.0   0.0   0.0   0.0
    ##  [4,]   0.0   0.0   0.0   0.0   0.0
    ##  [5,]   0.0   0.0   0.0   0.0   0.0
    ##  [6,]   0.0   0.0   0.0   0.0   0.0
    ##  [7,]   0.0   0.0   0.0   0.0   0.0
    ##  [8,]   0.0   0.0   0.0   0.0   0.0
    ##  [9,]   0.0   0.0   0.0   0.0   0.0
    ## [10,]   0.0   0.0   0.0   0.0   0.0
    ## [11,]   0.0   0.0   0.0   0.0   0.0
    ## [12,]   0.0   0.0   0.0   0.0   0.0
    ## [13,]   0.0   0.0   0.0   0.0   0.0
    ## [14,]   0.0   0.0   0.0   0.0   0.0
    ## [15,]   0.0   0.0   0.0   0.0   0.0
    ## [16,]   0.0   0.0   0.0   0.0   0.0
    ## [17,]   0.0   0.0   0.0   0.0   0.0
    ## [18,]   0.0   0.0   0.0   0.0   0.0
    ## [19,]   0.0   0.0   0.0   0.0   0.0
    ## [20,]   0.0   0.0   0.0   0.0   0.0
    ## [21,]   0.0   0.0   0.0   0.0   0.0
    ## [22,]   0.0   0.0   0.0   0.0   0.0
    ## [23,]   0.0   0.0   0.0   0.0   0.0
    ## [24,]   0.0   0.0   0.0   0.0   0.0
    ## [25,]   0.0   0.0   0.0   0.0   0.0
    ## [26,]   0.0   0.0   0.0   0.0   0.0
    ## [27,]   0.0   0.0   0.0   0.0   0.0
    ## [28,]   0.0   0.0   0.0   0.0   0.0
    ## [29,]   0.0   0.0   0.0   0.0   0.0
    ## [30,]   0.0   0.0   0.0   0.0   0.0
    ## [31,]   0.0   0.0   0.0   0.0   0.0
    ## [32,]   0.0   0.0   0.0   0.0   0.0
    ## [33,]   0.0   0.0   0.0   0.0   0.0
    ## [34,]   0.0   0.0   0.0   0.0   0.0
    ## [35,]   0.0   0.0   0.0   0.0   0.0
    ## [36,]   0.5   0.0   0.0   0.0   0.0
    ## [37,]   0.0   0.5   0.0   0.0   0.0
    ## [38,]   0.0   0.0   0.5   0.0   0.0
    ## [39,]   0.0   0.0   0.0   0.5   0.0
    ## [40,]   0.0   0.0   0.0   0.0   0.5

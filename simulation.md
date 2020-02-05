Simulation
================
Jungang Zou
2/4/2020

# Simulation Function:

input:

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

output:

  - beta : a tibble consist of true beta, lr beta, step beta, lasso
    beta, beta type
  - data : the simulated data
  - response : the simulated response data
  - c : c is specified or calculated
  - percentage: scaled percentage
  - lr\_model : the fitted linear regression model
  - step\_model : the fitted stepwise forward model
  - lasso\_model: the fitted lasso model

The relationship between c and beta: If c is specified and beta is NULL,
then the function will calculate different zone for beta, and simulate
beta from the calculated zone by percentage. If beta is specified, c
will be no use. And new c will be calculate in function.

``` r
# this function is to calculate lasso function by cv
lasso_function = function(data, response){
  lambda_best = glmnet::cv.glmnet(x = data, y = response, alpha = 1, nfolds = 10) #cv to find best lambda
  lambda = lambda_best$lambda.min #best lambda
  lasso_best_fit = glmnet::glmnet(x = data, y = response, alpha = 1, nlambda = 1, lambda = lambda) #fit the lasso model
  lasso_best_fit
}


# This is the main function of simulation
simulation <- function(n, p, percentage, c = 1, mu = 1, sigma = 0.5, beta = NULL, seed = 100, error_sigma = 1){
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
  
  # scale percentage
  strong = percentage[1] / sum(percentage)
  weak_corr = percentage[2] / sum(percentage)
  weak_no_corr = percentage[3] / sum(percentage)
  corr_feature = as.integer((strong + weak_corr) * p)
  no_corr_feature = p - corr_feature
  
  
  # simulate correlated data from multivariate Guassian distribution
  data = mvrnorm(n = n, mu = rep(mu, corr_feature), Sigma = matrix(rep(sigma, corr_feature * corr_feature), nrow = corr_feature, ncol = corr_feature))
  
  # simulate uncorrelated data
  for (i in 1:no_corr_feature) {
    no_corr = rnorm(n, mu, sigma)
    data = cbind(data, no_corr)
  }
  
  # set the variable names
  data_columns_names = paste("V", 1:p, sep = "")
  colnames(data) <- data_columns_names
  
  
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
  
  # collect output
  result = list(beta = beta_result, data = data, response = response, c = c, percentage = c(strong, weak_corr, weak_no_corr), lr_model = lr, step_model = step_lr, lasso_model = lasso)
  return(result)
}


a = simulation(10, 20, c(1,3, 2))
```

    ## Warning in simulation(10, 20, c(1, 3, 2)): #Features are more than #samples

    ## Start:  AIC=5.61
    ## response ~ data

    ## Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations
    ## per fold

``` r
a
```

    ## $beta
    ## # A tibble: 20 x 5
    ##    `true beta` `lr coefficient` `step coefficien… `lasso coefficien… type  
    ##          <dbl>            <dbl>             <dbl>              <dbl> <chr> 
    ##  1    0.991              2.37              2.37             1.09e+ 0 strong
    ##  2    0.788              0                 0                5.00e- 7 strong
    ##  3    0.860              0                 0                5.70e-15 strong
    ##  4    0.239              0                 0                0.       weak_…
    ##  5   -0.440              0                 0                0.       weak_…
    ##  6   -0.285              0                 0                0.       weak_…
    ##  7   -0.381              0                 0                0.       weak_…
    ##  8    0.464              0                 0                0.       weak_…
    ##  9   -0.463              0                 0                0.       weak_…
    ## 10   -0.214              0                 0                0.       weak_…
    ## 11    0.0840             0                 0                0.       weak_…
    ## 12   -0.350              0                 0                0.       weak_…
    ## 13   -0.157              0                 0                0.       weak_…
    ## 14   -0.508             -1.34             -1.34             0.       weak  
    ## 15   -0.0852             2.17              2.17             0.       weak  
    ## 16    0.528              0.787             0.787            0.       weak  
    ## 17    0.000690          -0.270            -0.270            0.       weak  
    ## 18    0.387              1.02              1.02             1.37e- 1 weak  
    ## 19   -0.345              3.43              3.43             0.       weak  
    ## 20    0.423             -0.0122           -0.0122           0.       weak  
    ## 
    ## $data
    ##              V1        V2        V3        V4        V5        V6
    ##  [1,] 1.3551036 1.3551036 1.3551036 1.3551036 1.3551036 1.3551036
    ##  [2,] 0.9069934 0.9069934 0.9069934 0.9069934 0.9069934 0.9069934
    ##  [3,] 1.0558028 1.0558028 1.0558028 1.0558028 1.0558028 1.0558028
    ##  [4,] 0.3729485 0.3729484 0.3729484 0.3729484 0.3729484 0.3729484
    ##  [5,] 0.9172888 0.9172888 0.9172888 0.9172888 0.9172888 0.9172888
    ##  [6,] 0.7746945 0.7746945 0.7746945 0.7746945 0.7746945 0.7746945
    ##  [7,] 1.4113881 1.4113881 1.4113881 1.4113881 1.4113881 1.4113881
    ##  [8,] 0.4947491 0.4947491 0.4947491 0.4947491 0.4947491 0.4947491
    ##  [9,] 1.5835465 1.5835465 1.5835465 1.5835465 1.5835465 1.5835465
    ## [10,] 1.2544611 1.2544609 1.2544609 1.2544609 1.2544609 1.2544609
    ##              V7        V8        V9       V10       V11       V12
    ##  [1,] 1.3551036 1.3551036 1.3551036 1.3551036 1.3551036 1.3551036
    ##  [2,] 0.9069934 0.9069934 0.9069934 0.9069934 0.9069934 0.9069934
    ##  [3,] 1.0558028 1.0558028 1.0558028 1.0558028 1.0558028 1.0558028
    ##  [4,] 0.3729484 0.3729484 0.3729484 0.3729484 0.3729484 0.3729484
    ##  [5,] 0.9172888 0.9172888 0.9172888 0.9172888 0.9172888 0.9172888
    ##  [6,] 0.7746945 0.7746945 0.7746945 0.7746945 0.7746945 0.7746945
    ##  [7,] 1.4113881 1.4113881 1.4113881 1.4113881 1.4113881 1.4113881
    ##  [8,] 0.4947491 0.4947491 0.4947491 0.4947491 0.4947491 0.4947491
    ##  [9,] 1.5835465 1.5835465 1.5835465 1.5835465 1.5835465 1.5835465
    ## [10,] 1.2544609 1.2544609 1.2544609 1.2544609 1.2544609 1.2544609
    ##             V13       V14         V15       V16       V17         V18
    ##  [1,] 1.3551036 0.7911028  0.79578749 1.5173432 0.8788653  1.31603703
    ##  [2,] 0.9069934 0.5748096 -0.06824693 1.8267516 1.0295157  1.10070676
    ##  [3,] 1.0558028 1.3445231  1.07841096 0.9910266 0.9113641  0.95446468
    ##  [4,] 0.3729484 0.7699019  1.33002445 0.9878983 1.3973401  1.14474206
    ##  [5,] 0.9172888 1.6740922  0.50908279 1.1251235 1.0033689  0.97265753
    ##  [6,] 0.7746945 1.2215357  0.44317815 0.8314377 0.6851049 -0.02092493
    ##  [7,] 1.4113881 0.9245369  0.78132616 0.9433231 0.8737551  1.17918462
    ##  [8,] 0.4947491 1.2277744  0.74194438 0.9505585 0.6547889  0.81369957
    ##  [9,] 1.5835465 0.9799227  1.20949800 1.1320434 1.1012711  1.63415442
    ## [10,] 1.2544609 1.2280605  1.06707772 1.0694918 1.4231907  2.08430016
    ##             V19       V20
    ##  [1,] 0.3801386 0.8834967
    ##  [2,] 1.2949369 0.8745916
    ##  [3,] 1.0620096 1.4769477
    ##  [4,] 0.7381461 0.8670137
    ##  [5,] 1.3101140 1.9476380
    ##  [6,] 1.3541108 0.7850046
    ##  [7,] 0.9534008 1.7877735
    ##  [8,] 0.8524016 1.0809706
    ##  [9,] 0.4570924 0.4572735
    ## [10,] 0.6875925 1.2884687
    ## 
    ## $response
    ##             [,1]
    ##  [1,]  1.1052821
    ##  [2,]  2.0343344
    ##  [3,]  1.6027811
    ##  [4,]  1.1605670
    ##  [5,]  1.9012464
    ##  [6,] -0.1095985
    ##  [7,]  3.0865284
    ##  [8,] -0.3465144
    ##  [9,]  3.8112358
    ## [10,]  1.7918566
    ## 
    ## $c
    ## [1] 1
    ## 
    ## $percentage
    ## [1] 0.1666667 0.5000000 0.3333333
    ## 
    ## $lr_model
    ## 
    ## Call:
    ## lm(formula = response ~ data)
    ## 
    ## Coefficients:
    ## (Intercept)       dataV1       dataV2       dataV3       dataV4  
    ##    -5.93592      2.37001           NA           NA           NA  
    ##      dataV5       dataV6       dataV7       dataV8       dataV9  
    ##          NA           NA           NA           NA           NA  
    ##     dataV10      dataV11      dataV12      dataV13      dataV14  
    ##          NA           NA           NA           NA     -1.33790  
    ##     dataV15      dataV16      dataV17      dataV18      dataV19  
    ##     2.16608      0.78656     -0.26978      1.02013      3.42670  
    ##     dataV20  
    ##    -0.01216  
    ## 
    ## 
    ## $step_model
    ## 
    ## Call:
    ## lm(formula = response ~ data)
    ## 
    ## Coefficients:
    ## (Intercept)       dataV1       dataV2       dataV3       dataV4  
    ##    -5.93592      2.37001           NA           NA           NA  
    ##      dataV5       dataV6       dataV7       dataV8       dataV9  
    ##          NA           NA           NA           NA           NA  
    ##     dataV10      dataV11      dataV12      dataV13      dataV14  
    ##          NA           NA           NA           NA     -1.33790  
    ##     dataV15      dataV16      dataV17      dataV18      dataV19  
    ##     2.16608      0.78656     -0.26978      1.02013      3.42670  
    ##     dataV20  
    ##    -0.01216  
    ## 
    ## 
    ## $lasso_model
    ## 
    ## Call:  glmnet::glmnet(x = data, y = response, alpha = 1, nlambda = 1,      lambda = lambda) 
    ## 
    ##   Df   %Dev Lambda
    ## 1  4 0.4325 0.4434

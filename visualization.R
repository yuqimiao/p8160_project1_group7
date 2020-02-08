# Simulation result analysis

library(tidyverse)

# read in saved simualtion results
beta_results = readRDS(file = "beta_estimates.rds")
fitted_results = readRDS(file = "fitted_values.rds")
measures_result =  readRDS(file = "measurement_scores.rds")


# Extract results of simulation setting 1:
##1- #100: sample size = 200, sigma(correlation) = 0.3
beta_df_n200_cor03 <-
  bind_rows(beta_results[1:100]) %>% 
  rename(idx_sim = idx) %>% 
  mutate(bias2 = (coefficient - true_beta)^2,
         bias = coefficient - true_beta)
beta_df_n200_cor03$idx_beta = rep(1:100, each = 2)
beta_df_n200_cor03 = beta_df_n200_cor03 %>% 
  select(idx_sim, idx_beta, everything())

## explore data
names(beta_df_n200_cor03)
unique(beta_df_n200_cor03$idx_sim)
unique(beta_df_n200_cor03$true_beta) %>% sort()
unique(beta_df_n200_cor03$type)

## number of different types of predictors
table(beta_results[[1]]$type, beta_results[[1]]$method)

## plots of beta estimates
beta_df_n200_cor03 %>% 
  group_by(idx_sim, method, type) %>% 
  summarise(beta_rmse = sqrt(mean(bias2)),
            beta_bias_mean = mean(bias)) %>% 
  ggplot(aes(x = method, y = beta_rmse, color = type)) +
    geom_boxplot()

beta_df_n200_cor03 %>% 
  ggplot(aes(x = method, y = bias, color = type)) +
    geom_boxplot()

beta_df_n200_cor03 %>% 
  ggplot(aes(x = method, y = bias)) +
    geom_boxplot()


# add simulation index 1-100 to fitted_results
for (i in 1:length(fitted_results)) {
  fitted_results[[i]]$idx_sim = i %% 100  
}

# predicted y 
y_df_n200_cor03 <-
  bind_rows(fitted_results[1:100]) %>%  
  gather(key = "method", value = "y_pred", step_fitted:lasso_fitted) %>% 
  mutate(loss_l2 = (y_pred - true_y)^2) %>% 
  select(idx_sim, everything()) 

y_df_n200_cor03 %>% 
  group_by(idx_sim, method) %>% 
  summarize(rmse = mean(loss_l2)) %>% 
  ggplot(aes(x = method, y = rmse)) +
    geom_boxplot()

y_df <-
  bind_rows(fitted_results) %>% 
  gather(key = "method", value = "y_pred", step_fitted:lasso_fitted) %>% 
  mutate(loss_l2 = (y_pred - true_y)^2,
         n = as.factor(n)) %>% 
  group_by(n, correlation, idx_sim, method) %>% 
  summarize(rmse = sqrt(mean(loss_l2)))

# RMSE of the outcome for all the simulation settings
y_df %>%  
  ggplot(aes(x = method, y = rmse, color = n)) +
    geom_boxplot() +
    facet_grid(~correlation)
  
  

# 
  

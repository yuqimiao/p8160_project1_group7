source("source.R")


# Save simulation results and fitted values
beta_results <- list()
fitted_results <- list()
sample_size = 5000
  temp_results <- simulation(time = 100, n = sample_size, p = 200, percentage = c(5, 1, 4), 
                             sigma = 1, type = "definite")
saveRDS(temp_results, file = paste("result_less_weak_corr_n", sample_size, ".rds", sep = ""))


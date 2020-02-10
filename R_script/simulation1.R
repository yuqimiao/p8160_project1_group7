source("source.R")


# Save simulation results and fitted values
beta_results <- list()
fitted_results <- list()
sample_size = 200
  temp_results <- simulation(time = 100, n = sample_size, p = 200, percentage = c(5, 3, 2), 
                             sigma = 0.5, type = "definite")
saveRDS(temp_results, file = paste("result_more_weak_corr_n", sample_size, ".rds", sep = ""))


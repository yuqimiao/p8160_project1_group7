source("source2.R")


# Save simulation results and fitted values
beta_results <- list()
fitted_results <- list()
sample_size = 1000
  temp_results <- simulation(time = 100, n = sample_size, p = 200, percentage = c(5, 3, 2), 
                             sigma = 1, type = "definite")
saveRDS(temp_results, file = paste("result_missing_weak_n", sample_size, ".rds", sep = ""))


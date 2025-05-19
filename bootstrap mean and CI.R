

bootstrap_mean_ci <- function(data, cat_var, num_var, W = "weight.adjusted",
                              nIter = 5000, confprob = 0.05) {
  
  # Get unique values of the categorical variable
  cat_levels <- unique(data[[cat_var]])
  
  # Create an empty data frame to store results
  results <- data.frame(category = character(),
                        mean = numeric(),
                        lower_ci = numeric(),
                        upper_ci = numeric(),
                        n = numeric())
  
  for (cat in cat_levels) {
    # Subset data for the current category
    cat_data <- data[data[[cat_var]] == cat, ]
    
    # Extract the numerical and weight variables
    values <- cat_data[[num_var]]
    weights <- cat_data[[W]]
    
    # record sample size
    n_cat <- nrow(cat_data)
    
    # Calculate weighted mean
    weighted_mean <- sum(values * weights) / sum(weights)
    
    # Perform bootstrapping
    boot_means <- replicate(nIter, {
      # Sample indices with replacement
      sample_indices <- sample(seq_along(values), replace = TRUE)
      # Calculate weighted mean for this bootstrap sample
      sum(values[sample_indices] * weights[sample_indices]) / sum(weights[sample_indices])
    })
    
    # Confidence interval
    ci <- quantile(boot_means, c(confprob / 2, 1 - confprob / 2))
    
    
    # Store the results
    results <- rbind(results, data.frame(category = cat,
                                         mean = mean(boot_means),
                                         lower_ci = ci[1],
                                         upper_ci = ci[2],
                                         n = n_cat))
  }
  
  return(results)
}







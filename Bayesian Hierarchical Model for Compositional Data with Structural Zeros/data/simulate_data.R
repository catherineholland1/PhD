###############################.
## SIMULATE DATA FOR JOURNAL ##
###############################.

## PACKAGES ####

library(dplyr)
library(readr)

##################################.

set.seed(123)

setwd("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/Statistical Modelling Journal/Forensic/Data")

##################################.

## PARAMETERS ##

n_items <- 100
n_rep <- 4
types <- 1:5

n_elements <- 10

# Create a data frame to hold results
simulated_data <- data.frame()

# Simulate function
simulate_compositions <- function(n = 1, d, zero_prob = 0.3) {
  
  # n = number of samples
  # d = number of components
  # zero_prob = probability a component is zero
  
  compositions <- matrix(0, nrow = n, ncol = d)
  
  for (i in 1:n) {
    vals <- runif(d, min = 0, max = 1)
    # Randomly set some components to zero
    zero_indices <- runif(d) < zero_prob
    vals[zero_indices] <- 0
    
    # If all zeros, resample to avoid division by zero
    while (sum(vals) == 0) {
      vals <- runif(d, min = 0, max = 1)
      zero_indices <- runif(d) < zero_prob
      vals[zero_indices] <- 0
    }
    
    # Normalize so sum = 1
    compositions[i, ] <- vals / sum(vals)
  }
  
  return(compositions)
}


## SIMULATE ##
for (item_id in 1:n_items) {
  item_type <- sample(types, 1)
  
  for (rep_id in 1:n_rep) {
    # Generate one compositional vector (1 x n_elements)
    comp <- simulate_compositions(n = 1, d = n_elements, zero_prob = 0.1)
    comp_vec <- as.numeric(comp)
    
    variables <- data.frame(setNames(as.list(comp_vec), paste0("variable_", 1:n_elements)))
    
    # Add new row to simulated_data
    simulated_data <- rbind(simulated_data, data.frame(
      type = item_type,
      item = item_id,
      replicate = rep_id,
      variables
    ))
  }
}


head(simulated_data)

rowSums(simulated_data[, 4:13])

sum(simulated_data[, 4:13] == 0)


##################################.

## SAVE ####

simulated_data %>%
  write_rds(., "simulated_compositional_data.rds")

##################################.
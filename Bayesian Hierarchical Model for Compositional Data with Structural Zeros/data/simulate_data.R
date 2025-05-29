###############################.
## SIMULATE DATA FOR JOURNAL ##
###############################.

## PACKAGES ####

library(dplyr)
library(readr)
library(MASS)
library(tidyverse)

##################################.

set.seed(123)

setwd("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/Statistical Modelling Journal/Forensic/Data")

##################################.

## PARAMETERS ##

n_items <- 500
n_rep <- 4
types <- 1:5

n_components <- 5

total = 100


## SIMULATE DATA ##

item_effects <- mvrnorm(n_items, mu = rep(0, n_components - 1), Sigma = diag(0.5, n_components - 1))

## Decide which components are non-zero per item ##
# Always keep variable_1 (index 1), others randomly dropped
non_zero_components_list <- map(1:n_items, function(i) {
  keep_idx <- sample(2:n_components, sample(1:(n_components - 1), 1))  # pick others to keep
  sort(unique(c(1, keep_idx)))  # always keep variable_1
})

## SIMULATE OBSERVATIONS ##
simulated_data <- map_df(1:n_items, function(item_id) {
  keep_idx <- non_zero_components_list[[item_id]]
  
  map_df(1:n_rep, function(obs_id) {
    noise <- rnorm(n_components - 1, 0, 0.3)
    log_ratios <- item_effects[item_id, ] + noise
    log_ratios <- c(log_ratios, 0)  # last component anchored
    
    masked_log_ratios <- rep(-Inf, n_components)
    masked_log_ratios[keep_idx] <- log_ratios[keep_idx]
    
    exp_vals <- exp(masked_log_ratios)
    props <- exp_vals / sum(exp_vals)
    comps <- round(props * total, 2)
    
    tibble(
      item = item_id,
      rep = obs_id,
      variable_1 = comps[1],
      variable_2 = comps[2],
      variable_3 = comps[3],
      variable_4 = comps[4],
      variable_5 = comps[5]
    )
  })
})


simulated_data <- simulated_data %>%
  group_by(item) %>%
  mutate(type = sample(1:5, 1, replace = T), .after = 'item') %>% 
  ungroup()


## CHECK ##

head(simulated_data)

rowSums((simulated_data %>% select(starts_with('variable'))))

sum((simulated_data %>% select(starts_with('variable'))) == 0)
sum((simulated_data %>% select(starts_with('variable'))) == 0) / (nrow(simulated_data)*n_components)

apply(simulated_data %>% select(starts_with('variable')), 2, function(x) sum(x == 0))
apply(simulated_data %>% select(starts_with('variable')), 2, function(x) sum(x == 0)) / nrow(simulated_data)


##################################.

## SAVE ####

simulated_data %>%
  write_rds(., "simulated_compositional_data.rds")

##################################.
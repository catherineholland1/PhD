############################################.
## HIDDEN MARKOV MODEL - GDM - SIMULATION ##
############################################.

## PACKAGES ####

library(dplyr)
library(tidyr)
library(ggplot2)


#######################.

## LOAD IN OUTPUT ####

## SET WORKING DIRECTORY ##

setwd("GDM")


## LOAD IN SAMPLES ####
load('output/mcmc_output-parallel.RData')


#######################.

## SIMULATION ####

output <- as_tibble(do.call('rbind', samples))


N_samples <- nrow(output)


#######################.

#### SIMULATION FUNCTION ####

simulate_hmm_count <- function(N, N_states, N_variants, N_countries,
                               y, 
                               p0, p,
                               nu, phi,
                               cluster_index) {
  
  
  # INITILISE TO STORE SIMULATION
  simulated_data <- array(0, dim = c(N, N_variants, N_countries))
  
  sum_y_minus_x <- array(0, dim = c(N, N_variants, N_countries))
  
  
  ## LOOP OVER EACH COUNTRY ##
  
  for (m in 1:N_countries) {
    
    ## FIRST VARIANT ##
    
    # INITIALISE STATE AND OBSERVATION
    state <- numeric(N)
    observation <- numeric(N)
    
    # SAMPLE INTIAL STATE
    state[1] <- sample(1:N_states, 1, prob = p0[1, , cluster_index[m]])
    
    # SAMPLE REMAINING STATE SEQUENCE
    for (i in 2:N) {
      
      # PROBABILITIES FOR EACH STATE
      state_probs <- p[state[i-1], , 1, cluster_index[m]]  
      
      # SAMPLE STATE
      state[i] <- sample(1:N_states, 1, prob = state_probs)
    }
    
    # SIMULATE COUNTS
    alpha <- sapply(1:N, function(i)nu[1, state[i], m] * phi[1, state[i], cluster_index[m]])
    beta <- sapply(1:N, function(i)phi[1, state[i], cluster_index[m]]*(1 - nu[1, state[i], m]))
    beta_probabilities <- rbeta(N, alpha, beta)
    observation <- rbinom(N, y[, m], beta_probabilities)
    
    # SAVE SIMULATED VALUES
    simulated_data[, 1, m] <- observation
    
    
    ## ALL OTHER VARIANTS ##
    
    for (v in 2:N_variants) {
      
      # INITIALISE STATE AND OBSERVATION
      state <- numeric(N)
      observation <- numeric(N)

      # SAMPLE INTIAL STATE
      state[1] <- sample(1:N_states, 1, prob = p0[1, , cluster_index[m]])
      
      # COMPUTE Y - X
      sum_y_minus_x[1, (v-1), m] <- y[1, m] - sum(simulated_data[1, 1:(v-1), m])
      
      # SAMPLE REMAINING STATE SEQUENCE
      for (i in 2:N) {
        
        # PROBABILITIES FOR EACH STATE
        state_probs <- p[state[i-1], , v, cluster_index[m]]  
        
        # SAMPLE STATE
        state[i] <- sample(1:N_states, 1, prob = state_probs)
        
        # COMPUTE Y - X
        sum_y_minus_x[i, (v-1), m] <- y[i, m] - sum(simulated_data[i, 1:(v-1), m])
        
      }
      
      # SIMULATE COUNTS
      alpha <- sapply(1:N,function(i)nu[v, state[i], m] * phi[v, state[i], cluster_index[m]])
      beta <- sapply(1:N,function(i)phi[v, state[i], cluster_index[m]]*(1 - nu[v, state[i], m]))
      beta_probabilities <- rbeta(N, alpha, beta)
      #observation <- rbinom(N, sum_y_minus_x[ , v, m], beta_probabilities)
      if(v>2){
        observation <- rbinom(N, y[, m]-apply(simulated_data[, 1:(v-1), m],1,sum), 
                              beta_probabilities)
      }else{
        observation <- rbinom(N, y[, m]-simulated_data[, 1, m], beta_probabilities)
      }
      
      # SAVE SIMULATED VALUES
      simulated_data[, v, m] <- observation
      
      }
    
  }
  
  return(simulated_data)
  
}


#######################.


#### DATA ####

x_array <- x
y_array <- y


### PARAMETERS ####

nu_samples <- array(NA, dim = c(N_variants, N_states, N_countries, N_samples))
for (i in 1:N_samples) {
  nu_row <- output[i, ] %>% select(starts_with("nu")) 
  nu_samples[, , , i] <- array(unlist(nu_row), dim = c(N_variants, N_states, N_countries))
}
nu_samples


phi_samples <- array(NA, dim = c(N_variants, N_states, N_clusters, N_samples))
for (i in 1:N_samples) {
  phi_row <- output[i, ] %>% select(starts_with("phi")) 
  phi_samples[, , , i] <- array(unlist(phi_row), dim = c(N_variants, N_states, N_clusters))
}
phi_samples


p0_samples <- array(NA, dim = c(N_variants, N_states, N_clusters, N_samples))
for (i in 1:N_samples) {
  p0_row <- output[i, ] %>% select(starts_with("p0")) 
  p0_samples[, , , i] <- array(unlist(p0_row), dim = c(N_variants, N_states, N_clusters))
}
p0_samples


p_samples <- array(NA, dim = c(N_states, N_states, N_variants, N_clusters, N_samples))
for (i in 1:N_samples){
  p_row <- output[i, ] %>% select(starts_with("p[")) 
  p_samples[, , , , i] <- array(unlist(p_row), dim = c(N_states, N_states, N_variants, N_clusters))
}
p_samples



### RUN SIMULATION ####

simulation_replicates <- array(NA, dim = c(N, N_variants, N_countries, N_samples))
for (i in 1:N_samples){
  simulation_replicates[, , , i] <- simulate_hmm_count(N = N,
                                                       N_states = N_states, 
                                                       N_variants = N_variants,
                                                       N_countries = N_countries,
                                                       cluster_index = cluster_index,
                                                       y = y_array,
                                                       p0 = p0_samples[, , , i], 
                                                       p = p_samples[, , , , i],
                                                       nu = nu_samples[, , , i], 
                                                       phi = phi_samples[, , , i])
}
simulation_replicates


#######################.


## SAVE REPLICATES ####

save(simulation_replicates, file = 'output/replicates_GDM.RData')


#######################.


## SAVE OUTPUT SUMMARY ####

gamma_samples <- array(NA, dim = c(N_variants, N_states, N_countries, N_samples))
for (i in 1:N_samples) {
  gamma_row <- output[i, ] %>% select(starts_with("gamma")) 
  gamma_samples[, , , i] <- array(unlist(gamma_row), dim = c(N_variants, N_states, N_countries))
}
gamma_samples


output_matrices <- list()
output_matrices[["p"]] <- p_samples %>% apply(., c(1,2,3), median)
output_matrices[["p0"]] <- p0_samples %>% apply(., c(1,2), median)
output_matrices[["gamma"]] <- gamma_samples %>% apply(., c(1,2,3), median)
output_matrices[["nu"]] <- nu_samples %>% apply(., c(1,2,3), median)
output_matrices[["phi"]] <- phi_samples %>% apply(., c(1,2), median)


sink("output/hmm_summary_output_matrices.txt")
print(output_matrices)
sink()


#######################.


## SUMMARY STATISTICS ####


### WINDOWS ####

generate_window_seq <- function(N, window_length) {
  window_seq <- list()
  for (i in 1:(N - (window_length-1))) {
    week_seq <- seq(i, i + (window_length-1))
    window_seq[[i]] <- data.frame(week = week_seq, window = rep(i, length(week_seq)))
  }
  do.call(rbind, window_seq)
}

window_seq <- generate_window_seq(N, window_length = 10)


#######################.


### WINDOWS TO ORIGINAL DATA ####

x_data_windows <- list()
for (index in 1:length(x)) {
  x_data_windows[[index]] <- x[[index]] %>% 
    mutate(country = countries[index]) %>%
    mutate(week = row_number()) %>% 
    right_join(., window_seq, by = "week", relationship = "many-to-many")
}


x_data_windows <- do.call(rbind, x_data_windows)


### COMPUTE SD & MEAN ####

x_summary <- x_data_windows %>%
  select(-week) %>%
  group_by(country, window) %>%
  summarise(across(c(alpha:omicron), sd, .names = "{.col}_sd")) %>%
  ungroup() %>%
  group_by(country) %>%
  summarise(across(c(alpha_sd:omicron_sd), mean, .names = "mean_{.col}"))


#######################.

### SAVE DATA SUMMARY WITH WINDOWS ####

save(x_summary, file = 'output/x_data_windows.RData')


#######################.

### WINDOWS TO SIMULATION ####

all_rep_data <- list()
for (replicate_index in 1:N_samples) {
  
  rep_data <- do.call(rbind, lapply(1:N_countries, function(i) {
    slice <- simulation_replicates[, , i, replicate_index]
    slice <- cbind(slice, country = rep(countries[i], N), week = 1:N)
    return(slice)
  })) %>%
    as.data.frame() %>%
    mutate_at(vars(starts_with("V"), week), as.numeric) %>%
    rename(alpha = "V1",
           beta = "V2",
           gamma = "V3",
           delta = "V4",
           omicron = "V5") %>% 
    right_join(., window_seq, by = "week", relationship = "many-to-many") %>%
    mutate(window = as.numeric(window))
  
  all_rep_data[[replicate_index]] <- rep_data
}

#all_rep_data <- do.call(rbind, all_rep_data)


### COMPUTE SD & MEAN ####

summary <- list(length = length(all_rep_data))
for (i in 1:N_samples) {
  summary[[i]] <- all_rep_data[[i]] %>%
    group_by(window, country) %>%
    summarise(across(c(alpha:omicron), sd, .names = "{.col}_sd")) %>%
    ungroup() %>%
    group_by(country) %>%
    summarise(across(c(alpha_sd:omicron_sd), mean, .names = "mean_{.col}"))
}


combined_summary <- bind_rows(summary) 


#######################.

### SAVE SUMMARY WITH WINDOWS ####

save(combined_summary, file = 'output/GDM_windows.RData')

#######################.

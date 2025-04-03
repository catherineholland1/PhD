###########################################################.
## RANDOM WALK TIME SERIES MODEL - SIMULATION - PARALLEL ##
###########################################################.

## PACKAGES ####

library(dplyr)
library(tidyr)
library(ggplot2)


#######################.

## LOAD IN OUTPUT ####

## SET WORKING DIRECTORY ##

setwd("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD R Code/Time Series - COVID-19/RW")


## LOAD IN SAMPLES ####
load('output/rw_mcmc_output-parallel.RData')


#######################.

## SIMULATION ####

output <- as_tibble(do.call('rbind', samples))

N_samples <- nrow(output)


#######################.

#### SIMULATION FUNCTION ####

simulate_rw <- function(N, N_variants, N_countries,
                        y, 
                        lambda_mu,
                        lambda_sigma,
                        phi,
                        sigma,
                        cluster_index) {
  
  
  # INITILISE TO STORE SIMULATION
  simulated_data <- array(0, dim = c(N, N_variants, N_countries))
  
  # INITIALISE RANDOM WALK #
  lambda <- array(NA, dim = c(N, N_variants, N_countries))
  eta <- array(NA, dim = c(N, N_variants, N_countries))
  
  
  ## RANDOM WALK ##
  
  for (m in 1:N_countries) {
    for (v in 1:N_variants) {
      
      lambda[1, v, m] <- rnorm(1, lambda_mu[v, cluster_index[m]], sd = lambda_sigma[v, cluster_index[m]])
      eta[1, v, m] <- 1 / (1 + exp(-lambda[1, v, m]))
      
      for (t in 2:N) {
        
        lambda[t, v, m] <- rnorm(1, lambda[(t-1), v, m], sd = sigma[v, m]) 
        eta[t, v, m] <- 1 / (1 + exp(-lambda[t, v, m]))
        
      }
    }
  }
  
  ## LOOP OVER EACH COUNTRY ##
  
  for (m in 1:N_countries) {
    
    ## FIRST VARIANT ##
    
    # INITIALISE OBSERVATION
    observation <- numeric(N)
    
    
    # SIMULATE COUNTS
    alpha <- sapply(1:N, function(i)eta[1:N, 1, m] * phi[1, cluster_index[m]])
    beta <- sapply(1:N, function(i)phi[1, cluster_index[m]]*(1 - eta[1:N, 1, m]))
    beta_probabilities <- rbeta(N, alpha, beta)
    observation <- rbinom(N, y[, m], beta_probabilities)
    
    # SAVE SIMULATED VALUES
    simulated_data[, 1, m] <- observation
    
    
    ## ALL OTHER VARIANTS ##
    
    for (v in 2:N_variants) {
      
      # INITIALISE OBSERVATION
      observation <- numeric(N)
      
      # SIMULATE COUNTS
      alpha <- sapply(1:N,function(i)eta[1:N, v, m] * phi[v, cluster_index[m]])
      beta <- sapply(1:N,function(i)phi[v, cluster_index[m]]*(1 - eta[1:N, v, m]))
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

x_array <- array(unlist(x), dim = c(N, N_variants, N_countries))
y_array <- array(unlist(y), dim = c(N, N_countries))


#### PARAMETERS ####

lambda_mu <- matrix(unlist(output[1, ] %>% select(starts_with("lambda_mu"))), nrow = N_variants, ncol = N_clusters, byrow = F)
lambda_sigma <- matrix(unlist(output[1, ] %>% select(starts_with("lambda_sigma"))), nrow = N_variants, ncol = N_clusters, byrow = F)


phi_samples <- array(NA, dim = c(N_variants, N_clusters, N_samples))
for (i in 1:N_samples) {
  phi_row <- output[i, ] %>% select(starts_with("phi")) 
  phi_samples[, , i] <- matrix(unlist(phi_row), nrow = N_variants, ncol = N_clusters, byrow = F)
}
phi_samples


sigma_samples <- array(NA, dim = c(N_variants, N_countries, N_samples))
for (i in 1:N_samples) {
  sigma_row <- output[i, ] %>% select(starts_with("sigma")) 
  sigma_samples[, , i] <- matrix(unlist(sigma_row), nrow = N_variants, N_countries, byrow = F)
}
sigma_samples



#### RUN SIMULATION ####

simulation_replicates <- array(NA, dim = c(N, N_variants, N_countries, N_samples))
for (i in 1:N_samples){
  simulation_replicates[, , , i] <- simulate_rw(N = N,
                                                N_variants = N_variants,
                                                N_countries = N_countries,
                                                y = y_array,
                                                lambda_mu = lambda_mu,
                                                lambda_sigma = lambda_sigma,
                                                phi = phi_samples[, , i],
                                                sigma = sigma_samples[ , , i],
                                                cluster_index = cluster_index)
}
simulation_replicates


#######################.

## SAVE REPLICATES ####

save(simulation_replicates, file = 'output/replicates_RW.RData')


#######################.


## SAVE OUTPUT SUMMARY ####

output_matrices <- list()
output_matrices[["lambda_mu"]] <- lambda_mu
output_matrices[["lambda_sigma"]] <- lambda_sigma
output_matrices[["sigma"]] <- sigma_samples %>% apply(., c(1,2), median)
output_matrices[["phi"]] <- phi_samples %>% apply(., c(1,2), median)


sink("output/rw_summary_output_matrices.txt")
print(output_matrices)
sink()


#######################.

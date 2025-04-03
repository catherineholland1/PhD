###################################.
## HIDDEN MARKOV MODEL - GDM-HMM ##
###################################.

# STATE 1 - DORMANT BEFORE OUTBREAK
# STATE 2 - ACTIVE - INCREASING
# STATE 3 - PEAK
# STATE 4 - ACTIVE - DECREASING
# STATE 5 - DORMANT AFTER OUTBREAK


#######################.

## PACKAGES ####

library(dplyr)
library(tidyr)
library(nimble)
library(coda)
library(mcmcplots)
library(ggplot2)
library(reshape2)
library(abind)
library(parallel)

#######################.


set.seed(123)


## SET WORKING DIRECTORY ##

setwd("GDM")

## GET BETA-BINOMIAL FUNCTION ##

source("forward_algorithm_count_nimblefunction.R")


#######################.

## DATA ####

# VARIANTS
variants <- c("alpha", 
              "beta",
              "gamma", 
              "delta",
              "omicron")

# COUNTRIES / CLUSTERS
countries_cluster <- cbind(countries = c("United Kingdom",
                                         "France",
                                         "Italy",
                                         "Germany",
                                         "Spain",
                                         "Australia",
                                         "Brazil",
                                         "Mexico",
                                         "South Africa",
                                         "Canada"), 
                           cluster = 3) %>%
  rbind(cbind(countries = c("Armenia",
                            "Azerbaijan",
                            "China",
                            "Dominica",
                            "Estonia",
                            "Fiji",
                            "India",
                            "Kenya",
                            "Qatar",
                            "United Arab Emirates"),
              cluster = 1)) %>%
  rbind(cbind(countries = c("Argentina",
                            "Bulgaria",
                            "Dominican Republic",
                            "Greece",
                            "Jamaica",
                            "Morocco",
                            "New Zealand",
                            "Pakistan",
                            "Peru",
                            "Philippines"),
              cluster = 2)) %>%
  as.data.frame() %>%
  arrange(countries)


cluster_index <- as.numeric(countries_cluster$cluster)



variant_data <- readRDS(file = "") %>%
  select_if(~ !is.numeric(.) || sum(.) != 0) 


data <- variant_data %>%
  filter(
    # COUNTRIES
    country %in% countries_cluster$countries,
    # VARIANTS
    variant %in% variants) %>%
  arrange(country, match(variant, variants))


data_long <- data %>%
  pivot_longer(starts_with("20"), names_to = "date", values_to = "n") %>%
  mutate(date = as.Date(date))


# Y = TOTAL COUNT

y <- data %>%
  filter(measure == "total") %>% 
  pivot_longer(., cols = starts_with("20"), names_to = "date", values_to = "y") %>% 
  select(date, country, y) %>% 
  group_by(country, date) %>% 
  summarise(y = max(y)) %>% 
  ungroup() %>%
  select(-date) %>%
  split(.$country) %>%
  lapply(., function(df) {
    select(df, -country)
  })
y


# X = VARIANT COUNT

x <- data %>%
  filter(measure == "count") %>% 
  pivot_longer(., cols = starts_with("20"), names_to = "date", values_to = "x") %>% 
  pivot_wider(., names_from = variant, values_from = x) %>%
  mutate(across(alpha:omicron, ~if_else(is.na(.x), 0, .x))) %>%
  select(date, country, alpha:omicron) %>%
  split(.$country) %>%
  lapply(., function(df) {
    select(df, -c(date, country))
  })
x


# DATES IN DATA
dates <- data_long %>%
  select(date) %>%
  distinct() %>%
  t() %>%
  as.Date()
dates


###################.

## CONSTANTS ####

# NUMBER OF STATES
N_states = 5

# NUMBER OF VARIANTS
N_variants = length(unique(variants))

# NUMBER OF COUNTRIES
N_countries = length(unique(countries_cluster$countries))

# NUMBER OF CLUSTER
N_clusters = length(unique(countries_cluster$cluster))

# NUMBER OF TIME POINTS
N = length(dates)


###################.

## DATA ARRAY ####

y_array = array(unlist(y), dim = c(N, N_countries))

x_array =  array(unlist(x), dim = c(N, N_variants, N_countries))


sum_y_minus_x <- array(NA, dim = c(N, N_variants, N_countries))

for (m in 1:N_countries) {
  
  for (v in 2:N_variants) {
    sum_y_minus_x[1, (v-1), m] <- y_array[1, m] - sum(x_array[1, 1:(v-1), m])
    for (i in 2:N) {
      sum_y_minus_x[i, (v-1), m] <- y_array[i, m] - sum(x_array[i, 1:(v-1),m])
    }
  }
}
sum_y_minus_x



#######################.

## NIMBLE MODEL ####

# NUMBER OF CHAINS
nchains = 4

# NUMBER OF CLUSTERS
this_cluster <- makeCluster(nchains)


### FUNCTION ####
run_MCMC <- function(seed,
                     x_data, y, sum_y_minus_x,
                     N, N_states, N_variants, N_clusters, N_countries,
                     cluster_index,
                     niter, nburnin, nthin) {
  
  library(nimble)
  
  ## GET FORWARD ALGORITHM BETA-BINOMIAL FUNCTION ##
  source("forward_algorithm_count_nimblefunction.R")
  
  ## MODEL CODE ##
  hmm_code <- nimbleCode({
    
    ## LOOP OVER VARIANTS
    for (v in 1:N_variants) {
      
      ## LOOP OVER CLUSTERS
      for (cl in 1:N_clusters) {
        
        # INITIAL PROBABILITY 
        p0[v, 1:N_states, cl] ~ ddirch(a_p0[1:N_states])
        
        
        # TRANSITION MATRIX 
        prob_1_to_2[v, cl] ~ dbeta(a1, b1)
        prob_2_to_3[v, cl] ~ dbeta(a2, b2)
        prob_3_to_4[v, cl] ~ dbeta(a3, b3)
        prob_4_to_5[v, cl] ~ dbeta(a4, b4)
        
        
        p[1, 1, v, cl] <- 1 - prob_1_to_2[v, cl]
        p[1, 2, v, cl] <- prob_1_to_2[v, cl]
        p[1, 3, v, cl] <- 0
        p[1, 4, v, cl] <- 0
        p[1, 5, v, cl] <- 0
        
        p[2, 1, v, cl] <- 0
        p[2, 2, v, cl] <- 1 - prob_2_to_3[v, cl]
        p[2, 3, v, cl] <- prob_2_to_3[v, cl]
        p[2, 4, v, cl] <- 0
        p[2, 5, v, cl] <- 0
        
        p[3, 1, v, cl] <- 0
        p[3, 2, v, cl] <- 0
        p[3, 3, v, cl] <- 1 - prob_3_to_4[v, cl]
        p[3, 4, v, cl] <- prob_3_to_4[v, cl]
        p[3, 5, v, cl] <- 0
        
        p[4, 1, v, cl] <- 0
        p[4, 2, v, cl] <- 0
        p[4, 3, v, cl] <- 0
        p[4, 4, v, cl] <- 1 - prob_4_to_5[v, cl]
        p[4, 5, v, cl] <- prob_4_to_5[v, cl]
        
        p[5, 1, v, cl] <- 0
        p[5, 2, v, cl] <- 0
        p[5, 3, v, cl] <- 0
        p[5, 4, v, cl] <- 0
        p[5, 5, v, cl] <- 1
        
        for (s in 1:N_states){
          
          phi[v, s, cl] ~ dgamma(2, 0.05)
          
        }
      }
      
      for (m in 1:N_countries) {
        
        gamma[v, 1, m] <- -10
        gamma[v, N_states, m] <- -10
        
        for (s in 2:(N_states - 1)) {
          
          gamma[v, s, m] ~ dnorm(0, sd = 1)
          
        }
        
        for (s in 1:N_states) {
          
          logit(nu[v, s, m]) <- gamma[v, s, m]
          
        }
        
        constraint[v, m] ~ dconstraint(gamma[v, 3, m] > max(gamma[v, 2, m], gamma[v, 4, m]) &
                                         gamma[v, 2, m] > max(gamma[v, 1, m], gamma[v, 5, m]) &
                                         gamma[v, 4, m] > max(gamma[v, 1, m], gamma[v, 5, m]))
        
      }
    }
    
    for (m in 1:N_countries) {
      
      # FIRST VARIANT
      x_data[1:N, 1, m] ~  dhmm_betabinomial(y = y[1:N, m],
                                             p0 = p0[1, 1:N_states, cluster_index[m]], 
                                             p = p[1:N_states, 1:N_states, 1, cluster_index[m]], 
                                             N = N, 
                                             N_states = N_states, 
                                             nu = nu[1, 1:N_states, m],
                                             phi = phi[1, 1:N_states, cluster_index[m]])
      
      
      # LOOP - ALL OTHER VARIANTS
      for (v in 2:N_variants) {
        
        x_data[1:N, v, m] ~  dhmm_betabinomial(y = sum_y_minus_x[1:N, (v-1), m],
                                               p0 = p0[v, 1:N_states, cluster_index[m]], 
                                               p = p[1:N_states, 1:N_states, v, cluster_index[m]], 
                                               N = N, 
                                               N_states = N_states,
                                               nu = nu[v, 1:N_states, m],
                                               phi = phi[v, 1:N_states, cluster_index[m]])
      }
    }
    
  })
  
  
  ## DATA ##
  data <- list(x_data = x_data,
               constraint = matrix(1, nrow = N_variants, ncol = N_countries))
  
  
  ## CONSTANTS ##
  constants <- list(y = y,
                    sum_y_minus_x = sum_y_minus_x,
                    
                    N = N,
                    N_states = N_states,
                    N_variants = N_variants,
                    N_countries = N_countries,
                    N_clusters = N_clusters,
                    
                    cluster_index = cluster_index,
                    
                    a_p0 = c(0.35, 0.1, 0.1, 0.1, 0.35) * 10,
                    
                    a1 = 1,
                    a2 = 1,
                    a3 = 1,
                    a4 = 1,
                    b1 = 1,
                    b2 = 1,
                    b3 = 1,
                    b4 = 1)
  
  
  ## INITIAL VALUES ##
  inits <- list(p0 = array(rep(matrix(c(0.6, 0.1, 0.1, 0.1, 0.1), byrow = T, nrow = N_variants, ncol = N_states), N_clusters), 
                           dim = c(N_variants, N_states, N_clusters)),       
                
                prob_1_to_2 = matrix(0.05, nrow = N_variants, ncol = N_clusters),
                prob_2_to_3 = matrix(0.05, nrow = N_variants, ncol = N_clusters),
                prob_3_to_4 = matrix(0.05, nrow = N_variants, ncol = N_clusters),
                prob_4_to_5 = matrix(0.05, nrow = N_variants, ncol = N_clusters),
                
                phi = array(0.2, dim = c(N_variants, N_states, N_clusters)),
                
                gamma = array(rep(matrix(c(0.1, 0.2, 0.3, 0.2, 0.1), byrow = T, nrow = N_variants, ncol = N_states), N_countries), 
                              dim = c(N_variants, N_states, N_countries)))
  
  
  ## BUILD MODEL ##
  model <- nimbleModel(code = hmm_code,
                       data = data,
                       constants = constants,
                       inits = inits)
  
  
  ## COMPILE MODEL ##
  compile_model <- compileNimble(model)
  

  ## CONFIGURE MCMC ##
  conf_hmm_model <- configureMCMC(compile_model, 
                                  monitors = c("p0",
                                               "p",
                                               "prob_1_to_2",
                                               "prob_2_to_3",
                                               "prob_3_to_4",
                                               "prob_4_to_5",
                                               "nu",
                                               "gamma",
                                               "phi"
                                  ),
                                  print = TRUE, 
                                  useConjugacy = FALSE,
                                  onlySlice = TRUE)
  
  
  ## BUILD MCMC ##
  hmm_model <- buildMCMC(conf_hmm_model)
  
  
  ## COMPILE MCMC ##
  compile_hmm_model <- compileNimble(hmm_model)
  
  
  ## RUN ##
  samples <- runMCMC(compile_hmm_model,
                     niter = niter,
                     nburnin = nburnin,
                     thin = nthin,
                     samplesAsCodaMCMC = T,
                     setSeed = seed*10000)
  
  return(samples)
  
}


## RUN MODEL ####
run_time <- system.time({samples <- parLapply(cl = this_cluster, 
                                              X = 1:nchains,
                                              fun = run_MCMC, 
                                              
                                              x_data = x_array, 
                                              y = y_array,
                                              sum_y_minus_x = sum_y_minus_x,
                                              
                                              N = N, 
                                              N_states = N_states,
                                              N_variants = N_variants, 
                                              N_clusters = N_clusters, 
                                              N_countries = N_countries,
                                              
                                              cluster_index = cluster_index,
                                              
                                              niter =   2000, 
                                              nburnin = 1000, 
                                              nthin =      1)})


# CLOSE CLUSTER
stopCluster(this_cluster)


# RUN TIME
run_time


#######################.

## OUTPUT ##

samples_mcmc <- as.mcmc.list(samples)

# COMBINE ALL CHAINS
output <- as_tibble(do.call('rbind', samples))


N_samples <- nrow(output)


## SAVE OUTPUT ####
save(samples, file = 'output/mcmc_output-parallel.RData')


### GELMAN ####

gelman.diag(samples, 
            transform = TRUE, 
            autoburnin = FALSE, 
            multivariate = FALSE)


#### SUBSETS ####

subset_p0 <- function(chain) {
  cols <- grep("^p0", colnames(chain), value = TRUE)
  return(chain[, cols])
}
gelman_p0 <- gelman.diag(lapply(samples, subset_p0), 
                         transform = TRUE, 
                         autoburnin = FALSE, 
                         multivariate = FALSE)
gelman_p0
gelman_p0[[1]][, "Point est."] %>% mean()
gelman_p0[[1]][, "Point est."] %>% median()


subset_p <- function(chain) {
  cols <- grep("^p\\[\\d+(,\\s*\\d+)*\\]$", colnames(as.matrix(chain)), value = TRUE)
  return(chain[, cols])
}
gelman_p <- gelman.diag(lapply(samples, subset_p), 
                         transform = TRUE, 
                         autoburnin = FALSE, 
                         multivariate = FALSE)
gelman_p
gelman_p[[1]][, "Point est."] %>% mean(., na.rm = T)
gelman_p[[1]][, "Point est."] %>% median(., na.rm = T)


subset_prob_1_to_2 <- function(chain) {
  cols <- grep("^prob_1_to_2", colnames(chain), value = TRUE)
  return(chain[, cols])
}
gelman_prob_1_to_2 <- gelman.diag(lapply(samples, subset_prob_1_to_2), 
                        transform = TRUE, 
                        autoburnin = FALSE, 
                        multivariate = FALSE)
gelman_prob_1_to_2
gelman_prob_1_to_2[[1]][, "Point est."] %>% mean()
gelman_prob_1_to_2[[1]][, "Point est."] %>% median()


subset_prob_2_to_3 <- function(chain) {
  cols <- grep("^prob_2_to_3", colnames(chain), value = TRUE)
  return(chain[, cols])
}
gelman_prob_2_to_3 <- gelman.diag(lapply(samples, subset_prob_2_to_3), 
                                  transform = TRUE, 
                                  autoburnin = FALSE, 
                                  multivariate = FALSE)
gelman_prob_2_to_3
gelman_prob_2_to_3[[1]][, "Point est."] %>% mean()
gelman_prob_2_to_3[[1]][, "Point est."] %>% median()


subset_prob_3_to_4 <- function(chain) {
  cols <- grep("^prob_3_to_4", colnames(chain), value = TRUE)
  return(chain[, cols])
}
gelman_prob_3_to_4 <- gelman.diag(lapply(samples, subset_prob_3_to_4), 
                                  transform = TRUE, 
                                  autoburnin = FALSE, 
                                  multivariate = FALSE)
gelman_prob_3_to_4
gelman_prob_3_to_4[[1]][, "Point est."] %>% mean()
gelman_prob_3_to_4[[1]][, "Point est."] %>% median()


subset_prob_4_to_5 <- function(chain) {
  cols <- grep("^prob_4_to_5", colnames(chain), value = TRUE)
  return(chain[, cols])
}
gelman_prob_4_to_5 <- gelman.diag(lapply(samples, subset_prob_4_to_5), 
                                  transform = TRUE, 
                                  autoburnin = FALSE, 
                                  multivariate = FALSE)
gelman_prob_4_to_5
gelman_prob_4_to_5[[1]][, "Point est."] %>% mean()
gelman_prob_4_to_5[[1]][, "Point est."] %>% median()


subset_nu <- function(chain) {
  cols <- grep("^nu", colnames(chain), value = TRUE)
  return(chain[, cols])
}
gelman_nu <- gelman.diag(lapply(samples, subset_nu), 
                         transform = TRUE, 
                         autoburnin = FALSE, 
                         multivariate = FALSE)
gelman_nu
gelman_nu[[1]][, "Point est."] %>% mean(., na.rm = T)
gelman_nu[[1]][, "Point est."] %>% median(., na.rm = T)


subset_gamma <- function(chain) {
  cols <- grep("^gamma", colnames(chain), value = TRUE)
  return(chain[, cols])
}
gelman_gamma <- gelman.diag(lapply(samples, subset_gamma), 
                         transform = TRUE, 
                         autoburnin = FALSE, 
                         multivariate = FALSE)
gelman_gamma
gelman_gamma[[1]][, "Point est."] %>% mean(., na.rm = T)
gelman_gamma[[1]][, "Point est."] %>% median(., na.rm = T)


subset_phi <- function(chain) {
  cols <- grep("^phi", colnames(chain), value = TRUE)
  return(chain[, cols])
}
gelman_phi <- gelman.diag(lapply(samples, subset_phi), 
                         transform = TRUE, 
                         autoburnin = FALSE, 
                         multivariate = FALSE)
gelman_phi
gelman_phi[[1]][, "Point est."] %>% mean(., na.rm = T)
gelman_phi[[1]][, "Point est."] %>% median(., na.rm = T)



gelman_output <- list()
gelman_output[["p0"]] <- gelman_p0
gelman_output[["p"]] <- gelman_p
gelman_output[["prob_1_to_2"]] <- gelman_prob_1_to_2
gelman_output[["prob_2_to_2"]] <- gelman_prob_2_to_3
gelman_output[["prob_3_to_4"]] <- gelman_prob_3_to_4
gelman_output[["prob_4_to_5"]] <- gelman_prob_4_to_5
gelman_output[["nu"]] <- gelman_nu
gelman_output[["gamma"]] <- gelman_gamma
gelman_output[["phi"]] <- gelman_phi

sink("output/hmm_gelman-parallel.txt")
print(gelman_output)
sink()


median(c(gelman_p0[[1]][, "Point est."] %>% median(., na.rm=T),
         gelman_p[[1]][, "Point est."] %>% median(., na.rm=T),
         gelman_prob_1_to_2[[1]][, "Point est."] %>% median(., na.rm=T),
         gelman_prob_2_to_3[[1]][, "Point est."] %>% median(., na.rm=T),
         gelman_prob_3_to_4[[1]][, "Point est."] %>% median(., na.rm=T),
         gelman_prob_4_to_5[[1]][, "Point est."] %>% median(., na.rm=T),
         gelman_nu[[1]][, "Point est."] %>% median(., na.rm=T),
         gelman_gamma[[1]][, "Point est."] %>% median(., na.rm=T),
         gelman_phi[[1]][, "Point est."] %>% median(., na.rm=T)))


mean(c(mean(gelman_p0[[1]][, "Point est."]<=1.05, na.rm = TRUE),
       mean(gelman_p[[1]][, "Point est."]<=1.05, na.rm = TRUE),
       mean(gelman_prob_1_to_2[[1]][, "Point est."]<=1.05, na.rm = TRUE),
       mean(gelman_prob_2_to_3[[1]][, "Point est."]<=1.05, na.rm = TRUE),
       mean(gelman_prob_3_to_4[[1]][, "Point est."]<=1.05, na.rm = TRUE),
       mean(gelman_prob_4_to_5[[1]][, "Point est."]<=1.05, na.rm = TRUE),
       mean(gelman_nu[[1]][, "Point est."]<=1.05, na.rm = TRUE),
       mean(gelman_gamma[[1]][, "Point est."]<=1.05, na.rm = TRUE),
       mean(gelman_phi[[1]][, "Point est."]<=1.05, na.rm = TRUE)))

#######################.

## PLOT - OUTPUT ####

### NU ####

#### ALPHA ####

plot(samples_mcmc[,'nu[1, 1, 1]'], main = "nu, Alpha - State 1 - Australia")
plot(samples_mcmc[,'nu[1, 2, 1]'], main = "nu, Alpha - State 2 - Australia")
plot(samples_mcmc[,'nu[1, 3, 1]'], main = "nu, Alpha - State 3 - Australia")
plot(samples_mcmc[,'nu[1, 4, 1]'], main = "nu, Alpha - State 4 - Australia")
plot(samples_mcmc[,'nu[1, 5, 1]'], main = "nu, Alpha - State 5 - Australia")


#### BETA ####

plot(samples_mcmc[,'nu[2, 1, 1]'], main = "nu, Beta - State 1 - Australia")
plot(samples_mcmc[,'nu[2, 2, 1]'], main = "nu, Beta - State 2 - Australia")
plot(samples_mcmc[,'nu[2, 3, 1]'], main = "nu, Beta - State 3 - Australia")
plot(samples_mcmc[,'nu[2, 4, 1]'], main = "nu, Beta - State 4 - Australia")
plot(samples_mcmc[,'nu[2, 5, 1]'], main = "nu, Beta - State 5 - Australia")


#### GAMMA ####

plot(samples_mcmc[,'nu[3, 1, 1]'], main = "nu, Gamma - State 1 - Australia")
plot(samples_mcmc[,'nu[3, 2, 1]'], main = "nu, Gamma - State 2 - Australia")
plot(samples_mcmc[,'nu[3, 3, 1]'], main = "nu, Gamma - State 3 - Australia")
plot(samples_mcmc[,'nu[3, 4, 1]'], main = "nu, Gamma - State 4 - Australia")
plot(samples_mcmc[,'nu[3, 5, 1]'], main = "nu, Gamma - State 5 - Australia")


#### DELTA ####

plot(samples_mcmc[,'nu[4, 1, 1]'], main = "nu, Delta - State 1 - Australia")
plot(samples_mcmc[,'nu[4, 2, 1]'], main = "nu, Delta - State 2 - Australia")
plot(samples_mcmc[,'nu[4, 3, 1]'], main = "nu, Delta - State 3 - Australia")
plot(samples_mcmc[,'nu[4, 4, 1]'], main = "nu, Delta - State 4 - Australia")
plot(samples_mcmc[,'nu[4, 5, 1]'], main = "nu, Delta - State 5 - Australia")


#### OMICRON ####

plot(samples_mcmc[,'nu[5, 1, 1]'], main = "nu, Omicron - State 1 - Australia")
plot(samples_mcmc[,'nu[5, 2, 1]'], main = "nu, Omicron - State 2 - Australia")
plot(samples_mcmc[,'nu[5, 3, 1]'], main = "nu, Omicron - State 3 - Australia")
plot(samples_mcmc[,'nu[5, 4, 1]'], main = "nu, Omicron - State 4 - Australia")
plot(samples_mcmc[,'nu[5, 5, 1]'], main = "nu, Omicron - State 5 - Australia")


#### BOXPLOTS ####

##### ALPHA ####

alpha_state_1 <- matrix(NA, nrow = N_samples, ncol = N_countries)
for (i in 1:N_countries) alpha_state_1[, i] <- c(as.vector(samples[[1]][, paste0("nu[1, 1, ", i, "]")]), as.vector(samples[[2]][, paste0("nu[1, 1, ", i, "]")]))
alpha_state_1 <- alpha_state_1 %>% 
  as.data.frame(.) %>%
  setNames(1:N_countries) %>%
  pivot_longer(., cols = everything(),
               names_to = "country",
               values_to = "state")

alpha_state_2 <- matrix(NA, nrow = N_samples, ncol = N_countries)
for (i in 1:N_countries) alpha_state_2[, i] <- c(as.vector(samples[[1]][, paste0("nu[1, 2, ", i, "]")]), as.vector(samples[[2]][, paste0("nu[1, 2, ", i, "]")])) 
alpha_state_2 <- alpha_state_2 %>% 
  as.data.frame(.) %>%
  setNames(1:N_countries) %>%
  pivot_longer(., cols = everything(),
               names_to = "country",
               values_to = "state")

alpha_state_3 <- matrix(NA, nrow = N_samples, ncol = N_countries)
for (i in 1:N_countries) alpha_state_3[, i] <- c(as.vector(samples[[1]][, paste0("nu[1, 3, ", i, "]")]), as.vector(samples[[2]][, paste0("nu[1, 3, ", i, "]")]))  
alpha_state_3 <- alpha_state_3 %>% 
  as.data.frame(.) %>%
  setNames(1:N_countries) %>%
  pivot_longer(., cols = everything(),
               names_to = "country",
               values_to = "state")

alpha_state_4 <- matrix(NA, nrow = N_samples, ncol = N_countries)
for (i in 1:N_countries) alpha_state_4[, i] <- c(as.vector(samples[[1]][, paste0("nu[1, 4, ", i, "]")]), as.vector(samples[[2]][, paste0("nu[1, 4, ", i, "]")]))  
alpha_state_4 <- alpha_state_4 %>% 
  as.data.frame(.) %>%
  setNames(1:N_countries) %>%
  pivot_longer(., cols = everything(),
               names_to = "country",
               values_to = "state")

alpha_state_5 <- matrix(NA, nrow = N_samples, ncol = N_countries)
for (i in 1:N_countries) alpha_state_5[, i] <- c(as.vector(samples[[1]][, paste0("nu[1, 5, ", i, "]")]), as.vector(samples[[2]][, paste0("nu[1, 5, ", i, "]")])) 
alpha_state_5 <- alpha_state_5 %>% 
  as.data.frame(.) %>%
  setNames(1:N_countries) %>%
  pivot_longer(., cols = everything(),
               names_to = "country",
               values_to = "state")


nu_alpha <- array(NA, dim = c(N_samples, N_states, N_countries))
for (i in 1:N_countries) {
  for (j in 1:N_states) {
    nu_alpha[, j, i] <- get(paste0("alpha_state_", j)) %>% filter(country == i) %>% select(state) %>% unlist() 
  }
}


for (i in 1:N_countries) {
  plot(nu_alpha[, , i] %>%
         as.data.frame(.) %>%
         pivot_longer(., cols = everything(),
                      names_to = "state",
                      values_to = "nu") %>%
         ggplot(., aes(x = state, y = nu)) +
         geom_boxplot() +
         ggtitle(paste0(countries_cluster$countries[i], " - alpha")))
}


##### BETA ####

beta_state_1 <- matrix(NA, nrow = N_samples, ncol = N_countries)
for (i in 1:N_countries) beta_state_1[, i] <- c(as.vector(samples[[1]][, paste0("nu[2, 1, ", i, "]")]), as.vector(samples[[2]][, paste0("nu[2, 1, ", i, "]")]))
beta_state_1 <- beta_state_1 %>% 
  as.data.frame(.) %>%
  setNames(1:N_countries) %>%
  pivot_longer(., cols = everything(),
               names_to = "country",
               values_to = "state")

beta_state_2 <- matrix(NA, nrow = N_samples, ncol = N_countries)
for (i in 1:N_countries) beta_state_2[, i] <- c(as.vector(samples[[1]][, paste0("nu[2, 2, ", i, "]")]), as.vector(samples[[2]][, paste0("nu[2, 2, ", i, "]")])) 
beta_state_2 <- beta_state_2 %>% 
  as.data.frame(.) %>%
  setNames(1:N_countries) %>%
  pivot_longer(., cols = everything(),
               names_to = "country",
               values_to = "state")

beta_state_3 <- matrix(NA, nrow = N_samples, ncol = N_countries)
for (i in 1:N_countries) beta_state_3[, i] <- c(as.vector(samples[[1]][, paste0("nu[2, 3, ", i, "]")]), as.vector(samples[[2]][, paste0("nu[2, 3, ", i, "]")]))  
beta_state_3 <- beta_state_3 %>% 
  as.data.frame(.) %>%
  setNames(1:N_countries) %>%
  pivot_longer(., cols = everything(),
               names_to = "country",
               values_to = "state")

beta_state_4 <- matrix(NA, nrow = N_samples, ncol = N_countries)
for (i in 1:N_countries) beta_state_4[, i] <- c(as.vector(samples[[1]][, paste0("nu[2, 4, ", i, "]")]), as.vector(samples[[2]][, paste0("nu[2, 4, ", i, "]")]))  
beta_state_4 <- beta_state_4 %>% 
  as.data.frame(.) %>%
  setNames(1:N_countries) %>%
  pivot_longer(., cols = everything(),
               names_to = "country",
               values_to = "state")

beta_state_5 <- matrix(NA, nrow = N_samples, ncol = N_countries)
for (i in 1:N_countries) beta_state_5[, i] <- c(as.vector(samples[[1]][, paste0("nu[2, 5, ", i, "]")]), as.vector(samples[[2]][, paste0("nu[2, 5, ", i, "]")])) 
beta_state_5 <- beta_state_5 %>% 
  as.data.frame(.) %>%
  setNames(1:N_countries) %>%
  pivot_longer(., cols = everything(),
               names_to = "country",
               values_to = "state")


nu_beta <- array(NA, dim = c(N_samples, N_states, N_countries))
for (i in 1:N_countries) {
  for (j in 1:N_states) {
    nu_beta[, j, i] <- get(paste0("beta_state_", j)) %>% filter(country == i) %>% select(state) %>% unlist() 
  }
}


for (i in 1:N_countries) {
  plot(nu_beta[, , i] %>%
         as.data.frame(.) %>%
         pivot_longer(., cols = everything(),
                      names_to = "state",
                      values_to = "nu") %>%
         ggplot(., aes(x = state, y = nu)) +
         geom_boxplot() +
         ggtitle(paste0(countries_cluster$countries[i], " - beta")))
}



##### GAMMA ####

gamma_state_1 <- matrix(NA, nrow = N_samples, ncol = N_countries)
for (i in 1:N_countries) gamma_state_1[, i] <- c(as.vector(samples[[1]][, paste0("nu[3, 1, ", i, "]")]), as.vector(samples[[2]][, paste0("nu[3, 1, ", i, "]")]))
gamma_state_1 <- gamma_state_1 %>% 
  as.data.frame(.) %>%
  setNames(1:N_countries) %>%
  pivot_longer(., cols = everything(),
               names_to = "country",
               values_to = "state")

gamma_state_2 <- matrix(NA, nrow = N_samples, ncol = N_countries)
for (i in 1:N_countries) gamma_state_2[, i] <- c(as.vector(samples[[1]][, paste0("nu[3, 2, ", i, "]")]), as.vector(samples[[2]][, paste0("nu[3, 2, ", i, "]")])) 
gamma_state_2 <- gamma_state_2 %>% 
  as.data.frame(.) %>%
  setNames(1:N_countries) %>%
  pivot_longer(., cols = everything(),
               names_to = "country",
               values_to = "state")

gamma_state_3 <- matrix(NA, nrow = N_samples, ncol = N_countries)
for (i in 1:N_countries) gamma_state_3[, i] <- c(as.vector(samples[[1]][, paste0("nu[3, 3, ", i, "]")]), as.vector(samples[[2]][, paste0("nu[3, 3, ", i, "]")]))  
gamma_state_3 <- gamma_state_3 %>% 
  as.data.frame(.) %>%
  setNames(1:N_countries) %>%
  pivot_longer(., cols = everything(),
               names_to = "country",
               values_to = "state")

gamma_state_4 <- matrix(NA, nrow = N_samples, ncol = N_countries)
for (i in 1:N_countries) gamma_state_4[, i] <- c(as.vector(samples[[1]][, paste0("nu[3, 4, ", i, "]")]), as.vector(samples[[2]][, paste0("nu[3, 4, ", i, "]")]))  
gamma_state_4 <- gamma_state_4 %>% 
  as.data.frame(.) %>%
  setNames(1:N_countries) %>%
  pivot_longer(., cols = everything(),
               names_to = "country",
               values_to = "state")

gamma_state_5 <- matrix(NA, nrow = N_samples, ncol = N_countries)
for (i in 1:N_countries) gamma_state_5[, i] <- c(as.vector(samples[[1]][, paste0("nu[3, 5, ", i, "]")]), as.vector(samples[[2]][, paste0("nu[3, 5, ", i, "]")])) 
gamma_state_5 <- gamma_state_5 %>% 
  as.data.frame(.) %>%
  setNames(1:N_countries) %>%
  pivot_longer(., cols = everything(),
               names_to = "country",
               values_to = "state")


nu_gamma <- array(NA, dim = c(N_samples, N_states, N_countries))
for (i in 1:N_countries) {
  for (j in 1:N_states) {
    nu_gamma[, j, i] <- get(paste0("gamma_state_", j)) %>% filter(country == i) %>% select(state) %>% unlist() 
  }
}


for (i in 1:N_countries) {
  plot(nu_gamma[, , i] %>%
         as.data.frame(.) %>%
         pivot_longer(., cols = everything(),
                      names_to = "state",
                      values_to = "nu") %>%
         ggplot(., aes(x = state, y = nu)) +
         geom_boxplot() +
         ggtitle(paste0(countries_cluster$countries[i], " - gamma")))
}



##### DELTA ####

delta_state_1 <- matrix(NA, nrow = N_samples, ncol = N_countries)
for (i in 1:N_countries) delta_state_1[, i] <- c(as.vector(samples[[1]][, paste0("nu[4, 1, ", i, "]")]), as.vector(samples[[2]][, paste0("nu[4, 1, ", i, "]")]))
delta_state_1 <- delta_state_1 %>% 
  as.data.frame(.) %>%
  setNames(1:N_countries) %>%
  pivot_longer(., cols = everything(),
               names_to = "country",
               values_to = "state")

delta_state_2 <- matrix(NA, nrow = N_samples, ncol = N_countries)
for (i in 1:N_countries) delta_state_2[, i] <- c(as.vector(samples[[1]][, paste0("nu[4, 2, ", i, "]")]), as.vector(samples[[2]][, paste0("nu[4, 2, ", i, "]")])) 
delta_state_2 <- delta_state_2 %>% 
  as.data.frame(.) %>%
  setNames(1:N_countries) %>%
  pivot_longer(., cols = everything(),
               names_to = "country",
               values_to = "state")

delta_state_3 <- matrix(NA, nrow = N_samples, ncol = N_countries)
for (i in 1:N_countries) delta_state_3[, i] <- c(as.vector(samples[[1]][, paste0("nu[4, 3, ", i, "]")]), as.vector(samples[[2]][, paste0("nu[4, 3, ", i, "]")]))  
delta_state_3 <- delta_state_3 %>% 
  as.data.frame(.) %>%
  setNames(1:N_countries) %>%
  pivot_longer(., cols = everything(),
               names_to = "country",
               values_to = "state")

delta_state_4 <- matrix(NA, nrow = N_samples, ncol = N_countries)
for (i in 1:N_countries) delta_state_4[, i] <- c(as.vector(samples[[1]][, paste0("nu[4, 4, ", i, "]")]), as.vector(samples[[2]][, paste0("nu[4, 4, ", i, "]")]))  
delta_state_4 <- delta_state_4 %>% 
  as.data.frame(.) %>%
  setNames(1:N_countries) %>%
  pivot_longer(., cols = everything(),
               names_to = "country",
               values_to = "state")

delta_state_5 <- matrix(NA, nrow = N_samples, ncol = N_countries)
for (i in 1:N_countries) delta_state_5[, i] <- c(as.vector(samples[[1]][, paste0("nu[4, 5, ", i, "]")]), as.vector(samples[[2]][, paste0("nu[4, 5, ", i, "]")])) 
delta_state_5 <- delta_state_5 %>% 
  as.data.frame(.) %>%
  setNames(1:N_countries) %>%
  pivot_longer(., cols = everything(),
               names_to = "country",
               values_to = "state")


nu_delta <- array(NA, dim = c(N_samples, N_states, N_countries))
for (i in 1:N_countries) {
  for (j in 1:N_states) {
    nu_delta[, j, i] <- get(paste0("delta_state_", j)) %>% filter(country == i) %>% select(state) %>% unlist() 
  }
}


for (i in 1:N_countries) {
  plot(nu_delta[, , i] %>%
         as.data.frame(.) %>%
         pivot_longer(., cols = everything(),
                      names_to = "state",
                      values_to = "nu") %>%
         ggplot(., aes(x = state, y = nu)) +
         geom_boxplot() +
         ggtitle(paste0(countries_cluster$countries[i], " - delta")))
}



##### OMICRON ####

omicron_state_1 <- matrix(NA, nrow = N_samples, ncol = N_countries)
for (i in 1:N_countries) omicron_state_1[, i] <- c(as.vector(samples[[1]][, paste0("nu[5, 1, ", i, "]")]), as.vector(samples[[2]][, paste0("nu[5, 1, ", i, "]")]))
omicron_state_1 <- omicron_state_1 %>% 
  as.data.frame(.) %>%
  setNames(1:N_countries) %>%
  pivot_longer(., cols = everything(),
               names_to = "country",
               values_to = "state")

omicron_state_2 <- matrix(NA, nrow = N_samples, ncol = N_countries)
for (i in 1:N_countries) omicron_state_2[, i] <- c(as.vector(samples[[1]][, paste0("nu[5, 2, ", i, "]")]), as.vector(samples[[2]][, paste0("nu[5, 2, ", i, "]")])) 
omicron_state_2 <- omicron_state_2 %>% 
  as.data.frame(.) %>%
  setNames(1:N_countries) %>%
  pivot_longer(., cols = everything(),
               names_to = "country",
               values_to = "state")

omicron_state_3 <- matrix(NA, nrow = N_samples, ncol = N_countries)
for (i in 1:N_countries) omicron_state_3[, i] <- c(as.vector(samples[[1]][, paste0("nu[5, 3, ", i, "]")]), as.vector(samples[[2]][, paste0("nu[5, 3, ", i, "]")]))  
omicron_state_3 <- omicron_state_3 %>% 
  as.data.frame(.) %>%
  setNames(1:N_countries) %>%
  pivot_longer(., cols = everything(),
               names_to = "country",
               values_to = "state")

omicron_state_4 <- matrix(NA, nrow = N_samples, ncol = N_countries)
for (i in 1:N_countries) omicron_state_4[, i] <- c(as.vector(samples[[1]][, paste0("nu[5, 4, ", i, "]")]), as.vector(samples[[2]][, paste0("nu[5, 4, ", i, "]")]))  
omicron_state_4 <- omicron_state_4 %>% 
  as.data.frame(.) %>%
  setNames(1:N_countries) %>%
  pivot_longer(., cols = everything(),
               names_to = "country",
               values_to = "state")

omicron_state_5 <- matrix(NA, nrow = N_samples, ncol = N_countries)
for (i in 1:N_countries) omicron_state_5[, i] <- c(as.vector(samples[[1]][, paste0("nu[5, 5, ", i, "]")]), as.vector(samples[[2]][, paste0("nu[5, 5, ", i, "]")])) 
omicron_state_5 <- omicron_state_5 %>% 
  as.data.frame(.) %>%
  setNames(1:N_countries) %>%
  pivot_longer(., cols = everything(),
               names_to = "country",
               values_to = "state")


nu_omicron <- array(NA, dim = c(N_samples, N_states, N_countries))
for (i in 1:N_countries) {
  for (j in 1:N_states) {
    nu_omicron[, j, i] <- get(paste0("omicron_state_", j)) %>% filter(country == i) %>% select(state) %>% unlist() 
  }
}


for (i in 1:N_countries) {
  plot(nu_omicron[, , i] %>%
         as.data.frame(.) %>%
         pivot_longer(., cols = everything(),
                      names_to = "state",
                      values_to = "nu") %>%
         ggplot(., aes(x = state, y = nu)) +
         geom_boxplot() +
         ggtitle(paste0(countries_cluster$countries[i], " - omicron")))
}



##### COMBINE VARIANTS ####

nu_alpha_data <- abind(nu_alpha, along = 3) %>%
  melt(varnames = c("sample", "state", "country")) %>%
  mutate(variant = "alpha") %>%
  as.data.frame(.)

nu_beta_data <- abind(nu_beta, along = 3) %>%
  melt(varnames = c("sample", "state", "country")) %>%
  mutate(variant = "beta") %>%
  as.data.frame(.)

nu_gamma_data <- abind(nu_gamma, along = 3) %>%
  melt(varnames = c("sample", "state", "country")) %>%
  mutate(variant = "gamma") %>%
  as.data.frame(.)

nu_delta_data <- abind(nu_delta, along = 3) %>%
  melt(varnames = c("sample", "state", "country")) %>%
  mutate(variant = "delta") %>%
  as.data.frame(.)

nu_omicron_data <- abind(nu_omicron, along = 3) %>%
  melt(varnames = c("sample", "state", "country")) %>%
  mutate(variant = "omicron") %>%
  as.data.frame(.)


nu_plot <- nu_alpha_data %>%
  rbind(nu_beta_data,
        nu_gamma_data,
        nu_delta_data,
        nu_omicron_data) %>%
  mutate(variant = factor(variant, levels = variants)) %>%
  mutate(state = as.factor(state))


for (i in 1:N_countries) {
  plot(nu_plot %>%
         filter(country == i) %>%
         ggplot(., aes(x = state, y = value, fill = variant)) +
         geom_boxplot() + 
         facet_wrap(~variant, scales = "free_y") +
         scale_fill_brewer(palette = "Set2") +
         theme(axis.text.x=element_text(angle=60, hjust=1)) +
         ggtitle(countries_cluster$countries[i])
  )
}


## SAVE PLOT ##

pdf("output/count_model_nu.pdf", width=8,height=6)
for (i in 1:N_countries) {
  plot(nu_plot %>%
         filter(country == i) %>%
         ggplot(., aes(x = state, y = value, fill = variant)) +
         geom_boxplot() + 
         facet_wrap(~variant, scales = "free_y") +
         scale_fill_brewer(palette = "Set2") +
         theme(axis.text.x=element_text(angle=60, hjust=1)) +
         ggtitle(countries_cluster$countries[i])
  )
} 
dev.off()


pdf("output/GDM_boxplot_nu-cluster1.pdf")
nu_plot %>%
  filter(country == 29) %>%
  ggplot(., aes(x = state, y = value, fill = variant)) +
  geom_boxplot() + 
  facet_wrap(~variant) +
  scale_fill_brewer(name = "Variant", palette = "Set2") +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  labs(x = "State", y = expression(nu)) +
  theme(legend.position = "bottom")
dev.off()

pdf("output/GDM_boxplot_nu-cluster2.pdf")
nu_plot %>%
  filter(country == 22) %>%
  ggplot(., aes(x = state, y = value, fill = variant)) +
  geom_boxplot() + 
  facet_wrap(~variant) +
  scale_fill_brewer(name = "Variant", palette = "Set2") +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  labs(x = "State", y = expression(nu)) +
  theme(legend.position = "bottom")
dev.off()

pdf("output/GDM_boxplot_nu-cluster3.pdf")
nu_plot %>%
  filter(country == 30) %>%
  ggplot(., aes(x = state, y = value, fill = variant)) +
  geom_boxplot() + 
  facet_wrap(~variant) +
  scale_fill_brewer(name = "Variant", palette = "Set2") +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  labs(x = "State", y = expression(nu)) +
  theme(legend.position = "bottom")
dev.off()


pdf("output/GDM_boxplot_nu.pdf")
nu_plot %>%
  filter(country == 29) %>%
  mutate(cluster = "cluster 1 - UAE") %>%
  rbind(nu_plot %>% filter(country == 22) %>% mutate(cluster = "cluster 2 - NZ")) %>%
  rbind(nu_plot %>% filter(country == 30) %>% mutate(cluster = "cluster 3 - UK")) %>%
  ggplot(., aes(x = state, y = value, fill = variant)) +
  geom_boxplot() + 
  facet_grid(cluster~variant) +
  scale_fill_brewer(name = "Variant", palette = "Set2") +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  labs(x = "State", y = expression(nu)) +
  theme(legend.position = "bottom")
dev.off()




### PHI ####

#### ALPHA ####

plot(samples_mcmc[,'phi[1, 1, 1]'], main = "Alpha - State 1 - cluster 1")
plot(samples_mcmc[,'phi[1, 1, 2]'], main = "Alpha - State 1 - cluster 2")
plot(samples_mcmc[,'phi[1, 1, 3]'], main = "Alpha - State 1 - cluster 3")

plot(samples_mcmc[,'phi[1, 2, 1]'], main = "Alpha - State 2 - cluster 1")
plot(samples_mcmc[,'phi[1, 2, 2]'], main = "Alpha - State 2 - cluster 2")
plot(samples_mcmc[,'phi[1, 2, 3]'], main = "Alpha - State 2 - cluster 3")

plot(samples_mcmc[,'phi[1, 3, 1]'], main = "Alpha - State 3 - cluster 1")
plot(samples_mcmc[,'phi[1, 3, 2]'], main = "Alpha - State 3 - cluster 2")
plot(samples_mcmc[,'phi[1, 3, 3]'], main = "Alpha - State 3 - cluster 3")

plot(samples_mcmc[,'phi[1, 4, 1]'], main = "Alpha - State 4 - cluster 1")
plot(samples_mcmc[,'phi[1, 4, 2]'], main = "Alpha - State 4 - cluster 2")
plot(samples_mcmc[,'phi[1, 4, 3]'], main = "Alpha - State 4 - cluster 3")

plot(samples_mcmc[,'phi[1, 5, 1]'], main = "Alpha - State 5 - cluster 1")
plot(samples_mcmc[,'phi[1, 5, 2]'], main = "Alpha - State 5 - cluster 2")
plot(samples_mcmc[,'phi[1, 5, 3]'], main = "Alpha - State 5 - cluster 3")


#### BETA ####

plot(samples_mcmc[,'phi[2, 1, 1]'], main = "Beta - State 1 - cluster 1")
plot(samples_mcmc[,'phi[2, 1, 2]'], main = "Beta - State 1 - cluster 2")
plot(samples_mcmc[,'phi[2, 1, 3]'], main = "Beta - State 1 - cluster 3")

plot(samples_mcmc[,'phi[2, 2, 1]'], main = "Beta - State 2 - cluster 1")
plot(samples_mcmc[,'phi[2, 2, 2]'], main = "Beta - State 2 - cluster 2")
plot(samples_mcmc[,'phi[2, 2, 3]'], main = "Beta - State 2 - cluster 3")

plot(samples_mcmc[,'phi[2, 3, 1]'], main = "Beta - State 3 - cluster 1")
plot(samples_mcmc[,'phi[2, 3, 2]'], main = "Beta - State 3 - cluster 2")
plot(samples_mcmc[,'phi[2, 3, 3]'], main = "Beta - State 3 - cluster 3")

plot(samples_mcmc[,'phi[2, 4, 1]'], main = "Beta - State 4 - cluster 1")
plot(samples_mcmc[,'phi[2, 4, 2]'], main = "Beta - State 4 - cluster 2")
plot(samples_mcmc[,'phi[2, 4, 3]'], main = "Beta - State 4 - cluster 3")

plot(samples_mcmc[,'phi[2, 5, 1]'], main = "Beta - State 5 - cluster 1")
plot(samples_mcmc[,'phi[2, 5, 2]'], main = "Beta - State 5 - cluster 2")
plot(samples_mcmc[,'phi[2, 5, 3]'], main = "Beta - State 5 - cluster 3")


#### GAMMA ####

plot(samples_mcmc[,'phi[3, 1, 1]'], main = "Gamma - State 1 - cluster 1")
plot(samples_mcmc[,'phi[3, 1, 2]'], main = "Gamma - State 1 - cluster 2")
plot(samples_mcmc[,'phi[3, 1, 3]'], main = "Gamma - State 1 - cluster 3")

plot(samples_mcmc[,'phi[3, 2, 1]'], main = "Gamma - State 2 - cluster 1")
plot(samples_mcmc[,'phi[3, 1, 2]'], main = "Gamma - State 1 - cluster 2")
plot(samples_mcmc[,'phi[3, 1, 3]'], main = "Gamma - State 1 - cluster 3")

plot(samples_mcmc[,'phi[3, 3, 1]'], main = "Gamma - State 3 - cluster 1")
plot(samples_mcmc[,'phi[3, 3, 2]'], main = "Gamma - State 3 - cluster 2")
plot(samples_mcmc[,'phi[3, 3, 3]'], main = "Gamma - State 3 - cluster 3")

plot(samples_mcmc[,'phi[3, 4, 1]'], main = "Gamma - State 4 - cluster 1")
plot(samples_mcmc[,'phi[3, 4, 2]'], main = "Gamma - State 4 - cluster 2")
plot(samples_mcmc[,'phi[3, 4, 3]'], main = "Gamma - State 4 - cluster 3")

plot(samples_mcmc[,'phi[3, 5, 1]'], main = "Gamma - State 5 - cluster 1")
plot(samples_mcmc[,'phi[3, 5, 2]'], main = "Gamma - State 5 - cluster 2")
plot(samples_mcmc[,'phi[3, 5, 3]'], main = "Gamma - State 5 - cluster 3")


#### DELTA ####

plot(samples_mcmc[,'phi[4, 1, 1]'], main = "Delta - State 1 - cluster 1")
plot(samples_mcmc[,'phi[4, 1, 2]'], main = "Delta - State 1 - cluster 2")
plot(samples_mcmc[,'phi[4, 1, 3]'], main = "Delta - State 1 - cluster 3")

plot(samples_mcmc[,'phi[4, 2, 1]'], main = "Delta - State 2 - cluster 1")
plot(samples_mcmc[,'phi[4, 2, 2]'], main = "Delta - State 2 - cluster 2")
plot(samples_mcmc[,'phi[4, 2, 3]'], main = "Delta - State 2 - cluster 3")

plot(samples_mcmc[,'phi[4, 3, 1]'], main = "Delta - State 3 - cluster 1")
plot(samples_mcmc[,'phi[4, 3, 2]'], main = "Delta - State 3 - cluster 2")
plot(samples_mcmc[,'phi[4, 3, 3]'], main = "Delta - State 3 - cluster 3")

plot(samples_mcmc[,'phi[4, 4, 1]'], main = "Delta - State 4 - cluster 1")
plot(samples_mcmc[,'phi[4, 4, 2]'], main = "Delta - State 4 - cluster 2")
plot(samples_mcmc[,'phi[4, 4, 3]'], main = "Delta - State 4 - cluster 3")

plot(samples_mcmc[,'phi[4, 5, 1]'], main = "Delta - State 5 - cluster 1")
plot(samples_mcmc[,'phi[4, 5, 2]'], main = "Delta - State 5 - cluster 2")
plot(samples_mcmc[,'phi[4, 5, 3]'], main = "Delta - State 5 - cluster 3")


#### OMICRON ####

plot(samples_mcmc[,'phi[5, 1, 1]'], main = "Omicron - State 1 - cluster 1")
plot(samples_mcmc[,'phi[5, 1, 2]'], main = "Omicron - State 1 - cluster 2")
plot(samples_mcmc[,'phi[5, 1, 3]'], main = "Omicron - State 1 - cluster 3")

plot(samples_mcmc[,'phi[5, 2, 1]'], main = "Omicron - State 2 - cluster 1")
plot(samples_mcmc[,'phi[5, 2, 2]'], main = "Omicron - State 2 - cluster 2")
plot(samples_mcmc[,'phi[5, 2, 3]'], main = "Omicron - State 2 - cluster 3")

plot(samples_mcmc[,'phi[5, 3, 1]'], main = "Omicron - State 3 - cluster 1")
plot(samples_mcmc[,'phi[5, 3, 2]'], main = "Omicron - State 3 - cluster 2")
plot(samples_mcmc[,'phi[5, 3, 3]'], main = "Omicron - State 3 - cluster 3")

plot(samples_mcmc[,'phi[5, 4, 1]'], main = "Omicron - State 4 - cluster 1")
plot(samples_mcmc[,'phi[5, 4, 2]'], main = "Omicron - State 4 - cluster 2")
plot(samples_mcmc[,'phi[5, 4, 3]'], main = "Omicron - State 4 - cluster 3")

plot(samples_mcmc[,'phi[5, 5, 1]'], main = "Omicron - State 5 - cluster 1")
plot(samples_mcmc[,'phi[5, 5, 2]'], main = "Omicron - State 5 - cluster 2")
plot(samples_mcmc[,'phi[5, 5, 3]'], main = "Omicron - State 5 - cluster 3")


### INITIAL STATE ####

#### p0_1 ####

plot(samples_mcmc[,'p0[1, 1, 1]'], main = "p0 - Alpha - State 1 - cluster 1")
plot(samples_mcmc[,'p0[1, 1, 2]'], main = "p0 - Alpha - State 1 - cluster 2")
plot(samples_mcmc[,'p0[1, 1, 3]'], main = "p0 - Alpha - State 1 - cluster 3")

plot(samples_mcmc[,'p0[2, 1, 1]'], main = "p0 - Beta - State 1 - cluster 1")
plot(samples_mcmc[,'p0[2, 1, 2]'], main = "p0 - Beta - State 1 - cluster 2")
plot(samples_mcmc[,'p0[2, 1, 3]'], main = "p0 - Beta - State 1 - cluster 3")

plot(samples_mcmc[,'p0[3, 1, 1]'], main = "p0 - Gamma - State 1 - cluster 1")
plot(samples_mcmc[,'p0[3, 1, 2]'], main = "p0 - Gamma - State 1 - cluster 2")
plot(samples_mcmc[,'p0[3, 1, 3]'], main = "p0 - Gamma - State 1 - cluster 3")

plot(samples_mcmc[,'p0[4, 1, 1]'], main = "p0 - Delta - State 1 - cluster 1")
plot(samples_mcmc[,'p0[4, 1, 2]'], main = "p0 - Delta - State 1 - cluster 2")
plot(samples_mcmc[,'p0[4, 1, 3]'], main = "p0 - Delta - State 1 - cluster 3")

plot(samples_mcmc[,'p0[5, 1, 1]'], main = "p0 - Omicron - State 1 - cluster 1")
plot(samples_mcmc[,'p0[5, 1, 2]'], main = "p0 - Omicron - State 1 - cluster 2")
plot(samples_mcmc[,'p0[5, 1, 3]'], main = "p0 - Omicron - State 1 - cluster 3")


#### p0_2 ####

plot(samples_mcmc[,'p0[1, 2, 1]'], main = "p0 - Alpha - State 2 - cluster 1")
plot(samples_mcmc[,'p0[1, 2, 2]'], main = "p0 - Alpha - State 2 - cluster 2")
plot(samples_mcmc[,'p0[1, 2, 3]'], main = "p0 - Alpha - State 2 - cluster 3")

plot(samples_mcmc[,'p0[2, 2, 1]'], main = "p0 - Beta - State 2 - cluster 1")
plot(samples_mcmc[,'p0[2, 2, 2]'], main = "p0 - Beta - State 2 - cluster 2")
plot(samples_mcmc[,'p0[2, 2, 3]'], main = "p0 - Beta - State 2 - cluster 3")

plot(samples_mcmc[,'p0[3, 2, 1]'], main = "p0 - Gamma - State 2 - cluster 1")
plot(samples_mcmc[,'p0[3, 2, 2]'], main = "p0 - Gamma - State 2 - cluster 2")
plot(samples_mcmc[,'p0[3, 2, 3]'], main = "p0 - Gamma - State 2 - cluster 3")

plot(samples_mcmc[,'p0[4, 2, 1]'], main = "p0 - Delta - State 2 - cluster 1")
plot(samples_mcmc[,'p0[4, 2, 2]'], main = "p0 - Delta - State 2 - cluster 2")
plot(samples_mcmc[,'p0[4, 2, 3]'], main = "p0 - Delta - State 2 - cluster 3")

plot(samples_mcmc[,'p0[5, 2, 1]'], main = "p0 - Omicron - State 2 - cluster 1")
plot(samples_mcmc[,'p0[5, 2, 2]'], main = "p0 - Omicron - State 2 - cluster 2")
plot(samples_mcmc[,'p0[5, 2, 3]'], main = "p0 - Omicron - State 2 - cluster 3")


#### p0_3 ####

plot(samples_mcmc[,'p0[1, 3, 1]'], main = "p0 - Alpha - State 3 - cluster 1")
plot(samples_mcmc[,'p0[1, 3, 2]'], main = "p0 - Alpha - State 3 - cluster 2")
plot(samples_mcmc[,'p0[1, 3, 3]'], main = "p0 - Alpha - State 3 - cluster 3")

plot(samples_mcmc[,'p0[2, 3, 1]'], main = "p0 - Beta - State 3 - cluster 1")
plot(samples_mcmc[,'p0[2, 3, 2]'], main = "p0 - Beta - State 3 - cluster 2")
plot(samples_mcmc[,'p0[2, 3, 3]'], main = "p0 - Beta - State 3 - cluster 3")

plot(samples_mcmc[,'p0[3, 3, 1]'], main = "p0 - Gamma - State 3 - cluster 1")
plot(samples_mcmc[,'p0[3, 3, 2]'], main = "p0 - Gamma - State 3 - cluster 2")
plot(samples_mcmc[,'p0[3, 3, 3]'], main = "p0 - Gamma - State 3 - cluster 3")

plot(samples_mcmc[,'p0[4, 3, 1]'], main = "p0 - Delta - State 3 - cluster 1")
plot(samples_mcmc[,'p0[4, 3, 2]'], main = "p0 - Delta - State 3 - cluster 2")
plot(samples_mcmc[,'p0[4, 3, 3]'], main = "p0 - Delta - State 3 - cluster 3")

plot(samples_mcmc[,'p0[5, 3, 1]'], main = "p0 - Omicron - State 3 - cluster 1")
plot(samples_mcmc[,'p0[5, 3, 2]'], main = "p0 - Omicron - State 3 - cluster 2")
plot(samples_mcmc[,'p0[5, 3, 3]'], main = "p0 - Omicron - State 3 - cluster 3")


#### p0_4 ####

plot(samples_mcmc[,'p0[1, 4, 1]'], main = "p0 - Alpha - State 4 - cluster 1")
plot(samples_mcmc[,'p0[1, 4, 2]'], main = "p0 - Alpha - State 4 - cluster 2")
plot(samples_mcmc[,'p0[1, 4, 3]'], main = "p0 - Alpha - State 4 - cluster 3")

plot(samples_mcmc[,'p0[2, 4, 1]'], main = "p0 - Beta - State 4 - cluster 1")
plot(samples_mcmc[,'p0[2, 4, 2]'], main = "p0 - Beta - State 4 - cluster 2")
plot(samples_mcmc[,'p0[2, 4, 3]'], main = "p0 - Beta - State 4 - cluster 3")

plot(samples_mcmc[,'p0[3, 4, 1]'], main = "p0 - Gamma - State 4 - cluster 1")
plot(samples_mcmc[,'p0[3, 4, 2]'], main = "p0 - Gamma - State 4 - cluster 2")
plot(samples_mcmc[,'p0[3, 4, 3]'], main = "p0 - Gamma - State 4 - cluster 3")

plot(samples_mcmc[,'p0[4, 4, 1]'], main = "p0 - Delta - State 4 - cluster 1")
plot(samples_mcmc[,'p0[4, 4, 2]'], main = "p0 - Delta - State 4 - cluster 2")
plot(samples_mcmc[,'p0[4, 4, 3]'], main = "p0 - Delta - State 4 - cluster 3")

plot(samples_mcmc[,'p0[5, 4, 1]'], main = "p0 - Omicron - State 4 - cluster 1")
plot(samples_mcmc[,'p0[5, 4, 2]'], main = "p0 - Omicron - State 4 - cluster 2")
plot(samples_mcmc[,'p0[5, 4, 3]'], main = "p0 - Omicron - State 4 - cluster 3")


#### p0_5 ####

plot(samples_mcmc[,'p0[1, 5, 1]'], main = "p0 - Alpha - State 5 - cluster 1")
plot(samples_mcmc[,'p0[1, 5, 2]'], main = "p0 - Alpha - State 5 - cluster 2")
plot(samples_mcmc[,'p0[1, 5, 3]'], main = "p0 - Alpha - State 5 - cluster 3")

plot(samples_mcmc[,'p0[2, 5, 1]'], main = "p0 - Beta - State 5 - cluster 1")
plot(samples_mcmc[,'p0[2, 5, 2]'], main = "p0 - Beta - State 5 - cluster 2")
plot(samples_mcmc[,'p0[2, 5, 3]'], main = "p0 - Beta - State 5 - cluster 3")

plot(samples_mcmc[,'p0[3, 5, 1]'], main = "p0 - Gamma - State 5 - cluster 1")
plot(samples_mcmc[,'p0[3, 5, 2]'], main = "p0 - Gamma - State 5 - cluster 2")
plot(samples_mcmc[,'p0[3, 5, 3]'], main = "p0 - Gamma - State 5 - cluster 3")

plot(samples_mcmc[,'p0[4, 5, 1]'], main = "p0 - Delta - State 5 - cluster 1")
plot(samples_mcmc[,'p0[4, 5, 2]'], main = "p0 - Delta - State 5 - cluster 2")
plot(samples_mcmc[,'p0[4, 5, 3]'], main = "p0 - Delta - State 5 - cluster 3")

plot(samples_mcmc[,'p0[5, 5, 1]'], main = "p0 - Omicron - State 5 - cluster 1")
plot(samples_mcmc[,'p0[5, 5, 2]'], main = "p0 - Omicron - State 5 - cluster 2")
plot(samples_mcmc[,'p0[5, 5, 3]'], main = "p0 - Omicron - State 5 - cluster 3")



### TRANSITION MATRIX ####

#### STATE 1 ####

# P_11

plot(samples_mcmc[,'p[1, 1, 1, 1]'], main = "Alpha - Stay in State 1 - cluster 1")
plot(samples_mcmc[,'p[1, 1, 1, 2]'], main = "Alpha - Stay in State 1 - cluster 2")
plot(samples_mcmc[,'p[1, 1, 1, 3]'], main = "Alpha - Stay in State 1 - cluster 3")

plot(samples_mcmc[,'p[1, 1, 2, 1]'], main = "Beta - Stay in State 1 - cluster 1")
plot(samples_mcmc[,'p[1, 1, 2, 2]'], main = "Beta - Stay in State 1 - cluster 2")
plot(samples_mcmc[,'p[1, 1, 2, 3]'], main = "Beta - Stay in State 1 - cluster 3")

plot(samples_mcmc[,'p[1, 1, 3, 1]'], main = "Gamma - Stay in State 1 - cluster 1")
plot(samples_mcmc[,'p[1, 1, 3, 2]'], main = "Gamma - Stay in State 1 - cluster 2")
plot(samples_mcmc[,'p[1, 1, 3, 3]'], main = "Gamma - Stay in State 1 - cluster 3")

plot(samples_mcmc[,'p[1, 1, 4, 1]'], main = "Delta - Stay in State 1 - cluster 1")
plot(samples_mcmc[,'p[1, 1, 4, 2]'], main = "Delta - Stay in State 1 - cluster 2")
plot(samples_mcmc[,'p[1, 1, 4, 3]'], main = "Delta - Stay in State 1 - cluster 3")

plot(samples_mcmc[,'p[1, 1, 5, 1]'], main = "Omicron - Stay in State 1 - cluster 1")
plot(samples_mcmc[,'p[1, 1, 5, 2]'], main = "Omicron - Stay in State 1 - cluster 2")
plot(samples_mcmc[,'p[1, 1, 5, 3]'], main = "Omicron - Stay in State 1 - cluster 3")


# PROB_1_TO_2

plot(samples_mcmc[,'prob_1_to_2[1, 1]'], main = "Alpha - Move from State 1 to State 2 - cluster 1")
plot(samples_mcmc[,'prob_1_to_2[1, 2]'], main = "Alpha - Move from State 1 to State 2 - cluster 2")
plot(samples_mcmc[,'prob_1_to_2[1, 3]'], main = "Alpha - Move from State 1 to State 2 - cluster 3")

plot(samples_mcmc[,'prob_1_to_2[2, 1]'], main = "Beta - Move from State 1 to State 2 - cluster 1")
plot(samples_mcmc[,'prob_1_to_2[2, 2]'], main = "Beta - Move from State 1 to State 2 - cluster 2")
plot(samples_mcmc[,'prob_1_to_2[2, 3]'], main = "Beta - Move from State 1 to State 2 - cluster 3")

plot(samples_mcmc[,'prob_1_to_2[3, 1]'], main = "Gamma - Move from State 1 to State 2 - cluster 1")
plot(samples_mcmc[,'prob_1_to_2[3, 2]'], main = "Gamma - Move from State 1 to State 2 - cluster 2")
plot(samples_mcmc[,'prob_1_to_2[3, 3]'], main = "Gamma - Move from State 1 to State 2 - cluster 3")

plot(samples_mcmc[,'prob_1_to_2[4, 1]'], main = "Delta - Move from State 1 to State 2 - cluster 1")
plot(samples_mcmc[,'prob_1_to_2[4, 2]'], main = "Delta - Move from State 1 to State 2 - cluster 2")
plot(samples_mcmc[,'prob_1_to_2[4, 3]'], main = "Delta - Move from State 1 to State 2 - cluster 3")

plot(samples_mcmc[,'prob_1_to_2[5, 1]'], main = "Omicron - Move from State 1 to State 2 - cluster 1")
plot(samples_mcmc[,'prob_1_to_2[5, 2]'], main = "Omicron - Move from State 1 to State 2 - cluster 2")
plot(samples_mcmc[,'prob_1_to_2[5, 3]'], main = "Omicron - Move from State 1 to State 2 - cluster 3")


#### STATE 2 ####

# P_22

plot(samples_mcmc[,'p[2, 2, 1, 1]'], main = "Alpha - Stay in State 2 - cluster 1")
plot(samples_mcmc[,'p[2, 2, 1, 2]'], main = "Alpha - Stay in State 2 - cluster 2")
plot(samples_mcmc[,'p[2, 2, 1, 3]'], main = "Alpha - Stay in State 2 - cluster 3")

plot(samples_mcmc[,'p[2, 2, 2, 1]'], main = "Beta - Stay in State 2 - cluster 1")
plot(samples_mcmc[,'p[2, 2, 2, 2]'], main = "Beta - Stay in State 2 - cluster 2")
plot(samples_mcmc[,'p[2, 2, 2, 3]'], main = "Beta - Stay in State 2 - cluster 3")

plot(samples_mcmc[,'p[2, 2, 3, 1]'], main = "Gamma - Stay in State 2 - cluster 1")
plot(samples_mcmc[,'p[2, 2, 3, 2]'], main = "Gamma - Stay in State 2 - cluster 2")
plot(samples_mcmc[,'p[2, 2, 3, 3]'], main = "Gamma - Stay in State 2 - cluster 3")

plot(samples_mcmc[,'p[2, 2, 4, 1]'], main = "Delta - Stay in State 2 - cluster 1")
plot(samples_mcmc[,'p[2, 2, 4, 2]'], main = "Delta - Stay in State 2 - cluster 2")
plot(samples_mcmc[,'p[2, 2, 4, 3]'], main = "Delta - Stay in State 2 - cluster 3")

plot(samples_mcmc[,'p[2, 2, 5, 1]'], main = "Omicron - Stay in State 2 - cluster 1")
plot(samples_mcmc[,'p[2, 2, 5, 2]'], main = "Omicron - Stay in State 2 - cluster 2")
plot(samples_mcmc[,'p[2, 2, 5, 3]'], main = "Omicron - Stay in State 2 - cluster 3")


# PROB_2_TO_3

plot(samples_mcmc[,'prob_2_to_3[1, 1]'], main = "Alpha - Move from State 2 to State 3 - cluster 1")
plot(samples_mcmc[,'prob_2_to_3[1, 2]'], main = "Alpha - Move from State 2 to State 3 - cluster 2")
plot(samples_mcmc[,'prob_2_to_3[1, 3]'], main = "Alpha - Move from State 2 to State 3 - cluster 3")

plot(samples_mcmc[,'prob_2_to_3[2, 1]'], main = "Beta - Move from State 2 to State 3 - cluster 1")
plot(samples_mcmc[,'prob_2_to_3[2, 2]'], main = "Beta - Move from State 2 to State 3 - cluster 2")
plot(samples_mcmc[,'prob_2_to_3[2, 3]'], main = "Beta - Move from State 2 to State 3 - cluster 3")

plot(samples_mcmc[,'prob_2_to_3[3, 1]'], main = "Gamma - Move from State 2 to State 3 - cluster 1")
plot(samples_mcmc[,'prob_2_to_3[3, 2]'], main = "Gamma - Move from State 2 to State 3 - cluster 2")
plot(samples_mcmc[,'prob_2_to_3[3, 3]'], main = "Gamma - Move from State 2 to State 3 - cluster 3")

plot(samples_mcmc[,'prob_2_to_3[4, 1]'], main = "Delta - Move from State 2 to State 3 - cluster 1")
plot(samples_mcmc[,'prob_2_to_3[4, 2]'], main = "Delta - Move from State 2 to State 3 - cluster 2")
plot(samples_mcmc[,'prob_2_to_3[4, 3]'], main = "Delta - Move from State 2 to State 3 - cluster 3")

plot(samples_mcmc[,'prob_2_to_3[5, 1]'], main = "Omicron - Move from State 2 to State 3 - cluster 1")
plot(samples_mcmc[,'prob_2_to_3[5, 2]'], main = "Omicron - Move from State 2 to State 3 - cluster 2")
plot(samples_mcmc[,'prob_2_to_3[5, 3]'], main = "Omicron - Move from State 2 to State 3 - cluster 3")


#### STATE 3 ####

# P_33

plot(samples_mcmc[,'p[3, 3, 1, 1]'], main = "Alpha - Stay in State 3 - cluster 1")
plot(samples_mcmc[,'p[3, 3, 1, 2]'], main = "Alpha - Stay in State 3 - cluster 2")
plot(samples_mcmc[,'p[3, 3, 1, 3]'], main = "Alpha - Stay in State 3 - cluster 3")

plot(samples_mcmc[,'p[3, 3, 2, 1]'], main = "Beta - Stay in State 3 - cluster 1")
plot(samples_mcmc[,'p[3, 3, 2, 2]'], main = "Beta - Stay in State 3 - cluster 2")
plot(samples_mcmc[,'p[3, 3, 2, 3]'], main = "Beta - Stay in State 3 - cluster 3")

plot(samples_mcmc[,'p[3, 3, 3, 1]'], main = "Gamma - Stay in State 3 - cluster 1")
plot(samples_mcmc[,'p[3, 3, 3, 2]'], main = "Gamma - Stay in State 3 - cluster 2")
plot(samples_mcmc[,'p[3, 3, 3, 3]'], main = "Gamma - Stay in State 3 - cluster 3")

plot(samples_mcmc[,'p[3, 3, 4, 1]'], main = "Delta - Stay in State 3 - cluster 1")
plot(samples_mcmc[,'p[3, 3, 4, 2]'], main = "Delta - Stay in State 3 - cluster 2")
plot(samples_mcmc[,'p[3, 3, 4, 3]'], main = "Delta - Stay in State 3 - cluster 3")

plot(samples_mcmc[,'p[3, 3, 5, 1]'], main = "Omicron - Stay in State 3 - cluster 1")
plot(samples_mcmc[,'p[3, 3, 5, 2]'], main = "Omicron - Stay in State 3 - cluster 2")
plot(samples_mcmc[,'p[3, 3, 5, 3]'], main = "Omicron - Stay in State 3 - cluster 3")



# PROB_3_TO_4

plot(samples_mcmc[,'prob_3_to_4[1, 1]'], main = "Alpha - Move from State 3 to State 4 - cluster 1")
plot(samples_mcmc[,'prob_3_to_4[1, 2]'], main = "Alpha - Move from State 3 to State 4 - cluster 2")
plot(samples_mcmc[,'prob_3_to_4[1, 3]'], main = "Alpha - Move from State 3 to State 4 - cluster 3")

plot(samples_mcmc[,'prob_3_to_4[2, 1]'], main = "Beta - Move from State 3 to State 4 - cluster 1")
plot(samples_mcmc[,'prob_3_to_4[2, 2]'], main = "Beta - Move from State 3 to State 4 - cluster 2")
plot(samples_mcmc[,'prob_3_to_4[2, 3]'], main = "Beta - Move from State 3 to State 4 - cluster 3")

plot(samples_mcmc[,'prob_3_to_4[3, 1]'], main = "Gamma - Move from State 3 to State 4 - cluster 1")
plot(samples_mcmc[,'prob_3_to_4[3, 2]'], main = "Gamma - Move from State 3 to State 4 - cluster 2")
plot(samples_mcmc[,'prob_3_to_4[3, 3]'], main = "Gamma - Move from State 3 to State 4 - cluster 3")

plot(samples_mcmc[,'prob_3_to_4[4, 1]'], main = "Delta - Move from State 3 to State 4 - cluster 1")
plot(samples_mcmc[,'prob_3_to_4[4, 2]'], main = "Delta - Move from State 3 to State 4 - cluster 2")
plot(samples_mcmc[,'prob_3_to_4[4, 3]'], main = "Delta - Move from State 3 to State 4 - cluster 3")

plot(samples_mcmc[,'prob_3_to_4[5, 1]'], main = "Omicron - Move from State 3 to State 4 - cluster 1")
plot(samples_mcmc[,'prob_3_to_4[5, 2]'], main = "Omicron - Move from State 3 to State 4 - cluster 2")
plot(samples_mcmc[,'prob_3_to_4[5, 3]'], main = "Omicron - Move from State 3 to State 4 - cluster 3")


#### STATE 4 ####

# P_44

plot(samples_mcmc[,'p[4, 4, 1, 1]'], main = "Alpha - Stay in State 4 - cluster 1")
plot(samples_mcmc[,'p[4, 4, 1, 2]'], main = "Alpha - Stay in State 4 - cluster 2")
plot(samples_mcmc[,'p[4, 4, 1, 3]'], main = "Alpha - Stay in State 4 - cluster 3")

plot(samples_mcmc[,'p[4, 4, 2, 1]'], main = "Beta - Stay in State 4 - cluster 1")
plot(samples_mcmc[,'p[4, 4, 2, 2]'], main = "Beta - Stay in State 4 - cluster 2")
plot(samples_mcmc[,'p[4, 4, 2, 3]'], main = "Beta - Stay in State 4 - cluster 3")

plot(samples_mcmc[,'p[4, 4, 3, 1]'], main = "Gamma - Stay in State 4 - cluster 1")
plot(samples_mcmc[,'p[4, 4, 3, 2]'], main = "Gamma - Stay in State 4 - cluster 2")
plot(samples_mcmc[,'p[4, 4, 3, 3]'], main = "Gamma - Stay in State 4 - cluster 3")

plot(samples_mcmc[,'p[4, 4, 4, 1]'], main = "Delta - Stay in State 4 - cluster 1")
plot(samples_mcmc[,'p[4, 4, 4, 2]'], main = "Delta - Stay in State 4 - cluster 2")
plot(samples_mcmc[,'p[4, 4, 4, 3]'], main = "Delta - Stay in State 4 - cluster 3")

plot(samples_mcmc[,'p[4, 4, 5, 1]'], main = "Omicron - Stay in State 4 - cluster 1")
plot(samples_mcmc[,'p[4, 4, 5, 2]'], main = "Omicron - Stay in State 4 - cluster 2")
plot(samples_mcmc[,'p[4, 4, 5, 3]'], main = "Omicron - Stay in State 4 - cluster 3")


# PROB_4_TO_5

plot(samples_mcmc[,'prob_4_to_5[1, 1]'], main = "Alpha - Move from State 4 to State 5 - cluster 1")
plot(samples_mcmc[,'prob_4_to_5[1, 2]'], main = "Alpha - Move from State 4 to State 5 - cluster 2")
plot(samples_mcmc[,'prob_4_to_5[1, 3]'], main = "Alpha - Move from State 4 to State 5 - cluster 3")

plot(samples_mcmc[,'prob_4_to_5[2, 1]'], main = "Beta - Move from State 4 to State 5 - cluster 1")
plot(samples_mcmc[,'prob_4_to_5[2, 2]'], main = "Beta - Move from State 4 to State 5 - cluster 2")
plot(samples_mcmc[,'prob_4_to_5[2, 3]'], main = "Beta - Move from State 4 to State 5 - cluster 3")

plot(samples_mcmc[,'prob_4_to_5[3, 1]'], main = "Gamma - Move from State 4 to State 5 - cluster 1")
plot(samples_mcmc[,'prob_4_to_5[3, 2]'], main = "Gamma - Move from State 4 to State 5 - cluster 2")
plot(samples_mcmc[,'prob_4_to_5[3, 3]'], main = "Gamma - Move from State 4 to State 5 - cluster 3")

plot(samples_mcmc[,'prob_4_to_5[4, 1]'], main = "Delta - Move from State 4 to State 5 - cluster 1")
plot(samples_mcmc[,'prob_4_to_5[4, 2]'], main = "Delta - Move from State 4 to State 5 - cluster 2")
plot(samples_mcmc[,'prob_4_to_5[4, 3]'], main = "Delta - Move from State 4 to State 5 - cluster 3")

plot(samples_mcmc[,'prob_4_to_5[5, 1]'], main = "Omicron - Move from State 4 to State 5 - cluster 1")
plot(samples_mcmc[,'prob_4_to_5[5, 2]'], main = "Omicron - Move from State 4 to State 5 - cluster 2")
plot(samples_mcmc[,'prob_4_to_5[5, 3]'], main = "Omicron - Move from State 4 to State 5 - cluster 3")


#######################.


## TRACEPLOTS OUTPUT ####

### P0 ####

p0_state_1_cluster_3 <- samples_mcmc[, c(
  "p0[1, 1, 3]", 
  "p0[2, 1, 3]",
  "p0[3, 1, 3]",
  "p0[4, 1, 3]",
  "p0[5, 1, 3]")]
pdf("output/p0_state_1_cluster_3.pdf", width = 10, height = 6)
layout(matrix(c(1:5), 1, 5, byrow = TRUE))
traplot(p0_state_1_cluster_3, main = c("alpha", "beta", "gamma", "delta", "omicron"))
dev.off()


p0_state_2_cluster_3 <- samples_mcmc[, c(
  "p0[1, 2, 3]", 
  "p0[2, 2, 3]",
  "p0[3, 2, 3]",
  "p0[4, 2, 3]",
  "p0[5, 2, 3]")]
pdf("output/p0_state_2_cluster_3.pdf", width = 8, height = 6)
layout(matrix(c(1:5), 1, 5, byrow = TRUE))
traplot(p0_state_2_cluster_3, main = c("alpha", "beta", "gamma", "delta", "omicron"))
dev.off()


p0_state_3_cluster_3 <- samples_mcmc[, c(
  "p0[1, 3, 3]", 
  "p0[2, 3, 3]",
  "p0[3, 3, 3]",
  "p0[4, 3, 3]",
  "p0[5, 3, 3]")]
pdf("output/p0_state_3_cluster_3.pdf", width = 8, height = 6)
layout(matrix(c(1:5), 1, 5, byrow = TRUE))
traplot(p0_state_3_cluster_3, main = c("alpha", "beta", "gamma", "delta", "omicron"))
dev.off()


p0_state_4_cluster_3 <- samples_mcmc[, c(
  "p0[1, 4, 3]", 
  "p0[2, 4, 3]",
  "p0[3, 4, 3]",
  "p0[4, 4, 3]",
  "p0[5, 4, 3]")]
pdf("output/p0_state_4_cluster_3.pdf", width = 8, height = 6)
layout(matrix(c(1:5), 1, 5, byrow = TRUE))
traplot(p0_state_4_cluster_3, main = c("alpha", "beta", "gamma", "delta", "omicron"))
dev.off()


p0_state_5_cluster_3 <- samples_mcmc[, c(
  "p0[1, 5, 3]", 
  "p0[2, 5, 3]",
  "p0[3, 5, 3]",
  "p0[4, 5, 3]",
  "p0[5, 5, 3]")]
pdf("output/p0_state_5_cluster_3.pdf", width = 8, height = 6)
layout(matrix(c(1:5), 1, 5, byrow = TRUE))
traplot(p0_state_5_cluster_3, main = c("alpha", "beta", "gamma", "delta", "omicron"))
dev.off()


### PHI ####


phi_state_1_cluster_3 <- samples_mcmc[, c(
  "phi[1, 1, 3]", 
  "phi[2, 1, 3]",
  "phi[3, 1, 3]",
  "phi[4, 1, 3]",
  "phi[5, 1, 3]")]
pdf("output/phi_state_1_cluster_3.pdf", width = 10, height = 6)
layout(matrix(c(1:5), 1, 5, byrow = TRUE))
traplot(phi_state_1_cluster_3, main = c("alpha", "beta", "gamma", "delta", "omicron"))
dev.off()


phi_state_2_cluster_3 <- samples_mcmc[, c(
  "phi[1, 2, 3]", 
  "phi[2, 2, 3]",
  "phi[3, 2, 3]",
  "phi[4, 2, 3]",
  "phi[5, 2, 3]")]
pdf("output/phi_state_2_cluster_3.pdf", width = 8, height = 6)
layout(matrix(c(1:5), 1, 5, byrow = TRUE))
traplot(phi_state_2_cluster_3, main = c("alpha", "beta", "gamma", "delta", "omicron"))
dev.off()


phi_state_3_cluster_3 <- samples_mcmc[, c(
  "phi[1, 3, 3]", 
  "phi[2, 3, 3]",
  "phi[3, 3, 3]",
  "phi[4, 3, 3]",
  "phi[5, 3, 3]")]
pdf("output/phi_state_3_cluster_3.pdf", width = 8, height = 6)
layout(matrix(c(1:5), 1, 5, byrow = TRUE))
traplot(phi_state_3_cluster_3, main = c("alpha", "beta", "gamma", "delta", "omicron"))
dev.off()


phi_state_4_cluster_3 <- samples_mcmc[, c(
  "phi[1, 4, 3]", 
  "phi[2, 4, 3]",
  "phi[3, 4, 3]",
  "phi[4, 4, 3]",
  "phi[5, 4, 3]")]
pdf("output/phi_state_4_cluster_3.pdf", width = 8, height = 6)
layout(matrix(c(1:5), 1, 5, byrow = TRUE))
traplot(phi_state_4_cluster_3, main = c("alpha", "beta", "gamma", "delta", "omicron"))
dev.off()


phi_state_5_cluster_3 <- samples_mcmc[, c(
  "phi[1, 5, 3]", 
  "phi[2, 5, 3]",
  "phi[3, 5, 3]",
  "phi[4, 5, 3]",
  "phi[5, 5, 3]")]
pdf("output/phi_state_5_cluster_3.pdf", width = 8, height = 6)
layout(matrix(c(1:5), 1, 5, byrow = TRUE))
traplot(phi_state_5_cluster_3, main = c("alpha", "beta", "gamma", "delta", "omicron"))
dev.off()


### PROB 1 TO 2 ####

kappa_1_2_cluster_3 <- samples_mcmc[, c(
  "prob_1_to_2[1, 3]",
  "prob_1_to_2[2, 3]",
  "prob_1_to_2[3, 3]",
  "prob_1_to_2[4, 3]",
  "prob_1_to_2[5, 3]")]
pdf("output/kappa_1_2_cluster_3.pdf", width = 8, height = 6)
layout(matrix(c(1:5), 1, 5, byrow = TRUE))
traplot(kappa_1_2_cluster_3, main = c("alpha", "beta", "gamma", "delta", "omicron"))
dev.off()


### PROB 2 TO 3 ####

kappa_2_3_cluster_3 <- samples_mcmc[, c(
  "prob_2_to_3[1, 3]",
  "prob_2_to_3[2, 3]",
  "prob_2_to_3[3, 3]",
  "prob_2_to_3[4, 3]",
  "prob_2_to_3[5, 3]")]
pdf("output/kappa_2_3_cluster_3.pdf", width = 8, height = 6)
layout(matrix(c(1:5), 1, 5, byrow = TRUE))
traplot(kappa_2_3_cluster_3, main = c("alpha", "beta", "gamma", "delta", "omicron"))
dev.off()


### PROB 3 TO 4 ####

kappa_3_4_cluster_3 <- samples_mcmc[, c(
  "prob_3_to_4[1, 3]",
  "prob_3_to_4[2, 3]",
  "prob_3_to_4[3, 3]",
  "prob_3_to_4[4, 3]",
  "prob_3_to_4[5, 3]")]
pdf("output/kappa_3_4_cluster_3.pdf", width = 8, height = 6)
layout(matrix(c(1:5), 1, 5, byrow = TRUE))
traplot(kappa_3_4_cluster_3, main = c("alpha", "beta", "gamma", "delta", "omicron"))
dev.off()


### PROB 4 TO 5 ####

kappa_4_5_cluster_3 <- samples_mcmc[, c(
  "prob_4_to_5[1, 3]",
  "prob_4_to_5[2, 3]",
  "prob_4_to_5[3, 3]",
  "prob_4_to_5[4, 3]",
  "prob_4_to_5[5, 3]")]
pdf("output/kappa_4_5_cluster_3.pdf", width = 8, height = 6)
layout(matrix(c(1:5), 1, 5, byrow = TRUE))
traplot(kappa_4_5_cluster_3, main = c("alpha", "beta", "gamma", "delta", "omicron"))
dev.off()


#######################.

## TRANSITION PROBABILITY ####

p_df <- output%>%select(starts_with("p["))%>%
  unlist()%>%
  array(dim=c(nrow(output),5,5,5,3))%>%
  melt(varnames=c("sample","i","j","variant","cluster"))

kappa_df <- p_df%>%filter(i==j,i!=5,i!=1)%>%
  mutate(kappa=1/(1-value),state=i)%>%
  group_by(state,variant,cluster)%>%
  summarise(lower=quantile(kappa,0.025),
            median=median(kappa),
            upper=quantile(kappa,0.975))

pdf("output/transition_prob_plot.pdf", width = 8, height = 5)
kappa_df %>%
  mutate(variant = case_when(variant == 1 ~ "alpha",
                             variant == 2 ~ "beta",
                             variant == 3 ~ "gamma",
                             variant == 4 ~ "delta",
                             variant == 5 ~ "omicron"),
         cluster=as.factor(cluster),
         state = factor(case_match(state,
                                   1 ~ "Dormant State (State 1)",
                                   2 ~ "Increasing State (State 2)",
                                   3 ~ "Dominant State (State 3)",
                                   4 ~ "Decreasing State (State 4)"),
                        levels=c("Dormant State (State 1)",
                                 "Increasing State (State 2)",
                                 "Dominant State (State 3)",
                                 "Decreasing State (State 4)"))) %>%
  mutate(variant = factor(variant, levels = c("alpha",
                                              "beta",
                                              "gamma",
                                              "delta",
                                              "omicron",
                                              "VOI"))) %>%
  
  ggplot(.,aes(x=variant, col = cluster))+
  facet_wrap(~state)+
  #geom_errorbar(aes(ymin=lower,ymax=upper),alpha=0.5)+
  geom_point(aes(y=median,shape=cluster),size=2.5) +
  scale_color_brewer(palette = "Set2") +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))
  labs(col = "Cluster", 
       shape = "Cluster",
       x = NULL,
       y = "Expected persistence length (weeks)") +
  #ggtitle("Expected persistence lengths for the Active HMM States (States 2,3,4) by Variant and Cluster") +
  theme(legend.position = "bottom")
dev.off()


#######################.

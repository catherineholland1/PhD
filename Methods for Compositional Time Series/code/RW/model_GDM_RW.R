###################################.
## RANDOM WALK TIME SERIES MODEL ##
###################################.

## PACKAGES ####

library(dplyr)
library(tidyr)
library(nimble)
library(mcmcplots)
library(ggplot2)
library(parallel)

#######################.


set.seed(123)


## SET WORKING DIRECTORY ##

setwd("RW")

## GET BETA-BINOMIAL FUNCTION ##

source("beta-binomial_nimblefunction.R")


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
  split(.$country) %>%
  lapply(., function(df) {
    df %>%
      select(variant, starts_with("20")) %>% 
      t() %>% 
      as.data.frame() %>% 
      janitor::row_to_names(row_number = 1) %>%
      mutate(across(everything(), as.numeric)) %>%
      apply(., 1, max)
  })
y



# X = VARIANT COUNT
x <- data %>%
  filter(measure == "count") %>%
  split(.$country) %>%
  lapply(., function(df) {
    df %>%
      select(variant, starts_with("20")) %>% 
      t() %>% 
      as.data.frame() %>% 
      janitor::row_to_names(row_number = 1) %>%
      mutate(across(everything(), as.numeric))
  })
x


add_missing_variants <- function(df) {
  missing_variants <- setdiff(variants, colnames(df))
  if(length(missing_variants) > 0) {
    df[, missing_variants] <- 0
  }
  df <- select(df, alpha, beta, gamma, delta, omicron)
  return(df)
}
x <- lapply(x, add_missing_variants)
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

# NUMBER OF VARIANTS
N_variants = length(unique(variants))

# NUMBER OF COUNTRIES
N_countries = length(unique(countries_cluster$countries))

# NUMBER OF CLUSTER
N_clusters = length(unique(countries_cluster$cluster))

# NUMBER OF TIME POINTS
N = length(dates)


#######################.

## NIMBLE MODEL ####

# NUMBER OF CHAINS
nchains = 4

# NUMBER OF CLUSTERS
this_cluster <- makeCluster(nchains)


### FUNCTION ###
run_MCMC <- function(seed,
                     x_data, y, 
                     N, N_variants, N_clusters, N_countries,
                     cluster_index,
                     niter, nburnin, nthin) {
  
  library(nimble)
  
  ## GET BETA-BINOMIAL FUNCTION ##
  source("beta-binomial_nimblefunction.R")
  
  
  ## MODEL CODE ##
  rw_code <- nimbleCode({
    
    
    for(v in 1:N_variants){
      
      for (cl in 1:N_clusters) {
        
        phi[v, cl] ~ dgamma(2, 0.05)
        
        lambda_mu[v, cl] ~ dnorm(0, sd = 10)
        lambda_sigma[v, cl] ~ T(dnorm(0,sd=1), 0, )
        
      }
      
      for(m in 1:N_countries){
        
        sigma[v, m] ~ T(dnorm(0,sd=1), 0, )
      }
    }
    
    
    ## LOOP OVER COUNTRIES ##
    
    for (m in 1:N_countries) {
      
      ## LOOP OVER VARIANTS ##
      for (v in 1:N_variants) {
        
        ## RANDOM WALK ##
        
        lambda[1, v, m] ~ dnorm(lambda_mu[v, cluster_index[m]], sd = lambda_sigma[v, cluster_index[m]])
        logit(eta[1, v, m]) <- lambda[1, v, m]
        
        
        for (t in 2:N) {
          
          lambda[t, v, m] ~ dnorm(lambda[(t-1), v, m], sd = sigma[v, m]) 
          logit(eta[t, v, m]) <- lambda[t, v, m]
          
        }
        
      }
      
      ## LOOP OVER TIME POINTS ##
      
      for (t in 1:N) {
        
        ## FIRST VARIANT ##
        
        x_data[t, 1, m] ~ dbetabinomial(eta[t, 1, m], 
                                        phi[1, cluster_index[m]], 
                                        y[t, m])
        
        ## LOOP OVER ALL OTHER VARIANTS ##
        
        for (v in 2:N_variants) {
          
          x_data[t, v, m] ~ dbetabinomial(eta[t, v, m], 
                                          phi[v, cluster_index[m]], 
                                          y[t, m] - sum(x_data[t, 1:(v-1), m]))
          
        }
      }
    }
    
  })
  
  
  ## DATA ##
  data <- list(x_data = x_data)
  
  
  ## CONSTANTS ##
  constants <- list(y = y,
                    
                    N = N,
                    N_variants = N_variants,
                    N_countries = N_countries,
                    N_clusters = N_clusters,
                    
                    cluster_index = cluster_index)
  
  
  ## INITIAL VALUES ##
  inits <- list(lambda = array(rnorm(N*N_variants*N_countries,0,0.1), dim = c(N,N_variants,N_countries)),
                lambda_mu = matrix(0, nrow = N_variants, ncol = N_clusters),
                lambda_sigma = matrix(1, nrow = N_variants, ncol = N_clusters),
                
                phi = matrix(runif(N_variants,5,50), nrow = N_variants, ncol = N_clusters),
                sigma = matrix(runif(N_variants*N_countries,0,0.1),nrow=N_variants))
  
  
  ## BUILD MODEL ##
  rw_model <- nimbleModel(code = rw_code,
                          data = data,
                          constants = constants,
                          inits = inits)
  
  
  ## COMPILE MODEL ##
  compile_rw_model <- compileNimble(rw_model)
  
  
  ## CONFIGURE MCMC ##
  conf_rw_model <- configureMCMC(compile_rw_model, 
                                 monitors = c("eta",
                                              "phi",
                                              "lambda",
                                              "lambda_mu",
                                              "lambda_sigma",
                                              "sigma"),
                                 print = TRUE, 
                                 useConjugacy = TRUE)
  
  
  ## BUILD MCMC ##
  mcmc_rw_model <- buildMCMC(conf_rw_model)
  
  
  ## COMPILE MCMC ##
  compile_mcmc_rw_model <- compileNimble(mcmc_rw_model)
  
  
  ## RUN ##
  samples <- runMCMC(compile_mcmc_rw_model,
                     niter = niter,
                     nburnin = nburnin,
                     thin = nthin,
                     samplesAsCodaMCMC = T,
                     setSeed = seed*10000)
  
  return(samples)
  
}



## DATA ARRAY ##

y_array = array(unlist(y), dim = c(N, N_countries))

x_array <- array(unlist(x), dim = c(N, N_variants, N_countries))



## RUN MODEL ####
run_time <- system.time({samples <- parLapply(cl = this_cluster, 
                                              X = 1:nchains,
                                              fun = run_MCMC, 
                                              
                                              x_data = x_array, 
                                              y = y_array,
                                              
                                              N = N, 
                                              N_variants = N_variants, 
                                              N_clusters = N_clusters, 
                                              N_countries = N_countries,
                                              
                                              cluster_index = cluster_index,
                                              
                                              niter =   400000,
                                              nburnin = 300000,
                                              nthin =   100)})


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
save(samples, file = 'output/rw_mcmc_output-parallel.RData')


### GELMAN ####

subset_phi <- function(chain) {
  cols <- grep("^phi", colnames(chain), value = TRUE)
  return(chain[, cols])
}

gelman_phi <- gelman.diag(lapply(samples, subset_phi), 
                          transform = TRUE, 
                          autoburnin = FALSE, 
                          multivariate = FALSE)
gelman_phi
gelman_phi[[1]][, "Point est."] %>% mean()
gelman_phi[[1]][, "Point est."] %>% median()


subset_lambda_mu <- function(chain) {
  cols <- grep("^lambda_mu", colnames(chain), value = TRUE)
  return(chain[, cols])
}

gelman_lambda_mu <- gelman.diag(lapply(samples, subset_lambda_mu), 
                                transform = TRUE, 
                                autoburnin = FALSE, 
                                multivariate = FALSE)
gelman_lambda_mu
gelman_lambda_mu[[1]][, "Point est."] %>% mean()
gelman_lambda_mu[[1]][, "Point est."] %>% median()


subset_lambda_sigma <- function(chain) {
  cols <- grep("^lambda_sigma", colnames(chain), value = TRUE)
  return(chain[, cols])
}

gelman_lambda_sigma <- gelman.diag(lapply(samples, subset_lambda_sigma), 
                                   transform = TRUE, 
                                   autoburnin = FALSE, 
                                   multivariate = FALSE)
gelman_lambda_sigma
gelman_lambda_sigma[[1]][, "Point est."] %>% mean()
gelman_lambda_sigma[[1]][, "Point est."] %>% median()


subset_sigma <- function(chain) {
  cols <- grep("^sigma", colnames(chain), value = TRUE)
  return(chain[, cols])
}

gelman_sigma <- gelman.diag(lapply(samples, subset_sigma), 
                            transform = TRUE, 
                            autoburnin = FALSE, 
                            multivariate = FALSE)
gelman_sigma
gelman_sigma[[1]][, "Point est."] %>% mean()
gelman_sigma[[1]][, "Point est."] %>% median()


gelman_output <- list()
gelman_output[["phi"]] <- gelman_phi
gelman_output[["lambda_mu"]] <- gelman_lambda_mu
gelman_output[["lambda_sigma"]] <- gelman_lambda_sigma
gelman_output[["sigma"]] <- gelman_sigma


sink("output/rw_gelman-parallel.txt")
print(gelman_output)
sink()


median(c(gelman_phi[[1]][, "Point est."] %>% median(., na.rm=T),
         gelman_lambda_mu[[1]][, "Point est."] %>% median(., na.rm=T),
         gelman_lambda_sigma[[1]][, "Point est."] %>% median(., na.rm=T),
         gelman_sigma[[1]][, "Point est."] %>% median(., na.rm=T)))


mean(c(mean(gelman_phi[[1]][, "Point est."]<=1.05, na.rm = TRUE),
       mean(gelman_lambda_mu[[1]][, "Point est."]<=1.05, na.rm = TRUE),
       mean(gelman_lambda_sigma[[1]][, "Point est."]<=1.05, na.rm = TRUE),
       mean(gelman_sigma[[1]][, "Point est."]<=1.05, na.rm = TRUE)))


#######################.

## PLOT - OUTPUT ####

### PHI ####

plot(samples_mcmc[,'phi[1, 1]'], main = "phi - alpha - cluster 1")
plot(samples_mcmc[,'phi[1, 2]'], main = "phi - alpha - cluster 2")
plot(samples_mcmc[,'phi[1, 3]'], main = "phi - alpha - cluster 3")

plot(samples_mcmc[,'phi[2, 1]'], main = "phi - beta - cluster 1")
plot(samples_mcmc[,'phi[2, 2]'], main = "phi - beta - cluster 2")
plot(samples_mcmc[,'phi[2, 3]'], main = "phi - beta - cluster 3")

plot(samples_mcmc[,'phi[3, 1]'], main = "phi - gamma - cluster 1")
plot(samples_mcmc[,'phi[3, 2]'], main = "phi - gamma - cluster 2")
plot(samples_mcmc[,'phi[3, 3]'], main = "phi - gamma - cluster 3")

plot(samples_mcmc[,'phi[4, 1]'], main = "phi - delta - cluster 1")
plot(samples_mcmc[,'phi[4, 2]'], main = "phi - delta - cluster 2")
plot(samples_mcmc[,'phi[4, 3]'], main = "phi - delta - cluster 3")

plot(samples_mcmc[,'phi[5, 1]'], main = "phi - omicron - cluster 1")
plot(samples_mcmc[,'phi[5, 2]'], main = "phi - omicron - cluster 2")
plot(samples_mcmc[,'phi[5, 3]'], main = "phi - omicron - cluster 3")


### LAMBDA MU ####

plot(samples_mcmc[,'lambda_mu[1, 1]'], main = "lambda_mu - alpha - cluster 1")
plot(samples_mcmc[,'lambda_mu[1, 2]'], main = "lambda_mu - alpha - cluster 2")
plot(samples_mcmc[,'lambda_mu[1, 3]'], main = "lambda_mu - alpha - cluster 3")

plot(samples_mcmc[,'lambda_mu[2, 1]'], main = "lambda_mu - beta - cluster 1")
plot(samples_mcmc[,'lambda_mu[2, 2]'], main = "lambda_mu - beta - cluster 2")
plot(samples_mcmc[,'lambda_mu[2, 3]'], main = "lambda_mu - beta - cluster 3")

plot(samples_mcmc[,'lambda_mu[3, 1]'], main = "lambda_mu - gamma - cluster 1")
plot(samples_mcmc[,'lambda_mu[3, 2]'], main = "lambda_mu - gamma - cluster 2")
plot(samples_mcmc[,'lambda_mu[3, 3]'], main = "lambda_mu - gamma - cluster 3")

plot(samples_mcmc[,'lambda_mu[4, 1]'], main = "lambda_mu - delta - cluster 1")
plot(samples_mcmc[,'lambda_mu[4, 2]'], main = "lambda_mu - delta - cluster 2")
plot(samples_mcmc[,'lambda_mu[4, 3]'], main = "lambda_mu - delta - cluster 3")

plot(samples_mcmc[,'lambda_mu[5, 1]'], main = "lambda_mu - omicron - cluster 1")
plot(samples_mcmc[,'lambda_mu[5, 2]'], main = "lambda_mu - omicron - cluster 2")
plot(samples_mcmc[,'lambda_mu[5, 3]'], main = "lambda_mu - omicron - cluster 3")


### LAMBDA SIGMA ####

plot(samples_mcmc[,'lambda_sigma[1, 1]'], main = "lambda_sigma - alpha - cluster 1")
plot(samples_mcmc[,'lambda_sigma[1, 2]'], main = "lambda_sigma - alpha - cluster 2")
plot(samples_mcmc[,'lambda_sigma[1, 3]'], main = "lambda_sigma - alpha - cluster 3")

plot(samples_mcmc[,'lambda_sigma[2, 1]'], main = "lambda_sigma - beta - cluster 1")
plot(samples_mcmc[,'lambda_sigma[2, 2]'], main = "lambda_sigma - beta - cluster 2")
plot(samples_mcmc[,'lambda_sigma[2, 3]'], main = "lambda_sigma - beta - cluster 3")

plot(samples_mcmc[,'lambda_sigma[3, 1]'], main = "lambda_sigma - gamma - cluster 1")
plot(samples_mcmc[,'lambda_sigma[3, 2]'], main = "lambda_sigma - gamma - cluster 2")
plot(samples_mcmc[,'lambda_sigma[3, 3]'], main = "lambda_sigma - gamma - cluster 3")

plot(samples_mcmc[,'lambda_sigma[4, 1]'], main = "lambda_sigma - delta - cluster 1")
plot(samples_mcmc[,'lambda_sigma[4, 2]'], main = "lambda_sigma - delta - cluster 2")
plot(samples_mcmc[,'lambda_sigma[4, 3]'], main = "lambda_sigma - delta - cluster 3")

plot(samples_mcmc[,'lambda_sigma[5, 1]'], main = "lambda_sigma - omicron - cluster 1")
plot(samples_mcmc[,'lambda_sigma[5, 2]'], main = "lambda_sigma - omicron - cluster 2")
plot(samples_mcmc[,'lambda_sigma[5, 3]'], main = "lambda_sigma - omicron - cluster 3")


### SIGMA ####

plot(samples_mcmc[,'sigma[1, 3]'], main = "sigma - alpha, Australia")
plot(samples_mcmc[,'sigma[2, 3]'], main = "sigma - beta, Australia")
plot(samples_mcmc[,'sigma[3, 3]'], main = "sigma - gamma, Australia")
plot(samples_mcmc[,'sigma[4, 3]'], main = "sigma - delta, Australia")
plot(samples_mcmc[,'sigma[5, 3]'], main = "sigma - omicron, Australia")


#######################.

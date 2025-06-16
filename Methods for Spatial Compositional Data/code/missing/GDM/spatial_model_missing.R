############################################################.
## SPATIAL MULTIVARIATE GAM NIMBLE MODEL - MISSING VALUES ##
############################################################.


## PARALLEL CHAINS ##


## PACKAGES ####

library(dplyr)
library(tidyr)
library(nimble)
library(nimbleHMC)
library(coda)
library(mcmcplots)
library(ggplot2)
library(mgcv)
library(parallel)


#######################.

## SET WORKING DIRECTORY ##

setwd("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD R Code/Spatial - Trees/Missing Multivariate Model/GDM")

source("beta-binomial_nimblefunction.R")


#######################.

## DATA ####

tree_data <- readRDS("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD R Code/Spatial - Trees/Data/tree_count_100_data.rds") 


#### SPATIAL LOCATIONS ####

spatial_locations <- tree_data %>%
  select(X_coord, Y_coord)


### TREE TYPES ####

tree_types <- c("larch",
                "oak",
                "sitka_spruce",
                "sycamore")


tree_types_data <- tree_data %>% 
  select(X_coord,
         Y_coord,
         all_of(tree_types),
         total) 


###################.

## RANDOMLY SELECT ROWS ####

N <- nrow(tree_types_data)

N_fit <- 1000

set.seed(451810)
fit_rows <- sample(1:nrow(tree_types_data),N_fit,replace = FALSE)


###################.

## SETUP MISSING DATA ####

selected_data <- tree_types_data[fit_rows,]

missing_tree_data <- selected_data

n_missing <- rep(0,N_fit)


#### SELECT DATA - 1 MISSING ####

missing_1_row <- sample(1:length(fit_rows),N_fit/5,
                        replace = FALSE)
n_missing[missing_1_row] <- 1

#### SELECT DATA - 2 MISSING ####

missing_2_rows <- sample(setdiff(1:length(fit_rows),
                                 missing_1_row),
                         N_fit/5,
                         replace = FALSE)
n_missing[missing_2_rows] <- 2

#### SELECT DATA - 3 MISSING ####

missing_3_rows <- sample(setdiff(1:length(fit_rows),
                                 c(missing_1_row,missing_2_rows)),
                         N_fit/5,
                         replace = FALSE)
n_missing[missing_3_rows] <- 3


#### CREATE NA DATA ####

for(i in 1:N_fit){
  if(n_missing[i]>0){
    missing_tree_data[i,sample(3:6,n_missing[i],replace=FALSE)] <- NA
  }
}

missing_tree_data <- missing_tree_data %>%
  mutate(n_missing = n_missing) %>%
  mutate(across(c(larch:sycamore),
                ~ if_else(is.na(.), TRUE, FALSE),
                .names = "{.col}_missing"))


#### PLOT ####

missing_tree_data %>%
  pivot_longer(larch:sycamore, 
               names_to = "type",
               values_to = "count") %>% 
  mutate(type = factor(type, levels = tree_types)) %>%
  ggplot(., aes(x = X_coord, y = Y_coord, fill = count)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  facet_wrap(~type) +
  labs(x = "x coordinate", y = "y coordinate") +
  theme(legend.position = "bottom") +
  ggtitle("Heatmap for Tree Types - Count")


#### TOTAL COUNT ####

y = missing_tree_data$total


#### TREE COUNTS ####

x_missing = missing_tree_data %>%
  select(larch:sycamore)

x_actual = selected_data %>%
  select(larch:sycamore)


###################.

## CONSTANTS ####

#N = nrow(tree_types_data)

# NUMBER OF ROWS TO FIT
#N_fit <- length(fit_rows)
#N = N_fit

# NUMBER OF TREE TYPES
N_types = length(tree_types)

# NUMBER OF BASIS
N_basis = 400


#######################.

## GAM MODEL ####

larch_gam_bin <- gam(larch/total ~ s(X_coord, Y_coord, k=N_basis), # CHANGED
                     family = "binomial",
                     weights = missing_tree_data$total,
                     data = missing_tree_data)

Z <- predict(larch_gam_bin,
             missing_tree_data,
             type="lpmatrix")


larch_jagam_bin <- jagam(larch/total ~ s(X_coord, Y_coord, k=N_basis), # CHANGED
                         family = "binomial",
                         weights = tree_types_data[fit_rows,]$total,
                         data = tree_types_data[fit_rows,],
                         file = "larch.jags")

S1 <- larch_jagam_bin$jags.data$S1


#######################.

## NIMBLE MODEL ####

# NUMBER OF CORES
n_cores <- min(detectCores(), 100, 20)

# NUMBER OF CHAINS
nchains = 4

# NUMBER OF CLUSTERS
this_cluster <- makeCluster(nchains)

#parallel::clusterSetRNGStream(this_cluster, iseed = 12345)


### SEED ####

initial_seed <- round(runif(nchains, 1234, 8750874))
initial_seed

chain_seed <- round(runif(nchains, 5468, 5632751))
chain_seed


### MODEL FUNCTION ####

run_MCMC <- function(X,
                     initial_seed, chain_seed,
                     x_missing, y,
                     missing_data,
                     Z, S1,
                     N, N_types, N_basis,
                     niter, nburn, nthin) {
  
  library(nimble)
  
  source("beta-binomial_nimblefunction.R")
  
  #### MODEL CODE ####
  
  spatial_model_code <- nimbleCode({
    
    ## EACH TREE TYPE ##
    
    for (k in 1:N_types) {
      
      ## PRIORS ##
      
      for (i in 1:2) {
        
        lambda[i,k] ~ dgamma(0.05, 0.005)
        rho[i,k] <- log(lambda[i,k])
        
      }
      
      #K1[1:(N_basis - 1), 1:(N_basis-1), k] <- S1[1:(N_basis - 1), 1:(N_basis - 1)] * lambda[1,k] + S1[1:(N_basis - 1), N_basis:((N_basis - 1)*2)] * lambda[2,k]
      K1[1:(N_basis - 1), 1:(N_basis-1), k] <- S1[1:(N_basis - 1), 1:(N_basis - 1)] / (sigma[1,k]^2) + S1[1:(N_basis - 1), N_basis:((N_basis - 1)*2)] / (sigma[2,k]^2)
      
      
      b[1,k] ~ dnorm(0, sd = 10)
      b[2:N_basis,k] ~ dmnorm(zero[1:(N_basis-1)], K1[1:(N_basis - 1), 1:(N_basis-1), k])
      
      
      ## LINEAR PREDICTOR ##
      eta[1:N, k] <- Z[1:N, 1:N_basis] %*% b[1:N_basis, k]
      
      
      ## EXPECTED RESPONSE ##
      mu[1:N, k] <- expit(eta[1:N, k])
      
      ## MODEL COUNTS X ##
      mu_mean[k] <- mean(mu[1:N, k])
      
      for (i in 1:N) {
        
        log(phi[i, k]) <- psi[1, k] + 
          psi[2, k]*(mu[i, k]-mu_mean[k]) +
          psi[3, k]*(mu[i, k]-mu_mean[k])^2 +
          psi[4, k]*(mu[i, k]-mu_mean[k])^2
        
      }
      
      psi[1, k] ~ dnorm(2,sd=2)
      
      for(j in 2:4){
        psi[j, k] ~ dnorm(0,sd=1)
      }
      
    }
    
    for (i in 1:N) {
      
      # FIRST TYPE
      x[i, 1] ~ dbetabinomial(mu[i, 1], phi[i, 1], y[i])
      
      for (k in 2:N_types) {
        
        x[i, k] ~ dbetabinomial(mu[i, k], phi[i, k], y[i] - sum(x[i, 1:(k-1)]))
        
      }
    }
    
  })
  
  
  #### DATA ####
  data <- list(x = x_missing)
  
  
  #### CONSTANTS ####
  constants <- list(y = y,
                    Z = Z,
                    S1 = S1,
                    N = N,
                    N_basis = N_basis,
                    N_types = N_types,
                    zero = rep(0, (N_basis - 1)))
  
  
  #### INITIAL VALUES ####
  init_fun <- function(){
    
    x_inits <- matrix(NA, nrow = N, ncol = N_types) 
    for(i in 1:N){
      if(is.na(sum(x_missing[i,]))){
        n_missing <- sum(is.na(x_missing[i,]))
        x_inits[i,is.na(x_missing[i,])] <- rmulti(1,y-sum(x_missing[i,which(!is.na(x_missing[i,]))]),
                                                  prob=rep(1/(n_missing+1),n_missing+1))[1:n_missing]
      }
    }
    
    inits <- list(b = matrix(rnorm(N_basis, 0, 0.25), nrow = N_basis, ncol = N_types),
                  #lambda = matrix(rep(1, 2), nrow = 2, ncol = N_types),
                  sigma = matrix(rep(1,2), nrow = 2, ncol = N_types),
                  log_sigma = matrix(rep(1,2), nrow = 2, ncol = N_types),
                  psi = matrix(c(rnorm(1,2,sd=0.25),rnorm(3,0,0.25)), nrow = 4, ncol = N_types),
                  x = x_inits)
    return(inits)
  }
  
  set.seed(initial_seed[X])
  initial_values <- init_fun()
  
  
  #### BUILD MODEL ####
  model <- nimbleModel(code = spatial_model_code,
                       data = data,
                       constants = constants,
                       inits = init_fun(),
                       buildDerivs = TRUE)
  
  
  
  #### COMPILE MODEL ####
  compile_model <- compileNimble(model)
  
  compile_model$calculate()
  compile_x <- compile_model$x
  
  
  #### CONFIGURE MCMC ####
  conf_model <- configureMCMC(compile_model, 
                              monitors = c("x",
                                           "b",
                                           "log_sigma",
                                           "sigma",
                                           "eta",
                                           "mu",
                                           "phi",
                                           "psi"
                              ),
                              print = TRUE, 
                              useConjugacy = FALSE)
  
  
  #### SAMPLERS ####
  
  conf_model$removeSamplers(c('lambda', "psi", "b"))
  for(j in 1:N_types){
    conf_model$addSampler(target = paste0("b[1:",N_basis,",",j,"]"), type = "AF_slice")
  }
  for(j in 1:N_types){
    conf_model$addSampler(target = paste0("psi[1:",4,",",j,"]"), type = "AF_slice")
  }
  for (node in 1:length(model$expandNodeNames("lambda"))) {
    conf_model$addSampler(target = model$expandNodeNames("lambda")[node], type = "slice")
  }
  
  #### BUILD MCMC ####
  spatial_model <- buildMCMC(conf_model)
  
  
  #### COMPILE MCMC ####
  compile_spatial_model <- compileNimble(spatial_model,
                                         project = model,
                                         resetFunctions = TRUE)
  
  #### RUN MODEL ####
  samples <- runMCMC(compile_spatial_model, 
                     niter =   niter, 
                     nburnin = nburn,
                     thin =   nthin, 
                     samplesAsCodaMCMC = TRUE,
                     setSeed = chain_seed[X])
  
  #return(samples)
  return(list(initial_values = initial_values, samples = samples))
  
}



### RUN MODEL ####
run_time <- system.time({samples <- parLapply(cl = this_cluster, 
                                              X = 1:nchains,
                                              fun = run_MCMC, 
                                              
                                              initial_seed = initial_seed,
                                              chain_seed = chain_seed,
                                              
                                              x_missing = x_missing, 
                                              y = y,
                                              
                                              Z = Z,
                                              S1 = S1,
                                              
                                              N = N_fit, 
                                              N_types = N_types,
                                              N_basis = N_basis,
                                              
                                              niter = 10000, 
                                              nburn =  5000, 
                                              nthin =   10)})




# CLOSE CLUSTER
stopCluster(this_cluster)

run_time


#######################.

## OUTPUT ####

save(samples, file = paste0('output/spatial_model_missing_', Sys.Date(), '_', N_basis, '_20000.RData'))
#load('larch_mcmc_output.RData')


initial_values_list <- lapply(samples, function(chain) chain$initial_values)

samples_list <- lapply(samples, function(chain) chain$samples)


# COMBINE ALL CHAINS
output <- as_tibble(do.call('rbind', samples_list))


N_samples <- nrow(output)


samples_mcmc <- as.mcmc.list(samples_list)


### GELMAN ####

# gelman.diag(samples, 
#             transform = TRUE, 
#             autoburnin = FALSE, 
#             multivariate = FALSE)


subset_b <- function(chain) {
  cols <- grep("^b", colnames(chain), value = TRUE)
  return(chain[, cols])
}

gelman_b <- gelman.diag(lapply(samples_list, subset_b), 
                        transform = TRUE, 
                        autoburnin = FALSE, 
                        multivariate = FALSE)
gelman_b
gelman_b[[1]][, "Point est."] %>% mean()
gelman_b[[1]][, "Point est."] %>% median()

mean(gelman_b[[1]][, "Point est."]<=1.05, na.rm = TRUE)


c(gelman_b[[1]][1:N_basis, "Point est."] %>% mean(),
  gelman_b[[1]][N_basis:(N_basis*2), "Point est."] %>% mean(),
  gelman_b[[1]][(N_basis*2):(N_basis*3), "Point est."] %>% mean(),
  gelman_b[[1]][(N_basis*3):(N_basis*4), "Point est."] %>% mean()
)

c(gelman_b[[1]][1:N_basis, "Point est."] %>% median(),
  gelman_b[[1]][N_basis:(N_basis*2), "Point est."] %>% median(),
  gelman_b[[1]][(N_basis*2):(N_basis*3), "Point est."] %>% median(),
  gelman_b[[1]][(N_basis*3):(N_basis*4), "Point est."] %>% median()
)


subset_log_sigma <- function(chain) {
  cols <- grep("^log_sigma", colnames(chain), value = TRUE)
  return(chain[, cols])
}

gelman_log_sigma <- gelman.diag(lapply(samples_list, subset_log_sigma), 
                            transform = TRUE, 
                            autoburnin = FALSE, 
                            multivariate = FALSE)
gelman_log_sigma
gelman_log_sigma[[1]][, "Point est."] %>% mean()
gelman_log_sigma[[1]][, "Point est."] %>% median()

mean(gelman_log_sigma[[1]][, "Point est."]<=1.05, na.rm = TRUE)


subset_phi <- function(chain) {
  cols <- grep("^phi", colnames(chain), value = TRUE)
  return(chain[, cols])
}

gelman_phi <- gelman.diag(lapply(samples_list, subset_phi), 
                          transform = TRUE, 
                          autoburnin = FALSE, 
                          multivariate = FALSE)
gelman_phi
gelman_phi[[1]][, "Point est."] %>% mean()
gelman_phi[[1]][, "Point est."] %>% median()

mean(gelman_phi[[1]][, "Point est."]<=1.05, na.rm = TRUE)


subset_psi <- function(chain) {
  cols <- grep("^psi", colnames(chain), value = TRUE)
  return(chain[, cols])
}

gelman_psi <- gelman.diag(lapply(samples_list, subset_psi), 
                          transform = TRUE, 
                          autoburnin = FALSE, 
                          multivariate = FALSE)
gelman_psi
gelman_psi[[1]][, "Point est."] %>% mean()
gelman_psi[[1]][, "Point est."] %>% median()

mean(gelman_psi[[1]][, "Point est."]<=1.05, na.rm = TRUE)


median(c(gelman_b[[1]][, "Point est."] %>% median(),
         gelman_log_sigma[[1]][, "Point est."] %>% median(),
         gelman_phi[[1]][, "Point est."] %>% median(),
         gelman_psi[[1]][, "Point est."] %>% median()))

mean(c(mean(gelman_b[[1]][, "Point est."]<=1.05, na.rm = TRUE),
       mean(gelman_log_sigma[[1]][, "Point est."]<=1.05, na.rm = TRUE),
       mean(gelman_phi[[1]][, "Point est."]<=1.05, na.rm = TRUE),
       mean(gelman_psi[[1]][, "Point est."]<=1.05, na.rm = TRUE)))

#######################.

### TRACEPLOT ####

#### SIGMA ####

traplot(samples_mcmc[, "log_sigma[1, 1]"], main = "log_sigma_1, larch")
traplot(samples_mcmc[, "log_sigma[2, 1]"], main = "log_sigma_2, larch")

traplot(samples_mcmc[, "log_sigma[1, 2]"], main = "log_sigma_1, oak")
traplot(samples_mcmc[, "log_sigma[2, 2]"], main = "log_sigma_2, oak")

traplot(samples_mcmc[, "log_sigma[1, 3]"], main = "log_sigma_1, sitka spruce")
traplot(samples_mcmc[, "log_sigma[2, 3]"], main = "log_sigma_2, sitka spruce")

traplot(samples_mcmc[, "log_sigma[1, 4]"], main = "log_sigma_1, sycamore")
traplot(samples_mcmc[, "log_sigma[2, 4]"], main = "log_sigma_2, sycamore")


#### PHI #####

traplot(samples_mcmc[, "phi[1, 1]"], main = "phi_1, larch")
traplot(samples_mcmc[, "phi[2, 1]"], main = "phi_2, larch")

traplot(samples_mcmc[, "phi[1, 2]"], main = "phi_1, oak")
traplot(samples_mcmc[, "phi[2, 2]"], main = "phi_2, oak")

traplot(samples_mcmc[, "phi[1, 3]"], main = "phi_1, sitka spruce")
traplot(samples_mcmc[, "phi[2, 3]"], main = "phi_2, sitka spruce")

traplot(samples_mcmc[, "phi[1, 4]"], main = "phi_1, sycamore")
traplot(samples_mcmc[, "phi[2, 4]"], main = "phi_2, sycamore")


#### PSI ####

traplot(samples_mcmc[, "psi[1, 1]"], main = "psi_1, larch")
traplot(samples_mcmc[, "psi[2, 1]"], main = "psi_2, larch")

traplot(samples_mcmc[, "psi[1, 2]"], main = "psi_1, oak")
traplot(samples_mcmc[, "psi[2, 2]"], main = "psi_2, oak")

traplot(samples_mcmc[, "psi[1, 3]"], main = "psi_1, sitka spruce")
traplot(samples_mcmc[, "psi[2, 3]"], main = "psi_2, sitka spruce")

traplot(samples_mcmc[, "psi[1, 4]"], main = "psi_1, sycamore")
traplot(samples_mcmc[, "psi[2, 4]"], main = "psi_2, sycamore")


#### B ####

traplot(samples_mcmc[, "b[1, 1]"], main = "b_1, larch")
traplot(samples_mcmc[, "b[2, 1]"], main = "b_2, larch")
traplot(samples_mcmc[, "b[3, 1]"], main = "b_3, larch")
traplot(samples_mcmc[, "b[4, 1]"], main = "b_4, larch")
traplot(samples_mcmc[, "b[5, 1]"], main = "b_5, larch")

traplot(samples_mcmc[, "b[10, 1]"], main = "b_10, larch")

traplot(samples_mcmc[, "b[1, 2]"], main = "b_1, oak")
traplot(samples_mcmc[, "b[2, 2]"], main = "b_2, oak")
traplot(samples_mcmc[, "b[3, 2]"], main = "b_3, oak")
traplot(samples_mcmc[, "b[4, 2]"], main = "b_4, oak")
traplot(samples_mcmc[, "b[5, 2]"], main = "b_5, oak")

traplot(samples_mcmc[, "b[10, 2]"], main = "b_10, oak")

traplot(samples_mcmc[, "b[1, 3]"], main = "b_1, sitka spruce")
traplot(samples_mcmc[, "b[2, 3]"], main = "b_2, sitka spruce")
traplot(samples_mcmc[, "b[3, 3]"], main = "b_3, sitka spruce")
traplot(samples_mcmc[, "b[4, 3]"], main = "b_4, sitka spruce")
traplot(samples_mcmc[, "b[5, 3]"], main = "b_5, sitka spruce")

traplot(samples_mcmc[, "b[10, 3]"], main = "b_10, sitka spruce")

traplot(samples_mcmc[, "b[1, 4]"], main = "b_1, sycamore")
traplot(samples_mcmc[, "b[2, 4]"], main = "b_2, sycamore")
traplot(samples_mcmc[, "b[3, 4]"], main = "b_3, sycamore")
traplot(samples_mcmc[, "b[4, 4]"], main = "b_4, sycamore")
traplot(samples_mcmc[, "b[5, 4]"], main = "b_5, sycamore")

traplot(samples_mcmc[, "b[10, 4]"], main = "b_10, sycamore")


#######################.

## SPLINE VARIANCE ####

b_samples <- array(NA, dim = c(N_basis, N_types, N_samples))
for (i in 1:N_samples) {
  b_row <- output[i, ] %>% select(starts_with("b")) 
  b_samples[, , i] <- array(unlist(b_row), dim = c(N_basis, N_types))
}
b_samples


colMeans(apply(b_samples, c(1,2), var))


#######################.

## PREDICTION ####

x_samples <- array(NA, dim = c(N_fit, N_types, N_samples))
for (i in 1:N_samples) {
  x_row <- output[i, ] %>% select(starts_with("x")) 
  x_samples[, , i] <- array(unlist(x_row), dim = c(N_fit, N_types))
}
x_samples


#######################.

## SAVE PREDICTIONS ####

saveRDS(x_samples, "output/predictions_GDM.rds")


#######################.

## EXAMINE PREDICTIONS ####


predicted_comparison_data_median <- x_samples %>%
  apply(., c(1,2), median) %>%
  
  as.data.frame() %>%
  rename(larch = "V1",
         oak = "V2",
         sitka_spruce = "V3",
         sycamore = "V4") %>%
  cbind(missing_tree_data %>% select(c(X_coord, Y_coord,
                                       contains("missing")))) %>%
  pivot_longer(., cols = larch:sycamore, 
               names_to = "type", values_to = "predicted") %>%
  left_join((selected_data %>% pivot_longer(., cols = larch:sycamore, 
                                            names_to = "type", values_to = "actual") %>%
               select(-total)), by = c("X_coord", "Y_coord", "type"))


predicted_comparison_data_mean <- x_samples %>%
  apply(., c(1,2), mean) %>%
  
  as.data.frame() %>%
  rename(larch = "V1",
         oak = "V2",
         sitka_spruce = "V3",
         sycamore = "V4") %>%
  cbind(missing_tree_data %>% select(c(X_coord, Y_coord,
                                       contains("missing")))) %>%
  pivot_longer(., cols = larch:sycamore, 
               names_to = "type", values_to = "predicted") %>%
  left_join((selected_data %>% pivot_longer(., cols = larch:sycamore, 
                                            names_to = "type", values_to = "actual") %>%
               select(-total)), by = c("X_coord", "Y_coord", "type")) %>%
  
  mutate(predicted = round(predicted))


### PLOT ####

#### PRED VS ACTUAL ####

predicted_comparison_data_mean %>%
  
  filter((type == "larch" & larch_missing == TRUE) |
           (type == "oak" & oak_missing == TRUE) |
           (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
           (type == "sycamore" & sycamore_missing == TRUE)) %>%
  
  ggplot(., aes(x = actual, y = predicted)) +
  geom_point() +
  geom_abline(col = "red") +
  facet_wrap(~type*n_missing, ncol = max(n_missing), scales = "free_y") +
  ggtitle("GDM Model - Predicted vs Actual")


#### RESID VS PRED ####

predicted_comparison_data_mean %>%
  
  filter((type == "larch" & larch_missing == TRUE) |
           (type == "oak" & oak_missing == TRUE) |
           (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
           (type == "sycamore" & sycamore_missing == TRUE)) %>%
  
  mutate(residuals = actual - predicted) %>%
  
  ggplot(., aes(x = predicted, y = residuals)) +
  geom_point() +
  facet_wrap(~type*n_missing, ncol = max(n_missing)) +
  ggtitle("GDM Model - Residuals vs Predicted")


#### RESID VS ACTUAL ####

predicted_comparison_data_mean %>%
  
  filter((type == "larch" & larch_missing == TRUE) |
           (type == "oak" & oak_missing == TRUE) |
           (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
           (type == "sycamore" & sycamore_missing == TRUE)) %>%
  
  mutate(residuals = actual - predicted) %>%
  
  ggplot(., aes(x = actual, y = residuals)) +
  geom_point() +
  facet_wrap(~type*n_missing, ncol = max(n_missing)) +
  ggtitle("GDM Model - Residuals vs Actual")



#### HEATMAP ####

predicted_comparison_data_mean %>%
  
  filter((type == "larch" & larch_missing == TRUE) |
           (type == "oak" & oak_missing == TRUE) |
           (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
           (type == "sycamore" & sycamore_missing == TRUE)) %>%
  
  pivot_longer(., cols = actual:predicted, values_to = "count", names_to = "model") %>%
  mutate(type = factor(type, levels = tree_types)) %>%
  
  ggplot(., aes(x = X_coord, y = Y_coord, fill = count)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  facet_wrap(~type*model, ncol = 2) +
  labs(x = "x coordinate", y = "y coordinate") +
  theme(legend.position = "bottom") +
  ggtitle("Heatmap for Tree Types - Actual vs Predicted - GDM Model")


predicted_comparison_data_median %>%
  
  filter((type == "larch" & larch_missing == TRUE) |
           (type == "oak" & oak_missing == TRUE) |
           (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
           (type == "sycamore" & sycamore_missing == TRUE)) %>%
  
  pivot_longer(., cols = actual:predicted, values_to = "count", names_to = "model") %>%
  mutate(type = factor(type, levels = tree_types)) %>%
  
  ggplot(., aes(x = X_coord, y = Y_coord, fill = count)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  facet_wrap(~type*model, ncol = 2) +
  labs(x = "x coordinate", y = "y coordinate") +
  theme(legend.position = "bottom") +
  ggtitle("Heatmap for Tree Types - Actual vs Predicted")



#### DENSITY ####

predicted_comparison_data_mean %>%
  
  filter((type == "larch" & larch_missing == TRUE) |
           (type == "oak" & oak_missing == TRUE) |
           (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
           (type == "sycamore" & sycamore_missing == TRUE)) %>%
  
  ggplot(aes(x = log(predicted), fill = type)) +
  
  geom_density(aes(color = type), alpha = 0.8) +
  
  facet_wrap(~ type + n_missing) +
  
  geom_vline(aes(xintercept = log(mean(actual)), 
                 color = "Original Data"),
             linetype = "dashed", 
             linewidth = 1) +
  
  scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
  scale_fill_manual(name = "Tree Type", values = c("larch" = "#66C2A5",
                                                   "oak" = "#8DA0CB",
                                                   "sitka_spruce" = "#E78AC3",
                                                   "sycamore" = "#FFD92F")) +
  
  labs(title = "Density Plot of Predicted and Actual Counts - GDM Model",
       x = "Log of Predicted Counts", y = "Density") +
  theme(legend.position = "bottom")


predicted_comparison_data_median %>%
  
  filter((type == "larch" & larch_missing == TRUE) |
           (type == "oak" & oak_missing == TRUE) |
           (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
           (type == "sycamore" & sycamore_missing == TRUE)) %>%
  
  ggplot(aes(x = log(predicted), fill = type)) +
  
  geom_density(aes(color = type), alpha = 0.8) +
  
  facet_wrap(~ type + n_missing) +
  
  geom_vline(aes(xintercept = log(mean(actual)), 
                 color = "Original Data"),
             linetype = "dashed", 
             linewidth = 1) +
  
  scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
  scale_fill_manual(name = "Tree Type", values = c("larch" = "#66C2A5",
                                                   "oak" = "#8DA0CB",
                                                   "sitka_spruce" = "#E78AC3",
                                                   "sycamore" = "#FFD92F")) +
  
  labs(title = "Density Plot of Predicted and Actual Counts",
       x = "Log of Predicted Counts", y = "Density") +
  theme(legend.position = "bottom")


#### QUANTIFY ####

## MAE ##

predicted_comparison_data_mean %>%
  filter((type == "larch" & larch_missing == TRUE) |
           (type == "oak" & oak_missing == TRUE) |
           (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
           (type == "sycamore" & sycamore_missing == TRUE)) %>%
  group_by(type, n_missing) %>%
  summarise(mae = mean(abs(predicted - actual), na.rm = TRUE))


predicted_comparison_data_median %>%
  filter((type == "larch" & larch_missing == TRUE) |
           (type == "oak" & oak_missing == TRUE) |
           (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
           (type == "sycamore" & sycamore_missing == TRUE)) %>%
  group_by(type, n_missing) %>%
  summarise(mae = mean(abs(predicted - actual), na.rm = TRUE))


## RMSE ##

predicted_comparison_data_mean %>%
  filter((type == "larch" & larch_missing == TRUE) |
           (type == "oak" & oak_missing == TRUE) |
           (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
           (type == "sycamore" & sycamore_missing == TRUE)) %>%
  group_by(type, n_missing) %>%
  summarise(rmse = sqrt(mean(predicted - actual)^2))


predicted_comparison_data_median %>%
  filter((type == "larch" & larch_missing == TRUE) |
           (type == "oak" & oak_missing == TRUE) |
           (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
           (type == "sycamore" & sycamore_missing == TRUE)) %>%
  group_by(type, n_missing) %>%
  summarise(rmse = sqrt(mean(predicted - actual)^2))


#######################.
###########################################.
## SPATIAL MULTIVARIATE GAM NIMBLE MODEL ##
###########################################.

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

setwd("multi-tree/")

source("beta-binomial_nimblefunction.R")


#######################.

## DATA ####

tree_data <- readRDS("") 


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



tree_types_data %>%
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

y = tree_types_data$total


#### TREE COUNTS ####

x = tree_types_data %>%
  select(larch:sycamore)


###################.

## RANDOMLY SELECT ROWS ####

set.seed(451810)
fit_rows <- sample(1:nrow(tree_types_data),1000,replace = FALSE)


###################.

## CONSTANTS ####

N = nrow(tree_types_data)

# NUMBER OF ROWS TO FIT
N_fit <- length(fit_rows)

# NUMBER OF TREE TYPES
N_types = length(tree_types)

# NUMBER OF BASIS
N_basis = 400
#N_basis = 200


#######################.

## GAM MODEL ####

larch_gam_bin <- gam(larch/total ~ s(X_coord, Y_coord, k=N_basis), # CHANGED
                     family = "binomial",
                     weights = tree_types_data[fit_rows,]$total,
                     data = tree_types_data[fit_rows,])

Z <- predict(larch_gam_bin,
             tree_types_data,
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
                     x_data, y,
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
      
      K1[1:(N_basis - 1), 1:(N_basis-1), k] <- S1[1:(N_basis - 1), 1:(N_basis - 1)] * lambda[1,k] + S1[1:(N_basis - 1), N_basis:((N_basis - 1)*2)] * lambda[2,k]
      
      
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
  data <- list(x = x_data)
  
  
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
    inits <- list(b = matrix(rnorm(N_basis, 0, 0.25), nrow = N_basis, ncol = N_types),
                  lambda = matrix(rep(1, 2), nrow = 2, ncol = N_types),
                  psi = matrix(c(rnorm(1,2,sd=0.25),rnorm(3,0,0.25)), nrow = 4, ncol = N_types))
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
  
  #compile_model$calculate()
  
  
  #### CONFIGURE MCMC ####
  conf_model <- configureMCMC(compile_model, 
                              monitors = c("b",
                                           "lambda",
                                           "rho",
                                           "eta",
                                           "mu",
                                           "phi",
                                           "psi"
                              ),
                              print = TRUE, 
                              useConjugacy = FALSE)
  
  
  #### SAMPLERS ####
  
  conf_model$removeSamplers(c('lambda', "psi", "b"))
  #conf_model$addSampler(target = "psi", type = "AF_slice")
  #conf_model$addSampler(target = "b", type = "AF_slice")
  for(j in 1:N_types){
   conf_model$addSampler(target = paste0("b[1:",N_basis,",",j,"]"), type = "AF_slice")
  }
  for(j in 1:N_types){
    conf_model$addSampler(target = paste0("psi[1:",4,",",j,"]"), type = "AF_slice")
  }
  for (node in 1:length(model$expandNodeNames("lambda"))) {
    conf_model$addSampler(target = model$expandNodeNames("lambda")[node], type = "slice")
  }
  
  #browser()
  
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
                     setSeed = chain_seed[X]) # Chain seed
  
  #return(samples)
  return(list(initial_values = initial_values, samples = samples))
  
}


### RUN MODEL ####
run_time <- system.time({samples <- parLapply(cl = this_cluster, 
                                              X = 1:nchains,
                                              fun = run_MCMC, 
                                              
                                              initial_seed = initial_seed,
                                              chain_seed = chain_seed,
                                              
                                              x_data = x[fit_rows, ], 
                                              y = y[fit_rows],
                                              
                                              Z = Z[fit_rows, ],
                                              S1 = S1,
                                              
                                              N = N_fit, 
                                              N_types = N_types,
                                              N_basis = N_basis,
                                              
                                              niter = 6000, 
                                              nburn = 5000, 
                                              nthin =   1)})


run_time


# CLOSE CLUSTER
stopCluster(this_cluster)


#######################.

## OUTPUT ####

save(samples, file = paste0('output/spatial_model.RData'))



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


subset_lambda <- function(chain) {
  cols <- grep("^lambda", colnames(chain), value = TRUE)
  return(chain[, cols])
}

gelman_lambda <- gelman.diag(lapply(samples_list, subset_lambda), 
                             transform = TRUE, 
                             autoburnin = FALSE, 
                             multivariate = FALSE)
gelman_lambda
gelman_lambda[[1]][, "Point est."] %>% mean()
gelman_lambda[[1]][, "Point est."] %>% median()

mean(gelman_lambda[[1]][, "Point est."]<=1.05, na.rm = TRUE)


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
         gelman_lambda[[1]][, "Point est."] %>% median(),
         gelman_phi[[1]][, "Point est."] %>% median(),
         gelman_psi[[1]][, "Point est."] %>% median()))


mean(c(mean(gelman_b[[1]][, "Point est."]<=1.05, na.rm = TRUE),
       mean(gelman_lambda[[1]][, "Point est."]<=1.05, na.rm = TRUE),
       mean(gelman_phi[[1]][, "Point est."]<=1.05, na.rm = TRUE),
       mean(gelman_psi[[1]][, "Point est."]<=1.05, na.rm = TRUE)))

#######################.

### TRACEPLOT ####

#### LAMBDA ####

traplot(samples_mcmc[, "lambda[1, 1]"], main = "lambda_1, larch")
traplot(samples_mcmc[, "lambda[2, 1]"], main = "lambda_2, larch")

traplot(samples_mcmc[, "lambda[1, 2]"], main = "lambda_1, oak")
traplot(samples_mcmc[, "lambda[2, 2]"], main = "lambda_2, oak")

traplot(samples_mcmc[, "lambda[1, 3]"], main = "lambda_1, sitka spruce")
traplot(samples_mcmc[, "lambda[2, 3]"], main = "lambda_2, sitka spruce")

traplot(samples_mcmc[, "lambda[1, 4]"], main = "lambda_1, sycamore")
traplot(samples_mcmc[, "lambda[2, 4]"], main = "lambda_2, sycamore")


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


# plot(x = mu_posterior_median$mu_posterior_median, y = (larch_data$larch[train_index] / larch_data$total[train_index]))
# abline(a=0, b = 1, col = "red")


#######################.

## MODEL CHECKING ####

set.seed(123)

### SIMULATE NEW Y ####

simulate_tree_count <- function(N,
                                y, 
                                mu, phi) {
  
  
  # INITILISE TO STORE SIMULATION
  simulated_data <- array(0, dim = c(N, N_types))
  
  # SIMULATE COUNTS FIRST TYPE
  observation <- rbinom(N, 
                        y,
                        rbeta(N, 
                              mu[,1]*phi[,1],
                              (1-mu[,1])*phi[,1]))
  
  # SAVE SIMULATED VALUES
  simulated_data[, 1] <- observation
  
  # ALL OTHER TYPES 
  for (k in 2:N_types) {
    
    if (k>2){
      observation <- rbinom(N,
                            y - apply(simulated_data[, 1:(k-1)],1,sum),
                            rbeta(N, mu[,k]*phi[,k],
                                  (1-mu[,k])*phi[,k]))
    }else{
      observation <- rbinom(N, 
                            y - simulated_data[, 1], 
                            rbeta(N, mu[,k]*phi[,k],
                                  (1-mu[,k])*phi[,k]))
      
    }
    
    # SAVE SIMULATED VALUES
    simulated_data[, k] <- observation
    
  }
  
  
  return(simulated_data)
  
}

#### SIMULATION DATA ####

#b_samples <- output[, grep("^b", colnames(output))]
b_samples <- array(NA, dim = c(N_basis, N_types, N_samples))
for (i in 1:N_samples) {
  b_row <- output[i, ] %>% select(starts_with("b")) 
  b_samples[, , i] <- array(unlist(b_row), dim = c(N_basis, N_types))
}
b_samples

#psi_samples <- select(output,starts_with("psi"))%>%as.matrix()
psi_samples <- array(NA, dim = c(N_types, N_samples))
for (i in 1:N_samples) {
  psi_row <- output[i, ] %>% select(starts_with("psi")) 
  psi_samples[, i] <- array(unlist(psi_row), dim = c(N_types))
}
psi_samples

#mu_samples <- expit(as.matrix(t(b_samples))%*%t(Z[fit_rows, ]))
mu_samples <- array(NA, dim = c(N_samples, N_fit, N_types))
for (i in 1:N_types) {
  mu_samples[, , i] <- expit(as.matrix(t(b_samples[, i, ]))%*%t(Z[fit_rows,]))
}
mu_samples

phi_samples <- array(NA, dim = c(N_samples, N_fit, N_types))
for(i in 1:N_types) {
  for(j in 1:N_fit){
    mu_mean <- apply(mu_samples[, , i],1,mean)
    phi_samples[, j, i] <- exp(psi_samples[i, 1] +
                                 psi_samples[i,2]*(mu_samples[,j,i]-mu_mean) +
                                 psi_samples[i,3]*(mu_samples[,j,i]-mu_mean)^2 +
                                 psi_samples[i,4]*(mu_samples[,j,i]-mu_mean)^3)
  }
}
phi_samples


#### RUN SIMULATION ####

simulation_replicates <- array(NA, dim = c(N_fit, N_types, N_samples))
for (i in 1:N_samples){
  simulation_replicates[, , i] <- simulate_tree_count(N = N_fit,
                                                      y = y[fit_rows],
                                                      mu = mu_samples[i, ,], 
                                                      phi = phi_samples[i, ,])
}


simulation_replicates


#### SAVE REPLICATES ####


save(simulation_replicates, file = 'replicates_polynomial.RData')


### DENSITY ####

#### MEAN ####

for (i in 1:N_types) {
  
  plot(simulation_replicates[, i, ] %>%
         apply(.,2, mean) %>%
         as.data.frame() %>%
         rename(mean = ".") %>%
         ggplot(., aes(x = mean, fill = type)) +
         
         geom_density(aes(color = "Replicates", fill = "Replicates"), alpha = 0.8) +
         
         geom_vline(data = tree_types_data[fit_rows, ] %>%
                      select(tree_types[i]) %>%
                      mutate(mean =  apply(.,2,mean)), 
                    aes(xintercept = mean, color = "Original Data"), 
                    linetype = "dashed", linewidth = 1) +
         
         scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
         scale_fill_manual(name = "", values = c("Replicates" = "#66C2A5")) +
         labs(title = paste0("Posterior Density Plot of the Mean of ", stringr::str_to_title(stringr::str_replace(tree_types[i], "_", " ")))) +
         theme(legend.position = "bottom"))
}


pdf("output/posterior_density_mean.pdf")
for (i in 1:N_types) {
  
  plot(simulation_replicates[, i, ] %>%
         apply(.,2, mean) %>%
         as.data.frame() %>%
         rename(mean = ".") %>%
         ggplot(., aes(x = mean, fill = type)) +
         
         geom_density(aes(color = "Replicates", fill = "Replicates"), alpha = 0.8) +
         
         geom_vline(data = tree_types_data[fit_rows, ] %>%
                      select(tree_types[i]) %>%
                      mutate(mean =  apply(.,2,mean)), 
                    aes(xintercept = mean, color = "Original Data"), 
                    linetype = "dashed", linewidth = 1) +
         
         scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
         scale_fill_manual(name = "", values = c("Replicates" = "#66C2A5")) +
         labs(title = paste0("Posterior Density Plot of the Mean of ", stringr::str_to_title(stringr::str_replace(tree_types[i], "_", " ")))) +
         theme(legend.position = "bottom"))
}
dev.off()


for (i in 1:N_types) {
  pdf(paste0("output/posterior_density_mean-", tree_types[i],".pdf"))
  
  plot(simulation_replicates[, i, ] %>%
         apply(.,2, mean) %>%
         as.data.frame() %>%
         rename(mean = ".") %>%
         ggplot(., aes(x = mean, fill = type)) +
         
         geom_density(aes(color = "Replicates", fill = "Replicates"), alpha = 0.8) +
         
         geom_vline(data = tree_types_data[fit_rows, ] %>%
                      select(tree_types[i]) %>%
                      mutate(mean =  apply(.,2,mean)), 
                    aes(xintercept = mean, color = "Original Data"), 
                    linetype = "dashed", linewidth = 1) +
         
         scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
         scale_fill_manual(name = "", values = c("Replicates" = "#66C2A5")) +
         labs(title = paste0("Posterior Density Plot of the Mean of ", stringr::str_to_title(stringr::str_replace(tree_types[i], "_", " ")))) +
         theme(legend.position = "bottom"))
  dev.off()
}


#### STANDARD DEVIATION ####

for (i in 1:N_types) {
  
  plot(simulation_replicates[, i, ] %>%
         apply(.,2, sd) %>%
         as.data.frame() %>%
         rename(sd = ".") %>%
         ggplot(., aes(x = sd, fill = type)) +
         
         geom_density(aes(color = "Replicates", fill = "Replicates"), alpha = 0.8) +
         
         geom_vline(data = tree_types_data[fit_rows, ] %>%
                      select(tree_types[i]) %>%
                      mutate(sd =  apply(.,2,sd)), 
                    aes(xintercept = sd, color = "Original Data"), 
                    linetype = "dashed", linewidth = 1) +
         
         scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
         scale_fill_manual(name = "", values = c("Replicates" = "#66C2A5")) +
         labs(title = paste0("Posterior Density Plot of the Standard Deviation of ", stringr::str_to_title(stringr::str_replace(tree_types[i], "_", " ")))) +
         theme(legend.position = "bottom"))
}


pdf("output/posterior_density_sd.pdf")
for (i in 1:N_types) {
  
  plot(simulation_replicates[, i, ] %>%
         apply(.,2, sd) %>%
         as.data.frame() %>%
         rename(sd = ".") %>%
         ggplot(., aes(x = sd, fill = type)) +
         
         geom_density(aes(color = "Replicates", fill = "Replicates"), alpha = 0.8) +
         
         geom_vline(data = tree_types_data[fit_rows, ] %>%
                      select(tree_types[i]) %>%
                      mutate(sd =  apply(.,2,sd)), 
                    aes(xintercept = sd, color = "Original Data"), 
                    linetype = "dashed", linewidth = 1) +
         
         scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
         scale_fill_manual(name = "", values = c("Replicates" = "#66C2A5")) +
         labs(title = paste0("Posterior Density Plot of the Standard Deviation of ", stringr::str_to_title(stringr::str_replace(tree_types[i], "_", " ")))) +
         theme(legend.position = "bottom"))
}
dev.off()


for (i in 1:N_types) {
  pdf(paste0("output/posterior_density_sd-", tree_types[i],".pdf"))
  
  plot(simulation_replicates[, i, ] %>%
         apply(.,2, sd) %>%
         as.data.frame() %>%
         rename(sd = ".") %>%
         ggplot(., aes(x = sd, fill = type)) +
         
         geom_density(aes(color = "Replicates", fill = "Replicates"), alpha = 0.8) +
         
         geom_vline(data = tree_types_data[fit_rows, ] %>%
                      select(tree_types[i]) %>%
                      mutate(sd =  apply(.,2,sd)), 
                    aes(xintercept = sd, color = "Original Data"), 
                    linetype = "dashed", linewidth = 1) +
         
         scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
         scale_fill_manual(name = "", values = c("Replicates" = "#66C2A5")) +
         labs(title = paste0("Posterior Density Plot of the Standard Deviation of ", stringr::str_to_title(stringr::str_replace(tree_types[i], "_", " ")))) +
         theme(legend.position = "bottom"))
  dev.off()
}



### QUANTILE ####

data_quantiles <- tree_types_data %>%
  select(larch:sycamore) %>%
  pivot_longer(cols = everything(), values_to = "count", names_to = "type") %>%
  group_by(type) %>%
  summarise(quantile = seq(0.01, 0.99, by = 0.01),
            value = quantile(count, probs = seq(0.01, 0.99, by = 0.01))) %>%
  ungroup()


quantiles <- vector("list", length = N_types)
for (i in 1:N_types) {
  
  quantiles[[i]] <- simulation_replicates[,i,] %>%
    apply(., 2, function(x) quantile(x, probs = seq(0.01, 0.99, by = 0.01))) %>%
    as.data.frame(.) %>%
    tibble::rownames_to_column(var = "quantile") %>%
    pivot_longer(., cols = starts_with("V"), names_to = "sample", values_to = "x") %>%
    mutate(sample = gsub("^V", "", sample)) %>%
    mutate(quantile = gsub("%", "", quantile),
           quantile = as.numeric(quantile) / 100) %>%
    mutate(type = tree_types[i])
}

quantile_data <- bind_rows(quantiles)


quantile_summary <- quantile_data %>%
  group_by(type, quantile) %>%
  summarise(lower = quantile(x, 0.025),
            mean = mean(x),
            median = median(x),
            upper = quantile(x, 0.975))


for (i in 1:N_types) {
  
  plot(quantile_summary %>%
         filter(type == tree_types[i]) %>%
         
         ggplot(., aes(x = quantile, y = median)) +
         
         geom_point(aes(shape = "Replicates Median"), size = 2) +
         
         geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95% Interval"), alpha = 0.2) +
         
         geom_line(data = data_quantiles %>%
                     filter(type == tree_types[i]), 
                   aes(x = quantile, 
                       y = value, col = "Original Data")) +
         
         scale_shape_manual(name = "", values = c("Replicates Median" = 1)) +
         scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
         scale_fill_manual(name = "", values = c("95% Interval" = "#8DA0CB")) +
         theme(legend.position = "bottom") +
         labs(y = "value") +
         ggtitle(paste0("Posterior Quantile Plot of ", stringr::str_to_title(stringr::str_replace(tree_types[i], "_", " ")))) 
  )
}


pdf("output/posterior_quantile.pdf")
for (i in 1:N_types) {
  
  plot(quantile_summary %>%
         filter(type == tree_types[i]) %>%
         
         ggplot(., aes(x = quantile, y = median)) +
         
         geom_point(aes(shape = "Replicates Median"), size = 2) +
         
         geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95% Interval"), alpha = 0.2) +
         
         geom_line(data = data_quantiles %>%
                     filter(type == tree_types[i]), 
                   aes(x = quantile, 
                       y = value, col = "Original Data")) +
         
         scale_shape_manual(name = "", values = c("Replicates Median" = 1)) +
         scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
         scale_fill_manual(name = "", values = c("95% Interval" = "#8DA0CB")) +
         theme(legend.position = "bottom") +
         labs(y = "value") +
         ggtitle(paste0("Posterior Quantile Plot of ", stringr::str_to_title(stringr::str_replace(tree_types[i], "_", " ")))) 
  )
}
dev.off()


for (i in 1:N_types) {
  pdf(paste0("output/posterior_quantile-", tree_types[i],".pdf"))
  
  plot(quantile_summary %>%
         filter(type == tree_types[i]) %>%
         
         ggplot(., aes(x = quantile, y = median)) +
         
         geom_point(aes(shape = "Replicates Median"), size = 2) +
         
         geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95% Interval"), alpha = 0.2) +
         
         geom_line(data = data_quantiles %>%
                     filter(type == tree_types[i]), 
                   aes(x = quantile, 
                       y = value, col = "Original Data")) +
         
         scale_shape_manual(name = "", values = c("Replicates Median" = 1)) +
         scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
         scale_fill_manual(name = "", values = c("95% Interval" = "#8DA0CB")) +
         theme(legend.position = "bottom") +
         labs(y = "value") +
         ggtitle(paste0("Posterior Quantile Plot of ", stringr::str_to_title(stringr::str_replace(tree_types[i], "_", " ")))) 
  )
  dev.off()
}


#######################.

####################################.
## SPATIAL LARCH GAM NIMBLE MODEL ##
####################################.

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

setwd("larch")

source("beta-binomial_nimblefunction.R")


#######################.

## DATA ####

tree_data <- readRDS("") 


#### SPATIAL LOCATIONS ####

spatial_locations <- tree_data %>%
  select(X_coord, Y_coord)


### LARCH ####

larch_data <- tree_data %>% 
  select(X_coord,
         Y_coord,
         larch,
         total) 


larch_data %>%
  ggplot(., aes(x = X_coord, y = Y_coord, fill = larch)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  labs(x = "x coordinate", y = "y coordinate") +
  theme(legend.position = "bottom") +
  ggtitle("Heatmap for Larch - Count")



larch_data %>%
  mutate(proportion = larch / total) %>%
  ggplot(., aes(x = X_coord, y = Y_coord, fill = proportion)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  labs(x = "x coordinate", y = "y coordinate") +
  theme(legend.position = "bottom") +
  ggtitle("Heatmap for Larch - Proportion")



#### TOTAL COUNT ####

y = larch_data$total


#### TREE COUNT ####

x = larch_data$larch


###################.

## CONSTANTS ####

N = nrow(larch_data)

# NUMBER OF BASIS
N_basis = 400


#######################.

## TRAIN / TEST DATA ####

N_train <- ceiling(N*0.5) # CHANGED
N_test <- floor(N*0.5)

set.seed(453459) # CHANGED
train_index <- sort(sample(1:N,N_train,replace=FALSE))
test_index <- setdiff(1:N,train_index)


### SUMMARY ####

train_x = x[train_index]
train_y = y[train_index]

test_x = x[test_index]
test_y = y[test_index]


N_train = length(train_index)
N_test = length(test_index)


### SUMMARY ####

mean(train_x)
mean(test_x)

median(train_x)
median(test_x)

mean(train_y)
mean(test_y)


#######################.

larch_gam_bin <- gam(larch/total ~ s(X_coord, Y_coord, k=N_basis), # CHANGED
                     family = "binomial",
                     weights = larch_data$total[train_index],
                     data = larch_data[train_index,])

Z <- predict(larch_gam_bin,
             larch_data,
             type="lpmatrix")


larch_jagam_bin <- jagam(larch/total ~ s(X_coord, Y_coord, k=N_basis), # CHANGED
                         family = "binomial",
                         weights = larch_data$total,
                         data = larch_data,
                         file = "larch.jags")

# Z <- larch_jagam_bin$jags.data$X
# Z

S1 <- larch_jagam_bin$jags.data$S1
#S1



#######################.

## NIMBLE MODEL ####

# NUMBER OF CHAINS
nchains = 4

# NUMBER OF CLUSTERS
this_cluster <- makeCluster(nchains)


### MODEL FUNCTION ####

run_MCMC <- function(seed, 
                     x_data, y,
                     Z, S1,
                     N, N_basis,
                     niter, nburn, nthin) {
  
  library(nimble)
  
  source("beta-binomial_nimblefunction.R")
  
  
  #### MODEL CODE ####
  
  spatial_model_code <- nimbleCode({
    
    phi[1] ~ dgamma(2, 0.05)
    
    ## PRIORS ##
    
    for (i in 1:2) {
      lambda[i] ~ dgamma(0.05, 0.005)
      rho[i] <- log(lambda[i])
      
    }
    
    K1[1:(N_basis - 1), 1:(N_basis-1)] <- S1[1:(N_basis - 1), 1:(N_basis - 1)] * lambda[1]  + S1[1:(N_basis - 1), N_basis:((N_basis - 1)*2)] * lambda[2]
    
    
    b[1] ~ dnorm(0, sd = 10)
    b[2:N_basis] ~ dmnorm(zero[1:(N_basis-1)], K1[1:(N_basis - 1), 1:(N_basis-1)])
    
    
    ## LINEAR PREDICTOR ##
    eta[1:N] <- Z[1:N, 1:N_basis] %*% b[1:N_basis]
    
    
    ## EXPECTED RESPONSE ##
    mu[1:N] <- expit(eta[1:N])
    
    ## MODEL COUNTS X ##
    for (i in 1:N) {
      
      x_data[i] ~ dbetabinomial(mu[i],
                                phi[1],
                                y[i])
      
    }
    
  })
  
  
  #### DATA ####
  data <- list(x_data = x_data)
  
  
  #### CONSTANTS ####
  constants <- list(y = y,
                    Z = Z,
                    S1 = S1,
                    N = N,
                    N_basis = N_basis,
                    zero = rep(0, (N_basis - 1)))
  
  
  #### INITIAL VALUES ####
  init_fun <- function(){
    inits <- list(b = rnorm(N_basis, 0, 0.1),
                  lambda = rep(1, 2),
                  phi = runif(1,5,50))
    return(inits)
  }
  
  
  
  #### BUILD MODEL ####
  model <- nimbleModel(code = spatial_model_code,
                       data = data,
                       constants = constants,
                       inits = init_fun(),
                       buildDerivs = TRUE)
  
  
  
  #### COMPILE MODEL ####
  compile_model <- compileNimble(model)
  
  
  compile_model$calculate()
  
  
  #### CONFIGURE MCMC ####
  conf_model <- configureMCMC(compile_model, 
                              monitors = c("b",
                                           "lambda",
                                           "rho",
                                           "eta",
                                           "mu",
                                           "phi"
                              ),
                              print = TRUE, 
                              useConjugacy = FALSE)
  
  
  #### SAMPLERS ####
  
  conf_model$removeSamplers(c('lambda', "b"))
  conf_model$addSampler(target = "lambda[1]", type = "slice")
  conf_model$addSampler(target = "lambda[2]", type = "slice")
  conf_model$addSampler(target = "b", type = "AF_slice")
  
  
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
                     thin =    nthin, 
                     samplesAsCodaMCMC = TRUE,
                     setSeed = seed*10000)
  
  return(samples)
  
}


### RUN MODEL ####
run_time <- system.time({samples <- parLapply(cl = this_cluster, 
                                              X = 1:nchains,
                                              fun = run_MCMC, 
                                              
                                              x_data = train_x, 
                                              y = train_y,
                                              
                                              Z = Z[train_index, ],
                                              S1 = S1,
                                              
                                              N = N_train, 
                                              N_basis = N_basis,
                                              
                                              niter = 2000, 
                                              nburn = 1000, 
                                              nthin =    1)})


# CLOSE CLUSTER
stopCluster(this_cluster)



# RUN TIME
run_time
# 19449.38 seconds


#######################.

## OUTPUT ####

save(samples, file = paste0('output/larch_mcmc_output-fixed.RData'))



# COMBINE ALL CHAINS
output <- as_tibble(do.call('rbind', samples))


samples_mcmc <- as.mcmc.list(samples)


N_samples <- nrow(output)


### GELMAN ####


subset_b <- function(chain) {
  cols <- grep("^b", colnames(chain), value = TRUE)
  return(chain[, cols])
}

gelman_b <- gelman.diag(lapply(samples, subset_b), 
                        transform = TRUE, 
                        autoburnin = FALSE, 
                        multivariate = FALSE)
gelman_b
gelman_b[[1]][, "Point est."] %>% mean()
gelman_b[[1]][, "Point est."] %>% median()


mean(gelman_b[[1]][, "Point est."]<=1.05, na.rm = TRUE)


subset_lambda <- function(chain) {
  cols <- grep("^lambda", colnames(chain), value = TRUE)
  return(chain[, cols])
}

gelman_lambda <- gelman.diag(lapply(samples, subset_lambda), 
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

gelman_phi <- gelman.diag(lapply(samples, subset_phi), 
                          transform = TRUE, 
                          autoburnin = FALSE, 
                          multivariate = FALSE)
gelman_phi
gelman_phi[[1]][, "Point est."] %>% mean()
gelman_phi[[1]][, "Point est."] %>% median()

mean(gelman_phi[[1]][, "Point est."]<=1.05, na.rm = TRUE)


median(gelman_b[[1]][, "Point est."] %>% median(),
       gelman_lambda[[1]][, "Point est."] %>% median(),
       gelman_phi[[1]][, "Point est."] %>% median())

mean(mean(gelman_b[[1]][, "Point est."]<=1.05, na.rm = TRUE),
     mean(gelman_lambda[[1]][, "Point est."]<=1.05, na.rm = TRUE),
     mean(gelman_phi[[1]][, "Point est."]<=1.05, na.rm = TRUE))


#######################.

## SPLINE VARIANCE ####

mu_samples <- array(NA, dim = c(N_train, N_samples))
for (i in 1:N_samples) {
  mu_row <- output[i, ] %>% select(starts_with("mu")) 
  mu_samples[, i] <- array(unlist(mu_row), dim = c(N_train))
}
mu_samples


mean_mu <- apply(mu_samples, 1, mean)
sd(mean_mu) 

#######################.

## TRACEPLOT - OUTPUT ####

# LAMBDA #

traplot(samples_mcmc[, "lambda[1]"], main = "lambda_1")
traplot(samples_mcmc[, "lambda[2]"], main = "lambda_2")


# PHI #

traplot(samples_mcmc[, "phi[1]"], main = "phi_1")
traplot(samples_mcmc[, "phi[2]"], main = "phi_2")
traplot(samples_mcmc[, "phi[3]"], main = "phi_3")
traplot(samples_mcmc[, "phi[4]"], main = "phi_4")
traplot(samples_mcmc[, "phi[5]"], main = "phi_5")


# B #

traplot(samples_mcmc[, "b[1]"], main = "b_1")
traplot(samples_mcmc[, "b[2]"], main = "b_2")
traplot(samples_mcmc[, "b[3]"], main = "b_3")
traplot(samples_mcmc[, "b[4]"], main = "b_4")
traplot(samples_mcmc[, "b[5]"], main = "b_5")



#######################.


plot(x = mu_posterior_median$mu_posterior_median, y = (larch_data$larch[train_index] / larch_data$total[train_index]))
abline(a=0, b = 1, col = "red")


#######################.

## HEATMAP - OUTPUT ####

mu_samples <- output[, grep("mu", colnames(output))]

mu_posterior_median <- apply(mu_samples, 2, median) %>%
  as.data.frame() %>%
  rename(mu_posterior_median = ".") %>%
  cbind(spatial_locations[train_index, ])


mu_posterior_median %>%
  ggplot(., aes(x = X_coord, y = Y_coord, fill = mu_posterior_median)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  labs(x = "x coordinate", y = "y coordinate") +
  theme(legend.position = "bottom") +
  ggtitle("Heatmap for Posterior Median - Larch")


pdf("output/heatmap_larch-posterior-train_fixed.pdf")
mu_posterior_median %>%
  ggplot(., aes(x = X_coord, y = Y_coord, fill = mu_posterior_median)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  labs(x = "x coordinate", y = "y coordinate") +
  theme(legend.position = "bottom") +
  ggtitle("Heatmap for Posterior Median - Larch")
dev.off()


#######################.

## MODEL CHECKING ####

set.seed(123)

### SIMULATE NEW Y ####

simulate_tree_count <- function(N,
                                y, 
                                mu, phi) {
  
  
  # INITILISE TO STORE SIMULATION
  simulated_data <- array(0, dim = c(N))
  
  # SIMULATE COUNTS
  # alpha <- as.numeric(mu * phi)
  # beta <- as.numeric(phi*(1 - mu))
  # beta_probabilities <- rbeta(N, alpha, beta)
  # observation <- rbinom(N, y, beta_probabilities)
  
  observation <- rbinom(N, 
                        y,
                        rbeta(N, 
                              mu*phi,
                              (1-mu)*phi))
  
  # SAVE SIMULATED VALUES
  simulated_data <- observation
  
  
  return(simulated_data)
  
}


#### SIMULATION DATA ####

# mu_samples <- output[, grep("mu", colnames(output))]
# mu_samples
# 
# phi_samples <- output[, grep("phi", colnames(output))]
# phi_samples


b_samples <- output[, grep("^b", colnames(output))]

mu_samples <- expit(as.matrix(b_samples)%*%t(Z[train_index, ]))

phi_samples <- output[, grep("phi", colnames(output))]



#### RUN SIMULATION ####

simulation_replicates <- array(NA, dim = c(N_train, N_samples))
for (i in 1:N_samples){
  simulation_replicates[, i] <- simulate_tree_count(N = N_train,
                                                    y = train_y,
                                                    mu = as.numeric(mu_samples[i, ]), 
                                                    phi = as.numeric(phi_samples[i, ]))
}


simulation_replicates


#### SAVE REPLICATES ####


save(simulation_replicates, file = 'output/replicates_fixed.RData')


### DENSITY ####

#### MEAN ####

simulation_replicates %>%
  apply(.,2, mean) %>%
  as.data.frame() %>%
  rename(mean = ".") %>%
  ggplot(., aes(x = mean, fill = type)) +
  
  geom_density(aes(color = "Replicates", fill = "Replicates"), alpha = 0.8) +
  
  geom_vline(data = larch_data[train_index, ] %>%
               summarise(mean_count = mean(larch)), 
             aes(xintercept = mean_count, color = "Original Data"), 
             linetype = "dashed", linewidth = 1) +
  
  scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
  scale_fill_manual(name = "", values = c("Replicates" = "#66C2A5")) +
  labs(title = "Density Plot of the Mean of the Larch Replicates") +
  theme(legend.position = "bottom")


pdf("output/posterior_density_mean-fixed.pdf")
simulation_replicates %>%
  apply(.,2, mean) %>%
  as.data.frame() %>%
  rename(mean = ".") %>%
  ggplot(., aes(x = mean, fill = type)) +
  
  geom_density(aes(color = "Replicates", fill = "Replicates"), alpha = 0.8) +
  
  geom_vline(data = larch_data[train_index, ] %>%
               summarise(mean_count = mean(larch)), 
             aes(xintercept = mean_count, color = "Original Data"), 
             linetype = "dashed", linewidth = 1) +
  
  scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
  scale_fill_manual(name = "", values = c("Replicates" = "#66C2A5")) +
  labs(title = "Density Plot of the Mean of the Larch Replicates") +
  theme(legend.position = "bottom")
dev.off()


#### STANDARD DEVIATION ####

simulation_replicates %>%
  apply(.,2, sd) %>%
  as.data.frame() %>%
  rename(sd = ".") %>%
  ggplot(., aes(x = sd, fill = type)) +
  
  geom_density(aes(color = "Replicates", fill = "Replicates"), alpha = 0.8) +
  
  geom_vline(data = larch_data[train_index, ] %>%
               summarise(sd_count = sd(larch)), 
             aes(xintercept = sd_count, color = "Original Data"), 
             linetype = "dashed", linewidth = 1) +
  
  scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
  scale_fill_manual(name = "", values = c("Replicates" = "#66C2A5")) +
  labs(title = "Density Plot of the Standard Deviation of the Larch Replicates") +
  theme(legend.position = "bottom")


pdf("output/posterior_density_sd-fixed.pdf")
simulation_replicates %>%
  apply(.,2, sd) %>%
  as.data.frame() %>%
  rename(sd = ".") %>%
  ggplot(., aes(x = sd, fill = type)) +
  geom_density(aes(color = "Replicates", fill = "Replicates"), alpha = 0.8) +
  geom_vline(data = larch_data[train_index, ] %>%
               summarise(sd_count = sd(larch)), 
             aes(xintercept = sd_count, color = "Original Data"), 
             linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
  scale_fill_manual(name = "", values = c("Replicates" = "#66C2A5")) +
  labs(title = "Density Plot of the Standard Deviation of the Larch Replicates") +
  theme(legend.position = "bottom")
dev.off()


### QUANTILE ####

data_quantiles <- as.data.frame(quantile(larch_data$larch, probs = seq(0.01, 0.99, by = 0.01))) %>%
  tibble::rownames_to_column(var = "quantile") %>%
  rename_with(~"value", 2) %>%
  mutate(quantile = gsub("%", "", quantile),
         quantile = as.numeric(quantile) / 100)


quantiles <- simulation_replicates %>%
  apply(., 2, function(x) quantile(x, probs = seq(0.01, 0.99, by = 0.01))) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(var = "quantile") %>%
  pivot_longer(., cols = starts_with("V"), names_to = "sample", values_to = "x") %>%
  mutate(sample = gsub("^V", "", sample)) %>%
  mutate(quantile = gsub("%", "", quantile),
         quantile = as.numeric(quantile) / 100)


quantile_summary <- quantiles %>%
  group_by(quantile) %>%
  summarise(lower = quantile(x, 0.025),
            mean = mean(x),
            median = median(x),
            upper = quantile(x, 0.975))


quantile_summary %>%
  ggplot(., aes(x = quantile, y = median)) +
  
  geom_point(aes(shape = "Replicates Median"), size = 2) +
  
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95% Interval"), alpha = 0.2) +
  
  geom_line(data = data_quantiles, 
            aes(x = quantile, 
                y = value, col = "Original Data")) +
  
  scale_shape_manual(name = "", values = c("Replicates Median" = 1)) +
  scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
  scale_fill_manual(name = "", values = c("95% Interval" = "#8DA0CB")) +
  theme(legend.position = "bottom") +
  labs(y = "value") +
  ggtitle("Quantile Plot of the Larch Replicates") 


pdf("output/quantile_plot-fixed.pdf")
quantile_summary %>%
  ggplot(., aes(x = quantile, y = median)) +
  geom_point(aes(shape = "Replicates Median"), size = 2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95% Interval"), alpha = 0.2) +
  geom_line(data = data_quantiles, 
            aes(x = quantile, 
                y = value, col = "Original Data")) +
  
  scale_shape_manual(name = "", values = c("Replicates Median" = 1)) +
  scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
  scale_fill_manual(name = "", values = c("95% Interval" = "#8DA0CB")) +
  theme(legend.position = "bottom") +
  labs(y = "value") +
  ggtitle("Quantile Plot of the Larch Replicates") 
dev.off()


#######################.

## PREDICTION ####

# b_samples <- output[, grep("^b", colnames(output))]

# eta_test <- array(NA, dim = c(N_samples, N_test))
# for (i in 1:N_samples) eta_test[i, ] <- Z[test_index, ] %*% as.numeric(b_samples[i,])
# 
# mu_test <- 1 / (1 + exp(-eta_test))
# mu_test <- exp(eta_test) / (1 + exp(eta_test))


mu_test <- expit(as.matrix(b_samples)%*%t(Z[test_index, ]))


predicted_counts <- array(NA, dim = c(N_test, N_samples))
for (i in 1:N_samples) {
  predicted_counts[,i] <- simulate_tree_count(N = N_test,
                                              y = test_y,
                                              mu = mu_test[i, ], 
                                              phi = as.numeric(phi_samples[i,]))
}
predicted_counts


#### SAVE PREDICTED COUNTS ####

save(predicted_counts, file = 'output/predicted_counts_fixed.RData')


### DENSITY ####

#### MEAN ####

predicted_counts %>%
  apply(.,2, mean) %>%
  as.data.frame() %>%
  rename(mean = ".") %>%
  ggplot(., aes(x = mean, fill = type)) +
  
  geom_density(aes(color = "Replicates", fill = "Replicates"), alpha = 0.8) +
  
  geom_vline(data = larch_data[test_index, ] %>%
               summarise(mean_count = mean(larch)), 
             aes(xintercept = mean_count, color = "Original Data"), 
             linetype = "dashed", linewidth = 1) +
  
  scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
  scale_fill_manual(name = "", values = c("Replicates" = "#66C2A5")) +
  labs(title = "Posterior Predictive Density Plot of the Mean for Larch") +
  theme(legend.position = "bottom")


pdf("output/posterior_predictive_density_mean-fixed.pdf")
predicted_counts %>%
  apply(.,2, mean) %>%
  as.data.frame() %>%
  rename(mean = ".") %>%
  ggplot(., aes(x = mean, fill = type)) +
  
  geom_density(aes(color = "Replicates", fill = "Replicates"), alpha = 0.8) +
  
  geom_vline(data = larch_data[test_index, ] %>%
               summarise(mean_count = mean(larch)), 
             aes(xintercept = mean_count, color = "Original Data"), 
             linetype = "dashed", linewidth = 1) +
  
  scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
  scale_fill_manual(name = "", values = c("Replicates" = "#66C2A5")) +
  labs(title = "Posterior Predictive Density Plot of the Mean for Larch") +
  theme(legend.position = "bottom")
dev.off()


#### STANDARD DEVIATION ####

predicted_counts %>%
  apply(.,2, sd) %>%
  as.data.frame() %>%
  rename(sd = ".") %>%
  ggplot(., aes(x = sd, fill = type)) +
  
  geom_density(aes(color = "Replicates", fill = "Replicates"), alpha = 0.8) +
  
  geom_vline(data = larch_data[test_index, ] %>%
               summarise(sd_count = sd(larch)), 
             aes(xintercept = sd_count, color = "Original Data"), 
             linetype = "dashed", linewidth = 1) +
  
  scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
  scale_fill_manual(name = "", values = c("Replicates" = "#66C2A5")) +
  labs(title = "Posterior Predictive Density Plot of the SD for Larch") +
  theme(legend.position = "bottom")


pdf("output/posterior_predictive_density_sd-fixed.pdf")
predicted_counts %>%
  apply(.,2, sd) %>%
  as.data.frame() %>%
  rename(sd = ".") %>%
  ggplot(., aes(x = sd, fill = type)) +
  geom_density(aes(color = "Replicates", fill = "Replicates"), alpha = 0.8) +
  geom_vline(data = larch_data[test_index, ] %>%
               summarise(sd_count = sd(larch)), 
             aes(xintercept = sd_count, color = "Original Data"), 
             linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
  scale_fill_manual(name = "", values = c("Replicates" = "#66C2A5")) +
  labs(title = "Posterior Predictive Density Plot of the SD for Larch") +
  theme(legend.position = "bottom")
dev.off()


### QUANTILE ####

predicted_quantiles <- predicted_counts %>%
  apply(., 2, function(x) quantile(x, probs = seq(0.01, 0.99, by = 0.01))) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(var = "quantile") %>%
  pivot_longer(., cols = starts_with("V"), names_to = "sample", values_to = "x") %>%
  mutate(sample = gsub("^V", "", sample)) %>%
  mutate(quantile = gsub("%", "", quantile),
         quantile = as.numeric(quantile) / 100)


predicted_quantile_summary <- predicted_quantiles %>%
  group_by(quantile) %>%
  summarise(lower = quantile(x, 0.025),
            mean = mean(x),
            median = median(x),
            upper = quantile(x, 0.975))


predicted_quantile_summary %>%
  ggplot(., aes(x = quantile, y = median)) +
  
  geom_point(aes(shape = "Predicted Median"), size = 2) +
  
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95% Interval"), alpha = 0.2) +
  
  geom_line(data = data_quantiles, 
            aes(x = quantile, 
                y = value, col = "Original Data")) +
  
  scale_shape_manual(name = "", values = c("Predicted Median" = 1)) +
  scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
  scale_fill_manual(name = "", values = c("95% Interval" = "#8DA0CB")) +
  theme(legend.position = "bottom") +
  labs(y = "value") +
  ggtitle("Posterior Predictive Quantile Plot for Larch") 


pdf("output/posterior_predictive_quantile_plot-fixed.pdf")
predicted_quantile_summary %>%
  ggplot(., aes(x = quantile, y = median)) +
  geom_point(aes(shape = "Predicted Median"), size = 2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95% Interval"), alpha = 0.2) +
  geom_line(data = data_quantiles, 
            aes(x = quantile, 
                y = value, col = "Original Data")) +
  
  scale_shape_manual(name = "", values = c("Predicted Median" = 1)) +
  scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
  scale_fill_manual(name = "", values = c("95% Interval" = "#8DA0CB")) +
  theme(legend.position = "bottom") +
  labs(y = "value") +
  ggtitle("Posterior Predictive Quantile Plot for Larch") 
dev.off()


#######################.

## HEATMAP - OUTPUT - TRAIN / TEST ####

eta_output <- array(NA, dim = c(N_samples, N))
for (i in 1:N_samples) eta_output[i, ] <- Z %*% as.numeric(b_samples[i,])

#mu_output <- 1 / (1 + exp(-eta_output))
#mu_output <- exp(eta_output) / (1 + exp(eta_output))
mu_output <- expit(as.matrix(b_samples)%*%t(Z))


mu_output_median <- apply(mu_output, 2, median) %>%
  as.data.frame() %>%
  rename(median = ".") %>%
  cbind(spatial_locations)


mu_output_median %>%
  ggplot(., aes(x = X_coord, y = Y_coord, fill = median)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  labs(x = "x coordinate", y = "y coordinate") +
  theme(legend.position = "bottom") +
  ggtitle("Heatmap for Posterior Median - Larch")


pdf("output/heatmap_mean_surface-fixed.pdf")
mu_output_median %>%
  ggplot(., aes(x = X_coord, y = Y_coord, fill = median)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  labs(x = "x coordinate", y = "y coordinate") +
  theme(legend.position = "bottom") +
  ggtitle("Heatmap for Posterior Median - Larch")
dev.off()


#######################.

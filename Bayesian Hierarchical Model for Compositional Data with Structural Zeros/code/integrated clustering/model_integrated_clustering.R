#########################################################.
## BAYESIAN HIERARCHICAL MODEL - INTEGRATED CLUSTERING ##
#########################################################.

##################################.

## PACKAGES ####

library(coda)
library(R6)
library(nimble)
library(dplyr)
library(tidyr)
library(tibble)
library(parallel)
library(ggplot2)
#library(GGally)
library(mcmcplots)
library(reshape2)
#library(factoextra)
library(dendextend)
library(cluster)
library(mnormt)
library(class)


##################################.

## SET WORKING DIRECTORY ####

setwd("")

##################################.

## DATA ####

glass <- read.table("", header = T)

glass$Type <- as.factor(sapply(glass$Name, function(x) substr(x, 1, 1)))

##################################.

## TRAIN / TEST DATA ####

set.seed(453459) 

kfold = 5

N_item <- length(unique(glass$Item))

N_type <- length(unique(glass$Type))

glass_items <- unique(glass$Item)

resampled_items <- sample(glass_items, length(glass_items), replace = FALSE)

data_partition <- split(resampled_items, sort(resampled_items%%5))


## TEST DATA ##

# EACH TEST SET IS TEST SET 

for (i in 1:kfold) {
  assign(paste0("test_data_", i), 
         glass %>%
           filter(Item %in% data_partition[[i]]))
}


## TRAINING DATA ##

# EACH TRAINING SET IS REST OF SET NOT IN TEST SET 

for (i in 1:kfold) {
  assign(paste0("train_data_", i), 
         glass %>%
           filter(!(Item %in% data_partition[[i]])))
}


N_train = c(nrow(train_data_1), nrow(train_data_2), nrow(train_data_3), nrow(train_data_4), nrow(train_data_5))
N_train
N_test = c(nrow(test_data_1), nrow(test_data_2), nrow(test_data_3), nrow(test_data_4), nrow(test_data_5))
N_test

N_cluster = k = 5


for (kf in 1:kfold) {
  
  assign(paste0("data_", kf),
         get(paste0("train_data_", kf)) %>%
           mutate(Type = factor(Type, levels = c("b", "c", "h", "p", "w"))) %>%
           mutate(Type = as.character(Type)) %>%
           rbind(get(paste0("test_data_", kf)) %>% mutate(Type = NA)) %>%
           arrange(Item))
}


##################################.

# NUMBER OF CHAINS
nchains = 8

# NUMBER OF CLUSTERS
this_cluster <- makeCluster(nchains)


# FUNCTION
run_MCMC <- function(X, seed,
                     item_presence_absence,
                     N, 
                     z, p, t, I, J, s, 
                     N_cluster,
                     type_index, unique_item_index, item_type_index, piece_index, 
                     niter, nburnin, nthin) {
  
  library(nimble)
  library(dplyr)
  
  # MODEL CODE #
  code <- nimbleCode({ 
    
    
    # CLUSTERS #
    p_cluster[1:N_cluster] ~ ddirch(prior_cluster[1:N_cluster])
    
    for (i in 1:I) {
      item_cluster[i] ~ dcat(p_cluster[1:N_cluster])
    }
    
    for (cl in 1:N_cluster) {
      
      for (e in 1:p) {
        
        # PROBABILITY PRESENCE / ABSENCE
        q[cl, e] ~ dunif(0, 1)
        
      }
      
      # MEAN PROBABILITY PRESENCE / ABSENCE
      #mean_q[cl] <- mean(q[cl, 1:p])
      
      
      # TYPE #
      p_type[1:t, cl] ~ ddirch(prop_type[1:t]) 
      
      
      # PRIORS #
      
      Sigma_Prior[1:p, 1:p, cl] <- diag(p)/1000
      
      Psi[1:p, 1:p, cl] ~ dwish(df = d2, B[1:p, 1:p, cl])
      
      Psi_inv[1:p, 1:p, cl] <- inverse(Psi[1:p, 1:p, cl])
      
      Sigma[1:p, 1:p, cl] ~ dwish(df = p, Sigma_Prior[1:p, 1:p, cl])
      
      Sigma_inv[1:p, 1:p, cl] <- inverse(Sigma[1:p, 1:p, cl])
      
      for (k in 1:t){
        
        theta[1:p, k, cl] ~ dmnorm(mean = mean_zero[1:p], cov = inv_phi[1:p, 1:p, cl])
        
        Omega[1:p, 1:p, k, cl] ~ dwish(df = d1[k], A[1:p, 1:p, k, cl])
        
        Omega_inv[1:p, 1:p, k, cl] <- inverse(Omega[1:p, 1:p, k, cl])
        
      }
      
      for(i in 1:I){
        
        for (e in 1:p) {
          
          u[i, e] ~ dbern(q[item_cluster[i], e])
          
        }
        
        b[i, 1:p] ~ dmnorm(theta[1:p, item_type[i], cl], Omega_inv[1:p, 1:p, item_type[i], cl]) 
        
        for(j in 1:J){
          
          c[i, j, 1:p] ~ dmnorm(b[i, 1:p], Psi[1:p, 1:p, cl])
          
        }
      }
    }
    
    for (i in 1:I){
      item_type[i] ~ dcat(p_type[1:t, item_cluster[i]])
    }
    
    for(i in 1:N){
      z[i, 1:p] ~ dmnorm(mean = c[item[i], piece[i], 1:p], Sigma[1:p, 1:p, item_cluster[item[i]]])
    }
    
  })
  
  # DATA #
  data <- list(z = z,
               u = item_presence_absence,
               item_type = item_type_index
  )
  
  # CONSTANTS #
  constants <- list(p = p,
                    t = t,
                    I = I,
                    J = J,
                    
                    N = N,
                    
                    piece = piece_index,
                    #type = type_index,
                    item = unique_item_index,
                    
                    N_cluster = N_cluster,
                    
                    mean_zero = rep(0,p),
                    inv_phi = array(s * diag(1, p), dim = c(p, p, N_cluster)),
                    
                    A = array(diag(p)/1000, dim=c(p, p, t, N_cluster)),
                    B = array(diag(p)/1000, dim=c(p, p, N_cluster)),
                    d1 = rep(p, t),
                    d2 = p,
                    
                    prop_type = rep(1/t, t),
                    prior_cluster = rep(1/N_cluster, N_cluster)
  )
  
  # INITIAL VALUES #
  init_fun <- function(){
    
    item_type_inits <- rep(NA, I)
    for (i in 1:I) {
      if(is.na(item_type_index[i])) {item_type_inits[i] <- sample(1:t, 1, replace = T)}
    }
    
    inits <- list(theta = array(0.001, dim = c(p, t, N_cluster)),
                  Sigma = array(diag(1, p), dim = c(p, p, N_cluster)),
                  Omega = array(diag(1, p), dim = c(p, p, t, N_cluster)),
                  Psi = array(diag(1, p), dim = c(p, p, N_cluster)),
                  b = matrix(rnorm(I*p, 0, 0.1), nrow = I, ncol = p),
                  c = array(rnorm(I*J*p, 0, 0.1), dim = c(I, J, p)),
                  
                  p_type = matrix(rep(1/t, N_cluster), nrow = t, ncol = N_cluster),
                  item_type = item_type_inits,
                  
                  p_cluster = rep(1/N_cluster, N_cluster),
                  item_cluster = sample(1:N_cluster, I, replace = T),
                  
                  q = matrix(runif(N_cluster*p, 0, 1), nrow = N_cluster, ncol = p) %>%
                    as.data.frame() %>%
                    mutate(mean = rowMeans(.)) %>%
                    arrange(desc(mean)) %>%
                    select(-mean) %>%
                    as.matrix(),
                  mean_q = sort(runif(N_cluster, 0, 1), decreasing = TRUE)
    )
    return(inits)
  }
  
  
  # BUILD MODEL #
  model <- nimbleModel(code = code,
                       data = data,
                       constants = constants,
                       inits = init_fun())
  
  
  # COMPILE MODEL #
  compile_model <- compileNimble(model)
  compile_model$calculate()
  
  # CONFIGURE MCMC #
  confmcmc_model <- configureMCMC(compile_model, print = TRUE, 
                                  monitors = c("theta", 
                                               "Sigma", 
                                               "Sigma_inv",
                                               "Omega", 
                                               "Omega_inv",
                                               "Psi",
                                               "Psi_inv",
                                               "item_type",
                                               "p_type",
                                               "p_cluster",
                                               "item_cluster",
                                               "q",
                                               "mean_q",
                                               "u"),
                                  useConjugacy = FALSE)
  
  
  # BUILD MCMC #
  mcmc_model <- buildMCMC(confmcmc_model)
  
  # COMPILE MCMC #
  compile_mcmc_model <- compileNimble(mcmc_model)
  
  
  # RUN #
  samples <- runMCMC(compile_mcmc_model, 
                     niter = niter, 
                     nburnin = nburnin,
                     thin = nthin, 
                     nchains = 1,
                     samplesAsCodaMCMC = TRUE,
                     setSeed = seed[X])
  
  return(samples)
  
}


##################################.

## RUN KFOLD ## 

elements <- c("Na", "Mg", "Alu", "Si", "K", "Ca", "Fe")

kfold_results <- vector("list", kfold)

run_time <- system.time({
  
  for (kf in 1:kfold) {
    
    model_data <- get(paste0("data_", kf))
    
    model_data_indicator <- model_data %>%
      select(Item, O:Fe) %>%
      mutate(across(Na:Fe, ~ .x/O)) %>%
      select(-O) %>%
      group_by(Item) %>%
      summarise_all(mean) %>%
      mutate(across(Na:Fe, ~if_else(.==0, 0, 1))) %>%
      ungroup() %>%
      column_to_rownames(var = "Item")
    
    ##################################.
    
    ## MODEL DATA ####
    
    comps <- model_data %>%
      mutate(across(Na:Fe, ~ sqrt(.x/O))) %>% 
      select(Na:Fe, Item, Type, piece) %>%
      rename(Piece = piece)
    
    z = comps %>%
      select(Na:Fe) %>%
      select_if(colSums(.) != 0)
    
    
    ##################################.
    
    ### CONSTANTS ####
    
    p = length(elements)
    N = nrow(comps)
    t = N_type
    I = length(levels(as.factor(model_data$Item)))
    J = length(levels(as.factor(model_data$piece)))
    K = 3
    
    s = 1000
    
    
    type_index <- recode(comps$Type,
                         "b" = 1,
                         "c" = 2,
                         "h" = 3,
                         "p" = 4,
                         "w" = 5,
                         .default = NA_real_) 
    
    unique_items <- unique(comps$Item)
    unique_item_index <- sapply(comps$Item, function(x) which(unique_items==x))
    
    unique_item_df <- select(comps, Type, Item) %>% unique.data.frame()
    item_type_index = recode(unique_item_df$Type,
                             "b" = 1,
                             "c" = 2,
                             "h" = 3,
                             "p" = 4,
                             "w" = 5,
                             .default = NA_real_) 
    
    piece_index <- as.numeric(comps$Piece)
    
    
    ##################################.
    
    ## RUN FUNCTION ##
    output <- parLapply(cl = this_cluster, 
                        X = 1:nchains,
                        fun = run_MCMC, 
                        seed = round(runif(nchains, 5468, 5632751)),
                        
                        item_presence_absence = model_data_indicator,
                        
                        N = N, 
                        
                        z = z, 
                        p = p, 
                        t = t, 
                        I = I, 
                        J = J, 
                        s = s, 
                        
                        type_index = type_index,
                        unique_item_index = unique_item_index,
                        item_type_index = item_type_index,
                        piece_index = piece_index,
                        
                        N_cluster = N_cluster,
                        
                        niter =   400000, 
                        nburnin = 300000, 
                        nthin =      100)
    
    
    kfold_results[[kf]] <- list(output = output)
    
  }
})


# CLOSE CLUSTER
stopCluster(this_cluster)


# RUN TIME
run_time


# EXTRACT OUTPUT
mcmc_output <- lapply(kfold_results, function(x) {as.mcmc.list(x$output)})


## SAVE OUTPUT ##
save(mcmc_output, file=paste0('mcmc_output.RData'))
save.image(file=paste0('mcmc_workspace.RData'))


# COMBINE ALL CHAINS
samples <- lapply(mcmc_output, function(mcmc_list) {
  as_tibble(do.call('rbind', mcmc_list))
})


N_samples <- nrow(samples[[1]])


# NUMBER OF DRAWS
n_draws = nrow(samples[[1]])
n_draws_chain = n_draws / nchains


##################################.

## OUTPUT ####

subset_theta <- function(chain) {
  cols <- grep("^theta", colnames(chain), value = TRUE)
  return(chain[, cols])
}

gelman_theta <- gelman.diag(lapply(mcmc_output[[1]], subset_theta), 
                            transform = TRUE, 
                            autoburnin = FALSE, 
                            multivariate = FALSE)
gelman_theta
gelman_theta[[1]][, "Point est."] %>% mean()
gelman_theta[[1]][, "Point est."] %>% median()

mean(gelman_theta[[1]][, "Point est."]<=1.05, na.rm = TRUE)


subset_Sigma <- function(chain) {
  cols <- grep("^Sigma", colnames(chain), value = TRUE)
  return(chain[, cols])
}

gelman_Sigma <- gelman.diag(lapply(mcmc_output[[1]], subset_Sigma), 
                            transform = TRUE, 
                            autoburnin = FALSE, 
                            multivariate = FALSE)
gelman_Sigma
gelman_Sigma[[1]][, "Point est."] %>% mean()
gelman_Sigma[[1]][, "Point est."] %>% median()

mean(gelman_Sigma[[1]][, "Point est."]<=1.05, na.rm = TRUE)


subset_Omega <- function(chain) {
  cols <- grep("^Omega", colnames(chain), value = TRUE)
  return(chain[, cols])
}

gelman_Omega <- gelman.diag(lapply(mcmc_output[[1]], subset_Omega), 
                            transform = TRUE, 
                            autoburnin = FALSE, 
                            multivariate = FALSE)
gelman_Omega
gelman_Omega[[1]][, "Point est."] %>% mean()
gelman_Omega[[1]][, "Point est."] %>% median()

mean(gelman_Omega[[1]][, "Point est."]<=1.05, na.rm = TRUE)


subset_Psi <- function(chain) {
  cols <- grep("^Psi", colnames(chain), value = TRUE)
  return(chain[, cols])
}

gelman_Psi <- gelman.diag(lapply(mcmc_output[[1]], subset_Psi), 
                          transform = TRUE, 
                          autoburnin = FALSE, 
                          multivariate = FALSE)
gelman_Psi
gelman_Psi[[1]][, "Point est."] %>% mean()
gelman_Psi[[1]][, "Point est."] %>% median()

mean(gelman_Psi[[1]][, "Point est."]<=1.05, na.rm = TRUE)


subset_p_type <- function(chain) {
  cols <- grep("^p_type", colnames(chain), value = TRUE)
  return(chain[, cols])
}

gelman_p_type <- gelman.diag(lapply(mcmc_output[[1]], subset_p_type), 
                             transform = TRUE, 
                             autoburnin = FALSE, 
                             multivariate = FALSE)
gelman_p_type
gelman_p_type[[1]][, "Point est."] %>% mean()
gelman_p_type[[1]][, "Point est."] %>% median()

mean(gelman_p_type[[1]][, "Point est."]<=1.05, na.rm = TRUE)


subset_p_cluster <- function(chain) {
  cols <- grep("^p_cluster", colnames(chain), value = TRUE)
  return(chain[, cols])
}

gelman_p_cluster <- gelman.diag(lapply(mcmc_output[[1]], subset_p_cluster), 
                                transform = TRUE, 
                                autoburnin = FALSE, 
                                multivariate = FALSE)
gelman_p_cluster
gelman_p_cluster[[1]][, "Point est."] %>% mean()
gelman_p_cluster[[1]][, "Point est."] %>% median()

mean(gelman_p_cluster[[1]][, "Point est."]<=1.05, na.rm = TRUE)


subset_q <- function(chain) {
  cols <- grep("^q", colnames(chain), value = TRUE)
  return(chain[, cols])
}

gelman_q <- gelman.diag(lapply(mcmc_output[[1]], subset_q), 
                        transform = TRUE, 
                        autoburnin = FALSE, 
                        multivariate = FALSE)
gelman_q
gelman_q[[1]][, "Point est."] %>% mean()
gelman_q[[1]][, "Point est."] %>% median()

mean(gelman_q[[1]][, "Point est."]<=1.05, na.rm = TRUE)


subset_mean_q <- function(chain) {
  cols <- grep("^mean_q", colnames(chain), value = TRUE)
  return(chain[, cols])
}

gelman_mean_q <- gelman.diag(lapply(mcmc_output[[1]], subset_mean_q), 
                             transform = TRUE, 
                             autoburnin = FALSE, 
                             multivariate = FALSE)
gelman_mean_q
gelman_mean_q[[1]][, "Point est."] %>% mean()
gelman_mean_q[[1]][, "Point est."] %>% median()

mean(gelman_mean_q[[1]][, "Point est."]<=1.05, na.rm = TRUE)


median(c(gelman_theta[[1]][, "Point est."] %>% median(),
         gelman_Omega[[1]][, "Point est."] %>% median(),
         gelman_Sigma[[1]][, "Point est."] %>% median(),
         gelman_Psi[[1]][, "Point est."] %>% median(),
         gelman_p_type[[1]][, "Point est."] %>% median(),
         gelman_p_cluster[[1]][, "Point est."] %>% median(),
         gelman_q[[1]][, "Point est."] %>% median(),
         gelman_mean_q[[1]][, "Point est."] %>% median()))


mean(c(mean(gelman_theta[[1]][, "Point est."]<=1.05, na.rm = TRUE),
       mean(gelman_Omega[[1]][, "Point est."]<=1.05, na.rm = TRUE),
       mean(gelman_Sigma[[1]][, "Point est."]<=1.05, na.rm = TRUE),
       mean(gelman_Psi[[1]][, "Point est."]<=1.05, na.rm = TRUE),
       mean(gelman_p_type[[1]][, "Point est."]<=1.05, na.rm = TRUE),
       mean(gelman_p_cluster[[1]][, "Point est."]<=1.05, na.rm = TRUE),
       mean(gelman_q[[1]][, "Point est."]<=1.05, na.rm = TRUE),
       mean(gelman_mean_q[[1]][, "Point est."]<=1.05, na.rm = TRUE)))



## ITEM TYPE ##
subset_item_type <- function(chain) {
  cols <- grep("^item_type", colnames(chain), value = TRUE)
  return(chain[, cols])
}

gelman_item_type <- gelman.diag(lapply(mcmc_output[[1]], subset_item_type), 
                                transform = TRUE, 
                                autoburnin = FALSE, 
                                multivariate = FALSE)

gelman_item_type
gelman_item_type[[1]][, "Point est."] %>% mean(., na.rm = T)
gelman_item_type[[1]][, "Point est."] %>% median(., na.rm = T)
mean(gelman_item_type[[1]][, "Point est."]<=1.05, na.rm = TRUE)


for (kf in 1:kfold) {
  
  gelman_item_type <- gelman.diag(lapply(mcmc_output[[kf]], subset_item_type), 
                                  transform = TRUE, 
                                  autoburnin = FALSE, 
                                  multivariate = FALSE)
  
  assign(paste0("test_items_", kf), get(paste0("data_", kf)) %>%
           filter(is.na(Type)) %>%
           select(Item) %>%
           unique()
  )
  
  print(gelman_item_type[[1]] %>%
          as.data.frame() %>%
          mutate(Item = row_number()) %>%
          filter(Item %in% as.vector(unlist(get(paste0("test_items_", kf))))) %>%
          select(Item, gelman = `Point est.`) %>%
          summarise(mean = mean(gelman, na.rm = T),
                    median = median(gelman, na.rm = T),
                    mean_less_1.05 = mean(gelman<=1.05, na.rm = T))
  )
  
}



##################################.

## TRACEPLOTS ####

for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "theta[1, 1, 1]"], main = paste0("theta_1, Na, cluster 1 - kfold ", kf)))
for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "theta[2, 1, 1]"], main = paste0("theta_1, Mg, cluster 1 - kfold ", kf)))
for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "theta[3, 1, 1]"], main = paste0("theta_1, Alu, cluster 1 - kfold ", kf)))
for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "theta[4, 1, 1]"], main = paste0("theta_1, Si, cluster 1 - kfold ", kf)))
for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "theta[5, 1, 1]"], main = paste0("theta_1, K, cluster 1 - kfold ", kf)))
for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "theta[6, 1, 1]"], main = paste0("theta_1, Ca, cluster 1 - kfold ", kf)))
for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "theta[7, 1, 1]"], main = paste0("theta_1, Fe, cluster 1 - kfold ", kf)))


for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "p_type[1, 1]"], main = paste0("bulb, cluster 1 - kfold ", kf)))
for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "p_type[2, 1]"], main = paste0("car window, cluster 1 - kfold ", kf)))
for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "p_type[3, 1]"], main = paste0("headlamp, cluster 1 - kfold ", kf)))
for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "p_type[4, 1]"], main = paste0("container, cluster 1 - kfold ", kf)))
for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "p_type[5, 1]"], main = paste0("building window, cluster 1 - kfold ", kf)))


for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "q[1, 1]"]))
for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "q[1, 2]"]))
for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "q[1, 3]"]))
for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "q[1, 4]"]))
for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "q[1, 5]"]))
for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "q[1, 6]"]))
for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "q[1, 7]"]))

for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "q[2, 1]"]))
for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "q[2, 2]"]))
for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "q[2, 3]"]))
for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "q[2, 4]"]))
for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "q[2, 5]"]))
for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "q[2, 6]"]))
for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "q[2, 7]"]))


for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "p_type[1, 1]"]))
for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "p_type[1, 2]"]))
for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "p_type[1, 3]"]))
for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "p_type[1, 4]"]))
for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "p_type[1, 5]"]))

for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "p_cluster[1]"]))
for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "p_cluster[2]"]))
for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "p_cluster[3]"]))
for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "p_cluster[4]"]))
for (kf in 1:kfold) print(traplot(mcmc_output[[kf]][, "p_cluster[5]"], main = ""))

for (i in 1:N_cluster) {
  pdf(paste0("output/traceplot_p_cluster_", i,".pdf"), height = 2)
  layout(matrix(1:5, nrow = 1, ncol = 5)) 
  traplot(mcmc_output[[1]][, paste0("p_cluster[",i,"]")], main = "")
  dev.off()
}

##################################.

mode <- function(x) {
  uniq_x <- unique(x)
  uniq_x[which.max(tabulate(match(x, uniq_x)))]
}

##################################.

## CLUSTERS ####

cluster_samples <- array(NA, dim = c(I, N_samples, kfold))
for (kf in 1:kfold){
  for (i in 1:N_samples) {
    cluster_row <- samples[[kf]][i, ] %>% select(starts_with("item_cluster")) 
    cluster_samples[, i, kf] <- array(unlist(cluster_row), dim = c(I))
  }
}
cluster_samples

cluster_mode_kfold <- apply(cluster_samples, c(1,3), mode)
cluster_mode <- apply(cluster_samples, 1, mode)


## CHECK NUMBER OF CLUSTERS PER ITEM ##

cluster_mode_check <- cluster_mode_kfold %>%
  as.data.frame() %>%
  mutate(n_different = apply(., 1, function(row) length(unique(row))))

hist(cluster_mode_check$n_different,
     breaks = seq(0.5, max(cluster_mode_check$n_different) + 0.5, by = 1))
min(cluster_mode_check$n_different)
max(cluster_mode_check$n_different)



for (kf in 1:kfold) {
  
  assign(paste0("clustered_data_", kf),
         get(paste0("data_", kf)) %>%
           left_join(cluster_mode_kfold[, kf] %>%
                       as.data.frame() %>%
                       cbind(Item = unique(get(paste0("data_", kf))$Item)), by = "Item") %>%
           rename(cluster = ".")
  )
  
  assign(paste0("cluster_train_types_", kf), 
         get(paste0("clustered_data_", kf)) %>%
           filter(!is.na(Type)) %>%
           group_by(cluster, Type) %>%
           count(.) %>%
           mutate(n = n/12) %>%
           pivot_wider(names_from = cluster, values_from = n, values_fill = 0) %>%
           mutate(Type = factor(Type, levels = c("b", "c", "h", "p", "w"))) %>%
           arrange(Type) %>%
           ungroup() %>%
           
           filter(!is.na(Type)) %>%
           
           bind_rows(summarise(., Type = "total", across(where(is.numeric), sum, na.rm = TRUE)))
  )
  
  print(get(paste0("cluster_train_types_", kf)))
  
}



### CLUSTER MEANS ##

cluster_mean_samples <- array(NA, dim = c(N_cluster, N_samples, kfold))
for (kf in 1:kfold){
  for (i in 1:N_samples) {
    cluster_row <- samples[[kf]][i, ] %>% select(starts_with("mean_q")) 
    cluster_mean_samples[, i, kf] <- array(unlist(cluster_row), dim = c(N_cluster))
  }
}
cluster_mean_samples


cluster_mean_samples[, 1,]

cluster_means <- apply(cluster_mean_samples, c(1,3), mean)
cluster_means 



### CLUSTER PRESENCE / ABSENCE PROBABILITY ##

cluster_indicator_samples <- array(NA, dim = c(N_cluster, p, N_samples, kfold))
for (kf in 1:kfold){
  for (i in 1:N_samples) {
    cluster_row <- samples[[kf]][i, ] %>% select(starts_with("q")) 
    cluster_indicator_samples[, , i, kf] <- array(unlist(cluster_row), dim = c(N_cluster, p))
  }
}
cluster_indicator_samples


cluster_indicator <- apply(cluster_indicator_samples, c(1,2,4), mode)
cluster_indicator
dim(cluster_indicator)

cluster_indicator[,,1]
rowMeans(cluster_indicator[,,1])
cluster_indicator %>% apply(., c(1,3), mean)


##################################.

## PREDICTED TYPES ####

type_samples <- array(NA, dim = c(I, N_samples, kfold))
for (kf in 1:kfold){
  for (i in 1:N_samples) {
    type_row <- samples[[kf]][i, ] %>% select(starts_with("item_type")) 
    type_samples[, i, kf] <- array(unlist(type_row), dim = c(I))
  }
}
type_samples


saveRDS(type_samples, 'output/item_type_samples.rds')


type_samples[,,1] %>%
  as.data.frame() %>%
  mutate(item = row_number()) %>%
  filter(item %in% c(1:12)) %>%
  melt(., id.vars = "item", variable.name = "Item", value.name = "value") %>%
  select(-Item) %>%
  ggplot(., aes(x = value)) +
  geom_density() +
  facet_wrap(~item)



item_types_mode <- apply(type_samples, c(1,3), mode)


predicted_types_list <- vector("list", kfold)
predicted_types_list_test <- vector("list", kfold)
for (kf in 1:kfold) {
  
  predicted_types_list[[kf]] <- get(paste0("train_data_", kf)) %>%
    rbind(get(paste0("test_data_", kf))) %>%
    select(Item, Type) %>%
    unique() %>%
    mutate(Type = as.numeric(Type)) %>%
    cbind(classified = item_types_mode[, kf])
  
  predicted_types_list_test[[kf]] <- get(paste0("test_data_", kf)) %>%
    select(Item, Type) %>%
    unique() %>%
    cbind(classified = item_types_mode[unique(get(paste0("test_data_", kf))$Item), kf]) %>% 
    mutate(Type = as.numeric(Type))
  
} 
predicted_types <- do.call(rbind, predicted_types_list)
predicted_types_test <- do.call(rbind, predicted_types_list_test)


correct_classification_test <- predicted_types_test %>%
  group_by(Item, Type) %>%
  summarise(correct = if_else(classified == Type, T, F)) %>%
  group_by(Type) %>%
  summarise(
    correct = sum(correct == T),         
    total = n()
  ) %>%
  mutate(classification_rate = (correct / total) * 100) %>% 
  ungroup()
correct_classification_test

correct_classification_test %>%
  summarise(correct = sum(correct),
            total = sum(total)) %>%
  summarise(classification_rate = (correct / total) * 100)


## ACROSS SAMPLES ##

correct_classification_samples <- vector("list", kfold)

for (kf in 1:kfold) {
  
  correct_classification_matrix <- matrix(NA, nrow = N_samples, ncol = t)
  
  for (i in 1:N_samples) {
    
    test_items <- get(paste0("test_data_", kf)) %>%
      select(Item, Type) %>%
      unique() 
    
    classified_types <- type_samples[, i, kf]
    classified_types_test <- classified_types[as.numeric(test_items$Item)] 
    
    correct_test <- test_items %>%
      mutate(classified = classified_types[match(Item, seq_along(classified_types))]) %>%
      mutate(Type = as.numeric(Type)) %>%
      group_by(Type) %>%
      summarise(correct = sum(classified == Type, na.rm = TRUE), 
                total = n()) %>%
      mutate(classification_rate = (correct / total)*100)
    
    correct_classification_matrix[i, ] <- correct_test$classification_rate
  }
  
  correct_classification_samples[[kf]] <- as.data.frame(correct_classification_matrix) %>%
    rename(classification_rate = V1)
}

classification_rates_all <- bind_rows(correct_classification_samples, .id = "fold") %>%
  pivot_longer(., cols = c(classification_rate, V2, V3, V3, V4, V5), values_to = "value", names_to = "Type") %>%
  mutate(Type = case_when(Type == "V1" | Type == "classification_rate" ~ "bulb",
                          Type == "V2" ~ "car window",
                          Type == "V3" ~ "headlamp",
                          Type == "V4" ~ "container",
                          Type == "V5" ~ "building window")) %>%
  mutate(Type = factor(Type, levels = c("bulb", "car window", "headlamp", "container", "building window"))) %>%
  arrange(Type) %>%
  rename(classification_rate = value)


classification_rates_all %>%
  mutate(Type = factor(Type, levels = c("bulb", "car window", "headlamp", "container", "building window"))) %>%
  ggplot(., aes(x = classification_rate, fill = Type)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Type, scales = "free_y") +
  
  labs(
    x = "",
    y = ""
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "bottom")

pdf("output/classification_density.pdf", height = 4)
classification_rates_all %>%
  mutate(Type = case_when(Type == "car_window" ~ "car window",
                          Type == "building_window" ~ "building window",
                          TRUE ~ Type)) %>%
  mutate(Type = factor(Type, levels = c("bulb", "car window", "headlamp", "container", "building window"))) %>%
  ggplot(., aes(x = classification_rate, fill = Type)) +
  geom_density(alpha = 0.5, adjust = 1.5) +
  facet_wrap(~Type, scale = 'free_y', ncol = 3) +
  #facet_grid(~Type, scales = 'free', space = 'free') +
  #facet_grid2(~Type, scales='free', independent='y') +
  labs(
    x = "",
    y = ""
  ) +
  scale_fill_brewer(palette = "Set2", name = "Classified Glass Type") +
  theme(legend.position = "bottom",
        legend.margin = margin(t = -15))
dev.off()


classification_rates_all %>%
  group_by(Type) %>%
  summarise(classification_rate = mean(classification_rate))

classification_rates_all %>%
  summarise(classification_rate = mean(classification_rate))


##################################.

## CLASSIFICATION SCORES ##

source("code/classification_score_functions.R")

glass_types <- c("bulb", "car window", "headlamp", "container", "building window")


test_classification_table <- predicted_types_test %>%
  group_by(classified, Type) %>%
  reframe(n = n()) %>% 
  ungroup() %>%
  as.data.frame() %>%
  pivot_wider(., names_from = classified, values_from = n, values_fill = 0) %>%
  rename(bulb = `1`,
         `car window` = `2`,
         headlamp = `3`,
         container = `4`,
         `building window` = `5`) %>%
  arrange(Type) %>%
  select(-Type) %>%
  as.matrix()

rownames(test_classification_table) <- glass_types

test_classification_table 


# NUMBER OF INSTANCES #
n_instances = sum(test_classification_table) 

# NUMBER OF CLASSES #
n_class = nrow(test_classification_table) 

# CORRECT CLASSIFIED #
n_correct_classified = diag(test_classification_table)

# NUMBER OF INSTANCES PER CLASS #
n_instances_class = apply(test_classification_table, 2, sum) 

# NUMBER OF PREDICTIONS PER CLASS #
n_predictions = apply(test_classification_table, 1, sum) 

# INSTANCES OVER ACTUAL CLASS #
d_class = n_instances_class / n_instances

# INSTANCES OVER PREDICTED CLASS #
d_predictions = n_predictions / n_instances 


#### ACCURACY ####

accuracy <- sum(n_correct_classified) / n_instances
accuracy

#### % MISSCLASSIFIED ####

missclassified <- 100 - (sum(diag(test_classification_table)) / sum(test_classification_table) * 100)
missclassified


#### CONFUSION MATRIX ####

calculate_metrics <- function(glass_type, classification_table) {
  
  TP <- classification_table[glass_type, glass_type]
  FP <- sum(classification_table[glass_type, ]) - TP
  FN <- sum(classification_table[, glass_type]) - TP
  TN <- sum(classification_table) - (TP + FP + FN)
  
  return(matrix(c(TP, FP, FN, TN), nrow = 2, byrow = TRUE))
}

confusion_matrices <- list()
for (i in 1:t) {
  confusion_matrices[[i]] <- calculate_metrics(glass_types[[i]], test_classification_table)
}
confusion_matrices


#### PRECISION ####

precision = n_correct_classified / n_predictions
precision


#### RECALL ####

recall <- n_correct_classified / n_instances_class
recall


#### F1 SCORE ####

F1 <- (2 * precision * recall) / (precision + recall)
F1



#### MCC ####
# MATTHEWS CORRELATION COEFFICIENT #

MCC <- function(TP, FP, TN, FN) {
  MCC <- ((TP * TN) - (FP * FN)) / 
    sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
}

MCC_bulb <- MCC(TP = confusion_matrices[[1]][1,1], 
                FP = confusion_matrices[[1]][1,2],
                FN = confusion_matrices[[1]][2,1],
                TN = confusion_matrices[[1]][2,2])

MCC_car <- MCC(TP = confusion_matrices[[2]][1,1], 
               FP = confusion_matrices[[2]][1,2],
               FN = confusion_matrices[[2]][2,1],
               TN = confusion_matrices[[2]][2,2])

MCC_headlamp <- MCC(TP = confusion_matrices[[3]][1,1], 
                    FP = confusion_matrices[[3]][1,2],
                    FN = confusion_matrices[[3]][2,1],
                    TN = confusion_matrices[[3]][2,2])

MCC_container <- MCC(TP = confusion_matrices[[4]][1,1], 
                     FP = confusion_matrices[[4]][1,2],
                     FN = confusion_matrices[[4]][2,1],
                     TN = confusion_matrices[[4]][2,2])

MCC_building <- MCC(TP = confusion_matrices[[5]][1,1], 
                    FP = confusion_matrices[[5]][1,2],
                    FN = confusion_matrices[[5]][2,1],
                    TN = confusion_matrices[[5]][2,2])


MCC <- matrix(c(MCC_bulb,
                MCC_car,
                MCC_headlamp,
                MCC_container,
                MCC_building),
              ncol = 5,
              byrow = TRUE)  %>%
  'colnames<-'(c("bulb", "car window", "headlamp", "container", "building window"))
MCC


#### GOODMAN AND KRUSKAL TAU ####

GKtau(test_classification_table)

#### COHEN KAPPA ####

CohenKappa(test_classification_table)

#### THIEL U ####

TheilU(test_classification_table)


##################################.

## POSTERIOR PROBABILTIES - ITEM TYPE ####

test_item_index <- c(unique(test_data_1$Item),
                     unique(test_data_2$Item),
                     unique(test_data_3$Item),
                     unique(test_data_4$Item),
                     unique(test_data_5$Item))

test_item_type_index <- rbind(test_data_1 %>% select(Item, Type) %>% mutate(Type = as.numeric(Type)) %>% unique(),
                              test_data_2 %>% select(Item, Type) %>% mutate(Type = as.numeric(Type)) %>% unique(),
                              test_data_3 %>% select(Item, Type) %>% mutate(Type = as.numeric(Type)) %>% unique(),
                              test_data_4 %>% select(Item, Type) %>% mutate(Type = as.numeric(Type)) %>% unique(),
                              test_data_5 %>% select(Item, Type) %>% mutate(Type = as.numeric(Type)) %>% unique())


test_type_samples <- do.call(rbind, lapply(1:kfold, function(kf) {
  
  type_samples[unique(get(paste0("test_data_", kf))$Item), , kf]
  
}))

test_type_samples <- cbind(Item = test_item_index, test_type_samples)
dim(test_type_samples)


posterior_prob <- test_type_samples %>%
  as.data.frame() %>%
  group_by(Item) %>%
  reframe(
    prob_bulb = sum(c_across(starts_with("V")) == 1) / N_samples,
    prob_car = sum(c_across(starts_with("V")) == 2) / N_samples,
    prob_headlamp = sum(c_across(starts_with("V")) == 3) / N_samples,
    prob_container = sum(c_across(starts_with("V")) == 4) / N_samples,
    prob_building = sum(c_across(starts_with("V")) == 5) / N_samples
  )
posterior_prob

rowSums(posterior_prob[,-1])


posterior_prob_type <- posterior_prob %>%
  left_join(test_item_type_index, by = "Item")


### BRIER SCORE ####

brier_score_data <- posterior_prob_type %>%
  
  mutate(
    actual_bulb = ifelse(Type == 1, 1, 0),
    actual_car = ifelse(Type == 2, 1, 0),
    actual_headlamp = ifelse(Type == 3, 1, 0),
    actual_container = ifelse(Type == 4, 1, 0),
    actual_building = ifelse(Type == 5, 1, 0)
  ) %>%
  
  mutate(
    diff_bulb = (actual_bulb - prob_bulb)^2,
    diff_car = (actual_car - prob_car)^2,
    diff_headlamp = (actual_headlamp - prob_headlamp)^2,
    diff_container = (actual_container - prob_container)^2,
    diff_building = (actual_building - prob_building)^2
  )


brier_score <- brier_score_data %>%
  summarise(
    brier_bulb = mean(diff_bulb),
    brier_car = mean(diff_car),
    brier_headlamp = mean(diff_headlamp),
    brier_container = mean(diff_container),
    brier_building = mean(diff_building)
  )

brier_score


overall_brier_score <- brier_score_data %>%
  mutate(
    overall_diff = (actual_bulb - prob_bulb)^2 + 
      (actual_car - prob_car)^2 +
      (actual_headlamp - prob_headlamp)^2 +
      (actual_container - prob_container)^2 +
      (actual_building - prob_building)^2
  ) %>%
  summarise(overall_brier_score = mean(overall_diff))

overall_brier_score


### ECE ####

ece <- posterior_prob_type %>%
  mutate(predicted_prob = pmax(prob_bulb, prob_car, prob_headlamp, prob_container, prob_building)) %>%
  mutate(prob_bin = cut(predicted_prob, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE, labels = FALSE)) %>%
  group_by(prob_bin) %>%
  summarise(
    mean_prob = mean(predicted_prob),
    observed_frequency = mean(Type == 1),  # Adjust for the event type
    calibration_error = abs(mean_prob - observed_frequency)
  ) %>%
  summarise(mean_ece = mean(calibration_error))


ece


### PLOT ####

posterior_prob_type %>%
  pivot_longer(., cols = -c(Item, Type), names_to = "classified_type", values_to = "prob") %>%
  group_by(Item) %>%
  arrange(desc(prob)) %>% 
  slice(1:2) %>%
  ungroup() %>%
  mutate(Type = case_when(Type == 1 ~ "bulb",
                          Type == 2 ~ "car window",
                          Type == 3 ~ "headlamp",
                          Type == 4 ~ "container",
                          Type == 5 ~ "building window")) %>%
  mutate(classified_type = case_when(classified_type == "prob_bulb" ~ "bulb",
                                     classified_type == "prob_car" ~ "car window",
                                     classified_type == "prob_headlamp" ~ "headlamp",
                                     classified_type == "prob_container" ~ "container",
                                     classified_type == "prob_building" ~ "building window")) %>%
  mutate(Type = factor(Type , levels = c("bulb", "car window", "headlamp", "container", "building window"))) %>%
  mutate(classified_type = factor(classified_type , levels = c("bulb", "car window", "headlamp", "container", "building window"))) %>%
  group_by(Type,) %>% 
  mutate(id = dense_rank(Item)) %>% #
  ungroup() %>% 
  
  ggplot(aes(x = id, y = prob)) +
  geom_point(aes(color = factor(classified_type), shape = factor(classified_type)), size = 2.5) + 
  facet_wrap(~Type, scales = "free", ncol = 1) +
  scale_color_brewer(palette = "Set2") +
  labs(color = "Classified Glass Type", shape = "Classified Glass Type", y = "Posterior Probability") +
  theme(legend.position = "bottom",
        legend.box = "vertical")


pdf("output/classification_prob.pdf", height = 6)
posterior_prob_type %>%
  pivot_longer(., cols = -c(Item, Type), names_to = "classified_type", values_to = "prob") %>%
  group_by(Item) %>%
  arrange(desc(prob)) %>% 
  slice(1:2) %>%
  ungroup() %>%
  mutate(Type = case_when(Type == 1 ~ "bulb",
                          Type == 2 ~ "car window",
                          Type == 3 ~ "headlamp",
                          Type == 4 ~ "container",
                          Type == 5 ~ "building window")) %>%
  mutate(classified_type = case_when(classified_type == "prob_bulb" ~ "bulb",
                                     classified_type == "prob_car" ~ "car window",
                                     classified_type == "prob_headlamp" ~ "headlamp",
                                     classified_type == "prob_container" ~ "container",
                                     classified_type == "prob_building" ~ "building window")) %>%
  mutate(Type = factor(Type , levels = c("bulb", "car window", "headlamp", "container", "building window"))) %>%
  mutate(classified_type = factor(classified_type , levels = c("bulb", "car window", "headlamp", "container", "building window"))) %>%
  group_by(Type,) %>% 
  mutate(id = dense_rank(Item)) %>% #
  ungroup() %>% 
  
  ggplot(aes(x = id, y = prob)) +
  geom_point(aes(color = factor(classified_type), shape = factor(classified_type)), size = 2.5) + 
  facet_wrap(~Type, scales = "free_x", ncol = 2) +
  scale_color_brewer(palette = "Set2") +
  labs(color = "Classified Glass Type", shape = "Classified Glass Type", 
       y = "", 
       x = "") +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.margin = margin(t = -15))
dev.off()



##################################.

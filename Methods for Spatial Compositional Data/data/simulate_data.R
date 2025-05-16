#######################################.
## SIMULATE SPATIAL DATA FOR JOURNAL ##
#######################################.

## PACKAGES ####

library(dplyr)
library(mgcv)
library(compositions)
library(readr)

##################################.

set.seed(123)

setwd("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/Statistical Modelling Journal/Spatial/Data")

##################################.

## CREATE SQUARE SPATITAL GRID ####

n_side <- 32  
# 32 x 32 = 1024 points

grid <- expand.grid(
  x = seq(0, 1, length.out = n_side),
  y = seq(0, 1, length.out = n_side)
)

n_locations <- nrow(grid)
n_components <- 10  

##################################.

## SIMULATE SPATIAL SURFACES FOR EACH COMPONENT ####

fitted_surfaces <- list()
sim_data_matrix <- matrix(NA, nrow = n_locations, ncol = n_components)

for (i in 1:n_components) {
  # Generate smooth random surface (GAM with thin-plate spline)
  z <- rnorm(n_locations)  # add some noise
  gam_model <- gam(z ~ s(x, y, bs = "tp"), data = cbind(grid, z = z), family = gaussian())
  
  # Simulate one realization from the GAM
  sim_vals <- simulate(gam_model, nsim = 1)[, 1]
  
  # Store simulated raw values
  sim_data_matrix[, i] <- sim_vals
}

##################################.

## APPLY INV CLR TO CONVERT TO PROPORTIONS ####

proportions <- t(apply(sim_data_matrix, 1, function(row) {
  inv_clr <- exp(row) / sum(exp(row))
  return(inv_clr)
}))


##################################.

## CONSTRUCT DATA ####

simulated_data <- data.frame(
  grid,
  proportions
)

# Rename proportion columns
names(simulated_data)[-(1:2)] <- paste0("variable_", 1:n_components)

head(simulated_data)
rowSums(simulated_data[,-c(1:2)])

##################################.

## SAVE ####

simulated_data %>%
  write_rds(., "simulated_spatial_compositional_data.rds")

##################################.
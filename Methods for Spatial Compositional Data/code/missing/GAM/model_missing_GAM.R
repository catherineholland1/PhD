####################################.
## SPATIAL GAM MODEL - COMPARISON ##
####################################.

## PACKAGES ####

library(dplyr)
library(tidyr)
library(ggplot2)
library(mgcv)

#######################.

## SET WORKING DIRECTORY ##

setwd("GAM")


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


###################.

## CONSTANTS ####

# NUMBER OF TREE TYPES
N_types = length(tree_types)

# NUMBER OF BASIS
#N_basis = 200
N_basis = 400
#N_basis = 600


#######################.

## GAM MODELS ####

larch_gam <- gam(larch / total ~ s(X_coord, Y_coord, k = N_basis),
                 family = quasibinomial(link="logit"),
                 weights = missing_tree_data$total,
                 data = missing_tree_data,
                 method="REML")

oak_gam <- gam(oak / total ~ s(X_coord, Y_coord, k = N_basis),
               family = quasibinomial(link="logit"),
               weights = missing_tree_data$total,
               data = missing_tree_data,
               method="REML")

sitka_spruce_gam <- gam(sitka_spruce / total ~ s(X_coord, Y_coord, k = N_basis),
                        family = quasibinomial(link="logit"),
                        weights = missing_tree_data$total,
                        data = missing_tree_data,
                        method="REML")

sycamore_gam <- gam(sycamore / total ~ s(X_coord, Y_coord, k = N_basis),
                    family = quasibinomial(link="logit"),
                    weights = missing_tree_data$total,
                    data = missing_tree_data,
                    method="REML")

#######################.

## PREDICT ####

predicted_counts <- matrix(NA, nrow = N_fit, ncol = N_types)

for (i in 1:N_types) {
  
  
  predicted_proportions <- predict(get(paste0(tree_types[i], "_gam")),
                                   newdata = missing_tree_data,
                                   type = "response")
  
  predicted_counts[, i] <- round(predicted_proportions * missing_tree_data$total)

}

predicted_counts


### MANUALLY ####

larch_predict <- predict(larch_gam,
                         newdata = missing_tree_data,
                         type = "response")

oak_predict <- predict(oak_gam,
                       newdata = missing_tree_data,
                       type = "response")

sitka_spruce_predict <- predict(sitka_spruce_gam,
                                newdata = missing_tree_data,
                                type = "response")

sycamore_predict <- predict(sycamore_gam,
                            newdata = missing_tree_data,
                            type = "response")


predicted_counts <- cbind(larch = round(larch_predict * missing_tree_data$total),
                          oak = round(oak_predict * missing_tree_data$total),
                          sitka_spruce = round(sitka_spruce_predict * missing_tree_data$total),
                          sycamore = round(sycamore_predict * missing_tree_data$total))


#######################.

## SAVE PREDICTED COUNTS ####

saveRDS(predicted_counts, paste0("predictions_GAM-", N_basis, ".rds"))


#######################.

## EXAMINE PREDICTIONS ####

predicted_counts_data <- predicted_counts %>%
  as.data.frame() %>%
  rename(larch = "V1",
         oak = "V2",
         sitka_spruce = "V3",
         sycamore = "V4") %>%
  cbind(missing_tree_data %>% select(contains("coord"),
                                     contains("missing"))) %>%
  pivot_longer(., cols = larch:sycamore, 
               names_to = "type", values_to = "predicted") %>%
  left_join((selected_data %>% pivot_longer(., cols = larch:sycamore, 
                                            names_to = "type", values_to = "actual") %>%
               select(-total)), by = c("X_coord", "Y_coord", "type"))


### PLOT ####

#### PRED VS ACTUAL ####

predicted_counts_data %>%
  
  filter((type == "larch" & larch_missing == TRUE) |
           (type == "oak" & oak_missing == TRUE) |
           (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
           (type == "sycamore" & sycamore_missing == TRUE)) %>%
  
  ggplot(., aes(x = actual, y = predicted)) +
  geom_point() +
  geom_abline(col = "red") +
  facet_wrap(~type*n_missing, ncol = max(n_missing)) +
  ggtitle("GAM Model - Predicted vs Actual")


#### RESID VS PRED ####

predicted_counts_data %>%
  
  filter((type == "larch" & larch_missing == TRUE) |
           (type == "oak" & oak_missing == TRUE) |
           (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
           (type == "sycamore" & sycamore_missing == TRUE)) %>%
  
  mutate(residuals = actual - predicted) %>%
  
  ggplot(., aes(x = predicted, y = residuals)) +
  geom_point() +
  facet_wrap(~type*n_missing, ncol = max(n_missing)) +
  ggtitle("GAM Model - Residuals vs Predicted")


#### RESID VS ACTUAL ####

predicted_counts_data %>%
  
  filter((type == "larch" & larch_missing == TRUE) |
           (type == "oak" & oak_missing == TRUE) |
           (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
           (type == "sycamore" & sycamore_missing == TRUE)) %>%
  
  mutate(residuals = actual - predicted) %>%
  
  ggplot(., aes(x = actual, y = residuals)) +
  geom_point() +
  facet_wrap(~type*n_missing, ncol = max(n_missing)) +
  ggtitle("GAM Model - Residuals vs Actual")


#### HEATMAP ####

predicted_counts_data %>%
  
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
  ggtitle("Heatmap for Tree Types - Actual vs Predicted - GAM Model")


#### DENSITY ####

predicted_counts_data %>%
  
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
  
  labs(title = "Density Plot of Predicted and Actual Counts - GAM Model",
       x = "Log of Predicted Counts", y = "Density") +
  theme(legend.position = "bottom")


#### QUANTIFY ####

## MAE ##

predicted_counts_data %>%
  filter((type == "larch" & larch_missing == TRUE) |
           (type == "oak" & oak_missing == TRUE) |
           (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
           (type == "sycamore" & sycamore_missing == TRUE)) %>%
  group_by(type, n_missing) %>%
  summarise(mae = mean(abs(predicted - actual), na.rm = TRUE))


## RMSE ##

predicted_counts_data %>%
  filter((type == "larch" & larch_missing == TRUE) |
           (type == "oak" & oak_missing == TRUE) |
           (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
           (type == "sycamore" & sycamore_missing == TRUE)) %>%
  group_by(type, n_missing) %>%
  summarise(rmse = sqrt(mean(predicted - actual)^2))


#######################.

## SPLINE VARIANCE ####

spline_variances <- vector(length = N_types)

for (i in 1:N_types) {
  
  current_gam <- get(paste0(tree_types[i], "_gam"))
  
  gam_coefficients <- coef(current_gam)
  var_cov_matrix <- vcov(current_gam)
  
  spline_indices <- which(grepl("s\\(X_coord|s\\(Y_coord", names(gam_coefficients)))
  
  spline_se <- sqrt(diag(var_cov_matrix[spline_indices, spline_indices]))
  
  spline_variances_current <- spline_se^2
  
  spline_variances[i] <- mean(spline_variances_current)
}


spline_variances


#######################.

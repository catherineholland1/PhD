#####################################.
## SPATIAL LARCH OUTPUT COMAPRISON ##
#####################################.

## PACKAGES ####

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggh4x)

#######################.

## SET WORKING DIRECTORY ##

setwd("comparison")


#######################.

## DATA ####

tree_data <- readRDS("") 


### LARCH ####

larch_data <- tree_data %>% 
  select(X_coord,
         Y_coord,
         larch,
         total) 



N = nrow(larch_data)


## TRAIN / TEST DATA ####

N_train <- ceiling(N*0.5) # CHANGED
N_test <- floor(N*0.5)

set.seed(453459) # CHANGED
train_index <- sort(sample(1:N,N_train,replace=FALSE))
test_index <- setdiff(1:N,train_index)


#######################.

## TRAIN OUTPUT ####

## LOAD LARCH ##

load('larch/output/replicates_polynomial.RData')
larch_train <- simulation_replicates

rm(simulation_replicates)

## LOAD LARCH COMPARISON ##

load('larch/output/replicates_fixed.RData')
larch_comparison_train <- simulation_replicates

rm(simulation_replicates)



### COMBINE REPLICATES ####

train_output_mean <- rbind(mutate(as.data.frame(larch_train %>% apply(., 2, mean)), method = "Polynomial phi") %>%
                             rename(mean = "larch_train %>% apply(., 2, mean)"),
                           mutate(as.data.frame(larch_comparison_train %>% apply(., 2, mean)), method = "Fixed phi") %>%
                             rename(mean = "larch_comparison_train %>% apply(., 2, mean)"))

train_output_sd <- rbind(mutate(as.data.frame(larch_train %>% apply(., 2, sd)), method = "Polynomial phi") %>%
                             rename(sd = "larch_train %>% apply(., 2, sd)"),
                           mutate(as.data.frame(larch_comparison_train %>% apply(., 2, sd)), method = "Fixed phi") %>%
                             rename(sd = "larch_comparison_train %>% apply(., 2, sd)"))


### PLOT ####

#### DENSITY - MEAN ####

train_output_mean %>%
  
  mutate(method = factor(method, levels = c("Polynomial phi", "Fixed phi"))) %>%
  
  ggplot(., aes(x = mean, fill = method)) +
  
  geom_density(alpha = 0.5) +
  
  geom_vline(data = larch_data[train_index, ] %>%
               summarise(mean = mean(larch)), 
             aes(xintercept = mean, color = "Original Data"), 
             linetype = "dashed", linewidth = 1) +
  
  scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
  scale_fill_manual(name = "Method", values = c("Polynomial phi" = "#66C2A5",
                                          "Fixed phi" = "#E78AC3")) +
  
  labs() +
  theme(legend.position = "bottom")


pdf("output/train_density_mean.pdf")
train_output_mean %>%
  
  mutate(method = factor(method, levels = c("Polynomial phi", "Fixed phi"))) %>%
  
  ggplot(., aes(x = mean, fill = method)) +
  
  geom_density(alpha = 0.5) +
  
  geom_vline(data = larch_data[train_index, ] %>%
               summarise(mean = mean(larch)), 
             aes(xintercept = mean, color = "Original Data"), 
             linetype = "dashed", linewidth = 1) +
  
  scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
  scale_fill_manual(name = "Method", values = c("Polynomial phi" = "#66C2A5",
                                                "Fixed phi" = "#E78AC3")) +
  
  labs(title = "Posterior Density Plot of the Mean for Larch - Train") +
  theme(legend.position = "bottom")
dev.off()


#### DENSITY - SD ####

train_output_sd %>%
  mutate(method = factor(method, levels = c("Polynomial phi", "Fixed phi"))) %>%
  
  ggplot(., aes(x = sd, fill = method)) +
  
  geom_density(alpha = 0.5) +
  
  geom_vline(data = larch_data[train_index, ] %>%
               summarise(sd = sd(larch)), 
             aes(xintercept = sd, color = "Original Data"), 
             linetype = "dashed", linewidth = 1) +
  
  scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
  scale_fill_manual(name = "Method", values = c("Polynomial phi" = "#66C2A5",
                                                "Fixed phi" = "#E78AC3")) +
  
  labs(title = "Posterior Density Plot of the Standard Deviation for Larch - Train") +
  theme(legend.position = "bottom")


pdf("output/train_density_sd.pdf")
train_output_sd %>%
  mutate(method = factor(method, levels = c("Polynomial phi", "Fixed phi"))) %>%
  
  ggplot(., aes(x = sd, fill = method)) +
  
  geom_density(alpha = 0.5) +
  
  geom_vline(data = larch_data[train_index, ] %>%
               summarise(sd = sd(larch)), 
             aes(xintercept = sd, color = "Original Data"), 
             linetype = "dashed", linewidth = 1) +
  
  scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
  scale_fill_manual(name = "Method", values = c("Polynomial phi" = "#66C2A5",
                                                "Fixed phi" = "#E78AC3")) +
  
  labs(title = "Posterior Density Plot of the Standard Deviation for Larch - Train") +
  theme(legend.position = "bottom")
dev.off()


#### QUANTILE ####

train_data_quantiles <- as.data.frame(quantile(larch_data[train_index,]$larch, probs = seq(0.01, 0.99, by = 0.01))) %>%
  tibble::rownames_to_column(var = "quantile") %>%
  rename_with(~"value", 2) %>%
  mutate(quantile = gsub("%", "", quantile),
         quantile = as.numeric(quantile) / 100)


polynomial_train_quantiles <- larch_train %>%
  apply(., 2, function(x) quantile(x, probs = seq(0.01, 0.99, by = 0.01))) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(var = "quantile") %>%
  pivot_longer(., cols = starts_with("V"), names_to = "sample", values_to = "x") %>%
  mutate(sample = gsub("^V", "", sample)) %>%
  mutate(quantile = gsub("%", "", quantile),
         quantile = as.numeric(quantile) / 100)

fixed_train_quantiles <- larch_comparison_train %>%
  apply(., 2, function(x) quantile(x, probs = seq(0.01, 0.99, by = 0.01))) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(var = "quantile") %>%
  pivot_longer(., cols = starts_with("V"), names_to = "sample", values_to = "x") %>%
  mutate(sample = gsub("^V", "", sample)) %>%
  mutate(quantile = gsub("%", "", quantile),
         quantile = as.numeric(quantile) / 100)


train_quantile_summary <- rbind(polynomial_train_quantiles %>%
                                  group_by(quantile) %>%
                                  summarise(lower = quantile(x, 0.025),
                                            mean = mean(x),
                                            median = median(x),
                                            upper = quantile(x, 0.975)) %>%
                                  mutate(method = "Polynomial phi"),
                                fixed_train_quantiles %>%
                                  group_by(quantile) %>%
                                  summarise(lower = quantile(x, 0.025),
                                            mean = mean(x),
                                            median = median(x),
                                            upper = quantile(x, 0.975)) %>%
                                  mutate(method = "Fixed phi"))


train_quantile_summary %>%
  
  mutate(method = factor(method, levels = c("Polynomial phi", "Fixed phi"))) %>%
  
  ggplot(., aes(x = quantile, y = mean)) +
  
  geom_point(aes(col = method), size = 2, shape = 16) +
  
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95% Interval"), alpha = 0.5) +
  
  geom_line(data = train_data_quantiles, 
            aes(x = quantile, 
                y = value, col = "Original Data")) +
  
  facet_wrap(~method) +
  
  scale_color_manual(name = "", values = c("Original Data" = "#FC8D62",
                                           "Polynomial phi" = "#66C2A5",
                                           "Fixed phi" = "#E78AC3")) +
  scale_fill_manual(name = "", values = c("95% Interval" = "#8DA0CB")) +
  theme(legend.position = "bottom") +
  labs(y = "value") +
  ggtitle("Quantile Plot for Larch - Train") 


pdf("output/train_quantile.pdf", height = 3)
train_quantile_summary %>%
  
  mutate(method = factor(method, levels = c("Polynomial phi", "Fixed phi"))) %>%
  
  ggplot(., aes(x = quantile, y = mean)) +
  
  geom_point(aes(col = method), size = 2, shape = 16) +
  
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95% Interval"), alpha = 0.5) +
  
  geom_line(data = train_data_quantiles, 
            aes(x = quantile, 
                y = value, col = "Original Data")) +
  
  facet_wrap(~method) +
  
  scale_color_manual(name = "", values = c("Original Data" = "#FC8D62",
                                           "Polynomial phi" = "#66C2A5",
                                           "Fixed phi" = "#E78AC3")) +
  scale_fill_manual(name = "", values = c("95% Interval" = "#8DA0CB")) +
  theme(legend.position = "bottom",
        legend.margin = margin(t = -15)) +
  labs(y = "", x = "") 
dev.off()



### STATISTICAL SUMMARY ####

larch_train_data <- larch_data[train_index,]$larch


#### MAE ####

## MEAN ##

train_mae_polynomial <- mean(abs(apply(larch_train, 1, mean) - larch_train_data))
train_mae_fixed <- mean(abs(apply(larch_comparison_train, 1, mean) - larch_train_data))

c(train_mae_polynomial, train_mae_fixed)


## QUANTILE ##

train_mae_quantile_polynomial <- mean(abs(filter(train_quantile_summary,method=="Polynomial phi")$median - train_data_quantiles$value))
train_mae_quantile_fixed <- mean(abs(filter(train_quantile_summary,method=="Fixed phi")$median - train_data_quantiles$value))

c(train_mae_quantile_polynomial, train_mae_quantile_fixed)


#### RMSE ####

train_rmse_polynomial <- sqrt(mean((apply(larch_train, 1, mean) - larch_train_data)^2))
train_rmse_fixed <- sqrt(mean((apply(larch_comparison_train, 1, mean) - larch_train_data)^2))

c(train_rmse_polynomial, train_rmse_fixed)


#### PREDICITION COVERAGE INTERVALS ####

train_output_intervals_polynomial <- larch_train %>% 
  apply(., 1, function(x) quantile(x, probs = c(0.005, 0.025, 0.05, 0.075, 0.1, 0.5, 0.9, 0.925, 0.95, 0.975, 0.995))) %>% 
  t(.) %>%
  as.data.frame()

train_output_intervals_fixed <- larch_comparison_train %>% 
  apply(., 1, function(x) quantile(x, probs = c(0.005, 0.025, 0.05, 0.075, 0.1, 0.5, 0.9, 0.925, 0.95, 0.975, 0.995))) %>% 
  t(.) %>%
  as.data.frame()


train_output_intervals_check_polynomial <- train_output_intervals_polynomial %>%
  mutate(larch = larch_data[train_index, ]$larch) %>%
  mutate(interval_80 = larch >= `10%` & larch <= `90%`,
         interval_85 = larch >= `7.5%` & larch <= `92.5%`,
         interval_90 = larch >= `5%` & larch <= `95%`,
         interval_95 = larch >= `2.5%` & larch <= `97.5%`,
         interval_99 = larch >= `0.5%` & larch <= `99.5%`)

train_output_intervals_check_polynomial %>%
  summarise(across(starts_with("interval"), ~ mean(. == TRUE) * 100))


train_output_intervals_check_polynomial %>% 
  summarise(range_80 = mean(`90%` - `10%`),
            range_85 = mean(`92.5%` - `7.5%`),
            range_90 = mean(`95%` - `5%`),
            range_95 = mean(`97.5%` - `2.5%`),
            range_99 = mean(`99.5%` - `0.5%`))


train_output_intervals_check_fixed <- train_output_intervals_fixed %>%
  mutate(larch = larch_data[train_index, ]$larch) %>%
  mutate(interval_80 = larch >= `10%` & larch <= `90%`,
         interval_85 = larch >= `7.5%` & larch <= `92.5%`,
         interval_90 = larch >= `5%` & larch <= `95%`,
         interval_95 = larch >= `2.5%` & larch <= `97.5%`,
         interval_99 = larch >= `0.5%` & larch <= `99.5%`)

train_output_intervals_check_fixed %>%
  summarise(across(starts_with("interval"), ~ mean(. == TRUE) * 100))


train_output_intervals_check_fixed %>% 
  summarise(range_80 = mean(`90%` - `10%`),
            range_85 = mean(`92.5%` - `7.5%`),
            range_90 = mean(`95%` - `5%`),
            range_95 = mean(`97.5%` - `2.5%`),
            range_99 = mean(`99.5%` - `0.5%`))


#######################.

## TEST OUTPUT ####

## LOAD LARCH ##

load('~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD R Code/Spatial - Trees/Larch Model/output/predicted_counts_polynomial.RData')
larch_test <- predicted_counts

rm(predicted_counts)


## LOAD LARCH COMPARISON ##

load('~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD R Code/Spatial - Trees/Larch Model/output/predicted_counts_fixed.RData')
larch_comparison_test <- predicted_counts

rm(predicted_counts)


### COMBINE REPLICATES ####

test_output_mean <- rbind(mutate(as.data.frame(larch_test %>% apply(., 2, mean)), method = "Polynomial phi") %>%
                             rename(mean = "larch_test %>% apply(., 2, mean)"),
                           mutate(as.data.frame(larch_comparison_test %>% apply(., 2, mean)), method = "Fixed phi") %>%
                             rename(mean = "larch_comparison_test %>% apply(., 2, mean)"))

test_output_sd <- rbind(mutate(as.data.frame(larch_test %>% apply(., 2, sd)), method = "Polynomial phi") %>%
                           rename(sd = "larch_test %>% apply(., 2, sd)"),
                         mutate(as.data.frame(larch_comparison_test %>% apply(., 2, sd)), method = "Fixed phi") %>%
                           rename(sd = "larch_comparison_test %>% apply(., 2, sd)"))


### PLOT ####

#### DENSITY - MEAN ####

test_output_mean %>%
  
  mutate(method = factor(method, levels = c("Polynomial phi", "Fixed phi"))) %>%
  
  ggplot(., aes(x = mean, fill = method)) +
  
  geom_density(alpha = 0.5) +
  
  geom_vline(data = larch_data[test_index, ] %>%
               summarise(mean = mean(larch)), 
             aes(xintercept = mean, color = "Original Data"), 
             linetype = "dashed", linewidth = 1) +
  
  scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
  scale_fill_manual(name = "Method", values = c("Polynomial phi" = "#66C2A5",
                                                "Fixed phi" = "#E78AC3")) +
  
  labs(title = "Posterior Density Plot of the Mean for Larch - Test") +
  theme(legend.position = "bottom")


pdf("output/test_density_mean.pdf")
test_output_mean %>%
  
  mutate(method = factor(method, levels = c("Polynomial phi", "Fixed phi"))) %>%
  
  ggplot(., aes(x = mean, fill = method)) +
  
  geom_density(alpha = 0.5) +
  
  geom_vline(data = larch_data[test_index, ] %>%
               summarise(mean = mean(larch)), 
             aes(xintercept = mean, color = "Original Data"), 
             linetype = "dashed", linewidth = 1) +
  
  scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
  scale_fill_manual(name = "Method", values = c("Polynomial phi" = "#66C2A5",
                                                "Fixed phi" = "#E78AC3")) +
  
  labs(title = "Posterior Density Plot of the Mean for Larch - Test") +
  theme(legend.position = "bottom")
dev.off()


#### DENSITY - SD ####

test_output_sd %>%
  
  mutate(method = factor(method, levels = c("Polynomial phi", "Fixed phi"))) %>%
  
  ggplot(., aes(x = sd, fill = method)) +
  
  geom_density(alpha = 0.5) +
  
  geom_vline(data = larch_data[test_index, ] %>%
               summarise(sd = sd(larch)), 
             aes(xintercept = sd, color = "Original Data"), 
             linetype = "dashed", linewidth = 1) +
  
  scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
  scale_fill_manual(name = "Method", values = c("Polynomial phi" = "#66C2A5",
                                                "Fixed phi" = "#E78AC3")) +
  
  labs(title = "Posterior Density Plot of the Standard Deviation for Larch - Test") +
  theme(legend.position = "bottom")


pdf("output/test_density_sd.pdf")
test_output_sd %>%
  
  mutate(method = factor(method, levels = c("Polynomial phi", "Fixed phi"))) %>%
  
  ggplot(., aes(x = sd, fill = method)) +
  
  geom_density(alpha = 0.5) +
  
  geom_vline(data = larch_data[test_index, ] %>%
               summarise(sd = sd(larch)), 
             aes(xintercept = sd, color = "Original Data"), 
             linetype = "dashed", linewidth = 1) +
  
  scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
  scale_fill_manual(name = "Method", values = c("Polynomial phi" = "#66C2A5",
                                                "Fixed phi" = "#E78AC3")) +
  
  labs(title = "Posterior Density Plot of the Standard Deviation for Larch - Test") +
  theme(legend.position = "bottom")
dev.off()


#### QUANTILES ####

test_data_quantiles <- as.data.frame(quantile(larch_data[test_index,]$larch, probs = seq(0.01, 0.99, by = 0.01))) %>%
  tibble::rownames_to_column(var = "quantile") %>%
  rename_with(~"value", 2) %>%
  mutate(quantile = gsub("%", "", quantile),
         quantile = as.numeric(quantile) / 100)

polynomial_test_quantiles <- larch_test %>%
  apply(., 2, function(x) quantile(x, probs = seq(0.01, 0.99, by = 0.01))) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(var = "quantile") %>%
  pivot_longer(., cols = starts_with("V"), names_to = "sample", values_to = "x") %>%
  mutate(sample = gsub("^V", "", sample)) %>%
  mutate(quantile = gsub("%", "", quantile),
         quantile = as.numeric(quantile) / 100)

fixed_test_quantiles <- larch_comparison_test %>%
  apply(., 2, function(x) quantile(x, probs = seq(0.01, 0.99, by = 0.01))) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(var = "quantile") %>%
  pivot_longer(., cols = starts_with("V"), names_to = "sample", values_to = "x") %>%
  mutate(sample = gsub("^V", "", sample)) %>%
  mutate(quantile = gsub("%", "", quantile),
         quantile = as.numeric(quantile) / 100)


test_quantile_summary <- rbind(polynomial_test_quantiles %>%
                                  group_by(quantile) %>%
                                  summarise(lower = quantile(x, 0.025),
                                            mean = mean(x),
                                            median = median(x),
                                            upper = quantile(x, 0.975)) %>%
                                  mutate(method = "Polynomial phi"),
                                fixed_test_quantiles %>%
                                  group_by(quantile) %>%
                                  summarise(lower = quantile(x, 0.025),
                                            mean = mean(x),
                                            median = median(x),
                                            upper = quantile(x, 0.975)) %>%
                                  mutate(method = "Fixed phi"))


test_quantile_summary %>%
  
  mutate(method = factor(method, levels = c("Polynomial phi", "Fixed phi"))) %>%
  
  ggplot(., aes(x = quantile, y = mean)) +
  
  geom_point(aes(col = method), size = 2, shape = 16) +
  
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95% Interval"), alpha = 0.5) +
  
  geom_line(data = test_data_quantiles, 
            aes(x = quantile, 
                y = value, col = "Original Data")) +
  
  facet_wrap(~method) +
  
  scale_color_manual(name = "", values = c("Original Data" = "#FC8D62",
                                           "Polynomial phi" = "#66C2A5",
                                           "Fixed phi" = "#E78AC3")) +
  scale_fill_manual(name = "", values = c("95% Interval" = "#8DA0CB")) +
  theme(legend.position = "bottom") +
  labs(y = "value") +
  ggtitle("Quantile Plot for Larch - Test") 


pdf("output/test_quantile.pdf", height = 3)
test_quantile_summary %>%
  
  mutate(method = factor(method, levels = c("Polynomial phi", "Fixed phi"))) %>%
  
  ggplot(., aes(x = quantile, y = median)) +
  
  geom_point(aes(col = method), size = 2, shape = 16) +
  
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95% Interval"), alpha = 0.5) +
  
  geom_line(data = test_data_quantiles, 
            aes(x = quantile, 
                y = value, col = "Original Data")) +
  
  facet_wrap(~method) +
  
  scale_color_manual(name = "", values = c("Original Data" = "#FC8D62",
                                           "Polynomial phi" = "#66C2A5",
                                           "Fixed phi" = "#E78AC3")) +
  scale_fill_manual(name = "", values = c("95% Interval" = "#8DA0CB")) +
  theme(legend.position = "bottom",
        legend.margin = margin(t=-15)) +
  labs(y = "", x = "")
dev.off()


### STATISTICAL SUMMARY ####


larch_test_data <- larch_data[test_index,]$larch


#### MAE ####

## MEAN ##

test_mae_polynomial <- mean(abs(apply(larch_test, 1, mean) - larch_test_data))
test_mae_fixed <- mean(abs(apply(larch_comparison_test, 1, mean) - larch_test_data))

c(test_mae_polynomial, test_mae_fixed)


## QUANTILE ##

test_mae_quantile_polynomial <- mean(abs(filter(test_quantile_summary,method=="Polynomial phi")$median-test_data_quantiles$value))
test_mae_quantile_fixed <- mean(abs(filter(test_quantile_summary,method=="Fixed phi")$median-test_data_quantiles$value))

c(test_mae_quantile_polynomial, test_mae_quantile_fixed)


#### RMSE ####

test_rmse_polynomial <- sqrt(mean((apply(larch_test, 1, mean) - larch_test_data)^2))
test_rmse_fixed <- sqrt(mean((apply(larch_comparison_test, 1, mean) - larch_test_data)^2))

c(test_rmse_polynomial, test_rmse_fixed)


#### PREDICITION COVERAGE INTERVALS ####

test_output_intervals_polynomial <- larch_test %>% 
  apply(., 1, function(x) quantile(x, probs = c(0.005, 0.025, 0.05, 0.075, 0.1, 0.5, 0.9, 0.925, 0.95, 0.975, 0.995))) %>% 
  t(.) %>%
  as.data.frame()

test_output_intervals_fixed <- larch_comparison_test %>% 
  apply(., 1, function(x) quantile(x, probs = c(0.005, 0.025, 0.05, 0.075, 0.1, 0.5, 0.9, 0.925, 0.95, 0.975, 0.995))) %>% 
  t(.) %>%
  as.data.frame()


test_output_intervals_check_polynomial <- test_output_intervals_polynomial %>%
  mutate(larch = larch_data[test_index, ]$larch) %>%
  mutate(interval_80 = larch >= `10%` & larch <= `90%`,
         interval_85 = larch >= `7.5%` & larch <= `92.5%`,
         interval_90 = larch >= `5%` & larch <= `95%`,
         interval_95 = larch >= `2.5%` & larch <= `97.5%`,
         interval_99 = larch >= `0.5%` & larch <= `99.5%`)

test_output_intervals_check_polynomial %>%
  summarise(across(starts_with("interval"), ~ mean(. == TRUE) * 100))


test_output_intervals_check_polynomial %>% 
  summarise(range_80 = mean(`90%` - `10%`),
            range_85 = mean(`92.5%` - `7.5%`),
            range_90 = mean(`95%` - `5%`),
            range_95 = mean(`97.5%` - `2.5%`),
            range_99 = mean(`99.5%` - `0.5%`))


test_output_intervals_check_fixed <- test_output_intervals_fixed %>%
  mutate(larch = larch_data[test_index, ]$larch) %>%
  mutate(interval_80 = larch >= `10%` & larch <= `90%`,
         interval_85 = larch >= `7.5%` & larch <= `92.5%`,
         interval_90 = larch >= `5%` & larch <= `95%`,
         interval_95 = larch >= `2.5%` & larch <= `97.5%`,
         interval_99 = larch >= `0.5%` & larch <= `99.5%`)

test_output_intervals_check_fixed %>%
  summarise(across(starts_with("interval"), ~ mean(. == TRUE) * 100))


test_output_intervals_check_fixed %>% 
  summarise(range_80 = mean(`90%` - `10%`),
            range_85 = mean(`92.5%` - `7.5%`),
            range_90 = mean(`95%` - `5%`),
            range_95 = mean(`97.5%` - `2.5%`),
            range_99 = mean(`99.5%` - `0.5%`))
  

#######################.

## COMBINE TRAIN / TEST PLOT ####

### MEAN ####

pdf("output/density_mean.pdf", height = 3)
train_output_mean %>%
  mutate(cv = "train") %>%
  rbind(test_output_mean %>%
          mutate(cv = "test")) %>% 
  
  mutate(cv = factor(stringr::str_to_title(cv), levels = c("Train", "Test"))) %>%
  mutate(method = factor(method, levels = c("Polynomial phi", "Fixed phi"))) %>%
  
  ggplot(., aes(x = mean, fill = method)) +
  
  geom_density(alpha = 0.5) +
  
  facet_wrap(~cv) +
  
  geom_vline(data = bind_rows(
    larch_data %>% 
      filter(row_number() %in% train_index) %>% 
      summarise(mean = mean(larch)) %>% 
      mutate(cv = "train"),
    larch_data %>% 
      filter(row_number() %in% test_index) %>% 
      summarise(mean = mean(larch)) %>% 
      mutate(cv = "test")
  ) %>% mutate(cv = factor(stringr::str_to_title(cv), levels = c("Train", "Test"))),
  aes(xintercept = mean, color = "Original Data"), 
  linetype = "dashed", linewidth = 1) +
  
  scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
  scale_fill_manual(name = "Method", values = c("Polynomial phi" = "#66C2A5",
                                                "Fixed phi" = "#E78AC3")) +
  
  labs(x = "", y = "") +
  theme(legend.position = "bottom",
        legend.margin = margin(t = -15))
dev.off()


### STANDARD DEVIATION ####

pdf("output/density_sd.pdf", height = 3)
train_output_sd %>%
  mutate(cv = "train") %>%
  rbind(test_output_sd %>%
          mutate(cv = "test")) %>% 
  
  mutate(cv = factor(stringr::str_to_title(cv), levels = c("Train", "Test"))) %>%
  mutate(method = factor(method, levels = c("Polynomial phi", "Fixed phi"))) %>%
  
  ggplot(., aes(x = sd, fill = method)) +
  
  geom_density(alpha = 0.5) +
  
  facet_wrap(~cv) +
  
  geom_vline(data = bind_rows(
    larch_data %>% 
      filter(row_number() %in% train_index) %>% 
      summarise(sd = sd(larch)) %>% 
      mutate(cv = "train"),
    larch_data %>% 
      filter(row_number() %in% test_index) %>% 
      summarise(sd = sd(larch)) %>% 
      mutate(cv = "test")
  ) %>% mutate(cv = factor(stringr::str_to_title(cv), levels = c("Train", "Test"))),
  aes(xintercept = sd, color = "Original Data"), 
  linetype = "dashed", linewidth = 1) +
  
  scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
  scale_fill_manual(name = "Method", values = c("Polynomial phi" = "#66C2A5",
                                                "Fixed phi" = "#E78AC3")) +
  
  labs(x = "", y = "") +
  theme(legend.position = "bottom",
        legend.margin = margin(t = -15))
dev.off()


### QUANTILE ####

pdf('output/quantiles.pdf', height = 4)
train_quantile_summary %>%
  mutate(cv = 'Train') %>%
  rbind(test_quantile_summary %>% mutate(cv = 'Test')) %>%
  
  mutate(method = factor(method, levels = c("Polynomial phi", "Fixed phi"))) %>%
  mutate(cv = factor(cv, levels = c("Train", "Test"))) %>%
  
  ggplot(., aes(x = quantile, y = median)) +
  
  geom_point(aes(col = method), size = 2, shape = 16) +
  
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95% Interval"), alpha = 0.5) +
  
  geom_line(data = train_data_quantiles %>% mutate(cv = 'Train') %>% rbind(test_data_quantiles %>% mutate(cv = 'Test')) %>%
              mutate(cv = factor(cv, levels = c("Train", "Test"))), 
            aes(x = quantile, 
                y = value, col = "Original Data")) +
  
  facet_grid2(cv~method) +
  
  scale_color_manual(name = "", values = c("Original Data" = "#FC8D62",
                                           "Polynomial phi" = "#66C2A5",
                                           "Fixed phi" = "#E78AC3")) +
  scale_fill_manual(name = "", values = c("95% Interval" = "#8DA0CB")) +
  theme(legend.position = "bottom",
        legend.margin = margin(t=-15)) +
  labs(y = "", x = "")
dev.off()


#######################.

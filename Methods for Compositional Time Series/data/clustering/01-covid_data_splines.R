##########################.
## COVID DATA - SPLINES ##.
##########################.


## PACKAGES ##

library(dplyr)
library(tibble)
library(mgcv)
library(forcats)
library(tidyr)


#######################.

## GET DATA ####

data <- readRDS(file = "") %>%
  select_if(~ !is.numeric(.) || sum(.) != 0) %>%
  mutate(variant = fct_relevel(variant, c("alpha", "beta", "delta", "gamma", "omicron", "VOI", "all")))


## POPULATION DATA ##
population_data <- readRDS(file = "")


#######################.

## DATA MANIPULATION ####

# SQRT
sqrt_data <- data %>%
  mutate(across(starts_with("20"), ~ if_else(variant != "all", sqrt(.x), .x))) %>%
  filter(measure == "proportion" | measure == "total_covid_pop") %>%
  arrange(country, variant, measure)


#######################.

## SPLINES ####

spline_data = sqrt_data


## ALPHA ####

alpha_data <- spline_data %>%
  filter(variant == "alpha") %>%
  select(-c(variant, measure)) %>%
  pivot_longer(., cols = starts_with("20"), names_to = "date", values_to = "alpha") %>%
  mutate(date = as.Date(date))


# FIT SPLINE FOR EACH COUNTRY IN VARIANT
alpha_model <- vector("list", length = length(unique(alpha_data$country)))
for (i in 1:length(unique(alpha_data$country))) {
  
  alpha_model[[i]] <- alpha_data %>%
    filter(country == unique(alpha_data$country)[i]) %>%
    gam(alpha ~ s(as.numeric(date),
                k = 10,
                bs = "tp"
                ), 
      data = .)
}



# COEFFICIENTS

alpha_coeff <- vector("list", length(unique(alpha_data$country)))
for (i in 1:length(unique(alpha_data$country))) {
  alpha_coeff[[i]] <- unlist(alpha_model[[i]]["coefficients"])
}
names(alpha_coeff) <- unique(alpha_data$country)


alpha_output_data <- as.data.frame(alpha_coeff) %>%
  t(.) %>%
  `colnames<-`(c("alpha_intercept", "alpha_x1", "alpha_x2", "alpha_x3", "alpha_x4", "alpha_x5", "alpha_x6", "alpha_x7", "alpha_x8", "alpha_x9")) %>%
  as.data.frame() %>%
  rownames_to_column("country")


# RESIDUALS

alpha_dev_residuals <- vector("list", length(unique(alpha_data$country)))
for (i in 1:length(unique(alpha_data$country))) {
  alpha_dev_residuals[[i]] <- unlist(residuals(alpha_model[[i]], type = "deviance"))
}
alpha_dev_residuals <- unlist(alpha_dev_residuals)


# FITTED VALUES

alpha_fitted_values <- vector("list", length(unique(alpha_data$country)))
for (i in 1:length(unique(alpha_data$country))) {
  alpha_fitted_values[[i]] <- unlist(fitted.values(alpha_model[[i]]))
}
alpha_fitted_values <- unlist(alpha_fitted_values)


# TEST RESIDUALS

qqnorm(alpha_dev_residuals,
       main = "Alpha")
qqline(alpha_dev_residuals, col = "red")

hist(alpha_dev_residuals,
     main = "Alpha Deviance Residuals")
abline(v = 0, col = "red", lty = 2)

plot(alpha_fitted_values, 
     alpha_dev_residuals,
     main = "Alpha - Residuals vs. Fitted",
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, col = "red", lty = 2)



## BETA ####

beta_data <- spline_data %>%
  filter(variant == "beta") %>%
  select(-c(variant, measure)) %>%
  pivot_longer(., cols = starts_with("20"), names_to = "date", values_to = "beta") %>%
  mutate(date = as.Date(date))

# FIT SPLINE FOR EACH COUNTRY IN VARIANT
beta_model <- vector("list", length = length(unique(beta_data$country)))
for (i in 1:length(unique(beta_data$country))) {
  
  beta_model[[i]] <- beta_data %>%
    filter(country == unique(beta_data$country)[i]) %>%
    gam(beta ~ s(as.numeric(date),
                  k = 10,
                  bs = "tp"
    ), 
    data = .)
}


# COEFFICIENTS

beta_coeff <- vector("list", length(unique(beta_data$country)))
for (i in 1:length(unique(beta_data$country))) {
  beta_coeff[[i]] <- unlist(beta_model[[i]]["coefficients"])
}
names(beta_coeff) <- unique(beta_data$country)


beta_output_data <- as.data.frame(beta_coeff) %>%
  t(.) %>%
  `colnames<-`(c("beta_intercept", "beta_x1", "beta_x2", "beta_x3", "beta_x4", "beta_x5", "beta_x6", "beta_x7", "beta_x8", "beta_x9")) %>%
  as.data.frame() %>%
  rownames_to_column("country")


# RESIDUALS

beta_dev_residuals <- vector("list", length(unique(beta_data$country)))
for (i in 1:length(unique(beta_data$country))) {
  beta_dev_residuals[[i]] <- unlist(residuals(beta_model[[i]], type = "deviance"))
}
beta_dev_residuals <- unlist(beta_dev_residuals)


# FITTED VALUES

beta_fitted_values <- vector("list", length(unique(beta_data$country)))
for (i in 1:length(unique(beta_data$country))) {
  beta_fitted_values[[i]] <- unlist(fitted.values(beta_model[[i]]))
}
beta_fitted_values <- unlist(beta_fitted_values)


# TEST RESIDUALS

qqnorm(beta_dev_residuals,
       main = "Beta")
qqline(beta_dev_residuals, col = "red")

hist(beta_dev_residuals,
     main = "Beta Deviance Residuals")
abline(v = 0, col = "red", lty = 2)

plot(beta_fitted_values, 
     beta_dev_residuals,
     main = "Beta - Residuals vs. Fitted",
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, col = "red", lty = 2)


## DELTA ####

delta_data <- spline_data %>%
  filter(variant == "delta") %>%
  select(-c(variant, measure)) %>%
  pivot_longer(., cols = starts_with("20"), names_to = "date", values_to = "delta") %>%
  mutate(date = as.Date(date))


# FIT SPLINE FOR EACH COUNTRY IN VARIANT
delta_model <- vector("list", length = length(unique(delta_data$country)))
for (i in 1:length(unique(delta_data$country))) {
  
  delta_model[[i]] <- delta_data %>%
    filter(country == unique(delta_data$country)[i]) %>%
    gam(delta ~ s(as.numeric(date),
                  k = 10,
                  bs = "tp"
    ), 
    data = .)
}



# COEFFICIENTS

delta_coeff <- vector("list", length(unique(delta_data$country)))
for (i in 1:length(unique(delta_data$country))) {
  delta_coeff[[i]] <- unlist(delta_model[[i]]["coefficients"])
}
names(delta_coeff) <- unique(delta_data$country)


delta_output_data <- as.data.frame(delta_coeff) %>%
  t(.) %>%
  `colnames<-`(c("delta_intercept", "delta_x1", "delta_x2", "delta_x3", "delta_x4", "delta_x5", "delta_x6", "delta_x7", "delta_x8", "delta_x9")) %>%
  as.data.frame() %>%
  rownames_to_column("country")


# RESIDUALS

delta_dev_residuals <- vector("list", length(unique(delta_data$country)))
for (i in 1:length(unique(delta_data$country))) {
  delta_dev_residuals[[i]] <- unlist(residuals(delta_model[[i]], type = "deviance"))
}
delta_dev_residuals <- unlist(delta_dev_residuals)


# FITTED VALUES

delta_fitted_values <- vector("list", length(unique(delta_data$country)))
for (i in 1:length(unique(delta_data$country))) {
  delta_fitted_values[[i]] <- unlist(fitted.values(delta_model[[i]]))
}
delta_fitted_values <- unlist(delta_fitted_values)


# TEST RESIDUALS

qqnorm(delta_dev_residuals,
       main = "Delta")
qqline(delta_dev_residuals, col = "red")

hist(delta_dev_residuals,
     main = "Delta Deviance Residuals")
abline(v = 0, col = "red", lty = 2)

plot(delta_fitted_values, 
     delta_dev_residuals,
     main = "Delta - Residuals vs. Fitted",
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, col = "red", lty = 2)



## GAMMA ####

gamma_data <- spline_data %>%
  filter(variant == "gamma") %>%
  select(-c(variant, measure)) %>%
  pivot_longer(., cols = starts_with("20"), names_to = "date", values_to = "gamma") %>%
  mutate(date = as.Date(date))

# FIT SPLINE FOR EACH COUNTRY IN VARIANT
gamma_model <- vector("list", length = length(unique(gamma_data$country)))
for (i in 1:length(unique(gamma_data$country))) {
  
  gamma_model[[i]] <- gamma_data %>%
    filter(country == unique(gamma_data$country)[i]) %>%
    gam(gamma ~ s(as.numeric(date),
                  k = 10,
                  bs = "tp"
    ), 
    data = .)
}



# COEFFICIENTS

gamma_coeff <- vector("list", length(unique(gamma_data$country)))
for (i in 1:length(unique(gamma_data$country))) {
  gamma_coeff[[i]] <- unlist(gamma_model[[i]]["coefficients"])
}
names(gamma_coeff) <- unique(gamma_data$country)


gamma_output_data <- as.data.frame(gamma_coeff) %>%
  t(.) %>%
  `colnames<-`(c("gamma_intercept", "gamma_x1", "gamma_x2", "gamma_x3", "gamma_x4", "gamma_x5", "gamma_x6", "gamma_x7", "gamma_x8", "gamma_x9")) %>%
  as.data.frame() %>%
  rownames_to_column("country")


# RESIDUALS

gamma_dev_residuals <- vector("list", length(unique(gamma_data$country)))
for (i in 1:length(unique(gamma_data$country))) {
  gamma_dev_residuals[[i]] <- unlist(residuals(gamma_model[[i]], type = "deviance"))
}
gamma_dev_residuals <- unlist(gamma_dev_residuals)


# FITTED VALUES

gamma_fitted_values <- vector("list", length(unique(gamma_data$country)))
for (i in 1:length(unique(gamma_data$country))) {
  gamma_fitted_values[[i]] <- unlist(fitted.values(gamma_model[[i]]))
}
gamma_fitted_values <- unlist(gamma_fitted_values)


# TEST RESIDUALS

qqnorm(gamma_dev_residuals,
       main = "Gamma")
qqline(gamma_dev_residuals, col = "red")

hist(gamma_dev_residuals,
     main = "Gamma Deviance Residuals")
abline(v = 0, col = "red", lty = 2)

plot(gamma_fitted_values, 
     gamma_dev_residuals,
     main = "Gamma - Residuals vs. Fitted",
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, col = "red", lty = 2)


## OMICRON ####

omicron_data <- spline_data %>%
  filter(variant == "omicron") %>%
  select(-c(variant, measure)) %>%
  pivot_longer(., cols = starts_with("20"), names_to = "date", values_to = "omicron") %>%
  mutate(date = as.Date(date))


# FIT SPLINE FOR EACH COUNTRY IN VARIANT
omicron_model <- vector("list", length = length(unique(omicron_data$country)))
for (i in 1:length(unique(omicron_data$country))) {
  
  omicron_model[[i]] <- omicron_data %>%
    filter(country == unique(omicron_data$country)[i]) %>%
    gam(omicron ~ s(as.numeric(date),
                  k = 10,
                  bs = "tp"
    ), 
    data = .)
}


# COEFFICIENTS

omicron_coeff <- vector("list", length(unique(omicron_data$country)))
for (i in 1:length(unique(omicron_data$country))) {
  omicron_coeff[[i]] <- unlist(omicron_model[[i]]["coefficients"])
}
names(omicron_coeff) <- unique(omicron_data$country)


omicron_output_data <- as.data.frame(omicron_coeff) %>%
  t(.) %>%
  `colnames<-`(c("omicron_intercept", "omicron_x1", "omicron_x2", "omicron_x3", "omicron_x4", "omicron_x5", "omicron_x6", "omicron_x7", "omicron_x8", "omicron_x9")) %>%
  as.data.frame() %>%
  rownames_to_column("country")


# RESIDUALS

omicron_dev_residuals <- vector("list", length(unique(omicron_data$country)))
for (i in 1:length(unique(omicron_data$country))) {
  omicron_dev_residuals[[i]] <- unlist(residuals(omicron_model[[i]], type = "deviance"))
}
omicron_dev_residuals <- unlist(omicron_dev_residuals)


# FITTED VALUES

omicron_fitted_values <- vector("list", length(unique(omicron_data$country)))
for (i in 1:length(unique(omicron_data$country))) {
  omicron_fitted_values[[i]] <- unlist(fitted.values(omicron_model[[i]]))
}
omicron_fitted_values <- unlist(omicron_fitted_values)


# TEST RESIDUALS

qqnorm(omicron_dev_residuals,
       main = "Omicron")
qqline(omicron_dev_residuals, col = "red")

hist(omicron_dev_residuals,
     main = "Omicron Deviance Residuals")
abline(v = 0, col = "red", lty = 2)

plot(omicron_fitted_values, 
     omicron_dev_residuals,
     main = "Omicron - Residuals vs. Fitted",
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, col = "red", lty = 2)



## VOI ####

VOI_data <- spline_data %>%
  filter(variant == "VOI") %>%
  select(-c(variant, measure)) %>%
  pivot_longer(., cols = starts_with("20"), names_to = "date", values_to = "VOI") %>%
  mutate(date = as.Date(date))


# FIT SPLINE FOR EACH COUNTRY IN VARIANT
VOI_model <- vector("list", length = length(unique(VOI_data$country)))
for (i in 1:length(unique(VOI_data$country))) {
  
  VOI_model[[i]] <- VOI_data %>%
    filter(country == unique(VOI_data$country)[i]) %>%
    gam(VOI ~ s(as.numeric(date),
                  k = 10,
                  bs = "tp"
    ), 
    data = .)
}


# COEFFICIENTS

VOI_coeff <- vector("list", length(unique(VOI_data$country)))
for (i in 1:length(unique(VOI_data$country))) {
  VOI_coeff[[i]] <- unlist(VOI_model[[i]]["coefficients"])
}
names(VOI_coeff) <- unique(VOI_data$country)


VOI_output_data <- as.data.frame(VOI_coeff) %>%
  t(.) %>%
  `colnames<-`(c("VOI_intercept", "VOI_x1", "VOI_x2", "VOI_x3", "VOI_x4", "VOI_x5", "VOI_x6", "VOI_x7", "VOI_x8", "VOI_x9")) %>%
  as.data.frame() %>%
  rownames_to_column("country")


# RESIDUALS

VOI_dev_residuals <- vector("list", length(unique(VOI_data$country)))
for (i in 1:length(unique(VOI_data$country))) {
  VOI_dev_residuals[[i]] <- unlist(residuals(VOI_model[[i]], type = "deviance"))
}
VOI_dev_residuals <- unlist(VOI_dev_residuals)


# FITTED VALUES

VOI_fitted_values <- vector("list", length(unique(VOI_data$country)))
for (i in 1:length(unique(VOI_data$country))) {
  VOI_fitted_values[[i]] <- unlist(fitted.values(VOI_model[[i]]))
}
VOI_fitted_values <- unlist(VOI_fitted_values)


# TEST RESIDUALS

qqnorm(VOI_dev_residuals,
       main = "VOI")
qqline(VOI_dev_residuals, col = "red")

hist(VOI_dev_residuals,
     main = "VOI Deviance Residuals")
abline(v = 0, col = "red", lty = 2)

plot(VOI_fitted_values, 
     VOI_dev_residuals,
     main = "Residuals vs. Fitted",
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, col = "red", lty = 2)



## TOTAL ####

total_data <- spline_data %>%
  filter(variant == "all") %>%
  select(-c(variant, measure)) %>%
  pivot_longer(., cols = starts_with("20"), names_to = "date", values_to = "n") %>%
  mutate(date = as.Date(date)) %>%
  left_join(population_data, by = "country")


# FIT SPLINE FOR EACH COUNTRY IN VARIANT
total_model <- vector("list", length = length(unique(total_data$country)))
for (i in 1:length(unique(total_data$country))) {

  total_model[[i]] <- total_data %>%
    filter(country == unique(total_data$country)[i]) %>%
    gam(n ~ s(as.numeric(date),
              k = 10,
              bs = "tp") +
              offset(log(population)),
      data = .,
      family = nb)
}


# COEFFICIENTS

total_coeff <- vector("list", length(unique(total_data$country)))
for (i in 1:length(unique(total_data$country))) {
  total_coeff[[i]] <- unlist(total_model[[i]]["coefficients"])
}
names(total_coeff) <- unique(total_data$country)


total_output_data <- as.data.frame(total_coeff) %>%
  t(.) %>%
  `colnames<-`(c("total_intercept", "total_x1", "total_x2", "total_x3", "total_x4", "total_x5", "total_x6", "total_x7", "total_x8", "total_x9")) %>%
  as.data.frame() %>%
  rownames_to_column("country")


# RESIDUALS

total_dev_residuals <- vector("list", length(unique(total_data$country)))
for (i in 1:length(unique(total_data$country))) {
  total_dev_residuals[[i]] <- unlist(residuals(total_model[[i]], type = "deviance"))
}
total_dev_residuals <- unlist(total_dev_residuals)


# FITTED VALUES

total_fitted_values <- vector("list", length(unique(total_data$country)))
for (i in 1:length(unique(total_data$country))) {
  total_fitted_values[[i]] <- unlist(fitted.values(total_model[[i]]))
}
total_fitted_values <- unlist(total_fitted_values)


# TEST RESIDUALS

qqnorm(total_dev_residuals,
       main = "Total Cases by Population")
qqline(total_dev_residuals, col = "red")

hist(total_dev_residuals,
     main = "Total Cases by Population Deviance Residuals")
abline(v = 0, col = "red", lty = 2)

plot(total_fitted_values, 
     total_dev_residuals,
     main = "Residuals vs. Fitted",
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, col = "red", lty = 2)


## COMBINE DATA ####

output_data <- alpha_output_data %>%
  full_join(beta_output_data, by = "country") %>%
  full_join(delta_output_data, by = "country") %>%
  full_join(gamma_output_data, by = "country") %>%
  full_join(omicron_output_data, by = "country") %>%
  full_join(VOI_output_data, by = "country") %>%
  full_join(total_output_data, by = "country")

output_data$country = stringr::str_replace_all(output_data$country, "\\.", " ")


#######################.

## SAVE DATA ####

output_data %>%
  saveRDS(., file = "spline_coeff_data.RDS")


#######################.

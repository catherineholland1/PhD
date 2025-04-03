library(dplyr)
library(ggplot2)

# Define parameters for the DLM
sigma_lambda_high <- 10
sigma_alpha <- 0

sigma_lambda <- 0
sigma_alpha_high <- 10

n <- 1000  # Number of time points

# Simulate time series for the DLM
simulate_dlm <- function(n, sigma_lambda, sigma_alpha) {
  lambda <- numeric(n)
  alpha <- rnorm(n, mean = 0, sd = sigma_alpha)
  epsilon <- rnorm(n, mean = 0, sd = sigma_lambda)
  
  lambda[1] <- rnorm(1, mean = 0, sd = sigma_lambda)  # Initial value of lambda
  
  for (t in 2:n) {
    lambda[t] <- lambda[t - 1] + alpha[t - 1] + epsilon[t]
  }
  
  return(lambda)
}

# Simulate the time series using the DLM
dlm_lambda_high <- simulate_dlm(n, sigma_lambda_high, sigma_alpha)
dlm_alpha_high <- simulate_dlm(n, sigma_lambda, sigma_alpha_high)


# Create data frame for plotting
time <- 1:n
dlm_data <- data.frame(time = time, rbind(cbind(series = 1, lambda = dlm_lambda_high),
                                      cbind(series = 2, lambda = dlm_alpha_high))
)

dlm_data %>%
  ggplot(., aes(x = time, y = lambda)) +
  geom_line() +
  facet_wrap(~series) +
  labs(title = "Simulated Time Series of lambda_t from DLM",
       x = "Time",
       y = expression(lambda[t])) +
  theme_minimal()


pdf("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD R Code/Time Series - COVID-19/DLM/example_dlm_plot.pdf", height = 3)
dlm_data %>%
  mutate(series = as.factor(series)) %>%
  ggplot(aes(x = time, y = lambda)) +
  geom_line(aes(col = series)) +
  facet_wrap(~series) +
  scale_color_manual(name = "Time Series", values = c("1" = "#66C2A5",
                                                      "2" = "#8DA0CB")) +
  labs(x = "Time", y = expression(lambda)) +
  theme(legend.position = "bottom") 
dev.off()

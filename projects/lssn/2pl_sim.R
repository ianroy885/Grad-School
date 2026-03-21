# 2PL Simulation Exercise
library(tidyverse)
library(ltm)
# 1. GenmultcompView# 1. Generate data
# 2. Fit model
# 3. Explore estimators


n <- 30 # individuals
k <- 10 # questions
totals <- n*k


alpha <- runif(k, 0.5, 2) # discrimination parameter
beta <- runif(k, -2, 3) # difficulty parameter

theta_mu <- 0 # mean person ability
theta_sig <- 1 # sd person ability

theta <- rnorm(n, theta_mu, theta_sig) # generate 30 ability parameters

discrim_ability <- theta %*% t(alpha) # discrimination * ability matrix
intercepts <- matrix(rep(beta, n), nrow = n, byrow = T) # get intercepts per person per item

probs <- plogis(intercepts + discrim_ability) # prob with logistic distribution := 1/(1+exp(.))
results <- ifelse(runif(totals) < probs, 1, 0) # matrix of 0-1 responses 


tidy_data <- results %>%
  as_tibble() %>%
  mutate(person_id = row_number()) %>%
  gather(item, respnse, -person_id) %>%
  mutate(item = as.numeric(as.factor(item)))

# do irt
modeled <- ltm(results ~ z1, IRT.param = T)
coef(modeled)
summary(modeled)

# plot characteristic curves
plot(modeled, type = "ICC", legend=T) # item characteristic curves
plot(modeled, type = "IIC", legend=T) # item info curves
plot(modeled, type = "IIC", legend=T, items=0) # test info function
information(modeled, c(-3,3)) # amount of info between -3,3 ability

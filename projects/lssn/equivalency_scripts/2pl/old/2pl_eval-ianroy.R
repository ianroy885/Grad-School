source('./projects/lssn/2pl_sim_functions.R')

n <- 1000
items <- 30

theta_true <- rnorm(n, 0, 1) # ability
a_true <- rlnorm(items, 0.5, 0.25) # discrimination
b_true <- rnorm(items, 0, 1) # difficulty

# response matrix
X <- matrix(NA, n, items)

# generate data and insert into the response matrix
for (i in 1:items){
  
  p <- p2pl(theta_true, a_true[i], b_true[i]) # 500 abilities, whats P(x=1) across the dif discrim/difficulty
  
  X[,i] <- rbinom(n, 1, p) # for each individual, whats P(X=1|prob of endorsement) 
  
}



fit <- estimate_2pl(X)

# Compare estimates to truth
cor(fit$a, a_true)      # discrimination recovery
cor(fit$b, b_true)      # difficulty recovery
cor(fit$theta, theta_true)  # ability recovery

# Plot
par(mfrow = c(1, 3))
plot(a_true, fit$a, main = "Discrimination (a)", xlab = "True", ylab = "Est"); abline(0,1)
plot(b_true, fit$b, main = "Difficulty (b)",     xlab = "True", ylab = "Est"); abline(0,1)
plot(theta_true, fit$theta, main = "Ability (θ)", xlab = "True", ylab = "Est"); abline(0,1)

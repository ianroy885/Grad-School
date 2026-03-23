n <- 500
items <- 20

theta_true <- rnorm(n, 0, 1) # ability
a_true <- rlnorm(items, 0.5, 0.25) # discrimination
b_true <- rnorm(items, 0, 1) # difficulty

# 2PL prob
p2pl <- function(theta, a, b){
  
  1 / (1 + exp(-a * (theta - b)))
  
}

# response matrix
X <- matrix(NA, n, items)

for (i in seq_len(items)){
  
  p <- p2pl(theta_true, a_true[i], b_true[i])
  X[,i] <- rbinom(n, 1, p)
  
}

log_lik <- function(params, X, theta){
  
  n_items <- ncol(X)
  a <- exp(params[1:n_items]) # keeps a>0
  b <- params[(n_items+1):(2*n_items)] 
  
  ll <- 0
  
  for (i in 1:n_items){
    p <- p2pl(theta, a[i],b[i])
    p <- pmax(p, 1e-10); pmin(p, 1 - 1e-10) # constrains to be [1e-10, 1-1e-10] since log(0) and log(1) go to inf
    ll <- ll + sum(X[,i] * log(p) + (1- X[,i]) * log(1-p)) # log likelihood of bernoulli across both items and people
  }
  print(paste("L =", ll))
  return(-ll) # negative for minimization
  
}

estimate_2pl <- function(X, tol=1e-4, maxit=50){
  
  n_persons <- nrow(X)
  n_items <- ncol(X)
  
  #initalize
  theta <- scale(rowMeans(X))[,1] # take average ability of each person and center it
  a <- rep(1, n_items)
  b <- rep(0, n_items)
  
  for (iter in 1:maxit){
    
    a_old <- a; b_old <- b
    
    # 1. Update item parameters given thera
    init_params <- c(log(a), b) # log for unconstrained optimization
    opt <- optim(
      par = init_params,
      fn = log_lik,
      X = X,
      theta = theta,
      method = "BFGS",
      control = list(maxit = 200)
    )
    
    a <- exp(opt$par[1:n_items])
    b <- opt$par[(n_items + 1):(2*n_items)]
    
    # 2. Update theta given item parameters
    theta <- sapply(1:n_persons, function(j){
      obj <- function(th){
        p <- p2pl(th, a, b)
        p <- pmax(p, 1e-10); p <- pmin(p, 1 - 1e-10)
        -sum(X[j,] * log(p) + (1 - X[j, ]) * log(1 - p))
      }
      optimize(obj, interval = c(-6,6))$minimum # finds the most likely theta per person given the item parameters
    })
    
    # center theta
    theta <- (theta - mean(theta)) / sd(theta)
    
    #check for convergence
    delta <- max(abs(a - a_old), abs(b - b_old))
    cat(sprintf("Iter %d | max delta = %.5f\n", iter, delta))
    if (delta < tol) break
    
  }
  
  list(a = a, b = b, theta = theta)
  
}

fit <- estimate_2pl(X)

#compare estimates to truth
cor(fit$a, a_true)
cor(fit$b, b_true)
cor(fit$theta, theta_true)

# Plot
par(mfrow = c(1, 3))
plot(a_true, fit$a, main = "Discrimination (a)", xlab = "True", ylab = "Est"); abline(0,1)
plot(b_true, fit$b, main = "Difficulty (b)",     xlab = "True", ylab = "Est"); abline(0,1)
plot(theta_true, fit$theta, main = "Ability (θ)", xlab = "True", ylab = "Est"); abline(0,1)
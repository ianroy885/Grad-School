# functions for the 2pl

# 2PL model for endorsing an item - alter this to get to lssn
p2pl <- function(theta, a, b){
  
  1 / (1 + exp(-a * (theta - b)))
  
}


# return -LL(theta | X)
log_lik <- function(a, b, theta, X){
  
  n_items <- ncol(X)
  
  a <- exp(a) # keeps a>0
  
  # initialize log likelihood
  ll <- 0
  
  # for each item...
  for (i in 1:n_items){
    
    # grab the prob of endorsement 
    p <- p2pl(theta, a[i],b[i])
    
    # constrain it to be within [1e-10, 1-1e-10] to avoid log(0) and log(1)
    p <- pmin(pmax(p, 1e-10), 1 - 1e-10)
    
    # take log likelihood using bernoulli pdf across both items (outside sum) and people (inside sum)
    ll <- ll + sum(X[,i] * log(p) + (1- X[,i]) * log(1-p))
  }
  
  #print(paste("L =", ll))
  return(-ll) # negative for minimization
  
}



# alternate version for one item at time
ll_item <- function(params, x_i, theta){
  a <- exp(params[1])#params[1]
  b <- params[2]
  p <- p2pl(theta,a,b)
  p <- pmin(pmax(p,1e-10),1-1e-10)
  ll <- -sum(x_i * log(p) + (1- x_i) * log(1-p))
  
  if(!is.finite(ll)){
    cat("Non finite LL. a=",a, "b=",b, "bad p:", any(!is.finite(p)),"\n")
  }
  
  return(ll)
}


estimate_2pl <- function(X, tol=1e-6, maxit=50){
  
  # browser()
  
  n_persons <- nrow(X)
  n_items <- ncol(X)
  
  # initalize baselines for each parameter as a starting point
  theta <- scale(rowMeans(X))[,1] 
  a <- rep(1, n_items)
  b <- rep(0, n_items)
  
  # cat('a:',a, 'log(a):',log(a), 'b:',b, theta)
  
  for (iter in 1:maxit){
    
    a_old <- a 
    b_old <- b
    
    # 1. Update item parameters given theta
    for (item in 1:n_items){
      
      # optimiize to find a and b 
      opt <- optim(
        par = c(log(a[item]), b[item]),
        fn = ll_item,
        x_i = X[, item],
        theta = theta,
        method = "BFGS", # waht is the optimization functin
        control = list(maxit = 200)
      )
      
      a[item] <- exp(opt$par[1])
      b[item] <- opt$par[2]
      
    }
    
    # 2. Update theta given item parameters
    theta <- sapply(1:n_persons, function(j){
      
      # obj will use the a and b from optim to get an estimate for theta
      obj <- function(th){
        
        
        p <- p2pl(th, a, b)
        p <- pmin(pmax(p,1e-10),1-1e-10)
        
        # for each specific person return -LL(theta | a,b)
        ll <- -sum(X[j,] * log(p) + (1 - X[j, ]) * log(1 - p))
        
        # add N(0,1) prior penalty on theta to pull away from extreme responders (all right or all wrong)
        prior <- -dnorm(th, 0, 1, log = TRUE)
        
        return(ll + prior)
        
      }
      
      # return an estimate for theta for each person 
      # within a bound of -6,6 since covers probability mass of theta ~ N(0,1)
      # taking the minimum since the -LL is used 
      optimize(obj, interval = c(-6,6))$minimum 
    }
    )
    
    # center theta
    theta <- (theta - mean(theta)) / sd(theta)
    
    #check for convergence
    delta <- max(abs(a - a_old), abs(b - b_old))
    cat(sprintf("Iter %d | max delta = %.5f\n", iter, delta))
    if (delta < tol) break
    
  }
  
  list(a = a, b = b, theta = theta)
  
}

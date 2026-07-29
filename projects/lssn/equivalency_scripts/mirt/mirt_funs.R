# goal is to make a 2pl from scratch

# Step 1: Define Model, likelihood function
# Step 2: Get LL(theta | data) using data, prob of endorsement via your model
# Step 3: Get the params
# Step 4: 

########## mirt package utilizes fixed quadrature EM for exploratory models and 
# Metropolis Hastings Robbins Monro method for exploratory, confirmatory, and polytomous models

######################
# Up to date version # 
######################

# define model
m_twopl <- function(theta, beta, alpha){
    
  model <- t(alpha) %*% theta - beta # expanding to fit x amount of latent traits
  
  1 / (1 + exp(-model))
}

m_twopl_vectorized <- Vectorize(m_twopl, c('alpha','beta'))

ll_bern_item <- function(X_j, params, theta){
  
  alpha <- params[1]
  beta <- params[2]
  
  # returns vector of probabilities for each person endorsing an item to use for LL
  p <- m_twopl(theta, beta, alpha)
  
  # constrain it to avoid log(0) or log(1)
  p[p==0] <- 1e-10
  p[p==1] <- 1 - 1e-10
  
  ### now need to get LL(P|X)
  
  # bernoulli -LL
  ll_item <- -sum(X_j * log(p) - X_j * log(1-p) + 1)

  if (is.nan(ll_item)|is.infinite(ll_item)){
    browser()
  }
  
  return(ll_item)
  
}

ll_bern_theta <- function(X_i, alpha, beta, theta){
  
  # returns vector of probabilities for one person endorsing each item to use for LL
  p <- m_twopl(theta, beta, alpha)
  
  # constrain it to avoid log(0) or log(1)
  p[p==0] <- 1e-10
  p[p==1] <- 1 - 1e-10
  
  ### now need to get LL(P|X) for one person
  
  # bernoulli -LL
  ll_thetas <- -sum(X_i * log(p) - X_i * log(1-p) + 1)
  
  if (is.nan(ll_thetas)|is.infinite(ll_thetas)){
    browser()
  }
  
  return(ll_thetas)
  
  # # returns matrix of probabilities for each person endorsing an item to use for LL
  # p <- m_twopl_vectorized(theta, beta, alpha)
  
  # # initialize empty LL matrix
  # ll_mat <- matrix(nrow = nrow(X), ncol = ncol(X))
  # 
  # ### now need to get LL(P|X)
  # for (i in 1:nrow(X)){
  #   for (j in 1:ncol(X)){
  #     
  #     # bernoulli negative LL
  #     ll_mat[i,j] <- -(X[i,j] * log(p[i,j]) - X[i,j] * log(1-p[i,j]) + 1)
  #     
  #   }
  # }
  # 
  # ll_thetas <- rowSums(ll_mat)
  
  
  # return(ll_thetas)
  
}

get_params <- function(input_data, max_iter = 500, tolerance = 1e-6){
 
  # get base stuff
  n <- nrow(input_data)
  k <- ncol(input_data)
  
  # get original estimates
  thetas <- apply(input_data, 1, mean)
  alphas <- rep(1, k)
  betas <- rep(0, k)
  
  delta <- 1
  
  for (iter in 1:max_iter){
    
    alpha_old <- alphas
    beta_old <- betas
    
    # start optimization loop 
    for (j in 1:k){
      # browser()

      # get estimates for alpha and beta
      model_params <- optim(
        par = c(alphas[j], betas[j]), # want to grab the alpas and betas for that item
        fn = ll_bern_item, # minimize log likelihood
        X_j = input_data[,j], # one item at a time
        theta = thetas,
        method = "BFGS",
        control = list(maxit = 1000)
      )
      
      # grab updated alpha and beta params
      alphas[j] <- model_params$par[1]
      betas[j] <- model_params$par[2]
 
    }
    
    # use updated params to get better estimate of theta
    # i.e. get new P(theta or p|X)
    
    theta_old <- thetas
    
    for (i in 1:n){
      
      thetas[i] <- optim(
        par = thetas[i], # want to grab the thetas
        fn = ll_bern_theta, # minimize -log likelihood
        X = input_data[i,], # one person at a time
        alpha = alphas,
        beta = betas,
        method = "BFGS",
        control = list(maxit = 1000)
      )$par[1]
    }
    
    delta <- max(c((alphas - alpha_old),(betas - beta_old)))
    
    
    if (delta < tolerance) break
    
  }
  
  return(list(a = alphas, b = betas, theta = thetas))
  
  
}

# finally <- get_params(sim_data)

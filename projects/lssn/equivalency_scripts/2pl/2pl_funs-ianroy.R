# goal is to make a 2pl from scratch

# Step 1: Define Model, likelihood function
# Step 2: Get LL(theta | data) using data, prob of endorsement via your model
# Step 3: Get the params
# Step 4: 


# define model
twopl <- Vectorize(
  
  function(theta, beta, alpha){
    
    1 / (1 + exp(-alpha * (theta - beta)))
    
  }, c('alpha','beta')
)

# likelihood fun of bernoulli
log_likelihood_bern <- #Vectorize(
  function(k,p){
  
  k * log(p) - k * log(1-p) + 1 # not doing any sums since we are looking at one success/failure at a time
  
  }#, 'k'
#)

# finding likelihood of theta | data
theta_likelihood <- function(input_data, alpha, beta, theta){
  
  num_items <- ncol(input_data)
  num_people <- nrow(input_data)

  
  # finding probability of endorsement given data
  p_endorse <- twopl(
    theta = theta, 
    alpha = alpha,
    beta = beta
  )
  
  
  # initialize results storage
  ll_abilities_matrix <- matrix(nrow = num_people, ncol = num_items)
  ll_abilities_vector <- matrix(nrow = num_people, ncol = 1)
  
  for (k in 1:num_people){
    
    # grab the ith person
    current_person <- input_data[k,]
    
    for (p in 1:num_items){
      
      # individual log likelihood of single response given the prob
      ll_abilities_matrix[k,p] <- k * log(p_endorse[k,p]) - k * log(1-p_endorse[k,p]) + 1
      
    }
    
    # compute actual log likelihood by summing across individual ones per item
    ll_abilities_vector[k,1] <- sum(ll_abilities_matrix[k,])
    
  }

  return(ll_abilities_vector)
  
  
}

# so i have a vector of n ability log likelihoods, but i am unsure on how best to proceed. so plan on learning more about optimization first before i continue











# define model
twopl <- function(theta, beta, alpha){
    
    1 / (1 + exp(-alpha * (theta - beta)))
}

twopl_vectorized <- Vectorize(
  function(theta, beta, alpha){
  
  1 / (1 + exp(-alpha * (theta - beta)))
  }, c('alpha','beta')
)

ll_bern_item <- function(X_j, params, theta){
  
  alpha <- params[1]
  beta <- params[2]
  
  # returns vector of probabilities for each person endorsing an item to use for LL
  p <- twopl(theta, beta, alpha)
  
  ### now need to get LL(P|X)
  
  # bernoulli -LL
  ll_item <- -sum(X_j * log(p) - X_j * log(1-p) + 1)
  
  return(ll_item)
  
}

ll_bern_theta <- function(X_i, alpha, beta, theta){
  
  # returns vector of probabilities for one person endorsing each item to use for LL - wait but this is returning 200 probs
  p <- twopl(theta, beta, alpha)
  
  # constrain p 
  p <- pmin(pmax(p,1e-10),1-1e-10)
  
  ### now need to get LL(P|X) for one person
  
  # bernoulli -LL
  ll_thetas <- -sum(X_i * log(p) - X_i * log(1-p) + 1)
  
  if(is.nan(ll_thetas | is.infinite(ll_thetas))){
    cat("bad p:", p, "ll:",ll_thetas,"\n")
  }
  
  return(ll_thetas)
  
  # # returns matrix of probabilities for each person endorsing an item to use for LL
  # p <- twopl_vectorized(theta, beta, alpha)
  
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
  
    # start optimization loop 
    for (i in 1:k){
      
      beta_old <- betas
      alpha_old <- alphas
      
      # get estimates for alpha and beta
      model_params <- optim(
        par = c(alpha_old[k], beta_old[k]), # want to grab the alpas and betas for that item
        fn = ll_bern_item, # minimize log likelihood
        X_j = input_data[,k], # one item at a time
        theta = thetas,
        method = "BFGS",
        control = list(maxit = 200)
      )
      
      # grab updated alpha and beta params
      alphas[k] <- model_params$par[1]
      betas[k] <- model_params$par[2]
      
    }
    
    # use updated params to get better estimate of theta
    # i.e. get new P(theta or p|X)
    
    theta_old <- thetas
    thetas <- rep(NA, n)
    
    for (i in 1:n){
      
      thetas[i] <- optim(
        par = theta_old[i], # want to grab the thetas
        fn = ll_bern_theta, # minimize log likelihood
        X = input_data[i,], # one person at a time
        alpha = alphas,
        beta = betas,
        method = "BFGS",
        control = list(maxit = 200)
      )$par[1]
    }
    
    delta <- max(c((alphas - alpha_old),(betas - beta_old)))
    print(delta)
    
    if (delta < tolerance) break
    
  }
  
  return(list(a = alphas, b = betas, theta = thetas))
  
  
}

finally <- get_params(sim_data)


finally$theta

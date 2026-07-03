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


# goal is to make a 2pl from scratch

# Step 1: Define Model
# Step 2: Using the data, get likelihood of the alpha, beta params
# Step 3: Using the likelihood, use optim to get the ability estimates since first deriv of likelihood at 0 is the MLE
# Step 4: 


# define model
twopl <- Vectorize(
  
  function(theta, beta, alpha){
    
    1 / (1 + exp(-alpha * (theta - beta)))
    
  }, c('alpha','beta')
)

# likelihood fun of bernoulli
log_likelihood_bern <- Vectorize(
  function(k,p){
  
  k * sum(log(p)) - k * sum(log(1-p)) + 1
  
  }, 'k'
)

# finding likelihood of theta | data
theta_likelihood <- function(input_data, alpha, beta, theta){
  
  # dim
  num_items <- ncol(input_data)
  num_people <- nrow(input_data)
  lik_matrix <- matrix(
      data = NA, 
      nrow=num_people, 
      ncol=1
    )
  
  # finding probability of endorsement given data
  p_endorse <- twopl(
    theta = theta, 
    alpha = alpha,
    beta = beta
  )
  
  # find likelihood of the theta for each person across all items given the p(X=1)
  ll <- apply(input_data, 2, function(x){ # for each item
    
    apply(p_endorse, 1, function(y){ # across each row (person)
      
      log_likelihood_bern(k = x, p = y)
      
    })
    
  })
  
}

test <- theta_likelihood(input_data = sim_data, alpha = difficulties, beta = discrims, theta = abilities)
test
dim(test)

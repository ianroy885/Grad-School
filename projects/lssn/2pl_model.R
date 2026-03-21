# building a 2pl from scratch

p_correct <- function(alpha, beta, personal){
  
  odds <- exp(alpha * (personal - beta))
  log_odds <- odds / (1 + odds)
  
  return(log_odds)
  
}

# Expectation maximization using gaussian quadrature
#   1. Assume normal distribution of plausible parameters
#   2. Break of the distribution into chunks with QuadPoints
#   3. Each chunk of the distribution will have different likelihood weights
#   4. Given the data, what is the likelihood of each parameter, then grab that 
#     (ie derive the expected number of correct responses within each chunk, so map the data to an estimate that maximizes that)
get_likelihood <- function(Y, alphas, betas, QuadPoints, Weights, NumPatterns){
  
  # initialize vector to store likelihoods
  PL <- rep(0, NumPatterns)
  num_items <- ncol(Y)
  prob_correct <- c(rep(na, num_items))
  prob_response <- c(rep(na, num_items))
  likelihood <- matrix(nrow=NumPatterns, ncol=num_items)
  
  
  # Loop over groups of people where patterns of responses get grouped together
  # 1. across each person (grouped together for easier estimation, but think of as person)
  # 2. Across each latent ability estimate (gaussian chunked up)
  # 3. Across each item
  #       For the item response (correct/incorrect), what's the P(Y==response | ability, difficulty, discrimination)
  #       Then multiply together all of these probabilities to get an overall estimate for a person,ability likelihood
  for (l in 1:NumPatterns){
    
    for (k in 1:length(QuadPoints)){
      
      
      for (i in 1:num_items){
        
 
        item <- Y[l, i] # get the specific item/person response 
        prob_correct[i] <- p_correct(alphas[i], betas[i], QuadPoints[k]) # estimate probability of correct given the ability, difficulty, discrim 
        prob_response[i] <- ifelse(item==1, prediction, 1 - prediction) # fill in with what was the probability of the response given
        
        if (is.na(likelihood[l,k])){
          
          likelihood[l,k] <- prob_response[i] # store likelihood estimate of the response given 
          
        }
        else{
          
          likelihood[l,k] <- likelihood[l,k] * prediction[i] # multiply the probability of response to each item within the person,ability pair
          
        }
        
      }
      
    }
    
    PL[l] <- sum(likelihood[l,] * Weights) # for a group/person, what's the integral over their distribution of likelihoods ?
    
  }
  
  return(list(likelihood, PL))
  
}



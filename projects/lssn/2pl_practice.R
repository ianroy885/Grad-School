library(ltm)
library(psych)
library(Rlab)



###########################
# Generating the 2pl data # 
###########################


n <- 200
items <- 10

abilities <- rnorm(n, 0, 1)
difficulties <- rnorm(items, 0, 1)
discriminations <- rnorm(items, 0.5, 0.25)

# 2pl model to estimate probably of getting a 1 given ability and difficulty
pl2_model <- Vectorize(
  
  function(theta, beta, alpha){
    
    1 / (1+exp(-alpha*(theta - beta)))
    
  }, 
  
  vectorize.args= c('beta','alpha') # only want to vectorize beta and alpha since will iterate over ability 
  
)

# matrix of P(X=1|ability_i, difficulty_j) for each person for each item
p_success <- sapply(abilities, function(x)
{
  pl2_model(x, difficulties, discriminations)
}

)

# generate the damn data
irt_data <- apply(
  
  # using matrix of probabilities
  p_success, 
  
  # across each row (item)
  1, 
  
  # find whether they got a 1 or 0 
  function(x){
    rbinom(
      n=200, 
      size=1,
      prob=x)
    
  })

#####################
# Playing with 2pl  # 
#####################

pl2_model <- ltm(irt_data ~ z1)
summary(pl2_model)
plot(pl2_model, type=c('ICC'))
plot(pl2_model, type=c('IIC'))

# assessing fit
item.fit(pl2_model)

# parameter recovery
ability_recovery <- ltm::factor.scores(pl2_model)
summary(ability_recovery$coef)
plot(ability_recovery) # recovers distribution well

# unidimensionality
unidimTest(pl2_model, irt_data)

library(ltm)
library(psych)
library(Rlab)



#############################
# Generating the rasch data # 
#############################


n <- 200
items <- 10

abilities <- rnorm(n, 0, 1)
difficulties <- rnorm(items, 0, 1)

# 1pl model to estimate probably of getting a 1 given ability and difficulty
rasch_model <- Vectorize(
  
  function(theta, beta){
  
    1 / (1+exp(-(theta - beta)))
  
  }, 
  
  vectorize.args= 'beta' # only want to vectorize beta since this will be person-by-person
  
)

# matrix of P(X=1|ability_i, difficulty_j) for each person for each item
p_success <- sapply(abilities, function(x)
    {
      rasch_model(x, difficulties)
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

#######################
# Playing with Rasch  # 
#######################

pl1_rasch <- rasch(irt_data)
summary(pl1_rasch)
coef(pl1_rasch)[,1] - difficulties
plot(pl1_rasch, type=c('ICC')) # gives info on how difficult each question is
plot(pl1_rasch, type=c('IIC')) # gives info on what abilities we are getting most info from for each item

# testing the fit
item.fit(pl1_rasch, simulate.p.value=T) # lower p would mean misfit

# estimating laten abilities (parameter recovery baby!)
recovered_ability <- ltm::factor.scores(pl1_rasch)
summary(recovered_ability$score.dat$`Item 10`)
plot(recovered_ability) # looks to match our generator of rnorm(0,1)

# testing for unidimensionality
unidimTest(pl1_rasch, irt_data) # if we simulated data with a discrim param, we would see it saying multidimensional even though it's uni 

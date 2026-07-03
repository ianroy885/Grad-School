# testing out the 2pl model

# data setup
n <- 200
items <- 10
abilities <- rnorm(n, 0, 1)
difficulties <- rnorm(items, 0, 1)
discrims <- rnorm(items, 0, 1)

# get prob of endorsement for each item, person
p_endorse <- twopl(abilities, difficulties, discrims)

# generate data with bernoulli
sim_data <- apply(p_endorse, 2, function(x){
  rbinom(
    n=n,
    size=1,
    prob=x
  )
})



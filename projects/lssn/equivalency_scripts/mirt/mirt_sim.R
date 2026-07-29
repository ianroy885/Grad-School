# testing out the 2pl model
library(ggplot2)
library(cowplot)
source('./projects/lssn/equivalency_scripts/2pl/2pl_funs.R')
# data setup
n <- 200
items <- 10
abilities <- rnorm(n, 0, 1)
difficulties <- rnorm(items, 0, 1)
discrims <- rnorm(items, 0, 1)

# get prob of endorsement for each item, person
p_endorse <- m_twopl_vectorized(abilities, difficulties, discrims)

# generate data with bernoulli
sim_data <- apply(p_endorse, 2, function(x){
  rbinom(
    n=n,
    size=1,
    prob=x
  )
})

model_params <- get_params(sim_data)
comparisons <- list(a = data.frame(
                    model = model_params$a,
                    pop = discrims
                    ),
                   b = data.frame(
                     model = model_params$b,
                     pop = difficulties
                   ),
                   theta = data.frame(
                     model = model_params$theta,
                     pop = abilities
                   )
                   )


alpha_plot <- ggplot(data = comparisons$a) + geom_point(aes(x = model, y = pop)) + 
  geom_smooth(aes(x = model, y = pop), method="lm", se=F, col='red') + 
  xlab('Model Estimate') + ylab('Pop Value') + labs(title = 'Alphas')

betas_plot <- ggplot(data = comparisons$b) + geom_point(aes(x = model, y = pop)) + 
  geom_smooth(aes(x = model, y = pop), method="lm", se=F, col='red') + 
  xlab('Model Estimate') + ylab('Pop Value') + labs(title = 'Betas')

thetas_plot <- ggplot(data = comparisons$theta) + geom_point(aes(x = model, y = pop)) + 
  geom_smooth(aes(x = model, y = pop), method="lm", se=F, col='red') + 
  xlab('Model Estimate') + ylab('Pop Value') + labs(title = 'Thetas')

plot_grid(alpha_plot, betas_plot, thetas_plot, ncol=3)
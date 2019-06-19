library(tidyverse)
library(rstan)
library(loo)
library(plyr)
library(magrittr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Sampling parameters
n_cores = 3
n_chain = 3
iterations = 6000

# Load df
setwd("SNOD/")
load(file = "./data/snod.RData") 
set.seed(125)
#d %<>% filter(image %in% sample(unique(image),20),
#         subj %in% sample(unique(subj),20)) %>%
#  mutate(image = as.integer(factor(image)),
#         subj = as.integer(factor(subj)))

#d %<>% filter(subj %in% sample(unique(subj), 10), 
#              image %in% sample(unique(image), 10)) %>%
#  mutate(subj = as.integer(factor(subj)),
#         image = as.integer(factor(image)))

d %>% group_by(subj, image) %>% dplyr::count() %>% ungroup() -> nB
maxB <- max(nB$n)

d %>% select(subj, image) %>% unique() %>% group_by(subj) %>% dplyr::count() -> Nwords

nS <- length(unique(d$subj))
nI <- length(unique(d$image))

subj <- d$subj
image <- d$image

N <- nrow(d)

S <- rep(1:nS, times = Nwords$n)

# Matrix for keystrokes
y <- array(-1, c(length(S), maxB))
n = 0
for (i in 1:nS) {
  for(j in nB[nB$subj == i,]$image){
    n = n + 1
    subj_data <- filter(d, subj == i, image == j)
    y[n, 1:nB[nB$subj == i & nB$image == j,]$n] <- subj_data$IKI  
  }
}
head(y)

# Predictor


# Turn data into list
dat <- within( list(), {
  nS <- nS
  nI <- nI
  
  subj <- subj
  image <- image
  spell_H <- d$spell_H
  nKeys <- nB$n
  nK <- length(nB$n)
  maxB <- maxB
  
  y <- y
  N <- N
} );str(dat)


# Initialise start values
start <- 
  function(chain_id = 1){
    list(sigma = .1
         , phi = 0
         , alpha_raw = rep(0,2)
         , alpha_mu = 5
         , alpha_sigma = .1
         , b_spell = 0
         , spell_mu = 0
         , spell_sigma = .1
         , u = rep(0.1, dat$nS)
         , sigma_u = .1
         , w = rep(0.1, dat$nI)
         , sigma_w = .1
         , alpha = chain_id
    )
  }

start_ll <- lapply(1:n_chain, function(id) start(chain_id = id) )

# --------------
# Stan models ##
# --------------
#---- 
# Mixture of two gaussians with unequal variance
#---- 
# Load model
ar1 <- stan_model(file = "stanin/ARspell.stan")

# Parameters to omit in output
omit <- c("RE")

# Fit model
m <- sampling(ar1, 
              data = dat,
              init = start_ll,
              iter = iterations,
              warmup = iterations/2,
              chains = n_chain, 
              cores = n_cores,
              refresh = 250,
              save_warmup = FALSE, # Don't save the warmup
              include = FALSE, # Don't include the following parameters in the output
              pars = omit,
              thin = 1,
              seed = 81,
              control = list(max_treedepth = 16,
                             adapt_delta = 0.99,
                             stepsize = .001)
)


# Save model
saveRDS(m, 
        file = "stanout/ARspell.rda",
        compress = "xz")

# Traceplots
#names(m)
param <- c("alpha", "b_spell", "phi", "sigma", "alpha_mu", "alpha_sigma", "alpha_raw") 
summary(print(m, pars = param, probs = c(.025,.975)))
traceplot(m, param, inc_warmup = F)
#traceplot(m, "u", inc_warmup = F)






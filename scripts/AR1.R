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

# This script runs a multi-level mixed autoregression model.
# The response related effects are estimated separaletly for keystroke data for
# first keystroke [1] and withih words [2]. 
# Autoregression for within work keystrokes are captured in phi.
# Random intercets are included for subjs and images.

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


# Turn data into list
dat <- within( list(), {
  nS <- nS
  nI <- nI
  
  subj <- subj
  image <- image
  
  nKeys <- nB$n
  nK <- length(nB$n)
  maxB <- maxB
  
#  wordinit <- d$ispwordinit
  spRT <- d$spRT
  name_H <- d$name_H
  spell_H <- d$spell_H
  word_length <- d$word.length
  word_frequency <- d$mc_response_wdfrq
  letter_frq <- d$letter_frq
  dig_frq_pre <- d$dig_frq_pre
  dig_frq <- d$dig_frq
  dig_frq_cond <- d$dig_frq_cond
  
  y <- y
  N <- N
} );str(dat)


# Initialise start values
start <- 
  function(chain_id = 1){
    list(sigma = .1
         , phi = 0
         , b_name = rep(0, 2)
         , name_mu = 0
         , name_sigma = .1
         , b_spell = rep(0, 2)
         , spell_mu = 0
         , spell_sigma = .1
         , b_spoken = rep(0, 2)
         , spken_mu = 0
         , spoken_sigma = .1
         , b_wlen = rep(0, 2)
         , wlen_mu = 0
         , wlen_sigma = .1
         , b_freq = rep(0, 2)
         , freq_mu = 0
         , freq_sigma = .1
         , b_dig_freq_init = rep(0, 2)
         , dfi_mu = 0
         , dfi_sigma = .1
         , b_dig_freq_mid = rep(0, 2)
         , dfm_mu = 0
         , dfm_sigma = .1
         , b_dig_cond_prob = rep(0, 2)
         , dcp_mu = 0
         , dcp_sigma = .1
         , b_letter_freq = rep(0, 2)
         , letter_mu = 0
         , letter_sigma = .1
         , alpha_raw = rep(0,2)
         , alpha_mu = 5
         , alpha_sigma = .1
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
ar1 <- stan_model(file = "stanin/AR1.stan")

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
        file = "stanout/AR1.rda",
        compress = "xz")

# Traceplots

slopes <- names(m)[grepl("b_",names(m))]
param <- c("alpha", slopes, "phi", "sigma") 
summary(print(m, pars = param, probs = c(.025,.975)))
traceplot(m, param, inc_warmup = F)
#traceplot(m, "u", inc_warmup = F)






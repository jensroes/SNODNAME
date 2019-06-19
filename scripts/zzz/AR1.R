# Load packages
library(tidyverse)
library(rstan)
library(loo)
library(plyr)
library(magrittr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Sampling parameters
set.seed(125)
n_cores = 3
n_chain = 3
iterations = 10000

setwd("Frontline/")


# Load df
d <- read_csv("data/ct.csv") %>% 
  filter(component == "Consonants", 
        # component == "Tapping", 
         !is.na(age), 
         age >= 50, 
         session == 1,
         !is.na(IKI), 
         IKI > 0, 
         target == 1) %>%
  mutate(logfile = paste0(subj, logfile),
         subj = as.integer(factor(logfile)),
         bigram = as.integer(factor(bigram))) %>% 
  select(-logfile) %>%
  select(bigram, IKI, subj) #%>%
#  filter(subj %in% 1:10)

d %<>% group_by(subj) %>% dplyr::mutate(bid = 1:n())
d %>% group_by(subj) %>% dplyr::count() %>% ungroup() -> nB
nB %>% summarise(max = as.integer(max(n))) -> maxB
nS <- length(unique(d$subj))
subj <- d$subj
bigram <- d$bigram
N <- nrow(d)
y <- array(-1, c(nS, maxB))

for (i in 1:nS) {
  subj_data <- filter(d, subj == i)
  y[i, 1:nB[nB$subj == i,]$n] <- subj_data$IKI  
}

dat <- within( list(), {
  nS <- nS
  nB <- nB$n
  maxB <- maxB$max
  y <- y
  subj <- subj
  N <- N
} );str(dat)


# Initialise start values
start <- 
  function(chain_id = 1){
    list(sigma = .1
         , phi = 0
         , alpha_raw = 0
         , alpha_mu = 5
         , alpha_sigma = .1
         , u = rep(0.1, dat$nS)
         , sigma_u = 1
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
                             adapt_delta = 0.99)
)


# Save model
saveRDS(m, 
        file = "stanout/AR1.rda",
        compress = "xz")

# Traceplots
#names(m)
param <- c("alpha","phi", "sigma", "alpha_mu", "alpha_sigma", "alpha_raw") 
summary(print(m, pars = param, probs = c(.025,.975)))
traceplot(m, param, inc_warmup = F)
#traceplot(m, "u", inc_warmup = F)

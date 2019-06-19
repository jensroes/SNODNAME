library(loo)
library(tidyverse)
library(rstan)


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
  select(bigram, IKI, subj)

path <- "Frontline/stanout/"
(files <- dir(path))



m <- readRDS(paste0(path, files[4] ))
names(m)
pars <- c("alpha",  "phi", "sigma") #"sigmap_e", "sigma_e", "sigma_diff"
print(m, pars, probs = c(.025, .975))
traceplot(m, pars)
ll <- extract_log_lik(m, merge_chains = F) 
r_eff <- relative_eff(exp(ll)) 
loo(ll, r_eff = r_eff, cores =2)


#Theta <- extract(m, 'theta')
#Theta <- unlist(Theta, use.names=FALSE)
y_pred <- rstan::extract(m, 'y_tilde')
#y_pred <- unlist(y_pred, use.names=FALSE)

y_pred$y_tilde %>% t() %>% as_tibble() %>%
  select(paste0("V",1:500)) %>%
  gather(iteration, y_tilde) %>%
  mutate(y_tilde = exp(y_tilde)) %>%
  ggplot(aes(x = y_tilde, group = iteration)) +
  geom_density(alpha = .1, color = "grey") +
  geom_density(data = d, aes(x = IKI, group = NULL), color = "red") +
  scale_x_continuous(limits = c(0, 15000)) +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank())

y_pred$y_tilde %>% t() %>% as_tibble() %>%
  gather(iteration, y_tilde) -> tmp
which(tmp$y_tilde < 0)


m <- readRDS(paste0(path, files[2] ))
pars <- c("alpha", "delta", "phi", "theta", "sigma")
print(m, pars, probs = c(.025, .975))
traceplot(m, pars)
ll <- extract_log_lik(m, merge_chains = F) 
r_eff <- relative_eff(exp(log_lik)) 
loo(ll, r_eff = r_eff, cores = 2)


m <- readRDS(paste0(path, files[3] ))
pars <- c("alpha","phi", "sigma")
print(m, pars, probs = c(.025, .975))
traceplot(m, pars)
ll <- extract_log_lik(m, merge_chains = F) 
r_eff <- relative_eff(exp(log_lik)) 
loo(ll, r_eff, cores =2)


m <- readRDS(paste0(path, files[4] ))
pars <- c("alpha","phi", "theta", "sigma", "sigma_diff", "sigma_e", "sigmap_e")
print(m, pars, probs = c(.025, .975))
traceplot(m, pars)
ll <- extract_log_lik(m, merge_chains = F) 
r_eff <- relative_eff(exp(log_lik)) 
loo(ll, r_eff, cores =2)


m <- readRDS(paste0(path, files[5] ))
pars <- c("alpha","phi", "delta", "theta", "sigma", "sigma_diff", "sigma_e", "sigmap_e")
print(m, pars, probs = c(.025, .975))
traceplot(m, pars)
ll <- extract_log_lik(m, merge_chains = F) 
r_eff <- relative_eff(exp(log_lik)) 
loo(ll, r_eff, cores =2)


m <- readRDS(paste0(path, files[6] ))
pars <- c("alpha","phi", "gamma", "delta", "theta", "sigma", "sigma_diff", "sigma_e", "sigmap_e")
print(m, pars, probs = c(.025, .975))
traceplot(m, pars)
ll <- extract_log_lik(m, merge_chains = F) 
r_eff <- relative_eff(exp(log_lik)) 
loo(ll, r_eff, cores =2)


# The equal variance model doesn't make sense at the moment. Sigma can't be the same for both components the way it is specificed
m <- readRDS(paste0(path, files[7] ))
pars <- c("alpha","phi", "gamma", "delta", "theta", "sigma")
print(m, pars, probs = c(.025, .975))
traceplot(m, pars)
ll <- extract_log_lik(m, merge_chains = F) 
r_eff <- relative_eff(exp(log_lik)) 
loo(ll, r_eff, cores =2)


m <- readRDS(paste0(path, files[8] ))
pars <- c("alpha","phi", "gamma",  "sigma")
print(m, pars, probs = c(.025, .975))
traceplot(m, pars)
ll <- extract_log_lik(m, merge_chains = F) 
r_eff <- relative_eff(exp(log_lik)) 
loo(ll, r_eff, cores =2)

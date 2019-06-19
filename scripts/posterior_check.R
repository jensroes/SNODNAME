library(loo)
library(tidyverse)
library(rstan)

path <- "stanout/"
(files <- dir(path))

#m <- readRDS(paste0(path, files[5] ))
m <- readRDS(paste0(path, "ARspell2MoG2.rda" ))

#names(m)
pars <- c("alpha", "b_spell",  "delta", "theta",  "phi",   "sigma","sigmap_e", "sigma_e", "sigma_diff") #"
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



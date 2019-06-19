library(loo)
library(tidyverse)
library(rstan)
library(ggridges)
library(rethinking)

# Data
load(file = "./data/snod.RData") 
set.seed(125)
d %<>% filter(image %in% sample(unique(image),20),
              subj %in% sample(unique(subj),20)) %>%
  mutate(image = as.integer(factor(image)),
         subj = as.integer(factor(subj)))


# Model
path <- "stanout/"
(files <- dir(path))
#m <- readRDS(paste0(path, files[5] ))
m <- readRDS(paste0(path, "AR1.rda" ))

#names(m)
slopes <- names(m)[grepl("b_",names(m))]
param <- c("alpha", slopes, "phi", "sigma") 
summary(print(m, pars = param, probs = c(.025,.975)))
traceplot(m, param, inc_warmup = F)

ll <- extract_log_lik(m, merge_chains = F) 
r_eff <- relative_eff(exp(ll)) 
loo(ll, r_eff = r_eff, cores =2)


msamps <- rstan::extract(m,pars = param) %>% 
  as.data.frame() %>% as_tibble() 

write_csv(msamps, "results/posteriorAR1.csv")

# Plot predicted data

y_pred <- rstan::extract(m, 'y_tilde')
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

(IVlevs = unique(gsub(slopes,pattern = "\\[|\\]|1|2", replacement = "")))
(IVlevs = IVlevs[c(2,1,3:9)])
IVlabs = c("Spoken naming RT"
           ,"Name diversity (H)" 
           ,"Spell diversity (H)"
           ,"Word length"
           ,"Word frequency"
           ,"Digraph frequency (initial)"
           ,"Digraph frequency (mid)"
           ,"Digraph conditional probability"
           ,"Letter frequency")


samps <- rstan::extract(m,pars = param) %>% 
  as.data.frame() %>% as_tibble() %>%
  gather(Parameter, value) %>%
  mutate(Parameter = gsub('^\\.|\\.$', '', Parameter)) %>%
  separate(col = Parameter, into = c("Parameter", "Pos"), sep = "\\.") %>%
  mutate(Pos = ifelse(is.na(Pos), 3, Pos)) %>%
  group_by(Parameter) %>%
  dplyr::mutate(id = 1:n()) %>%
  spread(Parameter, value) %>%
  ungroup() %>%
  select(-id)



pop_effs = samps %>%
  rename(Intercept = alpha) %>%
  group_by(Pos) %>%
  mutate_at(vars(starts_with('b_')), list(~ (exp(Intercept + .) - exp(Intercept)))) %>%
  select(starts_with('b_'), Pos) %>%
  gather(IV,param,-Pos) %>% 
  ungroup() %>%
  filter(!is.na(param)) %>%
  mutate(Pos = as.integer(Pos)) %>%
  mutate(IV = factor(IV, levels = rev(IVlevs), labels = rev(IVlabs))
         ,RT_IKI = factor(Pos, levels = c('1','2')
                          , labels = c('Response onset latency (RT)', 'Keystroke latency (IKI)')))


# do plot ----

s = 18 #base size
m = .3 # label and point size multiplier

my_theme = theme_light(base_size = s) +
  theme(
    axis.text.x = element_text(size = s),
    #axis.text.x = element_blank(),
    axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)
                                ,size = s*.8),
    #axis.title.x = element_blank(),
    axis.text.y = element_text(size = s),
    #axis.title.y = element_text(size = s),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    #axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    #axis.line = element_line(),
    
    
    #legend.position = 'right',
    legend.position = c(.8,.7),
    plot.title = element_text(hjust = 0.5,
                              size = s*.9),
    #legend.text = element_text(size = s,
    #                           margin = margin(t = 10, r = 0, b = 10, l = 10)),
    legend.title = element_blank(),
    legend.spacing.x = unit(2, "char"),
    legend.key.size = unit(2, 'char'),
    #legend.spacing.y = unit(3, "char"),
    
    strip.text = element_text(size = s),
    
    panel.border=element_blank()
    #panel.grid = element_blank()
  )


#pop_effs %>%
#  ggplot(aes(y = IV, x = param)) +
#  geom_density_ridges(fill = 'dark grey') +
#  geom_vline(aes(xintercept = 0)
#             ,linetype = 'dashed') +
#  facet_wrap(~RT_IKI, scales ='free_x') +
#  my_theme +
#  labs(x = 'Effect (ms) of 1SD increase in predictor')

#ggsave('./plots/posterior_density.png', height = 16, width = 30, units = 'cm')

# do CI plot ----

pop_effs %>% 
  filter(!(Pos == 1 & IV %in% c("Digraph conditional probability", "Digraph frequency (mid)"))) %>%
  group_by(RT_IKI, IV) %>% 
  summarise(Median = median(param)
            ,lo = PI(param,.95)[1]
            ,hi = PI(param,.95)[2]) %>% 
  ggplot(aes(y = IV, x = Median)) +
  geom_vline(aes(xintercept = 0),linetype = 'dashed') +
  geom_errorbarh(aes(xmin = lo, xmax = hi, height = 0)
                 ,size = 2, colour = 'dark grey') +
  geom_point(size = 4, colour = 'dark red') +
  facet_wrap(~RT_IKI, scales ='free_x') +
  my_theme +
  labs(x = 'Effect (ms) of 1SD increase in predictor, 95% CrI')

ggsave('results/CIs.png', height = 16, width = 30, units = 'cm')

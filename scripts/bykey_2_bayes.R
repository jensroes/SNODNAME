# setup ----
library(tidyverse)
library(magrittr)
library(brms)


load("./data/snod_EN_keydata_2.Rdata");Df


Df %<>% filter(mc_response_wdlen == 1)

length(unique(Df$image)) #number of images out of 260 with single word modal responses

# within word effects ----

stdz <- function(x){(x-mean(x, na.rm = T))/sd(x, na.rm = T)}

Df2 = Df %>% 
  mutate(letter_frq =  stdz(log(letter_frq))
         ,dig_frq_pre = stdz(log(dig_frq))
         ,dig_frq = stdz(log(lag(dig_frq,1)))
         ,dig_frq_cond = stdz(log(lag(dig_frq_cond,1)))
         ,spRT = stdz(sp_RT)
         ,word.length = stdz(word.length)
         ,mc_response_wdfrq = stdz(log(mc_response_wdfrq))
         ,IKI = log(IKI)) %>%
  filter(!is.na(dig_frq)
         ,!is.na(dig_frq_pre)
         ,!is.na(letter_frq)) %>%  
  filter(is_mc_response) %>% 
  filter(IKI <= mean(IKI, na.rm = T) + sd(IKI, na.rm = T)*3)


Df2samp = Df2 %>% 
  filter(subno %in% sample(unique(Df2$subno),10)
         ,image %in% sample(unique(Df2$image),10))

bm1 = brm(IKI ~ 1 + spRT + name_H 
          + letter_frq + dig_frq_pre 
          + dig_frq + dig_frq_cond
          + spell_H + word.length + mc_response_wdfrq
          + (1|subno) + (1|image)
         , prior = set_prior('normal(0, 100)', class = 'b')
         , data = Df2
         , cores = 36)

save(bm1,file = './models/bm1.Rdata')


# prior_summary(bm1)
# summary(bm1)
# plot(bm1)
# plot(bm1, pars = '^b')
# stanplot(bm1, pars = c('spRT','name_H','letter_frq',
#                        'dig_frq_pre','dig_frq','dig_frq_cond',
#                        'spell_H','word.length','mc_response_wdfrq'))
# pp_check(bm1)
# marginal_effects(bm1)
# 
# pop_effs = posterior_samples(bm1) %>%
#   mutate(b_spRT = exp(b_Intercept + b_spRT) - exp(b_Intercept)
#          ,name_H = exp(b_Intercept + b_spRT) - exp(b_Intercept))

ps = posterior_samples(bm1)

write_csv(ps,file = './models/bm1_posterior_samples')



# word-initial effects ----

stdz <- function(x){(x-mean(x, na.rm = T))/sd(x, na.rm = T)}
Df3 = Df %>% 
  group_by(subno, image) %>% 
  mutate(M_dig_frq = ifelse(letter.id <= 2, NA, lag(dig_frq,1)) #ignore first digraph when getting mean for word
         ,M_dig_frq = mean(M_dig_frq, na.rm = T)
         ,M_dig_frq_cond = mean(lag(dig_frq_cond,1), na.rm = T)) %>% 
  ungroup() %>% 
  mutate(letter_frq =  stdz(log(letter_frq))
         ,dig_frq_pre = stdz(log(dig_frq))
         ,M_dig_frq = stdz(log(M_dig_frq))
         ,M_dig_frq_cond = stdz(log(M_dig_frq_cond))
         ,spRT = stdz(sp_RT)
         ,IKI = log(IKI)
         ,word.length = stdz(word.length)
         ,mc_response_wdfrq = stdz(log(mc_response_wdfrq))) %>%
  filter(!is.na(dig_frq_pre)
         ,!is.na(letter_frq)) %>%  
  filter(is_mc_response
         ,letter.id == 1) %>% 
  filter(IKI <= mean(IKI, na.rm = T) + sd(IKI, na.rm = T)*3)


Df3samp = Df3 %>% 
  filter(image %in% sample(unique(image),20),
         subno %in% sample(unique(subno),20))

# bivariate cors
round(Df3 %>% select(spRT, name_H, letter_frq, dig_frq_pre 
              , M_dig_frq, M_dig_frq_cond
              , spell_H, word.length, mc_response_wdfrq) %>% 
  cor(.),3)

bm2 = brm(IKI ~ 1 + spRT + name_H 
          + letter_frq + dig_frq_pre 
          + M_dig_frq + M_dig_frq_cond
          + spell_H + word.length + mc_response_wdfrq
          + (1|subno) + (1|image)
          , prior = set_prior('normal(0, 100)', class = 'b')
          , data = Df3
          , cores = 36)

save(bm2,file = './models/bm2.Rdata')



# stanplot(bm2, pars = c('spRT','name_H','letter_frq',
#                        'dig_frq_pre','M_dig_frq','M_dig_frq_cond',
#                        'spell_H','word.length','mc_response_wdfrq'))
# 
# 
ps2 = posterior_samples(bm2)
save(ps2,file = './models/bm1_posterior_samples')
# 
# plot(bm2)

# pop_effs = posterior_samples(bm2) %>%
#   rename(Intercept = b_Intercept) %>%
#   mutate_at(vars(starts_with('b_')), list(~ (exp(Intercept + .) - exp(Intercept)))) %>% 
#   select(starts_with('b_')) %>% 
#   gather(IV,param)

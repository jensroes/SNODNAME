# setup ----
library(tidyverse)
library(magrittr)


load("./data/snod_EN_keydata_2.Rdata");Df


Df %<>% filter(mc_response_wdlen == 1)

length(unique(Df$image)) #number of images out of 260 with single word modal responses

stdz <- function(x){(x-mean(x, na.rm = T))/sd(x, na.rm = T)}

Df2 = Df %>% 
  mutate(letter_frq =  stdz(log(letter_frq))
         ,dig_frq_pre = stdz(log(dig_frq))
         ,dig_frq = stdz(log(lag(dig_frq,1)))
         ,dig_frq_cond = stdz(log(lag(dig_frq_cond,1)))
         ,spRT = stdz(sp_RT)
         ,word.length = stdz(word.length)
         ,mc_response_wdfrq = stdz(log(mc_response_wdfrq))) %>%
  filter(!is.na(dig_frq)
         ,!is.na(dig_frq_pre)
         ,!is.na(letter_frq)) %>%  
  filter(is_mc_response) 

Df3 = Df %>% 
  group_by(subno, image) %>% 
  mutate(M_dig_frq = ifelse(letter.id <= 2, NA, lag(dig_frq,1)) #ignore first digraph when getting mean for word
         ,M_dig_frq = mean(M_dig_frq, na.rm = T)
         ,M_dig_frq_cond = mean(lag(dig_frq_cond,1), na.rm = T)) %>% 
  ungroup() %>% 
  mutate(letter_frq =  stdz(log(letter_frq))
         ,dig_frq_pre = stdz(log(dig_frq))
         ,dig_frq = stdz(log(M_dig_frq))
         ,dig_frq_cond = stdz(log(M_dig_frq_cond))
         ,spRT = stdz(sp_RT)
         ,word.length = stdz(word.length)
         ,mc_response_wdfrq = stdz(log(mc_response_wdfrq))) %>%
  filter(!is.na(dig_frq_pre)
         ,!is.na(letter_frq)) %>%  
  filter(is_mc_response
         ,letter.id == 1) 

Df2 %>% bind_rows(Df3) %>%
  arrange(subno, image, letter.id) %>%
  mutate(image = as.integer(factor(image)),
         subj = as.integer(factor(subno))) %>%
  select(image, word, subj, IKI, ispwordinit, letter.id, spRT, 
         name_H, spell_H, word.length, mc_response_wdfrq, dig_frq,
         dig_frq_pre, dig_frq_cond, letter_frq) -> d; d


save(d, file = "./data/snod.RData")


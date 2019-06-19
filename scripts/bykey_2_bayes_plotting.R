ps1 = read_csv(file = './models/bm1_posterior_samples.csv') %>% mutate(RT_IKI = 'IKI')
ps2 = read_csv(file = './models/bm2_posterior_samples.csv') %>% 
  mutate(RT_IKI = 'RT') %>%
  rename(b_dig_frq = b_M_dig_frq
         ,b_dig_frq_cond = b_M_dig_frq_cond)

ps = rbind(ps2,ps1)


IVlabs = c("Spoken naming RT"
            ,"Name diversity (H)" 
            ,"Spell diversity (H)"
            ,"Word length"
            ,"Word frequency"
            ,"Digraph frequency (initial)"
            ,"Digraph frequency (mid) [Mean]"
            ,"Digraph conditional probability [Mean]"
            ,"Letter frequency")

IVlevs = c("b_spRT"
           ,"b_name_H"
           ,"b_spell_H"
           ,"b_word.length"
           ,"b_mc_response_wdfrq"
           ,"b_dig_frq_pre"
           ,"b_dig_frq"
           ,"b_dig_frq_cond"
           ,"b_letter_frq")

pop_effs = ps %>%
  rename(Intercept = b_Intercept) %>%
  mutate_at(vars(starts_with('b_')), list(~ (exp(Intercept + .) - exp(Intercept)))) %>%
  select(starts_with('b_'), RT_IKI) %>%
  gather(IV,param, -RT_IKI) %>% 
  mutate(IV = factor(IV, levels = rev(IVlevs), labels = rev(IVlabs))
         ,RT_IKI = factor(RT_IKI, levels = c('RT','IKI')
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


library(ggridges)
pop_effs %>%
  ggplot(aes(y = IV, x = param)) +
  geom_density_ridges(fill = 'dark grey') +
  geom_vline(aes(xintercept = 0)
             ,linetype = 'dashed') +
  facet_wrap(~RT_IKI, scales ='free_x') +
  my_theme +
  labs(x = 'Effect (ms) of 1SD increase in predictor')

ggsave('./plots/posterior_density.png', height = 16, width = 30, units = 'cm')

# do CI plot ----

pop_effs %>% 
  group_by(RT_IKI, IV) %>% 
  summarise(Median = median(param)
            ,lo = quantile(param,.05)
            ,hi = quantile(param,.95)) %>% 
  ggplot(aes(y = IV, x = Median)) +
  geom_vline(aes(xintercept = 0),linetype = 'dashed') +
  geom_errorbarh(aes(xmin = lo, xmax = hi, height = 0)
                 ,size = 2, colour = 'dark grey') +
  geom_point(size = 4, colour = 'dark red') +
  facet_wrap(~RT_IKI, scales ='free_x') +
  my_theme +
  labs(x = 'Effect (ms) of 1SD increase in predictor, 95% CrI')

ggsave('./plots/CIs.png', height = 16, width = 30, units = 'cm')

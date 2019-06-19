library(loo)
library(tidyverse)
library(magrittr)

path <- "Frontline/stanout/"
path <- "stanout/"

(files <- dir(path, pattern = ".rda"))
ms <- gsub(files, pattern = ".rda", replacement = "")

for(i in 1:length(files)){
  m <- readRDS(paste0(path, files[i] ))
  log_lik <- extract_log_lik(m, merge_chains = F) 
  r_eff <- relative_eff(exp(log_lik)) 
  assign(paste0("loo_",ms[i]), loo(log_lik, r_eff = r_eff, cores = 2))#,  
  print(ms[i]); rm(list = "m")
}

(loos <- ls(pattern = "loo_"))
mc <- do.call(compare, lapply(loos[-c(10)], as.name))

mc %<>% as.data.frame() %>%
  round(2) %>%
  mutate(model=row.names(.)) %>%
  select(model, elpd_diff:se_elpd_loo) %>%
  remove_rownames();mc

file_out <- paste0(path, "loo_results.csv")
write_csv(mc, file_out)


# Load libraries
library(tidybayes)
library(tidyverse)
library(rstan)
library(loo)




#### Model A #####


samples_m_A_posterior <- readRDS(file="Models/samples_m_A_posterior.rds")


log_lik_A <- extract_log_lik(samples_m_A_posterior, merge_chains = FALSE)


r_eff_A <- relative_eff(exp(log_lik_A), cores = parallel::detectCores()) 

loo_A <- loo(log_lik_A, r_eff = r_eff_A, cores = parallel::detectCores())


# print(loo_A)




#### Model B #####


samples_m_B_posterior <- readRDS(file="Models/samples_m_B_posterior.rds")


log_lik_B <- extract_log_lik(samples_m_B_posterior, merge_chains = FALSE)


r_eff_B <- relative_eff(exp(log_lik_B), cores = parallel::detectCores()) 

loo_B <- loo(log_lik_B, r_eff = r_eff_B, cores = parallel::detectCores())


# print(loo_B)




#### Model C #####


samples_m_C_posterior <- readRDS(file="Models/samples_m_C_posterior.rds")


log_lik_C <- extract_log_lik(samples_m_C_posterior, merge_chains = FALSE)


r_eff_C <- relative_eff(exp(log_lik_C), cores = parallel::detectCores()) 

loo_C <- loo(log_lik_C, r_eff = r_eff_C, cores = parallel::detectCores())


# print(loo_C)




#### Model D #####


samples_m_D_posterior <- readRDS(file="Models/samples_m_D_posterior.rds")


log_lik_D <- extract_log_lik(samples_m_D_posterior, merge_chains = FALSE)


r_eff_D <- relative_eff(exp(log_lik_D), cores = parallel::detectCores()) 

loo_D <- loo(log_lik_D, r_eff = r_eff_D, cores = parallel::detectCores())


# print(loo_D)





######### Model Comparison #########



comp <- loo_compare(loo_A, loo_B, loo_C, loo_D)


print(comp) 


# As Data Frame

df_comp <- data.frame(comp)|>
  mutate(Model = c("D","C","B","A"))|>
  tibble()|>
  select(Model,
         `Expected Log Predictive Density` = elpd_loo,
         `ELPD Difference` = elpd_diff,
         `ELPD Difference Standard Error` = se_diff)


df_comp

saveRDS(df_comp, file = "Models/df_comp.rds")


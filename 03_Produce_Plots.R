
# Libraries
library(tidyverse)
library(tidybayes)

# Useful functions
f_inv_logit <-  function(x) {exp(x)/(1+exp(x))} # Inverse_logit


species <- c("Blackcap","Chiffchaff","Robin")


months <- c("October", "November", "December", "January", "February", "March", 
            "April")

winters <- 2006:2017


years <- 2007:2018





# Logistic Function Plot


seq_logistic_x <- seq(from = -6, to = 6, length.out = 100)  

seq_logistic_y <- map_dbl(seq_logistic_x, .f = function(x) f_inv_logit(x))


plt_logistic <- ggplot(data = tibble(x = seq_logistic_x,
                     y = seq_logistic_y),
       aes(x = x, y = y))+
  geom_line(alpha = 0.5, colour = "#440154",lwd=1.2)+
  geom_hline(yintercept = 0.5)+
  ylab("Output")+
  xlab("Input Values")+
  ggtitle("The Logistic Function")+
  theme_bw()


saveRDS(plt_logistic, file = "Plots/plt_logistic.rds")




# Logistic Output


seq_logistic_x <- rnorm(10000, mean = 0, sd = 10)  

seq_logistic_y <- map_dbl(seq_logistic_x, .f = function(x) f_inv_logit(x))


plt_logistic_output <- ggplot(data = tibble(x = seq_logistic_y),
                       aes(x, after_stat(density)))+
  geom_histogram(alpha = 0.5, fill = "#440154",lwd=1.2, bins = 30)+
  ylab("Probability")+
  xlab("Density")+
  ggtitle("Potential Pitfall of Priors on the Logistic Scale")+
  theme_bw()

plt_logistic_output

saveRDS(plt_logistic_output, file = "Plots/plt_logistic_output.rds")




###### Histograms for Prior Probability ######


df_prior <- expand_grid(species = species,
                        adult = 0:1,
                        month = 2,
                        winter = 2)|>
  mutate(adult = if_else(adult == 1, "Adult", "Juvenile"),
         i = 1:6)

samples_m_D_prior <- readRDS(file="Models/samples_m_D_prior.rds")


# Phi


draws_phi_prior <- gather_draws(samples_m_D_prior,
                                pars = phi_prior[i])|>
  left_join(df_prior, by = "i")


plt_phi_prior <- draws_phi_prior|>
  ggplot(aes(.value, after_stat(density)))+
  geom_histogram(alpha = 0.5, fill = "#440154", bins = 30)+
  facet_wrap(c("species", "adult"))+
  theme_bw()+
  xlab("Probability")+
  ylab("Density")+
  ggtitle("Histograms of Prior Distributions for Intra-Winter \nSurvival Probabilities")

plt_phi_prior

saveRDS(plt_phi_prior, file = "Plots/plt_phi_prior.rds")




# Psi


draws_psi_prior <- gather_draws(samples_m_D_prior,
                                pars = psi_prior[i])|>
  left_join(df_prior, by = "i")|>
  filter(adult == "Adult")


plt_psi_prior <- draws_psi_prior|>
  ggplot(aes(.value, after_stat(density)))+
  geom_histogram(alpha= 0.5, fill = "#440154", bins = 30)+
  facet_wrap(~species)+
  theme_bw()+
  xlab("Probability")+
  ylab("Density")+
  ggtitle("Histograms of Prior Distributions for Inter-Winter \nSurvival Probabilities")

plt_psi_prior

saveRDS(plt_psi_prior, file = "Plots/plt_psi_prior.rds")




# P


draws_p_prior <- gather_draws(samples_m_D_prior,
                                pars = p_prior[i])|>
  left_join(df_prior, by = "i")


plt_p_prior <- draws_p_prior|>
  ggplot(aes(.value, after_stat(density)))+
  geom_histogram(alpha = 0.5, fill = "#440154", bins = 30)+
  facet_wrap(c("species", "adult"))+
  theme_bw()+
  xlab("Probability")+
  ylab("Density")+
  ggtitle("Histograms of Prior Distributions for Capture Probabilities")

plt_p_prior

saveRDS(plt_p_prior, file = "Plots/plt_p_prior.rds")









###### Trace Plots #######


samples_m_D_posterior <- readRDS(file="Models/samples_m_D_posterior.rds")


# Species Phi

draws_spec_posterior <- gather_draws(samples_m_D_posterior,
                                     pars = a_spec_phi[i])|>
  mutate(species = species[i])


plt_trace_spec_phi <- draws_spec_posterior|>
  mutate(.chain = as.factor(.chain))|>
  rename(Chain = .chain)|>
  ggplot(aes(x = .iteration, y = .value, colour = Chain))+
  geom_line(alpha = 0.5)+
  facet_wrap(~species, scales = "free_y", nrow = 3)+
  ggtitle("Trace Plot for Species Parameters - Intra-Winter Probabilities")+
  xlab("Sample")+
  ylab("")+
  scale_color_viridis_d()+
  theme_bw()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

plt_trace_spec_phi


saveRDS(plt_trace_spec_phi, file = "Plots/plt_trace_spec_phi.rds")



# Month Phi

draws_month_posterior <- gather_draws(samples_m_D_posterior,
                                     pars = b_mon_phi[i])|>
  mutate(month = months[i])|>
  filter(month != "October")


plt_trace_month_phi <- draws_month_posterior|>
  mutate(.chain = as.factor(.chain),
         month = factor(month, levels = months))|>
  rename(Chain = .chain)|>
  ggplot(aes(x = .iteration, y = .value, colour = Chain))+
  geom_line(alpha = 0.5)+
  facet_wrap(~month, scales = "free_y", nrow = 6)+
  ggtitle("Trace Plot for Month Parameters - Intra-Winter Probabilities")+
  xlab("Sample")+
  ylab("")+
  scale_color_viridis_d()+
  theme_bw()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

plt_trace_month_phi


saveRDS(plt_trace_month_phi, file = "Plots/plt_trace_month_phi.rds")







################### Results #####################





samples_m_D_posterior <- readRDS(file="Models/samples_m_D_posterior.rds")


###### Intra-Winter Survival Results #######


#### Species ####


draws_spec_phi <- gather_draws(samples_m_D_posterior,
                                      pars = a_spec_phi[i])|>
  mutate(Species = species[i])


# Histogram
plt_spec_phi <- draws_spec_phi|>
  ggplot(aes(x = .value, after_stat(density), fill = Species))+
  geom_histogram(alpha = 0.5, bins = 30)+
  ggtitle("Histograms for Species Parameters - Intra-Winter Survival")+
  xlab("Parameter Value")+
  ylab("Density")+
  scale_fill_viridis_d()+
  theme_bw()



plt_spec_phi

saveRDS(plt_spec_phi, file = "Plots/plt_spec_phi.rds")


# Credible Intervals

plt_spec_phi_conf <- draws_spec_phi|>
  group_by(Species)|> 
  mean_qi(.value, .width=0.95)|> 
  ggplot(aes(y = Species, x = .value, 
             xmin = .lower, xmax = .upper)) +
  geom_pointinterval(alpha = 0.5, colour = "#440154")+
  ggtitle("Credible Intervals for Species Parameters \n- Intra-Winter Survival")+
  ylab("Species")+
  xlab("Parameter Value")+
  scale_colour_viridis_d()+
  theme_bw()+
  scale_y_discrete(limits=rev)


plt_spec_phi_conf


saveRDS(plt_spec_phi_conf, file = "Plots/plt_spec_phi_conf.rds")



#### Age ####



draws_age_phi <- gather_draws(samples_m_D_posterior,
                               pars = b_adult_phi)


# Histogram
plt_age_phi <- draws_age_phi|>
  ggplot(aes(x = .value, after_stat(density)))+
  geom_histogram(alpha = 0.5, bins = 30, fill = "#440154")+
  ggtitle("Histogram for Adult Parameter - Intra-Winter Survival")+
  xlab("Parameter Value")+
  ylab("Density")+
  scale_fill_viridis_d()+
  theme_bw()



plt_age_phi


saveRDS(plt_age_phi, file = "Plots/plt_age_phi.rds")


# Credible Intervals

plt_age_phi_conf <- draws_age_phi|>
  mean_qi(.value, .width=0.95)|> 
  ggplot(aes(x = .value, 
             xmin = .lower, xmax = .upper)) +
  geom_pointinterval(alpha = 0.5, colour = "#440154")+
  geom_vline(xintercept = 0, lty = 2)+
  ggtitle("Credible Interval for Adult Parameter \n- Intra-Winter Survival")+
  ylab("")+
  xlab("Parameter Value")+
  scale_colour_viridis_d()+
  theme_bw()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())


plt_age_phi_conf


saveRDS(plt_age_phi_conf, file = "Plots/plt_age_phi_conf.rds")


# Example probabilities


a_spec_phi <- spread_draws(samples_m_D_posterior,
                              pars = a_spec_phi[i])|>
  filter(i == 3)|> # robin
  pull(a_spec_phi)


b_adult_phi <- spread_draws(samples_m_D_posterior,
                           pars = b_adult_phi)|>
  pull(b_adult_phi)


b_mon_phi <- spread_draws(samples_m_D_posterior,
                           pars = b_mon_phi[i])|>
  filter(i == 4)|> # January
  pull(b_mon_phi)


b_wint_phi <- spread_draws(samples_m_D_posterior,
                          pars = b_wint_phi[i])|>
  filter(i == 3)|> # 2009
  pull(b_wint_phi)


draws_age_phi_prob = tibble(draws = 1:4000,
                            a_spec_phi = a_spec_phi,
                            b_adult_phi = b_adult_phi,
                            b_mon_phi = b_mon_phi,
                            b_wint_phi = b_wint_phi)|>
  mutate(Juvenile = f_inv_logit(a_spec_phi + b_mon_phi + b_wint_phi),
         Adult = f_inv_logit(a_spec_phi + b_adult_phi + b_mon_phi + b_wint_phi))|>
  select(draws, Juvenile, Adult)|>
  pivot_longer(c(Juvenile, Adult), names_to = "Age", values_to = "Probability")



# Credible Intervals

plt_age_phi_prob <- draws_age_phi_prob|>
  group_by(Age)|> 
  mean_qi(Probability, .width=0.95)|> 
  ggplot(aes(y = Age, x = Probability, 
             xmin = .lower, xmax = .upper)) +
  geom_pointinterval(alpha = 0.5, colour = "#440154")+
  ggtitle("Credible Intervals for Robin Intra-Winter \nSurvival Probabilities")+
  ylab("Age")+
  xlab("Probability of Survival")+
  scale_colour_viridis_d()+
  theme_bw()+
  xlim(c(0,1))+
  scale_y_discrete(limits=rev)


plt_age_phi_prob


saveRDS(plt_age_phi_prob, file = "Plots/plt_age_phi_prob.rds")



#### Month ####



draws_mon_phi <- gather_draws(samples_m_D_posterior,
                               pars = b_mon_phi[i])|>
  mutate(Month = months[i])|>
  filter(Month != "October")


# Credible Intervals

plt_mon_phi_conf <- draws_mon_phi|>
  group_by(Month)|> 
  mean_qi(.value, .width=0.95)|> 
  mutate(Month = factor(Month, levels = months))|>
  ggplot(aes(y = Month, x = .value, 
             xmin = .lower, xmax = .upper)) +
  geom_pointinterval(alpha = 0.5, colour = "#440154")+
  ggtitle("Month Parameters \n- Intra-Winter Survival")+
  ylab("Month")+
  xlab("Parameter Value")+
  scale_colour_viridis_d()+
  theme_bw()+
  scale_y_discrete(limits=rev)


plt_mon_phi_conf


saveRDS(plt_mon_phi_conf, file = "Plots/plt_mon_phi_conf.rds")




#### Winter ####



draws_wint_phi <- gather_draws(samples_m_D_posterior,
                              pars = b_wint_phi[i])|>
  mutate(Winter = winters[i])


# Credible Intervals

plt_wint_phi_conf <- draws_wint_phi|>
  group_by(Winter)|> 
  mean_qi(.value, .width=0.95)|> 
  ggplot(aes(y = Winter, x = .value, 
             xmin = .lower, xmax = .upper)) +
  geom_pointinterval(alpha = 0.5, colour = "#440154")+
  ggtitle("Winter Parameters \n- Intra-Winter Survival")+
  ylab("Winter")+
  xlab("Parameter Value")+
  scale_colour_viridis_d()+
  theme_bw()+
  scale_y_reverse(breaks = 2006:2017)


plt_wint_phi_conf


saveRDS(plt_wint_phi_conf, file = "Plots/plt_wint_phi_conf.rds")






###### Inter-Winter Parameters Results #######



#### Species ####


draws_spec_psi <- gather_draws(samples_m_D_posterior,
                               pars = a_spec_psi[i])|>
  mutate(Species = species[i])


# Histogram
plt_spec_psi <- draws_spec_psi|>
  ggplot(aes(x = .value, after_stat(density), fill = Species))+
  geom_histogram(alpha = 0.5, bins = 30)+
  ggtitle("Histograms for Species Parameters - Inter-Winter Survival")+
  xlab("Parameter Value")+
  ylab("Density")+
  scale_fill_viridis_d()+
  theme_bw()



plt_spec_psi

saveRDS(plt_spec_psi, file = "Plots/plt_spec_psi.rds")


# Credible Intervals

plt_spec_psi_conf <- draws_spec_psi|>
  group_by(Species)|> 
  mean_qi(.value, .width=0.95)|> 
  ggplot(aes(y = Species, x = .value, 
             xmin = .lower, xmax = .upper)) +
  geom_pointinterval(alpha = 0.5, colour = "#440154")+
  ggtitle("Credible Intervals for Species Parameters \n- Inter-Winter Survival")+
  ylab("Species")+
  xlab("Parameter Value")+
  scale_colour_viridis_d()+
  theme_bw()+
  scale_y_discrete(limits=rev)


plt_spec_psi_conf


saveRDS(plt_spec_psi_conf, file = "Plots/plt_spec_psi_conf.rds")





#### Year ####




draws_year_psi <- gather_draws(samples_m_D_posterior,
                               pars = b_wint_psi[i])|>
  mutate(Year = years[i])


# Credible Intervals

plt_year_psi_conf <- draws_year_psi|>
  group_by(Year)|> 
  mean_qi(.value, .width=0.95)|> 
  ggplot(aes(y = Year, x = .value, 
             xmin = .lower, xmax = .upper)) +
  geom_pointinterval(alpha = 0.5, colour = "#440154")+
  ggtitle("Credible Intervals for Year Parameters \n- Inter-Winter Survival")+
  ylab("Year")+
  xlab("Parameter Value")+
  scale_colour_viridis_d()+
  theme_bw()+
  scale_y_reverse(breaks = 2007:2018)


plt_year_psi_conf


saveRDS(plt_year_psi_conf, file = "Plots/plt_year_psi_conf.rds")










###### Capture Probabilities Results #######




#### Species ####


draws_spec_p <- gather_draws(samples_m_D_posterior,
                               pars = a_spec_p[i])|>
  mutate(Species = species[i])


# Histogram
plt_spec_p <- draws_spec_p|>
  ggplot(aes(x = .value, after_stat(density), fill = Species))+
  geom_histogram(alpha = 0.5, bins = 30)+
  ggtitle("Histograms for Species Capture Parameters")+
  xlab("Parameter Value")+
  ylab("Density")+
  scale_fill_viridis_d()+
  theme_bw()



plt_spec_p

saveRDS(plt_spec_p, file = "Plots/plt_spec_p.rds")


# Credible Intervals

plt_spec_p_conf <- draws_spec_p|>
  group_by(Species)|> 
  mean_qi(.value, .width=0.95)|> 
  ggplot(aes(y = Species, x = .value, 
             xmin = .lower, xmax = .upper)) +
  geom_pointinterval(alpha = 0.5, colour = "#440154")+
  ggtitle("Credible Intervals for Species \nCapture Parameters")+
  ylab("Species")+
  xlab("Parameter Value")+
  scale_colour_viridis_d()+
  theme_bw()+
  scale_y_discrete(limits=rev)


plt_spec_p_conf


saveRDS(plt_spec_p_conf, file = "Plots/plt_spec_p_conf.rds")



#### Age ####



draws_age_p <- gather_draws(samples_m_D_posterior,
                              pars = b_adult_p)


# Histogram
plt_age_p <- draws_age_p|>
  ggplot(aes(x = .value, after_stat(density)))+
  geom_histogram(alpha = 0.5, bins = 30, fill = "#440154")+
  ggtitle("Histogram for Adult Capture Parameter")+
  xlab("Parameter Value")+
  ylab("Density")+
  scale_fill_viridis_d()+
  theme_bw()



plt_age_p


saveRDS(plt_age_p, file = "Plots/plt_age_p.rds")


# Credible Intervals

plt_age_p_conf <- draws_age_p|>
  mean_qi(.value, .width=0.95)|> 
  ggplot(aes(x = .value, 
             xmin = .lower, xmax = .upper)) +
  geom_pointinterval(alpha = 0.5, colour = "#440154")+
  geom_vline(xintercept = 0, lty = 2)+
  ggtitle("Credible Interval for Adult Parameter \nCapture Parameter")+
  ylab("")+
  xlab("Parameter Value")+
  scale_colour_viridis_d()+
  xlim(c(-2,2))+
  theme_bw()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())


plt_age_p_conf


saveRDS(plt_age_p_conf, file = "Plots/plt_age_p_conf.rds")


# Example probabilities


b_mon_p <- spread_draws(samples_m_D_posterior,
                        pars = b_mon_p[i])|>
  filter(i == 4)|> # January
  pull(b_mon_p)


b_wint_p <- spread_draws(samples_m_D_posterior,
                         pars = b_wint_p[i])|>
  filter(i == 3)|> # 2009
  pull(b_wint_p)

example_prob_p <- spread_draws(samples_m_D_posterior,
                           pars = a_spec_p[i],
                         b_adult_p)|>
  ungroup()|>
  mutate(Species = species[i],
         b_mon_p = rep(b_mon_p,3),
         b_wint_p = rep(b_wint_p,3))


draws_age_p_prob <-  example_prob_p|>
  mutate(Juvenile = f_inv_logit(a_spec_p + b_mon_p + b_wint_p),
         Adult = f_inv_logit(a_spec_p + b_adult_p + b_mon_p + b_wint_p))|>
  select(.draw, Species, Juvenile, Adult)|>
  pivot_longer(c(Juvenile, Adult), names_to = "Age", values_to = "Probability")



# Credible Intervals

plt_age_p_prob <- draws_age_p_prob|>
  group_by(Age, Species)|> 
  mean_qi(Probability, .width=0.95)|> 
  ggplot(aes(y = Age, x = Probability, 
             xmin = .lower, xmax = .upper)) +
  geom_pointinterval(alpha = 0.5, colour = "#440154")+
  facet_wrap(~Species)+
  ggtitle("Capture Probabilities by Species and Age - January 2009")+
  ylab("Age")+
  xlab("Probability of Survival")+
  scale_colour_viridis_d()+
  xlim(0,0.05)+
  theme_bw()+
  scale_y_discrete(limits=rev)+ 
  theme(axis.text.x = element_text(angle = 30))


plt_age_p_prob


saveRDS(plt_age_p_prob, file = "Plots/plt_age_p_prob.rds")



#### Month ####



draws_mon_p <- gather_draws(samples_m_D_posterior,
                              pars = b_mon_p[i])|>
  mutate(Month = months[i])|>
  filter(Month != "October")


# Credible Intervals

plt_mon_p_conf <- draws_mon_p|>
  group_by(Month)|> 
  mean_qi(.value, .width=0.95)|> 
  mutate(Month = factor(Month, levels = months))|>
  ggplot(aes(y = Month, x = .value, 
             xmin = .lower, xmax = .upper)) +
  geom_pointinterval(alpha = 0.5, colour = "#440154")+
  ggtitle("Credible Intervals for Month \nCapture Parameters")+
  ylab("Month")+
  xlab("Parameter Value")+
  scale_colour_viridis_d()+
  theme_bw()+
  scale_y_discrete(limits=rev)


plt_mon_p_conf


saveRDS(plt_mon_p_conf, file = "Plots/plt_mon_p_conf.rds")




#### Winter ####



draws_wint_p <- gather_draws(samples_m_D_posterior,
                               pars = b_wint_p[i])|>
  mutate(Winter = winters[i])


# Credible Intervals

plt_wint_p_conf <- draws_wint_p|>
  group_by(Winter)|> 
  mean_qi(.value, .width=0.95)|> 
  ggplot(aes(y = Winter, x = .value, 
             xmin = .lower, xmax = .upper)) +
  geom_pointinterval(alpha = 0.5, colour = "#440154")+
  ggtitle("Credible Intervals for Winter \nCapture Parameters")+
  ylab("Winter")+
  xlab("Parameter Value")+
  scale_colour_viridis_d()+
  theme_bw()+
  scale_y_reverse(breaks = 2006:2017)


plt_wint_p_conf


saveRDS(plt_wint_p_conf, file = "Plots/plt_wint_p_conf.rds")









###### Missing Data Results #######


missing_sample <- sample(1:145, 6)


draws_miss_p <- gather_draws(samples_m_D_posterior,
                             pars = p_adult[i])|>
  filter(i %in% missing_sample)


# Histogram
plt_miss_p <- draws_miss_p|>
  ggplot(aes(x = .value, after_stat(density)))+
  geom_histogram(alpha = 0.5, fill = "#440154", bins = 30)+
  facet_wrap(~i)+
  ggtitle("Histograms for a Sample of \nMissing Adult Probabilities")+
  xlab("Probability")+
  ylab("Density")+
  scale_fill_viridis_d()+
  theme_bw()+ 
  theme(axis.text.x = element_text(angle = 90))


plt_miss_p


saveRDS(plt_miss_p, file = "Plots/plt_miss_p.rds")



#### Random Effects ####



draws_miss_re <- spread_draws(samples_m_D_posterior,
                             pars = alpha_miss, beta_miss)|>
  mutate(Probability = rbeta(4000, alpha_miss, beta_miss))



# Histogram
plt_miss_re <- draws_miss_re|>
  ggplot(aes(x = Probability, after_stat(density)))+
  geom_histogram(alpha = 0.5, fill = "#440154", bins = 30)+
  ggtitle("Probability Implied by \nRandom Effects Parameters")+
  xlab("Probability")+
  ylab("Density")+
  scale_fill_viridis_d()+
  theme_bw()


plt_miss_re


saveRDS(plt_miss_re, file = "Plots/plt_miss_re.rds")



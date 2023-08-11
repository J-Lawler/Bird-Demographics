
########## Contents ##########

# Preamble

# Data Preparation

# Annual Model - Loops

## Simulation

## Model Code


# Annual Model - Data Frame Format

## Simulation

## Model Code


########## Preamble ##########

# Load libraries
library(tidybayes)
library(tidyverse)
library(rstan)
library(arrow)


species <- c("Blackcap","Chiffchaff","Robin")

months <- c("October", "November", "December", "January", "February", "March", 
            "April")

# Useful functions
f_inv_logit <-  function(x) {exp(x)/(1+exp(x))} # Inverse_logit

# Load in data
df_model <- read_csv("Data_Clean/df_model.csv")





########## Data Preparation ##########

# f(i) - the first observation of individual i
# l(i) - the last observation of individual i
# T - the final capture event i.e. April 2019

f_prep_for_stan <- function(df){
  
  # Species list
  species <- c("Blackcap","Chiffchaff","Robin")
  
  # Label Months - October 1, April 7
  df <- df|>
    mutate(Month = factor(Month, levels = months),
           month_label = as.integer(Month))
  
  # Label Winter
  df <- df|>
    mutate(winter_label = Winter-min(df$Winter)+1)
  
  # Label ages
  df <- df|>
    mutate(adult = if_else(Age==4,1,Age-1)) # all juveniles set to 0, adults to 1, unknowns to 2
  
  # Create vector of number of total observations for each individual
  total_obs <- df|>
    group_by(Individual)|>
    summarise(Count = sum(Present))|>
    pull(Count)
  
  # Add cumulative column to data frame
  df <- df|>
    group_by(Individual)|>
    mutate(Cumulative = cumsum(Present))|>
    ungroup()
  
  # Data frame of f(i) + 1 to l(i) for each individual
  df_obs <-  df|>
    filter(Cumulative != 0, # remove rows before f(i)
           !(Cumulative == 1 & Present == 1), # remove f(i)
           !(Cumulative == total_obs[Individual] & Present ==0)) # remove after l(i)
  
  # Data frame of l(i) to T
  df_chi <- df|>
    filter(Cumulative == total_obs[Individual]) # remove before l(i)
  
  ## Count species and individuals
  n_species <- length(unique(df$Species))
  n_indiv <- length(unique(df$Individual))
  
  # Mark final event T
  T_locations <- if_else(df_chi$Year == 2018 & df_chi$Month == "April",
                         1,0)
  
  # Mark last capture event for each individual
  l_locations <- which(df_chi$Present==1)
  
  # Mark Octobers
  october <- if_else(df_obs$Month == "October",1,0)
  october_chi <- if_else(df_chi$Month == "October",1,0)
  
  # Convert species to integers
  species_n <- map_dbl(.x = df_obs$Species, 
                       .f = function(x) which(x == species))
  species_n_chi <- map_dbl(.x = df_chi$Species, 
                       .f = function(x) which(x == species))
  
  ### Missingness ###
  
  # Labels of individuals with missing age
  indiv_miss_labels <- df|>filter(adult==2)|>pull(Individual)|>unique()
  
  # Number of individuals missing
  n_missing <- length(indiv_miss_labels)
  
  # For each row with missing data, which individual is concerned?
  # Vectors length n or n_chi: 0 if not missing,
  # the appropriate position in indiv_miss_labels if missing

  ind_missing <- rep(0,nrow(df_obs))
  
  for(i in 1:nrow(df_obs)){
    
    if(df_obs$adult[i]==2){
      
      ind_missing[i] = which(df_obs$Individual[i] == indiv_miss_labels)
      
    }
    
  }
  
  ind_missing_chi <- rep(0,nrow(df_chi)) # same process for df_chi
  
  for(i in 1:nrow(df_chi)){
    
    if(df_chi$adult[i]==2){
      
      ind_missing_chi[i] = which(df_chi$Individual[i] == indiv_miss_labels)
      
    }
    
  }
  
  df_obs$october <- october
  df_chi$october_chi <- october_chi
  df_obs$ind_missing <- ind_missing
  df_chi$ind_missing_chi <- ind_missing_chi
  
  
  ### Priors ###
  
  df_prior <- expand_grid(species = 1:3,
                          adult = 0:1,
                          month = 2,
                          winter = 2)
  
  ### Create list ###
  stan_list <- list(n_species = n_species,
                   n_indiv = n_indiv,
                   n_months = length(unique(df$Month)),
                   n_winters = length(unique(df$Winter)),
                   #n_events = length(unique(df$Event)),
                   n_missing = n_missing,
                   n_prior = nrow(df_prior),
                   prior = 0,
                   posterior = 0,
                   loo_calcs = 0,
                   n = nrow(df_obs),
                   indiv = df_obs$Individual,
                   species = species_n,
                   adult = df_obs$adult,
                   month = df_obs$month_label,
                   winter = df_obs$winter_label,
                   #event = df_obs$Event,
                   present = df_obs$Present,
                   october = df_obs$october,
                   ind_missing = df_obs$ind_missing,
                   n_chi = nrow(df_chi),
                   #indiv_chi = df_chi$Individual,
                   species_chi = species_n_chi,
                   adult_chi = df_chi$adult,
                   month_chi = df_chi$month_label,
                   winter_chi = df_chi$winter_label,
                   #event_chi = df_chi$Event,
                   present_chi = df_chi$Present,
                   october_chi = df_chi$october_chi,
                   ind_missing_chi = df_chi$ind_missing_chi,
                   T = T_locations,
                   l = l_locations,
                   species_prior = df_prior$species,
                   adult_prior = df_prior$adult,
                   month_prior = df_prior$month,
                   winter_prior = df_prior$winter)
  
  # Return objects
  output <- list("df_obs" = df_obs,
                 "df_chi" = df_chi,
                 "df_prior" = df_prior,
                 "stan_list" = stan_list,
                 "indiv_miss_labels" = indiv_miss_labels)
  
  return(output)
}

 


########## Model A - Data Frame Format ##########



##### Model Code #####


code_m_A <- 
  "data{
  
  //
  int<lower=0> n_species;
  int<lower=0> n_indiv;
  int<lower=0> n_months;
  int<lower=0> n_winters;
  int<lower=0> n_missing;
  int<lower=0> n_prior;
  
  // Switches
  int<lower=0, upper = 1> prior;
  int<lower=0, upper = 1> posterior;
  int<lower=0, upper = 1> loo_calcs;
  
  // Observed Calculations, f(i) + 1 to l(i)
  int<lower=0> n;
  array[n] int<lower=1, upper = n_species> species;
  vector[n] adult;
  array[n] int<lower=1, upper = n_months> month;
  array[n] int<lower=1, upper = n_winters> winter;
  array[n] int<lower=1, upper = n_indiv> indiv;
  vector[n] present;
  vector[n] october; // Indicator for October
  array[n] int<lower=0, upper = n_missing> ind_missing;
  
  // Chi Calculations, l(i) to T
  int<lower=0> n_chi;
  array[n_chi] int<lower=1, upper = n_species> species_chi;
  vector[n_chi] adult_chi;
  array[n_chi] int<lower=1, upper = n_months> month_chi;
  array[n_chi] int<lower=1, upper = n_winters> winter_chi;
  vector[n_chi] present_chi;
  vector[n_chi] october_chi; // Indicator for October
  array[n_chi] int T; // 1 if time T, 0 otherwise
  array[n_chi] int<lower=0, upper = n_missing> ind_missing_chi;
  array[n_indiv] int l; // identifies row of final observation

  // Data for Generating Prior Distributions of Phi, Psi, P
  array[n_prior] int<lower=1, upper = n_species> species_prior;
  vector[n_prior] adult_prior;
  array[n_prior] int<lower=1, upper = n_months> month_prior;
  array[n_prior] int<lower=1, upper = n_winters> winter_prior;
  
}
parameters{
  
  // Phi - Winter monthly survival
  real a_phi; // intercept

  
  // Psi - Summer survival 
  real a_psi; // intercept
  
  
  // P - Capture probabilities
  real a_p; // intercept

  
  // Probability of adult for missing individuals
  //vector<lower = 0, upper = 1>[n_missing] p_adult;
  //real<lower = 0> alpha_miss;
  //real<lower = 0> beta_miss;
  
 
}
model{

  // Priors
  
  // Phi
  a_phi ~ normal(0,1.5);
  
  
  // Psi
  a_psi ~ normal(0,1.5);
  
  
  // P
  a_p ~ normal(0,1.5);

  
  // Missingness Priors
  //p_adult ~ beta(alpha_miss,beta_miss);
  
  
  // Missingness RE Priors
  
  //alpha_miss ~ normal(1,0.2);
  //beta_miss ~ normal(1,0.2);
  
  
  
  // Posterior
  if(posterior == 1){
    
    
  // Linear predictors
  
  vector[n] lin_phi = a_phi * ones_vector(n);
                      
  vector[n] lin_psi = a_psi * ones_vector(n);
  
  vector[n] lin_pi = october.*lin_psi + (1-october).*lin_phi;
  
  vector[n] lin_p = a_p * ones_vector(n);
  
  // Probabilities from predictors
  
  vector[n] pi = inv_logit(lin_pi);
  vector[n] p = inv_logit(lin_p);
  

  // Likelihood
  
  target += log(pi) + present.*log(p) + (1-present).*log1m(p);
  
  
  // Chi Parameter
  
  
  vector[n_chi] lin_phi_chi = a_phi * ones_vector(n_chi);
                      
  vector[n_chi] lin_psi_chi = a_psi * ones_vector(n_chi);
  
  vector[n_chi] lin_pi_chi = october_chi.*lin_psi_chi 
                              + (1-october_chi).*lin_phi_chi;
  
  vector[n_chi] lin_p_chi = a_p * ones_vector(n_chi);
  
  // Probabilities from Predictors
  vector[n_chi] pi_chi = inv_logit(lin_pi_chi);
  vector[n_chi] p_chi = inv_logit(lin_p_chi);
  
  
  vector[n_chi] chi;
  
  for (t in 1:n_chi){
    
    int s = n_chi - t + 1; // working backwards from end of data frame
    
    if (T[s] != 1) { 
    
      chi[s] = (1 - pi_chi[s+1]) 
                          + pi_chi[s+1]*(1-p_chi[s+1])*chi[s+1];
    
    } else { 
    
      chi[s] = 1;
    
    }  
  
  
  }
    
  target += log(chi[l]);
  
  } // end of posterior conditional section 
  
}
generated quantities{

// Probabilities for prior predictive simulation

real phi_prior;
real psi_prior;
real p_prior;

if(prior == 1){

  // Phi
  
  phi_prior = inv_logit(a_phi);
                              
  // Psi

  psi_prior = inv_logit(a_psi);
  
  // P

  p_prior = inv_logit(a_p);

}


// Log Likelihood for Model Comparison


vector[n_indiv] log_lik = zeros_vector(n_indiv);

if(loo_calcs == 1){
 
 
  // Linear predictors
  
  vector[n] lin_phi = a_phi * ones_vector(n);
                      
  vector[n] lin_psi = a_psi * ones_vector(n);
  
  vector[n] lin_pi = october.*lin_psi + (1-october).*lin_phi;
  
  vector[n] lin_p = a_p * ones_vector(n);
  
  // Probabilities from predictors
  
  vector[n] pi = inv_logit(lin_pi);
  vector[n] p = inv_logit(lin_p);
  

  // Likelihood
  
  vector[n] log_lik_n = log(pi) + present.*log(p) + (1-present).*log1m(p);
  
  for(i in 1:n){
  
  log_lik[indiv[i]] = log_lik[indiv[i]] + log_lik_n[i];
  
  }
  
  // Chi Parameter
  
 
  vector[n_chi] lin_phi_chi = a_phi * ones_vector(n_chi);
                      
  vector[n_chi] lin_psi_chi = a_psi * ones_vector(n_chi);
  
  vector[n_chi] lin_pi_chi = october_chi.*lin_psi_chi 
                              + (1-october_chi).*lin_phi_chi;
  
  vector[n_chi] lin_p_chi = a_p * ones_vector(n_chi);
  
  // Probabilities from Predictors
  vector[n_chi] pi_chi = inv_logit(lin_pi_chi);
  vector[n_chi] p_chi = inv_logit(lin_p_chi);
  
  // Missing assigned adult
  
  vector[n_chi] chi;
  
  for (t in 1:n_chi){
    
    int s = n_chi - t + 1; // working backwards from end of data frame
    
    if (T[s] != 1) { 
    
      chi[s] = (1 - pi_chi[s+1]) 
                          + pi_chi[s+1]*(1-p_chi[s+1])*chi[s+1];
    
    } else { 
    
      chi[s] = 1;
    
    }  
  
  
  }
    
  log_lik = log_lik + log(chi[l]);
  
}


}"


# Create list for model
outputs_A <- f_prep_for_stan(df_model)




##### Prior Predictive #####


prior_list <- outputs_A$stan_list

# Permit generation of prior predictive
prior_list$prior <- 1


if(!file.exists("Models/samples_m_A_prior.rds")){
  
  # Compiles the model
  stan_m_A <- stan_model(model_name = "stan_m_A",model_code=code_m_A)
  
  # Sample
  samples_m_A_prior <-  sampling(stan_m_A,
                                 data = prior_list, 
                                 chains=1, iter = 10000,
                                 cores = parallel::detectCores(),
                                 thin = 10,
                                 control = list(adapt_delta = 0.95),
                                 pars = "log_lik",
                                 include = FALSE)
  
  # Save file
  saveRDS(samples_m_A_prior,file="Models/samples_m_A_prior.rds")
}else{
  samples_m_A_prior <- readRDS(file="Models/samples_m_A_prior.rds")
}




##### Plot Prior Predictive #####


## Draws

# Phi
gather_draws(samples_m_A_prior, 
             phi_prior)|>
  ggplot(aes(x = .value))+
  geom_histogram(fill="#440154", bins = 30)+
  ggtitle("Prior Predictive Distribution - Phi Probabilities")+
  ylab("Count")+
  xlab("Value")+
  theme_bw()


# Psi
gather_draws(samples_m_A_prior, 
             psi_prior)|>
  ggplot(aes(x = .value))+
  geom_histogram(fill="#440154", bins = 30)+
  ggtitle("Prior Predictive Distribution - Psi Probabilities")+
  ylab("Count")+
  xlab("Value")+
  theme_bw()


# P
gather_draws(samples_m_A_prior, 
             p_prior)|>
  ggplot(aes(x = .value))+
  geom_histogram(fill="#440154", bins = 30)+
  ggtitle("Prior Predictive Distribution - Capture Probabilities")+
  ylab("Count")+
  xlab("Value")+
  theme_bw()


# Species
gather_draws(samples_m_A_prior, 
             a_phi)|>
  ggplot(aes(x = .value))+
  geom_histogram(fill="#440154", bins = 30)+
  ggtitle("Prior Predictive Distribution - Alpha Parameter")+
  ylab("Count")+
  xlab("Value")+
  theme_bw()

##### Posterior Predictive #####




posterior_list <- outputs_A$stan_list

# Permit generation of posterior
posterior_list$posterior <- 1
posterior_list$loo_calcs <- 1


if(!file.exists("Models/samples_m_A_posterior.rds")){
  
  # Compiles the model
  stan_m_A <- stan_model(model_name = "stan_m_A",model_code=code_m_A)
  
  start_time <- Sys.time()
  # Sample
  samples_m_A_posterior <-  sampling(stan_m_A,
                                 data = posterior_list, 
                                 chains=4, iter = 2000,
                                 cores = parallel::detectCores(),
                                 #thin = 10,
                                 control = list(adapt_delta = 0.95),
                                 pars = c("phi_prior", "psi_prior", "p_prior"),
                                 include = FALSE)
  
  end_time <- Sys.time()
  
  end_time-start_time
  
  # Save file
  saveRDS(samples_m_A_posterior,file="Models/samples_m_A_posterior.rds")
}else{
  samples_m_A_posterior <- readRDS(file="Models/samples_m_A_posterior.rds")
}






##### Plot Posterior Predictive #####


## Draws

# Phi
gather_draws(samples_m_A_posterior, 
             a_phi)|>
  mutate(phi = f_inv_logit(.value))|>
  ggplot(aes(x = phi))+
  geom_histogram(fill="#440154", bins = 30)+
  ggtitle("Prior Predictive Distribution - Phi")+
  ylab("Count")+
  xlab("Value")+
  theme_bw()


# Phi
gather_draws(samples_m_A_posterior, 
             a_psi)|>
  mutate(psi = f_inv_logit(.value))|>
  ggplot(aes(x = psi))+
  geom_histogram(fill="#440154", bins = 30)+
  ggtitle("Prior Predictive Distribution - Psi")+
  ylab("Count")+
  xlab("Value")+
  theme_bw()


# P
gather_draws(samples_m_A_posterior, 
             a_p)|>
  mutate(p = f_inv_logit(.value))|>
  ggplot(aes(x = p))+
  geom_histogram(fill="#440154", bins = 30)+
  ggtitle("Prior Predictive Distribution - P")+
  ylab("Count")+
  xlab("Value")+
  theme_bw()






#######################################



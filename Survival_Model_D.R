
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
  T_locations <- if_else(df_chi$Year == 2019 & df_chi$Month == "April",
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
  
  ### Create list ###
  stan_list <- list(n_species = n_species,
                   n_indiv = n_indiv,
                   n_months = length(unique(df$Month)),
                   n_winters = length(unique(df$Winter)),
                   n_events = length(unique(df$Event)),
                   n_missing = n_missing,
                   n = nrow(df_obs),
                   #indiv = df_obs$Individual,
                   species = species_n,
                   adult = df_obs$adult,
                   month = df_obs$month_label,
                   winter = df_obs$winter_label,
                   event = df_obs$Event,
                   present = df_obs$Present,
                   october = df_obs$october,
                   ind_missing = df_obs$ind_missing,
                   n_chi = nrow(df_chi),
                   #indiv_chi = df_chi$Individual,
                   species_chi = species_n_chi,
                   adult_chi = df_chi$adult,
                   month_chi = df_chi$month_label,
                   winter_chi = df_chi$winter_label,
                   event_chi = df_chi$Event,
                   present_chi = df_chi$Present,
                   october_chi = df_chi$october_chi,
                   ind_missing_chi = df_chi$ind_missing_chi,
                   T = T_locations,
                   l = l_locations)
  
  # Return objects
  output <- list("df_obs" = df_obs,
                 "df_chi" = df_chi, 
                 "stan_list" = stan_list,
                 "indiv_miss_labels" = indiv_miss_labels)
  
  return(output)
}

 


########## Annual Model - Data Frame Format ##########



##### Model Code #####


code_m_D <- 
  "data{
  
  //
  int<lower=0> n_species;
  int<lower=0> n_indiv;
  int<lower=0> n_months;
  int<lower=0> n_winters;
  int<lower=0> n_events;
  int<lower=0> n_missing;
  
  // Observed Calculations, f(i) + 1 to l(i)
  int<lower=0> n;
  array[n] int<lower=1, upper = n_species> species;
  vector[n] adult;
  array[n] int<lower=1, upper = n_months> month;
  array[n] int<lower=1, upper = n_winters> winter;
  array[n] int<lower=1, upper = n_events> event;
  vector[n] present;
  vector[n] october; // Indicator for October
  array[n] int<lower=0, upper = n_missing> ind_missing;
  
  // Chi Calculations, l(i) to T
  int<lower=0> n_chi;
  array[n_chi] int<lower=1, upper = n_species> species_chi;
  vector[n_chi] adult_chi;
  array[n_chi] int<lower=1, upper = n_months> month_chi;
  array[n_chi] int<lower=1, upper = n_winters> winter_chi;
  array[n_chi] int<lower=1, upper = n_events> event_chi;
  vector[n_chi] present_chi;
  vector[n_chi] october_chi; // Indicator for October
  array[n_chi] int T; // 1 if time T, 0 otherwise
  array[n_chi] int<lower=0, upper = n_missing> ind_missing_chi;
  array[n_indiv] int l; // identifies row of final observation

}
parameters{
  
  // Phi - Winter monthly survival
  vector[n_species] a_spec_phi; // intercept
  real b_adult_phi; // adult parameter
  vector[n_months] b_mon_phi; // month parameter
  vector[n_winters] b_wint_phi; // winter year parameter

  
  // Psi - Summer survival 
  vector[n_species] a_spec_psi; // intercept
  vector[n_winters] b_wint_psi; // winter year parameter
  
  
  // P - Capture probabilities
  vector[n_species] a_spec_p; // intercept
  real b_adult_p; // adult parameter
  vector[n_events] b_event_p; // event parameter
  
  
  // Probability of adult for missing individuals
  vector<lower = 0, upper = 1>[n_missing] p_adult;
  real<lower = 0> alpha_miss;
  real<lower = 0> beta_miss;
  
  // Random Effects Parameters
  
  // Phi
  real mu_spec_phi;
  real<lower=0> sig_spec_phi;
  
  real mu_mon_phi;
  real<lower=0> sig_mon_phi;
  
  real mu_wint_phi;
  real<lower=0> sig_wint_phi;
  
  
  // Psi
  
  real mu_spec_psi;
  real<lower=0> sig_spec_psi;
  
  real mu_wint_psi;
  real<lower=0> sig_wint_psi;
  
  // P
  
  real mu_spec_p;
  real<lower=0> sig_spec_p;
  
  real mu_event_p;
  real<lower=0> sig_event_p;
    
}
model{

  // Priors
  
  // Phi
  a_spec_phi ~ normal(mu_spec_phi,sig_spec_phi);
  b_adult_phi ~ normal(0,0.5);
  b_mon_phi ~ normal(mu_mon_phi,sig_mon_phi);
  b_wint_phi ~ normal(mu_wint_phi,sig_wint_phi);
  
  
  // Psi
  a_spec_psi ~ normal(mu_spec_psi,sig_spec_psi);
  b_wint_psi ~ normal(mu_wint_psi,sig_wint_psi);
  
  
  // P
  a_spec_p ~ normal(mu_spec_p,sig_spec_p);
  b_adult_p ~ normal(0,0.5);
  b_event_p ~ normal(mu_event_p,sig_event_p);
  
  // Missingness Priors
  p_adult ~ beta(alpha_miss,beta_miss);
  
  // Random Effects Priors
  
  // Phi
  mu_spec_phi ~ normal(0,1.0);
  sig_spec_phi ~ exponential(1);
  
  mu_mon_phi ~ normal(0,0.5);
  sig_mon_phi ~ exponential(1);
  
  mu_wint_phi ~ normal(0,0.5);
  sig_wint_phi ~ exponential(1);
  
  
  // Psi
  
  mu_spec_psi ~ normal(0,1.0);
  sig_spec_psi ~ exponential(1);
  
  mu_wint_psi ~ normal(0,0.5);
  sig_wint_psi ~ exponential(1);
  
  // P
  
  mu_spec_p ~ normal(0,1.0);
  sig_spec_p ~ exponential(1);
  
  mu_event_p ~ normal(0,0.5);
  sig_event_p ~ exponential(1);
  
  // Missingness RE Priors
  
  alpha_miss ~ exponential(1);
  beta_miss ~ exponential(1);
  
  // Missingness - Create adult vector including missingness
  vector[n] vec_p_adult =  zeros_vector(n);
  
  for(i in 1:n){
  
    if(adult[i]!=2){
    
      vec_p_adult[i] = adult[i];
    
    } else{
    
      vec_p_adult[i] = p_adult[ind_missing[i]];

    }
  }
  
  // Linear predictors
  
  vector[n] lin_phi = a_spec_phi[species] + b_adult_phi*vec_p_adult 
                      + b_mon_phi[month] + b_wint_phi[winter];
                      
  vector[n] lin_psi = a_spec_psi[species] + b_wint_psi[winter];
  
  vector[n] lin_pi = october.*lin_psi + (1-october).*lin_phi;
  
  vector[n] lin_p = a_spec_p[species] + b_adult_p*vec_p_adult 
                      + b_event_p[event];
  
  // Probabilities from predictors
  
  vector[n] pi = inv_logit(lin_pi);
  vector[n] p = inv_logit(lin_p);
  
  
  // Likelihood

  target += log(pi) + present.*log(p) + (1-present).*log1m(p);
  
  // Chi Parameter
  
  // Missingness - Create adult vector including missingness
  vector[n_chi] vec_p_adult_chi =  zeros_vector(n_chi);
  
  for(i in 1:n_chi){
  
    if(adult_chi[i]!=2){
    
      vec_p_adult_chi[i] = adult_chi[i];
    
    } else{
    
      vec_p_adult_chi[i] = p_adult[ind_missing_chi[i]];

    }
  }
  
  vector[n_chi] lin_phi_chi = a_spec_phi[species_chi] + b_adult_phi*vec_p_adult_chi 
                      + b_mon_phi[month_chi] + b_wint_phi[winter_chi];
                      
  vector[n_chi] lin_psi_chi = a_spec_psi[species_chi] + b_wint_psi[winter_chi];
  
  vector[n_chi] lin_pi_chi = october_chi.*lin_psi_chi 
                              + (1-october_chi).*lin_phi_chi;
  
  vector[n_chi] lin_p_chi = a_spec_p[species_chi] + b_adult_p*vec_p_adult_chi 
                      + b_event_p[event_chi];
  
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
  
  target += log(chi[l]);
  
}"


# Compiles the model
stan_m_D <- stan_model(model_name = "stan_m_D",model_code=code_m_D)


# Create list for model
outputs_D <- f_prep_for_stan(df_model)


# Sample
samples_m_D <-  sampling(stan_m_D,
                         data = outputs_D$stan_list, 
                         chains=1, iter = 1000,
                         cores = parallel::detectCores())







#######################################


# Prepare data

df_sim_stan <- df_sim|> 
  filter(total_count_obs[Individual]>1)|>
  group_by(Individual)|>
  mutate(Cumulative = cumsum(Present))|>
  filter(Cumulative != 0)|>
  ungroup()

# Relabel individuals
relabel <- rep(0, nrow(df_sim_stan))

for(i in 1:nrow(df_sim_stan)){
  
  relabel[i] = which(df_sim_stan$Individual[i]==unique(df_sim_stan$Individual))  
  
}

# df_sim_stan now includes for each individual every row from the first time
# they are observed to the final event T
df_sim_stan <- df_sim_stan|>
  mutate(Individual = relabel)


# Total Counts for Filtered and Relabelled Individuals
total_count_relabel <- df_sim_stan|>
  group_by(Individual)|>
  summarise(Count = sum(Present))|>
  pull(Count)

# Observed Calculations, f(i) + 1 to l(i)
df_sim_stan_obs <- df_sim_stan|>
  filter(!(Cumulative == 1 & Present == 1), # remove f(i)
         !(Cumulative == total_count_relabel[Individual] & Present ==0)) # remove after l(i)
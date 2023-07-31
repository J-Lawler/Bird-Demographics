
########## Preamble ##########

# Load libraries
library(VIM)
library(tidyverse)


months <- c("January", "February", "March", "April", "May", "June", "July", "August",
            "September", "October", "November", "December")

species <- c("Blackcap","Chiffchaff","Robin")


########## Contents ##########

# Data Cleaning

  ## Blackcap

  ## Chiffchaff

  ## Robin

  ## Combined


# Create Annual Data


# Data Exploration


########## Data Cleaning ##########




##### Blackcap


df_blackcap_raw <- read_csv2("Data_Raw/CR.blackcap_FixRing.csv")


df_blackcap <- df_blackcap_raw|>
  mutate(Individual = 1:nrow(df_blackcap_raw))|>
  pivot_longer(-c(Individual, age_at_ringing), names_to = "Date", values_to = "Present")|>
  mutate(Winter = as.numeric(str_sub(Date,1,4)),
         Month_Num = as.numeric(str_sub(Date,5,6)),
         Year = if_else(Month_Num < 6, Winter + 1, Winter), 
         Month = months[Month_Num])|>
  select(Year, Winter, Month, Month_Num, Individual,
         Age_Ring = age_at_ringing, Present)|>
  arrange(Individual, Year, Month_Num)


write_csv(df_blackcap, file = "Data_Clean/df_blackcap.csv")

# Number of individual blackcaps
num_blackcap <- max(df_blackcap$Individual)





##### chiffchaff


df_chiffchaff_raw <- read_csv2("Data_Raw/CR.chifchaf_FixRing.csv")


df_chiffchaff <- df_chiffchaff_raw|>
  mutate(Individual = 1:nrow(df_chiffchaff_raw))|>
  pivot_longer(-c(Individual, age_at_ringing), names_to = "Date", values_to = "Present")|>
  mutate(Winter = as.numeric(str_sub(Date,1,4)),
         Month_Num = as.numeric(str_sub(Date,5,6)),
         Year = if_else(Month_Num < 6, Winter + 1, Winter), 
         Month = months[Month_Num],
         Individual = Individual + num_blackcap)|> # relabel individuals
  select(Year, Winter, Month, Month_Num, Individual,
         Age_Ring = age_at_ringing, Present)|>
  arrange(Individual, Year, Month_Num)


write_csv(df_chiffchaff, file = "Data_Clean/df_chiffchaff.csv")

num_blackcaps_chiffchaffs <- max(df_chiffchaff$Individual)


##### Robin


df_robin_raw <- read_csv2("Data_Raw/CR.robin_FixRing.csv")


df_robin <- df_robin_raw|>
  mutate(Individual = 1:nrow(df_robin_raw))|>
  pivot_longer(-c(Individual, age_at_ringing), names_to = "Date", values_to = "Present")|>
  mutate(Winter = as.numeric(str_sub(Date,1,4)),
         Month_Num = as.numeric(str_sub(Date,5,6)),
         Year = if_else(Month_Num < 6, Winter + 1, Winter),
         Month = months[Month_Num],
         Individual = Individual + num_blackcaps_chiffchaffs)|> # relabel individuals
  select(Year, Winter, Month, Month_Num, Individual,
         Age_Ring = age_at_ringing, Present)|>
  arrange(Individual, Year, Month_Num)



write_csv(df_robin, file = "Data_Clean/df_robin.csv")





##### Combined

df_combined <- bind_rows(df_blackcap, df_chiffchaff, df_robin,
                         .id = "Species")|>
  mutate(Species = species[as.numeric(Species)])



write_csv(df_combined, file = "Data_Clean/df_combined.csv")




########## Age Covariate ##########


# Read in combined file, convert to annual (including age covariate)
# Calculative cumulative sightings 
df_age <- read_csv(file = "Data_Clean/df_combined.csv")|>
  group_by(Species, Individual, Winter, Age_Ring)|>
  summarise(Present = sum(Present))|>
  mutate(Present = if_else(Present>=1,1,0))|>
  select(Species, Winter, Individual, Present, Age_Ring)|>
  ungroup()|>
  group_by(Individual)|>
  mutate(Cumulative = cumsum(Present))|>
  ungroup()


# Here:
# 1 - Juvenile
# 2 - Adult
# 3 - Unknown
# 4 - NA
f_age <- function(cumulative,present,age_orig){
  
  ages <- c("juvenile", "adult", "Unknown")
  
  if(cumulative==0){
    
    return(4)
    
  } else if(cumulative==1 & present == 1){
    
    return(which(age_orig==ages))
    
  } else{
    
    return(2)
    
  }
}


age_current <- rep(0,nrow(df_age))

for(i in 1:length(age_current)){
  
  age_current[i] <- f_age(cumulative = df_age$Cumulative[i],
                          present = df_age$Present[i],
                          age_orig = df_age$Age_Ring[i])
  
}


saveRDS(age_current, file = "Data_Clean/age_current.RDS")


########## Annual Data ##########

# This data condenses each 7 month winter period into a single observation

# If an individual is observed at any point over winter, they get a 1. Else 0.

# If they are observed multiple times, they just get a 1.


df_combined <- read_csv("Data_Clean/df_combined.csv")
age_current <- readRDS(file = "Data_Clean/age_current.RDS")

df_annual <- df_combined|>
  group_by(Species, Individual, Winter)|>
  summarise(Present = sum(Present))|>
  mutate(Present = if_else(Present>=1,1,0))|>
  select(Species, Winter, Individual, Present)|>
  ungroup()|>
  mutate(Age = age_current)


write_csv(df_annual, file = "Data_Clean/df_annual.csv")




########## Weather Covariate ##########






########## Monthly data with covariates ##########

df_current_age <- read_csv(file = "Data_Clean/df_annual.csv")|>
  select(Individual, Winter, Age)


df_combined <- read_csv("Data_Clean/df_combined.csv")

df_events <- df_combined|>
  select(Year,Month)|>
  unique()

df_events <- df_events|>
  mutate(Event = 1:nrow(df_events))

df_model <- df_combined|>
  left_join(df_current_age,
            by = c("Individual" = "Individual",
                   "Winter" = "Winter"))|>
  left_join(df_events,
            by = c("Year" = "Year",
                   "Month" = "Month"))


write_csv(df_model, file = "Data_Clean/df_model.csv")


########## Data Exploration ##########


# Plot of observed count each month, grouped by year, faceted by species

# Should change plot to start in October, go through April(?) and stop
# remove implicit interpolation through the summer months
df_combined|>
  mutate(Year = factor(Year, levels = unique(df_combined$Year)))|>
  group_by(Species,Year, Month_Num)|>
  summarise(Count = sum(Present))|>
  ggplot(aes(x = Month_Num, y = Count, colour = Year))+
  geom_line()+
  facet_wrap(~Species)+
  scale_x_continuous(labels = function(x) months[x])+
  scale_color_viridis_d()+
  theme_bw()+
  xlab("Month")+
  ggtitle("Observation Count by Month")+
  theme(axis.text.x = element_text(angle = 90))




# Plot of how many birds are seen across multiple winters

df_annual <- read_csv("Data_Clean/df_annual.csv")


df_annual|>
  group_by(Individual)|>
  summarise(Present = sum(Present))|>
  ggplot(aes(x = Present))+
  geom_histogram(fill = "#440154")+
  theme_bw()+
  xlab("Number of Winters Observed")+
  ylab("Count of Individuals")+
  ggtitle("How Many Birds are Seen Over Multiple Winters?")






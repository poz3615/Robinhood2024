library(dplyr)
library(tidyr)
library(ggplot2)
library(rjags)
library(R2jags)
library(readxl)
library(truncnorm)
library(scoringRules)
library(viridis)
library(loo)
library(patchwork)


######################################### Data Cleaning #########################################
# Mosquito development rate in days, a lot small so filtered dataset to >1 to try to get the hang of it
# Read data, get correct columns, correct metric
mosquito.data <- read.csv("development-rate-observations.clean.csv")
mosquito.data <- mosquito.data[,c("Trait", "originaltraitunit","interactor1", "Temp", "interactor1number", 
                                          "originalerrorpos", "originalerrorunit")]
mosquito.data <- mosquito.data |> 
  mutate(
    Trait = ifelse(`originaltraitunit` == "hours", Trait / 24, Trait),
    originalerrorpos = ifelse(`originaltraitunit` == "hours", originalerrorpos / 24, originalerrorpos),
    `originaltraitunit` = "days"  # Standardizing the unit for all rows
  )
mosquito.data <- mosquito.data %>%
  filter(!(is.na(originalerrorpos) & interactor1number != 1)) # Remove rows where OriginalErrorPos is NA and Interactor1Number is not 1



library(dplyr)
# Small- 0-2.55
# Medium- 2.56-3.05
# Large - 3.06-Inf
mosquito.data <- mosquito.data %>%
  mutate(wing.size = case_when(
    interactor1 == "Aedes aegypti" ~ "medium",
    interactor1 == "Aedes albopictus" ~ "small",
    interactor1 == "Aedes japonicus japonicus" ~ "large",
    interactor1 == "Aedes atropalpus" ~ "large",
    interactor1 == "Aedes krombeini" ~ "medium",
    interactor1 == "Aedes notoscriptus" ~ "medium",
    interactor1 == "Aedes camptorhynchus" ~ "large",
    TRUE ~ NA_character_
  ))

mosquito.data <- mosquito.data %>%
  mutate(climate = case_when(
    interactor1 == "Aedes aegypti" ~ "tropical",
    interactor1 == "Aedes albopictus" ~ "tropical",
    interactor1 == "Aedes japonicus japonicus" ~ "temperate",
    interactor1 == "Aedes atropalpus" ~ "temperate",
    interactor1 == "Aedes krombeini" ~ "tropical",
    interactor1 == "Aedes notoscriptus" ~ "subtropical",
    interactor1 == "Aedes camptorhynchus" ~ "subtropical",
    TRUE ~ NA_character_
  ))

mosquito.data$climate <- as.factor(mosquito.data$climate)
mosquito.data$wing.size <- as.factor(mosquito.data$wing.size)
# Plotted for each species
ggplot(mosquito.data, aes(y = Trait, x = Temp)) +
  geom_point() +
  facet_wrap(~ interactor1) +
  labs(title = "Temperature by Trait for Each Factor Level",
       x = "Temperature",
       y = "Development Time") +
  theme_minimal()

# Create a new expanded dataset with pseudo-distribution data
m.expanded_df1 <- data.frame()

# Get the data that is NOT individual
m.filtered_df1 <- mosquito.data |> filter(interactor1number != 1)
m.filtered_df1 <- m.filtered_df1[m.filtered_df1$originalerrorunit %in% c("se", "sem", "SEM", "sd"), ]
# Loop over each row in the dataset of mean data
for (i in 1:nrow(m.filtered_df1)) {
  # Extract values for each row
  # Mean
  mean_value <- m.filtered_df1$Trait[i]
  
  # Standard error times the sqrt(n) to get sd
  if (m.filtered_df1$originalerrorunit[i] %in% c("se", "sem", "SEM")){
    sd_value <- m.filtered_df1$originalerrorpos[i] * sqrt(m.filtered_df1$interactor1number[i])
  }
  else{
    sd_value <- m.filtered_df1$originalerrorpos[i]
  }
  # Species
  interactor1 <- m.filtered_df1$interactor1[i]
  # Temp
  interactor1temp <- m.filtered_df1$Temp[i]
  # Sample size
  num_points <- m.filtered_df1$interactor1number[i]
  # Climate
  climate <- m.filtered_df1$climate[i]
  # Size
  wing.size <- m.filtered_df1$wing.size[i]
  # Generate sample size of random values from truncated normal (0, Inf) distribution
  generated_values <- rtruncnorm(num_points, a=0, b=Inf, mean = mean_value, sd = sd_value)
  # Create a temporary dataframe for the generated values
  temp_df <- data.frame(
    # Simulated values
    Trait = generated_values,
    # Carry over sd, species, temperature, climate, and size
    sd = sd_value,
    Species = interactor1,
    Temp = interactor1temp,
    climate = climate,
    wing.size = wing.size
  )
  # Append the generated rows to the expanded dataset
  m.expanded_df1 <- rbind(m.expanded_df1, temp_df)
}

# Take the individual data
m.filtered_df <- mosquito.data |> filter(interactor1number == 1)
# Summarize it, grouping by species and temperature, to get the mean, sd, and sample size
m.summarized_df <- m.filtered_df %>%
  group_by(interactor1, Temp, climate, wing.size) %>%
  summarize(
    mean_value = mean(Trait, na.rm = TRUE),
    sd_value = sd(Trait, na.rm = TRUE),
    count = n(),
    .groups = "drop"
  )
# Resample n data points with appropriate mean and standard deviation
m.expanded_df <- m.summarized_df %>%
  rowwise() %>%
  do(data.frame(
    Species = .$interactor1,
    Temp = .$Temp,
    Trait = rtruncnorm(.$count, a = 0, b = Inf, mean = .$mean_value, sd = .$sd_value),
    climate = .$climate,
    wing.size = .$wing.size  
  ))

# Combine the new mean data and the new individual data, selecting 
mosquito.pd <- rbind(m.expanded_df1[, c(1, 3:6)], m.expanded_df)
# Get rid of NAs
mosquito.pd <- na.omit(mosquito.pd)
# Drop unused levels of species
mosquito.pd$Species <- as.factor(mosquito.pd$Species)
mosquito.pd$Species <- droplevels(mosquito.pd$Species)


# Plot new data
ggplot(mosquito.pd, aes(y = Trait, x = Temp)) +
  geom_point() +
  facet_wrap(~ Species) +
  labs(title = "Temperature by Trait for Each Factor Level",
       x = "Temperature",
       y = "Development Time") +
  theme_minimal()



# Species number for hierarchical model
m.species.pd <- as.factor(mosquito.pd$Species)
# Number of species
m.n.species.pd <- length(unique(mosquito.pd$Species))

# Climate number for hierarchical model
m.climate.pd <- as.factor(mosquito.pd$climate)
# Number of climates
m.n.climate.pd <- length(unique(mosquito.pd$climate))

# Size number for hierarchical model
m.size.pd <- as.factor(mosquito.pd$Size)
# Small: 0-3
# Medium: 3.01-5
# Large: 5.01+
# Number of sizes
m.n.size.pd <- length(unique(mosquito.pd$wing.size))

m.N.obs.pd <- length(mosquito.pd$Trait)

######################################### By Full Data #########################################
sink("aedes.one.data.txt")
cat("
model {
  # Priors
   a ~ dexp(0.1)             # Width
  Topt ~ dnorm(30, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  c ~ dexp(0.001)
  
  # Likelihood
  for (i in 1:m.N.obs.pd) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], tau1) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], tau1)
  }}
", file = "aedes.one.data.txt")


# Parameters to monitor
parameters <- c("a", "c", "tau1", "Topt", "log_lik")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Initial Values
inits = vector('list', nc)
GenInits = function() {
  a <- rexp(1, 0.1)
  c <- rexp(1, 0.001)
  tau1 <- rexp(1, 0.01)
  Topt <- rnorm(1, 30, 1/0.001)
  list(
    a = a,             
    c = c,
    tau1 = tau1,
    Topt = Topt
  )
}
for(i in 1:nc){
  inits[[i]] = GenInits()
}

# Bundle all data in a list for JAGS

jag.data<-list(trait = mosquito.pd$Trait, m.N.obs.pd = m.N.obs.pd, temp = mosquito.pd$Temp)

m.one.data.pd.fit <- jags(data=jag.data, parameters.to.save=parameters, 
                       model.file="aedes.one.data.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                       n.iter=ni, DIC=T, working.directory=getwd())
m.one.data.pd.fit.mcmc <- as.mcmc(m.one.data.pd.fit) ## makes an "mcmc" object
closeAllConnections()
m.one.data.pd.fit$BUGSoutput$DIC

# Extract log-likelihood matrix (samples × observations)
m.log.lik.full.pd <- m.one.data.pd.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
m.waic.full.pd <- waic(m.log.lik.full.pd)$estimates["waic", "Estimate"]


# Get posterior means
m.summary_fit.pd <- summary(m.one.data.pd.fit.mcmc)$statistics
m.a_mean.pd <- m.summary_fit.pd["a", "Mean"]
m.c_mean.pd <- m.summary_fit.pd["c", "Mean"]
m.means_Topt.pd <- m.summary_fit.pd["Topt", "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)

# Get predicted values
m.full.data.pd <- function (x){
  m.a_mean.pd * (x - m.means_Topt.pd)^2 + m.c_mean.pd
}


# Optionally, convert to a data frame for easier plotting
m.fitted_data.one.data.pd <- data.frame(Temperature = rep(m.temp_seq, m.n.species.pd), 
                                     Fitted = m.full.data(m.temp_seq),
                                     Species = factor(rep(c("Aedes aegypti", "Aedes albopictus",
                                                            "Aedes atropalpus", "Aedes camptorhynchus",
                                                            "Aedes japonicus japonicus",
                                                            "Aedes krombeini", "Aedes notoscriptus"), each = length(m.temp_seq))))

m.plot_data.one.data.pd <- data.frame(
  Temperature = mosquito.pd$Temp,
  Trait = mosquito.pd$Trait,
  Species = as.factor(mosquito.pd$interactor1)
)

# Plot observed data and fitted curves with facet wrap by species
m.p.full.pd <- ggplot() +
  geom_point(data = m.plot_data.one.data.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.one.data.pd, aes(x = Temperature, y = Fitted), color = "blue", size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) 
ggsave("mosquito.full.plot.png", width = 8, height = 5)


# Merge the observed and predicted data for comparison
m.merged_data.one.data.pd <- m.fitted_data.one.data.pd %>%
  left_join(plot_data.one.data, by = c("Temperature", "Species"))


# Calculate RMSE for each species
m.mse_by_species.one.data.pd <- m.merged_data.one.data.pd %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.pd <- m.merged_data.one.data.pd %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(m.mse_by_species.one.data.pd,
     m.fitted_data.one.data.pd,
     m.merged_data.one.data.pd,
     m.one.data.pd.fit$BUGSoutput$DIC,
     m.overall_rmse.pd,
     m.waic.full.pd,
     m.p.full.pd,
     file="m.one.data.pd.fit.RData")

######################################### By Species Individual Aegypti #########################################
a.aegypti.pd <- mosquito.pd[mosquito.pd$interactor1 == "Aedes aegypti", ]

m.N.obs.pd.ind <- length(a.aegypti.pd$Trait)

sink("aedes.individual.txt")
cat("
model {
  # Priors
   a ~ dexp(0.1)             # Width
  Topt ~ dnorm(30, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  c ~ dexp(0.001)
  
  # Likelihood
  for (i in 1:m.N.obs.pd.ind) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], tau1) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], tau1)
  }}
", file = "aedes.individual.txt")



# Parameters to monitor
parameters <- c("a", "c", "tau1", "Topt", "log_lik")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Initial Values
inits = vector('list', nc)
GenInits = function() {
  a <- rexp(1, 0.1)
  c <- rexp(1, 0.001)
  tau1 <- rexp(1, 0.01)
  Topt <- rnorm(1, 30, 1/0.001)
  list(
    a = a,             
    c = c,
    tau1 = tau1,
    Topt = Topt
  )
}
for(i in 1:nc){
  inits[[i]] = GenInits()
}

# Bundle all data in a list for JAGS

jag.data<-list(trait = a.aegypti.pd$Trait, m.N.obs.pd.ind = m.N.obs.pd.ind, temp = a.aegypti.pd$Temp)

a.aegypti.pd.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                         model.file="aedes.individual.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                         n.iter=ni, DIC=T, working.directory=getwd())
a.aegypti.pd.fit.mcmc <- as.mcmc(a.aegypti.pd.fit) ## makes an "mcmc" object
closeAllConnections()
a.aegypti.pd.fit$BUGSoutput$DIC


# Get posterior means
m.summary_fit.a.aegypti.pd <- summary(a.aegypti.pd.fit.mcmc)$statistics
m.a_mean.pd <- m.summary_fit.a.aegypti.pd["a", "Mean"]
m.c_mean.pd <- m.summary_fit.a.aegypti.pd["c", "Mean"]
m.means_Topt.pd <- m.summary_fit.a.aegypti.pd["Topt", "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)

# Get predicted values
a.aegypti.fitted.pd <- function (x){
  m.a_mean.pd * (x - m.means_Topt.pd)^2 + m.c_mean.pd
}

# Data frame of temperature, fitted value, species
m.fitted_data.a.aegypti.pd <- data.frame(Temperature = rep(m.temp_seq),
                                         Fitted = a.aegypti.fitted.pd(m.temp_seq))
# Original data with temperature, trait, and species as numeric
m.plot_data.a.aegypti.pd <- data.frame(
  Temperature = a.aegypti.pd$Temp,
  Trait = a.aegypti.pd$Trait
)

# Plot observed data and fitted curves 
m.p.a.aegypti.pd <- ggplot() +
  geom_point(data = m.plot_data.a.aegypti.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.a.aegypti.pd, aes(x = Temperature, y = Fitted), color = "seagreen4", size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  theme(legend.position = "none")
m.p.a.aegypti.pd

# Merge the observed and predicted data for comparison
m.merged_data.a.aegypti.pd <- m.fitted_data.a.aegypti.pd %>%
  left_join(a.aegypti, by = c("Temperature" = "Temp"))

# Calculate RMSE for each species
m.rmse.a.aegypti.pd <- m.merged_data.a.aegypti.pd %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE),
    SE = (Trait - Fitted)^2
  )

# Extract log-likelihood matrix (samples × observations)
m.log.lik.a.aegypti.pd <- a.aegypti.pd.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
m.waic.a.aegypti.pd <- waic(m.log.lik.a.aegypti.pd)$estimates["waic", "Estimate"]



save(m.rmse.a.aegypti.pd,
     m.fitted_data.a.aegypti.pd,
     m.merged_data.a.aegypti.pd,
     a.aegypti.pd.fit$BUGSoutput$DIC,
     m.waic.a.aegypti.pd,
     m.p.a.aegypti.pd,
     file="m.a.aegypti.pd.fit.RData")

######################################### By Species Individual Albopictus #########################################
a.albopictus.pd <- mosquito.pd[mosquito.pd$interactor1 == "Aedes albopictus", ]

m.N.obs.pd.ind <- length(a.albopictus.pd$Trait)

sink("aedes.individual.txt")
cat("
model {
  # Priors
   a ~ dexp(0.1)             # Width
  Topt ~ dnorm(30, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  c ~ dexp(0.001)
  
  # Likelihood
  for (i in 1:m.N.obs.pd.ind) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], tau1) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], tau1)
  }}
", file = "aedes.individual.txt")



# Parameters to monitor
parameters <- c("a", "c", "tau1", "Topt", "log_lik")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Initial Values
inits = vector('list', nc)
GenInits = function() {
  a <- rexp(1, 0.1)
  c <- rexp(1, 0.001)
  tau1 <- rexp(1, 0.01)
  Topt <- rnorm(1, 30, 1/0.001)
  list(
    a = a,             
    c = c,
    tau1 = tau1,
    Topt = Topt
  )
}
for(i in 1:nc){
  inits[[i]] = GenInits()
}

# Bundle all data in a list for JAGS

jag.data<-list(trait = a.albopictus.pd$Trait, m.N.obs.pd.ind = m.N.obs.pd.ind, temp = a.albopictus.pd$Temp)

a.albopictus.pd.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                            model.file="aedes.individual.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                            n.iter=ni, DIC=T, working.directory=getwd())
a.albopictus.pd.fit.mcmc <- as.mcmc(a.albopictus.pd.fit) ## makes an "mcmc" object
closeAllConnections()
a.albopictus.pd.fit$BUGSoutput$DIC


# Get posterior means
m.summary_fit.a.albopictus.pd <- summary(a.albopictus.pd.fit.mcmc)$statistics
m.a_mean.pd <- m.summary_fit.a.albopictus.pd["a", "Mean"]
m.c_mean.pd <- m.summary_fit.a.albopictus.pd["c", "Mean"]
m.means_Topt.pd <- m.summary_fit.a.albopictus.pd["Topt", "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)

# Get predicted values
a.albopictus.fitted.pd <- function (x){
  m.a_mean.pd * (x - m.means_Topt.pd)^2 + m.c_mean.pd
}

# Data frame of temperature, fitted value, species
m.fitted_data.a.albopictus.pd <- data.frame(Temperature = rep(m.temp_seq),
                                            Fitted = a.albopictus.fitted.pd(m.temp_seq))
# Original data with temperature, trait, and species as numeric
m.plot_data.a.albopictus.pd <- data.frame(
  Temperature = a.albopictus.pd$Temp,
  Trait = a.albopictus.pd$Trait
)

# Plot observed data and fitted curves 
m.p.a.albopictus.pd <- ggplot() +
  geom_point(data = m.plot_data.a.albopictus.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.a.albopictus.pd, aes(x = Temperature, y = Fitted), color = "seagreen4", size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  theme(legend.position = "none")
m.p.a.albopictus.pd

# Merge the observed and predicted data for comparison
m.merged_data.a.albopictus.pd <- m.fitted_data.a.albopictus.pd %>%
  left_join(a.albopictus, by = c("Temperature" = "Temp"))

# Calculate RMSE for each species
m.rmse.a.albopictus.pd <- m.merged_data.a.albopictus.pd %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE),
    SE = (Trait - Fitted)^2
  )

# Extract log-likelihood matrix (samples × observations)
m.log.lik.a.albopictus.pd <- a.albopictus.pd.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
m.waic.a.albopictus.pd <- waic(m.log.lik.a.albopictus.pd)$estimates["waic", "Estimate"]


save(m.rmse.a.albopictus.pd,
     m.fitted_data.a.albopictus.pd,
     m.merged_data.a.albopictus.pd,
     a.albopictus.pd.fit$BUGSoutput$DIC,
     m.waic.a.albopictus.pd,
     m.p.a.albopictus.pd,
     file="m.a.albopictus.pd.fit.RData")
######################################### By Species Individual Atropalpus #########################################
a.atropalpus.pd <- mosquito.pd[mosquito.pd$interactor1 == "Aedes atropalpus", ]

m.N.obs.pd.ind <- length(a.atropalpus.pd$Trait)

sink("aedes.individual.txt")
cat("
model {
  # Priors
   a ~ dexp(0.1)             # Width
  Topt ~ dnorm(30, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  c ~ dexp(0.001)
  
  # Likelihood
  for (i in 1:m.N.obs.pd.ind) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], tau1) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], tau1)
  }}
", file = "aedes.individual.txt")


# Parameters to monitor
parameters <- c("a", "c", "tau1", "Topt", "log_lik")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Initial Values
inits = vector('list', nc)
GenInits = function() {
  a <- rexp(1, 0.1)
  c <- rexp(1, 0.001)
  tau1 <- rexp(1, 0.01)
  Topt <- rnorm(1, 30, 1/0.001)
  list(
    a = a,             
    c = c,
    tau1 = tau1,
    Topt = Topt
  )
}
for(i in 1:nc){
  inits[[i]] = GenInits()
}

# Bundle all data in a list for JAGS

jag.data<-list(trait = a.atropalpus.pd$Trait, m.N.obs.pd.ind = m.N.obs.pd.ind, temp = a.atropalpus.pd$Temp)

a.atropalpus.pd.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                            model.file="aedes.individual.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                            n.iter=ni, DIC=T, working.directory=getwd())
a.atropalpus.pd.fit.mcmc <- as.mcmc(a.atropalpus.pd.fit) ## makes an "mcmc" object
closeAllConnections()
a.atropalpus.pd.fit$BUGSoutput$DIC


# Get posterior means
m.summary_fit.a.atropalpus.pd <- summary(a.atropalpus.pd.fit.mcmc)$statistics
m.a_mean.pd <- m.summary_fit.a.atropalpus.pd["a", "Mean"]
m.c_mean.pd <- m.summary_fit.a.atropalpus.pd["c", "Mean"]
m.means_Topt.pd <- m.summary_fit.a.atropalpus.pd["Topt", "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)

# Get predicted values
a.atropalpus.fitted.pd <- function (x){
  m.a_mean.pd * (x - m.means_Topt.pd)^2 + m.c_mean.pd
}

# Data frame of temperature, fitted value, species
m.fitted_data.a.atropalpus.pd <- data.frame(Temperature = rep(m.temp_seq),
                                            Fitted = a.atropalpus.fitted.pd(m.temp_seq))
# Original data with temperature, trait, and species as numeric
m.plot_data.a.atropalpus.pd <- data.frame(
  Temperature = a.atropalpus.pd$Temp,
  Trait = a.atropalpus.pd$Trait
)

# Plot observed data and fitted curves 
m.p.a.atropalpus.pd <- ggplot() +
  geom_point(data = m.plot_data.a.atropalpus.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.a.atropalpus.pd, aes(x = Temperature, y = Fitted), color = "seagreen4", size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  theme(legend.position = "none")
m.p.a.atropalpus.pd

# Merge the observed and predicted data for comparison
m.merged_data.a.atropalpus.pd <- m.fitted_data.a.atropalpus.pd %>%
  left_join(a.atropalpus, by = c("Temperature" = "Temp"))

# Calculate RMSE for each species
m.rmse.a.atropalpus.pd <- m.merged_data.a.atropalpus.pd %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE),
    SE = (Trait - Fitted)^2
  )

# Extract log-likelihood matrix (samples × observations)
m.log.lik.a.atropalpus.pd <- a.atropalpus.pd.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
m.waic.a.atropalpus.pd <- waic(m.log.lik.a.atropalpus.pd)$estimates["waic", "Estimate"]
# 

save(m.rmse.a.atropalpus.pd,
     m.fitted_data.a.atropalpus.pd,
     m.merged_data.a.atropalpus.pd,
     a.atropalpus.pd.fit$BUGSoutput$DIC,
     m.waic.a.atropalpus.pd,
     m.p.a.atropalpus.pd,
     file="m.a.atropalpus.pd.fit.RData")
######################################### By Species Individual Camptorhynchus #########################################
a.camptorhynchus.pd <- mosquito.pd[mosquito.pd$interactor1 == "Aedes camptorhynchus", ]


m.N.obs.pd.ind <- length(a.camptorhynchus.pd$Trait)

sink("aedes.individual.txt")
cat("
model {
  # Priors
   a ~ dexp(0.1)             # Width
  Topt ~ dnorm(30, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  c ~ dexp(0.001)
  
  # Likelihood
  for (i in 1:m.N.obs.pd.ind) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], tau1) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], tau1)
  }}
", file = "aedes.individual.txt")



# Parameters to monitor
parameters <- c("a", "c", "tau1", "Topt", "log_lik")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Initial Values
inits = vector('list', nc)
GenInits = function() {
  a <- rexp(1, 0.1)
  c <- rexp(1, 0.001)
  tau1 <- rexp(1, 0.01)
  Topt <- rnorm(1, 30, 1/0.001)
  list(
    a = a,             
    c = c,
    tau1 = tau1,
    Topt = Topt
  )
}
for(i in 1:nc){
  inits[[i]] = GenInits()
}

# Bundle all data in a list for JAGS

jag.data<-list(trait = a.camptorhynchus.pd$Trait, m.N.obs.pd.ind = m.N.obs.pd.ind, temp = a.camptorhynchus.pd$Temp)

a.camptorhynchus.pd.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                                model.file="aedes.individual.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                                n.iter=ni, DIC=T, working.directory=getwd())
a.camptorhynchus.pd.fit.mcmc <- as.mcmc(a.camptorhynchus.pd.fit) ## makes an "mcmc" object
closeAllConnections()
a.camptorhynchus.pd.fit$BUGSoutput$DIC


# Get posterior means
m.summary_fit.a.camptorhynchus.pd <- summary(a.camptorhynchus.pd.fit.mcmc)$statistics
m.a_mean.pd <- m.summary_fit.a.camptorhynchus.pd["a", "Mean"]
m.c_mean.pd <- m.summary_fit.a.camptorhynchus.pd["c", "Mean"]
m.means_Topt.pd <- m.summary_fit.a.camptorhynchus.pd["Topt", "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)

# Get predicted values
a.camptorhynchus.fitted.pd <- function (x){
  m.a_mean.pd * (x - m.means_Topt.pd)^2 + m.c_mean.pd
}

# Data frame of temperature, fitted value, species
m.fitted_data.a.camptorhynchus.pd <- data.frame(Temperature = rep(m.temp_seq),
                                                Fitted = a.camptorhynchus.fitted.pd(m.temp_seq))
# Original data with temperature, trait, and species as numeric
m.plot_data.a.camptorhynchus.pd <- data.frame(
  Temperature = a.camptorhynchus.pd$Temp,
  Trait = a.camptorhynchus.pd$Trait
)

# Plot observed data and fitted curves 
m.p.a.camptorhynchus.pd <- ggplot() +
  geom_point(data = m.plot_data.a.camptorhynchus.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.a.camptorhynchus.pd, aes(x = Temperature, y = Fitted), color = "seagreen4", size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  theme(legend.position = "none")
m.p.a.camptorhynchus.pd

# Merge the observed and predicted data for comparison
m.merged_data.a.camptorhynchus.pd <- m.fitted_data.a.camptorhynchus.pd %>%
  left_join(a.camptorhynchus, by = c("Temperature" = "Temp"))

# Calculate RMSE for each species
m.rmse.a.camptorhynchus.pd <- m.merged_data.a.camptorhynchus.pd %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE),
    SE = (Trait - Fitted)^2
  )

# Extract log-likelihood matrix (samples × observations)
m.log.lik.a.camptorhynchus.pd <- a.camptorhynchus.pd.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
m.waic.a.camptorhynchus.pd <- waic(m.log.lik.a.camptorhynchus.pd)$estimates["waic", "Estimate"]



save(m.rmse.a.camptorhynchus.pd,
     m.fitted_data.a.camptorhynchus.pd,
     m.merged_data.a.camptorhynchus.pd,
     a.camptorhynchus.pd.fit$BUGSoutput$DIC,
     m.waic.a.camptorhynchus.pd,
     m.p.a.camptorhynchus.pd,
     file="m.a.camptorhynchus.pd.fit.RData")
######################################### By Species Individual Japonicus Japonicus #########################################
a.jj.pd <- mosquito.pd[mosquito.pd$interactor1 == "Aedes japonicus japonicus", ]


m.N.obs.pd.ind <- length(a.jj.pd$Trait)

sink("aedes.individual.txt")
cat("
model {
  # Priors
   a ~ dexp(0.1)             # Width
  Topt ~ dnorm(30, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  c ~ dexp(0.001)
  
  # Likelihood
  for (i in 1:m.N.obs.pd.ind) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], tau1) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], tau1)
  }}
", file = "aedes.individual.txt")


# Parameters to monitor
parameters <- c("a", "c", "tau1", "Topt", "log_lik")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Initial Values
inits = vector('list', nc)
GenInits = function() {
  a <- rexp(1, 0.1)
  c <- rexp(1, 0.001)
  tau1 <- rexp(1, 0.01)
  Topt <- rnorm(1, 30, 1/0.001)
  list(
    a = a,             
    c = c,
    tau1 = tau1,
    Topt = Topt
  )
}
for(i in 1:nc){
  inits[[i]] = GenInits()
}

# Bundle all data in a list for JAGS

jag.data<-list(trait = a.jj.pd$Trait, m.N.obs.pd.ind = m.N.obs.pd.ind, temp = a.jj.pd$Temp)

a.jj.pd.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                    model.file="aedes.individual.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                    n.iter=ni, DIC=T, working.directory=getwd())
a.jj.pd.fit.mcmc <- as.mcmc(a.jj.pd.fit) ## makes an "mcmc" object
closeAllConnections()
a.jj.pd.fit$BUGSoutput$DIC


# Get posterior means
m.summary_fit.a.jj.pd <- summary(a.jj.pd.fit.mcmc)$statistics
m.a_mean.pd <- m.summary_fit.a.jj.pd["a", "Mean"]
m.c_mean.pd <- m.summary_fit.a.jj.pd["c", "Mean"]
m.means_Topt.pd <- m.summary_fit.a.jj.pd["Topt", "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)

# Get predicted values
a.jj.fitted.pd <- function (x){
  m.a_mean.pd * (x - m.means_Topt.pd)^2 + m.c_mean.pd
}

# Data frame of temperature, fitted value, species
m.fitted_data.a.jj.pd <- data.frame(Temperature = rep(m.temp_seq),
                                    Fitted = a.jj.fitted.pd(m.temp_seq))
# Original data with temperature, trait, and species as numeric
m.plot_data.a.jj.pd <- data.frame(
  Temperature = a.jj.pd$Temp,
  Trait = a.jj.pd$Trait
)

# Plot observed data and fitted curves 
m.p.a.jj.pd <- ggplot() +
  geom_point(data = m.plot_data.a.jj.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.a.jj.pd, aes(x = Temperature, y = Fitted), color = "seagreen4", size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  theme(legend.position = "none")
m.p.a.jj.pd

# Merge the observed and predicted data for comparison
m.merged_data.a.jj.pd <- m.fitted_data.a.jj.pd %>%
  left_join(a.jj, by = c("Temperature" = "Temp"))

# Calculate RMSE for each species
m.rmse.a.jj.pd <- m.merged_data.a.jj.pd %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE),
    SE = (Trait - Fitted)^2
  )

# Extract log-likelihood matrix (samples × observations)
m.log.lik.a.jj.pd <- a.jj.pd.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
m.waic.a.jj.pd <- waic(m.log.lik.a.jj.pd)$estimates["waic", "Estimate"]



save(m.rmse.a.jj.pd,
     m.fitted_data.a.jj.pd,
     m.merged_data.a.jj.pd,
     a.jj.pd.fit$BUGSoutput$DIC,
     m.waic.a.jj.pd,
     m.p.a.jj.pd,
     file="m.a.jj.pd.fit.RData")
######################################### By Species Individual Krombeini #########################################
a.krombeini.pd <- mosquito.pd[mosquito.pd$interactor1 == "Aedes krombeini", ]

m.N.obs.pd.ind <- length(a.krombeini.pd$Trait)

sink("aedes.individual.txt")
cat("
model {
  # Priors
   a ~ dexp(0.1)             # Width
  Topt ~ dnorm(30, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  c ~ dexp(0.001)
  
  # Likelihood
  for (i in 1:m.N.obs.pd.ind) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], tau1) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], tau1)
  }}
", file = "aedes.individual.txt")


# Parameters to monitor
parameters <- c("a", "c", "tau1", "Topt", "log_lik")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Initial Values
inits = vector('list', nc)
GenInits = function() {
  a <- rexp(1, 0.1)
  c <- rexp(1, 0.001)
  tau1 <- rexp(1, 0.01)
  Topt <- rnorm(1, 30, 1/0.001)
  list(
    a = a,             
    c = c,
    tau1 = tau1,
    Topt = Topt
  )
}
for(i in 1:nc){
  inits[[i]] = GenInits()
}

# Bundle all data in a list for JAGS

jag.data<-list(trait = a.krombeini.pd$Trait, m.N.obs.pd.ind = m.N.obs.pd.ind, temp = a.krombeini.pd$Temp)

a.krombeini.pd.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                           model.file="aedes.individual.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                           n.iter=ni, DIC=T, working.directory=getwd())
a.krombeini.pd.fit.mcmc <- as.mcmc(a.krombeini.pd.fit) ## makes an "mcmc" object
closeAllConnections()
a.krombeini.pd.fit$BUGSoutput$DIC


# Get posterior means
m.summary_fit.a.krombeini.pd <- summary(a.krombeini.pd.fit.mcmc)$statistics
m.a_mean.pd <- m.summary_fit.a.krombeini.pd["a", "Mean"]
m.c_mean.pd <- m.summary_fit.a.krombeini.pd["c", "Mean"]
m.means_Topt.pd <- m.summary_fit.a.krombeini.pd["Topt", "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)

# Get predicted values
a.krombeini.fitted.pd <- function (x){
  m.a_mean.pd * (x - m.means_Topt.pd)^2 + m.c_mean.pd
}

# Data frame of temperature, fitted value, species
m.fitted_data.a.krombeini.pd <- data.frame(Temperature = rep(m.temp_seq),
                                           Fitted = a.krombeini.fitted.pd(m.temp_seq))
# Original data with temperature, trait, and species as numeric
m.plot_data.a.krombeini.pd <- data.frame(
  Temperature = a.krombeini.pd$Temp,
  Trait = a.krombeini.pd$Trait
)

# Plot observed data and fitted curves 
m.p.a.krombeini.pd <- ggplot() +
  geom_point(data = m.plot_data.a.krombeini.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.a.krombeini.pd, aes(x = Temperature, y = Fitted), color = "seagreen4", size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  theme(legend.position = "none")
m.p.a.krombeini.pd

# Merge the observed and predicted data for comparison
m.merged_data.a.krombeini.pd <- m.fitted_data.a.krombeini.pd %>%
  left_join(a.krombeini, by = c("Temperature" = "Temp"))

# Calculate RMSE for each species
m.rmse.a.krombeini.pd <- m.merged_data.a.krombeini.pd %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE),
    SE = (Trait - Fitted)^2
  )

# Extract log-likelihood matrix (samples × observations)
m.log.lik.a.krombeini.pd <- a.krombeini.pd.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
m.waic.a.krombeini.pd <- waic(m.log.lik.a.krombeini.pd)$estimates["waic", "Estimate"]



save(m.rmse.a.krombeini.pd,
     m.fitted_data.a.krombeini.pd,
     m.merged_data.a.krombeini.pd,
     a.krombeini.pd.fit$BUGSoutput$DIC,
     m.waic.a.krombeini.pd,
     m.p.a.krombeini.pd,
     file="m.a.krombeini.pd.fit.RData")
######################################### By Species Individual Notoscriptus #########################################
a.notoscriptus.pd <- mosquito.pd[mosquito.pd$interactor1 == "Aedes notoscriptus", ]

m.N.obs.pd.ind <- length(a.notoscriptus.pd$Trait)

sink("aedes.individual.txt")
cat("
model {
  # Priors
   a ~ dexp(0.1)             # Width
  Topt ~ dnorm(30, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  c ~ dexp(0.001)
  
  # Likelihood
  for (i in 1:m.N.obs.pd.ind) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], tau1) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], tau1)
  }}
", file = "aedes.individual.txt")


# Parameters to monitor
parameters <- c("a", "c", "tau1", "Topt", "log_lik")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Initial Values
inits = vector('list', nc)
GenInits = function() {
  a <- rexp(1, 0.1)
  c <- rexp(1, 0.001)
  tau1 <- rexp(1, 0.01)
  Topt <- rnorm(1, 30, 1/0.001)
  list(
    a = a,             
    c = c,
    tau1 = tau1,
    Topt = Topt
  )
}
for(i in 1:nc){
  inits[[i]] = GenInits()
}

# Bundle all data in a list for JAGS

jag.data<-list(trait = a.notoscriptus.pd$Trait, m.N.obs.pd.ind = m.N.obs.pd.ind, temp = a.notoscriptus.pd$Temp)

a.notoscriptus.pd.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                              model.file="aedes.individual.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                              n.iter=ni, DIC=T, working.directory=getwd())
a.notoscriptus.pd.fit.mcmc <- as.mcmc(a.notoscriptus.pd.fit) ## makes an "mcmc" object
closeAllConnections()
a.notoscriptus.pd.fit$BUGSoutput$DIC


# Get posterior means
m.summary_fit.a.notoscriptus.pd <- summary(a.notoscriptus.pd.fit.mcmc)$statistics
m.a_mean.pd <- m.summary_fit.a.notoscriptus.pd["a", "Mean"]
m.c_mean.pd <- m.summary_fit.a.notoscriptus.pd["c", "Mean"]
m.means_Topt.pd <- m.summary_fit.a.notoscriptus.pd["Topt", "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)

# Get predicted values
a.notoscriptus.fitted.pd <- function (x){
  m.a_mean.pd * (x - m.means_Topt.pd)^2 + m.c_mean.pd
}

# Data frame of temperature, fitted value, species
m.fitted_data.a.notoscriptus.pd <- data.frame(Temperature = rep(m.temp_seq),
                                              Fitted = a.notoscriptus.fitted.pd(m.temp_seq))
# Original data with temperature, trait, and species as numeric
m.plot_data.a.notoscriptus.pd <- data.frame(
  Temperature = a.notoscriptus.pd$Temp,
  Trait = a.notoscriptus.pd$Trait
)

# Plot observed data and fitted curves 
m.p.a.notoscriptus.pd <- ggplot() +
  geom_point(data = m.plot_data.a.notoscriptus.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.a.notoscriptus.pd, aes(x = Temperature, y = Fitted), color = "seagreen4", size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  theme(legend.position = "none")
m.p.a.notoscriptus.pd

# Merge the observed and predicted data for comparison
m.merged_data.a.notoscriptus.pd <- m.fitted_data.a.notoscriptus.pd %>%
  left_join(a.notoscriptus, by = c("Temperature" = "Temp"))

# Calculate RMSE for each species
m.rmse.a.notoscriptus.pd <- m.merged_data.a.notoscriptus.pd %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE),
    SE = (Trait - Fitted)^2
  )

# Extract log-likelihood matrix (samples × observations)
m.log.lik.a.notoscriptus.pd <- a.notoscriptus.pd.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
m.waic.a.notoscriptus.pd <- waic(m.log.lik.a.notoscriptus.pd)$estimates["waic", "Estimate"]


save(m.rmse.a.notoscriptus.pd,
     m.fitted_data.a.notoscriptus.pd,
     m.merged_data.a.notoscriptus.pd,
     a.notoscriptus.pd.fit$BUGSoutput$DIC,
     m.waic.a.notoscriptus.pd,
     m.p.a.notoscriptus.pd,
     file="m.a.notoscriptus.pd.fit.RData")

######################################### By Species Individual Summary #########################################
# Define RMSE values for each species
m.se.individual.pd <- c(m.rmse.a.aegypti.pd$SE, m.rmse.a.albopictus.pd$SE, m.rmse.a.atropalpus.pd$SE,
                        m.rmse.a.camptorhynchus.pd$SE, m.rmse.a.jj.pd$SE,
                        m.rmse.a.krombeini.pd$SE, m.rmse.a.notoscriptus.pd$SE)

# Compute weighted average RMSE
m.overall_rmse.species.ind.pd <- sqrt(mean(m.se.individual.pd, na.rm = TRUE))

m.DIC.species.individual.pd <- a.aegypti.pd.fit$BUGSoutput$DIC+ a.albopictus.pd.fit$BUGSoutput$DIC +
  a.atropalpus.pd.fit$BUGSoutput$DIC + a.camptorhynchus.pd.fit$BUGSoutput$DIC + a.jj.pd.fit$BUGSoutput$DIC +
  a.krombeini.pd.fit$BUGSoutput$DIC + a.notoscriptus.pd.fit$BUGSoutput$DIC

m.waic.species.ind.pd <- m.waic.a.aegypti.pd + m.waic.a.albopictus.pd +
  m.waic.a.atropalpus.pd + m.waic.a.camptorhynchus.pd +
  m.waic.a.jj.pd + m.waic.a.krombeini.pd + m.waic.a.notoscriptus.pd

# Arrange plots into a 3x3 grid
m.grid.ind.species.pd <- grid.arrange(
  m.p.a.aegypti.pd + labs(title = "Aedes aegypti", x = " ", y = "Development Time in Days"),
  m.p.a.albopictus.pd + labs(title = "Aedes albopictus", x = " ", y = " "),
  m.p.a.atropalpus.pd + labs(title = "Aedes atropalpus", x = " ", y = " "),
  m.p.a.camptorhynchus.pd + labs(title = "Aedes camptorhynchus", x = " ", y = "Development Time in Days"),
  m.p.a.jj.pd + labs(title = "Aedes japonicus japonicus", x = " ", y = " "),
  m.p.a.krombeini.pd + labs(title = "Aedes krombeini", x = " ", y = " "),
  m.p.a.notoscriptus.pd + labs(title = "Aedes notoscriptus", y = "Development Time in Days"),
  nrow = 3, ncol = 3
)
ggsave("mosquito.species.ind.plot.png", m.grid.ind.species.pd, width = 8.5, height = 6)

######################################### By Species Hierarchical #########################################
sink("aedes.by.species.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Species-specific 
  mu.c ~ dexp(10)             # Mean for species effects c
  mu.Topt ~ dnorm(30, 0.1)             # Mean for species effects Topt
  tau.Topt ~ dexp(0.01)             # Tau for species effect Topt
  mu.a ~ dexp(0.01)            # Mean for species effects Topt

  for (j in 1:m.n.species.pd) {
    mu.species[j] ~ dexp(mu.c)       # Species-specific low point
    mu.species.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Species-specific temp
    mu.species.a[j] ~ dexp(mu.a)    # Species-specific temp
  }
  # Likelihood
  for (i in 1:m.N.obs.pd) {
    mu[i] <- mu.species.a[m.species.pd[i]] * (mu.species.Topt[m.species.pd[i]] - temp[i])^2 + mu.species[m.species.pd[i]]
    trait[i] ~ dnorm(mu[i], tau1) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], tau1)
  }}
", file = "aedes.by.species.txt")


# Parameters to monitor
parameters <- c("mu.a", "mu.c", "mu.species", "log_lik",
                "tau1", "mu.Topt", "tau.Topt", "mu.species.Topt", "mu.species.a")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Initial Values
inits = vector('list', nc)
GenInits = function() {
  mu.a <- rexp(1, 0.01)
  mu.c <- rexp(1, 10)
  tau1 <- rexp(1, 0.01)
  mu.Topt <- rnorm(1, 30, 1/0.1)
  tau.Topt <- rexp(1, 0.01)
  list(
    mu.a = mu.a,             
    mu.c = mu.c,
    tau1 = tau1,
    mu.Topt = mu.Topt,
    tau.Topt = tau.Topt
  )
}
for(i in 1:nc){
  inits[[i]] = GenInits()
}

# Bundle all data in a list for JAGS
#jag.data<-list(trait = mosquito.data$Trait, N.obs = N.obs, temp = mosquito.data$Temp)

jag.data <- list(trait = mosquito.pd$Trait, m.N.obs.pd = m.N.obs.pd, temp = mosquito.pd$Temp, 
                 m.n.species.pd = m.n.species.pd, m.species.pd = m.species.pd)

m.species.pd.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                      model.file="aedes.by.species.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                      n.iter=ni, DIC=T, working.directory=getwd())
m.species.pd.fit.mcmc <- as.mcmc(m.species.pd.fit) ## makes an "mcmc" object
closeAllConnections()
m.species.pd.fit$BUGSoutput$DIC

# Extract log-likelihood matrix (samples × observations)
m.log.lik.species.pd <- m.species.pd.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
m.waic.species.pd <- waic(m.log.lik.species.pd)$estimates["waic", "Estimate"]

# Get posterior means
m.summary_fit.species.pd <- summary(m.species.pd.fit.mcmc)$statistics
m.a_mean_species.pd <- m.summary_fit.species.pd[grep("^mu.species.a\\[", rownames(m.summary_fit.species.pd)), "Mean"]
m.c_mean_species.pd <- m.summary_fit.species.pd[grep("^mu.species\\[", rownames(m.summary_fit.species.pd)), "Mean"]
m.Topt_mean_species.pd <- m.summary_fit.species.pd[grep("^mu.species.Topt\\[", rownames(m.summary_fit.species.pd)), "Mean"]


# Generate temperature sequence 
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)

# Get predicted values
m.spec1.pd <- function (x){
  m.a_mean_species.pd[1] * (x - m.Topt_mean_species.pd[1])^2 + m.c_mean_species.pd[1]
}
m.spec2.pd <- function (x){
  m.a_mean_species.pd[2] * (x - m.Topt_mean_species.pd[2])^2 + m.c_mean_species.pd[2]
}
m.spec3.pd <- function (x){
  m.a_mean_species.pd[3] * (x - m.Topt_mean_species.pd[3])^2 + m.c_mean_species.pd[3]
}
m.spec4.pd <- function (x){
  m.a_mean_species.pd[4] * (x - m.Topt_mean_species.pd[4])^2 + m.c_mean_species.pd[4]
}
m.spec5.pd <- function (x){
  m.a_mean_species.pd[5] * (x - m.Topt_mean_species.pd[5])^2 + m.c_mean_species.pd[5]
}
m.spec6.pd <- function (x){
  m.a_mean_species.pd[6] * (x - m.Topt_mean_species.pd[6])^2 + m.c_mean_species.pd[6]
}
m.spec7.pd <- function (x){
  m.a_mean_species.pd[7] * (x - m.Topt_mean_species.pd[7])^2 + m.c_mean_species.pd[7]
}

# Optionally, convert to a data frame for easier plotting
m.fitted_data.species.pd <- data.frame(Temperature = rep(m.temp_seq, m.n.species.pd), 
                                    Fitted = c(m.spec1.pd(m.temp_seq),
                                               m.spec2.pd(m.temp_seq),
                                               m.spec3.pd(m.temp_seq),
                                               m.spec4.pd(m.temp_seq),
                                               m.spec5.pd(m.temp_seq),
                                               m.spec6.pd(m.temp_seq),
                                               m.spec7.pd(m.temp_seq)),
                                    Species = factor(rep(c("Aedes aegypti", "Aedes albopictus",
                                                           "Aedes atropalpus", "Aedes camptorhynchus",
                                                           "Aedes japonicus japonicus",
                                                           "Aedes krombeini", "Aedes notoscriptus"), each = length(m.temp_seq))))

m.plot_data.species.pd <- data.frame(
  Temperature = mosquito.pd$Temp,
  Trait = mosquito.pd$Trait,
  Species = as.factor(mosquito.pd$interactor1)
)

# Plot observed data and fitted curves with facet wrap by species
m.p.species.pd <- ggplot() +
  geom_point(data = m.plot_data.species.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.species.pd, aes(x = Temperature, y = Fitted, color = Species), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~ Species) 


# Merge the observed and predicted data for comparison
m.merged_data.species.pd <- m.fitted_data.species.pd %>%
  left_join(m.plot_data.species, by = c("Temperature", "Species"))


# Calculate RMSE for each species
m.mse_by_species.species.pd <- m.merged_data.species.pd %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.species.pd <- m.merged_data.species.pd %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 



save(m.mse_by_species.species.pd,
     m.fitted_data.species.pd,
     m.merged_data.species.pd,
     m.waic.species.pd,
     m.species.pd.fit$BUGSoutput$DIC,
     m.overall_rmse.species.pd,
     m.p.species.pd,
     file="m.species.pd.fit.RData")

######################################### By Species Hierarchical Different Variances #########################################
sink("aedes.by.species.tau.txt")
cat("
model {
  # Priors
  # Species-specific 
  mu.c ~ dexp(10)             # Mean for species effects c
  tau1 ~ dexp(0.01)           # Tau for species effect  
  mu.Topt ~ dnorm(30, 0.1)             # Mean for species effects Topt
  tau.Topt ~ dexp(0.01)             # Tau for species effect Topt
  mu.a ~ dexp(0.01)            # Mean for species effects Topt

  for (j in 1:m.n.species.pd) {
    mu.species[j] ~ dexp(mu.c)       # Species-specific low point
    mu.species.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Species-specific temp
    mu.species.a[j] ~ dexp(mu.a)    # Species-specific temp
    species.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:m.N.obs.pd) {
    mu[i] <- mu.species.a[m.species.pd[i]] * (mu.species.Topt[m.species.pd[i]] - temp[i])^2 + mu.species[m.species.pd[i]]
    trait[i] ~ dnorm(mu[i], species.tau[m.species.pd[i]]) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], species.tau[m.species.pd[i]])
  }}
", file = "aedes.by.species.tau.txt")


# Parameters to monitor
parameters <- c("mu.a", "mu.c", "mu.species", "species.tau", "log_lik",
                "tau1", "mu.Topt", "tau.Topt", "mu.species.Topt", "mu.species.a")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Initial Values
inits = vector('list', nc)
GenInits = function() {
  mu.a <- rexp(1, 0.01)
  mu.c <- rexp(1, 10)
  tau1 <- rexp(1, 0.01)
  mu.Topt <- rnorm(1, 30, 1/0.1)
  tau.Topt <- rexp(1, 0.01)
  list(
    mu.a = mu.a,             
    mu.c = mu.c,
    tau1 = tau1,
    mu.Topt = mu.Topt,
    tau.Topt = tau.Topt
  )
}
for(i in 1:nc){
  inits[[i]] = GenInits()
}

# Bundle all data in a list for JAGS
jag.data <- list(trait = mosquito.pd$Trait, m.N.obs.pd = m.N.obs.pd, temp = mosquito.pd$Temp, 
                 m.n.species.pd = m.n.species.pd, m.species.pd = m.species.pd)

m.species.pd.tau.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                          model.file="aedes.by.species.tau.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                          n.iter=ni, DIC=T, working.directory=getwd())
m.species.pd.tau.fit.mcmc <- as.mcmc(m.species.pd.tau.fit) ## makes an "mcmc" object
closeAllConnections()
m.species.pd.tau.fit$BUGSoutput$DIC

# Extract log-likelihood matrix (samples × observations)
m.log.lik.species.tau.pd <- m.species.pd.tau.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
m.waic.species.tau.pd <- waic(m.log.lik.species.tau.pd)$estimates["waic", "Estimate"]

# Get posterior means
m.summary_fit.species.tau.pd <- summary(m.species.pd.tau.fit.mcmc)$statistics
m.a_mean_species.tau.pd <- m.summary_fit.species.tau.pd[grep("^mu.species.a\\[", rownames(m.summary_fit.species.tau.pd)), "Mean"]
m.c_mean_species.tau.pd <- m.summary_fit.species.tau.pd[grep("^mu.species\\[", rownames(m.summary_fit.species.tau.pd)), "Mean"]
m.Topt_mean_species.tau.pd <- m.summary_fit.species.tau.pd[grep("^mu.species.Topt\\[", rownames(m.summary_fit.species.tau.pd)), "Mean"]


# Generate temperature sequence 
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)

# Get predicted values
m.spec1.tau.pd <- function (x){
  m.a_mean_species.tau.pd[1] * (x - m.Topt_mean_species.tau.pd[1])^2 + m.c_mean_species.tau.pd[1]
}
m.spec2.tau.pd <- function (x){
  m.a_mean_species.tau.pd[2] * (x - m.Topt_mean_species.tau.pd[2])^2 + m.c_mean_species.tau.pd[2]
}
m.spec3.tau.pd <- function (x){
  m.a_mean_species.tau.pd[3] * (x - m.Topt_mean_species.tau.pd[3])^2 + m.c_mean_species.tau.pd[3]
}
m.spec4.tau.pd <- function (x){
  m.a_mean_species.tau.pd[4] * (x - m.Topt_mean_species.tau.pd[4])^2 + m.c_mean_species.tau.pd[4]
}
m.spec5.tau.pd <- function (x){
  m.a_mean_species.tau.pd[5] * (x - m.Topt_mean_species.tau.pd[5])^2 + m.c_mean_species.tau.pd[5]
}
m.spec6.tau.pd <- function (x){
  m.a_mean_species.tau.pd[6] * (x - m.Topt_mean_species.tau.pd[6])^2 + m.c_mean_species.tau.pd[6]
}
m.spec7.tau.pd <- function (x){
  m.a_mean_species.tau.pd[7] * (x - m.Topt_mean_species.tau.pd[7])^2 + m.c_mean_species.tau.pd[7]
}

# Optionally, convert to a data frame for easier plotting
m.fitted_data.species.tau.pd <- data.frame(Temperature = rep(m.temp_seq, m.n.species.pd), 
                                        Fitted = c(m.spec1.tau.pd(m.temp_seq),
                                                   m.spec2.tau.pd(m.temp_seq),
                                                   m.spec3.tau.pd(m.temp_seq),
                                                   m.spec4.tau.pd(m.temp_seq),
                                                   m.spec5.tau.pd(m.temp_seq),
                                                   m.spec6.tau.pd(m.temp_seq),
                                                   m.spec7.tau.pd(m.temp_seq)),
                                        Species = factor(rep(c("Aedes aegypti", "Aedes albopictus",
                                                               "Aedes atropalpus", "Aedes camptorhynchus",
                                                               "Aedes japonicus japonicus",
                                                               "Aedes krombeini", "Aedes notoscriptus"), each = length(m.temp_seq))))

m.plot_data.species.pd <- data.frame(
  Temperature = mosquito.pd$Temp,
  Trait = mosquito.pd$Trait,
  Species = as.factor(mosquito.pd$interactor1)
)

# Plot observed data and fitted curves with facet wrap by species
m.p.species.tau.pd <- ggplot() +
  geom_point(data = m.plot_data.species.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.species.tau.pd, aes(x = Temperature, y = Fitted, color = Species), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~ Species) 
#ggsave("mosquito.species.tau.plot.png", width = 8, height = 5)

# Merge the observed and predicted data for comparison
m.merged_data.species.tau.pd <- m.fitted_data.species.tau.pd %>%
  left_join(m.plot_data.species, by = c("Temperature", "Species"))


# Calculate RMSE for each species
m.mse_by_species.species.tau.pd <- m.merged_data.species.tau.pd %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.species.tau.pd <- m.merged_data.species.tau.pd %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(m.mse_by_species.species.tau.pd,
     m.fitted_data.species.tau.pd,
     m.merged_data.species.tau.pd,
     m.waic.species.tau.pd,
     m.species.pd.tau.fit$BUGSouput$DIC,
     m.overall_rmse.species.tau.pd,
     m.p.species.tau.pd,
     file="m.species.pd.tau.fit.RData")

######################################### By Climate #########################################
m.climate.pd <- as.factor(mosquito.pd$climate)
m.n.climate.pd <- length(unique(mosquito.pd$climate))
sink("aedes.by.climate.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Climate-specific 
  mu.c ~ dexp(10)             # Mean for climate effects c
  mu.Topt ~ dnorm(30, 0.1)             # Mean for climate effects Topt
  tau.Topt ~ dexp(0.01)             # Tau for climate effect Topt
  mu.a ~ dexp(0.01)            # Mean for climate effects Topt

  for (j in 1:m.n.climate.pd) {
    mu.climate[j] ~ dexp(mu.c)       # Climate-specific low point
    mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific temp
    mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific temp
    climate.tau[j] ~ dexp(tau1)
  }

  # Likelihood
  for (i in 1:m.N.obs.pd) {
    mu[i] <- mu.climate.a[m.climate.pd[i]] * (mu.climate.Topt[m.climate.pd[i]] - temp[i])^2 + mu.climate[m.climate.pd[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[m.climate.pd[i]]) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], climate.tau[m.climate.pd[i]])
  }}
", file = "aedes.by.climate.txt")

# Parameters to monitor
parameters <- c("mu.a", "mu.c", "mu.climate", "climate.tau", "log_lik",
                "tau1", "mu.Topt", "tau.Topt", "mu.climate.Topt", "mu.climate.a")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Initial Values
inits = vector('list', nc)
GenInits = function() {
  mu.a <- rexp(1, 0.01)
  mu.c <- rexp(1, 10)
  tau1 <- rexp(1, 0.01)
  mu.Topt <- rnorm(1, 30, 1/0.1)
  tau.Topt <- rexp(1, 0.01)
  list(
    mu.a = mu.a,             
    mu.c = mu.c,
    tau1 = tau1,
    mu.Topt = mu.Topt,
    tau.Topt = tau.Topt
  )
}
for(i in 1:nc){
  inits[[i]] = GenInits()
}

# Bundle all data in a list for JAGS
jag.data<-list(trait = mosquito.pd$Trait, m.N.obs.pd = m.N.obs.pd, temp = mosquito.pd$Temp, 
               m.n.climate.pd = m.n.climate.pd, m.climate.pd = m.climate.pd)

m.climate.pd.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                      model.file="aedes.by.climate.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                      n.iter=ni, DIC=T, working.directory=getwd())
m.climate.pd.fit.mcmc <- as.mcmc(m.climate.pd.fit) ## makes an "mcmc" object
closeAllConnections()
m.climate.pd.fit$BUGSoutput$DIC

# Extract log-likelihood matrix (samples × observations)
m.log.lik.climate.pd <- m.climate.pd.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
m.waic.climate.pd <- waic(m.log.lik.climate.pd)$estimates["waic", "Estimate"]

# Get posterior means
m.summary_fit.climate.pd <- summary(m.climate.pd.fit.mcmc)$statistics
m.a_mean_climate.pd <- m.summary_fit.climate.pd[grep("^mu.climate.a\\[", rownames(m.summary_fit.climate.pd)), "Mean"]
m.c_mean_climate.pd <- m.summary_fit.climate.pd[grep("^mu.climate\\[", rownames(m.summary_fit.climate.pd)), "Mean"]
m.Topt_mean_climate.pd <- m.summary_fit.climate.pd[grep("^mu.climate.Topt\\[", rownames(m.summary_fit.climate.pd)), "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)

# Get predicted values
m.clim1.pd <- function (x){
  m.a_mean_climate.pd[1] * (x - m.Topt_mean_climate.pd[1])^2 + m.c_mean_climate.pd[1]
}
m.clim2.pd <- function (x){
  m.a_mean_climate.pd[2] * (x - m.Topt_mean_climate.pd[2])^2 + m.c_mean_climate.pd[2]
}
m.clim3.pd <- function (x){
  m.a_mean_climate.pd[3] * (x - m.Topt_mean_climate.pd[3])^2 + m.c_mean_climate.pd[3]
}


# Optionally, convert to a data frame for easier plotting
m.fitted_data.climate.pd <- data.frame(Temperature = rep(m.temp_seq, m.n.climate), 
                                    Fitted = c(m.clim1.pd(m.temp_seq),
                                               m.clim2.pd(m.temp_seq),
                                               m.clim3.pd(m.temp_seq)),
                                    Climate = factor(rep(c("subtropical", "temperate", "tropical"), each = length(m.temp_seq))))

m.plot_data.climate.pd <- data.frame(
  Temperature = mosquito.pd$Temp,
  Trait = mosquito.pd$Trait,
  Climate = as.factor(mosquito.pd$climate),
  Species = as.factor(mosquito.pd$interactor1)
)

m.species.pd_to_climate <- unique(m.plot_data.climate.pd[, c("Species", "Climate")])  # Ensure species matches correct climate
m.species.pd_to_climate$Climate <- as.factor(m.species.pd_to_climate$Climate)
# Merge fitted climate curves with species-to-climate mapping
m.fitted_data.climate.pd <- data.frame(
  Temperature = rep(m.temp_seq, m.n.climate),
  Fitted = c(m.clim1.pd(m.temp_seq),
             m.clim2.pd(m.temp_seq),
             m.clim3.pd(m.temp_seq)),
  Climate = factor(rep(c("subtropical", "temperate", "tropical"), each = length(m.temp_seq)))
) %>%
  left_join(m.species.pd_to_climate, by = "Climate")

m.plot_data.climate.pd$Climate <- factor(m.plot_data.climate.pd$Climate, levels = c("temperate", "subtropical", "tropical"))
m.fitted_data.climate.pd$Climate <- factor(m.fitted_data.climate.pd$Climate, levels = c("temperate", "subtropical", "tropical"))
# Plot observed data and fitted curves with facet wrap by species
m.p.climate.pd <- ggplot() +
  geom_point(data = m.plot_data.climate.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.climate.pd, aes(x = Temperature, y = Fitted, color = Climate), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) 

# Merge the observed and predicted data for comparison
m.merged_data.climate.pd <- m.fitted_data.climate.pd %>%
  left_join(m.plot_data.climate, by = c("Temperature", "Climate", "Species"))


# Calculate RMSE for each climate
m.mse_by_species.climate.pd <- m.merged_data.climate.pd |> 
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.climate.pd <- m.merged_data.climate.pd |> 
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 



save(m.mse_by_species.climate.pd,
     m.fitted_data.climate.pd,
     m.merged_data.climate.pd,
     m.p.climate.pd,
     m.climate.pd.fit$BUGSouput$DIC,
     m.waic.climate.pd,
     m.overall_rmse.climate.pd,
     file="m.climate.pd.fit.RData")


######################################### By Wing Size Category #########################################
m.size.pd <- as.factor(mosquito.pd$wing.size)
m.n.size.pd <- length(unique(mosquito.pd$wing.size))

sink("aedes.by.size.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Size-specific 
  mu.c ~ dexp(10)             # Mean for size effects c
  mu.Topt ~ dnorm(30, 0.1)             # Mean for size effects Topt
  tau.Topt ~ dexp(0.01)             # Tau for size effect Topt
  mu.a ~ dexp(0.01)            # Mean for size effects Topt

  for (j in 1:m.n.size.pd) {
    mu.size[j] ~ dexp(mu.c)       # Size-specific low point
    mu.size.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Size-specific temp
    mu.size.a[j] ~ dexp(mu.a)    # Size-specific temp
    size.tau[j] ~ dexp(tau1)
  }

  # Likelihood
  for (i in 1:m.N.obs.pd) {
    mu[i] <- mu.size.a[m.size.pd[i]] * (mu.size.Topt[m.size.pd[i]] - temp[i])^2 + mu.size[m.size.pd[i]]
    trait[i] ~ dnorm(mu[i], size.tau[m.size.pd[i]]) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], size.tau[m.size.pd[i]])
  }}
", file = "aedes.by.size.txt")

# Parameters to monitor
parameters <- c("mu.a", "mu.c", "mu.size", "size.tau", "log_lik",
                "tau1", "mu.Topt", "tau.Topt", "mu.size.Topt", "mu.size.a")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Initial Values
inits = vector('list', nc)
GenInits = function() {
  mu.a <- rexp(1, 0.01)
  mu.c <- rexp(1, 10)
  tau1 <- rexp(1, 0.01)
  mu.Topt <- rnorm(1, 30, 1/0.1)
  tau.Topt <- rexp(1, 0.01)
  list(
    mu.a = mu.a,             
    mu.c = mu.c,
    tau1 = tau1,
    mu.Topt = mu.Topt,
    tau.Topt = tau.Topt
  )
}
for(i in 1:nc){
  inits[[i]] = GenInits()
}

# Bundle all data in a list for JAGS
jag.data<-list(trait = mosquito.pd$Trait, m.N.obs.pd = m.N.obs.pd, temp = mosquito.pd$Temp, 
               m.n.size.pd = m.n.size.pd, m.size.pd = m.size.pd)

m.size.pd.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                   model.file="aedes.by.size.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                   n.iter=ni, DIC=T, working.directory=getwd())
m.size.pd.fit.mcmc <- as.mcmc(m.size.pd.fit) ## makes an "mcmc" object
closeAllConnections()
m.size.pd.fit$BUGSoutput$DIC

# Extract log-likelihood matrix (samples × observations)
m.log.lik.size.pd <- m.size.pd.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
m.waic.size.pd <- waic(m.log.lik.size.pd)$estimates["waic", "Estimate"]

# Get posterior means
m.summary_fit.size.pd <- summary(m.size.pd.fit.mcmc)$statistics
m.a_mean_size.pd <- m.summary_fit.size.pd[grep("^mu.size.a\\[", rownames(m.summary_fit.size.pd)), "Mean"]
m.c_mean_size.pd <- m.summary_fit.size.pd[grep("^mu.size\\[", rownames(m.summary_fit.size.pd)), "Mean"]
m.Topt_mean_size.pd <- m.summary_fit.size.pd[grep("^mu.size.Topt\\[", rownames(m.summary_fit.size.pd)), "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)

# Get predicted values
m.size.pd1 <- function (x){
  m.a_mean_size.pd[1] * (x - m.Topt_mean_size.pd[1])^2 + m.c_mean_size.pd[1]
}
m.size.pd2 <- function (x){
  m.a_mean_size.pd[2] * (x - m.Topt_mean_size.pd[2])^2 + m.c_mean_size.pd[2]
}
m.size.pd3 <- function (x){
  m.a_mean_size.pd[3] * (x - m.Topt_mean_size.pd[3])^2 + m.c_mean_size.pd[3]
}


# Optionally, convert to a data frame for easier plotting
m.fitted_data.size.pd <- data.frame(Temperature = rep(m.temp_seq, m.n.size.pd), 
                                 Fitted = c(m.size.pd1(m.temp_seq),
                                            m.size.pd2(m.temp_seq),
                                            m.size.pd3(m.temp_seq)),
                                 Size = factor(rep(c("large", "medium", "small"), each = length(m.temp_seq))))

m.plot_data.size.pd <- data.frame(
  Temperature = mosquito.pd$Temp,
  Trait = mosquito.pd$Trait,
  Size = as.factor(mosquito.pd$wing.size),
  Species = as.factor(mosquito.pd$interactor1)
)

m.species.pd_to_size <- unique(m.plot_data.size.pd[, c("Species", "Size")])  # Ensure species matches correct climate
m.species.pd_to_size$Size <- as.factor(m.species.pd_to_size$Size)
# Merge fitted climate curves with species-to-climate mapping
m.fitted_data.size.pd <- data.frame(
  Temperature = rep(m.temp_seq, m.n.size.pd),
  Fitted = c(m.size.pd1(m.temp_seq),
             m.size.pd2(m.temp_seq),
             m.size.pd3(m.temp_seq)),
  Size = factor(rep(c("large", "medium", "small"), each = length(m.temp_seq)))
) %>%
  left_join(m.species.pd_to_size, by = "Size")

m.plot_data.size.pd$Size <- factor(m.plot_data.size.pd$Size, levels = c("small", "medium", "large"))
m.fitted_data.size.pd$Size <- factor(m.fitted_data.size.pd$Size, levels = c("small", "medium", "large"))
# Plot observed data and fitted curves with facet wrap by species
m.p.size.pd <- ggplot() +
  geom_point(data = m.plot_data.size.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.size.pd, aes(x = Temperature, y = Fitted, color = Size), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) 


# Merge the observed and predicted data for comparison
m.merged_data.size.pd <- m.fitted_data.size.pd %>%
  left_join(m.plot_data.size, by = c("Temperature", "Size", "Species"))


# Calculate RMSE for each species
m.mse_by_species.size.pd <- m.merged_data.size.pd %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.size.pd <- m.merged_data.size.pd %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 



save(m.mse_by_species.size.pd,
     m.fitted_data.size.pd,
     m.merged_data.size.pd,
     m.p.size.pd,
     m.size.pd.fit$BUGSoutput$DIC,
     m.waic.size.pd,
     m.overall_rmse.size.pd,
     file="m.size.pd.fit.RData")

######################################### By Wing Size Category Removed #########################################
mosquito.pd1 <- mosquito.pd[!(mosquito.pd$interactor1 %in% c("Aedes notoscriptus", "Aedes krombeini")), ]
m.size.pd <- as.factor(mosquito.pd1$wing.size)
m.n.size.pd <- length(unique(mosquito.pd1$wing.size))
m.N.obs.pd.rm <- length(mosquito.pd1$Trait)

sink("aedes.by.size.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Size-specific 
  mu.c ~ dexp(10)             # Mean for size effects c
  mu.Topt ~ dnorm(30, 0.1)             # Mean for size effects Topt
  tau.Topt ~ dexp(0.01)             # Tau for size effect Topt
  mu.a ~ dexp(0.01)            # Mean for size effects Topt

  for (j in 1:m.n.size.pd) {
    mu.size[j] ~ dexp(mu.c)       # Size-specific low point
    mu.size.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Size-specific temp
    mu.size.a[j] ~ dexp(mu.a)    # Size-specific temp
    size.tau[j] ~ dexp(tau1)
  }

  # Likelihood
  for (i in 1:m.N.obs.pd.rm) {
    mu[i] <- mu.size.a[m.size.pd[i]] * (mu.size.Topt[m.size.pd[i]] - temp[i])^2 + mu.size[m.size.pd[i]]
    trait[i] ~ dnorm(mu[i], size.tau[m.size.pd[i]]) T(0, ) 
  }}
", file = "aedes.by.size.txt")

# Parameters to monitor
parameters <- c("mu.a", "mu.c", "mu.size", "size.tau",
                "tau1", "mu.Topt", "tau.Topt", "mu.size.Topt", "mu.size.a")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Initial Values
inits = vector('list', nc)
GenInits = function() {
  mu.a <- rexp(1, 0.01)
  mu.c <- rexp(1, 10)
  tau1 <- rexp(1, 0.01)
  mu.Topt <- rnorm(1, 30, 1/0.1)
  tau.Topt <- rexp(1, 0.01)
  list(
    mu.a = mu.a,             
    mu.c = mu.c,
    tau1 = tau1,
    mu.Topt = mu.Topt,
    tau.Topt = tau.Topt
  )
}
for(i in 1:nc){
  inits[[i]] = GenInits()
}

# Bundle all data in a list for JAGS
jag.data<-list(trait = mosquito.pd1$Trait, m.N.obs.pd.rm = m.N.obs.pd.rm, temp = mosquito.pd1$Temp, 
               m.n.size.pd = m.n.size.pd, m.size.pd = m.size.pd)

m.size.pd.rm.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                      model.file="aedes.by.size.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                      n.iter=ni, DIC=T, working.directory=getwd())
m.size.pd.rm.fit.mcmc <- as.mcmc(m.size.pd.rm.fit) ## makes an "mcmc" object
closeAllConnections()
m.size.pd.rm.fit$BUGSoutput$DIC

# Get posterior means
m.summary_fit.size.rm.pd <- summary(m.size.pd.rm.fit.mcmc)$statistics
m.a_mean_size.rm.pd <- m.summary_fit.size.rm.pd[grep("^mu.size.a\\[", rownames(m.summary_fit.size.rm.pd)), "Mean"]
m.c_mean_size.rm.pd <- m.summary_fit.size.rm.pd[grep("^mu.size\\[", rownames(m.summary_fit.size.rm.pd)), "Mean"]
m.Topt_mean_size.rm.pd <- m.summary_fit.size.rm.pd[grep("^mu.size.Topt\\[", rownames(m.summary_fit.size.rm.pd)), "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)

# Get predicted values
m.size.pd1.rm <- function (x){
  m.a_mean_size.rm.pd[1] * (x - m.Topt_mean_size.rm.pd[1])^2 + m.c_mean_size.rm.pd[1]
}
m.size.pd2.rm <- function (x){
  m.a_mean_size.rm.pd[2] * (x - m.Topt_mean_size.rm.pd[2])^2 + m.c_mean_size.rm.pd[2]
}
m.size.pd3.rm <- function (x){
  m.a_mean_size.rm.pd[3] * (x - m.Topt_mean_size.rm.pd[3])^2 + m.c_mean_size.rm.pd[3]
}


# Optionally, convert to a data frame for easier plotting
m.fitted_data.size.rm.pd <- data.frame(Temperature = rep(m.temp_seq, m.n.size.pd), 
                                    Fitted = c(m.size.pd1.rm(m.temp_seq),
                                               m.size.pd2.rm(m.temp_seq),
                                               m.size.pd3.rm(m.temp_seq)),
                                    Size = factor(rep(c("large", "medium", "small"), each = length(m.temp_seq))))

m.plot_data.size.rm.pd <- data.frame(
  Temperature = mosquito.pd$Temp,
  Trait = mosquito.pd$Trait,
  Size = as.factor(mosquito.pd$wing.size),
  Species = as.factor(mosquito.pd$interactor1)
)

m.species.pd_to_size.rm <- unique(m.plot_data.size.pd[, c("Species", "Size")])  # Ensure species matches correct size
m.species.pd_to_size.rm$Size <- as.factor(m.species.pd_to_size.rm$Size)
# Merge fitted size curves with species-to-size mapping
m.fitted_data.size.rm.pd <- data.frame(
  Temperature = rep(m.temp_seq, m.n.size.pd),
  Fitted = c(m.size.pd1.rm(m.temp_seq),
             m.size.pd2.rm(m.temp_seq),
             m.size.pd3.rm(m.temp_seq)),
  Size = factor(rep(c("large", "medium", "small"), each = length(m.temp_seq)))
) %>%
  left_join(m.species.pd_to_size.rm, by = "Size")

m.plot_data.size.rm.pd$Size <- factor(m.plot_data.size.rm.pd$Size, levels = c("small", "medium", "large"))
m.fitted_data.size.rm.pd$Size <- factor(m.fitted_data.size.rm.pd$Size, levels = c("small", "medium", "large"))
# Plot observed data and fitted curves with facet wrap by species
m.p.size.rm.pd <- ggplot() +
  geom_point(data = m.plot_data.size.rm.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.size.rm.pd, aes(x = Temperature, y = Fitted, color = Size), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) 


# Merge the observed and predicted data for comparison
m.merged_data.size.rm.pd <- m.fitted_data.size.rm.pd %>%
  left_join(m.plot_data.size.rm, by = c("Temperature", "Size", "Species"))


# Calculate RMSE for each species
m.mse_by_species.size.rm.pd <- m.merged_data.size.rm.pd %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.size.rm.pd <- m.merged_data.size.rm.pd %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 



save(m.mse_by_species.size.rm.pd,
     m.fitted_data.size.rm.pd,
     m.merged_data.size.rm.pd,
     #m.size.pd.rm.fit.mcmc,
     #m.size.pd.rm.fit,
     m.p.size.rm.pd,
     file="m.size.pd.rm.fit.RData")


######################################### By Wing Size and Climate #########################################
m.climate.pd <- as.factor(mosquito.pd$climate)
m.n.climate.pd <- length(unique(mosquito.pd$climate))
m.size.pd <- as.factor(mosquito.pd$wing.size)
m.n.size.pd <- length(unique(mosquito.pd$wing.size))
sink("aedes.by.climate.size.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Size-specific 
  mu.c ~ dexp(10)             # Mean for size effects c
  # Climate-specific
  mu.Topt ~ dnorm(30, 0.1)             # Mean for climate effects Topt
  tau.Topt ~ dexp(0.01)             # Tau for species climate Topt
  mu.a ~ dexp(0.01)            # Mean for climate effects Topt

  for (j in 1:m.n.climate.pd) {
    mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific temp
    mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific temp
  }
    for (k in 1:m.n.size.pd) {
    mu.size[k] ~ dexp(mu.c)       # Size-specific low point
    size.tau[k] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:m.N.obs.pd) {
    mu[i] <- mu.climate.a[m.climate.pd[i]] * (mu.climate.Topt[m.climate.pd[i]] - temp[i])^2 + mu.size[m.size.pd[i]]
    trait[i] ~ dnorm(mu[i], size.tau[m.size.pd[i]])  T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], size.tau[m.size.pd[i]])
  }}
", file = "aedes.by.climate.size.txt")


# Parameters to monitor
parameters <- c("mu.a", "mu.c", "mu.size", "size.tau", "log_lik",
                "tau1", "mu.Topt", "tau.Topt", "mu.climate.Topt", "mu.climate.a")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Initial Values
inits = vector('list', nc)
GenInits = function() {
  mu.a <- rexp(1, 0.01)
  mu.c <- rexp(1, 10)
  tau1 <- rexp(1, 0.01)
  mu.Topt <- rnorm(1, 30, 1/0.1)
  tau.Topt <- rexp(1, 0.01)
  list(
    mu.a = mu.a,             
    mu.c = mu.c,
    tau1 = tau1,
    mu.Topt = mu.Topt,
    tau.Topt = tau.Topt
  )
}
for(i in 1:nc){
  inits[[i]] = GenInits()
}

# Bundle all data in a list for JAGS
jag.data<-list(trait = mosquito.pd$Trait, m.N.obs.pd = m.N.obs.pd, temp = mosquito.pd$Temp, 
               m.n.climate.pd = m.n.climate.pd, m.climate.pd = m.climate.pd, m.n.size.pd = m.n.size.pd, m.size.pd = m.size.pd)

m.climate.pd.size.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                           model.file="aedes.by.climate.size.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                           n.iter=ni, DIC=T, working.directory=getwd())
m.climate.pd.size.fit.mcmc <- as.mcmc(m.climate.pd.size.fit) ## makes an "mcmc" object
closeAllConnections()
m.climate.pd.size.fit$BUGSoutput$DIC

# Extract log-likelihood matrix (samples × observations)
m.log.lik.climate.size.pd <- m.climate.pd.size.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
m.waic.climate.size.pd <- waic(m.log.lik.climate.size.pd)$estimates["waic", "Estimate"]

# Get posterior means
m.summary_fit.climate.size.pd <- summary(m.climate.pd.size.fit.mcmc)$statistics
m.a_mean_climate.size.pd <- m.summary_fit.climate.size.pd[grep("^mu.climate.a\\[", rownames(m.summary_fit.climate.size.pd)), "Mean"]
m.c_mean_climate.size.pd <- m.summary_fit.climate.size.pd[grep("^mu.size\\[", rownames(m.summary_fit.climate.size.pd)), "Mean"]
m.Topt_mean_climate.size.pd <- m.summary_fit.climate.size.pd[grep("^mu.climate.Topt\\[", rownames(m.summary_fit.climate.size.pd)), "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)

# Get predicted values
m.clim1size1.pd <- function (x){
  m.a_mean_climate.size.pd[1] * (x - m.Topt_mean_climate.size.pd[1])^2 + m.c_mean_climate.size.pd[1]
}
m.clim2size1.pd <- function (x){
  m.a_mean_climate.size.pd[2] * (x - m.Topt_mean_climate.size.pd[2])^2 + m.c_mean_climate.size.pd[1]
}
m.clim3size1.pd <- function (x){
  m.a_mean_climate.size.pd[3] * (x - m.Topt_mean_climate.size.pd[3])^2 + m.c_mean_climate.size.pd[1]
}
m.clim1size2.pd <- function (x){
  m.a_mean_climate.size.pd[1] * (x - m.Topt_mean_climate.size.pd[1])^2 + m.c_mean_climate.size.pd[2]
}
m.clim2size2.pd <- function (x){
  m.a_mean_climate.size.pd[2] * (x - m.Topt_mean_climate.size.pd[2])^2 + m.c_mean_climate.size.pd[2]
}
m.clim3size2.pd <- function (x){
  m.a_mean_climate.size.pd[3] * (x - m.Topt_mean_climate.size.pd[3])^2 + m.c_mean_climate.size.pd[2]
}
m.clim1size3.pd <- function (x){
  m.a_mean_climate.size.pd[1] * (x - m.Topt_mean_climate.size.pd[1])^2 + m.c_mean_climate.size.pd[3]
}
m.clim2size3.pd <- function (x){
  m.a_mean_climate.size.pd[2] * (x - m.Topt_mean_climate.size.pd[2])^2 + m.c_mean_climate.size.pd[3]
}
m.clim3size3.pd <- function (x){
  m.a_mean_climate.size.pd[3] * (x - m.Topt_mean_climate.size.pd[3])^2 + m.c_mean_climate.size.pd[3]
}

# Optionally, convert to a data frame for easier plotting
m.fitted_data.climate.size.pd <- data.frame(Temperature = rep(m.temp_seq, m.n.climate.pd * m.n.size.pd), 
                                         Fitted = c(m.clim1size1.pd(m.temp_seq),
                                                    m.clim2size1.pd(m.temp_seq),
                                                    m.clim3size1.pd(m.temp_seq),
                                                    m.clim1size2.pd(m.temp_seq),
                                                    m.clim2size2.pd(m.temp_seq),
                                                    m.clim3size2.pd(m.temp_seq),
                                                    m.clim1size3.pd(m.temp_seq),
                                                    m.clim2size3.pd(m.temp_seq),
                                                    m.clim3size3.pd(m.temp_seq)),
                                         Climate = factor(rep(rep(c("subtropical", "temperate", "tropical"), each = length(m.temp_seq)), 3)),
                                         Size = factor(rep(c("large", "medium", "small"), each = length(m.temp_seq) * 3)))

m.plot_data.climate.size.pd <- data.frame(
  Temperature = mosquito.pd$Temp,
  Trait = mosquito.pd$Trait,
  Climate = as.factor(mosquito.pd$climate),
  Species = as.factor(mosquito.pd$interactor1),
  Size = as.factor(mosquito.pd$wing.size)
)

m.species.pd_to_climate.size <- unique(m.plot_data.climate.siz.pde[, c("Species", "Climate", "Size")])
m.species.pd_to_climate.size$Climate <- as.factor(m.species.pd_to_climate.size$Climate)
m.species.pd_to_climate.size$Size <- as.factor(m.species.pd_to_climate.size$Size)

# Merge fitted climate curves with species-to-climate mapping and filter for matching species
m.fitted_data.climate.size.pd <- data.frame(
  Temperature = rep(m.temp_seq, m.n.climate.pd * m.n.size.pd), 
  Fitted = c(m.clim1size1.pd(m.temp_seq),
             m.clim2size1.pd(m.temp_seq),
             m.clim3size1.pd(m.temp_seq),
             m.clim1size2.pd(m.temp_seq),
             m.clim2size2.pd(m.temp_seq),
             m.clim3size2.pd(m.temp_seq),
             m.clim1size3.pd(m.temp_seq),
             m.clim2size3.pd(m.temp_seq),
             m.clim3size3.pd(m.temp_seq)),
  Climate = factor(rep(rep(c("subtropical", "temperate", "tropical"), each = length(m.temp_seq)), 3)),
  Size = factor(rep(c("large", "medium", "small"), each = length(m.temp_seq) * 3))) |> 
  left_join(m.species.pd_to_climate.size, by = c("Climate", "Size")) |> 
  filter(!is.na(Species))  # Keep only rows where Species is non-missing

# Combine Climate and Size for coloring
m.fitted_data.climate.size.pd$Combo <- with(m.fitted_data.climate.size.pd, paste(Climate, Size, sep = ", "))
m.fitted_data.climate.size.pd$Combo <- factor(m.fitted_data.climate.size.pd$Combo, levels = c("temperate, large", 
                                                                                        "subtropical, medium", "subtropical, large",
                                                                                        "tropical, small", "tropical, medium"))
# Plot observed data and fitted curves with facet wrap by species
m.p.climate.size.pd <- ggplot() +
  geom_point(data = m.plot_data.climate.size.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.climate.size.pd, aes(x = Temperature, y = Fitted, color = Combo), size = 1) +
  labs(
    x = "Temperature",
    y = "Development Time in Days",
    color = "Climate Size Combo"
  ) +
  theme_minimal() +
  facet_wrap(~ Species)
#ggsave("mosquito.climate.size.plot.png", width = 8, height = 5)

# Merge the observed and predicted data for comparison
m.merged_data.climate.size.pd <- m.fitted_data.climate.size.pd %>%
  left_join(m.plot_data.climate.size, by = c("Temperature", "Climate", "Species"))


# Calculate RMSE for each species
m.mse_by_species.climate.size.pd <- m.merged_data.climate.size.pd %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.climate.size.pd <- m.merged_data.climate.size.pd %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 



save(m.mse_by_species.climate.size.pd,
     m.fitted_data.climate.size.pd,
     m.merged_data.climate.size.pd,
     m.p.climate.size.pd,
     m.climate.pd.size.fit$BUGSoutput$DIC,
     m.waic.climate.size.pd,
     m.overall_rmse.climate.size.pd,
     file="m.climate.pd.size.fit.RData")

######################################### RMSE, DIC #########################################
# df.aedes.pd <- data.frame(
#   Method = c("Full Data", "By Species Individual", "By Species Hierarchical", 
#              "By Species Hierarchical, Different Variances",
#              "By Climate", "By Wing Size",
#              "By Wing Size and Climate"),
#   DIC = c(m.one.data.pd.fit$BUGSoutput$DIC, m.DIC.species.individual.pd, 
#           m.species.pd.fit$BUGSoutput$DIC,
#           m.species.pd.tau.fit$BUGSoutput$DIC, m.climate.pd.fit$BUGSoutput$DIC,
#           m.size.pd.fit$BUGSoutput$DIC, 
#           m.climate.pd.size.fit$BUGSoutput$DIC),
#   wAIC = c(m.waic.full.pd, m.waic.species.ind.pd, 
#            m.waic.species.pd, m.waic.species.tau.pd,
#            m.waic.climate.pd, m.waic.size.pd,
#            m.waic.climate.size.pd),
#   Total.RMSE = as.numeric(c(m.overall_rmse.pd, m.overall_rmse.species.ind.pd,
#                             m.overall_rmse.species.pd,
#                             m.overall_rmse.species.tau.pd, m.overall_rmse.climate.pd,
#                             m.overall_rmse.size.pd,
#                             m.overall_rmse.climate.size.pd))
# )
# 
# df.aedes.pd
# 
# df.aedes.pd$Method <- reorder(df.aedes.pd$Method, -df.aedes.pd$Total.RMSE)
# ggplot(df.aedes.pd, aes(x = Method, y = Total.RMSE, fill = Method)) +
#   geom_bar(stat = "identity", position = "dodge", color = "black") +
#   labs(title = "RMSE by Method", x = "Method",
#        y = "Root Mean Square Error (RMSE)") +
#   coord_cartesian(ylim = c(3.5, 3.9)) +
#   theme_minimal() +
#   scale_fill_viridis_d(option = "magma") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         legend.position = "none")
# ggsave("mosquito.RMSE.png", width = 8, height = 5)
# 
# 
# df.aedes.pd$Method <- reorder(df.aedes.pd$Method, -df.aedes.pd$DIC)
# ggplot(df.aedes.pd, aes(x = Method, y = DIC, fill = Method)) +
#   geom_bar(stat = "identity", position = "dodge", color = "black") +
#   labs(title = "DIC by Method", x = "Method",
#        y = "Deviance Information Criterion (DIC)") +
#   theme_minimal() +
#   coord_cartesian(ylim = c(5700, 6100)) +
#   scale_fill_viridis_d(option = "magma") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         legend.position = "none")
# ggsave("mosquito.DIC.png", width = 8, height = 5)
# 
# save(df.aedes.pd, file = "Dataset.Mosquito.RData")


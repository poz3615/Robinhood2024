library(dplyr)
library(tidyr)
library(ggplot2)
library(rjags)
library(R2jags)
library(readxl)
library(truncnorm)
library(scoringRules)


######################################### Data Cleaning #########################################
# Mosquito development rate in days, a lot small so filtered dataset to >1 to try to get the hang of it
# Read data, get correct columns, correct metric
test <- read.csv("development-rate-observations.csv")
test <- test[,c("Trait", "originaltraitunit","interactor1", "Temp")]
test$Trait <- ifelse(test$originaltraitunit == "hours", test$Trait/24, test$Trait)
test$Trait <- ifelse(test$originaltraitunit %in% c("days-1", "day-1"), 1/test$Trait, test$Trait)
#test <- test[test$originaltraitunit %in% c("days","day"),]
test <- test[test$Trait > 0.1 & test$Trait < 100, ]
N.obs <- length(test$Trait)

library(dplyr)
# Small- 0-2.55
# Medium- 2.56-3.05
# Large - 3.06-Inf
test <- test %>%
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

test <- test %>%
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

test$climate <- as.factor(test$climate)
test$wing.size <- as.factor(test$wing.size)
# Plotted for each species
ggplot(test, aes(y = Trait, x = Temp, col = climate)) +
  geom_point() +
  facet_wrap(~ interactor1) +
  labs(title = "Temperature by Trait for Each Factor Level",
       x = "Temperature",
       y = "Oviposition") +
  theme_minimal()


# Making vector for each species 
test$interactor1 <- as.factor(test$interactor1)
species <- as.numeric(test$interactor1)
n.species <- length(unique(test$interactor1))

# Model!
# Used an exponential prior on a to keep it positive
# Topt normal prior
# Tau normal
# Drew species specific mu and tau
# Got the mu for each species
# This mu then is used as the "c" parameter in the quadratic, which is the bottom of the curve

######################################### By Species #########################################
sink("aedes.by.species.txt")
cat("
model {
  # Priors
#  a ~ dexp(0.001)             # Width
  #Topt ~ dnorm(20, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Species-specific 
  mu.c ~ dnorm(0, 0.1)             # Mean for species effects c
  tau.c ~ dexp(0.01)             # Tau for species effect c
  
  mu.Topt ~ dnorm(30, 0.1)             # Mean for species effects Topt
  tau.Topt ~ dexp(0.01)             # Tau for species effect Topt
  
  mu.a ~ dexp(0.01)            # Mean for species effects Topt

  for (j in 1:n.species) {
    mu.species[j] ~ dnorm(mu.c, tau.c)       # Species-specific low point
    mu.species.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Species-specific temp
    mu.species.a[j] ~ dexp(mu.a)    # Species-specific temp
  }

  # Likelihood
  for (i in 1:N.obs) {
    mu[i] <- mu.species.a[species[i]] * (mu.species.Topt[species[i]] - temp[i])^2 + mu.species[species[i]]
    trait[i] ~ dnorm(mu[i], tau1)  
  }}
", file = "aedes.by.species.txt")

# exponential asmpytote which is not zero basically just +c at the end determining the minimum development time

inits <- function() {
  list(
    mu.a = 1,             
    mu.c = 0,          
    tau.c = 1,
    tau1 = 1,
    mu.Topt = 30,
    tau.Topt = 1
  )
}

# Parameters to monitor
parameters <- c("mu.a", "mu.c", "tau.c", "mu.species", 
                "tau1", "mu.Topt", "tau.Topt", "mu.species.Topt", "mu.species.a")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Bundle all data in a list for JAGS
#jag.data<-list(trait = test$Trait, N.obs = N.obs, temp = test$Temp)

jag.data<-list(trait = test$Trait, N.obs = N.obs, temp = test$Temp, 
               n.species = n.species, species = species)

species.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
               model.file="aedes.by.species.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
               n.iter=ni, DIC=T, working.directory=getwd())
species.fit.mcmc <- as.mcmc(species.fit) ## makes an "mcmc" object
closeAllConnections()
species.fit$BUGSoutput$DIC


# Get posterior means
summary_fit.species <- summary(species.fit.mcmc)$statistics
a_mean_species <- summary_fit.species[grep("^mu.species.a\\[", rownames(summary_fit.species)), "Mean"]
c_mean_species <- summary_fit.species[grep("^mu.species\\[", rownames(summary_fit.species)), "Mean"]
Topt_mean_species <- summary_fit.species[grep("^mu.species.Topt\\[", rownames(summary_fit.species)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

# Get predicted values
spec1 <- function (x){
  a_mean_species[1] * (x - Topt_mean_species[1])^2 + c_mean_species[1]
}
spec2 <- function (x){
  a_mean_species[2] * (x - Topt_mean_species[2])^2 + c_mean_species[2]
}
spec3 <- function (x){
  a_mean_species[3] * (x - Topt_mean_species[3])^2 + c_mean_species[3]
}
spec4 <- function (x){
  a_mean_species[4] * (x - Topt_mean_species[4])^2 + c_mean_species[4]
}
spec5 <- function (x){
  a_mean_species[5] * (x - Topt_mean_species[5])^2 + c_mean_species[5]
}
spec6 <- function (x){
  a_mean_species[6] * (x - Topt_mean_species[6])^2 + c_mean_species[6]
}
spec7 <- function (x){
  a_mean_species[7] * (x - Topt_mean_species[7])^2 + c_mean_species[7]
}

# Optionally, convert to a data frame for easier plotting
fitted_data.species <- data.frame(Temperature = rep(temp_seq, n.species), 
                          Fitted = c(spec1(temp_seq),
                                     spec2(temp_seq),
                                     spec3(temp_seq),
                                     spec4(temp_seq),
                                     spec5(temp_seq),
                                     spec6(temp_seq),
                                     spec7(temp_seq)),
                          Species = factor(rep(1:n.species, each = length(temp_seq))))

plot_data.species <- data.frame(
  Temperature = test$Temp,
  Trait = test$Trait,
  Species = as.numeric(test$interactor1)
)

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.species, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.species, aes(x = Temperature, y = Fitted, color = Species), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Species",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 



# Ensure test and fitted_data are aligned correctly
# Make sure the Species column in both datasets uses the same format (factor or numeric)
fitted_data.species$Species <- as.numeric(fitted_data.species$Species)

# Merge the observed and predicted data for comparison
merged_data.species <- fitted_data.species %>%
  left_join(plot_data.species, by = c("Temperature", "Species"))


# Calculate RMSE for each species
rmse_by_species.species <- merged_data.species %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_species.species)

mean(rmse_by_species.species$RMSE)

######################################### By Full Data #########################################
sink("aedes.one.data.txt")
cat("
model {
  # Priors
   a ~ dexp(0.1)             # Width
  Topt ~ dnorm(30, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  c ~ dnorm(10, 0.001)
  
  # Likelihood
  for (i in 1:N.obs) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], tau1)  
  }}
", file = "aedes.one.data.txt")


inits <- function() {
  list(
    a = 0.1,             
    c = 8,
    tau1 = 1,
    Topt = 30
  )
}

# Parameters to monitor
parameters <- c("a", "c", "tau1", "Topt")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Bundle all data in a list for JAGS
#jag.data<-list(trait = test$Trait, N.obs = N.obs, temp = test$Temp)

jag.data<-list(trait = test$Trait, N.obs = N.obs, temp = test$Temp)

one.data.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
               model.file="aedes.one.data.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
               n.iter=ni, DIC=T, working.directory=getwd())
one.data.fit.mcmc <- as.mcmc(one.data.fit) ## makes an "mcmc" object
closeAllConnections()
one.data.fit$BUGSoutput$DIC

# Get posterior means
summary_fit <- summary(one.data.fit.mcmc)$statistics
a_mean <- summary_fit["a", "Mean"]
c_mean <- summary_fit["c", "Mean"]
means_Topt <- summary_fit["Topt", "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

# Get predicted values
full.data <- function (x){
  a_mean * (x - means_Topt)^2 + c_mean
}


# Optionally, convert to a data frame for easier plotting
fitted_data.one.data <- data.frame(Temperature = rep(temp_seq, n.species), 
                          Fitted = full.data(temp_seq),
                          Species = factor(rep(1:n.species, each = length(temp_seq))))

plot_data.one.data <- data.frame(
  Temperature = test$Temp,
  Trait = test$Trait,
  Species = as.numeric(test$interactor1)
)

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.one.data, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.one.data, aes(x = Temperature, y = Fitted, color = Species), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Species",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 

# Ensure test and fitted_data are aligned correctly
# Make sure the Species column in both datasets uses the same format (factor or numeric)
fitted_data.one.data$Species <- as.numeric(fitted_data.one.data$Species)

# Merge the observed and predicted data for comparison
merged_data.one.data <- fitted_data.one.data %>%
  left_join(plot_data.one.data, by = c("Temperature", "Species"))


# Calculate RMSE for each species
rmse_by_species.one.data <- merged_data.one.data %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_species.one.data)
mean(rmse_by_species.one.data$RMSE)


######################################### By Climate #########################################

climate <- as.numeric(test$climate)
n.climate <- length(unique(test$climate))
sink("aedes.by.climate.txt")
cat("
model {
  # Priors
#  a ~ dexp(0.001)             # Width
  #Topt ~ dnorm(20, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Species-specific 
  mu.c ~ dnorm(0, 0.1)             # Mean for species effects c
  tau.c ~ dexp(0.01)             # Tau for species effect c
  
  mu.Topt ~ dnorm(30, 0.1)             # Mean for species effects Topt
  tau.Topt ~ dexp(0.01)             # Tau for species effect Topt
  
  mu.a ~ dexp(0.01)            # Mean for species effects Topt

  for (j in 1:n.climate) {
    mu.climate[j] ~ dnorm(mu.c, tau.c)       # Species-specific low point
    mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Species-specific temp
    mu.climate.a[j] ~ dexp(mu.a)    # Species-specific temp
  }

  # Likelihood
  for (i in 1:N.obs) {
    mu[i] <- mu.climate.a[climate[i]] * (mu.climate.Topt[climate[i]] - temp[i])^2 + mu.climate[climate[i]]
    trait[i] ~ dnorm(mu[i], tau1)  
  }}
", file = "aedes.by.climate.txt")

# exponential asmpytote which is not zero basically just +c at the end determining the minimum development time

inits <- function() {
  list(
    mu.a = 1,             
    mu.c = 0,          
    tau.c = 1,
    tau1 = 1,
    mu.Topt = 30,
    tau.Topt = 1
  )
}

# Parameters to monitor
parameters <- c("mu.a", "mu.c", "tau.c", "mu.climate", 
                "tau1", "mu.Topt", "tau.Topt", "mu.climate.Topt", "mu.climate.a")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Bundle all data in a list for JAGS
#jag.data<-list(trait = test$Trait, N.obs = N.obs, temp = test$Temp)

jag.data<-list(trait = test$Trait, N.obs = N.obs, temp = test$Temp, 
               n.climate = n.climate, climate = climate)

climate.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
               model.file="aedes.by.climate.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
               n.iter=ni, DIC=T, working.directory=getwd())
climate.fit.mcmc <- as.mcmc(climate.fit) ## makes an "mcmc" object
closeAllConnections()
climate.fit$BUGSoutput$DIC

# Get posterior means
summary_fit.climate <- summary(climate.fit.mcmc)$statistics
a_mean_climate <- summary_fit.climate[grep("^mu.climate.a\\[", rownames(summary_fit.climate)), "Mean"]
c_mean_climate <- summary_fit.climate[grep("^mu.climate\\[", rownames(summary_fit.climate)), "Mean"]
Topt_mean_climate <- summary_fit.climate[grep("^mu.climate.Topt\\[", rownames(summary_fit.climate)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

# Get predicted values
clim1 <- function (x){
  a_mean_climate[1] * (x - Topt_mean_climate[1])^2 + c_mean_climate[1]
}
clim2 <- function (x){
  a_mean_climate[2] * (x - Topt_mean_climate[2])^2 + c_mean_climate[2]
}
clim3 <- function (x){
  a_mean_climate[3] * (x - Topt_mean_climate[3])^2 + c_mean_climate[3]
}


# Optionally, convert to a data frame for easier plotting
fitted_data.climate <- data.frame(Temperature = rep(temp_seq, n.climate), 
                          Fitted = c(clim1(temp_seq),
                                     clim2(temp_seq),
                                     clim3(temp_seq)),
                          Climate = factor(rep(1:n.climate, each = length(temp_seq))))

plot_data.climate <- data.frame(
  Temperature = test$Temp,
  Trait = test$Trait,
  Climate = as.numeric(test$climate),
  Species = as.numeric(test$interactor1)
)

species_to_climate <- unique(plot_data.climate[, c("Species", "Climate")])  # Ensure species matches correct climate
species_to_climate$Climate <- as.factor(species_to_climate$Climate)
# Merge fitted climate curves with species-to-climate mapping
fitted_data.climate <- data.frame(
  Temperature = rep(temp_seq, n.climate),
  Fitted = c(clim1(temp_seq),
             clim2(temp_seq),
             clim3(temp_seq)),
  Climate = factor(rep(1:n.climate, each = length(temp_seq)))
) %>%
  left_join(species_to_climate, by = "Climate")

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.climate, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.climate, aes(x = Temperature, y = Fitted, color = Climate), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Species",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 



# Ensure test and fitted_data are aligned correctly
# Make sure the Species column in both datasets uses the same format (factor or numeric)
fitted_data.climate$Climate <- as.numeric(fitted_data.climate$Climate)

# Merge the observed and predicted data for comparison
merged_data.climate <- fitted_data.climate %>%
  left_join(plot_data.climate, by = c("Temperature", "Climate", "Species"))


# Calculate RMSE for each species
rmse_by_species.climate <- merged_data.climate %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_species.climate)

mean(rmse_by_species.climate$RMSE)





######################################### By Wing Size Category #########################################
size <- as.numeric(test$wing.size)
n.size <- length(unique(test$wing.size))

sink("aedes.by.size.txt")
cat("
model {
  # Priors
#  a ~ dexp(0.001)             # Width
  #Topt ~ dnorm(20, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Species-specific 
  mu.c ~ dnorm(0, 0.1)             # Mean for species effects c
  tau.c ~ dexp(0.01)             # Tau for species effect c
  
  mu.Topt ~ dnorm(30, 0.1)             # Mean for species effects Topt
  tau.Topt ~ dexp(0.01)             # Tau for species effect Topt
  
  mu.a ~ dexp(0.01)            # Mean for species effects Topt

  for (j in 1:n.size) {
    mu.size[j] ~ dnorm(mu.c, tau.c)       # Species-specific low point
    mu.size.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Species-specific temp
    mu.size.a[j] ~ dexp(mu.a)    # Species-specific temp
  }

  # Likelihood
  for (i in 1:N.obs) {
    mu[i] <- mu.size.a[size[i]] * (mu.size.Topt[size[i]] - temp[i])^2 + mu.size[size[i]]
    trait[i] ~ dnorm(mu[i], tau1)  
  }}
", file = "aedes.by.size.txt")

# exponential asmpytote which is not zero basically just +c at the end determining the minimum development time

inits <- function() {
  list(
    mu.a = 1,             
    mu.c = 0,          
    tau.c = 1,
    tau1 = 1,
    mu.Topt = 30,
    tau.Topt = 1
  )
}

# Parameters to monitor
parameters <- c("mu.a", "mu.c", "tau.c", "mu.size", 
                "tau1", "mu.Topt", "tau.Topt", "mu.size.Topt", "mu.size.a")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Bundle all data in a list for JAGS
#jag.data<-list(trait = test$Trait, N.obs = N.obs, temp = test$Temp)

jag.data<-list(trait = test$Trait, N.obs = N.obs, temp = test$Temp, 
               n.size = n.size, size = size)

size.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
               model.file="aedes.by.size.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
               n.iter=ni, DIC=T, working.directory=getwd())
size.fit.mcmc <- as.mcmc(size.fit) ## makes an "mcmc" object
closeAllConnections()
size.fit$BUGSoutput$DIC


# Get posterior means
summary_fit.size <- summary(size.fit.mcmc)$statistics
a_mean_size <- summary_fit.size[grep("^mu.size.a\\[", rownames(summary_fit.size)), "Mean"]
c_mean_size <- summary_fit.size[grep("^mu.size\\[", rownames(summary_fit.size)), "Mean"]
Topt_mean_size <- summary_fit.size[grep("^mu.size.Topt\\[", rownames(summary_fit.size)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

# Get predicted values
size1 <- function (x){
  a_mean_size[1] * (x - Topt_mean_size[1])^2 + c_mean_size[1]
}
size2 <- function (x){
  a_mean_size[2] * (x - Topt_mean_size[2])^2 + c_mean_size[2]
}
size3 <- function (x){
  a_mean_size[3] * (x - Topt_mean_size[3])^2 + c_mean_size[3]
}


# Optionally, convert to a data frame for easier plotting
fitted_data.size <- data.frame(Temperature = rep(temp_seq, n.size), 
                          Fitted = c(size1(temp_seq),
                                     size2(temp_seq),
                                     size3(temp_seq)),
                          Size = factor(rep(1:n.size, each = length(temp_seq))))

plot_data.size <- data.frame(
  Temperature = test$Temp,
  Trait = test$Trait,
  Size = as.numeric(test$wing.size),
  Species = as.numeric(test$interactor1)
)

species_to_size <- unique(plot_data.size[, c("Species", "Size")])  # Ensure species matches correct climate
species_to_size$Size <- as.factor(species_to_size$Size)
# Merge fitted climate curves with species-to-climate mapping
fitted_data.size <- data.frame(
  Temperature = rep(temp_seq, n.size),
  Fitted = c(size1(temp_seq),
             size2(temp_seq),
             size3(temp_seq)),
  Size = factor(rep(1:n.size, each = length(temp_seq)))
) %>%
  left_join(species_to_size, by = "Size")

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.size, aes(x = Temperature, y = Fitted, color = Size), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Size",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 



# Ensure test and fitted_data are aligned correctly
# Make sure the Species column in both datasets uses the same format (factor or numeric)
fitted_data.size$Size <- as.numeric(fitted_data.size$Size)

# Merge the observed and predicted data for comparison
merged_data.size <- fitted_data.size %>%
  left_join(plot_data.size, by = c("Temperature", "Size", "Species"))


# Calculate RMSE for each species
rmse_by_species.size <- merged_data.size %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_species.size)

mean(rmse_by_species.size$RMSE)




######################################### By Species- Take temperature out #########################################
aedes <- test[!(test$Temp %in% c(22, 25, 27, 34)), ]
N.obs.no25 <- length(aedes$Trait)
species.no25 <- as.numeric(aedes$interactor1)
n.species.no25 <- length(unique(aedes$interactor1))
sink("aedes.by.species.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Species-specific 
  mu.c ~ dnorm(0, 0.1)             # Mean for species effects c
  tau.c ~ dexp(0.01)             # Tau for species effect c
  
  mu.Topt ~ dnorm(30, 0.1)             # Mean for species effects Topt
  tau.Topt ~ dexp(0.01)             # Tau for species effect Topt
  
  mu.a ~ dexp(0.01)            # Mean for species effects Topt

  for (j in 1:n.species.no25) {
    mu.species[j] ~ dnorm(mu.c, tau.c)       # Species-specific low point
    mu.species.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Species-specific temp
    mu.species.a[j] ~ dexp(mu.a)    # Species-specific temp
  }

  # Likelihood
  for (i in 1:N.obs.no25) {
    mu[i] <- mu.species.a[species.no25[i]] * (mu.species.Topt[species.no25[i]] - temp[i])^2 + mu.species[species.no25[i]]
    trait[i] ~ dnorm(mu[i], tau1)  
  }}
", file = "aedes.by.species.txt")

# exponential asmpytote which is not zero basically just +c at the end determining the minimum development time

inits <- function() {
  list(
    mu.a = 1,             
    mu.c = 0,          
    tau.c = 1,
    tau1 = 1,
    mu.Topt = 30,
    tau.Topt = 1
  )
}

# Parameters to monitor
parameters <- c("mu.a", "mu.c", "tau.c", "mu.species", 
                "tau1", "mu.Topt", "tau.Topt", "mu.species.Topt", "mu.species.a")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Bundle all data in a list for JAGS
#jag.data<-list(trait = test$Trait, N.obs = N.obs, temp = test$Temp)

jag.data <- list(trait = aedes$Trait, N.obs.no25 = N.obs.no25, temp = aedes$Temp, 
               n.species.no25 = n.species.no25, species.no25 = species.no25)

species.fit.no25 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                    model.file="aedes.by.species.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                    n.iter=ni, DIC=T, working.directory=getwd())
species.fit.no25.mcmc <- as.mcmc(species.fit.no25) ## makes an "mcmc" object
closeAllConnections()
species.fit.no25$BUGSoutput$DIC


# Get posterior means
summary_fit.species.no25 <- summary(species.fit.no25.mcmc)$statistics
a_mean_species.no25 <- summary_fit.species.no25[grep("^mu.species.a\\[", rownames(summary_fit.species.no25)), "Mean"]
c_mean_species.no25 <- summary_fit.species.no25[grep("^mu.species\\[", rownames(summary_fit.species.no25)), "Mean"]
Topt_mean_species.no25 <- summary_fit.species.no25[grep("^mu.species.Topt\\[", rownames(summary_fit.species.no25)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

# Get predicted values
spec1.no25 <- function (x){
  a_mean_species.no25[1] * (x - Topt_mean_species.no25[1])^2 + c_mean_species.no25[1]
}
spec2.no25 <- function (x){
  a_mean_species.no25[2] * (x - Topt_mean_species.no25[2])^2 + c_mean_species.no25[2]
}
spec3.no25 <- function (x){
  a_mean_species.no25[3] * (x - Topt_mean_species.no25[3])^2 + c_mean_species.no25[3]
}
spec4.no25 <- function (x){
  a_mean_species.no25[4] * (x - Topt_mean_species.no25[4])^2 + c_mean_species.no25[4]
}
spec5.no25 <- function (x){
  a_mean_species.no25[5] * (x - Topt_mean_species.no25[5])^2 + c_mean_species.no25[5]
}
spec6.no25 <- function (x){
  a_mean_species.no25[6] * (x - Topt_mean_species.no25[6])^2 + c_mean_species.no25[6]
}
spec7.no25 <- function (x){
  a_mean_species.no25[7] * (x - Topt_mean_species.no25[7])^2 + c_mean_species.no25[7]
}

# Optionally, convert to a data frame for easier plotting
fitted_data.species.no25 <- data.frame(Temperature = rep(temp_seq, n.species.no25), 
                                  Fitted = c(spec1.no25(temp_seq),
                                             spec2.no25(temp_seq),
                                             spec3.no25(temp_seq),
                                             spec4.no25(temp_seq),
                                             spec5.no25(temp_seq),
                                             spec6.no25(temp_seq),
                                             spec7.no25(temp_seq)),
                                  Species = factor(rep(1:n.species.no25, each = length(temp_seq))))

plot_data.species <- data.frame(
  Temperature = test$Temp,
  Trait = test$Trait,
  Species = as.numeric(test$interactor1)
)

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.species, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.species.no25, aes(x = Temperature, y = Fitted, color = Species), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Species",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 



# Ensure test and fitted_data are aligned correctly
# Make sure the Species column in both datasets uses the same format (factor or numeric)
fitted_data.species.no25$Species <- as.numeric(fitted_data.species.no25$Species)

# Merge the observed and predicted data for comparison
merged_data.species.no25 <- fitted_data.species.no25 %>%
  left_join(plot_data.species, by = c("Temperature", "Species"))


# Calculate RMSE for each species
rmse_by_species.species.no25 <- merged_data.species.no25 %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_species.species.no25)

mean(rmse_by_species.species.no25$RMSE)

######################################### By Full Data- Take temperature out #########################################
aedes <- test[!(test$Temp %in% c(22, 25, 27, 34)), ]
N.obs.no25 <- length(aedes$Trait)
sink("aedes.one.data.txt")
cat("
model {
  # Priors
   a ~ dexp(0.1)             # Width
  Topt ~ dnorm(30, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  c ~ dnorm(10, 0.001)
  
  # Likelihood
  for (i in 1:N.obs.no25) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], tau1)  
  }}
", file = "aedes.one.data.txt")


inits <- function() {
  list(
    a = 0.1,             
    c = 8,
    tau1 = 1,
    Topt = 30
  )
}

# Parameters to monitor
parameters <- c("a", "c", "tau1", "Topt")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Bundle all data in a list for JAGS
#jag.data<-list(trait = test$Trait, N.obs = N.obs, temp = test$Temp)

jag.data <- list(trait = aedes$Trait, N.obs.no25 = N.obs.no25, temp = aedes$Temp)

one.data.no25.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                     model.file="aedes.one.data.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                     n.iter=ni, DIC=T, working.directory=getwd())
one.data.no25.fit.mcmc <- as.mcmc(one.data.no25.fit) ## makes an "mcmc" object
closeAllConnections()
one.data.no25.fit$BUGSoutput$DIC

# Get posterior means
summary_fit.no25 <- summary(one.data.no25.fit.mcmc)$statistics
a_mean.no25 <- summary_fit.no25["a", "Mean"]
c_mean.no25 <- summary_fit.no25["c", "Mean"]
means_Topt.no25 <- summary_fit.no25["Topt", "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

# Get predicted values
full.data.no25 <- function (x){
  a_mean.no25 * (x - means_Topt.no25)^2 + c_mean.no25
}


# Optionally, convert to a data frame for easier plotting
fitted_data.one.data.no25 <- data.frame(Temperature = rep(temp_seq, n.species.no25), 
                                   Fitted = full.data.no25(temp_seq),
                                   Species = factor(rep(1:n.species.no25, each = length(temp_seq))))

plot_data.one.data <- data.frame(
  Temperature = test$Temp,
  Trait = test$Trait,
  Species = as.numeric(test$interactor1)
)

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.one.data, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.one.data.no25, aes(x = Temperature, y = Fitted, color = Species), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Species",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 

# Ensure test and fitted_data are aligned correctly
# Make sure the Species column in both datasets uses the same format (factor or numeric)
fitted_data.one.data.no25$Species <- as.numeric(fitted_data.one.data.no25$Species)

# Merge the observed and predicted data for comparison
merged_data.one.data.no25 <- fitted_data.one.data.no25 %>%
  left_join(plot_data.one.data, by = c("Temperature", "Species"))


# Calculate RMSE for each species
rmse_by_species.one.data.no25 <- merged_data.one.data.no25 %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_species.one.data.no25)
mean(rmse_by_species.one.data.no25$RMSE)


######################################### By Climate- Take temperature out #########################################
aedes <- test[!(test$Temp %in% c(22, 25, 27, 34)), ]
N.obs.no25 <- length(aedes$Trait)
climate.no25 <- as.numeric(aedes$climate)
n.climate.no25 <- length(unique(aedes$climate))
sink("aedes.by.climate.txt")
cat("
model {
  # Priors
#  a ~ dexp(0.001)             # Width
  #Topt ~ dnorm(20, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Species-specific 
  mu.c ~ dnorm(0, 0.1)             # Mean for species effects c
  tau.c ~ dexp(0.01)             # Tau for species effect c
  
  mu.Topt ~ dnorm(30, 0.1)             # Mean for species effects Topt
  tau.Topt ~ dexp(0.01)             # Tau for species effect Topt
  
  mu.a ~ dexp(0.01)            # Mean for species effects Topt

  for (j in 1:n.climate.no25) {
    mu.climate[j] ~ dnorm(mu.c, tau.c)       # Species-specific low point
    mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Species-specific temp
    mu.climate.a[j] ~ dexp(mu.a)    # Species-specific temp
  }

  # Likelihood
  for (i in 1:N.obs.no25) {
    mu[i] <- mu.climate.a[climate.no25[i]] * (mu.climate.Topt[climate.no25[i]] - temp[i])^2 + mu.climate[climate.no25[i]]
    trait[i] ~ dnorm(mu[i], tau1)  
  }}
", file = "aedes.by.climate.txt")

# exponential asmpytote which is not zero basically just +c at the end determining the minimum development time

inits <- function() {
  list(
    mu.a = 1,             
    mu.c = 0,          
    tau.c = 1,
    tau1 = 1,
    mu.Topt = 30,
    tau.Topt = 1
  )
}

# Parameters to monitor
parameters <- c("mu.a", "mu.c", "tau.c", "mu.climate", 
                "tau1", "mu.Topt", "tau.Topt", "mu.climate.Topt", "mu.climate.a")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Bundle all data in a list for JAGS
#jag.data<-list(trait = test$Trait, N.obs = N.obs, temp = test$Temp)

jag.data <- list(trait = aedes$Trait, N.obs.no25 = N.obs.no25, temp = aedes$Temp, 
               n.climate.no25 = n.climate.no25, climate.no25 = climate.no25)

climate.no25.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                    model.file="aedes.by.climate.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                    n.iter=ni, DIC=T, working.directory=getwd())
climate.no25.fit.mcmc <- as.mcmc(climate.no25.fit) ## makes an "mcmc" object
closeAllConnections()
climate.no25.fit$BUGSoutput$DIC

# Get posterior means
summary_fit.climate.no25 <- summary(climate.no25.fit.mcmc)$statistics
a_mean_climate.no25 <- summary_fit.climate.no25[grep("^mu.climate.a\\[", rownames(summary_fit.climate.no25)), "Mean"]
c_mean_climate.no25 <- summary_fit.climate.no25[grep("^mu.climate\\[", rownames(summary_fit.climate.no25)), "Mean"]
Topt_mean_climate.no25 <- summary_fit.climate.no25[grep("^mu.climate.Topt\\[", rownames(summary_fit.climate.no25)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

# Get predicted values
clim1.no25 <- function (x){
  a_mean_climate.no25[1] * (x - Topt_mean_climate.no25[1])^2 + c_mean_climate.no25[1]
}
clim2.no25 <- function (x){
  a_mean_climate.no25[2] * (x - Topt_mean_climate.no25[2])^2 + c_mean_climate.no25[2]
}
clim3.no25 <- function (x){
  a_mean_climate.no25[3] * (x - Topt_mean_climate.no25[3])^2 + c_mean_climate.no25[3]
}


# Optionally, convert to a data frame for easier plotting
fitted_data.climate.no25 <- data.frame(Temperature = rep(temp_seq, n.climate.no25), 
                                  Fitted = c(clim1.no25(temp_seq),
                                             clim2.no25(temp_seq),
                                             clim3.no25(temp_seq)),
                                  Climate = factor(rep(1:n.climate.no25, each = length(temp_seq))))

plot_data.climate <- data.frame(
  Temperature = test$Temp,
  Trait = test$Trait,
  Climate = as.numeric(test$climate),
  Species = as.numeric(test$interactor1)
)
plot_data.climate.no25 <- data.frame(
  Temperature = aedes$Temp,
  Trait = aedes$Trait,
  Climate = as.numeric(aedes$climate),
  Species = as.numeric(aedes$interactor1)
)

species_to_climate <- unique(plot_data.climate[, c("Species", "Climate")])  # Ensure species matches correct climate
species_to_climate$Climate <- as.factor(species_to_climate$Climate)
# Merge fitted climate curves with species-to-climate mapping
fitted_data.climate.no25 <- data.frame(
  Temperature = rep(temp_seq, n.climate.no25),
  Fitted = c(clim1.no25(temp_seq),
             clim2.no25(temp_seq),
             clim3.no25(temp_seq)),
  Climate = factor(rep(1:n.climate.no25, each = length(temp_seq)))
) %>%
  left_join(species_to_climate, by = "Climate")

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.climate, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.climate.no25, aes(x = Temperature, y = Fitted, color = Climate), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Species",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 



# Ensure test and fitted_data are aligned correctly
# Make sure the Species column in both datasets uses the same format (factor or numeric)
fitted_data.climate.no25$Climate <- as.numeric(fitted_data.climate.no25$Climate)

# Merge the observed and predicted data for comparison
merged_data.climate.no25 <- fitted_data.climate.no25 %>%
  left_join(plot_data.climate, by = c("Temperature", "Climate", "Species"))


# Calculate RMSE for each species
rmse_by_species.climate.no25 <- merged_data.climate.no25 %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_species.climate.no25)

mean(rmse_by_species.climate.no25$RMSE)

######################################### By Wing Size Category- Take temperature out #########################################
aedes <- test[!(test$Temp %in% c(22, 25, 27, 34)), ]
N.obs.no25 <- length(aedes$Trait)
size.no25 <- as.numeric(aedes$wing.size)
n.size.no25 <- length(unique(aedes$wing.size))
sink("aedes.by.size.txt")
cat("
model {
  # Priors
#  a ~ dexp(0.001)             # Width
  #Topt ~ dnorm(20, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Species-specific 
  mu.c ~ dnorm(0, 0.1)             # Mean for species effects c
  tau.c ~ dexp(0.01)             # Tau for species effect c
  
  mu.Topt ~ dnorm(30, 0.1)             # Mean for species effects Topt
  tau.Topt ~ dexp(0.01)             # Tau for species effect Topt
  
  mu.a ~ dexp(0.01)            # Mean for species effects Topt

  for (j in 1:n.size.no25) {
    mu.size[j] ~ dnorm(mu.c, tau.c)       # Species-specific low point
    mu.size.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Species-specific temp
    mu.size.a[j] ~ dexp(mu.a)    # Species-specific temp
  }

  # Likelihood
  for (i in 1:N.obs.no25) {
    mu[i] <- mu.size.a[size.no25[i]] * (mu.size.Topt[size.no25[i]] - temp[i])^2 + mu.size[size.no25[i]]
    trait[i] ~ dnorm(mu[i], tau1)  
  }}
", file = "aedes.by.size.txt")

# exponential asmpytote which is not zero basically just +c at the end determining the minimum development time

inits <- function() {
  list(
    mu.a = 1,             
    mu.c = 0,          
    tau.c = 1,
    tau1 = 1,
    mu.Topt = 30,
    tau.Topt = 1
  )
}

# Parameters to monitor
parameters <- c("mu.a", "mu.c", "tau.c", "mu.size", 
                "tau1", "mu.Topt", "tau.Topt", "mu.size.Topt", "mu.size.a")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Bundle all data in a list for JAGS
#jag.data<-list(trait = test$Trait, N.obs = N.obs, temp = test$Temp)

jag.data <- list(trait = test$Trait, N.obs.no25 = N.obs.no25, temp = test$Temp, 
               n.size.no25 = n.size.no25, size.no25 = size.no25)

size.no25.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                 model.file="aedes.by.size.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                 n.iter=ni, DIC=T, working.directory=getwd())
size.no25.fit.mcmc <- as.mcmc(size.no25.fit) ## makes an "mcmc" object
closeAllConnections()
size.no25.fit$BUGSoutput$DIC


# Get posterior means
summary_fit.size.no25 <- summary(size.no25.fit.mcmc)$statistics
a_mean_size.no25 <- summary_fit.size.no25[grep("^mu.size.a\\[", rownames(summary_fit.size.no25)), "Mean"]
c_mean_size.no25 <- summary_fit.size.no25[grep("^mu.size\\[", rownames(summary_fit.size.no25)), "Mean"]
Topt_mean_size.no25 <- summary_fit.size.no25[grep("^mu.size.Topt\\[", rownames(summary_fit.size.no25)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

# Get predicted values
size1.no25 <- function (x){
  a_mean_size.no25[1] * (x - Topt_mean_size.no25[1])^2 + c_mean_size.no25[1]
}
size2.no25 <- function (x){
  a_mean_size.no25[2] * (x - Topt_mean_size.no25[2])^2 + c_mean_size.no25[2]
}
size3.no25 <- function (x){
  a_mean_size.no25[3] * (x - Topt_mean_size.no25[3])^2 + c_mean_size.no25[3]
}


# Optionally, convert to a data frame for easier plotting
fitted_data.size.no25 <- data.frame(Temperature = rep(temp_seq, n.size.no25), 
                               Fitted = c(size1.no25(temp_seq),
                                          size2.no25(temp_seq),
                                          size3.no25(temp_seq)),
                               Size = factor(rep(1:n.size.no25, each = length(temp_seq))))

plot_data.size <- data.frame(
  Temperature = test$Temp,
  Trait = test$Trait,
  Size = as.numeric(test$wing.size),
  Species = as.numeric(test$interactor1)
)

species_to_size <- unique(plot_data.size[, c("Species", "Size")])  # Ensure species matches correct climate
species_to_size$Size <- as.factor(species_to_size$Size)
# Merge fitted climate curves with species-to-climate mapping
fitted_data.size.no25 <- data.frame(
  Temperature = rep(temp_seq, n.size.no25),
  Fitted = c(size1.no25(temp_seq),
             size2.no25(temp_seq),
             size3.no25(temp_seq)),
  Size = factor(rep(1:n.size.no25, each = length(temp_seq)))
) %>%
  left_join(species_to_size, by = "Size")

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.size.no25, aes(x = Temperature, y = Fitted, color = Size), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Size",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 



# Ensure test and fitted_data are aligned correctly
# Make sure the Species column in both datasets uses the same format (factor or numeric)
fitted_data.size.no25$Size <- as.numeric(fitted_data.size.no25$Size)

# Merge the observed and predicted data for comparison
merged_data.size.no25 <- fitted_data.size.no25 %>%
  left_join(plot_data.size, by = c("Temperature", "Size", "Species"))


# Calculate RMSE for each species
rmse_by_species.size.no25 <- merged_data.size.no25 %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_species.size.no25)

mean(rmse_by_species.size.no25$RMSE)




######################################### RMSE, DIC- Take temperature out #########################################
table.aedes.no25 <- matrix(0, nrow = 4, ncol = 2)
row.names(table.aedes.no25) <- c("Full Data", "By Species", "By Climate", "By Wing Size")
colnames(table.aedes.no25) <- c("DIC", "RMSE")
table.aedes.no25[1, 1] <- one.data.no25.fit$BUGSoutput$DIC
table.aedes.no25[2, 1] <- species.fit.no25$BUGSoutput$DIC
table.aedes.no25[3, 1] <- climate.no25.fit$BUGSoutput$DIC
table.aedes.no25[4, 1] <- size.no25.fit$BUGSoutput$DIC
table.aedes.no25[1, 2] <- mean(rmse_by_species.one.data.no25$RMSE)
table.aedes.no25[2, 2] <- mean(rmse_by_species.species.no25$RMSE)
table.aedes.no25[3, 2] <- mean(rmse_by_species.climate.no25$RMSE)
table.aedes.no25[4, 2] <- mean(rmse_by_species.size.no25$RMSE)

table.aedes.no25







######################################### By Full Data- Out of sample prediction #########################################
aedes.rm <- test[test$interactor1 != "Aedes japonicus japonicus", ]
N.obs.rm <- length(aedes.rm$Trait)
sink("aedes.one.data.txt")
cat("
model {
  # Priors
   a ~ dexp(0.1)             # Width
  Topt ~ dnorm(30, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  c ~ dnorm(10, 0.001)
  
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], tau1)  
  }}
", file = "aedes.one.data.txt")


inits <- function() {
  list(
    a = 0.1,             
    c = 8,
    tau1 = 1,
    Topt = 30
  )
}

# Parameters to monitor
parameters <- c("a", "c", "tau1", "Topt")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Bundle all data in a list for JAGS
#jag.data<-list(trait = test$Trait, N.obs = N.obs, temp = test$Temp)

jag.data <- list(trait = aedes.rm$Trait, N.obs.rm = N.obs.rm, temp = aedes.rm$Temp)

one.data.jj.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                     model.file="aedes.one.data.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                     n.iter=ni, DIC=T, working.directory=getwd())
one.data.jj.fit.mcmc <- as.mcmc(one.data.jj.fit) ## makes an "mcmc" object
closeAllConnections()
one.data.rm.fit$BUGSoutput$DIC

# Get posterior means
summary_fit.rm <- summary(one.data.jj.fit.mcmc)$statistics
a_mean.rm <- summary_fit.rm["a", "Mean"]
c_mean.rm <- summary_fit.rm["c", "Mean"]
means_Topt.rm <- summary_fit.rm["Topt", "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

# Get predicted values
full.data.rm <- function (x){
  a_mean.rm * (x - means_Topt.rm)^2 + c_mean.rm
}


# Optionally, convert to a data frame for easier plotting
fitted_data.one.data.jj <- data.frame(Temperature = rep(temp_seq, n.species), 
                                   Fitted = full.data.rm(temp_seq),
                                   Species = factor(rep(1:n.species, each = length(temp_seq))))

plot_data.one.data <- data.frame(
  Temperature = test$Temp,
  Trait = test$Trait,
  Species = as.numeric(test$interactor1)
)

# Plot observed data and fitted curves with facet wrap by species
p.jj <- ggplot() +
  geom_point(data = plot_data.one.data, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.one.data.jj, aes(x = Temperature, y = Fitted), col = "blue", size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Species",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 

# Ensure test and fitted_data are aligned correctly
# Make sure the Species column in both datasets uses the same format (factor or numeric)
fitted_data.one.data.jj$Species <- as.numeric(fitted_data.one.data.jj$Species)

# Merge the observed and predicted data for comparison
merged_data.one.data.jj<- fitted_data.one.data.jj %>%
  left_join(plot_data.one.data, by = c("Temperature", "Species"))


# Calculate RMSE for each species
rmse_by_species.one.data.jj <- merged_data.one.data.jj %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_species.one.data.jj)
mean(rmse_by_species.one.data.jj$RMSE)

save(rmse_by_species.one.data.jj,
     fitted_data.one.data.jj,
     merged_data.one.data.jj,
     one.data.jj.fit.mcmc,
     one.data.jj.fit,
     p.jj,
     file="jj.full.data.RData")


######################################### By Climate- Out of sample prediction #########################################
aedes.rm <- test[test$interactor1 != "Aedes notoscriptus", ]
climate.rm <- as.factor(aedes.rm$climate)
n.climate.rm <- length(unique(aedes.rm$climate))
N.obs.rm <- length(aedes.rm$Trait)
sink("aedes.by.climate.txt")
cat("
model {
  # Priors
#  a ~ dexp(0.001)             # Width
  #Topt ~ dnorm(20, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Species-specific 
  mu.c ~ dnorm(0, 0.1)             # Mean for species effects c
  tau.c ~ dexp(0.01)             # Tau for species effect c
  
  mu.Topt ~ dnorm(30, 0.1)             # Mean for species effects Topt
  tau.Topt ~ dexp(0.01)             # Tau for species effect Topt
  
  mu.a ~ dexp(0.01)            # Mean for species effects Topt

  for (j in 1:n.climate.rm) {
    mu.climate[j] ~ dnorm(mu.c, tau.c)       # Species-specific low point
    mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Species-specific temp
    mu.climate.a[j] ~ dexp(mu.a)    # Species-specific temp
  }

  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.climate.a[climate.rm[i]] * (mu.climate.Topt[climate.rm[i]] - temp[i])^2 + mu.climate[climate.rm[i]]
    trait[i] ~ dnorm(mu[i], tau1)  
  }}
", file = "aedes.by.climate.txt")

# exponential asmpytote which is not zero basically just +c at the end determining the minimum development time

inits <- function() {
  list(
    mu.a = 1,             
    mu.c = 0,          
    tau.c = 1,
    tau1 = 1,
    mu.Topt = 30,
    tau.Topt = 1
  )
}

# Parameters to monitor
parameters <- c("mu.a", "mu.c", "tau.c", "mu.climate", 
                "tau1", "mu.Topt", "tau.Topt", "mu.climate.Topt", "mu.climate.a")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Bundle all data in a list for JAGS
#jag.data<-list(trait = test$Trait, N.obs = N.obs, temp = test$Temp)

jag.data <- list(trait = aedes.rm$Trait, N.obs.rm = N.obs.rm, temp = aedes.rm$Temp, 
               n.climate.rm = n.climate.rm, climate.rm = climate.rm)

climate.notoscriptus.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                    model.file="aedes.by.climate.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                    n.iter=ni, DIC=T, working.directory=getwd())
climate.notoscriptus.fit.mcmc <- as.mcmc(climate.notoscriptus.fit) ## makes an "mcmc" object
closeAllConnections()
climate.rm.fit$BUGSoutput$DIC

# Get posterior means
summary_fit.rm.climate <- summary(climate.notoscriptus.fit.mcmc)$statistics
a_mean_climate.rm <- summary_fit.rm.climate[grep("^mu.climate.a\\[", rownames(summary_fit.rm.climate)), "Mean"]
c_mean_climate.rm <- summary_fit.rm.climate[grep("^mu.climate\\[", rownames(summary_fit.rm.climate)), "Mean"]
Topt_mean_climate.rm <- summary_fit.rm.climate[grep("^mu.climate.Topt\\[", rownames(summary_fit.rm.climate)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

# Get predicted values
clim1.rm <- function (x){
  a_mean_climate.rm[1] * (x - Topt_mean_climate.rm[1])^2 + c_mean_climate.rm[1]
}
clim2.rm <- function (x){
  a_mean_climate.rm[2] * (x - Topt_mean_climate.rm[2])^2 + c_mean_climate.rm[2]
}
clim3.rm <- function (x){
  a_mean_climate.rm[3] * (x - Topt_mean_climate.rm[3])^2 + c_mean_climate.rm[3]
}


# Optionally, convert to a data frame for easier plotting
fitted_data.climate.notoscriptus <- data.frame(Temperature = rep(temp_seq, n.climate.rm), 
                                  Fitted = c(clim1.rm(temp_seq),
                                             clim2.rm(temp_seq),
                                             clim3.rm(temp_seq)),
                                  Climate = factor(rep(c("subtropical", "temperate", "tropical"), each = length(temp_seq))))

plot_data.climate <- data.frame(
  Temperature = test$Temp,
  Trait = test$Trait,
  Climate = as.factor(test$climate),
  Species = as.factor(test$interactor1)
)

species_to_climate <- unique(plot_data.climate[, c("Species", "Climate")])  # Ensure species matches correct climate
species_to_climate$Climate <- as.factor(species_to_climate$Climate)
# Merge fitted climate curves with species-to-climate mapping
fitted_data.climate.notoscriptus <- data.frame(
  Temperature = rep(temp_seq, n.climate.rm),
  Fitted = c(clim1.rm(temp_seq),
             clim2.rm(temp_seq),
             clim3.rm(temp_seq)),
  Climate = factor(rep(c("subtropical", "temperate", "tropical"), each = length(temp_seq)))
) %>%
  left_join(species_to_climate, by = "Climate")

# Plot observed data and fitted curves with facet wrap by species
p.climate.notoscriptus <- ggplot() +
  geom_point(data = plot_data.climate, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.climate.notoscriptus, aes(x = Temperature, y = Fitted, color = Climate), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Species",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 



# Ensure test and fitted_data are aligned correctly
# Make sure the Species column in both datasets uses the same format (factor or numeric)
#fitted_data.climate.notoscriptus$Climate <- as.numeric(fitted_data.climate.notoscriptus$Climate)

# Merge the observed and predicted data for comparison
merged_data.climate.notoscriptus <- fitted_data.climate.notoscriptus %>%
  left_join(plot_data.climate, by = c("Temperature", "Climate", "Species"))


# Calculate RMSE for each species
rmse_by_species.climate.notoscriptus <- merged_data.climate.notoscriptus %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_species.climate.notoscriptus)

mean(rmse_by_species.climate.notoscriptus$RMSE)

save(rmse_by_species.climate.notoscriptus,
     fitted_data.climate.notoscriptus,
     merged_data.climate.notoscriptus,
     climate.notoscriptus.fit.mcmc,
     climate.notoscriptus.fit,
     p.climate.notoscriptus,
     file="notoscriptus.climate.RData")


######################################### By Wing Size Category- Out of sample prediction #########################################
aedes.rm <- test[test$interactor1 != "Aedes notoscriptus", ]
size.rm <- as.numeric(aedes.rm$wing.size)
n.size.rm <- length(unique(aedes.rm$wing.size))
N.obs.rm <- length(aedes.rm$Trait)

sink("aedes.by.size.txt")
cat("
model {
  # Priors
#  a ~ dexp(0.001)             # Width
  #Topt ~ dnorm(20, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Species-specific 
  mu.c ~ dnorm(0, 0.1)             # Mean for species effects c
  tau.c ~ dexp(0.01)             # Tau for species effect c
  
  mu.Topt ~ dnorm(30, 0.1)             # Mean for species effects Topt
  tau.Topt ~ dexp(0.01)             # Tau for species effect Topt
  
  mu.a ~ dexp(0.01)            # Mean for species effects Topt

  for (j in 1:n.size.rm) {
    mu.size[j] ~ dnorm(mu.c, tau.c)       # Species-specific low point
    mu.size.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Species-specific temp
    mu.size.a[j] ~ dexp(mu.a)    # Species-specific temp
  }

  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.size.a[size.rm[i]] * (mu.size.Topt[size.rm[i]] - temp[i])^2 + mu.size[size.rm[i]]
    trait[i] ~ dnorm(mu[i], tau1)  
  }}
", file = "aedes.by.size.txt")

# exponential asmpytote which is not zero basically just +c at the end determining the minimum development time

inits <- function() {
  list(
    mu.a = 1,             
    mu.c = 0,          
    tau.c = 1,
    tau1 = 1,
    mu.Topt = 30,
    tau.Topt = 1
  )
}

# Parameters to monitor
parameters <- c("mu.a", "mu.c", "tau.c", "mu.size", 
                "tau1", "mu.Topt", "tau.Topt", "mu.size.Topt", "mu.size.a")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Bundle all data in a list for JAGS
#jag.data<-list(trait = test$Trait, N.obs = N.obs, temp = test$Temp)

jag.data <- list(trait = aedes.rm$Trait, N.obs.rm = N.obs.rm, temp = aedes.rm$Temp, 
               n.size.rm = n.size.rm, size.rm = size.rm)

size.rm.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                 model.file="aedes.by.size.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                 n.iter=ni, DIC=T, working.directory=getwd())
size.rm.fit.mcmc <- as.mcmc(size.rm.fit) ## makes an "mcmc" object
closeAllConnections()
size.rm.fit$BUGSoutput$DIC


# Get posterior means
summary_fit.size.rm <- summary(size.rm.fit.mcmc)$statistics
a_mean_size.rm <- summary_fit.size.rm[grep("^mu.size.a\\[", rownames(summary_fit.size.rm)), "Mean"]
c_mean_size.rm <- summary_fit.size.rm[grep("^mu.size\\[", rownames(summary_fit.size.rm)), "Mean"]
Topt_mean_size.rm <- summary_fit.size.rm[grep("^mu.size.Topt\\[", rownames(summary_fit.size.rm)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

# Get predicted values
size1.rm <- function (x){
  a_mean_size.rm[1] * (x - Topt_mean_size.rm[1])^2 + c_mean_size.rm[1]
}
size2.rm <- function (x){
  a_mean_size.rm[2] * (x - Topt_mean_size.rm[2])^2 + c_mean_size.rm[2]
}
size3.rm <- function (x){
  a_mean_size.rm[3] * (x - Topt_mean_size.rm[3])^2 + c_mean_size.rm[3]
}


# Optionally, convert to a data frame for easier plotting
fitted_data.size.rm <- data.frame(Temperature = rep(temp_seq, n.size.rm), 
                               Fitted = c(size1.rm(temp_seq),
                                          size2.rm(temp_seq),
                                          size3.rm(temp_seq)),
                               Size = factor(rep(1:n.size.rm, each = length(temp_seq))))

plot_data.size <- data.frame(
  Temperature = test$Temp,
  Trait = test$Trait,
  Size = as.numeric(test$wing.size),
  Species = as.numeric(test$interactor1)
)

species_to_size <- unique(plot_data.size[, c("Species", "Size")])  # Ensure species matches correct climate
species_to_size$Size <- as.factor(species_to_size$Size)
# Merge fitted climate curves with species-to-climate mapping
fitted_data.size.rm <- data.frame(
  Temperature = rep(temp_seq, n.size.rm),
  Fitted = c(size1.rm(temp_seq),
             size2.rm(temp_seq),
             size3.rm(temp_seq)),
  Size = factor(rep(1:n.size.rm, each = length(temp_seq)))
) %>%
  left_join(species_to_size, by = "Size")

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.size.rm, aes(x = Temperature, y = Fitted, color = Size), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Size",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 



# Ensure test and fitted_data are aligned correctly
# Make sure the Species column in both datasets uses the same format (factor or numeric)
fitted_data.size.rm$Size <- as.numeric(fitted_data.size.rm$Size)

# Merge the observed and predicted data for comparison
merged_data.size.rm <- fitted_data.size.rm %>%
  left_join(plot_data.size, by = c("Temperature", "Size", "Species"))


# Calculate RMSE for each species
rmse_by_species.size.rm <- merged_data.size.rm %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_species.size.rm)

mean(rmse_by_species.size.rm$RMSE)





######################################### By Wing Size and Climate #########################################
climate <- as.numeric(test$climate)
n.climate <- length(unique(test$climate))
size <- as.numeric(test$wing.size)
n.size <- length(unique(test$wing.size))
sink("aedes.by.climate.size.txt")
cat("
model {
  # Priors
#  a ~ dexp(0.001)             # Width
  #Topt ~ dnorm(20, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Species-specific 
  mu.c ~ dnorm(0, 0.1)             # Mean for species effects c
  tau.c ~ dexp(0.01)             # Tau for species effect c
  
  mu.Topt ~ dnorm(30, 0.1)             # Mean for species effects Topt
  tau.Topt ~ dexp(0.01)             # Tau for species effect Topt
  
  mu.a ~ dexp(0.01)            # Mean for species effects Topt

  for (j in 1:n.climate) {
    mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific temp
    mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific temp
  }
    for (k in 1:n.size) {
    mu.size[k] ~ dnorm(mu.c, tau.c)       # Size-specific low point
  }
  # Likelihood
  for (i in 1:N.obs) {
    mu[i] <- mu.climate.a[climate[i]] * (mu.climate.Topt[climate[i]] - temp[i])^2 + mu.size[size[i]]
    trait[i] ~ dnorm(mu[i], tau1)  
  }}
", file = "aedes.by.climate.size.txt")

inits <- function() {
  list(
    mu.a = 1,             
    mu.c = 0,          
    tau.c = 1,
    tau1 = 1,
    mu.Topt = 30,
    tau.Topt = 1
  )
}

# Parameters to monitor
parameters <- c("mu.a", "mu.c", "tau.c", "mu.size", 
                "tau1", "mu.Topt", "tau.Topt", "mu.climate.Topt", "mu.climate.a")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Bundle all data in a list for JAGS
#jag.data<-list(trait = test$Trait, N.obs = N.obs, temp = test$Temp)

jag.data<-list(trait = test$Trait, N.obs = N.obs, temp = test$Temp, 
               n.climate = n.climate, climate = climate, n.size = n.size, size = size)

climate.size.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                    model.file="aedes.by.climate.size.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                    n.iter=ni, DIC=T, working.directory=getwd())
climate.size.fit.mcmc <- as.mcmc(climate.size.fit) ## makes an "mcmc" object
closeAllConnections()
climate.size.fit$BUGSoutput$DIC

# Get posterior means
summary_fit.climate.size <- summary(climate.size.fit.mcmc)$statistics
a_mean_climate.size <- summary_fit.climate.size[grep("^mu.climate.a\\[", rownames(summary_fit.climate.size)), "Mean"]
c_mean_climate.size <- summary_fit.climate.size[grep("^mu.size\\[", rownames(summary_fit.climate.size)), "Mean"]
Topt_mean_climate.size <- summary_fit.climate.size[grep("^mu.climate.Topt\\[", rownames(summary_fit.climate.size)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

# Get predicted values
clim1size1 <- function (x){
  a_mean_climate.size[1] * (x - Topt_mean_climate.size[1])^2 + c_mean_climate.size[1]
}
clim2size1 <- function (x){
  a_mean_climate.size[2] * (x - Topt_mean_climate.size[2])^2 + c_mean_climate.size[1]
}
clim3size1 <- function (x){
  a_mean_climate.size[3] * (x - Topt_mean_climate.size[3])^2 + c_mean_climate.size[1]
}
clim1size2 <- function (x){
  a_mean_climate.size[1] * (x - Topt_mean_climate.size[1])^2 + c_mean_climate.size[2]
}
clim2size2 <- function (x){
  a_mean_climate.size[2] * (x - Topt_mean_climate.size[2])^2 + c_mean_climate.size[2]
}
clim3size2 <- function (x){
  a_mean_climate.size[3] * (x - Topt_mean_climate.size[3])^2 + c_mean_climate.size[2]
}
clim1size3 <- function (x){
  a_mean_climate.size[1] * (x - Topt_mean_climate.size[1])^2 + c_mean_climate.size[3]
}
clim2size3 <- function (x){
  a_mean_climate.size[2] * (x - Topt_mean_climate.size[2])^2 + c_mean_climate.size[3]
}
clim3size3 <- function (x){
  a_mean_climate.size[3] * (x - Topt_mean_climate.size[3])^2 + c_mean_climate.size[3]
}

# Optionally, convert to a data frame for easier plotting
fitted_data.climate.size <- data.frame(Temperature = rep(temp_seq, n.climate*3), 
                                  Fitted = c(clim1size1(temp_seq),
                                             clim2size1(temp_seq),
                                             clim3size1(temp_seq),
                                             clim1size2(temp_seq),
                                             clim2size2(temp_seq),
                                             clim3size2(temp_seq),
                                             clim1size3(temp_seq),
                                             clim2size3(temp_seq),
                                             clim3size3(temp_seq)),
                                  Climate = factor(rep(rep(1:n.climate, each = length(temp_seq)), 3)),
                                  Size = factor(rep(1:n.size, each = length(temp_seq) * 3)))

plot_data.climate.size <- data.frame(
  Temperature = test$Temp,
  Trait = test$Trait,
  Climate = as.numeric(test$climate),
  Species = as.numeric(test$interactor1),
  Size = as.numeric(test$wing.size)
)

species_to_climate.size <- unique(plot_data.climate.size[, c("Species", "Climate", "Size")])
species_to_climate.size$Climate <- as.factor(species_to_climate.size$Climate)
species_to_climate.size$Size <- as.factor(species_to_climate.size$Size)

# Merge fitted climate curves with species-to-climate mapping and filter for matching species
fitted_data.climate.size <- data.frame(
  Temperature = rep(temp_seq, n.climate * n.size), 
  Fitted = c(clim1size1(temp_seq),
             clim2size1(temp_seq),
             clim3size1(temp_seq),
             clim1size2(temp_seq),
             clim2size2(temp_seq),
             clim3size2(temp_seq),
             clim1size3(temp_seq),
             clim2size3(temp_seq),
             clim3size3(temp_seq)),
  Climate = factor(rep(rep(1:n.climate, each = length(temp_seq)), n.size)),
  Size = factor(rep(1:n.size, each = length(temp_seq) * n.climate))
) |> 
  left_join(species_to_climate.size, by = c("Climate", "Size")) |> 
  filter(!is.na(Species))  # Keep only rows where Species is non-missing

# Combine Climate and Size for coloring
fitted_data.climate.size$Combo <- with(fitted_data.climate.size, paste(Climate, Size, sep = ", "))

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.climate.size, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.climate.size, aes(x = Temperature, y = Fitted, color = Combo), size = 1) +
  labs(
    title = "Observed Data and Quadratic Fit by Species",
    x = "Temperature",
    y = "Trait Value",
    color = "Climate_Size Combo"
  ) +
  theme_minimal() +
  facet_wrap(~ Species)


# Ensure test and fitted_data are aligned correctly
# Make sure the Species column in both datasets uses the same format (factor or numeric)
fitted_data.climate.size$Climate <- as.numeric(fitted_data.climate.size$Climate)
fitted_data.climate.size$Size <- as.numeric(fitted_data.climate.size$Size)

# Merge the observed and predicted data for comparison
merged_data.climate.size <- fitted_data.climate.size %>%
  left_join(plot_data.climate.size, by = c("Temperature", "Climate", "Species"))


# Calculate RMSE for each species
rmse_by_species.climate.size <- merged_data.climate.size %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_species.climate.size)

mean(rmse_by_species.climate.size$RMSE)


######################################### RMSE, DIC #########################################
table.aedes <- matrix(0, nrow = 5, ncol = 2)
row.names(table.aedes) <- c("Full Data", "By Species", "By Climate", "By Wing Size", "By Wing Size and Climate")
colnames(table.aedes) <- c("DIC", "RMSE")
table.aedes[1, 1] <- one.data.fit$BUGSoutput$DIC
table.aedes[2, 1] <- species.fit$BUGSoutput$DIC
table.aedes[3, 1] <- climate.fit$BUGSoutput$DIC
table.aedes[4, 1] <- size.fit$BUGSoutput$DIC
table.aedes[5, 1] <- climate.size.fit$BUGSoutput$DIC
table.aedes[1, 2] <- mean(rmse_by_species.one.data$RMSE)
table.aedes[2, 2] <- mean(rmse_by_species.species$RMSE)
table.aedes[3, 2] <- mean(rmse_by_species.climate$RMSE)
table.aedes[4, 2] <- mean(rmse_by_species.size$RMSE)
table.aedes[5, 2] <- mean(rmse_by_species.climate.size$RMSE)

table.aedes









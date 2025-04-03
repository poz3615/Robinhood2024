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


################################### Data Cleaning ###################################
# Set seed for consistent data generating
set.seed(52)
# Read data in
tick.raw <- read.csv("Extra.Data.Tick.csv")
# Convert all characteristics to factors
tick.raw$Species <- as.factor(tick.raw$Species)
tick.raw$Host <- as.factor(tick.raw$Host)
tick.raw$Genus <- as.factor(tick.raw$Genus)
tick.raw$Climate <- as.factor(tick.raw$Climate)
tick.raw$Continent <- as.factor(tick.raw$Continent)
tick.raw$Size <- as.factor(tick.raw$Size)

# Plot by species
ggplot(tick.raw, aes(y = Trait, x = Temp)) +
  geom_point() +
  facet_wrap(~ Species) +
  labs(title = "Temperature by Trait for Each Factor Level",
       x = "Temperature",
       y = "Oviposition") +
  theme_minimal()

# Removing the species that aren't perfectly cleaned for now
tick.raw <- tick.raw[!(tick.raw$Species %in% c("Ixodes anatis", "Rhipicephalus (Boophilus) annulatus", "Haemaphysalis longicornis")), ]
# Drop the factor levels of species that aren't used anymore
tick.raw$Species <- droplevels(tick.raw$Species)

# Species number for hierarchical model
species.raw <- as.factor(tick.raw$Species)
# Number of species
n.species.raw <- length(unique(tick.raw$Species))

# Genus number for hierarchical model
genus.raw <- as.factor(tick.raw$Genus)
# Number of genus
n.genus.raw <- length(unique(tick.raw$Genus))

# Host number for hierarchical model
host.raw <- as.factor(tick.raw$Host)
# Number of hosts
n.host.raw <- length(unique(tick.raw$Host))

# Climate number for hierarchical model
climate.raw <- as.factor(tick.raw$Climate)
# Number of climates
n.climate.raw <- length(unique(tick.raw$Climate))

# Continent number for hierarchical model
continent.raw <- as.factor(tick.raw$Continent)
# Number of continents
n.continent.raw <- length(unique(tick.raw$Continent))

# Size number for hierarchical model
size.raw <- as.factor(tick.raw$Size)
# Small: 0-3
# Medium: 3.01-5
# Large: 5.01+
# Number of sizes
n.size.raw <- length(unique(tick.raw$Size))

N.obs.raw <- length(tick.raw$Trait)
################################### Full Data ###################################
# Model using full data
sink("full_data_robin.txt")
cat("
model {
  # Priors
   a ~ dexp(0.1)             # Width
  Topt ~ dnorm(30, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  c ~ dexp(0.001)   # Oviposition rate at Topt
  
  # Likelihood
  for (i in 1:N.obs.raw) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], tau1) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], tau1)
  }}
", file = "full_data_robin.txt")


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
jag.data <- list(trait = tick.raw$Trait, N.obs.raw = N.obs.raw, 
                 temp = tick.raw$Temp)

# Fit fill model
full.data.fit.raw <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                      model.file = "full_data_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                      n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
full.data.fit.raw.mcmc <- as.mcmc(full.data.fit.raw) 
closeAllConnections()
#DIC
full.data.fit.raw$BUGSoutput$DIC

# Extract log-likelihood matrix (samples × observations)
log.lik.full.raw <- full.data.fit.raw$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.full.raw <- waic(log.lik.full.raw)$estimates["waic", "Estimate"]


# Get posterior means
summary_fit.raw <- summary(full.data.fit.raw.mcmc)$statistics
a_mean.raw <- summary_fit.raw["a", "Mean"]
c_mean.raw <- summary_fit.raw["c", "Mean"]
Topt_mean.raw <- summary_fit.raw["Topt", "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp), max(tick.raw$Temp), by = 0.1)

# Get predicted values
full.data.raw <- function (x){
  a_mean.raw * (x - Topt_mean.raw)^2 + c_mean.raw
}


# Data frame of temperature, fitted value, species
fitted_data.full.data.raw <- data.frame(Temperature = rep(temp_seq, n.species.raw), 
                                    Fitted = full.data.raw(temp_seq),
                                    Species = factor(rep(c("Amblyomma lepidum", "Dermacentor andersoni",
                                                           "Dermacentor nitens", "Haemaphysalis leporispalustris",
                                                           "Hyalomma aegyptium", "Hyalomma dromedarii",
                                                           "Hyalomma impeltatum", "Hyalomma lusitanicum",
                                                           "Hyalomma schulzei"), each = length(temp_seq))))
# Original data with temperature, trait, and species as numeric
plot_data.full.raw <- data.frame(
  Temperature = tick.raw$Temp,
  Trait = tick.raw$Trait,
  Species = as.factor(tick.raw$Species)
)

# Plot observed data and fitted curves with facet wrap by species
full.raw <- ggplot() +
  geom_point(data = plot_data.full.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.full.data.raw, aes(x = Temperature, y = Fitted), color = "lightslateblue", size = 1) +
  labs(title = "Observed Data and Quadratic Fit Using Full Data",
       x = "Temperature", y = "Oviposition Time in Days") +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~ Species) 
full.raw
#ggsave("tick.full.plot.png", width = 8, heigh = 5)

# Merge the observed and predicted data for comparison
merged_data.full.raw <- fitted_data.full.data.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species"))

# Calculate RMSE for each species
mse_by_full_data.raw <- merged_data.full.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.raw <- merged_data.full.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 


save(mse_by_full_data.raw,
     overall_rmse.raw,
     waic.full.raw,
     full.raw,
     full.data.fit.raw$BUGSoutput$DIC,
     file="full.data.raw.RData")


################################### By Species Hierarchical ###################################
# Model using species effect
sink("species_robin.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Species-specific 
  mu.c ~ dexp(10)          # Mean for species effects c (vertex y)
  mu.Topt ~ dnorm(30, 0.1)          # Mean for species effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for species effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for species effects a (width)

  for (j in 1:n.species.raw) {
    mu.species.c[j] ~ dexp(mu.c)       # Species-specific vertex y
    mu.species.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Species-specific vertex x
    mu.species.a[j] ~ dexp(mu.a)    # Species-specific width
  }
  # Likelihood
  for (i in 1:N.obs.raw) {
    mu[i] <- mu.species.a[species.raw[i]] * (mu.species.Topt[species.raw[i]] - temp[i])^2 + mu.species.c[species.raw[i]]
    trait[i] ~ dnorm(mu[i], tau1) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], tau1)
  }}
", file = "species_robin.txt")


# Parameters to monitor
parameters <- c("mu.species.a", "mu.species.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt",
                "tau1", "mu.species.Topt", "log_lik")


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
jag.data <- list(trait = tick.raw$Trait, N.obs.raw = N.obs.raw, 
                 temp = tick.raw$Temp, 
                 n.species.raw = n.species.raw, species.raw = species.raw)
# Fit model
species.raw.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                    model.file = "species_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                    n.iter = ni, DIC = T, working.directory = getwd())
# Make an mcmc object
species.raw.fit.mcmc <- as.mcmc(species.raw.fit) 
closeAllConnections()
species.raw.fit$BUGSoutput$DIC

# Extract log-likelihood matrix (samples × observations)
log.lik.species.raw <- species.raw.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.species.raw <- waic(log.lik.species.raw)$estimates["waic", "Estimate"]



# Get posterior means for each species effect
summary_fit.species.raw <- summary(species.raw.fit.mcmc)$statistics
a_species_mean.raw <- summary_fit.species.raw[grep("^mu.species.a\\[", rownames(summary_fit.species.raw)), "Mean"]
c_species_mean.raw <- summary_fit.species.raw[grep("^mu.species.c\\[", rownames(summary_fit.species.raw)), "Mean"]
Topt_species_mean.raw <- summary_fit.species.raw[grep("^mu.species.Topt\\[", rownames(summary_fit.species.raw)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp) ,
                max(tick.raw$Temp) , by = 0.1)

# Get predicted values for each species using species effect parameters
species1.raw <- function (x){
  a_species_mean.raw[1] * (x - Topt_species_mean.raw[1])^2 + c_species_mean.raw[1]
}
species2.raw <- function (x){
  a_species_mean.raw[2] * (x - Topt_species_mean.raw[2])^2 + c_species_mean.raw[2]
}
species3.raw <- function (x){
  a_species_mean.raw[3] * (x - Topt_species_mean.raw[3])^2 + c_species_mean.raw[3]
}
species4.raw <- function (x){
  a_species_mean.raw[4] * (x - Topt_species_mean.raw[4])^2 + c_species_mean.raw[4]
}
species5.raw <- function (x){
  a_species_mean.raw[5] * (x - Topt_species_mean.raw[5])^2 + c_species_mean.raw[5]
}
species6.raw <- function (x){
  a_species_mean.raw[6] * (x - Topt_species_mean.raw[6])^2 + c_species_mean.raw[6]
}
species7.raw <- function (x){
  a_species_mean.raw[7] * (x - Topt_species_mean.raw[7])^2 + c_species_mean.raw[7]
}
species8.raw <- function (x){
  a_species_mean.raw[8] * (x - Topt_species_mean.raw[8])^2 + c_species_mean.raw[8]
}
species9.raw <- function (x){
  a_species_mean.raw[9] * (x - Topt_species_mean.raw[9])^2 + c_species_mean.raw[9]
}

# Make data frame of temperature, fitted value for each species, and species
fitted_data.species.raw <- data.frame(Temperature = rep(temp_seq, n.species.raw), 
                                  Fitted = c(species1.raw(temp_seq),
                                             species2.raw(temp_seq),
                                             species3.raw(temp_seq),
                                             species4.raw(temp_seq),
                                             species5.raw(temp_seq),
                                             species6.raw(temp_seq),
                                             species7.raw(temp_seq),
                                             species8.raw(temp_seq),
                                             species9.raw(temp_seq)),
                                  Species = factor(rep(c("Amblyomma lepidum", "Dermacentor andersoni",
                                                         "Dermacentor nitens", "Haemaphysalis leporispalustris",
                                                         "Hyalomma aegyptium", "Hyalomma dromedarii",
                                                         "Hyalomma impeltatum", "Hyalomma lusitanicum",
                                                         "Hyalomma schulzei"), each = length(temp_seq))))
# Data frame of original temperature, trait, and numeric species
plot_data.species.raw <- data.frame(
  Temperature = tick.raw$Temp,
  Trait = tick.raw$Trait,
  Species = as.factor(tick.raw$Species)
)

# Plot observed data and fitted curves with facet wrap by species
species.raw.plot <- ggplot() +
  geom_point(data = plot_data.species.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.species.raw, aes(x = Temperature, y = Fitted, color = Species), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Species",
       x = "Temperature", y = "Ovipositon Time in Days") +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~ Species) 
species.raw.plot
#ggsave("tick.species.plot.png", width = 8, height = 5)

# Merge the observed and predicted data for comparison
merged_data.species.raw <- fitted_data.species.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species"))

# Calculate RMSE for each species
mse_by_species.raw <- merged_data.species.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )


# Calculate overall RMSE (no grouping)
overall_rmse.species.raw <- merged_data.species.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 


save(mse_by_species.raw,
     fitted_data.species.raw,
     merged_data.species.raw,
     species.raw.fit.mcmc,
     species.raw.fit,
     log.lik.species.raw,
     file="species.raw.fit.RData")

################################### By Species Hierarchical Different Variance ###################################
# Model using species effect
sink("species_tau_robin.txt")
cat("
model {
  # Priors
  # Species-specific 
  mu.c ~ dexp(10)          # Mean for species effects c (vertex y)
  mu.Topt ~ dnorm(30, 0.1)          # Mean for species effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for species effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for species effects a (width)
  tau1 ~ dexp(0.01)           # Tau for species effect

  for (j in 1:n.species.raw) {
    mu.species.c[j] ~ dexp(mu.c)       # Species-specific vertex y
    mu.species.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Species-specific vertex x
    mu.species.a[j] ~ dexp(mu.a)    # Species-specific width
    species.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.raw) {
    mu[i] <- mu.species.a[species.raw[i]] * (mu.species.Topt[species.raw[i]] - temp[i])^2 + mu.species.c[species.raw[i]]
    trait[i] ~ dnorm(mu[i], species.tau[species.raw[i]]) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], species.tau[species.raw[i]])
  }}
", file = "species_tau_robin.txt")


# Parameters to monitor
parameters <- c("mu.species.a", "mu.species.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt",
                "tau1", "mu.species.Topt", "species.tau", "log_lik")


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
jag.data <- list(trait = tick.raw$Trait, N.obs.raw = N.obs.raw, 
                 temp = tick.raw$Temp, 
                 n.species.raw = n.species.raw, species.raw = species.raw)
# Fit model
species.tau.raw.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                        model.file = "species_tau_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                        n.iter = ni, DIC = T, working.directory = getwd())
# Make an mcmc object
species.tau.raw.fit.mcmc <- as.mcmc(species.tau.raw.fit) 
closeAllConnections()
species.tau.raw.fit$BUGSoutput$DIC

# Extract log-likelihood matrix (samples × observations)
log.lik.species.tau.raw <- species.tau.raw.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.species.tau.raw <- waic(log.lik.species.tau.raw)$estimates["waic", "Estimate"]


# Get posterior means for each species effect
summary_fit.species.tau.raw <- summary(species.tau.raw.fit.mcmc)$statistics
a_species.tau_mean.raw <- summary_fit.species.tau.raw[grep("^mu.species.a\\[", rownames(summary_fit.species.tau.raw)), "Mean"]
c_species.tau_mean.raw <- summary_fit.species.tau.raw[grep("^mu.species.c\\[", rownames(summary_fit.species.tau.raw)), "Mean"]
Topt_species.tau_mean.raw <- summary_fit.species.tau.raw[grep("^mu.species.Topt\\[", rownames(summary_fit.species.tau.raw)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp) ,
                max(tick.raw$Temp) , by = 0.1)

# Get predicted values for each species using species effect parameters
species1.tau.raw <- function (x){
  a_species.tau_mean.raw[1] * (x - Topt_species.tau_mean.raw[1])^2 + c_species.tau_mean.raw[1]
}
species2.tau.raw <- function (x){
  a_species.tau_mean.raw[2] * (x - Topt_species.tau_mean.raw[2])^2 + c_species.tau_mean.raw[2]
}
species3.tau.raw <- function (x){
  a_species.tau_mean.raw[3] * (x - Topt_species.tau_mean.raw[3])^2 + c_species.tau_mean.raw[3]
}
species4.tau.raw <- function (x){
  a_species.tau_mean.raw[4] * (x - Topt_species.tau_mean.raw[4])^2 + c_species.tau_mean.raw[4]
}
species5.tau.raw <- function (x){
  a_species.tau_mean.raw[5] * (x - Topt_species.tau_mean.raw[5])^2 + c_species.tau_mean.raw[5]
}
species6.tau.raw <- function (x){
  a_species.tau_mean.raw[6] * (x - Topt_species.tau_mean.raw[6])^2 + c_species.tau_mean.raw[6]
}
species7.tau.raw <- function (x){
  a_species.tau_mean.raw[7] * (x - Topt_species.tau_mean.raw[7])^2 + c_species.tau_mean.raw[7]
}
species8.tau.raw <- function (x){
  a_species.tau_mean.raw[8] * (x - Topt_species.tau_mean.raw[8])^2 + c_species.tau_mean.raw[8]
}
species9.tau.raw <- function (x){
  a_species.tau_mean.raw[9] * (x - Topt_species.tau_mean.raw[9])^2 + c_species.tau_mean.raw[9]
}

# Make data frame of temperature, fitted value for each species, and species
fitted_data.species.tau.raw <- data.frame(Temperature = rep(temp_seq, n.species.raw), 
                                      Fitted = c(species1.tau.raw(temp_seq),
                                                 species2.tau.raw(temp_seq),
                                                 species3.tau.raw(temp_seq),
                                                 species4.tau.raw(temp_seq),
                                                 species5.tau.raw(temp_seq),
                                                 species6.tau.raw(temp_seq),
                                                 species7.tau.raw(temp_seq),
                                                 species8.tau.raw(temp_seq),
                                                 species9.tau.raw(temp_seq)),
                                      Species = factor(rep(c("Amblyomma lepidum", "Dermacentor andersoni",
                                                             "Dermacentor nitens", "Haemaphysalis leporispalustris",
                                                             "Hyalomma aegyptium", "Hyalomma dromedarii",
                                                             "Hyalomma impeltatum", "Hyalomma lusitanicum",
                                                             "Hyalomma schulzei"), each = length(temp_seq))))
# Data frame of original temperature, trait, and numeric species
plot_data.species.raw <- data.frame(
  Temperature = tick.raw$Temp,
  Trait = tick.raw$Trait,
  Species = as.factor(tick.raw$Species)
)

# Plot observed data and fitted curves with facet wrap by species
species.tau.raw.plot <- ggplot() +
  geom_point(data = plot_data.species.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.species.tau.raw, aes(x = Temperature, y = Fitted, color = Species), size = 1) +
  labs(title = "Observed Data and Quadratic Fit By Species",
       subtitle = "Different Variances",
       x = "Temperature", y = "Oviposition Time in Days") +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~ Species) 
#ggsave("tick.species.tau.plot.png", width = 8, heigh = 5)
species.tau.raw.plot


# Merge the observed and predicted data for comparison
merged_data.species.tau.raw <- fitted_data.species.tau.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species"))

# Calculate RMSE for each species
mse_by_species.tau.raw <- merged_data.species.tau.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.species.tau.raw <- merged_data.species.tau.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_species.tau.raw,
     fitted_data.species.tau.raw,
     merged_data.species.tau.raw,
     species.tau.raw.fit.mcmc,
     species.tau.raw.fit,
     log.lik.species.tau.raw,
     file="species.tau.raw.fit.RData")

################################### By Species Hierarchical Different Var, Same Params ###################################
# Model using species effect
sink("species_tau.param_robin.txt")
cat("
model {
  # Priors
   a ~ dexp(0.1)             # Width
  Topt ~ dnorm(30, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for species effect
  c ~ dexp(0.001)   # Oviposition rate at Topt

  for (j in 1:n.species.raw) {
    species.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.raw) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], species.tau[species.raw[i]]) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], species.tau[species.raw[i]])
  }}
", file = "species_tau.param_robin.txt")


# Parameters to monitor
parameters <- c("a", "c", "tau1", "Topt", "species.tau", "log_lik")


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
jag.data <- list(trait = tick.raw$Trait, N.obs.raw = N.obs.raw, 
                 temp = tick.raw$Temp, 
                 n.species.raw = n.species.raw, species.raw = species.raw)
# Fit model
species.tau.param.raw.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                              model.file = "species_tau.param_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                              n.iter = ni, DIC = T, working.directory = getwd())
# Make an mcmc object
species.tau.param.raw.fit.mcmc <- as.mcmc(species.tau.param.raw.fit) 
closeAllConnections()
species.tau.param.raw.fit$BUGSoutput$DIC

# Extract log-likelihood matrix (samples × observations)
log.lik.species.tau.param.raw <- species.tau.param.raw.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.species.tau.param.raw <- waic(log.lik.species.tau.param.raw)$estimates["waic", "Estimate"]


# Get posterior means for each species effect
summary_fit.species.tau.param.raw <- summary(species.tau.param.raw.fit.mcmc)$statistics
a_species.tau.param_mean.raw <- summary_fit.species.tau.param.raw["a", "Mean"]
c_species.tau.param_mean.raw <- summary_fit.species.tau.param.raw["c", "Mean"]
Topt_species.tau.param_mean.raw <- summary_fit.species.tau.param.raw["Topt", "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp) ,
                max(tick.raw$Temp) , by = 0.1)

# Get predicted values
species.tau.param.raw <- function (x){
  a_species.tau.param_mean.raw * (x - Topt_species.tau.param_mean.raw)^2 + c_species.tau.param_mean.raw
}

# Make data frame of temperature, fitted value for each species, and species
fitted_data.species.tau.param.raw <- data.frame(Temperature = rep(temp_seq, n.species.raw), 
                                            Fitted = species.tau.param.raw(temp_seq),
                                            Species = factor(rep(c("Amblyomma lepidum", "Dermacentor andersoni",
                                                                   "Dermacentor nitens", "Haemaphysalis leporispalustris",
                                                                   "Hyalomma aegyptium", "Hyalomma dromedarii",
                                                                   "Hyalomma impeltatum", "Hyalomma lusitanicum",
                                                                   "Hyalomma schulzei"), each = length(temp_seq))))
# Data frame of original temperature, trait, and numeric species
plot_data.species.raw <- data.frame(
  Temperature = tick.raw$Temp,
  Trait = tick.raw$Trait,
  Species = as.factor(tick.raw$Species)
)

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.species.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.species.tau.param.raw, aes(x = Temperature, y = Fitted, color = Species), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Species",
       x = "Temperature", y = "Oviposition Time in Days") +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~ Species) 

# Merge the observed and predicted data for comparison
merged_data.species.tau.param.raw <- fitted_data.species.tau.param.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species"))

# Calculate RMSE for each species
mse_by_species.tau.param.raw <- merged_data.species.tau.param.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.species.tau.param.raw <- merged_data.species.tau.param.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_species.tau.param.raw,
     fitted_data.species.tau.param.raw,
     merged_data.species.tau.param.raw,
     species.tau.param.raw.fit.mcmc,
     species.tau.param.raw.fit,
     log.lik.species.tau.param.raw,
     file="species.tau.param.raw.fit.RData")

################################### By Species Individual- Lepidum ###################################
a.lepidum.raw <- tick.raw[tick.raw$Species == "Amblyomma lepidum", ]

N.obs.ind.raw <- length(a.lepidum.raw$Trait)

# Model using individual species data
sink("species_individual_robin.txt")
cat("
model {
  # Priors
   a ~ dexp(0.1)             # Width
  Topt ~ dnorm(30, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  c ~ dexp(0.001)      # Oviposition rate at Topt
  
  # Likelihood
  for (i in 1:N.obs.ind.raw) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], tau1) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], tau1)
  }}
", file = "species_individual_robin.txt")


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
jag.data <- list(trait = a.lepidum.raw$Trait, N.obs.ind.raw = N.obs.ind.raw, 
                 temp = a.lepidum.raw$Temp)

# Fit fill model
a.lepidum.raw.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                       model.file = "species_individual_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                       n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
a.lepidum.raw.fit.mcmc <- as.mcmc(a.lepidum.raw.fit) 
closeAllConnections()
#DIC
a.lepidum.raw.fit$BUGSoutput$DIC

# Get posterior means
summary_fit.a.lepidum.raw <- summary(a.lepidum.raw.fit.mcmc)$statistics
a_mean.raw <- summary_fit.a.lepidum.raw["a", "Mean"]
c_mean.raw <- summary_fit.a.lepidum.raw["c", "Mean"]
Topt_mean.raw <- summary_fit.a.lepidum.raw["Topt", "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp), max(tick.raw$Temp), by = 0.1)
# 
# # Get predicted values
h.schulzei.fitted.raw <- function (x){
  a_mean.raw * (x - Topt_mean.raw)^2 + c_mean.raw
}


# Data frame of temperature, fitted value, species
fitted_data.a.lepidum.raw <- data.frame(Temperature = rep(temp_seq),
                                     Fitted = h.schulzei.fitted.raw(temp_seq))
# Original data with temperature, trait, and species as numeric
plot_data.a.lepidum.raw <- data.frame(
  Temperature = a.lepidum.raw$Temp,
  Trait = a.lepidum.raw$Trait
)

# Plot observed data and fitted curves 
p.a.lepidum.raw <- ggplot() +
  geom_point(data = plot_data.a.lepidum.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.a.lepidum.raw, aes(x = Temperature, y = Fitted), color = "seagreen1", size = 1) +
  labs(title = "Observed Data and Quadratic Fit",
       x = "Temperature", y = "Oviposition Time in Days") +
  theme_minimal() +
  theme(legend.position = "none")
p.a.lepidum.raw

# Merge the observed and predicted data for comparison
merged_data.a.lepidum.raw <- fitted_data.a.lepidum.raw %>%
  left_join(a.lepidum.raw, by = c("Temperature" = "Temp"))

# Calculate RMSE for each species
rmse.a.lepidum.raw <- merged_data.a.lepidum.raw %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE),
    SE = (Trait - Fitted)^2 
  )


# Extract log-likelihood matrix (samples × observations)
log.lik.a.lepidum.raw <- a.lepidum.raw.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.a.lepidum.raw <- waic(log.lik.a.lepidum.raw)$estimates["waic", "Estimate"]


save(rmse.a.lepidum.raw,
     fitted_data.a.lepidum.raw,
     merged_data.a.lepidum.raw,
     a.lepidum.raw.fit.mcmc,
     a.lepidum.raw.fit,
     p.a.lepidum.raw,
     waic.a.lepidum.raw,
     file="a.lepidum.raw.fit.RData")

################################### By Species Individual- Andersoni ###################################
d.andersoni.raw <- tick.raw[tick.raw$Species == "Dermacentor andersoni", ]


N.obs.ind.raw <- length(d.andersoni.raw$Trait)

# Model using individual species data
sink("species_individual_robin.txt")
cat("
model {
  # Priors
   a ~ dexp(0.1)             # Width
  Topt ~ dnorm(30, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  c ~ dexp(0.001)      # Oviposition rate at Topt
  
  # Likelihood
  for (i in 1:N.obs.ind.raw) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], tau1) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], tau1)
  }}
", file = "species_individual_robin.txt")


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
jag.data <- list(trait = d.andersoni.raw$Trait, N.obs.ind.raw = N.obs.ind.raw, 
                 temp = d.andersoni.raw$Temp)

# Fit fill model
d.andersoni.raw.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                           model.file = "species_individual_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                           n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
d.andersoni.raw.fit.mcmc <- as.mcmc(d.andersoni.raw.fit) 
closeAllConnections()
#DIC
d.andersoni.raw.fit$BUGSoutput$DIC

# Get posterior means
summary_fit.d.andersoni.raw <- summary(d.andersoni.raw.fit.mcmc)$statistics
a_mean.raw <- summary_fit.d.andersoni.raw["a", "Mean"]
c_mean.raw <- summary_fit.d.andersoni.raw["c", "Mean"]
Topt_mean.raw <- summary_fit.d.andersoni.raw["Topt", "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp), max(tick.raw$Temp), by = 0.1)
# 
# # Get predicted values
h.schulzei.fitted.raw <- function (x){
  a_mean.raw * (x - Topt_mean.raw)^2 + c_mean.raw
}


# Data frame of temperature, fitted value, species
fitted_data.d.andersoni.raw <- data.frame(Temperature = rep(temp_seq),
                                         Fitted = h.schulzei.fitted.raw(temp_seq))
# Original data with temperature, trait, and species as numeric
plot_data.d.andersoni.raw <- data.frame(
  Temperature = d.andersoni.raw$Temp,
  Trait = d.andersoni.raw$Trait
)

# Plot observed data and fitted curves 
p.d.andersoni.raw <- ggplot() +
  geom_point(data = plot_data.d.andersoni.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.d.andersoni.raw, aes(x = Temperature, y = Fitted), color = "seagreen1", size = 1) +
  labs(title = "Observed Data and Quadratic Fit",
       x = "Temperature", y = "Oviposition Time in Days") +
  theme_minimal() +
  theme(legend.position = "none")
p.d.andersoni.raw

# Merge the observed and predicted data for comparison
merged_data.d.andersoni.raw <- fitted_data.d.andersoni.raw %>%
  left_join(d.andersoni.raw, by = c("Temperature" = "Temp"))

# Calculate RMSE for each species
rmse.d.andersoni.raw <- merged_data.d.andersoni.raw %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE),
    SE = (Trait - Fitted)^2 
  )


# Extract log-likelihood matrix (samples × observations)
log.lik.d.andersoni.raw <- d.andersoni.raw.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.d.andersoni.raw <- waic(log.lik.d.andersoni.raw)$estimates["waic", "Estimate"]


save(rmse.d.andersoni.raw,
     fitted_data.d.andersoni.raw,
     merged_data.d.andersoni.raw,
     d.andersoni.raw.fit.mcmc,
     d.andersoni.raw.fit,
     p.d.andersoni.raw,
     waic.d.andersoni.raw,
     file="d.andersoni.raw.fit.RData")

################################### By Species Individual- Nitens ###################################
d.nitens.raw <- tick.raw[tick.raw$Species == "Dermacentor nitens", ]

N.obs.ind.raw <- length(d.nitens.raw$Trait)

# Model using individual species data
sink("species_individual_robin.txt")
cat("
model {
  # Priors
   a ~ dexp(0.1)             # Width
  Topt ~ dnorm(30, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  c ~ dexp(0.001)      # Oviposition rate at Topt
  
  # Likelihood
  for (i in 1:N.obs.ind.raw) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], tau1) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], tau1)
  }}
", file = "species_individual_robin.txt")


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
jag.data <- list(trait = d.nitens.raw$Trait, N.obs.ind.raw = N.obs.ind.raw, 
                 temp = d.nitens.raw$Temp)

# Fit fill model
d.nitens.raw.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                           model.file = "species_individual_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                           n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
d.nitens.raw.fit.mcmc <- as.mcmc(d.nitens.raw.fit) 
closeAllConnections()
#DIC
d.nitens.raw.fit$BUGSoutput$DIC

# Get posterior means
summary_fit.d.nitens.raw <- summary(d.nitens.raw.fit.mcmc)$statistics
a_mean.raw <- summary_fit.d.nitens.raw["a", "Mean"]
c_mean.raw <- summary_fit.d.nitens.raw["c", "Mean"]
Topt_mean.raw <- summary_fit.d.nitens.raw["Topt", "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp), max(tick.raw$Temp), by = 0.1)
# 
# # Get predicted values
h.schulzei.fitted.raw <- function (x){
  a_mean.raw * (x - Topt_mean.raw)^2 + c_mean.raw
}


# Data frame of temperature, fitted value, species
fitted_data.d.nitens.raw <- data.frame(Temperature = rep(temp_seq),
                                         Fitted = h.schulzei.fitted.raw(temp_seq))
# Original data with temperature, trait, and species as numeric
plot_data.d.nitens.raw <- data.frame(
  Temperature = d.nitens.raw$Temp,
  Trait = d.nitens.raw$Trait
)

# Plot observed data and fitted curves 
p.d.nitens.raw <- ggplot() +
  geom_point(data = plot_data.d.nitens.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.d.nitens.raw, aes(x = Temperature, y = Fitted), color = "seagreen1", size = 1) +
  labs(title = "Observed Data and Quadratic Fit",
       x = "Temperature", y = "Oviposition Time in Days") +
  theme_minimal() +
  theme(legend.position = "none")
p.d.nitens.raw

# Merge the observed and predicted data for comparison
merged_data.d.nitens.raw <- fitted_data.d.nitens.raw %>%
  left_join(d.nitens.raw, by = c("Temperature" = "Temp"))

# Calculate RMSE for each species
rmse.d.nitens.raw <- merged_data.d.nitens.raw %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE),
    SE = (Trait - Fitted)^2 
  )


# Extract log-likelihood matrix (samples × observations)
log.lik.d.nitens.raw <- d.nitens.raw.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.d.nitens.raw <- waic(log.lik.d.nitens.raw)$estimates["waic", "Estimate"]


save(rmse.d.nitens.raw,
     fitted_data.d.nitens.raw,
     merged_data.d.nitens.raw,
     d.nitens.raw.fit.mcmc,
     d.nitens.raw.fit,
     p.d.nitens.raw,
     waic.d.nitens.raw,
     file="d.nitens.raw.fit.RData")

################################### By Species Individual- Leporispalustris ###################################
h.leporispalustris.raw <- tick.raw[tick.raw$Species == "Haemaphysalis leporispalustris", ]

N.obs.ind.raw <- length(h.leporispalustris.raw$Trait)

# Model using individual species data
sink("species_individual_robin.txt")
cat("
model {
  # Priors
   a ~ dexp(0.1)             # Width
  Topt ~ dnorm(30, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  c ~ dexp(0.001)      # Oviposition rate at Topt
  
  # Likelihood
  for (i in 1:N.obs.ind.raw) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], tau1) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], tau1)
  }}
", file = "species_individual_robin.txt")


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
jag.data <- list(trait = h.leporispalustris.raw$Trait, N.obs.ind.raw = N.obs.ind.raw, 
                 temp = h.leporispalustris.raw$Temp)

# Fit fill model
h.leporispalustris.raw.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                           model.file = "species_individual_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                           n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
h.leporispalustris.raw.fit.mcmc <- as.mcmc(h.leporispalustris.raw.fit) 
closeAllConnections()
#DIC
h.leporispalustris.raw.fit$BUGSoutput$DIC

# Get posterior means
summary_fit.h.leporispalustris.raw <- summary(h.leporispalustris.raw.fit.mcmc)$statistics
a_mean.raw <- summary_fit.h.leporispalustris.raw["a", "Mean"]
c_mean.raw <- summary_fit.h.leporispalustris.raw["c", "Mean"]
Topt_mean.raw <- summary_fit.h.leporispalustris.raw["Topt", "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp), max(tick.raw$Temp), by = 0.1)
# 
# # Get predicted values
h.schulzei.fitted.raw <- function (x){
  a_mean.raw * (x - Topt_mean.raw)^2 + c_mean.raw
}


# Data frame of temperature, fitted value, species
fitted_data.h.leporispalustris.raw <- data.frame(Temperature = rep(temp_seq),
                                         Fitted = h.schulzei.fitted.raw(temp_seq))
# Original data with temperature, trait, and species as numeric
plot_data.h.leporispalustris.raw <- data.frame(
  Temperature = h.leporispalustris.raw$Temp,
  Trait = h.leporispalustris.raw$Trait
)

# Plot observed data and fitted curves 
p.h.leporispalustris.raw <- ggplot() +
  geom_point(data = plot_data.h.leporispalustris.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.h.leporispalustris.raw, aes(x = Temperature, y = Fitted), color = "seagreen1", size = 1) +
  labs(title = "Observed Data and Quadratic Fit",
       x = "Temperature", y = "Oviposition Time in Days") +
  theme_minimal() +
  theme(legend.position = "none")
p.h.leporispalustris.raw

# Merge the observed and predicted data for comparison
merged_data.h.leporispalustris.raw <- fitted_data.h.leporispalustris.raw %>%
  left_join(h.leporispalustris.raw, by = c("Temperature" = "Temp"))

# Calculate RMSE for each species
rmse.h.leporispalustris.raw <- merged_data.h.leporispalustris.raw %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE),
    SE = (Trait - Fitted)^2 
  )


# Extract log-likelihood matrix (samples × observations)
log.lik.h.leporispalustris.raw <- h.leporispalustris.raw.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.h.leporispalustris.raw <- waic(log.lik.h.leporispalustris.raw)$estimates["waic", "Estimate"]


save(rmse.h.leporispalustris.raw,
     fitted_data.h.leporispalustris.raw,
     merged_data.h.leporispalustris.raw,
     h.leporispalustris.raw.fit.mcmc,
     h.leporispalustris.raw.fit,
     p.h.leporispalustris.raw,
     waic.h.leporispalustris.raw,
     file="h.leporispalustris.raw.fit.RData")

################################### By Species Individual- Aegyptium ###################################
h.aegyptium.raw <- tick.raw[tick.raw$Species == "Hyalomma aegyptium", ]

N.obs.ind.raw <- length(h.aegyptium.raw$Trait)

# Model using individual species data
sink("species_individual_robin.txt")
cat("
model {
  # Priors
   a ~ dexp(0.1)             # Width
  Topt ~ dnorm(30, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  c ~ dexp(0.001)      # Oviposition rate at Topt
  
  # Likelihood
  for (i in 1:N.obs.ind.raw) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], tau1) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], tau1)
  }}
", file = "species_individual_robin.txt")


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
jag.data <- list(trait = h.aegyptium.raw$Trait, N.obs.ind.raw = N.obs.ind.raw, 
                 temp = h.aegyptium.raw$Temp)

# Fit fill model
h.aegyptium.raw.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                           model.file = "species_individual_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                           n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
h.aegyptium.raw.fit.mcmc <- as.mcmc(h.aegyptium.raw.fit) 
closeAllConnections()
#DIC
h.aegyptium.raw.fit$BUGSoutput$DIC

# Get posterior means
summary_fit.h.aegyptium.raw <- summary(h.aegyptium.raw.fit.mcmc)$statistics
a_mean.raw <- summary_fit.h.aegyptium.raw["a", "Mean"]
c_mean.raw <- summary_fit.h.aegyptium.raw["c", "Mean"]
Topt_mean.raw <- summary_fit.h.aegyptium.raw["Topt", "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp), max(tick.raw$Temp), by = 0.1)
# 
# # Get predicted values
h.schulzei.fitted.raw <- function (x){
  a_mean.raw * (x - Topt_mean.raw)^2 + c_mean.raw
}


# Data frame of temperature, fitted value, species
fitted_data.h.aegyptium.raw <- data.frame(Temperature = rep(temp_seq),
                                         Fitted = h.schulzei.fitted.raw(temp_seq))
# Original data with temperature, trait, and species as numeric
plot_data.h.aegyptium.raw <- data.frame(
  Temperature = h.aegyptium.raw$Temp,
  Trait = h.aegyptium.raw$Trait
)

# Plot observed data and fitted curves 
p.h.aegyptium.raw <- ggplot() +
  geom_point(data = plot_data.h.aegyptium.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.h.aegyptium.raw, aes(x = Temperature, y = Fitted), color = "seagreen1", size = 1) +
  labs(title = "Observed Data and Quadratic Fit",
       x = "Temperature", y = "Oviposition Time in Days") +
  theme_minimal() +
  theme(legend.position = "none")
p.h.aegyptium.raw

# Merge the observed and predicted data for comparison
merged_data.h.aegyptium.raw <- fitted_data.h.aegyptium.raw %>%
  left_join(h.aegyptium.raw, by = c("Temperature" = "Temp"))

# Calculate RMSE for each species
rmse.h.aegyptium.raw <- merged_data.h.aegyptium.raw %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE),
    SE = (Trait - Fitted)^2 
  )


# Extract log-likelihood matrix (samples × observations)
log.lik.h.aegyptium.raw <- h.aegyptium.raw.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.h.aegyptium.raw <- waic(log.lik.h.aegyptium.raw)$estimates["waic", "Estimate"]

save(rmse.h.aegyptium.raw,
     fitted_data.h.aegyptium.raw,
     merged_data.h.aegyptium.raw,
     h.aegyptium.raw.fit.mcmc,
     h.aegyptium.raw.fit,
     p.h.aegyptium.raw,
     waic.h.aegyptium.raw,
     file="h.aegyptium.raw.fit.RData")

################################### By Species Individual- Dromdarii ###################################
h.dromedarii.raw <- tick.raw[tick.raw$Species == "Hyalomma dromedarii", ]

N.obs.ind.raw <- length(h.dromedarii.raw$Trait)

# Model using individual species data
sink("species_individual_robin.txt")
cat("
model {
  # Priors
   a ~ dexp(0.1)             # Width
  Topt ~ dnorm(30, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  c ~ dexp(0.001)      # Oviposition rate at Topt
  
  # Likelihood
  for (i in 1:N.obs.ind.raw) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], tau1) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], tau1)
  }}
", file = "species_individual_robin.txt")


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
jag.data <- list(trait = h.dromedarii.raw$Trait, N.obs.ind.raw = N.obs.ind.raw, 
                 temp = h.dromedarii.raw$Temp)

# Fit fill model
h.dromedarii.raw.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                           model.file = "species_individual_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                           n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
h.dromedarii.raw.fit.mcmc <- as.mcmc(h.dromedarii.raw.fit) 
closeAllConnections()
#DIC
h.dromedarii.raw.fit$BUGSoutput$DIC

# Get posterior means
summary_fit.h.dromedarii.raw <- summary(h.dromedarii.raw.fit.mcmc)$statistics
a_mean.raw <- summary_fit.h.dromedarii.raw["a", "Mean"]
c_mean.raw <- summary_fit.h.dromedarii.raw["c", "Mean"]
Topt_mean.raw <- summary_fit.h.dromedarii.raw["Topt", "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp), max(tick.raw$Temp), by = 0.1)
# 
# # Get predicted values
h.schulzei.fitted.raw <- function (x){
  a_mean.raw * (x - Topt_mean.raw)^2 + c_mean.raw
}


# Data frame of temperature, fitted value, species
fitted_data.h.dromedarii.raw <- data.frame(Temperature = rep(temp_seq),
                                         Fitted = h.schulzei.fitted.raw(temp_seq))
# Original data with temperature, trait, and species as numeric
plot_data.h.dromedarii.raw <- data.frame(
  Temperature = h.dromedarii.raw$Temp,
  Trait = h.dromedarii.raw$Trait
)

# Plot observed data and fitted curves 
p.h.dromedarii.raw <- ggplot() +
  geom_point(data = plot_data.h.dromedarii.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.h.dromedarii.raw, aes(x = Temperature, y = Fitted), color = "seagreen1", size = 1) +
  labs(title = "Observed Data and Quadratic Fit",
       x = "Temperature", y = "Oviposition Time in Days") +
  theme_minimal() +
  theme(legend.position = "none")
p.h.dromedarii.raw

# Merge the observed and predicted data for comparison
merged_data.h.dromedarii.raw <- fitted_data.h.dromedarii.raw %>%
  left_join(h.dromedarii.raw, by = c("Temperature" = "Temp"))

# Calculate RMSE for each species
rmse.h.dromedarii.raw <- merged_data.h.dromedarii.raw %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE),
    SE = (Trait - Fitted)^2 
  )


# Extract log-likelihood matrix (samples × observations)
log.lik.h.dromedarii.raw <- h.dromedarii.raw.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.h.dromedarii.raw <- waic(log.lik.h.dromedarii.raw)$estimates["waic", "Estimate"]

save(rmse.h.dromedarii.raw,
     fitted_data.h.dromedarii.raw,
     merged_data.h.dromedarii.raw,
     h.dromedarii.raw.fit.mcmc,
     h.dromedarii.raw.fit,
     p.h.dromedarii.raw,
     waic.h.dromedarii.raw,
     file="h.dromedarii.raw.fit.RData")

################################### By Species Individual- Impeltatum ###################################
h.impeltatum.raw <- tick.raw[tick.raw$Species == "Hyalomma impeltatum", ]

N.obs.ind.raw <- length(h.impeltatum.raw$Trait)

# Model using individual species data
sink("species_individual_robin.txt")
cat("
model {
  # Priors
   a ~ dexp(0.1)             # Width
  Topt ~ dnorm(30, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  c ~ dexp(0.001)      # Oviposition rate at Topt
  
  # Likelihood
  for (i in 1:N.obs.ind.raw) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], tau1) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], tau1)
  }}
", file = "species_individual_robin.txt")


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
jag.data <- list(trait = h.impeltatum.raw$Trait, N.obs.ind.raw = N.obs.ind.raw, 
                 temp = h.impeltatum.raw$Temp)

# Fit fill model
h.impeltatum.raw.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                           model.file = "species_individual_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                           n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
h.impeltatum.raw.fit.mcmc <- as.mcmc(h.impeltatum.raw.fit) 
closeAllConnections()
#DIC
h.impeltatum.raw.fit$BUGSoutput$DIC

# Get posterior means
summary_fit.h.impeltatum.raw <- summary(h.impeltatum.raw.fit.mcmc)$statistics
a_mean.raw <- summary_fit.h.impeltatum.raw["a", "Mean"]
c_mean.raw <- summary_fit.h.impeltatum.raw["c", "Mean"]
Topt_mean.raw <- summary_fit.h.impeltatum.raw["Topt", "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp), max(tick.raw$Temp), by = 0.1)
# 
# # Get predicted values
h.schulzei.fitted.raw <- function (x){
  a_mean.raw * (x - Topt_mean.raw)^2 + c_mean.raw
}


# Data frame of temperature, fitted value, species
fitted_data.h.impeltatum.raw <- data.frame(Temperature = rep(temp_seq),
                                         Fitted = h.schulzei.fitted.raw(temp_seq))
# Original data with temperature, trait, and species as numeric
plot_data.h.impeltatum.raw <- data.frame(
  Temperature = h.impeltatum.raw$Temp,
  Trait = h.impeltatum.raw$Trait
)

# Plot observed data and fitted curves 
p.h.impeltatum.raw <- ggplot() +
  geom_point(data = plot_data.h.impeltatum.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.h.impeltatum.raw, aes(x = Temperature, y = Fitted), color = "seagreen1", size = 1) +
  labs(title = "Observed Data and Quadratic Fit",
       x = "Temperature", y = "Oviposition Time in Days") +
  theme_minimal() +
  theme(legend.position = "none")
p.h.impeltatum.raw

# Merge the observed and predicted data for comparison
merged_data.h.impeltatum.raw <- fitted_data.h.impeltatum.raw %>%
  left_join(h.impeltatum.raw, by = c("Temperature" = "Temp"))

# Calculate RMSE for each species
rmse.h.impeltatum.raw <- merged_data.h.impeltatum.raw %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE),
    SE = (Trait - Fitted)^2 
  )


# Extract log-likelihood matrix (samples × observations)
log.lik.h.impeltatum.raw <- h.impeltatum.raw.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.h.impeltatum.raw <- waic(log.lik.h.impeltatum.raw)$estimates["waic", "Estimate"]


save(rmse.h.impeltatum.raw,
     fitted_data.h.impeltatum.raw,
     merged_data.h.impeltatum.raw,
     h.impeltatum.raw.fit.mcmc,
     h.impeltatum.raw.fit,
     p.h.impeltatum.raw,
     waic.h.impeltatum.raw,
     file="h.impeltatum.raw.fit.RData")

################################### By Species Individual- Lusitanicum ###################################
h.lusitanicum.raw <- tick.raw[tick.raw$Species == "Hyalomma lusitanicum", ]

N.obs.ind.raw <- length(h.lusitanicum.raw$Trait)

# Model using individual species data
sink("species_individual_robin.txt")
cat("
model {
  # Priors
   a ~ dexp(0.1)             # Width
  Topt ~ dnorm(30, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  c ~ dexp(0.001)      # Oviposition rate at Topt
  
  # Likelihood
  for (i in 1:N.obs.ind.raw) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], tau1) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], tau1)
  }}
", file = "species_individual_robin.txt")


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
jag.data <- list(trait = h.lusitanicum.raw$Trait, N.obs.ind.raw = N.obs.ind.raw, 
                 temp = h.lusitanicum.raw$Temp)

# Fit fill model
h.lusitanicum.raw.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                           model.file = "species_individual_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                           n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
h.lusitanicum.raw.fit.mcmc <- as.mcmc(h.lusitanicum.raw.fit) 
closeAllConnections()
#DIC
h.lusitanicum.raw.fit$BUGSoutput$DIC

# Get posterior means
summary_fit.h.lusitanicum.raw <- summary(h.lusitanicum.raw.fit.mcmc)$statistics
a_mean.raw <- summary_fit.h.lusitanicum.raw["a", "Mean"]
c_mean.raw <- summary_fit.h.lusitanicum.raw["c", "Mean"]
Topt_mean.raw <- summary_fit.h.lusitanicum.raw["Topt", "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp), max(tick.raw$Temp), by = 0.1)
# 
# # Get predicted values
h.schulzei.fitted.raw <- function (x){
  a_mean.raw * (x - Topt_mean.raw)^2 + c_mean.raw
}


# Data frame of temperature, fitted value, species
fitted_data.h.lusitanicum.raw <- data.frame(Temperature = rep(temp_seq),
                                         Fitted = h.schulzei.fitted.raw(temp_seq))
# Original data with temperature, trait, and species as numeric
plot_data.h.lusitanicum.raw <- data.frame(
  Temperature = h.lusitanicum.raw$Temp,
  Trait = h.lusitanicum.raw$Trait
)

# Plot observed data and fitted curves 
p.h.lusitanicum.raw <- ggplot() +
  geom_point(data = plot_data.h.lusitanicum.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.h.lusitanicum.raw, aes(x = Temperature, y = Fitted), color = "seagreen1", size = 1) +
  labs(title = "Observed Data and Quadratic Fit",
       x = "Temperature", y = "Oviposition Time in Days") +
  theme_minimal() +
  theme(legend.position = "none")
p.h.lusitanicum.raw

# Merge the observed and predicted data for comparison
merged_data.h.lusitanicum.raw <- fitted_data.h.lusitanicum.raw %>%
  left_join(h.lusitanicum.raw, by = c("Temperature" = "Temp"))

# Calculate RMSE for each species
rmse.h.lusitanicum.raw <- merged_data.h.lusitanicum.raw %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE),
    SE = (Trait - Fitted)^2 
  )


# Extract log-likelihood matrix (samples × observations)
log.lik.h.lusitanicum.raw <- h.lusitanicum.raw.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.h.lusitanicum.raw <- waic(log.lik.h.lusitanicum.raw)$estimates["waic", "Estimate"]

save(rmse.h.lusitanicum.raw,
     fitted_data.h.lusitanicum.raw,
     merged_data.h.lusitanicum.raw,
     h.lusitanicum.raw.fit.mcmc,
     h.lusitanicum.raw.fit,
     p.h.lusitanicum.raw,
     waic.h.lusitanicum.raw,
     file="h.lusitanicum.raw.fit.RData")

################################### By Species Individual- Schulzei ###################################
h.schulzei.raw <- tick.raw[tick.raw$Species == "Hyalomma schulzei", ]

N.obs.ind.raw <- length(h.schulzei.raw$Trait)

# Model using individual species data
sink("species_individual_robin.txt")
cat("
model {
  # Priors
   a ~ dexp(0.1)             # Width
  Topt ~ dnorm(30, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  c ~ dexp(0.001)      # Oviposition rate at Topt
  
  # Likelihood
  for (i in 1:N.obs.ind.raw) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], tau1) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], tau1)
  }}
", file = "species_individual_robin.txt")


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
jag.data <- list(trait = h.schulzei.raw$Trait, N.obs.ind.raw = N.obs.ind.raw, 
                 temp = h.schulzei.raw$Temp)

# Fit fill model
h.schulzei.raw.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                           model.file = "species_individual_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                           n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
h.schulzei.raw.fit.mcmc <- as.mcmc(h.schulzei.raw.fit) 
closeAllConnections()
#DIC
h.schulzei.raw.fit$BUGSoutput$DIC

# Get posterior means
summary_fit.h.schulzei.raw <- summary(h.schulzei.raw.fit.mcmc)$statistics
a_mean.raw <- summary_fit.h.schulzei.raw["a", "Mean"]
c_mean.raw <- summary_fit.h.schulzei.raw["c", "Mean"]
Topt_mean.raw <- summary_fit.h.schulzei.raw["Topt", "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp), max(tick.raw$Temp), by = 0.1)
# 
# # Get predicted values
h.schulzei.fitted.raw <- function (x){
  a_mean.raw * (x - Topt_mean.raw)^2 + c_mean.raw
}


# Data frame of temperature, fitted value, species
fitted_data.h.schulzei.raw <- data.frame(Temperature = rep(temp_seq),
                                         Fitted = h.schulzei.fitted.raw(temp_seq))
# Original data with temperature, trait, and species as numeric
plot_data.h.schulzei.raw <- data.frame(
  Temperature = h.schulzei.raw$Temp,
  Trait = h.schulzei.raw$Trait
)

# Plot observed data and fitted curves 
p.h.schulzei.raw <- ggplot() +
  geom_point(data = plot_data.h.schulzei.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.h.schulzei.raw, aes(x = Temperature, y = Fitted), color = "seagreen1", size = 1) +
  labs(title = "Observed Data and Quadratic Fit",
       x = "Temperature", y = "Oviposition Time in Days") +
  theme_minimal() +
  theme(legend.position = "none")
p.h.schulzei.raw

# Merge the observed and predicted data for comparison
merged_data.h.schulzei.raw <- fitted_data.h.schulzei.raw %>%
  left_join(h.schulzei.raw, by = c("Temperature" = "Temp"))

# Calculate RMSE for each species
rmse.h.schulzei.raw <- merged_data.h.schulzei.raw %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE),
    SE = (Trait - Fitted)^2 
  )


# Extract log-likelihood matrix (samples × observations)
log.lik.h.schulzei.raw <- h.schulzei.raw.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.h.schulzei.raw <- waic(log.lik.h.schulzei.raw)$estimates["waic", "Estimate"]


save(rmse.h.schulzei.raw,
     fitted_data.h.schulzei.raw,
     merged_data.h.schulzei.raw,
     h.schulzei.raw.fit.mcmc,
     h.schulzei.raw.fit,
     p.h.schulzei.raw,
     waic.h.schulzei.raw,
     file="h.schulzei.raw.fit.RData")

################################### By Species Individual- Summary ###################################
waic.species.ind.raw <- waic(log.lik.a.lepidum.raw)$estimates["waic", "Estimate"] + 
  waic(log.lik.d.andersoni.raw)$estimates["waic", "Estimate"] +
  waic(log.lik.d.nitens.raw)$estimates["waic", "Estimate"] +
  waic(log.lik.h.leporispalustris.raw)$estimates["waic", "Estimate"] + 
  waic(log.lik.h.aegyptium.raw)$estimates["waic", "Estimate"] +
  waic(log.lik.h.dromedarii.raw)$estimates["waic", "Estimate"] + 
  waic(log.lik.h.impeltatum.raw)$estimates["waic", "Estimate"] +
  waic(log.lik.h.lusitanicum.raw)$estimates["waic", "Estimate"] +
  waic(log.lik.h.schulzei.raw)$estimates["waic", "Estimate"]

DIC.species.individual.raw <- a.lepidum.raw.fit$BUGSoutput$DIC + d.andersoni.raw.fit$BUGSoutput$DIC + 
  d.nitens.raw.fit$BUGSoutput$DIC +
  h.leporispalustris.raw.fit$BUGSoutput$DIC + h.aegyptium.raw.fit$BUGSoutput$DIC + 
  h.dromedarii.raw.fit$BUGSoutput$DIC +
  h.impeltatum.raw.fit$BUGSoutput$DIC + h.lusitanicum.raw.fit$BUGSoutput$DIC + 
  h.schulzei.raw.fit$BUGSoutput$DIC

# Define RMSE values for each species
se.individual.raw <- c(rmse.a.lepidum.raw$SE, rmse.d.andersoni.raw$SE, rmse.d.nitens.raw$SE,
                       rmse.h.leporispalustris.raw$SE, rmse.h.aegyptium.raw$SE,
                       rmse.h.dromedarii.raw$SE, rmse.h.impeltatum.raw$SE,
                       rmse.h.lusitanicum.raw$SE, rmse.h.schulzei.raw$SE)

# Compute weighted average RMSE
overall_rmse.species.ind.raw <- sqrt(mean(se.individual.raw, na.rm = TRUE))

# Arrange plots into a 3x3 grid
grid.ind.species.raw <- grid.arrange(
  p.a.lepidum.raw + labs(title = "Amblyomma lepidum", x = " ", y = "Oviposition in Days"),
  p.d.andersoni.raw + labs(title = "Dermacentor andersoni", x = " ", y = " "),
  p.d.nitens.raw + labs(title = "Dermacentor nitens", x = " ", y = " "),
  p.h.leporispalustris.raw + labs(title = "Haemaphysalis leporispalustris", x = " ", y = "Oviposition in Days"),
  p.h.aegyptium.raw + labs(title = "Hyalomma aegyptium", x = " ", y = " "),
  p.h.dromedarii.raw + labs(title = "Hyalomma dromedarii", x = " ", y = " "),
  p.h.impeltatum.raw + labs(title = "Hyalomma impeltatum", y = "Oviposition in Days"),
  p.h.lusitanicum.raw + labs(title = "Hyalomma lusitanicum", y = " "),
  p.h.schulzei.raw + labs(title = "Hyalomma schulzei", y = " "),
  nrow = 3, ncol = 3
)
#ggsave("tick.species.ind.plot.png", grid.ind.species, width = 8, height = 5)


################################### By Genus ###################################
# Model using genus effect
sink("genus_robin.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Genus-specific 
  mu.c ~ dexp(10)          # Mean for genus effects c (vertex y)
  mu.Topt ~ dnorm(30, 0.1)          # Mean for genus effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for genus effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for genus effects a (width)

  for (j in 1:n.genus.raw) {
    mu.genus.c[j] ~ dexp(mu.c)       # Genus-specific vertex y
    mu.genus.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Genus-specific vertex x
    mu.genus.a[j] ~ dexp(mu.a)    # Genus-specific width
    genus.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.raw) {
    mu[i] <- mu.genus.a[genus.raw[i]] * (mu.genus.Topt[genus.raw[i]] - temp[i])^2 + mu.genus.c[genus.raw[i]]
    trait[i] ~ dnorm(mu[i], genus.tau[genus.raw[i]]) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], genus.tau[genus.raw[i]])
  }}
", file = "genus_robin.txt")


# Parameters to monitor
parameters <- c("mu.genus.a", "mu.genus.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt", "genus.tau",
                "tau1", "mu.genus.Topt", "log_lik")


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
jag.data <- list(trait = tick.raw$Trait, N.obs.raw = N.obs.raw, 
                 temp = tick.raw$Temp, 
                 n.genus.raw = n.genus.raw, genus.raw = genus.raw)
# Fit model
genus.raw.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                  model.file = "genus_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                  n.iter = ni, DIC = T, working.directory = getwd())
# Make an mcmc object
genus.raw.fit.mcmc <- as.mcmc(genus.raw.fit)
closeAllConnections()
genus.raw.fit$BUGSoutput$DIC

# Extract log-likelihood matrix (samples × observations)
log.lik.genus.raw <- genus.raw.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.genus.raw <- waic(log.lik.genus.raw)$estimates["waic", "Estimate"]




# Get posterior means using parameters from genus effect
summary_fit.genus.raw <- summary(genus.raw.fit.mcmc)$statistics
a_genus_mean.raw <- summary_fit.genus.raw[grep("^mu.genus.a\\[", rownames(summary_fit.genus.raw)), "Mean"]
c_genus_mean.raw <- summary_fit.genus.raw[grep("^mu.genus.c\\[", rownames(summary_fit.genus.raw)), "Mean"]
Topt_genus_mean.raw <- summary_fit.genus.raw[grep("^mu.genus.Topt\\[", rownames(summary_fit.genus.raw)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp) ,
                max(tick.raw$Temp) , by = 0.1)

# Get predicted values for each genus
genus1.raw <- function (x){
  a_genus_mean.raw[1] * (x - Topt_genus_mean.raw[1])^2 + c_genus_mean.raw[1]
}
genus2.raw <- function (x){
  a_genus_mean.raw[2] * (x - Topt_genus_mean.raw[2])^2 + c_genus_mean.raw[2]
}
genus3.raw <- function (x){
  a_genus_mean.raw[3] * (x - Topt_genus_mean.raw[3])^2 + c_genus_mean.raw[3]
}
genus4.raw <- function (x){
  a_genus_mean.raw[4] * (x - Topt_genus_mean.raw[4])^2 + c_genus_mean.raw[4]
}

# Data frame of temperature, fitted values, and genus
fitted_data.genus.raw <- data.frame(Temperature = rep(temp_seq, n.genus.raw), 
                                Fitted = c(genus1.raw(temp_seq),
                                           genus2.raw(temp_seq),
                                           genus3.raw(temp_seq),
                                           genus4.raw(temp_seq)),
                                Genus = factor(rep(c("Amblyomma", "Dermacentor",
                                                     "Haemaphysalis", "Hyalomma"), each = length(temp_seq))))
# Original temperature, trait, numeric genus, and numeric species to plot by species
plot_data.genus.raw <- data.frame(
  Temperature = tick.raw$Temp,
  Trait = tick.raw$Trait,
  Genus = as.factor(tick.raw$Genus),
  Species = as.factor(tick.raw$Species)
)

# Get unique combos of species and genus to align the correct curve on the plot
species_to_genus.raw <- unique(plot_data.genus.raw[, c("Species", "Genus")])
species_to_genus.raw$Genus <- as.factor(species_to_genus.raw$Genus)
# Merge fitted genus curves with species-to-genus mapping
fitted_data.genus.raw <- data.frame(
  Temperature = rep(temp_seq, n.genus.raw),
  Fitted = c(genus1.raw(temp_seq),
             genus2.raw(temp_seq),
             genus3.raw(temp_seq),
             genus4.raw(temp_seq)),
  Genus = factor(rep(c("Amblyomma", "Dermacentor",
                       "Haemaphysalis", "Hyalomma"), each = length(temp_seq)))) %>%
  left_join(species_to_genus.raw, by = "Genus")

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.genus.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.genus.raw, aes(x = Temperature, y = Fitted, color = Genus), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Genus",
       x = "Temperature", y = "Oviposition Time in Days") +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~ Species) 

# Merge the observed and fitted data for comparison
merged_data.genus.raw <- fitted_data.genus.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species", "Genus" = "Genus"))

# Calculate RMSE for each species
mse_by_genus.raw <- merged_data.genus.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.genus.raw <- merged_data.genus.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 



save(mse_by_genus.raw,
     fitted_data.genus.raw,
     merged_data.genus.raw,
     genus.raw.fit.mcmc,
     genus.raw.fit,
     log.lik.genus.raw,
     file="genus.raw.fit.RData")

################################### By Climate ###################################
# Model using climate effect
sink("climate_robin.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Climate-specific 
  mu.c ~ dexp(10)          # Mean for climate effect c (vertex y)
  mu.Topt ~ dnorm(30, 0.1)          # Mean for climate effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for climate effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for climate effects a (width)

  for (j in 1:n.climate.raw) {
    mu.climate.c[j] ~ dexp(mu.c)       # Climate-specific vertex y
    mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific vertex x
    mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific width
    climate.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.raw) {
    mu[i] <- mu.climate.a[climate.raw[i]] * (mu.climate.Topt[climate.raw[i]] - temp[i])^2 + mu.climate.c[climate.raw[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[climate.raw[i]]) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], climate.tau[climate.raw[i]])
  }}
", file = "climate_robin.txt")


# Parameters to monitor
parameters <- c("mu.climate.a", "mu.climate.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt", "climate.tau",
                "tau1", "mu.climate.Topt", "log_lik")


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
jag.data <- list(trait = tick.raw$Trait, N.obs.raw = N.obs.raw, 
                 temp = tick.raw$Temp, 
                 n.climate.raw = n.climate.raw, climate.raw = climate.raw)
# Fit model
climate.raw.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                    model.file = "climate_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                    n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
climate.raw.fit.mcmc <- as.mcmc(climate.raw.fit) 
closeAllConnections()
climate.raw.fit$BUGSoutput$DIC


# Extract log-likelihood matrix (samples × observations)
log.lik.climate.raw <- climate.raw.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.climate.raw <- waic(log.lik.climate.raw)$estimates["waic", "Estimate"]

# Get posterior means for parameters for each climate group
summary_fit.climate.raw <- summary(climate.raw.fit.mcmc)$statistics
a_climate_mean.raw <- summary_fit.climate.raw[grep("^mu.climate.a\\[", rownames(summary_fit.climate.raw)), "Mean"]
c_climate_mean.raw <- summary_fit.climate.raw[grep("^mu.climate.c\\[", rownames(summary_fit.climate.raw)), "Mean"]
Topt_climate_mean.raw <- summary_fit.climate.raw[grep("^mu.climate.Topt\\[", rownames(summary_fit.climate.raw)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp) ,
                max(tick.raw$Temp) , by = 0.1)

# Get predicted values for each climate group
climate1.raw <- function (x){
  a_climate_mean.raw[1] * (x - Topt_climate_mean.raw[1])^2 + c_climate_mean.raw[1]
}
climate2.raw <- function (x){
  a_climate_mean.raw[2] * (x - Topt_climate_mean.raw[2])^2 + c_climate_mean.raw[2]
}
climate3.raw <- function (x){
  a_climate_mean.raw[3] * (x - Topt_climate_mean.raw[3])^2 + c_climate_mean.raw[3]
}
climate4.raw <- function (x){
  a_climate_mean.raw[4] * (x - Topt_climate_mean.raw[4])^2 + c_climate_mean.raw[4]
}

# Data frame of temperature, fitted values, and climate
fitted_data.climate.raw <- data.frame(Temperature = rep(temp_seq, n.climate.raw), 
                                  Fitted = c(climate1.raw(temp_seq),
                                             climate2.raw(temp_seq),
                                             climate3.raw(temp_seq),
                                             climate4.raw(temp_seq)),
                                  Climate = factor(rep(c("Mixed", "Subtropical", 
                                                         "Temperate", "Tropical"), each = length(temp_seq))))
# Original data with temperature, trait, numeric climate, and numeric species
plot_data.climate.raw <- data.frame(
  Temperature = tick.raw$Temp,
  Trait = tick.raw$Trait,
  Climate = as.factor(tick.raw$Climate),
  Species = as.factor(tick.raw$Species)
)

# Get unique combinations of species and climate
species_to_climate.raw <- unique(plot_data.climate.raw[, c("Species", "Climate")])  
species_to_climate.raw$Climate <- as.factor(species_to_climate.raw$Climate)
# Merge fitted climate curves with species-to-climate mapping
fitted_data.climate.raw <- data.frame(
  Temperature = rep(temp_seq, n.climate.raw),
  Fitted = c(climate1.raw(temp_seq),
             climate2.raw(temp_seq),
             climate3.raw(temp_seq),
             climate4.raw(temp_seq)),
  Climate = factor(rep(c("Mixed", "Subtropical", 
                         "Temperate", "Tropical"), each = length(temp_seq)))) %>%
  left_join(species_to_climate.raw, by = "Climate")

plot_data.climate.raw$Climate <- factor(plot_data.climate.raw$Climate, levels = c("Mixed", "Temperate", "Subtropical", "Tropical"))
fitted_data.climate.raw$Climate <- factor(fitted_data.climate.raw$Climate, levels = c("Mixed", "Temperate", "Subtropical", "Tropical"))
# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.climate.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.climate.raw, aes(x = Temperature, y = Fitted, color = Climate), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Climate",
       x = "Temperature", y = "Oviposition Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) 

# Merge the observed and predicted data for comparison
merged_data.climate.raw <- fitted_data.climate.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species", "Climate" = "Climate"))

# Calculate RMSE for each species
mse_by_climate.raw <- merged_data.climate.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.climate.raw <- merged_data.climate.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 


save(mse_by_climate.raw,
     fitted_data.climate.raw,
     merged_data.climate.raw,
     climate.raw.fit.mcmc,
     climate.raw.fit,
     log.lik.climate.raw,
     file="climate.raw.fit.RData")

################################### By Host ###################################
# Model by typical number of host species
sink("host_robin.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Host-specific 
  mu.c ~ dexp(10)          # Mean for host effects c (vertex y)
  mu.Topt ~ dnorm(30, 0.1)          # Mean for host effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for host effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for host effects a (width)
  
  for (j in 1:n.host.raw) {
    mu.host.c[j] ~ dexp(mu.c)       # Host-specific vertex y
    mu.host.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Host-specific vertex x
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
    host.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.raw) {
    mu[i] <- mu.host.a[host.raw[i]] * (mu.host.Topt[host.raw[i]] - temp[i])^2 + mu.host.c[host.raw[i]]
    trait[i] ~ dnorm(mu[i], host.tau[host.raw[i]]) T(0, ) 
    log_lik[i] <- logdensity.norm(trait[i], mu[i], host.tau[host.raw[i]])
  }}
", file = "host_robin.txt")


# Parameters to monitor
parameters <- c("mu.host.a", "mu.host.c", "mu.c", "mu.a",
                "mu.Topt", "tau.Topt", "host.tau",
                "tau1", "mu.host.Topt", "log_lik") 


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
jag.data <- list(trait = tick.raw$Trait, N.obs.raw = N.obs.raw, 
                 temp = tick.raw$Temp, 
                 n.host.raw = n.host.raw, host.raw = host.raw)
# Fit model
host.raw.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                 model.file = "host_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                 n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.raw.fit.mcmc <- as.mcmc(host.raw.fit) 
closeAllConnections()
host.raw.fit$BUGSoutput$DIC


# Extract log-likelihood matrix (samples × observations)
log.lik.host.raw <- host.raw.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.host.raw <- waic(log.lik.host.raw)$estimates["waic", "Estimate"]

# Get posterior means for parameters for each host group
summary_fit.host.raw <- summary(host.raw.fit.mcmc)$statistics
a_host_mean.raw <- summary_fit.host.raw[grep("^mu.host.a\\[", rownames(summary_fit.host.raw)), "Mean"]
c_host_mean.raw <- summary_fit.host.raw[grep("^mu.host.c\\[", rownames(summary_fit.host.raw)), "Mean"]
Topt_host_mean.raw <- summary_fit.host.raw[grep("^mu.host.Topt\\[", rownames(summary_fit.host.raw)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp) ,
                max(tick.raw$Temp) , by = 0.1)

# Get predicted values for each host group
host1.raw <- function (x){
  a_host_mean.raw[1] * (x - Topt_host_mean.raw[1])^2 + c_host_mean.raw[1]
}
host2.raw <- function (x){
  a_host_mean.raw[2] * (x - Topt_host_mean.raw[2])^2 + c_host_mean.raw[2]
}
host3.raw <- function (x){
  a_host_mean.raw[3] * (x - Topt_host_mean.raw[3])^2 + c_host_mean.raw[3]
}


# Data frame with temperature, fitted value, and host group
fitted_data.host.raw <- data.frame(Temperature = rep(temp_seq, n.host.raw), 
                               Fitted = c(host1.raw(temp_seq),
                                          host2.raw(temp_seq),
                                          host3.raw(temp_seq)),
                               Host = factor(rep(c("One", "Three", "Two"), each = length(temp_seq))))
# Original data with temperature, trait, numeric host group, numeric species
plot_data.host.raw <- data.frame(
  Temperature = tick.raw$Temp,
  Trait = tick.raw$Trait,
  Host = as.factor(tick.raw$Host),
  Species = as.factor(tick.raw$Species)
)

# Get unique combinations of species and host for correct plotting
species_to_host.raw <- unique(plot_data.host.raw[, c("Species", "Host")])  
species_to_host.raw$Host <- as.factor(species_to_host.raw$Host)
# Merge fitted host curves with species-to-host mapping
fitted_data.host.raw <- data.frame(
  Temperature = rep(temp_seq, n.host.raw),
  Fitted = c(host1.raw(temp_seq),
             host2.raw(temp_seq),
             host3.raw(temp_seq)),
  Host = factor(rep(c("One", "Three", "Two"), each = length(temp_seq)))
) %>%
  left_join(species_to_host.raw, by = "Host")

plot_data.host.raw$Host <- factor(plot_data.host.raw$Host, levels = c("One", "Two", "Three"))
fitted_data.host.raw$Host <- factor(fitted_data.host.raw$Host, levels = c("One", "Two", "Three"))
# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.host.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.raw, aes(x = Temperature, y = Fitted, color = Host), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Host",
       x = "Temperature", y = "Oviposition Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) 
#ggsave("tick.host.plot.png", width = 8, height = 5)


# Merge the observed and predicted data for comparison
merged_data.host.raw <- fitted_data.host.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host"))

# Calculate RMSE for each species
mse_by_host.raw <- merged_data.host.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.host.raw <- merged_data.host.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 


save(mse_by_host.raw,
     fitted_data.host.raw,
     merged_data.host.raw,
     host.raw.fit.mcmc,
     host.raw.fit,
     log.lik.host.raw,
     file="host.raw.fit.RData")


################################### By Size ###################################
# Model by size group
sink("size_robin.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Host-specific 
  mu.c ~ dexp(10)          # Mean for size effects c (vertex y)
  mu.Topt ~ dnorm(30, 0.1)          # Mean for size effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for size effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for size effects a (width)

  for (j in 1:n.size.raw) {
    mu.size.c[j] ~ dexp(mu.c)       # Size-specific vertex y
    mu.size.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Size-specific vertex x
    mu.size.a[j] ~ dexp(mu.a)    # Size-specific width
    size.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.raw) {
    mu[i] <- mu.size.a[size.raw[i]] * (mu.size.Topt[size.raw[i]] - temp[i])^2 + mu.size.c[size.raw[i]]
    trait[i] ~ dnorm(mu[i], size.tau[size.raw[i]]) T(0, )  
    log_lik[i] <- logdensity.norm(trait[i], mu[i], size.tau[size.raw[i]])
  }}
", file = "size_robin.txt")


# Parameters to monitor
parameters <- c("mu.size.a", "mu.size.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt", "size.tau",
                "tau1", "mu.size.Topt", "log_lik")


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
jag.data <- list(trait = tick.raw$Trait, N.obs.raw = N.obs.raw, 
                 temp = tick.raw$Temp, 
                 n.size.raw = n.size.raw, size.raw = size.raw)
# Fit model
size.raw.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                 model.file = "size_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                 n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
size.raw.fit.mcmc <- as.mcmc(size.raw.fit)
closeAllConnections()
size.raw.fit$BUGSoutput$DIC


# Extract log-likelihood matrix (samples × observations)
log.lik.size.raw <- size.raw.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.size.raw <- waic(log.lik.size.raw)$estimates["waic", "Estimate"]

# Get posterior means for parameters for each size 
summary_fit.size.raw <- summary(size.raw.fit.mcmc)$statistics
a_size_mean.raw <- summary_fit.size.raw[grep("^mu.size.a\\[", rownames(summary_fit.size.raw)), "Mean"]
c_size_mean.raw <- summary_fit.size.raw[grep("^mu.size.c\\[", rownames(summary_fit.size.raw)), "Mean"]
Topt_size_mean.raw <- summary_fit.size.raw[grep("^mu.size.Topt\\[", rownames(summary_fit.size.raw)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp) ,
                max(tick.raw$Temp) , by = 0.1)

# Get predicted values for each size category
size1.raw <- function (x){
  a_size_mean.raw[1] * (x - Topt_size_mean.raw[1])^2 + c_size_mean.raw[1]
}
size2.raw <- function (x){
  a_size_mean.raw[2] * (x - Topt_size_mean.raw[2])^2 + c_size_mean.raw[2]
}
size3.raw <- function (x){
  a_size_mean.raw[3] * (x - Topt_size_mean.raw[3])^2 + c_size_mean.raw[3]
}


# Make data frame with temperature, fitted value, and size category
fitted_data.size.raw <- data.frame(Temperature = rep(temp_seq, n.size.raw), 
                               Fitted = c(size1.raw(temp_seq),
                                          size1.raw(temp_seq),
                                          size3.raw(temp_seq)),
                               Size = factor(rep(c("Large", "Medium", "Small"), each = length(temp_seq))))
# Original data with temperature, trait, numeric size, and numeric species
plot_data.size.raw <- data.frame(
  Temperature = tick.raw$Temp,
  Trait = tick.raw$Trait,
  Size = as.factor(tick.raw$Size),
  Species = as.factor(tick.raw$Species)
)

# Get unique combinations of species and size for correct fitted lines on plots
species_to_size.raw <- unique(plot_data.size.raw[, c("Species", "Size")])  
species_to_size.raw$Size <- as.factor(species_to_size.raw$Size)
# Merge fitted host curves with species-to-size mapping
fitted_data.size.raw <- data.frame(
  Temperature = rep(temp_seq, n.size.raw),
  Fitted = c(size1.raw(temp_seq),
             size2.raw(temp_seq),
             size3.raw(temp_seq)),
  Size = factor(rep(c("Large", "Medium", "Small"), each = length(temp_seq)))
) %>%
  left_join(species_to_size.raw, by = "Size")

plot_data.size.raw$Size <- factor(plot_data.size.raw$Size, levels = c("Small", "Medium", "Large"))
fitted_data.size.raw$Size <- factor(fitted_data.size.raw$Size, levels = c("Small", "Medium", "Large"))
# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.size.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.size.raw, aes(x = Temperature, y = Fitted, color = Size), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Size",
       x = "Temperature", y = "Oviposition Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) 


# Merge the observed and predicted data for comparison
merged_data.size.raw <- fitted_data.size.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species", "Size" = "Size"))

# Calculate RMSE for each species
mse_by_size.raw <- merged_data.size.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.size.raw <- merged_data.size.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_size.raw,
     fitted_data.size.raw,
     merged_data.size.raw,
     size.raw.fit.mcmc,
     size.raw.fit,
     log.lik.size.raw,
     file="size.raw.fit.RData")

################################### By Continent ###################################
# Model by number of continents tick is found on
sink("continent_robin.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Continent-specific 
  mu.c ~ dexp(10)          # Mean for continent effects c (vertex y)
  mu.Topt ~ dnorm(30, 0.1)          # Mean for continent effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for continent effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for continent effects a (width)
 

  for (j in 1:n.continent.raw) {
    mu.continent.c[j] ~ dexp(mu.c)       # Continent-specific vertex y
    mu.continent.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Continent-specific vertex x
    mu.continent.a[j] ~ dexp(mu.a)    # Continent-specific width
    continent.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.raw) {
    mu[i] <- mu.continent.a[continent.raw[i]] * (mu.continent.Topt[continent.raw[i]] - temp[i])^2 + mu.continent.c[continent.raw[i]]
    trait[i] ~ dnorm(mu[i], continent.tau[continent.raw[i]]) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], continent.tau[continent.raw[i]])
  }}
", file = "continent_robin.txt")


# Parameters to monitor
parameters <- c("mu.continent.a", "mu.continent.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt", "continent.tau",
                "tau1", "mu.continent.Topt", "log_lik")


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
jag.data <- list(trait = tick.raw$Trait, N.obs.raw = N.obs.raw, 
                 temp = tick.raw$Temp, 
                 n.continent.raw = n.continent.raw, continent.raw = continent.raw)
# Fit model
continent.raw.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                      model.file = "continent_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                      n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
continent.raw.fit.mcmc <- as.mcmc(continent.raw.fit)
closeAllConnections()
continent.raw.fit$BUGSoutput$DIC


# Extract log-likelihood matrix (samples × observations)
log.lik.continent.raw <- continent.raw.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.continent.raw <- waic(log.lik.continent.raw)$estimates["waic", "Estimate"]

# Get posterior means for parameters for different continent groups
summary_fit.continent.raw <- summary(continent.raw.fit.mcmc)$statistics
a_continent_mean.raw <- summary_fit.continent.raw[grep("^mu.continent.a\\[", rownames(summary_fit.continent.raw)), "Mean"]
c_continent_mean.raw <- summary_fit.continent.raw[grep("^mu.continent.c\\[", rownames(summary_fit.continent.raw)), "Mean"]
Topt_continent_mean.raw <- summary_fit.continent.raw[grep("^mu.continent.Topt\\[", rownames(summary_fit.continent.raw)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp) ,
                max(tick.raw$Temp) , by = 0.1)

# Get predicted values for each continent group
continent1.raw <- function (x){
  a_continent_mean.raw[1] * (x - Topt_continent_mean.raw[1])^2 + c_continent_mean.raw[1]
}
continent2.raw <- function (x){
  a_continent_mean.raw[2] * (x - Topt_continent_mean.raw[2])^2 + c_continent_mean.raw[2]
}
continent3.raw <- function (x){
  a_continent_mean.raw[3] * (x - Topt_continent_mean.raw[3])^2 + c_continent_mean.raw[3]
}

# Data frame of temperature, fitted value, and continent group
fitted_data.continent.raw <- data.frame(Temperature = rep(temp_seq, n.continent.raw), 
                                    Fitted = c(continent1.raw(temp_seq),
                                               continent2.raw(temp_seq),
                                               continent3.raw(temp_seq)),
                                    Continent = factor(rep(c("More than two", "One", "Two"), each = length(temp_seq))))
# Original data with temperature, trait, numeric continent, numeric species
plot_data.continent.raw <- data.frame(
  Temperature = tick.raw$Temp,
  Trait = tick.raw$Trait,
  Continent = as.factor(tick.raw$Continent),
  Species = as.factor(tick.raw$Species)
)
# Get unique combinations of species and continent for correct plotting
species_to_continent.raw <- unique(plot_data.continent.raw[, c("Species", "Continent")]) 
species_to_continent.raw$Continent <- as.factor(species_to_continent.raw$Continent)
# Merge fitted continent curves with species-to-continent mapping
fitted_data.continent.raw <- data.frame(
  Temperature = rep(temp_seq, n.continent.raw),
  Fitted = c(continent1.raw(temp_seq),
             continent2.raw(temp_seq),
             continent3.raw(temp_seq)),
  Continent = factor(rep(c("More than two", "One", "Two"), each = length(temp_seq)))
) %>%
  left_join(species_to_continent.raw, by = "Continent")

plot_data.continent.raw$Continent <- factor(plot_data.continent.raw$Continent, levels = c("One", "Two", "More than two"))
fitted_data.continent.raw$Continent <- factor(fitted_data.continent.raw$Continent, levels = c("One", "Two", "More than two"))
# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.continent.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.continent.raw, aes(x = Temperature, y = Fitted, color = Continent), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Continent",
       x = "Temperature", y = "Oviposition Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) 

# Merge the observed and predicted data for comparison
merged_data.continent.raw <- fitted_data.continent.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species", "Continent" = "Continent"))

# Calculate RMSE for each species
mse_by_continent.raw <- merged_data.continent.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.continent.raw <- merged_data.continent.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_continent.raw,
     fitted_data.continent.raw,
     merged_data.continent.raw,
     continent.raw.fit.mcmc,
     continent.raw.fit,
     log.lik.continent.raw,
     file="continent.raw.fit.RData")

################################### By Host and Climate ###################################
# Model by typical number of host species
sink("host.climate_robin.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Host-specific 
  mu.c ~ dexp(10)          # Mean for host effects c (vertex y)
  mu.a ~ dexp(0.01)            # Mean for host effects a (width)
  # Climate-specific
  mu.Topt ~ dnorm(30, 0.1)          # Mean for climate effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for climate effect Topt (vertex x)
  
  for (j in 1:n.host.raw) {
    mu.host.c[j] ~ dexp(mu.c)       # Host-specific vertex y
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
  }
  for (k in 1:n.climate.raw) {
    mu.climate.Topt[k] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific vertex x
    climate.tau[k] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.raw) {
    mu[i] <- mu.host.a[host.raw[i]] * (mu.climate.Topt[climate.raw[i]] - temp[i])^2 + mu.host.c[host.raw[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[climate.raw[i]]) T(0, ) 
    log_lik[i] <- logdensity.norm(trait[i], mu[i], climate.tau[climate.raw[i]])
  }}
", file = "host.climate_robin.txt")


# Parameters to monitor
parameters <- c("mu.host.a", "mu.host.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt", "climate.tau",
                "tau1", "mu.climate.Topt", "log_lik")


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
jag.data <- list(trait = tick.raw$Trait, N.obs.raw = N.obs.raw, 
                 temp = tick.raw$Temp, 
                 n.host.raw = n.host.raw, host.raw = host.raw, 
                 n.climate.raw = n.climate.raw, climate.raw = climate.raw)
# Fit model
host.climate.raw.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                         model.file = "host.climate_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                         n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.climate.raw.fit.mcmc <- as.mcmc(host.climate.raw.fit) 
closeAllConnections()
host.climate.raw.fit$BUGSoutput$DIC


# Extract log-likelihood matrix (samples × observations)
log.lik.host.climate.raw <- host.climate.raw.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.host.climate.raw <- waic(log.lik.host.climate.raw)$estimates["waic", "Estimate"]

# Get posterior means for parameters for each host group
summary_fit.host.climate.raw <- summary(host.climate.raw.fit.mcmc)$statistics
a_host.climate_mean.raw <- summary_fit.host.climate.raw[grep("^mu.host.a\\[", rownames(summary_fit.host.climate.raw)), "Mean"]
c_host.climate_mean.raw <- summary_fit.host.climate.raw[grep("^mu.host.c\\[", rownames(summary_fit.host.climate.raw)), "Mean"]
Topt_host.climate_mean.raw <- summary_fit.host.climate.raw[grep("^mu.climate.Topt\\[", rownames(summary_fit.host.climate.raw)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp) ,
                max(tick.raw$Temp) , by = 0.1)

# Get predicted values for each host group
host1clim1.raw <- function (x){
  a_host.climate_mean.raw[1] * (x - Topt_host.climate_mean.raw[1])^2 + c_host.climate_mean.raw[1]
}
host2clim1.raw <- function (x){
  a_host.climate_mean.raw[2] * (x - Topt_host.climate_mean.raw[1])^2 + c_host.climate_mean.raw[2]
}
host3clim1.raw <- function (x){
  a_host.climate_mean.raw[3] * (x - Topt_host.climate_mean.raw[1])^2 + c_host.climate_mean.raw[3]
}
host1clim2.raw <- function (x){
  a_host.climate_mean.raw[1] * (x - Topt_host.climate_mean.raw[2])^2 + c_host.climate_mean.raw[1]
}
host2clim2.raw <- function (x){
  a_host.climate_mean.raw[2] * (x - Topt_host.climate_mean.raw[2])^2 + c_host.climate_mean.raw[2]
}
host3clim2.raw <- function (x){
  a_host.climate_mean.raw[3] * (x - Topt_host.climate_mean.raw[2])^2 + c_host.climate_mean.raw[3]
}

host1clim3.raw <- function (x){
  a_host.climate_mean.raw[1] * (x - Topt_host.climate_mean.raw[3])^2 + c_host.climate_mean.raw[1]
}
host2clim3.raw <- function (x){
  a_host.climate_mean.raw[2] * (x - Topt_host.climate_mean.raw[3])^2 + c_host.climate_mean.raw[2]
}
host3clim3.raw <- function (x){
  a_host.climate_mean.raw[3] * (x - Topt_host.climate_mean.raw[3])^2 + c_host.climate_mean.raw[3]
}

host1clim4.raw <- function (x){
  a_host.climate_mean.raw[1] * (x - Topt_host.climate_mean.raw[4])^2 + c_host.climate_mean.raw[1]
}
host2clim4.raw <- function (x){
  a_host.climate_mean.raw[2] * (x - Topt_host.climate_mean.raw[4])^2 + c_host.climate_mean.raw[2]
}
host3clim4.raw <- function (x){
  a_host.climate_mean.raw[3] * (x - Topt_host.climate_mean.raw[4])^2 + c_host.climate_mean.raw[3]
}



# Data frame with temperature, fitted value, and host group
fitted_data.host.climate.raw <- data.frame(Temperature = rep(temp_seq, n.host.raw * n.climate.raw), 
                                       Fitted = c(host1clim1.raw(temp_seq),
                                                  host2clim1.raw(temp_seq),
                                                  host3clim1.raw(temp_seq),
                                                  host1clim2.raw(temp_seq),
                                                  host2clim2.raw(temp_seq),
                                                  host3clim2.raw(temp_seq),
                                                  host1clim3.raw(temp_seq),
                                                  host2clim3.raw(temp_seq),
                                                  host3clim3.raw(temp_seq),
                                                  host1clim4.raw(temp_seq),
                                                  host2clim4.raw(temp_seq),
                                                  host3clim4.raw(temp_seq)),
                                       Host = factor(rep(rep(c("One", "Three", "Two"), each = length(temp_seq)), n.climate.raw)),
                                       Climate = factor(rep(c("Mixed", "Subtropical", "Temperate", "Tropical"),
                                                            each = length(temp_seq) * n.size.raw)))
# Original data with temperature, trait, numeric host group, numeric species
plot_data.host.climate.raw <- data.frame(
  Temperature = tick.raw$Temp,
  Trait = tick.raw$Trait,
  Host = as.factor(tick.raw$Host),
  Species = as.factor(tick.raw$Species),
  Climate = as.factor(tick.raw$Climate)
)

# Get unique combinations of species and host for correct plotting
species_to_host.climate.raw <- unique(plot_data.host.climate.raw[, c("Species", "Host", "Climate")])  
species_to_host.climate.raw$Host <- as.factor(species_to_host.climate.raw$Host)
species_to_host.climate.raw$Climate <- as.factor(species_to_host.climate.raw$Climate)

# Merge fitted host curves with species-to-host mapping
fitted_data.host.climate.raw <- data.frame(Temperature = rep(temp_seq, n.host.raw * n.climate.raw), 
                                       Fitted = c(host1clim1.raw(temp_seq),
                                                  host2clim1.raw(temp_seq),
                                                  host3clim1.raw(temp_seq),
                                                  host1clim2.raw(temp_seq),
                                                  host2clim2.raw(temp_seq),
                                                  host3clim2.raw(temp_seq),
                                                  host1clim3.raw(temp_seq),
                                                  host2clim3.raw(temp_seq),
                                                  host3clim3.raw(temp_seq),
                                                  host1clim4.raw(temp_seq),
                                                  host2clim4.raw(temp_seq),
                                                  host3clim4.raw(temp_seq)),
                                       Host = factor(rep(rep(c("One", "Three", "Two"), each = length(temp_seq)), n.climate.raw)),
                                       Climate = factor(rep(c("Mixed", "Subtropical", "Temperate", "Tropical"),
                                                            each = length(temp_seq) * n.size.raw))) %>%
  left_join(species_to_host.climate.raw, by = c("Host", "Climate"))|> 
  filter(!is.na(Species)) 

plot_data.host.climate.raw$Host <- factor(plot_data.host.climate.raw$Host, levels = c("One", "Two", "Three"))
fitted_data.host.climate.raw$Host <- factor(fitted_data.host.climate.raw$Host, levels = c("One", "Two", "Three"))

# Combine Climate and Size for coloring
fitted_data.host.climate.raw$Combo <- with(fitted_data.host.climate.raw, paste(Host, Climate, sep = ", "))
fitted_data.host.climate.raw$Combo <- factor(fitted_data.host.climate.raw$Combo, levels = c("One, Tropical", 
                                                                                    "Two, Subtropical", "Three, Mixed",
                                                                                    "Three, Subtropical", "Three, Temperate",
                                                                                    "Three, Tropical"))

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.host.climate.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.climate.raw, aes(x = Temperature, y = Fitted, color = Combo), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Host and Climate",
       x = "Temperature", y = "Oviposition Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) 

# Merge the observed and predicted data for comparison
merged_data.host.climate.raw <- fitted_data.host.climate.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host", "Climate" = "Climate"))

# Calculate RMSE for each species
mse_by_host.climate.raw <- merged_data.host.climate.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.host.climate.raw <- merged_data.host.climate.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 


save(mse_by_host.climate.raw,
     fitted_data.host.climate.raw,
     merged_data.host.climate.raw,
     host.climate.raw.fit.mcmc,
     host.climate.raw.fit,
     log.lik.host.climate.raw,
     file="host.climate.raw.fit.RData")

################################### By Host and Genus ###################################
# Model by typical number of host species
sink("host.genus_robin.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Host-specific 
   mu.a ~ dexp(0.01)            # Mean for host effects a (width)
   # Genus-specific
  mu.c ~ dexp(10)          # Mean for genus effects c (vertex y)
  mu.Topt ~ dnorm(30, 0.1)          # Mean for genus effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for genus effect Topt (vertex x)
 

  for (j in 1:n.host.raw) {
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
    host.tau[j] ~ dexp(tau1)
  }
  for (k in 1:n.genus.raw) {
    mu.genus.c[k] ~ dexp(mu.c)       # Genus-specific vertex y
    mu.genus.Topt[k] ~ dnorm(mu.Topt, tau.Topt)    # Genus-specific vertex x
  }
  # Likelihood
  for (i in 1:N.obs.raw) {
    mu[i] <- mu.host.a[host.raw[i]] * (mu.genus.Topt[genus.raw[i]] - temp[i])^2 + mu.genus.c[genus.raw[i]]
    trait[i] ~ dnorm(mu[i], host.tau[host.raw[i]]) T(0, )  
    log_lik[i] <- logdensity.norm(trait[i], mu[i], host.tau[host.raw[i]])
  }}
", file = "host.genus_robin.txt")


# Parameters to monitor
parameters <- c("mu.host.a", "mu.genus.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt", "host.tau",
                "tau1", "mu.genus.Topt", "log_lik")


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
jag.data <- list(trait = tick.raw$Trait, N.obs.raw = N.obs.raw, 
                 temp = tick.raw$Temp, 
                 n.host.raw = n.host.raw, host.raw = host.raw, 
                 n.genus.raw = n.genus.raw, genus.raw = genus.raw)
# Fit model
host.genus.raw.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                       model.file = "host.genus_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                       n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.genus.raw.fit.mcmc <- as.mcmc(host.genus.raw.fit) 
closeAllConnections()
host.genus.raw.fit$BUGSoutput$DIC


# Extract log-likelihood matrix (samples × observations)
log.lik.host.genus.raw <- host.genus.raw.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.host.genus.raw <- waic(log.lik.host.genus.raw)$estimates["waic", "Estimate"]

# Get posterior means for parameters for each host group
summary_fit.host.genus.raw <- summary(host.genus.raw.fit.mcmc)$statistics
a_host.genus_mean.raw <- summary_fit.host.genus.raw[grep("^mu.host.a\\[", rownames(summary_fit.host.genus.raw)), "Mean"]
c_host.genus_mean.raw <- summary_fit.host.genus.raw[grep("^mu.genus.c\\[", rownames(summary_fit.host.genus.raw)), "Mean"]
Topt_host.genus_mean.raw <- summary_fit.host.genus.raw[grep("^mu.genus.Topt\\[", rownames(summary_fit.host.genus.raw)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp) ,
                max(tick.raw$Temp) , by = 0.1)

# Get predicted values for each host group
host1genus1.raw <- function (x){
  a_host.genus_mean.raw[1] * (x - Topt_host.genus_mean.raw[1])^2 + c_host.genus_mean.raw[1]
}
host2genus1.raw <- function (x){
  a_host.genus_mean.raw[2] * (x - Topt_host.genus_mean.raw[1])^2 + c_host.genus_mean.raw[1]
}
host3genus1.raw <- function (x){
  a_host.genus_mean.raw[3] * (x - Topt_host.genus_mean.raw[1])^2 + c_host.genus_mean.raw[1]
}
host1genus2.raw <- function (x){
  a_host.genus_mean.raw[1] * (x - Topt_host.genus_mean.raw[2])^2 + c_host.genus_mean.raw[2]
}
host2genus2.raw <- function (x){
  a_host.genus_mean.raw[2] * (x - Topt_host.genus_mean.raw[2])^2 + c_host.genus_mean.raw[2]
}
host3genus2.raw <- function (x){
  a_host.genus_mean.raw[3] * (x - Topt_host.genus_mean.raw[2])^2 + c_host.genus_mean.raw[2]
}

host1genus3.raw <- function (x){
  a_host.genus_mean.raw[1] * (x - Topt_host.genus_mean.raw[3])^2 + c_host.genus_mean.raw[3]
}
host2genus3.raw <- function (x){
  a_host.genus_mean.raw[2] * (x - Topt_host.genus_mean.raw[3])^2 + c_host.genus_mean.raw[3]
}
host3genus3.raw <- function (x){
  a_host.genus_mean.raw[3] * (x - Topt_host.genus_mean.raw[3])^2 + c_host.genus_mean.raw[3]
}

host1genus4.raw <- function (x){
  a_host.genus_mean.raw[1] * (x - Topt_host.genus_mean.raw[4])^2 + c_host.genus_mean.raw[4]
}
host2genus4.raw <- function (x){
  a_host.genus_mean.raw[2] * (x - Topt_host.genus_mean.raw[4])^2 + c_host.genus_mean.raw[4]
}
host3genus4.raw <- function (x){
  a_host.genus_mean.raw[3] * (x - Topt_host.genus_mean.raw[4])^2 + c_host.genus_mean.raw[4]
}



# Data frame with temperature, fitted value, and host group
fitted_data.host.genus.raw <- data.frame(Temperature = rep(temp_seq, n.host.raw * n.genus.raw), 
                                     Fitted = c(host1genus1.raw(temp_seq),
                                                host2genus1.raw(temp_seq),
                                                host3genus1.raw(temp_seq),
                                                host1genus2.raw(temp_seq),
                                                host2genus2.raw(temp_seq),
                                                host3genus2.raw(temp_seq),
                                                host1genus3.raw(temp_seq),
                                                host2genus3.raw(temp_seq),
                                                host3genus3.raw(temp_seq),
                                                host1genus4.raw(temp_seq),
                                                host2genus4.raw(temp_seq),
                                                host3genus4.raw(temp_seq)),
                                     Host = factor(rep(rep(c("One", "Three", "Two"), each = length(temp_seq)), n.genus.raw)),
                                     Genus = factor(rep(c("Amblyomma", "Dermacentor", 
                                                          "Haemaphysalis", "Hyalomma"),
                                                        each = length(temp_seq) * n.host.raw)))
# Original data with temperature, trait, numeric host group, numeric species
plot_data.host.genus.raw <- data.frame(
  Temperature = tick.raw$Temp,
  Trait = tick.raw$Trait,
  Host = as.factor(tick.raw$Host),
  Species = as.factor(tick.raw$Species),
  Genus = as.factor(tick.raw$Genus)
)

# Get unique combinations of species and host for correct plotting
species_to_host.genus.raw <- unique(plot_data.host.genus.raw[, c("Species", "Host", "Genus")])  
species_to_host.genus.raw$Host <- as.factor(species_to_host.genus.raw$Host)
species_to_host.genus.raw$Genus <- as.factor(species_to_host.genus.raw$Genus)

# Merge fitted host curves with species-to-host mapping
fitted_data.host.genus.raw <- data.frame(Temperature = rep(temp_seq, n.host.raw * n.genus.raw), 
                                     Fitted = c(host1genus1.raw(temp_seq),
                                                host2genus1.raw(temp_seq),
                                                host3genus1.raw(temp_seq),
                                                host1genus2.raw(temp_seq),
                                                host2genus2.raw(temp_seq),
                                                host3genus2.raw(temp_seq),
                                                host1genus3.raw(temp_seq),
                                                host2genus3.raw(temp_seq),
                                                host3genus3.raw(temp_seq),
                                                host1genus4.raw(temp_seq),
                                                host2genus4.raw(temp_seq),
                                                host3genus4.raw(temp_seq)),
                                     Host = factor(rep(rep(c("One", "Three", "Two"), each = length(temp_seq)), n.genus.raw)),
                                     Genus = factor(rep(c("Amblyomma", "Dermacentor", 
                                                          "Haemaphysalis", "Hyalomma"),
                                                        each = length(temp_seq) * n.host.raw))) %>%
  left_join(species_to_host.genus.raw, by = c("Host", "Genus"))|> 
  filter(!is.na(Species)) 

plot_data.host.genus.raw$Host <- factor(plot_data.host.genus.raw$Host, levels = c("One", "Two", "Three"))
fitted_data.host.genus.raw$Host <- factor(fitted_data.host.genus.raw$Host, levels = c("One", "Two", "Three"))

# Combine genus and Size for coloring
fitted_data.host.genus.raw$Combo <- with(fitted_data.host.genus.raw, paste(Host, Genus, sep = ", "))
fitted_data.host.genus.raw$Combo <- factor(fitted_data.host.genus.raw$Combo, levels = c("One, Dermacentor", 
                                                                                "Two, Hyalomma", "Three, Amblyomma",
                                                                                "Three, Dermacentor", "Three, Haemaphysalis",
                                                                                "Three, Hyalomma"))

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.host.genus.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.genus.raw, aes(x = Temperature, y = Fitted, color = Combo), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Host and Genus",
       x = "Temperature", y = "Oviposition Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) 

# Merge the observed and predicted data for comparison
merged_data.host.genus.raw <- fitted_data.host.genus.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host", "Genus" = "Genus"))

# Calculate RMSE for each species
mse_by_host.genus.raw <- merged_data.host.genus.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.host.genus.raw <- merged_data.host.genus.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE)))


save(mse_by_host.genus.raw,
     fitted_data.host.genus.raw,
     merged_data.host.genus.raw,
     host.genus.raw.fit.mcmc,
     host.genus.raw.fit,
     log.lik.host.genus.raw,
     file="host.genus.raw.fit.RData")

################################### By Genus and Climate ###################################
# Model by typical number of host species
sink("genus.climate_robin.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # genus-specific 
  mu.c ~ dexp(10)          # Mean for genus effects c (vertex y)
   mu.Topt ~ dnorm(30, 0.1)          # Mean for genus effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for genus effect Topt (vertex x)
  # Climate-specific
  mu.a ~ dexp(0.01)            # Mean for climate effects a (width)

  for (j in 1:n.genus.raw) {
    mu.genus.c[j] ~ dexp(mu.c)       # genus-specific vertex y
    mu.genus.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # genus-specific vertex x
  }
  for (k in 1:n.climate.raw) {
    mu.climate.a[k] ~ dexp(mu.a)    # Cimate-specific width
    climate.tau[k] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.raw) {
    mu[i] <- mu.climate.a[climate.raw[i]] * (mu.genus.Topt[genus.raw[i]] - temp[i])^2 + mu.genus.c[genus.raw[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[climate.raw[i]]) T(0, ) 
    log_lik[i] <- logdensity.norm(trait[i], mu[i], climate.tau[climate.raw[i]])
  }}
", file = "genus.climate_robin.txt")


# Parameters to monitor
parameters <- c("mu.climate.a", "mu.genus.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt", "climate.tau",
                "tau1", "mu.genus.Topt", "log_lik")


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
jag.data <- list(trait = tick.raw$Trait, N.obs.raw = N.obs.raw, 
                 temp = tick.raw$Temp, 
                 n.genus.raw = n.genus.raw, genus.raw = genus.raw, 
                 n.climate.raw = n.climate.raw, climate.raw = climate.raw)
# Fit model
genus.climate.raw.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                          model.file = "genus.climate_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                          n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
genus.climate.raw.fit.mcmc <- as.mcmc(genus.climate.raw.fit) 
closeAllConnections()
genus.climate.raw.fit$BUGSoutput$DIC


# Extract log-likelihood matrix (samples × observations)
log.lik.genus.climate.raw <- genus.climate.raw.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.genus.climate.raw <- waic(log.lik.genus.climate.raw)$estimates["waic", "Estimate"]


# Get posterior means for parameters for each genus group
summary_fit.genus.climate.raw <- summary(genus.climate.raw.fit.mcmc)$statistics
a_genus.climate_mean.raw <- summary_fit.genus.climate.raw[grep("^mu.climate.a\\[", rownames(summary_fit.genus.climate.raw)), "Mean"]
c_genus.climate_mean.raw <- summary_fit.genus.climate.raw[grep("^mu.genus.c\\[", rownames(summary_fit.genus.climate.raw)), "Mean"]
Topt_genus.climate_mean.raw <- summary_fit.genus.climate.raw[grep("^mu.genus.Topt\\[", rownames(summary_fit.genus.climate.raw)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp) ,
                max(tick.raw$Temp) , by = 0.1)

# Get predicted values for each genus group
genus1clim1.raw <- function (x){
  a_genus.climate_mean.raw[1] * (x - Topt_genus.climate_mean.raw[1])^2 + c_genus.climate_mean.raw[1]
}
genus2clim1.raw <- function (x){
  a_genus.climate_mean.raw[1] * (x - Topt_genus.climate_mean.raw[2])^2 + c_genus.climate_mean.raw[2]
}
genus3clim1.raw <- function (x){
  a_genus.climate_mean.raw[1] * (x - Topt_genus.climate_mean.raw[3])^2 + c_genus.climate_mean.raw[3]
}
genus4clim1.raw <- function (x){
  a_genus.climate_mean.raw[1] * (x - Topt_genus.climate_mean.raw[4])^2 + c_genus.climate_mean.raw[4]
}
genus1clim2.raw <- function (x){
  a_genus.climate_mean.raw[2] * (x - Topt_genus.climate_mean.raw[1])^2 + c_genus.climate_mean.raw[1]
}
genus2clim2.raw <- function (x){
  a_genus.climate_mean.raw[2] * (x - Topt_genus.climate_mean.raw[2])^2 + c_genus.climate_mean.raw[2]
}
genus3clim2.raw <- function (x){
  a_genus.climate_mean.raw[2] * (x - Topt_genus.climate_mean.raw[3])^2 + c_genus.climate_mean.raw[3]
}
genus4clim2.raw <- function (x){
  a_genus.climate_mean.raw[2] * (x - Topt_genus.climate_mean.raw[4])^2 + c_genus.climate_mean.raw[4]
}
genus1clim3.raw <- function (x){
  a_genus.climate_mean.raw[3] * (x - Topt_genus.climate_mean.raw[1])^2 + c_genus.climate_mean.raw[1]
}
genus2clim3.raw <- function (x){
  a_genus.climate_mean.raw[3] * (x - Topt_genus.climate_mean.raw[2])^2 + c_genus.climate_mean.raw[2]
}
genus3clim3.raw <- function (x){
  a_genus.climate_mean.raw[3] * (x - Topt_genus.climate_mean.raw[3])^2 + c_genus.climate_mean.raw[3]
}
genus4clim3.raw <- function (x){
  a_genus.climate_mean.raw[3] * (x - Topt_genus.climate_mean.raw[4])^2 + c_genus.climate_mean.raw[4]
}
genus1clim4.raw <- function (x){
  a_genus.climate_mean.raw[4] * (x - Topt_genus.climate_mean.raw[1])^2 + c_genus.climate_mean.raw[1]
}
genus2clim4.raw <- function (x){
  a_genus.climate_mean.raw[4] * (x - Topt_genus.climate_mean.raw[2])^2 + c_genus.climate_mean.raw[2]
}
genus3clim4.raw <- function (x){
  a_genus.climate_mean.raw[4] * (x - Topt_genus.climate_mean.raw[3])^2 + c_genus.climate_mean.raw[3]
}
genus4clim4.raw <- function (x){
  a_genus.climate_mean.raw[4] * (x - Topt_genus.climate_mean.raw[4])^2 + c_genus.climate_mean.raw[4]
}



# Data frame with temperature, fitted value, and genus group
fitted_data.genus.climate.raw <- data.frame(Temperature = rep(temp_seq, n.genus.raw * n.climate.raw), 
                                        Fitted = c(genus1clim1.raw(temp_seq),
                                                   genus2clim1.raw(temp_seq),
                                                   genus3clim1.raw(temp_seq),
                                                   genus4clim1.raw(temp_seq),
                                                   genus1clim2.raw(temp_seq),
                                                   genus2clim2.raw(temp_seq),
                                                   genus3clim2.raw(temp_seq),
                                                   genus4clim2.raw(temp_seq),
                                                   genus1clim3.raw(temp_seq),
                                                   genus2clim3.raw(temp_seq),
                                                   genus3clim3.raw(temp_seq),
                                                   genus4clim3.raw(temp_seq),
                                                   genus1clim4.raw(temp_seq),
                                                   genus2clim4.raw(temp_seq),
                                                   genus3clim4.raw(temp_seq),
                                                   genus4clim4.raw(temp_seq)),
                                        Genus = factor(rep(rep(c("Amblyomma", "Dermacentor", 
                                                                 "Haemaphysalis", "Hyalomma"), each = length(temp_seq)), n.climate.raw)),
                                        Climate = factor(rep(c("Mixed", "Subtropical", "Temperate", "Tropical"),
                                                             each = length(temp_seq) * n.genus.raw)))

# Original data with temperature, trait, numeric genus group, numeric species
plot_data.genus.climate.raw <- data.frame(
  Temperature = tick.raw$Temp,
  Trait = tick.raw$Trait,
  Genus = as.factor(tick.raw$Genus),
  Species = as.factor(tick.raw$Species),
  Climate = as.factor(tick.raw$Climate)
)

# Get unique combinations of species and genus for correct plotting
species_to_genus.climate.raw <- unique(plot_data.genus.climate.raw[, c("Species", "Genus", "Climate")])  
species_to_genus.climate.raw$Genus <- as.factor(species_to_genus.climate.raw$Genus)
species_to_genus.climate.raw$Climate <- as.factor(species_to_genus.climate.raw$Climate)

# Merge fitted genus curves with species-to-genus mapping
fitted_data.genus.climate.raw <- data.frame(Temperature = rep(temp_seq, n.genus.raw * n.climate.raw), 
                                        Fitted = c(genus1clim1.raw(temp_seq),
                                                   genus2clim1.raw(temp_seq),
                                                   genus3clim1.raw(temp_seq),
                                                   genus4clim1.raw(temp_seq),
                                                   genus1clim2.raw(temp_seq),
                                                   genus2clim2.raw(temp_seq),
                                                   genus3clim2.raw(temp_seq),
                                                   genus4clim2.raw(temp_seq),
                                                   genus1clim3.raw(temp_seq),
                                                   genus2clim3.raw(temp_seq),
                                                   genus3clim3.raw(temp_seq),
                                                   genus4clim3.raw(temp_seq),
                                                   genus1clim4.raw(temp_seq),
                                                   genus2clim4.raw(temp_seq),
                                                   genus3clim4.raw(temp_seq),
                                                   genus4clim4.raw(temp_seq)),
                                        Genus = factor(rep(rep(c("Amblyomma", "Dermacentor", 
                                                                 "Haemaphysalis", "Hyalomma"), each = length(temp_seq)), n.climate.raw)),
                                        Climate = factor(rep(c("Mixed", "Subtropical", "Temperate", "Tropical"),
                                                             each = length(temp_seq) * n.genus.raw))) %>%
  left_join(species_to_genus.climate.raw, by = c("Genus", "Climate"))|> 
  filter(!is.na(Species)) 


# Combine Climate and Size for coloring
fitted_data.genus.climate.raw$Combo <- with(fitted_data.genus.climate.raw, paste(Genus, Climate, sep = ", "))


# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.genus.climate.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.genus.climate.raw, aes(x = Temperature, y = Fitted, color = Combo), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Genus and Climate",
       x = "Temperature", y = "Ovipostion Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) 
#ggsave("tick.genus.climate.plot.png", width = 8, height = 5)

# Merge the observed and predicted data for comparison
merged_data.genus.climate.raw <- fitted_data.genus.climate.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species", 
                         "Genus" = "Genus", "Climate" = "Climate"))

# Calculate RMSE for each species
mse_by_genus.climate.raw <- merged_data.genus.climate.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.genus.climate.raw <- merged_data.genus.climate.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 


save(mse_by_genus.climate.raw,
     fitted_data.genus.climate.raw,
     merged_data.genus.climate.raw,
     genus.climate.raw.fit.mcmc,
     genus.climate.raw.fit,
     log.lik.genus.climate.raw,
     file="genus.climate.raw.fit.RData")

################################### By Host, Climate, and Size ###################################
# Model by typical number of host species
sink("host.climate.size_robin.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Host-specific 
    mu.a ~ dexp(0.01)            # Mean for host effects a (width)
  # Size-specific
    mu.c ~ dexp(10)          # Mean for host effects c (vertex y)
  # Climate-specific
  mu.Topt ~ dnorm(30, 0.1)          # Mean for cliamte effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for climate effect Topt (vertex x)

  for (j in 1:n.host.raw) {
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
  }
  for (k in 1:n.climate.raw) {
    mu.climate.Topt[k] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific vertex x
    climate.tau[k] ~ dexp(tau1)
  }
  for (s in 1:n.size.raw) {
    mu.size.c[s] ~ dexp(mu.c)       # Size-specific vertex y
  }
  # Likelihood
  for (i in 1:N.obs.raw) {
    mu[i] <- mu.host.a[host.raw[i]] * (mu.climate.Topt[climate.raw[i]] - temp[i])^2 + mu.size.c[size.raw[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[climate.raw[i]]) T(0, ) 
    log_lik[i] <- logdensity.norm(trait[i], mu[i], climate.tau[climate.raw[i]])
  }}
", file = "host.climate.size_robin.txt")


# Parameters to monitor
parameters <- c("mu.host.a", "mu.size.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt", "climate.tau",
                "tau1", "mu.climate.Topt", "log_lik")


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
jag.data <- list(trait = tick.raw$Trait, N.obs.raw = N.obs.raw, 
                 temp = tick.raw$Temp, 
                 n.host.raw = n.host.raw, host.raw = host.raw, 
                 n.climate.raw = n.climate.raw, climate.raw = climate.raw,
                 n.size.raw = n.size.raw, size.raw = size.raw)
# Fit model
h.c.s.raw.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                  model.file = "host.climate.size_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                  n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
h.c.s.raw.fit.mcmc <- as.mcmc(h.c.s.raw.fit) 
closeAllConnections()
h.c.s.raw.fit$BUGSoutput$DIC


# Extract log-likelihood matrix (samples × observations)
log.lik.h.c.s.raw <- h.c.s.raw.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.h.c.s.raw <- waic(log.lik.h.c.s.raw)$estimates["waic", "Estimate"]

# Get posterior means for parameters for each host group
summary_fit.h.c.s.raw <- summary(h.c.s.raw.fit.mcmc)$statistics
a_h.c.s_mean.raw <- summary_fit.h.c.s.raw[grep("^mu.host.a\\[", rownames(summary_fit.h.c.s.raw)), "Mean"]
c_h.c.s_mean.raw <- summary_fit.h.c.s.raw[grep("^mu.size.c\\[", rownames(summary_fit.h.c.s.raw)), "Mean"]
Topt_h.c.s_mean.raw <- summary_fit.h.c.s.raw[grep("^mu.climate.Topt\\[", rownames(summary_fit.h.c.s.raw)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp) ,
                max(tick.raw$Temp) , by = 0.1)

# Get predicted values for each  group
for (host in 1:3) {
  for (clim in 1:4) {
    for (size in 1:3) {
      func_name <- paste0("host", host, "clim", clim, "size", size, ".raw")
      assign(func_name, 
             local({
               h <- host
               c <- clim
               s <- size
               force(h); force(c); force(s)  # Ensure values are captured
               
               function(x) {
                 a_h.c.s_mean.raw[h] * (x - Topt_h.c.s_mean.raw[c])^2 + c_h.c.s_mean.raw[s]
               }
             }))
    }
  }
}

# Initialize empty lists to store values
fitted_list <- list()
host_labels <- c("One", "Three", "Two")
climate_labels <- c("Mixed", "Subtropical", "Temperate", "Tropical")
size_labels <- c("Large", "Medium", "Small")

# Loop through all host, climate, and size combinations
for (host in 1:3) {
  for (clim in 1:4) {
    for (size in 1:3) {
      # Name the function for appropriate host, climate, size combo
      func_name <- paste0("host", host, "clim", clim, "size", size, ".raw")
      # Get function by name and apply to temp_seq
      fitted_values <- get(func_name)(temp_seq)
      # Store fitted values
      fitted_list[[func_name]] <- data.frame(
        Temperature = temp_seq,
        Fitted = fitted_values,
        Host = factor(rep(host_labels[host], length(temp_seq))),
        Climate = factor(rep(climate_labels[clim], length(temp_seq))),
        Size = factor(rep(size_labels[size], length(temp_seq)))
      )
    }
  }
}
# Combine all individual data frames into one
fitted_data.h.c.s.raw <- do.call(rbind, fitted_list)

# Original data with temperature, trait, factor host group, species, climate, size
plot_data.h.c.s.raw <- data.frame(
  Temperature = tick.raw$Temp,
  Trait = tick.raw$Trait,
  Host = as.factor(tick.raw$Host),
  Species = as.factor(tick.raw$Species),
  Climate = as.factor(tick.raw$Climate),
  Size = as.factor(tick.raw$Size)
)

# Get unique combinations of species and host for correct plotting
species_to_h.c.s.raw <- unique(plot_data.h.c.s.raw[, c("Species", "Host", "Climate", "Size")])  
species_to_h.c.s.raw$Host <- as.factor(species_to_h.c.s.raw$Host)
species_to_h.c.s.raw$Climate <- as.factor(species_to_h.c.s.raw$Climate)
species_to_h.c.s.raw$Size <- as.factor(species_to_h.c.s.raw$Size)

# Merge fitted host curves with species-to-host mapping
fitted_data.h.c.s.raw <- fitted_data.h.c.s.raw |>  
  left_join(species_to_h.c.s.raw, by = c("Host", "Climate", "Size"))|> 
  filter(!is.na(Species)) 

# Order factor levels
plot_data.h.c.s.raw$Host <- factor(plot_data.h.c.s.raw$Host, levels = c("One", "Two", "Three"))
fitted_data.h.c.s.raw$Host <- factor(fitted_data.h.c.s.raw$Host, levels = c("One", "Two", "Three"))
plot_data.h.c.s.raw$Size <- factor(plot_data.h.c.s.raw$Size, levels = c("Small", "Medium", "Large"))
fitted_data.h.c.s.raw$Size <- factor(fitted_data.h.c.s.raw$Size, levels = c("Small", "Medium", "Large"))

# Combine Host, Climate, and Size for coloring
fitted_data.h.c.s.raw$Combo <- with(fitted_data.h.c.s.raw, paste(Host, Climate, Size, sep = ", "))
fitted_data.h.c.s.raw$Combo <- factor(fitted_data.h.c.s.raw$Combo, levels = c("One, Tropical, Medium",
                                                                      "Two, Subtropical, Medium",
                                                                      "Two, Subtropical, Large",
                                                                      "Three, Mixed, Small",
                                                                      "Three, Mixed, Large",
                                                                      "Three, Subtropical, Small",
                                                                      "Three, Subtropical, Large",
                                                                      "Three, Temperate, Large",
                                                                      "Three, Tropical, Large"))

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.h.c.s.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.h.c.s.raw, aes(x = Temperature, y = Fitted, color = Combo), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Host, Climate, and Size",
       x = "Temperature", y = "Oviposition Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) 

# Merge the observed and predicted data for comparison
merged_data.h.c.s.raw <- fitted_data.h.c.s.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species", 
                         "Host" = "Host", "Climate" = "Climate",
                         "Size" = "Size"))

# Calculate RMSE for each species
mse_by_h.c.s.raw <- merged_data.h.c.s.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.h.c.s.raw <- merged_data.h.c.s.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 


save(mse_by_h.c.s.raw,
     fitted_data.h.c.s.raw,
     merged_data.h.c.s.raw,
     h.c.s.raw.fit.mcmc,
     h.c.s.raw.fit,
     log.lik.h.c.s.raw,
     file="host.climate.size.raw.fit.RData")

################################### By Host, Climate, and Genus ###################################
# Model by typical number of host species
sink("host.climate.genus_robin.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Host-specific 
   mu.a ~ dexp(0.01)            # Mean for host effects a (width)
  # Genus-specific
    mu.c ~ dexp(10)          # Mean for host effects c (vertex y)
  # Climate-specific
  mu.Topt ~ dnorm(30, 0.1)          # Mean for cliamte effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for climate effect Topt (vertex x)

  for (j in 1:n.host.raw) {
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
  }
  for (k in 1:n.climate.raw) {
    mu.climate.Topt[k] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific vertex x
    climate.tau[k] ~ dexp(tau1)
  }
  for (s in 1:n.genus.raw) {
    mu.genus.c[s] ~ dexp(mu.c)       # genus-specific vertex y
  }
  # Likelihood
  for (i in 1:N.obs.raw) {
    mu[i] <- mu.host.a[host.raw[i]] * (mu.climate.Topt[climate.raw[i]] - temp[i])^2 + mu.genus.c[genus.raw[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[climate.raw[i]]) T(0, ) 
    log_lik[i] <- logdensity.norm(trait[i], mu[i], climate.tau[climate.raw[i]])
  }}
", file = "host.climate.genus_robin.txt")



# Parameters to monitor
parameters <- c("mu.host.a", "mu.genus.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt", "climate.tau",
                "tau1", "mu.climate.Topt", "log_lik")


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
jag.data <- list(trait = tick.raw$Trait, N.obs.raw = N.obs.raw, 
                 temp = tick.raw$Temp, 
                 n.host.raw = n.host.raw, host.raw = host.raw, 
                 n.climate.raw = n.climate.raw, climate.raw = climate.raw,
                 n.genus.raw = n.genus.raw, genus.raw = genus.raw)
# Fit model
h.c.g.raw.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                  model.file = "host.climate.genus_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                  n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
h.c.g.raw.fit.mcmc <- as.mcmc(h.c.g.raw.fit) 
closeAllConnections()
h.c.g.raw.fit$BUGSoutput$DIC


# Extract log-likelihood matrix (samples × observations)
log.lik.h.c.g.raw <- h.c.g.raw.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.h.c.g.raw <- waic(log.lik.h.c.g.raw)$estimates["waic", "Estimate"]

# Get posterior means for parameters for each host group
summary_fit.h.c.g.raw <- summary(h.c.g.raw.fit.mcmc)$statistics
a_h.c.g_mean.raw <- summary_fit.h.c.g.raw[grep("^mu.host.a\\[", rownames(summary_fit.h.c.g.raw)), "Mean"]
c_h.c.g_mean.raw <- summary_fit.h.c.g.raw[grep("^mu.genus.c\\[", rownames(summary_fit.h.c.g.raw)), "Mean"]
Topt_h.c.g_mean.raw <- summary_fit.h.c.g.raw[grep("^mu.climate.Topt\\[", rownames(summary_fit.h.c.g.raw)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp) ,
                max(tick.raw$Temp) , by = 0.1)

# Get predicted values for each  group
for (host in 1:3) {
  for (clim in 1:4) {
    for (genus in 1:4) {
      func_name <- paste0("host", host, "clim", clim, "genus", genus, ".raw")
      assign(func_name, 
             local({
               h <- host
               c <- clim
               g <- genus
               force(h); force(c); force(g)  # Ensure values are captured
               
               function(x) {
                 a_h.c.g_mean.raw[h] * (x - Topt_h.c.g_mean.raw[c])^2 + c_h.c.g_mean.raw[g]
               }
             }))
    }
  }
}

# Initialize empty lists to store values
fitted_list <- list()
host_labels <- c("One", "Three", "Two")
climate_labels <- c("Mixed", "Subtropical", "Temperate", "Tropical")
genus_labels <- c("Amblyomma", "Dermacentor", "Haemaphysalis", "Hyalomma")

# Loop through all host, climate, and size combinations
for (host in 1:3) {
  for (clim in 1:4) {
    for (genus in 1:4) {
      # Name the function for appropriate host, climate, size combo
      func_name <- paste0("host", host, "clim", clim, "genus", genus, ".raw")
      # Get function by name and apply to temp_seq
      fitted_values <- get(func_name)(temp_seq)
      # Store fitted values
      fitted_list[[func_name]] <- data.frame(
        Temperature = temp_seq,
        Fitted = fitted_values,
        Host = factor(rep(host_labels[host], length(temp_seq))),
        Climate = factor(rep(climate_labels[clim], length(temp_seq))),
        Genus = factor(rep(genus_labels[genus], length(temp_seq)))
      )
    }
  }
}
# Combine all individual data frames into one
fitted_data.h.c.g.raw <- do.call(rbind, fitted_list)

# Original data with temperature, trait, factor host group, species, climate, size
plot_data.h.c.g.raw <- data.frame(
  Temperature = tick.raw$Temp,
  Trait = tick.raw$Trait,
  Host = as.factor(tick.raw$Host),
  Species = as.factor(tick.raw$Species),
  Climate = as.factor(tick.raw$Climate),
  Genus = as.factor(tick.raw$Genus)
)

# Get unique combinations of species and host for correct plotting
species_to_h.c.g.raw <- unique(plot_data.h.c.g.raw[, c("Species", "Host", "Climate", "Genus")])  
species_to_h.c.g.raw$Host <- as.factor(species_to_h.c.g.raw$Host)
species_to_h.c.g.raw$Climate <- as.factor(species_to_h.c.g.raw$Climate)
species_to_h.c.g.raw$Genus <- as.factor(species_to_h.c.g.raw$Genus)

# Merge fitted host curves with species-to-host mapping
fitted_data.h.c.g.raw <- fitted_data.h.c.g.raw |>  
  left_join(species_to_h.c.g.raw, by = c("Host", "Climate", "Genus"))|> 
  filter(!is.na(Species)) 

# Order factor levels
plot_data.h.c.g.raw$Host <- factor(plot_data.h.c.s.raw$Host, levels = c("One", "Two", "Three"))
fitted_data.h.c.g.raw$Host <- factor(fitted_data.h.c.s.raw$Host, levels = c("One", "Two", "Three"))

# Combine Host, Climate, and Size for coloring
fitted_data.h.c.g.raw$Combo <- with(fitted_data.h.c.g.raw, paste(Host, Climate, Genus, sep = ", "))

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.h.c.g.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.h.c.g.raw, aes(x = Temperature, y = Fitted, color = Combo), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Host, Climate, and Genus",
       x = "Temperature", y = "Oviposition in Days") +
  theme_minimal() +
  facet_wrap(~ Species) 
#ggsave("tick.h.c.g.plot.png", width = 8, height = 5)

# Merge the observed and predicted data for comparison
merged_data.h.c.g.raw <- fitted_data.h.c.g.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species", 
                         "Host" = "Host", "Climate" = "Climate",
                         "Genus" = "Genus"))

# Calculate RMSE for each species
mse_by_h.c.g.raw <- merged_data.h.c.g.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.h.c.g.raw <- merged_data.h.c.g.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 


save(mse_by_h.c.g.raw,
     fitted_data.h.c.g.raw,
     merged_data.h.c.g.raw,
     h.c.g.raw.fit.mcmc,
     h.c.g.raw.fit,
     log.lik.h.c.g.raw,
     file="host.climate.genus.raw.fit.RData")

################################### Trace, Pairs, Histograms ###################################
# Make mcmc object a data frame, choosing columns with variables of interest
df.quad.priors <- as.data.frame(quad.fit.mcmc[[1]][,c(2, 3, 43, 44, 45, 46)])
# Make mcmc object a data frame for each parameter to compare difference from effects
df.quad.a <- as.data.frame(quad.fit.mcmc[[1]][,c(4:16)])
df.quad.Topt <- as.data.frame(quad.fit.mcmc[[1]][,c(17:29)])
df.quad.c <- as.data.frame(quad.fit.mcmc[[1]][,c(30:42)])
# Pairs plots
pairs(df.quad.priors)
pairs(df.quad.a)
pairs(df.quad.Topt)
pairs(df.quad.c)
# Trace plots
plot(quad.fit.mcmc)

# Mcmc object into list
samp.quad <- as.mcmc.list(quad.fit.mcmc)
# Parameter sets to plot a, Topt, and c
param_sets <- list(
  c(4:15),
  c(16:27),
  c(28:39),
  c(2, 3, 40, 41, 42, 43)
)

# Loop over each parameter set
for (set in param_sets) {
  par(mfrow = c(3,5))
  for (i in set) {
    hist(samp.quad[[1]][, i], xlab = paste("Parameter", i), 
         main = paste("Posterior Samples of Parameter", i))
  }
}


################################### RMSE, DIC ###################################
# Table of average RMSE across species and DIC of the model
table.tick.raw <- data.frame(
  Method = factor(rep(c("Full Data", "By Species Individual", "By Species Hierarchical", 
                        "By Species Hierarchical Different Variances", "By Species Different Variances, Same Parameters", 
                        "By Genus", 
                        "By Climate", "By Host", "By Size", "By Continents",
                        "By Host and Climate", "By Host and Genus", 
                        "By Genus and Climate", "By Host, Climate, and Size",
                        "By Host, Climate, and Genus"), each = 9)),
  DIC = rep(c(full.data.fit.raw$BUGSoutput$DIC, DIC.species.individual.raw, species.raw.fit$BUGSoutput$DIC,
              species.tau.raw.fit$BUGSoutput$DIC, species.tau.param.raw.fit$BUGSoutput$DIC,
              genus.raw.fit$BUGSoutput$DIC,
              climate.raw.fit$BUGSoutput$DIC, host.raw.fit$BUGSoutput$DIC, size.raw.fit$BUGSoutput$DIC, 
              continent.raw.fit$BUGSoutput$DIC,
              host.climate.raw.fit$BUGSoutput$DIC, host.genus.raw.fit$BUGSoutput$DIC,
              genus.climate.raw.fit$BUGSoutput$DIC,
              h.c.s.raw.fit$BUGSoutput$DIC, h.c.g.raw.fit$BUGSoutput$DIC), each = 9),
  wAIC = rep(c(waic.full.raw, waic.species.ind.raw, waic.species.raw,
               waic.species.tau.raw, waic.species.tau.param.raw,
               waic.genus.raw, waic.climate.raw, waic.host.raw,
               waic.size.raw, waic.continent.raw,
               waic.host.climate.raw, waic.host.genus.raw,
               waic.genus.climate.raw, 
               waic.h.c.s.raw, waic.h.c.g.raw), each = 9),
  Species = factor(rep(c("Amblyomma lepidum", "Dermacentor andersoni",
                         "Dermacentor nitens", "Haemaphysalis leporispalustris",
                         "Hyalomma aegyptium", "Hyalomma dromedarii",
                         "Hyalomma impeltatum", "Hyalomma lusitanicum",
                         "Hyalomma schulzei"), times = 15)),
  # MSE = c(mse_by_full_data.raw$MSE, as.numeric(rmse.a.lepidum.raw$MSE), as.numeric(rmse.d.andersoni.raw$MSE),
  #         as.numeric(rmse.d.nitens.raw$MSE), 
  #         as.numeric(rmse.h.leporispalustris.raw$MSE), as.numeric(rmse.h.aegyptium.raw$MSE), 
  #         as.numeric(rmse.h.dromedarii.raw$MSE), as.numeric(rmse.h.impeltatum.raw$MSE), 
  #         as.numeric(rmse.h.lusitanicum.raw$MSE), as.numeric(rmse.h.schulzei.raw$MSE), 
  #         mse_by_species.raw$MSE, mse_by_species.tau.raw$MSE, 
  #         mse_by_species.tau.param.raw$MSE, mse_by_genus.raw$MSE,
  #         mse_by_climate.raw$MSE, mse_by_host.raw$MSE, mse_by_size.raw$MSE,
  #         mse_by_continent.raw$MSE, mse_by_host.climate.raw$MSE,
  #         mse_by_host.genus.raw$MSE, mse_by_genus.climate.raw$MSE,
  #         mse_by_h.c.s.raw$MSE, mse_by_h.c.g.raw$MSE),
  Total.RMSE = as.numeric(rep(c(overall_rmse.raw, overall_rmse.species.ind.raw,
                     overall_rmse.species.raw, overall_rmse.species.tau.raw, 
                     overall_rmse.species.tau.param.raw, overall_rmse.genus.raw, 
                     overall_rmse.climate.raw, overall_rmse.host.raw, overall_rmse.size.raw, 
                     overall_rmse.continent.raw, overall_rmse.host.climate.raw, 
                     overall_rmse.host.genus.raw,
                     overall_rmse.genus.climate.raw,
                     overall_rmse.h.c.s.raw, overall_rmse.h.c.g.raw), each = 9))
)

table.tick.avg.raw <- data.frame(
  Method = factor(c("Full Data", "By Species Individual", "By Species Hierarchical", 
                    "By Species Hierarchical Different Variances", 
                    "By Species Different Variances, Same Parameters", 
                    "By Genus", 
                    "By Climate", "By Host", "By Size", "By Continents",
                    "By Host and Climate", "By Host and Genus", 
                    "By Genus and Climate", "By Host, Climate, and Size",
                    "By Host, Climate, and Genus")),
  DIC = c(full.data.fit.raw$BUGSoutput$DIC, DIC.species.individual.raw, species.raw.fit$BUGSoutput$DIC,
          species.tau.raw.fit$BUGSoutput$DIC, species.tau.param.raw.fit$BUGSoutput$DIC,
          genus.raw.fit$BUGSoutput$DIC,
          climate.raw.fit$BUGSoutput$DIC, host.raw.fit$BUGSoutput$DIC, size.raw.fit$BUGSoutput$DIC, 
          continent.raw.fit$BUGSoutput$DIC,
          host.climate.raw.fit$BUGSoutput$DIC, host.genus.raw.fit$BUGSoutput$DIC,
          genus.climate.raw.fit$BUGSoutput$DIC,
          h.c.s.raw.fit$BUGSoutput$DIC, h.c.g.raw.fit$BUGSoutput$DIC),
  wAIC = c(waic.full.raw, waic.species.ind.raw, waic.species.raw,
           waic.species.tau.raw, waic.species.tau.param.raw,
           waic.genus.raw, waic.climate.raw, waic.host.raw,
           waic.size.raw, waic.continent.raw,
           waic.host.climate.raw, waic.host.genus.raw,
           waic.genus.climate.raw, 
           waic.h.c.s.raw, waic.h.c.g.raw),
  Total.RMSE = round(as.numeric(c(overall_rmse.raw, overall_rmse.species.ind.raw,
                                  overall_rmse.species.raw, overall_rmse.species.tau.raw, 
                                  overall_rmse.species.tau.param.raw, overall_rmse.genus.raw, 
                                  overall_rmse.climate.raw, overall_rmse.host.raw, overall_rmse.size.raw, 
                                  overall_rmse.continent.raw, overall_rmse.host.climate.raw, 
                                  overall_rmse.host.genus.raw,
                                  overall_rmse.genus.climate.raw,
                                  overall_rmse.h.c.s.raw, overall_rmse.h.c.g.raw)), 3)
)

table.tick.avg.raw$IC <- table.tick.avg.raw$DIC + table.tick.avg.raw$wAIC
# Unhelpful bar plot with species on x-axis, RMSE on y, colored by method
table.tick.raw$Method <- reorder(table.tick.raw$Method, -table.tick.raw$Total.RMSE)
ggplot(table.tick.raw, aes(x = Species, y = Total.RMSE, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  labs(title = "RMSE by Species and Method", x = "Species",
       y = "Root Mean Square Error (RMSE)", fill = "Method") +
  theme_minimal() +
  scale_fill_viridis_d(option = "mako") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

table.tick.avg.raw$Method <- reorder(table.tick.avg.raw$Method, -table.tick.avg.raw$Total.RMSE)
ggplot(table.tick.avg.raw, aes(x = Method, y = Total.RMSE, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  #labs(title = "DIC by Model", x = "Model",
  #    y = "Deviance Information Criterion (DIC)") +
  labs(title = "RMSE by Model", x = "Model",
       y = "Root Mean Square Error (RMSE)") +
  #coord_cartesian(ylim = c(6.5, 12)) +
  theme_minimal() +
  scale_fill_viridis_d(option = "mako") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
#ggsave("tick.rm.rawSE.png", width = 8, height = 5)


table.tick.raw$Method <- reorder(table.tick.raw$Method, -table.tick.raw$DIC)
ggplot(table.tick.raw, aes(x = Method, y = DIC, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  labs(title = "RMSE by Species and Method", x = "Method",
       y = "Root Mean Square Error (RMSE)") +
  theme_minimal() +
  #coord_cartesian(ylim = c(9900, 10700)) +
  scale_fill_viridis_d(option = "mako") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

save(table.tick.raw,
     table.tick.avg.raw, file = "Dataset.Tick.Raw.RData")


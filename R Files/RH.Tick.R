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
tick <- read.csv("Extra.Data.Tick.csv")
# Convert all characteristics to factors
tick$Species <- as.factor(tick$Species)
tick$Host <- as.factor(tick$Host)
tick$Genus <- as.factor(tick$Genus)
tick$Climate <- as.factor(tick$Climate)
tick$Continent <- as.factor(tick$Continent)
tick$Size <- as.factor(tick$Size)

# Plot by species
ggplot(tick, aes(y = Trait, x = Temp)) +
  geom_point() +
  facet_wrap(~ Species) +
  labs(title = "Temperature by Trait for Each Factor Level",
       x = "Temperature",
       y = "Oviposition") +
  theme_minimal()

# Removing the species that aren't perfectly cleaned for now
tick <- tick[!(tick$Species %in% c("Ixodes anatis", "Rhipicephalus (Boophilus) annulatus", "Haemaphysalis longicornis")), ]
# Drop the factor levels of species that aren't used anymore
tick$Species <- droplevels(tick$Species)

# Create a new expanded dataset with pseudo-distribution data
expanded_df1 <- data.frame()

# Get the data that is NOT individual
filtered_df1 <- tick |> filter(Count != 1)

# Loop over each row in the dataset of mean data
for (i in 1:nrow(filtered_df1)) {
  # Extract values for each row
  # Mean
  mean_value <- filtered_df1$Trait[i]
  # Standard error times the sqrt(n) to get sd
  if (filtered_df1$Error.Type[i] %in% c("se", "SE")){
    sd_value <- filtered_df1$Error[i] * sqrt(filtered_df1$Count[i])
  }
  else{
    sd_value <- filtered_df1$Error[i]
  }
  # Species
  interactor1 <- filtered_df1$Species[i]
  # Genus 
  Genus <- filtered_df1$Genus[i]
  # Temp
  interactor1temp <- filtered_df1$Temp[i]
  # Sample size
  num_points <- filtered_df1$Count[i]
  # Host #
  Host <- filtered_df1$Host[i]
  # Climate
  Climate <- filtered_df1$Climate[i]
  # Number of Continents
  Continent <- filtered_df1$Continent[i]
  # Size
  Size <- filtered_df1$Size[i]
  # Generate sample size of random values from truncated normal (0, Inf) distribution
  generated_values <- rtruncnorm(num_points, a=0, b=Inf, mean = mean_value, sd = sd_value)
  # Create a temporary dataframe for the generated values
  temp_df <- data.frame(
    # Simulated values
    Trait = generated_values,
    # Carry over sd, species, temperature, genus, host, climate, continent, and size
    sd = sd_value,
    Species = interactor1,
    Temp = interactor1temp,
    Genus = Genus,
    Host = Host,
    Climate = Climate,
    Continent = Continent,
    Size = Size
  )
  # Append the generated rows to the expanded dataset
  expanded_df1 <- rbind(expanded_df1, temp_df)
}

# Take the individual data
filtered_df <- tick |> filter(Count == 1)
# Summarize it, grouping by species and temperature, to get the mean, sd, and sample size
summarized_df <- filtered_df %>%
  group_by(Genus, Species, Temp, Host, Climate, Continent, Size) %>%
  summarize(
    mean_value = mean(Trait, na.rm = TRUE),
    sd_value = sd(Trait, na.rm = TRUE),
    count = n(),
    .groups = "drop"
  )
# Resample n data points with appropriate mean and standard deviation
expanded_df <- summarized_df %>%
  rowwise() %>%
  do(data.frame(
    Species = .$Species,
    Temp = .$Temp,
    Trait = rnorm(.$count, mean = .$mean_value, sd = .$sd_value),
    Genus = .$Genus,
    Host = .$Host,
    Climate = .$Climate,
    Continent = .$Continent,
    Size = .$Size  
  ))

# Combine the new mean data and the new individual data, selecting 
tick.abr.new <- rbind(expanded_df1[, c(1, 3:9)], expanded_df)
# Get rid of NAs
tick.abr.new <- na.omit(tick.abr.new)
# Drop unused levels of species
tick.abr.new$Species <- droplevels(tick.abr.new$Species)

# Plot new data
ggplot(tick.abr.new, aes(y = Trait, x = Temp)) +
  geom_point() +
  facet_wrap(~ Species) +
  labs(title = "Temperature by Trait for Each Factor Level",
       x = "Temperature",
       y = "Oviposition") +
  theme_minimal()


# Species number for hierarchical model
species <- as.factor(tick.abr.new$Species)
# Number of species
n.species <- length(unique(tick.abr.new$Species))

# Genus number for hierarchical model
genus <- as.factor(tick.abr.new$Genus)
# Number of genus
n.genus <- length(unique(tick.abr.new$Genus))

# Host number for hierarchical model
host <- as.factor(tick.abr.new$Host)
# Number of hosts
n.host <- length(unique(tick.abr.new$Host))

# Climate number for hierarchical model
climate <- as.factor(tick.abr.new$Climate)
# Number of climates
n.climate <- length(unique(tick.abr.new$Climate))

# Continent number for hierarchical model
continent <- as.factor(tick.abr.new$Continent)
# Number of continents
n.continent <- length(unique(tick.abr.new$Continent))

# Size number for hierarchical model
size <- as.factor(tick.abr.new$Size)
# Small: 0-3
# Medium: 3.01-5
# Large: 5.01+
# Number of sizes
n.size <- length(unique(tick.abr.new$Size))

N.obs <- length(tick.abr.new$Trait)
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
  for (i in 1:N.obs) {
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
jag.data <- list(trait = tick.abr.new$Trait, N.obs = N.obs, 
                 temp = tick.abr.new$Temp)

# Fit fill model
full.data.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                    model.file = "full_data_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                    n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
full.data.fit.mcmc <- as.mcmc(full.data.fit) 
closeAllConnections()
#DIC
full.data.fit$BUGSoutput$DIC

# Extract log-likelihood matrix (samples × observations)
log.lik.full <- full.data.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.full <- waic(log.lik.full)$estimates["waic", "Estimate"]


# Get posterior means
summary_fit <- summary(full.data.fit.mcmc)$statistics
a_mean <- summary_fit["a", "Mean"]
c_mean <- summary_fit["c", "Mean"]
Topt_mean <- summary_fit["Topt", "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp), max(tick.abr.new$Temp), by = 0.1)

# Get predicted values
full.data <- function (x){
  a_mean * (x - Topt_mean)^2 + c_mean
}


# Data frame of temperature, fitted value, species
fitted_data.full.data <- data.frame(Temperature = rep(temp_seq, n.species), 
                                   Fitted = full.data(temp_seq),
                                   Species = factor(rep(c("Amblyomma lepidum", "Dermacentor andersoni",
                                                          "Dermacentor nitens", "Haemaphysalis leporispalustris",
                                                          "Hyalomma aegyptium", "Hyalomma dromedarii",
                                                          "Hyalomma impeltatum", "Hyalomma lusitanicum",
                                                          "Hyalomma schulzei"), each = length(temp_seq))))
# Original data with temperature, trait, and species as numeric
plot_data.full <- data.frame(
  Temperature = tick.abr.new$Temp,
  Trait = tick.abr.new$Trait,
  Species = as.factor(tick.abr.new$Species)
)

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.full, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.full.data, aes(x = Temperature, y = Fitted), color = "lightslateblue", size = 1) +
  labs(x = "Temperature", y = "Oviposition Time in Days") +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~ Species) 
ggsave("tick.full.plot.png", width = 8, height = 5)

# Merge the observed and predicted data for comparison
merged_data.full <- fitted_data.full.data %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species"))

# Calculate RMSE for each species
mse_by_full_data <- merged_data.full %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse <- merged_data.full %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 


# save(mse_by_full_data,
#      overall_rmse,
#      fitted_data.full.data,
#      merged_data.full,
#      full.data.fit.mcmc,
#      full.data.fit,
#      log.lik.full,
#      file="full.data.fit.RData")


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

  for (j in 1:n.species) {
    mu.species.c[j] ~ dexp(mu.c)       # Species-specific vertex y
    mu.species.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Species-specific vertex x
    mu.species.a[j] ~ dexp(mu.a)    # Species-specific width
  }
  # Likelihood
  for (i in 1:N.obs) {
    mu[i] <- mu.species.a[species[i]] * (mu.species.Topt[species[i]] - temp[i])^2 + mu.species.c[species[i]]
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
jag.data <- list(trait = tick.abr.new$Trait, N.obs = N.obs, 
               temp = tick.abr.new$Temp, 
               n.species = n.species, species = species)
# Fit model
species.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                 model.file = "species_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                 n.iter = ni, DIC = T, working.directory = getwd())
# Make an mcmc object
species.fit.mcmc <- as.mcmc(species.fit) 
closeAllConnections()
species.fit$BUGSoutput$DIC

# Extract log-likelihood matrix (samples × observations)
log.lik.species <- species.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.species <- waic(log.lik.species)$estimates["waic", "Estimate"]



# Get posterior means for each species effect
summary_fit.species <- summary(species.fit.mcmc)$statistics
a_species_mean <- summary_fit.species[grep("^mu.species.a\\[", rownames(summary_fit.species)), "Mean"]
c_species_mean <- summary_fit.species[grep("^mu.species.c\\[", rownames(summary_fit.species)), "Mean"]
Topt_species_mean <- summary_fit.species[grep("^mu.species.Topt\\[", rownames(summary_fit.species)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp) ,
                max(tick.abr.new$Temp) , by = 0.1)

# Get predicted values for each species using species effect parameters
species1 <- function (x){
  a_species_mean[1] * (x - Topt_species_mean[1])^2 + c_species_mean[1]
}
species2 <- function (x){
  a_species_mean[2] * (x - Topt_species_mean[2])^2 + c_species_mean[2]
}
species3 <- function (x){
  a_species_mean[3] * (x - Topt_species_mean[3])^2 + c_species_mean[3]
}
species4 <- function (x){
  a_species_mean[4] * (x - Topt_species_mean[4])^2 + c_species_mean[4]
}
species5 <- function (x){
  a_species_mean[5] * (x - Topt_species_mean[5])^2 + c_species_mean[5]
}
species6 <- function (x){
  a_species_mean[6] * (x - Topt_species_mean[6])^2 + c_species_mean[6]
}
species7 <- function (x){
  a_species_mean[7] * (x - Topt_species_mean[7])^2 + c_species_mean[7]
}
species8 <- function (x){
  a_species_mean[8] * (x - Topt_species_mean[8])^2 + c_species_mean[8]
}
species9 <- function (x){
  a_species_mean[9] * (x - Topt_species_mean[9])^2 + c_species_mean[9]
}

# Make data frame of temperature, fitted value for each species, and species
fitted_data.species <- data.frame(Temperature = rep(temp_seq, n.species), 
                                  Fitted = c(species1(temp_seq),
                                             species2(temp_seq),
                                             species3(temp_seq),
                                             species4(temp_seq),
                                             species5(temp_seq),
                                             species6(temp_seq),
                                             species7(temp_seq),
                                             species8(temp_seq),
                                             species9(temp_seq)),
                                  Species = factor(rep(c("Amblyomma lepidum", "Dermacentor andersoni",
                                                         "Dermacentor nitens", "Haemaphysalis leporispalustris",
                                                         "Hyalomma aegyptium", "Hyalomma dromedarii",
                                                         "Hyalomma impeltatum", "Hyalomma lusitanicum",
                                                         "Hyalomma schulzei"), each = length(temp_seq))))
# Data frame of original temperature, trait, and numeric species
plot_data.species <- data.frame(
  Temperature = tick.abr.new$Temp,
  Trait = tick.abr.new$Trait,
  Species = as.factor(tick.abr.new$Species)
)

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.species, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.species, aes(x = Temperature, y = Fitted, color = Species), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Species",
       x = "Temperature", y = "Ovipositon in Days") +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~ Species) 
ggsave("tick.species.plot.png", width = 8, height = 5)

# Merge the observed and predicted data for comparison
merged_data.species <- fitted_data.species %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species"))

# Calculate RMSE for each species
mse_by_species <- merged_data.species %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )


# Calculate overall RMSE (no grouping)
overall_rmse.species <- merged_data.species %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

# 
# save(mse_by_species,
#      fitted_data.species,
#      merged_data.species,
#      species.fit.mcmc,
#      species.fit,
#      log.lik.species,
#      file="species.fit.RData")

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

  for (j in 1:n.species) {
    mu.species.c[j] ~ dexp(mu.c)       # Species-specific vertex y
    mu.species.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Species-specific vertex x
    mu.species.a[j] ~ dexp(mu.a)    # Species-specific width
    species.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs) {
    mu[i] <- mu.species.a[species[i]] * (mu.species.Topt[species[i]] - temp[i])^2 + mu.species.c[species[i]]
    trait[i] ~ dnorm(mu[i], species.tau[species[i]]) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], species.tau[species[i]])
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
jag.data <- list(trait = tick.abr.new$Trait, N.obs = N.obs, 
                 temp = tick.abr.new$Temp, 
                 n.species = n.species, species = species)
# Fit model
species.tau.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                    model.file = "species_tau_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                    n.iter = ni, DIC = T, working.directory = getwd())
# Make an mcmc object
species.tau.fit.mcmc <- as.mcmc(species.tau.fit) 
closeAllConnections()
species.tau.fit$BUGSoutput$DIC

# Extract log-likelihood matrix (samples × observations)
log.lik.species.tau <- species.tau.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.species.tau <- waic(log.lik.species.tau)$estimates["waic", "Estimate"]


# Get posterior means for each species effect
summary_fit.species.tau <- summary(species.tau.fit.mcmc)$statistics
a_species.tau_mean <- summary_fit.species.tau[grep("^mu.species.a\\[", rownames(summary_fit.species.tau)), "Mean"]
c_species.tau_mean <- summary_fit.species.tau[grep("^mu.species.c\\[", rownames(summary_fit.species.tau)), "Mean"]
Topt_species.tau_mean <- summary_fit.species.tau[grep("^mu.species.Topt\\[", rownames(summary_fit.species.tau)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp) ,
                max(tick.abr.new$Temp) , by = 0.1)

# Get predicted values for each species using species effect parameters
species1.tau <- function (x){
  a_species.tau_mean[1] * (x - Topt_species.tau_mean[1])^2 + c_species.tau_mean[1]
}
species2.tau <- function (x){
  a_species.tau_mean[2] * (x - Topt_species.tau_mean[2])^2 + c_species.tau_mean[2]
}
species3.tau <- function (x){
  a_species.tau_mean[3] * (x - Topt_species.tau_mean[3])^2 + c_species.tau_mean[3]
}
species4.tau <- function (x){
  a_species.tau_mean[4] * (x - Topt_species.tau_mean[4])^2 + c_species.tau_mean[4]
}
species5.tau <- function (x){
  a_species.tau_mean[5] * (x - Topt_species.tau_mean[5])^2 + c_species.tau_mean[5]
}
species6.tau <- function (x){
  a_species.tau_mean[6] * (x - Topt_species.tau_mean[6])^2 + c_species.tau_mean[6]
}
species7.tau <- function (x){
  a_species.tau_mean[7] * (x - Topt_species.tau_mean[7])^2 + c_species.tau_mean[7]
}
species8.tau <- function (x){
  a_species.tau_mean[8] * (x - Topt_species.tau_mean[8])^2 + c_species.tau_mean[8]
}
species9.tau <- function (x){
  a_species.tau_mean[9] * (x - Topt_species.tau_mean[9])^2 + c_species.tau_mean[9]
}

# Make data frame of temperature, fitted value for each species, and species
fitted_data.species.tau <- data.frame(Temperature = rep(temp_seq, n.species), 
                                  Fitted = c(species1.tau(temp_seq),
                                             species2.tau(temp_seq),
                                             species3.tau(temp_seq),
                                             species4.tau(temp_seq),
                                             species5.tau(temp_seq),
                                             species6.tau(temp_seq),
                                             species7.tau(temp_seq),
                                             species8.tau(temp_seq),
                                             species9.tau(temp_seq)),
                                  Species = factor(rep(c("Amblyomma lepidum", "Dermacentor andersoni",
                                                         "Dermacentor nitens", "Haemaphysalis leporispalustris",
                                                         "Hyalomma aegyptium", "Hyalomma dromedarii",
                                                         "Hyalomma impeltatum", "Hyalomma lusitanicum",
                                                         "Hyalomma schulzei"), each = length(temp_seq))))
# Data frame of original temperature, trait, and numeric species
plot_data.species <- data.frame(
  Temperature = tick.abr.new$Temp,
  Trait = tick.abr.new$Trait,
  Species = as.factor(tick.abr.new$Species)
)

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.species, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.species.tau, aes(x = Temperature, y = Fitted, color = Species), size = 1) +
  labs(x = "Temperature", y = "Oviposition Time in Days") +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~ Species) 
ggsave("tick.species.tau.plot.png", width = 8, height = 5)


# Merge the observed and predicted data for comparison
merged_data.species.tau <- fitted_data.species.tau %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species"))

# Calculate RMSE for each species
mse_by_species.tau <- merged_data.species.tau %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.species.tau <- merged_data.species.tau %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

# save(mse_by_species.tau,
#      fitted_data.species.tau,
#      merged_data.species.tau,
#      species.tau.fit.mcmc,
#      species.tau.fit,
#      log.lik.species.tau,
#      file="species.tau.fit.RData")

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

  for (j in 1:n.species) {
    species.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], species.tau[species[i]]) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], species.tau[species[i]])
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
jag.data <- list(trait = tick.abr.new$Trait, N.obs = N.obs, 
                 temp = tick.abr.new$Temp, 
                 n.species = n.species, species = species)
# Fit model
species.tau.param.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                        model.file = "species_tau.param_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                        n.iter = ni, DIC = T, working.directory = getwd())
# Make an mcmc object
species.tau.param.fit.mcmc <- as.mcmc(species.tau.param.fit) 
closeAllConnections()
species.tau.param.fit$BUGSoutput$DIC

# Extract log-likelihood matrix (samples × observations)
log.lik.species.tau.param <- species.tau.param.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.species.tau.param <- waic(log.lik.species.tau.param)$estimates["waic", "Estimate"]


# Get posterior means for each species effect
summary_fit.species.tau.param <- summary(species.tau.param.fit.mcmc)$statistics
a_species.tau.param_mean <- summary_fit.species.tau.param["a", "Mean"]
c_species.tau.param_mean <- summary_fit.species.tau.param["c", "Mean"]
Topt_species.tau.param_mean <- summary_fit.species.tau.param["Topt", "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp) ,
                max(tick.abr.new$Temp) , by = 0.1)

# Get predicted values
species.tau.param <- function (x){
  a_species.tau.param_mean * (x - Topt_species.tau.param_mean)^2 + c_species.tau.param_mean
}

# Make data frame of temperature, fitted value for each species, and species
fitted_data.species.tau.param <- data.frame(Temperature = rep(temp_seq, n.species), 
                                            Fitted = species.tau.param(temp_seq),
                                            Species = factor(rep(c("Amblyomma lepidum", "Dermacentor andersoni",
                                                                   "Dermacentor nitens", "Haemaphysalis leporispalustris",
                                                                   "Hyalomma aegyptium", "Hyalomma dromedarii",
                                                                   "Hyalomma impeltatum", "Hyalomma lusitanicum",
                                                                   "Hyalomma schulzei"), each = length(temp_seq))))
# Data frame of original temperature, trait, and numeric species
plot_data.species <- data.frame(
  Temperature = tick.abr.new$Temp,
  Trait = tick.abr.new$Trait,
  Species = as.factor(tick.abr.new$Species)
)

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.species, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.species.tau.param, aes(x = Temperature, y = Fitted, color = Species), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Species",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~ Species) 

# Merge the observed and predicted data for comparison
merged_data.species.tau.param <- fitted_data.species.tau.param %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species"))

# Calculate RMSE for each species
mse_by_species.tau.param <- merged_data.species.tau.param %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.species.tau.param <- merged_data.species.tau.param %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

# save(mse_by_species.tau.param,
#      fitted_data.species.tau.param,
#      merged_data.species.tau.param,
#      species.tau.param.fit.mcmc,
#      species.tau.param.fit,
#      log.lik.species.tau.param,
#      file="species.tau.param.fit.RData")

################################### By Species Individual- Lepidum ###################################
a.lepidum <- tick.abr.new[tick.abr.new$Species == "Amblyomma lepidum", ]

N.obs.ind <- length(a.lepidum$Trait)

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
  for (i in 1:N.obs.ind) {
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
jag.data <- list(trait = a.lepidum$Trait, N.obs.ind = N.obs.ind, 
                 temp = a.lepidum$Temp)

# Fit fill model
a.lepidum.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                       model.file = "species_individual_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                       n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
a.lepidum.fit.mcmc <- as.mcmc(a.lepidum.fit) 
closeAllConnections()
#DIC
a.lepidum.fit$BUGSoutput$DIC

# Get posterior means
summary_fit.a.lepidum <- summary(a.lepidum.fit.mcmc)$statistics
a_mean <- summary_fit.a.lepidum["a", "Mean"]
c_mean <- summary_fit.a.lepidum["c", "Mean"]
Topt_mean <- summary_fit.a.lepidum["Topt", "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp), max(tick.abr.new$Temp), by = 0.1)
# 
# # Get predicted values
a.lepidum.fitted <- function (x){
  a_mean * (x - Topt_mean)^2 + c_mean
}


# Data frame of temperature, fitted value, species
fitted_data.a.lepidum <- data.frame(Temperature = rep(temp_seq),
                                     Fitted = a.lepidum.fitted(temp_seq))
# Original data with temperature, trait, and species as numeric
plot_data.a.lepidum <- data.frame(
  Temperature = a.lepidum$Temp,
  Trait = a.lepidum$Trait
)

# Plot observed data and fitted curves 
p.a.lepidum <- ggplot() +
  geom_point(data = plot_data.a.lepidum, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.a.lepidum, aes(x = Temperature, y = Fitted), color = "seagreen1", size = 1) +
  labs(title = "Observed Data and Quadratic Fit",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  theme(legend.position = "none")
p.a.lepidum

# Merge the observed and predicted data for comparison
merged_data.a.lepidum <- fitted_data.a.lepidum %>%
  left_join(a.lepidum, by = c("Temperature" = "Temp"))

# Calculate RMSE for each species
rmse.a.lepidum <- merged_data.a.lepidum %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE),
    SE = (Trait - Fitted)^2 
  )


# Extract log-likelihood matrix (samples × observations)
log.lik.a.lepidum <- a.lepidum.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.a.lepidum <- waic(log.lik.a.lepidum)$estimates["waic", "Estimate"]


save(rmse.a.lepidum,
     fitted_data.a.lepidum,
     merged_data.a.lepidum,
     a.lepidum.fit.mcmc,
     a.lepidum.fit,
     p.a.lepidum,
     waic.a.lepidum,
     file="a.lepidum.fit.RData")

################################### By Species Individual- Andersoni ###################################
d.andersoni <- tick.abr.new[tick.abr.new$Species == "Dermacentor andersoni", ]


N.obs.ind <- length(d.andersoni$Trait)

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
  for (i in 1:N.obs.ind) {
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
jag.data <- list(trait = d.andersoni$Trait, N.obs.ind = N.obs.ind, 
                 temp = d.andersoni$Temp)

# Fit fill model
d.andersoni.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                       model.file = "species_individual_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                       n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
d.andersoni.fit.mcmc <- as.mcmc(d.andersoni.fit) 
closeAllConnections()
#DIC
d.andersoni.fit$BUGSoutput$DIC

# Get posterior means
summary_fit.d.andersoni <- summary(d.andersoni.fit.mcmc)$statistics
a_mean <- summary_fit.d.andersoni["a", "Mean"]
c_mean <- summary_fit.d.andersoni["c", "Mean"]
Topt_mean <- summary_fit.d.andersoni["Topt", "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp), max(tick.abr.new$Temp), by = 0.1)
# 
# # Get predicted values
d.andersoni.fitted <- function (x){
  a_mean * (x - Topt_mean)^2 + c_mean
}


# Data frame of temperature, fitted value, species
fitted_data.d.andersoni <- data.frame(Temperature = rep(temp_seq),
                                     Fitted = d.andersoni.fitted(temp_seq))
# Original data with temperature, trait, and species as numeric
plot_data.d.andersoni <- data.frame(
  Temperature = d.andersoni$Temp,
  Trait = d.andersoni$Trait
)

# Plot observed data and fitted curves 
p.d.andersoni <- ggplot() +
  geom_point(data = plot_data.d.andersoni, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.d.andersoni, aes(x = Temperature, y = Fitted), color = "seagreen1", size = 1) +
  labs(title = "Observed Data and Quadratic Fit",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  theme(legend.position = "none")
p.d.andersoni

# Merge the observed and predicted data for comparison
merged_data.d.andersoni <- fitted_data.d.andersoni %>%
  left_join(d.andersoni, by = c("Temperature" = "Temp"))

# Calculate RMSE for each species
rmse.d.andersoni <- merged_data.d.andersoni %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE),
    SE = (Trait - Fitted)^2 
  )


# Extract log-likelihood matrix (samples × observations)
log.lik.d.andersoni <- d.andersoni.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.d.andersoni <- waic(log.lik.d.andersoni)$estimates["waic", "Estimate"]

save(rmse.d.andersoni,
     fitted_data.d.andersoni,
     merged_data.d.andersoni,
     d.andersoni.fit.mcmc,
     d.andersoni.fit,
     p.d.andersoni,
     waic.d.andersoni,
     file="d.andersoni.fit.RData")

################################### By Species Individual- Nitens ###################################
d.nitens <- tick.abr.new[tick.abr.new$Species == "Dermacentor nitens", ]


N.obs.ind <- length(d.nitens$Trait)

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
  for (i in 1:N.obs.ind) {
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
jag.data <- list(trait = d.nitens$Trait, N.obs.ind = N.obs.ind, 
                 temp = d.nitens$Temp)

# Fit fill model
d.nitens.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                       model.file = "species_individual_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                       n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
d.nitens.fit.mcmc <- as.mcmc(d.nitens.fit) 
closeAllConnections()
#DIC
d.nitens.fit$BUGSoutput$DIC

# Get posterior means
summary_fit.d.nitens <- summary(d.nitens.fit.mcmc)$statistics
a_mean <- summary_fit.d.nitens["a", "Mean"]
c_mean <- summary_fit.d.nitens["c", "Mean"]
Topt_mean <- summary_fit.d.nitens["Topt", "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp), max(tick.abr.new$Temp), by = 0.1)
# 
# # Get predicted values
d.nitens.fitted <- function (x){
  a_mean * (x - Topt_mean)^2 + c_mean
}


# Data frame of temperature, fitted value, species
fitted_data.d.nitens <- data.frame(Temperature = rep(temp_seq),
                                     Fitted = d.nitens.fitted(temp_seq))
# Original data with temperature, trait, and species as numeric
plot_data.d.nitens <- data.frame(
  Temperature = d.nitens$Temp,
  Trait = d.nitens$Trait
)

# Plot observed data and fitted curves 
p.d.nitens <- ggplot() +
  geom_point(data = plot_data.d.nitens, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.d.nitens, aes(x = Temperature, y = Fitted), color = "seagreen1", size = 1) +
  labs(title = "Observed Data and Quadratic Fit",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  theme(legend.position = "none")
p.d.nitens

# Merge the observed and predicted data for comparison
merged_data.d.nitens <- fitted_data.d.nitens %>%
  left_join(d.nitens, by = c("Temperature" = "Temp"))

# Calculate RMSE for each species
rmse.d.nitens <- merged_data.d.nitens %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE),
    SE = (Trait - Fitted)^2 
  )


# Extract log-likelihood matrix (samples × observations)
log.lik.d.nitens <- d.nitens.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.d.nitens <- waic(log.lik.d.nitens)$estimates["waic", "Estimate"]


save(rmse.d.nitens,
     fitted_data.d.nitens,
     merged_data.d.nitens,
     d.nitens.fit.mcmc,
     d.nitens.fit,
     p.d.nitens,
     waic.d.nitens,
     file="d.nitens.fit.RData")

################################### By Species Individual- Leporispalustris ###################################
h.leporispalustris <- tick.abr.new[tick.abr.new$Species == "Haemaphysalis leporispalustris", ]

N.obs.ind <- length(h.leporispalustris$Trait)

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
  for (i in 1:N.obs.ind) {
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
jag.data <- list(trait = h.leporispalustris$Trait, N.obs.ind = N.obs.ind, 
                 temp = h.leporispalustris$Temp)

# Fit fill model
h.leporispalustris.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                       model.file = "species_individual_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                       n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
h.leporispalustris.fit.mcmc <- as.mcmc(h.leporispalustris.fit) 
closeAllConnections()
#DIC
h.leporispalustris.fit$BUGSoutput$DIC

# Get posterior means
summary_fit.h.leporispalustris <- summary(h.leporispalustris.fit.mcmc)$statistics
a_mean <- summary_fit.h.leporispalustris["a", "Mean"]
c_mean <- summary_fit.h.leporispalustris["c", "Mean"]
Topt_mean <- summary_fit.h.leporispalustris["Topt", "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp), max(tick.abr.new$Temp), by = 0.1)
# 
# # Get predicted values
h.leporispalustris.fitted <- function (x){
  a_mean * (x - Topt_mean)^2 + c_mean
}


# Data frame of temperature, fitted value, species
fitted_data.h.leporispalustris <- data.frame(Temperature = rep(temp_seq),
                                     Fitted = h.leporispalustris.fitted(temp_seq))
# Original data with temperature, trait, and species as numeric
plot_data.h.leporispalustris <- data.frame(
  Temperature = h.leporispalustris$Temp,
  Trait = h.leporispalustris$Trait
)

# Plot observed data and fitted curves 
p.h.leporispalustris <- ggplot() +
  geom_point(data = plot_data.h.leporispalustris, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.h.leporispalustris, aes(x = Temperature, y = Fitted), color = "seagreen1", size = 1) +
  labs(title = "Observed Data and Quadratic Fit",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  theme(legend.position = "none")
p.h.leporispalustris

# Merge the observed and predicted data for comparison
merged_data.h.leporispalustris <- fitted_data.h.leporispalustris %>%
  left_join(h.leporispalustris, by = c("Temperature" = "Temp"))

# Calculate RMSE for each species
rmse.h.leporispalustris <- merged_data.h.leporispalustris %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE),
    SE = (Trait - Fitted)^2 
  )


# Extract log-likelihood matrix (samples × observations)
log.lik.h.leporispalustris <- h.leporispalustris.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.h.leporispalustris <- waic(log.lik.h.leporispalustris)$estimates["waic", "Estimate"]

save(rmse.h.leporispalustris,
     fitted_data.h.leporispalustris,
     merged_data.h.leporispalustris,
     h.leporispalustris.fit.mcmc,
     h.leporispalustris.fit,
     p.h.leporispalustris,
     waic.h.leporispalustris,
     file="h.leporispalustris.fit.RData")

################################### By Species Individual- Aegyptium ###################################
h.aegyptium <- tick.abr.new[tick.abr.new$Species == "Hyalomma aegyptium", ]


N.obs.ind <- length(h.aegyptium$Trait)

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
  for (i in 1:N.obs.ind) {
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
jag.data <- list(trait = h.aegyptium$Trait, N.obs.ind = N.obs.ind, 
                 temp = h.aegyptium$Temp)

# Fit fill model
h.aegyptium.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                       model.file = "species_individual_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                       n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
h.aegyptium.fit.mcmc <- as.mcmc(h.aegyptium.fit) 
closeAllConnections()
#DIC
h.aegyptium.fit$BUGSoutput$DIC

# Get posterior means
summary_fit.h.aegyptium <- summary(h.aegyptium.fit.mcmc)$statistics
a_mean <- summary_fit.h.aegyptium["a", "Mean"]
c_mean <- summary_fit.h.aegyptium["c", "Mean"]
Topt_mean <- summary_fit.h.aegyptium["Topt", "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp), max(tick.abr.new$Temp), by = 0.1)
# 
# # Get predicted values
h.aegyptium.fitted <- function (x){
  a_mean * (x - Topt_mean)^2 + c_mean
}


# Data frame of temperature, fitted value, species
fitted_data.h.aegyptium <- data.frame(Temperature = rep(temp_seq),
                                     Fitted = h.aegyptium.fitted(temp_seq))
# Original data with temperature, trait, and species as numeric
plot_data.h.aegyptium <- data.frame(
  Temperature = h.aegyptium$Temp,
  Trait = h.aegyptium$Trait
)

# Plot observed data and fitted curves 
p.h.aegyptium <- ggplot() +
  geom_point(data = plot_data.h.aegyptium, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.h.aegyptium, aes(x = Temperature, y = Fitted), color = "seagreen1", size = 1) +
  labs(title = "Observed Data and Quadratic Fit",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  theme(legend.position = "none")
p.h.aegyptium

# Merge the observed and predicted data for comparison
merged_data.h.aegyptium <- fitted_data.h.aegyptium %>%
  left_join(h.aegyptium, by = c("Temperature" = "Temp"))

# Calculate RMSE for each species
rmse.h.aegyptium <- merged_data.h.aegyptium %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE),
    SE = (Trait - Fitted)^2 
  )


# Extract log-likelihood matrix (samples × observations)
log.lik.h.aegyptium <- h.aegyptium.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.h.aegyptium <- waic(log.lik.h.aegyptium)$estimates["waic", "Estimate"]

save(rmse.h.aegyptium,
     fitted_data.h.aegyptium,
     merged_data.h.aegyptium,
     h.aegyptium.fit.mcmc,
     h.aegyptium.fit,
     p.h.aegyptium,
     waic.h.aegyptium,
     file="h.aegyptium.fit.RData")

################################### By Species Individual- Dromedarii ###################################
h.dromedarii <- tick.abr.new[tick.abr.new$Species == "Hyalomma dromedarii", ]

N.obs.ind <- length(h.dromedarii$Trait)

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
  for (i in 1:N.obs.ind) {
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
jag.data <- list(trait = h.dromedarii$Trait, N.obs.ind = N.obs.ind, 
                 temp = h.dromedarii$Temp)

# Fit fill model
h.dromedarii.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                       model.file = "species_individual_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                       n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
h.dromedarii.fit.mcmc <- as.mcmc(h.dromedarii.fit) 
closeAllConnections()
#DIC
h.dromedarii.fit$BUGSoutput$DIC

# Get posterior means
summary_fit.h.dromedarii <- summary(h.dromedarii.fit.mcmc)$statistics
a_mean <- summary_fit.h.dromedarii["a", "Mean"]
c_mean <- summary_fit.h.dromedarii["c", "Mean"]
Topt_mean <- summary_fit.h.dromedarii["Topt", "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp), max(tick.abr.new$Temp), by = 0.1)
# 
# # Get predicted values
h.dromedarii.fitted <- function (x){
  a_mean * (x - Topt_mean)^2 + c_mean
}


# Data frame of temperature, fitted value, species
fitted_data.h.dromedarii <- data.frame(Temperature = rep(temp_seq),
                                     Fitted = h.dromedarii.fitted(temp_seq))
# Original data with temperature, trait, and species as numeric
plot_data.h.dromedarii <- data.frame(
  Temperature = h.dromedarii$Temp,
  Trait = h.dromedarii$Trait
)

# Plot observed data and fitted curves 
p.h.dromedarii <- ggplot() +
  geom_point(data = plot_data.h.dromedarii, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.h.dromedarii, aes(x = Temperature, y = Fitted), color = "seagreen1", size = 1) +
  labs(title = "Observed Data and Quadratic Fit",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  theme(legend.position = "none")
p.h.dromedarii

# Merge the observed and predicted data for comparison
merged_data.h.dromedarii <- fitted_data.h.dromedarii %>%
  left_join(h.dromedarii, by = c("Temperature" = "Temp"))

# Calculate RMSE for each species
rmse.h.dromedarii <- merged_data.h.dromedarii %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE),
    SE = (Trait - Fitted)^2 
  )


# Extract log-likelihood matrix (samples × observations)
log.lik.h.dromedarii <- h.dromedarii.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.h.dromedarii <- waic(log.lik.h.dromedarii)$estimates["waic", "Estimate"]


save(rmse.h.dromedarii,
     fitted_data.h.dromedarii,
     merged_data.h.dromedarii,
     h.dromedarii.fit.mcmc,
     h.dromedarii.fit,
     p.h.dromedarii,
     waic.h.dromedarii,
     file="h.dromedarii.fit.RData")

################################### By Species Individual- Impeltatum ###################################
h.impeltatum <- tick.abr.new[tick.abr.new$Species == "Hyalomma impeltatum", ]


N.obs.ind <- length(h.impeltatum$Trait)

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
  for (i in 1:N.obs.ind) {
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
jag.data <- list(trait = h.impeltatum$Trait, N.obs.ind = N.obs.ind, 
                 temp = h.impeltatum$Temp)

# Fit fill model
h.impeltatum.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                       model.file = "species_individual_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                       n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
h.impeltatum.fit.mcmc <- as.mcmc(h.impeltatum.fit) 
closeAllConnections()
#DIC
h.impeltatum.fit$BUGSoutput$DIC

# Get posterior means
summary_fit.h.impeltatum <- summary(h.impeltatum.fit.mcmc)$statistics
a_mean <- summary_fit.h.impeltatum["a", "Mean"]
c_mean <- summary_fit.h.impeltatum["c", "Mean"]
Topt_mean <- summary_fit.h.impeltatum["Topt", "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp), max(tick.abr.new$Temp), by = 0.1)
# 
# # Get predicted values
h.impeltatum.fitted <- function (x){
  a_mean * (x - Topt_mean)^2 + c_mean
}


# Data frame of temperature, fitted value, species
fitted_data.h.impeltatum <- data.frame(Temperature = rep(temp_seq),
                                     Fitted = h.impeltatum.fitted(temp_seq))
# Original data with temperature, trait, and species as numeric
plot_data.h.impeltatum <- data.frame(
  Temperature = h.impeltatum$Temp,
  Trait = h.impeltatum$Trait
)

# Plot observed data and fitted curves 
p.h.impeltatum <- ggplot() +
  geom_point(data = plot_data.h.impeltatum, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.h.impeltatum, aes(x = Temperature, y = Fitted), color = "seagreen1", size = 1) +
  labs(title = "Observed Data and Quadratic Fit",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  theme(legend.position = "none")
p.h.impeltatum

# Merge the observed and predicted data for comparison
merged_data.h.impeltatum <- fitted_data.h.impeltatum %>%
  left_join(h.impeltatum, by = c("Temperature" = "Temp"))

# Calculate RMSE for each species
rmse.h.impeltatum <- merged_data.h.impeltatum %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE),
    SE = (Trait - Fitted)^2 
  )


# Extract log-likelihood matrix (samples × observations)
log.lik.h.impeltatum <- h.impeltatum.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.h.impeltatum <- waic(log.lik.h.impeltatum)$estimates["waic", "Estimate"]

save(rmse.h.impeltatum,
     fitted_data.h.impeltatum,
     merged_data.h.impeltatum,
     h.impeltatum.fit.mcmc,
     h.impeltatum.fit,
     p.h.impeltatum,
     waic.h.impeltatum,
     file="h.impeltatum.fit.RData")

################################### By Species Individual- Lusitanicum ###################################
h.lusitanicum <- tick.abr.new[tick.abr.new$Species == "Hyalomma lusitanicum", ]

N.obs.ind <- length(h.lusitanicum$Trait)

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
  for (i in 1:N.obs.ind) {
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
jag.data <- list(trait = h.lusitanicum$Trait, N.obs.ind = N.obs.ind, 
                 temp = h.lusitanicum$Temp)

# Fit fill model
h.lusitanicum.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                       model.file = "species_individual_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                       n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
h.lusitanicum.fit.mcmc <- as.mcmc(h.lusitanicum.fit) 
closeAllConnections()
#DIC
h.lusitanicum.fit$BUGSoutput$DIC

# Get posterior means
summary_fit.h.lusitanicum <- summary(h.lusitanicum.fit.mcmc)$statistics
a_mean <- summary_fit.h.lusitanicum["a", "Mean"]
c_mean <- summary_fit.h.lusitanicum["c", "Mean"]
Topt_mean <- summary_fit.h.lusitanicum["Topt", "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp), max(tick.abr.new$Temp), by = 0.1)
# 
# # Get predicted values
h.lusitanicum.fitted <- function (x){
  a_mean * (x - Topt_mean)^2 + c_mean
}


# Data frame of temperature, fitted value, species
fitted_data.h.lusitanicum <- data.frame(Temperature = rep(temp_seq),
                                     Fitted = h.lusitanicum.fitted(temp_seq))
# Original data with temperature, trait, and species as numeric
plot_data.h.lusitanicum <- data.frame(
  Temperature = h.lusitanicum$Temp,
  Trait = h.lusitanicum$Trait
)

# Plot observed data and fitted curves 
p.h.lusitanicum <- ggplot() +
  geom_point(data = plot_data.h.lusitanicum, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.h.lusitanicum, aes(x = Temperature, y = Fitted), color = "seagreen1", size = 1) +
  labs(title = "Observed Data and Quadratic Fit",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  theme(legend.position = "none")
p.h.lusitanicum

# Merge the observed and predicted data for comparison
merged_data.h.lusitanicum <- fitted_data.h.lusitanicum %>%
  left_join(h.lusitanicum, by = c("Temperature" = "Temp"))

# Calculate RMSE for each species
rmse.h.lusitanicum <- merged_data.h.lusitanicum %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE),
    SE = (Trait - Fitted)^2 
  )


# Extract log-likelihood matrix (samples × observations)
log.lik.h.lusitanicum <- h.lusitanicum.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.h.lusitanicum <- waic(log.lik.h.lusitanicum)$estimates["waic", "Estimate"]


save(rmse.h.lusitanicum,
     fitted_data.h.lusitanicum,
     merged_data.h.lusitanicum,
     h.lusitanicum.fit.mcmc,
     h.lusitanicum.fit,
     p.h.lusitanicum,
     waic.h.lusitanicum,
     file="h.lusitanicum.fit.RData")


################################### By Species Individual- Schulzei ###################################
h.schulzei <- tick.abr.new[tick.abr.new$Species == "Hyalomma schulzei", ]

N.obs.ind <- length(h.schulzei$Trait)

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
  for (i in 1:N.obs.ind) {
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
jag.data <- list(trait = h.schulzei$Trait, N.obs.ind = N.obs.ind, 
                 temp = h.schulzei$Temp)

# Fit fill model
h.schulzei.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                      model.file = "species_individual_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                      n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
h.schulzei.fit.mcmc <- as.mcmc(h.schulzei.fit) 
closeAllConnections()
#DIC
h.schulzei.fit$BUGSoutput$DIC

# Get posterior means
summary_fit.h.schulzei <- summary(h.schulzei.fit.mcmc)$statistics
a_mean <- summary_fit.h.schulzei["a", "Mean"]
c_mean <- summary_fit.h.schulzei["c", "Mean"]
Topt_mean <- summary_fit.h.schulzei["Topt", "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp), max(tick.abr.new$Temp), by = 0.1)
# 
# # Get predicted values
h.schulzei.fitted <- function (x){
  a_mean * (x - Topt_mean)^2 + c_mean
}


# Data frame of temperature, fitted value, species
fitted_data.h.schulzei <- data.frame(Temperature = rep(temp_seq),
                                    Fitted = h.schulzei.fitted(temp_seq))
# Original data with temperature, trait, and species as numeric
plot_data.h.schulzei <- data.frame(
  Temperature = h.schulzei$Temp,
  Trait = h.schulzei$Trait
)

# Plot observed data and fitted curves 
p.h.schulzei <- ggplot() +
  geom_point(data = plot_data.h.schulzei, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.h.schulzei, aes(x = Temperature, y = Fitted), color = "seagreen1", size = 1) +
  labs(title = "Observed Data and Quadratic Fit",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  theme(legend.position = "none")
p.h.schulzei

# Merge the observed and predicted data for comparison
merged_data.h.schulzei <- fitted_data.h.schulzei %>%
  left_join(h.schulzei, by = c("Temperature" = "Temp"))

# Calculate RMSE for each species
rmse.h.schulzei <- merged_data.h.schulzei %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE),
    SE = (Trait - Fitted)^2 
  )


# Extract log-likelihood matrix (samples × observations)
log.lik.h.schulzei <- h.schulzei.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.h.schulzei <- waic(log.lik.h.schulzei)$estimates["waic", "Estimate"]


save(rmse.h.schulzei,
     fitted_data.h.schulzei,
     merged_data.h.schulzei,
     h.schulzei.fit.mcmc,
     h.schulzei.fit,
     p.h.schulzei,
     waic.h.schulzei,
     file="h.schulzei.fit.RData")

################################### By Species Individual- Summary ###################################

waic.species.ind <- waic(log.lik.a.lepidum)$estimates["waic", "Estimate"] + 
  waic(log.lik.d.andersoni)$estimates["waic", "Estimate"] +
  waic(log.lik.d.nitens)$estimates["waic", "Estimate"] +
  waic(log.lik.h.leporispalustris)$estimates["waic", "Estimate"] + 
  waic(log.lik.h.aegyptium)$estimates["waic", "Estimate"] +
  waic(log.lik.h.dromedarii)$estimates["waic", "Estimate"] + 
  waic(log.lik.h.impeltatum)$estimates["waic", "Estimate"] +
  waic(log.lik.h.lusitanicum)$estimates["waic", "Estimate"] +
  waic(log.lik.h.schulzei)$estimates["waic", "Estimate"]

DIC.species.individual <- a.lepidum.fit$BUGSoutput$DIC + d.andersoni.fit$BUGSoutput$DIC + d.nitens.fit$BUGSoutput$DIC +
  h.leporispalustris.fit$BUGSoutput$DIC + h.aegyptium.fit$BUGSoutput$DIC + h.dromedarii.fit$BUGSoutput$DIC +
  h.impeltatum.fit$BUGSoutput$DIC + h.lusitanicum.fit$BUGSoutput$DIC + h.schulzei.fit$BUGSoutput$DIC

# Define RMSE values for each species
se.individual <- c(rmse.a.lepidum$SE, rmse.d.andersoni$SE, rmse.d.nitens$SE,
                   rmse.h.leporispalustris$SE, rmse.h.aegyptium$SE,
                   rmse.h.dromedarii$SE, rmse.h.impeltatum$SE,
                   rmse.h.lusitanicum$SE, rmse.h.schulzei$SE)

# Compute weighted average RMSE
overall_rmse.species.ind <- sqrt(mean(se.individual, na.rm = TRUE))

# Arrange plots into a 3x3 grid
grid.ind.species <- grid.arrange(
  p.a.lepidum + labs(title = "Amblyomma lepidum", x = " ", y = "Oviposition Time in Days") + 
    theme(axis.title.y = element_text(size = 8), plot.title = element_text(size = 10)) +
    coord_cartesian(ylim = c(0,100)),
  p.d.andersoni + labs(title = "Dermacentor andersoni", x = " ", y = " ") + 
    theme(axis.title.y = element_text(size = 8), plot.title = element_text(size = 10)) +
    coord_cartesian(ylim = c(0,100)),
  p.d.nitens + labs(title = "Dermacentor nitens", x = " ", y = " ") + 
    theme(axis.title.y = element_text(size = 8), plot.title = element_text(size = 10)) +
    coord_cartesian(ylim = c(0,100)),
  p.h.leporispalustris + labs(title = "Haemaphysalis leporispalustris", x = " ", y = "Oviposition Time in Days") + 
    theme(axis.title.y = element_text(size = 8), plot.title = element_text(size = 10)) +
    coord_cartesian(ylim = c(0,100)),
  p.h.aegyptium + labs(title = "Hyalomma aegyptium", x = " ", y = " ") + 
    theme(axis.title.y = element_text(size = 8), plot.title = element_text(size = 10)) +
    coord_cartesian(ylim = c(0,100)),
  p.h.dromedarii + labs(title = "Hyalomma dromedarii", x = " ", y = " ") + 
    theme(axis.title.y = element_text(size = 8), plot.title = element_text(size = 10)) +
    coord_cartesian(ylim = c(0,100)),
  p.h.impeltatum + labs(title = "Hyalomma impeltatum", y = "Oviposition Time in Days") + 
    theme(axis.title.y = element_text(size = 8), plot.title = element_text(size = 10)) +
    coord_cartesian(ylim = c(0,100)),
  p.h.lusitanicum + labs(title = "Hyalomma lusitanicum", y = " ") + 
    theme(axis.title.y = element_text(size = 8), plot.title = element_text(size = 10)) +
    coord_cartesian(ylim = c(0,100)),
  p.h.schulzei + labs(title = "Hyalomma schulzei", y = " ") + 
    theme(axis.title.y = element_text(size = 8), plot.title = element_text(size = 10)) +
    coord_cartesian(ylim = c(0,100)),
  nrow = 3, ncol = 3
)
ggsave("tick.species.ind.plot.png", grid.ind.species, width = 8, height = 5)

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

  for (j in 1:n.genus) {
    mu.genus.c[j] ~ dexp(mu.c)       # Genus-specific vertex y
    mu.genus.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Genus-specific vertex x
    mu.genus.a[j] ~ dexp(mu.a)    # Genus-specific width
    genus.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs) {
    mu[i] <- mu.genus.a[genus[i]] * (mu.genus.Topt[genus[i]] - temp[i])^2 + mu.genus.c[genus[i]]
    trait[i] ~ dnorm(mu[i], genus.tau[genus[i]]) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], genus.tau[genus[i]])
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
jag.data <- list(trait = tick.abr.new$Trait, N.obs = N.obs, 
                 temp = tick.abr.new$Temp, 
                 n.genus = n.genus, genus = genus)
# Fit model
genus.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                  model.file = "genus_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                  n.iter = ni, DIC = T, working.directory = getwd())
# Make an mcmc object
genus.fit.mcmc <- as.mcmc(genus.fit)
closeAllConnections()
genus.fit$BUGSoutput$DIC

# Extract log-likelihood matrix (samples × observations)
log.lik.genus <- genus.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.genus <- waic(log.lik.genus)$estimates["waic", "Estimate"]




# Get posterior means using parameters from genus effect
summary_fit.genus <- summary(genus.fit.mcmc)$statistics
a_genus_mean <- summary_fit.genus[grep("^mu.genus.a\\[", rownames(summary_fit.genus)), "Mean"]
c_genus_mean <- summary_fit.genus[grep("^mu.genus.c\\[", rownames(summary_fit.genus)), "Mean"]
Topt_genus_mean <- summary_fit.genus[grep("^mu.genus.Topt\\[", rownames(summary_fit.genus)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp) ,
                max(tick.abr.new$Temp) , by = 0.1)

# Get predicted values for each genus
genus1 <- function (x){
  a_genus_mean[1] * (x - Topt_genus_mean[1])^2 + c_genus_mean[1]
}
genus2 <- function (x){
  a_genus_mean[2] * (x - Topt_genus_mean[2])^2 + c_genus_mean[2]
}
genus3 <- function (x){
  a_genus_mean[3] * (x - Topt_genus_mean[3])^2 + c_genus_mean[3]
}
genus4 <- function (x){
  a_genus_mean[4] * (x - Topt_genus_mean[4])^2 + c_genus_mean[4]
}

# Data frame of temperature, fitted values, and genus
fitted_data.genus <- data.frame(Temperature = rep(temp_seq, n.genus), 
                                Fitted = c(genus1(temp_seq),
                                           genus2(temp_seq),
                                           genus3(temp_seq),
                                           genus4(temp_seq)),
                                Genus = factor(rep(c("Amblyomma", "Dermacentor",
                                                     "Haemaphysalis", "Hyalomma"), each = length(temp_seq))))
# Original temperature, trait, numeric genus, and numeric species to plot by species
plot_data.genus <- data.frame(
  Temperature = tick.abr.new$Temp,
  Trait = tick.abr.new$Trait,
  Genus = as.factor(tick.abr.new$Genus),
  Species = as.factor(tick.abr.new$Species)
)

# Get unique combos of species and genus to align the correct curve on the plot
species_to_genus <- unique(plot_data.genus[, c("Species", "Genus")])
species_to_genus$Genus <- as.factor(species_to_genus$Genus)
# Merge fitted genus curves with species-to-genus mapping
fitted_data.genus <- data.frame(
  Temperature = rep(temp_seq, n.genus),
  Fitted = c(genus1(temp_seq),
             genus2(temp_seq),
             genus3(temp_seq),
             genus4(temp_seq)),
  Genus = factor(rep(c("Amblyomma", "Dermacentor",
                       "Haemaphysalis", "Hyalomma"), each = length(temp_seq)))) %>%
  left_join(species_to_genus, by = "Genus")

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.genus, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.genus, aes(x = Temperature, y = Fitted, color = Genus), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Genus",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~ Species) 

# Merge the observed and fitted data for comparison
merged_data.genus <- fitted_data.genus %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Genus" = "Genus"))

# Calculate RMSE for each species
mse_by_genus <- merged_data.genus %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.genus <- merged_data.genus %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 



# save(mse_by_genus,
#      fitted_data.genus,
#      merged_data.genus,
#      genus.fit.mcmc,
#      genus.fit,
#      log.lik.genus,
#      file="genus.fit.RData")

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

  for (j in 1:n.climate) {
    mu.climate.c[j] ~ dexp(mu.c)       # Climate-specific vertex y
    mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific vertex x
    mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific width
    climate.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs) {
    mu[i] <- mu.climate.a[climate[i]] * (mu.climate.Topt[climate[i]] - temp[i])^2 + mu.climate.c[climate[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[climate[i]]) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], climate.tau[climate[i]])
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
jag.data <- list(trait = tick.abr.new$Trait, N.obs = N.obs, 
                 temp = tick.abr.new$Temp, 
                 n.climate = n.climate, climate = climate)
# Fit model
climate.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                    model.file = "climate_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                    n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
climate.fit.mcmc <- as.mcmc(climate.fit) 
closeAllConnections()
climate.fit$BUGSoutput$DIC


# Extract log-likelihood matrix (samples × observations)
log.lik.climate <- climate.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.climate <- waic(log.lik.climate)$estimates["waic", "Estimate"]

# Get posterior means for parameters for each climate group
summary_fit.climate <- summary(climate.fit.mcmc)$statistics
a_climate_mean <- summary_fit.climate[grep("^mu.climate.a\\[", rownames(summary_fit.climate)), "Mean"]
c_climate_mean <- summary_fit.climate[grep("^mu.climate.c\\[", rownames(summary_fit.climate)), "Mean"]
Topt_climate_mean <- summary_fit.climate[grep("^mu.climate.Topt\\[", rownames(summary_fit.climate)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp) ,
                max(tick.abr.new$Temp) , by = 0.1)

# Get predicted values for each climate group
climate1 <- function (x){
  a_climate_mean[1] * (x - Topt_climate_mean[1])^2 + c_climate_mean[1]
}
climate2 <- function (x){
  a_climate_mean[2] * (x - Topt_climate_mean[2])^2 + c_climate_mean[2]
}
climate3 <- function (x){
  a_climate_mean[3] * (x - Topt_climate_mean[3])^2 + c_climate_mean[3]
}
climate4 <- function (x){
  a_climate_mean[4] * (x - Topt_climate_mean[4])^2 + c_climate_mean[4]
}

# Data frame of temperature, fitted values, and climate
fitted_data.climate <- data.frame(Temperature = rep(temp_seq, n.climate), 
                                  Fitted = c(climate1(temp_seq),
                                             climate2(temp_seq),
                                             climate3(temp_seq),
                                             climate4(temp_seq)),
                                  Climate = factor(rep(c("Mixed", "Subtropical", 
                                                         "Temperate", "Tropical"), each = length(temp_seq))))
# Original data with temperature, trait, numeric climate, and numeric species
plot_data.climate <- data.frame(
  Temperature = tick.abr.new$Temp,
  Trait = tick.abr.new$Trait,
  Climate = as.factor(tick.abr.new$Climate),
  Species = as.factor(tick.abr.new$Species)
)

# Get unique combinations of species and climate
species_to_climate <- unique(plot_data.climate[, c("Species", "Climate")])  
species_to_climate$Climate <- as.factor(species_to_climate$Climate)
# Merge fitted climate curves with species-to-climate mapping
fitted_data.climate <- data.frame(
  Temperature = rep(temp_seq, n.climate),
  Fitted = c(climate1(temp_seq),
             climate2(temp_seq),
             climate3(temp_seq),
             climate4(temp_seq)),
  Climate = factor(rep(c("Mixed", "Subtropical", 
                         "Temperate", "Tropical"), each = length(temp_seq)))) %>%
  left_join(species_to_climate, by = "Climate")

plot_data.climate$Climate <- factor(plot_data.climate$Climate, levels = c("Mixed", "Temperate", "Subtropical", "Tropical"))
fitted_data.climate$Climate <- factor(fitted_data.climate$Climate, levels = c("Mixed", "Temperate", "Subtropical", "Tropical"))
# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.climate, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.climate, aes(x = Temperature, y = Fitted, color = Climate), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Climate",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 

# Merge the observed and predicted data for comparison
merged_data.climate <- fitted_data.climate %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Climate" = "Climate"))

# Calculate RMSE for each species
mse_by_climate <- merged_data.climate %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.climate <- merged_data.climate %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

# 
# save(mse_by_climate,
#      fitted_data.climate,
#      merged_data.climate,
#      climate.fit.mcmc,
#      climate.fit,
#      log.lik.climate,
#      file="climate.fit.RData")

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
  
  for (j in 1:n.host) {
    mu.host.c[j] ~ dexp(mu.c)       # Host-specific vertex y
    mu.host.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Host-specific vertex x
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
    host.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs) {
    mu[i] <- mu.host.a[host[i]] * (mu.host.Topt[host[i]] - temp[i])^2 + mu.host.c[host[i]]
    trait[i] ~ dnorm(mu[i], host.tau[host[i]]) T(0, ) 
    log_lik[i] <- logdensity.norm(trait[i], mu[i], host.tau[host[i]])
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
jag.data <- list(trait = tick.abr.new$Trait, N.obs = N.obs, 
                 temp = tick.abr.new$Temp, 
                 n.host = n.host, host = host)
# Fit model
host.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                  model.file = "host_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                  n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.fit.mcmc <- as.mcmc(host.fit) 
closeAllConnections()
host.fit$BUGSoutput$DIC


# Extract log-likelihood matrix (samples × observations)
log.lik.host <- host.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.host <- waic(log.lik.host)$estimates["waic", "Estimate"]

# Get posterior means for parameters for each host group
summary_fit.host <- summary(host.fit.mcmc)$statistics
a_host_mean <- summary_fit.host[grep("^mu.host.a\\[", rownames(summary_fit.host)), "Mean"]
c_host_mean <- summary_fit.host[grep("^mu.host.c\\[", rownames(summary_fit.host)), "Mean"]
Topt_host_mean <- summary_fit.host[grep("^mu.host.Topt\\[", rownames(summary_fit.host)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp) ,
                max(tick.abr.new$Temp) , by = 0.1)

# Get predicted values for each host group
host1 <- function (x){
  a_host_mean[1] * (x - Topt_host_mean[1])^2 + c_host_mean[1]
}
host2 <- function (x){
  a_host_mean[2] * (x - Topt_host_mean[2])^2 + c_host_mean[2]
}
host3 <- function (x){
  a_host_mean[3] * (x - Topt_host_mean[3])^2 + c_host_mean[3]
}


# Data frame with temperature, fitted value, and host group
fitted_data.host <- data.frame(Temperature = rep(temp_seq, n.host), 
                                Fitted = c(host1(temp_seq),
                                           host2(temp_seq),
                                           host3(temp_seq)),
                                Host = factor(rep(c("One", "Three", "Two"), each = length(temp_seq))))
# Original data with temperature, trait, numeric host group, numeric species
plot_data.host <- data.frame(
  Temperature = tick.abr.new$Temp,
  Trait = tick.abr.new$Trait,
  Host = as.factor(tick.abr.new$Host),
  Species = as.factor(tick.abr.new$Species)
)

# Get unique combinations of species and host for correct plotting
species_to_host <- unique(plot_data.host[, c("Species", "Host")])  
species_to_host$Host <- as.factor(species_to_host$Host)
# Merge fitted host curves with species-to-host mapping
fitted_data.host <- data.frame(
  Temperature = rep(temp_seq, n.host),
  Fitted = c(host1(temp_seq),
             host2(temp_seq),
             host3(temp_seq)),
  Host = factor(rep(c("One", "Three", "Two"), each = length(temp_seq)))
) %>%
  left_join(species_to_host, by = "Host")

plot_data.host$Host <- factor(plot_data.host$Host, levels = c("One", "Two", "Three"))
fitted_data.host$Host <- factor(fitted_data.host$Host, levels = c("One", "Two", "Three"))
# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.host, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host, aes(x = Temperature, y = Fitted, color = Host), size = 1) +
  labs(x = "Temperature", y = "Oviposition in Days") +
  theme_minimal() +
  facet_wrap(~ Species) 
ggsave("tick.host.plot.png", width = 8, height = 5)


# Merge the observed and predicted data for comparison
merged_data.host <- fitted_data.host %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host"))

# Calculate RMSE for each species
mse_by_host <- merged_data.host %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.host <- merged_data.host %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 


# save(mse_by_host,
#      fitted_data.host,
#      merged_data.host,
#      host.fit.mcmc,
#      host.fit,
#      log.lik.host,
#      file="host.fit.RData")
# 

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

  for (j in 1:n.size) {
    mu.size.c[j] ~ dexp(mu.c)       # Size-specific vertex y
    mu.size.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Size-specific vertex x
    mu.size.a[j] ~ dexp(mu.a)    # Size-specific width
    size.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs) {
    mu[i] <- mu.size.a[size[i]] * (mu.size.Topt[size[i]] - temp[i])^2 + mu.size.c[size[i]]
    trait[i] ~ dnorm(mu[i], size.tau[size[i]]) T(0, )  
    log_lik[i] <- logdensity.norm(trait[i], mu[i], size.tau[size[i]])
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
jag.data <- list(trait = tick.abr.new$Trait, N.obs = N.obs, 
                 temp = tick.abr.new$Temp, 
                 n.size = n.size, size = size)
# Fit model
size.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                 model.file = "size_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                 n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
size.fit.mcmc <- as.mcmc(size.fit)
closeAllConnections()
size.fit$BUGSoutput$DIC


# Extract log-likelihood matrix (samples × observations)
log.lik.size <- size.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.size <- waic(log.lik.size)$estimates["waic", "Estimate"]

# Get posterior means for parameters for each size 
summary_fit.size <- summary(size.fit.mcmc)$statistics
a_size_mean <- summary_fit.size[grep("^mu.size.a\\[", rownames(summary_fit.size)), "Mean"]
c_size_mean <- summary_fit.size[grep("^mu.size.c\\[", rownames(summary_fit.size)), "Mean"]
Topt_size_mean <- summary_fit.size[grep("^mu.size.Topt\\[", rownames(summary_fit.size)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp) ,
                max(tick.abr.new$Temp) , by = 0.1)

# Get predicted values for each size category
size1 <- function (x){
  a_size_mean[1] * (x - Topt_size_mean[1])^2 + c_size_mean[1]
}
size2 <- function (x){
  a_size_mean[2] * (x - Topt_size_mean[2])^2 + c_size_mean[2]
}
size3 <- function (x){
  a_size_mean[3] * (x - Topt_size_mean[3])^2 + c_size_mean[3]
}


# Make data frame with temperature, fitted value, and size category
fitted_data.size <- data.frame(Temperature = rep(temp_seq, n.size), 
                               Fitted = c(size1(temp_seq),
                                          size1(temp_seq),
                                          size3(temp_seq)),
                               Size = factor(rep(c("Large", "Medium", "Small"), each = length(temp_seq))))
# Original data with temperature, trait, numeric size, and numeric species
plot_data.size <- data.frame(
  Temperature = tick.abr.new$Temp,
  Trait = tick.abr.new$Trait,
  Size = as.factor(tick.abr.new$Size),
  Species = as.factor(tick.abr.new$Species)
)

# Get unique combinations of species and size for correct fitted lines on plots
species_to_size <- unique(plot_data.size[, c("Species", "Size")])  
species_to_size$Size <- as.factor(species_to_size$Size)
# Merge fitted host curves with species-to-size mapping
fitted_data.size <- data.frame(
  Temperature = rep(temp_seq, n.size),
  Fitted = c(size1(temp_seq),
             size2(temp_seq),
             size3(temp_seq)),
  Size = factor(rep(c("Large", "Medium", "Small"), each = length(temp_seq)))
) %>%
  left_join(species_to_size, by = "Size")

plot_data.size$Size <- factor(plot_data.size$Size, levels = c("Small", "Medium", "Large"))
fitted_data.size$Size <- factor(fitted_data.size$Size, levels = c("Small", "Medium", "Large"))
# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.size, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.size, aes(x = Temperature, y = Fitted, color = Size), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Size",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 


# Merge the observed and predicted data for comparison
merged_data.size <- fitted_data.size %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Size" = "Size"))

# Calculate RMSE for each species
mse_by_size <- merged_data.size %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.size <- merged_data.size %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 
# 
# save(mse_by_size,
#      fitted_data.size,
#      merged_data.size,
#      size.fit.mcmc,
#      size.fit,
#      log.lik.size,
#      file="size.fit.RData")

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
 

  for (j in 1:n.continent) {
    mu.continent.c[j] ~ dexp(mu.c)       # Continent-specific vertex y
    mu.continent.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Continent-specific vertex x
    mu.continent.a[j] ~ dexp(mu.a)    # Continent-specific width
    continent.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs) {
    mu[i] <- mu.continent.a[continent[i]] * (mu.continent.Topt[continent[i]] - temp[i])^2 + mu.continent.c[continent[i]]
    trait[i] ~ dnorm(mu[i], continent.tau[continent[i]]) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], continent.tau[continent[i]])
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
jag.data <- list(trait = tick.abr.new$Trait, N.obs = N.obs, 
                 temp = tick.abr.new$Temp, 
                 n.continent = n.continent, continent = continent)
# Fit model
continent.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                 model.file = "continent_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                 n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
continent.fit.mcmc <- as.mcmc(continent.fit)
closeAllConnections()
continent.fit$BUGSoutput$DIC


# Extract log-likelihood matrix (samples × observations)
log.lik.continent <- continent.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.continent <- waic(log.lik.continent)$estimates["waic", "Estimate"]

# Get posterior means for parameters for different continent groups
summary_fit.continent <- summary(continent.fit.mcmc)$statistics
a_continent_mean <- summary_fit.continent[grep("^mu.continent.a\\[", rownames(summary_fit.continent)), "Mean"]
c_continent_mean <- summary_fit.continent[grep("^mu.continent.c\\[", rownames(summary_fit.continent)), "Mean"]
Topt_continent_mean <- summary_fit.continent[grep("^mu.continent.Topt\\[", rownames(summary_fit.continent)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp) ,
                max(tick.abr.new$Temp) , by = 0.1)

# Get predicted values for each continent group
continent1 <- function (x){
  a_continent_mean[1] * (x - Topt_continent_mean[1])^2 + c_continent_mean[1]
}
continent2 <- function (x){
  a_continent_mean[2] * (x - Topt_continent_mean[2])^2 + c_continent_mean[2]
}
continent3 <- function (x){
  a_continent_mean[3] * (x - Topt_continent_mean[3])^2 + c_continent_mean[3]
}

# Data frame of temperature, fitted value, and continent group
fitted_data.continent <- data.frame(Temperature = rep(temp_seq, n.continent), 
                               Fitted = c(continent1(temp_seq),
                                          continent2(temp_seq),
                                          continent3(temp_seq)),
                               Continent = factor(rep(c("More than two", "One", "Two"), each = length(temp_seq))))
# Original data with temperature, trait, numeric continent, numeric species
plot_data.continent <- data.frame(
  Temperature = tick.abr.new$Temp,
  Trait = tick.abr.new$Trait,
  Continent = as.factor(tick.abr.new$Continent),
  Species = as.factor(tick.abr.new$Species)
)
# Get unique combinations of species and continent for correct plotting
species_to_continent <- unique(plot_data.continent[, c("Species", "Continent")]) 
species_to_continent$Continent <- as.factor(species_to_continent$Continent)
# Merge fitted continent curves with species-to-continent mapping
fitted_data.continent <- data.frame(
  Temperature = rep(temp_seq, n.continent),
  Fitted = c(continent1(temp_seq),
             continent2(temp_seq),
             continent3(temp_seq)),
  Continent = factor(rep(c("More than two", "One", "Two"), each = length(temp_seq)))
) %>%
  left_join(species_to_continent, by = "Continent")

plot_data.continent$Continent <- factor(plot_data.continent$Continent, levels = c("One", "Two", "More than two"))
fitted_data.continent$Continent <- factor(fitted_data.continent$Continent, levels = c("One", "Two", "More than two"))
# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.continent, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.continent, aes(x = Temperature, y = Fitted, color = Continent), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Continent",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 

# Merge the observed and predicted data for comparison
merged_data.continent <- fitted_data.continent %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Continent" = "Continent"))

# Calculate RMSE for each species
mse_by_continent <- merged_data.continent %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.continent <- merged_data.continent %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

# save(mse_by_continent,
#      fitted_data.continent,
#      merged_data.continent,
#      continent.fit.mcmc,
#      continent.fit,
#      log.lik.continent,
#      file="continent.fit.RData")

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
  
  for (j in 1:n.host) {
    mu.host.c[j] ~ dexp(mu.c)       # Host-specific vertex y
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
  }
  for (k in 1:n.climate) {
    mu.climate.Topt[k] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific vertex x
    climate.tau[k] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs) {
    mu[i] <- mu.host.a[host[i]] * (mu.climate.Topt[climate[i]] - temp[i])^2 + mu.host.c[host[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[climate[i]]) T(0, ) 
    log_lik[i] <- logdensity.norm(trait[i], mu[i], climate.tau[climate[i]])
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
jag.data <- list(trait = tick.abr.new$Trait, N.obs = N.obs, 
                 temp = tick.abr.new$Temp, 
                 n.host = n.host, host = host, n.climate = n.climate, climate = climate)
# Fit model
host.climate.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                 model.file = "host.climate_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                 n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.climate.fit.mcmc <- as.mcmc(host.climate.fit) 
closeAllConnections()
host.climate.fit$BUGSoutput$DIC


# Extract log-likelihood matrix (samples × observations)
log.lik.host.climate <- host.climate.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.host.climate <- waic(log.lik.host.climate)$estimates["waic", "Estimate"]

# Get posterior means for parameters for each host group
summary_fit.host.climate <- summary(host.climate.fit.mcmc)$statistics
a_host.climate_mean <- summary_fit.host.climate[grep("^mu.host.a\\[", rownames(summary_fit.host.climate)), "Mean"]
c_host.climate_mean <- summary_fit.host.climate[grep("^mu.host.c\\[", rownames(summary_fit.host.climate)), "Mean"]
Topt_host.climate_mean <- summary_fit.host.climate[grep("^mu.climate.Topt\\[", rownames(summary_fit.host.climate)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp) ,
                max(tick.abr.new$Temp) , by = 0.1)

# Get predicted values for each host group
host1clim1 <- function (x){
  a_host.climate_mean[1] * (x - Topt_host.climate_mean[1])^2 + c_host.climate_mean[1]
}
host2clim1 <- function (x){
  a_host.climate_mean[2] * (x - Topt_host.climate_mean[1])^2 + c_host.climate_mean[2]
}
host3clim1 <- function (x){
  a_host.climate_mean[3] * (x - Topt_host.climate_mean[1])^2 + c_host.climate_mean[3]
}
host1clim2 <- function (x){
  a_host.climate_mean[1] * (x - Topt_host.climate_mean[2])^2 + c_host.climate_mean[1]
}
host2clim2 <- function (x){
  a_host.climate_mean[2] * (x - Topt_host.climate_mean[2])^2 + c_host.climate_mean[2]
}
host3clim2 <- function (x){
  a_host.climate_mean[3] * (x - Topt_host.climate_mean[2])^2 + c_host.climate_mean[3]
}

host1clim3 <- function (x){
  a_host.climate_mean[1] * (x - Topt_host.climate_mean[3])^2 + c_host.climate_mean[1]
}
host2clim3 <- function (x){
  a_host.climate_mean[2] * (x - Topt_host.climate_mean[3])^2 + c_host.climate_mean[2]
}
host3clim3 <- function (x){
  a_host.climate_mean[3] * (x - Topt_host.climate_mean[3])^2 + c_host.climate_mean[3]
}

host1clim4 <- function (x){
  a_host.climate_mean[1] * (x - Topt_host.climate_mean[4])^2 + c_host.climate_mean[1]
}
host2clim4 <- function (x){
  a_host.climate_mean[2] * (x - Topt_host.climate_mean[4])^2 + c_host.climate_mean[2]
}
host3clim4 <- function (x){
  a_host.climate_mean[3] * (x - Topt_host.climate_mean[4])^2 + c_host.climate_mean[3]
}



# Data frame with temperature, fitted value, and host group
fitted_data.host.climate <- data.frame(Temperature = rep(temp_seq, n.host * n.climate), 
                               Fitted = c(host1clim1(temp_seq),
                                          host2clim1(temp_seq),
                                          host3clim1(temp_seq),
                                          host1clim2(temp_seq),
                                          host2clim2(temp_seq),
                                          host3clim2(temp_seq),
                                          host1clim3(temp_seq),
                                          host2clim3(temp_seq),
                                          host3clim3(temp_seq),
                                          host1clim4(temp_seq),
                                          host2clim4(temp_seq),
                                          host3clim4(temp_seq)),
                               Host = factor(rep(rep(c("One", "Three", "Two"), each = length(temp_seq)), n.climate)),
                               Climate = factor(rep(c("Mixed", "Subtropical", "Temperate", "Tropical"),
                                                    each = length(temp_seq) * n.size)))
# Original data with temperature, trait, numeric host group, numeric species
plot_data.host.climate <- data.frame(
  Temperature = tick.abr.new$Temp,
  Trait = tick.abr.new$Trait,
  Host = as.factor(tick.abr.new$Host),
  Species = as.factor(tick.abr.new$Species),
  Climate = as.factor(tick.abr.new$Climate)
)

# Get unique combinations of species and host for correct plotting
species_to_host.climate <- unique(plot_data.host.climate[, c("Species", "Host", "Climate")])  
species_to_host.climate$Host <- as.factor(species_to_host.climate$Host)
species_to_host.climate$Climate <- as.factor(species_to_host.climate$Climate)

# Merge fitted host curves with species-to-host mapping
fitted_data.host.climate <- data.frame(Temperature = rep(temp_seq, n.host * n.climate), 
                                       Fitted = c(host1clim1(temp_seq),
                                                  host2clim1(temp_seq),
                                                  host3clim1(temp_seq),
                                                  host1clim2(temp_seq),
                                                  host2clim2(temp_seq),
                                                  host3clim2(temp_seq),
                                                  host1clim3(temp_seq),
                                                  host2clim3(temp_seq),
                                                  host3clim3(temp_seq),
                                                  host1clim4(temp_seq),
                                                  host2clim4(temp_seq),
                                                  host3clim4(temp_seq)),
                                       Host = factor(rep(rep(c("One", "Three", "Two"), each = length(temp_seq)), n.climate)),
                                       Climate = factor(rep(c("Mixed", "Subtropical", "Temperate", "Tropical"),
                                                            each = length(temp_seq) * n.size))) %>%
  left_join(species_to_host.climate, by = c("Host", "Climate"))|> 
  filter(!is.na(Species)) 

plot_data.host.climate$Host <- factor(plot_data.host.climate$Host, levels = c("One", "Two", "Three"))
fitted_data.host.climate$Host <- factor(fitted_data.host.climate$Host, levels = c("One", "Two", "Three"))

# Combine Climate and Size for coloring
fitted_data.host.climate$Combo <- with(fitted_data.host.climate, paste(Host, Climate, sep = ", "))
fitted_data.host.climate$Combo <- factor(fitted_data.host.climate$Combo, levels = c("One, Tropical", 
                                                                                        "Two, Subtropical", "Three, Mixed",
                                                                                        "Three, Subtropical", "Three, Temperate",
                                                                                    "Three, Tropical"))

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.host.climate, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.climate, aes(x = Temperature, y = Fitted, color = Combo), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Host and Climate",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 

# Merge the observed and predicted data for comparison
merged_data.host.climate <- fitted_data.host.climate %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host", "Climate" = "Climate"))

# Calculate RMSE for each species
mse_by_host.climate <- merged_data.host.climate %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.host.climate <- merged_data.host.climate %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 


# save(mse_by_host.climate,
#      fitted_data.host.climate,
#      merged_data.host.climate,
#      host.climate.fit.mcmc,
#      host.climate.fit,
#      log.lik.host.climate,
#      file="host.climate.fit.RData")

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
 

  for (j in 1:n.host) {
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
    host.tau[j] ~ dexp(tau1)
  }
  for (k in 1:n.genus) {
    mu.genus.c[k] ~ dexp(mu.c)       # Genus-specific vertex y
    mu.genus.Topt[k] ~ dnorm(mu.Topt, tau.Topt)    # Genus-specific vertex x
  }
  # Likelihood
  for (i in 1:N.obs) {
    mu[i] <- mu.host.a[host[i]] * (mu.genus.Topt[genus[i]] - temp[i])^2 + mu.genus.c[genus[i]]
    trait[i] ~ dnorm(mu[i], host.tau[host[i]]) T(0, )  
    log_lik[i] <- logdensity.norm(trait[i], mu[i], host.tau[host[i]])
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
jag.data <- list(trait = tick.abr.new$Trait, N.obs = N.obs, 
                 temp = tick.abr.new$Temp, 
                 n.host = n.host, host = host, n.genus = n.genus, genus = genus)
# Fit model
host.genus.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                         model.file = "host.genus_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                         n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.genus.fit.mcmc <- as.mcmc(host.genus.fit) 
closeAllConnections()
host.genus.fit$BUGSoutput$DIC


# Extract log-likelihood matrix (samples × observations)
log.lik.host.genus <- host.genus.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.host.genus <- waic(log.lik.host.genus)$estimates["waic", "Estimate"]

# Get posterior means for parameters for each host group
summary_fit.host.genus <- summary(host.genus.fit.mcmc)$statistics
a_host.genus_mean <- summary_fit.host.genus[grep("^mu.host.a\\[", rownames(summary_fit.host.genus)), "Mean"]
c_host.genus_mean <- summary_fit.host.genus[grep("^mu.genus.c\\[", rownames(summary_fit.host.genus)), "Mean"]
Topt_host.genus_mean <- summary_fit.host.genus[grep("^mu.genus.Topt\\[", rownames(summary_fit.host.genus)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp) ,
                max(tick.abr.new$Temp) , by = 0.1)

# Get predicted values for each host group
host1genus1 <- function (x){
  a_host.genus_mean[1] * (x - Topt_host.genus_mean[1])^2 + c_host.genus_mean[1]
}
host2genus1 <- function (x){
  a_host.genus_mean[2] * (x - Topt_host.genus_mean[1])^2 + c_host.genus_mean[1]
}
host3genus1 <- function (x){
  a_host.genus_mean[3] * (x - Topt_host.genus_mean[1])^2 + c_host.genus_mean[1]
}
host1genus2 <- function (x){
  a_host.genus_mean[1] * (x - Topt_host.genus_mean[2])^2 + c_host.genus_mean[2]
}
host2genus2 <- function (x){
  a_host.genus_mean[2] * (x - Topt_host.genus_mean[2])^2 + c_host.genus_mean[2]
}
host3genus2 <- function (x){
  a_host.genus_mean[3] * (x - Topt_host.genus_mean[2])^2 + c_host.genus_mean[2]
}

host1genus3 <- function (x){
  a_host.genus_mean[1] * (x - Topt_host.genus_mean[3])^2 + c_host.genus_mean[3]
}
host2genus3 <- function (x){
  a_host.genus_mean[2] * (x - Topt_host.genus_mean[3])^2 + c_host.genus_mean[3]
}
host3genus3 <- function (x){
  a_host.genus_mean[3] * (x - Topt_host.genus_mean[3])^2 + c_host.genus_mean[3]
}

host1genus4 <- function (x){
  a_host.genus_mean[1] * (x - Topt_host.genus_mean[4])^2 + c_host.genus_mean[4]
}
host2genus4 <- function (x){
  a_host.genus_mean[2] * (x - Topt_host.genus_mean[4])^2 + c_host.genus_mean[4]
}
host3genus4 <- function (x){
  a_host.genus_mean[3] * (x - Topt_host.genus_mean[4])^2 + c_host.genus_mean[4]
}



# Data frame with temperature, fitted value, and host group
fitted_data.host.genus <- data.frame(Temperature = rep(temp_seq, n.host * n.genus), 
                                       Fitted = c(host1genus1(temp_seq),
                                                  host2genus1(temp_seq),
                                                  host3genus1(temp_seq),
                                                  host1genus2(temp_seq),
                                                  host2genus2(temp_seq),
                                                  host3genus2(temp_seq),
                                                  host1genus3(temp_seq),
                                                  host2genus3(temp_seq),
                                                  host3genus3(temp_seq),
                                                  host1genus4(temp_seq),
                                                  host2genus4(temp_seq),
                                                  host3genus4(temp_seq)),
                                       Host = factor(rep(rep(c("One", "Three", "Two"), each = length(temp_seq)), n.genus)),
                                       Genus = factor(rep(c("Amblyomma", "Dermacentor", 
                                                            "Haemaphysalis", "Hyalomma"),
                                                            each = length(temp_seq) * n.host)))
# Original data with temperature, trait, numeric host group, numeric species
plot_data.host.genus <- data.frame(
  Temperature = tick.abr.new$Temp,
  Trait = tick.abr.new$Trait,
  Host = as.factor(tick.abr.new$Host),
  Species = as.factor(tick.abr.new$Species),
  Genus = as.factor(tick.abr.new$Genus)
)

# Get unique combinations of species and host for correct plotting
species_to_host.genus <- unique(plot_data.host.genus[, c("Species", "Host", "Genus")])  
species_to_host.genus$Host <- as.factor(species_to_host.genus$Host)
species_to_host.genus$Genus <- as.factor(species_to_host.genus$Genus)

# Merge fitted host curves with species-to-host mapping
fitted_data.host.genus <- data.frame(Temperature = rep(temp_seq, n.host * n.genus), 
                                       Fitted = c(host1genus1(temp_seq),
                                                  host2genus1(temp_seq),
                                                  host3genus1(temp_seq),
                                                  host1genus2(temp_seq),
                                                  host2genus2(temp_seq),
                                                  host3genus2(temp_seq),
                                                  host1genus3(temp_seq),
                                                  host2genus3(temp_seq),
                                                  host3genus3(temp_seq),
                                                  host1genus4(temp_seq),
                                                  host2genus4(temp_seq),
                                                  host3genus4(temp_seq)),
                                       Host = factor(rep(rep(c("One", "Three", "Two"), each = length(temp_seq)), n.genus)),
                                       Genus = factor(rep(c("Amblyomma", "Dermacentor", 
                                                            "Haemaphysalis", "Hyalomma"),
                                                            each = length(temp_seq) * n.host))) %>%
  left_join(species_to_host.genus, by = c("Host", "Genus"))|> 
  filter(!is.na(Species)) 

plot_data.host.genus$Host <- factor(plot_data.host.genus$Host, levels = c("One", "Two", "Three"))
fitted_data.host.genus$Host <- factor(fitted_data.host.genus$Host, levels = c("One", "Two", "Three"))

# Combine genus and Size for coloring
fitted_data.host.genus$Combo <- with(fitted_data.host.genus, paste(Host, Genus, sep = ", "))
fitted_data.host.genus$Combo <- factor(fitted_data.host.genus$Combo, levels = c("One, Dermacentor", 
                                                                                    "Two, Hyalomma", "Three, Amblyomma",
                                                                                    "Three, Dermacentor", "Three, Haemaphysalis",
                                                                                    "Three, Hyalomma"))

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.host.genus, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.genus, aes(x = Temperature, y = Fitted, color = Combo), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Host and Genus",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 

# Merge the observed and predicted data for comparison
merged_data.host.genus <- fitted_data.host.genus %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host", "Genus" = "Genus"))

# Calculate RMSE for each species
mse_by_host.genus <- merged_data.host.genus %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.host.genus <- merged_data.host.genus %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE)))

# 
# save(mse_by_host.genus,
#      fitted_data.host.genus,
#      merged_data.host.genus,
#      host.genus.fit.mcmc,
#      host.genus.fit,
#      log.lik.host.genus,
#      file="host.genus.fit.RData")

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

  for (j in 1:n.genus) {
    mu.genus.c[j] ~ dexp(mu.c)       # genus-specific vertex y
    mu.genus.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # genus-specific vertex x
  }
  for (k in 1:n.climate) {
    mu.climate.a[k] ~ dexp(mu.a)    # Cimate-specific width
    climate.tau[k] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs) {
    mu[i] <- mu.climate.a[climate[i]] * (mu.genus.Topt[genus[i]] - temp[i])^2 + mu.genus.c[genus[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[climate[i]]) T(0, ) 
    log_lik[i] <- logdensity.norm(trait[i], mu[i], climate.tau[climate[i]])
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
jag.data <- list(trait = tick.abr.new$Trait, N.obs = N.obs, 
                 temp = tick.abr.new$Temp, 
                 n.genus = n.genus, genus = genus, n.climate = n.climate, climate = climate)
# Fit model
genus.climate.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                         model.file = "genus.climate_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                         n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
genus.climate.fit.mcmc <- as.mcmc(genus.climate.fit) 
closeAllConnections()
genus.climate.fit$BUGSoutput$DIC


# Extract log-likelihood matrix (samples × observations)
log.lik.genus.climate <- genus.climate.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.genus.climate <- waic(log.lik.genus.climate)$estimates["waic", "Estimate"]


# Get posterior means for parameters for each genus group
summary_fit.genus.climate <- summary(genus.climate.fit.mcmc)$statistics
a_genus.climate_mean <- summary_fit.genus.climate[grep("^mu.climate.a\\[", rownames(summary_fit.genus.climate)), "Mean"]
c_genus.climate_mean <- summary_fit.genus.climate[grep("^mu.genus.c\\[", rownames(summary_fit.genus.climate)), "Mean"]
Topt_genus.climate_mean <- summary_fit.genus.climate[grep("^mu.genus.Topt\\[", rownames(summary_fit.genus.climate)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp) ,
                max(tick.abr.new$Temp) , by = 0.1)

# Get predicted values for each genus group
genus1clim1 <- function (x){
  a_genus.climate_mean[1] * (x - Topt_genus.climate_mean[1])^2 + c_genus.climate_mean[1]
}
genus2clim1 <- function (x){
  a_genus.climate_mean[1] * (x - Topt_genus.climate_mean[2])^2 + c_genus.climate_mean[2]
}
genus3clim1 <- function (x){
  a_genus.climate_mean[1] * (x - Topt_genus.climate_mean[3])^2 + c_genus.climate_mean[3]
}
genus4clim1 <- function (x){
  a_genus.climate_mean[1] * (x - Topt_genus.climate_mean[4])^2 + c_genus.climate_mean[4]
}
genus1clim2 <- function (x){
  a_genus.climate_mean[2] * (x - Topt_genus.climate_mean[1])^2 + c_genus.climate_mean[1]
}
genus2clim2 <- function (x){
  a_genus.climate_mean[2] * (x - Topt_genus.climate_mean[2])^2 + c_genus.climate_mean[2]
}
genus3clim2 <- function (x){
  a_genus.climate_mean[2] * (x - Topt_genus.climate_mean[3])^2 + c_genus.climate_mean[3]
}
genus4clim2 <- function (x){
  a_genus.climate_mean[2] * (x - Topt_genus.climate_mean[4])^2 + c_genus.climate_mean[4]
}
genus1clim3 <- function (x){
  a_genus.climate_mean[3] * (x - Topt_genus.climate_mean[1])^2 + c_genus.climate_mean[1]
}
genus2clim3 <- function (x){
  a_genus.climate_mean[3] * (x - Topt_genus.climate_mean[2])^2 + c_genus.climate_mean[2]
}
genus3clim3 <- function (x){
  a_genus.climate_mean[3] * (x - Topt_genus.climate_mean[3])^2 + c_genus.climate_mean[3]
}
genus4clim3 <- function (x){
  a_genus.climate_mean[3] * (x - Topt_genus.climate_mean[4])^2 + c_genus.climate_mean[4]
}
genus1clim4 <- function (x){
  a_genus.climate_mean[4] * (x - Topt_genus.climate_mean[1])^2 + c_genus.climate_mean[1]
}
genus2clim4 <- function (x){
  a_genus.climate_mean[4] * (x - Topt_genus.climate_mean[2])^2 + c_genus.climate_mean[2]
}
genus3clim4 <- function (x){
  a_genus.climate_mean[4] * (x - Topt_genus.climate_mean[3])^2 + c_genus.climate_mean[3]
}
genus4clim4 <- function (x){
  a_genus.climate_mean[4] * (x - Topt_genus.climate_mean[4])^2 + c_genus.climate_mean[4]
}



# Data frame with temperature, fitted value, and genus group
fitted_data.genus.climate <- data.frame(Temperature = rep(temp_seq, n.genus * n.climate), 
                                       Fitted = c(genus1clim1(temp_seq),
                                                  genus2clim1(temp_seq),
                                                  genus3clim1(temp_seq),
                                                  genus4clim1(temp_seq),
                                                  genus1clim2(temp_seq),
                                                  genus2clim2(temp_seq),
                                                  genus3clim2(temp_seq),
                                                  genus4clim2(temp_seq),
                                                  genus1clim3(temp_seq),
                                                  genus2clim3(temp_seq),
                                                  genus3clim3(temp_seq),
                                                  genus4clim3(temp_seq),
                                                  genus1clim4(temp_seq),
                                                  genus2clim4(temp_seq),
                                                  genus3clim4(temp_seq),
                                                  genus4clim4(temp_seq)),
                                       Genus = factor(rep(rep(c("Amblyomma", "Dermacentor", 
                                                                "Haemaphysalis", "Hyalomma"), each = length(temp_seq)), n.climate)),
                                       Climate = factor(rep(c("Mixed", "Subtropical", "Temperate", "Tropical"),
                                                            each = length(temp_seq) * n.genus)))

# Original data with temperature, trait, numeric genus group, numeric species
plot_data.genus.climate <- data.frame(
  Temperature = tick.abr.new$Temp,
  Trait = tick.abr.new$Trait,
  Genus = as.factor(tick.abr.new$Genus),
  Species = as.factor(tick.abr.new$Species),
  Climate = as.factor(tick.abr.new$Climate)
)

# Get unique combinations of species and genus for correct plotting
species_to_genus.climate <- unique(plot_data.genus.climate[, c("Species", "Genus", "Climate")])  
species_to_genus.climate$Genus <- as.factor(species_to_genus.climate$Genus)
species_to_genus.climate$Climate <- as.factor(species_to_genus.climate$Climate)

# Merge fitted genus curves with species-to-genus mapping
fitted_data.genus.climate <- data.frame(Temperature = rep(temp_seq, n.genus * n.climate), 
                                       Fitted = c(genus1clim1(temp_seq),
                                                  genus2clim1(temp_seq),
                                                  genus3clim1(temp_seq),
                                                  genus4clim1(temp_seq),
                                                  genus1clim2(temp_seq),
                                                  genus2clim2(temp_seq),
                                                  genus3clim2(temp_seq),
                                                  genus4clim2(temp_seq),
                                                  genus1clim3(temp_seq),
                                                  genus2clim3(temp_seq),
                                                  genus3clim3(temp_seq),
                                                  genus4clim3(temp_seq),
                                                  genus1clim4(temp_seq),
                                                  genus2clim4(temp_seq),
                                                  genus3clim4(temp_seq),
                                                  genus4clim4(temp_seq)),
                                       Genus = factor(rep(rep(c("Amblyomma", "Dermacentor", 
                                                                "Haemaphysalis", "Hyalomma"), each = length(temp_seq)), n.climate)),
                                       Climate = factor(rep(c("Mixed", "Subtropical", "Temperate", "Tropical"),
                                                            each = length(temp_seq) * n.genus))) %>%
  left_join(species_to_genus.climate, by = c("Genus", "Climate"))|> 
  filter(!is.na(Species)) 


# Combine Climate and Size for coloring
fitted_data.genus.climate$Combo <- with(fitted_data.genus.climate, paste(Genus, Climate, sep = ", "))


# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.genus.climate, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.genus.climate, aes(x = Temperature, y = Fitted, color = Combo), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Genus and Climate",
       x = "Temperature", y = "Ovipostion in Days") +
  theme_minimal() +
  facet_wrap(~ Species) 
ggsave("tick.genus.climate.plot.png", width = 8, height = 5)

# Merge the observed and predicted data for comparison
merged_data.genus.climate <- fitted_data.genus.climate %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", 
                         "Genus" = "Genus", "Climate" = "Climate"))

# Calculate RMSE for each species
mse_by_genus.climate <- merged_data.genus.climate %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.genus.climate <- merged_data.genus.climate %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

# 
# save(mse_by_genus.climate,
#      fitted_data.genus.climate,
#      merged_data.genus.climate,
#      genus.climate.fit.mcmc,
#      genus.climate.fit,
#      log.lik.genus.climate,
#      file="genus.climate.fit.RData")

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

  for (j in 1:n.host) {
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
  }
  for (k in 1:n.climate) {
    mu.climate.Topt[k] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific vertex x
    climate.tau[k] ~ dexp(tau1)
  }
  for (s in 1:n.size) {
    mu.size.c[s] ~ dexp(mu.c)       # Size-specific vertex y
  }
  # Likelihood
  for (i in 1:N.obs) {
    mu[i] <- mu.host.a[host[i]] * (mu.climate.Topt[climate[i]] - temp[i])^2 + mu.size.c[size[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[climate[i]]) T(0, ) 
    log_lik[i] <- logdensity.norm(trait[i], mu[i], climate.tau[climate[i]])
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
jag.data <- list(trait = tick.abr.new$Trait, N.obs = N.obs, 
                 temp = tick.abr.new$Temp, 
                 n.host = n.host, host = host, n.climate = n.climate, climate = climate,
                 n.size = n.size, size = size)
# Fit model
h.c.s.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                         model.file = "host.climate.size_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                         n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
h.c.s.fit.mcmc <- as.mcmc(h.c.s.fit) 
closeAllConnections()
h.c.s.fit$BUGSoutput$DIC


# Extract log-likelihood matrix (samples × observations)
log.lik.h.c.s <- h.c.s.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.h.c.s <- waic(log.lik.h.c.s)$estimates["waic", "Estimate"]

# Get posterior means for parameters for each host group
summary_fit.h.c.s <- summary(h.c.s.fit.mcmc)$statistics
a_h.c.s_mean <- summary_fit.h.c.s[grep("^mu.host.a\\[", rownames(summary_fit.h.c.s)), "Mean"]
c_h.c.s_mean <- summary_fit.h.c.s[grep("^mu.size.c\\[", rownames(summary_fit.h.c.s)), "Mean"]
Topt_h.c.s_mean <- summary_fit.h.c.s[grep("^mu.climate.Topt\\[", rownames(summary_fit.h.c.s)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp) ,
                max(tick.abr.new$Temp) , by = 0.1)

# Get predicted values for each  group
for (host in 1:3) {
  for (clim in 1:4) {
    for (size in 1:3) {
      func_name <- paste0("host", host, "clim", clim, "size", size)
      assign(func_name, 
             local({
               h <- host
               c <- clim
               s <- size
               force(h); force(c); force(s)  # Ensure values are captured
               
               function(x) {
                 a_h.c.s_mean[h] * (x - Topt_h.c.s_mean[c])^2 + c_h.c.s_mean[s]
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
      func_name <- paste0("host", host, "clim", clim, "size", size)
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
fitted_data.h.c.s <- do.call(rbind, fitted_list)

# Original data with temperature, trait, factor host group, species, climate, size
plot_data.h.c.s <- data.frame(
  Temperature = tick.abr.new$Temp,
  Trait = tick.abr.new$Trait,
  Host = as.factor(tick.abr.new$Host),
  Species = as.factor(tick.abr.new$Species),
  Climate = as.factor(tick.abr.new$Climate),
  Size = as.factor(tick.abr.new$Size)
)

# Get unique combinations of species and host for correct plotting
species_to_h.c.s <- unique(plot_data.h.c.s[, c("Species", "Host", "Climate", "Size")])  
species_to_h.c.s$Host <- as.factor(species_to_h.c.s$Host)
species_to_h.c.s$Climate <- as.factor(species_to_h.c.s$Climate)
species_to_h.c.s$Size <- as.factor(species_to_h.c.s$Size)

# Merge fitted host curves with species-to-host mapping
fitted_data.h.c.s <- fitted_data.h.c.s |>  
  left_join(species_to_h.c.s, by = c("Host", "Climate", "Size"))|> 
  filter(!is.na(Species)) 

# Order factor levels
plot_data.h.c.s$Host <- factor(plot_data.h.c.s$Host, levels = c("One", "Two", "Three"))
fitted_data.h.c.s$Host <- factor(fitted_data.h.c.s$Host, levels = c("One", "Two", "Three"))
plot_data.h.c.s$Size <- factor(plot_data.h.c.s$Size, levels = c("Small", "Medium", "Large"))
fitted_data.h.c.s$Size <- factor(fitted_data.h.c.s$Size, levels = c("Small", "Medium", "Large"))

# Combine Host, Climate, and Size for coloring
fitted_data.h.c.s$Combo <- with(fitted_data.h.c.s, paste(Host, Climate, Size, sep = ", "))
fitted_data.h.c.s$Combo <- factor(fitted_data.h.c.s$Combo, levels = c("One, Tropical, Medium",
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
  geom_point(data = plot_data.h.c.s, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.h.c.s, aes(x = Temperature, y = Fitted, color = Combo), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Host, Climate, and Size",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 

# Merge the observed and predicted data for comparison
merged_data.h.c.s <- fitted_data.h.c.s %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", 
                         "Host" = "Host", "Climate" = "Climate",
                         "Size" = "Size"))

# Calculate RMSE for each species
mse_by_h.c.s <- merged_data.h.c.s %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.h.c.s <- merged_data.h.c.s %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 


# save(mse_by_h.c.s,
#      fitted_data.h.c.s,
#      merged_data.h.c.s,
#      h.c.s.fit.mcmc,
#      h.c.s.fit,
#      log.lik.h.c.s,
#      file="host.climate.size.fit.RData")

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

  for (j in 1:n.host) {
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
  }
  for (k in 1:n.climate) {
    mu.climate.Topt[k] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific vertex x
    climate.tau[k] ~ dexp(tau1)
  }
  for (s in 1:n.genus) {
    mu.genus.c[s] ~ dexp(mu.c)       # genus-specific vertex y
  }
  # Likelihood
  for (i in 1:N.obs) {
    mu[i] <- mu.host.a[host[i]] * (mu.climate.Topt[climate[i]] - temp[i])^2 + mu.genus.c[genus[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[climate[i]]) T(0, ) 
    log_lik[i] <- logdensity.norm(trait[i], mu[i], climate.tau[climate[i]])
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
jag.data <- list(trait = tick.abr.new$Trait, N.obs = N.obs, 
                 temp = tick.abr.new$Temp, 
                 n.host = n.host, host = host, n.climate = n.climate, climate = climate,
                 n.genus = n.genus, genus = genus)
# Fit model
h.c.g.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                  model.file = "host.climate.genus_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                  n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
h.c.g.fit.mcmc <- as.mcmc(h.c.g.fit) 
closeAllConnections()
h.c.g.fit$BUGSoutput$DIC


# Extract log-likelihood matrix (samples × observations)
log.lik.h.c.g <- h.c.g.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.h.c.g <- waic(log.lik.h.c.g)$estimates["waic", "Estimate"]

# Get posterior means for parameters for each host group
summary_fit.h.c.g <- summary(h.c.g.fit.mcmc)$statistics
a_h.c.g_mean <- summary_fit.h.c.g[grep("^mu.host.a\\[", rownames(summary_fit.h.c.g)), "Mean"]
c_h.c.g_mean <- summary_fit.h.c.g[grep("^mu.genus.c\\[", rownames(summary_fit.h.c.g)), "Mean"]
Topt_h.c.g_mean <- summary_fit.h.c.g[grep("^mu.climate.Topt\\[", rownames(summary_fit.h.c.g)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp) ,
                max(tick.abr.new$Temp) , by = 0.1)

# Get predicted values for each  group
for (host in 1:3) {
  for (clim in 1:4) {
    for (genus in 1:4) {
      func_name <- paste0("host", host, "clim", clim, "genus", genus)
      assign(func_name, 
             local({
               h <- host
               c <- clim
               g <- genus
               force(h); force(c); force(g)  # Ensure values are captured
               
               function(x) {
                 a_h.c.g_mean[h] * (x - Topt_h.c.g_mean[c])^2 + c_h.c.g_mean[g]
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
      func_name <- paste0("host", host, "clim", clim, "genus", genus)
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
fitted_data.h.c.g <- do.call(rbind, fitted_list)

# Original data with temperature, trait, factor host group, species, climate, size
plot_data.h.c.g <- data.frame(
  Temperature = tick.abr.new$Temp,
  Trait = tick.abr.new$Trait,
  Host = as.factor(tick.abr.new$Host),
  Species = as.factor(tick.abr.new$Species),
  Climate = as.factor(tick.abr.new$Climate),
  Genus = as.factor(tick.abr.new$Genus)
)

# Get unique combinations of species and host for correct plotting
species_to_h.c.g <- unique(plot_data.h.c.g[, c("Species", "Host", "Climate", "Genus")])  
species_to_h.c.g$Host <- as.factor(species_to_h.c.g$Host)
species_to_h.c.g$Climate <- as.factor(species_to_h.c.g$Climate)
species_to_h.c.g$Genus <- as.factor(species_to_h.c.g$Genus)

# Merge fitted host curves with species-to-host mapping
fitted_data.h.c.g <- fitted_data.h.c.g |>  
  left_join(species_to_h.c.g, by = c("Host", "Climate", "Genus"))|> 
  filter(!is.na(Species)) 

# Order factor levels
plot_data.h.c.g$Host <- factor(plot_data.h.c.s$Host, levels = c("One", "Two", "Three"))
fitted_data.h.c.g$Host <- factor(fitted_data.h.c.s$Host, levels = c("One", "Two", "Three"))

# Combine Host, Climate, and Size for coloring
fitted_data.h.c.g$Combo <- with(fitted_data.h.c.g, paste(Host, Climate, Genus, sep = ", "))

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.h.c.g, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.h.c.g, aes(x = Temperature, y = Fitted, color = Combo), size = 1) +
  labs(x = "Temperature", y = "Oviposition in Days") +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~ Species) 
ggsave("tick.h.c.g.plot.png", width = 6.5, height = 5)

# Merge the observed and predicted data for comparison
merged_data.h.c.g <- fitted_data.h.c.g %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", 
                         "Host" = "Host", "Climate" = "Climate",
                         "Genus" = "Genus"))

# Calculate RMSE for each species
mse_by_h.c.g <- merged_data.h.c.g %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.h.c.g <- merged_data.h.c.g %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 


# save(mse_by_h.c.g,
#      fitted_data.h.c.g,
#      merged_data.h.c.g,
#      h.c.g.fit.mcmc,
#      h.c.g.fit,
#      log.lik.h.c.g,
#      file="host.climate.genus.fit.RData")

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
table.tick <- data.frame(
  Method = factor(rep(c("Full Data", "By Species Individual", "By Species Hierarchical", 
                        "By Species Hierarchical Different Variances", "By Species Different Variances, Same Parameters", 
                        "By Genus", 
                        "By Climate", "By Host", "By Size", "By Continents",
                        "By Host and Climate", "By Host and Genus", 
                        "By Genus and Climate", "By Host, Climate, and Size",
                        "By Host, Climate, and Genus"), each = 9)),
  DIC = rep(c(full.data.fit$BUGSoutput$DIC, DIC.species.individual, species.fit$BUGSoutput$DIC,
              species.tau.fit$BUGSoutput$DIC, species.tau.param.fit$BUGSoutput$DIC,
              genus.fit$BUGSoutput$DIC,
              climate.fit$BUGSoutput$DIC, host.fit$BUGSoutput$DIC, size.fit$BUGSoutput$DIC, 
              continent.fit$BUGSoutput$DIC,
              host.climate.fit$BUGSoutput$DIC, host.genus.fit$BUGSoutput$DIC,
              genus.climate.fit$BUGSoutput$DIC,
              h.c.s.fit$BUGSoutput$DIC, h.c.g.fit$BUGSoutput$DIC), each = 9),
  wAIC = rep(c(waic.full, waic.species.ind, waic.species,
           waic.species.tau, waic.species.tau.param,
           waic.genus, waic.climate, waic.host,
           waic.size, waic.continent,
           waic.host.climate, waic.host.genus,
           waic.genus.climate, 
           waic.h.c.s, waic.h.c.g), each = 9),
  Species = factor(rep(c("Amblyomma lepidum", "Dermacentor andersoni",
                         "Dermacentor nitens", "Haemaphysalis leporispalustris",
                         "Hyalomma aegyptium", "Hyalomma dromedarii",
                         "Hyalomma impeltatum", "Hyalomma lusitanicum",
                         "Hyalomma schulzei"), times = 15)),
  # MSE = c(mse_by_full_data$RMSE, as.numeric(rmse.a.lepidum$MSE), as.numeric(rmse.d.andersoni$MSE),
  #          as.numeric(rmse.d.nitens$MSE), 
  #          as.numeric(rmse.h.leporispalustris$MSE), as.numeric(rmse.h.aegyptium$MSE), 
  #          as.numeric(rmse.h.dromedarii$MSE), as.numeric(rmse.h.impeltatum$MSE), 
  #          as.numeric(rmse.h.lusitanicum$MSE), as.numeric(rmse.h.schulzei$MSE), 
  #          mse_by_species$MSE, mse_by_species.tau$MSE, 
  #          mse_by_species.tau.param$MSE, mse_by_genus$MSE,
  #          mse_by_climate$MSE, mse_by_host$MSE, mse_by_size$MSE,
  #          mse_by_continent$MSE, mse_by_host.climate$MSE,
  #          mse_by_host.genus$MSE, mse_by_genus.climate$MSE,
  #          mse_by_h.c.s$MSE, mse_by_h.c.g$MSE),
  Total.RMSE = as.numeric(rep(c(overall_rmse, overall_rmse.species.ind,
                       overall_rmse.species, overall_rmse.species.tau, 
                       overall_rmse.species.tau.param, overall_rmse.genus, 
                       overall_rmse.climate, overall_rmse.host, overall_rmse.size, 
                       overall_rmse.continent, overall_rmse.host.climate, 
                       overall_rmse.host.genus,
                       overall_rmse.genus.climate,
                       overall_rmse.h.c.s, overall_rmse.h.c.g), each = 9))
)

table.tick.avg <- data.frame(
  Method = factor(c("Full Data", "By Species Individual", "By Species Hierarchical", 
                    "By Species Hierarchical Different Variances", 
                    "By Genus", 
             "By Climate", "By Host", "By Size", "By Continents",
             "By Host and Climate", "By Host and Genus", 
             "By Genus and Climate", "By Host, Climate, and Size",
             "By Host, Climate, and Genus")),
  DIC = c(full.data.fit$BUGSoutput$DIC, DIC.species.individual,
          species.fit$BUGSoutput$DIC, species.tau.fit$BUGSoutput$DIC, 
          genus.fit$BUGSoutput$DIC,
          climate.fit$BUGSoutput$DIC, host.fit$BUGSoutput$DIC, size.fit$BUGSoutput$DIC, 
          continent.fit$BUGSoutput$DIC,
          host.climate.fit$BUGSoutput$DIC, host.genus.fit$BUGSoutput$DIC,
          genus.climate.fit$BUGSoutput$DIC,
          h.c.s.fit$BUGSoutput$DIC,
          h.c.g.fit$BUGSoutput$DIC),
  wAIC = c(waic.full, waic.species.ind, waic.species,
           waic.species.tau, 
           waic.genus, waic.climate, waic.host,
           waic.size, waic.continent,
           waic.host.climate, waic.host.genus,
           waic.genus.climate, 
           waic.h.c.s, waic.h.c.g),
  Total.RMSE = round(as.numeric(c(overall_rmse, overall_rmse.species.ind,
                 overall_rmse.species, overall_rmse.species.tau, 
                overall_rmse.genus, 
                 overall_rmse.climate, overall_rmse.host, overall_rmse.size, 
                 overall_rmse.continent, overall_rmse.host.climate,
                 overall_rmse.host.genus,
                 overall_rmse.genus.climate,
                 overall_rmse.h.c.s, overall_rmse.h.c.g)), 3)
)

table.tick.avg$IC <- table.tick.avg$DIC + table.tick.avg$wAIC

table.tick.avg$Method <- reorder(table.tick.avg$Method, -table.tick.avg$Total.RMSE)
ggplot(table.tick.avg, aes(x = Method, y = Total.RMSE, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  #labs(title = "DIC by Model", x = "Model",
  #    y = "Deviance Information Criterion (DIC)") +
  labs(title = "RMSE by Method", x = "Method",
       y = "Root Mean Square Error (RMSE)") +
  coord_cartesian(ylim = c(6.5, 12)) +
  theme_minimal() +
  scale_fill_viridis_d(option = "mako") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
ggsave("tick.RMSE.png", width = 8.7, height = 5.5)


table.tick.avg$Method <- reorder(table.tick.avg$Method, -table.tick.avg$DIC)
ggplot(table.tick.avg, aes(x = Method, y = DIC, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  labs(title = "DIC by Method", x = "Method",
       y = "Deviance Information Criterion (DIC)") +
  theme_minimal() +
  coord_cartesian(ylim = c(9900, 10700)) +
  scale_fill_viridis_d(option = "mako") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
ggsave("tick.DIC.png", width = 8.7, height = 5.5)


save(table.tick,
     table.tick.avg, file = "Dataset.Tick.RData")

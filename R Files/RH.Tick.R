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
  geom_point(size = 3) +
  facet_wrap(~ Species) +
  labs(x = "Temperature",
       y = "Oviposition Time in Days") +
  theme_minimal()+ 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )  
#ggsave("Tick.Data.png", width = 10, height = 10)


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
ni <- 40000 # number of iterations in each chain
nb <- 20000 # number of 'burn in' iterations to discard
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
# # Get rid of log liks
# full.data.mcmc <- lapply(as.mcmc.list(full.data.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(full.data.mcmc), autoburnin = FALSE)
# DIC
full.data.fit$BUGSoutput$DIC

# # Get samples of all chains
# full.conv <- do.call(rbind, lapply(full.data.fit.mcmc, function(chain) as.data.frame(chain[, c("a", "Topt", "c")])))
# # Save pairs plots
# png('tick.full.pairs.png', height = 5, width = 5, units = "in", res = 250)
# pairs(full.conv)
# dev.off()
# # Save trace plots
# pdf('tick.full.trace.pdf')
# plot(full.data.fit.mcmc[, c("a", "Topt", "c", "tau1")])
# dev.off()

# full.mcmc.list <- as.mcmc.list(full.data.fit.mcmc)
# combined_samples.a <- unlist(lapply(full.mcmc.list, function(chain) chain[, "a"]))
# combined_samples.topt <- unlist(lapply(full.mcmc.list, function(chain) chain[, "Topt"]))
# combined_samples.c <- unlist(lapply(full.mcmc.list, function(chain) chain[, "c"]))
# png('tick.full.hist.png', height = 3, width = 9, units = "in", res = 250)
# par(mfrow = c(1,3))
# hist(combined_samples.a, xlab = "a", main = "Posterior Samples of a")
# hist(combined_samples.topt, xlab = "Topt", main = "Posterior Samples of Topt")
# hist(combined_samples.c, xlab = "c", main = "Posterior Samples of c")
# dev.off()

# # Extract log-likelihood matrix (samples × observations)
# log.lik.full <- full.data.fit$BUGSoutput$sims.list$log_lik
# # Compute WAIC
# waic.full <- waic(log.lik.full)$estimates["waic", "Estimate"]

# Get posterior samples
full.data.samp <- full.data.fit$BUGSoutput$sims.list
# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp), max(tick.abr.new$Temp), by = 0.1)
# Evaluate each posterior sample
full.eval <- matrix(NA, nrow = length(full.data.samp$a), ncol = length(temp_seq))
for (i in 1:length(full.data.samp$a)) {
  a_i <- full.data.samp$a[i]
  Topt_i <- full.data.samp$Topt[i]
  c_i <- full.data.samp$c[i]
  full.eval[i, ] <- a_i * (Topt_i - temp_seq)^2 + c_i
}
# Get the mean to plot
mean_curve.full <- apply(full.eval, 2, mean)

# Data frame of temperature, fitted value, species
fitted_data.full.data <- data.frame(Temperature = rep(temp_seq, n.species), 
                                   Fitted = mean_curve.full,
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
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )
ggsave("tick.full.plot.png", width = 9, height = 6)

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


save(mse_by_full_data,
     overall_rmse,
     fitted_data.full.data,
     merged_data.full,
     full.data.fit.mcmc,
     full.data.fit,
     log.lik.full,
     file="full.data.fit.RData")


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

# # Get rid of log liks
# species.mcmc <- lapply(as.mcmc.list(species.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(species.mcmc), autoburnin = FALSE)

# 
# # Extract log-likelihood matrix (samples × observations)
# log.lik.species <- species.fit$BUGSoutput$sims.list$log_lik
# # Compute WAIC
# waic.species <- waic(log.lik.species)$estimates["waic", "Estimate"]



# Get posterior samples
species.samp <- as.matrix(species.fit.mcmc)
# Get column indices
a_species <- grep("^mu\\.species\\.a\\[", colnames(species.samp))
Topt_species <- grep("^mu\\.species\\.Topt\\[", colnames(species.samp))
c_species <- grep("^mu\\.species\\.c\\[", colnames(species.samp))
# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp), max(tick.abr.new$Temp), by = 0.1)
# Prepare a list to store data frames
species_dfs <- list()
for (j in 1:n.species) {
  # Retrieve parameter vectors for species j
  a_samps <- species.samp[, a_species[j]]
  Topt_samps <- species.samp[, Topt_species[j]]
  c_samps <- species.samp[, c_species[j]]
  # For each temperature in temp_seq, evaluate across ALL posterior samples
  species_eval <- sapply(temp_seq, function(temp) {
    a_samps * (Topt_samps - temp)^2 + c_samps
  })
  # Compute pointwise posterior means and HDIs
  mean_curve.species <- apply(species_eval, 2, mean)
  species_dfs[[j]] <- mean_curve.species
}
# Combine all species
species_dfs <- unlist(species_dfs)
# Make data frame of temperature, fitted value for each species, and species
fitted_data.species <- data.frame(Temperature = rep(temp_seq, n.species), 
                                  Fitted = species_dfs,
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
  labs(x = "Temperature", y = "Ovipositon Time in Days") +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  ) 
ggsave("tick.species.plot.png", width = 9, height = 6)

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


save(mse_by_species,
     fitted_data.species,
     merged_data.species,
     species.fit.mcmc,
     species.fit,
     log.lik.species,
     file="species.fit.RData")

################################### By Species Hierarchical Different Variance ###################################
# Model using species effect
sink("species_tau_robin.txt")
cat("
model {
  # Priors
  # Species-specific 
  # mu.c ~ dexp(10)          # Mean for species effects c (vertex y)
  # mu.Topt ~ dnorm(30, 0.1)          # Mean for species effects Topt (vertex x)
  # tau.Topt ~ dexp(0.01)             # Tau for species effect Topt (vertex x)
  # mu.a ~ dexp(0.01)            # Mean for species effects a (width)
  mu.a ~ dunif(0, 2000)
  mu.Topt ~ dexp(0.001)
  mu.c ~ dgamma(10, 1)
  tau1 ~ dexp(0.01)           # Tau for species effect

  for (j in 1:n.species) {
    # mu.species.c[j] ~ dexp(mu.c)       # Species-specific vertex y
    # mu.species.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Species-specific vertex x
    # mu.species.a[j] ~ dexp(mu.a)    # Species-specific width
    mu.species.c[j] ~ dgamma(mu.c, 1)
    mu.species.Topt[j] ~ dexp(mu.Topt)
    mu.species.a[j] ~ dunif(0, mu.a)
    species.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs) {
    #mu[i] <- mu.species.a[species[i]] * (mu.species.Topt[species[i]] - temp[i])^2 + mu.species.c[species[i]]
    mu[i] <- mu.species.a[species[i]] * exp(- mu.species.Topt[species[i]] * temp[i]) + mu.species.c[species[i]]
    trait[i] ~ dnorm(mu[i], species.tau[species[i]]) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], species.tau[species[i]])
  }}
", file = "species_tau_robin.txt")


# Parameters to monitor
parameters <- c("mu.species.a", "mu.species.c", "mu.a", "mu.c",
                "mu.Topt", #"tau.Topt",
                "tau1", "mu.species.Topt", "species.tau", "log_lik")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Initial Values
inits = vector('list', nc)
GenInits = function() {
  # mu.a <- rexp(1, 0.01)
  # mu.c <- rexp(1, 10)
   tau1 <- rexp(1, 0.01)
  # mu.Topt <- rnorm(1, 30, 1/0.1)
  # tau.Topt <- rexp(1, 0.01)
  mu.a <- runif(1, 0, 2000)
  mu.c <- rgamma(1, 10, 1)
  mu.Topt <- rexp(1, 0.01)
  list(
    mu.a = mu.a,             
    mu.c = mu.c,
    tau1 = tau1,
    mu.Topt = mu.Topt
    #tau.Topt = tau.Topt
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
# 
# # Get rid of log liks
# species.tau.mcmc <- lapply(as.mcmc.list(species.tau.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(species.tau.mcmc), autoburnin = FALSE)

# # Extract log-likelihood matrix (samples × observations)
# log.lik.species.tau <- species.tau.fit$BUGSoutput$sims.list$log_lik
# # Compute WAIC
# waic.species.tau <- waic(log.lik.species.tau)$estimates["waic", "Estimate"]

# Get posterior samples
species.tau.samp <- as.matrix(species.tau.fit.mcmc)
# Get column indices
a_species.tau <- grep("^mu\\.species\\.a\\[", colnames(species.tau.samp))
Topt_species.tau <- grep("^mu\\.species\\.Topt\\[", colnames(species.tau.samp))
c_species.tau <- grep("^mu\\.species\\.c\\[", colnames(species.tau.samp))
# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp), max(tick.abr.new$Temp), by = 0.1)
# Prepare a list to store data frames
species.tau_dfs <- list()
for (j in 1:n.species) {
  # Retrieve parameter vectors for species j
  a_samps <- species.tau.samp[, a_species.tau[j]]
  Topt_samps <- species.tau.samp[, Topt_species.tau[j]]
  c_samps <- species.tau.samp[, c_species.tau[j]]
  # For each temperature in temp_seq, evaluate across ALL posterior samples
  species.tau_eval <- sapply(temp_seq, function(temp) {
    a_samps * (Topt_samps - temp)^2 + c_samps
  })
  # Compute pointwise posterior means and HDIs
  mean_curve.species.tau <- apply(species.tau_eval, 2, mean)
  species.tau_dfs[[j]] <- mean_curve.species.tau
}
# Combine all species
species.tau_dfs <- unlist(species.tau_dfs)

# Make data frame of temperature, fitted value for each species, and species
fitted_data.species.tau <- data.frame(Temperature = rep(temp_seq, n.species), 
                                  Fitted = species.tau_dfs,
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
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  ) 
ggsave("tick.species.tau.plot.png", width = 9, height = 6)


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

save(mse_by_species.tau,
     fitted_data.species.tau,
     merged_data.species.tau,
     species.tau.fit.mcmc,
     species.tau.fit,
     log.lik.species.tau,
     file="species.tau.fit.RData")

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

# # Get rid of log liks
# species.tau.param.mcmc <- lapply(as.mcmc.list(species.tau.param.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(species.tau.param.mcmc), autoburnin = FALSE)

# Extract log-likelihood matrix (samples × observations)
log.lik.species.tau.param <- species.tau.param.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.species.tau.param <- waic(log.lik.species.tau.param)$estimates["waic", "Estimate"]

# Get posterior samples
species.tau.param.samp <- species.tau.param.fit$BUGSoutput$sims.list
# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp), max(tick.abr.new$Temp), by = 0.1)
# Evaluate each posterior sample
species.tau.param.eval <- matrix(NA, nrow = length(species.tau.param.samp$a), ncol = length(temp_seq))
for (i in 1:length(species.tau.param.samp$a)) {
  a_i <- species.tau.param.samp$a[i]
  Topt_i <- species.tau.param.samp$Topt[i]
  c_i <- species.tau.param.samp$c[i]
  species.tau.param.eval[i, ] <- a_i * (Topt_i - temp_seq)^2 + c_i
}
# Get the mean to plot
mean_curve.species.tau.param <- apply(species.tau.param.eval, 2, mean)

# Make data frame of temperature, fitted value for each species, and species
fitted_data.species.tau.param <- data.frame(Temperature = rep(temp_seq, n.species), 
                                            Fitted = mean_curve.species.tau.param,
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
  geom_line(data = fitted_data.species.tau.param, aes(x = Temperature, y = Fitted), color = "forestgreen", size = 1) +
  labs(x = "Temperature", y = "Oviposition Time in Days") +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  ) 
ggsave("tick.species.tau.param.plot.png", width = 9, height = 6)

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

save(mse_by_species.tau.param,
     fitted_data.species.tau.param,
     merged_data.species.tau.param,
     species.tau.param.fit.mcmc,
     species.tau.param.fit,
     log.lik.species.tau.param,
     file="species.tau.param.fit.RData")

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


# Get rid of log liks
# lepidum.mcmc <- lapply(as.mcmc.list(a.lepidum.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(lepidum.mcmc), autoburnin = FALSE)

# Get posterior samples
a.lepidum.samp <- a.lepidum.fit$BUGSoutput$sims.list
# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp), max(tick.abr.new$Temp), by = 0.1)
# Evaluate each posterior sample
a.lepidum.eval <- matrix(NA, nrow = length(a.lepidum.samp$a), ncol = length(temp_seq))
for (i in 1:length(a.lepidum.samp$a)) {
  a_i <- a.lepidum.samp$a[i]
  Topt_i <- a.lepidum.samp$Topt[i]
  c_i <- a.lepidum.samp$c[i]
  a.lepidum.eval[i, ] <- a_i * (Topt_i - temp_seq)^2 + c_i
}
# Get the mean to plot
mean_curve.a.lepidum <- apply(a.lepidum.eval, 2, mean)


# Data frame of temperature, fitted value, species
fitted_data.a.lepidum <- data.frame(Temperature = rep(temp_seq),
                                     Fitted = mean_curve.a.lepidum)
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
     #waic.a.lepidum,
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

# Get rid of log liks
# andersoni.mcmc <- lapply(as.mcmc.list(d.andersoni.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(andersoni.mcmc), autoburnin = FALSE)

# Get posterior samples
d.andersoni.samp <- d.andersoni.fit$BUGSoutput$sims.list
# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp), max(tick.abr.new$Temp), by = 0.1)
# Evaluate each posterior sample
d.andersoni.eval <- matrix(NA, nrow = length(d.andersoni.samp$a), ncol = length(temp_seq))
for (i in 1:length(d.andersoni.samp$a)) {
  a_i <- d.andersoni.samp$a[i]
  Topt_i <- d.andersoni.samp$Topt[i]
  c_i <- d.andersoni.samp$c[i]
  d.andersoni.eval[i, ] <- a_i * (Topt_i - temp_seq)^2 + c_i
}
# Get the mean to plot
mean_curve.d.andersoni <- apply(d.andersoni.eval, 2, mean)

# Data frame of temperature, fitted value, species
fitted_data.d.andersoni <- data.frame(Temperature = rep(temp_seq),
                                     Fitted = mean_curve.d.andersoni)

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
     #waic.d.andersoni,
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

# Get rid of log liks
# nitens.mcmc <- lapply(as.mcmc.list(d.nitens.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(nitens.mcmc), autoburnin = FALSE)

# Get posterior samples
d.nitens.samp <- d.nitens.fit$BUGSoutput$sims.list
# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp), max(tick.abr.new$Temp), by = 0.1)
# Evaluate each posterior sample
d.nitens.eval <- matrix(NA, nrow = length(d.nitens.samp$a), ncol = length(temp_seq))
for (i in 1:length(d.nitens.samp$a)) {
  a_i <- d.nitens.samp$a[i]
  Topt_i <- d.nitens.samp$Topt[i]
  c_i <- d.nitens.samp$c[i]
  d.nitens.eval[i, ] <- a_i * (Topt_i - temp_seq)^2 + c_i
}
# Get the mean to plot
mean_curve.d.nitens <- apply(d.nitens.eval, 2, mean)

# Data frame of temperature, fitted value, species
fitted_data.d.nitens <- data.frame(Temperature = rep(temp_seq),
                                     Fitted = mean_curve.d.nitens)
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
     #waic.d.nitens,
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

# Get rid of log liks
# leporispalustris.mcmc <- lapply(as.mcmc.list(h.leporispalustris.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(leporispalustris.mcmc), autoburnin = FALSE)

# Get posterior samples
h.leporispalustris.samp <- h.leporispalustris.fit$BUGSoutput$sims.list
# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp), max(tick.abr.new$Temp), by = 0.1)
# Evaluate each posterior sample
h.leporispalustris.eval <- matrix(NA, nrow = length(h.leporispalustris.samp$a), ncol = length(temp_seq))
for (i in 1:length(h.leporispalustris.samp$a)) {
  a_i <- h.leporispalustris.samp$a[i]
  Topt_i <- h.leporispalustris.samp$Topt[i]
  c_i <- h.leporispalustris.samp$c[i]
  h.leporispalustris.eval[i, ] <- a_i * (Topt_i - temp_seq)^2 + c_i
}
# Get the mean to plot
mean_curve.h.leporispalustris <- apply(h.leporispalustris.eval, 2, mean)


# Data frame of temperature, fitted value, species
fitted_data.h.leporispalustris <- data.frame(Temperature = rep(temp_seq),
                                     Fitted = mean_curve.h.leporispalustris)
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
     #waic.h.leporispalustris,
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

# Get rid of log liks
# aegyptium.mcmc <- lapply(as.mcmc.list(h.aegyptium.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(aegyptium.mcmc), autoburnin = FALSE)

# Get posterior samples
h.aegyptium.samp <- h.aegyptium.fit$BUGSoutput$sims.list
# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp), max(tick.abr.new$Temp), by = 0.1)
# Evaluate each posterior sample
h.aegyptium.eval <- matrix(NA, nrow = length(h.aegyptium.samp$a), ncol = length(temp_seq))
for (i in 1:length(h.aegyptium.samp$a)) {
  a_i <- h.aegyptium.samp$a[i]
  Topt_i <- h.aegyptium.samp$Topt[i]
  c_i <- h.aegyptium.samp$c[i]
  h.aegyptium.eval[i, ] <- a_i * (Topt_i - temp_seq)^2 + c_i
}
# Get the mean to plot
mean_curve.h.aegyptium <- apply(h.aegyptium.eval, 2, mean)

# Data frame of temperature, fitted value, species
fitted_data.h.aegyptium <- data.frame(Temperature = rep(temp_seq),
                                     Fitted = mean_curve.h.aegyptium)
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
     #waic.h.aegyptium,
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

# Get rid of log liks
# dromedarii.mcmc <- lapply(as.mcmc.list(h.dromedarii.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(dromedarii.mcmc), autoburnin = FALSE)

# Get posterior samples
h.dromedarii.samp <- h.dromedarii.fit$BUGSoutput$sims.list
# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp), max(tick.abr.new$Temp), by = 0.1)
# Evaluate each posterior sample
h.dromedarii.eval <- matrix(NA, nrow = length(h.dromedarii.samp$a), ncol = length(temp_seq))
for (i in 1:length(h.dromedarii.samp$a)) {
  a_i <- h.dromedarii.samp$a[i]
  Topt_i <- h.dromedarii.samp$Topt[i]
  c_i <- h.dromedarii.samp$c[i]
  h.dromedarii.eval[i, ] <- a_i * (Topt_i - temp_seq)^2 + c_i
}
# Get the mean to plot
mean_curve.h.dromedarii <- apply(h.dromedarii.eval, 2, mean)

# Data frame of temperature, fitted value, species
fitted_data.h.dromedarii <- data.frame(Temperature = rep(temp_seq),
                                     Fitted = mean_curve.h.dromedarii)
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
     #waic.h.dromedarii,
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

# Get rid of log liks
# impeltatum.mcmc <- lapply(as.mcmc.list(h.impeltatum.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(impeltatum.mcmc), autoburnin = FALSE)

# Get posterior samples
h.impeltatum.samp <- h.impeltatum.fit$BUGSoutput$sims.list
# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp), max(tick.abr.new$Temp), by = 0.1)
# Evaluate each posterior sample
h.impeltatum.eval <- matrix(NA, nrow = length(h.impeltatum.samp$a), ncol = length(temp_seq))
for (i in 1:length(h.impeltatum.samp$a)) {
  a_i <- h.impeltatum.samp$a[i]
  Topt_i <- h.impeltatum.samp$Topt[i]
  c_i <- h.impeltatum.samp$c[i]
  h.impeltatum.eval[i, ] <- a_i * (Topt_i - temp_seq)^2 + c_i
}
# Get the mean to plot
mean_curve.h.impeltatum <- apply(h.impeltatum.eval, 2, mean)


# Data frame of temperature, fitted value, species
fitted_data.h.impeltatum <- data.frame(Temperature = rep(temp_seq),
                                     Fitted = mean_curve.h.impeltatum)
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
     #waic.h.impeltatum,
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

# Get rid of log liks
# lusitanicum.mcmc <- lapply(as.mcmc.list(h.lusitanicum.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(lusitanicum.mcmc), autoburnin = FALSE)

# Get posterior samples
h.lusitanicum.samp <- h.lusitanicum.fit$BUGSoutput$sims.list
# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp), max(tick.abr.new$Temp), by = 0.1)
# Evaluate each posterior sample
h.lusitanicum.eval <- matrix(NA, nrow = length(h.lusitanicum.samp$a), ncol = length(temp_seq))
for (i in 1:length(h.lusitanicum.samp$a)) {
  a_i <- h.lusitanicum.samp$a[i]
  Topt_i <- h.lusitanicum.samp$Topt[i]
  c_i <- h.lusitanicum.samp$c[i]
  h.lusitanicum.eval[i, ] <- a_i * (Topt_i - temp_seq)^2 + c_i
}
# Get the mean to plot
mean_curve.h.lusitanicum <- apply(h.lusitanicum.eval, 2, mean)



# Data frame of temperature, fitted value, species
fitted_data.h.lusitanicum <- data.frame(Temperature = rep(temp_seq),
                                     Fitted = mean_curve.h.lusitanicum)
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
     #waic.h.lusitanicum,
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

# Get rid of log liks
# schulzei.mcmc <- lapply(as.mcmc.list(h.schulzei.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(schulzei.mcmc), autoburnin = FALSE)

# Get posterior samples
h.schulzei.samp <- h.schulzei.fit$BUGSoutput$sims.list
# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp), max(tick.abr.new$Temp), by = 0.1)
# Evaluate each posterior sample
h.schulzei.eval <- matrix(NA, nrow = length(h.schulzei.samp$a), ncol = length(temp_seq))
for (i in 1:length(h.schulzei.samp$a)) {
  a_i <- h.schulzei.samp$a[i]
  Topt_i <- h.schulzei.samp$Topt[i]
  c_i <- h.schulzei.samp$c[i]
  h.schulzei.eval[i, ] <- a_i * (Topt_i - temp_seq)^2 + c_i
}
# Get the mean to plot
mean_curve.h.schulzei <- apply(h.schulzei.eval, 2, mean)


# Data frame of temperature, fitted value, species
fitted_data.h.schulzei <- data.frame(Temperature = rep(temp_seq),
                                    Fitted = mean_curve.h.schulzei)
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
     #waic.h.schulzei,
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

# Arrange plots into a 2x3 grid (your original)
big_grid.tick <- grid.arrange(
    p.a.lepidum + labs(title = "Amblyomma lepidum") + 
      theme(plot.title = element_text(size = 12, hjust = 0.5), axis.text.x = element_blank(),
            axis.text.y = element_text(size = 12),
            axis.title.y = element_blank(), axis.title.x = element_blank()) +
      coord_cartesian(ylim = c(0,100)),
    p.d.andersoni + labs(title = "Dermacentor andersoni") + 
      theme(plot.title = element_text(size = 12, hjust = 0.5), axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(), axis.title.x = element_blank()) +
      coord_cartesian(ylim = c(0,100)),
    p.d.nitens + labs(title = "Dermacentor nitens") + 
      theme(plot.title = element_text(size = 12, hjust = 0.5), axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(), axis.title.x = element_blank()) +
      coord_cartesian(ylim = c(0,100)),
    p.h.leporispalustris + labs(title = "Haemaphysalis leporispalustris") + 
      theme(plot.title = element_text(size = 12, hjust = 0.5), axis.text.x = element_blank(),
            axis.text.y = element_text(size = 12),
            axis.title.y = element_blank(), axis.title.x = element_blank()) +
      coord_cartesian(ylim = c(0,100)),
    p.h.aegyptium + labs(title = "Hyalomma aegyptium") + 
      theme(plot.title = element_text(size = 12, hjust = 0.5), axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(), axis.title.x = element_blank()) +
      coord_cartesian(ylim = c(0,100)),
    p.h.dromedarii + labs(title = "Hyalomma dromedarii") + 
      theme(plot.title = element_text(size = 12, hjust = 0.5), axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(), axis.title.x = element_blank()) +
      coord_cartesian(ylim = c(0,100)),
    p.h.impeltatum + labs(title = "Hyalomma impeltatum") + 
      theme(plot.title = element_text(size = 12, hjust = 0.5), axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            axis.title.y = element_blank(), axis.title.x = element_blank()) +
      coord_cartesian(ylim = c(0,100)),
    p.h.lusitanicum + labs(title = "Hyalomma lusitanicum") + 
      theme(plot.title = element_text(size = 12, hjust = 0.5), axis.text.x = element_text(size = 12),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(), axis.title.x = element_blank()) +
      coord_cartesian(ylim = c(0,100)),
    p.h.schulzei + labs(title = "Hyalomma schulzei") + 
      theme(plot.title = element_text(size = 12, hjust = 0.5), axis.text.x = element_text(size = 12),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(), axis.title.x = element_blank()) +
      coord_cartesian(ylim = c(0,100)),
    nrow = 3, ncol = 3
  )

# Create y-axis label grob (rotated)
y_label.tick <- textGrob("Oviposition Time in Days", rot = 90, gp = gpar(fontsize = 12))

# Combine y-axis label and grid side-by-side with some margin for y-label
ind.species.pd<- arrangeGrob(
  y_label.tick, big_grid.tick,
  ncol = 2,
  widths = unit.c(unit(0.5, "cm"), unit(1, "npc") - unit(0.5, "cm"))
)

# Create x-axis label grob
x_label <- textGrob("Temperature", gp = gpar(fontsize = 12))

# Add x-axis label below with some margin
ind.species.pd <- arrangeGrob(
  ind.species.pd,
  x_label,
  ncol = 1,
  heights = unit.c(unit(1, "npc") - unit(0.6, "cm"), unit(0.6, "cm"))
)

ggsave("tick.species.ind.plot.png", ind.species.pd, width = 8, height = 5)

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

# Get rid of log liks
# genus.mcmc <- lapply(as.mcmc.list(genus.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(genus.mcmc), autoburnin = FALSE)

# Extract log-likelihood matrix (samples × observations)
log.lik.genus <- genus.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.genus <- waic(log.lik.genus)$estimates["waic", "Estimate"]




# Get posterior samples
genus.samp <- as.matrix(genus.fit.mcmc)
# Get column indices
a_genus <- grep("^mu\\.genus\\.a\\[", colnames(genus.samp))
Topt_genus <- grep("^mu\\.genus\\.Topt\\[", colnames(genus.samp))
c_genus <- grep("^mu\\.genus\\.c\\[", colnames(genus.samp))
# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp), max(tick.abr.new$Temp), by = 0.1)
# Prepare a list to store data frames
genus_dfs <- list()
for (j in 1:n.genus) {
  # Retrieve parameter vectors for species j
  a_samps <- genus.samp[, a_genus[j]]
  Topt_samps <- genus.samp[, Topt_genus[j]]
  c_samps <- genus.samp[, c_genus[j]]
  # For each temperature in temp_seq, evaluate across ALL posterior samples
  genus_eval <- sapply(temp_seq, function(temp) {
    a_samps * (Topt_samps - temp)^2 + c_samps
  })
  # Compute pointwise posterior means and HDIs
  mean_curve.genus <- apply(genus_eval, 2, mean)
  genus_dfs[[j]] <- mean_curve.genus
}
# Combine all species
genus_dfs <- unlist(genus_dfs)
# Data frame of temperature, fitted values, and genus
fitted_data.genus <- data.frame(Temperature = rep(temp_seq, n.genus), 
                                Fitted = genus_dfs,
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
  Fitted = genus_dfs,
  Genus = factor(rep(c("Amblyomma", "Dermacentor",
                       "Haemaphysalis", "Hyalomma"), each = length(temp_seq)))) %>%
  left_join(species_to_genus, by = "Genus")

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.genus, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.genus, aes(x = Temperature, y = Fitted, color = Genus), size = 1) +
  labs(x = "Temperature", y = "Oviposition Time in Days") +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  ) 
ggsave("tick.genus.plot.png", width = 9, height = 6)

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



save(mse_by_genus,
     fitted_data.genus,
     merged_data.genus,
     genus.fit.mcmc,
     genus.fit,
     log.lik.genus,
     file="genus.fit.RData")

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

# Get rid of log liks
# climate.mcmc <- lapply(as.mcmc.list(climate.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(climate.mcmc), autoburnin = FALSE)


# Extract log-likelihood matrix (samples × observations)
log.lik.climate <- climate.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.climate <- waic(log.lik.climate)$estimates["waic", "Estimate"]

# Get posterior samples
climate.samp <- as.matrix(climate.fit.mcmc)
# Get column indices
a_climate <- grep("^mu\\.climate\\.a\\[", colnames(climate.samp))
Topt_climate <- grep("^mu\\.climate\\.Topt\\[", colnames(climate.samp))
c_climate <- grep("^mu\\.climate\\.c\\[", colnames(climate.samp))
# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp), max(tick.abr.new$Temp), by = 0.1)
# Prepare a list to store data frames
climate_dfs <- list()
for (j in 1:n.climate) {
  # Retrieve parameter vectors for species j
  a_samps <- climate.samp[, a_climate[j]]
  Topt_samps <- climate.samp[, Topt_climate[j]]
  c_samps <- climate.samp[, c_climate[j]]
  # For each temperature in temp_seq, evaluate across ALL posterior samples
  climate_eval <- sapply(temp_seq, function(temp) {
    a_samps * (Topt_samps - temp)^2 + c_samps
  })
  # Compute pointwise posterior means and HDIs
  mean_curve.climate <- apply(climate_eval, 2, mean)
  climate_dfs[[j]] <- mean_curve.climate
}
# Combine all species
climate_dfs <- unlist(climate_dfs)

# Data frame of temperature, fitted values, and climate
fitted_data.climate <- data.frame(Temperature = rep(temp_seq, n.climate), 
                                  Fitted = climate_dfs,
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
  Fitted = climate_dfs,
  Climate = factor(rep(c("Mixed", "Subtropical", 
                         "Temperate", "Tropical"), each = length(temp_seq)))) %>%
  left_join(species_to_climate, by = "Climate")

plot_data.climate$Climate <- factor(plot_data.climate$Climate, levels = c("Mixed", "Temperate", "Subtropical", "Tropical"))
fitted_data.climate$Climate <- factor(fitted_data.climate$Climate, levels = c("Mixed", "Temperate", "Subtropical", "Tropical"))
# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.climate, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.climate, aes(x = Temperature, y = Fitted, color = Climate), size = 1) +
  labs(x = "Temperature", y = "Oviposition Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  ) 
ggsave("tick.climate.plot.png", width = 11, height = 6)

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


save(mse_by_climate,
     fitted_data.climate,
     merged_data.climate,
     climate.fit.mcmc,
     climate.fit,
     log.lik.climate,
     file="climate.fit.RData")

################################### By Host ###################################
# Model by typical number of host species
sink("host_robin.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Host-specific 
  # mu.c ~ dexp(10)          # Mean for host effects c (vertex y)
  # mu.Topt ~ dnorm(30, 0.1)          # Mean for host effects Topt (vertex x)
  # tau.Topt ~ dexp(0.01)             # Tau for host effect Topt (vertex x)
  # mu.a ~ dexp(0.01)            # Mean for host effects a (width)
  mu.a ~ dunif(0, 2000)
  mu.Topt ~ dexp(0.001)
  mu.c ~ dgamma(10, 1)
  
  for (j in 1:n.host) {
    # mu.host.c[j] ~ dexp(mu.c)       # Host-specific vertex y
    # mu.host.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Host-specific vertex x
    # mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
    host.tau[j] ~ dexp(tau1)
    mu.host.c[j] ~ dgamma(mu.c, 1)
    mu.host.Topt[j] ~ dexp(mu.Topt)
    mu.host.a[j] ~ dunif(0, mu.a)
  }
  # Likelihood
  for (i in 1:N.obs) {
    #mu[i] <- mu.host.a[host[i]] * (mu.host.Topt[host[i]] - temp[i])^2 + mu.host.c[host[i]]
    mu[i] <- mu.host.a[host[i]] * exp(-mu.host.Topt[host[i]] * temp[i]) + mu.host.c[host[i]]
    trait[i] ~ dnorm(mu[i], host.tau[host[i]]) T(0, ) 
    log_lik[i] <- logdensity.norm(trait[i], mu[i], host.tau[host[i]])
  }}
", file = "host_robin.txt")


# Parameters to monitor
parameters <- c("mu.host.a", "mu.host.c", "mu.c", "mu.a",
                "mu.Topt", "host.tau",
                "tau1", "mu.host.Topt", "log_lik") 


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Initial Values
inits = vector('list', nc)
GenInits = function() {
  # mu.a <- rexp(1, 0.01)
  # mu.c <- rexp(1, 10)
  tau1 <- rexp(1, 0.01)
  #mu.Topt <- rnorm(1, 30, 1/0.1)
  #tau.Topt <- rexp(1, 0.01)
  mu.a <- runif(1, 0, 2000)
  mu.Topt <- rexp(1, 0.01)
  mu.c <- dgamma(1, 10, 1)
  list(
    mu.a = mu.a,             
    mu.c = mu.c,
    tau1 = tau1,
    mu.Topt = mu.Topt
    #tau.Topt = tau.Topt
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

# Get rid of log liks
# host.mcmc <- lapply(as.mcmc.list(host.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(host.mcmc), autoburnin = FALSE)


# Extract log-likelihood matrix (samples × observations)
log.lik.host <- host.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.host <- waic(log.lik.host)$estimates["waic", "Estimate"]

# Get posterior samples
host.samp <- as.matrix(host.fit.mcmc)
# Get column indices
a_host <- grep("^mu\\.host\\.a\\[", colnames(host.samp))
Topt_host <- grep("^mu\\.host\\.Topt\\[", colnames(host.samp))
c_host <- grep("^mu\\.host\\.c\\[", colnames(host.samp))
# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp), max(tick.abr.new$Temp), by = 0.1)
# Prepare a list to store data frames
host_dfs <- list()
for (j in 1:n.host) {
  # Retrieve parameter vectors for species j
  a_samps <- host.samp[, a_host[j]]
  Topt_samps <- host.samp[, Topt_host[j]]
  c_samps <- host.samp[, c_host[j]]
  # For each temperature in temp_seq, evaluate across ALL posterior samples
  host_eval <- sapply(temp_seq, function(temp) {
    a_samps * (Topt_samps - temp)^2 + c_samps
  })
  # Compute pointwise posterior means and HDIs
  mean_curve.host <- apply(host_eval, 2, mean)
  host_dfs[[j]] <- mean_curve.host
}
# Combine all species
host_dfs <- unlist(host_dfs)


# Data frame with temperature, fitted value, and host group
fitted_data.host <- data.frame(Temperature = rep(temp_seq, n.host), 
                                Fitted = host_dfs,
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
  Fitted = host_dfs,
  Host = factor(rep(c("One", "Three", "Two"), each = length(temp_seq)))
) %>%
  left_join(species_to_host, by = "Host")

plot_data.host$Host <- factor(plot_data.host$Host, levels = c("One", "Two", "Three"))
fitted_data.host$Host <- factor(fitted_data.host$Host, levels = c("One", "Two", "Three"))
# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.host, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host, aes(x = Temperature, y = Fitted, color = Host), size = 1) +
  labs(x = "Temperature", y = "Oviposition Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  ) 
ggsave("tick.host.plot.png", width = 11, height = 6)


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


save(mse_by_host,
     fitted_data.host,
     merged_data.host,
     host.fit.mcmc,
     host.fit,
     log.lik.host,
     file="host.fit.RData")
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

# Get rid of log liks
# size.mcmc <- lapply(as.mcmc.list(size.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(size.mcmc), autoburnin = FALSE)


# Extract log-likelihood matrix (samples × observations)
log.lik.size <- size.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.size <- waic(log.lik.size)$estimates["waic", "Estimate"]

# Get posterior samples
size.samp <- as.matrix(size.fit.mcmc)
# Get column indices
a_size <- grep("^mu\\.size\\.a\\[", colnames(size.samp))
Topt_size <- grep("^mu\\.size\\.Topt\\[", colnames(size.samp))
c_size <- grep("^mu\\.size\\.c\\[", colnames(size.samp))
# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp), max(tick.abr.new$Temp), by = 0.1)
# Prepare a list to store data frames
size_dfs <- list()
for (j in 1:n.size) {
  # Retrieve parameter vectors for species j
  a_samps <- size.samp[, a_size[j]]
  Topt_samps <- size.samp[, Topt_size[j]]
  c_samps <- size.samp[, c_size[j]]
  # For each temperature in temp_seq, evaluate across ALL posterior samples
  size_eval <- sapply(temp_seq, function(temp) {
    a_samps * (Topt_samps - temp)^2 + c_samps
  })
  # Compute pointwise posterior means and HDIs
  mean_curve.size <- apply(size_eval, 2, mean)
  size_dfs[[j]] <- mean_curve.size
}
# Combine all species
size_dfs <- unlist(size_dfs)

# Make data frame with temperature, fitted value, and size category
fitted_data.size <- data.frame(Temperature = rep(temp_seq, n.size), 
                               Fitted = size_dfs,
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
  Fitted = size_dfs,
  Size = factor(rep(c("Large", "Medium", "Small"), each = length(temp_seq)))
) %>%
  left_join(species_to_size, by = "Size")

plot_data.size$Size <- factor(plot_data.size$Size, levels = c("Small", "Medium", "Large"))
fitted_data.size$Size <- factor(fitted_data.size$Size, levels = c("Small", "Medium", "Large"))
# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.size, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.size, aes(x = Temperature, y = Fitted, color = Size), size = 1) +
  labs(x = "Temperature", y = "Oviposition Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  ) 
ggsave("tick.size.plot.png", width = 11, height = 6)


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

save(mse_by_size,
     fitted_data.size,
     merged_data.size,
     size.fit.mcmc,
     size.fit,
     log.lik.size,
     file="size.fit.RData")

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

# Get rid of log liks
# continent.mcmc <- lapply(as.mcmc.list(continent.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(continent.mcmc), autoburnin = FALSE)


# Extract log-likelihood matrix (samples × observations)
log.lik.continent <- continent.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.continent <- waic(log.lik.continent)$estimates["waic", "Estimate"]

# Get posterior samples
continent.samp <- as.matrix(continent.fit.mcmc)
# Get column indices
a_continent <- grep("^mu\\.continent\\.a\\[", colnames(continent.samp))
Topt_continent <- grep("^mu\\.continent\\.Topt\\[", colnames(continent.samp))
c_continent <- grep("^mu\\.continent\\.c\\[", colnames(continent.samp))
# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp), max(tick.abr.new$Temp), by = 0.1)
# Prepare a list to store data frames
continent_dfs <- list()
for (j in 1:n.continent) {
  # Retrieve parameter vectors for species j
  a_samps <- continent.samp[, a_continent[j]]
  Topt_samps <- continent.samp[, Topt_continent[j]]
  c_samps <- continent.samp[, c_continent[j]]
  # For each temperature in temp_seq, evaluate across ALL posterior samples
  continent_eval <- sapply(temp_seq, function(temp) {
    a_samps * (Topt_samps - temp)^2 + c_samps
  })
  # Compute pointwise posterior means and HDIs
  mean_curve.continent <- apply(continent_eval, 2, mean)
  continent_dfs[[j]] <- mean_curve.continent
}
# Combine all species
continent_dfs <- unlist(continent_dfs)

# Data frame of temperature, fitted value, and continent group
fitted_data.continent <- data.frame(Temperature = rep(temp_seq, n.continent), 
                               Fitted = continent_dfs,
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
  Fitted = continent_dfs,
  Continent = factor(rep(c("More than two", "One", "Two"), each = length(temp_seq)))
) %>%
  left_join(species_to_continent, by = "Continent")

plot_data.continent$Continent <- factor(plot_data.continent$Continent, levels = c("One", "Two", "More than two"))
fitted_data.continent$Continent <- factor(fitted_data.continent$Continent, levels = c("One", "Two", "More than two"))
# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.continent, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.continent, aes(x = Temperature, y = Fitted, color = Continent), size = 1) +
  labs(x = "Temperature", y = "Oviposition Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )
ggsave("tick.continent.plot.png", width = 11, height = 6)

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

save(mse_by_continent,
     fitted_data.continent,
     merged_data.continent,
     continent.fit.mcmc,
     continent.fit,
     log.lik.continent,
     file="continent.fit.RData")

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

# Get rid of log liks
# host.climate.mcmc <- lapply(as.mcmc.list(host.climate.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(host.climate.mcmc), autoburnin = FALSE)



# Extract log-likelihood matrix (samples × observations)
log.lik.host.climate <- host.climate.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.host.climate <- waic(log.lik.host.climate)$estimates["waic", "Estimate"]

# Get posterior samples
host.climate.samp <- as.matrix(host.climate.fit.mcmc)
host.climate_vector <- numeric()
# Get column indices
a_host.climate <- grep("^mu\\.host\\.a\\[", colnames(host.climate.samp))
Topt_host.climate <- grep("^mu\\.climate\\.Topt\\[", colnames(host.climate.samp))
c_host.climate <- grep("^mu\\.host\\.c\\[", colnames(host.climate.samp))
# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp), max(tick.abr.new$Temp), by = 0.1)
# Prepare a list to store data frames
host.climate_dfs <- list()
for(k in 1:n.host){
  for (j in 1:n.climate) {
    # Retrieve parameter vectors for species j
    a_samps <- host.climate.samp[, a_host.climate[k]]
    Topt_samps <- host.climate.samp[, Topt_host.climate[j]]
    c_samps <- host.climate.samp[, c_host.climate[k]]
    # For each temperature in temp_seq, evaluate across ALL posterior samples
    host.climate_eval <- sapply(temp_seq, function(temp) {
      a_samps * (Topt_samps - temp)^2 + c_samps
    })
    # Compute pointwise posterior means and HDIs
    mean_curve <- colMeans(host.climate_eval) # length(temp_seq)
    host.climate_vector <- c(host.climate_vector, mean_curve)
  }}




# Data frame with temperature, fitted value, and host group
fitted_data.host.climate <- data.frame(Temperature = rep(temp_seq, n.host * n.climate), 
                               Fitted = host.climate_vector,
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
                                       Fitted = host.climate_vector,
                                       Host = factor(rep(rep(c("One", "Three", "Two"), each = length(temp_seq)), n.climate)),
                                       Climate = factor(rep(c("Mixed", "Subtropical", "Temperate", "Tropical"),
                                                            each = length(temp_seq) * n.size))) %>%
  left_join(species_to_host.climate, by = c("Host", "Climate"))|> 
  filter(!is.na(Species)) 

plot_data.host.climate$Host <- factor(plot_data.host.climate$Host, levels = c("One", "Two", "Three"))
fitted_data.host.climate$Host <- factor(fitted_data.host.climate$Host, levels = c("One", "Two", "Three"))

# Combine Climate and Size for coloring
fitted_data.host.climate$`Host, Climate` <- with(fitted_data.host.climate, paste(Host, Climate, sep = ", "))
fitted_data.host.climate$`Host, Climate` <- factor(fitted_data.host.climate$`Host, Climate`, levels = c("One, Tropical", 
                                                                                        "Two, Subtropical", "Three, Mixed",
                                                                                        "Three, Subtropical", "Three, Temperate",
                                                                                    "Three, Tropical"))

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.host.climate, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.climate, aes(x = Temperature, y = Fitted, color = `Host, Climate`), size = 1) +
  labs(x = "Temperature", y = "Oviposition Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  ) 
ggsave("tick.host.climate.plot.png", width = 11.25, height = 6)

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


save(mse_by_host.climate,
     fitted_data.host.climate,
     merged_data.host.climate,
     host.climate.fit.mcmc,
     host.climate.fit,
     log.lik.host.climate,
     file="host.climate.fit.RData")

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

# Get rid of log liks
# host.genus.mcmc <- lapply(as.mcmc.list(host.genus.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(host.genus.mcmc), autoburnin = FALSE)


# Extract log-likelihood matrix (samples × observations)
log.lik.host.genus <- host.genus.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.host.genus <- waic(log.lik.host.genus)$estimates["waic", "Estimate"]
# Get posterior samples
host.genus.samp <- as.matrix(host.genus.fit.mcmc)
host.genus_vector <- numeric()
# Get column indices
a_host.genus <- grep("^mu\\.host\\.a\\[", colnames(host.genus.samp))
Topt_host.genus <- grep("^mu\\.genus\\.Topt\\[", colnames(host.genus.samp))
c_host.genus <- grep("^mu\\.genus\\.c\\[", colnames(host.genus.samp))
# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp), max(tick.abr.new$Temp), by = 0.1)
# Prepare a list to store data frames
for(k in 1:n.host){
  for (j in 1:n.genus) {
    # Retrieve parameter vectors for species j
    a_samps <- host.genus.samp[, a_host.genus[k]]
    Topt_samps <- host.genus.samp[, Topt_host.genus[j]]
    c_samps <- host.genus.samp[, c_host.genus[j]]
    # For each temperature in temp_seq, evaluate across ALL posterior samples
    host.genus_eval <- sapply(temp_seq, function(temp) {
      a_samps * (Topt_samps - temp)^2 + c_samps
    })
    # Compute pointwise posterior means and HDIs
    mean_curve <- colMeans(host.genus_eval) # length(temp_seq)
    host.genus_vector <- c(host.genus_vector, mean_curve)
  }}



# Data frame with temperature, fitted value, and host group
fitted_data.host.genus <- data.frame(Temperature = rep(temp_seq, n.host * n.genus), 
                                     Fitted = host.genus_vector,
                                     Host = factor(rep(c("One", "Three", "Two"), each = length(temp_seq) * n.genus)),
                                     Genus = factor(rep(rep(c("Amblyomma", "Dermacentor", 
                                                          "Haemaphysalis", "Hyalomma"),
                                                        each = length(temp_seq)), n.host)))
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
                                     Fitted = host.genus_vector,
                                     Host = factor(rep(c("One", "Three", "Two"), each = length(temp_seq)* n.genus)),
                                     Genus = factor(rep(rep(c("Amblyomma", "Dermacentor", 
                                                          "Haemaphysalis", "Hyalomma"),
                                                        each = length(temp_seq)), n.host))) %>%
  left_join(species_to_host.genus, by = c("Host", "Genus"))|> 
  filter(!is.na(Species)) 
plot_data.host.genus$Host <- factor(plot_data.host.genus$Host, levels = c("One", "Two", "Three"))
fitted_data.host.genus$Host <- factor(fitted_data.host.genus$Host, levels = c("One", "Two", "Three"))

# Combine genus and Size for coloring
fitted_data.host.genus$`Host, Genus` <- with(fitted_data.host.genus, paste(Host, Genus, sep = ", "))
fitted_data.host.genus$`Host, Genus` <- factor(fitted_data.host.genus$`Host, Genus`, levels = c("One, Dermacentor", 
                                                                                    "Two, Hyalomma", "Three, Amblyomma",
                                                                                    "Three, Dermacentor", "Three, Haemaphysalis",
                                                                                    "Three, Hyalomma"))

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.host.genus, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.genus, aes(x = Temperature, y = Fitted, color = `Host, Genus`), size = 1) +
  labs(x = "Temperature", y = "Oviposition Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  ) 
ggsave("tick.host.genus.plot.png", width = 11.5, height = 6)

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


save(mse_by_host.genus,
     fitted_data.host.genus,
     merged_data.host.genus,
     host.genus.fit.mcmc,
     host.genus.fit,
     log.lik.host.genus,
     file="host.genus.fit.RData")

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

# Get rid of log liks
# genus.climate.mcmc <- lapply(as.mcmc.list(genus.climate.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(genus.climate.mcmc), autoburnin = FALSE)


# Extract log-likelihood matrix (samples × observations)
log.lik.genus.climate <- genus.climate.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.genus.climate <- waic(log.lik.genus.climate)$estimates["waic", "Estimate"]

# Get posterior samples
genus.climate.samp <- as.matrix(genus.climate.fit.mcmc)
genus.climate_vector <- numeric()
# Get column indices
a_genus.climate <- grep("^mu\\.climate\\.a\\[", colnames(genus.climate.samp))
Topt_genus.climate <- grep("^mu\\.genus\\.Topt\\[", colnames(genus.climate.samp))
c_genus.climate <- grep("^mu\\.genus\\.c\\[", colnames(genus.climate.samp))
# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp), max(tick.abr.new$Temp), by = 0.1)
# Prepare a list to store data frames
for(k in 1:n.climate){
  for (j in 1:n.genus) {
    # Retrieve parameter vectors for species j
    a_samps <- genus.climate.samp[, a_genus.climate[k]]
    Topt_samps <- genus.climate.samp[, Topt_genus.climate[j]]
    c_samps <- genus.climate.samp[, c_genus.climate[j]]
    # For each temperature in temp_seq, evaluate across ALL posterior samples
    genus.climate_eval <- sapply(temp_seq, function(temp) {
      a_samps * (Topt_samps - temp)^2 + c_samps
    })
    # Compute pointwise posterior means and HDIs
    mean_curve <- colMeans(genus.climate_eval) # length(temp_seq)
    genus.climate_vector <- c(genus.climate_vector, mean_curve)
  }}

# Data frame with temperature, fitted value, and genus group
fitted_data.genus.climate <- data.frame(Temperature = rep(temp_seq, n.genus * n.climate), 
                                       Fitted = genus.climate_vector,
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
                                       Fitted = genus.climate_vector,
                                       Genus = factor(rep(rep(c("Amblyomma", "Dermacentor", 
                                                                "Haemaphysalis", "Hyalomma"), each = length(temp_seq)), n.climate)),
                                       Climate = factor(rep(c("Mixed", "Subtropical", "Temperate", "Tropical"),
                                                            each = length(temp_seq) * n.genus))) %>%
  left_join(species_to_genus.climate, by = c("Genus", "Climate"))|> 
  filter(!is.na(Species)) 


# Combine Climate and Size for coloring
fitted_data.genus.climate$`Genus, Climate` <- with(fitted_data.genus.climate, paste(Genus, Climate, sep = ", "))


# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.genus.climate, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.genus.climate, aes(x = Temperature, y = Fitted, color = `Genus, Climate`), size = 1) +
  labs(x = "Temperature", y = "Ovipostion Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )  
ggsave("tick.genus.climate.plot.png", width = 11.75, height = 6.25)

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


save(mse_by_genus.climate,
     fitted_data.genus.climate,
     merged_data.genus.climate,
     genus.climate.fit.mcmc,
     genus.climate.fit,
     log.lik.genus.climate,
     file="genus.climate.fit.RData")

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

# # Get rid of log liks
# h.c.s.mcmc <- lapply(as.mcmc.list(h.c.s.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(h.c.s.mcmc), autoburnin = FALSE)


# Extract log-likelihood matrix (samples × observations)
log.lik.h.c.s <- h.c.s.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.h.c.s <- waic(log.lik.h.c.s)$estimates["waic", "Estimate"]


# Get posterior samples
h.c.s.samp <- as.matrix(h.c.s.fit.mcmc)
h.c.s_vector <- numeric()
# Get column indices
a_h.c.s <- grep("^mu\\.host\\.a\\[", colnames(h.c.s.samp))
Topt_h.c.s <- grep("^mu\\.climate\\.Topt\\[", colnames(h.c.s.samp))
c_h.c.s <- grep("^mu\\.size\\.c\\[", colnames(h.c.s.samp))
# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp), max(tick.abr.new$Temp), by = 0.1)
# Prepare a list to store data frames
for(h in 1:n.host){
  for(k in 1:n.climate){
    for (j in 1:n.size) {
      # Retrieve parameter vectors for species j
      a_samps <- h.c.s.samp[, a_h.c.s[h]]
      Topt_samps <- h.c.s.samp[, Topt_h.c.s[k]]
      c_samps <- h.c.s.samp[, c_h.c.s[j]]
      # For each temperature in temp_seq, evaluate across ALL posterior samples
      h.c.s_eval <- sapply(temp_seq, function(temp) {
        a_samps * (Topt_samps - temp)^2 + c_samps
      })
      # Compute pointwise posterior means and HDIs
      mean_curve <- colMeans(h.c.s_eval) # length(temp_seq)
      h.c.s_vector <- c(h.c.s_vector, mean_curve)
    }}
}

fitted_data.h.c.s <- data.frame(
  Temperature = rep(temp_seq, n.host * n.climate * n.size),
  Fitted = h.c.s_vector,
  Host = factor(rep(c("One", "Three", "Two"), each = length(temp_seq) * n.climate * n.size)),
  Climate = factor(rep(rep(c("Mixed", "Subtropical", "Temperate", "Tropical"),
                       each = length(temp_seq) * n.host), n.size)),
  Size = factor(rep(rep(rep(c("Large", "Medium", "Small"), each = length(temp_seq)), n.host), n.climate))
)

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
fitted_data.h.c.s$`Host, Climate, Size` <- with(fitted_data.h.c.s, paste(Host, Climate, Size, sep = ", "))
fitted_data.h.c.s$`Host, Climate, Size` <- factor(fitted_data.h.c.s$`Host, Climate, Size`, levels = c("One, Tropical, Medium",
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
  geom_line(data = fitted_data.h.c.s, aes(x = Temperature, y = Fitted, color = `Host, Climate, Size`), size = 1) +
  labs(x = "Temperature", y = "Oviposition Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  ) 
ggsave("tick.h.c.s.plot.png", width = 12, height = 7)

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


save(mse_by_h.c.s,
     fitted_data.h.c.s,
     merged_data.h.c.s,
     h.c.s.fit.mcmc,
     h.c.s.fit,
     log.lik.h.c.s,
     file="host.climate.size.fit.RData")

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

# # Get rid of log liks
# h.c.g.mcmc <- lapply(as.mcmc.list(h.c.g.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(h.c.g.mcmc), autoburnin = FALSE)


# Extract log-likelihood matrix (samples × observations)
log.lik.h.c.g <- h.c.g.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
waic.h.c.g <- waic(log.lik.h.c.g)$estimates["waic", "Estimate"]

# Get posterior samples
h.c.g.samp <- as.matrix(h.c.g.fit.mcmc)
h.c.g_vector <- numeric()
# Get column indices
a_h.c.g <- grep("^mu\\.host\\.a\\[", colnames(h.c.g.samp))
Topt_h.c.g <- grep("^mu\\.climate\\.Topt\\[", colnames(h.c.g.samp))
c_h.c.g <- grep("^mu\\.genus\\.c\\[", colnames(h.c.g.samp))
# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp), max(tick.abr.new$Temp), by = 0.1)
# Prepare a list to store data frames
for(h in 1:n.host){
  for(k in 1:n.climate){
    for (j in 1:n.genus) {
      # Retrieve parameter vectors for species j
      a_samps <- h.c.g.samp[, a_h.c.g[h]]
      Topt_samps <- h.c.g.samp[, Topt_h.c.g[k]]
      c_samps <- h.c.g.samp[, c_h.c.g[j]]
      # For each temperature in temp_seq, evaluate across ALL posterior samples
      h.c.g_eval <- sapply(temp_seq, function(temp) {
        a_samps * (Topt_samps - temp)^2 + c_samps
      })
      # Compute pointwise posterior means and HDIs
      mean_curve <- colMeans(h.c.g_eval) # length(temp_seq)
      h.c.g_vector <- c(h.c.g_vector, mean_curve)
    }}
}

fitted_data.h.c.g <- data.frame(
  Temperature = rep(temp_seq, times = n.host * n.climate * n.genus),
  Fitted = h.c.g_vector,
  Host = factor(rep(c("One", "Three", "Two"),
                    each = length(temp_seq) * n.climate * n.genus)),
  Climate = factor(rep(rep(c("Mixed", "Subtropical", "Temperate", "Tropical"),
                           each = length(temp_seq) * n.genus),
                       times = n.host)),
  Genus = factor(rep(rep(c("Amblyomma", "Dermacentor", "Haemaphysalis", "Hyalomma"),
                         each = length(temp_seq)),
                     times = n.climate * n.host))
)

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
plot_data.h.c.g$Host <- factor(plot_data.h.c.g$Host, levels = c("One", "Two", "Three"))
fitted_data.h.c.g$Host <- factor(fitted_data.h.c.g$Host, levels = c("One", "Two", "Three"))

# Combine Host, Climate, and Size for coloring
fitted_data.h.c.g$`Host, Climate, Genus` <- with(fitted_data.h.c.g, paste(Host, Climate, Genus, sep = ", "))

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.h.c.g, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.h.c.g, aes(x = Temperature, y = Fitted, color = `Host, Climate, Genus`), size = 1) +
  labs(x = "Temperature", y = "Oviposition Time in Days") +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )  
ggsave("tick.h.c.g.plot.png", width = 9, height = 6)

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


save(mse_by_h.c.g,
     fitted_data.h.c.g,
     merged_data.h.c.g,
     h.c.g.fit.mcmc,
     h.c.g.fit,
     log.lik.h.c.g,
     file="host.climate.genus.fit.RData")

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
  Method = factor(rep(c("Aggregated Fit", "By Species Individual", "By Species Hierarchical", 
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
  Method = factor(c("Aggregated Fit", "By Species Individual", "By Species Hierarchical", 
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
table.tick.avg$`Delta DIC` <- table.tick.avg$DIC - DIC.species.individual

method_levels <- table.tick.avg$Method[order(-table.tick.avg$DIC)]
table.tick.avg$Method <- factor(table.tick.avg$Method, levels = method_levels)

ggplot(table.tick.avg, aes(x = `Delta DIC`, y = Method, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  labs(x = expression(Delta~DIC), y = "Method") +
  scale_fill_viridis_d(option = "mako") +
  theme_minimal() +
  #coord_cartesian(xlim = c(9900, 10700)) +
  theme(
    axis.title = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.position = "none"
  )
ggsave("tick.DIC.png", width = 9, height = 7)


ggplot(table.tick.avg, aes(x = Total.RMSE, y = Method, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  labs(x = "Root Mean Square Error (RMSE)", y = "Method") +
  scale_fill_viridis_d(option = "mako") +
  theme_minimal() +
  coord_cartesian(xlim = c(6, 12)) +
  theme(
    axis.title = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.position = "none"
  )

ggsave("tick.RMSE.png", width = 9, height = 7)


table.tick.avg$Method <- reorder(table.tick.avg$Method, -table.tick.avg$wAIC)
ggplot(table.tick.avg, aes(x = Method, y = wAIC, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  labs(title = "wAIC by Method", x = "Method",
       y = "Widely Applicable Information Criterion (wAIC)") +
  theme_minimal() +
  coord_cartesian(ylim = c(9900, 11000)) +
  scale_fill_viridis_d(option = "mako") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
ggsave("tick.wAIC.png", width = 8.7, height = 5.5)


save(table.tick,
     table.tick.avg, file = "Dataset.Tick.RData")



estimate.topt.df <- data.frame(
  Model = factor(rep(c("Aggregated Fit", "By Species Individual", "By Species Hierarchical", 
                   "By Species Hierarchical Different Variances", 
                   "By Genus", 
                   "By Climate", "By Host", "By Size", "By Continents",
                   "By Host and Climate", "By Host and Genus", 
                   "By Genus and Climate", "By Host, Climate, and Size",
                   "By Host, Climate, and Genus"), each = 2)),
  Species = as.factor(rep(c("Amblyomma lepidum", "Dermacentor andersoni"), times = 14)),
  Lower.HDI = c(hdi(full.data.fit$BUGSoutput$sims.list$Topt)[1], 
                hdi(full.data.fit$BUGSoutput$sims.list$Topt)[1],
                hdi(a.lepidum.fit$BUGSoutput$sims.list$Topt)[1],
                hdi(d.andersoni.fit$BUGSoutput$sims.list$Topt)[1],
                hdi(species.fit$BUGSoutput$sims.list$mu.species.Topt[, 1])[1], 
                hdi(species.fit$BUGSoutput$sims.list$mu.species.Topt[, 2])[1],
                hdi(species.tau.fit$BUGSoutput$sims.list$mu.species.Topt[, 1])[1],
                hdi(species.tau.fit$BUGSoutput$sims.list$mu.species.Topt[, 2])[1],
                hdi(genus.fit$BUGSoutput$sims.list$mu.genus.Topt[, 1])[1],
                hdi(genus.fit$BUGSoutput$sims.list$mu.genus.Topt[, 2])[1],
                hdi(climate.fit$BUGSoutput$sims.list$mu.climate.Topt[, 4])[1],
                hdi(climate.fit$BUGSoutput$sims.list$mu.climate.Topt[, 3])[1],
                rep(hdi(host.fit$BUGSoutput$sims.list$mu.host.Topt[, 2])[1], 2),
                rep(hdi(size.fit$BUGSoutput$sims.list$mu.size.Topt[, 1])[1], 2),
                rep(hdi(continent.fit$BUGSoutput$sims.list$mu.continent.Topt[, 2])[1], 2),
                hdi(host.climate.fit$BUGSoutput$sims.list$mu.climate.Topt[, 4])[1],
                hdi(host.climate.fit$BUGSoutput$sims.list$mu.climate.Topt[, 3])[1],
                hdi(host.genus.fit$BUGSoutput$sims.list$mu.genus.Topt[, 1])[1],
                hdi(host.genus.fit$BUGSoutput$sims.list$mu.genus.Topt[, 2])[1],
                hdi(genus.climate.fit$BUGSoutput$sims.list$mu.genus.Topt[, 1])[1],
                hdi(genus.climate.fit$BUGSoutput$sims.list$mu.genus.Topt[, 2])[1],
                hdi(h.c.s.fit$BUGSoutput$sims.list$mu.climate.Topt[, 4])[1],
                hdi(h.c.s.fit$BUGSoutput$sims.list$mu.climate.Topt[, 3])[1],
                hdi(h.c.g.fit$BUGSoutput$sims.list$mu.climate.Topt[, 4])[1],
                hdi(h.c.g.fit$BUGSoutput$sims.list$mu.climate.Topt[, 3])[1]),
  Upper.HDI = c(hdi(full.data.fit$BUGSoutput$sims.list$Topt)[2], 
                hdi(full.data.fit$BUGSoutput$sims.list$Topt)[2],
                hdi(a.lepidum.fit$BUGSoutput$sims.list$Topt)[2],
                hdi(d.andersoni.fit$BUGSoutput$sims.list$Topt)[2],
                hdi(species.fit$BUGSoutput$sims.list$mu.species.Topt[, 1])[2], 
                hdi(species.fit$BUGSoutput$sims.list$mu.species.Topt[, 2])[2],
                hdi(species.tau.fit$BUGSoutput$sims.list$mu.species.Topt[, 1])[2],
                hdi(species.tau.fit$BUGSoutput$sims.list$mu.species.Topt[, 2])[2],
                hdi(genus.fit$BUGSoutput$sims.list$mu.genus.Topt[, 1])[2],
                hdi(genus.fit$BUGSoutput$sims.list$mu.genus.Topt[, 2])[2],
                hdi(climate.fit$BUGSoutput$sims.list$mu.climate.Topt[, 4])[2],
                hdi(climate.fit$BUGSoutput$sims.list$mu.climate.Topt[, 3])[2],
                rep(hdi(host.fit$BUGSoutput$sims.list$mu.host.Topt[, 2])[2], 2),
                rep(hdi(size.fit$BUGSoutput$sims.list$mu.size.Topt[, 1])[2], 2),
                rep(hdi(continent.fit$BUGSoutput$sims.list$mu.continent.Topt[, 2])[2], 2),
                hdi(host.climate.fit$BUGSoutput$sims.list$mu.climate.Topt[, 4])[2],
                hdi(host.climate.fit$BUGSoutput$sims.list$mu.climate.Topt[, 3])[2],
                hdi(host.genus.fit$BUGSoutput$sims.list$mu.genus.Topt[, 1])[2],
                hdi(host.genus.fit$BUGSoutput$sims.list$mu.genus.Topt[, 2])[2],
                hdi(genus.climate.fit$BUGSoutput$sims.list$mu.genus.Topt[, 1])[2],
                hdi(genus.climate.fit$BUGSoutput$sims.list$mu.genus.Topt[, 2])[2],
                hdi(h.c.s.fit$BUGSoutput$sims.list$mu.climate.Topt[, 4])[2],
                hdi(h.c.s.fit$BUGSoutput$sims.list$mu.climate.Topt[, 3])[2],
                hdi(h.c.g.fit$BUGSoutput$sims.list$mu.climate.Topt[, 4])[2],
                hdi(h.c.g.fit$BUGSoutput$sims.list$mu.climate.Topt[, 3])[2]),
  Median = c(median(full.data.fit$BUGSoutput$sims.list$Topt), 
             median(full.data.fit$BUGSoutput$sims.list$Topt),
             median(a.lepidum.fit$BUGSoutput$sims.list$Topt),
             median(d.andersoni.fit$BUGSoutput$sims.list$Topt),
             median(species.fit$BUGSoutput$sims.list$mu.species.Topt[, 1]), 
             median(species.fit$BUGSoutput$sims.list$mu.species.Topt[, 2]),
             median(species.tau.fit$BUGSoutput$sims.list$mu.species.Topt[, 1]),
             median(species.tau.fit$BUGSoutput$sims.list$mu.species.Topt[, 2]),
             median(genus.fit$BUGSoutput$sims.list$mu.genus.Topt[, 1]),
             median(genus.fit$BUGSoutput$sims.list$mu.genus.Topt[, 2]),
             median(climate.fit$BUGSoutput$sims.list$mu.climate.Topt[, 4]),
             median(climate.fit$BUGSoutput$sims.list$mu.climate.Topt[, 3]),
             rep(median(host.fit$BUGSoutput$sims.list$mu.host.Topt[, 2]), 2),
             rep(median(size.fit$BUGSoutput$sims.list$mu.size.Topt[, 1]), 2),
             rep(median(continent.fit$BUGSoutput$sims.list$mu.continent.Topt[, 2]), 2),
             median(host.climate.fit$BUGSoutput$sims.list$mu.climate.Topt[, 4]),
             median(host.climate.fit$BUGSoutput$sims.list$mu.climate.Topt[, 3]),
             median(host.genus.fit$BUGSoutput$sims.list$mu.genus.Topt[, 1]),
             median(host.genus.fit$BUGSoutput$sims.list$mu.genus.Topt[, 2]),
             median(genus.climate.fit$BUGSoutput$sims.list$mu.genus.Topt[, 1]),
             median(genus.climate.fit$BUGSoutput$sims.list$mu.genus.Topt[, 2]),
             median(h.c.s.fit$BUGSoutput$sims.list$mu.climate.Topt[, 4]),
             median(h.c.s.fit$BUGSoutput$sims.list$mu.climate.Topt[, 3]),
             median(h.c.g.fit$BUGSoutput$sims.list$mu.climate.Topt[, 4]),
             median(h.c.g.fit$BUGSoutput$sims.list$mu.climate.Topt[, 3])))

estimate.topt.df$Highlight <- ifelse(estimate.topt.df$Model == "By Species Individual", "Highlight", "Normal")
estimate.topt.df$Model <- factor(
  estimate.topt.df$Model,
  levels = rev(unique(estimate.topt.df$Model))  # keeps the order as-is
)
# Extract HDI bounds for "By Species Individual"
hdi_bounds <- estimate.topt.df %>%
  filter(Model == "By Species Individual") %>%
  select(Species, Lower.HDI, Upper.HDI)

ggplot(estimate.topt.df, aes(x = Model, y = Median, color = (Model == "By Species Individual"))) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = Lower.HDI, ymax = Upper.HDI), 
                position = position_dodge(width = 0.5), width = 0.2) +
  geom_hline(data = hdi_bounds, aes(yintercept = Lower.HDI), linetype = "dashed", color = "red") +
  geom_hline(data = hdi_bounds, aes(yintercept = Upper.HDI), linetype = "dashed", color = "red") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  labs(
    x = "Model",
    y = "Optimal Temperature"
  ) +
  coord_flip() +
  facet_wrap(~ Species) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 16),
    legend.position = "none",
    panel.spacing = unit(4, "lines")
  )
ggsave("Topt.Ovi.png", height = 8, width = 14)

estimate.c.df <- data.frame(
  Model = factor(rep(c("Aggregated Fit", "By Species Individual", "By Species Hierarchical", 
                       "By Species Hierarchical Different Variances", 
                       "By Genus", 
                       "By Climate", "By Host", "By Size", "By Continents",
                       "By Host and Climate", "By Host and Genus", 
                       "By Genus and Climate", "By Host, Climate, and Size",
                       "By Host, Climate, and Genus"), each = 2)),
  Species = as.factor(rep(c("Amblyomma lepidum", "Dermacentor andersoni"), times = 14)),
  Lower.HDI = c(hdi(full.data.fit$BUGSoutput$sims.list$c)[1], 
                hdi(full.data.fit$BUGSoutput$sims.list$c)[1],
                hdi(a.lepidum.fit$BUGSoutput$sims.list$c)[1],
                hdi(d.andersoni.fit$BUGSoutput$sims.list$c)[1],
                hdi(species.fit$BUGSoutput$sims.list$mu.species.c[, 1])[1], 
                hdi(species.fit$BUGSoutput$sims.list$mu.species.c[, 2])[1],
                hdi(species.tau.fit$BUGSoutput$sims.list$mu.species.c[, 1])[1],
                hdi(species.tau.fit$BUGSoutput$sims.list$mu.species.c[, 2])[1],
                hdi(genus.fit$BUGSoutput$sims.list$mu.genus.c[, 1])[1],
                hdi(genus.fit$BUGSoutput$sims.list$mu.genus.c[, 2])[1],
                hdi(climate.fit$BUGSoutput$sims.list$mu.climate.c[, 4])[1],
                hdi(climate.fit$BUGSoutput$sims.list$mu.climate.c[, 3])[1],
                rep(hdi(host.fit$BUGSoutput$sims.list$mu.host.c[, 2])[1], 2),
                rep(hdi(size.fit$BUGSoutput$sims.list$mu.size.c[, 1])[1], 2),
                rep(hdi(continent.fit$BUGSoutput$sims.list$mu.continent.c[, 2])[1], 2),
                hdi(host.climate.fit$BUGSoutput$sims.list$mu.host.c[, 2])[1],
                hdi(host.climate.fit$BUGSoutput$sims.list$mu.host.c[, 2])[1],
                hdi(host.genus.fit$BUGSoutput$sims.list$mu.genus.c[, 1])[1],
                hdi(host.genus.fit$BUGSoutput$sims.list$mu.genus.c[, 2])[1],
                hdi(genus.climate.fit$BUGSoutput$sims.list$mu.genus.c[, 1])[1],
                hdi(genus.climate.fit$BUGSoutput$sims.list$mu.genus.c[, 2])[1],
                hdi(h.c.s.fit$BUGSoutput$sims.list$mu.size.c[, 1])[1],
                hdi(h.c.s.fit$BUGSoutput$sims.list$mu.size.c[, 1])[1],
                hdi(h.c.g.fit$BUGSoutput$sims.list$mu.genus.c[, 1])[1],
                hdi(h.c.g.fit$BUGSoutput$sims.list$mu.genus.c[, 2])[1]),
  Upper.HDI = c(hdi(full.data.fit$BUGSoutput$sims.list$c)[2], 
                hdi(full.data.fit$BUGSoutput$sims.list$c)[2],
                hdi(a.lepidum.fit$BUGSoutput$sims.list$c)[2],
                hdi(d.andersoni.fit$BUGSoutput$sims.list$c)[2],
                hdi(species.fit$BUGSoutput$sims.list$mu.species.c[, 1])[2], 
                hdi(species.fit$BUGSoutput$sims.list$mu.species.c[, 2])[2],
                hdi(species.tau.fit$BUGSoutput$sims.list$mu.species.c[, 1])[2],
                hdi(species.tau.fit$BUGSoutput$sims.list$mu.species.c[, 2])[2],
                hdi(genus.fit$BUGSoutput$sims.list$mu.genus.c[, 1])[2],
                hdi(genus.fit$BUGSoutput$sims.list$mu.genus.c[, 2])[2],
                hdi(climate.fit$BUGSoutput$sims.list$mu.climate.c[, 4])[2],
                hdi(climate.fit$BUGSoutput$sims.list$mu.climate.c[, 3])[2],
                rep(hdi(host.fit$BUGSoutput$sims.list$mu.host.c[, 2])[2], 2),
                rep(hdi(size.fit$BUGSoutput$sims.list$mu.size.c[, 1])[2], 2),
                rep(hdi(continent.fit$BUGSoutput$sims.list$mu.continent.c[, 2])[2], 2),
                hdi(host.climate.fit$BUGSoutput$sims.list$mu.host.c[, 2])[2],
                hdi(host.climate.fit$BUGSoutput$sims.list$mu.host.c[, 2])[2],
                hdi(host.genus.fit$BUGSoutput$sims.list$mu.genus.c[, 1])[2],
                hdi(host.genus.fit$BUGSoutput$sims.list$mu.genus.c[, 2])[2],
                hdi(genus.climate.fit$BUGSoutput$sims.list$mu.genus.c[, 1])[2],
                hdi(genus.climate.fit$BUGSoutput$sims.list$mu.genus.c[, 2])[2],
                hdi(h.c.s.fit$BUGSoutput$sims.list$mu.size.c[, 1])[2],
                hdi(h.c.s.fit$BUGSoutput$sims.list$mu.size.c[, 1])[2],
                hdi(h.c.g.fit$BUGSoutput$sims.list$mu.genus.c[, 1])[2],
                hdi(h.c.g.fit$BUGSoutput$sims.list$mu.genus.c[, 2])[2]),
  Median = c(median(full.data.fit$BUGSoutput$sims.list$c), 
             median(full.data.fit$BUGSoutput$sims.list$c),
             median(a.lepidum.fit$BUGSoutput$sims.list$c),
             median(d.andersoni.fit$BUGSoutput$sims.list$c),
             median(species.fit$BUGSoutput$sims.list$mu.species.c[, 1]), 
             median(species.fit$BUGSoutput$sims.list$mu.species.c[, 2]),
             median(species.tau.fit$BUGSoutput$sims.list$mu.species.c[, 1]),
             median(species.tau.fit$BUGSoutput$sims.list$mu.species.c[, 2]),
             median(genus.fit$BUGSoutput$sims.list$mu.genus.c[, 1]),
             median(genus.fit$BUGSoutput$sims.list$mu.genus.c[, 2]),
             median(climate.fit$BUGSoutput$sims.list$mu.climate.c[, 4]),
             median(climate.fit$BUGSoutput$sims.list$mu.climate.c[, 3]),
             rep(median(host.fit$BUGSoutput$sims.list$mu.host.c[, 2]), 2),
             rep(median(size.fit$BUGSoutput$sims.list$mu.size.c[, 1]), 2),
             rep(median(continent.fit$BUGSoutput$sims.list$mu.continent.c[, 2]), 2),
             median(host.climate.fit$BUGSoutput$sims.list$mu.host.c[, 2]),
             median(host.climate.fit$BUGSoutput$sims.list$mu.host.c[, 2]),
             median(host.genus.fit$BUGSoutput$sims.list$mu.genus.c[, 1]),
             median(host.genus.fit$BUGSoutput$sims.list$mu.genus.c[, 2]),
             median(genus.climate.fit$BUGSoutput$sims.list$mu.genus.c[, 1]),
             median(genus.climate.fit$BUGSoutput$sims.list$mu.genus.c[, 2]),
             median(h.c.s.fit$BUGSoutput$sims.list$mu.size.c[, 1]),
             median(h.c.s.fit$BUGSoutput$sims.list$mu.size.c[, 1]),
             median(h.c.g.fit$BUGSoutput$sims.list$mu.genus.c[, 1]),
             median(h.c.g.fit$BUGSoutput$sims.list$mu.genus.c[, 2])))

estimate.c.df$Highlight <- ifelse(estimate.c.df$Model == "By Species Individual", "Highlight", "Normal")
estimate.c.df$Model <- factor(
  estimate.c.df$Model,
  levels = rev(unique(estimate.c.df$Model))  # keeps the order as-is
)
# Extract HDI bounds for "By Species Individual"
hdi_bounds <- estimate.c.df %>%
  filter(Model == "By Species Individual") %>%
  select(Species, Lower.HDI, Upper.HDI)

ggplot(estimate.c.df, aes(x = Model, y = Median, color = (Model == "By Species Individual"))) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = Lower.HDI, ymax = Upper.HDI), 
                position = position_dodge(width = 0.5), width = 0.2) +
  geom_hline(data = hdi_bounds, aes(yintercept = Lower.HDI), linetype = "dashed", color = "red") +
  geom_hline(data = hdi_bounds, aes(yintercept = Upper.HDI), linetype = "dashed", color = "red") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  labs(
    x = "Model",
    y = "Minimum Oviposition"
  ) +
  coord_flip() +
  facet_wrap(~ Species) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 16),
    legend.position = "none",
    panel.spacing = unit(4, "lines")
  )
ggsave("Min.Ovi.png", height = 8, width = 14)

save(estimate.c.df, estimate.topt.df, file = "hdi.plot.tick.RData")

estimate.topt.rm.df <- data.frame(
  Model = factor(rep(c("By Genus", 
                       "By Host", 
                       "By Host and Climate",
                       "By Host, Climate, and Genus"), each = 8)),
  Species = as.factor(rep(c("Hyalomma aegyptium", "Hyalomma dromedarii", 
                            "Hyalomma lusitanicum", "Hyalomma schulzei"), times = 8)),
  Type = as.factor(rep(rep(c("Full", "Removed"), each = 4), 4)),
  Lower.HDI = c(rep(hdi(genus.fit$BUGSoutput$sims.list$mu.genus.Topt[, 4])[1], 4),
                hdi(genus.aegyptium.fit$BUGSoutput$sims.list$mu.genus.Topt[, 4])[1],
                hdi(genus.dromedarii.fit$BUGSoutput$sims.list$mu.genus.Topt[, 4])[1],
                hdi(genus.lusitanicum.fit$BUGSoutput$sims.list$mu.genus.Topt[, 4])[1],
                hdi(genus.schulzei.fit$BUGSoutput$sims.list$mu.genus.Topt[, 4])[1],
                rep(c(hdi(host.fit$BUGSoutput$sims.list$mu.host.Topt[, 2])[1],
                      hdi(host.fit$BUGSoutput$sims.list$mu.host.Topt[, 3])[1]), times = 2),
                hdi(host.aegyptium.fit$BUGSoutput$sims.list$mu.host.Topt[, 2])[1],
                hdi(host.dromedarii.fit$BUGSoutput$sims.list$mu.host.Topt[, 3])[1],
                hdi(host.lusitanicum.fit$BUGSoutput$sims.list$mu.host.Topt[, 2])[1],
                hdi(host.schulzei.fit$BUGSoutput$sims.list$mu.host.Topt[, 3])[1],
                rep(hdi(host.climate.fit$BUGSoutput$sims.list$mu.climate.Topt[, 2])[1], 4),
                hdi(host.climate.aegyptium.fit$BUGSoutput$sims.list$mu.climate.Topt[, 2])[1],
                hdi(host.climate.dromedarii.fit$BUGSoutput$sims.list$mu.climate.Topt[, 2])[1],
                hdi(host.climate.lusitanicum.fit$BUGSoutput$sims.list$mu.climate.Topt[, 2])[1],
                hdi(host.climate.schulzei.fit$BUGSoutput$sims.list$mu.climate.Topt[, 2])[1],
                rep(hdi(h.c.g.fit$BUGSoutput$sims.list$mu.climate.Topt[, 2])[1], 4),
                hdi(h.c.g.aegyptium.fit$BUGSoutput$sims.list$mu.climate.Topt[, 2])[1],
                hdi(h.c.g.dromedarii.fit$BUGSoutput$sims.list$mu.climate.Topt[, 2])[1],
                hdi(h.c.g.lusitanicum.fit$BUGSoutput$sims.list$mu.climate.Topt[, 2])[1],
                hdi(h.c.g.schulzei.fit$BUGSoutput$sims.list$mu.climate.Topt[, 2])[1]),
  Upper.HDI = c(rep(hdi(genus.fit$BUGSoutput$sims.list$mu.genus.Topt[, 4])[2], 4),
                hdi(genus.aegyptium.fit$BUGSoutput$sims.list$mu.genus.Topt[, 4])[2],
                hdi(genus.dromedarii.fit$BUGSoutput$sims.list$mu.genus.Topt[, 4])[2],
                hdi(genus.lusitanicum.fit$BUGSoutput$sims.list$mu.genus.Topt[, 4])[2],
                hdi(genus.schulzei.fit$BUGSoutput$sims.list$mu.genus.Topt[, 4])[2],
                rep(c(hdi(host.fit$BUGSoutput$sims.list$mu.host.Topt[, 2])[2],
                      hdi(host.fit$BUGSoutput$sims.list$mu.host.Topt[, 3])[2]), times = 2),
                hdi(host.aegyptium.fit$BUGSoutput$sims.list$mu.host.Topt[, 2])[2],
                hdi(host.dromedarii.fit$BUGSoutput$sims.list$mu.host.Topt[, 3])[2],
                hdi(host.lusitanicum.fit$BUGSoutput$sims.list$mu.host.Topt[, 2])[2],
                hdi(host.schulzei.fit$BUGSoutput$sims.list$mu.host.Topt[, 3])[2],
                rep(hdi(host.climate.fit$BUGSoutput$sims.list$mu.climate.Topt[, 2])[2], 4),
                hdi(host.climate.aegyptium.fit$BUGSoutput$sims.list$mu.climate.Topt[, 2])[2],
                hdi(host.climate.dromedarii.fit$BUGSoutput$sims.list$mu.climate.Topt[, 2])[2],
                hdi(host.climate.lusitanicum.fit$BUGSoutput$sims.list$mu.climate.Topt[, 2])[2],
                hdi(host.climate.schulzei.fit$BUGSoutput$sims.list$mu.climate.Topt[, 2])[2],
                rep(hdi(h.c.g.fit$BUGSoutput$sims.list$mu.climate.Topt[, 2])[2], 4),
                hdi(h.c.g.aegyptium.fit$BUGSoutput$sims.list$mu.climate.Topt[, 2])[2],
                hdi(h.c.g.dromedarii.fit$BUGSoutput$sims.list$mu.climate.Topt[, 2])[2],
                hdi(h.c.g.lusitanicum.fit$BUGSoutput$sims.list$mu.climate.Topt[, 2])[2],
                hdi(h.c.g.schulzei.fit$BUGSoutput$sims.list$mu.climate.Topt[, 2])[2]),
  Median = c(rep(median(genus.fit$BUGSoutput$sims.list$mu.genus.Topt[, 4]), 4),
             median(genus.aegyptium.fit$BUGSoutput$sims.list$mu.genus.Topt[, 4]),
             median(genus.dromedarii.fit$BUGSoutput$sims.list$mu.genus.Topt[, 4]),
             median(genus.lusitanicum.fit$BUGSoutput$sims.list$mu.genus.Topt[, 4]),
             median(genus.schulzei.fit$BUGSoutput$sims.list$mu.genus.Topt[, 4]),
             rep(c(median(host.fit$BUGSoutput$sims.list$mu.host.Topt[, 2]),
                   median(host.fit$BUGSoutput$sims.list$mu.host.Topt[, 3])), times = 2),
             median(host.aegyptium.fit$BUGSoutput$sims.list$mu.host.Topt[, 2]),
             median(host.dromedarii.fit$BUGSoutput$sims.list$mu.host.Topt[, 3]),
             median(host.lusitanicum.fit$BUGSoutput$sims.list$mu.host.Topt[, 2]),
             median(host.schulzei.fit$BUGSoutput$sims.list$mu.host.Topt[, 3]),
             rep(median(host.climate.fit$BUGSoutput$sims.list$mu.climate.Topt[, 2]), 4),
             median(host.climate.aegyptium.fit$BUGSoutput$sims.list$mu.climate.Topt[, 2]),
             median(host.climate.dromedarii.fit$BUGSoutput$sims.list$mu.climate.Topt[, 2]),
             median(host.climate.lusitanicum.fit$BUGSoutput$sims.list$mu.climate.Topt[, 2]),
             median(host.climate.schulzei.fit$BUGSoutput$sims.list$mu.climate.Topt[, 2]),
             rep(median(h.c.g.fit$BUGSoutput$sims.list$mu.climate.Topt[, 2]), 4),
             median(h.c.g.aegyptium.fit$BUGSoutput$sims.list$mu.climate.Topt[, 2]),
             median(h.c.g.dromedarii.fit$BUGSoutput$sims.list$mu.climate.Topt[, 2]),
             median(h.c.g.lusitanicum.fit$BUGSoutput$sims.list$mu.climate.Topt[, 2]),
             median(h.c.g.schulzei.fit$BUGSoutput$sims.list$mu.climate.Topt[, 2])))

estimate.topt.rm.df$Model <- factor(
  estimate.topt.rm.df$Model,
  levels = rev(unique(estimate.topt.rm.df$Model))  # keeps the order as-is
)

ggplot(estimate.topt.rm.df, aes(x = Model, y = Median, color = Type)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = Lower.HDI, ymax = Upper.HDI), 
                position = position_dodge(width = 0.5), width = 0.2) +
  scale_color_manual(
    values = c("Removed" = "blue", "Full" = "black"),
    breaks = c("Removed", "Full")  # order in legend
  ) +
  labs(
    x = "Model",
    y = "Optimal Temperature"
  ) +
  coord_flip() +
  facet_wrap(~ Species) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.title = element_blank(),
    panel.spacing = unit(3, "lines")
  )

ggsave("Topt.Ovi.RM.png", height = 7, width = 10)


estimate.c.rm.df <- data.frame(
  Model = factor(rep(c("By Genus", 
                       "By Host", 
                       "By Host and Climate",
                       "By Host, Climate, and Genus"), each = 8)),
  Species = as.factor(rep(c("Hyalomma aegyptium", "Hyalomma dromedarii", 
                            "Hyalomma lusitanicum", "Hyalomma schulzei"), times = 8)),
  Type = as.factor(rep(rep(c("Full", "Removed"), each = 4), 4)),
  Lower.HDI = c(rep(hdi(genus.fit$BUGSoutput$sims.list$mu.genus.c[, 4])[1], 4),
                hdi(genus.aegyptium.fit$BUGSoutput$sims.list$mu.genus.c[, 4])[1],
                hdi(genus.dromedarii.fit$BUGSoutput$sims.list$mu.genus.c[, 4])[1],
                hdi(genus.lusitanicum.fit$BUGSoutput$sims.list$mu.genus.c[, 4])[1],
                hdi(genus.schulzei.fit$BUGSoutput$sims.list$mu.genus.c[, 4])[1],
                rep(c(hdi(host.fit$BUGSoutput$sims.list$mu.host.c[, 2])[1],
                      hdi(host.fit$BUGSoutput$sims.list$mu.host.c[, 3])[1]), times = 2),
                hdi(host.aegyptium.fit$BUGSoutput$sims.list$mu.host.c[, 2])[1],
                hdi(host.dromedarii.fit$BUGSoutput$sims.list$mu.host.c[, 3])[1],
                hdi(host.lusitanicum.fit$BUGSoutput$sims.list$mu.host.c[, 2])[1],
                hdi(host.schulzei.fit$BUGSoutput$sims.list$mu.host.c[, 3])[1],
                rep(c(hdi(host.climate.fit$BUGSoutput$sims.list$mu.host.c[, 2])[1],
                      hdi(host.climate.fit$BUGSoutput$sims.list$mu.host.c[, 3])[1]), times = 2),
                hdi(host.climate.aegyptium.fit$BUGSoutput$sims.list$mu.host.c[, 2])[1],
                hdi(host.climate.dromedarii.fit$BUGSoutput$sims.list$mu.host.c[, 3])[1],
                hdi(host.climate.lusitanicum.fit$BUGSoutput$sims.list$mu.host.c[, 2])[1],
                hdi(host.climate.schulzei.fit$BUGSoutput$sims.list$mu.host.c[, 3])[1],
                rep(hdi(h.c.g.fit$BUGSoutput$sims.list$mu.genus.c[, 4])[1], 4),
                hdi(h.c.g.aegyptium.fit$BUGSoutput$sims.list$mu.genus.c[, 4])[1],
                hdi(h.c.g.dromedarii.fit$BUGSoutput$sims.list$mu.genus.c[, 4])[1],
                hdi(h.c.g.lusitanicum.fit$BUGSoutput$sims.list$mu.genus.c[, 4])[1],
                hdi(h.c.g.schulzei.fit$BUGSoutput$sims.list$mu.genus.c[, 4])[1]),
  Upper.HDI = c(rep(hdi(genus.fit$BUGSoutput$sims.list$mu.genus.c[, 4])[2], 4),
                hdi(genus.aegyptium.fit$BUGSoutput$sims.list$mu.genus.c[, 4])[2],
                hdi(genus.dromedarii.fit$BUGSoutput$sims.list$mu.genus.c[, 4])[2],
                hdi(genus.lusitanicum.fit$BUGSoutput$sims.list$mu.genus.c[, 4])[2],
                hdi(genus.schulzei.fit$BUGSoutput$sims.list$mu.genus.c[, 4])[2],
                rep(c(hdi(host.fit$BUGSoutput$sims.list$mu.host.c[, 2])[2],
                      hdi(host.fit$BUGSoutput$sims.list$mu.host.c[, 3])[2]), times = 2),
                hdi(host.aegyptium.fit$BUGSoutput$sims.list$mu.host.c[, 2])[2],
                hdi(host.dromedarii.fit$BUGSoutput$sims.list$mu.host.c[, 3])[2],
                hdi(host.lusitanicum.fit$BUGSoutput$sims.list$mu.host.c[, 2])[2],
                hdi(host.schulzei.fit$BUGSoutput$sims.list$mu.host.c[, 3])[2],
                rep(c(hdi(host.climate.fit$BUGSoutput$sims.list$mu.host.c[, 2])[2],
                      hdi(host.climate.fit$BUGSoutput$sims.list$mu.host.c[, 3])[2]), times = 2),
                hdi(host.climate.aegyptium.fit$BUGSoutput$sims.list$mu.host.c[, 2])[2],
                hdi(host.climate.dromedarii.fit$BUGSoutput$sims.list$mu.host.c[, 3])[2],
                hdi(host.climate.lusitanicum.fit$BUGSoutput$sims.list$mu.host.c[, 2])[2],
                hdi(host.climate.schulzei.fit$BUGSoutput$sims.list$mu.host.c[, 3])[2],
                rep(hdi(h.c.g.fit$BUGSoutput$sims.list$mu.genus.c[, 4])[2], 4),
                hdi(h.c.g.aegyptium.fit$BUGSoutput$sims.list$mu.genus.c[, 4])[2],
                hdi(h.c.g.dromedarii.fit$BUGSoutput$sims.list$mu.genus.c[, 4])[2],
                hdi(h.c.g.lusitanicum.fit$BUGSoutput$sims.list$mu.genus.c[, 4])[2],
                hdi(h.c.g.schulzei.fit$BUGSoutput$sims.list$mu.genus.c[, 4])[2]),
  Median = c(rep(median(genus.fit$BUGSoutput$sims.list$mu.genus.c[, 4]), 4),
             median(genus.aegyptium.fit$BUGSoutput$sims.list$mu.genus.c[, 4]),
             median(genus.dromedarii.fit$BUGSoutput$sims.list$mu.genus.c[, 4]),
             median(genus.lusitanicum.fit$BUGSoutput$sims.list$mu.genus.c[, 4]),
             median(genus.schulzei.fit$BUGSoutput$sims.list$mu.genus.c[, 4]),
             rep(c(median(host.fit$BUGSoutput$sims.list$mu.host.c[, 2]),
                   median(host.fit$BUGSoutput$sims.list$mu.host.c[, 3])), times = 2),
             median(host.aegyptium.fit$BUGSoutput$sims.list$mu.host.c[, 2]),
             median(host.dromedarii.fit$BUGSoutput$sims.list$mu.host.c[, 3]),
             median(host.lusitanicum.fit$BUGSoutput$sims.list$mu.host.c[, 2]),
             median(host.schulzei.fit$BUGSoutput$sims.list$mu.host.c[, 3]),
             rep(c(median(host.climate.fit$BUGSoutput$sims.list$mu.host.c[, 2]),
                   median(host.climate.fit$BUGSoutput$sims.list$mu.host.c[, 3])), times = 2),
             median(host.climate.aegyptium.fit$BUGSoutput$sims.list$mu.host.c[, 2]),
             median(host.climate.dromedarii.fit$BUGSoutput$sims.list$mu.host.c[, 3]),
             median(host.climate.lusitanicum.fit$BUGSoutput$sims.list$mu.host.c[, 2]),
             median(host.climate.schulzei.fit$BUGSoutput$sims.list$mu.host.c[, 3]),
             rep(median(h.c.g.fit$BUGSoutput$sims.list$mu.genus.c[, 4]), 4),
             median(h.c.g.aegyptium.fit$BUGSoutput$sims.list$mu.genus.c[, 4]),
             median(h.c.g.dromedarii.fit$BUGSoutput$sims.list$mu.genus.c[, 4]),
             median(h.c.g.lusitanicum.fit$BUGSoutput$sims.list$mu.genus.c[, 4]),
             median(h.c.g.schulzei.fit$BUGSoutput$sims.list$mu.genus.c[, 4])))

estimate.c.rm.df$Model <- factor(
  estimate.c.rm.df$Model,
  levels = rev(unique(estimate.c.rm.df$Model))  # keeps the order as-is
)
ggplot(estimate.c.rm.df, aes(x = Model, y = Median, color = Type)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = Lower.HDI, ymax = Upper.HDI), 
                position = position_dodge(width = 0.5), width = 0.2) +
  scale_color_manual(
    values = c("Removed" = "blue", "Full" = "black"),
    breaks = c("Removed", "Full")  # order in legend
  ) +
  labs(
    x = "Model",
    y = "Minimum Oviposition"
  ) +
  coord_flip() +
  facet_wrap(~ Species) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.title = element_blank(),
    panel.spacing = unit(3, "lines")
  )

ggsave("Min.Ovi.RM.png", height = 7, width = 10)
save(estimate.c.rm.df, estimate.topt.rm.df, file = "hdi.removed.tick.RData")


################################### Plots ###################################
fitted_data.a.lepidum$Species <- "Amblyomma lepidum"
fitted_data.a.lepidum$Fit <- "Individual Data"
fitted_data.d.andersoni$Species <- "Dermacentor andersoni"
fitted_data.d.andersoni$Fit <- "Individual Data"
fitted_data.d.nitens$Species <- "Dermacentor nitens"
fitted_data.d.nitens$Fit <- "Individual Data"
fitted_data.h.leporispalustris$Species <- "Haemaphysalis leporispalustris"
fitted_data.h.leporispalustris$Fit <- "Individual Data"
fitted_data.h.aegyptium$Species <- "Hyalomma aegyptium"
fitted_data.h.aegyptium$Fit <- "Individual Data"
fitted_data.h.dromedarii$Species <- "Hyalomma dromedarii"
fitted_data.h.dromedarii$Fit <- "Individual Data"
fitted_data.h.impeltatum$Species <- "Hyalomma impeltatum"
fitted_data.h.impeltatum$Fit <- "Individual Data"
fitted_data.h.lusitanicum$Species <- "Hyalomma lusitanicum"
fitted_data.h.lusitanicum$Fit <- "Individual Data"
fitted_data.h.schulzei$Species <- "Hyalomma schulzei"
fitted_data.h.schulzei$Fit <- "Individual Data"
fitted_data.full.data$Fit <- "Full Data"
fitted.full.ind <- rbind(fitted_data.full.data,
                         fitted_data.a.lepidum,
                         fitted_data.d.andersoni,
                         fitted_data.d.nitens,
                         fitted_data.h.leporispalustris,
                         fitted_data.h.aegyptium,
                         fitted_data.h.dromedarii,
                         fitted_data.h.impeltatum,
                         fitted_data.h.lusitanicum,
                         fitted_data.h.schulzei)

ggplot() +
  geom_point(data = plot_data.full, aes(x = Temperature, y = Trait), color = "grey52", alpha = 0.5) +
  geom_line(data = fitted.full.ind, aes(x = Temperature, y = Fitted, color = Fit), size = 1) +
  labs(x = "Temperature", y = "Oviposition Time in Days") +
  scale_color_viridis_d(option = "C") +
  theme_minimal() +
  #theme(legend.position = "none") +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )
ggsave("tick.full.ind.plot.png", width = 11, height = 6)


fitted_data.a.lepidum$Genus <- "Individual Data"
fitted_data.d.andersoni$Genus <- "Individual Data"
fitted_data.d.nitens$Genus <- "Individual Data"
fitted_data.h.leporispalustris$Genus <- "Individual Data"
fitted_data.h.aegyptium$Genus <- "Individual Data"
fitted_data.h.dromedarii$Genus <- "Individual Data"
fitted_data.h.impeltatum$Genus <- "Individual Data"
fitted_data.h.lusitanicum$Genus <- "Individual Data"
fitted_data.h.schulzei$Genus <- "Individual Data"
fitted_data.genus$Fit <- "Genus"
fitted.genus.ind <- rbind(fitted_data.genus,
                         fitted_data.a.lepidum,
                         fitted_data.d.andersoni,
                         fitted_data.d.nitens,
                         fitted_data.h.leporispalustris,
                         fitted_data.h.aegyptium,
                         fitted_data.h.dromedarii,
                         fitted_data.h.impeltatum,
                         fitted_data.h.lusitanicum,
                         fitted_data.h.schulzei)

# Separate Individual vs. Genus lines
individual_lines <- fitted.genus.ind %>% filter(Genus == "Individual Data")
genus_lines <- fitted.genus.ind %>% filter(Genus != "Individual Data")
ggplot() +
  geom_point(data = plot_data.genus, aes(x = Temperature, y = Trait), color = "grey52", alpha = 0.5) +
  geom_line(data = individual_lines, aes(x = Temperature, y = Fitted, color = Genus), size = 1) +
  geom_line(data = genus_lines, aes(x = Temperature, y = Fitted, color = Genus), size = 1) +
  labs(x = "Temperature", y = "Oviposition Time in Days") +
  scale_color_viridis_d(option = "C") +
  theme_minimal() +
  facet_wrap(~ Species) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )
ggsave("tick.genus.ind.plot.png", width = 11, height = 6)


fitted_data.a.lepidum$Host <- "Individual Data"
fitted_data.d.andersoni$Host <- "Individual Data"
fitted_data.d.nitens$Host <- "Individual Data"
fitted_data.h.leporispalustris$Host <- "Individual Data"
fitted_data.h.aegyptium$Host <- "Individual Data"
fitted_data.h.dromedarii$Host <- "Individual Data"
fitted_data.h.impeltatum$Host <- "Individual Data"
fitted_data.h.lusitanicum$Host <- "Individual Data"
fitted_data.h.schulzei$Host <- "Individual Data"
fitted_data.host$Fit <- "Host"
fitted.host.ind <- rbind(fitted_data.host,
                          fitted_data.a.lepidum[, -5],
                          fitted_data.d.andersoni[, -5],
                          fitted_data.d.nitens[, -5],
                          fitted_data.h.leporispalustris[, -5],
                          fitted_data.h.aegyptium[, -5],
                          fitted_data.h.dromedarii[, -5],
                          fitted_data.h.impeltatum[, -5],
                          fitted_data.h.lusitanicum[, -5],
                          fitted_data.h.schulzei[, -5])

# Separate Individual vs. Genus lines
individual_lines.host <- fitted.host.ind %>% filter(Host == "Individual Data")
host_lines <- fitted.host.ind %>% filter(Host != "Individual Data")
ggplot() +
  geom_point(data = plot_data.host, aes(x = Temperature, y = Trait), color = "grey52", alpha = 0.5) +
  geom_line(data = individual_lines.host, aes(x = Temperature, y = Fitted, color = Host), size = 1) +
  geom_line(data = host_lines, aes(x = Temperature, y = Fitted, color = Host), size = 1) +
  labs(x = "Temperature", y = "Oviposition Time in Days") +
  scale_color_viridis_d(option = "C") +
  theme_minimal() +
  facet_wrap(~ Species) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )
ggsave("tick.host.ind.plot.png", width = 11, height = 6)



fitted_data.h.c.g.plot <- fitted_data.h.c.g
fitted_data.h.c.g.plot$`Host, Climate, Genus` <- "Full Fit"
fitted.hcg.ind <- rbind(fitted_data.h.c.g.plot,
                        fitted_data.h.c.g.schulzei)
fitted.hcg.ind$`Host, Climate, Genus` <- factor(
  fitted.hcg.ind$`Host, Climate, Genus`,
  levels = c("Full Fit", "One, Tropical, Dermacentor", 
             "Two, Subtropical, Hyalomma", "Three, Mixed, Haemaphysalis",
             "Three, Mixed, Hyalomma",
             "Three, Temperate, Dermacentor",
             "Three, Subtropical, Hyalomma",
             "Three, Tropical, Amblyomma")
)
ggplot() +
  geom_point(data = plot_data.h.c.g, aes(x = Temperature, y = Trait), color = "grey52", alpha = 0.5) +
  #geom_line(data = fitted.hcg.ind, aes(x = Temperature, y = Fitted, color = `Host, Climate, Genus`), linewidth = 1) +
  # Full fit line (separate legend via linetype)
  geom_line(data = subset(fitted.hcg.ind, `Host, Climate, Genus` == "Full Fit"),
            aes(x = Temperature, y = Fitted, linetype = "Full Fit"), color = "red", size = 1) +
  # Other fits (separate color legend)
  geom_line(data = subset(fitted.hcg.ind, `Host, Climate, Genus` != "Full Fit"),
            aes(x = Temperature, y = Fitted, color = `Host, Climate, Genus`), size = 1) +
  # Legends
  scale_linetype_manual(name = " ", values = c("Full Fit" = "solid")) +
  scale_color_viridis_d(name = "Host, Climate, Genus", option = "C") +
  labs(x = "Temperature", y = "Oviposition Time in Days") +
  theme_minimal() +
  #theme(legend.position = "none") +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )  
ggsave("tick.hcg.ind.plot.png", width = 13, height = 6)


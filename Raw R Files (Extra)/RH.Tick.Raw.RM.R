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
################################### By Genus- Andersoni ###################################
tick.rm.raw <- tick.raw[tick.raw$Species != "Dermacentor andersoni", ]
genus.rm.raw <- as.factor(tick.rm.raw$Genus)
n.genus.rm.raw <- length(unique(tick.rm.raw$Genus))
N.obs.rm.raw <- length(tick.rm.raw$Trait)
unique(tick.rm.raw$Species)
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

  for (j in 1:n.genus.rm.raw) {
    mu.genus.c[j] ~ dexp(mu.c)       # Genus-specific vertex y
    mu.genus.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Genus-specific vertex x
    mu.genus.a[j] ~ dexp(mu.a)    # Genus-specific width
    genus.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm.raw) {
    mu[i] <- mu.genus.a[genus.rm.raw[i]] * (mu.genus.Topt[genus.rm.raw[i]] - temp[i])^2 + mu.genus.c[genus.rm.raw[i]]
    trait[i] ~ dnorm(mu[i], genus.tau[genus.rm.raw[i]]) T(0, ) 
  }}
", file = "genus_robin.txt")


# Parameters to monitor
parameters <- c("mu.genus.a", "mu.genus.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt", "genus.tau",
                "tau1", "mu.genus.Topt")


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
jag.data <- list(trait = tick.rm.raw$Trait, N.obs.rm.raw = N.obs.rm.raw, 
                 temp = tick.rm.raw$Temp, 
                 n.genus.rm.raw = n.genus.rm.raw, genus.rm.raw = genus.rm.raw)
# Fit model
genus.raw.andersoni.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                                model.file = "genus_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                                n.iter = ni, DIC = T, working.directory = getwd())
# Make an mcmc object
genus.raw.andersoni.fit.mcmc <- as.mcmc(genus.raw.andersoni.fit)
closeAllConnections()

# Get posterior means using parameters from genus effect
summary_fit.genus.rm.raw <- summary(genus.raw.andersoni.fit.mcmc)$statistics
a_genus_mean.rm.raw <- summary_fit.genus.rm.raw[grep("^mu.genus.a\\[", rownames(summary_fit.genus.rm.raw)), "Mean"]
c_genus_mean.rm.raw <- summary_fit.genus.rm.raw[grep("^mu.genus.c\\[", rownames(summary_fit.genus.rm.raw)), "Mean"]
Topt_genus_mean.rm.raw <- summary_fit.genus.rm.raw[grep("^mu.genus.Topt\\[", rownames(summary_fit.genus.rm.raw)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp) ,
                max(tick.raw$Temp) , by = 0.1)

# Get predicted values for each genus
genus1.rm.raw <- function (x){
  a_genus_mean.rm.raw[1] * (x - Topt_genus_mean.rm.raw[1])^2 + c_genus_mean.rm.raw[1]
}
genus2.rm.raw <- function (x){
  a_genus_mean.rm.raw[2] * (x - Topt_genus_mean.rm.raw[2])^2 + c_genus_mean.rm.raw[2]
}
genus3.rm.raw <- function (x){
  a_genus_mean.rm.raw[3] * (x - Topt_genus_mean.rm.raw[3])^2 + c_genus_mean.rm.raw[3]
}
genus4.rm.raw <- function (x){
  a_genus_mean.rm.raw[4] * (x - Topt_genus_mean.rm.raw[4])^2 + c_genus_mean.rm.raw[4]
}

# Data frame of temperature, fitted values, and genus
fitted_data.genus.andersoni.raw <- data.frame(Temperature = rep(temp_seq, n.genus.rm.raw), 
                                              Fitted = c(genus1.rm.raw(temp_seq),
                                                         genus2.rm.raw(temp_seq),
                                                         genus3.rm.raw(temp_seq),
                                                         genus4.rm.raw(temp_seq)),
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
fitted_data.genus.andersoni.raw <- data.frame(
  Temperature = rep(temp_seq, n.genus.rm.raw),
  Fitted = c(genus1.rm.raw(temp_seq),
             genus2.rm.raw(temp_seq),
             genus3.rm.raw(temp_seq),
             genus4.rm.raw(temp_seq)),
  Genus = factor(rep(c("Amblyomma", "Dermacentor",
                       "Haemaphysalis", "Hyalomma"), each = length(temp_seq)))
) %>%
  left_join(species_to_genus.raw, by = "Genus")

# Plot observed data and fitted curves with facet wrap by species
p.genus.andersoni.raw <- ggplot() +
  geom_point(data = plot_data.genus.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.genus.andersoni.raw, aes(x = Temperature, y = Fitted, color = Genus), size = 1) +
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
p.genus.andersoni.raw

# Merge the observed and fitted data for comparison
merged_data.genus.andersoni.raw <- fitted_data.genus.andersoni.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species", "Genus" = "Genus"))

# Calculate MSE for each species
mse_by_genus.andersoni.raw <- merged_data.genus.andersoni.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.genus.andersoni.raw <- merged_data.genus.andersoni.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_genus.andersoni.raw,
     fitted_data.genus.andersoni.raw,
     merged_data.genus.andersoni.raw,
     genus.raw.andersoni.fit.mcmc,
     genus.raw.andersoni.fit,
     p.genus.andersoni.raw,
     file="andersoni.raw.genus.RData")


################################### By Genus- Nitens ###################################
tick.rm.raw <- tick.raw[tick.raw$Species != "Dermacentor nitens", ]
genus.rm.raw <- as.factor(tick.rm.raw$Genus)
n.genus.rm.raw <- length(unique(tick.rm.raw$Genus))
N.obs.rm.raw <- length(tick.rm.raw$Trait)
unique(tick.rm.raw$Species)
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

  for (j in 1:n.genus.rm.raw) {
    mu.genus.c[j] ~ dexp(mu.c)       # Genus-specific vertex y
    mu.genus.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Genus-specific vertex x
    mu.genus.a[j] ~ dexp(mu.a)    # Genus-specific width
    genus.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm.raw) {
    mu[i] <- mu.genus.a[genus.rm.raw[i]] * (mu.genus.Topt[genus.rm.raw[i]] - temp[i])^2 + mu.genus.c[genus.rm.raw[i]]
    trait[i] ~ dnorm(mu[i], genus.tau[genus.rm.raw[i]]) T(0, ) 
  }}
", file = "genus_robin.txt")


# Parameters to monitor
parameters <- c("mu.genus.a", "mu.genus.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt", "genus.tau",
                "tau1", "mu.genus.Topt")


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
jag.data <- list(trait = tick.rm.raw$Trait, N.obs.rm.raw = N.obs.rm.raw, 
                 temp = tick.rm.raw$Temp, 
                 n.genus.rm.raw = n.genus.rm.raw, genus.rm.raw = genus.rm.raw)
# Fit model
genus.raw.nitens.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                                model.file = "genus_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                                n.iter = ni, DIC = T, working.directory = getwd())
# Make an mcmc object
genus.raw.nitens.fit.mcmc <- as.mcmc(genus.raw.nitens.fit)
closeAllConnections()

# Get posterior means using parameters from genus effect
summary_fit.genus.rm.raw <- summary(genus.raw.nitens.fit.mcmc)$statistics
a_genus_mean.rm.raw <- summary_fit.genus.rm.raw[grep("^mu.genus.a\\[", rownames(summary_fit.genus.rm.raw)), "Mean"]
c_genus_mean.rm.raw <- summary_fit.genus.rm.raw[grep("^mu.genus.c\\[", rownames(summary_fit.genus.rm.raw)), "Mean"]
Topt_genus_mean.rm.raw <- summary_fit.genus.rm.raw[grep("^mu.genus.Topt\\[", rownames(summary_fit.genus.rm.raw)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp) ,
                max(tick.raw$Temp) , by = 0.1)

# Get predicted values for each genus
genus1.rm.raw <- function (x){
  a_genus_mean.rm.raw[1] * (x - Topt_genus_mean.rm.raw[1])^2 + c_genus_mean.rm.raw[1]
}
genus2.rm.raw <- function (x){
  a_genus_mean.rm.raw[2] * (x - Topt_genus_mean.rm.raw[2])^2 + c_genus_mean.rm.raw[2]
}
genus3.rm.raw <- function (x){
  a_genus_mean.rm.raw[3] * (x - Topt_genus_mean.rm.raw[3])^2 + c_genus_mean.rm.raw[3]
}
genus4.rm.raw <- function (x){
  a_genus_mean.rm.raw[4] * (x - Topt_genus_mean.rm.raw[4])^2 + c_genus_mean.rm.raw[4]
}

# Data frame of temperature, fitted values, and genus
fitted_data.genus.nitens.raw <- data.frame(Temperature = rep(temp_seq, n.genus.rm.raw), 
                                              Fitted = c(genus1.rm.raw(temp_seq),
                                                         genus2.rm.raw(temp_seq),
                                                         genus3.rm.raw(temp_seq),
                                                         genus4.rm.raw(temp_seq)),
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
fitted_data.genus.nitens.raw <- data.frame(
  Temperature = rep(temp_seq, n.genus.rm.raw),
  Fitted = c(genus1.rm.raw(temp_seq),
             genus2.rm.raw(temp_seq),
             genus3.rm.raw(temp_seq),
             genus4.rm.raw(temp_seq)),
  Genus = factor(rep(c("Amblyomma", "Dermacentor",
                       "Haemaphysalis", "Hyalomma"), each = length(temp_seq)))
) %>%
  left_join(species_to_genus.raw, by = "Genus")

# Plot observed data and fitted curves with facet wrap by species
p.genus.nitens.raw <- ggplot() +
  geom_point(data = plot_data.genus.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.genus.nitens.raw, aes(x = Temperature, y = Fitted, color = Genus), size = 1) +
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
p.genus.nitens.raw

# Merge the observed and fitted data for comparison
merged_data.genus.nitens.raw <- fitted_data.genus.nitens.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species", "Genus" = "Genus"))

# Calculate MSE for each species
mse_by_genus.nitens.raw <- merged_data.genus.nitens.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.genus.nitens.raw <- merged_data.genus.nitens.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_genus.nitens.raw,
     fitted_data.genus.nitens.raw,
     merged_data.genus.nitens.raw,
     genus.raw.nitens.fit.mcmc,
     genus.raw.nitens.fit,
     p.genus.nitens.raw,
     file="nitens.raw.genus.RData")


################################### By Genus- Aegyptium ###################################
tick.rm.raw <- tick.raw[tick.raw$Species != "Hyalomma aegyptium", ]
genus.rm.raw <- as.factor(tick.rm.raw$Genus)
n.genus.rm.raw <- length(unique(tick.rm.raw$Genus))
N.obs.rm.raw <- length(tick.rm.raw$Trait)
unique(tick.rm.raw$Species)
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

  for (j in 1:n.genus.rm.raw) {
    mu.genus.c[j] ~ dexp(mu.c)       # Genus-specific vertex y
    mu.genus.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Genus-specific vertex x
    mu.genus.a[j] ~ dexp(mu.a)    # Genus-specific width
    genus.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm.raw) {
    mu[i] <- mu.genus.a[genus.rm.raw[i]] * (mu.genus.Topt[genus.rm.raw[i]] - temp[i])^2 + mu.genus.c[genus.rm.raw[i]]
    trait[i] ~ dnorm(mu[i], genus.tau[genus.rm.raw[i]]) T(0, ) 
  }}
", file = "genus_robin.txt")


# Parameters to monitor
parameters <- c("mu.genus.a", "mu.genus.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt", "genus.tau",
                "tau1", "mu.genus.Topt")


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
jag.data <- list(trait = tick.rm.raw$Trait, N.obs.rm.raw = N.obs.rm.raw, 
                 temp = tick.rm.raw$Temp, 
                 n.genus.rm.raw = n.genus.rm.raw, genus.rm.raw = genus.rm.raw)
# Fit model
genus.raw.aegyptium.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                                model.file = "genus_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                                n.iter = ni, DIC = T, working.directory = getwd())
# Make an mcmc object
genus.raw.aegyptium.fit.mcmc <- as.mcmc(genus.raw.aegyptium.fit)
closeAllConnections()

# Get posterior means using parameters from genus effect
summary_fit.genus.rm.raw <- summary(genus.raw.aegyptium.fit.mcmc)$statistics
a_genus_mean.rm.raw <- summary_fit.genus.rm.raw[grep("^mu.genus.a\\[", rownames(summary_fit.genus.rm.raw)), "Mean"]
c_genus_mean.rm.raw <- summary_fit.genus.rm.raw[grep("^mu.genus.c\\[", rownames(summary_fit.genus.rm.raw)), "Mean"]
Topt_genus_mean.rm.raw <- summary_fit.genus.rm.raw[grep("^mu.genus.Topt\\[", rownames(summary_fit.genus.rm.raw)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp) ,
                max(tick.raw$Temp) , by = 0.1)

# Get predicted values for each genus
genus1.rm.raw <- function (x){
  a_genus_mean.rm.raw[1] * (x - Topt_genus_mean.rm.raw[1])^2 + c_genus_mean.rm.raw[1]
}
genus2.rm.raw <- function (x){
  a_genus_mean.rm.raw[2] * (x - Topt_genus_mean.rm.raw[2])^2 + c_genus_mean.rm.raw[2]
}
genus3.rm.raw <- function (x){
  a_genus_mean.rm.raw[3] * (x - Topt_genus_mean.rm.raw[3])^2 + c_genus_mean.rm.raw[3]
}
genus4.rm.raw <- function (x){
  a_genus_mean.rm.raw[4] * (x - Topt_genus_mean.rm.raw[4])^2 + c_genus_mean.rm.raw[4]
}

# Data frame of temperature, fitted values, and genus
fitted_data.genus.aegyptium.raw <- data.frame(Temperature = rep(temp_seq, n.genus.rm.raw), 
                                              Fitted = c(genus1.rm.raw(temp_seq),
                                                         genus2.rm.raw(temp_seq),
                                                         genus3.rm.raw(temp_seq),
                                                         genus4.rm.raw(temp_seq)),
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
fitted_data.genus.aegyptium.raw <- data.frame(
  Temperature = rep(temp_seq, n.genus.rm.raw),
  Fitted = c(genus1.rm.raw(temp_seq),
             genus2.rm.raw(temp_seq),
             genus3.rm.raw(temp_seq),
             genus4.rm.raw(temp_seq)),
  Genus = factor(rep(c("Amblyomma", "Dermacentor",
                       "Haemaphysalis", "Hyalomma"), each = length(temp_seq)))
) %>%
  left_join(species_to_genus.raw, by = "Genus")

# Plot observed data and fitted curves with facet wrap by species
p.genus.aegyptium.raw <- ggplot() +
  geom_point(data = plot_data.genus.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.genus.aegyptium.raw, aes(x = Temperature, y = Fitted, color = Genus), size = 1) +
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
p.genus.aegyptium.raw

# Merge the observed and fitted data for comparison
merged_data.genus.aegyptium.raw <- fitted_data.genus.aegyptium.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species", "Genus" = "Genus"))

# Calculate MSE for each species
mse_by_genus.aegyptium.raw <- merged_data.genus.aegyptium.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.genus.aegyptium.raw <- merged_data.genus.aegyptium.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_genus.aegyptium.raw,
     fitted_data.genus.aegyptium.raw,
     merged_data.genus.aegyptium.raw,
     genus.raw.aegyptium.fit.mcmc,
     genus.raw.aegyptium.fit,
     p.genus.aegyptium.raw,
     file="aegyptium.raw.genus.RData")


################################### By Genus- Dromedarii ###################################
tick.rm.raw <- tick.raw[tick.raw$Species != "Hyalomma dromedarii", ]
genus.rm.raw <- as.factor(tick.rm.raw$Genus)
n.genus.rm.raw <- length(unique(tick.rm.raw$Genus))
N.obs.rm.raw <- length(tick.rm.raw$Trait)
unique(tick.rm.raw$Species)
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

  for (j in 1:n.genus.rm.raw) {
    mu.genus.c[j] ~ dexp(mu.c)       # Genus-specific vertex y
    mu.genus.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Genus-specific vertex x
    mu.genus.a[j] ~ dexp(mu.a)    # Genus-specific width
    genus.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm.raw) {
    mu[i] <- mu.genus.a[genus.rm.raw[i]] * (mu.genus.Topt[genus.rm.raw[i]] - temp[i])^2 + mu.genus.c[genus.rm.raw[i]]
    trait[i] ~ dnorm(mu[i], genus.tau[genus.rm.raw[i]]) T(0, ) 
  }}
", file = "genus_robin.txt")


# Parameters to monitor
parameters <- c("mu.genus.a", "mu.genus.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt", "genus.tau",
                "tau1", "mu.genus.Topt")


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
jag.data <- list(trait = tick.rm.raw$Trait, N.obs.rm.raw = N.obs.rm.raw, 
                 temp = tick.rm.raw$Temp, 
                 n.genus.rm.raw = n.genus.rm.raw, genus.rm.raw = genus.rm.raw)
# Fit model
genus.raw.dromedarii.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                                model.file = "genus_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                                n.iter = ni, DIC = T, working.directory = getwd())
# Make an mcmc object
genus.raw.dromedarii.fit.mcmc <- as.mcmc(genus.raw.dromedarii.fit)
closeAllConnections()

# Get posterior means using parameters from genus effect
summary_fit.genus.rm.raw <- summary(genus.raw.dromedarii.fit.mcmc)$statistics
a_genus_mean.rm.raw <- summary_fit.genus.rm.raw[grep("^mu.genus.a\\[", rownames(summary_fit.genus.rm.raw)), "Mean"]
c_genus_mean.rm.raw <- summary_fit.genus.rm.raw[grep("^mu.genus.c\\[", rownames(summary_fit.genus.rm.raw)), "Mean"]
Topt_genus_mean.rm.raw <- summary_fit.genus.rm.raw[grep("^mu.genus.Topt\\[", rownames(summary_fit.genus.rm.raw)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp) ,
                max(tick.raw$Temp) , by = 0.1)

# Get predicted values for each genus
genus1.rm.raw <- function (x){
  a_genus_mean.rm.raw[1] * (x - Topt_genus_mean.rm.raw[1])^2 + c_genus_mean.rm.raw[1]
}
genus2.rm.raw <- function (x){
  a_genus_mean.rm.raw[2] * (x - Topt_genus_mean.rm.raw[2])^2 + c_genus_mean.rm.raw[2]
}
genus3.rm.raw <- function (x){
  a_genus_mean.rm.raw[3] * (x - Topt_genus_mean.rm.raw[3])^2 + c_genus_mean.rm.raw[3]
}
genus4.rm.raw <- function (x){
  a_genus_mean.rm.raw[4] * (x - Topt_genus_mean.rm.raw[4])^2 + c_genus_mean.rm.raw[4]
}

# Data frame of temperature, fitted values, and genus
fitted_data.genus.dromedarii.raw <- data.frame(Temperature = rep(temp_seq, n.genus.rm.raw), 
                                              Fitted = c(genus1.rm.raw(temp_seq),
                                                         genus2.rm.raw(temp_seq),
                                                         genus3.rm.raw(temp_seq),
                                                         genus4.rm.raw(temp_seq)),
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
fitted_data.genus.dromedarii.raw <- data.frame(
  Temperature = rep(temp_seq, n.genus.rm.raw),
  Fitted = c(genus1.rm.raw(temp_seq),
             genus2.rm.raw(temp_seq),
             genus3.rm.raw(temp_seq),
             genus4.rm.raw(temp_seq)),
  Genus = factor(rep(c("Amblyomma", "Dermacentor",
                       "Haemaphysalis", "Hyalomma"), each = length(temp_seq)))
) %>%
  left_join(species_to_genus.raw, by = "Genus")

# Plot observed data and fitted curves with facet wrap by species
p.genus.dromedarii.raw <- ggplot() +
  geom_point(data = plot_data.genus.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.genus.dromedarii.raw, aes(x = Temperature, y = Fitted, color = Genus), size = 1) +
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
p.genus.dromedarii.raw

# Merge the observed and fitted data for comparison
merged_data.genus.dromedarii.raw <- fitted_data.genus.dromedarii.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species", "Genus" = "Genus"))

# Calculate MSE for each species
mse_by_genus.dromedarii.raw <- merged_data.genus.dromedarii.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.genus.dromedarii.raw <- merged_data.genus.dromedarii.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_genus.dromedarii.raw,
     fitted_data.genus.dromedarii.raw,
     merged_data.genus.dromedarii.raw,
     genus.raw.dromedarii.fit.mcmc,
     genus.raw.dromedarii.fit,
     p.genus.dromedarii.raw,
     file="dromedarii.raw.genus.RData")


################################### By Genus- Impeltatum ###################################
tick.rm.raw <- tick.raw[tick.raw$Species != "Hyalomma impeltatum", ]
genus.rm.raw <- as.factor(tick.rm.raw$Genus)
n.genus.rm.raw <- length(unique(tick.rm.raw$Genus))
N.obs.rm.raw <- length(tick.rm.raw$Trait)
unique(tick.rm.raw$Species)
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

  for (j in 1:n.genus.rm.raw) {
    mu.genus.c[j] ~ dexp(mu.c)       # Genus-specific vertex y
    mu.genus.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Genus-specific vertex x
    mu.genus.a[j] ~ dexp(mu.a)    # Genus-specific width
    genus.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm.raw) {
    mu[i] <- mu.genus.a[genus.rm.raw[i]] * (mu.genus.Topt[genus.rm.raw[i]] - temp[i])^2 + mu.genus.c[genus.rm.raw[i]]
    trait[i] ~ dnorm(mu[i], genus.tau[genus.rm.raw[i]]) T(0, ) 
  }}
", file = "genus_robin.txt")


# Parameters to monitor
parameters <- c("mu.genus.a", "mu.genus.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt", "genus.tau",
                "tau1", "mu.genus.Topt")


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
jag.data <- list(trait = tick.rm.raw$Trait, N.obs.rm.raw = N.obs.rm.raw, 
                 temp = tick.rm.raw$Temp, 
                 n.genus.rm.raw = n.genus.rm.raw, genus.rm.raw = genus.rm.raw)
# Fit model
genus.raw.impeltatum.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                                model.file = "genus_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                                n.iter = ni, DIC = T, working.directory = getwd())
# Make an mcmc object
genus.raw.impeltatum.fit.mcmc <- as.mcmc(genus.raw.impeltatum.fit)
closeAllConnections()

# Get posterior means using parameters from genus effect
summary_fit.genus.rm.raw <- summary(genus.raw.impeltatum.fit.mcmc)$statistics
a_genus_mean.rm.raw <- summary_fit.genus.rm.raw[grep("^mu.genus.a\\[", rownames(summary_fit.genus.rm.raw)), "Mean"]
c_genus_mean.rm.raw <- summary_fit.genus.rm.raw[grep("^mu.genus.c\\[", rownames(summary_fit.genus.rm.raw)), "Mean"]
Topt_genus_mean.rm.raw <- summary_fit.genus.rm.raw[grep("^mu.genus.Topt\\[", rownames(summary_fit.genus.rm.raw)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp) ,
                max(tick.raw$Temp) , by = 0.1)

# Get predicted values for each genus
genus1.rm.raw <- function (x){
  a_genus_mean.rm.raw[1] * (x - Topt_genus_mean.rm.raw[1])^2 + c_genus_mean.rm.raw[1]
}
genus2.rm.raw <- function (x){
  a_genus_mean.rm.raw[2] * (x - Topt_genus_mean.rm.raw[2])^2 + c_genus_mean.rm.raw[2]
}
genus3.rm.raw <- function (x){
  a_genus_mean.rm.raw[3] * (x - Topt_genus_mean.rm.raw[3])^2 + c_genus_mean.rm.raw[3]
}
genus4.rm.raw <- function (x){
  a_genus_mean.rm.raw[4] * (x - Topt_genus_mean.rm.raw[4])^2 + c_genus_mean.rm.raw[4]
}

# Data frame of temperature, fitted values, and genus
fitted_data.genus.impeltatum.raw <- data.frame(Temperature = rep(temp_seq, n.genus.rm.raw), 
                                              Fitted = c(genus1.rm.raw(temp_seq),
                                                         genus2.rm.raw(temp_seq),
                                                         genus3.rm.raw(temp_seq),
                                                         genus4.rm.raw(temp_seq)),
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
fitted_data.genus.impeltatum.raw <- data.frame(
  Temperature = rep(temp_seq, n.genus.rm.raw),
  Fitted = c(genus1.rm.raw(temp_seq),
             genus2.rm.raw(temp_seq),
             genus3.rm.raw(temp_seq),
             genus4.rm.raw(temp_seq)),
  Genus = factor(rep(c("Amblyomma", "Dermacentor",
                       "Haemaphysalis", "Hyalomma"), each = length(temp_seq)))
) %>%
  left_join(species_to_genus.raw, by = "Genus")

# Plot observed data and fitted curves with facet wrap by species
p.genus.impeltatum.raw <- ggplot() +
  geom_point(data = plot_data.genus.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.genus.impeltatum.raw, aes(x = Temperature, y = Fitted, color = Genus), size = 1) +
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
p.genus.impeltatum.raw

# Merge the observed and fitted data for comparison
merged_data.genus.impeltatum.raw <- fitted_data.genus.impeltatum.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species", "Genus" = "Genus"))

# Calculate MSE for each species
mse_by_genus.impeltatum.raw <- merged_data.genus.impeltatum.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.genus.impeltatum.raw <- merged_data.genus.impeltatum.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_genus.impeltatum.raw,
     fitted_data.genus.impeltatum.raw,
     merged_data.genus.impeltatum.raw,
     genus.raw.impeltatum.fit.mcmc,
     genus.raw.impeltatum.fit,
     p.genus.impeltatum.raw,
     file="impeltatum.raw.genus.RData")


################################### By Genus- Lusitanicum ###################################
tick.rm.raw <- tick.raw[tick.raw$Species != "Hyalomma lusitanicum", ]
genus.rm.raw <- as.factor(tick.rm.raw$Genus)
n.genus.rm.raw <- length(unique(tick.rm.raw$Genus))
N.obs.rm.raw <- length(tick.rm.raw$Trait)
unique(tick.rm.raw$Species)
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

  for (j in 1:n.genus.rm.raw) {
    mu.genus.c[j] ~ dexp(mu.c)       # Genus-specific vertex y
    mu.genus.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Genus-specific vertex x
    mu.genus.a[j] ~ dexp(mu.a)    # Genus-specific width
    genus.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm.raw) {
    mu[i] <- mu.genus.a[genus.rm.raw[i]] * (mu.genus.Topt[genus.rm.raw[i]] - temp[i])^2 + mu.genus.c[genus.rm.raw[i]]
    trait[i] ~ dnorm(mu[i], genus.tau[genus.rm.raw[i]]) T(0, ) 
  }}
", file = "genus_robin.txt")


# Parameters to monitor
parameters <- c("mu.genus.a", "mu.genus.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt", "genus.tau",
                "tau1", "mu.genus.Topt")


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
jag.data <- list(trait = tick.rm.raw$Trait, N.obs.rm.raw = N.obs.rm.raw, 
                 temp = tick.rm.raw$Temp, 
                 n.genus.rm.raw = n.genus.rm.raw, genus.rm.raw = genus.rm.raw)
# Fit model
genus.raw.lusitanicum.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                                model.file = "genus_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                                n.iter = ni, DIC = T, working.directory = getwd())
# Make an mcmc object
genus.raw.lusitanicum.fit.mcmc <- as.mcmc(genus.raw.lusitanicum.fit)
closeAllConnections()

# Get posterior means using parameters from genus effect
summary_fit.genus.rm.raw <- summary(genus.raw.lusitanicum.fit.mcmc)$statistics
a_genus_mean.rm.raw <- summary_fit.genus.rm.raw[grep("^mu.genus.a\\[", rownames(summary_fit.genus.rm.raw)), "Mean"]
c_genus_mean.rm.raw <- summary_fit.genus.rm.raw[grep("^mu.genus.c\\[", rownames(summary_fit.genus.rm.raw)), "Mean"]
Topt_genus_mean.rm.raw <- summary_fit.genus.rm.raw[grep("^mu.genus.Topt\\[", rownames(summary_fit.genus.rm.raw)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp) ,
                max(tick.raw$Temp) , by = 0.1)

# Get predicted values for each genus
genus1.rm.raw <- function (x){
  a_genus_mean.rm.raw[1] * (x - Topt_genus_mean.rm.raw[1])^2 + c_genus_mean.rm.raw[1]
}
genus2.rm.raw <- function (x){
  a_genus_mean.rm.raw[2] * (x - Topt_genus_mean.rm.raw[2])^2 + c_genus_mean.rm.raw[2]
}
genus3.rm.raw <- function (x){
  a_genus_mean.rm.raw[3] * (x - Topt_genus_mean.rm.raw[3])^2 + c_genus_mean.rm.raw[3]
}
genus4.rm.raw <- function (x){
  a_genus_mean.rm.raw[4] * (x - Topt_genus_mean.rm.raw[4])^2 + c_genus_mean.rm.raw[4]
}

# Data frame of temperature, fitted values, and genus
fitted_data.genus.lusitanicum.raw <- data.frame(Temperature = rep(temp_seq, n.genus.rm.raw), 
                                              Fitted = c(genus1.rm.raw(temp_seq),
                                                         genus2.rm.raw(temp_seq),
                                                         genus3.rm.raw(temp_seq),
                                                         genus4.rm.raw(temp_seq)),
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
fitted_data.genus.lusitanicum.raw <- data.frame(
  Temperature = rep(temp_seq, n.genus.rm.raw),
  Fitted = c(genus1.rm.raw(temp_seq),
             genus2.rm.raw(temp_seq),
             genus3.rm.raw(temp_seq),
             genus4.rm.raw(temp_seq)),
  Genus = factor(rep(c("Amblyomma", "Dermacentor",
                       "Haemaphysalis", "Hyalomma"), each = length(temp_seq)))
) %>%
  left_join(species_to_genus.raw, by = "Genus")

# Plot observed data and fitted curves with facet wrap by species
p.genus.lusitanicum.raw <- ggplot() +
  geom_point(data = plot_data.genus.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.genus.lusitanicum.raw, aes(x = Temperature, y = Fitted, color = Genus), size = 1) +
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
p.genus.lusitanicum.raw

# Merge the observed and fitted data for comparison
merged_data.genus.lusitanicum.raw <- fitted_data.genus.lusitanicum.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species", "Genus" = "Genus"))

# Calculate MSE for each species
mse_by_genus.lusitanicum.raw <- merged_data.genus.lusitanicum.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.genus.lusitanicum.raw <- merged_data.genus.lusitanicum.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_genus.lusitanicum.raw,
     fitted_data.genus.lusitanicum.raw,
     merged_data.genus.lusitanicum.raw,
     genus.raw.lusitanicum.fit.mcmc,
     genus.raw.lusitanicum.fit,
     p.genus.lusitanicum.raw,
     file="lusitanicum.raw.genus.RData")


################################### By Genus- Schulzei ###################################
tick.rm.raw <- tick.raw[tick.raw$Species != "Hyalomma schulzei", ]
genus.rm.raw <- as.factor(tick.rm.raw$Genus)
n.genus.rm.raw <- length(unique(tick.rm.raw$Genus))
N.obs.rm.raw <- length(tick.rm.raw$Trait)
unique(tick.rm.raw$Species)
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

  for (j in 1:n.genus.rm.raw) {
    mu.genus.c[j] ~ dexp(mu.c)       # Genus-specific vertex y
    mu.genus.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Genus-specific vertex x
    mu.genus.a[j] ~ dexp(mu.a)    # Genus-specific width
    genus.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm.raw) {
    mu[i] <- mu.genus.a[genus.rm.raw[i]] * (mu.genus.Topt[genus.rm.raw[i]] - temp[i])^2 + mu.genus.c[genus.rm.raw[i]]
    trait[i] ~ dnorm(mu[i], genus.tau[genus.rm.raw[i]]) T(0, ) 
  }}
", file = "genus_robin.txt")


# Parameters to monitor
parameters <- c("mu.genus.a", "mu.genus.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt", "genus.tau",
                "tau1", "mu.genus.Topt")


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
jag.data <- list(trait = tick.rm.raw$Trait, N.obs.rm.raw = N.obs.rm.raw, 
                 temp = tick.rm.raw$Temp, 
                 n.genus.rm.raw = n.genus.rm.raw, genus.rm.raw = genus.rm.raw)
# Fit model
genus.raw.schulzei.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                                model.file = "genus_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                                n.iter = ni, DIC = T, working.directory = getwd())
# Make an mcmc object
genus.raw.schulzei.fit.mcmc <- as.mcmc(genus.raw.schulzei.fit)
closeAllConnections()

# Get posterior means using parameters from genus effect
summary_fit.genus.rm.raw <- summary(genus.raw.schulzei.fit.mcmc)$statistics
a_genus_mean.rm.raw <- summary_fit.genus.rm.raw[grep("^mu.genus.a\\[", rownames(summary_fit.genus.rm.raw)), "Mean"]
c_genus_mean.rm.raw <- summary_fit.genus.rm.raw[grep("^mu.genus.c\\[", rownames(summary_fit.genus.rm.raw)), "Mean"]
Topt_genus_mean.rm.raw <- summary_fit.genus.rm.raw[grep("^mu.genus.Topt\\[", rownames(summary_fit.genus.rm.raw)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp) ,
                max(tick.raw$Temp) , by = 0.1)

# Get predicted values for each genus
genus1.rm.raw <- function (x){
  a_genus_mean.rm.raw[1] * (x - Topt_genus_mean.rm.raw[1])^2 + c_genus_mean.rm.raw[1]
}
genus2.rm.raw <- function (x){
  a_genus_mean.rm.raw[2] * (x - Topt_genus_mean.rm.raw[2])^2 + c_genus_mean.rm.raw[2]
}
genus3.rm.raw <- function (x){
  a_genus_mean.rm.raw[3] * (x - Topt_genus_mean.rm.raw[3])^2 + c_genus_mean.rm.raw[3]
}
genus4.rm.raw <- function (x){
  a_genus_mean.rm.raw[4] * (x - Topt_genus_mean.rm.raw[4])^2 + c_genus_mean.rm.raw[4]
}

# Data frame of temperature, fitted values, and genus
fitted_data.genus.schulzei.raw <- data.frame(Temperature = rep(temp_seq, n.genus.rm.raw), 
                                              Fitted = c(genus1.rm.raw(temp_seq),
                                                         genus2.rm.raw(temp_seq),
                                                         genus3.rm.raw(temp_seq),
                                                         genus4.rm.raw(temp_seq)),
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
fitted_data.genus.schulzei.raw <- data.frame(
  Temperature = rep(temp_seq, n.genus.rm.raw),
  Fitted = c(genus1.rm.raw(temp_seq),
             genus2.rm.raw(temp_seq),
             genus3.rm.raw(temp_seq),
             genus4.rm.raw(temp_seq)),
  Genus = factor(rep(c("Amblyomma", "Dermacentor",
                       "Haemaphysalis", "Hyalomma"), each = length(temp_seq)))
) %>%
  left_join(species_to_genus.raw, by = "Genus")

# Plot observed data and fitted curves with facet wrap by species
p.genus.schulzei.raw <- ggplot() +
  geom_point(data = plot_data.genus.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.genus.schulzei.raw, aes(x = Temperature, y = Fitted, color = Genus), size = 1) +
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
p.genus.schulzei.raw

# Merge the observed and fitted data for comparison
merged_data.genus.schulzei.raw <- fitted_data.genus.schulzei.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species", "Genus" = "Genus"))

# Calculate MSE for each species
mse_by_genus.schulzei.raw <- merged_data.genus.schulzei.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.genus.schulzei.raw <- merged_data.genus.schulzei.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_genus.schulzei.raw,
     fitted_data.genus.schulzei.raw,
     merged_data.genus.schulzei.raw,
     genus.raw.schulzei.fit.mcmc,
     genus.raw.schulzei.fit,
     p.genus.schulzei.raw,
     file="schulzei.raw.genus.RData")



################################### By Host - Lepidum ###################################
tick.rm.raw <- tick.raw[tick.raw$Species != "Amblyomma lepidum", ]
host.rm.raw <- as.factor(tick.rm.raw$Host)
n.host.rm.raw <- length(unique(tick.rm.raw$Host))
N.obs.rm.raw <- length(tick.rm.raw$Trait)
unique(tick.rm.raw$Species)
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

  for (j in 1:n.host.rm.raw) {
    mu.host.c[j] ~ dexp(mu.c)       # Host-specific vertex y
    mu.host.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Host-specific vertex x
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
    host.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm.raw) {
    mu[i] <- mu.host.a[host.rm.raw[i]] * (mu.host.Topt[host.rm.raw[i]] - temp[i])^2 + mu.host.c[host.rm.raw[i]]
    trait[i] ~ dnorm(mu[i], host.tau[host.rm.raw[i]]) T(0, )   
  }}
", file = "host_robin.txt")


# Parameters to monitor
parameters <- c("mu.host.a", "mu.host.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt", "host.tau",
                "tau1", "mu.host.Topt")


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
jag.data <- list(trait = tick.rm.raw$Trait, N.obs.rm.raw = N.obs.rm.raw, 
                 temp = tick.rm.raw$Temp, 
                 n.host.rm.raw = n.host.rm.raw, host.rm.raw = host.rm.raw)
# Fit model
host.raw.lepidum.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                             model.file = "host_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                             n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.raw.lepidum.fit.mcmc <- as.mcmc(host.raw.lepidum.fit) 
closeAllConnections()


# Get posterior means for parameters for each host group
summary_fit.host.rm.raw <- summary(host.raw.lepidum.fit.mcmc)$statistics
a_host_mean.rm.raw <- summary_fit.host.rm.raw[grep("^mu.host.a\\[", rownames(summary_fit.host.rm.raw)), "Mean"]
c_host_mean.rm.raw <- summary_fit.host.rm.raw[grep("^mu.host.c\\[", rownames(summary_fit.host.rm.raw)), "Mean"]
Topt_host_mean.rm.raw <- summary_fit.host.rm.raw[grep("^mu.host.Topt\\[", rownames(summary_fit.host.rm.raw)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp) ,
                max(tick.raw$Temp) , by = 0.1)

# Get predicted values for each host group
host1.rm.raw <- function (x){
  a_host_mean.rm.raw[1] * (x - Topt_host_mean.rm.raw[1])^2 + c_host_mean.rm.raw[1]
}
host2.rm.raw <- function (x){
  a_host_mean.rm.raw[2] * (x - Topt_host_mean.rm.raw[2])^2 + c_host_mean.rm.raw[2]
}
host3.rm.raw <- function (x){
  a_host_mean.rm.raw[3] * (x - Topt_host_mean.rm.raw[3])^2 + c_host_mean.rm.raw[3]
}


# Data frame with temperature, fitted value, and host group
fitted_data.host.lepidum.raw <- data.frame(Temperature = rep(temp_seq, n.host.rm.raw), 
                                           Fitted = c(host1.rm.raw(temp_seq),
                                                      host2.rm.raw(temp_seq),
                                                      host3.rm.raw(temp_seq)),
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
fitted_data.host.lepidum.raw <- data.frame(
  Temperature = rep(temp_seq, n.host.rm.raw),
  Fitted = c(host1.rm.raw(temp_seq),
             host2.rm.raw(temp_seq),
             host3.rm.raw(temp_seq)),
  Host = factor(rep(c("One", "Three", "Two"), each = length(temp_seq)))
) %>%
  left_join(species_to_host.raw, by = "Host")

plot_data.host.raw$Host <- factor(plot_data.host.raw$Host, levels = c("One", "Two", "Three"))
fitted_data.host.lepidum.raw$Host <- factor(fitted_data.host.lepidum.raw$Host, levels = c("One", "Two", "Three"))
# Plot observed data and fitted curves with facet wrap by species
p.host.lepidum.raw <- ggplot() +
  geom_point(data = plot_data.host.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.lepidum.raw, aes(x = Temperature, y = Fitted, color = Host), size = 1) +
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
p.host.lepidum.raw

# Merge the observed and predicted data for comparison
merged_data.host.lepidum.raw <- fitted_data.host.lepidum.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host"))

# Calculate RMSE for each species
mse_by_host.lepidum.raw <- merged_data.host.lepidum.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.host.lepidum.raw <- merged_data.host.lepidum.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_host.lepidum.raw,
     fitted_data.host.lepidum.raw,
     merged_data.host.lepidum.raw,
     host.raw.lepidum.fit.mcmc,
     host.raw.lepidum.fit,
     p.host.lepidum.raw,
     file="lepidum.raw.host.RData")

################################### By Host - Andersoni ###################################
tick.rm.raw <- tick.raw[tick.raw$Species != "Dermacentor andersoni", ]
host.rm.raw <- as.factor(tick.rm.raw$Host)
n.host.rm.raw <- length(unique(tick.rm.raw$Host))
N.obs.rm.raw <- length(tick.rm.raw$Trait)
unique(tick.rm.raw$Species)
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

  for (j in 1:n.host.rm.raw) {
    mu.host.c[j] ~ dexp(mu.c)       # Host-specific vertex y
    mu.host.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Host-specific vertex x
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
    host.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm.raw) {
    mu[i] <- mu.host.a[host.rm.raw[i]] * (mu.host.Topt[host.rm.raw[i]] - temp[i])^2 + mu.host.c[host.rm.raw[i]]
    trait[i] ~ dnorm(mu[i], host.tau[host.rm.raw[i]]) T(0, )   
  }}
", file = "host_robin.txt")


# Parameters to monitor
parameters <- c("mu.host.a", "mu.host.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt", "host.tau",
                "tau1", "mu.host.Topt")


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
jag.data <- list(trait = tick.rm.raw$Trait, N.obs.rm.raw = N.obs.rm.raw, 
                 temp = tick.rm.raw$Temp, 
                 n.host.rm.raw = n.host.rm.raw, host.rm.raw = host.rm.raw)
# Fit model
host.raw.andersoni.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                             model.file = "host_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                             n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.raw.andersoni.fit.mcmc <- as.mcmc(host.raw.andersoni.fit) 
closeAllConnections()


# Get posterior means for parameters for each host group
summary_fit.host.rm.raw <- summary(host.raw.andersoni.fit.mcmc)$statistics
a_host_mean.rm.raw <- summary_fit.host.rm.raw[grep("^mu.host.a\\[", rownames(summary_fit.host.rm.raw)), "Mean"]
c_host_mean.rm.raw <- summary_fit.host.rm.raw[grep("^mu.host.c\\[", rownames(summary_fit.host.rm.raw)), "Mean"]
Topt_host_mean.rm.raw <- summary_fit.host.rm.raw[grep("^mu.host.Topt\\[", rownames(summary_fit.host.rm.raw)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp) ,
                max(tick.raw$Temp) , by = 0.1)

# Get predicted values for each host group
host1.rm.raw <- function (x){
  a_host_mean.rm.raw[1] * (x - Topt_host_mean.rm.raw[1])^2 + c_host_mean.rm.raw[1]
}
host2.rm.raw <- function (x){
  a_host_mean.rm.raw[2] * (x - Topt_host_mean.rm.raw[2])^2 + c_host_mean.rm.raw[2]
}
host3.rm.raw <- function (x){
  a_host_mean.rm.raw[3] * (x - Topt_host_mean.rm.raw[3])^2 + c_host_mean.rm.raw[3]
}


# Data frame with temperature, fitted value, and host group
fitted_data.host.andersoni.raw <- data.frame(Temperature = rep(temp_seq, n.host.rm.raw), 
                                           Fitted = c(host1.rm.raw(temp_seq),
                                                      host2.rm.raw(temp_seq),
                                                      host3.rm.raw(temp_seq)),
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
fitted_data.host.andersoni.raw <- data.frame(
  Temperature = rep(temp_seq, n.host.rm.raw),
  Fitted = c(host1.rm.raw(temp_seq),
             host2.rm.raw(temp_seq),
             host3.rm.raw(temp_seq)),
  Host = factor(rep(c("One", "Three", "Two"), each = length(temp_seq)))
) %>%
  left_join(species_to_host.raw, by = "Host")

plot_data.host.raw$Host <- factor(plot_data.host.raw$Host, levels = c("One", "Two", "Three"))
fitted_data.host.andersoni.raw$Host <- factor(fitted_data.host.andersoni.raw$Host, levels = c("One", "Two", "Three"))
# Plot observed data and fitted curves with facet wrap by species
p.host.andersoni.raw <- ggplot() +
  geom_point(data = plot_data.host.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.andersoni.raw, aes(x = Temperature, y = Fitted, color = Host), size = 1) +
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
p.host.andersoni.raw

# Merge the observed and predicted data for comparison
merged_data.host.andersoni.raw <- fitted_data.host.andersoni.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host"))

# Calculate RMSE for each species
mse_by_host.andersoni.raw <- merged_data.host.andersoni.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.host.andersoni.raw <- merged_data.host.andersoni.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_host.andersoni.raw,
     fitted_data.host.andersoni.raw,
     merged_data.host.andersoni.raw,
     host.raw.andersoni.fit.mcmc,
     host.raw.andersoni.fit,
     p.host.andersoni.raw,
     file="andersoni.raw.host.RData")

################################### By Host - Leporispalustris ###################################
tick.rm.raw <- tick.raw[tick.raw$Species != "Haemaphysalis leporispalustris", ]
host.rm.raw <- as.factor(tick.rm.raw$Host)
n.host.rm.raw <- length(unique(tick.rm.raw$Host))
N.obs.rm.raw <- length(tick.rm.raw$Trait)
unique(tick.rm.raw$Species)
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

  for (j in 1:n.host.rm.raw) {
    mu.host.c[j] ~ dexp(mu.c)       # Host-specific vertex y
    mu.host.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Host-specific vertex x
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
    host.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm.raw) {
    mu[i] <- mu.host.a[host.rm.raw[i]] * (mu.host.Topt[host.rm.raw[i]] - temp[i])^2 + mu.host.c[host.rm.raw[i]]
    trait[i] ~ dnorm(mu[i], host.tau[host.rm.raw[i]]) T(0, )   
  }}
", file = "host_robin.txt")


# Parameters to monitor
parameters <- c("mu.host.a", "mu.host.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt", "host.tau",
                "tau1", "mu.host.Topt")


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
jag.data <- list(trait = tick.rm.raw$Trait, N.obs.rm.raw = N.obs.rm.raw, 
                 temp = tick.rm.raw$Temp, 
                 n.host.rm.raw = n.host.rm.raw, host.rm.raw = host.rm.raw)
# Fit model
host.raw.leporispalustris.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                             model.file = "host_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                             n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.raw.leporispalustris.fit.mcmc <- as.mcmc(host.raw.leporispalustris.fit) 
closeAllConnections()


# Get posterior means for parameters for each host group
summary_fit.host.rm.raw <- summary(host.raw.leporispalustris.fit.mcmc)$statistics
a_host_mean.rm.raw <- summary_fit.host.rm.raw[grep("^mu.host.a\\[", rownames(summary_fit.host.rm.raw)), "Mean"]
c_host_mean.rm.raw <- summary_fit.host.rm.raw[grep("^mu.host.c\\[", rownames(summary_fit.host.rm.raw)), "Mean"]
Topt_host_mean.rm.raw <- summary_fit.host.rm.raw[grep("^mu.host.Topt\\[", rownames(summary_fit.host.rm.raw)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp) ,
                max(tick.raw$Temp) , by = 0.1)

# Get predicted values for each host group
host1.rm.raw <- function (x){
  a_host_mean.rm.raw[1] * (x - Topt_host_mean.rm.raw[1])^2 + c_host_mean.rm.raw[1]
}
host2.rm.raw <- function (x){
  a_host_mean.rm.raw[2] * (x - Topt_host_mean.rm.raw[2])^2 + c_host_mean.rm.raw[2]
}
host3.rm.raw <- function (x){
  a_host_mean.rm.raw[3] * (x - Topt_host_mean.rm.raw[3])^2 + c_host_mean.rm.raw[3]
}


# Data frame with temperature, fitted value, and host group
fitted_data.host.leporispalustris.raw <- data.frame(Temperature = rep(temp_seq, n.host.rm.raw), 
                                           Fitted = c(host1.rm.raw(temp_seq),
                                                      host2.rm.raw(temp_seq),
                                                      host3.rm.raw(temp_seq)),
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
fitted_data.host.leporispalustris.raw <- data.frame(
  Temperature = rep(temp_seq, n.host.rm.raw),
  Fitted = c(host1.rm.raw(temp_seq),
             host2.rm.raw(temp_seq),
             host3.rm.raw(temp_seq)),
  Host = factor(rep(c("One", "Three", "Two"), each = length(temp_seq)))
) %>%
  left_join(species_to_host.raw, by = "Host")

plot_data.host.raw$Host <- factor(plot_data.host.raw$Host, levels = c("One", "Two", "Three"))
fitted_data.host.leporispalustris.raw$Host <- factor(fitted_data.host.leporispalustris.raw$Host, levels = c("One", "Two", "Three"))
# Plot observed data and fitted curves with facet wrap by species
p.host.leporispalustris.raw <- ggplot() +
  geom_point(data = plot_data.host.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.leporispalustris.raw, aes(x = Temperature, y = Fitted, color = Host), size = 1) +
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
p.host.leporispalustris.raw

# Merge the observed and predicted data for comparison
merged_data.host.leporispalustris.raw <- fitted_data.host.leporispalustris.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host"))

# Calculate RMSE for each species
mse_by_host.leporispalustris.raw <- merged_data.host.leporispalustris.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.host.leporispalustris.raw <- merged_data.host.leporispalustris.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_host.leporispalustris.raw,
     fitted_data.host.leporispalustris.raw,
     merged_data.host.leporispalustris.raw,
     host.raw.leporispalustris.fit.mcmc,
     host.raw.leporispalustris.fit,
     p.host.leporispalustris.raw,
     file="leporispalustris.raw.host.RData")

################################### By Host - Aegyptium ###################################
tick.rm.raw <- tick.raw[tick.raw$Species != "Hyalomma aegyptium", ]
host.rm.raw <- as.factor(tick.rm.raw$Host)
n.host.rm.raw <- length(unique(tick.rm.raw$Host))
N.obs.rm.raw <- length(tick.rm.raw$Trait)
unique(tick.rm.raw$Species)
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

  for (j in 1:n.host.rm.raw) {
    mu.host.c[j] ~ dexp(mu.c)       # Host-specific vertex y
    mu.host.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Host-specific vertex x
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
    host.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm.raw) {
    mu[i] <- mu.host.a[host.rm.raw[i]] * (mu.host.Topt[host.rm.raw[i]] - temp[i])^2 + mu.host.c[host.rm.raw[i]]
    trait[i] ~ dnorm(mu[i], host.tau[host.rm.raw[i]]) T(0, )   
  }}
", file = "host_robin.txt")


# Parameters to monitor
parameters <- c("mu.host.a", "mu.host.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt", "host.tau",
                "tau1", "mu.host.Topt")


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
jag.data <- list(trait = tick.rm.raw$Trait, N.obs.rm.raw = N.obs.rm.raw, 
                 temp = tick.rm.raw$Temp, 
                 n.host.rm.raw = n.host.rm.raw, host.rm.raw = host.rm.raw)
# Fit model
host.raw.aegyptium.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                             model.file = "host_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                             n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.raw.aegyptium.fit.mcmc <- as.mcmc(host.raw.aegyptium.fit) 
closeAllConnections()


# Get posterior means for parameters for each host group
summary_fit.host.rm.raw <- summary(host.raw.aegyptium.fit.mcmc)$statistics
a_host_mean.rm.raw <- summary_fit.host.rm.raw[grep("^mu.host.a\\[", rownames(summary_fit.host.rm.raw)), "Mean"]
c_host_mean.rm.raw <- summary_fit.host.rm.raw[grep("^mu.host.c\\[", rownames(summary_fit.host.rm.raw)), "Mean"]
Topt_host_mean.rm.raw <- summary_fit.host.rm.raw[grep("^mu.host.Topt\\[", rownames(summary_fit.host.rm.raw)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp) ,
                max(tick.raw$Temp) , by = 0.1)

# Get predicted values for each host group
host1.rm.raw <- function (x){
  a_host_mean.rm.raw[1] * (x - Topt_host_mean.rm.raw[1])^2 + c_host_mean.rm.raw[1]
}
host2.rm.raw <- function (x){
  a_host_mean.rm.raw[2] * (x - Topt_host_mean.rm.raw[2])^2 + c_host_mean.rm.raw[2]
}
host3.rm.raw <- function (x){
  a_host_mean.rm.raw[3] * (x - Topt_host_mean.rm.raw[3])^2 + c_host_mean.rm.raw[3]
}


# Data frame with temperature, fitted value, and host group
fitted_data.host.aegyptium.raw <- data.frame(Temperature = rep(temp_seq, n.host.rm.raw), 
                                           Fitted = c(host1.rm.raw(temp_seq),
                                                      host2.rm.raw(temp_seq),
                                                      host3.rm.raw(temp_seq)),
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
fitted_data.host.aegyptium.raw <- data.frame(
  Temperature = rep(temp_seq, n.host.rm.raw),
  Fitted = c(host1.rm.raw(temp_seq),
             host2.rm.raw(temp_seq),
             host3.rm.raw(temp_seq)),
  Host = factor(rep(c("One", "Three", "Two"), each = length(temp_seq)))
) %>%
  left_join(species_to_host.raw, by = "Host")

plot_data.host.raw$Host <- factor(plot_data.host.raw$Host, levels = c("One", "Two", "Three"))
fitted_data.host.aegyptium.raw$Host <- factor(fitted_data.host.aegyptium.raw$Host, levels = c("One", "Two", "Three"))
# Plot observed data and fitted curves with facet wrap by species
p.host.aegyptium.raw <- ggplot() +
  geom_point(data = plot_data.host.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.aegyptium.raw, aes(x = Temperature, y = Fitted, color = Host), size = 1) +
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
p.host.aegyptium.raw

# Merge the observed and predicted data for comparison
merged_data.host.aegyptium.raw <- fitted_data.host.aegyptium.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host"))

# Calculate RMSE for each species
mse_by_host.aegyptium.raw <- merged_data.host.aegyptium.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.host.aegyptium.raw <- merged_data.host.aegyptium.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_host.aegyptium.raw,
     fitted_data.host.aegyptium.raw,
     merged_data.host.aegyptium.raw,
     host.raw.aegyptium.fit.mcmc,
     host.raw.aegyptium.fit,
     p.host.aegyptium.raw,
     file="aegyptium.raw.host.RData")

################################### By Host - Dromedarii ###################################
tick.rm.raw <- tick.raw[tick.raw$Species != "Hyalomma dromedarii", ]
host.rm.raw <- as.factor(tick.rm.raw$Host)
n.host.rm.raw <- length(unique(tick.rm.raw$Host))
N.obs.rm.raw <- length(tick.rm.raw$Trait)
unique(tick.rm.raw$Species)
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

  for (j in 1:n.host.rm.raw) {
    mu.host.c[j] ~ dexp(mu.c)       # Host-specific vertex y
    mu.host.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Host-specific vertex x
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
    host.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm.raw) {
    mu[i] <- mu.host.a[host.rm.raw[i]] * (mu.host.Topt[host.rm.raw[i]] - temp[i])^2 + mu.host.c[host.rm.raw[i]]
    trait[i] ~ dnorm(mu[i], host.tau[host.rm.raw[i]]) T(0, )   
  }}
", file = "host_robin.txt")


# Parameters to monitor
parameters <- c("mu.host.a", "mu.host.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt", "host.tau",
                "tau1", "mu.host.Topt")


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
jag.data <- list(trait = tick.rm.raw$Trait, N.obs.rm.raw = N.obs.rm.raw, 
                 temp = tick.rm.raw$Temp, 
                 n.host.rm.raw = n.host.rm.raw, host.rm.raw = host.rm.raw)
# Fit model
host.raw.dromedarii.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                             model.file = "host_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                             n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.raw.dromedarii.fit.mcmc <- as.mcmc(host.raw.dromedarii.fit) 
closeAllConnections()


# Get posterior means for parameters for each host group
summary_fit.host.rm.raw <- summary(host.raw.dromedarii.fit.mcmc)$statistics
a_host_mean.rm.raw <- summary_fit.host.rm.raw[grep("^mu.host.a\\[", rownames(summary_fit.host.rm.raw)), "Mean"]
c_host_mean.rm.raw <- summary_fit.host.rm.raw[grep("^mu.host.c\\[", rownames(summary_fit.host.rm.raw)), "Mean"]
Topt_host_mean.rm.raw <- summary_fit.host.rm.raw[grep("^mu.host.Topt\\[", rownames(summary_fit.host.rm.raw)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp) ,
                max(tick.raw$Temp) , by = 0.1)

# Get predicted values for each host group
host1.rm.raw <- function (x){
  a_host_mean.rm.raw[1] * (x - Topt_host_mean.rm.raw[1])^2 + c_host_mean.rm.raw[1]
}
host2.rm.raw <- function (x){
  a_host_mean.rm.raw[2] * (x - Topt_host_mean.rm.raw[2])^2 + c_host_mean.rm.raw[2]
}
host3.rm.raw <- function (x){
  a_host_mean.rm.raw[3] * (x - Topt_host_mean.rm.raw[3])^2 + c_host_mean.rm.raw[3]
}


# Data frame with temperature, fitted value, and host group
fitted_data.host.dromedarii.raw <- data.frame(Temperature = rep(temp_seq, n.host.rm.raw), 
                                           Fitted = c(host1.rm.raw(temp_seq),
                                                      host2.rm.raw(temp_seq),
                                                      host3.rm.raw(temp_seq)),
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
fitted_data.host.dromedarii.raw <- data.frame(
  Temperature = rep(temp_seq, n.host.rm.raw),
  Fitted = c(host1.rm.raw(temp_seq),
             host2.rm.raw(temp_seq),
             host3.rm.raw(temp_seq)),
  Host = factor(rep(c("One", "Three", "Two"), each = length(temp_seq)))
) %>%
  left_join(species_to_host.raw, by = "Host")

plot_data.host.raw$Host <- factor(plot_data.host.raw$Host, levels = c("One", "Two", "Three"))
fitted_data.host.dromedarii.raw$Host <- factor(fitted_data.host.dromedarii.raw$Host, levels = c("One", "Two", "Three"))
# Plot observed data and fitted curves with facet wrap by species
p.host.dromedarii.raw <- ggplot() +
  geom_point(data = plot_data.host.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.dromedarii.raw, aes(x = Temperature, y = Fitted, color = Host), size = 1) +
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
p.host.dromedarii.raw

# Merge the observed and predicted data for comparison
merged_data.host.dromedarii.raw <- fitted_data.host.dromedarii.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host"))

# Calculate RMSE for each species
mse_by_host.dromedarii.raw <- merged_data.host.dromedarii.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.host.dromedarii.raw <- merged_data.host.dromedarii.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_host.dromedarii.raw,
     fitted_data.host.dromedarii.raw,
     merged_data.host.dromedarii.raw,
     host.raw.dromedarii.fit.mcmc,
     host.raw.dromedarii.fit,
     p.host.dromedarii.raw,
     file="dromedarii.raw.host.RData")

################################### By Host - Impeltatum ###################################
tick.rm.raw <- tick.raw[tick.raw$Species != "Hyalomma impeltatum", ]
host.rm.raw <- as.factor(tick.rm.raw$Host)
n.host.rm.raw <- length(unique(tick.rm.raw$Host))
N.obs.rm.raw <- length(tick.rm.raw$Trait)
unique(tick.rm.raw$Species)
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

  for (j in 1:n.host.rm.raw) {
    mu.host.c[j] ~ dexp(mu.c)       # Host-specific vertex y
    mu.host.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Host-specific vertex x
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
    host.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm.raw) {
    mu[i] <- mu.host.a[host.rm.raw[i]] * (mu.host.Topt[host.rm.raw[i]] - temp[i])^2 + mu.host.c[host.rm.raw[i]]
    trait[i] ~ dnorm(mu[i], host.tau[host.rm.raw[i]]) T(0, )   
  }}
", file = "host_robin.txt")


# Parameters to monitor
parameters <- c("mu.host.a", "mu.host.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt", "host.tau",
                "tau1", "mu.host.Topt")


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
jag.data <- list(trait = tick.rm.raw$Trait, N.obs.rm.raw = N.obs.rm.raw, 
                 temp = tick.rm.raw$Temp, 
                 n.host.rm.raw = n.host.rm.raw, host.rm.raw = host.rm.raw)
# Fit model
host.raw.impeltatum.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                             model.file = "host_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                             n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.raw.impeltatum.fit.mcmc <- as.mcmc(host.raw.impeltatum.fit) 
closeAllConnections()


# Get posterior means for parameters for each host group
summary_fit.host.rm.raw <- summary(host.raw.impeltatum.fit.mcmc)$statistics
a_host_mean.rm.raw <- summary_fit.host.rm.raw[grep("^mu.host.a\\[", rownames(summary_fit.host.rm.raw)), "Mean"]
c_host_mean.rm.raw <- summary_fit.host.rm.raw[grep("^mu.host.c\\[", rownames(summary_fit.host.rm.raw)), "Mean"]
Topt_host_mean.rm.raw <- summary_fit.host.rm.raw[grep("^mu.host.Topt\\[", rownames(summary_fit.host.rm.raw)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp) ,
                max(tick.raw$Temp) , by = 0.1)

# Get predicted values for each host group
host1.rm.raw <- function (x){
  a_host_mean.rm.raw[1] * (x - Topt_host_mean.rm.raw[1])^2 + c_host_mean.rm.raw[1]
}
host2.rm.raw <- function (x){
  a_host_mean.rm.raw[2] * (x - Topt_host_mean.rm.raw[2])^2 + c_host_mean.rm.raw[2]
}
host3.rm.raw <- function (x){
  a_host_mean.rm.raw[3] * (x - Topt_host_mean.rm.raw[3])^2 + c_host_mean.rm.raw[3]
}


# Data frame with temperature, fitted value, and host group
fitted_data.host.impeltatum.raw <- data.frame(Temperature = rep(temp_seq, n.host.rm.raw), 
                                           Fitted = c(host1.rm.raw(temp_seq),
                                                      host2.rm.raw(temp_seq),
                                                      host3.rm.raw(temp_seq)),
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
fitted_data.host.impeltatum.raw <- data.frame(
  Temperature = rep(temp_seq, n.host.rm.raw),
  Fitted = c(host1.rm.raw(temp_seq),
             host2.rm.raw(temp_seq),
             host3.rm.raw(temp_seq)),
  Host = factor(rep(c("One", "Three", "Two"), each = length(temp_seq)))
) %>%
  left_join(species_to_host.raw, by = "Host")

plot_data.host.raw$Host <- factor(plot_data.host.raw$Host, levels = c("One", "Two", "Three"))
fitted_data.host.impeltatum.raw$Host <- factor(fitted_data.host.impeltatum.raw$Host, levels = c("One", "Two", "Three"))
# Plot observed data and fitted curves with facet wrap by species
p.host.impeltatum.raw <- ggplot() +
  geom_point(data = plot_data.host.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.impeltatum.raw, aes(x = Temperature, y = Fitted, color = Host), size = 1) +
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
p.host.impeltatum.raw

# Merge the observed and predicted data for comparison
merged_data.host.impeltatum.raw <- fitted_data.host.impeltatum.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host"))

# Calculate RMSE for each species
mse_by_host.impeltatum.raw <- merged_data.host.impeltatum.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.host.impeltatum.raw <- merged_data.host.impeltatum.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_host.impeltatum.raw,
     fitted_data.host.impeltatum.raw,
     merged_data.host.impeltatum.raw,
     host.raw.impeltatum.fit.mcmc,
     host.raw.impeltatum.fit,
     p.host.impeltatum.raw,
     file="impeltatum.raw.host.RData")

################################### By Host - Lusitanicum ###################################
tick.rm.raw <- tick.raw[tick.raw$Species != "Hyalomma lusitanicum", ]
host.rm.raw <- as.factor(tick.rm.raw$Host)
n.host.rm.raw <- length(unique(tick.rm.raw$Host))
N.obs.rm.raw <- length(tick.rm.raw$Trait)
unique(tick.rm.raw$Species)
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

  for (j in 1:n.host.rm.raw) {
    mu.host.c[j] ~ dexp(mu.c)       # Host-specific vertex y
    mu.host.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Host-specific vertex x
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
    host.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm.raw) {
    mu[i] <- mu.host.a[host.rm.raw[i]] * (mu.host.Topt[host.rm.raw[i]] - temp[i])^2 + mu.host.c[host.rm.raw[i]]
    trait[i] ~ dnorm(mu[i], host.tau[host.rm.raw[i]]) T(0, )   
  }}
", file = "host_robin.txt")


# Parameters to monitor
parameters <- c("mu.host.a", "mu.host.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt", "host.tau",
                "tau1", "mu.host.Topt")


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
jag.data <- list(trait = tick.rm.raw$Trait, N.obs.rm.raw = N.obs.rm.raw, 
                 temp = tick.rm.raw$Temp, 
                 n.host.rm.raw = n.host.rm.raw, host.rm.raw = host.rm.raw)
# Fit model
host.raw.lusitanicum.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                             model.file = "host_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                             n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.raw.lusitanicum.fit.mcmc <- as.mcmc(host.raw.lusitanicum.fit) 
closeAllConnections()


# Get posterior means for parameters for each host group
summary_fit.host.rm.raw <- summary(host.raw.lusitanicum.fit.mcmc)$statistics
a_host_mean.rm.raw <- summary_fit.host.rm.raw[grep("^mu.host.a\\[", rownames(summary_fit.host.rm.raw)), "Mean"]
c_host_mean.rm.raw <- summary_fit.host.rm.raw[grep("^mu.host.c\\[", rownames(summary_fit.host.rm.raw)), "Mean"]
Topt_host_mean.rm.raw <- summary_fit.host.rm.raw[grep("^mu.host.Topt\\[", rownames(summary_fit.host.rm.raw)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp) ,
                max(tick.raw$Temp) , by = 0.1)

# Get predicted values for each host group
host1.rm.raw <- function (x){
  a_host_mean.rm.raw[1] * (x - Topt_host_mean.rm.raw[1])^2 + c_host_mean.rm.raw[1]
}
host2.rm.raw <- function (x){
  a_host_mean.rm.raw[2] * (x - Topt_host_mean.rm.raw[2])^2 + c_host_mean.rm.raw[2]
}
host3.rm.raw <- function (x){
  a_host_mean.rm.raw[3] * (x - Topt_host_mean.rm.raw[3])^2 + c_host_mean.rm.raw[3]
}


# Data frame with temperature, fitted value, and host group
fitted_data.host.lusitanicum.raw <- data.frame(Temperature = rep(temp_seq, n.host.rm.raw), 
                                           Fitted = c(host1.rm.raw(temp_seq),
                                                      host2.rm.raw(temp_seq),
                                                      host3.rm.raw(temp_seq)),
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
fitted_data.host.lusitanicum.raw <- data.frame(
  Temperature = rep(temp_seq, n.host.rm.raw),
  Fitted = c(host1.rm.raw(temp_seq),
             host2.rm.raw(temp_seq),
             host3.rm.raw(temp_seq)),
  Host = factor(rep(c("One", "Three", "Two"), each = length(temp_seq)))
) %>%
  left_join(species_to_host.raw, by = "Host")

plot_data.host.raw$Host <- factor(plot_data.host.raw$Host, levels = c("One", "Two", "Three"))
fitted_data.host.lusitanicum.raw$Host <- factor(fitted_data.host.lusitanicum.raw$Host, levels = c("One", "Two", "Three"))
# Plot observed data and fitted curves with facet wrap by species
p.host.lusitanicum.raw <- ggplot() +
  geom_point(data = plot_data.host.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.lusitanicum.raw, aes(x = Temperature, y = Fitted, color = Host), size = 1) +
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
p.host.lusitanicum.raw

# Merge the observed and predicted data for comparison
merged_data.host.lusitanicum.raw <- fitted_data.host.lusitanicum.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host"))

# Calculate RMSE for each species
mse_by_host.lusitanicum.raw <- merged_data.host.lusitanicum.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.host.lusitanicum.raw <- merged_data.host.lusitanicum.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_host.lusitanicum.raw,
     fitted_data.host.lusitanicum.raw,
     merged_data.host.lusitanicum.raw,
     host.raw.lusitanicum.fit.mcmc,
     host.raw.lusitanicum.fit,
     p.host.lusitanicum.raw,
     file="lusitanicum.raw.host.RData")

################################### By Host - Schulzei ###################################
tick.rm.raw <- tick.raw[tick.raw$Species != "Hyalomma schulzei", ]
host.rm.raw <- as.factor(tick.rm.raw$Host)
n.host.rm.raw <- length(unique(tick.rm.raw$Host))
N.obs.rm.raw <- length(tick.rm.raw$Trait)
unique(tick.rm.raw$Species)
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

  for (j in 1:n.host.rm.raw) {
    mu.host.c[j] ~ dexp(mu.c)       # Host-specific vertex y
    mu.host.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Host-specific vertex x
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
    host.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm.raw) {
    mu[i] <- mu.host.a[host.rm.raw[i]] * (mu.host.Topt[host.rm.raw[i]] - temp[i])^2 + mu.host.c[host.rm.raw[i]]
    trait[i] ~ dnorm(mu[i], host.tau[host.rm.raw[i]]) T(0, )   
  }}
", file = "host_robin.txt")


# Parameters to monitor
parameters <- c("mu.host.a", "mu.host.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt", "host.tau",
                "tau1", "mu.host.Topt")


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
jag.data <- list(trait = tick.rm.raw$Trait, N.obs.rm.raw = N.obs.rm.raw, 
                 temp = tick.rm.raw$Temp, 
                 n.host.rm.raw = n.host.rm.raw, host.rm.raw = host.rm.raw)
# Fit model
host.raw.schulzei.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                             model.file = "host_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                             n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.raw.schulzei.fit.mcmc <- as.mcmc(host.raw.schulzei.fit) 
closeAllConnections()


# Get posterior means for parameters for each host group
summary_fit.host.rm.raw <- summary(host.raw.schulzei.fit.mcmc)$statistics
a_host_mean.rm.raw <- summary_fit.host.rm.raw[grep("^mu.host.a\\[", rownames(summary_fit.host.rm.raw)), "Mean"]
c_host_mean.rm.raw <- summary_fit.host.rm.raw[grep("^mu.host.c\\[", rownames(summary_fit.host.rm.raw)), "Mean"]
Topt_host_mean.rm.raw <- summary_fit.host.rm.raw[grep("^mu.host.Topt\\[", rownames(summary_fit.host.rm.raw)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.raw$Temp) ,
                max(tick.raw$Temp) , by = 0.1)

# Get predicted values for each host group
host1.rm.raw <- function (x){
  a_host_mean.rm.raw[1] * (x - Topt_host_mean.rm.raw[1])^2 + c_host_mean.rm.raw[1]
}
host2.rm.raw <- function (x){
  a_host_mean.rm.raw[2] * (x - Topt_host_mean.rm.raw[2])^2 + c_host_mean.rm.raw[2]
}
host3.rm.raw <- function (x){
  a_host_mean.rm.raw[3] * (x - Topt_host_mean.rm.raw[3])^2 + c_host_mean.rm.raw[3]
}


# Data frame with temperature, fitted value, and host group
fitted_data.host.schulzei.raw <- data.frame(Temperature = rep(temp_seq, n.host.rm.raw), 
                                           Fitted = c(host1.rm.raw(temp_seq),
                                                      host2.rm.raw(temp_seq),
                                                      host3.rm.raw(temp_seq)),
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
fitted_data.host.schulzei.raw <- data.frame(
  Temperature = rep(temp_seq, n.host.rm.raw),
  Fitted = c(host1.rm.raw(temp_seq),
             host2.rm.raw(temp_seq),
             host3.rm.raw(temp_seq)),
  Host = factor(rep(c("One", "Three", "Two"), each = length(temp_seq)))
) %>%
  left_join(species_to_host.raw, by = "Host")

plot_data.host.raw$Host <- factor(plot_data.host.raw$Host, levels = c("One", "Two", "Three"))
fitted_data.host.schulzei.raw$Host <- factor(fitted_data.host.schulzei.raw$Host, levels = c("One", "Two", "Three"))
# Plot observed data and fitted curves with facet wrap by species
p.host.schulzei.raw <- ggplot() +
  geom_point(data = plot_data.host.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.schulzei.raw, aes(x = Temperature, y = Fitted, color = Host), size = 1) +
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
p.host.schulzei.raw

# Merge the observed and predicted data for comparison
merged_data.host.schulzei.raw <- fitted_data.host.schulzei.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host"))

# Calculate RMSE for each species
mse_by_host.schulzei.raw <- merged_data.host.schulzei.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.host.schulzei.raw <- merged_data.host.schulzei.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_host.schulzei.raw,
     fitted_data.host.schulzei.raw,
     merged_data.host.schulzei.raw,
     host.raw.schulzei.fit.mcmc,
     host.raw.schulzei.fit,
     p.host.schulzei.raw,
     file="schulzei.raw.host.RData")


################################### By Host and Climate- Leporispalustris ###################################
tick.rm.raw <- tick.raw[tick.raw$Species != "Haemaphysalis leporispalustris", ]
host.rm.raw <- as.factor(tick.rm.raw$Host)
n.host.rm.raw <- length(unique(tick.rm.raw$Host))
climate.rm.raw <- as.factor(tick.rm.raw$Climate)
n.climate.rm.raw <- length(unique(tick.rm.raw$Climate))
N.obs.rm.raw <- length(tick.rm.raw$Trait)
unique(tick.rm.raw$Species)

sink("host.climate_robin.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Host-specific
  mu.a ~ dexp(0.01)            # Mean for host effects a (width)
  mu.c ~ dexp(10)          # Mean for host effects c (vertex y)
  # Climate-specific
  mu.Topt ~ dnorm(30, 0.1)          # Mean for climate effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for climate effect Topt (vertex x)

  for (j in 1:n.host.rm.raw) {
    mu.host.c[j] ~ dexp(mu.c)       # Host-specific vertex y
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
  }
  for (k in 1:n.climate.rm.raw) {
    mu.climate.Topt[k] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific vertex x
    climate.tau[k] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm.raw) {
    mu[i] <- mu.host.a[host.rm.raw[i]] * (mu.climate.Topt[climate.rm.raw[i]] - temp[i])^2 + mu.host.c[host.rm.raw[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[climate.rm.raw[i]]) T(0, )
  }}
", file = "host_robin.txt")


# Parameters to monitor
parameters <- c("mu.host.a", "mu.host.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt", "climate.tau",
                "tau1", "mu.climate.Topt")


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
jag.data <- list(trait = tick.rm.raw$Trait, N.obs.rm.raw = N.obs.rm.raw, 
                 temp = tick.rm.raw$Temp, 
                 n.host.rm.raw = n.host.rm.raw, host.rm.raw = host.rm.raw, 
                 n.climate.rm.raw = n.climate.rm.raw, climate.rm.raw = climate.rm.raw)
# Fit model
host.climate.raw.leporispalustris.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                                              model.file = "host_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                                              n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.climate.raw.leporispalustris.fit.mcmc <- as.mcmc(host.climate.raw.leporispalustris.fit) 
closeAllConnections()

# Get posterior means for parameters for each host group
summary_fit.host.climate.raw <- summary(host.climate.raw.leporispalustris.fit.mcmc)$statistics
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
fitted_data.host.climate.leporispalustris.raw <- data.frame(Temperature = rep(temp_seq, n.host.rm.raw * n.climate.rm.raw), 
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
                                                            Host = factor(rep(rep(c("One", "Three", "Two"), each = length(temp_seq)), n.climate.rm.raw)),
                                                            Climate = factor(rep(c("Mixed", "Subtropical", "Temperate", "Tropical"),
                                                                                 each = length(temp_seq) * n.host.rm.raw)))
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
fitted_data.host.climate.leporispalustris.raw <- data.frame(Temperature = rep(temp_seq, n.host.rm.raw * n.climate.rm.raw), 
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
                                                            Host = factor(rep(rep(c("One", "Three", "Two"), each = length(temp_seq)), n.climate.rm.raw)),
                                                            Climate = factor(rep(c("Mixed", "Subtropical", "Temperate", "Tropical"),
                                                                                 each = length(temp_seq) * n.host.rm.raw))) %>%
  left_join(species_to_host.climate.raw, by = c("Host", "Climate"))|> 
  filter(!is.na(Species)) 

plot_data.host.climate.raw$Host <- factor(plot_data.host.climate.raw$Host, levels = c("One", "Two", "Three"))
fitted_data.host.climate.leporispalustris.raw$Host <- factor(fitted_data.host.climate.leporispalustris.raw$Host, levels = c("One", "Two", "Three"))

# Combine Climate and Size for coloring
fitted_data.host.climate.leporispalustris.raw$`Host, Climate` <- with(fitted_data.host.climate.leporispalustris.raw, paste(Host, Climate, sep = ", "))
fitted_data.host.climate.leporispalustris.raw$`Host, Climate` <- factor(fitted_data.host.climate.leporispalustris.raw$`Host, Climate`, levels = c("One, Tropical", 
                                                                                                                              "Two, Subtropical", "Three, Mixed",
                                                                                                                              "Three, Subtropical", "Three, Temperate",
                                                                                                                              "Three, Tropical"))

# Plot observed data and fitted curves with facet wrap by species
p.host.climate.leporispalustris.raw <- ggplot() +
  geom_point(data = plot_data.host.climate.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.climate.leporispalustris.raw, aes(x = Temperature, y = Fitted, color = `Host, Climate`), size = 1) +
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
p.host.climate.leporispalustris.raw

# Merge the observed and predicted data for comparison
merged_data.host.climate.leporispalustris.raw <- fitted_data.host.climate.leporispalustris.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host", "Climate" = "Climate"))

# Calculate RMSE for each species
mse_by_host.climate.leporispalustris.raw <- merged_data.host.climate.leporispalustris.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.host.climate.leporispalustris.raw <- merged_data.host.climate.leporispalustris.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 


save(mse_by_host.climate.leporispalustris.raw,
     fitted_data.host.climate.leporispalustris.raw,
     merged_data.host.climate.leporispalustris.raw,
     host.climate.raw.leporispalustris.fit.mcmc,
     host.climate.raw.leporispalustris.fit,
     p.host.climate.leporispalustris.raw,
     file="host.climate.raw.leporispalustris.fit.RData")

################################### By Host and Climate- Aegyptium ###################################
tick.rm.raw <- tick.raw[tick.raw$Species != "Hyalomma aegyptium", ]
host.rm.raw <- as.factor(tick.rm.raw$Host)
n.host.rm.raw <- length(unique(tick.rm.raw$Host))
climate.rm.raw <- as.factor(tick.rm.raw$Climate)
n.climate.rm.raw <- length(unique(tick.rm.raw$Climate))
N.obs.rm.raw <- length(tick.rm.raw$Trait)
unique(tick.rm.raw$Species)

sink("host.climate_robin.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Host-specific
  mu.a ~ dexp(0.01)            # Mean for host effects a (width)
  mu.c ~ dexp(10)          # Mean for host effects c (vertex y)
  # Climate-specific
  mu.Topt ~ dnorm(30, 0.1)          # Mean for climate effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for climate effect Topt (vertex x)

  for (j in 1:n.host.rm.raw) {
    mu.host.c[j] ~ dexp(mu.c)       # Host-specific vertex y
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
  }
  for (k in 1:n.climate.rm.raw) {
    mu.climate.Topt[k] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific vertex x
    climate.tau[k] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm.raw) {
    mu[i] <- mu.host.a[host.rm.raw[i]] * (mu.climate.Topt[climate.rm.raw[i]] - temp[i])^2 + mu.host.c[host.rm.raw[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[climate.rm.raw[i]]) T(0, )
  }}
", file = "host_robin.txt")


# Parameters to monitor
parameters <- c("mu.host.a", "mu.host.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt", "climate.tau",
                "tau1", "mu.climate.Topt")


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
jag.data <- list(trait = tick.rm.raw$Trait, N.obs.rm.raw = N.obs.rm.raw, 
                 temp = tick.rm.raw$Temp, 
                 n.host.rm.raw = n.host.rm.raw, host.rm.raw = host.rm.raw, 
                 n.climate.rm.raw = n.climate.rm.raw, climate.rm.raw = climate.rm.raw)
# Fit model
host.climate.raw.aegyptium.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                                              model.file = "host_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                                              n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.climate.raw.aegyptium.fit.mcmc <- as.mcmc(host.climate.raw.aegyptium.fit) 
closeAllConnections()

# Get posterior means for parameters for each host group
summary_fit.host.climate.raw <- summary(host.climate.raw.aegyptium.fit.mcmc)$statistics
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
fitted_data.host.climate.aegyptium.raw <- data.frame(Temperature = rep(temp_seq, n.host.rm.raw * n.climate.rm.raw), 
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
                                                            Host = factor(rep(rep(c("One", "Three", "Two"), each = length(temp_seq)), n.climate.rm.raw)),
                                                            Climate = factor(rep(c("Mixed", "Subtropical", "Temperate", "Tropical"),
                                                                                 each = length(temp_seq) * n.host.rm.raw)))
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
fitted_data.host.climate.aegyptium.raw <- data.frame(Temperature = rep(temp_seq, n.host.rm.raw * n.climate.rm.raw), 
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
                                                            Host = factor(rep(rep(c("One", "Three", "Two"), each = length(temp_seq)), n.climate.rm.raw)),
                                                            Climate = factor(rep(c("Mixed", "Subtropical", "Temperate", "Tropical"),
                                                                                 each = length(temp_seq) * n.host.rm.raw))) %>%
  left_join(species_to_host.climate.raw, by = c("Host", "Climate"))|> 
  filter(!is.na(Species)) 

plot_data.host.climate.raw$Host <- factor(plot_data.host.climate.raw$Host, levels = c("One", "Two", "Three"))
fitted_data.host.climate.aegyptium.raw$Host <- factor(fitted_data.host.climate.aegyptium.raw$Host, levels = c("One", "Two", "Three"))

# Combine Climate and Size for coloring
fitted_data.host.climate.aegyptium.raw$`Host, Climate` <- with(fitted_data.host.climate.aegyptium.raw, paste(Host, Climate, sep = ", "))
fitted_data.host.climate.aegyptium.raw$`Host, Climate` <- factor(fitted_data.host.climate.aegyptium.raw$`Host, Climate`, levels = c("One, Tropical", 
                                                                                                                              "Two, Subtropical", "Three, Mixed",
                                                                                                                              "Three, Subtropical", "Three, Temperate",
                                                                                                                              "Three, Tropical"))

# Plot observed data and fitted curves with facet wrap by species
p.host.climate.aegyptium.raw <- ggplot() +
  geom_point(data = plot_data.host.climate.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.climate.aegyptium.raw, aes(x = Temperature, y = Fitted, color = `Host, Climate`), size = 1) +
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
p.host.climate.aegyptium.raw

# Merge the observed and predicted data for comparison
merged_data.host.climate.aegyptium.raw <- fitted_data.host.climate.aegyptium.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host", "Climate" = "Climate"))

# Calculate RMSE for each species
mse_by_host.climate.aegyptium.raw <- merged_data.host.climate.aegyptium.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.host.climate.aegyptium.raw <- merged_data.host.climate.aegyptium.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 


save(mse_by_host.climate.aegyptium.raw,
     fitted_data.host.climate.aegyptium.raw,
     merged_data.host.climate.aegyptium.raw,
     host.climate.raw.aegyptium.fit.mcmc,
     host.climate.raw.aegyptium.fit,
     p.host.climate.aegyptium.raw,
     file="host.climate.raw.aegyptium.fit.RData")

################################### By Host and Climate- Dromedarii ###################################
tick.rm.raw <- tick.raw[tick.raw$Species != "Hyalomma dromedarii", ]
host.rm.raw <- as.factor(tick.rm.raw$Host)
n.host.rm.raw <- length(unique(tick.rm.raw$Host))
climate.rm.raw <- as.factor(tick.rm.raw$Climate)
n.climate.rm.raw <- length(unique(tick.rm.raw$Climate))
N.obs.rm.raw <- length(tick.rm.raw$Trait)
unique(tick.rm.raw$Species)

sink("host.climate_robin.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Host-specific
  mu.a ~ dexp(0.01)            # Mean for host effects a (width)
  mu.c ~ dexp(10)          # Mean for host effects c (vertex y)
  # Climate-specific
  mu.Topt ~ dnorm(30, 0.1)          # Mean for climate effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for climate effect Topt (vertex x)

  for (j in 1:n.host.rm.raw) {
    mu.host.c[j] ~ dexp(mu.c)       # Host-specific vertex y
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
  }
  for (k in 1:n.climate.rm.raw) {
    mu.climate.Topt[k] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific vertex x
    climate.tau[k] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm.raw) {
    mu[i] <- mu.host.a[host.rm.raw[i]] * (mu.climate.Topt[climate.rm.raw[i]] - temp[i])^2 + mu.host.c[host.rm.raw[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[climate.rm.raw[i]]) T(0, )
  }}
", file = "host_robin.txt")


# Parameters to monitor
parameters <- c("mu.host.a", "mu.host.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt", "climate.tau",
                "tau1", "mu.climate.Topt")


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
jag.data <- list(trait = tick.rm.raw$Trait, N.obs.rm.raw = N.obs.rm.raw, 
                 temp = tick.rm.raw$Temp, 
                 n.host.rm.raw = n.host.rm.raw, host.rm.raw = host.rm.raw, 
                 n.climate.rm.raw = n.climate.rm.raw, climate.rm.raw = climate.rm.raw)
# Fit model
host.climate.raw.dromedarii.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                                              model.file = "host_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                                              n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.climate.raw.dromedarii.fit.mcmc <- as.mcmc(host.climate.raw.dromedarii.fit) 
closeAllConnections()

# Get posterior means for parameters for each host group
summary_fit.host.climate.raw <- summary(host.climate.raw.dromedarii.fit.mcmc)$statistics
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
fitted_data.host.climate.dromedarii.raw <- data.frame(Temperature = rep(temp_seq, n.host.rm.raw * n.climate.rm.raw), 
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
                                                            Host = factor(rep(rep(c("One", "Three", "Two"), each = length(temp_seq)), n.climate.rm.raw)),
                                                            Climate = factor(rep(c("Mixed", "Subtropical", "Temperate", "Tropical"),
                                                                                 each = length(temp_seq) * n.host.rm.raw)))
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
fitted_data.host.climate.dromedarii.raw <- data.frame(Temperature = rep(temp_seq, n.host.rm.raw * n.climate.rm.raw), 
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
                                                            Host = factor(rep(rep(c("One", "Three", "Two"), each = length(temp_seq)), n.climate.rm.raw)),
                                                            Climate = factor(rep(c("Mixed", "Subtropical", "Temperate", "Tropical"),
                                                                                 each = length(temp_seq) * n.host.rm.raw))) %>%
  left_join(species_to_host.climate.raw, by = c("Host", "Climate"))|> 
  filter(!is.na(Species)) 

plot_data.host.climate.raw$Host <- factor(plot_data.host.climate.raw$Host, levels = c("One", "Two", "Three"))
fitted_data.host.climate.dromedarii.raw$Host <- factor(fitted_data.host.climate.dromedarii.raw$Host, levels = c("One", "Two", "Three"))

# Combine Climate and Size for coloring
fitted_data.host.climate.dromedarii.raw$`Host, Climate` <- with(fitted_data.host.climate.dromedarii.raw, paste(Host, Climate, sep = ", "))
fitted_data.host.climate.dromedarii.raw$`Host, Climate` <- factor(fitted_data.host.climate.dromedarii.raw$`Host, Climate`, levels = c("One, Tropical", 
                                                                                                                              "Two, Subtropical", "Three, Mixed",
                                                                                                                              "Three, Subtropical", "Three, Temperate",
                                                                                                                              "Three, Tropical"))

# Plot observed data and fitted curves with facet wrap by species
p.host.climate.dromedarii.raw <- ggplot() +
  geom_point(data = plot_data.host.climate.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.climate.dromedarii.raw, aes(x = Temperature, y = Fitted, color = `Host, Climate`), size = 1) +
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
p.host.climate.dromedarii.raw

# Merge the observed and predicted data for comparison
merged_data.host.climate.dromedarii.raw <- fitted_data.host.climate.dromedarii.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host", "Climate" = "Climate"))

# Calculate RMSE for each species
mse_by_host.climate.dromedarii.raw <- merged_data.host.climate.dromedarii.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.host.climate.dromedarii.raw <- merged_data.host.climate.dromedarii.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 


save(mse_by_host.climate.dromedarii.raw,
     fitted_data.host.climate.dromedarii.raw,
     merged_data.host.climate.dromedarii.raw,
     host.climate.raw.dromedarii.fit.mcmc,
     host.climate.raw.dromedarii.fit,
     p.host.climate.dromedarii.raw,
     file="host.climate.raw.dromedarii.fit.RData")

################################### By Host and Climate- Impeltatum ###################################
tick.rm.raw <- tick.raw[tick.raw$Species != "Hyalomma impeltatum", ]
host.rm.raw <- as.factor(tick.rm.raw$Host)
n.host.rm.raw <- length(unique(tick.rm.raw$Host))
climate.rm.raw <- as.factor(tick.rm.raw$Climate)
n.climate.rm.raw <- length(unique(tick.rm.raw$Climate))
N.obs.rm.raw <- length(tick.rm.raw$Trait)
unique(tick.rm.raw$Species)

sink("host.climate_robin.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Host-specific
  mu.a ~ dexp(0.01)            # Mean for host effects a (width)
  mu.c ~ dexp(10)          # Mean for host effects c (vertex y)
  # Climate-specific
  mu.Topt ~ dnorm(30, 0.1)          # Mean for climate effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for climate effect Topt (vertex x)

  for (j in 1:n.host.rm.raw) {
    mu.host.c[j] ~ dexp(mu.c)       # Host-specific vertex y
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
  }
  for (k in 1:n.climate.rm.raw) {
    mu.climate.Topt[k] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific vertex x
    climate.tau[k] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm.raw) {
    mu[i] <- mu.host.a[host.rm.raw[i]] * (mu.climate.Topt[climate.rm.raw[i]] - temp[i])^2 + mu.host.c[host.rm.raw[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[climate.rm.raw[i]]) T(0, )
  }}
", file = "host_robin.txt")


# Parameters to monitor
parameters <- c("mu.host.a", "mu.host.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt", "climate.tau",
                "tau1", "mu.climate.Topt")


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
jag.data <- list(trait = tick.rm.raw$Trait, N.obs.rm.raw = N.obs.rm.raw, 
                 temp = tick.rm.raw$Temp, 
                 n.host.rm.raw = n.host.rm.raw, host.rm.raw = host.rm.raw, 
                 n.climate.rm.raw = n.climate.rm.raw, climate.rm.raw = climate.rm.raw)
# Fit model
host.climate.raw.impeltatum.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                                              model.file = "host_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                                              n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.climate.raw.impeltatum.fit.mcmc <- as.mcmc(host.climate.raw.impeltatum.fit) 
closeAllConnections()

# Get posterior means for parameters for each host group
summary_fit.host.climate.raw <- summary(host.climate.raw.impeltatum.fit.mcmc)$statistics
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
fitted_data.host.climate.impeltatum.raw <- data.frame(Temperature = rep(temp_seq, n.host.rm.raw * n.climate.rm.raw), 
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
                                                            Host = factor(rep(rep(c("One", "Three", "Two"), each = length(temp_seq)), n.climate.rm.raw)),
                                                            Climate = factor(rep(c("Mixed", "Subtropical", "Temperate", "Tropical"),
                                                                                 each = length(temp_seq) * n.host.rm.raw)))
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
fitted_data.host.climate.impeltatum.raw <- data.frame(Temperature = rep(temp_seq, n.host.rm.raw * n.climate.rm.raw), 
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
                                                            Host = factor(rep(rep(c("One", "Three", "Two"), each = length(temp_seq)), n.climate.rm.raw)),
                                                            Climate = factor(rep(c("Mixed", "Subtropical", "Temperate", "Tropical"),
                                                                                 each = length(temp_seq) * n.host.rm.raw))) %>%
  left_join(species_to_host.climate.raw, by = c("Host", "Climate"))|> 
  filter(!is.na(Species)) 

plot_data.host.climate.raw$Host <- factor(plot_data.host.climate.raw$Host, levels = c("One", "Two", "Three"))
fitted_data.host.climate.impeltatum.raw$Host <- factor(fitted_data.host.climate.impeltatum.raw$Host, levels = c("One", "Two", "Three"))

# Combine Climate and Size for coloring
fitted_data.host.climate.impeltatum.raw$`Host, Climate` <- with(fitted_data.host.climate.impeltatum.raw, paste(Host, Climate, sep = ", "))
fitted_data.host.climate.impeltatum.raw$`Host, Climate` <- factor(fitted_data.host.climate.impeltatum.raw$`Host, Climate`, levels = c("One, Tropical", 
                                                                                                                              "Two, Subtropical", "Three, Mixed",
                                                                                                                              "Three, Subtropical", "Three, Temperate",
                                                                                                                              "Three, Tropical"))

# Plot observed data and fitted curves with facet wrap by species
p.host.climate.impeltatum.raw <- ggplot() +
  geom_point(data = plot_data.host.climate.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.climate.impeltatum.raw, aes(x = Temperature, y = Fitted, color = `Host, Climate`), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Host and Climate",
       x = "Temperature", y = "Oviposition Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )   
p.host.climate.impeltatum.raw

# Merge the observed and predicted data for comparison
merged_data.host.climate.impeltatum.raw <- fitted_data.host.climate.impeltatum.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host", "Climate" = "Climate"))

# Calculate RMSE for each species
mse_by_host.climate.impeltatum.raw <- merged_data.host.climate.impeltatum.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.host.climate.impeltatum.raw <- merged_data.host.climate.impeltatum.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 


save(mse_by_host.climate.impeltatum.raw,
     fitted_data.host.climate.impeltatum.raw,
     merged_data.host.climate.impeltatum.raw,
     host.climate.raw.impeltatum.fit.mcmc,
     host.climate.raw.impeltatum.fit,
     p.host.climate.impeltatum.raw,
     file="host.climate.raw.impeltatum.fit.RData")

################################### By Host and Climate- Lusitanicum ###################################
tick.rm.raw <- tick.raw[tick.raw$Species != "Hyalomma lusitanicum", ]
host.rm.raw <- as.factor(tick.rm.raw$Host)
n.host.rm.raw <- length(unique(tick.rm.raw$Host))
climate.rm.raw <- as.factor(tick.rm.raw$Climate)
n.climate.rm.raw <- length(unique(tick.rm.raw$Climate))
N.obs.rm.raw <- length(tick.rm.raw$Trait)
unique(tick.rm.raw$Species)

sink("host.climate_robin.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Host-specific
  mu.a ~ dexp(0.01)            # Mean for host effects a (width)
  mu.c ~ dexp(10)          # Mean for host effects c (vertex y)
  # Climate-specific
  mu.Topt ~ dnorm(30, 0.1)          # Mean for climate effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for climate effect Topt (vertex x)

  for (j in 1:n.host.rm.raw) {
    mu.host.c[j] ~ dexp(mu.c)       # Host-specific vertex y
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
  }
  for (k in 1:n.climate.rm.raw) {
    mu.climate.Topt[k] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific vertex x
    climate.tau[k] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm.raw) {
    mu[i] <- mu.host.a[host.rm.raw[i]] * (mu.climate.Topt[climate.rm.raw[i]] - temp[i])^2 + mu.host.c[host.rm.raw[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[climate.rm.raw[i]]) T(0, )
  }}
", file = "host_robin.txt")


# Parameters to monitor
parameters <- c("mu.host.a", "mu.host.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt", "climate.tau",
                "tau1", "mu.climate.Topt")


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
jag.data <- list(trait = tick.rm.raw$Trait, N.obs.rm.raw = N.obs.rm.raw, 
                 temp = tick.rm.raw$Temp, 
                 n.host.rm.raw = n.host.rm.raw, host.rm.raw = host.rm.raw, 
                 n.climate.rm.raw = n.climate.rm.raw, climate.rm.raw = climate.rm.raw)
# Fit model
host.climate.raw.lusitanicum.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                                              model.file = "host_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                                              n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.climate.raw.lusitanicum.fit.mcmc <- as.mcmc(host.climate.raw.lusitanicum.fit) 
closeAllConnections()

# Get posterior means for parameters for each host group
summary_fit.host.climate.raw <- summary(host.climate.raw.lusitanicum.fit.mcmc)$statistics
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
fitted_data.host.climate.lusitanicum.raw <- data.frame(Temperature = rep(temp_seq, n.host.rm.raw * n.climate.rm.raw), 
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
                                                            Host = factor(rep(rep(c("One", "Three", "Two"), each = length(temp_seq)), n.climate.rm.raw)),
                                                            Climate = factor(rep(c("Mixed", "Subtropical", "Temperate", "Tropical"),
                                                                                 each = length(temp_seq) * n.host.rm.raw)))
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
fitted_data.host.climate.lusitanicum.raw <- data.frame(Temperature = rep(temp_seq, n.host.rm.raw * n.climate.rm.raw), 
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
                                                            Host = factor(rep(rep(c("One", "Three", "Two"), each = length(temp_seq)), n.climate.rm.raw)),
                                                            Climate = factor(rep(c("Mixed", "Subtropical", "Temperate", "Tropical"),
                                                                                 each = length(temp_seq) * n.host.rm.raw))) %>%
  left_join(species_to_host.climate.raw, by = c("Host", "Climate"))|> 
  filter(!is.na(Species)) 

plot_data.host.climate.raw$Host <- factor(plot_data.host.climate.raw$Host, levels = c("One", "Two", "Three"))
fitted_data.host.climate.lusitanicum.raw$Host <- factor(fitted_data.host.climate.lusitanicum.raw$Host, levels = c("One", "Two", "Three"))

# Combine Climate and Size for coloring
fitted_data.host.climate.lusitanicum.raw$`Host, Climate` <- with(fitted_data.host.climate.lusitanicum.raw, paste(Host, Climate, sep = ", "))
fitted_data.host.climate.lusitanicum.raw$`Host, Climate` <- factor(fitted_data.host.climate.lusitanicum.raw$`Host, Climate`, levels = c("One, Tropical", 
                                                                                                                              "Two, Subtropical", "Three, Mixed",
                                                                                                                              "Three, Subtropical", "Three, Temperate",
                                                                                                                              "Three, Tropical"))

# Plot observed data and fitted curves with facet wrap by species
p.host.climate.lusitanicum.raw <- ggplot() +
  geom_point(data = plot_data.host.climate.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.climate.lusitanicum.raw, aes(x = Temperature, y = Fitted, color = `Host, Climate`), size = 1) +
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
p.host.climate.lusitanicum.raw

# Merge the observed and predicted data for comparison
merged_data.host.climate.lusitanicum.raw <- fitted_data.host.climate.lusitanicum.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host", "Climate" = "Climate"))

# Calculate RMSE for each species
mse_by_host.climate.lusitanicum.raw <- merged_data.host.climate.lusitanicum.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.host.climate.lusitanicum.raw <- merged_data.host.climate.lusitanicum.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 


save(mse_by_host.climate.lusitanicum.raw,
     fitted_data.host.climate.lusitanicum.raw,
     merged_data.host.climate.lusitanicum.raw,
     host.climate.raw.lusitanicum.fit.mcmc,
     host.climate.raw.lusitanicum.fit,
     p.host.climate.lusitanicum.raw,
     file="host.climate.raw.lusitanicum.fit.RData")

################################### By Host and Climate- Schulzei ###################################
tick.rm.raw <- tick.raw[tick.raw$Species != "Hyalomma schulzei", ]
host.rm.raw <- as.factor(tick.rm.raw$Host)
n.host.rm.raw <- length(unique(tick.rm.raw$Host))
climate.rm.raw <- as.factor(tick.rm.raw$Climate)
n.climate.rm.raw <- length(unique(tick.rm.raw$Climate))
N.obs.rm.raw <- length(tick.rm.raw$Trait)
unique(tick.rm.raw$Species)

sink("host.climate_robin.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Host-specific
  mu.a ~ dexp(0.01)            # Mean for host effects a (width)
  mu.c ~ dexp(10)          # Mean for host effects c (vertex y)
  # Climate-specific
  mu.Topt ~ dnorm(30, 0.1)          # Mean for climate effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for climate effect Topt (vertex x)

  for (j in 1:n.host.rm.raw) {
    mu.host.c[j] ~ dexp(mu.c)       # Host-specific vertex y
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
  }
  for (k in 1:n.climate.rm.raw) {
    mu.climate.Topt[k] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific vertex x
    climate.tau[k] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm.raw) {
    mu[i] <- mu.host.a[host.rm.raw[i]] * (mu.climate.Topt[climate.rm.raw[i]] - temp[i])^2 + mu.host.c[host.rm.raw[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[climate.rm.raw[i]]) T(0, )
  }}
", file = "host_robin.txt")


# Parameters to monitor
parameters <- c("mu.host.a", "mu.host.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt", "climate.tau",
                "tau1", "mu.climate.Topt")


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
jag.data <- list(trait = tick.rm.raw$Trait, N.obs.rm.raw = N.obs.rm.raw, 
                 temp = tick.rm.raw$Temp, 
                 n.host.rm.raw = n.host.rm.raw, host.rm.raw = host.rm.raw, 
                 n.climate.rm.raw = n.climate.rm.raw, climate.rm.raw = climate.rm.raw)
# Fit model
host.climate.raw.schulzei.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                                              model.file = "host_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                                              n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.climate.raw.schulzei.fit.mcmc <- as.mcmc(host.climate.raw.schulzei.fit) 
closeAllConnections()

# Get posterior means for parameters for each host group
summary_fit.host.climate.raw <- summary(host.climate.raw.schulzei.fit.mcmc)$statistics
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
fitted_data.host.climate.schulzei.raw <- data.frame(Temperature = rep(temp_seq, n.host.rm.raw * n.climate.rm.raw), 
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
                                                            Host = factor(rep(rep(c("One", "Three", "Two"), each = length(temp_seq)), n.climate.rm.raw)),
                                                            Climate = factor(rep(c("Mixed", "Subtropical", "Temperate", "Tropical"),
                                                                                 each = length(temp_seq) * n.host.rm.raw)))
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
fitted_data.host.climate.schulzei.raw <- data.frame(Temperature = rep(temp_seq, n.host.rm.raw * n.climate.rm.raw), 
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
                                                            Host = factor(rep(rep(c("One", "Three", "Two"), each = length(temp_seq)), n.climate.rm.raw)),
                                                            Climate = factor(rep(c("Mixed", "Subtropical", "Temperate", "Tropical"),
                                                                                 each = length(temp_seq) * n.host.rm.raw))) %>%
  left_join(species_to_host.climate.raw, by = c("Host", "Climate"))|> 
  filter(!is.na(Species)) 

plot_data.host.climate.raw$Host <- factor(plot_data.host.climate.raw$Host, levels = c("One", "Two", "Three"))
fitted_data.host.climate.schulzei.raw$Host <- factor(fitted_data.host.climate.schulzei.raw$Host, levels = c("One", "Two", "Three"))

# Combine Climate and Size for coloring
fitted_data.host.climate.schulzei.raw$`Host, Climate` <- with(fitted_data.host.climate.schulzei.raw, paste(Host, Climate, sep = ", "))
fitted_data.host.climate.schulzei.raw$`Host, Climate` <- factor(fitted_data.host.climate.schulzei.raw$`Host, Climate`, levels = c("One, Tropical", 
                                                                                                                              "Two, Subtropical", "Three, Mixed",
                                                                                                                              "Three, Subtropical", "Three, Temperate",
                                                                                                                              "Three, Tropical"))

# Plot observed data and fitted curves with facet wrap by species
p.host.climate.schulzei.raw <- ggplot() +
  geom_point(data = plot_data.host.climate.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.climate.schulzei.raw, aes(x = Temperature, y = Fitted, color = `Host, Climate`), size = 1) +
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
p.host.climate.schulzei.raw

# Merge the observed and predicted data for comparison
merged_data.host.climate.schulzei.raw <- fitted_data.host.climate.schulzei.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host", "Climate" = "Climate"))

# Calculate RMSE for each species
mse_by_host.climate.schulzei.raw <- merged_data.host.climate.schulzei.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.host.climate.schulzei.raw <- merged_data.host.climate.schulzei.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 


save(mse_by_host.climate.schulzei.raw,
     fitted_data.host.climate.schulzei.raw,
     merged_data.host.climate.schulzei.raw,
     host.climate.raw.schulzei.fit.mcmc,
     host.climate.raw.schulzei.fit,
     p.host.climate.schulzei.raw,
     file="host.climate.raw.schulzei.fit.RData")

################################### By Host, Genus, and Climate- Aegyptium ###################################
# Model by typical number of host species
tick.rm.raw <- tick.raw[tick.raw$Species != "Hyalomma aegyptium", ]
climate.rm.raw <- as.factor(tick.rm.raw$Climate)
n.climate.rm.raw <- length(unique(tick.rm.raw$Climate))
genus.rm.raw <- as.factor(tick.rm.raw$Genus)
n.genus.rm.raw <- length(unique(tick.rm.raw$Genus))
host.rm.raw <- as.factor(tick.rm.raw$Host)
n.host.rm.raw <- length(unique(tick.rm.raw$Host))
N.obs.rm.raw <- length(tick.rm.raw$Trait)
unique(tick.rm.raw$Species)
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

  for (j in 1:n.host.rm.raw) {
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
  }
  for (k in 1:n.climate.rm.raw) {
    mu.climate.Topt[k] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific vertex x
    climate.tau[k] ~ dexp(tau1)
  }
  for (s in 1:n.genus.rm.raw) {
    mu.genus.c[s] ~ dexp(mu.c)       # genus-specific vertex y
  }
  # Likelihood
  for (i in 1:N.obs.rm.raw) {
    mu[i] <- mu.host.a[host.rm.raw[i]] * (mu.climate.Topt[climate.rm.raw[i]] - temp[i])^2 + mu.genus.c[genus.rm.raw[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[climate.rm.raw[i]]) T(0, ) 
  }}
", file = "host.climate.genus_robin.txt")


# Parameters to monitor
parameters <- c("mu.host.a", "mu.genus.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt", "climate.tau",
                "tau1", "mu.climate.Topt")


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
jag.data <- list(trait = tick.rm.raw$Trait, N.obs.rm.raw = N.obs.rm.raw, 
                 temp = tick.rm.raw$Temp, 
                 n.host.rm.raw = n.host.rm.raw, host.rm.raw = host.rm.raw, 
                 n.climate.rm.raw = n.climate.rm.raw, climate.rm.raw = climate.rm.raw,
                 n.genus.rm.raw = n.genus.rm.raw, genus.rm.raw = genus.rm.raw)
# Fit model
h.c.g.raw.aegyptium.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                                model.file = "host.climate.genus_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                                n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
h.c.g.raw.aegyptium.fit.mcmc <- as.mcmc(h.c.g.raw.aegyptium.fit) 
closeAllConnections()

# Get posterior means for parameters for each host group
summary_fit.h.c.g.raw <- summary(h.c.g.raw.aegyptium.fit.mcmc)$statistics
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
fitted_data.h.c.g.aegyptium.raw <- do.call(rbind, fitted_list)

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
fitted_data.h.c.g.aegyptium.raw <- fitted_data.h.c.g.aegyptium.raw |>  
  left_join(species_to_h.c.g.raw, by = c("Host", "Climate", "Genus"))|> 
  filter(!is.na(Species)) 

# Order factor levels
plot_data.h.c.g.raw$Host <- factor(plot_data.h.c.g.raw$Host, levels = c("One", "Two", "Three"))
fitted_data.h.c.g.aegyptium.raw$Host <- factor(fitted_data.h.c.g.aegyptium.raw$Host, levels = c("One", "Two", "Three"))

# Combine Host, Climate, and Size for coloring
fitted_data.h.c.g.aegyptium.raw$`Host, Climate, Genus` <- with(fitted_data.h.c.g.aegyptium.raw, paste(Host, Climate, Genus, sep = ", "))

# Plot observed data and fitted curves with facet wrap by species
p.h.c.g.aegyptium.raw <- ggplot() +
  geom_point(data = plot_data.h.c.g.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.h.c.g.aegyptium.raw, aes(x = Temperature, y = Fitted, color = `Host, Climate, Genus`), size = 1) +
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
p.h.c.g.aegyptium.raw

# Merge the observed and predicted data for comparison
merged_data.h.c.g.aegyptium.raw <- fitted_data.h.c.g.aegyptium.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species", 
                             "Host" = "Host", "Climate" = "Climate",
                             "Genus" = "Genus"))

# Calculate RMSE for each species
mse_by_h.c.g.aegyptium.raw <- merged_data.h.c.g.aegyptium.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.h.c.g.aegyptium.raw <- merged_data.h.c.g.aegyptium.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 


save(mse_by_h.c.g.aegyptium.raw,
     fitted_data.h.c.g.aegyptium.raw,
     merged_data.h.c.g.aegyptium.raw,
     h.c.g.raw.aegyptium.fit.mcmc,
     h.c.g.raw.aegyptium.fit,
     p.h.c.g.aegyptium.raw,
     file="host.climate.genus.raw.aegyptium.fit.RData")

################################### By Host, Genus, and Climate- Dromedarii ###################################
# Model by typical number of host species
tick.rm.raw <- tick.raw[tick.raw$Species != "Hyalomma dromedarii", ]
climate.rm.raw <- as.factor(tick.rm.raw$Climate)
n.climate.rm.raw <- length(unique(tick.rm.raw$Climate))
genus.rm.raw <- as.factor(tick.rm.raw$Genus)
n.genus.rm.raw <- length(unique(tick.rm.raw$Genus))
host.rm.raw <- as.factor(tick.rm.raw$Host)
n.host.rm.raw <- length(unique(tick.rm.raw$Host))
N.obs.rm.raw <- length(tick.rm.raw$Trait)
unique(tick.rm.raw$Species)
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

  for (j in 1:n.host.rm.raw) {
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
  }
  for (k in 1:n.climate.rm.raw) {
    mu.climate.Topt[k] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific vertex x
    climate.tau[k] ~ dexp(tau1)
  }
  for (s in 1:n.genus.rm.raw) {
    mu.genus.c[s] ~ dexp(mu.c)       # genus-specific vertex y
  }
  # Likelihood
  for (i in 1:N.obs.rm.raw) {
    mu[i] <- mu.host.a[host.rm.raw[i]] * (mu.climate.Topt[climate.rm.raw[i]] - temp[i])^2 + mu.genus.c[genus.rm.raw[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[climate.rm.raw[i]]) T(0, ) 
  }}
", file = "host.climate.genus_robin.txt")


# Parameters to monitor
parameters <- c("mu.host.a", "mu.genus.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt", "climate.tau",
                "tau1", "mu.climate.Topt")


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
jag.data <- list(trait = tick.rm.raw$Trait, N.obs.rm.raw = N.obs.rm.raw, 
                 temp = tick.rm.raw$Temp, 
                 n.host.rm.raw = n.host.rm.raw, host.rm.raw = host.rm.raw, 
                 n.climate.rm.raw = n.climate.rm.raw, climate.rm.raw = climate.rm.raw,
                 n.genus.rm.raw = n.genus.rm.raw, genus.rm.raw = genus.rm.raw)
# Fit model
h.c.g.raw.dromedarii.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                                model.file = "host.climate.genus_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                                n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
h.c.g.raw.dromedarii.fit.mcmc <- as.mcmc(h.c.g.raw.dromedarii.fit) 
closeAllConnections()

# Get posterior means for parameters for each host group
summary_fit.h.c.g.raw <- summary(h.c.g.raw.dromedarii.fit.mcmc)$statistics
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
fitted_data.h.c.g.dromedarii.raw <- do.call(rbind, fitted_list)

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
fitted_data.h.c.g.dromedarii.raw <- fitted_data.h.c.g.dromedarii.raw |>  
  left_join(species_to_h.c.g.raw, by = c("Host", "Climate", "Genus"))|> 
  filter(!is.na(Species)) 

# Order factor levels
plot_data.h.c.g.raw$Host <- factor(plot_data.h.c.g.raw$Host, levels = c("One", "Two", "Three"))
fitted_data.h.c.g.dromedarii.raw$Host <- factor(fitted_data.h.c.g.dromedarii.raw$Host, levels = c("One", "Two", "Three"))

# Combine Host, Climate, and Size for coloring
fitted_data.h.c.g.dromedarii.raw$`Host, Climate, Genus` <- with(fitted_data.h.c.g.dromedarii.raw, paste(Host, Climate, Genus, sep = ", "))

# Plot observed data and fitted curves with facet wrap by species
p.h.c.g.dromedarii.raw <- ggplot() +
  geom_point(data = plot_data.h.c.g.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.h.c.g.dromedarii.raw, aes(x = Temperature, y = Fitted, color = `Host, Climate, Genus`), size = 1) +
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
p.h.c.g.dromedarii.raw

# Merge the observed and predicted data for comparison
merged_data.h.c.g.dromedarii.raw <- fitted_data.h.c.g.dromedarii.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species", 
                             "Host" = "Host", "Climate" = "Climate",
                             "Genus" = "Genus"))

# Calculate RMSE for each species
mse_by_h.c.g.dromedarii.raw <- merged_data.h.c.g.dromedarii.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.h.c.g.dromedarii.raw <- merged_data.h.c.g.dromedarii.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 


save(mse_by_h.c.g.dromedarii.raw,
     fitted_data.h.c.g.dromedarii.raw,
     merged_data.h.c.g.dromedarii.raw,
     h.c.g.raw.dromedarii.fit.mcmc,
     h.c.g.raw.dromedarii.fit,
     p.h.c.g.dromedarii.raw,
     file="host.climate.genus.raw.dromedarii.fit.RData")

################################### By Host, Genus, and Climate- Lusitanicum ###################################
# Model by typical number of host species
tick.rm.raw <- tick.raw[tick.raw$Species != "Hyalomma lusitanicum", ]
climate.rm.raw <- as.factor(tick.rm.raw$Climate)
n.climate.rm.raw <- length(unique(tick.rm.raw$Climate))
genus.rm.raw <- as.factor(tick.rm.raw$Genus)
n.genus.rm.raw <- length(unique(tick.rm.raw$Genus))
host.rm.raw <- as.factor(tick.rm.raw$Host)
n.host.rm.raw <- length(unique(tick.rm.raw$Host))
N.obs.rm.raw <- length(tick.rm.raw$Trait)
unique(tick.rm.raw$Species)
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

  for (j in 1:n.host.rm.raw) {
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
  }
  for (k in 1:n.climate.rm.raw) {
    mu.climate.Topt[k] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific vertex x
    climate.tau[k] ~ dexp(tau1)
  }
  for (s in 1:n.genus.rm.raw) {
    mu.genus.c[s] ~ dexp(mu.c)       # genus-specific vertex y
  }
  # Likelihood
  for (i in 1:N.obs.rm.raw) {
    mu[i] <- mu.host.a[host.rm.raw[i]] * (mu.climate.Topt[climate.rm.raw[i]] - temp[i])^2 + mu.genus.c[genus.rm.raw[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[climate.rm.raw[i]]) T(0, ) 
  }}
", file = "host.climate.genus_robin.txt")


# Parameters to monitor
parameters <- c("mu.host.a", "mu.genus.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt", "climate.tau",
                "tau1", "mu.climate.Topt")


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
jag.data <- list(trait = tick.rm.raw$Trait, N.obs.rm.raw = N.obs.rm.raw, 
                 temp = tick.rm.raw$Temp, 
                 n.host.rm.raw = n.host.rm.raw, host.rm.raw = host.rm.raw, 
                 n.climate.rm.raw = n.climate.rm.raw, climate.rm.raw = climate.rm.raw,
                 n.genus.rm.raw = n.genus.rm.raw, genus.rm.raw = genus.rm.raw)
# Fit model
h.c.g.raw.lusitanicum.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                                model.file = "host.climate.genus_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                                n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
h.c.g.raw.lusitanicum.fit.mcmc <- as.mcmc(h.c.g.raw.lusitanicum.fit) 
closeAllConnections()

# Get posterior means for parameters for each host group
summary_fit.h.c.g.raw <- summary(h.c.g.raw.lusitanicum.fit.mcmc)$statistics
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
fitted_data.h.c.g.lusitanicum.raw <- do.call(rbind, fitted_list)

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
fitted_data.h.c.g.lusitanicum.raw <- fitted_data.h.c.g.lusitanicum.raw |>  
  left_join(species_to_h.c.g.raw, by = c("Host", "Climate", "Genus"))|> 
  filter(!is.na(Species)) 

# Order factor levels
plot_data.h.c.g.raw$Host <- factor(plot_data.h.c.g.raw$Host, levels = c("One", "Two", "Three"))
fitted_data.h.c.g.lusitanicum.raw$Host <- factor(fitted_data.h.c.g.lusitanicum.raw$Host, levels = c("One", "Two", "Three"))

# Combine Host, Climate, and Size for coloring
fitted_data.h.c.g.lusitanicum.raw$`Host, Climate, Genus` <- with(fitted_data.h.c.g.lusitanicum.raw, paste(Host, Climate, Genus, sep = ", "))

# Plot observed data and fitted curves with facet wrap by species
p.h.c.g.lusitanicum.raw <- ggplot() +
  geom_point(data = plot_data.h.c.g.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.h.c.g.lusitanicum.raw, aes(x = Temperature, y = Fitted, color = `Host, Climate, Genus`), size = 1) +
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
p.h.c.g.lusitanicum.raw

# Merge the observed and predicted data for comparison
merged_data.h.c.g.lusitanicum.raw <- fitted_data.h.c.g.lusitanicum.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species", 
                             "Host" = "Host", "Climate" = "Climate",
                             "Genus" = "Genus"))

# Calculate RMSE for each species
mse_by_h.c.g.lusitanicum.raw <- merged_data.h.c.g.lusitanicum.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.h.c.g.lusitanicum.raw <- merged_data.h.c.g.lusitanicum.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 


save(mse_by_h.c.g.lusitanicum.raw,
     fitted_data.h.c.g.lusitanicum.raw,
     merged_data.h.c.g.lusitanicum.raw,
     h.c.g.raw.lusitanicum.fit.mcmc,
     h.c.g.raw.lusitanicum.fit,
     p.h.c.g.lusitanicum.raw,
     file="host.climate.genus.raw.lusitanicum.fit.RData")

################################### By Host, Genus, and Climate- Schulzei ###################################
# Model by typical number of host species
tick.rm.raw <- tick.raw[tick.raw$Species != "Hyalomma schulzei", ]
climate.rm.raw <- as.factor(tick.rm.raw$Climate)
n.climate.rm.raw <- length(unique(tick.rm.raw$Climate))
genus.rm.raw <- as.factor(tick.rm.raw$Genus)
n.genus.rm.raw <- length(unique(tick.rm.raw$Genus))
host.rm.raw <- as.factor(tick.rm.raw$Host)
n.host.rm.raw <- length(unique(tick.rm.raw$Host))
N.obs.rm.raw <- length(tick.rm.raw$Trait)
unique(tick.rm.raw$Species)
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

  for (j in 1:n.host.rm.raw) {
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
  }
  for (k in 1:n.climate.rm.raw) {
    mu.climate.Topt[k] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific vertex x
    climate.tau[k] ~ dexp(tau1)
  }
  for (s in 1:n.genus.rm.raw) {
    mu.genus.c[s] ~ dexp(mu.c)       # genus-specific vertex y
  }
  # Likelihood
  for (i in 1:N.obs.rm.raw) {
    mu[i] <- mu.host.a[host.rm.raw[i]] * (mu.climate.Topt[climate.rm.raw[i]] - temp[i])^2 + mu.genus.c[genus.rm.raw[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[climate.rm.raw[i]]) T(0, ) 
  }}
", file = "host.climate.genus_robin.txt")


# Parameters to monitor
parameters <- c("mu.host.a", "mu.genus.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt", "climate.tau",
                "tau1", "mu.climate.Topt")


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
jag.data <- list(trait = tick.rm.raw$Trait, N.obs.rm.raw = N.obs.rm.raw, 
                 temp = tick.rm.raw$Temp, 
                 n.host.rm.raw = n.host.rm.raw, host.rm.raw = host.rm.raw, 
                 n.climate.rm.raw = n.climate.rm.raw, climate.rm.raw = climate.rm.raw,
                 n.genus.rm.raw = n.genus.rm.raw, genus.rm.raw = genus.rm.raw)
# Fit model
h.c.g.raw.schulzei.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                                model.file = "host.climate.genus_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                                n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
h.c.g.raw.schulzei.fit.mcmc <- as.mcmc(h.c.g.raw.schulzei.fit) 
closeAllConnections()

# Get posterior means for parameters for each host group
summary_fit.h.c.g.raw <- summary(h.c.g.raw.schulzei.fit.mcmc)$statistics
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
fitted_data.h.c.g.schulzei.raw <- do.call(rbind, fitted_list)

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
fitted_data.h.c.g.schulzei.raw <- fitted_data.h.c.g.schulzei.raw |>  
  left_join(species_to_h.c.g.raw, by = c("Host", "Climate", "Genus"))|> 
  filter(!is.na(Species)) 

# Order factor levels
plot_data.h.c.g.raw$Host <- factor(plot_data.h.c.g.raw$Host, levels = c("One", "Two", "Three"))
fitted_data.h.c.g.schulzei.raw$Host <- factor(fitted_data.h.c.g.schulzei.raw$Host, levels = c("One", "Two", "Three"))

# Combine Host, Climate, and Size for coloring
fitted_data.h.c.g.schulzei.raw$`Host, Climate, Genus` <- with(fitted_data.h.c.g.schulzei.raw, paste(Host, Climate, Genus, sep = ", "))

# Plot observed data and fitted curves with facet wrap by species
p.h.c.g.schulzei.raw <- ggplot() +
  geom_point(data = plot_data.h.c.g.raw, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.h.c.g.schulzei.raw, aes(x = Temperature, y = Fitted, color = `Host, Climate, Genus`), size = 1) +
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
p.h.c.g.schulzei.raw

# Merge the observed and predicted data for comparison
merged_data.h.c.g.schulzei.raw <- fitted_data.h.c.g.schulzei.raw %>%
  left_join(tick.raw, by = c("Temperature" = "Temp", "Species" = "Species", 
                             "Host" = "Host", "Climate" = "Climate",
                             "Genus" = "Genus"))

# Calculate RMSE for each species
mse_by_h.c.g.schulzei.raw <- merged_data.h.c.g.schulzei.raw %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.h.c.g.schulzei.raw <- merged_data.h.c.g.schulzei.raw %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 


save(mse_by_h.c.g.schulzei.raw,
     fitted_data.h.c.g.schulzei.raw,
     merged_data.h.c.g.schulzei.raw,
     h.c.g.raw.schulzei.fit.mcmc,
     h.c.g.raw.schulzei.fit,
     p.h.c.g.schulzei.raw,
     file="host.climate.genus.raw.schulzei.fit.RData")


################################### RMSE ###################################
# Data frame
table.tick.rm.raw <- data.frame(
  Method = factor(c(rep("By Genus (RM)", 9 * 7), 
                    rep("By Host (RM)", 9 * 8), rep("By Host and Climate (RM)", 9 * 6), 
                    rep("By Host, Climate, and Genus (RM)", 9 * 4))),
  Species = factor(rep(c("Amblyomma lepidum", "Dermacentor andersoni",
                         "Dermacentor nitens", "Haemaphysalis leporispalustris",
                         "Hyalomma aegyptium", "Hyalomma dromedarii",
                         "Hyalomma impeltatum", "Hyalomma lusitanicum",
                         "Hyalomma schulzei"), times = 25)),
  # MSE = c(mse_by_genus.andersoni$MSE, mse_by_genus.nitens$MSE, 
  #         mse_by_genus.aegyptium$MSE, mse_by_genus.dromedarii$MSE, mse_by_genus.impeltatum$MSE,
  #         mse_by_genus.lusitanicum$MSE, mse_by_genus.schulzei$MSE, 
  #         mse_by_host.lepidum$MSE, mse_by_host.andersoni$MSE, mse_by_host.leporispalustris$MSE,
  #         mse_by_host.aegyptium$MSE, mse_by_host.dromedarii$MSE, mse_by_host.impeltatum$MSE,
  #         mse_by_host.lusitanicum$MSE, mse_by_host.schulzei$MSE,
  #         mse_by_host.climate.leporispalustris$MSE,
  #         mse_by_host.climate.aegyptium$MSE, mse_by_host.climate.dromedarii$MSE, 
  #         mse_by_host.climate.impeltatum$MSE,
  #         mse_by_host.climate.lusitanicum$MSE, mse_by_host.climate.schulzei$MSE,
  #         mse_by_h.c.g.aegyptium$MSE, mse_by_h.c.g$MSE,
  #         mse_by_h.c.g.lusitanicum$MSE, mse_by_h.c.g.schulzei$MSE),
  Total.RMSE = rep(as.numeric(c(overall_rmse.genus.andersoni.raw, overall_rmse.genus.nitens.raw,
                                overall_rmse.genus.aegyptium.raw, overall_rmse.genus.dromedarii.raw,
                                overall_rmse.genus.impeltatum.raw, overall_rmse.genus.lusitanicum.raw,
                                overall_rmse.genus.schulzei.raw,
                                overall_rmse.host.lepidum.raw, overall_rmse.host.andersoni.raw,
                                overall_rmse.host.leporispalustris.raw, overall_rmse.host.aegyptium.raw,
                                overall_rmse.host.dromedarii.raw, overall_rmse.host.impeltatum.raw,
                                overall_rmse.host.lusitanicum.raw, overall_rmse.host.schulzei.raw,
                                overall_rmse.host.climate.leporispalustris.raw,
                                overall_rmse.host.climate.aegyptium.raw, overall_rmse.host.climate.dromedarii.raw,
                                overall_rmse.host.climate.impeltatum.raw, overall_rmse.host.climate.lusitanicum.raw,
                                overall_rmse.host.climate.schulzei.raw,
                                overall_rmse.h.c.g.aegyptium.raw, overall_rmse.h.c.g.dromedarii.raw,
                                overall_rmse.h.c.g.lusitanicum.raw, overall_rmse.h.c.g.schulzei.raw)), each = 9),
  Removed = c(rep(c("Dermacentor andersoni",
                    "Dermacentor nitens",
                    "Hyalomma aegyptium", "Hyalomma dromedarii",
                    "Hyalomma impeltatum", "Hyalomma lusitanicum",
                    "Hyalomma schulzei"), each = 9),
              rep(c("Amblyomma lepidum", "Dermacentor andersoni",
                    "Haemaphysalis leporispalustris",
                    "Hyalomma aegyptium", "Hyalomma dromedarii",
                    "Hyalomma impeltatum", "Hyalomma lusitanicum",
                    "Hyalomma schulzei"), each = 9),
              rep(c("Haemaphysalis leporispalustris",
                    "Hyalomma aegyptium", "Hyalomma dromedarii",
                    "Hyalomma impeltatum", "Hyalomma lusitanicum",
                    "Hyalomma schulzei"), each = 9),
              rep(c("Hyalomma aegyptium",
                    "Hyalomma dromedarii",
                    "Hyalomma lusitanicum",
                    "Hyalomma schulzei"), each = 9))
)

tick.DIC.rm.raw <- data.frame(
  Method = rep(c("By Genus", "By Host",
                 "By Host and Climate",
                 "By Host, Genus, and Climate"), each = 4),
  Removed = rep(c("Hyalomma aegyptium", "Hyalomma dromedarii",
                  "Hyalomma lusitanicum", "Hyalomma schulzei"), times = 4),
  DIC = c(genus.raw.aegyptium.fit$BUGSoutput$DIC, genus.raw.dromedarii.fit$BUGSoutput$DIC,
          genus.raw.lusitanicum.fit$BUGSoutput$DIC, genus.raw.schulzei.fit$BUGSoutput$DIC,
          host.raw.aegyptium.fit$BUGSoutput$DIC, host.raw.dromedarii.fit$BUGSoutput$DIC,
          host.raw.lusitanicum.fit$BUGSoutput$DIC, host.raw.schulzei.fit$BUGSoutput$DIC,
          host.climate.raw.aegyptium.fit$BUGSoutput$DIC, host.climate.raw.dromedarii.fit$BUGSoutput$DIC,
          host.climate.raw.lusitanicum.fit$BUGSoutput$DIC, host.climate.raw.schulzei.fit$BUGSoutput$DIC,
          h.c.g.raw.aegyptium.fit$BUGSoutput$DIC, h.c.g.raw.dromedarii.fit$BUGSoutput$DIC,
          h.c.g.raw.lusitanicum.fit$BUGSoutput$DIC, h.c.g.raw.schulzei.fit$BUGSoutput$DIC)
)

table.tick.rm.raw |> 
  ggplot(aes(x = Removed, y = Total.RMSE, fill = Method))+
  geom_bar(stat = "identity", position = "dodge",
           color = "black")+
  theme_minimal() +
  scale_fill_viridis_d(option = "mako") +
  coord_cartesian(ylim = c(6, 11)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

tick.DIC.rm.raw |> 
  ggplot(aes(x = Removed, y = DIC, fill = Method))+
  geom_bar(stat = "identity", position = "dodge", color = "black")+
  theme_minimal() +
  scale_fill_viridis_d(option = "mako") +
  coord_cartesian(ylim = c(7500, 10600)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save(table.tick.rm.raw,
     tick.DIC.rm.raw, file = "Dataset.Tick.Raw.RM.RData")








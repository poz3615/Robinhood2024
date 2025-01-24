library(dplyr)
library(tidyr)
library(ggplot2)
library(rjags)
library(R2jags)
library(readxl)
library(truncnorm)
library(scoringRules)



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
  c ~ dnorm(10, 0.001)      # Oviposition rate at Topt
  
  # Likelihood
  for (i in 1:N.obs) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], tau1)  
  }}
", file = "full_data_robin.txt")

# Initial values
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
                                   Species = factor(rep(1:n.species, each = length(temp_seq))))
# Original data with temperature, trait, and species as numeric
plot_data.full <- data.frame(
  Temperature = tick.abr.new$Temp,
  Trait = tick.abr.new$Trait,
  Species = as.numeric(tick.abr.new$Species)
)

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.full, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.full.data, aes(x = Temperature, y = Fitted, color = Species), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Species",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 

# Ensure original data and fitted data are aligned correctly
# Make sure the species column in both datasets uses the same format (numeric)
fitted_data.full.data$Species <- as.numeric(fitted_data.full.data$Species)
tick$Species.number <- as.numeric(tick$Species)
# Merge the observed and predicted data for comparison
merged_data.full <- fitted_data.full.data %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species.number"))

# Calculate RMSE for each species
rmse_by_full_data <- merged_data.full %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_full_data)
# Average RMSE across species
mean(rmse_by_full_data$RMSE)


################################### By Species ###################################
# Model using species effect
sink("species_robin.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Species-specific 
  mu.c ~ dnorm(7, 0.1)          # Mean for species effects c (vertex y)
  tau.c ~ dexp(0.01)             # Tau for species effect c (vertex y)
  
  mu.Topt ~ dnorm(35, 0.1)          # Mean for species effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for species effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for species effects a (width)

  for (j in 1:n.species) {
    mu.species.c[j] ~ dnorm(mu.c, tau.c)       # Species-specific vertex y
    mu.species.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Species-specific vertex x
    mu.species.a[j] ~ dexp(mu.a)    # Species-specific width
  }
  # Likelihood
  for (i in 1:N.obs) {
    mu[i] <- mu.species.a[species[i]] * (mu.species.Topt[species[i]] - temp[i])^2 + mu.species.c[species[i]]
    trait[i] ~ dnorm(mu[i], tau1)  
  }}
", file = "species_robin.txt")

# Initial values (arbitrary)
inits <- function() {
  list(
    mu.a = 1,             
    mu.c = 0,          
    tau.c = 1,
    tau1 = 1,
    mu.Topt = 20,
    tau.Topt = 1
  )
}

# Parameters to monitor
parameters <- c("mu.species.a", "mu.species.c", "tau.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt",
                "tau1", "mu.species.Topt")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

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
# DIC
species.fit$BUGSoutput$DIC

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
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~ Species) 

# Ensure original data and fitted data are aligned correctly
# Make sure the species column in both datasets uses the same format (numeric)
#fitted_data.species$Species <- as.numeric(fitted_data.species$Species)
#tick$Species.number <- as.numeric(tick$Species)
# Merge the observed and predicted data for comparison
merged_data.species <- fitted_data.species %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species"))

# Calculate RMSE for each species
rmse_by_species <- merged_data.species %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_species)
# Average RMSE over species
mean(rmse_by_species$RMSE)




################################### By Genus ###################################
# Model using genus effect
sink("genus_robin.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Genus-specific 
  mu.c ~ dnorm(7, 0.1)          # Mean for genus effects c (vertex y)
  tau.c ~ dexp(0.01)             # Tau for genus effect c (vertex y)
  
  mu.Topt ~ dnorm(35, 0.1)          # Mean for genus effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for genus effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for genus effects a (width)

  for (j in 1:n.genus) {
    mu.genus.c[j] ~ dnorm(mu.c, tau.c)       # Genus-specific vertex y
    mu.genus.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Genus-specific vertex x
    mu.genus.a[j] ~ dexp(mu.a)    # Genus-specific width
  }
  # Likelihood
  for (i in 1:N.obs) {
    mu[i] <- mu.genus.a[genus[i]] * (mu.genus.Topt[genus[i]] - temp[i])^2 + mu.genus.c[genus[i]]
    trait[i] ~ dnorm(mu[i], tau1)  
  }}
", file = "genus_robin.txt")

# Initial values (arbitrary)
inits <- function() {
  list(
    mu.a = 1,             
    mu.c = 0,          
    tau.c = 1,
    tau1 = 1,
    mu.Topt = 20,
    tau.Topt = 1
  )
}

# Parameters to monitor
parameters <- c("mu.genus.a", "mu.genus.c", "tau.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt",
                "tau1", "mu.genus.Topt")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

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
# DIC
genus.fit$BUGSoutput$DIC

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

# Ensure original data and fitted data are aligned correctly
# Make sure the species column in both datasets uses the same format (numeric)
#fitted_data.genus$Genus <- as.numeric(fitted_data.genus$Genus)
#tick$Genus.number <- as.numeric(tick$Genus)
#tick$Species.number <- as.numeric(tick$Species)
# Merge the observed and fitted data for comparison
merged_data.genus <- fitted_data.genus %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Genus" = "Genus"))

# Calculate RMSE for each species
rmse_by_genus <- merged_data.genus %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_genus)
# Avergae RMSE across species
mean(rmse_by_genus$RMSE)


################################### By Climate ###################################
# Model using climate effect
sink("climate_robin.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Climate-specific 
  mu.c ~ dnorm(7, 0.1)          # Mean for climate effect c (vertex y)
  tau.c ~ dexp(0.01)             # Tau for climate effect c (vertex y)
  
  mu.Topt ~ dnorm(35, 0.1)          # Mean for climate effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for climate effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for climate effects a (width)

  for (j in 1:n.climate) {
    mu.climate.c[j] ~ dnorm(mu.c, tau.c)       # Climate-specific vertex y
    mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific vertex x
    mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific width
  }
  # Likelihood
  for (i in 1:N.obs) {
    mu[i] <- mu.climate.a[climate[i]] * (mu.climate.Topt[climate[i]] - temp[i])^2 + mu.climate.c[climate[i]]
    trait[i] ~ dnorm(mu[i], tau1)  
  }}
", file = "climate_robin.txt")

# Initial values (arbitrary)
inits <- function() {
  list(
    mu.a = 1,             
    mu.c = 0,          
    tau.c = 1,
    tau1 = 1,
    mu.Topt = 20,
    tau.Topt = 1
  )
}

# Parameters to monitor
parameters <- c("mu.climate.a", "mu.climate.c", "tau.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt",
                "tau1", "mu.climate.Topt")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

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
# DIC
climate.fit$BUGSoutput$DIC

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

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.climate, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.climate, aes(x = Temperature, y = Fitted, color = Climate), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Climate",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 

# Ensure test and fitted_data are aligned correctly
# Make sure the species column in both datasets uses the same format (numeric)
#fitted_data.climate$Climate <- as.numeric(fitted_data.climate$Climate)
#tick$Climate.number <- as.numeric(tick$Climate)
# Merge the observed and predicted data for comparison
merged_data.climate <- fitted_data.climate %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Climate" = "Climate"))

# Calculate RMSE for each species
rmse_by_climate <- merged_data.climate %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_climate)
# Average RMSE across species
mean(rmse_by_climate$RMSE)



################################### By Host ###################################
# Model by typical number of host species
sink("host_robin.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Host-specific 
  mu.c ~ dnorm(7, 0.1)          # Mean for host effects c (vertex y)
  tau.c ~ dexp(0.01)             # Tau for host effect c (vertex y)
  
  mu.Topt ~ dnorm(35, 0.1)          # Mean for host effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for host effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for host effects a (width)

  for (j in 1:n.host) {
    mu.host.c[j] ~ dnorm(mu.c, tau.c)       # Host-specific vertex y
    mu.host.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Host-specific vertex x
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
  }
  # Likelihood
  for (i in 1:N.obs) {
    mu[i] <- mu.host.a[host[i]] * (mu.host.Topt[host[i]] - temp[i])^2 + mu.host.c[host[i]]
    trait[i] ~ dnorm(mu[i], tau1)  
  }}
", file = "host_robin.txt")

# Initial values (arbitrary)
inits <- function() {
  list(
    mu.a = 1,             
    mu.c = 0,          
    tau.c = 1,
    tau1 = 1,
    mu.Topt = 20,
    tau.Topt = 1
  )
}

# Parameters to monitor
parameters <- c("mu.host.a", "mu.host.c", "tau.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt",
                "tau1", "mu.host.Topt")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

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
# DIC
host.fit$BUGSoutput$DIC

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

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.host, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host, aes(x = Temperature, y = Fitted, color = Host), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Host",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 

# Ensure test and fitted_data are aligned correctly
# Make sure the species column in both datasets uses the same format (numeric)
#fitted_data.host$Host <- as.numeric(fitted_data.host$Host)
#tick$Host.number <- as.numeric(tick$Host)
#tick$Species.number <- as.numeric(tick$Species)
# Merge the observed and predicted data for comparison
merged_data.host <- fitted_data.host %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host"))

# Calculate RMSE for each species
rmse_by_host <- merged_data.host %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_host)
# Average RMSE across species
mean(rmse_by_host$RMSE)

################################### By Size ###################################
# Model by size group
sink("size_robin.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Host-specific 
  mu.c ~ dnorm(7, 0.1)          # Mean for size effects c (vertex y)
  tau.c ~ dexp(0.01)             # Tau for size effect c (vertex y)
  
  mu.Topt ~ dnorm(35, 0.1)          # Mean for size effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for size effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for size effects a (width)

  for (j in 1:n.size) {
    mu.size.c[j] ~ dnorm(mu.c, tau.c)       # Size-specific vertex y
    mu.size.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Size-specific vertex x
    mu.size.a[j] ~ dexp(mu.a)    # Size-specific width
  }
  # Likelihood
  for (i in 1:N.obs) {
    mu[i] <- mu.size.a[size[i]] * (mu.size.Topt[size[i]] - temp[i])^2 + mu.size.c[size[i]]
    trait[i] ~ dnorm(mu[i], tau1)  
  }}
", file = "size_robin.txt")

# Initial values (arbitrary)
inits <- function() {
  list(
    mu.a = 1,             
    mu.c = 0,          
    tau.c = 1,
    tau1 = 1,
    mu.Topt = 20,
    tau.Topt = 1
  )
}

# Parameters to monitor
parameters <- c("mu.size.a", "mu.size.c", "tau.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt",
                "tau1", "mu.size.Topt")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

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
# DIC
size.fit$BUGSoutput$DIC

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

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.size, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.size, aes(x = Temperature, y = Fitted, color = Size), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Size",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 

# Ensure test and fitted data are aligned correctly
# Make sure the species column in both datasets uses the same format (numeric)
#fitted_data.size$Size <- as.numeric(fitted_data.size$Size)
#tick$Size.number <- as.numeric(tick$Size)
#tick$Species.number <- as.numeric(tick$Species)
# Merge the observed and predicted data for comparison
merged_data.size <- fitted_data.size %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Size" = "Size"))

# Calculate RMSE for each species
rmse_by_size <- merged_data.size %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_size)
# Average RMSE across species
mean(rmse_by_size$RMSE)

################################### By Continent ###################################
# Model by number of continents tick is found on
sink("continent_robin.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Continent-specific 
  mu.c ~ dnorm(7, 0.1)          # Mean for continent effects c (vertex y)
  tau.c ~ dexp(0.01)             # Tau for continent effect c (vertex y)
  
  mu.Topt ~ dnorm(35, 0.1)          # Mean for continent effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for continent effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for continent effects a (width)

  for (j in 1:n.continent) {
    mu.continent.c[j] ~ dnorm(mu.c, tau.c)       # Continent-specific vertex y
    mu.continent.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Continent-specific vertex x
    mu.continent.a[j] ~ dexp(mu.a)    # Continent-specific width
  }
  # Likelihood
  for (i in 1:N.obs) {
    mu[i] <- mu.continent.a[continent[i]] * (mu.continent.Topt[continent[i]] - temp[i])^2 + mu.continent.c[continent[i]]
    trait[i] ~ dnorm(mu[i], tau1)  
  }}
", file = "continent_robin.txt")

# Initial values (arbitrary)
inits <- function() {
  list(
    mu.a = 1,             
    mu.c = 0,          
    tau.c = 1,
    tau1 = 1,
    mu.Topt = 20,
    tau.Topt = 1
  )
}

# Parameters to monitor
parameters <- c("mu.continent.a", "mu.continent.c", "tau.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt",
                "tau1", "mu.continent.Topt")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

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
# DIC
continent.fit$BUGSoutput$DIC

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

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.continent, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.continent, aes(x = Temperature, y = Fitted, color = Continent), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Continent",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 

# Ensure test and fitted data are aligned correctly
# Make sure the species column in both datasets uses the same format (numeric)
#fitted_data.continent$Continent <- as.numeric(fitted_data.continent$Continent)
#tick$Continent.number <- as.numeric(tick$Continent)
#tick$Species.number <- as.numeric(tick$Species)
# Merge the observed and predicted data for comparison
merged_data.continent <- fitted_data.continent %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Continent" = "Continent"))

# Calculate RMSE for each species
rmse_by_continent <- merged_data.continent %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_continent)
# Average RMSE across species
mean(rmse_by_continent$RMSE)
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
table.tick <- matrix(0, nrow = 7, ncol = 2)
row.names(table.tick) <- c("Full Data", "By Species", "By Genus", 
                           "By Climate", "By Host", "By Size", "By Continents")
colnames(table.tick) <- c("DIC", "RMSE")
table.tick[1, 1] <- full.data.fit$BUGSoutput$DIC
table.tick[2, 1] <- species.fit$BUGSoutput$DIC
table.tick[3, 1] <- genus.fit$BUGSoutput$DIC
table.tick[4, 1] <- climate.fit$BUGSoutput$DIC
table.tick[5, 1] <- host.fit$BUGSoutput$DIC
table.tick[6, 1] <- size.fit$BUGSoutput$DIC
table.tick[7, 1] <- continent.fit$BUGSoutput$DIC
table.tick[1, 2] <- mean(rmse_by_full_data$RMSE)
table.tick[2, 2] <- mean(rmse_by_species$RMSE)
table.tick[3, 2] <- mean(rmse_by_genus$RMSE)
table.tick[4, 2] <- mean(rmse_by_climate$RMSE)
table.tick[5, 2] <- mean(rmse_by_host$RMSE)
table.tick[6, 2] <- mean(rmse_by_size$RMSE)
table.tick[7, 2] <- mean(rmse_by_continent$RMSE)

table.tick

################################### Take temperature out ###################################
# Take out 20, 25, 30
tick.notemp <- tick.abr.new[!(tick.abr.new$Temp %in% c(20, 25, 30)),]
tick.notemp$Species <- droplevels(tick.notemp$Species)
# New number of observations
N.obs.notemp <- length(tick.notemp$Trait)
# Species number for hierarchical model
species.notemp <- as.numeric(tick.notemp$Species)
# Number of species
n.species.notemp <- length(unique(tick.notemp$Species))

# Genus number for hierarchical model
genus.notemp <- as.numeric(tick.notemp$Genus)
# Number of genus
n.genus.notemp <- length(unique(tick.notemp$Genus))

# Host number for hierarchical model
host.notemp <- as.numeric(tick.notemp$Host)
# Number of hosts
n.host.notemp <- length(unique(tick.notemp$Host))

# Climate number for hierarchical model
climate.notemp <- as.numeric(tick.notemp$Climate)
# Number of climates
n.climate.notemp <- length(unique(tick.notemp$Climate))

# Continent number for hierarchical model
continent.notemp <- as.numeric(tick.notemp$Continent)
# Number of continents
n.continent.notemp <- length(unique(tick.notemp$Continent))

# Size number for hierarchical model
size.notemp <- as.numeric(tick.notemp$Size)
# Number of sizes
n.size.notemp <- length(unique(tick.notemp$Size))
################################### Full Data- Take temperature out ###################################
# Model using full data
sink("full_data_robin.txt")
cat("
model {
  # Priors
   a ~ dexp(0.1)             # Width
  Topt ~ dnorm(30, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  c ~ dnorm(10, 0.001)      # Oviposition rate at Topt
  
  # Likelihood
  for (i in 1:N.obs.notemp) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], tau1)  
  }}
", file = "full_data_robin.txt")

# Initial values
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
jag.data <- list(trait = tick.notemp$Trait, N.obs.notemp = N.obs.notemp, 
                 temp = tick.notemp$Temp)

# Fit fill model
full.data.notemp.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                      model.file = "full_data_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                      n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
full.data.notemp.fit.mcmc <- as.mcmc(full.data.notemp.fit) 
closeAllConnections()
#DIC
full.data.notemp.fit$BUGSoutput$DIC

# Get posterior means
summary_fit.notemp <- summary(full.data.notemp.fit.mcmc)$statistics
a_mean.notemp <- summary_fit.notemp["a", "Mean"]
c_mean.notemp <- summary_fit.notemp["c", "Mean"]
Topt_mean.notemp <- summary_fit.notemp["Topt", "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp), max(tick.abr.new$Temp), by = 0.1)

# Get predicted values
full.data <- function (x){
  a_mean.notemp * (x - Topt_mean.notemp)^2 + c_mean.notemp
}


# Data frame of temperature, fitted value, species
fitted_data.full.data.notemp <- data.frame(Temperature = rep(temp_seq, n.species.notemp), 
                                    Fitted = full.data(temp_seq),
                                    Species = factor(rep(1:n.species.notemp, each = length(temp_seq))))
# Original data with temperature, trait, and species as numeric
plot_data.full <- data.frame(
  Temperature = tick.abr.new$Temp,
  Trait = tick.abr.new$Trait,
  Species = as.numeric(tick.abr.new$Species)
)

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.full, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.full.data.notemp, aes(x = Temperature, y = Fitted, color = Species), size = 1) +
  labs(title = "Observed Data and Quadratic Fit",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 

# Ensure original data and fitted data are aligned correctly
# Make sure the species column in both datasets uses the same format (numeric)
fitted_data.full.data.notemp$Species <- as.numeric(fitted_data.full.data.notemp$Species)
tick$Species.number <- as.numeric(tick$Species)
# Merge the observed and predicted data for comparison
merged_data.full.notemp <- fitted_data.full.data.notemp %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species.number"))

# Calculate RMSE for each species
rmse_by_full_data.notemp <- merged_data.full.notemp %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_full_data.notemp)
# Average RMSE across species
mean(rmse_by_full_data.notemp$RMSE)


################################### By Species- Take temperature out ###################################
# Model using species effect
sink("species_robin.txt")
cat("
model {
  # Priors
  tau1 ~ dgamma(2, 2)           # Tau for trait
  
  # Species-specific 
  mu.c ~ dnorm(7, 0.1)          # Mean for species effects c (vertex y)
  tau.c ~ dexp(0.01)             # Tau for species effect c (vertex y)
  
  mu.Topt ~ dnorm(35, 0.1)          # Mean for species effects Topt (vertex x)
  tau.Topt ~ dgamma(2, 2)             # Tau for species effect Topt (vertex x)
  mu.a ~ dgamma(2, 2)            # Mean for species effects a (width)

  for (j in 1:n.species.notemp) {
    mu.species.c[j] ~ dnorm(mu.c, tau.c)       # Species-specific vertex y
    mu.species.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Species-specific vertex x
    mu.species.a[j] ~ dexp(mu.a)    # Species-specific width
  }
  # Likelihood
  for (i in 1:N.obs.notemp) {
    mu[i] <- mu.species.a[species.notemp[i]] * (mu.species.Topt[species.notemp[i]] - temp[i])^2 + mu.species.c[species.notemp[i]]
    trait[i] ~ dnorm(mu[i], tau1) 
  }}
", file = "species_robin.txt")

# Initial values (arbitrary)
inits <- function() {
  list(
    mu.a = 1,             
    mu.c = 0,          
    tau.c = 1,
    tau1 = 1,
    mu.Topt = 20,
    tau.Topt = 1
  )
}

# Parameters to monitor
parameters <- c("mu.species.a", "mu.species.c", "tau.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt",
                "tau1", "mu.species.Topt")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Bundle all data in a list for JAGS
jag.data <- list(trait = tick.notemp$Trait, N.obs.notemp = N.obs.notemp, 
                 temp = tick.notemp$Temp, 
                 n.species.notemp = n.species.notemp, species.notemp = species.notemp)
# Fit model
species.notemp.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                    model.file = "species_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                    n.iter = ni, DIC = T, working.directory = getwd())
# Make an mcmc object
species.notemp.fit.mcmc <- as.mcmc(species.notemp.fit) 
closeAllConnections()
# DIC
species.notemp.fit$BUGSoutput$DIC

# Get posterior means for each species effect
summary_fit.species.notemp <- summary(species.notemp.fit.mcmc)$statistics
a_species_mean.notemp <- summary_fit.species.notemp[grep("^mu.species.a\\[", rownames(summary_fit.species.notemp)), "Mean"]
c_species_mean.notemp <- summary_fit.species.notemp[grep("^mu.species.c\\[", rownames(summary_fit.species.notemp)), "Mean"]
Topt_species_mean.notemp <- summary_fit.species.notemp[grep("^mu.species.Topt\\[", rownames(summary_fit.species.notemp)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp) ,
                max(tick.abr.new$Temp) , by = 0.1)

# Get predicted values for each species using species effect parameters
species1.notemp <- function (x){
  a_species_mean.notemp[1] * (x - Topt_species_mean.notemp[1])^2 + c_species_mean.notemp[1]
}
species2.notemp <- function (x){
  a_species_mean.notemp[2] * (x - Topt_species_mean.notemp[2])^2 + c_species_mean.notemp[2]
}
species3.notemp <- function (x){
  a_species_mean.notemp[3] * (x - Topt_species_mean.notemp[3])^2 + c_species_mean.notemp[3]
}
species4.notemp <- function (x){
  a_species_mean.notemp[4] * (x - Topt_species_mean.notemp[4])^2 + c_species_mean.notemp[4]
}
species5.notemp <- function (x){
  a_species_mean.notemp[5] * (x - Topt_species_mean.notemp[5])^2 + c_species_mean.notemp[5]
}
species6.notemp <- function (x){
  a_species_mean.notemp[6] * (x - Topt_species_mean.notemp[6])^2 + c_species_mean.notemp[6]
}
species7.notemp <- function (x){
  a_species_mean.notemp[7] * (x - Topt_species_mean.notemp[7])^2 + c_species_mean.notemp[7]
}
species8.notemp <- function (x){
  a_species_mean.notemp[8] * (x - Topt_species_mean.notemp[8])^2 + c_species_mean.notemp[8]
}
species9.notemp <- function (x){
  a_species_mean.notemp[9] * (x - Topt_species_mean.notemp[9])^2 + c_species_mean.notemp[9]
}

# Make data frame of temperature, fitted value for each species, and species
fitted_data.species.notemp <- data.frame(Temperature = rep(temp_seq, n.species.notemp), 
                                  Fitted = c(species1.notemp(temp_seq),
                                             species2.notemp(temp_seq),
                                             species3.notemp(temp_seq),
                                             species4.notemp(temp_seq),
                                             species5.notemp(temp_seq),
                                             species6.notemp(temp_seq),
                                             species7.notemp(temp_seq),
                                             species8.notemp(temp_seq),
                                             species9.notemp(temp_seq)),
                                  Species = factor(rep(1:n.species.notemp, each = length(temp_seq))))
# Data frame of original temperature, trait, and numeric species
plot_data.species <- data.frame(
  Temperature = tick.abr.new$Temp,
  Trait = tick.abr.new$Trait,
  Species = as.numeric(tick.abr.new$Species)
)

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.species, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.species.notemp, aes(x = Temperature, y = Fitted, color = Species), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Species",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 

# Ensure original data and fitted data are aligned correctly
# Make sure the species column in both datasets uses the same format (numeric)
fitted_data.species.notemp$Species <- as.numeric(fitted_data.species.notemp$Species)
tick$Species.number <- as.numeric(tick$Species)
# Merge the observed and predicted data for comparison
merged_data.species.notemp <- fitted_data.species.notemp %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species.number"))

# Calculate RMSE for each species
rmse_by_species.notemp <- merged_data.species.notemp %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  ).notemp

# View the results
print(rmse_by_species.notemp)
# Average RMSE over species
mean(rmse_by_species.notemp$RMSE)




################################### By Genus- Take temperature out ###################################
# Model using genus effect
sink("genus_robin.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Genus-specific 
  mu.c ~ dnorm(7, 0.1)          # Mean for genus effects c (vertex y)
  tau.c ~ dexp(0.01)             # Tau for genus effect c (vertex y)
  
  mu.Topt ~ dnorm(35, 0.1)          # Mean for genus effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for genus effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for genus effects a (width)

  for (j in 1:n.genus.notemp) {
    mu.genus.c[j] ~ dnorm(mu.c, tau.c)       # Genus-specific vertex y
    mu.genus.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Genus-specific vertex x
    mu.genus.a[j] ~ dexp(mu.a)    # Genus-specific width
  }
  # Likelihood
  for (i in 1:N.obs.notemp) {
    mu[i] <- mu.genus.a[genus.notemp[i]] * (mu.genus.Topt[genus.notemp[i]] - temp[i])^2 + mu.genus.c[genus.notemp[i]]
    trait[i] ~ dnorm(mu[i], tau1)  
  }}
", file = "genus_robin.txt")

# Initial values (arbitrary)
inits <- function() {
  list(
    mu.a = 1,             
    mu.c = 0,          
    tau.c = 1,
    tau1 = 1,
    mu.Topt = 20,
    tau.Topt = 1
  )
}

# Parameters to monitor
parameters <- c("mu.genus.a", "mu.genus.c", "tau.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt",
                "tau1", "mu.genus.Topt")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Bundle all data in a list for JAGS
jag.data <- list(trait = tick.notemp$Trait, N.obs.notemp = N.obs.notemp, 
                 temp = tick.notemp$Temp, 
                 n.genus.notemp = n.genus.notemp, genus.notemp = genus.notemp)
# Fit model
genus.notemp.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                  model.file = "genus_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                  n.iter = ni, DIC = T, working.directory = getwd())
# Make an mcmc object
genus.notemp.fit.mcmc <- as.mcmc(genus.notemp.fit)
closeAllConnections()
# DIC
genus.notemp.fit$BUGSoutput$DIC

# Get posterior means using parameters from genus effect
summary_fit.genus.notemp <- summary(genus.notemp.fit.mcmc)$statistics
a_genus_mean.notemp <- summary_fit.genus.notemp[grep("^mu.genus.a\\[", rownames(summary_fit.genus.notemp)), "Mean"]
c_genus_mean.notemp <- summary_fit.genus.notemp[grep("^mu.genus.c\\[", rownames(summary_fit.genus.notemp)), "Mean"]
Topt_genus_mean.notemp <- summary_fit.genus.notemp[grep("^mu.genus.Topt\\[", rownames(summary_fit.genus.notemp)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp) ,
                max(tick.abr.new$Temp) , by = 0.1)

# Get predicted values for each genus
genus1.notemp <- function (x){
  a_genus_mean.notemp[1] * (x - Topt_genus_mean.notemp[1])^2 + c_genus_mean.notemp[1]
}
genus2.notemp <- function (x){
  a_genus_mean.notemp[2] * (x - Topt_genus_mean.notemp[2])^2 + c_genus_mean.notemp[2]
}
genus3.notemp <- function (x){
  a_genus_mean.notemp[3] * (x - Topt_genus_mean.notemp[3])^2 + c_genus_mean.notemp[3]
}
genus4.notemp <- function (x){
  a_genus_mean.notemp[4] * (x - Topt_genus_mean.notemp[4])^2 + c_genus_mean.notemp[4]
}

# Data frame of temperature, fitted values, and genus
fitted_data.genus.notemp <- data.frame(Temperature = rep(temp_seq, n.genus.notemp), 
                                Fitted = c(genus1.notemp(temp_seq),
                                           genus2.notemp(temp_seq),
                                           genus3.notemp(temp_seq),
                                           genus4.notemp(temp_seq)),
                                Genus = factor(rep(1:n.genus.notemp, each = length(temp_seq))))
# Original temperature, trait, numeric genus, and numeric species to plot by species
plot_data.genus <- data.frame(
  Temperature = tick.abr.new$Temp,
  Trait = tick.abr.new$Trait,
  Genus = as.numeric(tick.abr.new$Genus),
  Species = as.numeric(tick.abr.new$Species)
)

# Get unique combos of species and genus to align the correct curve on the plot
species_to_genus <- unique(plot_data.genus[, c("Species", "Genus")])
species_to_genus$Genus <- as.factor(species_to_genus$Genus)
# Merge fitted genus curves with species-to-genus mapping
fitted_data.genus.notemp <- data.frame(
  Temperature = rep(temp_seq, n.genus.notemp),
  Fitted = c(genus1.notemp(temp_seq),
             genus2.notemp(temp_seq),
             genus3.notemp(temp_seq),
             genus4.notemp(temp_seq)),
  Genus = factor(rep(1:n.genus.notemp, each = length(temp_seq)))
) %>%
  left_join(species_to_genus, by = "Genus")

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.genus, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.genus.notemp, aes(x = Temperature, y = Fitted, color = Genus), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Genus",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 

# Ensure original data and fitted data are aligned correctly
# Make sure the species column in both datasets uses the same format (numeric)
fitted_data.genus.notemp$Genus <- as.numeric(fitted_data.genus.notemp$Genus)
tick$Genus.number <- as.numeric(tick$Genus)
tick$Species.number <- as.numeric(tick$Species)
# Merge the observed and fitted data for comparison
merged_data.genus.notemp <- fitted_data.genus.notemp %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species.number", "Genus" = "Genus.number"))

# Calculate RMSE for each species
rmse_by_genus.notemp <- merged_data.genus.notemp %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_genus.notemp)
# Avergae RMSE across species
mean(rmse_by_genus.notemp$RMSE)


################################### By Climate- Take temperature out ###################################
# Model using climate effect
sink("climate_robin.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Climate-specific 
  mu.c ~ dnorm(7, 0.1)          # Mean for climate effect c (vertex y)
  tau.c ~ dexp(0.01)             # Tau for climate effect c (vertex y)
  
  mu.Topt ~ dnorm(35, 0.1)          # Mean for climate effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for climate effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for climate effects a (width)

  for (j in 1:n.climate.notemp) {
    mu.climate.c[j] ~ dnorm(mu.c, tau.c)       # Climate-specific vertex y
    mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific vertex x
    mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific width
  }
  # Likelihood
  for (i in 1:N.obs.notemp) {
    mu[i] <- mu.climate.a[climate.notemp[i]] * (mu.climate.Topt[climate.notemp[i]] - temp[i])^2 + mu.climate.c[climate.notemp[i]]
    trait[i] ~ dnorm(mu[i], tau1)  
  }}
", file = "climate_robin.txt")

# Initial values (arbitrary)
inits <- function() {
  list(
    mu.a = 1,             
    mu.c = 0,          
    tau.c = 1,
    tau1 = 1,
    mu.Topt = 20,
    tau.Topt = 1
  )
}

# Parameters to monitor
parameters <- c("mu.climate.a", "mu.climate.c", "tau.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt",
                "tau1", "mu.climate.Topt")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Bundle all data in a list for JAGS
jag.data <- list(trait = tick.notemp$Trait, N.obs.notemp = N.obs.notemp, 
                 temp = tick.notemp$Temp, 
                 n.climate.notemp = n.climate.notemp, climate.notemp = climate.notemp)
# Fit model
climate.notemp.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                    model.file = "climate_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                    n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
climate.notemp.fit.mcmc <- as.mcmc(climate.notemp.fit) 
closeAllConnections()
# DIC
climate.notemp.fit$BUGSoutput$DIC

# Get posterior means for parameters for each climate group
summary_fit.climate.notemp <- summary(climate.notemp.fit.mcmc)$statistics
a_climate_mean.notemp <- summary_fit.climate.notemp[grep("^mu.climate.a\\[", rownames(summary_fit.climate.notemp)), "Mean"]
c_climate_mean.notemp <- summary_fit.climate.notemp[grep("^mu.climate.c\\[", rownames(summary_fit.climate.notemp)), "Mean"]
Topt_climate_mean.notemp <- summary_fit.climate.notemp[grep("^mu.climate.Topt\\[", rownames(summary_fit.climate.notemp)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp) ,
                max(tick.abr.new$Temp) , by = 0.1)

# Get predicted values for each climate group
climate1.notemp <- function (x){
  a_climate_mean.notemp[1] * (x - Topt_climate_mean.notemp[1])^2 + c_climate_mean.notemp[1]
}
climate2.notemp <- function (x){
  a_climate_mean.notemp[2] * (x - Topt_climate_mean.notemp[2])^2 + c_climate_mean.notemp[2]
}
climate3.notemp <- function (x){
  a_climate_mean.notemp[3] * (x - Topt_climate_mean.notemp[3])^2 + c_climate_mean.notemp[3]
}
climate4.notemp <- function (x){
  a_climate_mean.notemp[4] * (x - Topt_climate_mean.notemp[4])^2 + c_climate_mean.notemp[4]
}
# Data frame of temperature, fitted values, and climate
fitted_data.climate.notemp <- data.frame(Temperature = rep(temp_seq, n.climate.notemp), 
                                  Fitted = c(climate1.notemp(temp_seq),
                                             climate2.notemp(temp_seq),
                                             climate3.notemp(temp_seq),
                                             climate4.notemp(temp_seq)),
                                  Climate = factor(rep(1:n.climate.notemp, each = length(temp_seq))))
# Original data with temperature, trait, numeric climate, and numeric species
plot_data.climate <- data.frame(
  Temperature = tick.abr.new$Temp,
  Trait = tick.abr.new$Trait,
  Climate = as.numeric(tick.abr.new$Climate),
  Species = as.numeric(tick.abr.new$Species)
)

# Get unique combinations of species and climate
species_to_climate <- unique(plot_data.climate[, c("Species", "Climate")])  
species_to_climate$Climate <- as.factor(species_to_climate$Climate)
# Merge fitted climate curves with species-to-climate mapping
fitted_data.climate.notemp <- data.frame(
  Temperature = rep(temp_seq, n.climate.notemp),
  Fitted = c(climate1.notemp(temp_seq),
             climate2.notemp(temp_seq),
             climate3.notemp(temp_seq),
             climate4.notemp(temp_seq)),
  Climate = factor(rep(1:n.climate.notemp, each = length(temp_seq)))
) %>%
  left_join(species_to_climate, by = "Climate")

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.climate, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.climate.notemp, aes(x = Temperature, y = Fitted, color = Climate), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Climate",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 

# Ensure test and fitted_data are aligned correctly
# Make sure the species column in both datasets uses the same format (numeric)
fitted_data.climate.notemp$Climate <- as.numeric(fitted_data.climate.notemp$Climate)
tick$Climate.number <- as.numeric(tick$Climate)
# Merge the observed and predicted data for comparison
merged_data.climate.notemp <- fitted_data.climate.notemp %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species.number", "Climate" = "Climate.number"))

# Calculate RMSE for each species
rmse_by_climate.notemp <- merged_data.climate.notemp %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_climate.notemp)
# Average RMSE across species
mean(rmse_by_climate.notemp$RMSE)



################################### By Host- Take temperature out ###################################
# Model by typical number of host species
sink("host_robin.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Host-specific 
  mu.c ~ dnorm(7, 0.1)          # Mean for host effects c (vertex y)
  tau.c ~ dexp(0.01)             # Tau for host effect c (vertex y)
  
  mu.Topt ~ dnorm(35, 0.1)          # Mean for host effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for host effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for host effects a (width)

  for (j in 1:n.host.notemp) {
    mu.host.c[j] ~ dnorm(mu.c, tau.c)       # Host-specific vertex y
    mu.host.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Host-specific vertex x
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
  }
  # Likelihood
  for (i in 1:N.obs.notemp) {
    mu[i] <- mu.host.a[host.notemp[i]] * (mu.host.Topt[host.notemp[i]] - temp[i])^2 + mu.host.c[host.notemp[i]]
    trait[i] ~ dnorm(mu[i], tau1)  
  }}
", file = "host_robin.txt")

# Initial values (arbitrary)
inits <- function() {
  list(
    mu.a = 1,             
    mu.c = 0,          
    tau.c = 1,
    tau1 = 1,
    mu.Topt = 20,
    tau.Topt = 1
  )
}

# Parameters to monitor
parameters <- c("mu.host.a", "mu.host.c", "tau.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt",
                "tau1", "mu.host.Topt")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Bundle all data in a list for JAGS
jag.data <- list(trait = tick.notemp$Trait, N.obs.notemp = N.obs.notemp, 
                 temp = tick.notemp$Temp, 
                 n.host.notemp = n.host.notemp, host.notemp = host.notemp)
# Fit model
host.notemp.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                 model.file = "host_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                 n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.notemp.fit.mcmc <- as.mcmc(host.notemp.fit) 
closeAllConnections()
# DIC
host.notemp.fit$BUGSoutput$DIC

# Get posterior means for parameters for each host group
summary_fit.host.notemp <- summary(host.notemp.fit.mcmc)$statistics
a_host_mean.notemp <- summary_fit.host.notemp[grep("^mu.host.a\\[", rownames(summary_fit.host.notemp)), "Mean"]
c_host_mean.notemp <- summary_fit.host.notemp[grep("^mu.host.c\\[", rownames(summary_fit.host.notemp)), "Mean"]
Topt_host_mean.notemp <- summary_fit.host.notemp[grep("^mu.host.Topt\\[", rownames(summary_fit.host.notemp)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp) ,
                max(tick.abr.new$Temp) , by = 0.1)

# Get predicted values for each host group
host1.notemp <- function (x){
  a_host_mean.notemp[1] * (x - Topt_host_mean.notemp[1])^2 + c_host_mean.notemp[1]
}
host2.notemp <- function (x){
  a_host_mean.notemp[2] * (x - Topt_host_mean.notemp[2])^2 + c_host_mean.notemp[2]
}

host3.notemp <- function (x){
  a_host_mean.notemp[3] * (x - Topt_host_mean.notemp[3])^2 + c_host_mean.notemp[3]
}

# Data frame with temperature, fitted value, and host group
fitted_data.host.notemp <- data.frame(Temperature = rep(temp_seq, n.host.notemp), 
                               Fitted = c(host1.notemp(temp_seq),
                                          host2.notemp(temp_seq)),
                               Host = factor(rep(1:n.host.notemp, each = length(temp_seq))))
# Original data with temperature, trait, numeric host group, numeric species
plot_data.host <- data.frame(
  Temperature = tick.abr.new$Temp,
  Trait = tick.abr.new$Trait,
  Host = as.numeric(tick.abr.new$Host),
  Species = as.numeric(tick.abr.new$Species)
)

# Get unique combinations of species and host for correct plotting
species_to_host <- unique(plot_data.host[, c("Species", "Host")])  
species_to_host$Host <- as.factor(species_to_host$Host)
# Merge fitted host curves with species-to-host mapping
fitted_data.host.notemp <- data.frame(
  Temperature = rep(temp_seq, n.host.notemp),
  Fitted = c(host1.notemp(temp_seq),
             host2.notemp(temp_seq)),
  Host = factor(rep(1:n.host.notemp, each = length(temp_seq)))
) %>%
  left_join(species_to_host, by = "Host")

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.host, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.notemp, aes(x = Temperature, y = Fitted, color = Host), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Host",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 

# Ensure test and fitted_data are aligned correctly
# Make sure the species column in both datasets uses the same format (numeric)
fitted_data.host.notemp$Host <- as.numeric(fitted_data.host.notemp$Host)
tick$Host.number <- as.numeric(tick$Host)
tick$Species.number <- as.numeric(tick$Species)
# Merge the observed and predicted data for comparison
merged_data.host.notemp <- fitted_data.host.notemp %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species.number", "Host" = "Host.number"))

# Calculate RMSE for each species
rmse_by_host.notemp <- merged_data.host.notemp %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_host.notemp)
# Average RMSE across species
mean(rmse_by_host.notemp$RMSE)

################################### By Size- Take temperature out ###################################
# Model by size group
sink("size_robin.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Host-specific 
  mu.c ~ dnorm(7, 0.1)          # Mean for size effects c (vertex y)
  tau.c ~ dexp(0.01)             # Tau for size effect c (vertex y)
  
  mu.Topt ~ dnorm(35, 0.1)          # Mean for size effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for size effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for size effects a (width)

  for (j in 1:n.size.notemp) {
    mu.size.c[j] ~ dnorm(mu.c, tau.c)       # Size-specific vertex y
    mu.size.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Size-specific vertex x
    mu.size.a[j] ~ dexp(mu.a)    # Size-specific width
  }
  # Likelihood
  for (i in 1:N.obs.notemp) {
    mu[i] <- mu.size.a[size.notemp[i]] * (mu.size.Topt[size.notemp[i]] - temp[i])^2 + mu.size.c[size.notemp[i]]
    trait[i] ~ dnorm(mu[i], tau1)  
  }}
", file = "size_robin.txt")

# Initial values (arbitrary)
inits <- function() {
  list(
    mu.a = 1,             
    mu.c = 0,          
    tau.c = 1,
    tau1 = 1,
    mu.Topt = 20,
    tau.Topt = 1
  )
}

# Parameters to monitor
parameters <- c("mu.size.a", "mu.size.c", "tau.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt",
                "tau1", "mu.size.Topt")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Bundle all data in a list for JAGS
jag.data <- list(trait = tick.notemp$Trait, N.obs.notemp = N.obs.notemp, 
                 temp = tick.notemp$Temp, 
                 n.size.notemp = n.size.notemp, size.notemp = size.notemp)
# Fit model
size.notemp.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                 model.file = "size_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                 n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
size.notemp.fit.mcmc <- as.mcmc(size.notemp.fit)
closeAllConnections()
# DIC
size.notemp.fit$BUGSoutput$DIC

# Get posterior means for parameters for each size 
summary_fit.size.notemp <- summary(size.notemp.fit.mcmc)$statistics
a_size_mean.notemp <- summary_fit.size.notemp[grep("^mu.size.a\\[", rownames(summary_fit.size.notemp)), "Mean"]
c_size_mean.notemp <- summary_fit.size.notemp[grep("^mu.size.c\\[", rownames(summary_fit.size.notemp)), "Mean"]
Topt_size_mean.notemp <- summary_fit.size.notemp[grep("^mu.size.Topt\\[", rownames(summary_fit.size.notemp)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp) ,
                max(tick.abr.new$Temp) , by = 0.1)

# Get predicted values for each size category
size1.notemp <- function (x){
  a_size_mean.notemp[1] * (x - Topt_size_mean.notemp[1])^2 + c_size_mean.notemp[1]
}
size2.notemp <- function (x){
  a_size_mean.notemp[2] * (x - Topt_size_mean.notemp[2])^2 + c_size_mean.notemp[2]
}
size3.notemp <- function (x){
  a_size_mean.notemp[3] * (x - Topt_size_mean.notemp[3])^2 + c_size_mean.notemp[3]
}


# Make data frame with temperature, fitted value, and size category
fitted_data.size.notemp <- data.frame(Temperature = rep(temp_seq, n.size.notemp), 
                               Fitted = c(size1.notemp(temp_seq),
                                          size2.notemp(temp_seq),
                                          size3.notemp(temp_seq)),
                               Size = factor(rep(1:n.size.notemp, each = length(temp_seq))))
# Original data with temperature, trait, numeric size, and numeric species
plot_data.size <- data.frame(
  Temperature = tick.abr.new$Temp,
  Trait = tick.abr.new$Trait,
  Size = as.numeric(tick.abr.new$Size),
  Species = as.numeric(tick.abr.new$Species)
)

# Get unique combinations of species and size for correct fitted lines on plots
species_to_size <- unique(plot_data.size[, c("Species", "Size")])  
species_to_size$Size <- as.factor(species_to_size$Size)
# Merge fitted host curves with species-to-size mapping
fitted_data.size.notemp <- data.frame(
  Temperature = rep(temp_seq, n.size.notemp),
  Fitted = c(size1.notemp(temp_seq),
             size2.notemp(temp_seq),
             size3.notemp(temp_seq)),
  Size = factor(rep(1:n.size.notemp, each = length(temp_seq)))
) %>%
  left_join(species_to_size, by = "Size")

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.size, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.size.notemp, aes(x = Temperature, y = Fitted, color = Size), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Size",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 

# Ensure test and fitted data are aligned correctly
# Make sure the species column in both datasets uses the same format (numeric)
fitted_data.size.notemp$Size <- as.numeric(fitted_data.size.notemp$Size)
tick$Size.number <- as.numeric(tick$Size)
tick$Species.number <- as.numeric(tick$Species)
# Merge the observed and predicted data for comparison
merged_data.size.notemp <- fitted_data.size.notemp %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species.number", "Size" = "Size.number"))

# Calculate RMSE for each species
rmse_by_size.notemp <- merged_data.size.notemp %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_size.notemp)
# Average RMSE across species
mean(rmse_by_size.notemp$RMSE)

################################### By Continent- Take temperature out ###################################
# Model by number of continents tick is found on
sink("continent_robin.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Continent-specific 
  mu.c ~ dnorm(7, 0.1)          # Mean for continent effects c (vertex y)
  tau.c ~ dexp(0.01)             # Tau for continent effect c (vertex y)
  
  mu.Topt ~ dnorm(35, 0.1)          # Mean for continent effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for continent effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for continent effects a (width)

  for (j in 1:n.continent.notemp) {
    mu.continent.c[j] ~ dnorm(mu.c, tau.c)       # Continent-specific vertex y
    mu.continent.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Continent-specific vertex x
    mu.continent.a[j] ~ dexp(mu.a)    # Continent-specific width
  }
  # Likelihood
  for (i in 1:N.obs.notemp) {
    mu[i] <- mu.continent.a[continent.notemp[i]] * (mu.continent.Topt[continent.notemp[i]] - temp[i])^2 + mu.continent.c[continent.notemp[i]]
    trait[i] ~ dnorm(mu[i], tau1)  
  }}
", file = "continent_robin.txt")

# Initial values (arbitrary)
inits <- function() {
  list(
    mu.a = 1,             
    mu.c = 7,          
    tau.c = 1,
    tau1 = 1,
    mu.Topt = 35,
    tau.Topt = 1
  )
}

# Parameters to monitor
parameters <- c("mu.continent.a", "mu.continent.c", "tau.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt",
                "tau1", "mu.continent.Topt")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Bundle all data in a list for JAGS
jag.data <- list(trait = tick.notemp$Trait, N.obs.notemp = N.obs.notemp, 
                 temp = tick.notemp$Temp, 
                 n.continent.notemp = n.continent.notemp, continent.notemp = continent.notemp)
# Fit model
continent.notemp.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                      model.file = "continent_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                      n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
continent.notemp.fit.mcmc <- as.mcmc(continent.notemp.fit)
closeAllConnections()
# DIC
continent.notemp.fit$BUGSoutput$DIC

# Get posterior means for parameters for different continent groups
summary_fit.continent.notemp <- summary(continent.notemp.fit.mcmc)$statistics
a_continent_mean.notemp <- summary_fit.continent.notemp[grep("^mu.continent.a\\[", rownames(summary_fit.continent.notemp)), "Mean"]
c_continent_mean.notemp <- summary_fit.continent.notemp[grep("^mu.continent.c\\[", rownames(summary_fit.continent.notemp)), "Mean"]
Topt_continent_mean.notemp <- summary_fit.continent.notemp[grep("^mu.continent.Topt\\[", rownames(summary_fit.continent.notemp)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp) ,
                max(tick.abr.new$Temp) , by = 0.1)

# Get predicted values for each continent group
continent1.notemp <- function (x){
  a_continent_mean.notemp[1] * (x - Topt_continent_mean.notemp[1])^2 + c_continent_mean.notemp[1]
}
continent2.notemp <- function (x){
  a_continent_mean.notemp[2] * (x - Topt_continent_mean.notemp[2])^2 + c_continent_mean.notemp[2]
}
continent3.notemp <- function (x){
  a_continent_mean.notemp[3] * (x - Topt_continent_mean.notemp[3])^2 + c_continent_mean.notemp[3]
}

# Data frame of temperature, fitted value, and continent group
fitted_data.continent.notemp <- data.frame(Temperature = rep(temp_seq, n.continent.notemp), 
                                    Fitted = c(continent1.notemp(temp_seq),
                                               continent2.notemp(temp_seq),
                                               continent3.notemp(temp_seq)),
                                    Continent = factor(rep(1:n.continent.notemp, each = length(temp_seq))))
# Original data with temperature, trait, numeric continent, numeric species
plot_data.continent <- data.frame(
  Temperature = tick.abr.new$Temp,
  Trait = tick.abr.new$Trait,
  Continent = as.numeric(tick.abr.new$Continent),
  Species = as.numeric(tick.abr.new$Species)
)
# Get unique combinations of species and continent for correct plotting
species_to_continent <- unique(plot_data.continent[, c("Species", "Continent")]) 
species_to_continent$Continent <- as.factor(species_to_continent$Continent)
# Merge fitted continent curves with species-to-continent mapping
fitted_data.continent.notemp <- data.frame(
  Temperature = rep(temp_seq, n.continent.notemp),
  Fitted = c(continent1.notemp(temp_seq),
             continent2.notemp(temp_seq),
             continent3.notemp(temp_seq)),
  Continent = factor(rep(1:n.continent.notemp, each = length(temp_seq)))
) %>%
  left_join(species_to_continent, by = "Continent")

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.continent, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.continent.notemp, aes(x = Temperature, y = Fitted, color = Continent), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Continent",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 

# Ensure test and fitted data are aligned correctly
# Make sure the species column in both datasets uses the same format (numeric)
fitted_data.continent.notemp$Continent <- as.numeric(fitted_data.continent.notemp$Continent)
tick$Continent.number <- as.numeric(tick$Continent)
tick$Species.number <- as.numeric(tick$Species)
# Merge the observed and predicted data for comparison
merged_data.continent.notemp <- fitted_data.continent.notemp %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species.number", "Continent" = "Continent.number"))

# Calculate RMSE for each species
rmse_by_continent.notemp <- merged_data.continent.notemp %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_continent.notemp)
# Average RMSE across species
mean(rmse_by_continent.notemp$RMSE)
################################### RMSE, DIC- Take temperature out ###################################
# Table of average RMSE across species and DIC of the model
table.tick.notemp <- matrix(0, nrow = 7, ncol = 2)
row.names(table.tick.notemp) <- c("Full Data", "By Species", "By Genus", 
                           "By Climate", "By Host", "By Size", "By Continents")
colnames(table.tick.notemp) <- c("DIC", "RMSE")
table.tick.notemp[1, 1] <- full.data.notemp.fit$BUGSoutput$DIC
table.tick.notemp[2, 1] <- species.notemp.fit$BUGSoutput$DIC
table.tick.notemp[3, 1] <- genus.notemp.fit$BUGSoutput$DIC
table.tick.notemp[4, 1] <- climate.notemp.fit$BUGSoutput$DIC
table.tick.notemp[5, 1] <- host.notemp.fit$BUGSoutput$DIC
table.tick.notemp[6, 1] <- size.notemp.fit$BUGSoutput$DIC
table.tick.notemp[7, 1] <- continent.notemp.fit$BUGSoutput$DIC
table.tick.notemp[1, 2] <- mean(rmse_by_full_data.notemp$RMSE)
table.tick.notemp[2, 2] <- mean(rmse_by_species.notemp$RMSE)
table.tick.notemp[3, 2] <- mean(rmse_by_genus.notemp$RMSE)
table.tick.notemp[4, 2] <- mean(rmse_by_climate.notemp$RMSE)
table.tick.notemp[5, 2] <- mean(rmse_by_host.notemp$RMSE)
table.tick.notemp[6, 2] <- mean(rmse_by_size.notemp$RMSE)
table.tick.notemp[7, 2] <- mean(rmse_by_continent.notemp$RMSE)

table.tick.notemp

# plot.tick <- data.frame(
#   Method = c("Full Data", "By Species", "By Genus", 
#              "By Climate", "By Host", "By Size", "By Continents"),
#   RMSE = c(mean(rmse_by_full_data$RMSE), mean(rmse_by_species$RMSE),
#            mean(rmse_by_genus$RMSE), mean(rmse_by_climate$RMSE),
#            mean(rmse_by_host$RMSE), mean(rmse_by_size$RMSE),
#            mean(rmse_by_continent$RMSE))
# )
# 
# plot.tick.notemp <- data.frame(
#   Method = c("Full Data", "By Species", "By Genus", 
#              "By Climate", "By Host", "By Size", "By Continents"),
#   RMSE = c(mean(rmse_by_full_data.notemp$RMSE), mean(rmse_by_species.notemp$RMSE),
#            mean(rmse_by_genus.notemp$RMSE), mean(rmse_by_climate.notemp$RMSE),
#            mean(rmse_by_host.notemp$RMSE), mean(rmse_by_size.notemp$RMSE),
#            mean(rmse_by_continent.notemp$RMSE))
# )
# 
# combined.plot <- bind_rows(
#   plot.tick %>% mutate(Table = "Original"),
#   plot.tick.notemp %>% mutate(Table = "No Temp")
# )
# 
# # Create a bar plot with grouping
# ggplot(combined.plot, aes(x = Method, y = RMSE, fill = Table)) +
#   geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black") +
#   labs(
#     title = "Paired RMSE Barplot for Methods",
#     x = "Method",
#     y = "RMSE",
#     fill = "Table"
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )



################################### Full Data- Out of sample prediction ###################################
tick.rm <- tick.abr.new[tick.abr.new$Species != "Hyalomma schulzei", ]
N.obs.rm <- length(tick.rm$Trait)

# Model using full data
sink("full_data_robin.txt")
cat("
model {
  # Priors
   a ~ dexp(0.1)             # Width
  Topt ~ dnorm(30, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  c ~ dnorm(10, 0.001)      # Oviposition rate at Topt
  
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], tau1)  
  }}
", file = "full_data_robin.txt")

# Initial values
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
jag.data <- list(trait = tick.rm$Trait, N.obs.rm = N.obs.rm, 
                 temp = tick.rm$Temp)

# Fit fill model
full.data.schulzei.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                      model.file = "full_data_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                      n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
full.data.schulzei.fit.mcmc <- as.mcmc(full.data.schulzei.fit) 
closeAllConnections()
#DIC
full.data.rm.fit$BUGSoutput$DIC

# Get posterior means
summary_fit.rm <- summary(full.data.schulzei.fit.mcmc)$statistics
a_mean.rm <- summary_fit.rm["a", "Mean"]
c_mean.rm <- summary_fit.rm["c", "Mean"]
Topt_mean.rm <- summary_fit.rm["Topt", "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp), max(tick.abr.new$Temp), by = 0.1)

# Get predicted values
full.data.rm <- function (x){
  a_mean.rm * (x - Topt_mean.rm)^2 + c_mean.rm
}


# Data frame of temperature, fitted value, species
fitted_data.full.data.schulzei <- data.frame(Temperature = rep(temp_seq, n.species), 
                                    Fitted = full.data.rm(temp_seq),
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
p.schulzei <- ggplot() +
  geom_point(data = plot_data.full, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.full.data.schulzei, aes(x = Temperature, y = Fitted), color = "purple", size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Species",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 

# Ensure original data and fitted data are aligned correctly
# Make sure the species column in both datasets uses the same format (numeric)
#fitted_data.full.data.rm$Species <- as.numeric(fitted_data.full.data.rm$Species)
#tick$Species.number <- as.numeric(tick$Species)
# Merge the observed and predicted data for comparison
merged_data.full.schulzei <- fitted_data.full.data.schulzei %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species"))

# Calculate RMSE for each species
rmse_by_full_data.schulzei <- merged_data.full.schulzei %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_full_data.schulzei)
# Average RMSE across species
mean(rmse_by_full_data.schulzei$RMSE)

save(rmse_by_full_data.schulzei,
     fitted_data.full.data.schulzei,
     merged_data.full.schulzei,
     full.data.schulzei.fit.mcmc,
     full.data.schulzei.fit,
     p.schulzei,
     file="schulzei.full.data.RData")

################################### By Genus- Out of sample prediction ###################################
tick.rm <- tick.abr.new[tick.abr.new$Species != "Hyalomma schulzei", ]
genus.rm <- as.factor(tick.rm$Genus)
n.genus.rm <- length(unique(tick.rm$Genus))
N.obs.rm <- length(tick.rm$Trait)
unique(tick.rm$Species)
# Model using genus effect
sink("genus_robin.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Genus-specific 
  mu.c ~ dnorm(7, 0.1)          # Mean for genus effects c (vertex y)
  tau.c ~ dexp(0.01)             # Tau for genus effect c (vertex y)
  
  mu.Topt ~ dnorm(35, 0.1)          # Mean for genus effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for genus effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for genus effects a (width)

  for (j in 1:n.genus.rm) {
    mu.genus.c[j] ~ dnorm(mu.c, tau.c)       # Genus-specific vertex y
    mu.genus.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Genus-specific vertex x
    mu.genus.a[j] ~ dexp(mu.a)    # Genus-specific width
  }
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.genus.a[genus.rm[i]] * (mu.genus.Topt[genus.rm[i]] - temp[i])^2 + mu.genus.c[genus.rm[i]]
    trait[i] ~ dnorm(mu[i], tau1)  
  }}
", file = "genus_robin.txt")

# Initial values (arbitrary)
inits <- function() {
  list(
    mu.a = 1,             
    mu.c = 0,          
    tau.c = 1,
    tau1 = 1,
    mu.Topt = 20,
    tau.Topt = 1
  )
}

# Parameters to monitor
parameters <- c("mu.genus.a", "mu.genus.c", "tau.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt",
                "tau1", "mu.genus.Topt")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Bundle all data in a list for JAGS
jag.data <- list(trait = tick.rm$Trait, N.obs.rm = N.obs.rm, 
                 temp = tick.rm$Temp, 
                 n.genus.rm = n.genus.rm, genus.rm = genus.rm)
# Fit model
genus.schulzei.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                  model.file = "genus_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                  n.iter = ni, DIC = T, working.directory = getwd())
# Make an mcmc object
genus.schulzei.fit.mcmc <- as.mcmc(genus.schulzei.fit)
closeAllConnections()
# DIC
genus.rm.fit$BUGSoutput$DIC

# Get posterior means using parameters from genus effect
summary_fit.genus.rm <- summary(genus.schulzei.fit.mcmc)$statistics
a_genus_mean.rm <- summary_fit.genus.rm[grep("^mu.genus.a\\[", rownames(summary_fit.genus.rm)), "Mean"]
c_genus_mean.rm <- summary_fit.genus.rm[grep("^mu.genus.c\\[", rownames(summary_fit.genus.rm)), "Mean"]
Topt_genus_mean.rm <- summary_fit.genus.rm[grep("^mu.genus.Topt\\[", rownames(summary_fit.genus.rm)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp) ,
                max(tick.abr.new$Temp) , by = 0.1)

# Get predicted values for each genus
genus1.rm <- function (x){
  a_genus_mean.rm[1] * (x - Topt_genus_mean.rm[1])^2 + c_genus_mean.rm[1]
}
genus2.rm <- function (x){
  a_genus_mean.rm[2] * (x - Topt_genus_mean.rm[2])^2 + c_genus_mean.rm[2]
}
genus3.rm <- function (x){
  a_genus_mean.rm[3] * (x - Topt_genus_mean.rm[3])^2 + c_genus_mean.rm[3]
}
genus4.rm <- function (x){
  a_genus_mean.rm[4] * (x - Topt_genus_mean.rm[4])^2 + c_genus_mean.rm[4]
}

# Data frame of temperature, fitted values, and genus
fitted_data.genus.schulzei <- data.frame(Temperature = rep(temp_seq, n.genus.rm), 
                                Fitted = c(genus1.rm(temp_seq),
                                           genus2.rm(temp_seq),
                                           genus3.rm(temp_seq),
                                           genus4.rm(temp_seq)),
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
fitted_data.genus.schulzei <- data.frame(
  Temperature = rep(temp_seq, n.genus.rm),
  Fitted = c(genus1.rm(temp_seq),
             genus2.rm(temp_seq),
             genus3.rm(temp_seq),
             genus4.rm(temp_seq)),
  Genus = factor(rep(c("Amblyomma", "Dermacentor",
                       "Haemaphysalis", "Hyalomma"), each = length(temp_seq)))
) %>%
  left_join(species_to_genus, by = "Genus")

# Plot observed data and fitted curves with facet wrap by species
p.genus.schulzei <- ggplot() +
  geom_point(data = plot_data.genus, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.genus.schulzei, aes(x = Temperature, y = Fitted, color = Genus), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Genus",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~ Species) 

# Ensure original data and fitted data are aligned correctly
# Make sure the species column in both datasets uses the same format (numeric)
#fitted_data.genus.rm$Genus <- as.numeric(fitted_data.genus.rm$Genus)
#tick$Genus.number <- as.numeric(tick$Genus)
#tick$Species.number <- as.numeric(tick$Species)
# Merge the observed and fitted data for comparison
merged_data.genus.schulzei <- fitted_data.genus.schulzei %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Genus" = "Genus"))

# Calculate RMSE for each species
rmse_by_genus.schulzei <- merged_data.genus.schulzei %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_genus.schulzei)
# Avergae RMSE across species
mean(rmse_by_genus.schulzei$RMSE)

save(rmse_by_genus.schulzei,
     fitted_data.genus.schulzei,
     merged_data.genus.schulzei,
     genus.schulzei.fit.mcmc,
     genus.schulzei.fit,
     p.genus.schulzei,
     file="schulzei.genus.RData")


################################### By Climate- Out of sample prediction ###################################
tick.rm <- tick.abr.new[tick.abr.new$Species != "Hyalomma schulzei", ]
climate.rm <- as.factor(tick.rm$Climate)
n.climate.rm <- length(unique(tick.rm$Climate))
N.obs.rm <- length(tick.rm$Trait)
unique(tick.rm$Species)
# Model using climate effect
sink("climate_robin.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Climate-specific 
  mu.c ~ dnorm(7, 0.1)          # Mean for climate effect c (vertex y)
  tau.c ~ dexp(0.01)             # Tau for climate effect c (vertex y)
  
  mu.Topt ~ dnorm(35, 0.1)          # Mean for climate effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for climate effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for climate effects a (width)

  for (j in 1:n.climate.rm) {
    mu.climate.c[j] ~ dnorm(mu.c, tau.c)       # Climate-specific vertex y
    mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific vertex x
    mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific width
  }
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.climate.a[climate.rm[i]] * (mu.climate.Topt[climate.rm[i]] - temp[i])^2 + mu.climate.c[climate.rm[i]]
    trait[i] ~ dnorm(mu[i], tau1)  
  }}
", file = "climate_robin.txt")

# Initial values (arbitrary)
inits <- function() {
  list(
    mu.a = 1,             
    mu.c = 0,          
    tau.c = 1,
    tau1 = 1,
    mu.Topt = 20,
    tau.Topt = 1
  )
}

# Parameters to monitor
parameters <- c("mu.climate.a", "mu.climate.c", "tau.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt",
                "tau1", "mu.climate.Topt")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Bundle all data in a list for JAGS
jag.data <- list(trait = tick.rm$Trait, N.obs.rm = N.obs.rm, 
                 temp = tick.rm$Temp, 
                 n.climate.rm = n.climate.rm, climate.rm = climate.rm)
# Fit model
climate.schulzei.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                    model.file = "climate_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                    n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
climate.schulzei.fit.mcmc <- as.mcmc(climate.schulzei.fit) 
closeAllConnections()
# DIC
climate.rm.fit$BUGSoutput$DIC

# Get posterior means for parameters for each climate group
summary_fit.climate.rm <- summary(climate.schulzei.fit.mcmc)$statistics
a_climate_mean.rm <- summary_fit.climate.rm[grep("^mu.climate.a\\[", rownames(summary_fit.climate.rm)), "Mean"]
c_climate_mean.rm <- summary_fit.climate.rm[grep("^mu.climate.c\\[", rownames(summary_fit.climate.rm)), "Mean"]
Topt_climate_mean.rm <- summary_fit.climate.rm[grep("^mu.climate.Topt\\[", rownames(summary_fit.climate.rm)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp) ,
                max(tick.abr.new$Temp) , by = 0.1)

# Get predicted values for each climate group
climate1.rm <- function (x){
  a_climate_mean.rm[1] * (x - Topt_climate_mean.rm[1])^2 + c_climate_mean.rm[1]
}
climate2.rm <- function (x){
  a_climate_mean.rm[2] * (x - Topt_climate_mean.rm[2])^2 + c_climate_mean.rm[2]
}
climate3.rm <- function (x){
  a_climate_mean.rm[3] * (x - Topt_climate_mean.rm[3])^2 + c_climate_mean.rm[3]
}
climate4.rm <- function (x){
  a_climate_mean.rm[4] * (x - Topt_climate_mean.rm[4])^2 + c_climate_mean.rm[4]
}

# Data frame of temperature, fitted values, and climate
fitted_data.climate.schulzei <- data.frame(Temperature = rep(temp_seq, n.climate.rm), 
                                  Fitted = c(climate1.rm(temp_seq),
                                             climate2.rm(temp_seq),
                                             climate3.rm(temp_seq),
                                             climate4.rm(temp_seq)),
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
fitted_data.climate.schulzei <- data.frame(
  Temperature = rep(temp_seq, n.climate.rm),
  Fitted = c(climate1.rm(temp_seq),
             climate2.rm(temp_seq),
             climate3.rm(temp_seq),
             climate4.rm(temp_seq)),
  Climate = factor(rep(c("Mixed", "Subtropical",
                         "Temperate", "Tropical"), each = length(temp_seq)))
) %>%
  left_join(species_to_climate, by = "Climate")

# Plot observed data and fitted curves with facet wrap by species
p.climate.schulzei <- ggplot() +
  geom_point(data = plot_data.climate, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.climate.schulzei, aes(x = Temperature, y = Fitted, color = Climate), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Climate",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 
p.climate.schulzei

# Ensure test and fitted_data are aligned correctly
# Make sure the species column in both datasets uses the same format (numeric)
#fitted_data.climate.schulzei$Climate <- as.numeric(fitted_data.climate.rm$Climate)
#tick$Climate.number <- as.numeric(tick$Climate)
# Merge the observed and predicted data for comparison
merged_data.climate.schulzei <- fitted_data.climate.schulzei %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Climate" = "Climate"))

# Calculate RMSE for each species
rmse_by_climate.schulzei <- merged_data.climate.schulzei %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_climate.schulzei)
# Average RMSE across species
mean(rmse_by_climate.schulzei$RMSE)

save(rmse_by_climate.schulzei,
     fitted_data.climate.schulzei,
     merged_data.climate.schulzei,
     climate.schulzei.fit.mcmc,
     climate.schulzei.fit,
     p.climate.schulzei,
     file="schulzei.climate.RData")


################################### By Host - Out of sample prediction ###################################
tick.rm <- tick.abr.new[tick.abr.new$Species != "Hyalomma schulzei", ]
host.rm <- as.factor(tick.rm$Host)
n.host.rm <- length(unique(tick.rm$Host))
N.obs.rm <- length(tick.rm$Trait)
unique(tick.rm$Species)
# Model by typical number of host species
sink("host_robin.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Host-specific 
  mu.c ~ dnorm(7, 0.1)          # Mean for host effects c (vertex y)
  tau.c ~ dexp(0.01)             # Tau for host effect c (vertex y)
  
  mu.Topt ~ dnorm(35, 0.1)          # Mean for host effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for host effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for host effects a (width)

  for (j in 1:n.host.rm) {
    mu.host.c[j] ~ dnorm(mu.c, tau.c)       # Host-specific vertex y
    mu.host.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Host-specific vertex x
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
  }
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.host.a[host.rm[i]] * (mu.host.Topt[host.rm[i]] - temp[i])^2 + mu.host.c[host.rm[i]]
    trait[i] ~ dnorm(mu[i], tau1)  
  }}
", file = "host_robin.txt")

# Initial values (arbitrary)
inits <- function() {
  list(
    mu.a = 1,             
    mu.c = 0,          
    tau.c = 1,
    tau1 = 1,
    mu.Topt = 20,
    tau.Topt = 1
  )
}

# Parameters to monitor
parameters <- c("mu.host.a", "mu.host.c", "tau.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt",
                "tau1", "mu.host.Topt")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Bundle all data in a list for JAGS
jag.data <- list(trait = tick.rm$Trait, N.obs.rm = N.obs.rm, 
                 temp = tick.rm$Temp, 
                 n.host.rm = n.host.rm, host.rm = host.rm)
# Fit model
host.schulzei.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                 model.file = "host_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                 n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.schulzei.fit.mcmc <- as.mcmc(host.schulzei.fit) 
closeAllConnections()
# DIC
host.rm.fit$BUGSoutput$DIC

# Get posterior means for parameters for each host group
summary_fit.host.rm <- summary(host.schulzei.fit.mcmc)$statistics
a_host_mean.rm <- summary_fit.host.rm[grep("^mu.host.a\\[", rownames(summary_fit.host.rm)), "Mean"]
c_host_mean.rm <- summary_fit.host.rm[grep("^mu.host.c\\[", rownames(summary_fit.host.rm)), "Mean"]
Topt_host_mean.rm <- summary_fit.host.rm[grep("^mu.host.Topt\\[", rownames(summary_fit.host.rm)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp) ,
                max(tick.abr.new$Temp) , by = 0.1)

# Get predicted values for each host group
host1.rm <- function (x){
  a_host_mean.rm[1] * (x - Topt_host_mean.rm[1])^2 + c_host_mean.rm[1]
}
host2.rm <- function (x){
  a_host_mean.rm[2] * (x - Topt_host_mean.rm[2])^2 + c_host_mean.rm[2]
}
host3.rm <- function (x){
  a_host_mean.rm[3] * (x - Topt_host_mean.rm[3])^2 + c_host_mean.rm[3]
}


# Data frame with temperature, fitted value, and host group
fitted_data.host.schulzei <- data.frame(Temperature = rep(temp_seq, n.host.rm), 
                               Fitted = c(host1.rm(temp_seq),
                                          host2.rm(temp_seq),
                                          host3.rm(temp_seq)),
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
fitted_data.host.schulzei <- data.frame(
  Temperature = rep(temp_seq, n.host.rm),
  Fitted = c(host1.rm(temp_seq),
             host2.rm(temp_seq),
             host3.rm(temp_seq)),
  Host = factor(rep(c("One", "Three", "Two"), each = length(temp_seq)))
) %>%
  left_join(species_to_host, by = "Host")

# Plot observed data and fitted curves with facet wrap by species
p.host.schulzei <- ggplot() +
  geom_point(data = plot_data.host, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.schulzei, aes(x = Temperature, y = Fitted, color = Host), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Host",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 
p.host.schulzei

# Ensure test and fitted_data are aligned correctly
# Make sure the species column in both datasets uses the same format (numeric)
#fitted_data.host.rm$Host <- as.numeric(fitted_data.host.rm$Host)
#tick$Host.number <- as.numeric(tick$Host)
#tick$Species.number <- as.numeric(tick$Species)
# Merge the observed and predicted data for comparison
merged_data.host.schulzei <- fitted_data.host.schulzei %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host"))

# Calculate RMSE for each species
rmse_by_host.schulzei <- merged_data.host.schulzei %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_host.schulzei)
# Average RMSE across species
mean(rmse_by_host.schulzei$RMSE)

save(rmse_by_host.schulzei,
     fitted_data.host.schulzei,
     merged_data.host.schulzei,
     host.schulzei.fit.mcmc,
     host.schulzei.fit,
     p.host.schulzei,
     file="schulzei.host.RData")

################################### By Size- Out of sample prediction ###################################
tick.rm <- tick.abr.new[tick.abr.new$Species != "Hyalomma schulzei", ]
size.rm <- as.numeric(tick.rm$Size)
n.size.rm <- length(unique(tick.rm$Size))
N.obs.rm <- length(tick.rm$Trait)

# Model by size group
sink("size_robin.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Host-specific 
  mu.c ~ dnorm(7, 0.1)          # Mean for size effects c (vertex y)
  tau.c ~ dexp(0.01)             # Tau for size effect c (vertex y)
  
  mu.Topt ~ dnorm(35, 0.1)          # Mean for size effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for size effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for size effects a (width)

  for (j in 1:n.size.rm) {
    mu.size.c[j] ~ dnorm(mu.c, tau.c)       # Size-specific vertex y
    mu.size.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Size-specific vertex x
    mu.size.a[j] ~ dexp(mu.a)    # Size-specific width
  }
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.size.a[size.rm[i]] * (mu.size.Topt[size.rm[i]] - temp[i])^2 + mu.size.c[size.rm[i]]
    trait[i] ~ dnorm(mu[i], tau1)  
  }}
", file = "size_robin.txt")

# Initial values (arbitrary)
inits <- function() {
  list(
    mu.a = 1,             
    mu.c = 0,          
    tau.c = 1,
    tau1 = 1,
    mu.Topt = 20,
    tau.Topt = 1
  )
}

# Parameters to monitor
parameters <- c("mu.size.a", "mu.size.c", "tau.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt",
                "tau1", "mu.size.Topt")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Bundle all data in a list for JAGS
jag.data <- list(trait = tick.rm$Trait, N.obs.rm = N.obs.rm, 
                 temp = tick.rm$Temp, 
                 n.size.rm = n.size.rm, size.rm = size.rm)
# Fit model
size.rm.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                 model.file = "size_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                 n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
size.rm.fit.mcmc <- as.mcmc(size.rm.fit)
closeAllConnections()
# DIC
size.rm.fit$BUGSoutput$DIC

# Get posterior means for parameters for each size 
summary_fit.size.rm <- summary(size.rm.fit.mcmc)$statistics
a_size_mean.rm <- summary_fit.size.rm[grep("^mu.size.a\\[", rownames(summary_fit.size.rm)), "Mean"]
c_size_mean.rm <- summary_fit.size.rm[grep("^mu.size.c\\[", rownames(summary_fit.size.rm)), "Mean"]
Topt_size_mean.rm <- summary_fit.size.rm[grep("^mu.size.Topt\\[", rownames(summary_fit.size.rm)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp) ,
                max(tick.abr.new$Temp) , by = 0.1)

# Get predicted values for each size category
size1.rm <- function (x){
  a_size_mean.rm[1] * (x - Topt_size_mean.rm[1])^2 + c_size_mean.rm[1]
}
size2.rm <- function (x){
  a_size_mean.rm[2] * (x - Topt_size_mean.rm[2])^2 + c_size_mean.rm[2]
}
size3.rm <- function (x){
  a_size_mean.rm[3] * (x - Topt_size_mean.rm[3])^2 + c_size_mean.rm[3]
}


# Make data frame with temperature, fitted value, and size category
fitted_data.size.rm <- data.frame(Temperature = rep(temp_seq, n.size.rm), 
                               Fitted = c(size1.rm(temp_seq),
                                          size1.rm(temp_seq),
                                          size3.rm(temp_seq)),
                               Size = factor(rep(1:n.size.rm, each = length(temp_seq))))
# Original data with temperature, trait, numeric size, and numeric species
plot_data.size <- data.frame(
  Temperature = tick.abr.new$Temp,
  Trait = tick.abr.new$Trait,
  Size = as.numeric(tick.abr.new$Size),
  Species = as.numeric(tick.abr.new$Species)
)

# Get unique combinations of species and size for correct fitted lines on plots
species_to_size <- unique(plot_data.size[, c("Species", "Size")])  
species_to_size$Size <- as.factor(species_to_size$Size)
# Merge fitted host curves with species-to-size mapping
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
  geom_point(data = plot_data.size, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.size.rm, aes(x = Temperature, y = Fitted, color = Size), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Size",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 

# Ensure test and fitted data are aligned correctly
# Make sure the species column in both datasets uses the same format (numeric)
fitted_data.size.rm$Size <- as.numeric(fitted_data.size.rm$Size)
tick$Size.number <- as.numeric(tick$Size)
tick$Species.number <- as.numeric(tick$Species)
# Merge the observed and predicted data for comparison
merged_data.size.rm <- fitted_data.size.rm %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species.number", "Size" = "Size.number"))

# Calculate RMSE for each species
rmse_by_size.rm <- merged_data.size.rm %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_size.rm)
# Average RMSE across species
mean(rmse_by_size.rm$RMSE)
################################### By Continent- Out of sample prediction ###################################
tick.rm <- tick.abr.new[tick.abr.new$Species != "Hyalomma schulzei", ]
continent.rm <- as.numeric(tick.rm$Continent)
n.continent.rm <- length(unique(tick.rm$Continent))
N.obs.rm <- length(tick.rm$Trait)

# Model by number of continents tick is found on
sink("continent_robin.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Continent-specific 
  mu.c ~ dnorm(7, 0.1)          # Mean for continent effects c (vertex y)
  tau.c ~ dexp(0.01)             # Tau for continent effect c (vertex y)
  
  mu.Topt ~ dnorm(35, 0.1)          # Mean for continent effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for continent effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for continent effects a (width)

  for (j in 1:n.continent.rm) {
    mu.continent.c[j] ~ dnorm(mu.c, tau.c)       # Continent-specific vertex y
    mu.continent.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Continent-specific vertex x
    mu.continent.a[j] ~ dexp(mu.a)    # Continent-specific width
  }
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.continent.a[continent.rm[i]] * (mu.continent.Topt[continent.rm[i]] - temp[i])^2 + mu.continent.c[continent.rm[i]]
    trait[i] ~ dnorm(mu[i], tau1)  
  }}
", file = "continent_robin.txt")

# Initial values (arbitrary)
inits <- function() {
  list(
    mu.a = 1,             
    mu.c = 0,          
    tau.c = 1,
    tau1 = 1,
    mu.Topt = 20,
    tau.Topt = 1
  )
}

# Parameters to monitor
parameters <- c("mu.continent.a", "mu.continent.c", "tau.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt",
                "tau1", "mu.continent.Topt")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Bundle all data in a list for JAGS
jag.data <- list(trait = tick.rm$Trait, N.obs.rm = N.obs.rm, 
                 temp = tick.rm$Temp, 
                 n.continent.rm = n.continent.rm, continent.rm = continent.rm)
# Fit model
continent.rm.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                      model.file = "continent_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                      n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
continent.rm.fit.mcmc <- as.mcmc(continent.rm.fit)
closeAllConnections()
# DIC
continent.rm.fit$BUGSoutput$DIC

# Get posterior means for parameters for different continent groups
summary_fit.continent.rm <- summary(continent.rm.fit.mcmc)$statistics
a_continent_mean.rm <- summary_fit.continent.rm[grep("^mu.continent.a\\[", rownames(summary_fit.continent.rm)), "Mean"]
c_continent_mean.rm <- summary_fit.continent.rm[grep("^mu.continent.c\\[", rownames(summary_fit.continent.rm)), "Mean"]
Topt_continent_mean.rm <- summary_fit.continent.rm[grep("^mu.continent.Topt\\[", rownames(summary_fit.continent.rm)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp) ,
                max(tick.abr.new$Temp) , by = 0.1)

# Get predicted values for each continent group
continent1.rm <- function (x){
  a_continent_mean.rm[1] * (x - Topt_continent_mean.rm[1])^2 + c_continent_mean.rm[1]
}
continent2.rm <- function (x){
  a_continent_mean.rm[2] * (x - Topt_continent_mean.rm[2])^2 + c_continent_mean.rm[2]
}
continent3.rm <- function (x){
  a_continent_mean.rm[3] * (x - Topt_continent_mean.rm[3])^2 + c_continent_mean.rm[3]
}

# Data frame of temperature, fitted value, and continent group
fitted_data.continent.rm <- data.frame(Temperature = rep(temp_seq, n.continent.rm), 
                                    Fitted = c(continent1.rm(temp_seq),
                                               continent2.rm(temp_seq),
                                               continent3.rm(temp_seq)),
                                    Continent = factor(rep(1:n.continent.rm, each = length(temp_seq))))
# Original data with temperature, trait, numeric continent, numeric species
plot_data.continent <- data.frame(
  Temperature = tick.abr.new$Temp,
  Trait = tick.abr.new$Trait,
  Continent = as.numeric(tick.abr.new$Continent),
  Species = as.numeric(tick.abr.new$Species)
)
# Get unique combinations of species and continent for correct plotting
species_to_continent <- unique(plot_data.continent[, c("Species", "Continent")]) 
species_to_continent$Continent <- as.factor(species_to_continent$Continent)
# Merge fitted continent curves with species-to-continent mapping
fitted_data.continent.rm <- data.frame(
  Temperature = rep(temp_seq, n.continent.rm),
  Fitted = c(continent1.rm(temp_seq),
             continent2.rm(temp_seq),
             continent3.rm(temp_seq)),
  Continent = factor(rep(1:n.continent.rm, each = length(temp_seq)))
) %>%
  left_join(species_to_continent, by = "Continent")

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data.continent, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.continent.rm, aes(x = Temperature, y = Fitted, color = Continent), size = 1) +
  labs(title = "Observed Data and Quadratic Fit by Continent",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 

# Ensure test and fitted data are aligned correctly
# Make sure the species column in both datasets uses the same format (numeric)
fitted_data.continent.rm$Continent <- as.numeric(fitted_data.continent.rm$Continent)
tick$Continent.number <- as.numeric(tick$Continent)
tick$Species.number <- as.numeric(tick$Species)
# Merge the observed and predicted data for comparison
merged_data.continent.rm <- fitted_data.continent.rm %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species.number", "Continent" = "Continent.number"))

# Calculate RMSE for each species
rmse_by_continent.rm <- merged_data.continent.rm %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_continent.rm)
# Average RMSE across species
mean(rmse_by_continent.rm$RMSE)

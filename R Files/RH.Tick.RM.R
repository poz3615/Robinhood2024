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
################################### By Genus- Andersoni ###################################
tick.rm <- tick.abr.new[tick.abr.new$Species != "Dermacentor andersoni", ]
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
  mu.c ~ dexp(10)          # Mean for genus effects c (vertex y)
  mu.Topt ~ dnorm(30, 0.1)          # Mean for genus effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for genus effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for genus effects a (width)

  for (j in 1:n.genus.rm) {
    mu.genus.c[j] ~ dexp(mu.c)       # Genus-specific vertex y
    mu.genus.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Genus-specific vertex x
    mu.genus.a[j] ~ dexp(mu.a)    # Genus-specific width
    genus.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.genus.a[genus.rm[i]] * (mu.genus.Topt[genus.rm[i]] - temp[i])^2 + mu.genus.c[genus.rm[i]]
    trait[i] ~ dnorm(mu[i], genus.tau[genus.rm[i]]) T(0, ) 
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
jag.data <- list(trait = tick.rm$Trait, N.obs.rm = N.obs.rm, 
                 temp = tick.rm$Temp, 
                 n.genus.rm = n.genus.rm, genus.rm = genus.rm)
# Fit model
genus.andersoni.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                            model.file = "genus_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                            n.iter = ni, DIC = T, working.directory = getwd())
# Make an mcmc object
genus.andersoni.fit.mcmc <- as.mcmc(genus.andersoni.fit)
closeAllConnections()

# Get posterior means using parameters from genus effect
summary_fit.genus.rm <- summary(genus.andersoni.fit.mcmc)$statistics
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
fitted_data.genus.andersoni <- data.frame(Temperature = rep(temp_seq, n.genus.rm), 
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
fitted_data.genus.andersoni <- data.frame(
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
p.genus.andersoni <- ggplot() +
  geom_point(data = plot_data.genus, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.genus.andersoni, aes(x = Temperature, y = Fitted, color = Genus), size = 1) +
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
p.genus.andersoni

# Merge the observed and fitted data for comparison
merged_data.genus.andersoni <- fitted_data.genus.andersoni %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Genus" = "Genus"))

# Calculate MSE for each species
mse_by_genus.andersoni <- merged_data.genus.andersoni %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.genus.andersoni <- merged_data.genus.andersoni %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_genus.andersoni,
     fitted_data.genus.andersoni,
     merged_data.genus.andersoni,
     genus.andersoni.fit.mcmc,
     genus.andersoni.fit,
     p.genus.andersoni,
     file="andersoni.genus.RData")


################################### By Genus- Nitens ###################################
tick.rm <- tick.abr.new[tick.abr.new$Species != "Dermacentor nitens", ]
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
  mu.c ~ dexp(10)          # Mean for genus effects c (vertex y)
  mu.Topt ~ dnorm(30, 0.1)          # Mean for genus effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for genus effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for genus effects a (width)

  for (j in 1:n.genus.rm) {
    mu.genus.c[j] ~ dexp(mu.c)       # Genus-specific vertex y
    mu.genus.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Genus-specific vertex x
    mu.genus.a[j] ~ dexp(mu.a)    # Genus-specific width
    genus.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.genus.a[genus.rm[i]] * (mu.genus.Topt[genus.rm[i]] - temp[i])^2 + mu.genus.c[genus.rm[i]]
    trait[i] ~ dnorm(mu[i], genus.tau[genus.rm[i]]) T(0, ) 
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
jag.data <- list(trait = tick.rm$Trait, N.obs.rm = N.obs.rm, 
                 temp = tick.rm$Temp, 
                 n.genus.rm = n.genus.rm, genus.rm = genus.rm)
# Fit model
genus.nitens.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                         model.file = "genus_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                         n.iter = ni, DIC = T, working.directory = getwd())
# Make an mcmc object
genus.nitens.fit.mcmc <- as.mcmc(genus.nitens.fit)
closeAllConnections()

# Get posterior means using parameters from genus effect
summary_fit.genus.rm <- summary(genus.nitens.fit.mcmc)$statistics
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
fitted_data.genus.nitens <- data.frame(Temperature = rep(temp_seq, n.genus.rm), 
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
fitted_data.genus.nitens <- data.frame(
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
p.genus.nitens <- ggplot() +
  geom_point(data = plot_data.genus, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.genus.nitens, aes(x = Temperature, y = Fitted, color = Genus), size = 1) +
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
p.genus.nitens

# Merge the observed and fitted data for comparison
merged_data.genus.nitens <- fitted_data.genus.nitens %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Genus" = "Genus"))

# Calculate MSE for each species
mse_by_genus.nitens <- merged_data.genus.nitens %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.genus.nitens <- merged_data.genus.nitens %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_genus.nitens,
     fitted_data.genus.nitens,
     merged_data.genus.nitens,
     genus.nitens.fit.mcmc,
     genus.nitens.fit,
     p.genus.nitens,
     file="nitens.genus.RData")

################################### By Genus- Aegyptium ###################################
tick.rm <- tick.abr.new[tick.abr.new$Species != "Hyalomma aegyptium", ]
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
  mu.c ~ dexp(10)          # Mean for genus effects c (vertex y)
  mu.Topt ~ dnorm(30, 0.1)          # Mean for genus effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for genus effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for genus effects a (width)

  for (j in 1:n.genus.rm) {
    mu.genus.c[j] ~ dexp(mu.c)       # Genus-specific vertex y
    mu.genus.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Genus-specific vertex x
    mu.genus.a[j] ~ dexp(mu.a)    # Genus-specific width
    genus.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.genus.a[genus.rm[i]] * (mu.genus.Topt[genus.rm[i]] - temp[i])^2 + mu.genus.c[genus.rm[i]]
    trait[i] ~ dnorm(mu[i], genus.tau[genus.rm[i]]) T(0, ) 
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
jag.data <- list(trait = tick.rm$Trait, N.obs.rm = N.obs.rm, 
                 temp = tick.rm$Temp, 
                 n.genus.rm = n.genus.rm, genus.rm = genus.rm)
# Fit model
genus.aegyptium.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                            model.file = "genus_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                            n.iter = ni, DIC = T, working.directory = getwd())
# Make an mcmc object
genus.aegyptium.fit.mcmc <- as.mcmc(genus.aegyptium.fit)
closeAllConnections()

# Get posterior means using parameters from genus effect
summary_fit.genus.rm <- summary(genus.aegyptium.fit.mcmc)$statistics
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
fitted_data.genus.aegyptium <- data.frame(Temperature = rep(temp_seq, n.genus.rm), 
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
fitted_data.genus.aegyptium <- data.frame(
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
p.genus.aegyptium <- ggplot() +
  geom_point(data = plot_data.genus, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.genus.aegyptium, aes(x = Temperature, y = Fitted, color = Genus), size = 1) +
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
p.genus.aegyptium

# Merge the observed and fitted data for comparison
merged_data.genus.aegyptium <- fitted_data.genus.aegyptium %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Genus" = "Genus"))

# Calculate MSE for each species
mse_by_genus.aegyptium <- merged_data.genus.aegyptium %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.genus.aegyptium <- merged_data.genus.aegyptium %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_genus.aegyptium,
     fitted_data.genus.aegyptium,
     merged_data.genus.aegyptium,
     genus.aegyptium.fit.mcmc,
     genus.aegyptium.fit,
     p.genus.aegyptium,
     file="aegyptium.genus.RData")

################################### By Genus- Dromedarii ###################################
tick.rm <- tick.abr.new[tick.abr.new$Species != "Hyalomma dromedarii", ]
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
  mu.c ~ dexp(10)          # Mean for genus effects c (vertex y)
  mu.Topt ~ dnorm(30, 0.1)          # Mean for genus effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for genus effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for genus effects a (width)

  for (j in 1:n.genus.rm) {
    mu.genus.c[j] ~ dexp(mu.c)       # Genus-specific vertex y
    mu.genus.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Genus-specific vertex x
    mu.genus.a[j] ~ dexp(mu.a)    # Genus-specific width
    genus.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.genus.a[genus.rm[i]] * (mu.genus.Topt[genus.rm[i]] - temp[i])^2 + mu.genus.c[genus.rm[i]]
    trait[i] ~ dnorm(mu[i], genus.tau[genus.rm[i]]) T(0, ) 
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
jag.data <- list(trait = tick.rm$Trait, N.obs.rm = N.obs.rm, 
                 temp = tick.rm$Temp, 
                 n.genus.rm = n.genus.rm, genus.rm = genus.rm)
# Fit model
genus.dromedarii.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                             model.file = "genus_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                             n.iter = ni, DIC = T, working.directory = getwd())
# Make an mcmc object
genus.dromedarii.fit.mcmc <- as.mcmc(genus.dromedarii.fit)
closeAllConnections()

# Get posterior means using parameters from genus effect
summary_fit.genus.rm <- summary(genus.dromedarii.fit.mcmc)$statistics
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
fitted_data.genus.dromedarii <- data.frame(Temperature = rep(temp_seq, n.genus.rm), 
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
fitted_data.genus.dromedarii <- data.frame(
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
p.genus.dromedarii <- ggplot() +
  geom_point(data = plot_data.genus, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.genus.dromedarii, aes(x = Temperature, y = Fitted, color = Genus), size = 1) +
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
p.genus.dromedarii

# Merge the observed and fitted data for comparison
merged_data.genus.dromedarii <- fitted_data.genus.dromedarii %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Genus" = "Genus"))

# Calculate MSE for each species
mse_by_genus.dromedarii <- merged_data.genus.dromedarii %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.genus.dromedarii <- merged_data.genus.dromedarii %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_genus.dromedarii,
     fitted_data.genus.dromedarii,
     merged_data.genus.dromedarii,
     genus.dromedarii.fit.mcmc,
     genus.dromedarii.fit,
     p.genus.dromedarii,
     file="dromedarii.genus.RData")

################################### By Genus- Impeltatum ###################################
tick.rm <- tick.abr.new[tick.abr.new$Species != "Hyalomma impeltatum", ]
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
  mu.c ~ dexp(10)          # Mean for genus effects c (vertex y)
  mu.Topt ~ dnorm(30, 0.1)          # Mean for genus effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for genus effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for genus effects a (width)

  for (j in 1:n.genus.rm) {
    mu.genus.c[j] ~ dexp(mu.c)       # Genus-specific vertex y
    mu.genus.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Genus-specific vertex x
    mu.genus.a[j] ~ dexp(mu.a)    # Genus-specific width
    genus.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.genus.a[genus.rm[i]] * (mu.genus.Topt[genus.rm[i]] - temp[i])^2 + mu.genus.c[genus.rm[i]]
    trait[i] ~ dnorm(mu[i], genus.tau[genus.rm[i]]) T(0, ) 
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
jag.data <- list(trait = tick.rm$Trait, N.obs.rm = N.obs.rm, 
                 temp = tick.rm$Temp, 
                 n.genus.rm = n.genus.rm, genus.rm = genus.rm)
# Fit model
genus.impeltatum.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                             model.file = "genus_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                             n.iter = ni, DIC = T, working.directory = getwd())
# Make an mcmc object
genus.impeltatum.fit.mcmc <- as.mcmc(genus.impeltatum.fit)
closeAllConnections()

# Get posterior means using parameters from genus effect
summary_fit.genus.rm <- summary(genus.impeltatum.fit.mcmc)$statistics
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
fitted_data.genus.impeltatum <- data.frame(Temperature = rep(temp_seq, n.genus.rm), 
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
fitted_data.genus.impeltatum <- data.frame(
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
p.genus.impeltatum <- ggplot() +
  geom_point(data = plot_data.genus, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.genus.impeltatum, aes(x = Temperature, y = Fitted, color = Genus), size = 1) +
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
p.genus.impeltatum

# Merge the observed and fitted data for comparison
merged_data.genus.impeltatum <- fitted_data.genus.impeltatum %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Genus" = "Genus"))

# Calculate MSE for each species
mse_by_genus.impeltatum <- merged_data.genus.impeltatum %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.genus.impeltatum <- merged_data.genus.impeltatum %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_genus.impeltatum,
     fitted_data.genus.impeltatum,
     merged_data.genus.impeltatum,
     genus.impeltatum.fit.mcmc,
     genus.impeltatum.fit,
     p.genus.impeltatum,
     file="impeltatum.genus.RData")

################################### By Genus- Lusitanicum ###################################
tick.rm <- tick.abr.new[tick.abr.new$Species != "Hyalomma lusitanicum", ]
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
  mu.c ~ dexp(10)          # Mean for genus effects c (vertex y)
  mu.Topt ~ dnorm(30, 0.1)          # Mean for genus effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for genus effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for genus effects a (width)

  for (j in 1:n.genus.rm) {
    mu.genus.c[j] ~ dexp(mu.c)       # Genus-specific vertex y
    mu.genus.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Genus-specific vertex x
    mu.genus.a[j] ~ dexp(mu.a)    # Genus-specific width
    genus.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.genus.a[genus.rm[i]] * (mu.genus.Topt[genus.rm[i]] - temp[i])^2 + mu.genus.c[genus.rm[i]]
    trait[i] ~ dnorm(mu[i], genus.tau[genus.rm[i]]) T(0, ) 
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
jag.data <- list(trait = tick.rm$Trait, N.obs.rm = N.obs.rm, 
                 temp = tick.rm$Temp, 
                 n.genus.rm = n.genus.rm, genus.rm = genus.rm)
# Fit model
genus.lusitanicum.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                              model.file = "genus_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                              n.iter = ni, DIC = T, working.directory = getwd())
# Make an mcmc object
genus.lusitanicum.fit.mcmc <- as.mcmc(genus.lusitanicum.fit)
closeAllConnections()

# Get posterior means using parameters from genus effect
summary_fit.genus.rm <- summary(genus.lusitanicum.fit.mcmc)$statistics
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
fitted_data.genus.lusitanicum <- data.frame(Temperature = rep(temp_seq, n.genus.rm), 
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
fitted_data.genus.lusitanicum <- data.frame(
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
p.genus.lusitanicum <- ggplot() +
  geom_point(data = plot_data.genus, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.genus.lusitanicum, aes(x = Temperature, y = Fitted, color = Genus), size = 1) +
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
p.genus.lusitanicum

# Merge the observed and fitted data for comparison
merged_data.genus.lusitanicum <- fitted_data.genus.lusitanicum %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Genus" = "Genus"))

# Calculate MSE for each species
mse_by_genus.lusitanicum <- merged_data.genus.lusitanicum %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.genus.lusitanicum <- merged_data.genus.lusitanicum %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_genus.lusitanicum,
     fitted_data.genus.lusitanicum,
     merged_data.genus.lusitanicum,
     genus.lusitanicum.fit.mcmc,
     genus.lusitanicum.fit,
     p.genus.lusitanicum,
     file="lusitanicum.genus.RData")

################################### By Genus- Schulzei ###################################
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
  mu.c ~ dexp(10)          # Mean for genus effects c (vertex y)
  mu.Topt ~ dnorm(30, 0.1)          # Mean for genus effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for genus effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for genus effects a (width)

  for (j in 1:n.genus.rm) {
    mu.genus.c[j] ~ dexp(mu.c)       # Genus-specific vertex y
    mu.genus.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Genus-specific vertex x
    mu.genus.a[j] ~ dexp(mu.a)    # Genus-specific width
    genus.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.genus.a[genus.rm[i]] * (mu.genus.Topt[genus.rm[i]] - temp[i])^2 + mu.genus.c[genus.rm[i]]
    trait[i] ~ dnorm(mu[i], genus.tau[genus.rm[i]]) T(0, ) 
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
p.genus.schulzei

# Merge the observed and fitted data for comparison
merged_data.genus.schulzei <- fitted_data.genus.schulzei %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Genus" = "Genus"))

# Calculate MSE for each species
mse_by_genus.schulzei <- merged_data.genus.schulzei %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.genus.schulzei <- merged_data.genus.schulzei %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_genus.schulzei,
     fitted_data.genus.schulzei,
     merged_data.genus.schulzei,
     genus.schulzei.fit.mcmc,
     genus.schulzei.fit,
     p.genus.schulzei,
     file="schulzei.genus.RData")


################################### By Host- Lepidum ###################################
tick.rm <- tick.abr.new[tick.abr.new$Species != "Amblyomma lepidum", ]
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
  mu.c ~ dexp(10)          # Mean for host effects c (vertex y)
  mu.Topt ~ dnorm(30, 0.1)          # Mean for host effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for host effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for host effects a (width)

  for (j in 1:n.host.rm) {
    mu.host.c[j] ~ dexp(mu.c)       # Host-specific vertex y
    mu.host.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Host-specific vertex x
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
    host.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.host.a[host.rm[i]] * (mu.host.Topt[host.rm[i]] - temp[i])^2 + mu.host.c[host.rm[i]]
    trait[i] ~ dnorm(mu[i], host.tau[host.rm[i]]) T(0, ) 
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
jag.data <- list(trait = tick.rm$Trait, N.obs.rm = N.obs.rm, 
                 temp = tick.rm$Temp, 
                 n.host.rm = n.host.rm, host.rm = host.rm)
# Fit model
host.lepidum.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                         model.file = "host_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                         n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.lepidum.fit.mcmc <- as.mcmc(host.lepidum.fit) 
closeAllConnections()


# Get posterior means for parameters for each host group
summary_fit.host.rm <- summary(host.lepidum.fit.mcmc)$statistics
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
fitted_data.host.lepidum <- data.frame(Temperature = rep(temp_seq, n.host.rm), 
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
fitted_data.host.lepidum <- data.frame(
  Temperature = rep(temp_seq, n.host.rm),
  Fitted = c(host1.rm(temp_seq),
             host2.rm(temp_seq),
             host3.rm(temp_seq)),
  Host = factor(rep(c("One", "Three", "Two"), each = length(temp_seq)))
) %>%
  left_join(species_to_host, by = "Host")

plot_data.host$Host <- factor(plot_data.host$Host, levels = c("One", "Two", "Three"))
fitted_data.host.lepidum$Host <- factor(fitted_data.host.lepidum$Host, levels = c("One", "Two", "Three"))
# Plot observed data and fitted curves with facet wrap by species
p.host.lepidum <- ggplot() +
  geom_point(data = plot_data.host, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.lepidum, aes(x = Temperature, y = Fitted, color = Host), size = 1) +
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
p.host.lepidum

# Merge the observed and predicted data for comparison
merged_data.host.lepidum <- fitted_data.host.lepidum %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host"))

# Calculate RMSE for each species
mse_by_host.lepidum <- merged_data.host.lepidum %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.host.lepidum <- merged_data.host.lepidum %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_host.lepidum,
     fitted_data.host.lepidum,
     merged_data.host.lepidum,
     host.lepidum.fit.mcmc,
     host.lepidum.fit,
     p.host.lepidum,
     file="lepidum.host.RData")

################################### By Host- Andersoni ###################################
tick.rm <- tick.abr.new[tick.abr.new$Species != "Dermacentor andersoni", ]
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
  mu.c ~ dexp(10)          # Mean for host effects c (vertex y)
  mu.Topt ~ dnorm(30, 0.1)          # Mean for host effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for host effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for host effects a (width)

  for (j in 1:n.host.rm) {
    mu.host.c[j] ~ dexp(mu.c)       # Host-specific vertex y
    mu.host.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Host-specific vertex x
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
    host.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.host.a[host.rm[i]] * (mu.host.Topt[host.rm[i]] - temp[i])^2 + mu.host.c[host.rm[i]]
    trait[i] ~ dnorm(mu[i], host.tau[host.rm[i]]) T(0, ) 
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
jag.data <- list(trait = tick.rm$Trait, N.obs.rm = N.obs.rm, 
                 temp = tick.rm$Temp, 
                 n.host.rm = n.host.rm, host.rm = host.rm)
# Fit model
host.andersoni.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                           model.file = "host_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                           n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.andersoni.fit.mcmc <- as.mcmc(host.andersoni.fit) 
closeAllConnections()


# Get posterior means for parameters for each host group
summary_fit.host.rm <- summary(host.andersoni.fit.mcmc)$statistics
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
fitted_data.host.andersoni <- data.frame(Temperature = rep(temp_seq, n.host.rm), 
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
fitted_data.host.andersoni <- data.frame(
  Temperature = rep(temp_seq, n.host.rm),
  Fitted = c(host1.rm(temp_seq),
             host2.rm(temp_seq),
             host3.rm(temp_seq)),
  Host = factor(rep(c("One", "Three", "Two"), each = length(temp_seq)))
) %>%
  left_join(species_to_host, by = "Host")

plot_data.host$Host <- factor(plot_data.host$Host, levels = c("One", "Two", "Three"))
fitted_data.host.andersoni$Host <- factor(fitted_data.host.andersoni$Host, levels = c("One", "Two", "Three"))
# Plot observed data and fitted curves with facet wrap by species
p.host.andersoni <- ggplot() +
  geom_point(data = plot_data.host, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.andersoni, aes(x = Temperature, y = Fitted, color = Host), size = 1) +
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
p.host.andersoni

# Merge the observed and predicted data for comparison
merged_data.host.andersoni <- fitted_data.host.andersoni %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host"))

# Calculate RMSE for each species
mse_by_host.andersoni <- merged_data.host.andersoni %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.host.andersoni <- merged_data.host.andersoni %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_host.andersoni,
     fitted_data.host.andersoni,
     merged_data.host.andersoni,
     host.andersoni.fit.mcmc,
     host.andersoni.fit,
     p.host.andersoni,
     file="andersoni.host.RData")

################################### By Host- Leporispalustris ###################################
tick.rm <- tick.abr.new[tick.abr.new$Species != "Haemaphysalis leporispalustris", ]
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
  mu.c ~ dexp(10)          # Mean for host effects c (vertex y)
  mu.Topt ~ dnorm(30, 0.1)          # Mean for host effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for host effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for host effects a (width)

  for (j in 1:n.host.rm) {
    mu.host.c[j] ~ dexp(mu.c)       # Host-specific vertex y
    mu.host.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Host-specific vertex x
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
    host.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.host.a[host.rm[i]] * (mu.host.Topt[host.rm[i]] - temp[i])^2 + mu.host.c[host.rm[i]]
    trait[i] ~ dnorm(mu[i], host.tau[host.rm[i]]) T(0, ) 
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
jag.data <- list(trait = tick.rm$Trait, N.obs.rm = N.obs.rm, 
                 temp = tick.rm$Temp, 
                 n.host.rm = n.host.rm, host.rm = host.rm)
# Fit model
host.leporispalustris.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                                  model.file = "host_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                                  n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.leporispalustris.fit.mcmc <- as.mcmc(host.leporispalustris.fit) 
closeAllConnections()


# Get posterior means for parameters for each host group
summary_fit.host.rm <- summary(host.leporispalustris.fit.mcmc)$statistics
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
fitted_data.host.leporispalustris <- data.frame(Temperature = rep(temp_seq, n.host.rm), 
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
fitted_data.host.leporispalustris <- data.frame(
  Temperature = rep(temp_seq, n.host.rm),
  Fitted = c(host1.rm(temp_seq),
             host2.rm(temp_seq),
             host3.rm(temp_seq)),
  Host = factor(rep(c("One", "Three", "Two"), each = length(temp_seq)))
) %>%
  left_join(species_to_host, by = "Host")

plot_data.host$Host <- factor(plot_data.host$Host, levels = c("One", "Two", "Three"))
fitted_data.host.leporispalustris$Host <- factor(fitted_data.host.leporispalustris$Host, levels = c("One", "Two", "Three"))
# Plot observed data and fitted curves with facet wrap by species
p.host.leporispalustris <- ggplot() +
  geom_point(data = plot_data.host, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.leporispalustris, aes(x = Temperature, y = Fitted, color = Host), size = 1) +
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
p.host.leporispalustris

# Merge the observed and predicted data for comparison
merged_data.host.leporispalustris <- fitted_data.host.leporispalustris %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host"))

# Calculate RMSE for each species
mse_by_host.leporispalustris <- merged_data.host.leporispalustris %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.host.leporispalustris <- merged_data.host.leporispalustris %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_host.leporispalustris,
     fitted_data.host.leporispalustris,
     merged_data.host.leporispalustris,
     host.leporispalustris.fit.mcmc,
     host.leporispalustris.fit,
     p.host.leporispalustris,
     file="leporispalustris.host.RData")

################################### By Host- Aegyptium ###################################
tick.rm <- tick.abr.new[tick.abr.new$Species != "Hyalomma aegyptium", ]
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
  mu.c ~ dexp(10)          # Mean for host effects c (vertex y)
  mu.Topt ~ dnorm(30, 0.1)          # Mean for host effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for host effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for host effects a (width)

  for (j in 1:n.host.rm) {
    mu.host.c[j] ~ dexp(mu.c)       # Host-specific vertex y
    mu.host.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Host-specific vertex x
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
    host.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.host.a[host.rm[i]] * (mu.host.Topt[host.rm[i]] - temp[i])^2 + mu.host.c[host.rm[i]]
    trait[i] ~ dnorm(mu[i], host.tau[host.rm[i]]) T(0, ) 
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
jag.data <- list(trait = tick.rm$Trait, N.obs.rm = N.obs.rm, 
                 temp = tick.rm$Temp, 
                 n.host.rm = n.host.rm, host.rm = host.rm)
# Fit model
host.aegyptium.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                           model.file = "host_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                           n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.aegyptium.fit.mcmc <- as.mcmc(host.aegyptium.fit) 
closeAllConnections()


# Get posterior means for parameters for each host group
summary_fit.host.rm <- summary(host.aegyptium.fit.mcmc)$statistics
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
fitted_data.host.aegyptium <- data.frame(Temperature = rep(temp_seq, n.host.rm), 
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
fitted_data.host.aegyptium <- data.frame(
  Temperature = rep(temp_seq, n.host.rm),
  Fitted = c(host1.rm(temp_seq),
             host2.rm(temp_seq),
             host3.rm(temp_seq)),
  Host = factor(rep(c("One", "Three", "Two"), each = length(temp_seq)))
) %>%
  left_join(species_to_host, by = "Host")

plot_data.host$Host <- factor(plot_data.host$Host, levels = c("One", "Two", "Three"))
fitted_data.host.aegyptium$Host <- factor(fitted_data.host.aegyptium$Host, levels = c("One", "Two", "Three"))
# Plot observed data and fitted curves with facet wrap by species
p.host.aegyptium <- ggplot() +
  geom_point(data = plot_data.host, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.aegyptium, aes(x = Temperature, y = Fitted, color = Host), size = 1) +
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
p.host.aegyptium

# Merge the observed and predicted data for comparison
merged_data.host.aegyptium <- fitted_data.host.aegyptium %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host"))

# Calculate RMSE for each species
mse_by_host.aegyptium <- merged_data.host.aegyptium %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.host.aegyptium <- merged_data.host.aegyptium %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_host.aegyptium,
     fitted_data.host.aegyptium,
     merged_data.host.aegyptium,
     host.aegyptium.fit.mcmc,
     host.aegyptium.fit,
     p.host.aegyptium,
     file="aegyptium.host.RData")

################################### By Host- Dromedarii ###################################
tick.rm <- tick.abr.new[tick.abr.new$Species != "Hyalomma dromedarii", ]
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
  mu.c ~ dexp(10)          # Mean for host effects c (vertex y)
  mu.Topt ~ dnorm(30, 0.1)          # Mean for host effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for host effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for host effects a (width)

  for (j in 1:n.host.rm) {
    mu.host.c[j] ~ dexp(mu.c)       # Host-specific vertex y
    mu.host.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Host-specific vertex x
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
    host.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.host.a[host.rm[i]] * (mu.host.Topt[host.rm[i]] - temp[i])^2 + mu.host.c[host.rm[i]]
    trait[i] ~ dnorm(mu[i], host.tau[host.rm[i]]) T(0, ) 
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
jag.data <- list(trait = tick.rm$Trait, N.obs.rm = N.obs.rm, 
                 temp = tick.rm$Temp, 
                 n.host.rm = n.host.rm, host.rm = host.rm)
# Fit model
host.dromedarii.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                            model.file = "host_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                            n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.dromedarii.fit.mcmc <- as.mcmc(host.dromedarii.fit) 
closeAllConnections()


# Get posterior means for parameters for each host group
summary_fit.host.rm <- summary(host.dromedarii.fit.mcmc)$statistics
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
fitted_data.host.dromedarii <- data.frame(Temperature = rep(temp_seq, n.host.rm), 
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
fitted_data.host.dromedarii <- data.frame(
  Temperature = rep(temp_seq, n.host.rm),
  Fitted = c(host1.rm(temp_seq),
             host2.rm(temp_seq),
             host3.rm(temp_seq)),
  Host = factor(rep(c("One", "Three", "Two"), each = length(temp_seq)))
) %>%
  left_join(species_to_host, by = "Host")

plot_data.host$Host <- factor(plot_data.host$Host, levels = c("One", "Two", "Three"))
fitted_data.host.dromedarii$Host <- factor(fitted_data.host.dromedarii$Host, levels = c("One", "Two", "Three"))
# Plot observed data and fitted curves with facet wrap by species
p.host.dromedarii <- ggplot() +
  geom_point(data = plot_data.host, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.dromedarii, aes(x = Temperature, y = Fitted, color = Host), size = 1) +
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
p.host.dromedarii

# Merge the observed and predicted data for comparison
merged_data.host.dromedarii <- fitted_data.host.dromedarii %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host"))

# Calculate RMSE for each species
mse_by_host.dromedarii <- merged_data.host.dromedarii %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.host.dromedarii <- merged_data.host.dromedarii %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_host.dromedarii,
     fitted_data.host.dromedarii,
     merged_data.host.dromedarii,
     host.dromedarii.fit.mcmc,
     host.dromedarii.fit,
     p.host.dromedarii,
     file="dromedarii.host.RData")

################################### By Host- Impeltatum ###################################
tick.rm <- tick.abr.new[tick.abr.new$Species != "Hyalomma impeltatum", ]
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
  mu.c ~ dexp(10)          # Mean for host effects c (vertex y)
  mu.Topt ~ dnorm(30, 0.1)          # Mean for host effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for host effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for host effects a (width)

  for (j in 1:n.host.rm) {
    mu.host.c[j] ~ dexp(mu.c)       # Host-specific vertex y
    mu.host.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Host-specific vertex x
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
    host.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.host.a[host.rm[i]] * (mu.host.Topt[host.rm[i]] - temp[i])^2 + mu.host.c[host.rm[i]]
    trait[i] ~ dnorm(mu[i], host.tau[host.rm[i]]) T(0, ) 
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
jag.data <- list(trait = tick.rm$Trait, N.obs.rm = N.obs.rm, 
                 temp = tick.rm$Temp, 
                 n.host.rm = n.host.rm, host.rm = host.rm)
# Fit model
host.impeltatum.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                            model.file = "host_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                            n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.impeltatum.fit.mcmc <- as.mcmc(host.impeltatum.fit) 
closeAllConnections()


# Get posterior means for parameters for each host group
summary_fit.host.rm <- summary(host.impeltatum.fit.mcmc)$statistics
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
fitted_data.host.impeltatum <- data.frame(Temperature = rep(temp_seq, n.host.rm), 
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
fitted_data.host.impeltatum <- data.frame(
  Temperature = rep(temp_seq, n.host.rm),
  Fitted = c(host1.rm(temp_seq),
             host2.rm(temp_seq),
             host3.rm(temp_seq)),
  Host = factor(rep(c("One", "Three", "Two"), each = length(temp_seq)))
) %>%
  left_join(species_to_host, by = "Host")

plot_data.host$Host <- factor(plot_data.host$Host, levels = c("One", "Two", "Three"))
fitted_data.host.impeltatum$Host <- factor(fitted_data.host.impeltatum$Host, levels = c("One", "Two", "Three"))
# Plot observed data and fitted curves with facet wrap by species
p.host.impeltatum <- ggplot() +
  geom_point(data = plot_data.host, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.impeltatum, aes(x = Temperature, y = Fitted, color = Host), size = 1) +
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
p.host.impeltatum

# Merge the observed and predicted data for comparison
merged_data.host.impeltatum <- fitted_data.host.impeltatum %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host"))

# Calculate RMSE for each species
mse_by_host.impeltatum <- merged_data.host.impeltatum %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.host.impeltatum <- merged_data.host.impeltatum %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_host.impeltatum,
     fitted_data.host.impeltatum,
     merged_data.host.impeltatum,
     host.impeltatum.fit.mcmc,
     host.impeltatum.fit,
     p.host.impeltatum,
     file="impeltatum.host.RData")

################################### By Host- Lusitanicum ###################################
tick.rm <- tick.abr.new[tick.abr.new$Species != "Hyalomma lusitanicum", ]
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
  mu.c ~ dexp(10)          # Mean for host effects c (vertex y)
  mu.Topt ~ dnorm(30, 0.1)          # Mean for host effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for host effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for host effects a (width)

  for (j in 1:n.host.rm) {
    mu.host.c[j] ~ dexp(mu.c)       # Host-specific vertex y
    mu.host.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Host-specific vertex x
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
    host.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.host.a[host.rm[i]] * (mu.host.Topt[host.rm[i]] - temp[i])^2 + mu.host.c[host.rm[i]]
    trait[i] ~ dnorm(mu[i], host.tau[host.rm[i]]) T(0, ) 
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
jag.data <- list(trait = tick.rm$Trait, N.obs.rm = N.obs.rm, 
                 temp = tick.rm$Temp, 
                 n.host.rm = n.host.rm, host.rm = host.rm)
# Fit model
host.lusitanicum.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                             model.file = "host_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                             n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.lusitanicum.fit.mcmc <- as.mcmc(host.lusitanicum.fit) 
closeAllConnections()


# Get posterior means for parameters for each host group
summary_fit.host.rm <- summary(host.lusitanicum.fit.mcmc)$statistics
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
fitted_data.host.lusitanicum <- data.frame(Temperature = rep(temp_seq, n.host.rm), 
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
fitted_data.host.lusitanicum <- data.frame(
  Temperature = rep(temp_seq, n.host.rm),
  Fitted = c(host1.rm(temp_seq),
             host2.rm(temp_seq),
             host3.rm(temp_seq)),
  Host = factor(rep(c("One", "Three", "Two"), each = length(temp_seq)))
) %>%
  left_join(species_to_host, by = "Host")

plot_data.host$Host <- factor(plot_data.host$Host, levels = c("One", "Two", "Three"))
fitted_data.host.lusitanicum$Host <- factor(fitted_data.host.lusitanicum$Host, levels = c("One", "Two", "Three"))
# Plot observed data and fitted curves with facet wrap by species
p.host.lusitanicum <- ggplot() +
  geom_point(data = plot_data.host, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.lusitanicum, aes(x = Temperature, y = Fitted, color = Host), size = 1) +
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
p.host.lusitanicum

# Merge the observed and predicted data for comparison
merged_data.host.lusitanicum <- fitted_data.host.lusitanicum %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host"))

# Calculate RMSE for each species
mse_by_host.lusitanicum <- merged_data.host.lusitanicum %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.host.lusitanicum <- merged_data.host.lusitanicum %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_host.lusitanicum,
     fitted_data.host.lusitanicum,
     merged_data.host.lusitanicum,
     host.lusitanicum.fit.mcmc,
     host.lusitanicum.fit,
     p.host.lusitanicum,
     file="lusitanicum.host.RData")

################################### By Host- Schulzei ###################################
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
  mu.c ~ dexp(10)          # Mean for host effects c (vertex y)
  mu.Topt ~ dnorm(30, 0.1)          # Mean for host effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for host effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for host effects a (width)

  for (j in 1:n.host.rm) {
    mu.host.c[j] ~ dexp(mu.c)       # Host-specific vertex y
    mu.host.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Host-specific vertex x
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
    host.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.host.a[host.rm[i]] * (mu.host.Topt[host.rm[i]] - temp[i])^2 + mu.host.c[host.rm[i]]
    trait[i] ~ dnorm(mu[i], host.tau[host.rm[i]]) T(0, ) 
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

plot_data.host$Host <- factor(plot_data.host$Host, levels = c("One", "Two", "Three"))
fitted_data.host.schulzei$Host <- factor(fitted_data.host.schulzei$Host, levels = c("One", "Two", "Three"))
# Plot observed data and fitted curves with facet wrap by species
p.host.schulzei <- ggplot() +
  geom_point(data = plot_data.host, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.schulzei, aes(x = Temperature, y = Fitted, color = Host), size = 1) +
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
p.host.schulzei

# Merge the observed and predicted data for comparison
merged_data.host.schulzei <- fitted_data.host.schulzei %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host"))

# Calculate RMSE for each species
mse_by_host.schulzei <- merged_data.host.schulzei %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.host.schulzei <- merged_data.host.schulzei %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_host.schulzei,
     fitted_data.host.schulzei,
     merged_data.host.schulzei,
     host.schulzei.fit.mcmc,
     host.schulzei.fit,
     p.host.schulzei,
     file="schulzei.host.RData")

################################### By Host and Climate- Leporispalustris ###################################
tick.rm <- tick.abr.new[tick.abr.new$Species != "Haemaphysalis leporispalustris", ]
host.rm <- as.factor(tick.rm$Host)
n.host.rm <- length(unique(tick.rm$Host))
N.obs.rm <- length(tick.rm$Trait)
climate.rm <- as.factor(tick.rm$Climate)
n.climate.rm <- length(unique(tick.rm$Climate))
N.obs.rm <- length(tick.rm$Trait)
unique(tick.rm$Species)

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

  for (j in 1:n.host.rm) {
    mu.host.c[j] ~ dexp(mu.c)       # Host-specific vertex y
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
  }
  for (k in 1:n.climate.rm) {
    mu.climate.Topt[k] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific vertex x
    climate.tau[k] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.host.a[host.rm[i]] * (mu.climate.Topt[climate.rm[i]] - temp[i])^2 + mu.host.c[host.rm[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[climate.rm[i]]) T(0, )  
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
jag.data <- list(trait = tick.rm$Trait, N.obs.rm = N.obs.rm, 
                 temp = tick.rm$Temp, 
                 n.host.rm = n.host.rm, host.rm = host.rm, 
                 n.climate.rm = n.climate.rm, climate.rm = climate.rm)
# Fit model
host.climate.leporispalustris.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                                          model.file = "host_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                                          n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.climate.leporispalustris.fit.mcmc <- as.mcmc(host.climate.leporispalustris.fit) 
closeAllConnections()

# Get posterior means for parameters for each host group
summary_fit.host.climate <- summary(host.climate.leporispalustris.fit.mcmc)$statistics
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
fitted_data.host.climate.leporispalustris <- data.frame(Temperature = rep(temp_seq, n.host * n.climate), 
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
fitted_data.host.climate.leporispalustris <- data.frame(Temperature = rep(temp_seq, n.host * n.climate), 
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
fitted_data.host.climate.leporispalustris$Host <- factor(fitted_data.host.climate.leporispalustris$Host, levels = c("One", "Two", "Three"))

# Combine Climate and Size for coloring
fitted_data.host.climate.leporispalustris$`Host, Climate` <- with(fitted_data.host.climate.leporispalustris, paste(Host, Climate, sep = ", "))
fitted_data.host.climate.leporispalustris$`Host, Climate` <- factor(fitted_data.host.climate.leporispalustris$`Host, Climate`, levels = c("One, Tropical", 
                                                                                                                      "Two, Subtropical", "Three, Mixed",
                                                                                                                      "Three, Subtropical", "Three, Temperate",
                                                                                                                      "Three, Tropical"))

# Plot observed data and fitted curves with facet wrap by species
p.host.climate.leporispalustris <- ggplot() +
  geom_point(data = plot_data.host.climate, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.climate.leporispalustris, aes(x = Temperature, y = Fitted, color = `Host, Climate`), size = 1) +
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
p.host.climate.leporispalustris

# Merge the observed and predicted data for comparison
merged_data.host.climate.leporispalustris <- fitted_data.host.climate.leporispalustris %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host", "Climate" = "Climate"))

# Calculate RMSE for each species
mse_by_host.climate.leporispalustris <- merged_data.host.climate.leporispalustris %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.host.climate.leporispalustris <- merged_data.host.climate.leporispalustris %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 


save(mse_by_host.climate.leporispalustris,
     fitted_data.host.climate.leporispalustris,
     merged_data.host.climate.leporispalustris,
     host.climate.leporispalustris.fit.mcmc,
     host.climate.leporispalustris.fit,
     p.host.climate.leporispalustris,
     file="host.climate.leporispalustris.fit.RData")

################################### By Host and Climate- Aegyptium ###################################
tick.rm <- tick.abr.new[tick.abr.new$Species != "Hyalomma aegyptium", ]
host.rm <- as.factor(tick.rm$Host)
n.host.rm <- length(unique(tick.rm$Host))
N.obs.rm <- length(tick.rm$Trait)
climate.rm <- as.factor(tick.rm$Climate)
n.climate.rm <- length(unique(tick.rm$Climate))
N.obs.rm <- length(tick.rm$Trait)
unique(tick.rm$Species)

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

  for (j in 1:n.host.rm) {
    mu.host.c[j] ~ dexp(mu.c)       # Host-specific vertex y
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
  }
  for (k in 1:n.climate.rm) {
    mu.climate.Topt[k] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific vertex x
    climate.tau[k] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.host.a[host.rm[i]] * (mu.climate.Topt[climate.rm[i]] - temp[i])^2 + mu.host.c[host.rm[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[climate.rm[i]]) T(0, )  
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
jag.data <- list(trait = tick.rm$Trait, N.obs.rm = N.obs.rm, 
                 temp = tick.rm$Temp, 
                 n.host.rm = n.host.rm, host.rm = host.rm, 
                 n.climate.rm = n.climate.rm, climate.rm = climate.rm)
# Fit model
host.climate.aegyptium.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                                   model.file = "host_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                                   n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.climate.aegyptium.fit.mcmc <- as.mcmc(host.climate.aegyptium.fit) 
closeAllConnections()

# Get posterior means for parameters for each host group
summary_fit.host.climate <- summary(host.climate.aegyptium.fit.mcmc)$statistics
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
fitted_data.host.climate.aegyptium <- data.frame(Temperature = rep(temp_seq, n.host * n.climate), 
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
fitted_data.host.climate.aegyptium <- data.frame(Temperature = rep(temp_seq, n.host * n.climate), 
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
fitted_data.host.climate.aegyptium$Host <- factor(fitted_data.host.climate.aegyptium$Host, levels = c("One", "Two", "Three"))

# Combine Climate and Size for coloring
fitted_data.host.climate.aegyptium$`Host, Climate` <- with(fitted_data.host.climate.aegyptium, paste(Host, Climate, sep = ", "))
fitted_data.host.climate.aegyptium$`Host, Climate` <- factor(fitted_data.host.climate.aegyptium$`Host, Climate`, levels = c("One, Tropical", 
                                                                                                        "Two, Subtropical", "Three, Mixed",
                                                                                                        "Three, Subtropical", "Three, Temperate",
                                                                                                        "Three, Tropical"))

# Plot observed data and fitted curves with facet wrap by species
p.host.climate.aegyptium <- ggplot() +
  geom_point(data = plot_data.host.climate, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.climate.aegyptium, aes(x = Temperature, y = Fitted, color = `Host, Climate`), size = 1) +
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
p.host.climate.aegyptium

# Merge the observed and predicted data for comparison
merged_data.host.climate.aegyptium <- fitted_data.host.climate.aegyptium %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host", "Climate" = "Climate"))

# Calculate RMSE for each species
mse_by_host.climate.aegyptium <- merged_data.host.climate.aegyptium %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.host.climate.aegyptium <- merged_data.host.climate.aegyptium %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 


save(mse_by_host.climate.aegyptium,
     fitted_data.host.climate.aegyptium,
     merged_data.host.climate.aegyptium,
     host.climate.aegyptium.fit.mcmc,
     host.climate.aegyptium.fit,
     p.host.climate.aegyptium,
     file="host.climate.aegyptium.fit.RData")

################################### By Host and Climate- Dromedarii ###################################
tick.rm <- tick.abr.new[tick.abr.new$Species != "Hyalomma dromedarii", ]
host.rm <- as.factor(tick.rm$Host)
n.host.rm <- length(unique(tick.rm$Host))
N.obs.rm <- length(tick.rm$Trait)
climate.rm <- as.factor(tick.rm$Climate)
n.climate.rm <- length(unique(tick.rm$Climate))
N.obs.rm <- length(tick.rm$Trait)
unique(tick.rm$Species)

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

  for (j in 1:n.host.rm) {
    mu.host.c[j] ~ dexp(mu.c)       # Host-specific vertex y
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
  }
  for (k in 1:n.climate.rm) {
    mu.climate.Topt[k] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific vertex x
    climate.tau[k] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.host.a[host.rm[i]] * (mu.climate.Topt[climate.rm[i]] - temp[i])^2 + mu.host.c[host.rm[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[climate.rm[i]]) T(0, )  
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
jag.data <- list(trait = tick.rm$Trait, N.obs.rm = N.obs.rm, 
                 temp = tick.rm$Temp, 
                 n.host.rm = n.host.rm, host.rm = host.rm, 
                 n.climate.rm = n.climate.rm, climate.rm = climate.rm)
# Fit model
host.climate.dromedarii.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                                    model.file = "host_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                                    n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.climate.dromedarii.fit.mcmc <- as.mcmc(host.climate.dromedarii.fit) 
closeAllConnections()

# Get posterior means for parameters for each host group
summary_fit.host.climate <- summary(host.climate.dromedarii.fit.mcmc)$statistics
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
fitted_data.host.climate.dromedarii <- data.frame(Temperature = rep(temp_seq, n.host * n.climate), 
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
fitted_data.host.climate.dromedarii <- data.frame(Temperature = rep(temp_seq, n.host * n.climate), 
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
fitted_data.host.climate.dromedarii$Host <- factor(fitted_data.host.climate.dromedarii$Host, levels = c("One", "Two", "Three"))

# Combine Climate and Size for coloring
fitted_data.host.climate.dromedarii$`Host, Climate` <- with(fitted_data.host.climate.dromedarii, paste(Host, Climate, sep = ", "))
fitted_data.host.climate.dromedarii$`Host, Climate` <- factor(fitted_data.host.climate.dromedarii$`Host, Climate`, levels = c("One, Tropical", 
                                                                                                          "Two, Subtropical", "Three, Mixed",
                                                                                                          "Three, Subtropical", "Three, Temperate",
                                                                                                          "Three, Tropical"))

# Plot observed data and fitted curves with facet wrap by species
p.host.climate.dromedarii <- ggplot() +
  geom_point(data = plot_data.host.climate, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.climate.dromedarii, aes(x = Temperature, y = Fitted, color = `Host, Climate`), size = 1) +
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
p.host.climate.dromedarii

# Merge the observed and predicted data for comparison
merged_data.host.climate.dromedarii <- fitted_data.host.climate.dromedarii %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host", "Climate" = "Climate"))

# Calculate RMSE for each species
mse_by_host.climate.dromedarii <- merged_data.host.climate.dromedarii %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.host.climate.dromedarii <- merged_data.host.climate.dromedarii %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 


save(mse_by_host.climate.dromedarii,
     fitted_data.host.climate.dromedarii,
     merged_data.host.climate.dromedarii,
     host.climate.dromedarii.fit.mcmc,
     host.climate.dromedarii.fit,
     p.host.climate.dromedarii,
     file="host.climate.dromedarii.fit.RData")

################################### By Host and Climate- Impeltatum ###################################
tick.rm <- tick.abr.new[tick.abr.new$Species != "Hyalomma impeltatum", ]
host.rm <- as.factor(tick.rm$Host)
n.host.rm <- length(unique(tick.rm$Host))
N.obs.rm <- length(tick.rm$Trait)
climate.rm <- as.factor(tick.rm$Climate)
n.climate.rm <- length(unique(tick.rm$Climate))
N.obs.rm <- length(tick.rm$Trait)
unique(tick.rm$Species)

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

  for (j in 1:n.host.rm) {
    mu.host.c[j] ~ dexp(mu.c)       # Host-specific vertex y
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
  }
  for (k in 1:n.climate.rm) {
    mu.climate.Topt[k] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific vertex x
    climate.tau[k] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.host.a[host.rm[i]] * (mu.climate.Topt[climate.rm[i]] - temp[i])^2 + mu.host.c[host.rm[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[climate.rm[i]]) T(0, )  
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
jag.data <- list(trait = tick.rm$Trait, N.obs.rm = N.obs.rm, 
                 temp = tick.rm$Temp, 
                 n.host.rm = n.host.rm, host.rm = host.rm, 
                 n.climate.rm = n.climate.rm, climate.rm = climate.rm)
# Fit model
host.climate.impeltatum.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                                    model.file = "host_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                                    n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.climate.impeltatum.fit.mcmc <- as.mcmc(host.climate.impeltatum.fit) 
closeAllConnections()

# Get posterior means for parameters for each host group
summary_fit.host.climate <- summary(host.climate.impeltatum.fit.mcmc)$statistics
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
fitted_data.host.climate.impeltatum <- data.frame(Temperature = rep(temp_seq, n.host * n.climate), 
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
fitted_data.host.climate.impeltatum <- data.frame(Temperature = rep(temp_seq, n.host * n.climate), 
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
fitted_data.host.climate.impeltatum$Host <- factor(fitted_data.host.climate.impeltatum$Host, levels = c("One", "Two", "Three"))

# Combine Climate and Size for coloring
fitted_data.host.climate.impeltatum$`Host, Climate` <- with(fitted_data.host.climate.impeltatum, paste(Host, Climate, sep = ", "))
fitted_data.host.climate.impeltatum$`Host, Climate` <- factor(fitted_data.host.climate.impeltatum$`Host, Climate`, levels = c("One, Tropical", 
                                                                                                          "Two, Subtropical", "Three, Mixed",
                                                                                                          "Three, Subtropical", "Three, Temperate",
                                                                                                          "Three, Tropical"))

# Plot observed data and fitted curves with facet wrap by species
p.host.climate.impeltatum <- ggplot() +
  geom_point(data = plot_data.host.climate, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.climate.impeltatum, aes(x = Temperature, y = Fitted, color = `Host, Climate`), size = 1) +
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
p.host.climate.impeltatum

# Merge the observed and predicted data for comparison
merged_data.host.climate.impeltatum <- fitted_data.host.climate.impeltatum %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host", "Climate" = "Climate"))

# Calculate RMSE for each species
mse_by_host.climate.impeltatum <- merged_data.host.climate.impeltatum %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.host.climate.impeltatum <- merged_data.host.climate.impeltatum %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 


save(mse_by_host.climate.impeltatum,
     fitted_data.host.climate.impeltatum,
     merged_data.host.climate.impeltatum,
     host.climate.impeltatum.fit.mcmc,
     host.climate.impeltatum.fit,
     p.host.climate.impeltatum,
     file="host.climate.impeltatum.fit.RData")

################################### By Host and Climate- Lusitanicum ###################################
tick.rm <- tick.abr.new[tick.abr.new$Species != "Hyalomma lusitanicum", ]
host.rm <- as.factor(tick.rm$Host)
n.host.rm <- length(unique(tick.rm$Host))
N.obs.rm <- length(tick.rm$Trait)
climate.rm <- as.factor(tick.rm$Climate)
n.climate.rm <- length(unique(tick.rm$Climate))
N.obs.rm <- length(tick.rm$Trait)
unique(tick.rm$Species)

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

  for (j in 1:n.host.rm) {
    mu.host.c[j] ~ dexp(mu.c)       # Host-specific vertex y
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
  }
  for (k in 1:n.climate.rm) {
    mu.climate.Topt[k] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific vertex x
    climate.tau[k] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.host.a[host.rm[i]] * (mu.climate.Topt[climate.rm[i]] - temp[i])^2 + mu.host.c[host.rm[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[climate.rm[i]]) T(0, )  
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
jag.data <- list(trait = tick.rm$Trait, N.obs.rm = N.obs.rm, 
                 temp = tick.rm$Temp, 
                 n.host.rm = n.host.rm, host.rm = host.rm, 
                 n.climate.rm = n.climate.rm, climate.rm = climate.rm)
# Fit model
host.climate.lusitanicum.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                                     model.file = "host_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                                     n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.climate.lusitanicum.fit.mcmc <- as.mcmc(host.climate.lusitanicum.fit) 
closeAllConnections()

# Get posterior means for parameters for each host group
summary_fit.host.climate <- summary(host.climate.lusitanicum.fit.mcmc)$statistics
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
fitted_data.host.climate.lusitanicum <- data.frame(Temperature = rep(temp_seq, n.host * n.climate), 
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
fitted_data.host.climate.lusitanicum <- data.frame(Temperature = rep(temp_seq, n.host * n.climate), 
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
fitted_data.host.climate.lusitanicum$Host <- factor(fitted_data.host.climate.lusitanicum$Host, levels = c("One", "Two", "Three"))

# Combine Climate and Size for coloring
fitted_data.host.climate.lusitanicum$`Host, Climate` <- with(fitted_data.host.climate.lusitanicum, paste(Host, Climate, sep = ", "))
fitted_data.host.climate.lusitanicum$`Host, Climate` <- factor(fitted_data.host.climate.lusitanicum$`Host, Climate`, levels = c("One, Tropical", 
                                                                                                            "Two, Subtropical", "Three, Mixed",
                                                                                                            "Three, Subtropical", "Three, Temperate",
                                                                                                            "Three, Tropical"))

# Plot observed data and fitted curves with facet wrap by species
p.host.climate.lusitanicum <- ggplot() +
  geom_point(data = plot_data.host.climate, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.climate.lusitanicum, aes(x = Temperature, y = Fitted, color = `Host, Climate`), size = 1) +
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
p.host.climate.lusitanicum

# Merge the observed and predicted data for comparison
merged_data.host.climate.lusitanicum <- fitted_data.host.climate.lusitanicum %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host", "Climate" = "Climate"))

# Calculate RMSE for each species
mse_by_host.climate.lusitanicum <- merged_data.host.climate.lusitanicum %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.host.climate.lusitanicum <- merged_data.host.climate.lusitanicum %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 


save(mse_by_host.climate.lusitanicum,
     fitted_data.host.climate.lusitanicum,
     merged_data.host.climate.lusitanicum,
     host.climate.lusitanicum.fit.mcmc,
     host.climate.lusitanicum.fit,
     p.host.climate.lusitanicum,
     file="host.climate.lusitanicum.fit.RData")

################################### By Host and Climate- Schulzei ###################################
tick.rm <- tick.abr.new[tick.abr.new$Species != "Hyalomma schulzei", ]
host.rm <- as.factor(tick.rm$Host)
n.host.rm <- length(unique(tick.rm$Host))
N.obs.rm <- length(tick.rm$Trait)
climate.rm <- as.factor(tick.rm$Climate)
n.climate.rm <- length(unique(tick.rm$Climate))
N.obs.rm <- length(tick.rm$Trait)
unique(tick.rm$Species)

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

  for (j in 1:n.host.rm) {
    mu.host.c[j] ~ dexp(mu.c)       # Host-specific vertex y
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
  }
  for (k in 1:n.climate.rm) {
    mu.climate.Topt[k] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific vertex x
    climate.tau[k] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.host.a[host.rm[i]] * (mu.climate.Topt[climate.rm[i]] - temp[i])^2 + mu.host.c[host.rm[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[climate.rm[i]]) T(0, )  
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
jag.data <- list(trait = tick.rm$Trait, N.obs.rm = N.obs.rm, 
                 temp = tick.rm$Temp, 
                 n.host.rm = n.host.rm, host.rm = host.rm, 
                 n.climate.rm = n.climate.rm, climate.rm = climate.rm)
# Fit model
host.climate.schulzei.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                                  model.file = "host_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                                  n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.climate.schulzei.fit.mcmc <- as.mcmc(host.climate.schulzei.fit) 
closeAllConnections()

# Get posterior means for parameters for each host group
summary_fit.host.climate <- summary(host.climate.schulzei.fit.mcmc)$statistics
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
fitted_data.host.climate.schulzei <- data.frame(Temperature = rep(temp_seq, n.host * n.climate), 
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
fitted_data.host.climate.schulzei <- data.frame(Temperature = rep(temp_seq, n.host * n.climate), 
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
fitted_data.host.climate.schulzei$Host <- factor(fitted_data.host.climate.schulzei$Host, levels = c("One", "Two", "Three"))

# Combine Climate and Size for coloring
fitted_data.host.climate.schulzei$`Host, Climate` <- with(fitted_data.host.climate.schulzei, paste(Host, Climate, sep = ", "))
fitted_data.host.climate.schulzei$`Host, Climate` <- factor(fitted_data.host.climate.schulzei$`Host, Climate`, levels = c("One, Tropical", 
                                                                                                      "Two, Subtropical", "Three, Mixed",
                                                                                                      "Three, Subtropical", "Three, Temperate",
                                                                                                      "Three, Tropical"))

# Plot observed data and fitted curves with facet wrap by species
p.host.climate.schulzei <- ggplot() +
  geom_point(data = plot_data.host.climate, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.climate.schulzei, aes(x = Temperature, y = Fitted, color = `Host, Climate`), size = 1) +
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
p.host.climate.schulzei

# Merge the observed and predicted data for comparison
merged_data.host.climate.schulzei <- fitted_data.host.climate.schulzei %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host", "Climate" = "Climate"))

# Calculate RMSE for each species
mse_by_host.climate.schulzei <- merged_data.host.climate.schulzei %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.host.climate.schulzei <- merged_data.host.climate.schulzei %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 


save(mse_by_host.climate.schulzei,
     fitted_data.host.climate.schulzei,
     merged_data.host.climate.schulzei,
     host.climate.schulzei.fit.mcmc,
     host.climate.schulzei.fit,
     p.host.climate.schulzei,
     file="host.climate.schulzei.fit.RData")

################################### By Host, Genus, and Climate- Aegyptium ###################################
# Model by typical number of host species
tick.rm <- tick.abr.new[tick.abr.new$Species != "Hyalomma aegyptium", ]
climate.rm <- as.factor(tick.rm$Climate)
n.climate.rm <- length(unique(tick.rm$Climate))
genus.rm <- as.factor(tick.rm$Genus)
n.genus.rm <- length(unique(tick.rm$Genus))
host.rm <- as.factor(tick.rm$Host)
n.host.rm <- length(unique(tick.rm$Host))
N.obs.rm <- length(tick.rm$Trait)
unique(tick.rm$Species)
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

  for (j in 1:n.host.rm) {
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
  }
  for (k in 1:n.climate.rm) {
    mu.climate.Topt[k] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific vertex x
    climate.tau[k] ~ dexp(tau1)
  }
  for (s in 1:n.genus.rm) {
    mu.genus.c[s] ~ dexp(mu.c)       # genus-specific vertex y
  }
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.host.a[host.rm[i]] * (mu.climate.Topt[climate.rm[i]] - temp[i])^2 + mu.genus.c[genus.rm[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[climate.rm[i]]) T(0, ) 
    log_lik[i] <- logdensity.norm(trait[i], mu[i], climate.tau[climate.rm[i]])
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
jag.data <- list(trait = tick.rm$Trait, N.obs.rm = N.obs.rm, 
                 temp = tick.rm$Temp, 
                 n.host.rm = n.host.rm, host.rm = host.rm, 
                 n.climate.rm = n.climate.rm, climate.rm = climate.rm,
                 n.genus.rm = n.genus.rm, genus.rm = genus.rm)
# Fit model
h.c.g.aegyptium.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                            model.file = "host.climate.genus_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                            n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
h.c.g.aegyptium.fit.mcmc <- as.mcmc(h.c.g.aegyptium.fit) 
closeAllConnections()

# Get posterior means for parameters for each host group
summary_fit.h.c.g <- summary(h.c.g.aegyptium.fit.mcmc)$statistics
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
fitted_data.h.c.g.aegyptium <- do.call(rbind, fitted_list)

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
fitted_data.h.c.g.aegyptium <- fitted_data.h.c.g.aegyptium |>  
  left_join(species_to_h.c.g, by = c("Host", "Climate", "Genus"))|> 
  filter(!is.na(Species)) 

# Order factor levels
plot_data.h.c.g$Host <- factor(plot_data.h.c.g$Host, levels = c("One", "Two", "Three"))
fitted_data.h.c.g.aegyptium$Host <- factor(fitted_data.h.c.g.aegyptium$Host, levels = c("One", "Two", "Three"))

# Combine Host, Climate, and Size for coloring
fitted_data.h.c.g.aegyptium$`Host, Climate, Genus` <- with(fitted_data.h.c.g.aegyptium, paste(Host, Climate, Genus, sep = ", "))

# Plot observed data and fitted curves with facet wrap by species
p.h.c.g.aegyptium <- ggplot() +
  geom_point(data = plot_data.h.c.g, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.h.c.g.aegyptium, aes(x = Temperature, y = Fitted, color = `Host, Climate, Genus`), size = 1) +
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
p.h.c.g.aegyptium
ggsave("tick.h.c.g.aegyptium.plot.png", plot = p.h.c.g.aegyptium, width = 8.5, height = 5)

# Merge the observed and predicted data for comparison
merged_data.h.c.g.aegyptium <- fitted_data.h.c.g.aegyptium %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", 
                         "Host" = "Host", "Climate" = "Climate",
                         "Genus" = "Genus"))

# Calculate RMSE for each species
mse_by_h.c.g.aegyptium <- merged_data.h.c.g.aegyptium %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.h.c.g.aegyptium <- merged_data.h.c.g.aegyptium %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 


save(mse_by_h.c.g.aegyptium,
     fitted_data.h.c.g.aegyptium,
     merged_data.h.c.g.aegyptium,
     h.c.g.aegyptium.fit.mcmc,
     h.c.g.aegyptium.fit,
     p.h.c.g.aegyptium,
     file="host.climate.genus.aegyptium.fit.RData")

################################### By Host, Genus, and Climate- Dromedarii ###################################
# Model by typical number of host species
tick.rm <- tick.abr.new[tick.abr.new$Species != "Hyalomma dromedarii", ]
climate.rm <- as.factor(tick.rm$Climate)
n.climate.rm <- length(unique(tick.rm$Climate))
genus.rm <- as.factor(tick.rm$Genus)
n.genus.rm <- length(unique(tick.rm$Genus))
host.rm <- as.factor(tick.rm$Host)
n.host.rm <- length(unique(tick.rm$Host))
N.obs.rm <- length(tick.rm$Trait)
unique(tick.rm$Species)
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

  for (j in 1:n.host.rm) {
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
  }
  for (k in 1:n.climate.rm) {
    mu.climate.Topt[k] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific vertex x
    climate.tau[k] ~ dexp(tau1)
  }
  for (s in 1:n.genus.rm) {
    mu.genus.c[s] ~ dexp(mu.c)       # genus-specific vertex y
  }
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.host.a[host.rm[i]] * (mu.climate.Topt[climate.rm[i]] - temp[i])^2 + mu.genus.c[genus.rm[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[climate.rm[i]]) T(0, ) 
    log_lik[i] <- logdensity.norm(trait[i], mu[i], climate.tau[climate.rm[i]])
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
jag.data <- list(trait = tick.rm$Trait, N.obs.rm = N.obs.rm, 
                 temp = tick.rm$Temp, 
                 n.host.rm = n.host.rm, host.rm = host.rm, 
                 n.climate.rm = n.climate.rm, climate.rm = climate.rm,
                 n.genus.rm = n.genus.rm, genus.rm = genus.rm)
# Fit model
h.c.g.dromedarii.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                             model.file = "host.climate.genus_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                             n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
h.c.g.dromedarii.fit.mcmc <- as.mcmc(h.c.g.dromedarii.fit) 
closeAllConnections()

# Get posterior means for parameters for each host group
summary_fit.h.c.g <- summary(h.c.g.dromedarii.fit.mcmc)$statistics
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
fitted_data.h.c.g.dromedarii <- do.call(rbind, fitted_list)

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
fitted_data.h.c.g.dromedarii <- fitted_data.h.c.g.dromedarii |>  
  left_join(species_to_h.c.g, by = c("Host", "Climate", "Genus"))|> 
  filter(!is.na(Species)) 

# Order factor levels
plot_data.h.c.g$Host <- factor(plot_data.h.c.g$Host, levels = c("One", "Two", "Three"))
fitted_data.h.c.g.dromedarii$Host <- factor(fitted_data.h.c.g.dromedarii$Host, levels = c("One", "Two", "Three"))

# Combine Host, Climate, and Size for coloring
fitted_data.h.c.g.dromedarii$`Host, Climate, Genus` <- with(fitted_data.h.c.g.dromedarii, paste(Host, Climate, Genus, sep = ", "))

# Plot observed data and fitted curves with facet wrap by species
p.h.c.g.dromedarii <- ggplot() +
  geom_point(data = plot_data.h.c.g, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.h.c.g.dromedarii, aes(x = Temperature, y = Fitted, color = `Host, Climate, Genus`), size = 1) +
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
p.h.c.g.dromedarii
ggsave("tick.h.c.g.dromedarii.plot.png", plot = p.h.c.g.dromedarii, width = 8.5, height = 5)

# Merge the observed and predicted data for comparison
merged_data.h.c.g.dromedarii <- fitted_data.h.c.g.dromedarii %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", 
                         "Host" = "Host", "Climate" = "Climate",
                         "Genus" = "Genus"))

# Calculate RMSE for each species
mse_by_h.c.g.dromedarii <- merged_data.h.c.g.dromedarii %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.h.c.g.dromedarii <- merged_data.h.c.g.dromedarii %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 


save(mse_by_h.c.g.dromedarii,
     fitted_data.h.c.g.dromedarii,
     merged_data.h.c.g.dromedarii,
     h.c.g.dromedarii.fit.mcmc,
     h.c.g.dromedarii.fit,
     p.h.c.g.dromedarii,
     file="host.climate.genus.dromedarii.fit.RData")

################################### By Host, Genus, and Climate- Lusitanicum ###################################
# Model by typical number of host species
tick.rm <- tick.abr.new[tick.abr.new$Species != "Hyalomma lusitanicum", ]
climate.rm <- as.factor(tick.rm$Climate)
n.climate.rm <- length(unique(tick.rm$Climate))
genus.rm <- as.factor(tick.rm$Genus)
n.genus.rm <- length(unique(tick.rm$Genus))
host.rm <- as.factor(tick.rm$Host)
n.host.rm <- length(unique(tick.rm$Host))
N.obs.rm <- length(tick.rm$Trait)
unique(tick.rm$Species)
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

  for (j in 1:n.host.rm) {
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
  }
  for (k in 1:n.climate.rm) {
    mu.climate.Topt[k] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific vertex x
    climate.tau[k] ~ dexp(tau1)
  }
  for (s in 1:n.genus.rm) {
    mu.genus.c[s] ~ dexp(mu.c)       # genus-specific vertex y
  }
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.host.a[host.rm[i]] * (mu.climate.Topt[climate.rm[i]] - temp[i])^2 + mu.genus.c[genus.rm[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[climate.rm[i]]) T(0, ) 
    log_lik[i] <- logdensity.norm(trait[i], mu[i], climate.tau[climate.rm[i]])
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
jag.data <- list(trait = tick.rm$Trait, N.obs.rm = N.obs.rm, 
                 temp = tick.rm$Temp, 
                 n.host.rm = n.host.rm, host.rm = host.rm, 
                 n.climate.rm = n.climate.rm, climate.rm = climate.rm,
                 n.genus.rm = n.genus.rm, genus.rm = genus.rm)
# Fit model
h.c.g.lusitanicum.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                              model.file = "host.climate.genus_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                              n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
h.c.g.lusitanicum.fit.mcmc <- as.mcmc(h.c.g.lusitanicum.fit) 
closeAllConnections()

# Get posterior means for parameters for each host group
summary_fit.h.c.g <- summary(h.c.g.lusitanicum.fit.mcmc)$statistics
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
fitted_data.h.c.g.lusitanicum <- do.call(rbind, fitted_list)

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
fitted_data.h.c.g.lusitanicum <- fitted_data.h.c.g.lusitanicum |>  
  left_join(species_to_h.c.g, by = c("Host", "Climate", "Genus"))|> 
  filter(!is.na(Species)) 

# Order factor levels
plot_data.h.c.g$Host <- factor(plot_data.h.c.g$Host, levels = c("One", "Two", "Three"))
fitted_data.h.c.g.lusitanicum$Host <- factor(fitted_data.h.c.g.lusitanicum$Host, levels = c("One", "Two", "Three"))

# Combine Host, Climate, and Size for coloring
fitted_data.h.c.g.lusitanicum$`Host, Climate, Genus` <- with(fitted_data.h.c.g.lusitanicum, paste(Host, Climate, Genus, sep = ", "))

# Plot observed data and fitted curves with facet wrap by species
p.h.c.g.lusitanicum <- ggplot() +
  geom_point(data = plot_data.h.c.g, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.h.c.g.lusitanicum, aes(x = Temperature, y = Fitted, color = `Host, Climate, Genus`), size = 1) +
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
p.h.c.g.lusitanicum
ggsave("tick.h.c.g.lusitanicum.plot.png", plot = p.h.c.g.lusitanicum, width = 8.5, height = 5)

# Merge the observed and predicted data for comparison
merged_data.h.c.g.lusitanicum <- fitted_data.h.c.g.lusitanicum %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", 
                         "Host" = "Host", "Climate" = "Climate",
                         "Genus" = "Genus"))

# Calculate RMSE for each species
mse_by_h.c.g.lusitanicum <- merged_data.h.c.g.lusitanicum %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.h.c.g.lusitanicum <- merged_data.h.c.g.lusitanicum %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 


save(mse_by_h.c.g.lusitanicum,
     fitted_data.h.c.g.lusitanicum,
     merged_data.h.c.g.lusitanicum,
     h.c.g.lusitanicum.fit.mcmc,
     h.c.g.lusitanicum.fit,
     p.h.c.g.lusitanicum,
     file="host.climate.genus.lusitanicum.fit.RData")

################################### By Host, Genus, and Climate- Schulzei ###################################
# Model by typical number of host species
tick.rm <- tick.abr.new[tick.abr.new$Species != "Hyalomma schulzei", ]
climate.rm <- as.factor(tick.rm$Climate)
n.climate.rm <- length(unique(tick.rm$Climate))
genus.rm <- as.factor(tick.rm$Genus)
n.genus.rm <- length(unique(tick.rm$Genus))
host.rm <- as.factor(tick.rm$Host)
n.host.rm <- length(unique(tick.rm$Host))
N.obs.rm <- length(tick.rm$Trait)
unique(tick.rm$Species)
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

  for (j in 1:n.host.rm) {
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
  }
  for (k in 1:n.climate.rm) {
    mu.climate.Topt[k] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific vertex x
    climate.tau[k] ~ dexp(tau1)
  }
  for (s in 1:n.genus.rm) {
    mu.genus.c[s] ~ dexp(mu.c)       # genus-specific vertex y
  }
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.host.a[host.rm[i]] * (mu.climate.Topt[climate.rm[i]] - temp[i])^2 + mu.genus.c[genus.rm[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[climate.rm[i]]) T(0, ) 
    log_lik[i] <- logdensity.norm(trait[i], mu[i], climate.tau[climate.rm[i]])
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
jag.data <- list(trait = tick.rm$Trait, N.obs.rm = N.obs.rm, 
                 temp = tick.rm$Temp, 
                 n.host.rm = n.host.rm, host.rm = host.rm, 
                 n.climate.rm = n.climate.rm, climate.rm = climate.rm,
                 n.genus.rm = n.genus.rm, genus.rm = genus.rm)
# Fit model
h.c.g.schulzei.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                           model.file = "host.climate.genus_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                           n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
h.c.g.schulzei.fit.mcmc <- as.mcmc(h.c.g.schulzei.fit) 
closeAllConnections()

# Get posterior means for parameters for each host group
summary_fit.h.c.g <- summary(h.c.g.schulzei.fit.mcmc)$statistics
a_h.c.g_mean <- summary_fit.h.c.g[grep("^mu.host.a\\[", rownames(summary_fit.h.c.g)), "Mean"]
c_h.c.g_mean <- summary_fit.h.c.g[grep("^mu.genus.c\\[", rownames(summary_fit.h.c.g)), "Mean"]
Topt_h.c.g_mean <- summary_fit.h.c.g[grep("^mu.climate.Topt\\[", rownames(summary_fit.h.c.g)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$Temp) ,
                max(tick.abr.new$Temp) , by = 0.1)

# Get posterior samples
h.c.g.samp <- as.matrix(h.c.g.schulzei.fit.mcmc)
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

fitted_data.h.c.g.schulzei <- data.frame(
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
fitted_data.h.c.g.schulzei <- fitted_data.h.c.g.schulzei |>  
  left_join(species_to_h.c.g, by = c("Host", "Climate", "Genus"))|> 
  filter(!is.na(Species)) 

# Order factor levels
plot_data.h.c.g$Host <- factor(plot_data.h.c.g$Host, levels = c("One", "Two", "Three"))
fitted_data.h.c.g.schulzei$Host <- factor(fitted_data.h.c.g.schulzei$Host, levels = c("One", "Two", "Three"))

# Combine Host, Climate, and Size for coloring
fitted_data.h.c.g.schulzei$`Host, Climate, Genus` <- with(fitted_data.h.c.g.schulzei, paste(Host, Climate, Genus, sep = ", "))

# Plot observed data and fitted curves with facet wrap by species
p.h.c.g.schulzei <- ggplot() +
  geom_point(data = plot_data.h.c.g, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.h.c.g.schulzei, aes(x = Temperature, y = Fitted, color = `Host, Climate, Genus`), size = 1) +
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
p.h.c.g.schulzei
ggsave("tick.h.c.g.schulzei.plot.png", plot = p.h.c.g.schulzei, width = 12.5, height = 6)

# Merge the observed and predicted data for comparison
merged_data.h.c.g.schulzei <- fitted_data.h.c.g.schulzei %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", 
                         "Host" = "Host", "Climate" = "Climate",
                         "Genus" = "Genus"))

# Calculate RMSE for each species
mse_by_h.c.g.schulzei <- merged_data.h.c.g.schulzei %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
overall_rmse.h.c.g.schulzei <- merged_data.h.c.g.schulzei %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 


save(mse_by_h.c.g.schulzei,
     fitted_data.h.c.g.schulzei,
     merged_data.h.c.g.schulzei,
     h.c.g.schulzei.fit.mcmc,
     h.c.g.schulzei.fit,
     p.h.c.g.schulzei,
     file="host.climate.genus.schulzei.fit.RData")

################################### RMSE ###################################
# Data frame
table.tick.rm <- data.frame(
  Method = factor(c(rep("By Genus (RM)", 9 * 7), 
                    rep("By Host (RM)", 9 * 8), rep("By Host and Climate (RM)", 9 * 6), 
                    rep("By Host, Climate, and Genus (RM)", 9 * 4))),
  Species = factor(rep(c("Amblyomma lepidum", "Dermacentor andersoni",
                         "Dermacentor nitens", "Haemaphysalis leporispalustris",
                         "Hyalomma aegyptium", "Hyalomma dromedarii",
                         "Hyalomma impeltatum", "Hyalomma lusitanicum",
                         "Hyalomma schulzei"), times = 25)),
  # MSE = c(mse_by_genus.andersoni$MSE, mse_by_genus.nitens$MSE, 
  #          mse_by_genus.aegyptium$MSE, mse_by_genus.dromedarii$MSE, mse_by_genus.impeltatum$MSE,
  #          mse_by_genus.lusitanicum$MSE, mse_by_genus.schulzei$MSE, 
  #          mse_by_host.lepidum$MSE, mse_by_host.andersoni$MSE, mse_by_host.leporispalustris$MSE,
  #          mse_by_host.aegyptium$MSE, mse_by_host.dromedarii$MSE, mse_by_host.impeltatum$MSE,
  #          mse_by_host.lusitanicum$MSE, mse_by_host.schulzei$MSE,
  #          mse_by_host.climate.leporispalustris$MSE,
  #          mse_by_host.climate.aegyptium$MSE, mse_by_host.climate.dromedarii$MSE, 
  #          mse_by_host.climate.impeltatum$MSE,
  #          mse_by_host.climate.lusitanicum$MSE, mse_by_host.climate.schulzei$MSE,
  #          mse_by_h.c.g.aegyptium$MSE, mse_by_h.c.g$MSE,
  #          mse_by_h.c.g.lusitanicum$MSE, mse_by_h.c.g.schulzei$MSE),
  Total.RMSE = rep(as.numeric(c(overall_rmse.genus.andersoni, overall_rmse.genus.nitens,
                                overall_rmse.genus.aegyptium, overall_rmse.genus.dromedarii,
                                overall_rmse.genus.impeltatum, overall_rmse.genus.lusitanicum,
                                overall_rmse.genus.schulzei,
                                overall_rmse.host.lepidum, overall_rmse.host.andersoni,
                                overall_rmse.host.leporispalustris, overall_rmse.host.aegyptium,
                                overall_rmse.host.dromedarii, overall_rmse.host.impeltatum,
                                overall_rmse.host.lusitanicum, overall_rmse.host.schulzei,
                                overall_rmse.host.climate.leporispalustris,
                                overall_rmse.host.climate.aegyptium, overall_rmse.host.climate.dromedarii,
                                overall_rmse.host.climate.impeltatum, overall_rmse.host.climate.lusitanicum,
                                overall_rmse.host.climate.schulzei,
                                overall_rmse.h.c.g.aegyptium, overall_rmse.h.c.g.dromedarii,
                                overall_rmse.h.c.g.lusitanicum, overall_rmse.h.c.g.schulzei)), each = 9),
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

tick.DIC.rm <- data.frame(
  Method = rep(c("By Genus", "By Host",
                 "By Host and Climate",
                 "By Host, Genus, and Climate"), each = 4),
  Removed = rep(c("Hyalomma aegyptium", "Hyalomma dromedarii",
                  "Hyalomma lusitanicum", "Hyalomma schulzei"), times = 4),
  DIC = c(genus.aegyptium.fit$BUGSoutput$DIC, genus.dromedarii.fit$BUGSoutput$DIC,
          genus.lusitanicum.fit$BUGSoutput$DIC, genus.schulzei.fit$BUGSoutput$DIC,
          host.aegyptium.fit$BUGSoutput$DIC, host.dromedarii.fit$BUGSoutput$DIC,
          host.lusitanicum.fit$BUGSoutput$DIC, host.schulzei.fit$BUGSoutput$DIC,
          host.climate.aegyptium.fit$BUGSoutput$DIC, host.climate.dromedarii.fit$BUGSoutput$DIC,
          host.climate.lusitanicum.fit$BUGSoutput$DIC, host.climate.schulzei.fit$BUGSoutput$DIC,
          h.c.g.aegyptium.fit$BUGSoutput$DIC, h.c.g.dromedarii.fit$BUGSoutput$DIC,
          h.c.g.lusitanicum.fit$BUGSoutput$DIC, h.c.g.schulzei.fit$BUGSoutput$DIC),
  RMSE = as.numeric(c(overall_rmse.genus.aegyptium, overall_rmse.genus.dromedarii,
                      overall_rmse.genus.lusitanicum, overall_rmse.genus.schulzei,
                      overall_rmse.host.aegyptium, overall_rmse.host.dromedarii, 
                      overall_rmse.host.lusitanicum, overall_rmse.host.schulzei,
                      overall_rmse.host.climate.aegyptium, overall_rmse.host.climate.dromedarii,
                      overall_rmse.host.climate.lusitanicum, overall_rmse.host.climate.schulzei,
                      overall_rmse.h.c.g.aegyptium, overall_rmse.h.c.g.dromedarii,
                      overall_rmse.h.c.g.lusitanicum, overall_rmse.h.c.g.schulzei))
)

tick.DIC.rm |> 
  ggplot(aes(x = Removed, y = RMSE, fill = Method))+
  geom_bar(stat = "identity", position = "dodge",
           color = "black")+
  theme_minimal() +
  scale_fill_viridis_d(option = "mako") +
  coord_cartesian(ylim = c(6, 11)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )
ggsave("tick.RMSE.RM.png", width = 10, height = 8)

tick.DIC.rm <- tick.DIC.rm %>%
  group_by(Removed) %>%
  mutate(Delta.DIC = DIC - min(DIC)) %>%
  ungroup()

tick.DIC.rm |> 
  ggplot(aes(x = Removed, y = Delta.DIC, fill = Method))+
  geom_bar(stat = "identity", position = "dodge", color = "black")+
  theme_minimal() +
  labs(y = expression(Delta~DIC)) +
  scale_fill_viridis_d(option = "mako") +
  #coord_cartesian(ylim = c(7500, 10600)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )
ggsave("tick.DIC.RM.png", width = 10, height = 8)

save(table.tick.rm,
     tick.DIC.rm, file = "Dataset.Tick.RM.RData")








###################################
###################################
###################################
################################### Full Data- Out of sample prediction ###################################
tick.rm <- tick.abr.new[tick.abr.new$Species != "Hyalomma schulzei", ]
N.obs.rm <- length(tick.rm$Trait)
unique(tick.rm$Species)
# Model using full data
sink("full_data_robin.txt")
cat("
model {
  # Priors
   a ~ dexp(0.1)             # Width
  Topt ~ dnorm(30, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  c ~ dexp(0.001)      # Oviposition rate at Topt
  
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], tau1) T(0, )  
  }}
", file = "full_data_robin.txt")


# Parameters to monitor
parameters <- c("a", "c", "tau1", "Topt")


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
jag.data <- list(trait = tick.rm$Trait, N.obs.rm = N.obs.rm, 
                 temp = tick.rm$Temp)

# Fit fill model
full.data.schulzei.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                               model.file = "full_data_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                               n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
full.data.schulzei.fit.mcmc <- as.mcmc(full.data.schulzei.fit) 
closeAllConnections()


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
p.schulzei

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
  mu.c ~ dexp(10)          # Mean for climate effect c (vertex y)
  mu.Topt ~ dnorm(30, 0.1)          # Mean for climate effects Topt (vertex x)
  tau.Topt ~ dexp(0.01)             # Tau for climate effect Topt (vertex x)
  mu.a ~ dexp(0.01)            # Mean for climate effects a (width)

  for (j in 1:n.climate.rm) {
    mu.climate.c[j] ~ dexp(mu.c)       # Climate-specific vertex y
    mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific vertex x
    mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific width
  }
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.climate.a[climate.rm[i]] * (mu.climate.Topt[climate.rm[i]] - temp[i])^2 + mu.climate.c[climate.rm[i]]
    trait[i] ~ dnorm(mu[i], tau1)  T(0, )
  }}
", file = "climate_robin.txt")


# Parameters to monitor
parameters <- c("mu.climate.a", "mu.climate.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt",
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

plot_data.climate$Climate <- factor(plot_data.climate$Climate, levels = c("Mixed", "Temperate", "Subtropical", "Tropical"))
fitted_data.climate.schulzei$Climate <- factor(fitted_data.climate.schulzei$Climate, levels = c("Mixed", "Temperate", "Subtropical", "Tropical"))
# Plot observed data and fitted curves with facet wrap by species
p.climate.schulzei <- ggplot() +
  geom_point(data = plot_data.climate, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.climate.schulzei, aes(x = Temperature, y = Fitted, color = Climate), size = 1) +
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
p.climate.schulzei


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


################################### By Size- Out of sample prediction ###################################
tick.rm <- tick.abr.new[tick.abr.new$Species != "Hyalomma schulzei", ]
size.rm <- as.factor(tick.rm$Size)
n.size.rm <- length(unique(tick.rm$Size))
N.obs.rm <- length(tick.rm$Trait)
unique(tick.rm$Species)
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

  for (j in 1:n.size.rm) {
    mu.size.c[j] ~ dexp(mu.c)       # Size-specific vertex y
    mu.size.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Size-specific vertex x
    mu.size.a[j] ~ dexp(mu.a)    # Size-specific width
  }
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.size.a[size.rm[i]] * (mu.size.Topt[size.rm[i]] - temp[i])^2 + mu.size.c[size.rm[i]]
    trait[i] ~ dnorm(mu[i], tau1) T(0, ) 
  }}
", file = "size_robin.txt")


# Parameters to monitor
parameters <- c("mu.size.a", "mu.size.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt",
                "tau1", "mu.size.Topt")


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
jag.data <- list(trait = tick.rm$Trait, N.obs.rm = N.obs.rm, 
                 temp = tick.rm$Temp, 
                 n.size.rm = n.size.rm, size.rm = size.rm)
# Fit model
size.schulzei.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                          model.file = "size_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                          n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
size.schulzei.fit.mcmc <- as.mcmc(size.schulzei.fit)
closeAllConnections()

# Get posterior means for parameters for each size 
summary_fit.size.rm <- summary(size.schulzei.fit.mcmc)$statistics
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
fitted_data.size.schulzei <- data.frame(Temperature = rep(temp_seq, n.size.rm), 
                                        Fitted = c(size1.rm(temp_seq),
                                                   size1.rm(temp_seq),
                                                   size3.rm(temp_seq)),
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
fitted_data.size.schulzei <- data.frame(
  Temperature = rep(temp_seq, n.size.rm),
  Fitted = c(size1.rm(temp_seq),
             size2.rm(temp_seq),
             size3.rm(temp_seq)),
  Size = factor(rep(c("Large", "Medium", "Small"), each = length(temp_seq)))
) %>%
  left_join(species_to_size, by = "Size")

plot_data.size$Size <- factor(plot_data.size$Size, levels = c("Small", "Medium", "Large"))
fitted_data.size.schulzei$Size <- factor(fitted_data.size.schulzei$Size, levels = c("Small", "Medium", "Large"))
# Plot observed data and fitted curves with facet wrap by species
p.size.schulzei <- ggplot() +
  geom_point(data = plot_data.size, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.size.schulzei, aes(x = Temperature, y = Fitted, color = Size), size = 1) +
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
p.size.schulzei

# Merge the observed and predicted data for comparison
merged_data.size.schulzei <- fitted_data.size.schulzei %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Size" = "Size"))

# Calculate RMSE for each species
rmse_by_size.schulzei <- merged_data.size.schulzei %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_size.schulzei)
# Average RMSE across species
mean(rmse_by_size.schulzei$RMSE)

save(rmse_by_size.schulzei,
     fitted_data.size.schulzei,
     merged_data.size.schulzei,
     size.schulzei.fit.mcmc,
     size.schulzei.fit,
     p.size.schulzei,
     file="schulzei.size.RData")

################################### By Continent- Out of sample prediction ###################################
tick.rm <- tick.abr.new[tick.abr.new$Species != "Hyalomma schulzei", ]
continent.rm <- as.factor(tick.rm$Continent)
n.continent.rm <- length(unique(tick.rm$Continent))
N.obs.rm <- length(tick.rm$Trait)
unique(tick.rm$Species)

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

  for (j in 1:n.continent.rm) {
    mu.continent.c[j] ~ dexp(mu.c)       # Continent-specific vertex y
    mu.continent.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Continent-specific vertex x
    mu.continent.a[j] ~ dexp(mu.a)    # Continent-specific width
  }
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.continent.a[continent.rm[i]] * (mu.continent.Topt[continent.rm[i]] - temp[i])^2 + mu.continent.c[continent.rm[i]]
    trait[i] ~ dnorm(mu[i], tau1) T(0, ) 
  }}
", file = "continent_robin.txt")


# Parameters to monitor
parameters <- c("mu.continent.a", "mu.continent.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt",
                "tau1", "mu.continent.Topt")


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
jag.data <- list(trait = tick.rm$Trait, N.obs.rm = N.obs.rm, 
                 temp = tick.rm$Temp, 
                 n.continent.rm = n.continent.rm, continent.rm = continent.rm)
# Fit model
continent.schulzei.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                               model.file = "continent_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                               n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
continent.schulzei.fit.mcmc <- as.mcmc(continent.schulzei.fit)
closeAllConnections()

# Get posterior means for parameters for different continent groups
summary_fit.continent.rm <- summary(continent.schulzei.fit.mcmc)$statistics
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
fitted_data.continent.schulzei <- data.frame(Temperature = rep(temp_seq, n.continent.rm), 
                                             Fitted = c(continent1.rm(temp_seq),
                                                        continent2.rm(temp_seq),
                                                        continent3.rm(temp_seq)),
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
fitted_data.continent.schulzei <- data.frame(
  Temperature = rep(temp_seq, n.continent.rm),
  Fitted = c(continent1.rm(temp_seq),
             continent2.rm(temp_seq),
             continent3.rm(temp_seq)),
  Continent = factor(rep(c("More than two", "One", "Two"), each = length(temp_seq)))
) %>%
  left_join(species_to_continent, by = "Continent")

plot_data.continent$Continent <- factor(plot_data.continent$Continent, levels = c("One", "Two", "More than two"))
fitted_data.continent.schulzei$Continent <- factor(fitted_data.continent.schulzei$Continent, levels = c("One", "Two", "More than two"))
# Plot observed data and fitted curves with facet wrap by species
p.continent.schulzei <- ggplot() +
  geom_point(data = plot_data.continent, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.continent.schulzei, aes(x = Temperature, y = Fitted, color = Continent), size = 1) +
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
p.continent.schulzei



# Merge the observed and predicted data for comparison
merged_data.continent.schulzei <- fitted_data.continent.schulzei %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Continent" = "Continent"))

# Calculate RMSE for each species
rmse_by_continent.schulzei <- merged_data.continent.schulzei %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_continent.schulzei)
# Average RMSE across species
mean(rmse_by_continent.schulzei$RMSE)

save(rmse_by_continent.schulzei,
     fitted_data.continent.schulzei,
     merged_data.continent.schulzei,
     continent.schulzei.fit.mcmc,
     continent.schulzei.fit,
     p.continent.schulzei,
     file="schulzei.continent.RData")

################################### By Host and Genus- Out of sample prediction ###################################
tick.rm <- tick.abr.new[tick.abr.new$Species != "Hyalomma schulzei", ]
host.rm <- as.factor(tick.rm$Host)
n.host.rm <- length(unique(tick.rm$Host))
genus.rm <- as.factor(tick.rm$Genus)
n.genus.rm <- length(unique(tick.rm$Genus))
N.obs.rm <- length(tick.rm$Trait)
unique(tick.rm$Species)

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

  for (j in 1:n.host.rm) {
    mu.host.a[j] ~ dexp(mu.a)    # Host-specific width
  }
  for (k in 1:n.genus.rm) {
    mu.genus.c[k] ~ dexp(mu.c)       # Genus-specific vertex y
    mu.genus.Topt[k] ~ dnorm(mu.Topt, tau.Topt)    # Genus-specific vertex x
  }
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.host.a[host.rm[i]] * (mu.genus.Topt[genus.rm[i]] - temp[i])^2 + mu.genus.c[genus.rm[i]]
    trait[i] ~ dnorm(mu[i], tau1) T(0, ) 
  }}
", file = "host.genus_robin.txt")

# Parameters to monitor
parameters <- c("mu.host.a", "mu.genus.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt",
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
jag.data <- list(trait = tick.rm$Trait, N.obs.rm = N.obs.rm, 
                 temp = tick.rm$Temp, 
                 n.host.rm = n.host.rm, host.rm = host.rm, n.genus.rm = n.genus.rm, genus.rm = genus.rm)
# Fit model
host.genus.schulzei.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                                model.file = "host.genus_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                                n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
host.genus.schulzei.fit.mcmc <- as.mcmc(host.genus.schulzei.fit) 
closeAllConnections()

# Get posterior means for parameters for each host group
summary_fit.host.genus <- summary(host.genus.schulzei.fit.mcmc)$statistics
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
fitted_data.host.genus.schulzei <- data.frame(Temperature = rep(temp_seq, n.host * n.genus), 
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
fitted_data.host.genus.schulzei <- data.frame(Temperature = rep(temp_seq, n.host * n.genus), 
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
fitted_data.host.genus.schulzei$Host <- factor(fitted_data.host.genus.schulzei$Host, levels = c("One", "Two", "Three"))

# Combine genus and Size for coloring
fitted_data.host.genus.schulzei$`Host, Genus` <- with(fitted_data.host.genus.schulzei, paste(Host, Genus, sep = ", "))
fitted_data.host.genus.schulzei$`Host, Genus` <- factor(fitted_data.host.genus.schulzei$`Host, Genus`, levels = c("One, Dermacentor", 
                                                                                                  "Two, Hyalomma", "Three, Amblyomma",
                                                                                                  "Three, Dermacentor", "Three, Haemaphysalis",
                                                                                                  "Three, Hyalomma"))

# Plot observed data and fitted curves with facet wrap by species
p.host.genus.schulzei <- ggplot() +
  geom_point(data = plot_data.host.genus, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.host.genus.schulzei, aes(x = Temperature, y = Fitted, color = `Host, Genus`), size = 1) +
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
p.host.genus.schulzei

# Merge the observed and predicted data for comparison
merged_data.host.genus.schulzei <- fitted_data.host.genus.schulzei %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", "Host" = "Host", "Genus" = "Genus"))

# Calculate RMSE for each species
rmse_by_host.genus.schulzei <- merged_data.host.genus.schulzei %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_host.genus.schulzei)
# Average RMSE across species
mean(rmse_by_host.genus.schulzei$RMSE)


save(rmse_by_host.genus.schulzei,
     fitted_data.host.genus.schulzei,
     merged_data.host.genus.schulzei,
     host.genus.schulzei.fit.mcmc,
     host.genus.schulzei.fit,
     p.host.genus.schulzei,
     file="host.genus.schulzei.fit.RData")

################################### By Genus and Climate- Out of sample prediction ###################################
tick.rm <- tick.abr.new[tick.abr.new$Species != "Hyalomma schulzei", ]
climate.rm <- as.factor(tick.rm$Climate)
n.climate.rm <- length(unique(tick.rm$Climate))
genus.rm <- as.factor(tick.rm$Genus)
n.genus.rm <- length(unique(tick.rm$Genus))
N.obs.rm <- length(tick.rm$Trait)
unique(tick.rm$Species)

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

  for (j in 1:n.genus.rm) {
    mu.genus.c[j] ~ dexp(mu.c)       # genus-specific vertex y
    mu.genus.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # genus-specific vertex x
  }
  for (k in 1:n.climate.rm) {
    mu.climate.a[k] ~ dexp(mu.a)    # Cimate-specific width
  }
  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.climate.a[climate.rm[i]] * (mu.genus.Topt[genus.rm[i]] - temp[i])^2 + mu.genus.c[genus.rm[i]]
    trait[i] ~ dnorm(mu[i], tau1) T(0, )  
  }}
", file = "genus.climate_robin.txt")

# Parameters to monitor
parameters <- c("mu.climate.a", "mu.genus.c", "mu.a", "mu.c",
                "mu.Topt", "tau.Topt",
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
jag.data <- list(trait = tick.rm$Trait, N.obs.rm = N.obs.rm, 
                 temp = tick.rm$Temp, 
                 n.genus.rm = n.genus.rm, genus.rm = genus.rm, 
                 n.climate.rm = n.climate.rm, climate.rm = climate.rm)
# Fit model
genus.climate.schulzei.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                                   model.file = "genus.climate_robin.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                                   n.iter = ni, DIC = T, working.directory = getwd())
# Make mcmc object
genus.climate.schulzei.fit.mcmc <- as.mcmc(genus.climate.schulzei.fit) 
closeAllConnections()


# Get posterior means for parameters for each genus group
summary_fit.genus.climate <- summary(genus.climate.schulzei.fit.mcmc)$statistics
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
fitted_data.genus.climate.schulzei <- data.frame(Temperature = rep(temp_seq, n.genus * n.climate), 
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
fitted_data.genus.climate.schulzei <- data.frame(Temperature = rep(temp_seq, n.genus * n.climate), 
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
fitted_data.genus.climate.schulzei$`Genus, Climate` <- with(fitted_data.genus.climate.schulzei, paste(Genus, Climate, sep = ", "))


# Plot observed data and fitted curves with facet wrap by species
p.genus.climate.schulzei <- ggplot() +
  geom_point(data = plot_data.genus.climate, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.genus.climate.schulzei, aes(x = Temperature, y = Fitted, color = `Genus, Climate`), size = 1) +
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
p.genus.climate.schulzei


# Merge the observed and predicted data for comparison
merged_data.genus.climate.schulzei <- fitted_data.genus.climate.schulzei %>%
  left_join(tick, by = c("Temperature" = "Temp", "Species" = "Species", 
                         "Genus" = "Genus", "Climate" = "Climate"))

# Calculate RMSE for each species
rmse_by_genus.climate.schulzei <- merged_data.genus.climate.schulzei %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_genus.climate.schulzei)
# Average RMSE across species
mean(rmse_by_genus.climate.schulzei$RMSE)



save(rmse_by_genus.climate.schulzei,
     fitted_data.genus.climate.schulzei,
     merged_data.genus.climate.schulzei,
     genus.climate.schulzei.fit.mcmc,
     genus.climate.schulzei.fit,
     p.genus.climate.schulzei,
     file="genus.climate.schulzei.fit.RData")


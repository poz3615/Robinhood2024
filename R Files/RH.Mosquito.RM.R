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
library(ggridges)
library(HDInterval)
library(gridExtra)


######################################### Data Cleaning #########################################
# Mosquito development rate in days, a lot Small so filtered dataset to >1 to try to get the hang of it
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
  mutate(Wing.Size = case_when(
    interactor1 == "Aedes aegypti" ~ "Medium",
    interactor1 == "Aedes albopictus" ~ "Small",
    interactor1 == "Aedes japonicus japonicus" ~ "Large",
    interactor1 == "Aedes atropalpus" ~ "Large",
    interactor1 == "Aedes krombeini" ~ "Medium",
    interactor1 == "Aedes notoscriptus" ~ "Medium",
    interactor1 == "Aedes camptorhynchus" ~ "Large",
    TRUE ~ NA_character_
  ))

mosquito.data <- mosquito.data %>%
  mutate(Climate = case_when(
    interactor1 == "Aedes aegypti" ~ "Tropical",
    interactor1 == "Aedes albopictus" ~ "Tropical",
    interactor1 == "Aedes japonicus japonicus" ~ "Temperate",
    interactor1 == "Aedes atropalpus" ~ "Temperate",
    interactor1 == "Aedes krombeini" ~ "Tropical",
    interactor1 == "Aedes notoscriptus" ~ "Subtropical",
    interactor1 == "Aedes camptorhynchus" ~ "Subtropical",
    TRUE ~ NA_character_
  ))

mosquito.data$Climate <- as.factor(mosquito.data$Climate)
mosquito.data$Wing.Size <- as.factor(mosquito.data$Wing.Size)
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
  # Sample Size
  num_points <- m.filtered_df1$interactor1number[i]
  # Climate
  Climate <- m.filtered_df1$Climate[i]
  # Size
  Wing.Size <- m.filtered_df1$Wing.Size[i]
  # Generate sample Size of random values from truncated normal (0, Inf) distribution
  generated_values <- rtruncnorm(num_points, a=0, b=Inf, mean = mean_value, sd = sd_value)
  # Create a temporary dataframe for the generated values
  temp_df <- data.frame(
    # Simulated values
    Trait = generated_values,
    # Carry over sd, species, temperature, Climate, and Size
    sd = sd_value,
    Species = interactor1,
    Temp = interactor1temp,
    Climate = Climate,
    Wing.Size = Wing.Size
  )
  # Append the generated rows to the expanded dataset
  m.expanded_df1 <- rbind(m.expanded_df1, temp_df)
}

# Take the individual data
m.filtered_df <- mosquito.data |> filter(interactor1number == 1)
# Summarize it, grouping by species and temperature, to get the mean, sd, and sample Size
m.summarized_df <- m.filtered_df %>%
  group_by(interactor1, Temp, Climate, Wing.Size) %>%
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
    Climate = .$Climate,
    Wing.Size = .$Wing.Size  
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
m.climate.pd <- as.factor(mosquito.pd$Climate)
# Number of Climates
m.n.climate.pd <- length(unique(mosquito.pd$Climate))

# Size number for hierarchical model
m.size.pd <- as.factor(mosquito.pd$Size)
# Small: 0-3
# Medium: 3.01-5
# Large: 5.01+
# Number of Sizes
m.n.size.pd <- length(unique(mosquito.pd$Wing.Size))

m.N.obs.pd <- length(mosquito.pd$Trait)

######################################### By Climate- Out of sample prediction Aegypti #########################################
aedes.rm.pd <- mosquito.pd[mosquito.pd$Species != "Aedes aegypti", ]
m.climate.pd.rm <- as.factor(aedes.rm.pd$Climate)
m.n.climate.rm.pd <- length(unique(aedes.rm.pd$Climate))
m.N.obs.pd.rm <- length(aedes.rm.pd$Trait)
unique(aedes.rm.pd$interactor1)
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait

  # # Climate-specific
  # mu.c ~ dexp(10)             # Mean for climate effects c
  # mu.Topt ~ dnorm(30, 0.1)             # Mean for climate effects Topt
  # tau.Topt ~ dexp(0.01)             # Tau for climate effect Topt
  # mu.a ~ dexp(0.01)            # Mean for climate effects Topt
  mu.c ~ dgamma(10, 1)
  mu.Topt ~ dexp(0.001)
 mu.a ~ dexp(10)
  mu.c.rate ~ dgamma(2, 1)

  for (j in 1:m.n.climate.rm.pd) {
    # mu.climate[j] ~ dexp(mu.c)       # Climate-specific low point
    # mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific temp
    # mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific temp
    climate.tau[j] ~ dexp(tau1)
    mu.climate[j] ~ dgamma(mu.c, mu.c.rate)
    mu.climate.Topt[j] ~ dexp(mu.Topt)
    mu.climate.a[j] ~ dexp(mu.a)
  }

  # Likelihood
  for (i in 1:m.N.obs.pd.rm) {
    mu[i] <- mu.climate.a[m.climate.pd.rm[i]] * exp(- mu.climate.Topt[m.climate.pd.rm[i]] * temp[i]) + mu.climate[m.climate.pd.rm[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[m.climate.pd.rm[i]]) T(0, )
  }}
", file = "aedes.by.climate.txt")


# Parameters to monitor
parameters <- c("mu.a", "mu.c", "mu.climate", "climate.tau", "mu.c.rate",
                "tau1", "mu.Topt", "mu.climate.Topt", "mu.climate.a")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 40000 # number of iterations in each chain
nb <- 20000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Initial Values

inits = vector('list', nc)
GenInits = function() {
  # mu.a <- rexp(1, 0.01)
  # mu.c <- rexp(1, 10)
  tau1 <- rexp(1, 0.01)
  # mu.Topt <- rnorm(1, 30, 1/0.1)
  mu.a <- rexp(1, 10)
  mu.Topt <- rexp(1, 0.001)
  mu.c <- rgamma(1, 10, 1)
  mu.c.rate <- rgamma(1, 2, 1)
  #tau.Topt <- rexp(1, 0.01)
  list(
    mu.a = mu.a,
    mu.c = mu.c,
    tau1 = tau1,
    mu.Topt = mu.Topt,
    mu.c.rate = mu.c.rate
    #tau.Topt = tau.Topt
  )
}
for(i in 1:nc){
  inits[[i]] = GenInits()
}

# Bundle all data in a list for JAGS
jag.data <- list(trait = aedes.rm.pd$Trait, m.N.obs.pd.rm = m.N.obs.pd.rm, temp = aedes.rm.pd$Temp,
                 m.n.climate.rm.pd = m.n.climate.rm.pd, m.climate.pd.rm = m.climate.pd.rm)

Climate.aegypti.pd.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters,
                               model.file="aedes.by.Climate.txt", n.thin=nt, n.chains=nc, n.burnin=nb,
                               n.iter=ni, DIC=T, working.directory=getwd())
Climate.aegypti.pd.fit.mcmc <- as.mcmc(Climate.aegypti.pd.fit) ## makes an "mcmc" object
closeAllConnections()

Climate.aegypti.pd.fit.DIC <- Climate.aegypti.pd.fit$BUGSoutput$DIC

# Get posterior samples
m.climate.samp.aegypti <- as.matrix(Climate.aegypti.pd.fit.mcmc)
# Get column indices
m.a_climate.aegypti <- grep("^mu\\.climate\\.a\\[", colnames(m.climate.samp.aegypti))
m.Topt_climate.aegypti <- grep("^mu\\.climate\\.Topt\\[", colnames(m.climate.samp.aegypti))
m.c_climate.aegypti <- grep("^mu\\.climate\\[", colnames(m.climate.samp.aegypti))
# Generate temperature sequence 
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)
# Prepare a list to store data frames
m.climate_dfs.aegypti <- list()
for (j in 1:m.n.climate.rm.pd) {
  # Retrieve parameter vectors for climate j
  m.a_samps <- m.climate.samp.aegypti[, m.a_climate.aegypti[j]]
  m.Topt_samps <- m.climate.samp.aegypti[, m.Topt_climate.aegypti[j]]
  m.c_samps <- m.climate.samp.aegypti[, m.c_climate.aegypti[j]]
  # For each temperature in temp_seq, evaluate across ALL posterior samples
  m.climate_eval.aegypti <- sapply(m.temp_seq, function(temp) {
    m.a_samps * exp(- m.Topt_samps * temp) + m.c_samps
  })
  # Compute pointwise posterior means and HDIs
  m.mean_curve.climate.aegypti <- apply(m.climate_eval.aegypti, 2, mean)
  m.climate_dfs.aegypti[[j]] <- m.mean_curve.climate.aegypti
}
# Combine all climate
m.climate_dfs.aegypti <- unlist(m.climate_dfs.aegypti)
# Optionally, convert to a data frame for easier plotting
fitted_data.Climate.aegypti.pd  <- data.frame(Temperature = rep(m.temp_seq, m.n.climate.rm.pd),
                                       Fitted = m.climate_dfs.aegypti,
                                       Climate = factor(rep(c("Subtropical", "Temperate", "Tropical"), each = length(m.temp_seq))))

m.plot_data.Climate.pd <- data.frame(
  Temperature = mosquito.pd$Temp,
  Trait = mosquito.pd$Trait,
  Climate = as.factor(mosquito.pd$Climate),
  Species = as.factor(mosquito.pd$Species)
)
m.plot_data.Climate <- data.frame(
  Temperature = mosquito.data$Temp,
  Trait = mosquito.data$Trait,
  Climate = as.factor(mosquito.data$Climate),
  Species = as.factor(mosquito.data$interactor1)
)

m.species.pd_to_Climate <- unique(m.plot_data.Climate.pd[, c("Species", "Climate")])  # Ensure species matches correct Climate
m.species.pd_to_Climate$Climate <- as.factor(m.species.pd_to_Climate$Climate)
# Merge fitted Climate curves with species-to-Climate mapping
fitted_data.Climate.aegypti.pd <- data.frame(Temperature = rep(m.temp_seq, m.n.climate.rm.pd),
                                             Fitted = m.climate_dfs.aegypti,
                                             Climate = factor(rep(c("Subtropical", "Temperate", "Tropical"), each = length(m.temp_seq)))) %>%
  left_join(m.species.pd_to_Climate, by = "Climate")

m.plot_data.Climate.pd$Climate <- factor(m.plot_data.Climate.pd$Climate, levels = c("Temperate", "Subtropical", "Tropical"))
fitted_data.Climate.aegypti.pd$Climate <- factor(fitted_data.Climate.aegypti.pd$Climate, levels = c("Temperate", "Subtropical", "Tropical"))
# Plot observed data and fitted curves with facet wrap by species
p.Climate.aegypti.pd <- ggplot() +
  geom_point(data = m.plot_data.Climate.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.Climate.aegypti.pd, aes(x = Temperature, y = Fitted, color = Climate), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label Size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )  
p.Climate.aegypti.pd


# Merge the observed and predicted data for comparison
merged_data.Climate.aegypti.pd <- fitted_data.Climate.aegypti.pd %>%
  left_join(m.plot_data.Climate, by = c("Temperature", "Climate", "Species"))


# Calculate RMSE for each Climate
m.mse_by_Climate.aegypti.pd <- merged_data.Climate.aegypti.pd |>
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.Climate.aegypti.pd <- merged_data.Climate.aegypti.pd |>
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE)))

# save(m.mse_by_Climate.aegypti.pd,
#      m.overall_rmse.Climate.aegypti.pd,
#      fitted_data.Climate.aegypti.pd,
#      merged_data.Climate.aegypti.pd,
#      Climate.aegypti.pd.fit.DIC,
#      Climate.aegypti.pd.fit,
#      p.Climate.aegypti.pd,
#      file="aegypti.pd.Climate.RData")

save(Climate.aegypti.pd.fit,
     m.climate.samp.aegypti,
     fitted_data.Climate.aegypti.pd,
     p.Climate.aegypti.pd,
     m.mse_by_Climate.aegypti.pd,
     m.overall_rmse.Climate.aegypti.pd,
     Climate.aegypti.pd.fit.DIC,
     file="exp.aegypti.pd.Climate.RData")


######################################### By Climate- Out of sample prediction Albopictus #########################################
aedes.rm.pd <- mosquito.pd[mosquito.pd$Species != "Aedes albopictus", ]
m.climate.pd.rm <- as.factor(aedes.rm.pd$Climate)
m.n.climate.rm.pd <- length(unique(aedes.rm.pd$Climate))
m.N.obs.pd.rm <- length(aedes.rm.pd$Trait)
unique(aedes.rm.pd$interactor1)
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait

  # # Climate-specific
  # mu.c ~ dexp(10)             # Mean for climate effects c
  # mu.Topt ~ dnorm(30, 0.1)             # Mean for climate effects Topt
  # tau.Topt ~ dexp(0.01)             # Tau for climate effect Topt
  # mu.a ~ dexp(0.01)            # Mean for climate effects Topt
  mu.c ~ dgamma(10, 1)
  mu.Topt ~ dexp(0.001)
 mu.a ~ dexp(10)
  mu.c.rate ~ dgamma(2, 1)

  for (j in 1:m.n.climate.rm.pd) {
    # mu.climate[j] ~ dexp(mu.c)       # Climate-specific low point
    # mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific temp
    # mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific temp
    climate.tau[j] ~ dexp(tau1)
    mu.climate[j] ~ dgamma(mu.c, mu.c.rate)
    mu.climate.Topt[j] ~ dexp(mu.Topt)
    mu.climate.a[j] ~ dexp(mu.a)
  }

  # Likelihood
  for (i in 1:m.N.obs.pd.rm) {
    mu[i] <- mu.climate.a[m.climate.pd.rm[i]] * exp(- mu.climate.Topt[m.climate.pd.rm[i]] * temp[i]) + mu.climate[m.climate.pd.rm[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[m.climate.pd.rm[i]]) T(0, )
  }}
", file = "aedes.by.climate.txt")

# Parameters to monitor
parameters <- c("mu.a", "mu.c", "mu.climate", "climate.tau",
                "tau1", "mu.Topt", "mu.c.rate", "mu.climate.Topt", "mu.climate.a")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 40000 # number of iterations in each chain
nb <- 20000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Initial Values
inits = vector('list', nc)
GenInits = function() {
  # mu.a <- rexp(1, 0.01)
  # mu.c <- rexp(1, 10)
  tau1 <- rexp(1, 0.01)
  # mu.Topt <- rnorm(1, 30, 1/0.1)
  mu.a <- rexp(1, 10)
  mu.Topt <- rexp(1, 0.001)
  mu.c <- rgamma(1, 10, 1)
  mu.c.rate <- rgamma(1, 2, 1)
  #tau.Topt <- rexp(1, 0.01)
  list(
    mu.a = mu.a,
    mu.c = mu.c,
    tau1 = tau1,
    mu.Topt = mu.Topt,
    mu.c.rate = mu.c.rate
    #tau.Topt = tau.Topt
  )
}
for(i in 1:nc){
  inits[[i]] = GenInits()
}

# Bundle all data in a list for JAGS
jag.data <- list(trait = aedes.rm.pd$Trait, m.N.obs.pd.rm = m.N.obs.pd.rm, temp = aedes.rm.pd$Temp,
                 m.n.climate.rm.pd = m.n.climate.rm.pd, m.climate.pd.rm = m.climate.pd.rm)

Climate.albopictus.pd.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters,
                                  model.file="aedes.by.Climate.txt", n.thin=nt, n.chains=nc, n.burnin=nb,
                                  n.iter=ni, DIC=T, working.directory=getwd())
Climate.albopictus.pd.fit.mcmc <- as.mcmc(Climate.albopictus.pd.fit) ## makes an "mcmc" object
closeAllConnections()

Climate.albopictus.pd.fit.DIC <- Climate.albopictus.pd.fit$BUGSoutput$DIC

# Get posterior samples
m.climate.samp.albopictus <- as.matrix(Climate.albopictus.pd.fit.mcmc)
# Get column indices
m.a_climate.albopictus <- grep("^mu\\.climate\\.a\\[", colnames(m.climate.samp.albopictus))
m.Topt_climate.albopictus <- grep("^mu\\.climate\\.Topt\\[", colnames(m.climate.samp.albopictus))
m.c_climate.albopictus <- grep("^mu\\.climate\\[", colnames(m.climate.samp.albopictus))
# Generate temperature sequence 
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)
# Prepare a list to store data frames
m.climate_dfs.albopictus <- list()
for (j in 1:m.n.climate.rm.pd) {
  # Retrieve parameter vectors for climate j
  m.a_samps <- m.climate.samp.albopictus[, m.a_climate.albopictus[j]]
  m.Topt_samps <- m.climate.samp.albopictus[, m.Topt_climate.albopictus[j]]
  m.c_samps <- m.climate.samp.albopictus[, m.c_climate.albopictus[j]]
  # For each temperature in temp_seq, evaluate across ALL posterior samples
  m.climate_eval.albopictus <- sapply(m.temp_seq, function(temp) {
    m.a_samps * exp(- m.Topt_samps * temp) + m.c_samps
  })
  # Compute pointwise posterior means and HDIs
  m.mean_curve.climate.albopictus <- apply(m.climate_eval.albopictus, 2, mean)
  m.climate_dfs.albopictus[[j]] <- m.mean_curve.climate.albopictus
}
# Combine all climate
m.climate_dfs.albopictus <- unlist(m.climate_dfs.albopictus)
# Optionally, convert to a data frame for easier plotting
fitted_data.Climate.albopictus.pd  <- data.frame(Temperature = rep(m.temp_seq, m.n.climate.rm.pd),
                                              Fitted = m.climate_dfs.albopictus,
                                              Climate = factor(rep(c("Subtropical", "Temperate", "Tropical"), each = length(m.temp_seq))))

m.plot_data.Climate.pd <- data.frame(
  Temperature = mosquito.pd$Temp,
  Trait = mosquito.pd$Trait,
  Climate = as.factor(mosquito.pd$Climate),
  Species = as.factor(mosquito.pd$Species)
)
m.plot_data.Climate <- data.frame(
  Temperature = mosquito.data$Temp,
  Trait = mosquito.data$Trait,
  Climate = as.factor(mosquito.data$Climate),
  Species = as.factor(mosquito.data$interactor1)
)

m.species.pd_to_Climate <- unique(m.plot_data.Climate.pd[, c("Species", "Climate")])  # Ensure species matches correct Climate
m.species.pd_to_Climate$Climate <- as.factor(m.species.pd_to_Climate$Climate)
# Merge fitted Climate curves with species-to-Climate mapping
fitted_data.Climate.albopictus.pd <- data.frame(Temperature = rep(m.temp_seq, m.n.climate.rm.pd),
                                             Fitted = m.climate_dfs.albopictus,
                                             Climate = factor(rep(c("Subtropical", "Temperate", "Tropical"), each = length(m.temp_seq)))) %>%
  left_join(m.species.pd_to_Climate, by = "Climate")

m.plot_data.Climate.pd$Climate <- factor(m.plot_data.Climate.pd$Climate, levels = c("Temperate", "Subtropical", "Tropical"))
fitted_data.Climate.albopictus.pd$Climate <- factor(fitted_data.Climate.albopictus.pd$Climate, levels = c("Temperate", "Subtropical", "Tropical"))
# Plot observed data and fitted curves with facet wrap by species
p.Climate.albopictus.pd <- ggplot() +
  geom_point(data = m.plot_data.Climate.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.Climate.albopictus.pd, aes(x = Temperature, y = Fitted, color = Climate), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label Size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )  
p.Climate.albopictus.pd

# Merge the observed and predicted data for comparison
merged_data.Climate.albopictus.pd <- fitted_data.Climate.albopictus.pd %>%
  left_join(m.plot_data.Climate, by = c("Temperature", "Climate", "Species"))


# Calculate RMSE for each Climate
m.mse_by_Climate.albopictus.pd <- merged_data.Climate.albopictus.pd |>
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.Climate.albopictus.pd <- merged_data.Climate.albopictus.pd |>
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE)))

# save(m.mse_by_Climate.albopictus.pd,
#      m.overall_rmse.Climate.albopictus.pd,
#      fitted_data.Climate.albopictus.pd,
#      merged_data.Climate.albopictus.pd,
#      Climate.albopictus.pd.fit,
#      Climate.albopictus.pd.fit.DIC,
#      p.Climate.albopictus.pd,
#      file="albopictus.pd.Climate.RData")

save(Climate.albopictus.pd.fit,
     m.climate.samp.albopictus,
     fitted_data.Climate.albopictus.pd,
     p.Climate.albopictus.pd,
     m.mse_by_Climate.albopictus.pd,
     m.overall_rmse.Climate.albopictus.pd,
     Climate.albopictus.pd.fit.DIC,
     file="exp.albopictus.pd.Climate.RData")

######################################### By Climate- Out of sample prediction Camptorhynchus #########################################
aedes.rm.pd <- mosquito.pd[mosquito.pd$Species != "Aedes camptorhynchus", ]
m.climate.pd.rm <- as.factor(aedes.rm.pd$Climate)
m.n.climate.rm.pd <- length(unique(aedes.rm.pd$Climate))
m.N.obs.pd.rm <- length(aedes.rm.pd$Trait)
unique(aedes.rm.pd$interactor1)
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait

  # # Climate-specific
  # mu.c ~ dexp(10)             # Mean for climate effects c
  # mu.Topt ~ dnorm(30, 0.1)             # Mean for climate effects Topt
  # tau.Topt ~ dexp(0.01)             # Tau for climate effect Topt
  # mu.a ~ dexp(0.01)            # Mean for climate effects Topt
  mu.c ~ dgamma(10, 1)
  mu.Topt ~ dexp(0.001)
 mu.a ~ dexp(10)
  mu.c.rate ~ dgamma(2, 1)

  for (j in 1:m.n.climate.rm.pd) {
    # mu.climate[j] ~ dexp(mu.c)       # Climate-specific low point
    # mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific temp
    # mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific temp
    climate.tau[j] ~ dexp(tau1)
    mu.climate[j] ~ dgamma(mu.c, mu.c.rate)
    mu.climate.Topt[j] ~ dexp(mu.Topt)
    mu.climate.a[j] ~ dexp(mu.a)
  }

  # Likelihood
  for (i in 1:m.N.obs.pd.rm) {
    mu[i] <- mu.climate.a[m.climate.pd.rm[i]] * exp(- mu.climate.Topt[m.climate.pd.rm[i]] * temp[i]) + mu.climate[m.climate.pd.rm[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[m.climate.pd.rm[i]]) T(0, )
  }}
", file = "aedes.by.climate.txt")

# Parameters to monitor
parameters <- c("mu.a", "mu.c", "mu.climate", "climate.tau",
                "tau1", "mu.Topt", "mu.c.rate", "mu.climate.Topt", "mu.climate.a")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 40000 # number of iterations in each chain
nb <- 20000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

inits = vector('list', nc)
GenInits = function() {
  # mu.a <- rexp(1, 0.01)
  # mu.c <- rexp(1, 10)
  tau1 <- rexp(1, 0.01)
  # mu.Topt <- rnorm(1, 30, 1/0.1)
  mu.a <- rexp(1, 10)
  mu.Topt <- rexp(1, 0.001)
  mu.c <- rgamma(1, 10, 1)
  mu.c.rate <- rgamma(1, 2, 1)
  #tau.Topt <- rexp(1, 0.01)
  list(
    mu.a = mu.a,
    mu.c = mu.c,
    tau1 = tau1,
    mu.Topt = mu.Topt,
    mu.c.rate = mu.c.rate
    #tau.Topt = tau.Topt
  )
}
for(i in 1:nc){
  inits[[i]] = GenInits()
}

# Bundle all data in a list for JAGS
jag.data <- list(trait = aedes.rm.pd$Trait, m.N.obs.pd.rm = m.N.obs.pd.rm, temp = aedes.rm.pd$Temp,
                 m.n.climate.rm.pd = m.n.climate.rm.pd, m.climate.pd.rm = m.climate.pd.rm)

Climate.camptorhynchus.pd.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters,
                                      model.file="aedes.by.Climate.txt", n.thin=nt, n.chains=nc, n.burnin=nb,
                                      n.iter=ni, DIC=T, working.directory=getwd())
Climate.camptorhynchus.pd.fit.mcmc <- as.mcmc(Climate.camptorhynchus.pd.fit) ## makes an "mcmc" object
closeAllConnections()

Climate.camptorhynchus.pd.fit.DIC <- Climate.camptorhynchus.pd.fit$BUGSoutput$DIC


# Get posterior samples
m.climate.samp.camptorhynchus <- as.matrix(Climate.camptorhynchus.pd.fit.mcmc)
# Get column indices
m.a_climate.camptorhynchus <- grep("^mu\\.climate\\.a\\[", colnames(m.climate.samp.camptorhynchus))
m.Topt_climate.camptorhynchus <- grep("^mu\\.climate\\.Topt\\[", colnames(m.climate.samp.camptorhynchus))
m.c_climate.camptorhynchus <- grep("^mu\\.climate\\[", colnames(m.climate.samp.camptorhynchus))
# Generate temperature sequence 
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)
# Prepare a list to store data frames
m.climate_dfs.camptorhynchus <- list()
for (j in 1:m.n.climate.rm.pd) {
  # Retrieve parameter vectors for climate j
  m.a_samps <- m.climate.samp.camptorhynchus[, m.a_climate.camptorhynchus[j]]
  m.Topt_samps <- m.climate.samp.camptorhynchus[, m.Topt_climate.camptorhynchus[j]]
  m.c_samps <- m.climate.samp.camptorhynchus[, m.c_climate.camptorhynchus[j]]
  # For each temperature in temp_seq, evaluate across ALL posterior samples
  m.climate_eval.camptorhynchus <- sapply(m.temp_seq, function(temp) {
    m.a_samps * exp(- m.Topt_samps * temp) + m.c_samps
  })
  # Compute pointwise posterior means and HDIs
  m.mean_curve.climate.camptorhynchus <- apply(m.climate_eval.camptorhynchus, 2, mean)
  m.climate_dfs.camptorhynchus[[j]] <- m.mean_curve.climate.camptorhynchus
}
# Combine all climate
m.climate_dfs.camptorhynchus <- unlist(m.climate_dfs.camptorhynchus)
# Optionally, convert to a data frame for easier plotting
fitted_data.Climate.camptorhynchus.pd  <- data.frame(Temperature = rep(m.temp_seq, m.n.climate.rm.pd),
                                                 Fitted = m.climate_dfs.camptorhynchus,
                                                 Climate = factor(rep(c("Subtropical", "Temperate", "Tropical"), each = length(m.temp_seq))))

m.plot_data.Climate.pd <- data.frame(
  Temperature = mosquito.pd$Temp,
  Trait = mosquito.pd$Trait,
  Climate = as.factor(mosquito.pd$Climate),
  Species = as.factor(mosquito.pd$Species)
)
m.plot_data.Climate <- data.frame(
  Temperature = mosquito.data$Temp,
  Trait = mosquito.data$Trait,
  Climate = as.factor(mosquito.data$Climate),
  Species = as.factor(mosquito.data$interactor1)
)

m.species.pd_to_Climate <- unique(m.plot_data.Climate.pd[, c("Species", "Climate")])  # Ensure species matches correct Climate
m.species.pd_to_Climate$Climate <- as.factor(m.species.pd_to_Climate$Climate)
# Merge fitted Climate curves with species-to-Climate mapping
fitted_data.Climate.camptorhynchus.pd <- data.frame(Temperature = rep(m.temp_seq, m.n.climate.rm.pd),
                                                Fitted = m.climate_dfs.camptorhynchus,
                                                Climate = factor(rep(c("Subtropical", "Temperate", "Tropical"), each = length(m.temp_seq)))) %>%
  left_join(m.species.pd_to_Climate, by = "Climate")

m.plot_data.Climate.pd$Climate <- factor(m.plot_data.Climate.pd$Climate, levels = c("Temperate", "Subtropical", "Tropical"))
fitted_data.Climate.camptorhynchus.pd$Climate <- factor(fitted_data.Climate.camptorhynchus.pd$Climate, levels = c("Temperate", "Subtropical", "Tropical"))
# Plot observed data and fitted curves with facet wrap by species
p.Climate.camptorhynchus.pd <- ggplot() +
  geom_point(data = m.plot_data.Climate.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.Climate.camptorhynchus.pd, aes(x = Temperature, y = Fitted, color = Climate), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label Size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )  
p.Climate.camptorhynchus.pd

# Merge the observed and predicted data for comparison
merged_data.Climate.camptorhynchus.pd <- fitted_data.Climate.camptorhynchus.pd %>%
  left_join(m.plot_data.Climate, by = c("Temperature", "Climate", "Species"))


# Calculate RMSE for each Climate
m.mse_by_Climate.camptorhynchus.pd <- merged_data.Climate.camptorhynchus.pd |>
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.Climate.camptorhynchus.pd <- merged_data.Climate.camptorhynchus.pd |>
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE)))

# save(m.mse_by_Climate.camptorhynchus.pd,
#      m.overall_rmse.Climate.camptorhynchus.pd,
#      fitted_data.Climate.camptorhynchus.pd,
#      merged_data.Climate.camptorhynchus.pd,
#      Climate.camptorhynchus.pd.fit,
#      Climate.camptorhynchus.pd.fit.DIC,
#      p.Climate.camptorhynchus.pd,
#      file="camptorhynchus.pd.Climate.RData")

save(Climate.camptorhynchus.pd.fit,
     m.climate.samp.camptorhynchus,
     fitted_data.Climate.camptorhynchus.pd,
     p.Climate.camptorhynchus.pd,
     m.mse_by_Climate.camptorhynchus.pd,
     m.overall_rmse.Climate.camptorhynchus.pd,
     Climate.camptorhynchus.pd.fit.DIC,
     file="exp.camptorhynchus.pd.Climate.RData")



######################################### By Climate- Out of sample prediction Krombeini #########################################
aedes.rm.pd <- mosquito.pd[mosquito.pd$Species != "Aedes krombeini", ]
m.climate.pd.rm <- as.factor(aedes.rm.pd$Climate)
m.n.climate.rm.pd <- length(unique(aedes.rm.pd$Climate))
m.N.obs.pd.rm <- length(aedes.rm.pd$Trait)
unique(aedes.rm.pd$interactor1)
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait

  # # Climate-specific
  # mu.c ~ dexp(10)             # Mean for climate effects c
  # mu.Topt ~ dnorm(30, 0.1)             # Mean for climate effects Topt
  # tau.Topt ~ dexp(0.01)             # Tau for climate effect Topt
  # mu.a ~ dexp(0.01)            # Mean for climate effects Topt
  mu.c ~ dgamma(10, 1)
  mu.Topt ~ dexp(0.001)
 mu.a ~ dexp(10)
  mu.c.rate ~ dgamma(2, 1)

  for (j in 1:m.n.climate.rm.pd) {
    # mu.climate[j] ~ dexp(mu.c)       # Climate-specific low point
    # mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific temp
    # mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific temp
    climate.tau[j] ~ dexp(tau1)
    mu.climate[j] ~ dgamma(mu.c, mu.c.rate)
    mu.climate.Topt[j] ~ dexp(mu.Topt)
    mu.climate.a[j] ~ dexp(mu.a)
  }

  # Likelihood
  for (i in 1:m.N.obs.pd.rm) {
    mu[i] <- mu.climate.a[m.climate.pd.rm[i]] * exp(- mu.climate.Topt[m.climate.pd.rm[i]] * temp[i]) + mu.climate[m.climate.pd.rm[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[m.climate.pd.rm[i]]) T(0, )
  }}
", file = "aedes.by.climate.txt")
# Parameters to monitor
parameters <- c("mu.a", "mu.c", "mu.climate", "climate.tau",
                "tau1", "mu.Topt", "mu.c.rate", "mu.climate.Topt", "mu.climate.a")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 40000 # number of iterations in each chain
nb <- 20000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Initial Values
inits = vector('list', nc)
GenInits = function() {
  # mu.a <- rexp(1, 0.01)
  # mu.c <- rexp(1, 10)
  tau1 <- rexp(1, 0.01)
  # mu.Topt <- rnorm(1, 30, 1/0.1)
  mu.a <- rexp(1, 10)
  mu.Topt <- rexp(1, 0.001)
  mu.c <- rgamma(1, 10, 1)
  mu.c.rate <- rgamma(1, 2, 1)
  #tau.Topt <- rexp(1, 0.01)
  list(
    mu.a = mu.a,
    mu.c = mu.c,
    tau1 = tau1,
    mu.Topt = mu.Topt,
    mu.c.rate = mu.c.rate
    #tau.Topt = tau.Topt
  )
}
for(i in 1:nc){
  inits[[i]] = GenInits()
}

# Bundle all data in a list for JAGS
jag.data <- list(trait = aedes.rm.pd$Trait, m.N.obs.pd.rm = m.N.obs.pd.rm, temp = aedes.rm.pd$Temp,
                 m.n.climate.rm.pd = m.n.climate.rm.pd, m.climate.pd.rm = m.climate.pd.rm)

Climate.krombeini.pd.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters,
                                 model.file="aedes.by.Climate.txt", n.thin=nt, n.chains=nc, n.burnin=nb,
                                 n.iter=ni, DIC=T, working.directory=getwd())
Climate.krombeini.pd.fit.mcmc <- as.mcmc(Climate.krombeini.pd.fit) ## makes an "mcmc" object
closeAllConnections()

Climate.krombeini.pd.fit.DIC <- Climate.krombeini.pd.fit$BUGSoutput$DIC

# Get posterior samples
m.climate.samp.krombeini <- as.matrix(Climate.krombeini.pd.fit.mcmc)
# Get column indices
m.a_climate.krombeini <- grep("^mu\\.climate\\.a\\[", colnames(m.climate.samp.krombeini))
m.Topt_climate.krombeini <- grep("^mu\\.climate\\.Topt\\[", colnames(m.climate.samp.krombeini))
m.c_climate.krombeini <- grep("^mu\\.climate\\[", colnames(m.climate.samp.krombeini))
# Generate temperature sequence 
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)
# Prepare a list to store data frames
m.climate_dfs.krombeini <- list()
for (j in 1:m.n.climate.rm.pd) {
  # Retrieve parameter vectors for climate j
  m.a_samps <- m.climate.samp.krombeini[, m.a_climate.krombeini[j]]
  m.Topt_samps <- m.climate.samp.krombeini[, m.Topt_climate.krombeini[j]]
  m.c_samps <- m.climate.samp.krombeini[, m.c_climate.krombeini[j]]
  # For each temperature in temp_seq, evaluate across ALL posterior samples
  m.climate_eval.krombeini <- sapply(m.temp_seq, function(temp) {
    m.a_samps * exp(- m.Topt_samps * temp) + m.c_samps
  })
  # Compute pointwise posterior means and HDIs
  m.mean_curve.climate.krombeini <- apply(m.climate_eval.krombeini, 2, mean)
  m.climate_dfs.krombeini[[j]] <- m.mean_curve.climate.krombeini
}
# Combine all climate
m.climate_dfs.krombeini <- unlist(m.climate_dfs.krombeini)
# Optionally, convert to a data frame for easier plotting
fitted_data.Climate.krombeini.pd  <- data.frame(Temperature = rep(m.temp_seq, m.n.climate.rm.pd),
                                                     Fitted = m.climate_dfs.krombeini,
                                                     Climate = factor(rep(c("Subtropical", "Temperate", "Tropical"), each = length(m.temp_seq))))

m.plot_data.Climate.pd <- data.frame(
  Temperature = mosquito.pd$Temp,
  Trait = mosquito.pd$Trait,
  Climate = as.factor(mosquito.pd$Climate),
  Species = as.factor(mosquito.pd$Species)
)
m.plot_data.Climate <- data.frame(
  Temperature = mosquito.data$Temp,
  Trait = mosquito.data$Trait,
  Climate = as.factor(mosquito.data$Climate),
  Species = as.factor(mosquito.data$interactor1)
)

m.species.pd_to_Climate <- unique(m.plot_data.Climate.pd[, c("Species", "Climate")])  # Ensure species matches correct Climate
m.species.pd_to_Climate$Climate <- as.factor(m.species.pd_to_Climate$Climate)
# Merge fitted Climate curves with species-to-Climate mapping
fitted_data.Climate.krombeini.pd <- data.frame(Temperature = rep(m.temp_seq, m.n.climate.rm.pd),
                                                    Fitted = m.climate_dfs.krombeini,
                                                    Climate = factor(rep(c("Subtropical", "Temperate", "Tropical"), each = length(m.temp_seq)))) %>%
  left_join(m.species.pd_to_Climate, by = "Climate")


m.plot_data.Climate.pd$Climate <- factor(m.plot_data.Climate.pd$Climate, levels = c("Temperate", "Subtropical", "Tropical"))
fitted_data.Climate.krombeini.pd$Climate <- factor(fitted_data.Climate.krombeini.pd$Climate, levels = c("Temperate", "Subtropical", "Tropical"))
# Plot observed data and fitted curves with facet wrap by species
p.Climate.krombeini.pd <- ggplot() +
  geom_point(data = m.plot_data.Climate.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.Climate.krombeini.pd, aes(x = Temperature, y = Fitted, color = Climate), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label Size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )  
p.Climate.krombeini.pd

# Merge the observed and predicted data for comparison
merged_data.Climate.krombeini.pd <- fitted_data.Climate.krombeini.pd %>%
  left_join(m.plot_data.Climate, by = c("Temperature", "Climate", "Species"))


# Calculate RMSE for each Climate
m.mse_by_Climate.krombeini.pd <- merged_data.Climate.krombeini.pd |>
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.Climate.krombeini.pd <- merged_data.Climate.krombeini.pd |>
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE)))

# save(m.mse_by_Climate.krombeini.pd,
#      m.overall_rmse.Climate.krombeini.pd,
#      fitted_data.Climate.krombeini.pd,
#      merged_data.Climate.krombeini.pd,
#      Climate.krombeini.pd.fit.DIC,
#      Climate.krombeini.pd.fit,
#      p.Climate.krombeini.pd,
#      file="krombeini.pd.Climate.RData")

save(Climate.krombeini.pd.fit,
     m.climate.samp.krombeini,
     fitted_data.Climate.krombeini.pd,
     p.Climate.krombeini.pd,
     m.mse_by_Climate.krombeini.pd,
     m.overall_rmse.Climate.krombeini.pd,
     Climate.krombeini.pd.fit.DIC,
     file="exp.krombeini.pd.Climate.RData")



######################################### By Climate- Out of sample prediction Notoscriptus #########################################
aedes.rm.pd <- mosquito.pd[mosquito.pd$Species != "Aedes notoscriptus", ]
m.climate.pd.rm <- as.factor(aedes.rm.pd$Climate)
m.n.climate.rm.pd <- length(unique(aedes.rm.pd$Climate))
m.N.obs.pd.rm <- length(aedes.rm.pd$Trait)
unique(aedes.rm.pd$interactor1)
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait

  # # Climate-specific
  # mu.c ~ dexp(10)             # Mean for climate effects c
  # mu.Topt ~ dnorm(30, 0.1)             # Mean for climate effects Topt
  # tau.Topt ~ dexp(0.01)             # Tau for climate effect Topt
  # mu.a ~ dexp(0.01)            # Mean for climate effects Topt
  mu.c ~ dgamma(10, 1)
  mu.Topt ~ dexp(0.001)
 mu.a ~ dexp(10)
  mu.c.rate ~ dgamma(2, 1)

  for (j in 1:m.n.climate.rm.pd) {
    # mu.climate[j] ~ dexp(mu.c)       # Climate-specific low point
    # mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific temp
    # mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific temp
    climate.tau[j] ~ dexp(tau1)
    mu.climate[j] ~ dgamma(mu.c, mu.c.rate)
    mu.climate.Topt[j] ~ dexp(mu.Topt)
    mu.climate.a[j] ~ dexp(mu.a)
  }

  # Likelihood
  for (i in 1:m.N.obs.pd.rm) {
    mu[i] <- mu.climate.a[m.climate.pd.rm[i]] * exp(- mu.climate.Topt[m.climate.pd.rm[i]] * temp[i]) + mu.climate[m.climate.pd.rm[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[m.climate.pd.rm[i]]) T(0, )
  }}
", file = "aedes.by.climate.txt")


# Parameters to monitor
parameters <- c("mu.a", "mu.c", "mu.climate", "climate.tau",
                "tau1", "mu.Topt", "mu.c.rate", "mu.climate.Topt", "mu.climate.a")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 40000 # number of iterations in each chain
nb <- 20000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Initial Values
inits = vector('list', nc)
GenInits = function() {
  # mu.a <- rexp(1, 0.01)
  # mu.c <- rexp(1, 10)
  tau1 <- rexp(1, 0.01)
  # mu.Topt <- rnorm(1, 30, 1/0.1)
  mu.a <- rexp(1, 10)
  mu.Topt <- rexp(1, 0.001)
  mu.c <- rgamma(1, 10, 1)
  mu.c.rate <- rgamma(1, 2, 1)
  #tau.Topt <- rexp(1, 0.01)
  list(
    mu.a = mu.a,
    mu.c = mu.c,
    tau1 = tau1,
    mu.Topt = mu.Topt,
    mu.c.rate = mu.c.rate
    #tau.Topt = tau.Topt
  )
}
for(i in 1:nc){
  inits[[i]] = GenInits()
}

# Bundle all data in a list for JAGS
jag.data <- list(trait = aedes.rm.pd$Trait, m.N.obs.pd.rm = m.N.obs.pd.rm, temp = aedes.rm.pd$Temp,
                 m.n.climate.rm.pd = m.n.climate.rm.pd, m.climate.pd.rm = m.climate.pd.rm)

Climate.notoscriptus.pd.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters,
                                    model.file="aedes.by.Climate.txt", n.thin=nt, n.chains=nc, n.burnin=nb,
                                    n.iter=ni, DIC=T, working.directory=getwd())
Climate.notoscriptus.pd.fit.mcmc <- as.mcmc(Climate.notoscriptus.pd.fit) ## makes an "mcmc" object
closeAllConnections()

Climate.notoscriptus.pd.fit.DIC <- Climate.notoscriptus.pd.fit$BUGSoutput$DIC

# Get posterior samples
m.climate.samp.notoscriptus <- as.matrix(Climate.notoscriptus.pd.fit.mcmc)
# Get column indices
m.a_climate.notoscriptus <- grep("^mu\\.climate\\.a\\[", colnames(m.climate.samp.notoscriptus))
m.Topt_climate.notoscriptus <- grep("^mu\\.climate\\.Topt\\[", colnames(m.climate.samp.notoscriptus))
m.c_climate.notoscriptus <- grep("^mu\\.climate\\[", colnames(m.climate.samp.notoscriptus))
# Generate temperature sequence 
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)
# Prepare a list to store data frames
m.climate_dfs.notoscriptus <- list()
for (j in 1:m.n.climate.rm.pd) {
  # Retrieve parameter vectors for climate j
  m.a_samps <- m.climate.samp.notoscriptus[, m.a_climate.notoscriptus[j]]
  m.Topt_samps <- m.climate.samp.notoscriptus[, m.Topt_climate.notoscriptus[j]]
  m.c_samps <- m.climate.samp.notoscriptus[, m.c_climate.notoscriptus[j]]
  # For each temperature in temp_seq, evaluate across ALL posterior samples
  m.climate_eval.notoscriptus <- sapply(m.temp_seq, function(temp) {
    m.a_samps * exp(- m.Topt_samps * temp) + m.c_samps
  })
  # Compute pointwise posterior means and HDIs
  m.mean_curve.climate.notoscriptus <- apply(m.climate_eval.notoscriptus, 2, mean)
  m.climate_dfs.notoscriptus[[j]] <- m.mean_curve.climate.notoscriptus
}
# Combine all climate
m.climate_dfs.notoscriptus <- unlist(m.climate_dfs.notoscriptus)
# Optionally, convert to a data frame for easier plotting
fitted_data.Climate.notoscriptus.pd  <- data.frame(Temperature = rep(m.temp_seq, m.n.climate.rm.pd),
                                                     Fitted = m.climate_dfs.notoscriptus,
                                                     Climate = factor(rep(c("Subtropical", "Temperate", "Tropical"), each = length(m.temp_seq))))

m.plot_data.Climate.pd <- data.frame(
  Temperature = mosquito.pd$Temp,
  Trait = mosquito.pd$Trait,
  Climate = as.factor(mosquito.pd$Climate),
  Species = as.factor(mosquito.pd$Species)
)
m.plot_data.Climate <- data.frame(
  Temperature = mosquito.data$Temp,
  Trait = mosquito.data$Trait,
  Climate = as.factor(mosquito.data$Climate),
  Species = as.factor(mosquito.data$interactor1)
)

m.species.pd_to_Climate <- unique(m.plot_data.Climate.pd[, c("Species", "Climate")])  # Ensure species matches correct Climate
m.species.pd_to_Climate$Climate <- as.factor(m.species.pd_to_Climate$Climate)
# Merge fitted Climate curves with species-to-Climate mapping
fitted_data.Climate.notoscriptus.pd <- data.frame(Temperature = rep(m.temp_seq, m.n.climate.rm.pd),
                                                    Fitted = m.climate_dfs.notoscriptus,
                                                    Climate = factor(rep(c("Subtropical", "Temperate", "Tropical"), each = length(m.temp_seq)))) %>%
  left_join(m.species.pd_to_Climate, by = "Climate")


m.plot_data.Climate.pd$Climate <- factor(m.plot_data.Climate.pd$Climate, levels = c("Temperate", "Subtropical", "Tropical"))
fitted_data.Climate.notoscriptus.pd$Climate <- factor(fitted_data.Climate.notoscriptus.pd$Climate, levels = c("Temperate", "Subtropical", "Tropical"))
# Plot observed data and fitted curves with facet wrap by species
p.Climate.notoscriptus.pd <- ggplot() +
  geom_point(data = m.plot_data.Climate.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.Climate.notoscriptus.pd, aes(x = Temperature, y = Fitted, color = Climate), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label Size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )  
p.Climate.notoscriptus.pd

# Merge the observed and predicted data for comparison
merged_data.Climate.notoscriptus.pd <- fitted_data.Climate.notoscriptus.pd %>%
  left_join(m.plot_data.Climate, by = c("Temperature", "Climate", "Species"))


# Calculate RMSE for each Climate
m.mse_by_Climate.notoscriptus.pd <- merged_data.Climate.notoscriptus.pd |>
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.Climate.notoscriptus.pd <- merged_data.Climate.notoscriptus.pd |>
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE)))

# save(m.mse_by_Climate.notoscriptus.pd,
#      m.overall_rmse.Climate.notoscriptus.pd,
#      fitted_data.Climate.notoscriptus.pd,
#      merged_data.Climate.notoscriptus.pd,
#      Climate.notoscriptus.pd.fit.DIC,
#      Climate.notoscriptus.pd.fit,
#      p.Climate.notoscriptus.pd,
#      file="notoscriptus.pd.Climate.RData")


save(Climate.notoscriptus.pd.fit,
     m.climate.samp.notoscriptus,
     fitted_data.Climate.notoscriptus.pd,
     p.Climate.notoscriptus.pd,
     m.mse_by_Climate.notoscriptus.pd,
     m.overall_rmse.Climate.notoscriptus.pd,
     Climate.notoscriptus.pd.fit.DIC,
     file="exp.notoscriptus.pd.Climate.RData")





######################################### By Wing Size Category- Out of sample prediction Aegypti #########################################
aedes.rm.pd <- mosquito.data[mosquito.data$interactor1 != "Aedes aegypti", ]
m.size.pd.rm <- as.factor(aedes.rm.pd$Wing.Size)
m.n.size.pd.rm <- length(unique(aedes.rm.pd$Wing.Size))
m.N.obs.pd.rm <- length(aedes.rm.pd$Trait)
unique(aedes.rm.pd$interactor1)
sink("aedes.by.size.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait

  # # Size-specific
  # mu.c ~ dexp(10)             # Mean for size effects c
  # mu.Topt ~ dnorm(30, 0.1)             # Mean for size effects Topt
  # tau.Topt ~ dexp(0.01)             # Tau for size effect Topt
  # mu.a ~ dexp(0.01)            # Mean for size effects Topt
  mu.c ~ dgamma(10, 1)
  mu.Topt ~ dexp(0.001)
 mu.a ~ dexp(10)
  mu.c.rate ~ dgamma(2, 1)

  for (j in 1:m.n.size.pd.rm) {
    # mu.size[j] ~ dexp(mu.c)       # Size-specific low point
    # mu.size.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Size-specific temp
    # mu.size.a[j] ~ dexp(mu.a)    # Size-specific temp
    size.tau[j] ~ dexp(tau1)
    mu.size[j] ~ dgamma(mu.c, mu.c.rate)
    mu.size.Topt[j] ~ dexp(mu.Topt)
    mu.size.a[j] ~ dexp(mu.a)
  }

  # Likelihood
  for (i in 1:m.N.obs.pd.rm) {
    mu[i] <- mu.size.a[m.size.pd.rm[i]] * exp(- mu.size.Topt[m.size.pd.rm[i]] * temp[i]) + mu.size[m.size.pd.rm[i]]
    trait[i] ~ dnorm(mu[i], size.tau[m.size.pd.rm[i]]) T(0, )
  }}
", file = "aedes.by.size.txt")


# Parameters to monitor
parameters <- c("mu.a", "mu.c", "mu.size", "size.tau",
                "tau1", "mu.Topt", "mu.c.rate", "mu.size.Topt", "mu.size.a")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 40000 # number of iterations in each chain
nb <- 20000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Initial Values
inits = vector('list', nc)
GenInits = function() {
  # mu.a <- rexp(1, 0.01)
  # mu.c <- rexp(1, 10)
  tau1 <- rexp(1, 0.01)
  # mu.Topt <- rnorm(1, 30, 1/0.1)
  mu.a <- rexp(1, 10)
  mu.Topt <- rexp(1, 0.001)
  mu.c <- rgamma(1, 10, 1)
  mu.c.rate <- rgamma(1, 2, 1)
  #tau.Topt <- rexp(1, 0.01)
  list(
    mu.a = mu.a,
    mu.c = mu.c,
    tau1 = tau1,
    mu.Topt = mu.Topt,
    mu.c.rate = mu.c.rate
    #tau.Topt = tau.Topt
  )
}
for(i in 1:nc){
  inits[[i]] = GenInits()
}

# Bundle all data in a list for JAGS
jag.data <- list(trait = aedes.rm.pd$Trait, m.N.obs.pd.rm = m.N.obs.pd.rm, temp = aedes.rm.pd$Temp,
                 m.n.size.pd.rm = m.n.size.pd.rm, m.size.pd.rm = m.size.pd.rm)

Size.aegypti.pd.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters,
                            model.file="aedes.by.size.txt", n.thin=nt, n.chains=nc, n.burnin=nb,
                            n.iter=ni, DIC=T, working.directory=getwd())
Size.aegypti.pd.fit.mcmc <- as.mcmc(Size.aegypti.pd.fit) ## makes an "mcmc" object
closeAllConnections()

Size.aegypti.pd.fit.DIC <- Size.aegypti.pd.fit$BUGSoutput$DIC


# Get posterior samples
m.size.samp.aegypti <- as.matrix(Size.aegypti.pd.fit.mcmc)
# Get column indices
m.a_size.aegypti <- grep("^mu\\.size\\.a\\[", colnames(m.size.samp.aegypti))
m.Topt_size.aegypti <- grep("^mu\\.size\\.Topt\\[", colnames(m.size.samp.aegypti))
m.c_size.aegypti <- grep("^mu\\.size\\[", colnames(m.size.samp.aegypti))
# Generate temperature sequence 
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)
# Prepare a list to store data frames
m.size_dfs.aegypti <- list()
for (j in 1:m.n.size.pd.rm) {
  # Retrieve parameter vectors for size j
  m.a_samps <- m.size.samp.aegypti[, m.a_size.aegypti[j]]
  m.Topt_samps <- m.size.samp.aegypti[, m.Topt_size.aegypti[j]]
  m.c_samps <- m.size.samp.aegypti[, m.c_size.aegypti[j]]
  # For each temperature in temp_seq, evaluate across ALL posterior samples
  m.size_eval.aegypti <- sapply(m.temp_seq, function(temp) {
    m.a_samps * exp(- m.Topt_samps * temp) + m.c_samps
  })
  # Compute pointwise posterior means and HDIs
  m.mean_curve.size.aegypti <- apply(m.size_eval.aegypti, 2, mean)
  m.size_dfs.aegypti[[j]] <- m.mean_curve.size.aegypti
}
# Combine all size
m.size_dfs.aegypti <- unlist(m.size_dfs.aegypti)
# Optionally, convert to a data frame for easier plotting
fitted_data.Size.aegypti.pd  <- data.frame(Temperature = rep(m.temp_seq, m.n.size.pd.rm),
                                              Fitted = m.size_dfs.aegypti,
                                              Size = factor(rep(c("Large", "Medium", "Small"), each = length(m.temp_seq))))


m.plot_data.Size.pd <- data.frame(
  Temperature = mosquito.pd$Temp,
  Trait = mosquito.pd$Trait,
  Size = as.factor(mosquito.pd$Wing.Size),
  Species = as.factor(mosquito.pd$Species)
)

m.plot_data.Size <- data.frame(
  Temperature = mosquito.data$Temp,
  Trait = mosquito.data$Trait,
  Size = as.factor(mosquito.data$Wing.Size),
  Species = as.factor(mosquito.data$interactor1)
)

m.species.pd_to_Size <- unique(m.plot_data.Size.pd[, c("Species", "Size")])  # Ensure species matches correct Size
m.species.pd_to_Size$Size <- as.factor(m.species.pd_to_Size$Size)
# Merge fitted Size curves with species-to-Size mapping
fitted_data.Size.aegypti.pd <- data.frame(Temperature = rep(m.temp_seq, m.n.size.pd.rm),
                                          Fitted = m.size_dfs.aegypti,
                                          Size = factor(rep(c("Large", "Medium", "Small"), each = length(m.temp_seq)))) %>%
  left_join(m.species.pd_to_Size, by = "Size")


m.plot_data.Size.pd$Size <- factor(m.plot_data.Size.pd$Size, levels = c("Small", "Medium", "Large"))
fitted_data.Size.aegypti.pd$Size <- factor(fitted_data.Size.aegypti.pd$Size, levels = c("Small", "Medium", "Large"))
# Plot observed data and fitted curves with facet wrap by species
p.Size.aegypti.pd <- ggplot() +
  geom_point(data = m.plot_data.Size.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.Size.aegypti.pd, aes(x = Temperature, y = Fitted, color = Size), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label Size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )  
p.Size.aegypti.pd


# Merge the observed and predicted data for comparison
merged_data.Size.aegypti.pd <- fitted_data.Size.aegypti.pd %>%
  left_join(m.plot_data.Size, by = c("Temperature", "Size", "Species"))


# Calculate RMSE for each species
mse_by_Size.aegypti.pd <- merged_data.Size.aegypti.pd %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.Size.aegypti.pd <- merged_data.Size.aegypti.pd |>
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE)))

# save(mse_by_Size.aegypti.pd,
#      m.overall_rmse.Size.aegypti.pd,
#      fitted_data.Size.aegypti.pd,
#      merged_data.Size.aegypti.pd,
#      Size.aegypti.pd.fit.DIC,
#      Size.aegypti.pd.fit,
#      p.Size.aegypti.pd,
#      file="aegypti.pd.Size.RData")

save(Size.aegypti.pd.fit,
     m.size.samp.aegypti,
     fitted_data.Size.aegypti.pd,
     p.Size.aegypti.pd,
     mse_by_Size.aegypti.pd,
     m.overall_rmse.Size.aegypti.pd,
     Size.aegypti.pd.fit.DIC,
     file="exp.aegypti.pd.size.RData")





######################################### By Wing Size Category- Out of sample prediction Camptorhynchus #########################################
aedes.rm.pd <- mosquito.pd[mosquito.pd$Species != "Aedes camptorhynchus", ]
m.size.pd.rm <- as.factor(aedes.rm.pd$Wing.Size)
m.n.size.pd.rm <- length(unique(aedes.rm.pd$Wing.Size))
m.N.obs.pd.rm <- length(aedes.rm.pd$Trait)
unique(aedes.rm.pd$interactor1)
sink("aedes.by.size.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait

  # # Size-specific
  # mu.c ~ dexp(10)             # Mean for size effects c
  # mu.Topt ~ dnorm(30, 0.1)             # Mean for size effects Topt
  # tau.Topt ~ dexp(0.01)             # Tau for size effect Topt
  # mu.a ~ dexp(0.01)            # Mean for size effects Topt
  mu.c ~ dgamma(10, 1)
  mu.Topt ~ dexp(0.001)
 mu.a ~ dexp(10)
  mu.c.rate ~ dgamma(2, 1)

  for (j in 1:m.n.size.pd.rm) {
    # mu.size[j] ~ dexp(mu.c)       # Size-specific low point
    # mu.size.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Size-specific temp
    # mu.size.a[j] ~ dexp(mu.a)    # Size-specific temp
    size.tau[j] ~ dexp(tau1)
    mu.size[j] ~ dgamma(mu.c, mu.c.rate)
    mu.size.Topt[j] ~ dexp(mu.Topt)
    mu.size.a[j] ~ dexp(mu.a)
  }

  # Likelihood
  for (i in 1:m.N.obs.pd.rm) {
    mu[i] <- mu.size.a[m.size.pd.rm[i]] * exp(- mu.size.Topt[m.size.pd.rm[i]] * temp[i]) + mu.size[m.size.pd.rm[i]]
    trait[i] ~ dnorm(mu[i], size.tau[m.size.pd.rm[i]]) T(0, )
  }}
", file = "aedes.by.size.txt")



# Parameters to monitor
parameters <- c("mu.a", "mu.c", "mu.size", "size.tau",
                "tau1", "mu.Topt", "mu.c.rate", "mu.size.Topt", "mu.size.a")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 40000 # number of iterations in each chain
nb <- 20000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Initial Values
inits = vector('list', nc)
GenInits = function() {
  # mu.a <- rexp(1, 0.01)
  # mu.c <- rexp(1, 10)
  tau1 <- rexp(1, 0.01)
  # mu.Topt <- rnorm(1, 30, 1/0.1)
  mu.a <- rexp(1, 10)
  mu.Topt <- rexp(1, 0.001)
  mu.c <- rgamma(1, 10, 1)
  mu.c.rate <- rgamma(1, 2, 1)
  #tau.Topt <- rexp(1, 0.01)
  list(
    mu.a = mu.a,
    mu.c = mu.c,
    tau1 = tau1,
    mu.Topt = mu.Topt,
    mu.c.rate = mu.c.rate
    #tau.Topt = tau.Topt
  )
}
for(i in 1:nc){
  inits[[i]] = GenInits()
}

# Bundle all data in a list for JAGS
jag.data <- list(trait = aedes.rm.pd$Trait, m.N.obs.pd.rm = m.N.obs.pd.rm, temp = aedes.rm.pd$Temp,
                 m.n.size.pd.rm = m.n.size.pd.rm, m.size.pd.rm = m.size.pd.rm)

Size.camptorhynchus.pd.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters,
                                   model.file="aedes.by.Size.txt", n.thin=nt, n.chains=nc, n.burnin=nb,
                                   n.iter=ni, DIC=T, working.directory=getwd())
Size.camptorhynchus.pd.fit.mcmc <- as.mcmc(Size.camptorhynchus.pd.fit) ## makes an "mcmc" object
closeAllConnections()

Size.camptorhynchus.pd.fit.DIC <- Size.camptorhynchus.pd.fit$BUGSoutput$DIC


# Get posterior samples
m.size.samp.camptorhynchus <- as.matrix(Size.camptorhynchus.pd.fit.mcmc)
# Get column indices
m.a_size.camptorhynchus <- grep("^mu\\.size\\.a\\[", colnames(m.size.samp.camptorhynchus))
m.Topt_size.camptorhynchus <- grep("^mu\\.size\\.Topt\\[", colnames(m.size.samp.camptorhynchus))
m.c_size.camptorhynchus <- grep("^mu\\.size\\[", colnames(m.size.samp.camptorhynchus))
# Generate temperature sequence 
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)
# Prepare a list to store data frames
m.size_dfs.camptorhynchus <- list()
for (j in 1:m.n.size.pd.rm) {
  # Retrieve parameter vectors for size j
  m.a_samps <- m.size.samp.camptorhynchus[, m.a_size.camptorhynchus[j]]
  m.Topt_samps <- m.size.samp.camptorhynchus[, m.Topt_size.camptorhynchus[j]]
  m.c_samps <- m.size.samp.camptorhynchus[, m.c_size.camptorhynchus[j]]
  # For each temperature in temp_seq, evaluate across ALL posterior samples
  m.size_eval.camptorhynchus <- sapply(m.temp_seq, function(temp) {
    m.a_samps * exp(- m.Topt_samps * temp) + m.c_samps
  })
  # Compute pointwise posterior means and HDIs
  m.mean_curve.size.camptorhynchus <- apply(m.size_eval.camptorhynchus, 2, mean)
  m.size_dfs.camptorhynchus[[j]] <- m.mean_curve.size.camptorhynchus
}
# Combine all size
m.size_dfs.camptorhynchus <- unlist(m.size_dfs.camptorhynchus)
# Optionally, convert to a data frame for easier plotting
fitted_data.Size.camptorhynchus.pd  <- data.frame(Temperature = rep(m.temp_seq, m.n.size.pd.rm),
                                           Fitted = m.size_dfs.camptorhynchus,
                                           Size = factor(rep(c("Large", "Medium", "Small"), each = length(m.temp_seq))))


m.plot_data.Size.pd <- data.frame(
  Temperature = mosquito.pd$Temp,
  Trait = mosquito.pd$Trait,
  Size = as.factor(mosquito.pd$Wing.Size),
  Species = as.factor(mosquito.pd$Species)
)

m.plot_data.Size <- data.frame(
  Temperature = mosquito.data$Temp,
  Trait = mosquito.data$Trait,
  Size = as.factor(mosquito.data$Wing.Size),
  Species = as.factor(mosquito.data$interactor1)
)

m.species.pd_to_Size <- unique(m.plot_data.Size.pd[, c("Species", "Size")])  # Ensure species matches correct Size
m.species.pd_to_Size$Size <- as.factor(m.species.pd_to_Size$Size)
# Merge fitted Size curves with species-to-Size mapping
fitted_data.Size.camptorhynchus.pd <- data.frame(Temperature = rep(m.temp_seq, m.n.size.pd.rm),
                                          Fitted = m.size_dfs.camptorhynchus,
                                          Size = factor(rep(c("Large", "Medium", "Small"), each = length(m.temp_seq)))) %>%
  left_join(m.species.pd_to_Size, by = "Size")


m.plot_data.Size.pd$Size <- factor(m.plot_data.Size.pd$Size, levels = c("Small", "Medium", "Large"))
fitted_data.Size.camptorhynchus.pd$Size <- factor(fitted_data.Size.camptorhynchus.pd$Size, levels = c("Small", "Medium", "Large"))
# Plot observed data and fitted curves with facet wrap by species
p.Size.camptorhynchus.pd <- ggplot() +
  geom_point(data = m.plot_data.Size.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.Size.camptorhynchus.pd, aes(x = Temperature, y = Fitted, color = Size), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label Size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )  
p.Size.camptorhynchus.pd


# Merge the observed and predicted data for comparison
merged_data.Size.camptorhynchus.pd <- fitted_data.Size.camptorhynchus.pd %>%
  left_join(m.plot_data.Size, by = c("Temperature", "Size", "Species"))


# Calculate RMSE for each species
mse_by_Size.camptorhynchus.pd <- merged_data.Size.camptorhynchus.pd %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.Size.camptorhynchus.pd <- merged_data.Size.camptorhynchus.pd |>
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE)))

# save(mse_by_Size.camptorhynchus.pd,
#      m.overall_rmse.Size.camptorhynchus.pd,
#      fitted_data.Size.camptorhynchus.pd,
#      merged_data.Size.camptorhynchus.pd,
#      Size.camptorhynchus.pd.fit.DIC,
#      Size.camptorhynchus.pd.fit,
#      p.Size.camptorhynchus.pd,
#      file="camptorhynchus.pd.Size.RData")

save(Size.camptorhynchus.pd.fit,
     m.size.samp.camptorhynchus,
     fitted_data.Size.camptorhynchus.pd,
     p.Size.camptorhynchus.pd,
     mse_by_Size.camptorhynchus.pd,
     m.overall_rmse.Size.camptorhynchus.pd,
     Size.camptorhynchus.pd.fit.DIC,
     file="exp.camptorhynchus.pd.size.RData")



######################################### By Wing Size Category- Out of sample prediction Japonicus Japonicus #########################################
aedes.rm.pd <- mosquito.pd[mosquito.pd$Species != "Aedes japonicus japonicus", ]
m.size.pd.rm <- as.factor(aedes.rm.pd$Wing.Size)
m.n.size.pd.rm <- length(unique(aedes.rm.pd$Wing.Size))
m.N.obs.pd.rm <- length(aedes.rm.pd$Trait)
unique(aedes.rm.pd$interactor1)
sink("aedes.by.size.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait

  # # Size-specific
  # mu.c ~ dexp(10)             # Mean for size effects c
  # mu.Topt ~ dnorm(30, 0.1)             # Mean for size effects Topt
  # tau.Topt ~ dexp(0.01)             # Tau for size effect Topt
  # mu.a ~ dexp(0.01)            # Mean for size effects Topt
  mu.c ~ dgamma(10, 1)
  mu.Topt ~ dexp(0.001)
 mu.a ~ dexp(10)
  mu.c.rate ~ dgamma(2, 1)

  for (j in 1:m.n.size.pd.rm) {
    # mu.size[j] ~ dexp(mu.c)       # Size-specific low point
    # mu.size.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Size-specific temp
    # mu.size.a[j] ~ dexp(mu.a)    # Size-specific temp
    size.tau[j] ~ dexp(tau1)
    mu.size[j] ~ dgamma(mu.c, mu.c.rate)
    mu.size.Topt[j] ~ dexp(mu.Topt)
    mu.size.a[j] ~ dexp(mu.a)
  }

  # Likelihood
  for (i in 1:m.N.obs.pd.rm) {
    mu[i] <- mu.size.a[m.size.pd.rm[i]] * exp(- mu.size.Topt[m.size.pd.rm[i]] * temp[i]) + mu.size[m.size.pd.rm[i]]
    trait[i] ~ dnorm(mu[i], size.tau[m.size.pd.rm[i]]) T(0, )
  }}
", file = "aedes.by.size.txt")

# Parameters to monitor
parameters <- c("mu.a", "mu.c", "mu.size", "size.tau",
                "tau1", "mu.Topt", "mu.c.rate", "mu.size.Topt", "mu.size.a")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 40000 # number of iterations in each chain
nb <- 20000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Initial Values
inits = vector('list', nc)
GenInits = function() {
  # mu.a <- rexp(1, 0.01)
  # mu.c <- rexp(1, 10)
  tau1 <- rexp(1, 0.01)
  # mu.Topt <- rnorm(1, 30, 1/0.1)
  mu.a <- rexp(1, 10)
  mu.Topt <- rexp(1, 0.001)
  mu.c <- rgamma(1, 10, 1)
  mu.c.rate <- rgamma(1, 2, 1)
  #tau.Topt <- rexp(1, 0.01)
  list(
    mu.a = mu.a,
    mu.c = mu.c,
    tau1 = tau1,
    mu.Topt = mu.Topt,
    mu.c.rate = mu.c.rate
    #tau.Topt = tau.Topt
  )
}
for(i in 1:nc){
  inits[[i]] = GenInits()
}

# Bundle all data in a list for JAGS
jag.data <- list(trait = aedes.rm.pd$Trait, m.N.obs.pd.rm = m.N.obs.pd.rm, temp = aedes.rm.pd$Temp,
                 m.n.size.pd.rm = m.n.size.pd.rm, m.size.pd.rm = m.size.pd.rm)

Size.jj.pd.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters,
                       model.file="aedes.by.Size.txt", n.thin=nt, n.chains=nc, n.burnin=nb,
                       n.iter=ni, DIC=T, working.directory=getwd())
Size.jj.pd.fit.mcmc <- as.mcmc(Size.jj.pd.fit) ## makes an "mcmc" object
closeAllConnections()

Size.jj.pd.fit.DIC <- Size.jj.pd.fit$BUGSoutput$DIC


# Get posterior samples
m.size.samp.jj <- as.matrix(Size.jj.pd.fit.mcmc)
# Get column indices
m.a_size.jj <- grep("^mu\\.size\\.a\\[", colnames(m.size.samp.jj))
m.Topt_size.jj <- grep("^mu\\.size\\.Topt\\[", colnames(m.size.samp.jj))
m.c_size.jj <- grep("^mu\\.size\\[", colnames(m.size.samp.jj))
# Generate temperature sequence 
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)
# Prepare a list to store data frames
m.size_dfs.jj <- list()
for (j in 1:m.n.size.pd.rm) {
  # Retrieve parameter vectors for size j
  m.a_samps <- m.size.samp.jj[, m.a_size.jj[j]]
  m.Topt_samps <- m.size.samp.jj[, m.Topt_size.jj[j]]
  m.c_samps <- m.size.samp.jj[, m.c_size.jj[j]]
  # For each temperature in temp_seq, evaluate across ALL posterior samples
  m.size_eval.jj <- sapply(m.temp_seq, function(temp) {
    m.a_samps * exp(- m.Topt_samps * temp) + m.c_samps
  })
  # Compute pointwise posterior means and HDIs
  m.mean_curve.size.jj <- apply(m.size_eval.jj, 2, mean)
  m.size_dfs.jj[[j]] <- m.mean_curve.size.jj
}
# Combine all size
m.size_dfs.jj <- unlist(m.size_dfs.jj)
# Optionally, convert to a data frame for easier plotting
fitted_data.Size.jj.pd  <- data.frame(Temperature = rep(m.temp_seq, m.n.size.pd.rm),
                                                  Fitted = m.size_dfs.jj,
                                                  Size = factor(rep(c("Large", "Medium", "Small"), each = length(m.temp_seq))))


m.plot_data.Size.pd <- data.frame(
  Temperature = mosquito.pd$Temp,
  Trait = mosquito.pd$Trait,
  Size = as.factor(mosquito.pd$Wing.Size),
  Species = as.factor(mosquito.pd$Species)
)

m.plot_data.Size <- data.frame(
  Temperature = mosquito.data$Temp,
  Trait = mosquito.data$Trait,
  Size = as.factor(mosquito.data$Wing.Size),
  Species = as.factor(mosquito.data$interactor1)
)

m.species.pd_to_Size <- unique(m.plot_data.Size.pd[, c("Species", "Size")])  # Ensure species matches correct Size
m.species.pd_to_Size$Size <- as.factor(m.species.pd_to_Size$Size)
# Merge fitted Size curves with species-to-Size mapping
fitted_data.Size.jj.pd <- data.frame(Temperature = rep(m.temp_seq, m.n.size.pd.rm),
                                                 Fitted = m.size_dfs.jj,
                                                 Size = factor(rep(c("Large", "Medium", "Small"), each = length(m.temp_seq)))) %>%
  left_join(m.species.pd_to_Size, by = "Size")

m.plot_data.Size.pd$Size <- factor(m.plot_data.Size.pd$Size, levels = c("Small", "Medium", "Large"))
fitted_data.Size.jj.pd$Size <- factor(fitted_data.Size.jj.pd$Size, levels = c("Small", "Medium", "Large"))
# Plot observed data and fitted curves with facet wrap by species
p.Size.jj.pd <- ggplot() +
  geom_point(data = m.plot_data.Size.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.Size.jj.pd, aes(x = Temperature, y = Fitted, color = Size), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label Size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )  
p.Size.jj.pd


# Merge the observed and predicted data for comparison
merged_data.Size.jj.pd <- fitted_data.Size.jj.pd %>%
  left_join(m.plot_data.Size, by = c("Temperature", "Size", "Species"))


# Calculate RMSE for each species
mse_by_Size.jj.pd <- merged_data.Size.jj.pd %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.Size.jj.pd <- merged_data.Size.jj.pd |>
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE)))

# save(mse_by_Size.jj.pd,
#      m.overall_rmse.Size.jj.pd,
#      fitted_data.Size.jj.pd,
#      merged_data.Size.jj.pd,
#      Size.jj.pd.fit.DIC,
#      Size.jj.pd.fit,
#      p.Size.jj.pd,
#      file="jj.pd.Size.RData")

save(Size.jj.pd.fit,
     m.size.samp.jj,
     fitted_data.Size.jj.pd,
     p.Size.jj.pd,
     mse_by_Size.jj.pd,
     m.overall_rmse.Size.jj.pd,
     Size.jj.pd.fit.DIC,
     file="exp.jj.pd.size.RData")


######################################### By Wing Size Category- Out of sample prediction Krombeini #########################################
aedes.rm.pd <- mosquito.pd[mosquito.pd$Species != "Aedes krombeini", ]
m.size.pd.rm <- as.factor(aedes.rm.pd$Wing.Size)
m.n.size.pd.rm <- length(unique(aedes.rm.pd$Wing.Size))
m.N.obs.pd.rm <- length(aedes.rm.pd$Trait)
unique(aedes.rm.pd$interactor1)
sink("aedes.by.size.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait

  # # Size-specific
  # mu.c ~ dexp(10)             # Mean for size effects c
  # mu.Topt ~ dnorm(30, 0.1)             # Mean for size effects Topt
  # tau.Topt ~ dexp(0.01)             # Tau for size effect Topt
  # mu.a ~ dexp(0.01)            # Mean for size effects Topt
  mu.c ~ dgamma(10, 1)
  mu.Topt ~ dexp(0.001)
 mu.a ~ dexp(10)
  mu.c.rate ~ dgamma(2, 1)

  for (j in 1:m.n.size.pd.rm) {
    # mu.size[j] ~ dexp(mu.c)       # Size-specific low point
    # mu.size.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Size-specific temp
    # mu.size.a[j] ~ dexp(mu.a)    # Size-specific temp
    size.tau[j] ~ dexp(tau1)
    mu.size[j] ~ dgamma(mu.c, mu.c.rate)
    mu.size.Topt[j] ~ dexp(mu.Topt)
    mu.size.a[j] ~ dexp(mu.a)
  }

  # Likelihood
  for (i in 1:m.N.obs.pd.rm) {
    mu[i] <- mu.size.a[m.size.pd.rm[i]] * exp(- mu.size.Topt[m.size.pd.rm[i]] * temp[i]) + mu.size[m.size.pd.rm[i]]
    trait[i] ~ dnorm(mu[i], size.tau[m.size.pd.rm[i]]) T(0, )
  }}
", file = "aedes.by.size.txt")


# Parameters to monitor
parameters <- c("mu.a", "mu.c", "mu.size", "size.tau",
                "tau1", "mu.Topt", "mu.c.rate", "mu.size.Topt", "mu.size.a")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 40000 # number of iterations in each chain
nb <- 20000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Initial Values
inits = vector('list', nc)
GenInits = function() {
  # mu.a <- rexp(1, 0.01)
  # mu.c <- rexp(1, 10)
  tau1 <- rexp(1, 0.01)
  # mu.Topt <- rnorm(1, 30, 1/0.1)
  mu.a <- rexp(1, 10)
  mu.Topt <- rexp(1, 0.001)
  mu.c <- rgamma(1, 10, 1)
  mu.c.rate <- rgamma(1, 2, 1)
  #tau.Topt <- rexp(1, 0.01)
  list(
    mu.a = mu.a,
    mu.c = mu.c,
    tau1 = tau1,
    mu.Topt = mu.Topt,
    mu.c.rate = mu.c.rate
    #tau.Topt = tau.Topt
  )
}
for(i in 1:nc){
  inits[[i]] = GenInits()
}

# Bundle all data in a list for JAGS
jag.data <- list(trait = aedes.rm.pd$Trait, m.N.obs.pd.rm = m.N.obs.pd.rm, temp = aedes.rm.pd$Temp,
                 m.n.size.pd.rm = m.n.size.pd.rm, m.size.pd.rm = m.size.pd.rm)

Size.krombeini.pd.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters,
                              model.file="aedes.by.Size.txt", n.thin=nt, n.chains=nc, n.burnin=nb,
                              n.iter=ni, DIC=T, working.directory=getwd())
Size.krombeini.pd.fit.mcmc <- as.mcmc(Size.krombeini.pd.fit) ## makes an "mcmc" object
closeAllConnections()

Size.krombeini.pd.fit.DIC <- Size.krombeini.pd.fit$BUGSoutput$DIC


# Get posterior samples
m.size.samp.krombeini <- as.matrix(Size.krombeini.pd.fit.mcmc)
# Get column indices
m.a_size.krombeini <- grep("^mu\\.size\\.a\\[", colnames(m.size.samp.krombeini))
m.Topt_size.krombeini <- grep("^mu\\.size\\.Topt\\[", colnames(m.size.samp.krombeini))
m.c_size.krombeini <- grep("^mu\\.size\\[", colnames(m.size.samp.krombeini))
# Generate temperature sequence 
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)
# Prepare a list to store data frames
m.size_dfs.krombeini <- list()
for (j in 1:m.n.size.pd.rm) {
  # Retrieve parameter vectors for size j
  m.a_samps <- m.size.samp.krombeini[, m.a_size.krombeini[j]]
  m.Topt_samps <- m.size.samp.krombeini[, m.Topt_size.krombeini[j]]
  m.c_samps <- m.size.samp.krombeini[, m.c_size.krombeini[j]]
  # For each temperature in temp_seq, evaluate across ALL posterior samples
  m.size_eval.krombeini <- sapply(m.temp_seq, function(temp) {
    m.a_samps * exp(- m.Topt_samps * temp) + m.c_samps
  })
  # Compute pointwise posterior means and HDIs
  m.mean_curve.size.krombeini <- apply(m.size_eval.krombeini, 2, mean)
  m.size_dfs.krombeini[[j]] <- m.mean_curve.size.krombeini
}
# Combine all size
m.size_dfs.krombeini <- unlist(m.size_dfs.krombeini)
# Optionally, convert to a data frame for easier plotting
fitted_data.Size.krombeini.pd  <- data.frame(Temperature = rep(m.temp_seq, m.n.size.pd.rm),
                                                  Fitted = m.size_dfs.krombeini,
                                                  Size = factor(rep(c("Large", "Medium", "Small"), each = length(m.temp_seq))))


m.plot_data.Size.pd <- data.frame(
  Temperature = mosquito.pd$Temp,
  Trait = mosquito.pd$Trait,
  Size = as.factor(mosquito.pd$Wing.Size),
  Species = as.factor(mosquito.pd$Species)
)

m.plot_data.Size <- data.frame(
  Temperature = mosquito.data$Temp,
  Trait = mosquito.data$Trait,
  Size = as.factor(mosquito.data$Wing.Size),
  Species = as.factor(mosquito.data$interactor1)
)

m.species.pd_to_Size <- unique(m.plot_data.Size.pd[, c("Species", "Size")])  # Ensure species matches correct Size
m.species.pd_to_Size$Size <- as.factor(m.species.pd_to_Size$Size)
# Merge fitted Size curves with species-to-Size mapping
fitted_data.Size.krombeini.pd <- data.frame(Temperature = rep(m.temp_seq, m.n.size.pd.rm),
                                                 Fitted = m.size_dfs.krombeini,
                                                 Size = factor(rep(c("Large", "Medium", "Small"), each = length(m.temp_seq)))) %>%
  left_join(m.species.pd_to_Size, by = "Size")

m.plot_data.Size.pd$Size <- factor(m.plot_data.Size.pd$Size, levels = c("Small", "Medium", "Large"))
fitted_data.Size.krombeini.pd$Size <- factor(fitted_data.Size.krombeini.pd$Size, levels = c("Small", "Medium", "Large"))
# Plot observed data and fitted curves with facet wrap by species
p.Size.krombeini.pd <- ggplot() +
  geom_point(data = m.plot_data.Size.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.Size.krombeini.pd, aes(x = Temperature, y = Fitted, color = Size), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label Size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )  
p.Size.krombeini.pd


# Merge the observed and predicted data for comparison
merged_data.Size.krombeini.pd <- fitted_data.Size.krombeini.pd %>%
  left_join(m.plot_data.Size, by = c("Temperature", "Size", "Species"))


# Calculate RMSE for each species
mse_by_Size.krombeini.pd <- merged_data.Size.krombeini.pd %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.Size.krombeini.pd <- merged_data.Size.krombeini.pd |>
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE)))

# save(mse_by_Size.krombeini.pd,
#      m.overall_rmse.Size.krombeini.pd,
#      fitted_data.Size.krombeini.pd,
#      merged_data.Size.krombeini.pd,
#      Size.krombeini.pd.fit.DIC,
#      Size.krombeini.pd.fit,
#      p.Size.krombeini.pd,
#      file="krombeini.pd.Size.RData")


save(Size.krombeini.pd.fit,
     m.size.samp.krombeini,
     fitted_data.Size.krombeini.pd,
     p.Size.krombeini.pd,
     mse_by_Size.krombeini.pd,
     m.overall_rmse.Size.krombeini.pd,
     Size.krombeini.pd.fit.DIC,
     file="exp.krombeini.pd.size.RData")


######################################### By Wing Size Category- Out of sample prediction Notoscriptus #########################################
aedes.rm.pd <- mosquito.pd[mosquito.pd$Species != "Aedes notoscriptus", ]
m.size.pd.rm <- as.factor(aedes.rm.pd$Wing.Size)
m.n.size.pd.rm <- length(unique(aedes.rm.pd$Wing.Size))
m.N.obs.pd.rm <- length(aedes.rm.pd$Trait)
unique(aedes.rm.pd$interactor1)
sink("aedes.by.size.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait

  # # Size-specific
  # mu.c ~ dexp(10)             # Mean for size effects c
  # mu.Topt ~ dnorm(30, 0.1)             # Mean for size effects Topt
  # tau.Topt ~ dexp(0.01)             # Tau for size effect Topt
  # mu.a ~ dexp(0.01)            # Mean for size effects Topt
  mu.c ~ dgamma(10, 1)
  mu.Topt ~ dexp(0.001)
 mu.a ~ dexp(10)
  mu.c.rate ~ dgamma(2, 1)

  for (j in 1:m.n.size.pd.rm) {
    # mu.size[j] ~ dexp(mu.c)       # Size-specific low point
    # mu.size.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Size-specific temp
    # mu.size.a[j] ~ dexp(mu.a)    # Size-specific temp
    size.tau[j] ~ dexp(tau1)
    mu.size[j] ~ dgamma(mu.c, mu.c.rate)
    mu.size.Topt[j] ~ dexp(mu.Topt)
    mu.size.a[j] ~ dexp(mu.a)
  }

  # Likelihood
  for (i in 1:m.N.obs.pd.rm) {
    mu[i] <- mu.size.a[m.size.pd.rm[i]] * exp(- mu.size.Topt[m.size.pd.rm[i]] * temp[i]) + mu.size[m.size.pd.rm[i]]
    trait[i] ~ dnorm(mu[i], size.tau[m.size.pd.rm[i]]) T(0, )
  }}
", file = "aedes.by.size.txt")



# Parameters to monitor
parameters <- c("mu.a", "mu.c", "mu.size", "size.tau",
                "tau1", "mu.Topt", "mu.c.rate", "mu.size.Topt", "mu.size.a")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 40000 # number of iterations in each chain
nb <- 20000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Initial Values
inits = vector('list', nc)
GenInits = function() {
  # mu.a <- rexp(1, 0.01)
  # mu.c <- rexp(1, 10)
  tau1 <- rexp(1, 0.01)
  # mu.Topt <- rnorm(1, 30, 1/0.1)
  mu.a <- rexp(1, 10)
  mu.Topt <- rexp(1, 0.001)
  mu.c <- rgamma(1, 10, 1)
  mu.c.rate <- rgamma(1, 2, 1)
  #tau.Topt <- rexp(1, 0.01)
  list(
    mu.a = mu.a,
    mu.c = mu.c,
    tau1 = tau1,
    mu.Topt = mu.Topt,
    mu.c.rate = mu.c.rate
    #tau.Topt = tau.Topt
  )
}
for(i in 1:nc){
  inits[[i]] = GenInits()
}

# Bundle all data in a list for JAGS
jag.data <- list(trait = aedes.rm.pd$Trait, m.N.obs.pd.rm = m.N.obs.pd.rm, temp = aedes.rm.pd$Temp,
                 m.n.size.pd.rm = m.n.size.pd.rm, m.size.pd.rm = m.size.pd.rm)

Size.notoscriptus.pd.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters,
                                 model.file="aedes.by.Size.txt", n.thin=nt, n.chains=nc, n.burnin=nb,
                                 n.iter=ni, DIC=T, working.directory=getwd())
Size.notoscriptus.pd.fit.mcmc <- as.mcmc(Size.notoscriptus.pd.fit) ## makes an "mcmc" object
closeAllConnections()

Size.notoscriptus.pd.fit.DIC <- Size.notoscriptus.pd.fit$BUGSoutput$DIC


# Get posterior samples
m.size.samp.notoscriptus <- as.matrix(Size.notoscriptus.pd.fit.mcmc)
# Get column indices
m.a_size.notoscriptus <- grep("^mu\\.size\\.a\\[", colnames(m.size.samp.notoscriptus))
m.Topt_size.notoscriptus <- grep("^mu\\.size\\.Topt\\[", colnames(m.size.samp.notoscriptus))
m.c_size.notoscriptus <- grep("^mu\\.size\\[", colnames(m.size.samp.notoscriptus))
# Generate temperature sequence 
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)
# Prepare a list to store data frames
m.size_dfs.notoscriptus <- list()
for (j in 1:m.n.size.pd.rm) {
  # Retrieve parameter vectors for size j
  m.a_samps <- m.size.samp.notoscriptus[, m.a_size.notoscriptus[j]]
  m.Topt_samps <- m.size.samp.notoscriptus[, m.Topt_size.notoscriptus[j]]
  m.c_samps <- m.size.samp.notoscriptus[, m.c_size.notoscriptus[j]]
  # For each temperature in temp_seq, evaluate across ALL posterior samples
  m.size_eval.notoscriptus <- sapply(m.temp_seq, function(temp) {
    m.a_samps * exp(- m.Topt_samps * temp) + m.c_samps
  })
  # Compute pointwise posterior means and HDIs
  m.mean_curve.size.notoscriptus <- apply(m.size_eval.notoscriptus, 2, mean)
  m.size_dfs.notoscriptus[[j]] <- m.mean_curve.size.notoscriptus
}
# Combine all size
m.size_dfs.notoscriptus <- unlist(m.size_dfs.notoscriptus)
# Optionally, convert to a data frame for easier plotting
fitted_data.Size.notoscriptus.pd  <- data.frame(Temperature = rep(m.temp_seq, m.n.size.pd.rm),
                                                  Fitted = m.size_dfs.notoscriptus,
                                                  Size = factor(rep(c("Large", "Medium", "Small"), each = length(m.temp_seq))))


m.plot_data.Size.pd <- data.frame(
  Temperature = mosquito.pd$Temp,
  Trait = mosquito.pd$Trait,
  Size = as.factor(mosquito.pd$Wing.Size),
  Species = as.factor(mosquito.pd$Species)
)

m.plot_data.Size <- data.frame(
  Temperature = mosquito.data$Temp,
  Trait = mosquito.data$Trait,
  Size = as.factor(mosquito.data$Wing.Size),
  Species = as.factor(mosquito.data$interactor1)
)

m.species.pd_to_Size <- unique(m.plot_data.Size.pd[, c("Species", "Size")])  # Ensure species matches correct Size
m.species.pd_to_Size$Size <- as.factor(m.species.pd_to_Size$Size)
# Merge fitted Size curves with species-to-Size mapping
fitted_data.Size.notoscriptus.pd <- data.frame(Temperature = rep(m.temp_seq, m.n.size.pd.rm),
                                                 Fitted = m.size_dfs.notoscriptus,
                                                 Size = factor(rep(c("Large", "Medium", "Small"), each = length(m.temp_seq)))) %>%
  left_join(m.species.pd_to_Size, by = "Size")

m.plot_data.Size.pd$Size <- factor(m.plot_data.Size.pd$Size, levels = c("Small", "Medium", "Large"))
fitted_data.Size.notoscriptus.pd$Size <- factor(fitted_data.Size.notoscriptus.pd$Size, levels = c("Small", "Medium", "Large"))
# Plot observed data and fitted curves with facet wrap by species
p.Size.notoscriptus.pd <- ggplot() +
  geom_point(data = m.plot_data.Size.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.Size.notoscriptus.pd, aes(x = Temperature, y = Fitted, color = Size), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label Size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )  
p.Size.notoscriptus.pd


# Merge the observed and predicted data for comparison
merged_data.Size.notoscriptus.pd <- fitted_data.Size.notoscriptus.pd %>%
  left_join(m.plot_data.Size, by = c("Temperature", "Size", "Species"))


# Calculate RMSE for each species
mse_by_Size.notoscriptus.pd <- merged_data.Size.notoscriptus.pd %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.Size.notoscriptus.pd <- merged_data.Size.notoscriptus.pd |>
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE)))

# save(mse_by_Size.notoscriptus.pd,
#      m.overall_rmse.Size.notoscriptus.pd,
#      fitted_data.Size.notoscriptus.pd,
#      merged_data.Size.notoscriptus.pd,
#      Size.notoscriptus.pd.fit.DIC,
#      Size.notoscriptus.pd.fit,
#      p.Size.notoscriptus.pd,
#      file="notoscriptus.pd.Size.RData")

save(Size.notoscriptus.pd.fit,
     m.size.samp.notoscriptus,
     fitted_data.Size.notoscriptus.pd,
     p.Size.notoscriptus.pd,
     mse_by_Size.notoscriptus.pd,
     m.overall_rmse.Size.notoscriptus.pd,
     Size.notoscriptus.pd.fit.DIC,
     file="exp.notoscriptus.pd.size.RData")



######################################### By Wing Size and Climate- Out of sample prediction Aegypti #########################################
aedes.rm.pd <- mosquito.pd[mosquito.pd$Species != "Aedes aegypti", ]
m.N.obs.pd.rm <- length(aedes.rm.pd$Trait)
m.climate.pd.rm <- as.factor(aedes.rm.pd$Climate)
m.n.climate.rm.pd <- length(unique(aedes.rm.pd$Climate))
m.size.pd.rm <- as.factor(aedes.rm.pd$Wing.Size)
m.n.size.pd.rm <- length(unique(aedes.rm.pd$Wing.Size))
unique(aedes.rm.pd$interactor1)

sink("aedes.by.climate.size.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait

  # # Size-specific
  # mu.c ~ dexp(10)             # Mean for size effects c
  # # Climate-specific
  # mu.Topt ~ dnorm(30, 0.1)             # Mean for climate effects Topt
  # tau.Topt ~ dexp(0.01)             # Tau for species climate Topt
  # mu.a ~ dexp(0.01)            # Mean for climate effects Topt
  mu.c ~ dgamma(10, 1)
  mu.Topt ~ dexp(0.001)
 mu.a ~ dexp(10)
  mu.c.rate ~ dgamma(2, 1)

  for (j in 1:m.n.climate.rm.pd) {
    # mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific temp
    # mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific temp
    mu.climate.Topt[j] ~ dexp(mu.Topt)
    mu.climate.a[j] ~ dexp(mu.a)
  }
    for (k in 1:m.n.size.pd.rm) {
    # mu.size[k] ~ dexp(mu.c)       # Size-specific low point
    size.tau[k] ~ dexp(tau1)
    mu.size[k] ~ dgamma(mu.c, mu.c.rate)
  }
  # Likelihood
  for (i in 1:m.N.obs.pd.rm) {
    mu[i] <- mu.climate.a[m.climate.pd.rm[i]] * exp(- mu.climate.Topt[m.climate.pd.rm[i]] * temp[i]) + mu.size[m.size.pd.rm[i]]
    trait[i] ~ dnorm(mu[i], size.tau[m.size.pd.rm[i]])  T(0, )
  }}
", file = "aedes.by.climate.size.txt")


# Parameters to monitor
parameters <- c("mu.a", "mu.c", "mu.size", "size.tau",
                "tau1", "mu.Topt", "mu.c.rate", "mu.climate.Topt", "mu.climate.a")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 40000 # number of iterations in each chain
nb <- 20000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Initial Values
inits = vector('list', nc)
GenInits = function() {
  # mu.a <- rexp(1, 0.01)
  # mu.c <- rexp(1, 10)
  tau1 <- rexp(1, 0.01)
  # mu.Topt <- rnorm(1, 30, 1/0.1)
  #tau.Topt <- rexp(1, 0.01)
  mu.a <- rexp(1, 10)
  mu.Topt <- rexp(1, 0.001)
  mu.c <- rgamma(1, 10, 1)
  mu.c.rate <- rgamma(1, 2, 1)
  list(
    mu.a = mu.a,
    mu.c = mu.c,
    tau1 = tau1,
    mu.Topt = mu.Topt,
    mu.c.rate = mu.c.rate
    #tau.Topt = tau.Topt
  )
}
for(i in 1:nc){
  inits[[i]] = GenInits()
}

# Bundle all data in a list for JAGS

jag.data<-list(trait = aedes.rm.pd$Trait, m.N.obs.pd.rm = m.N.obs.pd.rm, temp = aedes.rm.pd$Temp,
               m.n.climate.rm.pd = m.n.climate.rm.pd, m.climate.pd.rm = m.climate.pd.rm, m.n.size.pd.rm = m.n.size.pd.rm, m.size.pd.rm = m.size.pd.rm)

m.climate.pd.Size.aegypti.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters,
                                      model.file="aedes.by.climate.size.txt", n.thin=nt, n.chains=nc, n.burnin=nb,
                                      n.iter=ni, DIC=T, working.directory=getwd())
m.climate.pd.Size.aegypti.fit.mcmc <- as.mcmc(m.climate.pd.Size.aegypti.fit) ## makes an "mcmc" object
closeAllConnections()

m.climate.pd.Size.aegypti.fit.DIC <- m.climate.pd.Size.aegypti.fit$BUGSoutput$DIC

# Get posterior samples
m.climate.size.aegypti.samp <- as.matrix(m.climate.pd.Size.aegypti.fit.mcmc)
m.climate.size_vector.aegypti <- numeric()
# Get column indices
a_climate.size <- grep("^mu\\.climate\\.a\\[", colnames(m.climate.size.aegypti.samp))
Topt_climate.size <- grep("^mu\\.climate\\.Topt\\[", colnames(m.climate.size.aegypti.samp))
c_climate.size <- grep("^mu\\.size\\[", colnames(m.climate.size.aegypti.samp))
# Generate temperature sequence 
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)
# Prepare a list to store data frames
m.climate.size_dfs.aegypti <- list()
for(k in 1:m.n.climate.rm.pd){
  for (j in 1:m.n.size.pd.rm) {
    # Retrieve parameter vectors for species j
    m.a_samps <- m.climate.size.aegypti.samp[, a_climate.size[k]]
    m.Topt_samps <- m.climate.size.aegypti.samp[, Topt_climate.size[k]]
    m.c_samps <- m.climate.size.aegypti.samp[, c_climate.size[j]]
    # For each temperature in temp_seq, evaluate across ALL posterior samples
    m.climate.size_eval.aegypti <- sapply(m.temp_seq, function(temp) {
      m.a_samps * exp(- m.Topt_samps * temp) + m.c_samps
    })
    # Compute pointwise posterior means and HDIs
    mean_curve.climate.size.aegypti <- colMeans(m.climate.size_eval.aegypti) # length(temp_seq)
    m.climate.size_vector.aegypti <- c(m.climate.size_vector.aegypti, mean_curve.climate.size.aegypti)
  }}


# Data frame with temperature, fitted value, and host group
m.fitted_data.Climate.Size.aegypti.pd <- data.frame(Temperature = rep(m.temp_seq, m.n.size.pd.rm * m.n.climate.rm.pd), 
                                                      Fitted = m.climate.size_vector.aegypti,
                                                      Climate = factor(rep(c("Subtropical", "Temperate", "Tropical"), each = length(m.temp_seq) * m.n.size.pd)),
                                                      Size = factor(rep(rep(c("Large", "Medium", "Small"),
                                                                            each = length(m.temp_seq), m.n.climate.pd))))

m.plot_data.Climate.Size.pd <- data.frame(
  Temperature = mosquito.pd$Temp,
  Trait = mosquito.pd$Trait,
  Climate = as.factor(mosquito.pd$Climate),
  Species = as.factor(mosquito.pd$Species),
  Size = as.factor(mosquito.pd$Wing.Size)
)

m.plot_data.Climate.Size <- data.frame(
  Temperature = mosquito.data$Temp,
  Trait = mosquito.data$Trait,
  Climate = as.factor(mosquito.data$Climate),
  Species = as.factor(mosquito.data$interactor1),
  Size = as.factor(mosquito.data$Wing.Size)
)

m.species.pd_to_Climate.Size <- unique(m.plot_data.Climate.Size.pd[, c("Species", "Climate", "Size")])
m.species.pd_to_Climate.Size$Climate <- as.factor(m.species.pd_to_Climate.Size$Climate)
m.species.pd_to_Climate.Size$Size <- as.factor(m.species.pd_to_Climate.Size$Size)

# Merge fitted Climate curves with species-to-Climate mapping and filter for matching species
m.fitted_data.Climate.Size.aegypti.pd <- data.frame(Temperature = rep(m.temp_seq, m.n.size.pd.rm * m.n.climate.rm.pd), 
                                                      Fitted = m.climate.size_vector.aegypti,
                                                      Climate = factor(rep(c("Subtropical", "Temperate", "Tropical"), each = length(m.temp_seq) * m.n.size.pd)),
                                                      Size = factor(rep(rep(c("Large", "Medium", "Small"),
                                                                            each = length(m.temp_seq), m.n.climate.pd)))) |>
  left_join(m.species.pd_to_Climate.Size, by = c("Climate", "Size")) |>
  filter(!is.na(Species))  # Keep only rows where Species is non-missing


# Combine Climate and Size for coloring
m.fitted_data.Climate.Size.aegypti.pd$`Climate, Size` <- with(m.fitted_data.Climate.Size.aegypti.pd, paste(Climate, Size, sep = ", "))
m.fitted_data.Climate.Size.aegypti.pd$`Climate, Size` <- factor(m.fitted_data.Climate.Size.aegypti.pd$`Climate, Size`, levels = c("Temperate, Large",
                                                                                                              "Subtropical, Medium", "Subtropical, Large",
                                                                                                              "Tropical, Small", "Tropical, Medium"))
# Plot observed data and fitted curves with facet wrap by species
p.Climate.Size.aegypti.pd <- ggplot() +
  geom_point(data = m.plot_data.Climate.Size.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.Climate.Size.aegypti.pd, aes(x = Temperature, y = Fitted, color = `Climate, Size`), size = 1) +
  labs(x = "Temperature",
       y = "Development Time in Days",
  ) +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label Size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )  
p.Climate.Size.aegypti.pd

# Merge the observed and predicted data for comparison
m.merged_data.Climate.Size.aegypti.pd <- m.fitted_data.Climate.Size.aegypti.pd %>%
  left_join(m.plot_data.Climate.Size, by = c("Temperature", "Climate", "Species"))


# Calculate RMSE for each species
m.mse_by_Climate.Size.aegypti.pd <- m.merged_data.Climate.Size.aegypti.pd %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.Climate.Size.aegypti.pd <- m.merged_data.Climate.Size.aegypti.pd |>
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE)))

save(m.mse_by_Climate.Size.aegypti.pd,
     m.overall_rmse.Climate.Size.aegypti.pd,
     m.fitted_data.Climate.Size.aegypti.pd,
     m.merged_data.Climate.Size.aegypti.pd,
     m.climate.pd.Size.aegypti.fit.DIC,
     m.climate.pd.Size.aegypti.fit,
     p.Climate.Size.aegypti.pd,
     file="aegypti.Climate.Size.pd.RData")

rm(m.climate.pd.Size.aegypti.fit)
rm(m.climate.pd.Size.aegypti.fit.mcmc)

save(m.climate.pd.Size.aegypti.fit,
     m.climate.size.aegypti.samp,
     m.fitted_data.Climate.Size.aegypti.pd,
     p.Climate.Size.aegypti.pd,
     m.mse_by_Climate.Size.aegypti.pd,
     m.overall_rmse.Climate.Size.aegypti.pd,
     m.climate.pd.Size.aegypti.fit.DIC,
     file="exp.aegypti.pd.size.climate.RData")

######################################### By Wing Size and Climate- Out of sample prediction Krombeini #########################################
aedes.rm.pd <- mosquito.data[mosquito.data$interactor1 != "Aedes krombeini", ]
m.N.obs.pd.rm <- length(aedes.rm.pd$Trait)
m.climate.pd.rm <- as.factor(aedes.rm.pd$Climate)
m.n.climate.rm.pd <- length(unique(aedes.rm.pd$Climate))
m.size.pd.rm <- as.factor(aedes.rm.pd$Wing.Size)
m.n.size.pd.rm <- length(unique(aedes.rm.pd$Wing.Size))
unique(aedes.rm.pd$interactor1)

sink("aedes.by.climate.size.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait

  # # Size-specific
  # mu.c ~ dexp(10)             # Mean for size effects c
  # # Climate-specific
  # mu.Topt ~ dnorm(30, 0.1)             # Mean for climate effects Topt
  # tau.Topt ~ dexp(0.01)             # Tau for species climate Topt
  # mu.a ~ dexp(0.01)            # Mean for climate effects Topt
  mu.c ~ dgamma(10, 1)
  mu.Topt ~ dexp(0.001)
 mu.a ~ dexp(10)
  mu.c.rate ~ dgamma(2, 1)

  for (j in 1:m.n.climate.rm.pd) {
    # mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific temp
    # mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific temp
    mu.climate.Topt[j] ~ dexp(mu.Topt)
    mu.climate.a[j] ~ dexp(mu.a)
  }
    for (k in 1:m.n.size.pd.rm) {
    # mu.size[k] ~ dexp(mu.c)       # Size-specific low point
    size.tau[k] ~ dexp(tau1)
    mu.size[k] ~ dgamma(mu.c, mu.c.rate)
  }
  # Likelihood
  for (i in 1:m.N.obs.pd.rm) {
    mu[i] <- mu.climate.a[m.climate.pd.rm[i]] * exp(- mu.climate.Topt[m.climate.pd.rm[i]] * temp[i]) + mu.size[m.size.pd.rm[i]]
    trait[i] ~ dnorm(mu[i], size.tau[m.size.pd.rm[i]])  T(0, )
  }}
", file = "aedes.by.climate.size.txt")
# Parameters to monitor
parameters <- c("mu.a", "mu.c", "mu.size", "size.tau",
                "tau1", "mu.Topt", "mu.c.rate", "mu.climate.Topt", "mu.climate.a")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 40000 # number of iterations in each chain
nb <- 20000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Initial Values
inits = vector('list', nc)
GenInits = function() {
  # mu.a <- rexp(1, 0.01)
  # mu.c <- rexp(1, 10)
  tau1 <- rexp(1, 0.01)
  # mu.Topt <- rnorm(1, 30, 1/0.1)
  #tau.Topt <- rexp(1, 0.01)
  mu.a <- rexp(1, 10)
  mu.Topt <- rexp(1, 0.001)
  mu.c <- rgamma(1, 10, 1)
  mu.c.rate <- rgamma(1, 2, 1)
  list(
    mu.a = mu.a,
    mu.c = mu.c,
    tau1 = tau1,
    mu.Topt = mu.Topt,
    mu.c.rate = mu.c.rate
    #tau.Topt = tau.Topt
  )
}
for(i in 1:nc){
  inits[[i]] = GenInits()
}

# Bundle all data in a list for JAGS

jag.data<-list(trait = aedes.rm.pd$Trait, m.N.obs.pd.rm = m.N.obs.pd.rm, temp = aedes.rm.pd$Temp,
               m.n.climate.rm.pd = m.n.climate.rm.pd, m.climate.pd.rm = m.climate.pd.rm, m.n.size.pd.rm = m.n.size.pd.rm, m.size.pd.rm = m.size.pd.rm)

m.climate.pd.Size.krombeini.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters,
                                        model.file="aedes.by.Climate.Size.txt", n.thin=nt, n.chains=nc, n.burnin=nb,
                                        n.iter=ni, DIC=T, working.directory=getwd())
m.climate.pd.Size.krombeini.fit.mcmc <- as.mcmc(m.climate.pd.Size.krombeini.fit) ## makes an "mcmc" object
closeAllConnections()

m.climate.pd.Size.krombeini.fit.DIC <- m.climate.pd.Size.krombeini.fit$BUGSoutput$DIC

# Get posterior samples
m.climate.size.krombeini.samp <- as.matrix(m.climate.pd.Size.krombeini.fit.mcmc)
m.climate.size_vector.krombeini <- numeric()
# Get column indices
a_climate.size <- grep("^mu\\.climate\\.a\\[", colnames(m.climate.size.krombeini.samp))
Topt_climate.size <- grep("^mu\\.climate\\.Topt\\[", colnames(m.climate.size.krombeini.samp))
c_climate.size <- grep("^mu\\.size\\[", colnames(m.climate.size.krombeini.samp))
# Generate temperature sequence 
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)
# Prepare a list to store data frames
m.climate.size_dfs.krombeini <- list()
for(k in 1:m.n.climate.rm.pd){
  for (j in 1:m.n.size.pd.rm) {
    # Retrieve parameter vectors for species j
    m.a_samps <- m.climate.size.krombeini.samp[, a_climate.size[k]]
    m.Topt_samps <- m.climate.size.krombeini.samp[, Topt_climate.size[k]]
    m.c_samps <- m.climate.size.krombeini.samp[, c_climate.size[j]]
    # For each temperature in temp_seq, evaluate across ALL posterior samples
    m.climate.size_eval.krombeini <- sapply(m.temp_seq, function(temp) {
      m.a_samps * exp(- m.Topt_samps * temp) + m.c_samps
    })
    # Compute pointwise posterior means and HDIs
    mean_curve.climate.size.krombeini <- colMeans(m.climate.size_eval.krombeini) # length(temp_seq)
    m.climate.size_vector.krombeini <- c(m.climate.size_vector.krombeini, mean_curve.climate.size.krombeini)
  }}


# Data frame with temperature, fitted value, and host group
m.fitted_data.Climate.Size.krombeini.pd <- data.frame(Temperature = rep(m.temp_seq, m.n.size.pd.rm * m.n.climate.rm.pd), 
                                                    Fitted = m.climate.size_vector.krombeini,
                                                    Climate = factor(rep(c("Subtropical", "Temperate", "Tropical"), each = length(m.temp_seq) * m.n.size.pd)),
                                                    Size = factor(rep(rep(c("Large", "Medium", "Small"),
                                                                          each = length(m.temp_seq), m.n.climate.pd))))

m.plot_data.Climate.Size.pd <- data.frame(
  Temperature = mosquito.pd$Temp,
  Trait = mosquito.pd$Trait,
  Climate = as.factor(mosquito.pd$Climate),
  Species = as.factor(mosquito.pd$Species),
  Size = as.factor(mosquito.pd$Wing.Size)
)

m.plot_data.Climate.Size <- data.frame(
  Temperature = mosquito.data$Temp,
  Trait = mosquito.data$Trait,
  Climate = as.factor(mosquito.data$Climate),
  Species = as.factor(mosquito.data$interactor1),
  Size = as.factor(mosquito.data$Wing.Size)
)

m.species.pd_to_Climate.Size <- unique(m.plot_data.Climate.Size.pd[, c("Species", "Climate", "Size")])
m.species.pd_to_Climate.Size$Climate <- as.factor(m.species.pd_to_Climate.Size$Climate)
m.species.pd_to_Climate.Size$Size <- as.factor(m.species.pd_to_Climate.Size$Size)

# Merge fitted Climate curves with species-to-Climate mapping and filter for matching species
m.fitted_data.Climate.Size.krombeini.pd <- data.frame(Temperature = rep(m.temp_seq, m.n.size.pd.rm * m.n.climate.rm.pd), 
                                                    Fitted = m.climate.size_vector.krombeini,
                                                    Climate = factor(rep(c("Subtropical", "Temperate", "Tropical"), each = length(m.temp_seq) * m.n.size.pd)),
                                                    Size = factor(rep(rep(c("Large", "Medium", "Small"),
                                                                          each = length(m.temp_seq), m.n.climate.pd)))) |>
  left_join(m.species.pd_to_Climate.Size, by = c("Climate", "Size")) |>
  filter(!is.na(Species))  # Keep only rows where Species is non-missing

# Combine Climate and Size for coloring
m.fitted_data.Climate.Size.krombeini.pd$`Climate, Size` <- with(m.fitted_data.Climate.Size.krombeini.pd, paste(Climate, Size, sep = ", "))
m.fitted_data.Climate.Size.krombeini.pd$`Climate, Size` <- factor(m.fitted_data.Climate.Size.krombeini.pd$`Climate, Size`, levels = c("Temperate, Large",
                                                                                                                  "Subtropical, Medium", "Subtropical, Large",
                                                                                                                  "Tropical, Small", "Tropical, Medium"))
# Plot observed data and fitted curves with facet wrap by species
p.Climate.Size.krombeini.pd <- ggplot() +
  geom_point(data = m.plot_data.Climate.Size.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.Climate.Size.krombeini.pd, aes(x = Temperature, y = Fitted, color = `Climate, Size`), size = 1) +
  labs(x = "Temperature",
       y = "Development Time in Days"
  ) +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label Size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )  
p.Climate.Size.krombeini.pd

# Merge the observed and predicted data for comparison
m.merged_data.Climate.Size.krombeini.pd <- m.fitted_data.Climate.Size.krombeini.pd %>%
  left_join(m.plot_data.Climate.Size, by = c("Temperature", "Climate", "Species"))


# Calculate RMSE for each species
m.mse_by_Climate.Size.krombeini.pd <- m.merged_data.Climate.Size.krombeini.pd %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.Climate.Size.krombeini.pd <- m.merged_data.Climate.Size.krombeini.pd |>
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE)))

# save(m.mse_by_Climate.Size.krombeini.pd,
#      m.overall_rmse.Climate.Size.krombeini.pd,
#      m.fitted_data.Climate.Size.krombeini.pd,
#      m.merged_data.Climate.Size.krombeini.pd,
#      m.climate.pd.Size.krombeini.fit.DIC,
#      m.climate.pd.Size.krombeini.fit,
#      p.Climate.Size.krombeini.pd,
#      file="krombeini.Climate.Size.pd.RData")

save(m.climate.pd.Size.krombeini.fit,
     m.climate.size.krombeini.samp,
     m.fitted_data.Climate.Size.krombeini.pd,
     p.Climate.Size.krombeini.pd,
     m.mse_by_Climate.Size.krombeini.pd,
     m.overall_rmse.Climate.Size.krombeini.pd,
     m.climate.pd.Size.krombeini.fit.DIC,
     file="exp.krombeini.pd.size.climate.RData")



######################################### RMSE ###################################
# Data frame
table.mosquito.rm.pd <- data.frame(
  Method = factor(c(rep("By Climate (RM)", 6 * 5),
                    rep("By Wing Size (RM)", 6 * 5),
                    rep("By Wing Size and Climate (RM)", 6 * 2))),
  Species = factor(rep(c("Aedes aegypti", "Aedes albopictus",
                        "Aedes camptorhynchus",
                         "Aedes japonicus japonicus", "Aedes krombeini",
                         "Aedes notoscriptus"), times = 12)),
  MSE = c(m.mse_by_Climate.aegypti.pd, m.mse_by_Climate.albopictus.pd,
          m.mse_by_Climate.camptorhynchus.pd,
          m.mse_by_Climate.krombeini.pd, m.mse_by_Climate.notoscriptus.pd,
          mse_by_Size.aegypti.pd,
          mse_by_Size.camptorhynchus.pd, mse_by_Size.jj.pd,
          mse_by_Size.krombeini.pd, mse_by_Size.notoscriptus.pd,
          m.mse_by_Climate.Size.aegypti.pd,
          m.mse_by_Climate.Size.krombeini.pd),
  Total.RMSE = rep(as.numeric(c(m.overall_rmse.Climate.aegypti.pd, m.overall_rmse.Climate.albopictus.pd,
                                m.overall_rmse.Climate.camptorhynchus.pd,
                                m.overall_rmse.Climate.krombeini.pd, m.overall_rmse.Climate.notoscriptus.pd,
                                m.overall_rmse.Size.aegypti.pd,
                                m.overall_rmse.Size.camptorhynchus.pd, m.overall_rmse.Size.jj.pd,
                                m.overall_rmse.Size.krombeini.pd, m.overall_rmse.Size.notoscriptus.pd,
                                m.overall_rmse.Climate.Size.aegypti.pd,
                                m.overall_rmse.Climate.Size.krombeini.pd)), each = 6),
  Removed = c(rep(c("Aedes aegypti", "Aedes albopictus",
                    "Aedes camptorhynchus",
                    "Aedes krombeini",
                    "Aedes notoscriptus"), each = 6),
              rep(c("Aedes aegypti",
                    "Aedes camptorhynchus",
                    "Aedes japonicus japonicus", "Aedes krombeini",
                    "Aedes notoscriptus"), each = 6),
              rep(c("Aedes aegypti",
                    "Aedes krombeini"), each = 6))
)

mosquito.DIC.rm.pd <- data.frame(
  Method = rep(c("By Climate", "By Wing Size",
                 "By Climate and Wing Size"), each = 2),
  Removed = rep(c("Aedes aegypti",
                  "Aedes krombeini"), times = 3),
  DIC = c(Climate.aegypti.pd.fit.DIC,
          Climate.notoscriptus.pd.fit.DIC,
          Size.aegypti.pd.fit.DIC,
          Size.notoscriptus.pd.fit.DIC,
          m.climate.pd.Size.aegypti.fit.DIC,
          m.climate.pd.Size.krombeini.fit.DIC),
  RMSE = as.numeric(c(m.overall_rmse.Climate.aegypti.pd,
                              m.overall_rmse.Climate.krombeini.pd,
                              m.overall_rmse.Size.aegypti.pd,
                              m.overall_rmse.Size.krombeini.pd,
                              m.overall_rmse.Climate.Size.aegypti.pd,
                              m.overall_rmse.Climate.Size.krombeini.pd))
)
mosquito.DIC.rm.pd <- mosquito.DIC.rm.pd %>%
  group_by(Removed) %>%
  mutate(Delta.DIC = DIC - min(DIC)) %>%
  ungroup()

mosquito.DIC.rm.pd |>
  ggplot(aes(x = Removed, y = RMSE, fill = Method))+
  geom_bar(stat = "identity", position = "dodge",
           color = "black")+
  theme_minimal() +
  scale_fill_viridis_d(option = "magma") +
  #coord_cartesian(ylim = c(3.5, 5.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )
ggsave("mosquito.RMSE.RM.png", width = 10, height = 8)

mosquito.DIC.rm.pd |>
  ggplot(aes(x = Removed, y = Delta.DIC, fill = Method))+
  geom_bar(stat = "identity", position = "dodge", color = "black")+
  labs(y = expression(Delta~DIC)) +
  theme_minimal() +
  scale_fill_viridis_d(option = "magma") +
  #coord_cartesian(ylim = c(2000, 6200)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )
ggsave("mosquito.DIC.RM.png", width = 10, height = 8)

save(table.mosquito.rm.pd,
     mosquito.DIC.rm.pd, file = "Dataset.Mosquito.RM.RData")


#########################################
#########################################
#########################################
######################################### By Full Data- Out of sample prediction #########################################
# aedes.rm.pd <- mosquito.pd[mosquito.pd$interactor1 != "Aedes notoscriptus", ]
# m.N.obs.pd.rm <- length(aedes.rm.pd$Trait)
# unique(aedes.rm.pd$interactor1)
# sink("aedes.one.data.txt")
# cat("
# model {
#   # Priors
#    a ~ dexp(0.1)             # Width
#   Topt ~ dnorm(30, 0.001)             # Topt
#   tau1 ~ dexp(0.01)           # Tau for trait
#   c ~ dexp(0.001)
#   
#   # Likelihood
#   for (i in 1:m.N.obs.pd.rm) {
#     mu[i] <- a * (Topt - temp[i])^2 + c
#     trait[i] ~ dnorm(mu[i], tau1) T(0, ) 
#   }}
# ", file = "aedes.one.data.txt")
# 
# 
# 
# # Parameters to monitor
# parameters <- c("a", "c", "tau1", "Topt")
# 
# 
# # MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
# ni <- 40000 # number of iterations in each chain
# nb <- 20000 # number of 'burn in' iterations to discard
# nt <- 8 # thinning rate - jags saves every nt iterations in each chain
# nc <- 5 # number of chains

# Initial Values
# inits = vector('list', nc)
# GenInits = function() {
#   a <- rexp(1, 0.1)
#   c <- rexp(1, 0.001)
#   tau1 <- rexp(1, 0.01)
#   Topt <- rnorm(1, 30, 1/0.001)
#   list(
#     a = a,             
#     c = c,
#     tau1 = tau1,
#     Topt = Topt
#   )
# }
# for(i in 1:nc){
#   inits[[i]] = GenInits()
# }
# 
# # Bundle all data in a list for JAGS
# jag.data <- list(trait = aedes.rm.pd$Trait, m.N.obs.pd.rm = m.N.obs.pd.rm, temp = aedes.rm.pd$Temp)
# 
# one.data.notoscriptus.pd.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
#                                   model.file="aedes.one.data.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
#                                   n.iter=ni, DIC=T, working.directory=getwd())
# one.data.notoscriptus.pd.fit.mcmc <- as.mcmc(one.data.pd.notoscriptus.fit) ## makes an "mcmc" object
# closeAllConnections()
# 
# # Get posterior means
# summary_fit.rm.pd <- summary(one.data.notoscriptus.pd.fit.mcmc)$statistics
# a_mean.rm.pd <- summary_fit.rm.pd["a", "Mean"]
# c_mean.rm.pd <- summary_fit.rm.pd["c", "Mean"]
# means_Topt.rm.pd <- summary_fit.rm.pd["Topt", "Mean"]
# 
# # Generate temperature sequence 
# m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)
# 
# # Get predicted values
# full.data.rm.pd <- function (x){
#   a_mean.rm.pd * (x - means_Topt.rm.pd)^2 + c_mean.rm.pd
# }
# 
# 
# # Optionally, convert to a data frame for easier plotting
# fitted_data.one.data.notoscriptus.pd <- data.frame(Temperature = rep(m.temp_seq, m.n.species.pd), 
#                                                 Fitted = full.data.rm.pd(m.temp_seq),
#                                                 Species = factor(rep(mos.spec, each = length(m.temp_seq))))
# 
# plot_data.one.data.pd <- data.frame(
#   Temperature = mosquito.pd$Temp,
#   Trait = mosquito.pd$Trait,
#   Species = as.factor(mosquito.pd$interactor1)
# )
# 
# # Plot observed data and fitted curves with facet wrap by species
# p.notoscriptus.pd <- ggplot() +
#   geom_point(data = plot_data.one.data.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
#   geom_line(data = fitted_data.one.data.notoscriptus.pd, aes(x = Temperature, y = Fitted), col = "blue", size = 1) +
#   labs(x = "Temperature", y = "Development Time in Days") +
#   theme_minimal() +
#   facet_wrap(~ Species) 
# p.notoscriptus.pd
# 
# # Merge the observed and predicted data for comparison
# merged_data.one.data.notoscriptus.pd <- fitted_data.one.data.notoscriptus.pd %>%
#   left_join(plot_data.one.data, by = c("Temperature", "Species"))
# 
# 
# # Calculate RMSE for each species
# mse_by_species.one.data.notoscriptus.pd <- merged_data.one.data.notoscriptus.pd %>%
#   group_by(Species) %>%
#   summarize(
#     MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
#   )
# 
# 
# 
# save(mse_by_species.one.data.notoscriptus.pd,
#      fitted_data.one.data.notoscriptus.pd,
#      merged_data.one.data.notoscriptus.pd,
#      one.data.notoscriptus.pd.fit.mcmc,
#      one.data.notoscriptus.pd.fit,
#      p.notoscriptus.pd,
#      file="notoscriptus.pd.full.data.RData")
# 
# 
######################################### By Wing Size Category Removed- Out of sample prediction #########################################
# aedes.rm.pd <- mosquito.pd[!(mosquito.pd$interactor1 %in% c("Aedes notoscriptus", "Aedes krombeini", "Aedes japonicus japonicus")), ]
# Size.rm.pd <- as.factor(aedes.rm.pd$Wing.Size)
# n.Size.rm.pd <- length(unique(aedes.rm.pd$Wing.Size))
# N.obs.rm.pd <- length(aedes.rm.pd$Trait)
# unique(aedes.rm.pd$interactor1)
# sink("aedes.by.Size.txt")
# cat("
# model {
#   # Priors
#   tau1 ~ dexp(0.01)           # Tau for trait
#   
#   # Size-specific 
#   mu.c ~ dexp(10)             # Mean for Size effects c
#   mu.Topt ~ dnorm(30, 0.1)             # Mean for Size effects Topt
#   tau.Topt ~ dexp(0.01)             # Tau for Size effect Topt
#   mu.a ~ dexp(0.01)            # Mean for Size effects Topt
# 
#   for (j in 1:n.Size.rm.pd) {
#     mu.Size[j] ~ dexp(mu.c)       # Size-specific low point
#     mu.Size.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Size-specific temp
#     mu.Size.a[j] ~ dexp(mu.a)    # Size-specific temp
#   }
# 
#   # Likelihood
#   for (i in 1:N.obs.rm.pd) {
#     mu[i] <- mu.Size.a[Size.rm.pd[i]] * (mu.Size.Topt[Size.rm.pd[i]] - temp[i])^2 + mu.Size[Size.rm.pd[i]]
#     trait[i] ~ dnorm(mu[i], tau1) T(0, ) 
#   }}
# ", file = "aedes.by.Size.txt")
# 
# 
# # Parameters to monitor
# parameters <- c("mu.a", "mu.c", "mu.Size", 
#                 "tau1", "mu.Topt", "tau.Topt", "mu.Size.Topt", "mu.Size.a")
# 
# 
# # MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
# ni <- 40000 # number of iterations in each chain
# nb <- 20000 # number of 'burn in' iterations to discard
# nt <- 8 # thinning rate - jags saves every nt iterations in each chain
# nc <- 5 # number of chains
# 
# # Initial Values
# inits = vector('list', nc)
# GenInits = function() {
#   mu.a <- rexp(1, 0.01)
#   mu.c <- rexp(1, 10)
#   tau1 <- rexp(1, 0.01)
#   mu.Topt <- rnorm(1, 30, 1/0.1)
#   tau.Topt <- rexp(1, 0.01)
#   list(
#     mu.a = mu.a,             
#     mu.c = mu.c,
#     tau1 = tau1,
#     mu.Topt = mu.Topt,
#     tau.Topt = tau.Topt
#   )
# }
# for(i in 1:nc){
#   inits[[i]] = GenInits()
# }
# 
# # Bundle all data in a list for JAGS
# jag.data <- list(trait = aedes.rm.pd$Trait, N.obs.rm.pd = N.obs.rm.pd, temp = aedes.rm.pd$Temp, 
#                  n.Size.rm.pd = n.Size.rm.pd, Size.rm.pd = Size.rm.pd)
# 
# Size.rm.jj.pd.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
#                        model.file="aedes.by.Size.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
#                        n.iter=ni, DIC=T, working.directory=getwd())
# Size.rm.jj.pd.fit.mcmc <- as.mcmc(Size.rm.jj.pd.fit) ## makes an "mcmc" object
# closeAllConnections()
# 
# # Get posterior means
# summary_fit.Size.rm.pd <- summary(Size.rm.jj.pd.fit.mcmc)$statistics
# a_mean_Size.rm.pd <- summary_fit.Size.rm.pd[grep("^mu.Size.a\\[", rownames(summary_fit.Size.rm.pd)), "Mean"]
# c_mean_Size.rm.pd <- summary_fit.Size.rm.pd[grep("^mu.Size\\[", rownames(summary_fit.Size.rm.pd)), "Mean"]
# Topt_mean_Size.rm.pd <- summary_fit.Size.rm.pd[grep("^mu.Size.Topt\\[", rownames(summary_fit.Size.rm.pd)), "Mean"]
# 
# # Generate temperature sequence 
# m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)
# 
# # Get predicted values
# m.size.pd1.rm <- function (x){
#   a_mean_Size.rm.pd[1] * (x - Topt_mean_Size.rm.pd[1])^2 + c_mean_Size.rm.pd[1]
# }
# m.size.pd2.rm <- function (x){
#   a_mean_Size.rm.pd[2] * (x - Topt_mean_Size.rm.pd[2])^2 + c_mean_Size.rm.pd[2]
# }
# m.size.pd3.rm <- function (x){
#   a_mean_Size.rm.pd[3] * (x - Topt_mean_Size.rm.pd[3])^2 + c_mean_Size.rm.pd[3]
# }
# 
# 
# # Optionally, convert to a data frame for easier plotting
# fitted_data.Size.rm.jj.pd <- data.frame(Temperature = rep(m.temp_seq, n.Size.rm.pd), 
#                                      Fitted = c(m.size.pd1.rm(m.temp_seq),
#                                                 m.size.pd2.rm(m.temp_seq),
#                                                 m.size.pd3.rm(m.temp_seq)),
#                                      Size = factor(rep(c("Large", "Medium", "Small"), each = length(temp_seq))))
# 
# m.plot_data.Size.pd <- data.frame(
#   Temperature = mosquito.pd$Temp,
#   Trait = mosquito.pd$Trait,
#   Size = as.factor(mosquito.pd$Wing.Size),
#   Species = as.factor(mosquito.pd$interactor1)
# )
# 
# m.species.pd_to_Size <- unique(m.plot_data.Size.pd[, c("Species", "Size")])  # Ensure species matches correct Climate
# m.species.pd_to_Size$Size <- as.factor(m.species.pd_to_Size$Size)
# # Merge fitted Climate curves with species-to-Climate mapping
# fitted_data.Size.rm.jj.pd <- data.frame(
#   Temperature = rep(m.temp_seq, n.Size.rm.pd),
#   Fitted = c(m.size.pd1.rm(m.temp_seq),
#              m.size.pd2.rm(m.temp_seq),
#              m.size.pd3.rm(m.temp_seq)),
#   Size = factor(rep(c("Large", "Medium", "Small"), each = length(temp_seq)))
# ) %>%
#   left_join(m.species.pd_to_Size, by = "Size")
# 
# m.plot_data.Size.pd$Size <- factor(m.plot_data.Size.pd$Size, levels = c("Small", "Medium", "Large"))
# fitted_data.Size.rm.jj.pd$Size <- factor(fitted_data.Size.rm.jj.pd$Size, levels = c("Small", "Medium", "Large"))
# # Plot observed data and fitted curves with facet wrap by species
# p.Size.rm.jj.pd <- ggplot() +
#   geom_point(data = m.plot_data.Size.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
#   geom_line(data = fitted_data.Size.rm.jj.pd, aes(x = Temperature, y = Fitted, color = Size), size = 1) +
#   labs(x = "Temperature", y = "Development Time in Days") +
#   theme_minimal() +
#   facet_wrap(~ Species) 
# p.Size.rm.jj.pd
# 
# 
# # Merge the observed and predicted data for comparison
# merged_data.Size.rm.jj.pd <- fitted_data.Size.rm.jj.pd %>%
#   left_join(m.plot_data.Size, by = c("Temperature", "Size", "Species"))
# 
# 
# # Calculate RMSE for each species
# mse_by_species.Size.rm.jj.pd <- merged_data.Size.rm.jj.pd %>%
#   group_by(Species) %>%
#   summarize(
#     MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
#   )
# 
# 
# save(mse_by_species.Size.rm.jj.pd,
#      fitted_data.Size.rm.jj.pd,
#      merged_data.Size.rm.jj.pd,
#      Size.rm.jj.pd.fit.mcmc,
#      Size.rm.jj.pd.fit,
#      p.Size.rm.jj.pd,
#      file="jj.Size.pd.rm.RData")
# 
# 
# 
# 
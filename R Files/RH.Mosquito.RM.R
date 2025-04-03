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

######################################### By Climate- Out of sample prediction Aegypti #########################################
aedes.rm.pd <- mosquito.pd[mosquito.pd$interactor1 != "Aedes aegypti", ]
m.climate.pd.rm <- as.factor(aedes.rm.pd$climate)
m.n.climate.rm.pd <- length(unique(aedes.rm.pd$climate))
m.N.obs.pd.rm <- length(aedes.rm.pd$Trait)
unique(aedes.rm.pd$interactor1)
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

  for (j in 1:m.n.climate.rm.pd) {
    mu.climate[j] ~ dexp(mu.c)       # Climate-specific low point
    mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific temp
    mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific temp
    climate.tau[j] ~ dexp(tau1)
  }

  # Likelihood
  for (i in 1:m.N.obs.pd.rm) {
    mu[i] <- mu.climate.a[m.climate.pd.rm[i]] * (mu.climate.Topt[m.climate.pd.rm[i]] - temp[i])^2 + mu.climate[m.climate.pd.rm[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[m.climate.pd.rm[i]]) T(0, )
  }}
", file = "aedes.by.climate.txt")


# Parameters to monitor
parameters <- c("mu.a", "mu.c", "mu.climate", "climate.tau",
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
jag.data <- list(trait = aedes.rm.pd$Trait, m.N.obs.pd.rm = m.N.obs.pd.rm, temp = aedes.rm.pd$Temp,
                 m.n.climate.rm.pd = m.n.climate.rm.pd, m.climate.pd.rm = m.climate.pd.rm)

climate.aegypti.pd.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters,
                               model.file="aedes.by.climate.txt", n.thin=nt, n.chains=nc, n.burnin=nb,
                               n.iter=ni, DIC=T, working.directory=getwd())
climate.aegypti.pd.fit.mcmc <- as.mcmc(climate.aegypti.pd.fit) ## makes an "mcmc" object
closeAllConnections()

# Get posterior means
summary_fit.rm.climate.pd <- summary(climate.aegypti.pd.fit.mcmc)$statistics
a_mean_climate.rm.pd <- summary_fit.rm.climate.pd[grep("^mu.climate.a\\[", rownames(summary_fit.rm.climate.pd)), "Mean"]
c_mean_climate.rm.pd <- summary_fit.rm.climate.pd[grep("^mu.climate\\[", rownames(summary_fit.rm.climate.pd)), "Mean"]
Topt_mean_climate.rm.pd <- summary_fit.rm.climate.pd[grep("^mu.climate.Topt\\[", rownames(summary_fit.rm.climate.pd)), "Mean"]

# Generate temperature sequence
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)

# Get predicted values
clim1.rm.pd <- function (x){
  a_mean_climate.rm.pd[1] * (x - Topt_mean_climate.rm.pd[1])^2 + c_mean_climate.rm.pd[1]
}
clim2.rm.pd <- function (x){
  a_mean_climate.rm.pd[2] * (x - Topt_mean_climate.rm.pd[2])^2 + c_mean_climate.rm.pd[2]
}
clim3.rm.pd <- function (x){
  a_mean_climate.rm.pd[3] * (x - Topt_mean_climate.rm.pd[3])^2 + c_mean_climate.rm.pd[3]
}


# Optionally, convert to a data frame for easier plotting
fitted_data.climate.aegypti.pd <- data.frame(Temperature = rep(m.temp_seq, m.n.climate.rm),
                                             Fitted = c(clim1.rm.pd(m.temp_seq),
                                                        clim2.rm.pd(m.temp_seq),
                                                        clim3.rm.pd(m.temp_seq)),
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
fitted_data.climate.aegypti.pd <- data.frame(
  Temperature = rep(m.temp_seq, m.n.climate.rm),
  Fitted = c(clim1.rm.pd(m.temp_seq),
             clim2.rm.pd(m.temp_seq),
             clim3.rm.pd(m.temp_seq)),
  Climate = factor(rep(c("subtropical", "temperate", "tropical"), each = length(m.temp_seq)))
) %>%
  left_join(m.species.pd_to_climate, by = "Climate")

m.plot_data.climate.pd$Climate <- factor(m.plot_data.climate.pd$Climate, levels = c("temperate", "subtropical", "tropical"))
fitted_data.climate.aegypti.pd$Climate <- factor(fitted_data.climate.aegypti.pd$Climate, levels = c("temperate", "subtropical", "tropical"))
# Plot observed data and fitted curves with facet wrap by species
p.climate.aegypti.pd <- ggplot() +
  geom_point(data = m.plot_data.climate.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.climate.aegypti.pd, aes(x = Temperature, y = Fitted, color = Climate), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species)
p.climate.aegypti.pd

# Merge the observed and predicted data for comparison
merged_data.climate.aegypti.pd <- fitted_data.climate.aegypti.pd %>%
  left_join(m.plot_data.climate, by = c("Temperature", "Climate", "Species"))


# Calculate RMSE for each climate
m.mse_by_climate.aegypti.pd <- merged_data.climate.aegypti.pd |>
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.climate.aegypti.pd <- merged_data.climate.aegypti.pd |>
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE)))

save(m.mse_by_climate.aegypti.pd,
     m.overall_rmse.climate.aegypti.pd,
     fitted_data.climate.aegypti.pd,
     merged_data.climate.aegypti.pd,
     climate.aegypti.pd.fit$BUGSoutput$DIC,
     p.climate.aegypti.pd,
     file="aegypti.pd.climate.RData")


######################################### By Climate- Out of sample prediction Albopictus #########################################
aedes.rm.pd <- mosquito.pd[mosquito.pd$interactor1 != "Aedes albopictus", ]
m.climate.pd.rm <- as.factor(aedes.rm.pd$climate)
m.n.climate.rm.pd <- length(unique(aedes.rm.pd$climate))
m.N.obs.pd.rm <- length(aedes.rm.pd$Trait)
unique(aedes.rm.pd$interactor1)
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

  for (j in 1:m.n.climate.rm.pd) {
    mu.climate[j] ~ dexp(mu.c)       # Climate-specific low point
    mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific temp
    mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific temp
    climate.tau[j] ~ dexp(tau1)
  }

  # Likelihood
  for (i in 1:m.N.obs.pd.rm) {
    mu[i] <- mu.climate.a[m.climate.pd.rm[i]] * (mu.climate.Topt[m.climate.pd.rm[i]] - temp[i])^2 + mu.climate[m.climate.pd.rm[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[m.climate.pd.rm[i]]) T(0, )
  }}
", file = "aedes.by.climate.txt")

# Parameters to monitor
parameters <- c("mu.a", "mu.c", "mu.climate", "climate.tau",
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
jag.data <- list(trait = aedes.rm.pd$Trait, m.N.obs.pd.rm = m.N.obs.pd.rm, temp = aedes.rm.pd$Temp,
                 m.n.climate.rm.pd = m.n.climate.rm.pd, m.climate.pd.rm = m.climate.pd.rm)

climate.albopictus.pd.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters,
                                  model.file="aedes.by.climate.txt", n.thin=nt, n.chains=nc, n.burnin=nb,
                                  n.iter=ni, DIC=T, working.directory=getwd())
climate.albopictus.pd.fit.mcmc <- as.mcmc(climate.albopictus.pd.fit) ## makes an "mcmc" object
closeAllConnections()

# Get posterior means
summary_fit.rm.climate.pd <- summary(climate.albopictus.pd.fit.mcmc)$statistics
a_mean_climate.rm.pd <- summary_fit.rm.climate.pd[grep("^mu.climate.a\\[", rownames(summary_fit.rm.climate.pd)), "Mean"]
c_mean_climate.rm.pd <- summary_fit.rm.climate.pd[grep("^mu.climate\\[", rownames(summary_fit.rm.climate.pd)), "Mean"]
Topt_mean_climate.rm.pd <- summary_fit.rm.climate.pd[grep("^mu.climate.Topt\\[", rownames(summary_fit.rm.climate.pd)), "Mean"]

# Generate temperature sequence
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)

# Get predicted values
clim1.rm.pd <- function (x){
  a_mean_climate.rm.pd[1] * (x - Topt_mean_climate.rm.pd[1])^2 + c_mean_climate.rm.pd[1]
}
clim2.rm.pd <- function (x){
  a_mean_climate.rm.pd[2] * (x - Topt_mean_climate.rm.pd[2])^2 + c_mean_climate.rm.pd[2]
}
clim3.rm.pd <- function (x){
  a_mean_climate.rm.pd[3] * (x - Topt_mean_climate.rm.pd[3])^2 + c_mean_climate.rm.pd[3]
}


# Optionally, convert to a data frame for easier plotting
fitted_data.climate.albopictus.pd <- data.frame(Temperature = rep(m.temp_seq, m.n.climate.rm),
                                                Fitted = c(clim1.rm.pd(m.temp_seq),
                                                           clim2.rm.pd(m.temp_seq),
                                                           clim3.rm.pd(m.temp_seq)),
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
fitted_data.climate.albopictus.pd <- data.frame(
  Temperature = rep(m.temp_seq, m.n.climate.rm),
  Fitted = c(clim1.rm.pd(m.temp_seq),
             clim2.rm.pd(m.temp_seq),
             clim3.rm.pd(m.temp_seq)),
  Climate = factor(rep(c("subtropical", "temperate", "tropical"), each = length(m.temp_seq)))
) %>%
  left_join(m.species.pd_to_climate, by = "Climate")

m.plot_data.climate.pd$Climate <- factor(m.plot_data.climate.pd$Climate, levels = c("temperate", "subtropical", "tropical"))
fitted_data.climate.albopictus.pd$Climate <- factor(fitted_data.climate.albopictus.pd$Climate, levels = c("temperate", "subtropical", "tropical"))
# Plot observed data and fitted curves with facet wrap by species
p.climate.albopictus.pd <- ggplot() +
  geom_point(data = m.plot_data.climate.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.climate.albopictus.pd, aes(x = Temperature, y = Fitted, color = Climate), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species)
p.climate.albopictus.pd

# Merge the observed and predicted data for comparison
merged_data.climate.albopictus.pd <- fitted_data.climate.albopictus.pd %>%
  left_join(m.plot_data.climate, by = c("Temperature", "Climate", "Species"))


# Calculate RMSE for each climate
m.mse_by_climate.albopictus.pd <- merged_data.climate.albopictus.pd |>
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.climate.albopictus.pd <- merged_data.climate.albopictus.pd |>
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE)))

save(m.mse_by_climate.albopictus.pd,
     m.overall_rmse.climate.albopictus.pd,
     fitted_data.climate.albopictus.pd,
     merged_data.climate.albopictus.pd,
     climate.albopictus.pd.fit$BUGSoutput$DIC,
     p.climate.albopictus.pd,
     file="albopictus.pd.climate.RData")

######################################### By Climate- Out of sample prediction Camptorhynchus #########################################
aedes.rm.pd <- mosquito.pd[mosquito.pd$interactor1 != "Aedes camptorhynchus", ]
m.climate.pd.rm <- as.factor(aedes.rm.pd$climate)
m.n.climate.rm.pd <- length(unique(aedes.rm.pd$climate))
m.N.obs.pd.rm <- length(aedes.rm.pd$Trait)
unique(aedes.rm.pd$interactor1)
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

  for (j in 1:m.n.climate.rm.pd) {
    mu.climate[j] ~ dexp(mu.c)       # Climate-specific low point
    mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific temp
    mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific temp
    climate.tau[j] ~ dexp(tau1)
  }

  # Likelihood
  for (i in 1:m.N.obs.pd.rm) {
    mu[i] <- mu.climate.a[m.climate.pd.rm[i]] * (mu.climate.Topt[m.climate.pd.rm[i]] - temp[i])^2 + mu.climate[m.climate.pd.rm[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[m.climate.pd.rm[i]]) T(0, )
  }}
", file = "aedes.by.climate.txt")


# Parameters to monitor
parameters <- c("mu.a", "mu.c", "mu.climate", "climate.tau",
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
jag.data <- list(trait = aedes.rm.pd$Trait, m.N.obs.pd.rm = m.N.obs.pd.rm, temp = aedes.rm.pd$Temp,
                 m.n.climate.rm.pd = m.n.climate.rm.pd, m.climate.pd.rm = m.climate.pd.rm)

climate.camptorhynchus.pd.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters,
                                      model.file="aedes.by.climate.txt", n.thin=nt, n.chains=nc, n.burnin=nb,
                                      n.iter=ni, DIC=T, working.directory=getwd())
climate.camptorhynchus.pd.fit.mcmc <- as.mcmc(climate.camptorhynchus.pd.fit) ## makes an "mcmc" object
closeAllConnections()

# Get posterior means
summary_fit.rm.climate.pd <- summary(climate.camptorhynchus.pd.fit.mcmc)$statistics
a_mean_climate.rm.pd <- summary_fit.rm.climate.pd[grep("^mu.climate.a\\[", rownames(summary_fit.rm.climate.pd)), "Mean"]
c_mean_climate.rm.pd <- summary_fit.rm.climate.pd[grep("^mu.climate\\[", rownames(summary_fit.rm.climate.pd)), "Mean"]
Topt_mean_climate.rm.pd <- summary_fit.rm.climate.pd[grep("^mu.climate.Topt\\[", rownames(summary_fit.rm.climate.pd)), "Mean"]

# Generate temperature sequence
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)

# Get predicted values
clim1.rm.pd <- function (x){
  a_mean_climate.rm.pd[1] * (x - Topt_mean_climate.rm.pd[1])^2 + c_mean_climate.rm.pd[1]
}
clim2.rm.pd <- function (x){
  a_mean_climate.rm.pd[2] * (x - Topt_mean_climate.rm.pd[2])^2 + c_mean_climate.rm.pd[2]
}
clim3.rm.pd <- function (x){
  a_mean_climate.rm.pd[3] * (x - Topt_mean_climate.rm.pd[3])^2 + c_mean_climate.rm.pd[3]
}


# Optionally, convert to a data frame for easier plotting
fitted_data.climate.camptorhynchus.pd <- data.frame(Temperature = rep(m.temp_seq, m.n.climate.rm),
                                                    Fitted = c(clim1.rm.pd(m.temp_seq),
                                                               clim2.rm.pd(m.temp_seq),
                                                               clim3.rm.pd(m.temp_seq)),
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
fitted_data.climate.camptorhynchus.pd <- data.frame(
  Temperature = rep(m.temp_seq, m.n.climate.rm),
  Fitted = c(clim1.rm.pd(m.temp_seq),
             clim2.rm.pd(m.temp_seq),
             clim3.rm.pd(m.temp_seq)),
  Climate = factor(rep(c("subtropical", "temperate", "tropical"), each = length(m.temp_seq)))
) %>%
  left_join(m.species.pd_to_climate, by = "Climate")

m.plot_data.climate.pd$Climate <- factor(m.plot_data.climate.pd$Climate, levels = c("temperate", "subtropical", "tropical"))
fitted_data.climate.camptorhynchus.pd$Climate <- factor(fitted_data.climate.camptorhynchus.pd$Climate, levels = c("temperate", "subtropical", "tropical"))
# Plot observed data and fitted curves with facet wrap by species
p.climate.camptorhynchus.pd <- ggplot() +
  geom_point(data = m.plot_data.climate.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.climate.camptorhynchus.pd, aes(x = Temperature, y = Fitted, color = Climate), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species)
p.climate.camptorhynchus.pd

# Merge the observed and predicted data for comparison
merged_data.climate.camptorhynchus.pd <- fitted_data.climate.camptorhynchus.pd %>%
  left_join(m.plot_data.climate, by = c("Temperature", "Climate", "Species"))


# Calculate RMSE for each climate
m.mse_by_climate.camptorhynchus.pd <- merged_data.climate.camptorhynchus.pd |>
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.climate.camptorhynchus.pd <- merged_data.climate.camptorhynchus.pd |>
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE)))

save(m.mse_by_climate.camptorhynchus.pd,
     m.overall_rmse.climate.camptorhynchus.pd,
     fitted_data.climate.camptorhynchus.pd,
     merged_data.climate.camptorhynchus.pd,
     climate.camptorhynchus.pd.fit$BUGSoutput$DIC,
     p.climate.camptorhynchus.pd,
     file="camptorhynchus.pd.climate.RData")

######################################### By Climate- Out of sample prediction Krombeini #########################################
aedes.rm.pd <- mosquito.pd[mosquito.pd$interactor1 != "Aedes krombeini", ]
m.climate.pd.rm <- as.factor(aedes.rm.pd$climate)
m.n.climate.rm.pd <- length(unique(aedes.rm.pd$climate))
m.N.obs.pd.rm <- length(aedes.rm.pd$Trait)
unique(aedes.rm.pd$interactor1)
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

  for (j in 1:m.n.climate.rm.pd) {
    mu.climate[j] ~ dexp(mu.c)       # Climate-specific low point
    mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific temp
    mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific temp
    climate.tau[j] ~ dexp(tau1)
  }

  # Likelihood
  for (i in 1:m.N.obs.pd.rm) {
    mu[i] <- mu.climate.a[m.climate.pd.rm[i]] * (mu.climate.Topt[m.climate.pd.rm[i]] - temp[i])^2 + mu.climate[m.climate.pd.rm[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[m.climate.pd.rm[i]]) T(0, )
  }}
", file = "aedes.by.climate.txt")


# Parameters to monitor
parameters <- c("mu.a", "mu.c", "mu.climate", "climate.tau",
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
jag.data <- list(trait = aedes.rm.pd$Trait, m.N.obs.pd.rm = m.N.obs.pd.rm, temp = aedes.rm.pd$Temp,
                 m.n.climate.rm.pd = m.n.climate.rm.pd, m.climate.pd.rm = m.climate.pd.rm)

climate.krombeini.pd.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters,
                                 model.file="aedes.by.climate.txt", n.thin=nt, n.chains=nc, n.burnin=nb,
                                 n.iter=ni, DIC=T, working.directory=getwd())
climate.krombeini.pd.fit.mcmc <- as.mcmc(climate.krombeini.pd.fit) ## makes an "mcmc" object
closeAllConnections()

# Get posterior means
summary_fit.rm.climate.pd <- summary(climate.krombeini.pd.fit.mcmc)$statistics
a_mean_climate.rm.pd <- summary_fit.rm.climate.pd[grep("^mu.climate.a\\[", rownames(summary_fit.rm.climate.pd)), "Mean"]
c_mean_climate.rm.pd <- summary_fit.rm.climate.pd[grep("^mu.climate\\[", rownames(summary_fit.rm.climate.pd)), "Mean"]
Topt_mean_climate.rm.pd <- summary_fit.rm.climate.pd[grep("^mu.climate.Topt\\[", rownames(summary_fit.rm.climate.pd)), "Mean"]

# Generate temperature sequence
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)

# Get predicted values
clim1.rm.pd <- function (x){
  a_mean_climate.rm.pd[1] * (x - Topt_mean_climate.rm.pd[1])^2 + c_mean_climate.rm.pd[1]
}
clim2.rm.pd <- function (x){
  a_mean_climate.rm.pd[2] * (x - Topt_mean_climate.rm.pd[2])^2 + c_mean_climate.rm.pd[2]
}
clim3.rm.pd <- function (x){
  a_mean_climate.rm.pd[3] * (x - Topt_mean_climate.rm.pd[3])^2 + c_mean_climate.rm.pd[3]
}


# Optionally, convert to a data frame for easier plotting
fitted_data.climate.krombeini.pd <- data.frame(Temperature = rep(m.temp_seq, m.n.climate.rm),
                                               Fitted = c(clim1.rm.pd(m.temp_seq),
                                                          clim2.rm.pd(m.temp_seq),
                                                          clim3.rm.pd(m.temp_seq)),
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
fitted_data.climate.krombeini.pd <- data.frame(
  Temperature = rep(m.temp_seq, m.n.climate.rm),
  Fitted = c(clim1.rm.pd(m.temp_seq),
             clim2.rm.pd(m.temp_seq),
             clim3.rm.pd(m.temp_seq)),
  Climate = factor(rep(c("subtropical", "temperate", "tropical"), each = length(m.temp_seq)))
) %>%
  left_join(m.species.pd_to_climate, by = "Climate")

m.plot_data.climate.pd$Climate <- factor(m.plot_data.climate.pd$Climate, levels = c("temperate", "subtropical", "tropical"))
fitted_data.climate.krombeini.pd$Climate <- factor(fitted_data.climate.krombeini.pd$Climate, levels = c("temperate", "subtropical", "tropical"))
# Plot observed data and fitted curves with facet wrap by species
p.climate.krombeini.pd <- ggplot() +
  geom_point(data = m.plot_data.climate.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.climate.krombeini.pd, aes(x = Temperature, y = Fitted, color = Climate), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species)
p.climate.krombeini.pd

# Merge the observed and predicted data for comparison
merged_data.climate.krombeini.pd <- fitted_data.climate.krombeini.pd %>%
  left_join(m.plot_data.climate, by = c("Temperature", "Climate", "Species"))


# Calculate RMSE for each climate
m.mse_by_climate.krombeini.pd <- merged_data.climate.krombeini.pd |>
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.climate.krombeini.pd <- merged_data.climate.krombeini.pd |>
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE)))

save(m.mse_by_climate.krombeini.pd,
     m.overall_rmse.climate.krombeini.pd,
     fitted_data.climate.krombeini.pd,
     merged_data.climate.krombeini.pd,
     climate.krombeini.pd.fit$BUGSoutput$DIC,
     p.climate.krombeini.pd,
     file="krombeini.pd.climate.RData")

######################################### By Climate- Out of sample prediction Notoscriptus #########################################
aedes.rm.pd <- mosquito.pd[mosquito.pd$interactor1 != "Aedes notoscriptus", ]
m.climate.pd.rm <- as.factor(aedes.rm.pd$climate)
m.n.climate.rm.pd <- length(unique(aedes.rm.pd$climate))
m.N.obs.pd.rm <- length(aedes.rm.pd$Trait)
unique(aedes.rm.pd$interactor1)
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

  for (j in 1:m.n.climate.rm.pd) {
    mu.climate[j] ~ dexp(mu.c)       # Climate-specific low point
    mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific temp
    mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific temp
    climate.tau[j] ~ dexp(tau1)
  }

  # Likelihood
  for (i in 1:m.N.obs.pd.rm) {
    mu[i] <- mu.climate.a[m.climate.pd.rm[i]] * (mu.climate.Topt[m.climate.pd.rm[i]] - temp[i])^2 + mu.climate[m.climate.pd.rm[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[m.climate.pd.rm[i]]) T(0, )
  }}
", file = "aedes.by.climate.txt")


# Parameters to monitor
parameters <- c("mu.a", "mu.c", "mu.climate", "climate.tau",
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
jag.data <- list(trait = aedes.rm.pd$Trait, m.N.obs.pd.rm = m.N.obs.pd.rm, temp = aedes.rm.pd$Temp,
                 m.n.climate.rm.pd = m.n.climate.rm.pd, m.climate.pd.rm = m.climate.pd.rm)

climate.notoscriptus.pd.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters,
                                    model.file="aedes.by.climate.txt", n.thin=nt, n.chains=nc, n.burnin=nb,
                                    n.iter=ni, DIC=T, working.directory=getwd())
climate.notoscriptus.pd.fit.mcmc <- as.mcmc(climate.notoscriptus.pd.fit) ## makes an "mcmc" object
closeAllConnections()

# Get posterior means
summary_fit.rm.climate.pd <- summary(climate.notoscriptus.pd.fit.mcmc)$statistics
a_mean_climate.rm.pd <- summary_fit.rm.climate.pd[grep("^mu.climate.a\\[", rownames(summary_fit.rm.climate.pd)), "Mean"]
c_mean_climate.rm.pd <- summary_fit.rm.climate.pd[grep("^mu.climate\\[", rownames(summary_fit.rm.climate.pd)), "Mean"]
Topt_mean_climate.rm.pd <- summary_fit.rm.climate.pd[grep("^mu.climate.Topt\\[", rownames(summary_fit.rm.climate.pd)), "Mean"]

# Generate temperature sequence
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)

# Get predicted values
clim1.rm.pd <- function (x){
  a_mean_climate.rm.pd[1] * (x - Topt_mean_climate.rm.pd[1])^2 + c_mean_climate.rm.pd[1]
}
clim2.rm.pd <- function (x){
  a_mean_climate.rm.pd[2] * (x - Topt_mean_climate.rm.pd[2])^2 + c_mean_climate.rm.pd[2]
}
clim3.rm.pd <- function (x){
  a_mean_climate.rm.pd[3] * (x - Topt_mean_climate.rm.pd[3])^2 + c_mean_climate.rm.pd[3]
}


# Optionally, convert to a data frame for easier plotting
fitted_data.climate.notoscriptus.pd <- data.frame(Temperature = rep(m.temp_seq, m.n.climate.rm),
                                                  Fitted = c(clim1.rm.pd(m.temp_seq),
                                                             clim2.rm.pd(m.temp_seq),
                                                             clim3.rm.pd(m.temp_seq)),
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
fitted_data.climate.notoscriptus.pd <- data.frame(
  Temperature = rep(m.temp_seq, m.n.climate.rm),
  Fitted = c(clim1.rm.pd(m.temp_seq),
             clim2.rm.pd(m.temp_seq),
             clim3.rm.pd(m.temp_seq)),
  Climate = factor(rep(c("subtropical", "temperate", "tropical"), each = length(m.temp_seq)))
) %>%
  left_join(m.species.pd_to_climate, by = "Climate")

m.plot_data.climate.pd$Climate <- factor(m.plot_data.climate.pd$Climate, levels = c("temperate", "subtropical", "tropical"))
fitted_data.climate.notoscriptus.pd$Climate <- factor(fitted_data.climate.notoscriptus.pd$Climate, levels = c("temperate", "subtropical", "tropical"))
# Plot observed data and fitted curves with facet wrap by species
p.climate.notoscriptus.pd <- ggplot() +
  geom_point(data = m.plot_data.climate.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.climate.notoscriptus.pd, aes(x = Temperature, y = Fitted, color = Climate), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species)
p.climate.notoscriptus.pd

# Merge the observed and predicted data for comparison
merged_data.climate.notoscriptus.pd <- fitted_data.climate.notoscriptus.pd %>%
  left_join(m.plot_data.climate, by = c("Temperature", "Climate", "Species"))


# Calculate RMSE for each climate
m.mse_by_climate.notoscriptus.pd <- merged_data.climate.notoscriptus.pd |>
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.climate.notoscriptus.pd <- merged_data.climate.notoscriptus.pd |>
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE)))

save(m.mse_by_climate.notoscriptus.pd,
     m.overall_rmse.climate.notoscriptus.pd,
     fitted_data.climate.notoscriptus.pd,
     merged_data.climate.notoscriptus.pd,
     climate.notoscriptus.pd.fit$BUGSoutput$DIC,
     p.climate.notoscriptus.pd,
     file="notoscriptus.pd.climate.RData")


######################################### By Wing Size Category- Out of sample prediction Aegypti #########################################
aedes.rm.pd <- mosquito.pd[mosquito.pd$interactor1 != "Aedes aegypti", ]
m.size.pd.rm <- as.factor(aedes.rm.pd$wing.size)
m.n.size.pd.rm <- length(unique(aedes.rm.pd$wing.size))
m.N.obs.pd.rm <- length(aedes.rm.pd$Trait)
unique(aedes.rm.pd$interactor1)
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

  for (j in 1:m.n.size.pd.rm) {
    mu.size[j] ~ dexp(mu.c)       # Size-specific low point
    mu.size.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Size-specific temp
    mu.size.a[j] ~ dexp(mu.a)    # Size-specific temp
    size.tau[j] ~ dexp(tau1)
  }

  # Likelihood
  for (i in 1:m.N.obs.pd.rm) {
    mu[i] <- mu.size.a[m.size.pd.rm[i]] * (mu.size.Topt[m.size.pd.rm[i]] - temp[i])^2 + mu.size[m.size.pd.rm[i]]
    trait[i] ~ dnorm(mu[i], size.tau[m.size.pd.rm[i]]) T(0, )
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
jag.data <- list(trait = aedes.rm.pd$Trait, m.N.obs.pd.rm = m.N.obs.pd.rm, temp = aedes.rm.pd$Temp,
                 m.n.size.pd.rm = m.n.size.pd.rm, m.size.pd.rm = m.size.pd.rm)

size.aegypti.pd.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters,
                            model.file="aedes.by.size.txt", n.thin=nt, n.chains=nc, n.burnin=nb,
                            n.iter=ni, DIC=T, working.directory=getwd())
size.aegypti.pd.fit.mcmc <- as.mcmc(size.aegypti.pd.fit) ## makes an "mcmc" object
closeAllConnections()


# Get posterior means
summary_fit.size.rm.pd <- summary(size.aegypti.pd.fit.mcmc)$statistics
a_mean_size.rm.pd <- summary_fit.size.rm.pd[grep("^mu.size.a\\[", rownames(summary_fit.size.rm.pd)), "Mean"]
c_mean_size.rm.pd <- summary_fit.size.rm.pd[grep("^mu.size\\[", rownames(summary_fit.size.rm.pd)), "Mean"]
Topt_mean_size.rm.pd <- summary_fit.size.rm.pd[grep("^mu.size.Topt\\[", rownames(summary_fit.size.rm.pd)), "Mean"]

# Generate temperature sequence
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)

# Get predicted values
size1.rm.pd <- function (x){
  a_mean_size.rm.pd[1] * (x - Topt_mean_size.rm.pd[1])^2 + c_mean_size.rm.pd[1]
}
size2.rm.pd <- function (x){
  a_mean_size.rm.pd[2] * (x - Topt_mean_size.rm.pd[2])^2 + c_mean_size.rm.pd[2]
}
size3.rm.pd <- function (x){
  a_mean_size.rm.pd[3] * (x - Topt_mean_size.rm.pd[3])^2 + c_mean_size.rm.pd[3]
}


# Optionally, convert to a data frame for easier plotting
fitted_data.size.aegypti.pd <- data.frame(Temperature = rep(m.temp_seq, m.n.size.pd.rm),
                                          Fitted = c(size1.rm.pd(m.temp_seq),
                                                     size2.rm.pd(m.temp_seq),
                                                     size3.rm.pd(m.temp_seq)),
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
fitted_data.size.aegypti <- data.frame(
  Temperature = rep(m.temp_seq, m.n.size.pd.rm),
  Fitted = c(size1.rm.pd(m.temp_seq),
             size2.rm.pd(m.temp_seq),
             size3.rm.pd(m.temp_seq)),
  Size = factor(rep(c("large", "medium", "small"), each = length(m.temp_seq)))
) %>%
  left_join(m.species.pd_to_size, by = "Size")

m.plot_data.size.pd$Size <- factor(m.plot_data.size.pd$Size, levels = c("small", "medium", "large"))
fitted_data.size.aegypti.pd$Size <- factor(fitted_data.size.aegypti.pd$Size, levels = c("small", "medium", "large"))
# Plot observed data and fitted curves with facet wrap by species
p.size.aegypti.pd <- ggplot() +
  geom_point(data = m.plot_data.size.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.size.aegypti.pd, aes(x = Temperature, y = Fitted, color = Size), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species)
p.size.aegypti.pd


# Merge the observed and predicted data for comparison
merged_data.size.aegypti.pd <- fitted_data.size.aegypti.pd %>%
  left_join(m.plot_data.size, by = c("Temperature", "Size", "Species"))


# Calculate RMSE for each species
mse_by_size.aegypti.pd <- merged_data.size.aegypti.pd %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.size.aegypti.pd <- merged_data.size.aegypti.pd |>
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE)))

save(mse_by_size.aegypti.pd,
     m.overall_rmse.size.aegypti.pd,
     fitted_data.size.aegypti.pd,
     merged_data.size.aegypti.pd,
     size.aegypti.pd.fit$BUGSoutput$DIC,
     p.size.aegypti.pd,
     file="aegypti.pd.size.RData")


######################################### By Wing Size Category- Out of sample prediction Camptorhynchus #########################################
aedes.rm.pd <- mosquito.pd[mosquito.pd$interactor1 != "Aedes camptorhynchus", ]
m.size.pd.rm <- as.factor(aedes.rm.pd$wing.size)
m.n.size.pd.rm <- length(unique(aedes.rm.pd$wing.size))
m.N.obs.pd.rm <- length(aedes.rm.pd$Trait)
unique(aedes.rm.pd$interactor1)
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

  for (j in 1:m.n.size.pd.rm) {
    mu.size[j] ~ dexp(mu.c)       # Size-specific low point
    mu.size.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Size-specific temp
    mu.size.a[j] ~ dexp(mu.a)    # Size-specific temp
    size.tau[j] ~ dexp(tau1)
  }

  # Likelihood
  for (i in 1:m.N.obs.pd.rm) {
    mu[i] <- mu.size.a[m.size.pd.rm[i]] * (mu.size.Topt[m.size.pd.rm[i]] - temp[i])^2 + mu.size[m.size.pd.rm[i]]
    trait[i] ~ dnorm(mu[i], size.tau[m.size.pd.rm[i]]) T(0, )
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
jag.data <- list(trait = aedes.rm.pd$Trait, m.N.obs.pd.rm = m.N.obs.pd.rm, temp = aedes.rm.pd$Temp,
                 m.n.size.pd.rm = m.n.size.pd.rm, m.size.pd.rm = m.size.pd.rm)

size.camptorhynchus.pd.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters,
                                   model.file="aedes.by.size.txt", n.thin=nt, n.chains=nc, n.burnin=nb,
                                   n.iter=ni, DIC=T, working.directory=getwd())
size.camptorhynchus.pd.fit.mcmc <- as.mcmc(size.camptorhynchus.pd.fit) ## makes an "mcmc" object
closeAllConnections()


# Get posterior means
summary_fit.size.rm.pd <- summary(size.camptorhynchus.pd.fit.mcmc)$statistics
a_mean_size.rm.pd <- summary_fit.size.rm.pd[grep("^mu.size.a\\[", rownames(summary_fit.size.rm.pd)), "Mean"]
c_mean_size.rm.pd <- summary_fit.size.rm.pd[grep("^mu.size\\[", rownames(summary_fit.size.rm.pd)), "Mean"]
Topt_mean_size.rm.pd <- summary_fit.size.rm.pd[grep("^mu.size.Topt\\[", rownames(summary_fit.size.rm.pd)), "Mean"]

# Generate temperature sequence
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)

# Get predicted values
size1.rm.pd <- function (x){
  a_mean_size.rm.pd[1] * (x - Topt_mean_size.rm.pd[1])^2 + c_mean_size.rm.pd[1]
}
size2.rm.pd <- function (x){
  a_mean_size.rm.pd[2] * (x - Topt_mean_size.rm.pd[2])^2 + c_mean_size.rm.pd[2]
}
size3.rm.pd <- function (x){
  a_mean_size.rm.pd[3] * (x - Topt_mean_size.rm.pd[3])^2 + c_mean_size.rm.pd[3]
}


# Optionally, convert to a data frame for easier plotting
fitted_data.size.camptorhynchus.pd <- data.frame(Temperature = rep(m.temp_seq, m.n.size.pd.rm),
                                                 Fitted = c(size1.rm.pd(m.temp_seq),
                                                            size2.rm.pd(m.temp_seq),
                                                            size3.rm.pd(m.temp_seq)),
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
fitted_data.size.camptorhynchus <- data.frame(
  Temperature = rep(m.temp_seq, m.n.size.pd.rm),
  Fitted = c(size1.rm.pd(m.temp_seq),
             size2.rm.pd(m.temp_seq),
             size3.rm.pd(m.temp_seq)),
  Size = factor(rep(c("large", "medium", "small"), each = length(m.temp_seq)))
) %>%
  left_join(m.species.pd_to_size, by = "Size")

m.plot_data.size.pd$Size <- factor(m.plot_data.size.pd$Size, levels = c("small", "medium", "large"))
fitted_data.size.camptorhynchus.pd$Size <- factor(fitted_data.size.camptorhynchus.pd$Size, levels = c("small", "medium", "large"))
# Plot observed data and fitted curves with facet wrap by species
p.size.camptorhynchus.pd <- ggplot() +
  geom_point(data = m.plot_data.size.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.size.camptorhynchus.pd, aes(x = Temperature, y = Fitted, color = Size), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species)
p.size.camptorhynchus.pd


# Merge the observed and predicted data for comparison
merged_data.size.camptorhynchus.pd <- fitted_data.size.camptorhynchus.pd %>%
  left_join(m.plot_data.size, by = c("Temperature", "Size", "Species"))


# Calculate RMSE for each species
mse_by_size.camptorhynchus.pd <- merged_data.size.camptorhynchus.pd %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.size.camptorhynchus.pd <- merged_data.size.camptorhynchus.pd |>
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE)))

save(mse_by_size.camptorhynchus.pd,
     m.overall_rmse.size.camptorhynchus.pd,
     fitted_data.size.camptorhynchus.pd,
     merged_data.size.camptorhynchus.pd,
     size.camptorhynchus.pd.fit$BUGSoutput$DIC,
     p.size.camptorhynchus.pd,
     file="camptorhynchus.pd.size.RData")


######################################### By Wing Size Category- Out of sample prediction Japonicus Japonicus #########################################
aedes.rm.pd <- mosquito.pd[mosquito.pd$interactor1 != "Aedes japonicus japonicus", ]
m.size.pd.rm <- as.factor(aedes.rm.pd$wing.size)
m.n.size.pd.rm <- length(unique(aedes.rm.pd$wing.size))
m.N.obs.pd.rm <- length(aedes.rm.pd$Trait)
unique(aedes.rm.pd$interactor1)
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

  for (j in 1:m.n.size.pd.rm) {
    mu.size[j] ~ dexp(mu.c)       # Size-specific low point
    mu.size.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Size-specific temp
    mu.size.a[j] ~ dexp(mu.a)    # Size-specific temp
    size.tau[j] ~ dexp(tau1)
  }

  # Likelihood
  for (i in 1:m.N.obs.pd.rm) {
    mu[i] <- mu.size.a[m.size.pd.rm[i]] * (mu.size.Topt[m.size.pd.rm[i]] - temp[i])^2 + mu.size[m.size.pd.rm[i]]
    trait[i] ~ dnorm(mu[i], size.tau[m.size.pd.rm[i]]) T(0, )
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
jag.data <- list(trait = aedes.rm.pd$Trait, m.N.obs.pd.rm = m.N.obs.pd.rm, temp = aedes.rm.pd$Temp,
                 m.n.size.pd.rm = m.n.size.pd.rm, m.size.pd.rm = m.size.pd.rm)

size.jj.pd.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters,
                       model.file="aedes.by.size.txt", n.thin=nt, n.chains=nc, n.burnin=nb,
                       n.iter=ni, DIC=T, working.directory=getwd())
size.jj.pd.fit.mcmc <- as.mcmc(size.jj.pd.fit) ## makes an "mcmc" object
closeAllConnections()


# Get posterior means
summary_fit.size.rm.pd <- summary(size.jj.pd.fit.mcmc)$statistics
a_mean_size.rm.pd <- summary_fit.size.rm.pd[grep("^mu.size.a\\[", rownames(summary_fit.size.rm.pd)), "Mean"]
c_mean_size.rm.pd <- summary_fit.size.rm.pd[grep("^mu.size\\[", rownames(summary_fit.size.rm.pd)), "Mean"]
Topt_mean_size.rm.pd <- summary_fit.size.rm.pd[grep("^mu.size.Topt\\[", rownames(summary_fit.size.rm.pd)), "Mean"]

# Generate temperature sequence
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)

# Get predicted values
size1.rm.pd <- function (x){
  a_mean_size.rm.pd[1] * (x - Topt_mean_size.rm.pd[1])^2 + c_mean_size.rm.pd[1]
}
size2.rm.pd <- function (x){
  a_mean_size.rm.pd[2] * (x - Topt_mean_size.rm.pd[2])^2 + c_mean_size.rm.pd[2]
}
size3.rm.pd <- function (x){
  a_mean_size.rm.pd[3] * (x - Topt_mean_size.rm.pd[3])^2 + c_mean_size.rm.pd[3]
}


# Optionally, convert to a data frame for easier plotting
fitted_data.size.jj.pd <- data.frame(Temperature = rep(m.temp_seq, m.n.size.pd.rm),
                                     Fitted = c(size1.rm.pd(m.temp_seq),
                                                size2.rm.pd(m.temp_seq),
                                                size3.rm.pd(m.temp_seq)),
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
fitted_data.size.jj <- data.frame(
  Temperature = rep(m.temp_seq, m.n.size.pd.rm),
  Fitted = c(size1.rm.pd(m.temp_seq),
             size2.rm.pd(m.temp_seq),
             size3.rm.pd(m.temp_seq)),
  Size = factor(rep(c("large", "medium", "small"), each = length(m.temp_seq)))
) %>%
  left_join(m.species.pd_to_size, by = "Size")

m.plot_data.size.pd$Size <- factor(m.plot_data.size.pd$Size, levels = c("small", "medium", "large"))
fitted_data.size.jj.pd$Size <- factor(fitted_data.size.jj.pd$Size, levels = c("small", "medium", "large"))
# Plot observed data and fitted curves with facet wrap by species
p.size.jj.pd <- ggplot() +
  geom_point(data = m.plot_data.size.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.size.jj.pd, aes(x = Temperature, y = Fitted, color = Size), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species)
p.size.jj.pd


# Merge the observed and predicted data for comparison
merged_data.size.jj.pd <- fitted_data.size.jj.pd %>%
  left_join(m.plot_data.size, by = c("Temperature", "Size", "Species"))


# Calculate RMSE for each species
mse_by_size.jj.pd <- merged_data.size.jj.pd %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.size.jj.pd <- merged_data.size.jj.pd |>
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE)))

save(mse_by_size.jj.pd,
     m.overall_rmse.size.jj.pd,
     fitted_data.size.jj.pd,
     merged_data.size.jj.pd,
     size.jj.pd.fit$BUGSoutput$DIC,
     p.size.jj.pd,
     file="jj.pd.size.RData")


######################################### By Wing Size Category- Out of sample prediction Krombeini #########################################
aedes.rm.pd <- mosquito.pd[mosquito.pd$interactor1 != "Aedes krombeini", ]
m.size.pd.rm <- as.factor(aedes.rm.pd$wing.size)
m.n.size.pd.rm <- length(unique(aedes.rm.pd$wing.size))
m.N.obs.pd.rm <- length(aedes.rm.pd$Trait)
unique(aedes.rm.pd$interactor1)
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

  for (j in 1:m.n.size.pd.rm) {
    mu.size[j] ~ dexp(mu.c)       # Size-specific low point
    mu.size.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Size-specific temp
    mu.size.a[j] ~ dexp(mu.a)    # Size-specific temp
    size.tau[j] ~ dexp(tau1)
  }

  # Likelihood
  for (i in 1:m.N.obs.pd.rm) {
    mu[i] <- mu.size.a[m.size.pd.rm[i]] * (mu.size.Topt[m.size.pd.rm[i]] - temp[i])^2 + mu.size[m.size.pd.rm[i]]
    trait[i] ~ dnorm(mu[i], size.tau[m.size.pd.rm[i]]) T(0, )
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
jag.data <- list(trait = aedes.rm.pd$Trait, m.N.obs.pd.rm = m.N.obs.pd.rm, temp = aedes.rm.pd$Temp,
                 m.n.size.pd.rm = m.n.size.pd.rm, m.size.pd.rm = m.size.pd.rm)

size.krombeini.pd.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters,
                              model.file="aedes.by.size.txt", n.thin=nt, n.chains=nc, n.burnin=nb,
                              n.iter=ni, DIC=T, working.directory=getwd())
size.krombeini.pd.fit.mcmc <- as.mcmc(size.krombeini.pd.fit) ## makes an "mcmc" object
closeAllConnections()


# Get posterior means
summary_fit.size.rm.pd <- summary(size.krombeini.pd.fit.mcmc)$statistics
a_mean_size.rm.pd <- summary_fit.size.rm.pd[grep("^mu.size.a\\[", rownames(summary_fit.size.rm.pd)), "Mean"]
c_mean_size.rm.pd <- summary_fit.size.rm.pd[grep("^mu.size\\[", rownames(summary_fit.size.rm.pd)), "Mean"]
Topt_mean_size.rm.pd <- summary_fit.size.rm.pd[grep("^mu.size.Topt\\[", rownames(summary_fit.size.rm.pd)), "Mean"]

# Generate temperature sequence
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)

# Get predicted values
size1.rm.pd <- function (x){
  a_mean_size.rm.pd[1] * (x - Topt_mean_size.rm.pd[1])^2 + c_mean_size.rm.pd[1]
}
size2.rm.pd <- function (x){
  a_mean_size.rm.pd[2] * (x - Topt_mean_size.rm.pd[2])^2 + c_mean_size.rm.pd[2]
}
size3.rm.pd <- function (x){
  a_mean_size.rm.pd[3] * (x - Topt_mean_size.rm.pd[3])^2 + c_mean_size.rm.pd[3]
}


# Optionally, convert to a data frame for easier plotting
fitted_data.size.krombeini.pd <- data.frame(Temperature = rep(m.temp_seq, m.n.size.pd.rm),
                                            Fitted = c(size1.rm.pd(m.temp_seq),
                                                       size2.rm.pd(m.temp_seq),
                                                       size3.rm.pd(m.temp_seq)),
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
fitted_data.size.krombeini <- data.frame(
  Temperature = rep(m.temp_seq, m.n.size.pd.rm),
  Fitted = c(size1.rm.pd(m.temp_seq),
             size2.rm.pd(m.temp_seq),
             size3.rm.pd(m.temp_seq)),
  Size = factor(rep(c("large", "medium", "small"), each = length(m.temp_seq)))
) %>%
  left_join(m.species.pd_to_size, by = "Size")

m.plot_data.size.pd$Size <- factor(m.plot_data.size.pd$Size, levels = c("small", "medium", "large"))
fitted_data.size.krombeini.pd$Size <- factor(fitted_data.size.krombeini.pd$Size, levels = c("small", "medium", "large"))
# Plot observed data and fitted curves with facet wrap by species
p.size.krombeini.pd <- ggplot() +
  geom_point(data = m.plot_data.size.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.size.krombeini.pd, aes(x = Temperature, y = Fitted, color = Size), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species)
p.size.krombeini.pd


# Merge the observed and predicted data for comparison
merged_data.size.krombeini.pd <- fitted_data.size.krombeini.pd %>%
  left_join(m.plot_data.size, by = c("Temperature", "Size", "Species"))


# Calculate RMSE for each species
mse_by_size.krombeini.pd <- merged_data.size.krombeini.pd %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.size.krombeini.pd <- merged_data.size.krombeini.pd |>
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE)))

save(mse_by_size.krombeini.pd,
     m.overall_rmse.size.krombeini.pd,
     fitted_data.size.krombeini.pd,
     merged_data.size.krombeini.pd,
     size.krombeini.pd.fit$BUGSoutput$DIC,
     p.size.krombeini.pd,
     file="krombeini.pd.size.RData")


######################################### By Wing Size Category- Out of sample prediction Notoscriptus #########################################
aedes.rm.pd <- mosquito.pd[mosquito.pd$interactor1 != "Aedes notoscriptus", ]
m.size.pd.rm <- as.factor(aedes.rm.pd$wing.size)
m.n.size.pd.rm <- length(unique(aedes.rm.pd$wing.size))
m.N.obs.pd.rm <- length(aedes.rm.pd$Trait)
unique(aedes.rm.pd$interactor1)
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

  for (j in 1:m.n.size.pd.rm) {
    mu.size[j] ~ dexp(mu.c)       # Size-specific low point
    mu.size.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Size-specific temp
    mu.size.a[j] ~ dexp(mu.a)    # Size-specific temp
    size.tau[j] ~ dexp(tau1)
  }

  # Likelihood
  for (i in 1:m.N.obs.pd.rm) {
    mu[i] <- mu.size.a[m.size.pd.rm[i]] * (mu.size.Topt[m.size.pd.rm[i]] - temp[i])^2 + mu.size[m.size.pd.rm[i]]
    trait[i] ~ dnorm(mu[i], size.tau[m.size.pd.rm[i]]) T(0, )
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
jag.data <- list(trait = aedes.rm.pd$Trait, m.N.obs.pd.rm = m.N.obs.pd.rm, temp = aedes.rm.pd$Temp,
                 m.n.size.pd.rm = m.n.size.pd.rm, m.size.pd.rm = m.size.pd.rm)

size.notoscriptus.pd.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters,
                                 model.file="aedes.by.size.txt", n.thin=nt, n.chains=nc, n.burnin=nb,
                                 n.iter=ni, DIC=T, working.directory=getwd())
size.notoscriptus.pd.fit.mcmc <- as.mcmc(size.notoscriptus.pd.fit) ## makes an "mcmc" object
closeAllConnections()


# Get posterior means
summary_fit.size.rm.pd <- summary(size.notoscriptus.pd.fit.mcmc)$statistics
a_mean_size.rm.pd <- summary_fit.size.rm.pd[grep("^mu.size.a\\[", rownames(summary_fit.size.rm.pd)), "Mean"]
c_mean_size.rm.pd <- summary_fit.size.rm.pd[grep("^mu.size\\[", rownames(summary_fit.size.rm.pd)), "Mean"]
Topt_mean_size.rm.pd <- summary_fit.size.rm.pd[grep("^mu.size.Topt\\[", rownames(summary_fit.size.rm.pd)), "Mean"]

# Generate temperature sequence
m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)

# Get predicted values
size1.rm.pd <- function (x){
  a_mean_size.rm.pd[1] * (x - Topt_mean_size.rm.pd[1])^2 + c_mean_size.rm.pd[1]
}
size2.rm.pd <- function (x){
  a_mean_size.rm.pd[2] * (x - Topt_mean_size.rm.pd[2])^2 + c_mean_size.rm.pd[2]
}
size3.rm.pd <- function (x){
  a_mean_size.rm.pd[3] * (x - Topt_mean_size.rm.pd[3])^2 + c_mean_size.rm.pd[3]
}


# Optionally, convert to a data frame for easier plotting
fitted_data.size.notoscriptus.pd <- data.frame(Temperature = rep(m.temp_seq, m.n.size.pd.rm),
                                               Fitted = c(size1.rm.pd(m.temp_seq),
                                                          size2.rm.pd(m.temp_seq),
                                                          size3.rm.pd(m.temp_seq)),
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
fitted_data.size.notoscriptus <- data.frame(
  Temperature = rep(m.temp_seq, m.n.size.pd.rm),
  Fitted = c(size1.rm.pd(m.temp_seq),
             size2.rm.pd(m.temp_seq),
             size3.rm.pd(m.temp_seq)),
  Size = factor(rep(c("large", "medium", "small"), each = length(m.temp_seq)))
) %>%
  left_join(m.species.pd_to_size, by = "Size")

m.plot_data.size.pd$Size <- factor(m.plot_data.size.pd$Size, levels = c("small", "medium", "large"))
fitted_data.size.notoscriptus.pd$Size <- factor(fitted_data.size.notoscriptus.pd$Size, levels = c("small", "medium", "large"))
# Plot observed data and fitted curves with facet wrap by species
p.size.notoscriptus.pd <- ggplot() +
  geom_point(data = m.plot_data.size.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.size.notoscriptus.pd, aes(x = Temperature, y = Fitted, color = Size), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species)
p.size.notoscriptus.pd


# Merge the observed and predicted data for comparison
merged_data.size.notoscriptus.pd <- fitted_data.size.notoscriptus.pd %>%
  left_join(m.plot_data.size, by = c("Temperature", "Size", "Species"))


# Calculate RMSE for each species
mse_by_size.notoscriptus.pd <- merged_data.size.notoscriptus.pd %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.size.notoscriptus.pd <- merged_data.size.notoscriptus.pd |>
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE)))

save(mse_by_size.notoscriptus.pd,
     m.overall_rmse.size.notoscriptus.pd,
     fitted_data.size.notoscriptus.pd,
     merged_data.size.notoscriptus.pd,
     size.notoscriptus.pd.fit$BUGSoutput$DIC,
     p.size.notoscriptus.pd,
     file="notoscriptus.pd.size.RData")


######################################### By Wing Size and Climate- Out of sample prediction Aegypti #########################################
aedes.rm.pd <- mosquito.pd[mosquito.pd$interactor1 != "Aedes aegypti", ]
m.N.obs.pd.rm <- length(aedes.rm.pd$Trait)
m.climate.pd.rm <- as.factor(aedes.rm.pd$climate)
m.n.climate.rm.pd <- length(unique(aedes.rm.pd$climate))
m.size.pd.rm <- as.factor(aedes.rm.pd$wing.size)
m.n.size.pd.rm <- length(unique(aedes.rm.pd$wing.size))
unique(aedes.rm.pd$interactor1)

sink("aedes.by.climate.size.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait

  # Size-specific
  mu.c ~ dexp(10)             # Mean for size effects c
  # Climate-specific
  mu.Topt ~ dnorm(30, 0.1)             # Mean for climate effects Topt
  tau.Topt ~ dexp(0.01)             # Tau for climate effect Topt
  mu.a ~ dexp(0.01)            # Mean for climate effects Topt

  for (j in 1:m.n.climate.rm.pd) {
    mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific temp
    mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific temp
  }
    for (k in 1:m.n.size.pd.rm) {
    mu.size[k] ~ dexp(mu.c)       # Size-specific low point
    size.tau[k] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:m.N.obs.pd.rm) {
    mu[i] <- mu.climate.a[m.climate.pd.rm[i]] * (mu.climate.Topt[m.climate.pd.rm[i]] - temp[i])^2 + mu.size[m.size.pd.rm[i]]
    trait[i] ~ dnorm(mu[i], size.tau[m.size.pd.rm[i]])  T(0, )
  }}
", file = "aedes.by.climate.size.txt")


# Parameters to monitor
parameters <- c("mu.a", "mu.c", "mu.size", "size.tau",
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

jag.data<-list(trait = aedes.rm.pd$Trait, m.N.obs.pd.rm = m.N.obs.pd.rm, temp = aedes.rm.pd$Temp,
               m.n.climate.rm.pd = m.n.climate.rm.pd, m.climate.pd.rm = m.climate.pd.rm, m.n.size.pd.rm = m.n.size.pd.rm, m.size.pd.rm = m.size.pd.rm)

m.climate.pd.size.aegypti.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters,
                                      model.file="aedes.by.climate.size.txt", n.thin=nt, n.chains=nc, n.burnin=nb,
                                      n.iter=ni, DIC=T, working.directory=getwd())
m.climate.pd.size.aegypti.fit.mcmc <- as.mcmc(m.climate.pd.size.aegypti.fit) ## makes an "mcmc" object
closeAllConnections()

# Get posterior means
m.summary_fit.climate.size.pd <- summary(m.climate.pd.size.aegypti.fit.mcmc)$statistics
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
m.fitted_data.climate.size.aegypti.pd <- data.frame(Temperature = rep(m.temp_seq, m.n.climate.rm.pd*3),
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

m.species.pd_to_climate.size <- unique(m.plot_data.climate.size.pd[, c("Species", "Climate", "Size")])
m.species.pd_to_climate.size$Climate <- as.factor(m.species.pd_to_climate.size$Climate)
m.species.pd_to_climate.size$Size <- as.factor(m.species.pd_to_climate.size$Size)

# Merge fitted climate curves with species-to-climate mapping and filter for matching species
m.fitted_data.climate.size.aegypti.pd <- data.frame(
  Temperature = rep(m.temp_seq, m.n.climate.rm.pd * m.n.size.pd.rm),
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
m.fitted_data.climate.size.aegypti.pd$Combo <- with(m.fitted_data.climate.size.aegypti.pd, paste(Climate, Size, sep = ", "))
m.fitted_data.climate.size.aegypti.pd$Combo <- factor(m.fitted_data.climate.size.aegypti.pd$Combo, levels = c("temperate, large",
                                                                                                              "subtropical, medium", "subtropical, large",
                                                                                                              "tropical, small", "tropical, medium"))
# Plot observed data and fitted curves with facet wrap by species
p.climate.size.aegypti.pd <- ggplot() +
  geom_point(data = m.plot_data.climate.size.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.climate.size.aegypti.pd, aes(x = Temperature, y = Fitted, color = Combo), size = 1) +
  labs(x = "Temperature",
       y = "Development Time in Days",
       color = "Climate_Size Combo"
  ) +
  theme_minimal() +
  facet_wrap(~ Species)
p.climate.size.aegypti.pd

# Merge the observed and predicted data for comparison
m.merged_data.climate.size.aegypti.pd <- m.fitted_data.climate.size.aegypti.pd %>%
  left_join(m.plot_data.climate.size, by = c("Temperature", "Climate", "Species"))


# Calculate RMSE for each species
m.mse_by_climate.size.aegypti.pd <- m.merged_data.climate.size.aegypti.pd %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.climate.size.aegypti.pd <- m.merged_data.climate.size.aegypti.pd |>
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE)))

save(m.mse_by_climate.size.aegypti.pd,
     m.overall_rmse.climate.size.aegypti.pd,
     m.fitted_data.climate.size.aegypti.pd,
     m.merged_data.climate.size.aegypti.pd,
     m.climate.pd.size.aegypti.fit$BUGSoutput$DIC,
     p.climate.size.aegypti.pd,
     file="aegypti.climate.size.pd.RData")


######################################### By Wing Size and Climate- Out of sample prediction Krombeini #########################################
aedes.rm.pd <- mosquito.pd[mosquito.pd$interactor1 != "Aedes krombeini", ]
m.N.obs.pd.rm <- length(aedes.rm.pd$Trait)
m.climate.pd.rm <- as.factor(aedes.rm.pd$climate)
m.n.climate.rm.pd <- length(unique(aedes.rm.pd$climate))
m.size.pd.rm <- as.factor(aedes.rm.pd$wing.size)
m.n.size.pd.rm <- length(unique(aedes.rm.pd$wing.size))
unique(aedes.rm.pd$interactor1)

sink("aedes.by.climate.size.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait

  # Size-specific
  mu.c ~ dexp(10)             # Mean for size effects c
  # Climate-specific
  mu.Topt ~ dnorm(30, 0.1)             # Mean for climate effects Topt
  tau.Topt ~ dexp(0.01)             # Tau for climate effect Topt
  mu.a ~ dexp(0.01)            # Mean for climate effects Topt

  for (j in 1:m.n.climate.rm.pd) {
    mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific temp
    mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific temp
  }
    for (k in 1:m.n.size.pd.rm) {
    mu.size[k] ~ dexp(mu.c)       # Size-specific low point
    size.tau[k] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:m.N.obs.pd.rm) {
    mu[i] <- mu.climate.a[m.climate.pd.rm[i]] * (mu.climate.Topt[m.climate.pd.rm[i]] - temp[i])^2 + mu.size[m.size.pd.rm[i]]
    trait[i] ~ dnorm(mu[i], size.tau[m.size.pd.rm[i]])  T(0, )
  }}
", file = "aedes.by.climate.size.txt")

# Parameters to monitor
parameters <- c("mu.a", "mu.c", "mu.size", "size.tau",
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

jag.data<-list(trait = aedes.rm.pd$Trait, m.N.obs.pd.rm = m.N.obs.pd.rm, temp = aedes.rm.pd$Temp,
               m.n.climate.rm.pd = m.n.climate.rm.pd, m.climate.pd.rm = m.climate.pd.rm, m.n.size.pd.rm = m.n.size.pd.rm, m.size.pd.rm = m.size.pd.rm)

m.climate.pd.size.krombeini.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters,
                                        model.file="aedes.by.climate.size.txt", n.thin=nt, n.chains=nc, n.burnin=nb,
                                        n.iter=ni, DIC=T, working.directory=getwd())
m.climate.pd.size.krombeini.fit.mcmc <- as.mcmc(m.climate.pd.size.krombeini.fit) ## makes an "mcmc" object
closeAllConnections()

# Get posterior means
m.summary_fit.climate.size.pd <- summary(m.climate.pd.size.krombeini.fit.mcmc)$statistics
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
m.fitted_data.climate.size.krombeini.pd <- data.frame(Temperature = rep(m.temp_seq, m.n.climate.rm.pd*3),
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

m.species.pd_to_climate.size <- unique(m.plot_data.climate.size.pd[, c("Species", "Climate", "Size")])
m.species.pd_to_climate.size$Climate <- as.factor(m.species.pd_to_climate.size$Climate)
m.species.pd_to_climate.size$Size <- as.factor(m.species.pd_to_climate.size$Size)

# Merge fitted climate curves with species-to-climate mapping and filter for matching species
m.fitted_data.climate.size.krombeini.pd <- data.frame(
  Temperature = rep(m.temp_seq, m.n.climate.rm.pd * m.n.size.pd.rm),
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
m.fitted_data.climate.size.krombeini.pd$Combo <- with(m.fitted_data.climate.size.krombeini.pd, paste(Climate, Size, sep = ", "))
m.fitted_data.climate.size.krombeini.pd$Combo <- factor(m.fitted_data.climate.size.krombeini.pd$Combo, levels = c("temperate, large",
                                                                                                                  "subtropical, medium", "subtropical, large",
                                                                                                                  "tropical, small", "tropical, medium"))
# Plot observed data and fitted curves with facet wrap by species
p.climate.size.krombeini.pd <- ggplot() +
  geom_point(data = m.plot_data.climate.size.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.climate.size.krombeini.pd, aes(x = Temperature, y = Fitted, color = Combo), size = 1) +
  labs(x = "Temperature",
       y = "Development Time in Days",
       color = "Climate_Size Combo"
  ) +
  theme_minimal() +
  facet_wrap(~ Species)
p.climate.size.krombeini.pd

# Merge the observed and predicted data for comparison
m.merged_data.climate.size.krombeini.pd <- m.fitted_data.climate.size.krombeini.pd %>%
  left_join(m.plot_data.climate.size, by = c("Temperature", "Climate", "Species"))


# Calculate RMSE for each species
m.mse_by_climate.size.krombeini.pd <- m.merged_data.climate.size.krombeini.pd %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.climate.size.krombeini.pd <- m.merged_data.climate.size.krombeini.pd |>
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE)))

save(m.mse_by_climate.size.krombeini.pd,
     m.overall_rmse.climate.size.krombeini.pd,
     m.fitted_data.climate.size.krombeini.pd,
     m.merged_data.climate.size.krombeini.pd,
     m.climate.pd.size.krombeini.fit$BUGSoutput$DIC,
     p.climate.size.krombeini.pd,
     file="krombeini.climate.size.pd.RData")

######################################### RMSE ###################################
# # Data frame
# table.mosquito.rm.pd <- data.frame(
#   Method = factor(c(rep("By Climate (RM)", 7 * 7), 
#                     rep("By Wing Size (RM)", 7 * 6), 
#                     rep("By Wing Size and Climate (RM)", 7 * 4))),
#   Species = factor(rep(c("Aedes aegypti", "Aedes albopictus",
#                          "Aedes atropalpus", "Aedes camptorhynchus",
#                          "Aedes japonicus japonicus", "Aedes krombeini",
#                          "Aedes notoscriptus"), times = 17)),
#   MSE = c(m.mse_by_climate.aegypti.pd, m.mse_by_climate.albopictus.pd,
#           m.mse_by_climate.atropalpus.pd,
#           m.mse_by_climate.camptorhynchus.pd, m.mse_by_climate.jj.pd,
#           m.mse_by_climate.krombeini.pd, m.mse_by_climate.notoscriptus.pd,
#           mse_by_size.aegypti.pd, mse_by_size.atropalpus.pd,
#           mse_by_size.camptorhynchus.pd, mse_by_size.jj.pd,
#           mse_by_size.krombeini.pd, mse_by_size.notoscriptus.pd,
#           m.mse_by_climate.size.aegypti.pd,
#           m.mse_by_climate.size.atropalpus.pd,
#           m.mse_by_climate.size.jj.pd,
#           m.mse_by_climate.size.krombeini.pd),
#   Total.RMSE = rep(as.numeric(c(m.overall_rmse.climate.aegypti.pd, m.overall_rmse.climate.albopictus.pd,
#                                 m.overall_rmse.climate.atropalpus.pd,
#                                 m.overall_rmse.climate.camptorhynchus.pd ,m.overall_rmse.climate.jj.pd,
#                                 m.overall_rmse.climate.krombeini.pd, m.overall_rmse.climate.notoscriptus.pd,
#                                 m.overall_rmse.size.aegypti.pd, m.overall_rmse.size.atropalpus.pd,
#                                 m.overall_rmse.size.camptorhynchus.pd, m.overall_rmse.size.jj.pd,
#                                 m.overall_rmse.size.krombeini.pd, m.overall_rmse.size.notoscriptus.pd,
#                                 m.overall_rmse.climate.size.aegypti.pd,
#                                 m.overall_rmse.climate.size.atropalpus.pd,
#                                 m.overall_rmse.climate.size.jj.pd,
#                                 m.overall_rmse.climate.size.krombeini.pd)), each = 7),
#   Removed = c(rep(c("Aedes aegypti", "Aedes albopictus",
#                     "Aedes atropalpus", "Aedes camptorhynchus",
#                     "Aedes japonicus japonicus", "Aedes krombeini",
#                     "Aedes notoscriptus"), each = 7),
#               rep(c("Aedes aegypti",
#                     "Aedes atropalpus", "Aedes camptorhynchus",
#                     "Aedes japonicus japonicus", "Aedes krombeini",
#                     "Aedes notoscriptus"), each = 7),
#               rep(c("Aedes aegypti", 
#                     "Aedes atropalpus",
#                     "Aedes japonicus japonicus", 
#                     "Aedes krombeini"), each = 7))
# )
# 
# mosquito.DIC.rm.pd <- data.frame(
#   Method = rep(c("By Climate", "By Wing Size",
#                  "By Climate and Wing Size"), each = 4),
#   Removed = rep(c("Aedes aegypti", "Aedes atropalpus",
#                   "Aedes japonicus japonicus", "Aedes krombeini"), times = 3),
#   DIC = c(climate.aegypti.pd.fit$BUGSoutput$DIC, climate.atropalpus.pd.fit$BUGSoutput$DIC,
#           climate.jj.pd.fit$BUGSoutput$DIC, climate.notoscriptus.pd.fit$BUGSoutput$DIC,
#           size.aegypti.pd.fit$BUGSoutput$DIC, size.atropalpus.pd.fit$BUGSoutput$DIC,
#           size.jj.pd.fit$BUGSoutput$DIC, size.notoscriptus.pd.fit$BUGSoutput$DIC,
#           m.climate.pd.size.aegypti.fit$BUGSoutput$DIC, m.climate.pd.size.atropalpus.fit$BUGSoutput$DIC,
#           m.climate.pd.size.jj.fit$BUGSoutput$DIC, m.climate.pd.size.krombeini.fit$BUGSoutput$DIC)
# )
# 
# table.mosquito.rm.pd |> 
#   ggplot(aes(x = Removed, y = Total.RMSE, fill = Method))+
#   geom_bar(stat = "identity", position = "dodge",
#            color = "black")+
#   theme_minimal() +
#   scale_fill_viridis_d(option = "mako") +
#   coord_cartesian(ylim = c(3.5, 5.5)) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# mosquito.DIC.rm.pd |> 
#   ggplot(aes(x = Removed, y = DIC, fill = Method))+
#   geom_bar(stat = "identity", position = "dodge", color = "black")+
#   theme_minimal() +
#   scale_fill_viridis_d(option = "mako") +
#   coord_cartesian(ylim = c(2000, 6200)) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# save(table.mosquito.rm.pd,
#      mosquito.DIC.rm.pd, file = "Dataset.Mosquito.RM.RData")
# 

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
# ni <- 60000 # number of iterations in each chain
# nb <- 30000 # number of 'burn in' iterations to discard
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
# size.rm.pd <- as.factor(aedes.rm.pd$wing.size)
# n.size.rm.pd <- length(unique(aedes.rm.pd$wing.size))
# N.obs.rm.pd <- length(aedes.rm.pd$Trait)
# unique(aedes.rm.pd$interactor1)
# sink("aedes.by.size.txt")
# cat("
# model {
#   # Priors
#   tau1 ~ dexp(0.01)           # Tau for trait
#   
#   # Size-specific 
#   mu.c ~ dexp(10)             # Mean for size effects c
#   mu.Topt ~ dnorm(30, 0.1)             # Mean for size effects Topt
#   tau.Topt ~ dexp(0.01)             # Tau for size effect Topt
#   mu.a ~ dexp(0.01)            # Mean for size effects Topt
# 
#   for (j in 1:n.size.rm.pd) {
#     mu.size[j] ~ dexp(mu.c)       # Size-specific low point
#     mu.size.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Size-specific temp
#     mu.size.a[j] ~ dexp(mu.a)    # Size-specific temp
#   }
# 
#   # Likelihood
#   for (i in 1:N.obs.rm.pd) {
#     mu[i] <- mu.size.a[size.rm.pd[i]] * (mu.size.Topt[size.rm.pd[i]] - temp[i])^2 + mu.size[size.rm.pd[i]]
#     trait[i] ~ dnorm(mu[i], tau1) T(0, ) 
#   }}
# ", file = "aedes.by.size.txt")
# 
# 
# # Parameters to monitor
# parameters <- c("mu.a", "mu.c", "mu.size", 
#                 "tau1", "mu.Topt", "tau.Topt", "mu.size.Topt", "mu.size.a")
# 
# 
# # MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
# ni <- 60000 # number of iterations in each chain
# nb <- 30000 # number of 'burn in' iterations to discard
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
#                  n.size.rm.pd = n.size.rm.pd, size.rm.pd = size.rm.pd)
# 
# size.rm.jj.pd.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
#                        model.file="aedes.by.size.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
#                        n.iter=ni, DIC=T, working.directory=getwd())
# size.rm.jj.pd.fit.mcmc <- as.mcmc(size.rm.jj.pd.fit) ## makes an "mcmc" object
# closeAllConnections()
# 
# # Get posterior means
# summary_fit.size.rm.pd <- summary(size.rm.jj.pd.fit.mcmc)$statistics
# a_mean_size.rm.pd <- summary_fit.size.rm.pd[grep("^mu.size.a\\[", rownames(summary_fit.size.rm.pd)), "Mean"]
# c_mean_size.rm.pd <- summary_fit.size.rm.pd[grep("^mu.size\\[", rownames(summary_fit.size.rm.pd)), "Mean"]
# Topt_mean_size.rm.pd <- summary_fit.size.rm.pd[grep("^mu.size.Topt\\[", rownames(summary_fit.size.rm.pd)), "Mean"]
# 
# # Generate temperature sequence 
# m.temp_seq <- seq(min(mosquito.pd$Temp), max(mosquito.pd$Temp), by = 0.1)
# 
# # Get predicted values
# m.size.pd1.rm <- function (x){
#   a_mean_size.rm.pd[1] * (x - Topt_mean_size.rm.pd[1])^2 + c_mean_size.rm.pd[1]
# }
# m.size.pd2.rm <- function (x){
#   a_mean_size.rm.pd[2] * (x - Topt_mean_size.rm.pd[2])^2 + c_mean_size.rm.pd[2]
# }
# m.size.pd3.rm <- function (x){
#   a_mean_size.rm.pd[3] * (x - Topt_mean_size.rm.pd[3])^2 + c_mean_size.rm.pd[3]
# }
# 
# 
# # Optionally, convert to a data frame for easier plotting
# fitted_data.size.rm.jj.pd <- data.frame(Temperature = rep(m.temp_seq, n.size.rm.pd), 
#                                      Fitted = c(m.size.pd1.rm(m.temp_seq),
#                                                 m.size.pd2.rm(m.temp_seq),
#                                                 m.size.pd3.rm(m.temp_seq)),
#                                      Size = factor(rep(c("large", "medium", "small"), each = length(temp_seq))))
# 
# m.plot_data.size.pd <- data.frame(
#   Temperature = mosquito.pd$Temp,
#   Trait = mosquito.pd$Trait,
#   Size = as.factor(mosquito.pd$wing.size),
#   Species = as.factor(mosquito.pd$interactor1)
# )
# 
# m.species.pd_to_size <- unique(m.plot_data.size.pd[, c("Species", "Size")])  # Ensure species matches correct climate
# m.species.pd_to_size$Size <- as.factor(m.species.pd_to_size$Size)
# # Merge fitted climate curves with species-to-climate mapping
# fitted_data.size.rm.jj.pd <- data.frame(
#   Temperature = rep(m.temp_seq, n.size.rm.pd),
#   Fitted = c(m.size.pd1.rm(m.temp_seq),
#              m.size.pd2.rm(m.temp_seq),
#              m.size.pd3.rm(m.temp_seq)),
#   Size = factor(rep(c("large", "medium", "small"), each = length(temp_seq)))
# ) %>%
#   left_join(m.species.pd_to_size, by = "Size")
# 
# m.plot_data.size.pd$Size <- factor(m.plot_data.size.pd$Size, levels = c("small", "medium", "large"))
# fitted_data.size.rm.jj.pd$Size <- factor(fitted_data.size.rm.jj.pd$Size, levels = c("small", "medium", "large"))
# # Plot observed data and fitted curves with facet wrap by species
# p.size.rm.jj.pd <- ggplot() +
#   geom_point(data = m.plot_data.size.pd, aes(x = Temperature, y = Trait), alpha = 0.5) +
#   geom_line(data = fitted_data.size.rm.jj.pd, aes(x = Temperature, y = Fitted, color = Size), size = 1) +
#   labs(x = "Temperature", y = "Development Time in Days") +
#   theme_minimal() +
#   facet_wrap(~ Species) 
# p.size.rm.jj.pd
# 
# 
# # Merge the observed and predicted data for comparison
# merged_data.size.rm.jj.pd <- fitted_data.size.rm.jj.pd %>%
#   left_join(m.plot_data.size, by = c("Temperature", "Size", "Species"))
# 
# 
# # Calculate RMSE for each species
# mse_by_species.size.rm.jj.pd <- merged_data.size.rm.jj.pd %>%
#   group_by(Species) %>%
#   summarize(
#     MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
#   )
# 
# 
# save(mse_by_species.size.rm.jj.pd,
#      fitted_data.size.rm.jj.pd,
#      merged_data.size.rm.jj.pd,
#      size.rm.jj.pd.fit.mcmc,
#      size.rm.jj.pd.fit,
#      p.size.rm.jj.pd,
#      file="jj.size.pd.rm.RData")
# 
# 
# 
# 
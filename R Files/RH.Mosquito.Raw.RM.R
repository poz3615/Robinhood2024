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
test <- test[,c("Trait", "originaltraitunit","interactor1", "Temp", "interactor1number", 
                "originalerrorpos", "originalerrorunit")]
test$Trait <- ifelse(test$originaltraitunit == "hours", test$Trait/24, test$Trait)
test$Trait <- ifelse(test$originaltraitunit %in% c("days-1", "day-1"), 1/test$Trait, test$Trait)
#test <- test[test$originaltraitunit %in% c("days","day"),]
test <- test[test$Trait > 0.1 & test$Trait < 100, ]
m.N.obs <- length(test$Trait)

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
m.species <- as.factor(test$interactor1)
m.n.species <- length(unique(test$interactor1))

######################################### By Climate- Aegypti #########################################
aedes.rm <- test[test$interactor1 != "Aedes aegypti", ]
m.climate.rm <- as.factor(aedes.rm$climate)
m.n.climate.rm <- length(unique(aedes.rm$climate))
m.N.obs.rm <- length(aedes.rm$Trait)
unique(aedes.rm$interactor1)
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

  for (j in 1:m.n.climate.rm) {
    mu.climate[j] ~ dexp(mu.c)       # Climate-specific low point
    mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific temp
    mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific temp
    climate.tau[j] ~ dexp(tau1)
  }

  # Likelihood
  for (i in 1:m.N.obs.rm) {
    mu[i] <- mu.climate.a[m.climate.rm[i]] * (mu.climate.Topt[m.climate.rm[i]] - temp[i])^2 + mu.climate[m.climate.rm[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[m.climate.rm[i]]) T(0, ) 
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
#jag.data<-list(trait = test$Trait, N.obs = N.obs, temp = test$Temp)

jag.data <- list(trait = aedes.rm$Trait, m.N.obs.rm = m.N.obs.rm, temp = aedes.rm$Temp, 
                 m.n.climate.rm = m.n.climate.rm, m.climate.rm = m.climate.rm)

climate.aegypti.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                            model.file="aedes.by.climate.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                            n.iter=ni, DIC=T, working.directory=getwd())
climate.aegypti.fit.mcmc <- as.mcmc(climate.aegypti.fit) ## makes an "mcmc" object
closeAllConnections()

# Get posterior means
summary_fit.rm.climate <- summary(climate.aegypti.fit.mcmc)$statistics
a_mean_climate.rm <- summary_fit.rm.climate[grep("^mu.climate.a\\[", rownames(summary_fit.rm.climate)), "Mean"]
c_mean_climate.rm <- summary_fit.rm.climate[grep("^mu.climate\\[", rownames(summary_fit.rm.climate)), "Mean"]
Topt_mean_climate.rm <- summary_fit.rm.climate[grep("^mu.climate.Topt\\[", rownames(summary_fit.rm.climate)), "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

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
fitted_data.climate.aegypti <- data.frame(Temperature = rep(m.temp_seq, m.n.climate.rm), 
                                          Fitted = c(clim1.rm(m.temp_seq),
                                                     clim2.rm(m.temp_seq),
                                                     clim3.rm(m.temp_seq)),
                                          Climate = factor(rep(c("subtropical", "temperate", "tropical"), each = length(m.temp_seq))))

m.plot_data.climate <- data.frame(
  Temperature = test$Temp,
  Trait = test$Trait,
  Climate = as.factor(test$climate),
  Species = as.factor(test$interactor1)
)

m.species_to_climate <- unique(m.plot_data.climate[, c("Species", "Climate")])  # Ensure species matches correct climate
m.species_to_climate$Climate <- as.factor(m.species_to_climate$Climate)
# Merge fitted climate curves with species-to-climate mapping
fitted_data.climate.aegypti <- data.frame(
  Temperature = rep(m.temp_seq, m.n.climate.rm),
  Fitted = c(clim1.rm(m.temp_seq),
             clim2.rm(m.temp_seq),
             clim3.rm(m.temp_seq)),
  Climate = factor(rep(c("subtropical", "temperate", "tropical"), each = length(m.temp_seq)))
) %>%
  left_join(m.species_to_climate, by = "Climate")

m.plot_data.climate$Climate <- factor(m.plot_data.climate$Climate, levels = c("temperate", "subtropical", "tropical"))
fitted_data.climate.aegypti$Climate <- factor(fitted_data.climate.aegypti$Climate, levels = c("temperate", "subtropical", "tropical"))
# Plot observed data and fitted curves with facet wrap by species
p.climate.aegypti <- ggplot() +
  geom_point(data = m.plot_data.climate, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.climate.aegypti, aes(x = Temperature, y = Fitted, color = Climate), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )   
p.climate.aegypti

# Merge the observed and predicted data for comparison
merged_data.climate.aegypti <- fitted_data.climate.aegypti %>%
  left_join(m.plot_data.climate, by = c("Temperature", "Climate", "Species"))


# Calculate RMSE for each climate
m.mse_by_climate.aegypti <- merged_data.climate.aegypti |> 
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.climate.aegypti <- merged_data.climate.aegypti |> 
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(m.mse_by_climate.aegypti,
     m.overall_rmse.climate.aegypti,
     fitted_data.climate.aegypti,
     merged_data.climate.aegypti,
     climate.aegypti.fit.mcmc,
     climate.aegypti.fit,
     p.climate.aegypti,
     file="aegypti.climate.RData")


######################################### By Climate- Albopictus #########################################
aedes.rm <- test[test$interactor1 != "Aedes albopictus", ]
m.climate.rm <- as.factor(aedes.rm$climate)
m.n.climate.rm <- length(unique(aedes.rm$climate))
m.N.obs.rm <- length(aedes.rm$Trait)
unique(aedes.rm$interactor1)
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

  for (j in 1:m.n.climate.rm) {
    mu.climate[j] ~ dexp(mu.c)       # Climate-specific low point
    mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific temp
    mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific temp
    climate.tau[j] ~ dexp(tau1)
  }

  # Likelihood
  for (i in 1:m.N.obs.rm) {
    mu[i] <- mu.climate.a[m.climate.rm[i]] * (mu.climate.Topt[m.climate.rm[i]] - temp[i])^2 + mu.climate[m.climate.rm[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[m.climate.rm[i]]) T(0, ) 
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
#jag.data<-list(trait = test$Trait, N.obs = N.obs, temp = test$Temp)

jag.data <- list(trait = aedes.rm$Trait, m.N.obs.rm = m.N.obs.rm, temp = aedes.rm$Temp, 
                 m.n.climate.rm = m.n.climate.rm, m.climate.rm = m.climate.rm)

climate.albopictus.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                            model.file="aedes.by.climate.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                            n.iter=ni, DIC=T, working.directory=getwd())
climate.albopictus.fit.mcmc <- as.mcmc(climate.albopictus.fit) ## makes an "mcmc" object
closeAllConnections()

# Get posterior means
summary_fit.rm.climate <- summary(climate.albopictus.fit.mcmc)$statistics
a_mean_climate.rm <- summary_fit.rm.climate[grep("^mu.climate.a\\[", rownames(summary_fit.rm.climate)), "Mean"]
c_mean_climate.rm <- summary_fit.rm.climate[grep("^mu.climate\\[", rownames(summary_fit.rm.climate)), "Mean"]
Topt_mean_climate.rm <- summary_fit.rm.climate[grep("^mu.climate.Topt\\[", rownames(summary_fit.rm.climate)), "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

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
fitted_data.climate.albopictus <- data.frame(Temperature = rep(m.temp_seq, m.n.climate.rm), 
                                          Fitted = c(clim1.rm(m.temp_seq),
                                                     clim2.rm(m.temp_seq),
                                                     clim3.rm(m.temp_seq)),
                                          Climate = factor(rep(c("subtropical", "temperate", "tropical"), each = length(m.temp_seq))))

m.plot_data.climate <- data.frame(
  Temperature = test$Temp,
  Trait = test$Trait,
  Climate = as.factor(test$climate),
  Species = as.factor(test$interactor1)
)

m.species_to_climate <- unique(m.plot_data.climate[, c("Species", "Climate")])  # Ensure species matches correct climate
m.species_to_climate$Climate <- as.factor(m.species_to_climate$Climate)
# Merge fitted climate curves with species-to-climate mapping
fitted_data.climate.albopictus <- data.frame(
  Temperature = rep(m.temp_seq, m.n.climate.rm),
  Fitted = c(clim1.rm(m.temp_seq),
             clim2.rm(m.temp_seq),
             clim3.rm(m.temp_seq)),
  Climate = factor(rep(c("subtropical", "temperate", "tropical"), each = length(m.temp_seq)))
) %>%
  left_join(m.species_to_climate, by = "Climate")

m.plot_data.climate$Climate <- factor(m.plot_data.climate$Climate, levels = c("temperate", "subtropical", "tropical"))
fitted_data.climate.albopictus$Climate <- factor(fitted_data.climate.albopictus$Climate, levels = c("temperate", "subtropical", "tropical"))
# Plot observed data and fitted curves with facet wrap by species
p.climate.albopictus <- ggplot() +
  geom_point(data = m.plot_data.climate, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.climate.albopictus, aes(x = Temperature, y = Fitted, color = Climate), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )   
p.climate.albopictus

# Merge the observed and predicted data for comparison
merged_data.climate.albopictus <- fitted_data.climate.albopictus %>%
  left_join(m.plot_data.climate, by = c("Temperature", "Climate", "Species"))


# Calculate RMSE for each climate
m.mse_by_climate.albopictus <- merged_data.climate.albopictus |> 
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.climate.albopictus <- merged_data.climate.albopictus |> 
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(m.mse_by_climate.albopictus,
     m.overall_rmse.climate.albopictus,
     fitted_data.climate.albopictus,
     merged_data.climate.albopictus,
     climate.albopictus.fit.mcmc,
     climate.albopictus.fit,
     p.climate.albopictus,
     file="albopictus.climate.RData")


######################################### By Climate- Atropalpus #########################################
aedes.rm <- test[test$interactor1 != "Aedes atropalpus", ]
m.climate.rm <- as.factor(aedes.rm$climate)
m.n.climate.rm <- length(unique(aedes.rm$climate))
m.N.obs.rm <- length(aedes.rm$Trait)
unique(aedes.rm$interactor1)
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

  for (j in 1:m.n.climate.rm) {
    mu.climate[j] ~ dexp(mu.c)       # Climate-specific low point
    mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific temp
    mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific temp
    climate.tau[j] ~ dexp(tau1)
  }

  # Likelihood
  for (i in 1:m.N.obs.rm) {
    mu[i] <- mu.climate.a[m.climate.rm[i]] * (mu.climate.Topt[m.climate.rm[i]] - temp[i])^2 + mu.climate[m.climate.rm[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[m.climate.rm[i]]) T(0, ) 
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
#jag.data<-list(trait = test$Trait, N.obs = N.obs, temp = test$Temp)

jag.data <- list(trait = aedes.rm$Trait, m.N.obs.rm = m.N.obs.rm, temp = aedes.rm$Temp, 
                 m.n.climate.rm = m.n.climate.rm, m.climate.rm = m.climate.rm)

climate.atropalpus.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                            model.file="aedes.by.climate.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                            n.iter=ni, DIC=T, working.directory=getwd())
climate.atropalpus.fit.mcmc <- as.mcmc(climate.atropalpus.fit) ## makes an "mcmc" object
closeAllConnections()

# Get posterior means
summary_fit.rm.climate <- summary(climate.atropalpus.fit.mcmc)$statistics
a_mean_climate.rm <- summary_fit.rm.climate[grep("^mu.climate.a\\[", rownames(summary_fit.rm.climate)), "Mean"]
c_mean_climate.rm <- summary_fit.rm.climate[grep("^mu.climate\\[", rownames(summary_fit.rm.climate)), "Mean"]
Topt_mean_climate.rm <- summary_fit.rm.climate[grep("^mu.climate.Topt\\[", rownames(summary_fit.rm.climate)), "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

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
fitted_data.climate.atropalpus <- data.frame(Temperature = rep(m.temp_seq, m.n.climate.rm), 
                                          Fitted = c(clim1.rm(m.temp_seq),
                                                     clim2.rm(m.temp_seq),
                                                     clim3.rm(m.temp_seq)),
                                          Climate = factor(rep(c("subtropical", "temperate", "tropical"), each = length(m.temp_seq))))

m.plot_data.climate <- data.frame(
  Temperature = test$Temp,
  Trait = test$Trait,
  Climate = as.factor(test$climate),
  Species = as.factor(test$interactor1)
)

m.species_to_climate <- unique(m.plot_data.climate[, c("Species", "Climate")])  # Ensure species matches correct climate
m.species_to_climate$Climate <- as.factor(m.species_to_climate$Climate)
# Merge fitted climate curves with species-to-climate mapping
fitted_data.climate.atropalpus <- data.frame(
  Temperature = rep(m.temp_seq, m.n.climate.rm),
  Fitted = c(clim1.rm(m.temp_seq),
             clim2.rm(m.temp_seq),
             clim3.rm(m.temp_seq)),
  Climate = factor(rep(c("subtropical", "temperate", "tropical"), each = length(m.temp_seq)))
) %>%
  left_join(m.species_to_climate, by = "Climate")

m.plot_data.climate$Climate <- factor(m.plot_data.climate$Climate, levels = c("temperate", "subtropical", "tropical"))
fitted_data.climate.atropalpus$Climate <- factor(fitted_data.climate.atropalpus$Climate, levels = c("temperate", "subtropical", "tropical"))
# Plot observed data and fitted curves with facet wrap by species
p.climate.atropalpus <- ggplot() +
  geom_point(data = m.plot_data.climate, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.climate.atropalpus, aes(x = Temperature, y = Fitted, color = Climate), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )   
p.climate.atropalpus

# Merge the observed and predicted data for comparison
merged_data.climate.atropalpus <- fitted_data.climate.atropalpus %>%
  left_join(m.plot_data.climate, by = c("Temperature", "Climate", "Species"))


# Calculate RMSE for each climate
m.mse_by_climate.atropalpus <- merged_data.climate.atropalpus |> 
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.climate.atropalpus <- merged_data.climate.atropalpus |> 
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(m.mse_by_climate.atropalpus,
     m.overall_rmse.climate.atropalpus,
     fitted_data.climate.atropalpus,
     merged_data.climate.atropalpus,
     climate.atropalpus.fit.mcmc,
     climate.atropalpus.fit,
     p.climate.atropalpus,
     file="atropalpus.climate.RData")


######################################### By Climate- Camptorhynchus #########################################
aedes.rm <- test[test$interactor1 != "Aedes camptorhynchus", ]
m.climate.rm <- as.factor(aedes.rm$climate)
m.n.climate.rm <- length(unique(aedes.rm$climate))
m.N.obs.rm <- length(aedes.rm$Trait)
unique(aedes.rm$interactor1)
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

  for (j in 1:m.n.climate.rm) {
    mu.climate[j] ~ dexp(mu.c)       # Climate-specific low point
    mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific temp
    mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific temp
    climate.tau[j] ~ dexp(tau1)
  }

  # Likelihood
  for (i in 1:m.N.obs.rm) {
    mu[i] <- mu.climate.a[m.climate.rm[i]] * (mu.climate.Topt[m.climate.rm[i]] - temp[i])^2 + mu.climate[m.climate.rm[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[m.climate.rm[i]]) T(0, ) 
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
#jag.data<-list(trait = test$Trait, N.obs = N.obs, temp = test$Temp)

jag.data <- list(trait = aedes.rm$Trait, m.N.obs.rm = m.N.obs.rm, temp = aedes.rm$Temp, 
                 m.n.climate.rm = m.n.climate.rm, m.climate.rm = m.climate.rm)

climate.camptorhynchus.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                            model.file="aedes.by.climate.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                            n.iter=ni, DIC=T, working.directory=getwd())
climate.camptorhynchus.fit.mcmc <- as.mcmc(climate.camptorhynchus.fit) ## makes an "mcmc" object
closeAllConnections()

# Get posterior means
summary_fit.rm.climate <- summary(climate.camptorhynchus.fit.mcmc)$statistics
a_mean_climate.rm <- summary_fit.rm.climate[grep("^mu.climate.a\\[", rownames(summary_fit.rm.climate)), "Mean"]
c_mean_climate.rm <- summary_fit.rm.climate[grep("^mu.climate\\[", rownames(summary_fit.rm.climate)), "Mean"]
Topt_mean_climate.rm <- summary_fit.rm.climate[grep("^mu.climate.Topt\\[", rownames(summary_fit.rm.climate)), "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

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
fitted_data.climate.camptorhynchus <- data.frame(Temperature = rep(m.temp_seq, m.n.climate.rm), 
                                          Fitted = c(clim1.rm(m.temp_seq),
                                                     clim2.rm(m.temp_seq),
                                                     clim3.rm(m.temp_seq)),
                                          Climate = factor(rep(c("subtropical", "temperate", "tropical"), each = length(m.temp_seq))))

m.plot_data.climate <- data.frame(
  Temperature = test$Temp,
  Trait = test$Trait,
  Climate = as.factor(test$climate),
  Species = as.factor(test$interactor1)
)

m.species_to_climate <- unique(m.plot_data.climate[, c("Species", "Climate")])  # Ensure species matches correct climate
m.species_to_climate$Climate <- as.factor(m.species_to_climate$Climate)
# Merge fitted climate curves with species-to-climate mapping
fitted_data.climate.camptorhynchus <- data.frame(
  Temperature = rep(m.temp_seq, m.n.climate.rm),
  Fitted = c(clim1.rm(m.temp_seq),
             clim2.rm(m.temp_seq),
             clim3.rm(m.temp_seq)),
  Climate = factor(rep(c("subtropical", "temperate", "tropical"), each = length(m.temp_seq)))
) %>%
  left_join(m.species_to_climate, by = "Climate")

m.plot_data.climate$Climate <- factor(m.plot_data.climate$Climate, levels = c("temperate", "subtropical", "tropical"))
fitted_data.climate.camptorhynchus$Climate <- factor(fitted_data.climate.camptorhynchus$Climate, levels = c("temperate", "subtropical", "tropical"))
# Plot observed data and fitted curves with facet wrap by species
p.climate.camptorhynchus <- ggplot() +
  geom_point(data = m.plot_data.climate, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.climate.camptorhynchus, aes(x = Temperature, y = Fitted, color = Climate), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )   
p.climate.camptorhynchus

# Merge the observed and predicted data for comparison
merged_data.climate.camptorhynchus <- fitted_data.climate.camptorhynchus %>%
  left_join(m.plot_data.climate, by = c("Temperature", "Climate", "Species"))


# Calculate RMSE for each climate
m.mse_by_climate.camptorhynchus <- merged_data.climate.camptorhynchus |> 
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.climate.camptorhynchus <- merged_data.climate.camptorhynchus |> 
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(m.mse_by_climate.camptorhynchus,
     m.overall_rmse.climate.camptorhynchus,
     fitted_data.climate.camptorhynchus,
     merged_data.climate.camptorhynchus,
     climate.camptorhynchus.fit.mcmc,
     climate.camptorhynchus.fit,
     p.climate.camptorhynchus,
     file="camptorhynchus.climate.RData")


######################################### By Climate- Japonicus Japonicus #########################################
aedes.rm <- test[test$interactor1 != "Aedes Japonicus Japonicus", ]
m.climate.rm <- as.factor(aedes.rm$climate)
m.n.climate.rm <- length(unique(aedes.rm$climate))
m.N.obs.rm <- length(aedes.rm$Trait)
unique(aedes.rm$interactor1)
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

  for (j in 1:m.n.climate.rm) {
    mu.climate[j] ~ dexp(mu.c)       # Climate-specific low point
    mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific temp
    mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific temp
    climate.tau[j] ~ dexp(tau1)
  }

  # Likelihood
  for (i in 1:m.N.obs.rm) {
    mu[i] <- mu.climate.a[m.climate.rm[i]] * (mu.climate.Topt[m.climate.rm[i]] - temp[i])^2 + mu.climate[m.climate.rm[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[m.climate.rm[i]]) T(0, ) 
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
#jag.data<-list(trait = test$Trait, N.obs = N.obs, temp = test$Temp)

jag.data <- list(trait = aedes.rm$Trait, m.N.obs.rm = m.N.obs.rm, temp = aedes.rm$Temp, 
                 m.n.climate.rm = m.n.climate.rm, m.climate.rm = m.climate.rm)

climate.jj.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                            model.file="aedes.by.climate.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                            n.iter=ni, DIC=T, working.directory=getwd())
climate.jj.fit.mcmc <- as.mcmc(climate.jj.fit) ## makes an "mcmc" object
closeAllConnections()

# Get posterior means
summary_fit.rm.climate <- summary(climate.jj.fit.mcmc)$statistics
a_mean_climate.rm <- summary_fit.rm.climate[grep("^mu.climate.a\\[", rownames(summary_fit.rm.climate)), "Mean"]
c_mean_climate.rm <- summary_fit.rm.climate[grep("^mu.climate\\[", rownames(summary_fit.rm.climate)), "Mean"]
Topt_mean_climate.rm <- summary_fit.rm.climate[grep("^mu.climate.Topt\\[", rownames(summary_fit.rm.climate)), "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

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
fitted_data.climate.jj <- data.frame(Temperature = rep(m.temp_seq, m.n.climate.rm), 
                                          Fitted = c(clim1.rm(m.temp_seq),
                                                     clim2.rm(m.temp_seq),
                                                     clim3.rm(m.temp_seq)),
                                          Climate = factor(rep(c("subtropical", "temperate", "tropical"), each = length(m.temp_seq))))

m.plot_data.climate <- data.frame(
  Temperature = test$Temp,
  Trait = test$Trait,
  Climate = as.factor(test$climate),
  Species = as.factor(test$interactor1)
)

m.species_to_climate <- unique(m.plot_data.climate[, c("Species", "Climate")])  # Ensure species matches correct climate
m.species_to_climate$Climate <- as.factor(m.species_to_climate$Climate)
# Merge fitted climate curves with species-to-climate mapping
fitted_data.climate.jj <- data.frame(
  Temperature = rep(m.temp_seq, m.n.climate.rm),
  Fitted = c(clim1.rm(m.temp_seq),
             clim2.rm(m.temp_seq),
             clim3.rm(m.temp_seq)),
  Climate = factor(rep(c("subtropical", "temperate", "tropical"), each = length(m.temp_seq)))
) %>%
  left_join(m.species_to_climate, by = "Climate")

m.plot_data.climate$Climate <- factor(m.plot_data.climate$Climate, levels = c("temperate", "subtropical", "tropical"))
fitted_data.climate.jj$Climate <- factor(fitted_data.climate.jj$Climate, levels = c("temperate", "subtropical", "tropical"))
# Plot observed data and fitted curves with facet wrap by species
p.climate.jj <- ggplot() +
  geom_point(data = m.plot_data.climate, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.climate.jj, aes(x = Temperature, y = Fitted, color = Climate), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )   
p.climate.jj

# Merge the observed and predicted data for comparison
merged_data.climate.jj <- fitted_data.climate.jj %>%
  left_join(m.plot_data.climate, by = c("Temperature", "Climate", "Species"))


# Calculate RMSE for each climate
m.mse_by_climate.jj <- merged_data.climate.jj |> 
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.climate.jj <- merged_data.climate.jj |> 
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(m.mse_by_climate.jj,
     m.overall_rmse.climate.jj,
     fitted_data.climate.jj,
     merged_data.climate.jj,
     climate.jj.fit.mcmc,
     climate.jj.fit,
     p.climate.jj,
     file="jj.climate.RData")


######################################### By Climate- Krombeini #########################################
aedes.rm <- test[test$interactor1 != "Aedes krombeini", ]
m.climate.rm <- as.factor(aedes.rm$climate)
m.n.climate.rm <- length(unique(aedes.rm$climate))
m.N.obs.rm <- length(aedes.rm$Trait)
unique(aedes.rm$interactor1)
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

  for (j in 1:m.n.climate.rm) {
    mu.climate[j] ~ dexp(mu.c)       # Climate-specific low point
    mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific temp
    mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific temp
    climate.tau[j] ~ dexp(tau1)
  }

  # Likelihood
  for (i in 1:m.N.obs.rm) {
    mu[i] <- mu.climate.a[m.climate.rm[i]] * (mu.climate.Topt[m.climate.rm[i]] - temp[i])^2 + mu.climate[m.climate.rm[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[m.climate.rm[i]]) T(0, ) 
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
#jag.data<-list(trait = test$Trait, N.obs = N.obs, temp = test$Temp)

jag.data <- list(trait = aedes.rm$Trait, m.N.obs.rm = m.N.obs.rm, temp = aedes.rm$Temp, 
                 m.n.climate.rm = m.n.climate.rm, m.climate.rm = m.climate.rm)

climate.krombeini.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                            model.file="aedes.by.climate.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                            n.iter=ni, DIC=T, working.directory=getwd())
climate.krombeini.fit.mcmc <- as.mcmc(climate.krombeini.fit) ## makes an "mcmc" object
closeAllConnections()

# Get posterior means
summary_fit.rm.climate <- summary(climate.krombeini.fit.mcmc)$statistics
a_mean_climate.rm <- summary_fit.rm.climate[grep("^mu.climate.a\\[", rownames(summary_fit.rm.climate)), "Mean"]
c_mean_climate.rm <- summary_fit.rm.climate[grep("^mu.climate\\[", rownames(summary_fit.rm.climate)), "Mean"]
Topt_mean_climate.rm <- summary_fit.rm.climate[grep("^mu.climate.Topt\\[", rownames(summary_fit.rm.climate)), "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

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
fitted_data.climate.krombeini <- data.frame(Temperature = rep(m.temp_seq, m.n.climate.rm), 
                                          Fitted = c(clim1.rm(m.temp_seq),
                                                     clim2.rm(m.temp_seq),
                                                     clim3.rm(m.temp_seq)),
                                          Climate = factor(rep(c("subtropical", "temperate", "tropical"), each = length(m.temp_seq))))

m.plot_data.climate <- data.frame(
  Temperature = test$Temp,
  Trait = test$Trait,
  Climate = as.factor(test$climate),
  Species = as.factor(test$interactor1)
)

m.species_to_climate <- unique(m.plot_data.climate[, c("Species", "Climate")])  # Ensure species matches correct climate
m.species_to_climate$Climate <- as.factor(m.species_to_climate$Climate)
# Merge fitted climate curves with species-to-climate mapping
fitted_data.climate.krombeini <- data.frame(
  Temperature = rep(m.temp_seq, m.n.climate.rm),
  Fitted = c(clim1.rm(m.temp_seq),
             clim2.rm(m.temp_seq),
             clim3.rm(m.temp_seq)),
  Climate = factor(rep(c("subtropical", "temperate", "tropical"), each = length(m.temp_seq)))
) %>%
  left_join(m.species_to_climate, by = "Climate")

m.plot_data.climate$Climate <- factor(m.plot_data.climate$Climate, levels = c("temperate", "subtropical", "tropical"))
fitted_data.climate.krombeini$Climate <- factor(fitted_data.climate.krombeini$Climate, levels = c("temperate", "subtropical", "tropical"))
# Plot observed data and fitted curves with facet wrap by species
p.climate.krombeini <- ggplot() +
  geom_point(data = m.plot_data.climate, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.climate.krombeini, aes(x = Temperature, y = Fitted, color = Climate), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )   
p.climate.krombeini

# Merge the observed and predicted data for comparison
merged_data.climate.krombeini <- fitted_data.climate.krombeini %>%
  left_join(m.plot_data.climate, by = c("Temperature", "Climate", "Species"))


# Calculate RMSE for each climate
m.mse_by_climate.krombeini <- merged_data.climate.krombeini |> 
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.climate.krombeini <- merged_data.climate.krombeini |> 
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(m.mse_by_climate.krombeini,
     m.overall_rmse.climate.krombeini,
     fitted_data.climate.krombeini,
     merged_data.climate.krombeini,
     climate.krombeini.fit.mcmc,
     climate.krombeini.fit,
     p.climate.krombeini,
     file="krombeini.climate.RData")


######################################### By Climate- Notoscriptus #########################################
aedes.rm <- test[test$interactor1 != "Aedes notoscriptus", ]
m.climate.rm <- as.factor(aedes.rm$climate)
m.n.climate.rm <- length(unique(aedes.rm$climate))
m.N.obs.rm <- length(aedes.rm$Trait)
unique(aedes.rm$interactor1)
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

  for (j in 1:m.n.climate.rm) {
    mu.climate[j] ~ dexp(mu.c)       # Climate-specific low point
    mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific temp
    mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific temp
    climate.tau[j] ~ dexp(tau1)
  }

  # Likelihood
  for (i in 1:m.N.obs.rm) {
    mu[i] <- mu.climate.a[m.climate.rm[i]] * (mu.climate.Topt[m.climate.rm[i]] - temp[i])^2 + mu.climate[m.climate.rm[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[m.climate.rm[i]]) T(0, ) 
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
#jag.data<-list(trait = test$Trait, N.obs = N.obs, temp = test$Temp)

jag.data <- list(trait = aedes.rm$Trait, m.N.obs.rm = m.N.obs.rm, temp = aedes.rm$Temp, 
                 m.n.climate.rm = m.n.climate.rm, m.climate.rm = m.climate.rm)

climate.notoscriptus.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                            model.file="aedes.by.climate.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                            n.iter=ni, DIC=T, working.directory=getwd())
climate.notoscriptus.fit.mcmc <- as.mcmc(climate.notoscriptus.fit) ## makes an "mcmc" object
closeAllConnections()

# Get posterior means
summary_fit.rm.climate <- summary(climate.notoscriptus.fit.mcmc)$statistics
a_mean_climate.rm <- summary_fit.rm.climate[grep("^mu.climate.a\\[", rownames(summary_fit.rm.climate)), "Mean"]
c_mean_climate.rm <- summary_fit.rm.climate[grep("^mu.climate\\[", rownames(summary_fit.rm.climate)), "Mean"]
Topt_mean_climate.rm <- summary_fit.rm.climate[grep("^mu.climate.Topt\\[", rownames(summary_fit.rm.climate)), "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

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
fitted_data.climate.notoscriptus <- data.frame(Temperature = rep(m.temp_seq, m.n.climate.rm), 
                                          Fitted = c(clim1.rm(m.temp_seq),
                                                     clim2.rm(m.temp_seq),
                                                     clim3.rm(m.temp_seq)),
                                          Climate = factor(rep(c("subtropical", "temperate", "tropical"), each = length(m.temp_seq))))

m.plot_data.climate <- data.frame(
  Temperature = test$Temp,
  Trait = test$Trait,
  Climate = as.factor(test$climate),
  Species = as.factor(test$interactor1)
)

m.species_to_climate <- unique(m.plot_data.climate[, c("Species", "Climate")])  # Ensure species matches correct climate
m.species_to_climate$Climate <- as.factor(m.species_to_climate$Climate)
# Merge fitted climate curves with species-to-climate mapping
fitted_data.climate.notoscriptus <- data.frame(
  Temperature = rep(m.temp_seq, m.n.climate.rm),
  Fitted = c(clim1.rm(m.temp_seq),
             clim2.rm(m.temp_seq),
             clim3.rm(m.temp_seq)),
  Climate = factor(rep(c("subtropical", "temperate", "tropical"), each = length(m.temp_seq)))
) %>%
  left_join(m.species_to_climate, by = "Climate")

m.plot_data.climate$Climate <- factor(m.plot_data.climate$Climate, levels = c("temperate", "subtropical", "tropical"))
fitted_data.climate.notoscriptus$Climate <- factor(fitted_data.climate.notoscriptus$Climate, levels = c("temperate", "subtropical", "tropical"))
# Plot observed data and fitted curves with facet wrap by species
p.climate.notoscriptus <- ggplot() +
  geom_point(data = m.plot_data.climate, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.climate.notoscriptus, aes(x = Temperature, y = Fitted, color = Climate), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )   
p.climate.notoscriptus

# Merge the observed and predicted data for comparison
merged_data.climate.notoscriptus <- fitted_data.climate.notoscriptus %>%
  left_join(m.plot_data.climate, by = c("Temperature", "Climate", "Species"))


# Calculate RMSE for each climate
m.mse_by_climate.notoscriptus <- merged_data.climate.notoscriptus |> 
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.climate.notoscriptus <- merged_data.climate.notoscriptus |> 
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(m.mse_by_climate.notoscriptus,
     m.overall_rmse.climate.notoscriptus,
     fitted_data.climate.notoscriptus,
     merged_data.climate.notoscriptus,
     climate.notoscriptus.fit.mcmc,
     climate.notoscriptus.fit,
     p.climate.notoscriptus,
     file="notoscriptus.climate.RData")


######################################### By Wing Size Category- Aegypti #########################################
aedes.rm <- test[test$interactor1 != "Aedes aegypti", ]
m.size.rm <- as.factor(aedes.rm$wing.size)
m.n.size.rm <- length(unique(aedes.rm$wing.size))
m.N.obs.rm <- length(aedes.rm$Trait)
unique(aedes.rm$interactor1)
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

  for (j in 1:m.n.size.rm) {
    mu.size[j] ~ dexp(mu.c)       # Size-specific low point
    mu.size.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Size-specific temp
    mu.size.a[j] ~ dexp(mu.a)    # Size-specific temp
    size.tau[j] ~ dexp(tau1)
  }

  # Likelihood
  for (i in 1:m.N.obs.rm) {
    mu[i] <- mu.size.a[m.size.rm[i]] * (mu.size.Topt[m.size.rm[i]] - temp[i])^2 + mu.size[m.size.rm[i]]
    trait[i] ~ dnorm(mu[i], size.tau[m.size.rm[i]]) T(0, ) 
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

jag.data <- list(trait = aedes.rm$Trait, m.N.obs.rm = m.N.obs.rm, temp = aedes.rm$Temp, 
                 m.n.size.rm = m.n.size.rm, m.size.rm = m.size.rm)

size.aegypti.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                         model.file="aedes.by.size.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                         n.iter=ni, DIC=T, working.directory=getwd())
size.aegypti.fit.mcmc <- as.mcmc(size.aegypti.fit) ## makes an "mcmc" object
closeAllConnections()


# Get posterior means
summary_fit.size.rm <- summary(size.aegypti.fit.mcmc)$statistics
a_mean_size.rm <- summary_fit.size.rm[grep("^mu.size.a\\[", rownames(summary_fit.size.rm)), "Mean"]
c_mean_size.rm <- summary_fit.size.rm[grep("^mu.size\\[", rownames(summary_fit.size.rm)), "Mean"]
Topt_mean_size.rm <- summary_fit.size.rm[grep("^mu.size.Topt\\[", rownames(summary_fit.size.rm)), "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

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
fitted_data.size.aegypti <- data.frame(Temperature = rep(m.temp_seq, m.n.size.rm), 
                                       Fitted = c(size1.rm(m.temp_seq),
                                                  size2.rm(m.temp_seq),
                                                  size3.rm(m.temp_seq)),
                                       Size = factor(rep(c("large", "medium", "small"), each = length(m.temp_seq))))

m.plot_data.size <- data.frame(
  Temperature = test$Temp,
  Trait = test$Trait,
  Size = as.factor(test$wing.size),
  Species = as.factor(test$interactor1)
)

m.species_to_size <- unique(m.plot_data.size[, c("Species", "Size")])  # Ensure species matches correct climate
m.species_to_size$Size <- as.factor(m.species_to_size$Size)
# Merge fitted climate curves with species-to-climate mapping
fitted_data.size.aegypti <- data.frame(
  Temperature = rep(m.temp_seq, m.n.size.rm),
  Fitted = c(size1.rm(m.temp_seq),
             size2.rm(m.temp_seq),
             size3.rm(m.temp_seq)),
  Size = factor(rep(c("large", "medium", "small"), each = length(m.temp_seq)))
) %>%
  left_join(m.species_to_size, by = "Size")

m.plot_data.size$Size <- factor(m.plot_data.size$Size, levels = c("small", "medium", "large"))
fitted_data.size.aegypti$Size <- factor(fitted_data.size.aegypti$Size, levels = c("small", "medium", "large"))
# Plot observed data and fitted curves with facet wrap by species
p.size.aegypti <- ggplot() +
  geom_point(data = m.plot_data.size, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.size.aegypti, aes(x = Temperature, y = Fitted, color = Size), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )   
p.size.aegypti


# Merge the observed and predicted data for comparison
merged_data.size.aegypti <- fitted_data.size.aegypti %>%
  left_join(m.plot_data.size, by = c("Temperature", "Size", "Species"))


# Calculate RMSE for each species
mse_by_size.aegypti <- merged_data.size.aegypti %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.size.aegypti <- merged_data.size.aegypti |> 
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_size.aegypti,
     m.overall_rmse.size.aegypti,
     fitted_data.size.aegypti,
     merged_data.size.aegypti,
     size.aegypti.fit.mcmc,
     size.aegypti.fit,
     p.size.aegypti,
     file="aegypti.size.RData")


######################################### By Wing Size Category- Atropalpus #########################################
aedes.rm <- test[test$interactor1 != "Aedes atropalpus", ]
m.size.rm <- as.factor(aedes.rm$wing.size)
m.n.size.rm <- length(unique(aedes.rm$wing.size))
m.N.obs.rm <- length(aedes.rm$Trait)
unique(aedes.rm$interactor1)
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

  for (j in 1:m.n.size.rm) {
    mu.size[j] ~ dexp(mu.c)       # Size-specific low point
    mu.size.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Size-specific temp
    mu.size.a[j] ~ dexp(mu.a)    # Size-specific temp
    size.tau[j] ~ dexp(tau1)
  }

  # Likelihood
  for (i in 1:m.N.obs.rm) {
    mu[i] <- mu.size.a[m.size.rm[i]] * (mu.size.Topt[m.size.rm[i]] - temp[i])^2 + mu.size[m.size.rm[i]]
    trait[i] ~ dnorm(mu[i], size.tau[m.size.rm[i]]) T(0, ) 
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

jag.data <- list(trait = aedes.rm$Trait, m.N.obs.rm = m.N.obs.rm, temp = aedes.rm$Temp, 
                 m.n.size.rm = m.n.size.rm, m.size.rm = m.size.rm)

size.atropalpus.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                         model.file="aedes.by.size.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                         n.iter=ni, DIC=T, working.directory=getwd())
size.atropalpus.fit.mcmc <- as.mcmc(size.atropalpus.fit) ## makes an "mcmc" object
closeAllConnections()


# Get posterior means
summary_fit.size.rm <- summary(size.atropalpus.fit.mcmc)$statistics
a_mean_size.rm <- summary_fit.size.rm[grep("^mu.size.a\\[", rownames(summary_fit.size.rm)), "Mean"]
c_mean_size.rm <- summary_fit.size.rm[grep("^mu.size\\[", rownames(summary_fit.size.rm)), "Mean"]
Topt_mean_size.rm <- summary_fit.size.rm[grep("^mu.size.Topt\\[", rownames(summary_fit.size.rm)), "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

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
fitted_data.size.atropalpus <- data.frame(Temperature = rep(m.temp_seq, m.n.size.rm), 
                                       Fitted = c(size1.rm(m.temp_seq),
                                                  size2.rm(m.temp_seq),
                                                  size3.rm(m.temp_seq)),
                                       Size = factor(rep(c("large", "medium", "small"), each = length(m.temp_seq))))

m.plot_data.size <- data.frame(
  Temperature = test$Temp,
  Trait = test$Trait,
  Size = as.factor(test$wing.size),
  Species = as.factor(test$interactor1)
)

m.species_to_size <- unique(m.plot_data.size[, c("Species", "Size")])  # Ensure species matches correct climate
m.species_to_size$Size <- as.factor(m.species_to_size$Size)
# Merge fitted climate curves with species-to-climate mapping
fitted_data.size.atropalpus <- data.frame(
  Temperature = rep(m.temp_seq, m.n.size.rm),
  Fitted = c(size1.rm(m.temp_seq),
             size2.rm(m.temp_seq),
             size3.rm(m.temp_seq)),
  Size = factor(rep(c("large", "medium", "small"), each = length(m.temp_seq)))
) %>%
  left_join(m.species_to_size, by = "Size")

m.plot_data.size$Size <- factor(m.plot_data.size$Size, levels = c("small", "medium", "large"))
fitted_data.size.atropalpus$Size <- factor(fitted_data.size.atropalpus$Size, levels = c("small", "medium", "large"))
# Plot observed data and fitted curves with facet wrap by species
p.size.atropalpus <- ggplot() +
  geom_point(data = m.plot_data.size, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.size.atropalpus, aes(x = Temperature, y = Fitted, color = Size), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )   
p.size.atropalpus


# Merge the observed and predicted data for comparison
merged_data.size.atropalpus <- fitted_data.size.atropalpus %>%
  left_join(m.plot_data.size, by = c("Temperature", "Size", "Species"))


# Calculate RMSE for each species
mse_by_size.atropalpus <- merged_data.size.atropalpus %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.size.atropalpus <- merged_data.size.atropalpus |> 
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_size.atropalpus,
     m.overall_rmse.size.atropalpus,
     fitted_data.size.atropalpus,
     merged_data.size.atropalpus,
     size.atropalpus.fit.mcmc,
     size.atropalpus.fit,
     p.size.atropalpus,
     file="atropalpus.size.RData")


######################################### By Wing Size Category- Camptorhynchus #########################################
aedes.rm <- test[test$interactor1 != "Aedes camptorhynchus", ]
m.size.rm <- as.factor(aedes.rm$wing.size)
m.n.size.rm <- length(unique(aedes.rm$wing.size))
m.N.obs.rm <- length(aedes.rm$Trait)
unique(aedes.rm$interactor1)
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

  for (j in 1:m.n.size.rm) {
    mu.size[j] ~ dexp(mu.c)       # Size-specific low point
    mu.size.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Size-specific temp
    mu.size.a[j] ~ dexp(mu.a)    # Size-specific temp
    size.tau[j] ~ dexp(tau1)
  }

  # Likelihood
  for (i in 1:m.N.obs.rm) {
    mu[i] <- mu.size.a[m.size.rm[i]] * (mu.size.Topt[m.size.rm[i]] - temp[i])^2 + mu.size[m.size.rm[i]]
    trait[i] ~ dnorm(mu[i], size.tau[m.size.rm[i]]) T(0, ) 
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

jag.data <- list(trait = aedes.rm$Trait, m.N.obs.rm = m.N.obs.rm, temp = aedes.rm$Temp, 
                 m.n.size.rm = m.n.size.rm, m.size.rm = m.size.rm)

size.camptorhynchus.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                         model.file="aedes.by.size.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                         n.iter=ni, DIC=T, working.directory=getwd())
size.camptorhynchus.fit.mcmc <- as.mcmc(size.camptorhynchus.fit) ## makes an "mcmc" object
closeAllConnections()


# Get posterior means
summary_fit.size.rm <- summary(size.camptorhynchus.fit.mcmc)$statistics
a_mean_size.rm <- summary_fit.size.rm[grep("^mu.size.a\\[", rownames(summary_fit.size.rm)), "Mean"]
c_mean_size.rm <- summary_fit.size.rm[grep("^mu.size\\[", rownames(summary_fit.size.rm)), "Mean"]
Topt_mean_size.rm <- summary_fit.size.rm[grep("^mu.size.Topt\\[", rownames(summary_fit.size.rm)), "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

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
fitted_data.size.camptorhynchus <- data.frame(Temperature = rep(m.temp_seq, m.n.size.rm), 
                                       Fitted = c(size1.rm(m.temp_seq),
                                                  size2.rm(m.temp_seq),
                                                  size3.rm(m.temp_seq)),
                                       Size = factor(rep(c("large", "medium", "small"), each = length(m.temp_seq))))

m.plot_data.size <- data.frame(
  Temperature = test$Temp,
  Trait = test$Trait,
  Size = as.factor(test$wing.size),
  Species = as.factor(test$interactor1)
)

m.species_to_size <- unique(m.plot_data.size[, c("Species", "Size")])  # Ensure species matches correct climate
m.species_to_size$Size <- as.factor(m.species_to_size$Size)
# Merge fitted climate curves with species-to-climate mapping
fitted_data.size.camptorhynchus <- data.frame(
  Temperature = rep(m.temp_seq, m.n.size.rm),
  Fitted = c(size1.rm(m.temp_seq),
             size2.rm(m.temp_seq),
             size3.rm(m.temp_seq)),
  Size = factor(rep(c("large", "medium", "small"), each = length(m.temp_seq)))
) %>%
  left_join(m.species_to_size, by = "Size")

m.plot_data.size$Size <- factor(m.plot_data.size$Size, levels = c("small", "medium", "large"))
fitted_data.size.camptorhynchus$Size <- factor(fitted_data.size.camptorhynchus$Size, levels = c("small", "medium", "large"))
# Plot observed data and fitted curves with facet wrap by species
p.size.camptorhynchus <- ggplot() +
  geom_point(data = m.plot_data.size, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.size.camptorhynchus, aes(x = Temperature, y = Fitted, color = Size), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )   
p.size.camptorhynchus


# Merge the observed and predicted data for comparison
merged_data.size.camptorhynchus <- fitted_data.size.camptorhynchus %>%
  left_join(m.plot_data.size, by = c("Temperature", "Size", "Species"))


# Calculate RMSE for each species
mse_by_size.camptorhynchus <- merged_data.size.camptorhynchus %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.size.camptorhynchus <- merged_data.size.camptorhynchus |> 
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_size.camptorhynchus,
     m.overall_rmse.size.camptorhynchus,
     fitted_data.size.camptorhynchus,
     merged_data.size.camptorhynchus,
     size.camptorhynchus.fit.mcmc,
     size.camptorhynchus.fit,
     p.size.camptorhynchus,
     file="camptorhynchus.size.RData")


######################################### By Wing Size Category- Japonicus Japonicus #########################################
aedes.rm <- test[test$interactor1 != "Aedes japonicus japonicus", ]
m.size.rm <- as.factor(aedes.rm$wing.size)
m.n.size.rm <- length(unique(aedes.rm$wing.size))
m.N.obs.rm <- length(aedes.rm$Trait)
unique(aedes.rm$interactor1)
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

  for (j in 1:m.n.size.rm) {
    mu.size[j] ~ dexp(mu.c)       # Size-specific low point
    mu.size.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Size-specific temp
    mu.size.a[j] ~ dexp(mu.a)    # Size-specific temp
    size.tau[j] ~ dexp(tau1)
  }

  # Likelihood
  for (i in 1:m.N.obs.rm) {
    mu[i] <- mu.size.a[m.size.rm[i]] * (mu.size.Topt[m.size.rm[i]] - temp[i])^2 + mu.size[m.size.rm[i]]
    trait[i] ~ dnorm(mu[i], size.tau[m.size.rm[i]]) T(0, ) 
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

jag.data <- list(trait = aedes.rm$Trait, m.N.obs.rm = m.N.obs.rm, temp = aedes.rm$Temp, 
                 m.n.size.rm = m.n.size.rm, m.size.rm = m.size.rm)

size.jj.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                         model.file="aedes.by.size.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                         n.iter=ni, DIC=T, working.directory=getwd())
size.jj.fit.mcmc <- as.mcmc(size.jj.fit) ## makes an "mcmc" object
closeAllConnections()


# Get posterior means
summary_fit.size.rm <- summary(size.jj.fit.mcmc)$statistics
a_mean_size.rm <- summary_fit.size.rm[grep("^mu.size.a\\[", rownames(summary_fit.size.rm)), "Mean"]
c_mean_size.rm <- summary_fit.size.rm[grep("^mu.size\\[", rownames(summary_fit.size.rm)), "Mean"]
Topt_mean_size.rm <- summary_fit.size.rm[grep("^mu.size.Topt\\[", rownames(summary_fit.size.rm)), "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

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
fitted_data.size.jj <- data.frame(Temperature = rep(m.temp_seq, m.n.size.rm), 
                                       Fitted = c(size1.rm(m.temp_seq),
                                                  size2.rm(m.temp_seq),
                                                  size3.rm(m.temp_seq)),
                                       Size = factor(rep(c("large", "medium", "small"), each = length(m.temp_seq))))

m.plot_data.size <- data.frame(
  Temperature = test$Temp,
  Trait = test$Trait,
  Size = as.factor(test$wing.size),
  Species = as.factor(test$interactor1)
)

m.species_to_size <- unique(m.plot_data.size[, c("Species", "Size")])  # Ensure species matches correct climate
m.species_to_size$Size <- as.factor(m.species_to_size$Size)
# Merge fitted climate curves with species-to-climate mapping
fitted_data.size.jj <- data.frame(
  Temperature = rep(m.temp_seq, m.n.size.rm),
  Fitted = c(size1.rm(m.temp_seq),
             size2.rm(m.temp_seq),
             size3.rm(m.temp_seq)),
  Size = factor(rep(c("large", "medium", "small"), each = length(m.temp_seq)))
) %>%
  left_join(m.species_to_size, by = "Size")

m.plot_data.size$Size <- factor(m.plot_data.size$Size, levels = c("small", "medium", "large"))
fitted_data.size.jj$Size <- factor(fitted_data.size.jj$Size, levels = c("small", "medium", "large"))
# Plot observed data and fitted curves with facet wrap by species
p.size.jj <- ggplot() +
  geom_point(data = m.plot_data.size, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.size.jj, aes(x = Temperature, y = Fitted, color = Size), size = 1) +
  labs(x = "Temperature", y = "Develpoment Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )   
p.size.jj


# Merge the observed and predicted data for comparison
merged_data.size.jj <- fitted_data.size.jj %>%
  left_join(m.plot_data.size, by = c("Temperature", "Size", "Species"))


# Calculate RMSE for each species
mse_by_size.jj <- merged_data.size.jj %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.size.jj <- merged_data.size.jj |> 
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_size.jj,
     m.overall_rmse.size.jj,
     fitted_data.size.jj,
     merged_data.size.jj,
     size.jj.fit.mcmc,
     size.jj.fit,
     p.size.jj,
     file="jj.size.RData")


######################################### By Wing Size Category- Krombeini #########################################
aedes.rm <- test[test$interactor1 != "Aedes krombeini", ]
m.size.rm <- as.factor(aedes.rm$wing.size)
m.n.size.rm <- length(unique(aedes.rm$wing.size))
m.N.obs.rm <- length(aedes.rm$Trait)
unique(aedes.rm$interactor1)
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

  for (j in 1:m.n.size.rm) {
    mu.size[j] ~ dexp(mu.c)       # Size-specific low point
    mu.size.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Size-specific temp
    mu.size.a[j] ~ dexp(mu.a)    # Size-specific temp
    size.tau[j] ~ dexp(tau1)
  }

  # Likelihood
  for (i in 1:m.N.obs.rm) {
    mu[i] <- mu.size.a[m.size.rm[i]] * (mu.size.Topt[m.size.rm[i]] - temp[i])^2 + mu.size[m.size.rm[i]]
    trait[i] ~ dnorm(mu[i], size.tau[m.size.rm[i]]) T(0, ) 
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

jag.data <- list(trait = aedes.rm$Trait, m.N.obs.rm = m.N.obs.rm, temp = aedes.rm$Temp, 
                 m.n.size.rm = m.n.size.rm, m.size.rm = m.size.rm)

size.krombeini.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                         model.file="aedes.by.size.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                         n.iter=ni, DIC=T, working.directory=getwd())
size.krombeini.fit.mcmc <- as.mcmc(size.krombeini.fit) ## makes an "mcmc" object
closeAllConnections()


# Get posterior means
summary_fit.size.rm <- summary(size.krombeini.fit.mcmc)$statistics
a_mean_size.rm <- summary_fit.size.rm[grep("^mu.size.a\\[", rownames(summary_fit.size.rm)), "Mean"]
c_mean_size.rm <- summary_fit.size.rm[grep("^mu.size\\[", rownames(summary_fit.size.rm)), "Mean"]
Topt_mean_size.rm <- summary_fit.size.rm[grep("^mu.size.Topt\\[", rownames(summary_fit.size.rm)), "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

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
fitted_data.size.krombeini <- data.frame(Temperature = rep(m.temp_seq, m.n.size.rm), 
                                       Fitted = c(size1.rm(m.temp_seq),
                                                  size2.rm(m.temp_seq),
                                                  size3.rm(m.temp_seq)),
                                       Size = factor(rep(c("large", "medium", "small"), each = length(m.temp_seq))))

m.plot_data.size <- data.frame(
  Temperature = test$Temp,
  Trait = test$Trait,
  Size = as.factor(test$wing.size),
  Species = as.factor(test$interactor1)
)

m.species_to_size <- unique(m.plot_data.size[, c("Species", "Size")])  # Ensure species matches correct climate
m.species_to_size$Size <- as.factor(m.species_to_size$Size)
# Merge fitted climate curves with species-to-climate mapping
fitted_data.size.krombeini <- data.frame(
  Temperature = rep(m.temp_seq, m.n.size.rm),
  Fitted = c(size1.rm(m.temp_seq),
             size2.rm(m.temp_seq),
             size3.rm(m.temp_seq)),
  Size = factor(rep(c("large", "medium", "small"), each = length(m.temp_seq)))
) %>%
  left_join(m.species_to_size, by = "Size")

m.plot_data.size$Size <- factor(m.plot_data.size$Size, levels = c("small", "medium", "large"))
fitted_data.size.krombeini$Size <- factor(fitted_data.size.krombeini$Size, levels = c("small", "medium", "large"))
# Plot observed data and fitted curves with facet wrap by species
p.size.krombeini <- ggplot() +
  geom_point(data = m.plot_data.size, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.size.krombeini, aes(x = Temperature, y = Fitted, color = Size), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )   
p.size.krombeini


# Merge the observed and predicted data for comparison
merged_data.size.krombeini <- fitted_data.size.krombeini %>%
  left_join(m.plot_data.size, by = c("Temperature", "Size", "Species"))


# Calculate RMSE for each species
mse_by_size.krombeini <- merged_data.size.krombeini %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.size.krombeini <- merged_data.size.krombeini |> 
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_size.krombeini,
     m.overall_rmse.size.krombeini,
     fitted_data.size.krombeini,
     merged_data.size.krombeini,
     size.krombeini.fit.mcmc,
     size.krombeini.fit,
     p.size.krombeini,
     file="krombeini.size.RData")


######################################### By Wing Size Category- Notoscriptus #########################################
aedes.rm <- test[test$interactor1 != "Aedes notoscriptus", ]
m.size.rm <- as.factor(aedes.rm$wing.size)
m.n.size.rm <- length(unique(aedes.rm$wing.size))
m.N.obs.rm <- length(aedes.rm$Trait)
unique(aedes.rm$interactor1)
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

  for (j in 1:m.n.size.rm) {
    mu.size[j] ~ dexp(mu.c)       # Size-specific low point
    mu.size.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Size-specific temp
    mu.size.a[j] ~ dexp(mu.a)    # Size-specific temp
    size.tau[j] ~ dexp(tau1)
  }

  # Likelihood
  for (i in 1:m.N.obs.rm) {
    mu[i] <- mu.size.a[m.size.rm[i]] * (mu.size.Topt[m.size.rm[i]] - temp[i])^2 + mu.size[m.size.rm[i]]
    trait[i] ~ dnorm(mu[i], size.tau[m.size.rm[i]]) T(0, ) 
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

jag.data <- list(trait = aedes.rm$Trait, m.N.obs.rm = m.N.obs.rm, temp = aedes.rm$Temp, 
                 m.n.size.rm = m.n.size.rm, m.size.rm = m.size.rm)

size.notoscriptus.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                         model.file="aedes.by.size.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                         n.iter=ni, DIC=T, working.directory=getwd())
size.notoscriptus.fit.mcmc <- as.mcmc(size.notoscriptus.fit) ## makes an "mcmc" object
closeAllConnections()


# Get posterior means
summary_fit.size.rm <- summary(size.notoscriptus.fit.mcmc)$statistics
a_mean_size.rm <- summary_fit.size.rm[grep("^mu.size.a\\[", rownames(summary_fit.size.rm)), "Mean"]
c_mean_size.rm <- summary_fit.size.rm[grep("^mu.size\\[", rownames(summary_fit.size.rm)), "Mean"]
Topt_mean_size.rm <- summary_fit.size.rm[grep("^mu.size.Topt\\[", rownames(summary_fit.size.rm)), "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

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
fitted_data.size.notoscriptus <- data.frame(Temperature = rep(m.temp_seq, m.n.size.rm), 
                                       Fitted = c(size1.rm(m.temp_seq),
                                                  size2.rm(m.temp_seq),
                                                  size3.rm(m.temp_seq)),
                                       Size = factor(rep(c("large", "medium", "small"), each = length(m.temp_seq))))

m.plot_data.size <- data.frame(
  Temperature = test$Temp,
  Trait = test$Trait,
  Size = as.factor(test$wing.size),
  Species = as.factor(test$interactor1)
)

m.species_to_size <- unique(m.plot_data.size[, c("Species", "Size")])  # Ensure species matches correct climate
m.species_to_size$Size <- as.factor(m.species_to_size$Size)
# Merge fitted climate curves with species-to-climate mapping
fitted_data.size.notoscriptus <- data.frame(
  Temperature = rep(m.temp_seq, m.n.size.rm),
  Fitted = c(size1.rm(m.temp_seq),
             size2.rm(m.temp_seq),
             size3.rm(m.temp_seq)),
  Size = factor(rep(c("large", "medium", "small"), each = length(m.temp_seq)))
) %>%
  left_join(m.species_to_size, by = "Size")

m.plot_data.size$Size <- factor(m.plot_data.size$Size, levels = c("small", "medium", "large"))
fitted_data.size.notoscriptus$Size <- factor(fitted_data.size.notoscriptus$Size, levels = c("small", "medium", "large"))
# Plot observed data and fitted curves with facet wrap by species
p.size.notoscriptus <- ggplot() +
  geom_point(data = m.plot_data.size, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.size.notoscriptus, aes(x = Temperature, y = Fitted, color = Size), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )   
p.size.notoscriptus


# Merge the observed and predicted data for comparison
merged_data.size.notoscriptus <- fitted_data.size.notoscriptus %>%
  left_join(m.plot_data.size, by = c("Temperature", "Size", "Species"))


# Calculate RMSE for each species
mse_by_size.notoscriptus <- merged_data.size.notoscriptus %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.size.notoscriptus <- merged_data.size.notoscriptus |> 
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(mse_by_size.notoscriptus,
     m.overall_rmse.size.notoscriptus,
     fitted_data.size.notoscriptus,
     merged_data.size.notoscriptus,
     size.notoscriptus.fit.mcmc,
     size.notoscriptus.fit,
     p.size.notoscriptus,
     file="notoscriptus.size.RData")


######################################### By Wing Size and Climate- Aegypti #########################################
aedes.rm <- test[test$interactor1 != "Aedes aegypti", ]
m.N.obs.rm <- length(aedes.rm$Trait)
m.climate.rm <- as.factor(aedes.rm$climate)
m.n.climate.rm <- length(unique(aedes.rm$climate))
m.size.rm <- as.factor(aedes.rm$wing.size)
m.n.size.rm <- length(unique(aedes.rm$wing.size))
unique(aedes.rm$interactor1)

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

  for (j in 1:m.n.climate.rm) {
    mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific temp
    mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific temp
  }
    for (k in 1:m.n.size.rm) {
    mu.size[k] ~ dexp(mu.c)       # Size-specific low point
    size.tau[k] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:m.N.obs.rm) {
    mu[i] <- mu.climate.a[m.climate.rm[i]] * (mu.climate.Topt[m.climate.rm[i]] - temp[i])^2 + mu.size[m.size.rm[i]]
    trait[i] ~ dnorm(mu[i], size.tau[m.size.rm[i]])  T(0, )
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

jag.data<-list(trait = aedes.rm$Trait, m.N.obs.rm = m.N.obs.rm, temp = aedes.rm$Temp, 
               m.n.climate.rm = m.n.climate.rm, m.climate.rm = m.climate.rm, m.n.size.rm = m.n.size.rm, m.size.rm = m.size.rm)

m.climate.size.aegypti.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                                   model.file="aedes.by.climate.size.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                                   n.iter=ni, DIC=T, working.directory=getwd())
m.climate.size.aegypti.fit.mcmc <- as.mcmc(m.climate.size.aegypti.fit) ## makes an "mcmc" object
closeAllConnections()

# Get posterior means
m.summary_fit.climate.size <- summary(m.climate.size.aegypti.fit.mcmc)$statistics
m.a_mean_climate.size <- m.summary_fit.climate.size[grep("^mu.climate.a\\[", rownames(m.summary_fit.climate.size)), "Mean"]
m.c_mean_climate.size <- m.summary_fit.climate.size[grep("^mu.size\\[", rownames(m.summary_fit.climate.size)), "Mean"]
m.Topt_mean_climate.size <- m.summary_fit.climate.size[grep("^mu.climate.Topt\\[", rownames(m.summary_fit.climate.size)), "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

# Get predicted values
m.clim1size1 <- function (x){
  m.a_mean_climate.size[1] * (x - m.Topt_mean_climate.size[1])^2 + m.c_mean_climate.size[1]
}
m.clim2size1 <- function (x){
  m.a_mean_climate.size[2] * (x - m.Topt_mean_climate.size[2])^2 + m.c_mean_climate.size[1]
}
m.clim3size1 <- function (x){
  m.a_mean_climate.size[3] * (x - m.Topt_mean_climate.size[3])^2 + m.c_mean_climate.size[1]
}
m.clim1size2 <- function (x){
  m.a_mean_climate.size[1] * (x - m.Topt_mean_climate.size[1])^2 + m.c_mean_climate.size[2]
}
m.clim2size2 <- function (x){
  m.a_mean_climate.size[2] * (x - m.Topt_mean_climate.size[2])^2 + m.c_mean_climate.size[2]
}
m.clim3size2 <- function (x){
  m.a_mean_climate.size[3] * (x - m.Topt_mean_climate.size[3])^2 + m.c_mean_climate.size[2]
}
m.clim1size3 <- function (x){
  m.a_mean_climate.size[1] * (x - m.Topt_mean_climate.size[1])^2 + m.c_mean_climate.size[3]
}
m.clim2size3 <- function (x){
  m.a_mean_climate.size[2] * (x - m.Topt_mean_climate.size[2])^2 + m.c_mean_climate.size[3]
}
m.clim3size3 <- function (x){
  m.a_mean_climate.size[3] * (x - m.Topt_mean_climate.size[3])^2 + m.c_mean_climate.size[3]
}

# Optionally, convert to a data frame for easier plotting
m.fitted_data.climate.size.aegypti <- data.frame(Temperature = rep(m.temp_seq, m.n.climate.rm*3), 
                                                 Fitted = c(m.clim1size1(m.temp_seq),
                                                            m.clim2size1(m.temp_seq),
                                                            m.clim3size1(m.temp_seq),
                                                            m.clim1size2(m.temp_seq),
                                                            m.clim2size2(m.temp_seq),
                                                            m.clim3size2(m.temp_seq),
                                                            m.clim1size3(m.temp_seq),
                                                            m.clim2size3(m.temp_seq),
                                                            m.clim3size3(m.temp_seq)),
                                                 Climate = factor(rep(rep(c("subtropical", "temperate", "tropical"), each = length(m.temp_seq)), 3)),
                                                 Size = factor(rep(c("large", "medium", "small"), each = length(m.temp_seq) * 3)))

m.plot_data.climate.size <- data.frame(
  Temperature = test$Temp,
  Trait = test$Trait,
  Climate = as.factor(test$climate),
  Species = as.factor(test$interactor1),
  Size = as.factor(test$wing.size)
)

m.species_to_climate.size <- unique(m.plot_data.climate.size[, c("Species", "Climate", "Size")])
m.species_to_climate.size$Climate <- as.factor(m.species_to_climate.size$Climate)
m.species_to_climate.size$Size <- as.factor(m.species_to_climate.size$Size)

# Merge fitted climate curves with species-to-climate mapping and filter for matching species
m.fitted_data.climate.size.aegypti <- data.frame(
  Temperature = rep(m.temp_seq, m.n.climate.rm * m.n.size.rm), 
  Fitted = c(m.clim1size1(m.temp_seq),
             m.clim2size1(m.temp_seq),
             m.clim3size1(m.temp_seq),
             m.clim1size2(m.temp_seq),
             m.clim2size2(m.temp_seq),
             m.clim3size2(m.temp_seq),
             m.clim1size3(m.temp_seq),
             m.clim2size3(m.temp_seq),
             m.clim3size3(m.temp_seq)),
  Climate = factor(rep(rep(c("subtropical", "temperate", "tropical"), each = length(m.temp_seq)), 3)),
  Size = factor(rep(c("large", "medium", "small"), each = length(m.temp_seq) * 3))) |> 
  left_join(m.species_to_climate.size, by = c("Climate", "Size")) |> 
  filter(!is.na(Species))  # Keep only rows where Species is non-missing

# Combine Climate and Size for coloring
m.fitted_data.climate.size.aegypti$`Climate, Size` <- with(m.fitted_data.climate.size.aegypti, paste(Climate, Size, sep = ", "))
m.fitted_data.climate.size.aegypti$`Climate, Size` <- factor(m.fitted_data.climate.size.aegypti$`Climate, Size`, levels = c("temperate, large", 
                                                                                                        "subtropical, medium", "subtropical, large",
                                                                                                        "tropical, small", "tropical, medium"))
# Plot observed data and fitted curves with facet wrap by species
p.climate.size.aegypti <- ggplot() +
  geom_point(data = m.plot_data.climate.size, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.climate.size.aegypti, aes(x = Temperature, y = Fitted, color = `Climate, Size`), size = 1) +
  labs(x = "Temperature",
    y = "Development Time in Days"
  ) +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )  
p.climate.size.aegypti

# Merge the observed and predicted data for comparison
m.merged_data.climate.size.aegypti <- m.fitted_data.climate.size.aegypti %>%
  left_join(m.plot_data.climate.size, by = c("Temperature", "Climate", "Species"))


# Calculate RMSE for each species
m.mse_by_climate.size.aegypti <- m.merged_data.climate.size.aegypti %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.climate.size.aegypti <- m.merged_data.climate.size.aegypti |> 
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(m.mse_by_climate.size.aegypti,
     m.overall_rmse.climate.size.aegypti,
     m.fitted_data.climate.size.aegypti,
     m.merged_data.climate.size.aegypti,
     m.climate.size.aegypti.fit.mcmc,
     m.climate.size.aegypti.fit,
     p.climate.size.aegypti,
     file="aegypti.climate.size.RData")


######################################### By Wing Size and Climate- Atropalpus #########################################
aedes.rm <- test[test$interactor1 != "Aedes atropalpus", ]
m.N.obs.rm <- length(aedes.rm$Trait)
m.climate.rm <- as.factor(aedes.rm$climate)
m.n.climate.rm <- length(unique(aedes.rm$climate))
m.size.rm <- as.factor(aedes.rm$wing.size)
m.n.size.rm <- length(unique(aedes.rm$wing.size))
unique(aedes.rm$interactor1)

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

  for (j in 1:m.n.climate.rm) {
    mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific temp
    mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific temp
  }
    for (k in 1:m.n.size.rm) {
    mu.size[k] ~ dexp(mu.c)       # Size-specific low point
    size.tau[k] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:m.N.obs.rm) {
    mu[i] <- mu.climate.a[m.climate.rm[i]] * (mu.climate.Topt[m.climate.rm[i]] - temp[i])^2 + mu.size[m.size.rm[i]]
    trait[i] ~ dnorm(mu[i], size.tau[m.size.rm[i]])  T(0, )
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

jag.data<-list(trait = aedes.rm$Trait, m.N.obs.rm = m.N.obs.rm, temp = aedes.rm$Temp, 
               m.n.climate.rm = m.n.climate.rm, m.climate.rm = m.climate.rm, m.n.size.rm = m.n.size.rm, m.size.rm = m.size.rm)

m.climate.size.atropalpus.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                                   model.file="aedes.by.climate.size.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                                   n.iter=ni, DIC=T, working.directory=getwd())
m.climate.size.atropalpus.fit.mcmc <- as.mcmc(m.climate.size.atropalpus.fit) ## makes an "mcmc" object
closeAllConnections()

# Get posterior means
m.summary_fit.climate.size <- summary(m.climate.size.atropalpus.fit.mcmc)$statistics
m.a_mean_climate.size <- m.summary_fit.climate.size[grep("^mu.climate.a\\[", rownames(m.summary_fit.climate.size)), "Mean"]
m.c_mean_climate.size <- m.summary_fit.climate.size[grep("^mu.size\\[", rownames(m.summary_fit.climate.size)), "Mean"]
m.Topt_mean_climate.size <- m.summary_fit.climate.size[grep("^mu.climate.Topt\\[", rownames(m.summary_fit.climate.size)), "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

# Get predicted values
m.clim1size1 <- function (x){
  m.a_mean_climate.size[1] * (x - m.Topt_mean_climate.size[1])^2 + m.c_mean_climate.size[1]
}
m.clim2size1 <- function (x){
  m.a_mean_climate.size[2] * (x - m.Topt_mean_climate.size[2])^2 + m.c_mean_climate.size[1]
}
m.clim3size1 <- function (x){
  m.a_mean_climate.size[3] * (x - m.Topt_mean_climate.size[3])^2 + m.c_mean_climate.size[1]
}
m.clim1size2 <- function (x){
  m.a_mean_climate.size[1] * (x - m.Topt_mean_climate.size[1])^2 + m.c_mean_climate.size[2]
}
m.clim2size2 <- function (x){
  m.a_mean_climate.size[2] * (x - m.Topt_mean_climate.size[2])^2 + m.c_mean_climate.size[2]
}
m.clim3size2 <- function (x){
  m.a_mean_climate.size[3] * (x - m.Topt_mean_climate.size[3])^2 + m.c_mean_climate.size[2]
}
m.clim1size3 <- function (x){
  m.a_mean_climate.size[1] * (x - m.Topt_mean_climate.size[1])^2 + m.c_mean_climate.size[3]
}
m.clim2size3 <- function (x){
  m.a_mean_climate.size[2] * (x - m.Topt_mean_climate.size[2])^2 + m.c_mean_climate.size[3]
}
m.clim3size3 <- function (x){
  m.a_mean_climate.size[3] * (x - m.Topt_mean_climate.size[3])^2 + m.c_mean_climate.size[3]
}

# Optionally, convert to a data frame for easier plotting
m.fitted_data.climate.size.atropalpus <- data.frame(Temperature = rep(m.temp_seq, m.n.climate.rm*3), 
                                                 Fitted = c(m.clim1size1(m.temp_seq),
                                                            m.clim2size1(m.temp_seq),
                                                            m.clim3size1(m.temp_seq),
                                                            m.clim1size2(m.temp_seq),
                                                            m.clim2size2(m.temp_seq),
                                                            m.clim3size2(m.temp_seq),
                                                            m.clim1size3(m.temp_seq),
                                                            m.clim2size3(m.temp_seq),
                                                            m.clim3size3(m.temp_seq)),
                                                 Climate = factor(rep(rep(c("subtropical", "temperate", "tropical"), each = length(m.temp_seq)), 3)),
                                                 Size = factor(rep(c("large", "medium", "small"), each = length(m.temp_seq) * 3)))

m.plot_data.climate.size <- data.frame(
  Temperature = test$Temp,
  Trait = test$Trait,
  Climate = as.factor(test$climate),
  Species = as.factor(test$interactor1),
  Size = as.factor(test$wing.size)
)

m.species_to_climate.size <- unique(m.plot_data.climate.size[, c("Species", "Climate", "Size")])
m.species_to_climate.size$Climate <- as.factor(m.species_to_climate.size$Climate)
m.species_to_climate.size$Size <- as.factor(m.species_to_climate.size$Size)

# Merge fitted climate curves with species-to-climate mapping and filter for matching species
m.fitted_data.climate.size.atropalpus <- data.frame(
  Temperature = rep(m.temp_seq, m.n.climate.rm * m.n.size.rm), 
  Fitted = c(m.clim1size1(m.temp_seq),
             m.clim2size1(m.temp_seq),
             m.clim3size1(m.temp_seq),
             m.clim1size2(m.temp_seq),
             m.clim2size2(m.temp_seq),
             m.clim3size2(m.temp_seq),
             m.clim1size3(m.temp_seq),
             m.clim2size3(m.temp_seq),
             m.clim3size3(m.temp_seq)),
  Climate = factor(rep(rep(c("subtropical", "temperate", "tropical"), each = length(m.temp_seq)), 3)),
  Size = factor(rep(c("large", "medium", "small"), each = length(m.temp_seq) * 3))) |> 
  left_join(m.species_to_climate.size, by = c("Climate", "Size")) |> 
  filter(!is.na(Species))  # Keep only rows where Species is non-missing

# Combine Climate and Size for coloring
m.fitted_data.climate.size.atropalpus$`Climate, Size` <- with(m.fitted_data.climate.size.atropalpus, paste(Climate, Size, sep = ", "))
m.fitted_data.climate.size.atropalpus$`Climate, Size` <- factor(m.fitted_data.climate.size.atropalpus$`Climate, Size`, levels = c("temperate, large", 
                                                                                                        "subtropical, medium", "subtropical, large",
                                                                                                        "tropical, small", "tropical, medium"))
# Plot observed data and fitted curves with facet wrap by species
p.climate.size.atropalpus <- ggplot() +
  geom_point(data = m.plot_data.climate.size, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.climate.size.atropalpus, aes(x = Temperature, y = Fitted, color = `Climate, Size`), size = 1) +
  labs(x = "Temperature",
    y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )  
p.climate.size.atropalpus

# Merge the observed and predicted data for comparison
m.merged_data.climate.size.atropalpus <- m.fitted_data.climate.size.atropalpus %>%
  left_join(m.plot_data.climate.size, by = c("Temperature", "Climate", "Species"))


# Calculate RMSE for each species
m.mse_by_climate.size.atropalpus <- m.merged_data.climate.size.atropalpus %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.climate.size.atropalpus <- m.merged_data.climate.size.atropalpus |> 
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(m.mse_by_climate.size.atropalpus,
     m.overall_rmse.climate.size.atropalpus,
     m.fitted_data.climate.size.atropalpus,
     m.merged_data.climate.size.atropalpus,
     m.climate.size.atropalpus.fit.mcmc,
     m.climate.size.atropalpus.fit,
     p.climate.size.atropalpus,
     file="atropalpus.climate.size.RData")


######################################### By Wing Size and Climate- Japonicus Japonicus #########################################
aedes.rm <- test[test$interactor1 != "Aedes Japonicus Japonicus", ]
m.N.obs.rm <- length(aedes.rm$Trait)
m.climate.rm <- as.factor(aedes.rm$climate)
m.n.climate.rm <- length(unique(aedes.rm$climate))
m.size.rm <- as.factor(aedes.rm$wing.size)
m.n.size.rm <- length(unique(aedes.rm$wing.size))
unique(aedes.rm$interactor1)

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

  for (j in 1:m.n.climate.rm) {
    mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific temp
    mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific temp
  }
    for (k in 1:m.n.size.rm) {
    mu.size[k] ~ dexp(mu.c)       # Size-specific low point
    size.tau[k] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:m.N.obs.rm) {
    mu[i] <- mu.climate.a[m.climate.rm[i]] * (mu.climate.Topt[m.climate.rm[i]] - temp[i])^2 + mu.size[m.size.rm[i]]
    trait[i] ~ dnorm(mu[i], size.tau[m.size.rm[i]])  T(0, )
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

jag.data<-list(trait = aedes.rm$Trait, m.N.obs.rm = m.N.obs.rm, temp = aedes.rm$Temp, 
               m.n.climate.rm = m.n.climate.rm, m.climate.rm = m.climate.rm, m.n.size.rm = m.n.size.rm, m.size.rm = m.size.rm)

m.climate.size.jj.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                                   model.file="aedes.by.climate.size.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                                   n.iter=ni, DIC=T, working.directory=getwd())
m.climate.size.jj.fit.mcmc <- as.mcmc(m.climate.size.jj.fit) ## makes an "mcmc" object
closeAllConnections()

# Get posterior means
m.summary_fit.climate.size <- summary(m.climate.size.jj.fit.mcmc)$statistics
m.a_mean_climate.size <- m.summary_fit.climate.size[grep("^mu.climate.a\\[", rownames(m.summary_fit.climate.size)), "Mean"]
m.c_mean_climate.size <- m.summary_fit.climate.size[grep("^mu.size\\[", rownames(m.summary_fit.climate.size)), "Mean"]
m.Topt_mean_climate.size <- m.summary_fit.climate.size[grep("^mu.climate.Topt\\[", rownames(m.summary_fit.climate.size)), "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

# Get predicted values
m.clim1size1 <- function (x){
  m.a_mean_climate.size[1] * (x - m.Topt_mean_climate.size[1])^2 + m.c_mean_climate.size[1]
}
m.clim2size1 <- function (x){
  m.a_mean_climate.size[2] * (x - m.Topt_mean_climate.size[2])^2 + m.c_mean_climate.size[1]
}
m.clim3size1 <- function (x){
  m.a_mean_climate.size[3] * (x - m.Topt_mean_climate.size[3])^2 + m.c_mean_climate.size[1]
}
m.clim1size2 <- function (x){
  m.a_mean_climate.size[1] * (x - m.Topt_mean_climate.size[1])^2 + m.c_mean_climate.size[2]
}
m.clim2size2 <- function (x){
  m.a_mean_climate.size[2] * (x - m.Topt_mean_climate.size[2])^2 + m.c_mean_climate.size[2]
}
m.clim3size2 <- function (x){
  m.a_mean_climate.size[3] * (x - m.Topt_mean_climate.size[3])^2 + m.c_mean_climate.size[2]
}
m.clim1size3 <- function (x){
  m.a_mean_climate.size[1] * (x - m.Topt_mean_climate.size[1])^2 + m.c_mean_climate.size[3]
}
m.clim2size3 <- function (x){
  m.a_mean_climate.size[2] * (x - m.Topt_mean_climate.size[2])^2 + m.c_mean_climate.size[3]
}
m.clim3size3 <- function (x){
  m.a_mean_climate.size[3] * (x - m.Topt_mean_climate.size[3])^2 + m.c_mean_climate.size[3]
}

# Optionally, convert to a data frame for easier plotting
m.fitted_data.climate.size.jj <- data.frame(Temperature = rep(m.temp_seq, m.n.climate.rm*3), 
                                                 Fitted = c(m.clim1size1(m.temp_seq),
                                                            m.clim2size1(m.temp_seq),
                                                            m.clim3size1(m.temp_seq),
                                                            m.clim1size2(m.temp_seq),
                                                            m.clim2size2(m.temp_seq),
                                                            m.clim3size2(m.temp_seq),
                                                            m.clim1size3(m.temp_seq),
                                                            m.clim2size3(m.temp_seq),
                                                            m.clim3size3(m.temp_seq)),
                                                 Climate = factor(rep(rep(c("subtropical", "temperate", "tropical"), each = length(m.temp_seq)), 3)),
                                                 Size = factor(rep(c("large", "medium", "small"), each = length(m.temp_seq) * 3)))

m.plot_data.climate.size <- data.frame(
  Temperature = test$Temp,
  Trait = test$Trait,
  Climate = as.factor(test$climate),
  Species = as.factor(test$interactor1),
  Size = as.factor(test$wing.size)
)

m.species_to_climate.size <- unique(m.plot_data.climate.size[, c("Species", "Climate", "Size")])
m.species_to_climate.size$Climate <- as.factor(m.species_to_climate.size$Climate)
m.species_to_climate.size$Size <- as.factor(m.species_to_climate.size$Size)

# Merge fitted climate curves with species-to-climate mapping and filter for matching species
m.fitted_data.climate.size.jj <- data.frame(
  Temperature = rep(m.temp_seq, m.n.climate.rm * m.n.size.rm), 
  Fitted = c(m.clim1size1(m.temp_seq),
             m.clim2size1(m.temp_seq),
             m.clim3size1(m.temp_seq),
             m.clim1size2(m.temp_seq),
             m.clim2size2(m.temp_seq),
             m.clim3size2(m.temp_seq),
             m.clim1size3(m.temp_seq),
             m.clim2size3(m.temp_seq),
             m.clim3size3(m.temp_seq)),
  Climate = factor(rep(rep(c("subtropical", "temperate", "tropical"), each = length(m.temp_seq)), 3)),
  Size = factor(rep(c("large", "medium", "small"), each = length(m.temp_seq) * 3))) |> 
  left_join(m.species_to_climate.size, by = c("Climate", "Size")) |> 
  filter(!is.na(Species))  # Keep only rows where Species is non-missing

# Combine Climate and Size for coloring
m.fitted_data.climate.size.jj$`Climate, Size` <- with(m.fitted_data.climate.size.jj, paste(Climate, Size, sep = ", "))
m.fitted_data.climate.size.jj$`Climate, Size` <- factor(m.fitted_data.climate.size.jj$`Climate, Size`, levels = c("temperate, large", 
                                                                                                        "subtropical, medium", "subtropical, large",
                                                                                                        "tropical, small", "tropical, medium"))
# Plot observed data and fitted curves with facet wrap by species
p.climate.size.jj <- ggplot() +
  geom_point(data = m.plot_data.climate.size, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.climate.size.jj, aes(x = Temperature, y = Fitted, color = `Climate, Size`), size = 1) +
  labs(x = "Temperature",
    y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )  
p.climate.size.jj

# Merge the observed and predicted data for comparison
m.merged_data.climate.size.jj <- m.fitted_data.climate.size.jj %>%
  left_join(m.plot_data.climate.size, by = c("Temperature", "Climate", "Species"))


# Calculate RMSE for each species
m.mse_by_climate.size.jj <- m.merged_data.climate.size.jj %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.climate.size.jj <- m.merged_data.climate.size.jj |> 
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(m.mse_by_climate.size.jj,
     m.overall_rmse.climate.size.jj,
     m.fitted_data.climate.size.jj,
     m.merged_data.climate.size.jj,
     m.climate.size.jj.fit.mcmc,
     m.climate.size.jj.fit,
     p.climate.size.jj,
     file="jj.climate.size.RData")


######################################### By Wing Size and Climate- Krombeini #########################################
aedes.rm <- test[test$interactor1 != "Aedes krombeini", ]
m.N.obs.rm <- length(aedes.rm$Trait)
m.climate.rm <- as.factor(aedes.rm$climate)
m.n.climate.rm <- length(unique(aedes.rm$climate))
m.size.rm <- as.factor(aedes.rm$wing.size)
m.n.size.rm <- length(unique(aedes.rm$wing.size))
unique(aedes.rm$interactor1)

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

  for (j in 1:m.n.climate.rm) {
    mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific temp
    mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific temp
  }
    for (k in 1:m.n.size.rm) {
    mu.size[k] ~ dexp(mu.c)       # Size-specific low point
    size.tau[k] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:m.N.obs.rm) {
    mu[i] <- mu.climate.a[m.climate.rm[i]] * (mu.climate.Topt[m.climate.rm[i]] - temp[i])^2 + mu.size[m.size.rm[i]]
    trait[i] ~ dnorm(mu[i], size.tau[m.size.rm[i]])  T(0, )
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

jag.data<-list(trait = aedes.rm$Trait, m.N.obs.rm = m.N.obs.rm, temp = aedes.rm$Temp, 
               m.n.climate.rm = m.n.climate.rm, m.climate.rm = m.climate.rm, m.n.size.rm = m.n.size.rm, m.size.rm = m.size.rm)

m.climate.size.krombeini.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                                   model.file="aedes.by.climate.size.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                                   n.iter=ni, DIC=T, working.directory=getwd())
m.climate.size.krombeini.fit.mcmc <- as.mcmc(m.climate.size.krombeini.fit) ## makes an "mcmc" object
closeAllConnections()

# Get posterior means
m.summary_fit.climate.size <- summary(m.climate.size.krombeini.fit.mcmc)$statistics
m.a_mean_climate.size <- m.summary_fit.climate.size[grep("^mu.climate.a\\[", rownames(m.summary_fit.climate.size)), "Mean"]
m.c_mean_climate.size <- m.summary_fit.climate.size[grep("^mu.size\\[", rownames(m.summary_fit.climate.size)), "Mean"]
m.Topt_mean_climate.size <- m.summary_fit.climate.size[grep("^mu.climate.Topt\\[", rownames(m.summary_fit.climate.size)), "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

# Get predicted values
m.clim1size1 <- function (x){
  m.a_mean_climate.size[1] * (x - m.Topt_mean_climate.size[1])^2 + m.c_mean_climate.size[1]
}
m.clim2size1 <- function (x){
  m.a_mean_climate.size[2] * (x - m.Topt_mean_climate.size[2])^2 + m.c_mean_climate.size[1]
}
m.clim3size1 <- function (x){
  m.a_mean_climate.size[3] * (x - m.Topt_mean_climate.size[3])^2 + m.c_mean_climate.size[1]
}
m.clim1size2 <- function (x){
  m.a_mean_climate.size[1] * (x - m.Topt_mean_climate.size[1])^2 + m.c_mean_climate.size[2]
}
m.clim2size2 <- function (x){
  m.a_mean_climate.size[2] * (x - m.Topt_mean_climate.size[2])^2 + m.c_mean_climate.size[2]
}
m.clim3size2 <- function (x){
  m.a_mean_climate.size[3] * (x - m.Topt_mean_climate.size[3])^2 + m.c_mean_climate.size[2]
}
m.clim1size3 <- function (x){
  m.a_mean_climate.size[1] * (x - m.Topt_mean_climate.size[1])^2 + m.c_mean_climate.size[3]
}
m.clim2size3 <- function (x){
  m.a_mean_climate.size[2] * (x - m.Topt_mean_climate.size[2])^2 + m.c_mean_climate.size[3]
}
m.clim3size3 <- function (x){
  m.a_mean_climate.size[3] * (x - m.Topt_mean_climate.size[3])^2 + m.c_mean_climate.size[3]
}

# Optionally, convert to a data frame for easier plotting
m.fitted_data.climate.size.krombeini <- data.frame(Temperature = rep(m.temp_seq, m.n.climate.rm*3), 
                                                 Fitted = c(m.clim1size1(m.temp_seq),
                                                            m.clim2size1(m.temp_seq),
                                                            m.clim3size1(m.temp_seq),
                                                            m.clim1size2(m.temp_seq),
                                                            m.clim2size2(m.temp_seq),
                                                            m.clim3size2(m.temp_seq),
                                                            m.clim1size3(m.temp_seq),
                                                            m.clim2size3(m.temp_seq),
                                                            m.clim3size3(m.temp_seq)),
                                                 Climate = factor(rep(rep(c("subtropical", "temperate", "tropical"), each = length(m.temp_seq)), 3)),
                                                 Size = factor(rep(c("large", "medium", "small"), each = length(m.temp_seq) * 3)))

m.plot_data.climate.size <- data.frame(
  Temperature = test$Temp,
  Trait = test$Trait,
  Climate = as.factor(test$climate),
  Species = as.factor(test$interactor1),
  Size = as.factor(test$wing.size)
)

m.species_to_climate.size <- unique(m.plot_data.climate.size[, c("Species", "Climate", "Size")])
m.species_to_climate.size$Climate <- as.factor(m.species_to_climate.size$Climate)
m.species_to_climate.size$Size <- as.factor(m.species_to_climate.size$Size)

# Merge fitted climate curves with species-to-climate mapping and filter for matching species
m.fitted_data.climate.size.krombeini <- data.frame(
  Temperature = rep(m.temp_seq, m.n.climate.rm * m.n.size.rm), 
  Fitted = c(m.clim1size1(m.temp_seq),
             m.clim2size1(m.temp_seq),
             m.clim3size1(m.temp_seq),
             m.clim1size2(m.temp_seq),
             m.clim2size2(m.temp_seq),
             m.clim3size2(m.temp_seq),
             m.clim1size3(m.temp_seq),
             m.clim2size3(m.temp_seq),
             m.clim3size3(m.temp_seq)),
  Climate = factor(rep(rep(c("subtropical", "temperate", "tropical"), each = length(m.temp_seq)), 3)),
  Size = factor(rep(c("large", "medium", "small"), each = length(m.temp_seq) * 3))) |> 
  left_join(m.species_to_climate.size, by = c("Climate", "Size")) |> 
  filter(!is.na(Species))  # Keep only rows where Species is non-missing

# Combine Climate and Size for coloring
m.fitted_data.climate.size.krombeini$`Climate, Size` <- with(m.fitted_data.climate.size.krombeini, paste(Climate, Size, sep = ", "))
m.fitted_data.climate.size.krombeini$`Climate, Size` <- factor(m.fitted_data.climate.size.krombeini$`Climate, Size`, levels = c("temperate, large", 
                                                                                                        "subtropical, medium", "subtropical, large",
                                                                                                        "tropical, small", "tropical, medium"))
# Plot observed data and fitted curves with facet wrap by species
p.climate.size.krombeini <- ggplot() +
  geom_point(data = m.plot_data.climate.size, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.climate.size.krombeini, aes(x = Temperature, y = Fitted, color = `Climate, Size`), size = 1) +
  labs(x = "Temperature",
    y = "Development Time in Days") +
  theme_minimal() +
  facet_wrap(~ Species) + 
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )  
p.climate.size.krombeini

# Merge the observed and predicted data for comparison
m.merged_data.climate.size.krombeini <- m.fitted_data.climate.size.krombeini %>%
  left_join(m.plot_data.climate.size, by = c("Temperature", "Climate", "Species"))


# Calculate RMSE for each species
m.mse_by_climate.size.krombeini <- m.merged_data.climate.size.krombeini %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.climate.size.krombeini <- m.merged_data.climate.size.krombeini |> 
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(m.mse_by_climate.size.krombeini,
     m.overall_rmse.climate.size.krombeini,
     m.fitted_data.climate.size.krombeini,
     m.merged_data.climate.size.krombeini,
     m.climate.size.krombeini.fit.mcmc,
     m.climate.size.krombeini.fit,
     p.climate.size.krombeini,
     file="krombeini.climate.size.RData")



######################################### RMSE ###################################
# Data frame
table.mosquito.rm <- data.frame(
  Method = factor(c(rep("By Climate (RM)", 7 * 7), 
                    rep("By Wing Size (RM)", 7 * 6), 
                    rep("By Wing Size and Climate (RM)", 7 * 4))),
  Species = factor(rep(c("Aedes aegypti", "Aedes albopictus",
                         "Aedes atropalpus", "Aedes camptorhynchus",
                         "Aedes japonicus japonicus", "Aedes krombeini",
                         "Aedes notoscriptus"), times = 17)),
  MSE = c(m.mse_by_climate.aegypti, m.mse_by_climate.albopictus,
          m.mse_by_climate.atropalpus,
          m.mse_by_climate.camptorhynchus, m.mse_by_climate.jj,
          m.mse_by_climate.krombeini, m.mse_by_climate.notoscriptus,
          mse_by_size.aegypti, mse_by_size.atropalpus,
          mse_by_size.camptorhynchus, mse_by_size.jj,
          mse_by_size.krombeini, mse_by_size.notoscriptus,
          m.mse_by_climate.size.aegypti,
          m.mse_by_climate.size.atropalpus,
          m.mse_by_climate.size.jj,
          m.mse_by_climate.size.krombeini),
  Total.RMSE = rep(as.numeric(c(m.overall_rmse.climate.aegypti, m.overall_rmse.climate.albopictus,
                                m.overall_rmse.climate.atropalpus,
                                m.overall_rmse.climate.camptorhynchus ,m.overall_rmse.climate.jj,
                                m.overall_rmse.climate.krombeini, m.overall_rmse.climate.notoscriptus,
                                m.overall_rmse.size.aegypti, m.overall_rmse.size.atropalpus,
                                m.overall_rmse.size.camptorhynchus,m.overall_rmse.size.jj,
                                m.overall_rmse.size.krombeini, m.overall_rmse.size.notoscriptus,
                                m.overall_rmse.climate.size.aegypti,
                                m.overall_rmse.climate.size.atropalpus,
                                m.overall_rmse.climate.size.jj,
                                m.overall_rmse.climate.size.krombeini)), each = 7),
  Removed = c(rep(c("Aedes aegypti", "Aedes albopictus",
                    "Aedes atropalpus", "Aedes camptorhynchus",
                    "Aedes japonicus japonicus", "Aedes krombeini",
                    "Aedes notoscriptus"), each = 7),
              rep(c("Aedes aegypti",
                    "Aedes atropalpus", "Aedes camptorhynchus",
                    "Aedes japonicus japonicus", "Aedes krombeini",
                    "Aedes notoscriptus"), each = 7),
              rep(c("Aedes aegypti", 
                    "Aedes atropalpus",
                    "Aedes japonicus japonicus", 
                    "Aedes krombeini"), each = 7))
)

mosquito.DIC.rm <- data.frame(
  Method = rep(c("By Climate", "By Wing Size",
                 "By Climate and Wing Size"), each = 4),
  Removed = rep(c("Aedes aegypti", "Aedes atropalpus",
                  "Aedes japonicus japonicus", "Aedes krombeini"), times = 3),
  DIC = c(climate.aegypti.fit$BUGSoutput$DIC, climate.atropalpus.fit$BUGSoutput$DIC,
          climate.jj.fit$BUGSoutput$DIC, climate.notoscriptus.fit$BUGSoutput$DIC,
          size.aegypti.fit$BUGSoutput$DIC, size.atropalpus.fit$BUGSoutput$DIC,
          size.jj.fit$BUGSoutput$DIC, size.notoscriptus.fit$BUGSoutput$DIC,
          m.climate.size.aegypti.fit$BUGSoutput$DIC, m.climate.size.atropalpus.fit$BUGSoutput$DIC,
          m.climate.size.jj.fit$BUGSoutput$DIC, m.climate.size.krombeini.fit$BUGSoutput$DIC)
)

table.mosquito.rm |> 
  ggplot(aes(x = Removed, y = Total.RMSE, fill = Method))+
  geom_bar(stat = "identity", position = "dodge",
           color = "black")+
  theme_minimal() +
  scale_fill_viridis_d(option = "mako") +
  coord_cartesian(ylim = c(3.5, 5.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

mosquito.DIC.rm |> 
  ggplot(aes(x = Removed, y = DIC, fill = Method))+
  geom_bar(stat = "identity", position = "dodge", color = "black")+
  theme_minimal() +
  scale_fill_viridis_d(option = "mako") +
  coord_cartesian(ylim = c(2000, 6200)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save(table.mosquito.rm,
     mosquito.DIC.rm, file = "Dataset.Mosquito.RM.RData")








#########################################
#########################################
#########################################
######################################### By Wing Size Category Removed- Out of sample prediction #########################################
aedes.rm <- test[!(test$interactor1 %in% c("Aedes notoscriptus", "Aedes krombeini", "Aedes japonicus japonicus")), ]
size.rm <- as.factor(aedes.rm$wing.size)
n.size.rm <- length(unique(aedes.rm$wing.size))
N.obs.rm <- length(aedes.rm$Trait)
unique(aedes.rm$interactor1)
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

  for (j in 1:n.size.rm) {
    mu.size[j] ~ dexp(mu.c)       # Size-specific low point
    mu.size.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Size-specific temp
    mu.size.a[j] ~ dexp(mu.a)    # Size-specific temp
  }

  # Likelihood
  for (i in 1:N.obs.rm) {
    mu[i] <- mu.size.a[size.rm[i]] * (mu.size.Topt[size.rm[i]] - temp[i])^2 + mu.size[size.rm[i]]
    trait[i] ~ dnorm(mu[i], tau1) T(0, ) 
  }}
", file = "aedes.by.size.txt")


# Parameters to monitor
parameters <- c("mu.a", "mu.c", "mu.size", 
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
#jag.data<-list(trait = test$Trait, N.obs = N.obs, temp = test$Temp)

jag.data <- list(trait = aedes.rm$Trait, N.obs.rm = N.obs.rm, temp = aedes.rm$Temp, 
                 n.size.rm = n.size.rm, size.rm = size.rm)

size.rm.jj.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                       model.file="aedes.by.size.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                       n.iter=ni, DIC=T, working.directory=getwd())
size.rm.jj.fit.mcmc <- as.mcmc(size.rm.jj.fit) ## makes an "mcmc" object
closeAllConnections()

# Get posterior means
summary_fit.size.rm <- summary(size.rm.jj.fit.mcmc)$statistics
a_mean_size.rm <- summary_fit.size.rm[grep("^mu.size.a\\[", rownames(summary_fit.size.rm)), "Mean"]
c_mean_size.rm <- summary_fit.size.rm[grep("^mu.size\\[", rownames(summary_fit.size.rm)), "Mean"]
Topt_mean_size.rm <- summary_fit.size.rm[grep("^mu.size.Topt\\[", rownames(summary_fit.size.rm)), "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

# Get predicted values
m.size1.rm <- function (x){
  a_mean_size.rm[1] * (x - Topt_mean_size.rm[1])^2 + c_mean_size.rm[1]
}
m.size2.rm <- function (x){
  a_mean_size.rm[2] * (x - Topt_mean_size.rm[2])^2 + c_mean_size.rm[2]
}
m.size3.rm <- function (x){
  a_mean_size.rm[3] * (x - Topt_mean_size.rm[3])^2 + c_mean_size.rm[3]
}


# Optionally, convert to a data frame for easier plotting
fitted_data.size.rm.jj <- data.frame(Temperature = rep(m.temp_seq, n.size.rm), 
                                     Fitted = c(m.size1.rm(m.temp_seq),
                                                m.size2.rm(m.temp_seq),
                                                m.size3.rm(m.temp_seq)),
                                     Size = factor(rep(c("large", "medium", "small"), each = length(temp_seq))))

m.plot_data.size <- data.frame(
  Temperature = test$Temp,
  Trait = test$Trait,
  Size = as.factor(test$wing.size),
  Species = as.factor(test$interactor1)
)

m.species_to_size <- unique(m.plot_data.size[, c("Species", "Size")])  # Ensure species matches correct climate
m.species_to_size$Size <- as.factor(m.species_to_size$Size)
# Merge fitted climate curves with species-to-climate mapping
fitted_data.size.rm.jj <- data.frame(
  Temperature = rep(m.temp_seq, n.size.rm),
  Fitted = c(m.size1.rm(m.temp_seq),
             m.size2.rm(m.temp_seq),
             m.size3.rm(m.temp_seq)),
  Size = factor(rep(c("large", "medium", "small"), each = length(temp_seq)))
) %>%
  left_join(m.species_to_size, by = "Size")

m.plot_data.size$Size <- factor(m.plot_data.size$Size, levels = c("small", "medium", "large"))
fitted_data.size.rm.jj$Size <- factor(fitted_data.size.rm.jj$Size, levels = c("small", "medium", "large"))
# Plot observed data and fitted curves with facet wrap by species
p.size.rm.jj <- ggplot() +
  geom_point(data = m.plot_data.size, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.size.rm.jj, aes(x = Temperature, y = Fitted, color = Size), size = 1) +
  labs(title = "Development Time in Days by Size",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 
p.size.rm.jj


# Merge the observed and predicted data for comparison
merged_data.size.rm.jj <- fitted_data.size.rm.jj %>%
  left_join(m.plot_data.size, by = c("Temperature", "Size", "Species"))


# Calculate RMSE for each species
rmse_by_species.size.rm.jj <- merged_data.size.rm.jj %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_species.size.rm.jj)

mean(rmse_by_species.size.rm.jj$RMSE)

save(rmse_by_species.size.rm.jj,
     fitted_data.size.rm.jj,
     merged_data.size.rm.jj,
     size.rm.jj.fit.mcmc,
     size.rm.jj.fit,
     p.size.rm.jj,
     file="jj.size.rm.RData")




######################################### By Full Data- Out of sample prediction #########################################
aedes.rm <- test[test$interactor1 != "Aedes notoscriptus", ]
m.N.obs.rm <- length(aedes.rm$Trait)
unique(aedes.rm$interactor1)
sink("aedes.one.data.txt")
cat("
model {
  # Priors
   a ~ dexp(0.1)             # Width
  Topt ~ dnorm(30, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  c ~ dexp(0.001)
  
  # Likelihood
  for (i in 1:m.N.obs.rm) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], tau1) T(0, ) 
  }}
", file = "aedes.one.data.txt")

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
#jag.data<-list(trait = test$Trait, N.obs = N.obs, temp = test$Temp)

jag.data <- list(trait = aedes.rm$Trait, m.N.obs.rm = m.N.obs.rm, temp = aedes.rm$Temp)

one.data.notoscriptus.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                                  model.file="aedes.one.data.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                                  n.iter=ni, DIC=T, working.directory=getwd())
one.data.notoscriptus.fit.mcmc <- as.mcmc(one.data.notoscriptus.fit) ## makes an "mcmc" object
closeAllConnections()

# Get posterior means
summary_fit.rm <- summary(one.data.notoscriptus.fit.mcmc)$statistics
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
fitted_data.one.data.notoscriptus <- data.frame(Temperature = rep(m.temp_seq, m.n.species), 
                                                Fitted = full.data.rm(m.temp_seq),
                                                Species = factor(rep(mos.spec, each = length(m.temp_seq))))

plot_data.one.data <- data.frame(
  Temperature = test$Temp,
  Trait = test$Trait,
  Species = as.factor(test$interactor1)
)

# Plot observed data and fitted curves with facet wrap by species
p.notoscriptus <- ggplot() +
  geom_point(data = plot_data.one.data, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data.one.data.notoscriptus, aes(x = Temperature, y = Fitted), col = "blue", size = 1) +
  labs(title = "Development Time in Days by Species",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 
p.notoscriptus

# Merge the observed and predicted data for comparison
merged_data.one.data.notoscriptus<- fitted_data.one.data.notoscriptus %>%
  left_join(plot_data.one.data, by = c("Temperature", "Species"))


# Calculate RMSE for each species
rmse_by_species.one.data.notoscriptus <- merged_data.one.data.notoscriptus %>%
  group_by(Species) %>%
  summarize(
    RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))
  )

# View the results
print(rmse_by_species.one.data.notoscriptus)
mean(rmse_by_species.one.data.notoscriptus$RMSE)

save(rmse_by_species.one.data.notoscriptus,
     fitted_data.one.data.notoscriptus,
     merged_data.one.data.notoscriptus,
     one.data.notoscriptus.fit.mcmc,
     one.data.notoscriptus.fit,
     p.notoscriptus,
     file="notoscriptus.full.data.RData")


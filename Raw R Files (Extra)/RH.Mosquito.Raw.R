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
ggplot(test, aes(y = Trait, x = Temp)) +
  geom_point() +
  facet_wrap(~ interactor1) +
  labs(#title = "Temperature by Trait for Each Factor Level",
       x = "Temperature",
       y = "Development Time in Days") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),       # facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )  
ggsave("Hierarchical.m.data.png", width = 8, height = 8)


# Making vector for each species 
test$interactor1 <- as.factor(test$interactor1)
m.species <- as.factor(test$interactor1)
m.n.species <- length(unique(test$interactor1))

######################################### By Full Data #########################################
sink("aedes.one.data.txt")
cat("
model {
  # Priors
   a ~ dexp(0.1)             # Width
  Topt ~ dnorm(30, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  c ~ dexp(0.001)
  
  # Likelihood
  for (i in 1:m.N.obs) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], tau1) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], tau1)
  }}
", file = "aedes.one.data.txt")


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

jag.data<-list(trait = test$Trait, m.N.obs = m.N.obs, temp = test$Temp)

m.one.data.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                       model.file="aedes.one.data.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                       n.iter=ni, DIC=T, working.directory=getwd())
m.one.data.fit.mcmc <- as.mcmc(m.one.data.fit) ## makes an "mcmc" object
closeAllConnections()
m.one.data.fit$BUGSoutput$DIC

# Get rid of log liks
# m.one.data.mcmc <- lapply(as.mcmc.list(m.one.data.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(m.one.data.mcmc), autoburnin = FALSE)

# Extract log-likelihood matrix (samples × observations)
m.log.lik.full <- m.one.data.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
m.waic.full <- waic(m.log.lik.full)$estimates["waic", "Estimate"]


# Get posterior means
m.summary_fit <- summary(m.one.data.fit.mcmc)$statistics
m.a_mean <- m.summary_fit["a", "Mean"]
m.c_mean <- m.summary_fit["c", "Mean"]
m.means_Topt <- m.summary_fit["Topt", "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

# Get predicted values
m.full.data <- function (x){
  m.a_mean * (x - m.means_Topt)^2 + m.c_mean
}


# Optionally, convert to a data frame for easier plotting
m.fitted_data.one.data <- data.frame(Temperature = rep(m.temp_seq, m.n.species), 
                                     Fitted = m.full.data(m.temp_seq),
                                     Species = factor(rep(c("Aedes aegypti", "Aedes albopictus",
                                                            "Aedes atropalpus", "Aedes camptorhynchus",
                                                            "Aedes japonicus japonicus",
                                                            "Aedes krombeini", "Aedes notoscriptus"), each = length(m.temp_seq))))

m.plot_data.one.data <- data.frame(
  Temperature = test$Temp,
  Trait = test$Trait,
  Species = as.factor(test$interactor1)
)

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = m.plot_data.one.data, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.one.data, aes(x = Temperature, y = Fitted), color = "blue", size = 1) +
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
ggsave("mosquito.full.plot.png", width = 8, height = 5)


# Merge the observed and predicted data for comparison
m.merged_data.one.data <- m.fitted_data.one.data %>%
  left_join(plot_data.one.data, by = c("Temperature", "Species"))


# Calculate RMSE for each species
m.mse_by_species.one.data <- m.merged_data.one.data %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse <- m.merged_data.one.data %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(m.mse_by_species.one.data,
     m.fitted_data.one.data,
     m.merged_data.one.data,
     m.one.data.fit.mcmc,
     m.one.data.fit,
     m.log.lik.full,
     file="m.one.data.fit.RData")

######################################### By Species Individual- Aegypti #########################################
a.aegypti <- test[test$interactor1 == "Aedes aegypti", ]

m.N.obs.ind <- length(a.aegypti$Trait)

sink("aedes.individual.txt")
cat("
model {
  # Priors
   a ~ dexp(0.1)             # Width
  Topt ~ dnorm(30, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  c ~ dexp(0.001)
  
  # Likelihood
  for (i in 1:m.N.obs.ind) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], tau1) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], tau1)
  }}
", file = "aedes.individual.txt")


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

jag.data<-list(trait = a.aegypti$Trait, m.N.obs.ind = m.N.obs.ind, temp = a.aegypti$Temp)

a.aegypti.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                           model.file="aedes.individual.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                           n.iter=ni, DIC=T, working.directory=getwd())
a.aegypti.fit.mcmc <- as.mcmc(a.aegypti.fit) ## makes an "mcmc" object
closeAllConnections()
a.aegypti.fit$BUGSoutput$DIC

# Get rid of log liks
# aegypti.mcmc <- lapply(as.mcmc.list(a.aegypti.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(aegypti.mcmc), autoburnin = FALSE)


# Get posterior means
m.summary_fit.a.aegypti <- summary(a.aegypti.fit.mcmc)$statistics
m.a_mean <- m.summary_fit.a.aegypti["a", "Mean"]
m.c_mean <- m.summary_fit.a.aegypti["c", "Mean"]
m.means_Topt <- m.summary_fit.a.aegypti["Topt", "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

# Get predicted values
a.aegypti.fitted <- function (x){
  m.a_mean * (x - m.means_Topt)^2 + m.c_mean
}

# Data frame of temperature, fitted value, species
m.fitted_data.a.aegypti <- data.frame(Temperature = rep(m.temp_seq),
                                           Fitted = a.aegypti.fitted(m.temp_seq))
# Original data with temperature, trait, and species as numeric
m.plot_data.a.aegypti <- data.frame(
  Temperature = a.aegypti$Temp,
  Trait = a.aegypti$Trait
)

# Plot observed data and fitted curves 
m.p.a.aegypti <- ggplot() +
  geom_point(data = m.plot_data.a.aegypti, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.a.aegypti, aes(x = Temperature, y = Fitted), color = "seagreen4", size = 1) +
  labs(title = "Observed Data and Quadratic Fit",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  theme(legend.position = "none")
m.p.a.aegypti

# Merge the observed and predicted data for comparison
m.merged_data.a.aegypti <- m.fitted_data.a.aegypti %>%
  left_join(a.aegypti, by = c("Temperature" = "Temp"))

# Calculate RMSE for each species
m.rmse.a.aegypti <- m.merged_data.a.aegypti %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE),
    SE = (Trait - Fitted)^2
  )

# Extract log-likelihood matrix (samples × observations)
m.log.lik.a.aegypti <- a.aegypti.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
m.waic.a.aegypti <- waic(m.log.lik.a.aegypti)$estimates["waic", "Estimate"]


save(m.rmse.a.aegypti,
     m.fitted_data.a.aegypti,
     m.merged_data.a.aegypti,
     a.aegypti.fit.mcmc,
     a.aegypti.fit,
     m.waic.a.aegypti,
     m.p.a.aegypti,
     file="m.a.aegypti.fit.RData")

######################################### By Species Individual- Albopictus #########################################
a.albopictus <- test[test$interactor1 == "Aedes albopictus", ]

m.N.obs.ind <- length(a.albopictus$Trait)

sink("aedes.individual.txt")
cat("
model {
  # Priors
   a ~ dexp(0.1)             # Width
  Topt ~ dnorm(30, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  c ~ dexp(0.001)
  
  # Likelihood
  for (i in 1:m.N.obs.ind) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], tau1) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], tau1)
  }}
", file = "aedes.individual.txt")


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

jag.data<-list(trait = a.albopictus$Trait, m.N.obs.ind = m.N.obs.ind, temp = a.albopictus$Temp)

a.albopictus.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                           model.file="aedes.individual.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                           n.iter=ni, DIC=T, working.directory=getwd())
a.albopictus.fit.mcmc <- as.mcmc(a.albopictus.fit) ## makes an "mcmc" object
closeAllConnections()
a.albopictus.fit$BUGSoutput$DIC

# Get rid of log liks
# albopictus.mcmc <- lapply(as.mcmc.list(a.albopictus.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(albopictus.mcmc), autoburnin = FALSE)


# Get posterior means
m.summary_fit.a.albopictus <- summary(a.albopictus.fit.mcmc)$statistics
m.a_mean <- m.summary_fit.a.albopictus["a", "Mean"]
m.c_mean <- m.summary_fit.a.albopictus["c", "Mean"]
m.means_Topt <- m.summary_fit.a.albopictus["Topt", "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

# Get predicted values
a.albopictus.fitted <- function (x){
  m.a_mean * (x - m.means_Topt)^2 + m.c_mean
}

# Data frame of temperature, fitted value, species
m.fitted_data.a.albopictus <- data.frame(Temperature = rep(m.temp_seq),
                                           Fitted = a.albopictus.fitted(m.temp_seq))
# Original data with temperature, trait, and species as numeric
m.plot_data.a.albopictus <- data.frame(
  Temperature = a.albopictus$Temp,
  Trait = a.albopictus$Trait
)

# Plot observed data and fitted curves 
m.p.a.albopictus <- ggplot() +
  geom_point(data = m.plot_data.a.albopictus, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.a.albopictus, aes(x = Temperature, y = Fitted), color = "seagreen4", size = 1) +
  labs(title = "Observed Data and Quadratic Fit",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  theme(legend.position = "none")
m.p.a.albopictus

# Merge the observed and predicted data for comparison
m.merged_data.a.albopictus <- m.fitted_data.a.albopictus %>%
  left_join(a.albopictus, by = c("Temperature" = "Temp"))

# Calculate RMSE for each species
m.rmse.a.albopictus <- m.merged_data.a.albopictus %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE),
    SE = (Trait - Fitted)^2
  )

# Extract log-likelihood matrix (samples × observations)
m.log.lik.a.albopictus <- a.albopictus.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
m.waic.a.albopictus <- waic(m.log.lik.a.albopictus)$estimates["waic", "Estimate"]

save(m.rmse.a.albopictus,
     m.fitted_data.a.albopictus,
     m.merged_data.a.albopictus,
     a.albopictus.fit.mcmc,
     a.albopictus.fit,
     m.waic.a.albopictus,
     m.p.a.albopictus,
     file="m.a.albopictus.fit.RData")

######################################### By Species Individual- Atropalpus #########################################
a.atropalpus <- test[test$interactor1 == "Aedes atropalpus", ]


m.N.obs.ind <- length(a.atropalpus$Trait)

sink("aedes.individual.txt")
cat("
model {
  # Priors
   a ~ dexp(0.1)             # Width
  Topt ~ dnorm(30, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  c ~ dexp(0.001)
  
  # Likelihood
  for (i in 1:m.N.obs.ind) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], tau1) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], tau1)
  }}
", file = "aedes.individual.txt")


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

jag.data<-list(trait = a.atropalpus$Trait, m.N.obs.ind = m.N.obs.ind, temp = a.atropalpus$Temp)

a.atropalpus.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                           model.file="aedes.individual.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                           n.iter=ni, DIC=T, working.directory=getwd())
a.atropalpus.fit.mcmc <- as.mcmc(a.atropalpus.fit) ## makes an "mcmc" object
closeAllConnections()
a.atropalpus.fit$BUGSoutput$DIC

# Get rid of log liks
# atropalpus.mcmc <- lapply(as.mcmc.list(a.atropalpus.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(atropalpus.mcmc), autoburnin = FALSE)


# Get posterior means
m.summary_fit.a.atropalpus <- summary(a.atropalpus.fit.mcmc)$statistics
m.a_mean <- m.summary_fit.a.atropalpus["a", "Mean"]
m.c_mean <- m.summary_fit.a.atropalpus["c", "Mean"]
m.means_Topt <- m.summary_fit.a.atropalpus["Topt", "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

# Get predicted values
a.atropalpus.fitted <- function (x){
  m.a_mean * (x - m.means_Topt)^2 + m.c_mean
}

# Data frame of temperature, fitted value, species
m.fitted_data.a.atropalpus <- data.frame(Temperature = rep(m.temp_seq),
                                           Fitted = a.atropalpus.fitted(m.temp_seq))
# Original data with temperature, trait, and species as numeric
m.plot_data.a.atropalpus <- data.frame(
  Temperature = a.atropalpus$Temp,
  Trait = a.atropalpus$Trait
)

# Plot observed data and fitted curves 
m.p.a.atropalpus <- ggplot() +
  geom_point(data = m.plot_data.a.atropalpus, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.a.atropalpus, aes(x = Temperature, y = Fitted), color = "seagreen4", size = 1) +
  labs(title = "Observed Data and Quadratic Fit",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  theme(legend.position = "none")
m.p.a.atropalpus

# Merge the observed and predicted data for comparison
m.merged_data.a.atropalpus <- m.fitted_data.a.atropalpus %>%
  left_join(a.atropalpus, by = c("Temperature" = "Temp"))

# Calculate RMSE for each species
m.rmse.a.atropalpus <- m.merged_data.a.atropalpus %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE),
    SE = (Trait - Fitted)^2
  )

# Extract log-likelihood matrix (samples × observations)
m.log.lik.a.atropalpus <- a.atropalpus.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
m.waic.a.atropalpus <- waic(m.log.lik.a.atropalpus)$estimates["waic", "Estimate"]


save(m.rmse.a.atropalpus,
     m.fitted_data.a.atropalpus,
     m.merged_data.a.atropalpus,
     a.atropalpus.fit.mcmc,
     a.atropalpus.fit,
     m.waic.a.atropalpus,
     m.p.a.atropalpus,
     file="m.a.atropalpus.fit.RData")

######################################### By Species Individual- Camptorhynchus #########################################
a.camptorhynchus <- test[test$interactor1 == "Aedes camptorhynchus", ]


m.N.obs.ind <- length(a.camptorhynchus$Trait)

sink("aedes.individual.txt")
cat("
model {
  # Priors
   a ~ dexp(0.1)             # Width
  Topt ~ dnorm(30, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  c ~ dexp(0.001)
  
  # Likelihood
  for (i in 1:m.N.obs.ind) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], tau1) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], tau1)
  }}
", file = "aedes.individual.txt")


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

jag.data<-list(trait = a.camptorhynchus$Trait, m.N.obs.ind = m.N.obs.ind, temp = a.camptorhynchus$Temp)

a.camptorhynchus.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                           model.file="aedes.individual.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                           n.iter=ni, DIC=T, working.directory=getwd())
a.camptorhynchus.fit.mcmc <- as.mcmc(a.camptorhynchus.fit) ## makes an "mcmc" object
closeAllConnections()
a.camptorhynchus.fit$BUGSoutput$DIC

# Get rid of log liks
# camptorhynchus.mcmc <- lapply(as.mcmc.list(a.camptorhynchus.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(camptorhynchus.mcmc), autoburnin = FALSE)


# Get posterior means
m.summary_fit.a.camptorhynchus <- summary(a.camptorhynchus.fit.mcmc)$statistics
m.a_mean <- m.summary_fit.a.camptorhynchus["a", "Mean"]
m.c_mean <- m.summary_fit.a.camptorhynchus["c", "Mean"]
m.means_Topt <- m.summary_fit.a.camptorhynchus["Topt", "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

# Get predicted values
a.camptorhynchus.fitted <- function (x){
  m.a_mean * (x - m.means_Topt)^2 + m.c_mean
}

# Data frame of temperature, fitted value, species
m.fitted_data.a.camptorhynchus <- data.frame(Temperature = rep(m.temp_seq),
                                           Fitted = a.camptorhynchus.fitted(m.temp_seq))
# Original data with temperature, trait, and species as numeric
m.plot_data.a.camptorhynchus <- data.frame(
  Temperature = a.camptorhynchus$Temp,
  Trait = a.camptorhynchus$Trait
)

# Plot observed data and fitted curves 
m.p.a.camptorhynchus <- ggplot() +
  geom_point(data = m.plot_data.a.camptorhynchus, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.a.camptorhynchus, aes(x = Temperature, y = Fitted), color = "seagreen4", size = 1) +
  labs(title = "Observed Data and Quadratic Fit",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  theme(legend.position = "none")
m.p.a.camptorhynchus

# Merge the observed and predicted data for comparison
m.merged_data.a.camptorhynchus <- m.fitted_data.a.camptorhynchus %>%
  left_join(a.camptorhynchus, by = c("Temperature" = "Temp"))

# Calculate RMSE for each species
m.rmse.a.camptorhynchus <- m.merged_data.a.camptorhynchus %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE),
    SE = (Trait - Fitted)^2
  )

# Extract log-likelihood matrix (samples × observations)
m.log.lik.a.camptorhynchus <- a.camptorhynchus.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
m.waic.a.camptorhynchus <- waic(m.log.lik.a.camptorhynchus)$estimates["waic", "Estimate"]


save(m.rmse.a.camptorhynchus,
     m.fitted_data.a.camptorhynchus,
     m.merged_data.a.camptorhynchus,
     a.camptorhynchus.fit.mcmc,
     a.camptorhynchus.fit,
     m.waic.a.camptorhynchus,
     m.p.a.camptorhynchus,
     file="m.a.camptorhynchus.fit.RData")

######################################### By Species Individual- Japonicus Japonicus #########################################
a.jj <- test[test$interactor1 == "Aedes japonicus japonicus", ]


m.N.obs.ind <- length(a.jj$Trait)

sink("aedes.individual.txt")
cat("
model {
  # Priors
   a ~ dexp(0.1)             # Width
  Topt ~ dnorm(30, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  c ~ dexp(0.001)
  
  # Likelihood
  for (i in 1:m.N.obs.ind) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], tau1) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], tau1)
  }}
", file = "aedes.individual.txt")


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

jag.data<-list(trait = a.jj$Trait, m.N.obs.ind = m.N.obs.ind, temp = a.jj$Temp)

a.jj.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                           model.file="aedes.individual.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                           n.iter=ni, DIC=T, working.directory=getwd())
a.jj.fit.mcmc <- as.mcmc(a.jj.fit) ## makes an "mcmc" object
closeAllConnections()
a.jj.fit$BUGSoutput$DIC

# Get rid of log liks
# jj.mcmc <- lapply(as.mcmc.list(a.jj.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(jj.mcmc), autoburnin = FALSE)


# Get posterior means
m.summary_fit.a.jj <- summary(a.jj.fit.mcmc)$statistics
m.a_mean <- m.summary_fit.a.jj["a", "Mean"]
m.c_mean <- m.summary_fit.a.jj["c", "Mean"]
m.means_Topt <- m.summary_fit.a.jj["Topt", "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

# Get predicted values
a.jj.fitted <- function (x){
  m.a_mean * (x - m.means_Topt)^2 + m.c_mean
}

# Data frame of temperature, fitted value, species
m.fitted_data.a.jj <- data.frame(Temperature = rep(m.temp_seq),
                                           Fitted = a.jj.fitted(m.temp_seq))
# Original data with temperature, trait, and species as numeric
m.plot_data.a.jj <- data.frame(
  Temperature = a.jj$Temp,
  Trait = a.jj$Trait
)

# Plot observed data and fitted curves 
m.p.a.jj <- ggplot() +
  geom_point(data = m.plot_data.a.jj, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.a.jj, aes(x = Temperature, y = Fitted), color = "seagreen4", size = 1) +
  labs(title = "Observed Data and Quadratic Fit",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  theme(legend.position = "none")
m.p.a.jj

# Merge the observed and predicted data for comparison
m.merged_data.a.jj <- m.fitted_data.a.jj %>%
  left_join(a.jj, by = c("Temperature" = "Temp"))

# Calculate RMSE for each species
m.rmse.a.jj <- m.merged_data.a.jj %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE),
    SE = (Trait - Fitted)^2
  )

# Extract log-likelihood matrix (samples × observations)
m.log.lik.a.jj <- a.jj.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
m.waic.a.jj <- waic(m.log.lik.a.jj)$estimates["waic", "Estimate"]




save(m.rmse.a.jj,
     m.fitted_data.a.jj,
     m.merged_data.a.jj,
     a.jj.fit.mcmc,
     a.jj.fit,
     m.waic.a.jj,
     m.p.a.jj,
     file="m.a.jj.fit.RData")

######################################### By Species Individual- Krombeini #########################################
a.krombeini <- test[test$interactor1 == "Aedes krombeini", ]

m.N.obs.ind <- length(a.krombeini$Trait)

sink("aedes.individual.txt")
cat("
model {
  # Priors
   a ~ dexp(0.1)             # Width
  Topt ~ dnorm(30, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  c ~ dexp(0.001)
  
  # Likelihood
  for (i in 1:m.N.obs.ind) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], tau1) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], tau1)
  }}
", file = "aedes.individual.txt")


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

jag.data<-list(trait = a.krombeini$Trait, m.N.obs.ind = m.N.obs.ind, temp = a.krombeini$Temp)

a.krombeini.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                           model.file="aedes.individual.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                           n.iter=ni, DIC=T, working.directory=getwd())
a.krombeini.fit.mcmc <- as.mcmc(a.krombeini.fit) ## makes an "mcmc" object
closeAllConnections()
a.krombeini.fit$BUGSoutput$DIC

# Get rid of log liks
# krombeini.mcmc <- lapply(as.mcmc.list(a.krombeini.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(krombeini.mcmc), autoburnin = FALSE)


# Get posterior means
m.summary_fit.a.krombeini <- summary(a.krombeini.fit.mcmc)$statistics
m.a_mean <- m.summary_fit.a.krombeini["a", "Mean"]
m.c_mean <- m.summary_fit.a.krombeini["c", "Mean"]
m.means_Topt <- m.summary_fit.a.krombeini["Topt", "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

# Get predicted values
a.krombeini.fitted <- function (x){
  m.a_mean * (x - m.means_Topt)^2 + m.c_mean
}

# Data frame of temperature, fitted value, species
m.fitted_data.a.krombeini <- data.frame(Temperature = rep(m.temp_seq),
                                           Fitted = a.krombeini.fitted(m.temp_seq))
# Original data with temperature, trait, and species as numeric
m.plot_data.a.krombeini <- data.frame(
  Temperature = a.krombeini$Temp,
  Trait = a.krombeini$Trait
)

# Plot observed data and fitted curves 
m.p.a.krombeini <- ggplot() +
  geom_point(data = m.plot_data.a.krombeini, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.a.krombeini, aes(x = Temperature, y = Fitted), color = "seagreen4", size = 1) +
  labs(title = "Observed Data and Quadratic Fit",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  theme(legend.position = "none")
m.p.a.krombeini

# Merge the observed and predicted data for comparison
m.merged_data.a.krombeini <- m.fitted_data.a.krombeini %>%
  left_join(a.krombeini, by = c("Temperature" = "Temp"))

# Calculate RMSE for each species
m.rmse.a.krombeini <- m.merged_data.a.krombeini %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE),
    SE = (Trait - Fitted)^2
  )

# Extract log-likelihood matrix (samples × observations)
m.log.lik.a.krombeini <- a.krombeini.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
m.waic.a.krombeini <- waic(m.log.lik.a.krombeini)$estimates["waic", "Estimate"]


save(m.rmse.a.krombeini,
     m.fitted_data.a.krombeini,
     m.merged_data.a.krombeini,
     a.krombeini.fit.mcmc,
     a.krombeini.fit,
     m.waic.a.krombeini,
     m.p.a.krombeini,
     file="m.a.krombeini.fit.RData")

######################################### By Species Individual- Notoscriptus #########################################
a.notoscriptus <- test[test$interactor1 == "Aedes notoscriptus", ]

m.N.obs.ind <- length(a.notoscriptus$Trait)

sink("aedes.individual.txt")
cat("
model {
  # Priors
   a ~ dexp(0.1)             # Width
  Topt ~ dnorm(30, 0.001)             # Topt
  tau1 ~ dexp(0.01)           # Tau for trait
  c ~ dexp(0.001)
  
  # Likelihood
  for (i in 1:m.N.obs.ind) {
    mu[i] <- a * (Topt - temp[i])^2 + c
    trait[i] ~ dnorm(mu[i], tau1) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], tau1)
  }}
", file = "aedes.individual.txt")


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

jag.data<-list(trait = a.notoscriptus$Trait, m.N.obs.ind = m.N.obs.ind, temp = a.notoscriptus$Temp)

a.notoscriptus.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                           model.file="aedes.individual.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                           n.iter=ni, DIC=T, working.directory=getwd())
a.notoscriptus.fit.mcmc <- as.mcmc(a.notoscriptus.fit) ## makes an "mcmc" object
closeAllConnections()
a.notoscriptus.fit$BUGSoutput$DIC

# Get rid of log liks
# notoscriptus.mcmc <- lapply(as.mcmc.list(a.notoscriptus.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(notoscriptus.mcmc), autoburnin = FALSE)


# Get posterior means
m.summary_fit.a.notoscriptus <- summary(a.notoscriptus.fit.mcmc)$statistics
m.a_mean <- m.summary_fit.a.notoscriptus["a", "Mean"]
m.c_mean <- m.summary_fit.a.notoscriptus["c", "Mean"]
m.means_Topt <- m.summary_fit.a.notoscriptus["Topt", "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

# Get predicted values
a.notoscriptus.fitted <- function (x){
  m.a_mean * (x - m.means_Topt)^2 + m.c_mean
}

# Data frame of temperature, fitted value, species
m.fitted_data.a.notoscriptus <- data.frame(Temperature = rep(m.temp_seq),
                                           Fitted = a.notoscriptus.fitted(m.temp_seq))
# Original data with temperature, trait, and species as numeric
m.plot_data.a.notoscriptus <- data.frame(
  Temperature = a.notoscriptus$Temp,
  Trait = a.notoscriptus$Trait
)

# Plot observed data and fitted curves 
m.p.a.notoscriptus <- ggplot() +
  geom_point(data = m.plot_data.a.notoscriptus, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.a.notoscriptus, aes(x = Temperature, y = Fitted), color = "seagreen4", size = 1) +
  labs(title = "Observed Data and Quadratic Fit",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  theme(legend.position = "none")
m.p.a.notoscriptus

# Merge the observed and predicted data for comparison
m.merged_data.a.notoscriptus <- m.fitted_data.a.notoscriptus %>%
  left_join(a.notoscriptus, by = c("Temperature" = "Temp"))

# Calculate RMSE for each species
m.rmse.a.notoscriptus <- m.merged_data.a.notoscriptus %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE),
    SE = (Trait - Fitted)^2
  )

# Extract log-likelihood matrix (samples × observations)
m.log.lik.a.notoscriptus <- a.notoscriptus.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
m.waic.a.notoscriptus <- waic(m.log.lik.a.notoscriptus)$estimates["waic", "Estimate"]



save(m.rmse.a.notoscriptus,
     m.fitted_data.a.notoscriptus,
     m.merged_data.a.notoscriptus,
     a.notoscriptus.fit.mcmc,
     a.notoscriptus.fit,
     m.waic.a.notoscriptus,
     m.p.a.notoscriptus,
     file="m.a.notoscriptus.fit.RData")

######################################### By Species Individual- Summary #########################################
# Define RMSE values for each species
m.se.individual <- c(m.rmse.a.aegypti$SE, m.rmse.a.albopictus$SE, m.rmse.a.atropalpus$SE,
                   m.rmse.a.camptorhynchus$SE, m.rmse.a.jj$SE,
                   m.rmse.a.krombeini$SE, m.rmse.a.notoscriptus$SE)

# Compute weighted average RMSE
m.overall_rmse.species.ind <- sqrt(mean(m.se.individual, na.rm = TRUE))

m.DIC.species.individual <- a.aegypti.fit$BUGSoutput$DIC+ a.albopictus.fit$BUGSoutput$DIC +
  a.atropalpus.fit$BUGSoutput$DIC + a.camptorhynchus.fit$BUGSoutput$DIC + a.jj.fit$BUGSoutput$DIC +
  a.krombeini.fit$BUGSoutput$DIC + a.notoscriptus.fit$BUGSoutput$DIC

m.waic.species.ind <- m.waic.a.aegypti + m.waic.a.albopictus +
  m.waic.a.atropalpus + m.waic.a.camptorhynchus +
  m.waic.a.jj + m.waic.a.krombeini + m.waic.a.notoscriptus

# Arrange plots into a 3x3 grid
m.grid.ind.species <- grid.arrange(
  m.p.a.aegypti + labs(title = "Aedes aegypti", x = " ", y = "Development Time in Days"),
  m.p.a.albopictus + labs(title = "Aedes albopictus", x = " ", y = " "),
  m.p.a.atropalpus + labs(title = "Aedes atropalpus", x = " ", y = " "),
  m.p.a.camptorhynchus + labs(title = "Aedes camptorhynchus", x = " ", y = "Development Time in Days"),
  m.p.a.jj + labs(title = "Aedes japonicus japonicus", x = " ", y = " "),
  m.p.a.krombeini + labs(title = "Aedes krombeini", x = " ", y = " "),
  m.p.a.notoscriptus + labs(title = "Aedes notoscriptus", y = "Development Time in Days"),
  nrow = 3, ncol = 3
)
ggsave("mosquito.species.ind.plot.png", m.grid.ind.species, width = 8.5, height = 6)


######################################### By Species Hierarchical #########################################
sink("aedes.by.species.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Species-specific 
  mu.c ~ dexp(10)             # Mean for species effects c
  mu.Topt ~ dnorm(30, 0.1)             # Mean for species effects Topt
  tau.Topt ~ dexp(0.01)             # Tau for species effect Topt
  mu.a ~ dexp(0.01)            # Mean for species effects Topt

  for (j in 1:m.n.species) {
    mu.species[j] ~ dexp(mu.c)       # Species-specific low point
    mu.species.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Species-specific temp
    mu.species.a[j] ~ dexp(mu.a)    # Species-specific temp
  }
  # Likelihood
  for (i in 1:m.N.obs) {
    mu[i] <- mu.species.a[m.species[i]] * (mu.species.Topt[m.species[i]] - temp[i])^2 + mu.species[m.species[i]]
    trait[i] ~ dnorm(mu[i], tau1) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], tau1)
  }}
", file = "aedes.by.species.txt")


# Parameters to monitor
parameters <- c("mu.a", "mu.c", "mu.species", "log_lik",
                "tau1", "mu.Topt", "tau.Topt", "mu.species.Topt", "mu.species.a")


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

jag.data <- list(trait = test$Trait, m.N.obs = m.N.obs, temp = test$Temp, 
                 m.n.species = m.n.species, m.species = m.species)

m.species.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
               model.file="aedes.by.species.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
               n.iter=ni, DIC=T, working.directory=getwd())
m.species.fit.mcmc <- as.mcmc(m.species.fit) ## makes an "mcmc" object
closeAllConnections()
m.species.fit$BUGSoutput$DIC

# Get rid of log liks
# m.species.mcmc <- lapply(as.mcmc.list(m.species.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(m.species.mcmc), autoburnin = FALSE)

# Extract log-likelihood matrix (samples × observations)
m.log.lik.species <- m.species.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
m.waic.species <- waic(m.log.lik.species)$estimates["waic", "Estimate"]

# Get posterior means
m.summary_fit.species <- summary(m.species.fit.mcmc)$statistics
m.a_mean_species <- m.summary_fit.species[grep("^mu.species.a\\[", rownames(m.summary_fit.species)), "Mean"]
m.c_mean_species <- m.summary_fit.species[grep("^mu.species\\[", rownames(m.summary_fit.species)), "Mean"]
m.Topt_mean_species <- m.summary_fit.species[grep("^mu.species.Topt\\[", rownames(m.summary_fit.species)), "Mean"]


# Generate temperature sequence 
m.temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

# Get predicted values
m.spec1 <- function (x){
  m.a_mean_species[1] * (x - m.Topt_mean_species[1])^2 + m.c_mean_species[1]
}
m.spec2 <- function (x){
  m.a_mean_species[2] * (x - m.Topt_mean_species[2])^2 + m.c_mean_species[2]
}
m.spec3 <- function (x){
  m.a_mean_species[3] * (x - m.Topt_mean_species[3])^2 + m.c_mean_species[3]
}
m.spec4 <- function (x){
  m.a_mean_species[4] * (x - m.Topt_mean_species[4])^2 + m.c_mean_species[4]
}
m.spec5 <- function (x){
  m.a_mean_species[5] * (x - m.Topt_mean_species[5])^2 + m.c_mean_species[5]
}
m.spec6 <- function (x){
  m.a_mean_species[6] * (x - m.Topt_mean_species[6])^2 + m.c_mean_species[6]
}
m.spec7 <- function (x){
  m.a_mean_species[7] * (x - m.Topt_mean_species[7])^2 + m.c_mean_species[7]
}

# Optionally, convert to a data frame for easier plotting
m.fitted_data.species <- data.frame(Temperature = rep(m.temp_seq, m.n.species), 
                          Fitted = c(m.spec1(m.temp_seq),
                                     m.spec2(m.temp_seq),
                                     m.spec3(m.temp_seq),
                                     m.spec4(m.temp_seq),
                                     m.spec5(m.temp_seq),
                                     m.spec6(m.temp_seq),
                                     m.spec7(m.temp_seq)),
                          Species = factor(rep(c("Aedes aegypti", "Aedes albopictus",
                                                 "Aedes atropalpus", "Aedes camptorhynchus",
                                                 "Aedes japonicus japonicus",
                                                 "Aedes krombeini", "Aedes notoscriptus"), each = length(m.temp_seq))))

m.plot_data.species <- data.frame(
  Temperature = test$Temp,
  Trait = test$Trait,
  Species = as.factor(test$interactor1)
)

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = m.plot_data.species, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.species, aes(x = Temperature, y = Fitted, color = Species), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
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


# Merge the observed and predicted data for comparison
m.merged_data.species <- m.fitted_data.species %>%
  left_join(m.plot_data.species, by = c("Temperature", "Species"))


# Calculate RMSE for each species
m.mse_by_species.species <- m.merged_data.species %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.species <- m.merged_data.species %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 



save(m.mse_by_species.species,
     m.fitted_data.species,
     m.merged_data.species,
     m.species.fit.mcmc,
     m.species.fit,
     m.log.lik.species,
     file="m.species.fit.RData")

######################################### By Species Hierarchical Different Variances #########################################
sink("aedes.by.species.tau.txt")
cat("
model {
  # Priors
  # Species-specific 
  mu.c ~ dexp(10)             # Mean for species effects c
  tau1 ~ dexp(0.01)           # Tau for species effect  
  mu.Topt ~ dnorm(30, 0.1)             # Mean for species effects Topt
  tau.Topt ~ dexp(0.01)             # Tau for species effect Topt
  mu.a ~ dexp(0.01)            # Mean for species effects Topt

  for (j in 1:m.n.species) {
    mu.species[j] ~ dexp(mu.c)       # Species-specific low point
    mu.species.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Species-specific temp
    mu.species.a[j] ~ dexp(mu.a)    # Species-specific temp
    species.tau[j] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:m.N.obs) {
    mu[i] <- mu.species.a[m.species[i]] * (mu.species.Topt[m.species[i]] - temp[i])^2 + mu.species[m.species[i]]
    trait[i] ~ dnorm(mu[i], species.tau[m.species[i]]) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], species.tau[m.species[i]])
  }}
", file = "aedes.by.species.tau.txt")

# Parameters to monitor
parameters <- c("mu.a", "mu.c", "mu.species", "species.tau", "log_lik",
                "tau1", "mu.Topt", "tau.Topt", "mu.species.Topt", "mu.species.a")


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

jag.data <- list(trait = test$Trait, m.N.obs = m.N.obs, temp = test$Temp, 
                 m.n.species = m.n.species, m.species = m.species)

m.species.tau.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                      model.file="aedes.by.species.tau.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                      n.iter=ni, DIC=T, working.directory=getwd())
m.species.tau.fit.mcmc <- as.mcmc(m.species.tau.fit) ## makes an "mcmc" object
closeAllConnections()
m.species.tau.fit$BUGSoutput$DIC

# Get rid of log liks
# m.species.tau.mcmc <- lapply(as.mcmc.list(m.species.tau.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(m.species.tau.mcmc), autoburnin = FALSE)

# Extract log-likelihood matrix (samples × observations)
m.log.lik.species.tau <- m.species.tau.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
m.waic.species.tau <- waic(m.log.lik.species.tau)$estimates["waic", "Estimate"]

# Get posterior means
m.summary_fit.species.tau <- summary(m.species.tau.fit.mcmc)$statistics
m.a_mean_species.tau <- m.summary_fit.species.tau[grep("^mu.species.a\\[", rownames(m.summary_fit.species.tau)), "Mean"]
m.c_mean_species.tau <- m.summary_fit.species.tau[grep("^mu.species\\[", rownames(m.summary_fit.species.tau)), "Mean"]
m.Topt_mean_species.tau <- m.summary_fit.species.tau[grep("^mu.species.Topt\\[", rownames(m.summary_fit.species.tau)), "Mean"]


# Generate temperature sequence 
m.temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

# Get predicted values
m.spec1.tau <- function (x){
  m.a_mean_species.tau[1] * (x - m.Topt_mean_species.tau[1])^2 + m.c_mean_species.tau[1]
}
m.spec2.tau <- function (x){
  m.a_mean_species.tau[2] * (x - m.Topt_mean_species.tau[2])^2 + m.c_mean_species.tau[2]
}
m.spec3.tau <- function (x){
  m.a_mean_species.tau[3] * (x - m.Topt_mean_species.tau[3])^2 + m.c_mean_species.tau[3]
}
m.spec4.tau <- function (x){
  m.a_mean_species.tau[4] * (x - m.Topt_mean_species.tau[4])^2 + m.c_mean_species.tau[4]
}
m.spec5.tau <- function (x){
  m.a_mean_species.tau[5] * (x - m.Topt_mean_species.tau[5])^2 + m.c_mean_species.tau[5]
}
m.spec6.tau <- function (x){
  m.a_mean_species.tau[6] * (x - m.Topt_mean_species.tau[6])^2 + m.c_mean_species.tau[6]
}
m.spec7.tau <- function (x){
  m.a_mean_species.tau[7] * (x - m.Topt_mean_species.tau[7])^2 + m.c_mean_species.tau[7]
}

# Optionally, convert to a data frame for easier plotting
m.fitted_data.species.tau <- data.frame(Temperature = rep(m.temp_seq, m.n.species), 
                                    Fitted = c(m.spec1.tau(m.temp_seq),
                                               m.spec2.tau(m.temp_seq),
                                               m.spec3.tau(m.temp_seq),
                                               m.spec4.tau(m.temp_seq),
                                               m.spec5.tau(m.temp_seq),
                                               m.spec6.tau(m.temp_seq),
                                               m.spec7.tau(m.temp_seq)),
                                    Species = factor(rep(c("Aedes aegypti", "Aedes albopictus",
                                                           "Aedes atropalpus", "Aedes camptorhynchus",
                                                           "Aedes japonicus japonicus",
                                                           "Aedes krombeini", "Aedes notoscriptus"), each = length(m.temp_seq))))

m.plot_data.species <- data.frame(
  Temperature = test$Temp,
  Trait = test$Trait,
  Species = as.factor(test$interactor1)
)

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = m.plot_data.species, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.species.tau, aes(x = Temperature, y = Fitted, color = Species), size = 1) +
  labs(x = "Temperature", y = "Development Time in Days") +
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
ggsave("mosquito.species.tau.plot.png", width = 8, height = 5)

# Merge the observed and predicted data for comparison
m.merged_data.species.tau <- m.fitted_data.species.tau %>%
  left_join(m.plot_data.species, by = c("Temperature", "Species"))


# Calculate RMSE for each species
m.mse_by_species.species.tau <- m.merged_data.species.tau %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.species.tau <- m.merged_data.species.tau %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 

save(m.mse_by_species.species.tau,
     m.fitted_data.species.tau,
     m.merged_data.species.tau,
     m.species.tau.fit.mcmc,
     m.species.tau.fit,
     m.log.lik.species.tau,
     file="m.species.tau.fit.RData")

######################################### By Climate #########################################
m.climate <- as.factor(test$climate)
m.n.climate <- length(unique(test$climate))
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

  for (j in 1:m.n.climate) {
    mu.climate[j] ~ dexp(mu.c)       # Climate-specific low point
    mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific temp
    mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific temp
    climate.tau[j] ~ dexp(tau1)
  }

  # Likelihood
  for (i in 1:m.N.obs) {
    mu[i] <- mu.climate.a[m.climate[i]] * (mu.climate.Topt[m.climate[i]] - temp[i])^2 + mu.climate[m.climate[i]]
    trait[i] ~ dnorm(mu[i], climate.tau[m.climate[i]]) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], climate.tau[m.climate[i]])
  }}
", file = "aedes.by.climate.txt")


# Parameters to monitor
parameters <- c("mu.a", "mu.c", "mu.climate", "climate.tau", "log_lik",
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

jag.data<-list(trait = test$Trait, m.N.obs = m.N.obs, temp = test$Temp, 
               m.n.climate = m.n.climate, m.climate = m.climate)

m.climate.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
               model.file="aedes.by.climate.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
               n.iter=ni, DIC=T, working.directory=getwd())
m.climate.fit.mcmc <- as.mcmc(m.climate.fit) ## makes an "mcmc" object
closeAllConnections()
m.climate.fit$BUGSoutput$DIC

# Get rid of log liks
# m.climate.mcmc <- lapply(as.mcmc.list(m.climate.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(m.climate.mcmc), autoburnin = FALSE)

# Extract log-likelihood matrix (samples × observations)
m.log.lik.climate <- m.climate.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
m.waic.climate<- waic(m.log.lik.climate)$estimates["waic", "Estimate"]

# Get posterior means
m.summary_fit.climate <- summary(m.climate.fit.mcmc)$statistics
m.a_mean_climate <- m.summary_fit.climate[grep("^mu.climate.a\\[", rownames(m.summary_fit.climate)), "Mean"]
m.c_mean_climate <- m.summary_fit.climate[grep("^mu.climate\\[", rownames(m.summary_fit.climate)), "Mean"]
m.Topt_mean_climate <- m.summary_fit.climate[grep("^mu.climate.Topt\\[", rownames(m.summary_fit.climate)), "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

# Get predicted values
m.clim1 <- function (x){
  m.a_mean_climate[1] * (x - m.Topt_mean_climate[1])^2 + m.c_mean_climate[1]
}
m.clim2 <- function (x){
  m.a_mean_climate[2] * (x - m.Topt_mean_climate[2])^2 + m.c_mean_climate[2]
}
m.clim3 <- function (x){
  m.a_mean_climate[3] * (x - m.Topt_mean_climate[3])^2 + m.c_mean_climate[3]
}


# Optionally, convert to a data frame for easier plotting
m.fitted_data.climate <- data.frame(Temperature = rep(m.temp_seq, m.n.climate), 
                          Fitted = c(m.clim1(m.temp_seq),
                                     m.clim2(m.temp_seq),
                                     m.clim3(m.temp_seq)),
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
m.fitted_data.climate <- data.frame(
  Temperature = rep(m.temp_seq, m.n.climate),
  Fitted = c(m.clim1(m.temp_seq),
             m.clim2(m.temp_seq),
             m.clim3(m.temp_seq)),
  Climate = factor(rep(c("subtropical", "temperate", "tropical"), each = length(m.temp_seq)))
) %>%
  left_join(m.species_to_climate, by = "Climate")

m.plot_data.climate$Climate <- factor(m.plot_data.climate$Climate, levels = c("temperate", "subtropical", "tropical"))
m.fitted_data.climate$Climate <- factor(m.fitted_data.climate$Climate, levels = c("temperate", "subtropical", "tropical"))
# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = m.plot_data.climate, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.climate, aes(x = Temperature, y = Fitted, color = Climate), size = 1) +
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

# Merge the observed and predicted data for comparison
m.merged_data.climate <- m.fitted_data.climate %>%
  left_join(m.plot_data.climate, by = c("Temperature", "Climate", "Species"))


# Calculate RMSE for each climate
m.mse_by_species.climate <- m.merged_data.climate |> 
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.climate <- m.merged_data.climate |> 
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 



save(m.mse_by_species.climate,
     m.fitted_data.climate,
     m.merged_data.climate,
     m.climate.fit.mcmc,
     m.climate.fit,
     m.log.lik.climate,
     file="m.climate.fit.RData")


######################################### By Wing Size Category #########################################
m.size <- as.factor(test$wing.size)
m.n.size <- length(unique(test$wing.size))

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

  for (j in 1:m.n.size) {
    mu.size[j] ~ dexp(mu.c)       # Size-specific low point
    mu.size.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Size-specific temp
    mu.size.a[j] ~ dexp(mu.a)    # Size-specific temp
    size.tau[j] ~ dexp(tau1)
  }

  # Likelihood
  for (i in 1:m.N.obs) {
    mu[i] <- mu.size.a[m.size[i]] * (mu.size.Topt[m.size[i]] - temp[i])^2 + mu.size[m.size[i]]
    trait[i] ~ dnorm(mu[i], size.tau[m.size[i]]) T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], size.tau[m.size[i]])
  }}
", file = "aedes.by.size.txt")

# Parameters to monitor
parameters <- c("mu.a", "mu.c", "mu.size", "size.tau", "log_lik",
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

jag.data<-list(trait = test$Trait, m.N.obs = m.N.obs, temp = test$Temp, 
               m.n.size = m.n.size, m.size = m.size)

m.size.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
               model.file="aedes.by.size.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
               n.iter=ni, DIC=T, working.directory=getwd())
m.size.fit.mcmc <- as.mcmc(m.size.fit) ## makes an "mcmc" object
closeAllConnections()
m.size.fit$BUGSoutput$DIC

# Get rid of log liks
# m.size.mcmc <- lapply(as.mcmc.list(m.size.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(m.size.mcmc), autoburnin = FALSE)

# Extract log-likelihood matrix (samples × observations)
m.log.lik.size <- m.size.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
m.waic.size <- waic(m.log.lik.size)$estimates["waic", "Estimate"]

# Get posterior means
m.summary_fit.size <- summary(m.size.fit.mcmc)$statistics
m.a_mean_size <- m.summary_fit.size[grep("^mu.size.a\\[", rownames(m.summary_fit.size)), "Mean"]
m.c_mean_size <- m.summary_fit.size[grep("^mu.size\\[", rownames(m.summary_fit.size)), "Mean"]
m.Topt_mean_size <- m.summary_fit.size[grep("^mu.size.Topt\\[", rownames(m.summary_fit.size)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

# Get predicted values
m.size1 <- function (x){
  m.a_mean_size[1] * (x - m.Topt_mean_size[1])^2 + m.c_mean_size[1]
}
m.size2 <- function (x){
  m.a_mean_size[2] * (x - m.Topt_mean_size[2])^2 + m.c_mean_size[2]
}
m.size3 <- function (x){
  m.a_mean_size[3] * (x - m.Topt_mean_size[3])^2 + m.c_mean_size[3]
}


# Optionally, convert to a data frame for easier plotting
m.fitted_data.size <- data.frame(Temperature = rep(m.temp_seq, m.n.size), 
                          Fitted = c(m.size1(m.temp_seq),
                                     m.size2(m.temp_seq),
                                     m.size3(m.temp_seq)),
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
m.fitted_data.size <- data.frame(
  Temperature = rep(m.temp_seq, m.n.size),
  Fitted = c(m.size1(m.temp_seq),
             m.size2(m.temp_seq),
             m.size3(m.temp_seq)),
  Size = factor(rep(c("large", "medium", "small"), each = length(m.temp_seq)))
) %>%
  left_join(m.species_to_size, by = "Size")

m.plot_data.size$Size <- factor(m.plot_data.size$Size, levels = c("small", "medium", "large"))
m.fitted_data.size$Size <- factor(m.fitted_data.size$Size, levels = c("small", "medium", "large"))
# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = m.plot_data.size, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.size, aes(x = Temperature, y = Fitted, color = Size), size = 1) +
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


# Merge the observed and predicted data for comparison
m.merged_data.size <- m.fitted_data.size %>%
  left_join(m.plot_data.size, by = c("Temperature", "Size", "Species"))


# Calculate RMSE for each species
m.mse_by_species.size <- m.merged_data.size %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.size <- m.merged_data.size %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 



save(m.mse_by_species.size,
     m.fitted_data.size,
     m.merged_data.size,
     m.size.fit.mcmc,
     m.size.fit,
     m.log.lik.size,
     file="m.size.fit.RData")

######################################### By Wing Size Category Removed #########################################
test1 <- test[!(test$interactor1 %in% c("Aedes notoscriptus", "Aedes krombeini")), ]
m.size <- as.factor(test1$wing.size)
m.n.size <- length(unique(test1$wing.size))
m.N.obs.rm <- length(test1$Trait)

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

  for (j in 1:m.n.size) {
    mu.size[j] ~ dexp(mu.c)       # Size-specific low point
    mu.size.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Size-specific temp
    mu.size.a[j] ~ dexp(mu.a)    # Size-specific temp
    size.tau[j] ~ dexp(tau1)
  }

  # Likelihood
  for (i in 1:m.N.obs.rm) {
    mu[i] <- mu.size.a[m.size[i]] * (mu.size.Topt[m.size[i]] - temp[i])^2 + mu.size[m.size[i]]
    trait[i] ~ dnorm(mu[i], size.tau[m.size[i]]) T(0, ) 
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
#jag.data<-list(trait = test$Trait, N.obs = N.obs, temp = test$Temp)

jag.data<-list(trait = test1$Trait, m.N.obs.rm = m.N.obs.rm, temp = test1$Temp, 
               m.n.size = m.n.size, m.size = m.size)

m.size.rm.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                   model.file="aedes.by.size.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                   n.iter=ni, DIC=T, working.directory=getwd())
m.size.rm.fit.mcmc <- as.mcmc(m.size.rm.fit) ## makes an "mcmc" object
closeAllConnections()
m.size.rm.fit$BUGSoutput$DIC

# Get rid of log liks
# m.size.rm.mcmc <- lapply(as.mcmc.list(m.size.rm.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(m.size.rm.mcmc), autoburnin = FALSE)

# Get posterior means
m.summary_fit.size.rm <- summary(m.size.rm.fit.mcmc)$statistics
m.a_mean_size.rm <- m.summary_fit.size.rm[grep("^mu.size.a\\[", rownames(m.summary_fit.size.rm)), "Mean"]
m.c_mean_size.rm <- m.summary_fit.size.rm[grep("^mu.size\\[", rownames(m.summary_fit.size.rm)), "Mean"]
m.Topt_mean_size.rm <- m.summary_fit.size.rm[grep("^mu.size.Topt\\[", rownames(m.summary_fit.size.rm)), "Mean"]

# Generate temperature sequence 
m.temp_seq <- seq(min(test$Temp), max(test$Temp), by = 0.1)

# Get predicted values
m.size1.rm <- function (x){
  m.a_mean_size.rm[1] * (x - m.Topt_mean_size.rm[1])^2 + m.c_mean_size.rm[1]
}
m.size2.rm <- function (x){
  m.a_mean_size.rm[2] * (x - m.Topt_mean_size.rm[2])^2 + m.c_mean_size.rm[2]
}
m.size3.rm <- function (x){
  m.a_mean_size.rm[3] * (x - m.Topt_mean_size.rm[3])^2 + m.c_mean_size.rm[3]
}


# Optionally, convert to a data frame for easier plotting
m.fitted_data.size.rm <- data.frame(Temperature = rep(m.temp_seq, m.n.size), 
                                 Fitted = c(m.size1.rm(m.temp_seq),
                                            m.size2.rm(m.temp_seq),
                                            m.size3.rm(m.temp_seq)),
                                 Size = factor(rep(c("large", "medium", "small"), each = length(m.temp_seq))))

m.plot_data.size.rm <- data.frame(
  Temperature = test$Temp,
  Trait = test$Trait,
  Size = as.factor(test$wing.size),
  Species = as.factor(test$interactor1)
)

m.species_to_size.rm <- unique(m.plot_data.size[, c("Species", "Size")])  # Ensure species matches correct size
m.species_to_size.rm$Size <- as.factor(m.species_to_size.rm$Size)
# Merge fitted size curves with species-to-size mapping
m.fitted_data.size.rm <- data.frame(
  Temperature = rep(m.temp_seq, m.n.size),
  Fitted = c(m.size1.rm(m.temp_seq),
             m.size2.rm(m.temp_seq),
             m.size3.rm(m.temp_seq)),
  Size = factor(rep(c("large", "medium", "small"), each = length(m.temp_seq)))
) %>%
  left_join(m.species_to_size.rm, by = "Size")

m.plot_data.size.rm$Size <- factor(m.plot_data.size.rm$Size, levels = c("small", "medium", "large"))
m.fitted_data.size.rm$Size <- factor(m.fitted_data.size.rm$Size, levels = c("small", "medium", "large"))
# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = m.plot_data.size.rm, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.size.rm, aes(x = Temperature, y = Fitted, color = Size), size = 1) +
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


# Merge the observed and predicted data for comparison
m.merged_data.size.rm <- m.fitted_data.size.rm %>%
  left_join(m.plot_data.size.rm, by = c("Temperature", "Size", "Species"))


# Calculate RMSE for each species
m.mse_by_species.size.rm <- m.merged_data.size.rm %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.size.rm <- m.merged_data.size.rm %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 



save(m.mse_by_species.size.rm,
     m.fitted_data.size.rm,
     m.merged_data.size.rm,
     m.size.rm.fit.mcmc,
     m.size.rm.fit,
     file="m.size.rm.fit.RData")


######################################### By Wing Size and Climate #########################################
m.climate <- as.factor(test$climate)
m.n.climate <- length(unique(test$climate))
m.size <- as.factor(test$wing.size)
m.n.size <- length(unique(test$wing.size))
sink("aedes.by.climate.size.txt")
cat("
model {
  # Priors
  tau1 ~ dexp(0.01)           # Tau for trait
  
  # Size-specific 
  mu.c ~ dexp(10)             # Mean for size effects c
  # Climate-specific
  mu.Topt ~ dnorm(30, 0.1)             # Mean for climate effects Topt
  tau.Topt ~ dexp(0.01)             # Tau for species climate Topt
  mu.a ~ dexp(0.01)            # Mean for climate effects Topt

  for (j in 1:m.n.climate) {
    mu.climate.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Climate-specific temp
    mu.climate.a[j] ~ dexp(mu.a)    # Climate-specific temp
  }
    for (k in 1:m.n.size) {
    mu.size[k] ~ dexp(mu.c)       # Size-specific low point
    size.tau[k] ~ dexp(tau1)
  }
  # Likelihood
  for (i in 1:m.N.obs) {
    mu[i] <- mu.climate.a[m.climate[i]] * (mu.climate.Topt[m.climate[i]] - temp[i])^2 + mu.size[m.size[i]]
    trait[i] ~ dnorm(mu[i], size.tau[m.size[i]])  T(0, )
    log_lik[i] <- logdensity.norm(trait[i], mu[i], size.tau[m.size[i]])
  }}
", file = "aedes.by.climate.size.txt")


# Parameters to monitor
parameters <- c("mu.a", "mu.c", "mu.size", "size.tau", "log_lik",
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

jag.data<-list(trait = test$Trait, m.N.obs = m.N.obs, temp = test$Temp, 
               m.n.climate = m.n.climate, m.climate = m.climate, m.n.size = m.n.size, m.size = m.size)

m.climate.size.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                         model.file="aedes.by.climate.size.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                         n.iter=ni, DIC=T, working.directory=getwd())
m.climate.size.fit.mcmc <- as.mcmc(m.climate.size.fit) ## makes an "mcmc" object
closeAllConnections()
m.climate.size.fit$BUGSoutput$DIC

# Get rid of log liks
# m.climate.size.mcmc <- lapply(as.mcmc.list(m.climate.size.fit.mcmc), function(chain) {
#   chain[, !grepl("log_lik", colnames(chain))]
# })
# 
# # Run Gelman-Rubin diagnostic
# gelman.diag(as.mcmc.list(m.climate.size.mcmc), autoburnin = FALSE)

# Extract log-likelihood matrix (samples × observations)
m.log.lik.climate.size <- m.climate.size.fit$BUGSoutput$sims.list$log_lik
# Compute WAIC
m.waic.climate.size <- waic(m.log.lik.climate.size)$estimates["waic", "Estimate"]

# Get posterior means
m.summary_fit.climate.size <- summary(m.climate.size.fit.mcmc)$statistics
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
m.fitted_data.climate.size <- data.frame(Temperature = rep(m.temp_seq, m.n.climate * m.n.size), 
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
m.fitted_data.climate.size <- data.frame(
  Temperature = rep(m.temp_seq, m.n.climate * m.n.size), 
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
m.fitted_data.climate.size$`Climate, Size` <- with(m.fitted_data.climate.size, paste(Climate, Size, sep = ", "))
m.fitted_data.climate.size$`Climate, Size` <- factor(m.fitted_data.climate.size$`Climate, Size`, levels = c("temperate, large", 
                                                                                        "subtropical, medium", "subtropical, large",
                                                                                        "tropical, small", "tropical, medium"))
# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = m.plot_data.climate.size, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = m.fitted_data.climate.size, aes(x = Temperature, y = Fitted, color = `Climate, Size`), size = 1) +
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
ggsave("mosquito.climate.size.plot.png", width = 8, height = 5)

# Merge the observed and predicted data for comparison
m.merged_data.climate.size <- m.fitted_data.climate.size %>%
  left_join(m.plot_data.climate.size, by = c("Temperature", "Climate", "Species"))


# Calculate RMSE for each species
m.mse_by_species.climate.size <- m.merged_data.climate.size %>%
  group_by(Species) %>%
  summarize(
    MSE = mean((Trait - Fitted)^2, na.rm = TRUE)
  )

# Calculate overall RMSE (no grouping)
m.overall_rmse.climate.size <- m.merged_data.climate.size %>%
  summarize(RMSE = sqrt(mean((Trait - Fitted)^2, na.rm = TRUE))) 



save(m.mse_by_species.climate.size,
     m.fitted_data.climate.size,
     m.merged_data.climate.size,
     m.climate.size.fit.mcmc,
     m.climate.size.fit,
     m.log.lik.climate.size,
     file="m.climate.size.fit.RData")

######################################### RMSE, DIC #########################################
df.aedes <- data.frame(
  Method = c("Full Data", "By Species Individual", "By Species Hierarchical", 
             "By Species Hierarchical, Different Variances",
             "By Climate", "By Wing Size",
             "By Wing Size and Climate"),
  DIC = c(m.one.data.fit$BUGSoutput$DIC, m.DIC.species.individual, 
          m.species.fit$BUGSoutput$DIC,
          m.species.tau.fit$BUGSoutput$DIC, m.climate.fit$BUGSoutput$DIC,
          m.size.fit$BUGSoutput$DIC, 
          m.climate.size.fit$BUGSoutput$DIC),
  wAIC = c(m.waic.full, m.waic.species.ind, 
           m.waic.species, m.waic.species.tau,
           m.waic.climate, m.waic.size,
           m.waic.climate.size),
  Total.RMSE = as.numeric(c(m.overall_rmse, m.overall_rmse.species.ind,
                            m.overall_rmse.species,
                 m.overall_rmse.species.tau, m.overall_rmse.climate,
                 m.overall_rmse.size,
                 m.overall_rmse.climate.size))
)

df.aedes

df.aedes$Method <- reorder(df.aedes$Method, -df.aedes$Total.RMSE)
ggplot(df.aedes, aes(x = Method, y = Total.RMSE, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  labs(title = "RMSE by Method", x = "Method",
       y = "Root Mean Square Error (RMSE)") +
  coord_cartesian(ylim = c(3.5, 3.9)) +
  theme_minimal() +
  scale_fill_viridis_d(option = "magma") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
ggsave("mosquito.RMSE.png", width = 8, height = 5)


df.aedes$Method <- reorder(df.aedes$Method, -df.aedes$DIC)
ggplot(df.aedes, aes(x = Method, y = DIC, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  labs(title = "DIC by Method", x = "Method",
       y = "Deviance Information Criterion (DIC)") +
  theme_minimal() +
  coord_cartesian(ylim = c(5700, 6100)) +
  scale_fill_viridis_d(option = "magma") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
ggsave("mosquito.DIC.png", width = 8, height = 5)

save(df.aedes, file = "Dataset.Mosquito.RData")



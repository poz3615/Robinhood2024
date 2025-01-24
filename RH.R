library(dplyr)
library(tidyr)
library(ggplot2)
library(rjags)
library(R2jags)
library(readxl)
library(truncnorm)
library(scoringRules)

################################### Data Cleaning ###################################
# Read data in
#tick <- read_excel("ovipositionData_Robinhood.xlsx")
tick <- read.csv("Extra.Data.Tick.csv")
tick$Species <- as.factor(tick$Species)
# Choose value, error, species, temp, sample size
tick.abr <- tick[,c("originaltraitvalue", "originalerrorpos",
                    "originalerrorunit", "interactor1",
                    "interactor1temp", "interactor1number")]

# Make species a factor
tick.abr$interactor1 <- as.factor(tick.abr$interactor1)
# If sample size is NA (individual data), make it 1
tick.abr$interactor1number[is.na(tick.abr$interactor1number)] <- 1

# Plot
ggplot(tick, aes(y = Trait, x = Temp)) +
  geom_point() +
  facet_wrap(~ Species) +
  labs(title = "Temperature by Trait for Each Factor Level",
       x = "Temperature",
       y = "Oviposition") +
  theme_minimal()

# Create a new expanded dataset
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
  # Temp
  interactor1temp <- filtered_df1$Temp[i]
  # Sample size
  num_points <- filtered_df1$Count[i]
  # Generate 'num_points' random values from normal distribution
  generated_values <- rtruncnorm(num_points, a=0, b=Inf, mean = mean_value, sd = sd_value)
  # Create a temporary dataframe for the generated values
  temp_df <- data.frame(
    # Simulated values
    originaltraitvalue = generated_values,
    sd = sd_value,
    interactor1 = interactor1,
    interactor1temp = interactor1temp
  )
  # Append the generated rows to the expanded dataset
  expanded_df1 <- rbind(expanded_df1, temp_df)
}

# Now take the individual data
filtered_df <- tick |> filter(Count == 1)
# Summarize it, grouping by species and temperature, to get the mean, sd, and sample size
summarized_df <- filtered_df %>%
  group_by(Species, Temp) %>%
  summarize(
    mean_value = mean(Trait, na.rm = TRUE),
    sd_value = sd(Trait, na.rm = TRUE),
    count = n(),
    .groups = "drop"
  )

# Generate rnorm samples for each group based on sample size, mean, and standard deviation
expanded_df <- summarized_df %>%
  rowwise() %>%
  do(data.frame(
    interactor1 = .$Species,
    interactor1temp = .$Temp,
    originaltraitvalue = rnorm(.$count, mean = .$mean_value, sd = .$sd_value)
  ))

# Combine the new mean data and the new individual data
tick.abr.new <- rbind(expanded_df1[,c(1, 3, 4)], expanded_df)

# Plot new data
ggplot(tick.abr.new, aes(y = originaltraitvalue, x = interactor1temp)) +
  geom_point() +
  facet_wrap(~ interactor1) +
  labs(title = "Temperature by Trait for Each Factor Level",
       x = "Temperature",
       y = "Oviposition") +
  theme_minimal()

# Change species to a factor in new combined cleaned data
tick.abr.new$interactor1 <- as.factor(tick.abr.new$interactor1)
# Species number for hierarchical model
species <- as.numeric(tick.abr.new$interactor1)
# Number of species
n.species <- length(unique(tick.abr.new$interactor1))
# Get rid of NAs
tick.abr.new <- na.omit(tick.abr.new)

################################### Quadratic Model ###################################

sink("quad_robin.txt")
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
    mu.species[j] ~ dnorm(mu.c, tau.c)       # Species-specific vertex y
    mu.species.Topt[j] ~ dnorm(mu.Topt, tau.Topt)    # Species-specific vertex x
    mu.species.a[j] ~ dexp(mu.a)    # Species-specific width
  }

  # Likelihood
  for (i in 1:N.obs) {
    mu[i] <- mu.species.a[species[i]] * (mu.species.Topt[species[i]] - temp[i])^2 + mu.species[species[i]]
    trait[i] ~ dnorm(mu[i], tau1)  
  }}
", file = "quad_robin.txt")

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
parameters <- c("mu.a", "mu.c", "tau.c", "mu.species", 
                "tau1", "mu.Topt", "tau.Topt", "mu.species.Topt", "mu.species.a")


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 60000 # number of iterations in each chain
nb <- 30000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

# Bundle all data in a list for JAGS
jag.data<-list(trait = tick.abr.new$originaltraitvalue, N.obs = N.obs, 
               temp = tick.abr.new$interactor1temp, 
               n.species = n.species, species = species)

quad.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
               model.file="quad_robin.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
               n.iter=ni, DIC=T, working.directory=getwd())
quad.fit.mcmc <- as.mcmc(quad.fit) ## makes an "mcmc" object
closeAllConnections()

# Get posterior means
summary_fit <- summary(quad.fit.mcmc)$statistics
a_mean <- summary_fit[grep("^mu.species.a\\[", rownames(summary_fit)), "Mean"]
species_means <- summary_fit[grep("^mu.species\\[", rownames(summary_fit)), "Mean"]
species_means_Topt <- summary_fit[grep("^mu.species.Topt\\[", rownames(summary_fit)), "Mean"]

# Generate temperature sequence 
temp_seq <- seq(min(tick.abr.new$interactor1temp) ,
                max(tick.abr.new$interactor1temp) , length.out = 100)

# Get predicted values
spec1 <- function (x){
  a_mean[1] * (x - species_means_Topt[1])^2 + species_means[1]
}
spec2 <- function (x){
  a_mean[2] * (x - species_means_Topt[2])^2 + species_means[2]
}
spec3 <- function (x){
  a_mean[3] * (x - species_means_Topt[3])^2 + species_means[3]
}
spec4 <- function (x){
  a_mean[4] * (x - species_means_Topt[4])^2 + species_means[4]
}
spec5 <- function (x){
  a_mean[5] * (x - species_means_Topt[5])^2 + species_means[5]
}
spec6 <- function (x){
  a_mean[6] * (x - species_means_Topt[6])^2 + species_means[6]
}
spec7 <- function (x){
  a_mean[7] * (x - species_means_Topt[7])^2 + species_means[7]
}
spec8 <- function (x){
  a_mean[8] * (x - species_means_Topt[8])^2 + species_means[8]
}
spec9 <- function (x){
  a_mean[9] * (x - species_means_Topt[9])^2 + species_means[9]
}
spec10 <- function (x){
  a_mean[10] * (x - species_means_Topt[10])^2 + species_means[10]
}
spec11 <- function (x){
  a_mean[11] * (x - species_means_Topt[11])^2 + species_means[11]
}
spec12 <- function (x){
  a_mean[12] * (x - species_means_Topt[12])^2 + species_means[12]
}


# Convert to a data frame
fitted_data <- data.frame(Temperature = rep(temp_seq, 10 ),#n.species), 
                          Fitted = c(spec1(temp_seq),
                                     spec2(temp_seq),
                                     spec3(temp_seq),
                                     spec4(temp_seq),
                                     spec5(temp_seq),
                                     spec6(temp_seq),
                                     spec7(temp_seq),
                                     spec8(temp_seq),
                                     spec9(temp_seq),
                                     spec10(temp_seq)),
                                     #spec11(temp_seq),
                                     #spec12(temp_seq)),
                          Species = factor(rep(unique(tick.abr.new$interactor1), each = length(temp_seq))))
# True data points
plot_data <- data.frame(
  Temperature = tick.abr.new$interactor1temp,
  Trait = tick.abr.new$originaltraitvalue,
  Species = tick.abr.new$interactor1 #as.numeric(tick.abr.new$interactor1)
)

# Plot observed data and fitted curves with facet wrap by species
ggplot() +
  geom_point(data = plot_data, aes(x = Temperature, y = Trait), alpha = 0.5) +
  geom_line(data = fitted_data, aes(x = Temperature, y = Fitted, color = Species), linewidth = 1) +
  labs(title = "Observed Data and Quadratic Fit by Species",
       x = "Temperature", y = "Trait Value") +
  theme_minimal() +
  facet_wrap(~ Species) 


################################### Trace, Pairs, Histograms ###################################
df.quad.priors <- as.data.frame(quad.fit.mcmc[[1]][,c(2, 3, 43, 44, 45, 46)])
df.quad.a <- as.data.frame(quad.fit.mcmc[[1]][,c(4:16)])
df.quad.Topt <- as.data.frame(quad.fit.mcmc[[1]][,c(17:29)])
df.quad.c <- as.data.frame(quad.fit.mcmc[[1]][,c(30:42)])
pairs(df.quad.priors)
pairs(df.quad.a)
pairs(df.quad.Topt)
pairs(df.quad.c)
plot(quad.fit.mcmc)

samp.quad <- as.mcmc.list(quad.fit.mcmc)
# Parameter sets
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


################################### RMSE, CRPS, AIC ###################################
# DIC
dic <- quad.fit$BUGSoutput$DIC
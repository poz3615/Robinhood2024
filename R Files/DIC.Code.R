# Not divided by tau!!!!

################################### Ticks
################################### By Species ###################################
# Initialize posterior_samples matrix to get DIC for each species for each sample
n_samples <- nrow(species.fit$BUGSoutput$sims.list$mu.species.a)  # Number of MCMC iterations
# Species names
tick.spec <- c("Amblyomma lepidum", "Dermacentor andersoni",
               "Dermacentor nitens", "Haemaphysalis leporispalustris",
               "Hyalomma aegyptium", "Hyalomma dromedarii",
               "Hyalomma impeltatum", "Hyalomma lusitanicum",
               "Hyalomma schulzei")
# Initialize a matrix to store deviance samples (rows for posterior samples, columns for species)
deviance_samples.species <- matrix(NA, nrow = n_samples, ncol = length(tick.spec))

# Calculate deviance for each sample and each species
for (s in 1:length(tick.spec)) {  
  # Use numeric index for species
  species_name <- tick.spec[s]
  # Which observations are species s
  species_indices <- which(tick.abr.new$Species == species_name)
  # Loop over posterior samples
  for (i in 1:n_samples) {
    # Posterior samples for a, Topt, and c used to make mu
    mu_sample.species <- species.fit$BUGSoutput$sims.list$mu.species.a[i, s] *
      (species.fit$BUGSoutput$sims.list$mu.species.Topt[i, s] - tick.abr.new$Temp[species_indices])^2 +
      species.fit$BUGSoutput$sims.list$mu.species.c[i, s]
    
    # Deviance calculated as -2 * sum of log likelihood - normal with true trait, posterior trait, tau predictor
    deviance_samples.species[i, s] <- -2 * sum(dnorm(tick.abr.new$Trait[species_indices], mean = mu_sample.species,
                                                     sd = sqrt(1 / species.fit$BUGSoutput$sims.list$tau1[i]), log = TRUE))
  }
}

# Initialize Dhat vector to store deviance at posterior mean for each species
Dhat.species <- numeric(length(unique(tick.abr.new$Species)))
tau_mean.species <- mean(species.fit$BUGSoutput$sims.list$tau1)
# Compute Dhat for each species
for (s in 1:length(tick.spec)) {  
  # Use numeric index for species
  species_name <- tick.spec[s]
  # Which observations are species s
  species_indices <- which(tick.abr.new$Species == species_name)
  # Predicted response using posterior mean for parameters for that species
  # Using temperatures from that species
  mu_mean.species <- a_species_mean[s] * (Topt_species_mean[s] - tick.abr.new$Temp[species_indices])^2 + c_species_mean[s]
  Dhat.species[s] <- -2 * sum(dnorm(tick.abr.new$Trait[species_indices], 
                                    mean = mu_mean.species, sd = sqrt(1 / tau_mean.species), log = TRUE))
}

# Compute Dbar for each species (posterior mean of the deviance)
# Deviance samples is a number of samples by 7 matrix for the deviance for each species
# Taking the mean of each column
Dbar.species <- apply(deviance_samples.species, 2, mean) 

pD.species <- Dbar.species - Dhat.species

# Compute DIC for each species
DIC.species <- Dhat.species + 2 * pD.species
# Sum of the DICs to get total DIC
sum(DIC.species)

################################### By Genus ###################################

# Initialize posterior_samples matrix to get DIC for each genus for each sample
n_samples <- nrow(genus.fit$BUGSoutput$sims.list$mu.genus.a)  # Number of MCMC iterations
# Species names
tick.genus <- c("Amblyomma", "Dermacentor", 
                "Haemaphysalis", "Hyalomma")
# Initialize a matrix to store deviance samples (rows for posterior samples, columns for genus)
deviance_samples.genus <- matrix(NA, nrow = n_samples, ncol = length(tick.genus))

# Calculate deviance for each sample and each genus
for (s in 1:length(tick.genus)) {  
  # Use numeric index for species
  genus_name <- tick.genus[s]
  # Which observations are species s
  genus_indices <- which(tick.abr.new$Genus == genus_name)
  # Loop over posterior samples
  for (i in 1:n_samples) {
    # Posterior samples for a, Topt, and c used to make mu
    mu_sample.genus <- genus.fit$BUGSoutput$sims.list$mu.genus.a[i, s] *
      (genus.fit$BUGSoutput$sims.list$mu.genus.Topt[i, s] - tick.abr.new$Temp[genus_indices])^2 +
      genus.fit$BUGSoutput$sims.list$mu.genus.c[i, s]
    
    # Deviance calculated as -2 * sum of log likelihood - normal with true trait, posterior trait, tau predictor
    deviance_samples.genus[i, s] <- -2 * sum(dnorm(tick.abr.new$Trait[genus_indices], mean = mu_sample.genus,
                                                   sd = sqrt(1 / genus.fit$BUGSoutput$sims.list$tau1[i]), log = TRUE))
  }
}

# Initialize Dhat vector to store deviance at posterior mean for each genus
Dhat.genus <- numeric(length(unique(tick.abr.new$Genus)))
tau_mean.genus <- mean(genus.fit$BUGSoutput$sims.list$tau1)
# Compute Dhat for each genus
for (s in 1:length(tick.genus)) {  
  # Use numeric index for genus
  genus_name <- tick.genus[s]
  # Which observations are genus s
  genus_indices <- which(tick.abr.new$Genus == genus_name)
  # Predicted response using posterior mean for parameters for that genus
  # Using temperatures from that genus
  mu_mean.genus <- a_genus_mean[s] * (Topt_genus_mean[s] - tick.abr.new$Temp[genus_indices])^2 + c_genus_mean[s]
  Dhat.genus[s] <- -2 * sum(dnorm(tick.abr.new$Trait[genus_indices], 
                                  mean = mu_mean.genus, sd = sqrt(1 / tau_mean.genus), log = TRUE))
}

# Compute Dbar for each genus (posterior mean of the deviance)
# Deviance samples is a number of samples by 4 matrix for the deviance for each genus
# Taking the mean of each column
Dbar.genus <- apply(deviance_samples.genus, 2, mean) 

pD.genus <- Dbar.genus - Dhat.genus

# Compute DIC for each genus
DIC.genus <- Dhat.genus + 2 * pD.genus
# Sum of the DICs to get total DIC
sum(DIC.genus)


################################### By Climate ###################################

# Initialize posterior_samples matrix to get DIC for each climate for each sample
n_samples <- nrow(climate.fit$BUGSoutput$sims.list$mu.climate.a)  # Number of MCMC iterations
# Species names
tick.climate <- c("Mixed", "Subtropical", 
                  "Temperate", "Tropical")
# Initialize a matrix to store deviance samples (rows for posterior samples, columns for climate)
deviance_samples.climate <- matrix(NA, nrow = n_samples, ncol = length(tick.climate))

# Calculate deviance for each sample and each climate
for (s in 1:length(tick.climate)) {  
  # Use numeric index for species
  climate_name <- tick.climate[s]
  # Which observations are species s
  climate_indices <- which(tick.abr.new$Climate == climate_name)
  # Loop over posterior samples
  for (i in 1:n_samples) {
    # Posterior samples for a, Topt, and c used to make mu
    mu_sample.climate <- climate.fit$BUGSoutput$sims.list$mu.climate.a[i, s] *
      (climate.fit$BUGSoutput$sims.list$mu.climate.Topt[i, s] - tick.abr.new$Temp[climate_indices])^2 +
      climate.fit$BUGSoutput$sims.list$mu.climate.c[i, s]
    
    # Deviance calculated as -2 * sum of log likelihood - normal with true trait, posterior trait, tau predictor
    deviance_samples.climate[i, s] <- -2 * sum(dnorm(tick.abr.new$Trait[climate_indices], mean = mu_sample.climate,
                                                     sd = sqrt(1 / climate.fit$BUGSoutput$sims.list$tau1[i]), log = TRUE))
  }
}

# Initialize Dhat vector to store deviance at posterior mean for each climate
Dhat.climate <- numeric(length(unique(tick.abr.new$Climate)))
tau_mean.climate <- mean(climate.fit$BUGSoutput$sims.list$tau1)
# Compute Dhat for each climate
for (s in 1:length(tick.climate)) {  
  # Use numeric index for climate
  climate_name <- tick.climate[s]
  # Which observations are climate s
  climate_indices <- which(tick.abr.new$Climate == climate_name)
  # Predicted response using posterior mean for parameters for that climate
  # Using temperatures from that climate
  mu_mean.climate <- a_climate_mean[s] * (Topt_climate_mean[s] - tick.abr.new$Temp[climate_indices])^2 + c_climate_mean[s]
  Dhat.climate[s] <- -2 * sum(dnorm(tick.abr.new$Trait[climate_indices], 
                                    mean = mu_mean.climate, sd = sqrt(1 / tau_mean.climate), log = TRUE))
}

# Compute Dbar for each climate (posterior mean of the deviance)
# Deviance samples is a number of samples by 4 matrix for the deviance for each climate
# Taking the mean of each column
Dbar.climate <- apply(deviance_samples.climate, 2, mean) 

pD.climate <- Dbar.climate - Dhat.climate

# Compute DIC for each climate
DIC.climate <- Dhat.climate + 2 * pD.climate
# Sum of the DICs to get total DIC
sum(DIC.climate)

################################### By Host ###################################

# Initialize posterior_samples matrix to get DIC for each host for each sample
n_samples <- nrow(host.fit$BUGSoutput$sims.list$mu.host.a)  # Number of MCMC iterations
# Species names
tick.host <- c("One", "Three", "Two")
# Initialize a matrix to store deviance samples (rows for posterior samples, columns for host)
deviance_samples.host <- matrix(NA, nrow = n_samples, ncol = length(tick.host))

# Calculate deviance for each sample and each host
for (s in 1:length(tick.host)) {  
  # Use numeric index for species
  host_name <- tick.host[s]
  # Which observations are species s
  host_indices <- which(tick.abr.new$Host == host_name)
  # Loop over posterior samples
  for (i in 1:n_samples) {
    # Posterior samples for a, Topt, and c used to make mu
    mu_sample.host <- host.fit$BUGSoutput$sims.list$mu.host.a[i, s] *
      (host.fit$BUGSoutput$sims.list$mu.host.Topt[i, s] - tick.abr.new$Temp[host_indices])^2 +
      host.fit$BUGSoutput$sims.list$mu.host.c[i, s]
    
    # Deviance calculated as -2 * sum of log likelihood - normal with true trait, posterior trait, tau predictor
    deviance_samples.host[i, s] <- -2 * sum(dnorm(tick.abr.new$Trait[host_indices], mean = mu_sample.host,
                                                  sd = sqrt(1 / host.fit$BUGSoutput$sims.list$tau1[i]), log = TRUE))
  }
}

# Initialize Dhat vector to store deviance at posterior mean for each host
Dhat.host <- numeric(length(unique(tick.abr.new$Host)))
tau_mean.host <- mean(host.fit$BUGSoutput$sims.list$tau1)
# Compute Dhat for each host
for (s in 1:length(tick.host)) {  
  # Use numeric index for host
  host_name <- tick.host[s]
  # Which observations are host s
  host_indices <- which(tick.abr.new$Host == host_name)
  # Predicted response using posterior mean for parameters for that host
  # Using temperatures from that host
  mu_mean.host <- a_host_mean[s] * (Topt_host_mean[s] - tick.abr.new$Temp[host_indices])^2 + c_host_mean[s]
  Dhat.host[s] <- -2 * sum(dnorm(tick.abr.new$Trait[host_indices], 
                                 mean = mu_mean.host, sd = sqrt(1 / tau_mean.host), log = TRUE))
}

# Compute Dbar for each host (posterior mean of the deviance)
# Deviance samples is a number of samples by 4 matrix for the deviance for each host
# Taking the mean of each column
Dbar.host <- apply(deviance_samples.host, 2, mean) 

pD.host <- Dbar.host - Dhat.host

# Compute DIC for each host
DIC.host <- Dhat.host + 2 * pD.host
# Sum of the DICs to get total DIC
sum(DIC.host)

################################### By Size ###################################

# Initialize posterior_samples matrix to get DIC for each size for each sample
n_samples <- nrow(size.fit$BUGSoutput$sims.list$mu.size.a)  # Number of MCMC iterations
# Species names
tick.size <- c("Large", "Medium", "Small")
# Initialize a matrix to store deviance samples (rows for posterior samples, columns for size)
deviance_samples.size <- matrix(NA, nrow = n_samples, ncol = length(tick.size))

# Calculate deviance for each sample and each size
for (s in 1:length(tick.size)) {  
  # Use numeric index for species
  size_name <- tick.size[s]
  # Which observations are species s
  size_indices <- which(tick.abr.new$Size == size_name)
  # Loop over posterior samples
  for (i in 1:n_samples) {
    # Posterior samples for a, Topt, and c used to make mu
    mu_sample.size <- size.fit$BUGSoutput$sims.list$mu.size.a[i, s] *
      (size.fit$BUGSoutput$sims.list$mu.size.Topt[i, s] - tick.abr.new$Temp[size_indices])^2 +
      size.fit$BUGSoutput$sims.list$mu.size.c[i, s]
    
    # Deviance calculated as -2 * sum of log likelihood - normal with true trait, posterior trait, tau predictor
    deviance_samples.size[i, s] <- -2 * sum(dnorm(tick.abr.new$Trait[size_indices], mean = mu_sample.size,
                                                  sd = sqrt(1 / size.fit$BUGSoutput$sims.list$tau1[i]), log = TRUE))
  }
}

# Initialize Dhat vector to store deviance at posterior mean for each size
Dhat.size <- numeric(length(unique(tick.abr.new$Size)))
tau_mean.size <- mean(size.fit$BUGSoutput$sims.list$tau1)
# Compute Dhat for each size
for (s in 1:length(tick.size)) {  
  # Use numeric index for size
  size_name <- tick.size[s]
  # Which observations are size s
  size_indices <- which(tick.abr.new$Size == size_name)
  # Predicted response using posterior mean for parameters for that size
  # Using temperatures from that size
  mu_mean.size <- a_size_mean[s] * (Topt_size_mean[s] - tick.abr.new$Temp[size_indices])^2 + c_size_mean[s]
  Dhat.size[s] <- -2 * sum(dnorm(tick.abr.new$Trait[size_indices], 
                                 mean = mu_mean.size, sd = sqrt(1 / tau_mean.size), log = TRUE))
}

# Compute Dbar for each size (posterior mean of the deviance)
# Deviance samples is a number of samples by 4 matrix for the deviance for each size
# Taking the mean of each column
Dbar.size <- apply(deviance_samples.size, 2, mean) 

pD.size <- Dbar.size - Dhat.size

# Compute DIC for each size
DIC.size <- Dhat.size + 2 * pD.size
# Sum of the DICs to get total DIC
sum(DIC.size)



################################### By Continent ###################################

# Initialize posterior_samples matrix to get DIC for each continent for each sample
n_samples <- nrow(continent.fit$BUGSoutput$sims.list$mu.continent.a)  # Number of MCMC iterations
# Species names
tick.continent <- c("More than two", "One", "Two")
# Initialize a matrix to store deviance samples (rows for posterior samples, columns for continent)
deviance_samples.continent <- matrix(NA, nrow = n_samples, ncol = length(tick.continent))

# Calculate deviance for each sample and each continent
for (s in 1:length(tick.continent)) {  
  # Use numeric index for species
  continent_name <- tick.continent[s]
  # Which observations are species s
  continent_indices <- which(tick.abr.new$Continent == continent_name)
  # Loop over posterior samples
  for (i in 1:n_samples) {
    # Posterior samples for a, Topt, and c used to make mu
    mu_sample.continent <- continent.fit$BUGSoutput$sims.list$mu.continent.a[i, s] *
      (continent.fit$BUGSoutput$sims.list$mu.continent.Topt[i, s] - tick.abr.new$Temp[continent_indices])^2 +
      continent.fit$BUGSoutput$sims.list$mu.continent.c[i, s]
    
    # Deviance calculated as -2 * sum of log likelihood - normal with true trait, posterior trait, tau predictor
    deviance_samples.continent[i, s] <- -2 * sum(dnorm(tick.abr.new$Trait[continent_indices], mean = mu_sample.continent,
                                                       sd = sqrt(1 / continent.fit$BUGSoutput$sims.list$tau1[i]), log = TRUE))
  }
}

# Initialize Dhat vector to store deviance at posterior mean for each continent
Dhat.continent <- numeric(length(unique(tick.abr.new$Continent)))
tau_mean.continent <- mean(continent.fit$BUGSoutput$sims.list$tau1)
# Compute Dhat for each continent
for (s in 1:length(tick.continent)) {  
  # Use numeric index for continent
  continent_name <- tick.continent[s]
  # Which observations are continent s
  continent_indices <- which(tick.abr.new$Continent == continent_name)
  # Predicted response using posterior mean for parameters for that continent
  # Using temperatures from that continent
  mu_mean.continent <- a_continent_mean[s] * (Topt_continent_mean[s] - tick.abr.new$Temp[continent_indices])^2 + c_continent_mean[s]
  Dhat.continent[s] <- -2 * sum(dnorm(tick.abr.new$Trait[continent_indices], 
                                      mean = mu_mean.continent, sd = sqrt(1 / tau_mean.continent), log = TRUE))
}

# Compute Dbar for each continent (posterior mean of the deviance)
# Deviance samples is a number of samples by 4 matrix for the deviance for each continent
# Taking the mean of each column
Dbar.continent <- apply(deviance_samples.continent, 2, mean) 

pD.continent <- Dbar.continent - Dhat.continent

# Compute DIC for each continent
DIC.continent <- Dhat.continent + 2 * pD.continent
# Sum of the DICs to get total DIC
sum(DIC.continent)


################################### By Host and Climate ###################################

# Initialize posterior_samples matrix to get DIC for each species for each sample
n_samples <- nrow(host.climate.fit$BUGSoutput$sims.list$mu.host.a)  # Number of MCMC iterations
# Unique climate/size combinations
tick.host <- c("One", "Three", "Two")
tick.climate <- c("Mixed", "Subtropical", "Temperate", "Tropical")
hc.combos <- expand.grid(Host = tick.host, Climate = tick.climate)



# Initialize a matrix to store deviance samples (rows for posterior samples, columns for climate/size combo)
deviance_samples.host.climate <- matrix(NA, nrow = n_samples, ncol = nrow(hc.combos))

# Loop through combinations using sequential indexing
for (index in 1:nrow(hc.combos)) {
  # Use numeric index for size and climate
  host_name <- hc.combos$Host[index]
  h <- as.numeric(hc.combos$Host[index])
  climate_name <- hc.combos$Climate[index]
  c <- as.numeric(hc.combos$Climate[index])
  # Which observations are climate c, size s
  host.climate_indices <- which(tick.abr.new$Climate == climate_name & tick.abr.new$Host == host_name)
  # Loop over posterior samples
  for (i in 1:n_samples) {
    # Compute mu using the correct indices
    mu_sample.host.climate <- host.climate.fit$BUGSoutput$sims.list$mu.host.a[i, h] *
      (host.climate.fit$BUGSoutput$sims.list$mu.climate.Topt[i, c] - tick.abr.new$Temp[host.climate_indices])^2 +
      host.climate.fit$BUGSoutput$sims.list$mu.host.c[i, h]
    # Deviance calculated as -2 * sum of log likelihood - normal with true trait, posterior trait, tau predictor
    deviance_samples.host.climate[i, index] <- -2 * sum(dnorm(tick.abr.new$Trait[host.climate_indices], mean = mu_sample.host.climate,
                                                              sd = sqrt(1 / host.climate.fit$BUGSoutput$sims.list$tau1[i]), log = TRUE))
  }
}

# Initialize Dhat vector to store deviance at posterior mean for each climate.size
Dhat.host.climate <- numeric(nrow(unique(hc.combos)))
tau_mean.host.climate <- mean(host.climate.fit$BUGSoutput$sims.list$tau1)
# Compute Dhat for each climate.size
for (index in 1:nrow(hc.combos)) {
  # Indices for climate.host h
  host_name <- hc.combos$Host[index]
  climate_name <- hc.combos$Climate[index]
  host.climate_indices <- which(tick.abr.new$Climate == climate_name & tick.abr.new$Host == host_name)
  h <- as.numeric(hc.combos$Host[index])
  c <- as.numeric(hc.combos$Climate[index])
  # Predicted response using posterior mean for parameters for that climate.size
  # Using temperatures from that climate.size
  mu_mean.host.climate <- a_host.climate_mean[h] * (Topt_host.climate_mean[c] - tick.abr.new$Temp[host.climate_indices])^2 + c_host.climate_mean[h]
  Dhat.host.climate[index] <- -2 * sum(dnorm(tick.abr.new$Trait[host.climate_indices], mean = mu_mean.host.climate, sd = sqrt(1 / tau_mean.host.climate), log = TRUE))
}

# Compute Dbar for each climate.size (posterior mean of the deviance)
# Deviance samples is a number of samples by 3 matrix for the deviance for each climate.size
# Taking the mean of each column
Dbar.host.climate <- apply(deviance_samples.host.climate, 2, mean) 

pD.host.climate <- Dbar.host.climate - Dhat.host.climate

# Compute DIC for each climate.size
DIC.host.climate <- Dhat.host.climate + 2 * pD.host.climate
# Sum of the DICs to get total DIC
sum(DIC.host.climate)


################################### By Host and Genus ###################################

# Initialize posterior_samples matrix to get DIC for each species for each sample
n_samples <- nrow(host.genus.fit$BUGSoutput$sims.list$mu.host.a)  # Number of MCMC iterations
# Unique host/genus combinations
tick.host <- c("One", "Three", "Two")
tick.genus <- c("Amblyomma", "Dermacentor", 
                "Haemaphysalis", "Hyalomma")
hg.combos <- expand.grid(Host = tick.host, Genus = tick.genus)

tick.abr.new$Genus <- droplevels(tick.abr.new$Genus)

# Initialize a matrix to store deviance samples (rows for posterior samples, columns for climate/size combo)
deviance_samples.host.genus <- matrix(NA, nrow = n_samples, ncol = nrow(hg.combos))

# Loop through combinations using sequential indexing
for (index in 1:nrow(hg.combos)) {
  # Use numeric index for host and genus
  host_name <- hg.combos$Host[index]
  h <- as.numeric(hg.combos$Host[index])
  genus_name <- hg.combos$Genus[index]
  g <- as.numeric(hg.combos$Genus[index])
  # Which observations are genus g host h
  host.genus_indices <- which(tick.abr.new$Genus == genus_name & tick.abr.new$Host == host_name)
  # Loop over posterior samples
  for (i in 1:n_samples) {
    # Compute mu using the correct indices
    mu_sample.host.genus <- host.genus.fit$BUGSoutput$sims.list$mu.host.a[i, h] *
      (host.genus.fit$BUGSoutput$sims.list$mu.genus.Topt[i, g] - tick.abr.new$Temp[host.genus_indices])^2 +
      host.genus.fit$BUGSoutput$sims.list$mu.genus.c[i, g]
    # Deviance calculated as -2 * sum of log likelihood - normal with true trait, posterior trait, tau predictor
    deviance_samples.host.genus[i, index] <- -2 * sum(dnorm(tick.abr.new$Trait[host.genus_indices], mean = mu_sample.host.genus,
                                                            sd = sqrt(1 / host.genus.fit$BUGSoutput$sims.list$tau1[i]), log = TRUE))
  }
}

# Initialize Dhat vector to store deviance at posterior mean for each host/genus
Dhat.host.genus <- numeric(nrow(unique(hg.combos)))
tau_mean.host.genus <- mean(host.genus.fit$BUGSoutput$sims.list$tau1)
# Compute Dhat for each host/genus
for (index in 1:nrow(hg.combos)) {
  # Indices for genus.host h
  host_name <- hg.combos$Host[index]
  genus_name <- hg.combos$Genus[index]
  host.genus_indices <- which(tick.abr.new$Genus == genus_name & tick.abr.new$Host == host_name)
  h <- as.numeric(hg.combos$Host[index])
  g <- as.numeric(hg.combos$Genus[index])
  # Predicted response using posterior mean for parameters for that host/genus
  # Using temperatures from that host/genus
  mu_mean.host.genus <- a_host.genus_mean[h] * (Topt_host.genus_mean[g] - tick.abr.new$Temp[host.genus_indices])^2 + c_host.genus_mean[g]
  Dhat.host.genus[index] <- -2 * sum(dnorm(tick.abr.new$Trait[host.genus_indices], mean = mu_mean.host.genus, sd = sqrt(1 / tau_mean.host.genus), log = TRUE))
}

# Compute Dbar for each host/genus (posterior mean of the deviance)
# Deviance samples is a number of samples by 3 matrix for the deviance for each host/genus
# Taking the mean of each column
Dbar.host.genus <- apply(deviance_samples.host.genus, 2, mean) 

pD.host.genus <- Dbar.host.genus - Dhat.host.genus

# Compute DIC for each host/genus
DIC.host.genus <- Dhat.host.genus + 2 * pD.host.genus
# Sum of the DICs to get total DIC
sum(DIC.host.genus)

################################### By Genus and Climate ###################################

# Initialize posterior_samples matrix to get DIC for each species for each sample
n_samples <- nrow(genus.climate.fit$BUGSoutput$sims.list$mu.climate.a)  # Number of MCMC iterations
# Unique climate/size combinations
tick.genus <- c("Amblyomma", "Dermacentor", 
                "Haemaphysalis", "Hyalomma")
tick.climate <- c("Mixed", "Subtropical", "Temperate", "Tropical")
gc.combos <- expand.grid(Genus = tick.genus, Climate = tick.climate)



# Initialize a matrix to store deviance samples (rows for posterior samples, columns for genus/climate combo)
deviance_samples.genus.climate <- matrix(NA, nrow = n_samples, ncol = nrow(gc.combos))

# Loop through combinations using sequential indexing
for (index in 1:nrow(gc.combos)) {
  # Use numeric index for genus and climate
  genus_name <- gc.combos$Genus[index]
  g <- as.numeric(gc.combos$Genus[index])
  climate_name <- gc.combos$Climate[index]
  c <- as.numeric(gc.combos$Climate[index])
  # Which observations are genus g, size s
  genus.climate_indices <- which(tick.abr.new$Climate == climate_name & tick.abr.new$Genus == genus_name)
  # Loop over posterior samples
  for (i in 1:n_samples) {
    # Compute mu using the correct indices
    mu_sample.genus.climate <- genus.climate.fit$BUGSoutput$sims.list$mu.climate.a[i, c] *
      (genus.climate.fit$BUGSoutput$sims.list$mu.genus.Topt[i, g] - tick.abr.new$Temp[genus.climate_indices])^2 +
      genus.climate.fit$BUGSoutput$sims.list$mu.genus.c[i, g]
    # Deviance calculated as -2 * sum of log likelihood - normal with true trait, posterior trait, tau predictor
    deviance_samples.genus.climate[i, index] <- -2 * sum(dnorm(tick.abr.new$Trait[genus.climate_indices], mean = mu_sample.genus.climate,
                                                               sd = sqrt(1 / genus.climate.fit$BUGSoutput$sims.list$tau1[i]), log = TRUE))
  }
}

# Initialize Dhat vector to store deviance at posterior mean for each genus/host
Dhat.genus.climate <- numeric(nrow(unique(gc.combos)))
tau_mean.genus.climate <- mean(genus.climate.fit$BUGSoutput$sims.list$tau1)
# Compute Dhat for each genus/host
for (index in 1:nrow(gc.combos)) {
  # Indices for climate.host h
  genus_name <- gc.combos$Genus[index]
  climate_name <- gc.combos$Climate[index]
  genus.climate_indices <- which(tick.abr.new$Climate == climate_name & tick.abr.new$Genus == genus_name)
  g <- as.numeric(gc.combos$Genus[index])
  c <- as.numeric(gc.combos$Climate[index])
  # Predicted response using posterior mean for parameters for that genus/climate
  # Using temperatures from that genus/climate
  mu_mean.genus.climate <- a_genus.climate_mean[c] * (Topt_genus.climate_mean[g] - tick.abr.new$Temp[genus.climate_indices])^2 + c_genus.climate_mean[g]
  Dhat.genus.climate[index] <- -2 * sum(dnorm(tick.abr.new$Trait[genus.climate_indices], mean = mu_mean.genus.climate, sd = sqrt(1 / tau_mean.genus.climate), log = TRUE))
}

# Compute Dbar for each genus/climate (posterior mean of the deviance)
# Deviance samples is a number of samples by 3 matrix for the deviance for each genus/climate
# Taking the mean of each column
Dbar.genus.climate <- apply(deviance_samples.genus.climate, 2, mean) 

pD.genus.climate <- Dbar.genus.climate - Dhat.genus.climate

# Compute DIC for each genus/climate
DIC.genus.climate <- Dhat.genus.climate + 2 * pD.genus.climate
# Sum of the DICs to get total DIC
sum(DIC.genus.climate)



################################### By Host, Climate, and Size ###################################

# Initialize posterior_samples matrix to get DIC for each species for each sample
n_samples <- nrow(h.c.s.fit$BUGSoutput$sims.list$mu.host.a)  # Number of MCMC iterations
# Unique climate/size combinations
tick.host <- c("One", "Three", "Two")
tick.climate <- c("Mixed", "Subtropical", "Temperate", "Tropical")
tick.size <- c("Large", "Medium", "Small")
h.c.s.combos <- expand.grid(Host = tick.host, Climate = tick.climate,
                            Size = tick.size)



# Initialize a matrix to store deviance samples (rows for posterior samples, columns for climate/size combo)
deviance_samples.h.c.s <- matrix(NA, nrow = n_samples, ncol = nrow(h.c.s.combos))

# Loop through combinations using sequential indexing
for (index in 1:nrow(h.c.s.combos)) {
  # Use numeric index for host, size, and climate
  host_name <- h.c.s.combos$Host[index]
  h <- as.numeric(h.c.s.combos$Host[index])
  climate_name <- h.c.s.combos$Climate[index]
  c <- as.numeric(h.c.s.combos$Climate[index])
  size_name <- h.c.s.combos$Size[index]
  s <- as.numeric(h.c.s.combos$Size[index])
  # Which observations are climate c, size s, host h
  h.c.s_indices <- which(tick.abr.new$Climate == climate_name & 
                           tick.abr.new$Host == host_name & tick.abr.new$Size == size_name)
  # Loop over posterior samples
  for (i in 1:n_samples) {
    # Compute mu using the correct indices
    mu_sample.h.c.s <- h.c.s.fit$BUGSoutput$sims.list$mu.host.a[i, h] *
      (h.c.s.fit$BUGSoutput$sims.list$mu.climate.Topt[i, c] - tick.abr.new$Temp[h.c.s_indices])^2 +
      h.c.s.fit$BUGSoutput$sims.list$mu.size.c[i, s]
    # Deviance calculated as -2 * sum of log likelihood - normal with true trait, posterior trait, tau predictor
    deviance_samples.h.c.s[i, index] <- -2 * sum(dnorm(tick.abr.new$Trait[h.c.s_indices], mean = mu_sample.h.c.s,
                                                       sd = sqrt(1 / h.c.s.fit$BUGSoutput$sims.list$tau1[i]), log = TRUE))
  }
}

# Initialize Dhat vector to store deviance at posterior mean for each climate.size
Dhat.h.c.s <- numeric(nrow(unique(h.c.s.combos)))
tau_mean.h.c.s <- mean(h.c.s.fit$BUGSoutput$sims.list$tau1)
# Compute Dhat for each climate.size
for (index in 1:nrow(h.c.s.combos)) {
  # Indices for climate.host.size 
  host_name <- h.c.s.combos$Host[index]
  climate_name <- h.c.s.combos$Climate[index]
  size_name <- h.c.s.combos$Size[index]
  h.c.s_indices <- which(tick.abr.new$Climate == climate_name & 
                           tick.abr.new$Host == host_name & tick.abr.new$Size == size_name)
  h <- as.numeric(h.c.s.combos$Host[index])
  c <- as.numeric(h.c.s.combos$Climate[index])
  s <- as.numeric(h.c.s.combos$Size[index])
  # Predicted response using posterior mean for parameters for that climate.size.host
  # Using temperatures from that climate.size
  mu_mean.h.c.s <- a_h.c.s_mean[h] * (Topt_h.c.s_mean[c] - tick.abr.new$Temp[h.c.s_indices])^2 + c_h.c.s_mean[s]
  Dhat.h.c.s[index] <- -2 * sum(dnorm(tick.abr.new$Trait[h.c.s_indices], mean = mu_mean.h.c.s, sd = sqrt(1 / tau_mean.h.c.s), log = TRUE))
}

# Compute Dbar for each climate.size.host (posterior mean of the deviance)
# Deviance samples is a number of samples by 3 matrix for the deviance for each climate.size.host
# Taking the mean of each column
Dbar.h.c.s <- apply(deviance_samples.h.c.s, 2, mean) 

pD.h.c.s <- Dbar.h.c.s - Dhat.h.c.s

# Compute DIC for each climate.size
DIC.h.c.s <- Dhat.h.c.s + 2 * pD.h.c.s
# Sum of the DICs to get total DIC
sum(DIC.h.c.s)







###################################
###################################
###################################
################################### Mosquitoes
################################### Full Data ###################################
# Initialize posterior_samples matrix to get DIC for each species for each sample
m.n_samples <- nrow(m.one.data.fit$BUGSoutput$sims.list$a)  # Number of MCMC iterations
# Species names

m.deviance_samples <- matrix(NA, nrow = m.n_samples, ncol = 1)

# Calculate deviance for each sample and each species

for (i in 1:m.n_samples) {
  # Posterior samples for a, Topt, and c used to make mu
  m.mu_sample <- m.one.data.fit$BUGSoutput$sims.list$a[i,] *
    (m.one.data.fit$BUGSoutput$sims.list$Topt[i,] - test$Temp)^2 +
    m.one.data.fit$BUGSoutput$sims.list$c[i,]
  
  # Deviance calculated as -2 * sum of log likelihood - normal with true trait, posterior trait, tau predictor
  m.deviance_samples[i,] <- -2 * sum(dnorm(test$Trait, mean = m.mu_sample,
                                           sd = sqrt(1 / m.one.data.fit$BUGSoutput$sims.list$tau1[i]), log = TRUE))
}


# Initialize Dhat vector to store deviance at posterior mean for each species
m.Dhat <- c()
m.tau_mean <- mean(m.one.data.fit$BUGSoutput$sims.list$tau1)
# Compute Dhat for each species

m.mu_mean <- m.a_mean * (m.means_Topt - test$Temp)^2 + m.c_mean
m.Dhat <- -2 * sum(dnorm(test$Trait, 
                         mean = m.mu_mean, sd = sqrt(1 / m.tau_mean), log = TRUE))

# Compute Dbar for each species (posterior mean of the deviance)
# Deviance samples is a number of samples by 7 matrix for the deviance for each species
# Taking the mean of each column
m.Dbar <- mean(m.deviance_samples)

m.pD <- m.Dbar - m.Dhat

# Compute DIC for each species
m.DIC <- m.Dhat + 2 * m.pD
# Sum of the DICs to get total DIC
sum(m.DIC)

######################################### By Species #########################################

# Initialize posterior_samples matrix to get DIC for each species for each sample
m.n_samples <- nrow(m.species.fit$BUGSoutput$sims.list$mu.species.a)  # Number of MCMC iterations
# Species names
mos.spec <- c("Aedes aegypti", "Aedes albopictus",
              "Aedes atropalpus", "Aedes camptorhynchus",
              "Aedes japonicus japonicus", "Aedes krombeini", "Aedes notoscriptus")
# Initialize a matrix to store deviance samples (rows for posterior samples, columns for species)
m.deviance_samples.species <- matrix(NA, nrow = m.n_samples, ncol = length(mos.spec))

# Calculate deviance for each sample and each species
for (s in 1:length(mos.spec)) {  
  # Use numeric index for species
  m.species_name <- mos.spec[s]
  # Which observations are species s
  m.species_indices <- which(test$interactor1 == m.species_name)
  # Loop over posterior samples
  for (i in 1:m.n_samples) {
    # Posterior samples for a, Topt, and c used to make mu
    m.mu_sample.species <- m.species.fit$BUGSoutput$sims.list$mu.species.a[i, s] *
      (m.species.fit$BUGSoutput$sims.list$mu.species.Topt[i, s] - test$Temp[m.species_indices])^2 +
      m.species.fit$BUGSoutput$sims.list$mu.species[i, s]
    
    # Deviance calculated as -2 * sum of log likelihood - normal with true trait, posterior trait, tau predictor
    m.deviance_samples.species[i, s] <- -2 * sum(dnorm(test$Trait[m.species_indices], mean = m.mu_sample.species,
                                                       sd = sqrt(1 / m.species.fit$BUGSoutput$sims.list$tau1[i]), log = TRUE))
  }
}

# Initialize Dhat vector to store deviance at posterior mean for each species
m.Dhat.species <- numeric(length(unique(test$interactor1)))
m.tau_mean.species <- mean(m.species.fit$BUGSoutput$sims.list$tau1)
# Compute Dhat for each species
for (s in 1:length(mos.spec)) {  
  # Use numeric index for species
  m.species_name <- mos.spec[s]
  # Which observations are species s
  m.species_indices <- which(test$interactor1 == m.species_name)
  # Predicted response using posterior mean for parameters for that species
  # Using temperatures from that species
  m.mu_mean.species <- m.a_mean_species[s] * (m.Topt_mean_species[s] - test$Temp[m.species_indices])^2 + m.c_mean_species[s]
  m.Dhat.species[s] <- -2 * sum(dnorm(test$Trait[m.species_indices], 
                                      mean = m.mu_mean.species, sd = sqrt(1 / m.tau_mean.species), log = TRUE))
}

# Compute Dbar for each species (posterior mean of the deviance)
# Deviance samples is a number of samples by 7 matrix for the deviance for each species
# Taking the mean of each column
m.Dbar.species <- apply(m.deviance_samples.species, 2, mean) 

m.pD.species <- m.Dbar.species - m.Dhat.species

# Compute DIC for each species
m.DIC.species <- m.Dhat.species + 2 * m.pD.species
# Sum of the DICs to get total DIC
sum(m.DIC.species)


######################################### By Climate #########################################

# Initialize posterior_samples matrix to get DIC for each species for each sample
m.n_samples <- nrow(m.climate.fit$BUGSoutput$sims.list$mu.climate.a)  # Number of MCMC iterations
# Species names
mos.climate <- c("subtropical", "temperate", "tropical")
# Initialize a matrix to store deviance samples (rows for posterior samples, columns for climate)
m.deviance_samples.climate <- matrix(NA, nrow = m.n_samples, ncol = length(mos.climate))

# Calculate deviance for each sample and each climate
for (s in 1:length(mos.climate)) {  
  # Use numeric index for climate
  m.climate_name <- mos.climate[s]
  # Which observations are climate s
  m.climate_indices <- which(test$climate == m.climate_name)
  # Loop over posterior samples
  for (i in 1:n_samples) {
    # Posterior samples for a, Topt, and c used to make mu
    m.mu_sample.climate <- m.climate.fit$BUGSoutput$sims.list$mu.climate.a[i, s] *
      (m.climate.fit$BUGSoutput$sims.list$mu.climate.Topt[i, s] - test$Temp[m.climate_indices])^2 +
      m.climate.fit$BUGSoutput$sims.list$mu.climate[i, s]
    
    # Deviance calculated as -2 * sum of log likelihood - normal with true trait, posterior trait, tau predictor
    m.deviance_samples.climate[i, s] <- -2 * sum(dnorm(test$Trait[m.climate_indices], mean = m.mu_sample.climate,
                                                       sd = sqrt(1 / m.climate.fit$BUGSoutput$sims.list$tau1[i]), log = TRUE))
  }
}

# Initialize Dhat vector to store deviance at posterior mean for each climate
m.Dhat.climate <- numeric(length(unique(test$climate)))
m.tau_mean.climate <- mean(m.climate.fit$BUGSoutput$sims.list$tau1)
# Compute Dhat for each climate
for (s in 1:length(unique(test$climate))) {
  # Indices for climate s
  m.climate_indices <- which(test$climate == unique(test$climate)[s]) 
  # Predicted response using posterior mean for parameters for that climate
  # Using temperatures from that climate
  m.mu_mean.climate <- m.a_mean_climate[s] * (m.Topt_mean_climate[s] - test$Temp[m.climate_indices])^2 + m.c_mean_climate[s]
  m.Dhat.climate[s] <- -2 * sum(dnorm(test$Trait[m.climate_indices], mean = m.mu_mean.climate, sd = sqrt(1 / m.tau_mean.climate), log = TRUE))
}

# Compute Dbar for each climate (posterior mean of the deviance)
# Deviance samples is a number of samples by 3 matrix for the deviance for each climate
# Taking the mean of each column
m.Dbar.climate <- apply(m.deviance_samples.climate, 2, mean) 

m.pD.climate <- m.Dbar.climate - m.Dhat.climate

# Compute DIC for each climate
m.DIC.climate <- m.Dhat.climate + 2 * m.pD.climate
# Sum of the DICs to get total DIC
sum(m.DIC.climate)

######################################### By Wing Size Category #########################################

# Initialize deviance_samples matrix to get DIC for each species for each sample
m.n_samples <- nrow(m.size.fit$BUGSoutput$sims.list$mu.size.a)  # Number of MCMC iterations
# Species names
mos.size <- c("large", "medium", "small")
# Initialize a matrix to store deviance samples (rows for posterior samples, columns for size)
m.deviance_samples.size <- matrix(NA, nrow = m.n_samples, ncol = length(mos.size))

# Calculate deviance for each sample and each size
for (s in 1:length(mos.size)) {  
  # Use numeric index for size
  m.size_name <- mos.size[s]
  # Which observations are size s
  m.size_indices <- which(test$wing.size == m.size_name)
  # Loop over posterior samples
  for (i in 1:m.n_samples) {
    # Posterior samples for a, Topt, and c used to make mu
    m.mu_sample.size <- m.size.fit$BUGSoutput$sims.list$mu.size.a[i, s] *
      (m.size.fit$BUGSoutput$sims.list$mu.size.Topt[i, s] - test$Temp[m.size_indices])^2 +
      m.size.fit$BUGSoutput$sims.list$mu.size[i, s]
    
    # Deviance calculated as -2 * sum of log likelihood - normal with true trait, posterior trait, tau predictor
    m.deviance_samples.size[i, s] <- -2 * sum(dnorm(test$Trait[m.size_indices], mean = m.mu_sample.size,
                                                    sd = sqrt(1 / m.size.fit$BUGSoutput$sims.list$tau1[i]), log = TRUE))
  }
}

# Initialize Dhat vector to store deviance at posterior mean for each size
m.Dhat.size <- numeric(length(unique(test$wing.size)))
m.tau_mean.size <- mean(m.size.fit$BUGSoutput$sims.list$tau1)
# Compute Dhat for each size
for (s in 1:length(unique(test$wing.size))) {
  # Indices for size s
  m.size_indices <- which(test$wing.size == unique(test$wing.size)[s]) 
  # Predicted response using posterior mean for parameters for that size
  # Using temperatures from that size
  m.mu_mean.size <- m.a_mean_size[s] * (m.Topt_mean_size[s] - test$Temp[m.size_indices])^2 + m.c_mean_size[s]
  m.Dhat.size[s] <- -2 * sum(dnorm(test$Trait[m.size_indices], mean = m.mu_mean.size, sd = sqrt(1 / m.tau_mean.size), log = TRUE))
}

# Compute Dbar for each size (posterior mean of the deviance)
# Deviance samples is a number of samples by 3 matrix for the deviance for each size
# Taking the mean of each column
m.Dbar.size <- apply(m.deviance_samples.size, 2, mean) 

m.pD.size <- m.Dbar.size - m.Dhat.size

# Compute DIC for each size
m.DIC.size <- m.Dhat.size + 2 * m.pD.size
# Sum of the DICs to get total DIC
sum(m.DIC.size)

######################################### By Wing Size Category Removed #########################################

# Initialize deviance_samples matrix to get DIC for each species for each sample
m.n_samples <- nrow(m.size.rm.fit$BUGSoutput$sims.list$mu.size.a)  # Number of MCMC iterations
# Species names
mos.size <- c("large", "medium", "small")
# Initialize a matrix to store deviance samples (rows for posterior samples, columns for size)
m.deviance_samples.size.rm <- matrix(NA, nrow = m.n_samples, ncol = length(mos.size))

# Calculate deviance for each sample and each size
for (s in 1:length(mos.size)) {  
  # Use numeric index for size
  m.size.rm_name <- mos.size[s]
  # Which observations are size s
  m.size.rm_indices <- which(test1$wing.size == m.size.rm_name)
  # Loop over posterior samples
  for (i in 1:m.n_samples) {
    # Posterior samples for a, Topt, and c used to make mu
    m.mu_sample.size.rm <- m.size.rm.fit$BUGSoutput$sims.list$mu.size.a[i, s] *
      (m.size.rm.fit$BUGSoutput$sims.list$mu.size.Topt[i, s] - test1$Temp[m.size.rm_indices])^2 +
      m.size.rm.fit$BUGSoutput$sims.list$mu.size[i, s]
    
    # Deviance calculated as -2 * sum of log likelihood - normal with true trait, posterior trait, tau predictor
    m.deviance_samples.size.rm[i, s] <- -2 * sum(dnorm(test1$Trait[m.size.rm_indices], mean = m.mu_sample.size.rm,
                                                       sd = sqrt(1 / m.size.rm.fit$BUGSoutput$sims.list$tau1[i]), log = TRUE))
  }
}

# Initialize Dhat vector to store deviance at posterior mean for each size
m.Dhat.size.rm <- numeric(length(unique(test1$wing.size)))
m.tau_mean.size.rm <- mean(m.size.rm.fit$BUGSoutput$sims.list$tau1)
# Compute Dhat for each size
for (s in 1:length(unique(test1$wing.size))) {
  # Indices for size s
  m.size.rm_indices <- which(test1$wing.size == unique(test1$wing.size)[s]) 
  # Predicted response using posterior mean for parameters for that size
  # Using temperatures from that size
  m.mu_mean.size.rm <- m.a_mean_size.rm[s] * (m.Topt_mean_size.rm[s] - test1$Temp[m.size_indices])^2 + m.c_mean_size.rm[s]
  m.Dhat.size.rm[s] <- -2 * sum(dnorm(test1$Trait[m.size.rm_indices], mean = m.mu_mean.size.rm, sd = sqrt(1 / m.tau_mean.size.rm), log = TRUE))
}

# Compute Dbar for each size (posterior mean of the deviance)
# Deviance samples is a number of samples by 3 matrix for the deviance for each size
# Taking the mean of each column
m.Dbar.size.rm <- apply(m.deviance_samples.size.rm, 2, mean) 

m.pD.size.rm <- m.Dbar.size.rm - m.Dhat.size.rm

# Compute DIC for each size
m.DIC.size.rm <- m.Dhat.size.rm + 2 * m.pD.size.rm
# Sum of the DICs to get total DIC
sum(m.DIC.size.rm)

######################################### By Wing Size and Climate #########################################

# Initialize posterior_samples matrix to get DIC for each species for each sample
m.n_samples <- nrow(m.climate.size.fit$BUGSoutput$sims.list$mu.climate.a)  # Number of MCMC iterations
# Unique climate/size combinations
sc.combos <- expand.grid(size = mos.size, climate = mos.climate)

# Initialize a matrix to store deviance samples (rows for posterior samples, columns for climate/size combo)
m.deviance_samples.climate.size <- matrix(NA, nrow = m.n_samples, ncol = nrow(combos))

# Loop through combinations using sequential indexing
for (index in 1:nrow(sc.combos)) {
  # Use numeric index for size and climate
  m.size_name <- sc.combos$size[index]
  s <- as.numeric(sc.combos$size[index])
  m.climate_name <- sc.combos$climate[index]
  c <- as.numeric(sc.combos$climate[index])
  # Which observations are climate c, size s
  m.climate.size_indices <- which(test$climate == m.climate_name & test$wing.size == m.size_name)
  # Loop over posterior samples
  for (i in 1:m.n_samples) {
    # Compute mu using the correct indices
    m.mu_sample.climate.size <- m.climate.size.fit$BUGSoutput$sims.list$mu.climate.a[i, c] *
      (m.climate.size.fit$BUGSoutput$sims.list$mu.climate.Topt[i, c] - test$Temp[m.climate.size_indices])^2 +
      m.climate.size.fit$BUGSoutput$sims.list$mu.size[i, s]
    # Deviance calculated as -2 * sum of log likelihood - normal with true trait, posterior trait, tau predictor
    m.deviance_samples.climate.size[i, index] <- -2 * sum(dnorm(test$Trait[m.climate.size_indices], mean = m.mu_sample.climate.size,
                                                                sd = sqrt(1 / m.climate.size.fit$BUGSoutput$sims.list$tau1[i]), log = TRUE))
  }
}

# Initialize Dhat vector to store deviance at posterior mean for each climate.size
m.Dhat.climate.size <- numeric(nrow(unique(sc.combos)))
m.tau_mean.climate.size <- mean(m.climate.size.fit$BUGSoutput$sims.list$tau1)
# Compute Dhat for each climate.size
for (index in 1:nrow(sc.combos)) {
  # Indices for climate.size s
  m.size_name <- sc.combos$size[index]
  m.climate_name <- sc.combos$climate[index]
  m.climate.size_indices <- which(test$climate == m.climate_name & test$wing.size == m.size_name)
  s <- as.numeric(sc.combos$size[index])
  c <- as.numeric(sc.combos$climate[index])
  # Predicted response using posterior mean for parameters for that climate.size
  # Using temperatures from that climate.size
  m.mu_mean.climate.size <- m.a_mean_climate.size[c] * (m.Topt_mean_climate.size[c] - test$Temp[m.climate.size_indices])^2 + m.c_mean_climate.size[s]
  m.Dhat.climate.size[index] <- -2 * sum(dnorm(test$Trait[m.climate.size_indices], mean = m.mu_mean.climate.size, sd = sqrt(1 / m.tau_mean.climate.size), log = TRUE))
}

# Compute Dbar for each climate.size (posterior mean of the deviance)
# Deviance samples is a number of samples by 3 matrix for the deviance for each climate.size
# Taking the mean of each column
m.Dbar.climate.size <- apply(m.deviance_samples.climate.size, 2, mean) 

m.pD.climate.size <- m.Dbar.climate.size - m.Dhat.climate.size

# Compute DIC for each climate.size
m.DIC.climate.size <- m.Dhat.climate.size + 2 * m.pD.climate.size
# Sum of the DICs to get total DIC
sum(m.DIC.climate.size)


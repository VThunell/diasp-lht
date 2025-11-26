# Load libraries, install if needed
pkgs <- c("tidyverse", "tidybayes","bayesplot", "nimbleHMC")

# remotes::install_github("nimble-dev/nimble", ref = "devel", subdir = "packages/nimble")

if(length(setdiff(pkgs,rownames(installed.packages()))) > 0){
  install.packages(setdiff(pkgs, rownames(installed.packages())), dependencies = T)
}

invisible(lapply(pkgs, library, character.only = T))

fec.code <- nimbleCode({
  
  # likelihood
  for(i in 1:nobs){
    
    n.eggs[i] ~ dnorm(mu[i], sd = sigma)
    mu[i] <- a[su[i]]*length_mm[i]^b[su[i]]
    #lognormal?
  }
  
  # priors
  sigma ~ dexp(1/10000)
  
  for(j in 1:nloc){
    #a[j] ~ dnorm(1, .1) #poor convergence but lower AIC cf dlnomr
    #a[j] ~ dlnorm(1, .1) # better conv. but lower AIC than dgammma
    a[j] ~ dgamma(mean = 1, sd = 1)
    b[j] ~ dnorm(0, 1)
  }
  
})

data.fec <- salfec %>%
  mutate(su = as.integer(factor(spat.unit))) %>%
  drop_na(su)
  # mutate(sai_location.int = as.integer(factor(sai_location))) %>%
  # drop_na(sai_location)
nloc <- length(unique(data.fec$su))

#initial values generating function
inits <- function(){
  list(a = rgamma(nloc,1,1),
       b = rnorm(nloc,0,1),
       sigma = rlnorm(1,log(1000),0.1))}

# build model
fec.model <- nimbleModel(fec.code,
                        constants = list(nloc = nloc,
                                         nobs = nrow(data.fec),
                                         su = data.fec$su),
                        inits=inits(),
                        data = data.fec %>% select(n.eggs, length_mm),
                        buildDerivs = TRUE)

fec.model$simulate()
fec.model$calculate()

# configure hmc
fec.confhmc <- configureHMC(fec.model, 
                            monitors = c("a", "b","sigma","mu"),
                            enableWAIC = TRUE)

# build mcmc (use buidlHMC() when not using configureHMC())
fec.hmc <- buildMCMC(fec.confhmc)

# compile model
fec.c <- compileNimble(fec.model)

# compile mcmc  and specify the project model
fec.hmcc <- compileNimble(fec.hmc)

# hmc samples
fec.samples <- runMCMC(fec.hmcc, niter = 16000, nburnin = 10000, thin = 3, nchains = 2, WAIC=TRUE, samplesAsCodaMCMC = TRUE)

# Save samples
saveRDS(fec.samples, file = paste0(home,"/data/fec_samples_1104.RData"))

# calculate summary and save those for temporal biphasic model fec predictions
fec.dpars <- mcmc.list(lapply(fec.samples$samples, as.mcmc)) %>% 
  spread_draws(b[su],a[su]) %>%
  ungroup() %>%
  summarise(mufb = mean(b),  sdfb = sd(b), 
            ma = mean(a), va  = var(a),
            shfa  = ma^2 / va, rafa   = ma / va,
            .by = su)

saveRDS(fec.dpars, file = paste0(home,"/data/fec_dpars.RData"))


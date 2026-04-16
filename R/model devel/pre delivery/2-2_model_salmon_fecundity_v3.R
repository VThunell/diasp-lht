# testing negative binomial as in de Eyto (because of overdipsersion)? 

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
    
    n.eggs[i] ~ dnegbin(prob = p[i], size = r)
    p[i] <- r/(r+mu[i])
    log(mu[i]) <- a[su[i]] + b[su[i]] * length.lc[i]
    
  }
  
  # priors
  r ~ dgamma(shape = 1, rate = 1)
  
  for(j in 1:nsu){
    a[j] ~ dnorm(mu.a, sd = sd.a)
    b[j] ~ dnorm(mu.b, sd = sd.b)
  }
  
  mu.a ~ dnorm(0, 1)
  sd.a ~ dexp(0.5)
  mu.b ~ dnorm(1, 1)
  sd.b ~ dexp(0.5)

})


nsu <- length(unique(data.fec$su))

#initial values generating function
inits <- function(){
  list(a = rnorm(nsu,0,1),
       mu.a = rnorm(1,0,1),
       #sd.a = rlnorm(1,log(0.2),0.1),
       sd.a = rexp(1,0.5),
       b = rnorm(nsu,1,1),
       mu.b = rnorm(1,1,1),
       #sd.b = rlnorm(1,log(0.5),0.1),
       sd.b = rexp(1,0.5),
       r = rgamma(1,1,1))
}

# build model
fec.model <- nimbleModel(fec.code,
                         constants = list(nsu = nsu,
                                          nobs = nrow(data.fec),
                                          su = data.fec$su),
                         inits=inits(),
                         data = data.fec %>% select(n.eggs, length.lc),
                         buildDerivs = TRUE)

fec.model$simulate()
fec.model$calculate()

dataNodes <- fec.model$getNodeNames(dataOnly = TRUE)
parentNodes <- fec.model$getParents(dataNodes, stochOnly = TRUE) #all of these should be added to monitor below to recreate other model variabsLes...
stnodes <- fec.model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
allvars <- fec.model$getVarNames(nodes = stnodes)
mvars <- allvars[!(grepl("lifted",allvars))]  

# calculate vars to id NAs
for(i in 1:length(mvars)){
  print(paste0(mvars[i]," ",fec.model$calculate(mvars[i]) ))
}

# configure hmc
fec.confhmc <- configureHMC(fec.model,
                            monitors = c("a","b","mu.a","mu.b","sd.a","sd.b","p","r","mu"),
                            enableWAIC = TRUE)

# build mcmc (use buidlHMC() when not using configureHMC())
fec.hmc <- buildMCMC(fec.confhmc)

# compile model
fec.c <- compileNimble(fec.model)

# compile mcmc  and specify the project model
fec.hmcc <- compileNimble(fec.hmc)

# hmc samples
fec.samples <- runMCMC(fec.hmcc, niter = 5000, nburnin = 3000, nchains = 2, WAIC=TRUE, samplesAsCodaMCMC = TRUE)

# Save model and samples
saveRDS(fec.samples, file = paste0(home,"/data/fec_samples_",Sys.Date(),".RData"))

# calculate summary and save those for temporal biphasic model fec predictions
# fec.dpars <- mcmc.list(lapply(fec.samples$samples, as.mcmc)) %>% 
#   spread_draws(b[su],a[su]) %>%
#   ungroup() %>%
#   summarise(mufb = mean(b),  sdfb = sd(b), 
#             ma = mean(a), va  = var(a),
#             shfa  = ma^2 / va, rafa   = ma / va,
#             .by = su)
# 
# saveRDS(fec.dpars, file = paste0(home,"/data/fec_dpars.RData"))


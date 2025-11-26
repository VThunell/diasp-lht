# Load libraries, install if needed
pkgs <- c("tidyverse", "devtools","nimble")

if(length(setdiff(pkgs,rownames(installed.packages()))) > 0){
  install.packages(setdiff(pkgs, rownames(installed.packages())), dependencies = T)
}

invisible(lapply(pkgs, library, character.only = T))

## Biphasic 7 - rw 
biph7.code <- nimbleCode({
  
  # likelihood
  for(i in 1:nobs){
    length_mm[i] ~ dnorm(l.mu[growth.year[i], ind.id[i]], sd = sig.l)
  }
  
  # Estimate length at age
  for(i in 1:nind){
    
    #length prediction 0+
    l.mu[1,i] <- lb.mu[spat.unit[i]] +
      (1-step(1-(age[i]+0.0001)))*exp(g[spat.unit[i], hatch.year.f[i]]) + step(1-(age[i]+0.0001))*exp(g[spat.unit[i], hatch.year.f[i]])*doy.dec[i]
    
    #length prediction 1+ and older
    for(a in 2:(growth.year[i]+1)){
      l.mu[a,i] <- l.mu[(a-1),i] + 
        # (1-step(a-(age[i]+0.0001)))*(step(smo.age[i] - (a-0.999))*exp(g[spat.unit[i], (hatch.year.f[i]+a-1)]) + 
        #                                (1-step(smo.age[i] - (a-0.999)))*((exp(par[spat.unit[i],(hatch.year.f[i]+(a-1)),1]) - 
        #                                                                     l.mu[(a-1),i])*(1 - exp(-exp(par[spat.unit[i],(hatch.year.f[i]+(a-1)),2]))))) +
        # 
        # step(a-(age[i]+0.0001))*(step(smo.age[i] - (a-0.999))*exp(g[spat.unit[i], (hatch.year.f[i]+a-1)])*doy.dec[i] + 
        #                            (1-step(smo.age[i] - (a-0.999)))*((exp(par[spat.unit[i],(hatch.year.f[i]+(a-1)),1]) - 
        #                                                                 l.mu[(a-1),i])*(1 - exp(-exp(par[spat.unit[i],(hatch.year.f[i]+(a-1)),2])*doy.dec[i]))))
        (1-step(a-(age[i]+0.0001)))*(step(smo.age[i] - (a-0.999))*exp(g[spat.unit[i], (hatch.year.f[i]+a-1)]) + 
                                       (1-step(smo.age[i] - (a-0.999)))*((exp(par[(hatch.year.f[i]+(a-1)),spat.unit[i],1]) - 
                                                                            l.mu[(a-1),i])*(1 - exp(-exp(par[(hatch.year.f[i]+(a-1)),spat.unit[i],2]))))) +
        
        step(a-(age[i]+0.0001))*(step(smo.age[i] - (a-0.999))*exp(g[spat.unit[i], (hatch.year.f[i]+a-1)])*doy.dec[i] + 
                                   (1-step(smo.age[i] - (a-0.999)))*((exp(par[(hatch.year.f[i]+(a-1)),spat.unit[i],1]) - 
                                                                        l.mu[(a-1),i])*(1 - exp(-exp(par[(hatch.year.f[i]+(a-1)),spat.unit[i],2])*doy.dec[i]))))
    }
  }
  
  #par[gy,su,p]
  par[1:nyear,1:nsu,1] <- log_p.year[1:nyear,1:nsu]
  par[1:nyear,1:nsu,2] <- log_p.year[1:nyear,(nsu+1):npars]

  # estimate smo.age
  for(l in 1:nind){
    smo.age[l] ~ dcat(smp[1:nK,spat.unit[l]])
  }
  
  # Priors
  sig.l ~ dexp(1/100)
  
  # LKJ prior on correlation matrix, see NIMBLE manual p45.
  Ustar[1:npars,1:npars] ~ dlkj_corr_cholesky(1.3, npars) # eta = 1.3
  U[1:npars,1:npars] <- uppertri_mult_diag(Ustar[1:npars, 1:npars], sig.par[1:npars])
  
  # first rw row, npars = nsu*number of vb-pars
  for(k in 1:npars){
     log_p.year[1,k] ~ dnorm(mu.par[k],10)
  }
  
  for(j in minyear:maxyear){
    # Multivariate random walk in the log scale
    log_p.year[(j+1), 1:npars] ~ dmnorm(log_p.year[j,1:npars], cholesky = U[1:npars, 1:npars], prec_param = 0)
  }
  # priors for smolt age, lb.mu, g
  for(i in 1:nsu){
    smp[1:nK,i] ~ ddirch(alpha[1:nK])
    #lb.mu[i] ~ dnorm(17, sd = 5) # ~15 to 20 mm, https://doi.org/10.1111/j.1095-8649.2009.02497.x sd is guesstimate
    lb.mu[i] ~ dunif(12, 22) # ~15 to 20 mm, https://doi.org/10.1111/j.1095-8649.2009.02497.x sd is guesstimate
    for(j in minyear:maxyear){
      # g
      g[i,j] ~ dnorm(mean = g.lmean, sd = g.lsd) 
    }
  }
  
  g.lsd <- sqrt(log(1 + (10^2) / (50^2))) # sqrt of variance (= sd) of the lognormal distr. 50 and 10 is arbitrary atm.
  g.lmean <- log(50) - 0.5 * g.lsd^2  # mean of the lognorm distr.
 
  for(i in 1:nsu){
    mu.par[i] ~ dnorm(mean = l.lmean, sd = 0.10)
    sig.par[i] ~ dlnorm(log(l.lsd),10)  
  }
  for(i in (nsu+1):npars){
    mu.par[i] ~ dnorm(mean = k.lmean, sd = 0.10)
    sig.par[i] ~ dlnorm(log(k.lsd),10)  
  }
 
  # K & l values from fb  
  l.lsd <- sqrt(log(1 + (278^2) / (1378^2)))
  l.lmean <- log(1378) - 0.5 * sig.par[1]^2 # sd of the mean based on sd:s of LKJ
  k.lsd <- sqrt(log(1 + (0.28^2) / (0.43^2)))
  k.lmean <- log(0.43) - 0.5 * sig.par[2]^2  
  
})

# Function creating the Cholesky of the covar. matrix (p45 Nimble manual)
uppertri_mult_diag <- nimbleFunction(
  run = function(mat = double(2), vec = double(1)) {
    returnType(double(2))
    p <- length(vec)
    out <- matrix(nrow = p, ncol = p, init = FALSE)
    for(k in 1:p)
      out[ , k] <- mat[ , k] * vec[k]
    return(out)
    # turn off buildDerivs for the i index
  }, buildDerivs = list(run = list(ignore = c('k'))))

data <- sallaa %>%
  mutate(yday.may1 = if_else(is.na(date), 1, yday(date) - yday("2025-05-01")+1), 
         doy.dec = if_else(yday.may1 > 0, yday.may1, yday.may1+365)/365,
         tot.age.dec = if_else(yday.may1 > 0, tot.age+doy.dec, tot.age+doy.dec-1),
         smo.age = if_else(age.type == "both", juv.age, NA),
         sea.age = if_else(is.na(sea.age), 0, sea.age),
         growth.year = trunc(tot.age.dec)+1, # + 1 to index age 0 individuals
         hatch.year = year-tot.age,
         hatch.year.f = as.numeric(factor(hatch.year, levels = sort(unique(hatch.year)), labels = seq_along(sort(unique(hatch.year))))),
         year.f = as.integer(factor(year)),
         spat.unit = as.integer(factor(spat.unit)),
         ind.id = row_number())

# get max year by, needed as hatch.year is integer and the last year is not (simpler solution?).
maxyear <- max(data$hatch.year.f) + data %>%
  filter(year == max(year)) %>%
  filter(tot.age == max(tot.age)) %>%
  distinct(tot.age) %>% pull(tot.age)
minyear <- min(data$hatch.year.f)

nind <- nrow(data) # n individuals
nsu <- length(unique(data$spat.unit))
npars <- 2*nsu # n parameters in the LKJ dmnorm * number of spatial units 
nages = length(unique(data$growth.year))+1
nyear = length(minyear:maxyear)
nK = max(data$smo.age %>% na.omit())

inits <- function(){
  list(#l.mu = matrix(abs(rnorm(nages * nind,500,25)), nrow = nages, ncol = nind), # this is not correct, only inits of =<tot.age
    sig.l = abs(rnorm(1,250,25)),
    sig.par = c(rlnorm(nsu,0,1),
                rlnorm(nsu,0,1)),
    mu.par = c(rnorm(nsu,7,0.2),
               rnorm(nsu,-1,0.6)),
    lb.mu = abs(runif(nsu,,5)),
    g = matrix(rnorm(nsu * length(minyear:maxyear),4,.2), nrow = nsu, ncol = length(minyear:maxyear))
    # smo.age = round(runif(nind,0,6))
    #Ustar, #U, #par
  )
}

# build model
biph7.model <- nimbleModel(biph7.code,
                     constants = list(npars = npars,
                                      nobs = nrow(data),
                                      nind = nind, # change when >1 sample per ind.id
                                      nK = nK,
                                      alpha = c(1, 2, 3, 2, 1, 1)[1:nK],
                                      nsu = nsu,
                                      age = data$tot.age.dec,
                                      growth.year = data$growth.year,
                                      doy.dec = data$doy.dec,
                                      ind.id = data$ind.id,
                                      hatch.year.f = data$hatch.year.f,
                                      minyear = minyear,
                                      maxyear = maxyear,
                                      nyear = nyear,
                                      spat.unit = data$spat.unit),
                     inits = inits(),
                     data = data %>% select(length_mm))

biph7.model$simulate() 
biph7.model$calculate()    

# identify nodes to sample 
dataNodes <- biph7.model$getNodeNames(dataOnly = TRUE)
parentNodes <- biph7.model$getParents(dataNodes, stochOnly = TRUE) #all of these should be added to monitor below to recreate other model variables...
stnodes <- biph7.model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
allvars <- biph7.model$getVarNames(nodes = stnodes)
mvars <- allvars[!(grepl("lifted",allvars))]  
unique(gsub("[0-9]", "", parentNodes))
unique(gsub("[0-9]", "", stnodes))
unique(gsub("[0-9]", "", allvars))
# calculate vars
for(i in 1:length(mvars)){
  print(paste0(mvars[i]," ",biph7.model$calculate(mvars[i]) ))
}

biph7.model$log_p.year
# compile model
biph7.c <- compileNimble(biph7.model)

# configure and build mcmc
#biph7.confmcmc <- configureMCMC(biph7.c, monitors = c("sig.l", "Ustar","lb.mu","smp","g","sig.par","smo.age","mu.par","par", "l.mu", "log_p.year"),
                               # useConjugacy = FALSE)
biph7.confmcmc <- configureMCMC(biph7.c, monitors = c("sig.l","lb.mu","g","smo.age","log_p.year"),
                                useConjugacy = FALSE)

biph7.mcmc <- buildMCMC(biph7.confmcmc, project = biph7.model)

# compile mcmc
biph7.mcmcc <- compileNimble(biph7.mcmc, project = biph7.model)

# MCMC samples
#biph7.samples <- runMCMC(biph7.mcmcc, niter = 150000, nburnin = 130000, thin=10, WAIC=TRUE) 
#gc()

# MCMC Samples
#biph7.mcmcc$run(10000, time = TRUE, nburnin=5000)
biph7.mcmcc$run(100000, time = TRUE, nburnin=80000,thin = 10)

biph7.samples <- as.matrix(biph7.mcmcc$mvSamples)

# Save samples
saveRDS(biph7.samples, file = paste0(home,"/data/biph7_samples.RData"))


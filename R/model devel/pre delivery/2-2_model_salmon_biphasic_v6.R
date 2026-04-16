# Load libraries, install if needed
pkgs <- c("tidyverse", "devtools","nimble")

if(length(setdiff(pkgs,rownames(installed.packages()))) > 0){
  install.packages(setdiff(pkgs, rownames(installed.packages())), dependencies = T)
}

invisible(lapply(pkgs, library, character.only = T))

## Biphasic 6
biph6.code <- nimbleCode({
  
  # likelihood
  for(i in 1:nobs){
    length_mm[i] ~ dnorm(l.mu[growth.year[i], ind.id[i]], sd = sig.l)
  }
  
  # Estimate length at age
  for(i in 1:n.ind){
    
    #length prediction 0+
    l.mu[1,i] <- lb.mu[spat.unit[i]] +  
      (1-step(1-(age[i]+0.0001)))*exp(g[spat.unit[i], hatch.year.f[i]]) + step(1-(age[i]+0.0001))*exp(g[spat.unit[i], hatch.year.f[i]])*doy.dec[i]
    
    #length prediction 1+ and older
    for(a in 2:(growth.year[i]+1)){
      l.mu[a,i] <- l.mu[(a-1),i] + 
        (1-step(a-(age[i]+0.0001)))*(step(smo.age[i] - (a-0.999))*exp(g[spat.unit[i], (hatch.year.f[i]+a-1)]) + (1-step(smo.age[i] - (a-0.999)))*((exp(par[spat.unit[i],(hatch.year.f[i]+(a-1)),1]) - l.mu[(a-1),i])*(1 - exp(-exp(par[spat.unit[i],(hatch.year.f[i]+(a-1)),2]))))) +
        
        step(a-(age[i]+0.0001))*(step(smo.age[i] - (a-0.999))*exp(g[spat.unit[i], (hatch.year.f[i]+a-1)])*doy.dec[i] + (1-step(smo.age[i] - (a-0.999)))*((exp(par[spat.unit[i],(hatch.year.f[i]+(a-1)),1]) - l.mu[(a-1),i])*(1 - exp(-exp(par[spat.unit[i],(hatch.year.f[i]+(a-1)),2])*doy.dec[i]))))
    }
  }
  
  # estimate smo.age
  for(l in 1:n.ind){
    smo.age[l] ~ dcat(smp[1:nK,spat.unit[l]])
  }
  
  # Priors
  sig.l ~ dexp(1/100)
  
  # LKJ prior on correlation matrix, see NIMBLE manual p45.
  Ustar[1:npars,1:npars] ~ dlkj_corr_cholesky(1.3, npars) # eta = 1.3
  U[1:npars,1:npars] <- uppertri_mult_diag(Ustar[1:npars, 1:npars], sig.par[1:npars])
  
  # priors for smolt age and juvenile growth rate
  for(i in 1:nsu){
    smp[1:nK,i] ~ ddirch(alpha[1:nK])  
    lb.mu[i] ~ dnorm(17, sd = 5) # ~15 to 20 mm, https://doi.org/10.1111/j.1095-8649.2009.02497.x sd is guesstimate
    for(j in minyear:maxyear){
      par[i,j,1:npars] ~ dmnorm(mu.par[1:npars], cholesky = U[1:npars, 1:npars], prec_param = 0)
      g[i,j] ~ dnorm(mean = g.lmean, sd = g.lsd) 
    }
  }
  
  g.lsd <- sqrt(log(1 + (10^2) / (50^2))) # sqrt of variance (= sd) of the lognormal distr. 50 and 10 is arbitrary atm.
  g.lmean <- log(50) - 0.5 * g.lsd^2  # mean of the lognorm distr.
  
  mu.par[1] ~ dnorm(mean = l.lmean, sd = 0.10)
  mu.par[2] ~ dnorm(mean = k.lmean, sd = 0.10)
  
  sig.par[1] ~ dlnorm(log(l.lsd),10)  
  sig.par[2] ~ dlnorm(log(k.lsd),10)  
  
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

n.ind <- nrow(data) # n individuals
npars <- 2 # n parameters in the LKJ dmnorm 
nsu <- length(unique(data$spat.unit))
n.ages = length(unique(data$growth.year))+1
nK = max(data$smo.age %>% na.omit())
# get max year by, needed as hatch.year is integer and the last year is not (simpler solution?).
maxyear <- max(data$hatch.year.f) + data %>%
  filter(year == max(year)) %>%
  filter(tot.age == max(tot.age)) %>%
  distinct(tot.age) %>% pull(tot.age)
minyear = min(data$hatch.year.f)

inits <- function(){
  list(#l.mu = matrix(abs(rnorm(n.ages * n.ind,500,25)), nrow = n.ages, ncol = n.ind), # this is not correct, only inits of =<tot.age
    sig.l = abs(rnorm(1,250,25)),
    sig.par = c(rlnorm(1,0,1),
                rlnorm(1,0,1)),
    mu.par = c(rnorm(1,7,0.2),
               rnorm(1,-1,0.6)),
    lb.mu = abs(rnorm(nsu,20,5)),
    g = matrix(rnorm(nsu * length(minyear:maxyear),4,.2), nrow = nsu, ncol = length(minyear:maxyear))
    # smo.age = round(runif(n.ind,0,6))
    #Ustar, #U, #par
  )
}

# build model
biph6.model <- nimbleModel(biph6.code,
                     constants = list(npars = npars,
                                              nobs = nrow(data),
                                              n.ind = n.ind, # change when >1 sample per ind.id
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
                                              spat.unit = data$spat.unit),
                     inits = inits(),
                     data = data %>% select(length_mm))

biph6.model$simulate() 
biph6.model$calculate()    

# identify nodes to sample 
dataNodes <- biph6.model$getNodeNames(dataOnly = TRUE)
parentNodes <- biph6.model$getParents(dataNodes, stochOnly = TRUE) #all of these should be added to monitor below to recreate other model variables...
stnodes <- biph6.model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
allvars <- biph6.model$getVarNames(nodes = stnodes)
mvars <- allvars[!(grepl("lifted",allvars))]  

# calculate vars
for(i in 1:length(mvars)){
  print(paste0(mvars[i]," ",biph6.model$calculate(mvars[i]) ))
}

# compile model
biph6.c <- compileNimble(biph6.model)

# configure and build mcmc
biph6.confmcmc <- configureMCMC(biph6.c, monitors = c("sig.l", "Ustar","lb.mu","smp","g","sig.par","smo.age","mu.par","par","l.mu"), useConjugacy = FALSE)
# remove non parent nodes to reduce size
# use model object to simulate from the post.

biph6.mcmc <- buildMCMC(biph6.confmcmc, project = biph6.model)

# compile mcmc
biph6.mcmcc <- compileNimble(biph6.mcmc, project = biph6.model)

# MCMC Samples
biph6.mcmcc$run(100000, time = TRUE, nburnin=80000,thin = 10)

biph6.samples2 <- as.matrix(biph6.mcmcc$mvSamples)

# Save samples
#saveRDS(biph6.samples, file = paste0(home,"/data/biph6_samples1.RData"))


# Biphasic model 8:
# adds a covariate effect of sex on k and Linf
# fixes maxyear to not be too large
# add an rw prior to g
# estimates fecundity ?

# Load libraries, install if needed
pkgs <- c("tidyverse", "nimble")

if(length(setdiff(pkgs,rownames(installed.packages()))) > 0){
  install.packages(setdiff(pkgs, rownames(installed.packages())), dependencies = T)
}

invisible(lapply(pkgs, library, character.only = T))

## Biphasic 9
biph9.code <- nimbleCode({
  
  # likelihood
  for(i in 1:nind){
    length_mm[i] ~ dnorm(l.mu[growth.year[i], ind.id[i]], sd = sig.l)
    
    #length prediction 0+
    l.mu[1,i] <- lb.mu[su[i]] +
      (1-step(1-(age[i]+0.001)))*(exp(g[su[i]]) + step(1-(age[i]+0.001))*exp(g[su[i]])*doy.dec[i])
    
    #NOT SEX length prediction 1+ and older
    for(a in 2:(growth.year[i]+1)){
      l.mu[a,i] <- l.mu[(a-1),i] +
        (1-step(a-(age[i]+0.001)))*(step(smo.age[i] - (a-0.999))*exp(g[su[i]]) +
                                      (1-step(smo.age[i] - (a-0.999)))*(((exp(par[(hatch.year.f[i]+(a-1)),su[i],1])) -
                                                                           l.mu[(a-1),i])*(1 - exp(-(exp(par[(hatch.year.f[i]+(a-1)),su[i],2])))))) +
        
        step(a-(age[i]+0.001))*(step(smo.age[i] - (a-0.999))*exp(g[su[i]])*doy.dec[i] +
                                  (1-step(smo.age[i] - (a-0.999)))*(((exp(par[(hatch.year.f[i]+(a-1)),su[i],1])) -
                                                                       l.mu[(a-1),i])*(1 - exp(-(exp(par[(hatch.year.f[i]+(a-1)),su[i],2])*doy.dec[i])))))
      
    }
  }
  
  
  #par[gy,su,p]
  par[1:nyear,1:nsu,1] <- log_p.year[1:nyear,1:nsu]
  par[1:nyear,1:nsu,2] <- log_p.year[1:nyear,13:npars]
  
  # estimate smo.age
  for(l in 1:nind){
    smo.age[l] ~ dcat(smp[1:nK,su[l]])
  }
  
  # Priors
  sig.l ~ dexp(1/350)
  
  # LKJ prior on correlation matrix for vb-pars
  Ustar[1:npars,1:npars] ~ dlkj_corr_cholesky(1.3, npars) # eta = 1.3
  U[1:npars,1:npars] <- uppertri_mult_diag(Ustar[1:npars, 1:npars], sig.par[1:npars])
  
  # LKJ prior on correlation matrix for g
  # Ugstar[1:nsu,1:nsu] ~ dlkj_corr_cholesky(1.3, nsu) # eta = 1.3
  # Ug[1:nsu,1:nsu] <- uppertri_mult_diag(Ugstar[1:nsu, 1:nsu], g.sig.par[1:nsu])
  
  # first year of vb-par rw, npars = nsu*number of vb-pars
  for(k in 1:npars){
    log_p.year[1,k] ~ dnorm(mu.par[k],10)
  }
  # rw for vb-par rw years > 1
  for(j in minyear:maxyear){
    log_p.year[(j+1), 1:npars] ~ dmnorm(log_p.year[j,1:npars], cholesky = U[1:npars, 1:npars], prec_param = 0)
  }
  
  # priors for smolt age, lb.mu, g
  for(i in 1:nsu){
    g[i] ~ dnorm(g.par[i], sd = 1) # rw first year
    smp[1:nK,i] ~ ddirch(alpha[1:nK])
    lb.mu[i] ~ dunif(12, 22) # ~15 to 20 mm, DOI: 10.1111/j.1095-8649.2009.02497.x  sd is guesstimate
  }
  
  # rw for g rw years > 1
  # for(j in minyear:maxyear){
  #   g[(j+1),1:nsu] ~ dmnorm(g[j,1:nsu], cholesky = Ug[1:nsu, 1:nsu], prec_param = 0) # rw for years >1
  # }
  
  for(i in 1:nsu){
    mu.par[i] ~ dnorm(l.lmean, sd = 1)
    sig.par[i] ~ dlnorm(log(l.lsd), sdlog = 0.1)  
    g.par[i] ~ dnorm(g.lmean, sd = 1)
    g.sig.par[i] ~ dlnorm(log(g.lsd),sdlog = 0.1)  
  }
  for(i in (nsu+1):npars){
    mu.par[i] ~ dnorm(mean = k.lmean, sd = 0.5)
    sig.par[i] ~ dlnorm(log(k.lsd), sdlog = 0.1)
  }
  
  # g, k & linf values (g and linf from fb)
  g.lsd <- sqrt(log(1 + (5^2) / (50^2))) # sqrt of variance (= sd) of the lognormal distr. 50 and 10 is arbitrary atm.
  g.lmean <- log(50) - 0.5 * log(g.lsd)^2  # mean of the lognorm distr.
  l.lsd <- sqrt(log(1 + (278^2) / (1378^2)))
  l.lmean <- log(1378) - 0.5 * log(l.lsd)^2 # sd of the mean based on sd:s of LKJ
  k.lsd <- sqrt(log(1 + (0.28^2) / (0.43^2)))
  k.lmean <- log(0.43) - 0.5 * log(k.lsd)^2  
  
  
  
})

#Function creating the Cholesky of the covar. matrix (p45 Nimble manual)
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

set.seed(1)
data.biph <- sallaa %>%
  slice_sample(n = 50000) %>%
  mutate(yday.may1 = if_else(is.na(date), 1, yday(date) - yday("2025-05-01")+1), 
         doy.dec = if_else(yday.may1 > 0, yday.may1, yday.may1+365)/365,
         tot.age.dec = if_else(yday.may1 > 0, tot.age+doy.dec, tot.age+doy.dec-1),
         smo.age = if_else(age.type == "both", juv.age, NA),
         sea.age = if_else(is.na(sea.age), 0, sea.age),
         growth.year = trunc(tot.age.dec)+1, # + 1 to index age 0 individuals
         hatch.year = year-tot.age,
         hatch.year.f = as.numeric(factor(hatch.year, levels = sort(unique(hatch.year)), labels = seq_along(sort(unique(hatch.year))))),
         year.f = as.integer(factor(year)),
         sex = as.integer(factor(sex)),
         su = as.integer(factor(spat.unit))) %>%
  filter(!(is.na(smo.age))) %>%
  mutate(ind.id = row_number())

nK = max(data.biph$smo.age %>% na.omit())
nsu = length(unique(data.biph$su))
maxyear = max(data.biph$hatch.year.f + data.biph$growth.year)
minyear = min(data.biph$hatch.year.f)

consts = list(#nobs = nrow(data.biph),
  nind = nrow(data.biph), # n individuals change when >1 sample per ind.id
  nK = nK,
  #eta = 2,
  alpha = c(1, 2, 3, 2, 1, 1)[1:nK],
  nsu = nsu,
  npars = 2*nsu, # n parameters in the LKJ dmnorm * number of spatial units 
  age = data.biph$tot.age.dec,
  growth.year = data.biph$growth.year,
  doy.dec = data.biph$doy.dec,
  ind.id = data.biph$ind.id,
  hatch.year.f = data.biph$hatch.year.f,
  #shfa = fec.dpars$shfa,
  #rafa = fec.dpars$rafa,
  #mufb = fec.dpars$mufb,
  #sdfb = fec.dpars$sdfb,
  #sex = data.biph$sex,
  maxyear = maxyear,
  minyear = minyear,
  nyear = length(minyear:maxyear),
  su = data.biph$su)

#nages = length(unique(data.biph$growth.year))+1

inits <- function(){
  list(sig.l = rexp(1,1/150),
       lb.mu = runif(consts$nsu,13,19),
       Ustar  = diag(consts$npars),
       Ugstar = diag(consts$nsu),
       sig.par = c(rlnorm(consts$nsu,log(0.2),sdlog = 0.1),
                   rlnorm(consts$nsu,log(0.6),sdlog = 0.1)),
       mu.par = c(rnorm(consts$nsu,7,0.2),
                  rnorm(consts$nsu,-1,0.6)),
       g.sig.par = rlnorm(consts$nsu,0,0.1),
       g.par = rnorm(consts$nsu,4,0.2)
       #corZ = matrix(rnorm((consts$npars-1)*(consts$npars-1), 0, 0.1),
      #               nrow = consts$npars-1, ncol = consts$npars-1),
      # corY = runif((consts$npars-1), 0, 0.1)
       #bL = rnorm(consts$nsu, mean = 0, sd = .1), 
       #bk = rnorm(consts$nsu, mean = 0, sd = .1)
  )
}


# build model
t <- Sys.time()
biph9.model <- nimbleModel(biph9.code,
                           constants = consts,
                           inits = inits(),
                           data = data.biph %>% select(length_mm))

biph9.model$initializeInfo()
biph9.model$simulate()
biph9.model$calculate()
Sys.time() - t # with nimble LKJ and 50K ind, this takes 4.15489 min

# identify nodes to sample 
dataNodes <- biph9.model$getNodeNames(dataOnly = TRUE)
deterNodes <- biph9.model$getNodeNames(determOnly = TRUE)
parentNodes <- biph9.model$getParents(dataNodes, stochOnly = TRUE) #all of these should be added to monitor below to recreate other model variables...
#biph9.model$getVarNames(nodes = parentNodes)
stnodes <- biph9.model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
allvars <- biph9.model$getVarNames(nodes = stnodes)
mvars <- allvars[!(grepl("lifted",allvars))]  

# calculate vars
for(i in 1:length(mvars)){
  print(paste0(mvars[i]," ",biph9.model$calculate(mvars[i]) ))
}

# compile model
biph9.c <- compileNimble(biph9.model)

# configure and build mcmc
#biph9.confmcmc <- configureMCMC(biph9.c, monitors = c("sig.l", "Ustar","lb.mu","smp","g","sig.par","smo.age","mu.par","par", "l.mu", "log_p.year"),
# useConjugacy = FALSE)
biph9.confmcmc <- configureMCMC(biph9.c, monitors = c("sig.l","Ustar","lb.mu","g","log_p.year","smo.age"), enableWAIC = TRUE,#,"bL","bk"),
                                useConjugacy = FALSE)

t <- Sys.time()
biph9.mcmc <- buildMCMC(biph9.confmcmc, project = biph9.model)
Sys.time() - t

# compile mcmc
t <- Sys.time()
biph9.mcmcc <- compileNimble(biph9.mcmc, project = biph9.model, resetFunctions = TRUE)
Sys.time() - t

# MCMC samples

#gc()

# MCMC Samples
t <- Sys.time()
biph9.samples_bLKJ <- runMCMC(biph9.mcmcc, niter = 10000, nburnin = 7000, 
                              thin = 2, nchains = 2, samplesAsCodaMCMC = TRUE, WAIC = TRUE)
Sys.time() - t #Time difference of 2.289687 hours
biph9.samples_bLKJ$WAIC #397769.2
node.sub <- grep("log_p.year\\[", colnames(biph9.samples_bLKJ$samples[[1]]), value = TRUE)
biph9.samples_bLKJ$samples[, node.sub[], drop = FALSE] %>% gelman.diag()

gc()

#saveRDS(biph9.samples, file = paste0(home,"/data/biph9_samples_1025.RData"))
#biph9.samples <- as.matrix(biph9.mcmcc$mvSamples)

# Save samples
#saveRDS(biph9.samples, file = paste0(home,"/data/biph9_samples_20000_1017.RData"))



#### becky LKJ

## Biphasic 9
biph9.code <- nimbleCode({
  
  # likelihood
  for(i in 1:nind){
    length_mm[i] ~ dnorm(l.mu[growth.year[i], ind.id[i]], sd = sig.l)
    
    #length prediction 0+
    l.mu[1,i] <- lb.mu[su[i]] +
      (1-step(1-(age[i]+0.001)))*(exp(g[su[i]]) + step(1-(age[i]+0.001))*exp(g[su[i]])*doy.dec[i])
    
    #NOT SEX length prediction 1+ and older
    for(a in 2:(growth.year[i]+1)){
      l.mu[a,i] <- l.mu[(a-1),i] +
        (1-step(a-(age[i]+0.001)))*(step(smo.age[i] - (a-0.999))*exp(g[su[i]]) +
                                      (1-step(smo.age[i] - (a-0.999)))*(((exp(par[(hatch.year.f[i]+(a-1)),su[i],1])) -
                                                                           l.mu[(a-1),i])*(1 - exp(-(exp(par[(hatch.year.f[i]+(a-1)),su[i],2])))))) +
        
        step(a-(age[i]+0.001))*(step(smo.age[i] - (a-0.999))*exp(g[su[i]])*doy.dec[i] +
                                  (1-step(smo.age[i] - (a-0.999)))*(((exp(par[(hatch.year.f[i]+(a-1)),su[i],1])) -
                                                                       l.mu[(a-1),i])*(1 - exp(-(exp(par[(hatch.year.f[i]+(a-1)),su[i],2])*doy.dec[i])))))
      
    }
  }
  
  
  #par[gy,su,p]
  par[1:nyear,1:nsu,1] <- log_p.year[1:nyear,1:nsu]
  par[1:nyear,1:nsu,2] <- log_p.year[1:nyear,13:npars]
  
  # estimate smo.age
  for(l in 1:nind){
    smo.age[l] ~ dcat(smp[1:nK,su[l]])
  }
  
  # Priors
  sig.l ~ dexp(1/350)
  
  # LKJ prior on correlation matrix for vb-pars
  # Ustar[1:npars,1:npars] ~ dlkj_corr_cholesky(1.3, npars) # eta = 1.3
  # U[1:npars,1:npars] <- uppertri_mult_diag(Ustar[1:npars, 1:npars], sig.par[1:npars])
  
  # LKJ prior on correlation matrix for g
  #Ugstar[1:nsu,1:nsu] ~ dlkj_corr_cholesky(1.3, nsu) # eta = 1.3
  #Ug[1:nsu,1:nsu] <- uppertri_mult_diag(Ugstar[1:nsu, 1:nsu], g.sig.par[1:nsu])
  
  # first year of vb-par rw, npars = nsu*number of vb-pars
  for(k in 1:npars){
    log_p.year[1,k] ~ dnorm(mu.par[k],10)
  }
  # rw for vb-par rw years > 1
  for(j in minyear:maxyear){
    #log_p.year[(j+1), 1:npars] ~ dmnorm(log_p.year[j,1:npars], cholesky = U[1:npars, 1:npars], prec_param = 0)
    log_p.year[(j+1), 1:npars] ~ dmnorm(log_p.year[j,1:npars], prec = tau_p[1:npars,1:npars])
  }
  
  tau_p[1:npars,1:npars] <- inverse(sigma_p[1:npars, 1:npars])
  
  for(k in 1:npars){
    for(j in 1:npars){
      sigma_p[k,j] <- Rnew[k,j] * sig.par[k] * sig.par[j]
    }
  }
  
  # priors for smolt age, lb.mu, g
  for(i in 1:nsu){
    g[i] ~ dnorm(g.par[i], sd = 1) # rw first year
    smp[1:nK,i] ~ ddirch(alpha[1:nK])
    lb.mu[i] ~ dunif(12, 22) # ~15 to 20 mm, DOI: 10.1111/j.1095-8649.2009.02497.x  sd is guesstimate
  }
  
  # # rw for g rw years > 1
  # for(j in minyear:maxyear){
  #   g[(j+1),1:nsu] ~ dmnorm(g[j,1:nsu], cholesky = Ug[1:nsu, 1:nsu], prec_param = 0) # rw for years >1
  # }
  
  for(i in 1:nsu){
    mu.par[i] ~ dnorm(l.lmean, sd = 1)
    sig.par[i] ~ dlnorm(log(l.lsd), sdlog = 0.1)  
    g.par[i] ~ dnorm(g.lmean, sd = 1)
    g.sig.par[i] ~ dlnorm(log(g.lsd),sdlog = 0.1)  
  }
  for(i in (nsu+1):npars){
    mu.par[i] ~ dnorm(mean = k.lmean, sd = 0.5)
    sig.par[i] ~ dlnorm(log(k.lsd), sdlog = 0.1)
  }
  
  # g, k & linf values (g and linf from fb)
  g.lsd <- sqrt(log(1 + (5^2) / (50^2))) # sqrt of variance (= sd) of the lognormal distr. 50 and 10 is arbitrary atm.
  g.lmean <- log(50) - 0.5 * log(g.lsd)^2  # mean of the lognorm distr.
  l.lsd <- sqrt(log(1 + (278^2) / (1378^2)))
  l.lmean <- log(1378) - 0.5 * log(l.lsd)^2 # sd of the mean based on sd:s of LKJ
  k.lsd <- sqrt(log(1 + (0.28^2) / (0.43^2)))
  k.lmean <- log(0.43) - 0.5 * log(k.lsd)^2  
  
  # fecundity fa and fb
  # for(i in 1:nsu){
  #   fa[i] ~ dgamma(shape = shfa[i], rate = rafa[i])
  #   fb[i] ~ dnorm(mean = mufb[i], sd = sdfb[i])
  # }
  
  #Prior for correlation matrix (LKJ prior)
  phi[1]  <- eta + (npars - 2)/2
  corY[1] ~ dbeta(phi[1], phi[1])
  r12   <- 2 * corY[1] - 1
  ##
  R[1,1]     <- 1
  R[1,2]     <- r12
  R[2,2]     <- sqrt(1 - r12^2)

  R[2:npars,1]   <- 0

  for (m in 2:(npars-1)) {
    ## Draw beta random variable
    phi[m] <- phi[(m-1)] - 0.5
    corY[m] ~ dbeta(m / 2, phi[m])
    ## Draw uniformly on a hypersphere
    for (jj in 1:m) {
      corZ[m, jj] ~ dnorm(0, 1)
    }
    scZ[m, 1:m] <- corZ[m, 1:m] / sqrt(inprod(corZ[m, 1:m], corZ[m, 1:m]))
    R[1:m,(m+1)] <- sqrt(corY[m]) * scZ[m,1:m]
    R[(m+1),(m+1)] <- sqrt(1 - corY[m])
    for(jk in (m+1):npars){
      R[jk,m] <- 0
    }
  }  #m

  Rnew[1:npars,1:npars] <- t(R[1:npars,1:npars]) %*% R[1:npars,1:npars]
  
})

#Function creating the Cholesky of the covar. matrix (p45 Nimble manual)
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


nK = max(data.biph$smo.age %>% na.omit())
nsu = length(unique(data.biph$su))
maxyear = max(data.biph$hatch.year.f + data.biph$growth.year)
minyear = min(data.biph$hatch.year.f)

consts = list(#nobs = nrow(data.biph),
  nind = nrow(data.biph), # n individuals change when >1 sample per ind.id
  nK = nK,
  eta = 2,
  alpha = c(1, 2, 3, 2, 1, 1)[1:nK],
  nsu = nsu,
  npars = 2*nsu, # n parameters in the LKJ dmnorm * number of spatial units 
  age = data.biph$tot.age.dec,
  growth.year = data.biph$growth.year,
  doy.dec = data.biph$doy.dec,
  ind.id = data.biph$ind.id,
  hatch.year.f = data.biph$hatch.year.f,
  #shfa = fec.dpars$shfa,
  #rafa = fec.dpars$rafa,
  #mufb = fec.dpars$mufb,
  #sdfb = fec.dpars$sdfb,
  #sex = data.biph$sex,
  maxyear = maxyear,
  minyear = minyear,
  nyear = length(minyear:maxyear),
  su = data.biph$su)

#nages = length(unique(data.biph$growth.year))+1

inits <- function(){
  list(sig.l = rexp(1,1/150),
       lb.mu = runif(consts$nsu,13,19),
       #Ustar  = diag(consts$npars),
       Ugstar = diag(consts$nsu),
       sig.par = c(rlnorm(consts$nsu,log(0.2),sdlog = 0.1),
                   rlnorm(consts$nsu,log(0.6),sdlog = 0.1)),
       mu.par = c(rnorm(consts$nsu,7,0.2),
                  rnorm(consts$nsu,-1,0.6)),
       g.sig.par = rlnorm(consts$nsu,0,0.1),
       g.par = rnorm(consts$nsu,4,0.2),
       corZ = matrix(rnorm((consts$npars-1)*(consts$npars-1), 0, 0.1),
                     nrow = consts$npars-1, ncol = consts$npars-1),
       corY = runif((consts$npars-1), 0, 0.1)
       #bL = rnorm(consts$nsu, mean = 0, sd = .1), 
       #bk = rnorm(consts$nsu, mean = 0, sd = .1)
  )
}


# build model
t <- Sys.time()
biph9.model_b <- nimbleModel(biph9.code,
                           constants = consts,
                           inits = inits(),
                           data = data.biph %>% select(length_mm))

biph9.model_b$initializeInfo()
biph9.model_b$simulate()
biph9.model_b$calculate()
Sys.time() - t # with nimble LKJ and 50K ind, this takes 4.15489 min
# with becky LKJ and 50K ind, this takes 7.609782 mins

# values(biph9.model, "sig.l")
# biph9.model$calculate("l.mu")
# biph9.model$calculate("sig.l[1888]")
# biph9.model$simulate("sig.l")
# which(is.na(biph9.model$l.mu))
# biph9.model$l.mu[1888]
# stochNodes <- biph9.model$getNodeNames(stochOnly = TRUE)
# biph9.model$getNodeNames(stochOnly = TRUE, includeData = FALSE, topOnly = TRUE)

# identify nodes to sample 
dataNodes <- biph9.model_b$getNodeNames(dataOnly = TRUE)
deterNodes <- biph9.model_b$getNodeNames(determOnly = TRUE)
parentNodes <- biph9.model_b$getParents(dataNodes, stochOnly = TRUE) #all of these should be added to monitor below to recreate other model variables...
#biph9.model$getVarNames(nodes = parentNodes)
stnodes <- biph9.model_b$getNodeNames(stochOnly = TRUE, includeData = FALSE)
allvars <- biph9.model_b$getVarNames(nodes = stnodes)
mvars <- allvars[!(grepl("lifted",allvars))]  

# calculate vars
for(i in 1:length(mvars)){
  print(paste0(mvars[i]," ",biph9.model_b$calculate(mvars[i]) ))
}

# compile model
biph9.c_b <- compileNimble(biph9.model_b)

# configure and build mcmc
#biph9.confmcmc <- configureMCMC(biph9.c, monitors = c("sig.l", "Ustar","lb.mu","smp","g","sig.par","smo.age","mu.par","par", "l.mu", "log_p.year"),
# useConjugacy = FALSE)
biph9.confmcmc_b <- configureMCMC(biph9.c_b, monitors = c("sig.l", "R","lb.mu","g","log_p.year","smo.age"), enableWAIC = TRUE,#,"bL","bk"),
                                useConjugacy = FALSE)

t <- Sys.time()
biph9.mcmc_b <- buildMCMC(biph9.confmcmc_b, project = biph9.model_b)
Sys.time() - t

# compile mcmc
t <- Sys.time()
biph9.mcmcc_b <- compileNimble(biph9.mcmc_b, project = biph9.model_b, resetFunctions = TRUE)
Sys.time() - t

# MCMC samples

#gc()

# MCMC Samples
t <- Sys.time()
#biph9.mcmcc$run(200,nburnin=150)
biph9.samples_b2LKJ <- runMCMC(biph9.mcmcc_b, niter = 10000, nburnin = 7000, 
                              thin = 2, nchains = 2, samplesAsCodaMCMC = TRUE, WAIC = TRUE)
#biph9.samples_nLKJ <- runMCMC(biph9.mcmcc, niter = 5000, nburnin = 1500, nchains = 2, samplesAsCodaMCMC = TRUE, WAIC = TRUE)
#biph9.samples <- runMCMC(biph9.mcmcc, niter = 2000, nburnin = 1500, thin=10)
#biph9.samples <- runMCMC(biph9.mcmcc, niter = 200000, nburnin = 185000, thin=10)#, WAIC=TRUE) 
Sys.time() - t
biph9.samples_b2LKJ$WAIC
node.sub <- grep("log_p.year\\[", colnames(biph9.samples_b2LKJ$samples[[1]]), value = TRUE)
biph9.samples_b2LKJ$samples[, node.sub[], drop = FALSE] %>% gelman.diag()


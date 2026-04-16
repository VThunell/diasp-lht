# Biphasic model 8:
# adds a covariate effect of sex on k and Linf ?
# fixes maxyear to not be too large
# add an rw prior to g
# estimates fecundity ?

# Load libraries, install if needed
# pkgs <- c("tidyverse", "devtools","nimble")
# 
# if(length(setdiff(pkgs,rownames(installed.packages()))) > 0){
#   install.packages(setdiff(pkgs, rownames(installed.packages())), dependencies = T)
# }
# 
# invisible(lapply(pkgs, library, character.only = T))

## Biphasic 8
biph8.code <- nimbleCode({
  
  # likelihood
  for(i in 1:nobs){
    length_mm[i] ~ dnorm(l.mu[growth.year[i], ind.id[i]], sd = sig.l)
  }
  
  # Estimate length at age
  for(i in 1:nind){
    
    #length prediction 0+
    l.mu[1,i] <- lb.mu[su[i]] +
      (1-step(1-(age[i]+0.001)))*(exp(g[hatch.year.f[i], su[i]]) + step(1-(age[i]+0.001))*exp(g[hatch.year.f[i], su[i]])*doy.dec[i])
    
    #NOT SEX length prediction 1+ and older
    for(a in 2:(growth.year[i]+1)){
      l.mu[a,i] <- l.mu[(a-1),i] +
        (1-step(a-(age[i]+0.001)))*(step(smo.age[i] - (a-0.999))*exp(g[(hatch.year.f[i]+(a-1)), su[i]]) +
                                      (1-step(smo.age[i] - (a-0.999)))*(((exp(par[(hatch.year.f[i]+(a-1)),su[i],1])) -
                                                                           l.mu[(a-1),i])*(1 - exp(-(exp(par[(hatch.year.f[i]+(a-1)),su[i],2])))))) +

        step(a-(age[i]+0.001))*(step(smo.age[i] - (a-0.999))*exp(g[(hatch.year.f[i]+(a-1)), su[i]])*doy.dec[i] +
                                  (1-step(smo.age[i] - (a-0.999)))*(((exp(par[(hatch.year.f[i]+(a-1)),su[i],1])) -
                                                                          l.mu[(a-1),i])*(1 - exp(-(exp(par[(hatch.year.f[i]+(a-1)),su[i],2])*doy.dec[i])))))
      # W SEX length prediction 1+ and older
      # for(a in 2:(growth.year[i]+1)){
      #   l.mu[a,i] <- l.mu[(a-1),i] +
      #     (1-step(a-(age[i]+0.001)))*(step(smo.age[i] - (a-0.999))*exp(g[(hatch.year.f[i]+(a-1)), su[i]]) +
      #                                   (1-step(smo.age[i] - (a-0.999)))*(((exp(par[(hatch.year.f[i]+(a-1)),su[i],1]) + bL[su[i]]*sex[i]) -
      #                                                                        l.mu[(a-1),i])*(1 - exp(-(exp(par[(hatch.year.f[i]+(a-1)),su[i],2]) + bk[su[i]]*sex[i]))))) +
      # 
      #     step(a-(age[i]+0.001))*(step(smo.age[i] - (a-0.999))*exp(g[(hatch.year.f[i]+(a-1)), su[i]])*doy.dec[i] +
      #                               (1-step(smo.age[i] - (a-0.999)))*(((exp(par[(hatch.year.f[i]+(a-1)),su[i],1]) + bL[su[i]]*sex[i]) -
      #                                                                    l.mu[(a-1),i])*(1 - exp(-((exp(par[(hatch.year.f[i]+(a-1)),su[i],2])) + bk[su[i]]*sex[i])*doy.dec[i]))))
      # 

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
  Ugstar[1:nsu,1:nsu] ~ dlkj_corr_cholesky(1.3, nsu) # eta = 1.3
  Ug[1:nsu,1:nsu] <- uppertri_mult_diag(Ugstar[1:nsu, 1:nsu], g.sig.par[1:nsu])
  
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
    g[1,i] ~ dnorm(g.par[i], sd = 1) # rw first year
    bL[i] ~ dnorm(mean = 0, sd = 1) 
    bk[i] ~ dnorm(mean = 0, sd = 1) 
    smp[1:nK,i] ~ ddirch(alpha[1:nK])
    #lb.mu[i] ~ dnorm(17, sd = 5) 
    lb.mu[i] ~ dunif(12, 22) # ~15 to 20 mm, https://doi.org/10.1111/j.1095-8649.2009.02497.x sd is guesstimate
  }
  
  # rw for g rw years > 1
  for(j in minyear:maxyear){
    g[(j+1),1:nsu] ~ dmnorm(g[j,1:nsu], cholesky = Ug[1:nsu, 1:nsu], prec_param = 0) # rw for years >1
  }
  
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

data.biph <- sallaa %>%
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
         su = as.integer(factor(spat.unit)),
         ind.id = row_number())

maxyear <- max(data.biph$hatch.year.f + data.biph$growth.year)
minyear <- min(data.biph$hatch.year.f)

nind <- nrow(data.biph) # n individuals
nsu <- length(unique(data.biph$su))
npars <- 2*nsu # n parameters in the LKJ dmnorm * number of spatial units 
nages = length(unique(data.biph$growth.year))+1
nyear = length(minyear:maxyear)
nK = max(data.biph$smo.age %>% na.omit())

inits <- function(){
  list(sig.l = rexp(1,1/250),
       lb.mu = runif(nsu,13,19),
       Ustar  = diag(npars),
       Ugstar = diag(nsu),
       sig.par = c(rlnorm(nsu,log(0.2),sdlog = 0.1),
                   rlnorm(nsu,log(0.6),sdlog = 0.1)),
       mu.par = c(rnorm(nsu,7,0.2),
                  rnorm(nsu,-1,0.6)),
       mu.par = c(rnorm(nsu,0,0.1),
                  rnorm(nsu,0,0.1)),
       g.sig.par = rlnorm(nsu,0,0.1),
       g.par = rnorm(nsu,4,0.2),
       bL = rnorm(nsu, mean = 0, sd = .1), 
       bk = rnorm(nsu, mean = 0, sd = .1)
      )
  }

inits <- function(){
  list(sig.l = runif(1,0.1,1),
       lb.mu = runif(nsu,13,19),
       Ustar  = diag(npars),
       Ugstar = diag(nsu),
       sig.par = c(rlnorm(nsu,-1,sdlog = 0.1),
                   rlnorm(nsu,-1,sdlog = 0.1)),
       mu.par = c(rnorm(nsu,0,0.1),
                  rnorm(nsu,0,0.1)),
       mu.par = c(rnorm(nsu,0,0.1),
                  rnorm(nsu,0,0.1)),
       g.sig.par = rlnorm(nsu,-1,0.1),
       g.par = rnorm(nsu,0,0.1),
       bL = rnorm(nsu, mean = 0, sd = .1), 
       bk = rnorm(nsu, mean = 0, sd = .1)
  )
}

# build model
biph8.model <- nimbleModel(biph8.code,
                     constants = list(npars = npars,
                                      nobs = nrow(data.biph),
                                      nind = nind, # change when >1 sample per ind.id
                                      nK = nK,
                                      alpha = c(1, 2, 3, 2, 1, 1)[1:nK],
                                      nsu = nsu,
                                      age = data.biph$tot.age.dec,
                                      growth.year = data.biph$growth.year,
                                      doy.dec = data.biph$doy.dec,
                                      ind.id = data.biph$ind.id,
                                      hatch.year.f = data.biph$hatch.year.f,
                                      sex = data.biph$sex,
                                      minyear = minyear,
                                      maxyear = maxyear,
                                      nyear = nyear,
                                      su = data.biph$su),
                     inits = inits(),
                     data = data.biph %>% select(length_mm))

biph8.model$initializeInfo()
biph8.model$simulate()
biph8.model$calculate()
# values(biph8.model, "sig.l")
# biph8.model$calculate("l.mu")
# biph8.model$calculate("sig.l[1888]")
# biph8.model$simulate("sig.l")
# which(is.na(biph8.model$l.mu))
# biph8.model$l.mu[1888]
# stochNodes <- biph8.model$getNodeNames(stochOnly = TRUE)
# biph8.model$getNodeNames(stochOnly = TRUE, includeData = FALSE, topOnly = TRUE)

# identify nodes to sample 
dataNodes <- biph8.model$getNodeNames(dataOnly = TRUE)
deterNodes <- biph8.model$getNodeNames(determOnly = TRUE)
parentNodes <- biph8.model$getParents(dataNodes, stochOnly = TRUE) #all of these should be added to monitor below to recreate other model variables...
stnodes <- biph8.model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
allvars <- biph8.model$getVarNames(nodes = stnodes)
mvars <- allvars[!(grepl("lifted",allvars))]  

# calculate vars
for(i in 1:length(mvars)){
  print(paste0(mvars[i]," ",biph8.model$calculate(mvars[i]) ))
}

# for(i in 1:length(deterNodes)){
#   if(is.na(biph8.model$calculate(deterNodes[i])))
#      {print(paste0(deterNodes[i]," ",biph8.model$calculate(deterNodes[i])))}
# }
# 
# # simulate vars
# for(i in 1:length(dataNodes)){
#   biph8.model$simulate(mvars[i])
#   print(paste0(mvars[i],"",values(biph8.model, mvars[i])))
# }
# for(i in 1:length(deterNodes)){
#   biph8.model$simulate(deterNodes[i])
#   if(sum(is.na(values(biph8.model, deterNodes[i])))>0)
#   {print(deterNodes[i])}
# }



# compile model
biph8.c <- compileNimble(biph8.model)

# configure and build mcmc
#biph8.confmcmc <- configureMCMC(biph8.c, monitors = c("sig.l", "Ustar","lb.mu","smp","g","sig.par","smo.age","mu.par","par", "l.mu", "log_p.year"),
                               # useConjugacy = FALSE)
biph8.confmcmc <- configureMCMC(biph8.c, monitors = c("sig.l", "Ustar","U","lb.mu","g","smo.age","log_p.year"),#,"bL","bk"),
                                useConjugacy = FALSE)

biph8.mcmc <- buildMCMC(biph8.confmcmc, project = biph8.model)

# compile mcmc
biph8.mcmcc <- compileNimble(biph8.mcmc, project = biph8.model)

# MCMC samples
#biph8.samples <- runMCMC(biph8.mcmcc, niter = 150000, nburnin = 130000, thin=10, WAIC=TRUE) 
#gc()

# MCMC Samples
t <- Sys.time()
biph8.mcmcc$run(200,nburnin=150)
Sys.time() - t

#the low logProb (less than -1e12) looks very probable (using dnorm on lb.mu):
#data[c(34984,35507,36333,36933,46327,37188,14572,14562,14654),]
#gc()
#biph8.mcmcc$run(10000, time = TRUE, nburnin=5000)
#biph8.mcmcc$run(300000, time = TRUE, nburnin=270000,thin = 10)

#biph8.samples <- as.matrix(biph8.mcmcc$mvSamples)

# Save samples
#saveRDS(biph8.samples, file = paste0(home,"/data/biph8_samples_20000_1017.RData"))


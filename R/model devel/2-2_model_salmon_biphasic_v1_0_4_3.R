# Biphasic model 1.0.4:
# this is 1_0_4 which worked but using new data and including smo_age data
# Load libraries, install if needed
pkgs <- c("tidyverse", "nimble")

if(length(setdiff(pkgs,rownames(installed.packages()))) > 0){
  install.packages(setdiff(pkgs, rownames(installed.packages())), dependencies = T)
}

invisible(lapply(pkgs, library, character.only = T))

## Biphasic code
biph1_0_4.code <- nimbleCode({
  
  l_sd_m <- sqrt(log(1 + 278^2 / 1378^2))
  l_sd_f <- sqrt(log(1 + 21^2 / 171^2))
  l_mean_m <- log(1378) - 0.5 * l_sd_m^2 # sd of the mean based on sd:s of LKJ
  l_mean_f <- log(171) - 0.5 * l_sd_f^2 # sd of the mean based on sd:s of LKJ
  k_sd_m <- sqrt(log(1 + 0.28^2 / 0.43^2))
  k_sd_f <- sqrt(log(1 + 0.2^2 / 0.61^2))
  k_mean_m <- log(0.43) - 0.5 * k_sd_m^2
  k_mean_f <- log(0.61) - 0.5 * k_sd_f^2 
  o_mean_m <- l_mean_m + k_mean_m # addition as they are on log scale
  o_mean_f <- k_mean_f + k_mean_f
  # lh_sd <- sqrt(log(1 + 3^2 / 17^2))
  # lh_mean <- log(17) - 0.5 * lh_sd^2
  # lh_sd <- sqrt(log(1 + 3^2 / 17^2)) # creates high lh_mu values no good using dlnorm
  # lh_mean <- log(17) - 0.5 * lh_sd^2
  
  # https://onlinelibrary.wiley.com/doi/full/10.1111/j.1095-8649.2012.03370.x # using only wild origin
  # mean(c(152,167,195,174,154,164,152,167,195,174,166,171,171,176,173,187,235,146,147,152))
  # sd(c(152,167,195,174,154,164,152,167,195,174,166,171,171,176,173,187,235,146,147,152))
  # for juvenile K, using mean for white spotted charr in https://doi.org/10.1111/j.1095-8649.2001.tb00220.x
  # lhmu / sd DOI: 10.1111/j.1095-8649.2009.02497.x  sd is guesstimate
  
  par_mu_m[1] ~ dnorm(o_mean_m, sd = 0.3)
  par_mu_f[1] ~ dnorm(o_mean_f, sd = 0.3)
  par_sig_m[1] ~ dlnorm(log(0.3), sdlog = 0.1)
  par_sig_f[1] ~ dlnorm(log(0.3), sdlog = 0.1)  
  par_mu_m[2] ~ dnorm(mean = k_mean_m, sd = 0.3)
  par_mu_f[2] ~ dnorm(mean = k_mean_f, sd = 0.3)
  par_sig_m[2] ~ dlnorm(log(k_sd_m), sdlog = 0.1)
  par_sig_f[2] ~ dlnorm(log(k_sd_f), sdlog = 0.1)

  #Prior for correlation matrix (LKJ prior)
  phi[1]  <- eta + (npars - 2)/2
  corY[1] ~ dbeta(phi[1], phi[1])
  r12   <- 2 * corY[1] - 1
  ##
  R[1,1]     <- 1
  R[1,2]     <- r12
  R[2,2]     <- sqrt(1 - r12^2)
  
  R[2,1]   <- 0
  
  Rnew[1:npars,1:npars] <- t(R[1:npars,1:npars]) %*% R[1:npars,1:npars]
  
  # freshwater
  phi_f[1]  <- eta + (npars - 2)/2
  corY_f[1] ~ dbeta(phi_f[1], phi_f[1])
  r12_f   <- 2 * corY_f[1] - 1
  ##
  R_f[1,1]     <- 1
  R_f[1,2]     <- r12_f
  R_f[2,2]     <- sqrt(1 - r12^2)
  
  R_f[2,1]   <- 0
  
  Rnew_f[1:npars,1:npars] <- t(R_f[1:npars,1:npars]) %*% R_f[1:npars,1:npars]
  
  for(k in 1:npars){
    for(j in 1:npars){
      sigma_p_m[k,j] <- Rnew[k,j] * par_sig_m[k] * par_sig_m[j]
      sigma_p_f[k,j] <- Rnew_f[k,j] * par_sig_f[k] * par_sig_f[j]
    }
  }
  
  tau_p_m[1:npars,1:npars] <- inverse(sigma_p_m[1:npars, 1:npars])
  tau_p_f[1:npars,1:npars] <- inverse(sigma_p_f[1:npars, 1:npars])
  
  # vB pars spatial MVN npars = nsu*2, number of vb-pars
  for(k in 1:nsu){
    par_m[k,1:npars] ~ dmnorm(par_mu_m[1:npars], prec = tau_p_m[1:npars,1:npars])
    par_f[k,1:npars] ~ dmnorm(par_mu_f[1:npars], prec = tau_p_f[1:npars,1:npars])
  }
  
  # priors for smolt age, lh.mu
  for(i in 1:nsu){
    smp[1:nC,i] ~ ddirch(alpha[1:nC])
    lh_mu[i] ~ dnorm(17, sd = 3) 
    #lh_mu[i] ~ dlnorm(lh_mean, sdlog = lh_sd) 
  }
  
  # estimate smo_age
  for(l in 1:nind){
    smo_age[l] ~ dcat(smp[1:nC,su[l]])
  }
  
  sig_l ~ dexp(1/350)
  
  # likelihood
  for(i in 1:nind){
    length_mm[i] ~ dnorm(l_mu[growth_year[i], ind_id[i]], sd = sig_l)
    
    l_mu[1,i] <- lh_mu[su[i]] +
      (1-step(1-(age[i]+0.001)))*((exp(par_f[su[i],1])/exp(par_f[su[i],2])) - lh_mu[su[i]])*(1 - exp(-exp(par_f[su[i],2]))) +
      step(1-(age[i]+0.001))*((exp(par_f[su[i],1])/exp(par_f[su[i],2])) - lh_mu[su[i]])*(1 - exp(-exp(par_f[su[i],2])))*doy_dec[i]
    
    for(a in 2:(growth_year[i]+1)){
      l_mu[a,i] <- l_mu[(a-1),i] +
        (1-step(a-(age[i]+0.001)))*(step(smo_age[i] - (a-0.999))*(( (exp(par_f[su[i],1])/exp(par_f[su[i],2])) - l_mu[(a-1),i])*(1 - exp(-exp(par_f[su[i],2])))) +
                                      (1-step(smo_age[i] - (a-0.999)))*(((exp(par_m[su[i],1])/exp(par_m[su[i],2])) - l_mu[(a-1),i])*(1 - exp(-(exp(par_m[su[i],2])))))) +
        
        step(a-(age[i]+0.001))* (step(smo_age[i] - (a-0.999))*(( (exp(par_f[su[i],1])/exp(par_f[su[i],2])) - l_mu[(a-1),i])*(1 - exp(-exp(par_f[su[i],2]))))*doy_dec[i] +
                                   (1-step(smo_age[i] - (a-0.999)))*(((exp(par_m[su[i],1])/exp(par_m[su[i],2])) - l_mu[(a-1),i])*(1 - exp(-(exp(par_m[su[i],2])))))*doy_dec[i]) 
      
    }
    # calculate Linf
    # Linf[su[i]] <- exp(par_m[su[i],1] - par_m[su[i],2])
  }
  
})

biph_step <- nimbleFunction(
  run = function(l_prev = double(0), L_f = double(0), L_m = double(0), k_f = double(0), k_m = double(0), doy_dec = double(0),
                 j = integer(0), age = double(0), smo_age = double(0)) {

    if(j > age) { # year of catch
      if(j <= smo_age) {
        return(l_prev + (L_f - l_prev)*(1 - exp(-k_f*doy_dec)))
      } else {
        return(l_prev + (L_m - l_prev)*(1 - exp(-k_m*doy_dec)))
      }
    } else { # full year
      if(j <= smo_age) {
        return(l_prev + (L_f - l_prev)*(1 - exp(-k_f)))
      } else {
        return(l_prev + (L_m - l_prev)*(1 - exp(-k_m)))
      }
    }
    returnType(double(0))
  }, buildDerivs = TRUE) # ignore indexes might fix som issues


nC = max(data.biph$smo_age %>% na.omit())
nsu = length(unique(data.biph$spat.unit))

consts = list(nind = nrow(data.biph),
              nC = nC,
              eta = 2,
              alpha = c(2, 3, 2, 1, 1, 1)[1:nC],
              nsu = nsu,
              npars = 2,
              age = data.biph$tot.age.dec,
              growth_year = data.biph$growth.year,
              doy_dec = data.biph$doy.dec,
              ind_id = data.biph$ind.id,
              su = data.biph$su)

initsmo_age = data.biph %>% 
  mutate(sa = if_else(is.na(smo_age), 2, NA)) %>% pull(sa)
inits <- function(){
  list(sig_l = rexp(1,1/150),
       lh_mu = runif(consts$nsu,10,20),
       par_sig_m = c(rlnorm(consts$npars,log(0.1),sdlog = 0.1)),
       par_sig_f = c(rlnorm(consts$npars,log(0.1),sdlog = 0.1)),
       par_mu_f = c(rnorm(1,5,0.2),
                    (rnorm(1,-1,0.2))),
       par_mu_f = c(rnorm(1,3,0.2),
                   (rnorm(1,-1,0.2))),
       smp = matrix(rep(1/nC,nC*nsu),nC,nsu),
       smo_age = initsmo_age#,
       #l_mu = matrix(runif(consts$nind*12,50,1000), nrow= 12, ncol = consts$nind)
       )
  }

# build model
t <- Sys.time()
biph1_0_4.model <- nimbleModel(biph1_0_4.code,
                             constants = consts,
                             inits = inits(),
                             data = data.biph %>% select(length_mm, smo_age),
                             buildDerivs = TRUE
                             #,calculate = FALSE
                           )


biph1_0_4.model$simulate()
biph1_0_4.model$calculate()
biph1_0_4.model$initializeInfo()
Sys.time() - t 
# identify nodes to sample 
dataNodes <- biph1_0_4.model$getNodeNames(dataOnly = TRUE)
deterNodes <- biph1_0_4.model$getNodeNames(determOnly = TRUE)
parentNodes <- biph1_0_4.model$getParents(dataNodes, stochOnly = TRUE) #all of these should be added to monitor below to recreate other model variables...
stnodes <- biph1_0_4.model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
allvars <- biph1_0_4.model$getVarNames(nodes = stnodes)
mvars <- allvars[!(grepl("lifted",allvars))]  

# calculate vars
for(i in 1:length(mvars)){
  print(paste0(mvars[i]," ",biph1_0_4.model$calculate(mvars[i]) ))
}

# compile model
biph1_0_4.c <- compileNimble(biph1_0_4.model)

# configure and build mcmc
#monits <- c(mvars, "l_mu")
monits <- c(mvars)#, "length_mm")

biph1_0_4.confmcmc <- configureHMC(biph1_0_4.c, monitors = monits, enableWAIC = TRUE)

t <- Sys.time()
biph1_0_4.mcmc <- buildMCMC(biph1_0_4.confmcmc, project = biph1_0_4.model)
Sys.time() - t

# compile mcmc
t <- Sys.time()
biph1_0_4.mcmcc <- compileNimble(biph1_0_4.mcmc, project = biph1_0_4.model)
Sys.time() - t

# MCMC Samples
t <- Sys.time()
biph1_0_4.samples <- runMCMC(biph1_0_4.mcmcc, niter = 2000, nburnin = 1000, nchains = 1, WAIC = TRUE)
#biph.samples <- runMCMC(biph1_0_4.mcmcc, niter = 100000, nburnin = 98500, nchains = 1, WAIC = TRUE)
Sys.time() - t

biph1_0_4.samples$WAIC #1524671

# Save samples
saveRDS(biph1_0_4.samples, file = paste0(home,"/data/biph1_0_4_samples_",Sys.Date(),".RData"))

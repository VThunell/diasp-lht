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
  
  # likelihood
  for(i in 1:nind){
    length_mm[i] ~ dnorm(l_mu[growth_year[i], ind_id[i]], sd = sig_l)
    
    l_mu[1,i] <- lb_mu[su[i]] +
      (1-step(1-(age[i]+0.001)))*(exp(par_f[su[i],1]) - lb_mu[su[i]])*(1 - exp(-exp(par_f[su[i],2]))) +
         step(1-(age[i]+0.001))*(exp(par_f[su[i],1]) - lb_mu[su[i]])*(1 - exp(-exp(par_f[su[i],2])))*doy_dec[i]
    
    for(a in 2:(growth_year[i]+1)){
      l_mu[a,i] <- l_mu[(a-1),i] +
        (1-step(a-(age[i]+0.001)))*(step(smo_age[i] - (a-0.999))*((exp(par_f[su[i],1]) - l_mu[(a-1),i])*(1 - exp(-exp(par_f[su[i],2])))) +
                                      (1-step(smo_age[i] - (a-0.999)))*(((exp(par_m[su[i],1])/exp(par_m[su[i],2])) - l_mu[(a-1),i])*(1 - exp(-(exp(par_m[su[i],2])))))) +
        
        step(a-(age[i]+0.001))* (step(smo_age[i] - (a-0.999))*((exp(par_f[su[i],1]) - l_mu[(a-1),i])*(1 - exp(-exp(par_f[su[i],2]))))*doy_dec[i] +
                                   (1-step(smo_age[i] - (a-0.999)))*(((exp(par_m[su[i],1])/exp(par_m[su[i],2])) - l_mu[(a-1),i])*(1 - exp(-(exp(par_m[su[i],2])))))*doy_dec[i]) 
    
    }
    # calculate Linf
    #Linf[su[i]] <- exp(par_m[su[i],1] - par_m[su[i],2])
  }
  
  # estimate smo_age
  for(l in 1:nind){
    smo_age[l] ~ dcat(smp[1:nK,su[l]])
  }
  
  # Priors
  sig_l ~ dexp(1/350)
  
  # vB pars spatial MVN npars = nsu*2, number of vb-pars
  for(k in 1:nsu){
    par_m[k,1:npars] ~ dmnorm(mu_par_m[1:npars], prec = tau_p_m[1:npars,1:npars])
    par_f[k,1:npars] ~ dmnorm(mu_par_f[1:npars], prec = tau_p_f[1:npars,1:npars])
  }
  
  tau_p_m[1:npars,1:npars] <- inverse(sigma_p_m[1:npars, 1:npars])
  tau_p_f[1:npars,1:npars] <- inverse(sigma_p_f[1:npars, 1:npars])
  
  for(k in 1:npars){
    for(j in 1:npars){
      sigma_p_m[k,j] <- Rnew[k,j] * sig_par_m[k] * sig_par_m[j]
      sigma_p_f[k,j] <- Rnew_f[k,j] * sig_par_f[k] * sig_par_f[j]
    }
  }
  
  # priors for smolt age, lb.mu, g
  for(i in 1:nsu){
    smp[1:nK,i] ~ ddirch(alpha[1:nK])
    lb_mu[i] ~ dnorm(17, sd = 2) 
  }
  
  mu_par_m[1] ~ dnorm(o_lmean_m, sd = 1)
  # mu_par_m[1] ~ dnorm(l_lmean_m, sd = 1)
  mu_par_f[1] ~ dnorm(l_lmean_f, sd = 1)
  sig_par_m[1] ~ dlnorm(log(l_lsd_f), sdlog = 0.1)  
  #sig_par_m[1] ~ dlnorm(log(l_lsd_m), sdlog = 0.1)  
  sig_par_f[1] ~ dlnorm(log(l_lsd_f), sdlog = 0.1)  
  mu_par_m[2] ~ dnorm(mean = k_lmean_m, sd = 0.5)
  mu_par_f[2] ~ dnorm(mean = k_lmean_f, sd = 0.5)
  sig_par_m[2] ~ dlnorm(log(k_lsd_m), sdlog = 0.1)
  sig_par_f[2] ~ dlnorm(log(k_lsd_f), sdlog = 0.1)
  
  # k & linf values (g and linf from fb)
  l_lsd_m <- sqrt(log(1 + (278^2) / (1378^2)))
  l_lsd_f <- sqrt(log(1 + (30^2) / (150^2)))
  l_lmean_m <- log(1378) - 0.5 * log(l_lsd_m)^2 # sd of the mean based on sd:s of LKJ
  l_lmean_f <- log(150) - 0.5 * log(l_lsd_f)^2 # sd of the mean based on sd:s of LKJ
  k_lsd_m <- sqrt(log(1 + (0.28^2) / (0.43^2)))
  k_lsd_f <- sqrt(log(1 + (0.1^2) / (0.3^2)))
  k_lmean_m <- log(0.43) - 0.5 * log(k_lsd_m)^2
  k_lmean_f <- log(0.3) - 0.5 * log(k_lsd_f)^2 
  o_lmean_m <- l_lmean_m + k_lmean_m
  
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

  
})

set.seed(1)
data.biph_m <- sallaa %>%
  drop_na(length_mm) %>%
  filter(age.type == "both",
         sea.age > 0) %>%
  slice_sample(n = 50000)

data.biph <- sallaa %>%
  filter(age.type == "juve.only") %>% 
  bind_rows(data.biph_m) %>%
  mutate(yday.may1 = if_else(is.na(date), 1, yday(date) - yday("2025-05-01")+1), 
         doy.dec = if_else(yday.may1 > 0, yday.may1, yday.may1+365)/365,
         tot.age.dec = if_else(yday.may1 > 0, tot.age+doy.dec, tot.age+doy.dec-1),
         smo_age = if_else(age.type == "both" | is_smo == 1, juv.age, NA),
         sea.age = if_else(is.na(sea.age), 0, sea.age),
         growth.year = trunc(tot.age.dec)+1, # + 1 to index age 0 individuals
         hatch.year = year-tot.age,
         hatch.year.f = as.numeric(factor(hatch.year, levels = sort(unique(hatch.year)), labels = seq_along(sort(unique(hatch.year))))),
         year.f = as.integer(factor(year)),
         sex = as.integer(factor(sex)),
         su = as.integer(factor(spat.unit)),
         basin = if_else(str_detect(spat.unit, "AU"),"BS","NA"),
         bs = as.integer(factor(basin))) %>%
  #filter(!(is.na(smo_age))) %>%
  mutate(ind.id = row_number()) 

nK = max(data.biph$smo_age %>% na.omit())
nsu = length(unique(data.biph$spat.unit))

consts = list(nind = nrow(data.biph), # n individuals change when >1 sample per ind.id
              nK = nK,
              eta = 2,
              alpha = c(1, 2, 3, 2, 1, 1)[1:nK],
              nsu = nsu,
              npars = 2,#*nsu, # n parameters in the LKJ dmnorm * number of spatial units 
              age = data.biph$tot.age.dec,# wy tot.age.dec and not tot.age??
              growth_year = data.biph$growth.year,
              doy_dec = data.biph$doy.dec,
              ind_id = data.biph$ind.id,
              su = data.biph$su)

inits <- function(){
  list(sig_l = rexp(1,1/150),
       lb_mu = runif(consts$nsu,13,19),
       sig_par_m = c(rlnorm(consts$npars,log(0.1),sdlog = 0.1)),
       sig_par_f = c(rlnorm(consts$npars,log(0.1),sdlog = 0.1)),
       mu_par_m = c(rnorm(consts$npars,6,0.2)),
       mu_par_f = c(rnorm(1,4,0.2),
                   (rnorm(1,-1,0.2)))
      )
}

# build model
t <- Sys.time()
biph1_0_4.model <- nimbleModel(biph1_0_4.code,
                             constants = consts,
                             inits = inits(),
                             data = data.biph %>% select(length_mm, smo_age),
                             buildDerivs = TRUE
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
monits <- c(mvars)

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
biph1_0_4.samples <- runMCMC(biph1_0_4.mcmcc, niter = 4000, nburnin = 2500, nchains = 1, samplesAsCodaMCMC = TRUE, WAIC = TRUE)
#biph.samples <- runMCMC(biph1_0_4.mcmcc, niter = 100000, nburnin = 98500, nchains = 1, WAIC = TRUE)
Sys.time() - t

# biph1_0.samples$WAIC 311294
biph.samples$WAIC #311301.6

# Save samples
saveRDS(biph.samples, file = paste0(home,"/data/biph1_0_4_samples_",Sys.Date(),".RData"))

# Biphasic model 1.0_2:
# Add basin level and update based on improvements in the spattemporal model. fixing smolt age and priors

# Load libraries, install if needed
pkgs <- c("tidyverse", "nimble")

if(length(setdiff(pkgs,rownames(installed.packages()))) > 0){
  install.packages(setdiff(pkgs, rownames(installed.packages())), dependencies = T)
}

invisible(lapply(pkgs, library, character.only = T))

## Biphasic code
biph.code <- nimbleCode({
  
  # k & linf values (g and linf from fb)
  #l_sd_m <- sqrt(log(1 + (278^2) / (1378^2))) # inforamtive
  l_sd_m <- sqrt(log(1 + (178^2) / (1378^2))) # lower than informative to increase shrinkage
  #l_sd_f <- sqrt(log(1 + (45^2) / (185^2)))  # informative 
  l_sd_f <- sqrt(log(1 + (25^2) / (185^2))) # lower than informative to increase shrinkage
  l_mean_m <- log(1378) - 0.5 * l_sd_m^2 # sd of the mean based on sd:s of LKJ
  l_mean_f <- log(185) - 0.5 * l_sd_f^2 # sd of the mean based on sd:s of LKJ
  #k_sd_m <- sqrt(log(1 + (0.28^2) / (0.43^2))) # inforamtive 
  k_sd_m <- sqrt(log(1 + (0.14^2) / (0.43^2))) # lower than informative to increase shrinkage
  k_sd_f <- sqrt(log(1 + (0.1^2) / (0.3^2)))
  k_mean_m <- log(0.43) - 0.5 * k_sd_m^2
  k_mean_f <- log(0.3) - 0.5 * k_sd_f^2 
  lb_sd <- sqrt(log(1 + (3^2) / (17^2)))
  lb_mean <- log(17) - 0.5 * (lb_sd^2)
  
  # https://onlinelibrary.wiley.com/doi/full/10.1111/j.1095-8649.2012.03370.x
  # mean(c(152,167,195,174,154,164,160,215,152,167,195,174,166,211,163,172,207,223,171,171,176,235,173,160,215,187,235,22,146,147,192,304,239,152,198,280,207,223))
  # sd(c(152,167,195,174,154,164,160,215,152,167,195,174,166,211,163,172,207,223,171,171,176,235,173,160,215,187,235,22,146,147,192,304,239,152,198,280,207,223))
  # lbmu / sd DOI: 10.1111/j.1095-8649.2009.02497.x  sd is guesstimate
  
  # mean and sd of freshwater growth parameters
  par_mu_f[1] ~ dnorm(l_mean_f, sd = 0.3)
  par_mu_f[2] ~ dnorm(k_mean_f, sd = 0.3)
  par_sig_f[1] ~ dlnorm(log(l_sd_f), sdlog = 0.5)
  par_sig_f[2] ~ dlnorm(log(k_sd_f), sdlog = 0.5)
  par_mu_m[1] ~ dnorm(l_mean_m, sd = 0.3)
  par_mu_m[2] ~ dnorm(k_mean_m, sd = 0.3)
  par_sig_m[1] ~ dlnorm(log(l_sd_m), sdlog = 0.5)
  par_sig_m[2] ~ dlnorm(log(k_sd_m), sdlog = 0.5)
  # l_mean_m_bs[1] ~ dnorm(l_mean_m, sd = 0.3)
  # l_mean_m_bs[2] ~ dnorm(l_mean_m, sd = 0.3)
  # k_mean_m_bs[1] ~ dnorm(k_mean_m, sd = 0.3)
  # k_mean_m_bs[2] ~ dnorm(k_mean_m, sd = 0.3)
  # l_sd_m_bs[1] ~ dlnorm(log(l_sd_m), sdlog = 0.5)
  # l_sd_m_bs[2] ~ dlnorm(log(l_sd_m), sdlog = 0.5)
  # k_sd_m_bs[1] ~ dlnorm(log(k_sd_m), sdlog = 0.5)
  # k_sd_m_bs[2] ~ dlnorm(log(k_sd_m), sdlog = 0.5)
  
  # Baltic and Atlantic basins
  #for(i in 1:nsu){
  # for(i in 1:nbs){
  #   par_mu_m[1,i] ~ dnorm(l_mean_m_bs[i], sd = 0.2)
  #   par_mu_m[1,(i+2)] ~ dnorm(l_mean_m_bs[i], sd = 0.2)
  #   par_mu_m[2,i] ~ dnorm(k_mean_m_bs[subs[i]], sd = 0.2)
  #   par_mu_m[2,(i+2)] ~ dnorm(k_mean_m_bs[subs[i]], sd = 0.2)
  #   par_sig_m[1,i] ~ dlnorm(log(l_sd_m_bs[subs[i]]), sdlog = 0.2)
  #   par_sig_m[1,(i+2)] ~ dlnorm(log(l_sd_m_bs[subs[i]]), sdlog = 0.2)
  #   par_sig_m[2,i] ~ dlnorm(log(k_sd_m_bs[subs[i]]), sdlog = 0.2)
  #   par_sig_m[2,(i+2)] ~ dlnorm(log(k_sd_m_bs[subs[i]]), sdlog = 0.2)
  # }
  
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
    #par_m[k,1:npars] ~ dmnorm(par_mu_m[1:npars,subs[k]], prec = tau_p_m[1:npars,1:npars])
    par_m[k,1:npars] ~ dmnorm(par_mu_m[1:npars], prec = tau_p_m[1:npars,1:npars])
    par_f[k,1:npars] ~ dmnorm(par_mu_f[1:npars], prec = tau_p_f[1:npars,1:npars])
  }
  
  # priors for smolt age, lb.mu, g
  for(i in 1:nsu){
    smp[1:nK,i] ~ ddirch(alpha[1:nK])
    lb_mu[i] ~ dlnorm(lb_mean, sdlog = lb_sd)
  }
  
  # Priors
  sig_l ~ dexp(1/350)
  
  # likelihood
  for(i in 1:nind){
    length_mm[i] ~ dnorm(l_mu[growth_year[i], ind_id[i]], sd = sig_l)
    
    l_mu[1,i] <- lb_mu[su[i]] +
      (1-step(1-(age[i]+0.001)))*(exp(par_f[su[i],1]) - lb_mu[su[i]])*(1 - exp(-exp(par_f[su[i],2]))) +
      step(1-(age[i]+0.001)) *(exp(par_f[su[i],1]) - lb_mu[su[i]])*(1 - exp(-exp(par_f[su[i],2])))*doy_dec[i]
    
    for(a in 2:(growth_year[i]+1)){
      l_mu[a,i] <- l_mu[(a-1),i] +
        (1-step(a-(age[i]+0.001)))*(step(smo_age[i] - (a-0.999))*((exp(par_f[su[i],1]) - l_mu[(a-1),i])*(1 - exp(-exp(par_f[su[i],2])))) +
                                      (1-step(smo_age[i] - (a-0.999)))*(((exp(par_m[su[i],1])) - l_mu[(a-1),i])*(1 - exp(-(exp(par_m[su[i],2])))))) +
        
        step(a-(age[i]+0.001))* (step(smo_age[i] - (a-0.999))*((exp(par_f[su[i],1]) - l_mu[(a-1),i])*(1 - exp(-exp(par_f[su[i],2]))))*doy_dec[i] +
                                   (1-step(smo_age[i] - (a-0.999)))*(((exp(par_m[su[i],1])) - l_mu[(a-1),i])*(1 - exp(-(exp(par_m[su[i],2])))))*doy_dec[i]) 
      
    }
  }
  
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
  filter(age.type == "both") %>%
  slice_sample(n = 30000)

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
  #filter(!(is.na(smo.age))) %>%
  mutate(ind.id = row_number()) 

nK = max(data.biph$smo_age %>% na.omit())
nsu = length(unique(data.biph$spat.unit))
#subs = data.biph %>% distinct(su,bs) %>% arrange(su)

consts = list(nind = nrow(data.biph), # n individuals change when >1 sample per ind.id
              nK = nK,
              eta = 2,
              alpha = c(1, 2, 3, 2, 1, 1)[1:nK],
              nsu = nsu,
              npars = 2,#*nsu, # n parameters in the LKJ dmnorm * number of spatial units 
              age = data.biph$tot.age.dec,
              growth_year = data.biph$growth.year,
              doy_dec = data.biph$doy.dec,
              #subs = subs$bs,
              ind_id = data.biph$ind.id,
              su = data.biph$su)

inits <- function(){
  list(sig_l = rexp(1,1/150),
       lb_mu = runif(consts$nsu,13,19),
       par_sig_m = c(rlnorm(consts$npars,log(0.1),sdlog = 0.1)),
       par_sig_f = c(rlnorm(consts$npars,log(0.1),sdlog = 0.1)),
       par_mu_m = c(rnorm(consts$npars,7,0.2)),
       par_mu_f = c(rnorm(1,4,0.2),
                   (rnorm(1,-1,0.2))),
       smo_age = sample(2:4,consts$nind, replace = TRUE)
      )
}

# build model
t <- Sys.time()
biph.model <- nimbleModel(biph.code,
                          constants = consts,
                          inits = inits(),
                          data = data.biph %>% select(length_mm,smo_age),
                          buildDerivs = TRUE
                          )

biph.model$simulate()
biph.model$calculate()
biph.model$initializeInfo()

Sys.time() - t 
# identify nodes to sample 
dataNodes <- biph.model$getNodeNames(dataOnly = TRUE)
deterNodes <- biph.model$getNodeNames(determOnly = TRUE)
parentNodes <- biph.model$getParents(dataNodes, stochOnly = TRUE) #all of these should be added to monitor below to recreate other model variables...
stnodes <- biph.model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
allvars <- biph.model$getVarNames(nodes = stnodes)
mvars <- allvars[!(grepl("lifted",allvars))]  

# calculate vars
for(i in 1:length(mvars)){
  print(paste0(mvars[i]," ",biph.model$calculate(mvars[i]) ))
}

# compile model
biph.c <- compileNimble(biph.model)

# configure and build mcmc
monits <- c(mvars, "l_mu")

biph.confmcmc <- configureHMC(biph.c, monitors = monits, enableWAIC = TRUE,#,"bL","bk"),
                                useConjugacy = FALSE)

t <- Sys.time()
biph.mcmc <- buildMCMC(biph.confmcmc, project = biph.model)
Sys.time() - t

# compile mcmc
t <- Sys.time()
biph.mcmcc <- compileNimble(biph.mcmc, project = biph.model)
Sys.time() - t

# MCMC Samples
t <- Sys.time()
biph.samples <- runMCMC(biph.mcmcc, niter = 6000, nburnin = 4500, nchains = 1, WAIC = TRUE)
Sys.time() - t

biph.samples$WAIC # 906068.9
# Save samples
saveRDS(biph.samples, file = paste0(home,"/data/biph1_0_1_samples_",Sys.Date(),".RData"))

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
  l_sd_m <- sqrt(log(1 + 278^2 / 1378^2)) # informative
  l_sd_f <- sqrt(log(1 + 21^2 / 171^2))  # informative 
  l_mean_m <- log(1378) - 0.5 * l_sd_m^2 # sd of the mean based on sd:s of LKJ
  l_mean_f <- log(171) - 0.5 * l_sd_f^2 # sd of the mean based on sd:s of LKJ
  k_sd_m <- sqrt(log(1 + 0.28^2 / 0.43^2))
  k_sd_f <- sqrt(log(1 + 0.28^2 / 0.61^2))
  k_mean_m <- log(0.43) - 0.5 * k_sd_m^2
  k_mean_f <- log(0.61) - 0.5 * k_sd_f^2
  lb_sd <- sqrt(log(1 + 3^2 / 17^2))
  lb_mean <- log(17) - 0.5 * lb_sd^2
  
  # https://onlinelibrary.wiley.com/doi/full/10.1111/j.1095-8649.2012.03370.x # using only wild origin
  # mean(c(152,167,195,174,154,164,152,167,195,174,166,171,171,176,173,187,235,146,147,152))
  # sd(c(152,167,195,174,154,164,152,167,195,174,166,171,171,176,173,187,235,146,147,152))
  # for juvenile K, using mean for white spotted charr in https://doi.org/10.1111/j.1095-8649.2001.tb00220.x
  # lbmu / sd DOI: 10.1111/j.1095-8649.2009.02497.x  sd is guesstimate
  
  # mean and sd of freshwater growth parameters
  par_mu_f[1] ~ dnorm(l_mean_f, sd = 0.1)
  par_mu_f[2] ~ dnorm(k_mean_f, sd = 0.1)
  par_sig_f[1] ~ dlnorm(log(l_sd_f), sdlog = 0.3)
  par_sig_f[2] ~ dlnorm(log(k_sd_f), sdlog = 0.3)
  par_mu_m[1] ~ dnorm(l_mean_m, sd = 0.1)
  par_mu_m[2] ~ dnorm(k_mean_m, sd = 0.1)
  par_sig_m[1] ~ dlnorm(log(l_sd_m), sdlog = 0.3) #try 0.1 for all par_sig sd:s
  par_sig_m[2] ~ dlnorm(log(k_sd_m), sdlog = 0.3)
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
  
  # estimate smo_age
  for(l in 1:nind){
    smo_age[l] ~ dcat(smp[1:nK,su[l]])
  }
  # Priors
  sig_l ~ dexp(1/300)
  
  # likelihood
  for(i in 1:nind){
    length_mm[i] ~ dnorm(l_mu[gy[i], ind_id[i]], sd = sig_l)
    
    l_mu[1,i] <- growth_step(
      l_prev  = lb_mu[su[i]],
      L_f     = exp(par_f[su[i],1]),
      L_m     = exp(par_m[su[i],1]),   
      k_f     = exp(par_f[su[i],2]),
      k_m     = exp(par_m[su[i],2]),
      doy_dec = doy_dec[i],
      j       = 1,
      age     = age[i],
      smo_age = smo_age[i]
      )
      
    # ) lb_mu[su[i]] +
    #   (1-step(1-(age[i]+0.001)))*(exp(par_f[su[i],1]) - lb_mu[su[i]])*(1 - exp(-exp(par_f[su[i],2]))) +
    #   step(1-(age[i]+0.001)) *(exp(par_f[su[i],1]) - lb_mu[su[i]])*(1 - exp(-exp(par_f[su[i],2])))*doy_dec[i]
    
    for(j in 2:(gy[i]+1)){
      l_mu[j, i] <- growth_step(
        l_prev  = l_mu[(j-1), i],
        L_f     = exp(par_f[su[i],1]),
        L_m     = exp(par_m[su[i],1]),   
        k_f     = exp(par_f[su[i],2]),
        k_m     = exp(par_m[su[i],2]),
        doy_dec = doy_dec[i],
        j       = j,
        age     = age[i],
        smo_age = smo_age[i]
        )
      }
      # l_mu[j,i] <- l_mu[(a-1),i] +
      #   (1-step(a-(age[i]+0.001)))*(step(smo_age[i] - (a-0.999))*((exp(par_f[su[i],1]) - l_mu[(a-1),i])*(1 - exp(-exp(par_f[su[i],2])))) +
      #                                 (1-step(smo_age[i] - (a-0.999)))*(((exp(par_m[su[i],1])) - l_mu[(a-1),i])*(1 - exp(-(exp(par_m[su[i],2])))))) +
      #   
      #   step(a-(age[i]+0.001))* (step(smo_age[i] - (a-0.999))*((exp(par_f[su[i],1]) - l_mu[(a-1),i])*(1 - exp(-exp(par_f[su[i],2]))))*doy_dec[i] +
      #                              (1-step(smo_age[i] - (a-0.999)))*(((exp(par_m[su[i],1])) - l_mu[(a-1),i])*(1 - exp(-(exp(par_m[su[i],2])))))*doy_dec[i]) 
      
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

# l_mu[1,i] <- growth_step(
#   l_prev  = lb_mu[su[i]],
#   L_f     = exp(par_f[su[i],1]),
#   L_m     = exp(par_m[su[i],1]),   
#   k_f     = exp(par_f[su[i],2]),
#   k_m     = exp(par_m[su[i],2]),
#   doy_dec = doy_dec[i],
#   j       = 1,
#   age     = age[i],
#   smo_age = smo_age[i]
# )

growth_step <- nimbleFunction(
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
    }, buildDerivs = TRUE)


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
  #filter(!(is.na(smo.age))) %>%
  mutate(ind.id = row_number()) 

nK = max(data.biph$smo_age %>% na.omit())
nsu = length(unique(data.biph$spat.unit))
#subs = data.biph %>% distinct(su,bs) %>% arrange(su)

consts = list(nind = nrow(data.biph), # n individuals change when >1 sample per ind.id
              nK = nK,
              eta = 2,
              alpha = c(2, 3, 2, 1, 1, 1)[1:nK],
              nsu = nsu,
              npars = 2,#*nsu, # n parameters in the LKJ dmnorm * number of spatial units 
              age = data.biph$tot.age.dec,
              gy = data.biph$growth.year,
              doy_dec = data.biph$doy.dec,
              #subs = subs$bs,
              ind_id = data.biph$ind.id,
              su = data.biph$su)

initsmo_age = data.biph %>% 
  mutate(sa = if_else(is.na(smo_age), 2, NA)) %>% pull(sa)
inits <- function(){
  list(sig_l = rexp(1,1/150),
       lb_mu = runif(consts$nsu,5,10),
       par_sig_m = c(rlnorm(consts$npars,log(0.1),sdlog = 0.1)),
       par_sig_f = c(rlnorm(consts$npars,log(0.1),sdlog = 0.1)),
       par_mu_m = c(rnorm(consts$npars,10,0.2)),
       par_mu_f = c(rnorm(1,4,0.2),
                   (rnorm(1,-1,0.2))),
       smo_age = initsmo_age,
       l_mu = matrix(runif(consts$nind*12, 50,1000), nrow= 12, ncol = consts$nind))
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
biph.c <- compileNimble(biph.model, showCompilerOutput = TRUE)

# configure and build mcmc
monits <- c(mvars[-c(13,10)])
#monits <- c(mvars[-c(13)], "l_mu")

#biph.confmcmc <- configureHMC(biph.c, monitors = monits, enableWAIC = TRUE)#,#,"bL","bk"),
#                                useConjugacy = FALSE, delta = 0.6)

biph.confmcmc <- configureMCMC(biph.c, monitors = monits, enableWAIC = TRUE)#,#,"bL","bk"),
#                                useConjugacy = FALSE, delta = 0.6)

t <- Sys.time()
biph.mcmc <- buildMCMC(biph.confmcmc, project = biph.model)
Sys.time() - t

# compile mcmc
biph.mcmcc <- compileNimble(biph.mcmc, project = biph.model)
t <- Sys.time()
Sys.time() - t

# MCMC Samples
t <- Sys.time()
#biph.samples <- runMCMC(biph.mcmcc, niter = 6000, nburnin = 4500, nchains = 1, WAIC = TRUE)
biph.samples <- runMCMC(biph.mcmcc, niter = 100000, nburnin = 98500, nchains = 1, WAIC = TRUE)
Sys.time() - t

biph.samples$WAIC # 906068.9
# Save samples
saveRDS(biph.samples, file = paste0(home,"/data/biph1_0_3_samples_",Sys.Date(),".RData"))

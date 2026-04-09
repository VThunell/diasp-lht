# Biphasic model 1.1:
# spatiotemporal version of v1.0

# Load libraries, install if needed
pkgs <- c("tidyverse", "nimble")

if(length(setdiff(pkgs,rownames(installed.packages()))) > 0){
  install.packages(setdiff(pkgs, rownames(installed.packages())), dependencies = T)
}

invisible(lapply(pkgs, library, character.only = T))

## Biphasic 
biph1_1.code <- nimbleCode({
  
  # likelihood
  for(i in 1:nind){
    length_mm[i] ~ dnorm(l_mu[growth_year[i], ind_id[i]], sd = sig_l)
    
    # length prediction 0+
    l_mu[1,i] <- lb_mu[su[i]] +
      (1-step(1-(age[i]+0.001)))*(exp(par_f[su[i],1]) - lb_mu[su[i]])*(1 - exp(-exp(par_f[su[i],2]))) +
      step(1-(age[i]+0.001))*(exp(par_f[su[i],1]) - lb_mu[su[i]])*(1 - exp(-exp(par_f[su[i],2])*doy_dec[i]))
    
    # length prediction 1+ and older
    for(j in 2:(growth_year[i]+1)){
      l_mu[j,i] <- l_mu[(j-1),i] +
        (1-step(j-(age[i]+0.001)))*(step(smo_age[i] - (j-0.999))*((exp(par_f[su[i],1]) - l_mu[(j-1),i])*(1 - exp(-exp(par_f[su[i],2])))) +
                                      (1-step(smo_age[i] - (j-0.999)))*(((exp(lf[(h_year[i]+(j-1)),su[i],1])) - l_mu[(j-1),i])*(1 - exp(-(exp(lf[(h_year[i]+(j-1)),su[i],2])))))) +
                                      
        step(j-(age[i]+0.001))*(step(smo_age[i] - (j-0.999))*((exp(par_f[su[i],1]) - l_mu[(j-1),i])*(1 - exp(-exp(par_f[su[i],2])*doy_dec[i]))) +
                                      (1-step(smo_age[i] - (j-0.999)))*(((exp(lf[(h_year[i]+(j-1)),su[i],1])) - l_mu[(j-1),i])*(1 - exp(-(exp(lf[(h_year[i]+(j-1)),su[i],2])*doy_dec[i])))))
        
      # l_mu[1,i] <- lb_mu[su[i]] +
      #   (1-step(1-(age[i]+0.001)))*(L_f[su[i],1] - lb_mu[su[i]])*(1 - exp(-k_f[su[i],2])) +
      #   step(1-(age[i]+0.001))*(L_f[su[i],1] - lb_mu[su[i]])*(1 - exp(-k_f[su[i],2]))*doy_dec[i]
      # 
      # # this (step(smo_age[i] - (j-0.999))) is 0 when j>smo.age. if age = smo.age, we model the growth year as marine
      # #NOT SEX length prediction 1+ and older
      # for(j in 2:(growth_year[i]+1)){
      #   l_mu[j,i] <- l_mu[(j-1),i] +
      #     (1-step(j-(age[i]+0.001)))*(step(smo_age[i] - (j-0.999))*((L_f[su[i],1] - l_mu[(j-1),i])*(1 - exp(-k_f[su[i],2]))) +
      #                                   (1-step(smo_age[i] - (j-0.999)))*((L_m[su[i],1] - l_mu[(j-1),i])*(1 - exp(-k_m[su[i],2])))) +
      #     
      #     step(j-(age[i]+0.001))*(step(smo_age[i] - (j-0.999))*((L_f[su[i],1] - l_mu[(j-1),i])*(1 - exp(-k_f[su[i],2]*doy_dec[i]))) +
      #                               (1-step(smo_age[i] - (j-0.999)))*((L_m[su[i],1] - l_mu[(j-1),i])*(1 - exp(-k_m[su[i],2]*doy_dec[i]))))
      #   
      #   L_f[su[i],1] <- exp(par_f[su[i],1])
      #   k_f[su[i],2] <- exp(par_f[su[i],2])
      #   L_m[su[i],1] <- exp(m_pl[(h_year[i]+(j-1)),su[i],1])
      #   k_m[su[i],2] <- exp(m_pl[(h_year[i]+(j-1)),su[i],2])
    }
  }
  
  # estimate smo_age
  for(l in 1:nind){
    smo_age[l] ~ dcat(smp[1:nK,su[l]])
  }
  
  # Priors
  sig_l ~ dexp(1/350)

  lf[1:nyear,1:nsu,1] <- par_m[1:nyear,1:nsu]
  lf[1:nyear,1:nsu,2] <- par_m[1:nyear,(nsu+1):nmpars]

  # first year of vb-par rw
  for(j in 1:nsu){
    par_m[1,j] ~ dnorm(mu_par_m[j], sd = l_lsd_m)
    par_m[1,(nsu+j)] ~ dnorm(mu_par_m[(nsu+j)], sd = k_lsd_m)
    
    par_f[j,1:nfpars] ~ dmnorm(mu_par_f[1:nfpars], prec = tau_p_f[1:nfpars,1:nfpars])
  }
  
  # transition of vB pars from year i to i+1 (1 for each sex)
  for(i in 1:(nyear-1)){
    par_m[(i+1),1:nmpars] ~ dmnorm(par_m[i,1:nmpars], prec = tau_p_m[1:nmpars,1:nmpars])
    #par_m[(i+1),(nsu+1):nmpars] <- par_m[i,1:nmpars]# for no temp var in linf
    #par_m[1:nsu,2,(i+1)] ~ dmnorm(par_m[1:nsu,2,i], prec = tau_p[1:nsu,1:nsu]) # for temp var in linf
  }
  
  tau_p_m[1:nmpars,1:nmpars] <- inverse(sigma_p_m[1:nmpars, 1:nmpars])
  tau_p_f[1:nfpars,1:nfpars] <- inverse(sigma_p_f[1:nfpars, 1:nfpars])
  
  for(k in 1:nmpars){
    for(j in 1:nmpars){
      sigma_p_m[k,j] <- Rnew[k,j] * sig_par_m[k] * sig_par_m[j]
    }
  }
  
  for(k in 1:nfpars){
    for(j in 1:nfpars){
      sigma_p_f[k,j] <- Rnew_f[k,j] * sig_par_f[k] * sig_par_f[j]
    }
  }
  
  # priors for smolt age, lb.mu, g
  for(i in 1:nsu){
    #bL[i] ~ dnorm(mean = 0, sd = 1) 
    #bk[i] ~ dnorm(mean = 0, sd = 1) 
    smp[1:nK,i] ~ ddirch(alpha[1:nK])
    lb_mu[i] ~ dnorm(17, sd = 3) 
    #lb.mu[i] ~ dunif(12, 22) # ~15 to 20 mm, DOI: 10.1111/j.1095-8649.2009.02497.x  sd is guesstimate
  }
  
  for(i in 1:nsu){
    #mu_par_m[i] <- l_lmean_m 
    mu_par_m[i] ~ dnorm(l_lmean_m, sd = 1)
    mu_par_m[(nsu+i)] ~ dnorm(mean = k_lmean_m, sd = 0.5)  
    #sig_par_m[i] <- log(l_lsd_m)
    sig_par_m[i] ~ dlnorm(log(l_lsd_m), sdlog = 0.1)
    sig_par_m[(nsu+i)] ~ dlnorm(log(k_lsd_m), sdlog = 0.1)
  }
  
  mu_par_f[1] ~ dnorm(l_lmean_f, sd = 1)  
  mu_par_f[2] ~ dnorm(mean = k_lmean_f, sd = 0.5)
  sig_par_f[1] ~ dlnorm(log(l_lsd_f), sdlog = 0.1)  
  sig_par_f[2] ~ dlnorm(log(k_lsd_f), sdlog = 0.1)
  
  # k & linf values (g and linf from fb)
  l_lsd_m <- sqrt(log(1 + (278^2) / (1378^2)))
  l_lsd_f <- sqrt(log(1 + (45^2) / (185^2)))
  l_lmean_m <- log(1378) - 0.5 * log(l_lsd_m)^2
  l_lmean_f <- log(185) - 0.5 * log(l_lsd_f)^2  
  k_lsd_m <- sqrt(log(1 + (0.28^2) / (0.43^2)))
  k_lsd_f <- sqrt(log(1 + (0.15^2) / (0.31^2)))
  k_lmean_m <- log(0.43) - 0.5 * log(k_lsd_m)^2
  k_lmean_f <- log(0.43) - 0.5 * log(k_lsd_f)^2 
  
  # https://onlinelibrary.wiley.com/doi/full/10.1111/j.1095-8649.2012.03370.x
  # mean(c(152,167,195,174,154,164,160,215,152,167,195,174,166,211,163,172,207,223,171,171,176,235,173,160,215,187,235,22,146,147,192,304,239,152,198,280,207,223))
  # sd(c(152,167,195,174,154,164,160,215,152,167,195,174,166,211,163,172,207,223,171,171,176,235,173,160,215,187,235,22,146,147,192,304,239,152,198,280,207,223))
  
  
  #Prior for correlation matrices (LKJ prior)
  # freshwater
  phi_f[1]  <- eta + (nfpars - 2)/2
  corY_f[1] ~ dbeta(phi_f[1], phi_f[1])
  r12_f   <- 2 * corY_f[1] - 1
  ##
  R_f[1,1]     <- 1
  R_f[1,2]     <- r12_f
  R_f[2,2]     <- sqrt(1 - r12^2)
  
  R_f[2,1]   <- 0
  
  Rnew_f[1:nfpars,1:nfpars] <- t(R_f[1:nfpars,1:nfpars]) %*% R_f[1:nfpars,1:nfpars]

  # marine 
  phi[1]  <- eta + (nmpars - 2)/2
  corY[1] ~ dbeta(phi[1], phi[1])
  r12   <- 2 * corY[1] - 1
  ##
  R[1,1]     <- 1
  R[1,2]     <- r12
  R[2,2]     <- sqrt(1 - r12^2)
  
  #R[2,1]   <- 0
  
  R[2:nmpars,1]  <- 0
  
  for (m in 2:(nmpars-1)) {
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
    for(jk in (m+1):nmpars){
      R[jk,m] <- 0
    }
  }  #m
  
  Rnew[1:nmpars,1:nmpars] <- t(R[1:nmpars,1:nmpars]) %*% R[1:nmpars,1:nmpars]
  
})

#Function creating the Cholesky of the covar. matrix (p45 Nimble manual)
# uppertri_mult_diag <- nimbleFunction(
#   run = function(mat = double(2), vec = double(1)) {
#     returnType(double(2))
#     p <- length(vec)
#     out <- matrix(nrow = p, ncol = p, init = FALSE)
#     for(k in 1:p)
#       out[ , k] <- mat[ , k] * vec[k]
#     return(out)
#     # turn off buildDerivs for the i index
#   }, buildDerivs = list(run = list(ignore = c('k'))))

set.seed(1)
data.biph_m <- sallaa %>%
  filter(age.type == "both") %>%
  slice_sample(n = 30000) 

unique(data.biph$spat.unit)
data.biph <- sallaa %>%
  filter(age.type == "juve.only") %>% 
  bind_rows(data.biph_m) %>%
  # combine Irish units and see if that helps the model
  mutate(spat.unit = if_else(spat.unit %in% c("Ireland-west","Ireland-south"), "Ireland",spat.unit)) %>%
  mutate(yday.may1 = if_else(is.na(date), 1, yday(date) - yday("2025-05-01")+1), 
         doy.dec = if_else(yday.may1 > 0, yday.may1, yday.may1+365)/365,
         tot.age.dec = if_else(yday.may1 > 0, tot.age+doy.dec, tot.age+doy.dec-1),
         smo.age = if_else(age.type == "both" | is_smo == 1, juv.age, NA),
         sea.age = if_else(is.na(sea.age), 0, sea.age),
         growth.year = trunc(tot.age.dec)+1, # + 1 to index age 0 individuals
         hatch.year = year-tot.age,
         hatch.year.f = as.numeric(factor(hatch.year, levels = sort(unique(hatch.year)), labels = seq_along(sort(unique(hatch.year))))),
         year.f = as.integer(factor(year)),
         sex = as.integer(factor(sex)),
         su = as.integer(factor(spat.unit))) %>%
  #filter(!(is.na(smo.age))) %>%
  mutate(ind.id = row_number()) 
  
nK = max(data.biph$smo.age %>% na.omit())
nsu = length(unique(data.biph$spat.unit))
maxyear = max(data.biph$hatch.year.f + data.biph$growth.year)
minyear = min(data.biph$hatch.year.f)
nyear = length(minyear:maxyear)

consts = list(nind = nrow(data.biph), # n individuals change when >1 sample per ind.id
              nK = nK,
              eta = 2,
              alpha = c(1, 2, 3, 2, 1, 1)[1:nK],
              nsu = nsu,
              nfpars = 2, 
              nmpars = nsu*2,
              age = data.biph$tot.age.dec,
              growth_year = data.biph$growth.year,
              doy_dec = data.biph$doy.dec,
              ind_id = data.biph$ind.id,
              h_year = data.biph$hatch.year.f,
              #shfa = fec.dpars$shfa,
              #rafa = fec.dpars$rafa,
              #mufb = fec.dpars$mufb,
              #sdfb = fec.dpars$sdfb,
              #sex = data.biph$sex,
              nyear = nyear,            
              su = data.biph$su)

inits <- function(){
  list(sig_l = rexp(1,1/150),
       lb_mu = runif(consts$nsu,13,19),
       sig_par_m = c(rlnorm(consts$nsu,log(0.2),sdlog = 0.1),
                     rlnorm(consts$nsu,log(0.6),sdlog = 0.1)),
       sig_par_f = c(rlnorm(1,log(0.2),sdlog = 0.1),
                     rlnorm(1,log(0.6),sdlog = 0.1)),
       mu_par_m = c(rnorm(consts$nsu,7,0.2),
                    rnorm(consts$nsu,-1,0.4)),
       mu_par_f = c(rnorm(1,2,0.2),
                    rnorm(1,-1,0.6)),
       par_m = matrix(rlnorm(consts$nmpars,log(0.2),sdlog = 0.1)*rlnorm(consts$nmpars,log(0.6),sdlog = 0.1),
                      nrow = consts$nyear, ncol = consts$nmpars),
       par_f = matrix(rlnorm(consts$nsu,log(0.2),sdlog = 0.1)*rlnorm(consts$nsu,log(0.6),sdlog = 0.1),
                      nrow = consts$nsu, ncol = consts$nfpars),
       corZ = matrix(rnorm((consts$nmpars-1)*(consts$nmpars-1), 0, 0.1),
                      nrow = consts$nmpars-1, ncol = consts$nmpars-1),
       # corY_f = matrix(rnorm((consts$nfpars-1)*(consts$nfpars-1), 0, 0.1),
       #                 nrow = consts$nfpars-1, ncol = consts$nfpars-1),
       corY = runif((consts$nmpars-1), 0, 0.1),
       corY_f = runif((consts$nfpars-1), 0, 0.1),
       smp = matrix(rep(rdirch(length(consts$alpha),consts$alpha),13),
                    nrow = length(consts$alpha), ncol = consts$nsu),
       smo_age = sample(1:5,consts$nind,replace = TRUE)
       #bL = rnorm(consts$nsu, mean = 0, sd = .1), 
       #bk = rnorm(consts$nsu, mean = 0, sd = .1)
      )
}

# build model
t <- Sys.time()
biph1_1.model <- nimbleModel(biph1_1.code,
                             constants = consts,
                             inits = inits(),
                             data = data.biph %>% select(length_mm),
                             buildDerivs = TRUE
                           )

biph1_1.model$simulate()
biph1_1.model$calculate()
biph1_1.model$initializeInfo()
Sys.time() - t 
# identify nodes to sample 
dataNodes <- biph1_1.model$getNodeNames(dataOnly = TRUE)
deterNodes <- biph1_1.model$getNodeNames(determOnly = TRUE)
parentNodes <- biph1_1.model$getParents(dataNodes, stochOnly = TRUE) #all of these should be added to monitor below to recreate other model variables...
stnodes <- biph1_1.model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
allvars <- biph1_1.model$getVarNames(nodes = stnodes)
mvars <- allvars[!(grepl("lifted",allvars))]  

# calculate vars
for(i in 1:length(mvars)){
  print(paste0(mvars[i]," ",biph1_1.model$calculate(mvars[i]) ))
}

# compile model
biph1_1.c <- compileNimble(biph1_1.model)

# configure and build mcmc
monits <- c(mvars, "l_mu")

biph1_1.confmcmc <- configureHMC(biph1_1.c, monitors = monits, enableWAIC = TRUE,#,"bL","bk"),
                                useConjugacy = FALSE)

t <- Sys.time()
biph1_1.mcmc <- buildMCMC(biph1_1.confmcmc, project = biph1_1.model)
Sys.time() - t

# compile mcmc
t <- Sys.time()
biph1_1.mcmcc <- compileNimble(biph1_1.mcmc, project = biph1_1.model)
Sys.time() - t

# MCMC Samples
t <- Sys.time()
biph1_1.samples <- runMCMC(biph1_1.mcmcc, niter = 7000, nburnin = 4000, thin = 2, nchains = 1, samplesAsCodaMCMC = TRUE, WAIC = TRUE)
Sys.time() - t
biph1_1.samples$WAIC

# biph1_1.samples$WAIC 312954.1
# biph1_1.samples$WAIC 339280.8
# node.sub <- grep("corY", colnames(biph1_1.samples$samples[[1]]), value = TRUE) 
# biph1_1.samples$samples[, node.sub[], drop = FALSE] %>% gelman.diag() 

# sig_l"     "lb_mu"     "corZ"      "corY_f"    "corY"      "smp"       "sig_par_m" "mu_par_m"  "sig_par_f"
# [10] "mu_par_f"  "smo_age"   "par_m"     "par_f"     "length_mm" "l_mu"     
node.sub <- grep("^par_f", colnames(biph1_1.samples$samples), value = TRUE)
biph1_1.samples$samples[, node.sub[sample(1:length(node.sub), 20)], drop = FALSE] %>% gelman.diag() 
biph1_1.samples$samples[, node.sub[], drop = FALSE] %>% gelman.diag() 

biph1_1.samples$samples[, node.sub[sample(1:length(node.sub), 20)], drop = FALSE] %>% mcmc_trace() 
biph1_1.samples$samples[, node.sub[], drop = FALSE] %>%  effectiveSize()
biph1_1.samples$samples[, node.sub[14:27], drop = FALSE] %>% mcmc_trace()
#gc()

#saveRDS(biph9.samples, file = paste0(home,"/data/biph9_samples_1025.RData"))
#biph9.samples <- as.matrix(biph9.mcmcc$mvSamples)

# Save samples
saveRDS(biph1_1.samples, file = paste0(home,"/data/biph1_1_samples_",Sys.Date(),".RData"))

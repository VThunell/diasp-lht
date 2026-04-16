# Biphasic model 10
# test nimbleFnction for l.mu, works fine until buildMCMC.

# Load libraries, install if needed
pkgs <- c("tidyverse", "nimble")

if(length(setdiff(pkgs,rownames(installed.packages()))) > 0){
  install.packages(setdiff(pkgs, rownames(installed.packages())), dependencies = T)
}

invisible(lapply(pkgs, library, character.only = T))

## Biphasic 10
biph10.code <- nimbleCode({
  
  # likelihood
  for(i in 1:nobs){
    length_mm[i] ~ dnorm(l.mu[growth.year[i], ind.id[i]], sd = sig.l)
  }
  
  l.mu[1:(nyear+1), 1:nind] <- calc_l_mu(
    lb_mu = lb.mu[1:nsu],
    g = g[1:nyear, 1:nsu],
    par = par[1:nyear, 1:nsu, 1:2],
    su = su[1:nind],
    age = age[1:nind],
    smo_age = smo.age[1:nind],
    hatch_year_f = hatch.year.f[1:nind],
    growth_year = growth.year[1:nind],
    doy_dec = doy.dec[1:nind]
  )
  # # Estimate length at age
  # for(i in 1:nind){
  #   
  #   #length prediction 0+
  #   l.mu[1,i] <- lb.mu[su[i]] +
  #     (1-step(1-(age[i]+0.001)))*(exp(g[hatch.year.f[i], su[i]]) + step(1-(age[i]+0.001))*exp(g[hatch.year.f[i], su[i]])*doy.dec[i])
  #   
  #   #NOT SEX length prediction 1+ and older
  #   for(a in 2:(growth.year[i]+1)){
  #     l.mu[a,i] <- l.mu[(a-1),i] +
  #       (1-step(a-(age[i]+0.001)))*(step(smo.age[i] - (a-0.999))*exp(g[(hatch.year.f[i]+(a-1)), su[i]]) +
  #                                     (1-step(smo.age[i] - (a-0.999)))*(((exp(par[(hatch.year.f[i]+(a-1)),su[i],1])) -
  #                                                                          l.mu[(a-1),i])*(1 - exp(-(exp(par[(hatch.year.f[i]+(a-1)),su[i],2])))))) +
  #       
  #       step(a-(age[i]+0.001))*(step(smo.age[i] - (a-0.999))*exp(g[(hatch.year.f[i]+(a-1)), su[i]])*doy.dec[i] +
  #                                 (1-step(smo.age[i] - (a-0.999)))*(((exp(par[(hatch.year.f[i]+(a-1)),su[i],1])) -
  #                                                                      l.mu[(a-1),i])*(1 - exp(-(exp(par[(hatch.year.f[i]+(a-1)),su[i],2])*doy.dec[i])))))
  #     
  #     f.mu[a,i] <- step(smo.age[i] - (a-0.999))*(l.mu[a,i]*fa[su[i]]^fb[su[i]])
  #   }
  # }
  # 
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
    smp[1:nK,i] ~ ddirch(alpha[1:nK])
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
  
  # fecundity fa and fb
  # for(i in 1:nsu){
  #   fa[i] ~ dgamma(shape = shfa[i], rate = rafa[i])
  #   fb[i] ~ dnorm(mean = mufb[i], sd = sdfb[i])
  # }
  
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



# Function for calculating l.mu and f.mu 
calc_l_mu <- nimbleFunction(
  run = function(
    lb_mu = double(1),g = double(2),par = double(3),su = double(1),age = double(1),smo_age = double(1),
    hatch_year_f = double(1),growth_year = double(1),doy_dec = double(1)) {
    returnType(double(2))
    
    nind <- length(su)
    nyears <- dim(g)[1]
    
    # initialize l_mu
    l_mu <- matrix(0.0, nrow = nyears + 1, ncol = nind)
    
    for (i in 1:nind) {
      su_i <- su[i]
      hyf_i <- hatch_year_f[i]
      ag_i <- age[i]
      smo_i <- smo_age[i]
      gy_i <- growth_year[i]
      doy_i <- doy_dec[i]
      
      # initial 0+ prediction
      l_prev <- lb_mu[su_i] +
        (1.0 - step(1.0 - (ag_i + 0.001))) *
        (exp(g[hyf_i, su_i]) + step(1.0 - (ag_i + 0.001)) * exp(g[hyf_i, su_i]) * doy_i)
      
      l_mu[1, i] <- l_prev
      
      # loop over years
      for (a in 1:nyears) {
        # ensure we only update within growth window
        use_year <- step(a - hyf_i) * step((hyf_i + gy_i) - a)
        if (use_year == 1.0) {
          
          pre_smolt <- step(smo_i - (a - hyf_i - 0.999))
          post_smolt <- 1.0 - pre_smolt
          before_age <- 1.0 - step((a - hyf_i) - (ag_i + 0.001))
          
          g_exp <- exp(g[a, su_i])
          linf_exp <- exp(par[a, su_i, 1])
          k_exp <- exp(par[a, su_i, 2])
          
          delta_len <- before_age * (pre_smolt * g_exp + post_smolt * ((linf_exp - l_prev) * (1.0 - exp(-k_exp)))) +
            (1.0 - before_age) * (pre_smolt * g_exp * doy_i + post_smolt * ((linf_exp - l_prev) * (1.0 - exp(-k_exp * doy_i))))
          
          l_new <- l_prev + delta_len
          l_prev <- l_new
          l_mu[a + 1, i] <- l_new
          
        } else {
          # keep previous value if outside growth window
          l_mu[a + 1, i] <- l_prev
        }
      } # end loop over a
    } # end loop over i
    
    return(l_mu)
  }
)

data.biph <- sallaa %>%
  mutate(yday.may1 = if_else(is.na(date), 1, yday(date) - yday("2025-05-01")+1), 
         doy.dec = if_else(yday.may1 > 0, yday.may1, yday.may1+365)/365,
         tot.age.dec = if_else(yday.may1 > 0, tot.age+doy.dec, tot.age+doy.dec-1),
         smo.age = if_else(age.type == "both", juv.age, NA),
         sea.age = if_else(is.na(sea.age), 0, sea.age),
         growth.year = as.integer(trunc(tot.age.dec)+1), # + 1 to index age 0 individuals
         hatch.year = year-tot.age,
         hatch.year.f = as.integer(factor(hatch.year, levels = sort(unique(hatch.year)), labels = seq_along(sort(unique(hatch.year))))),
         year.f = as.integer(factor(year)),
         sex = as.integer(factor(sex)),
         su = as.integer(factor(spat.unit)),
         ind.id = row_number())

nK = max(data.biph$smo.age %>% na.omit())
nsu = length(unique(data.biph$su))
maxyear = max(data.biph$hatch.year.f + data.biph$growth.year)
minyear = min(data.biph$hatch.year.f)

consts = list(nobs = nrow(data.biph),
              nind = nrow(data.biph), # n individuals change when >1 sample per ind.id
              nK = nK,
              alpha = c(1, 2, 3, 2, 1, 1)[1:nK],
              nsu = nsu,
              npars = 2*nsu, # n parameters in the LKJ dmnorm * number of spatial units 
              age = data.biph$tot.age.dec,
              growth.year = data.biph$growth.year,
              doy.dec = data.biph$doy.dec,
              ind.id = data.biph$ind.id,
              hatch.year.f = data.biph$hatch.year.f,
              # shfa = fec.dpars$shfa,
              # rafa = fec.dpars$rafa,
              # mufb = fec.dpars$mufb,
              # sdfb = fec.dpars$sdfb,
              #sex = data.biph$sex,
              maxyear = maxyear,
              minyear = minyear,
              nyear = length(minyear:maxyear),
              su = data.biph$su)

inits <- function(){
  list(sig.l = rexp(1,1/250),
       lb.mu = runif(consts$nsu,13,19),
       Ustar  = diag(consts$npars),
       Ugstar = diag(consts$nsu),
       sig.par = c(rlnorm(consts$nsu,log(0.2),sdlog = 0.1),
                   rlnorm(consts$nsu,log(0.6),sdlog = 0.1)),
       mu.par = c(rnorm(consts$nsu,7,0.2),
                  rnorm(consts$nsu,-1,0.6)),
       mu.par = c(rnorm(consts$nsu,0,0.1),
                  rnorm(consts$nsu,0,0.1)),
       g.sig.par = rlnorm(consts$nsu,0,0.1),
       g.par = rnorm(consts$nsu,4,0.2)
      )
  }

calc_l_mu.c <- compileNimble(calc_l_mu)

# build model
biph10.model <- nimbleModel(biph10.code,
                     constants = consts,
                     inits = inits(),
                     data = data.biph %>% select(length_mm))


biph10.model$initializeInfo()
t <- Sys.time()
biph10.model$simulate()
biph10.model$calculate()
Sys.time() - t

# identify nodes to sample 
dataNodes <- biph10.model$getNodeNames(dataOnly = TRUE)
deterNodes <- biph10.model$getNodeNames(determOnly = TRUE)
parentNodes <- biph10.model$getParents(dataNodes, stochOnly = TRUE) #all of these should be added to monitor below to recreate other model variables...
biph10.model$getVarNames(nodes = parentNodes)
stnodes <- biph10.model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
allvars <- biph10.model$getVarNames(nodes = stnodes)
mvars <- allvars[!(grepl("lifted",allvars))]  

# calculate vars
for(i in 1:length(mvars)){
  print(paste0(mvars[i]," ",biph10.model$calculate(mvars[i]) ))
}

# compile model
biph10.c <- compileNimble(biph10.model)

# configure and build mcmc
biph10.confmcmc <- configureMCMC(biph10.c, monitors = c("sig.l", "Ustar","Ugstar","lb.mu","g","log_p.year","smo.age"),useConjugacy = FALSE)
t <- Sys.time()
biph10.mcmc <- buildMCMC(biph10.confmcmc) # BUILDING THIS TAKES >8 hours
Sys.time() - t
# compile mcmc
biph10.mcmcc <- compileNimble(biph10.mcmc, project = biph10.model)

# MCMC Samples
t <- Sys.time()
biph10.mcmcc$run(200,nburnin=150)
#biph10.samples <- runMCMC(biph10.mcmcc, niter = 250000, nburnin = y35000, thin=10)#, WAIC=TRUE) 
Sys.time() - t

saveRDS(biph10.samples, file = paste0(home,"/data/biph10_samples_1025.RData"))
#biph10.samples <- as.matrix(biph10.mcmcc$mvSamples)

# Save samples
#saveRDS(biph10.samples, file = paste0(home,"/data/biph10_samples_20000_1017.RData"))


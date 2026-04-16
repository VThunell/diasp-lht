# Biphasic model 11:
# simplify biph 9

# Load libraries, install if needed
pkgs <- c("tidyverse", "nimble")

if(length(setdiff(pkgs,rownames(installed.packages()))) > 0){
  install.packages(setdiff(pkgs, rownames(installed.packages())), dependencies = T)
}

invisible(lapply(pkgs, library, character.only = T))

## Biphasic 11
biph0.code <- nimbleCode({
  
  # likelihood
  for(i in 1:nobs){
    length_mm[i] ~ dnorm(l.mu[su[i],smo.age[i],coh[i],age[i]], sd = sig.l)
  }
      
  for(j in 1:nsu) {
    for(k in 1:nK) {
      for(l in 1:ncoh) {

       l.vec[1:max.age[l]] <- grow_one(lb.mu[j], g[j], log_p.year, j, l,
                                max.age[l], max_age_global, k, nsu)

        # assign only up to cohort-specific max.age
        for(m in 1:max.age[l]) {
          l.mu[j,k,l,m] <- l.vec[m]
        }
      }
    }
  }
  
  # for(j in 1:nsu) {
  #   for(k in 1:nK) {
  #     for(l in 1:ncoh) {
  #       #max_age_valid[j,k,l] <- min(max.age + 1, NY - l + 1)
  #       #l.mu[j,k,l,1:max_age_valid[j,k,l]] <- grow_one(lb.mu[j],g[j],par,j,l,max.age,NY,k)
  #       # run recursion for entire fixed-length vector
  #       l.vec <- grow_one(lb.mu[j], g[j], par, j, l,
  #                         max.age.coh[l], max_age_global, k)
  #       
  #       # assign only valid ages
  #       for(m in 1:max.age.coh[l]) {
  #         l.mu[j,k,l,m] <- l.vec[m]
  #         
  #         #l.mu[j,k,l,1] <- lb.mu[j] + exp(g[j])
  #         #for(m in 2:(max.age + 1)) {
  #         #for(m in 2:(max.age + 1)) {
  #         #l.mu[j,k,l,m] <- l.mu[j,k,l,(m-1)]
  #         #is_pre_smolt <- step(k - m)   # 1 if a < k else 0
  #         # pre-smolt linear increment
  #         #Lin_growth <- l.mu[j,k,l,m-1] + exp(g[j])
  #         # VB post-smolt: careful: we need an appropriate par indexing -> example:
  #         #Linf_now <- exp(par[t + (a-1), s, 1])   # ensure t+(a-1) mapping valid
  #         #K_now    <- exp(par[t + (a-1), s, 2])
  #         #VB_pred  <- Linf_now - (Linf_now - l.mu[s,k,t,a-1]) * exp(-K_now * dt)
  #         #VB_pred  <- l.mu[j,k,l,m-1]*1.5
  #         
  #         # l.mu[j,k,l,m] <-
  #         #   (l.mu[j,k,l,m-1] + step(k - m)*exp(g[j])) +
  #         # (1 - step(k - m)) * l.mu[j,k,l,m-1]*1.1#*test[ncoh+m]
  #         
  #         # spatiotemporal rw 
  #         # l.mu[j,k,l,m] <-
  #         #    l.mu[j,k,l,m-1] + step(k - m)*exp(g[j]) +
  #         #   (1-step(k - m)) * (((exp(par[l+(m-1),j,1])) - l.mu[j,k,l,m-1])*(1 - exp(-exp(par[(l+(m-1)),j,2]))))
  #         
  #         # spatial only
  #         # l.mu[j,k,l,m] <- 
  #         #   l.mu[j,k,l,m-1] + step(k - m)*exp(g[j]) + 
  #         #   (1-step(k - m)) * (((exp(par[j,1])) - l.mu[j,k,l,m-1])*(1 - exp(-exp(par[j,2]))))
  #       }
  #     }
  #   }
  # }
 
  ### --- cell-level smolt-age composition ---
  # for(c in 1:Ncells) {
  #   p_k[c, 1:nK] ~ ddirch(alpha[1:nK])   # Dirichlet prior
  # }
  # 
  # ### --- Observations WITH known k (use direct mapping) ---
  # for(ii in 1:n_known){
  #   i <- idx_obs_known[ii]
  #   kobs <- k_obs[i]   # provided in data
  #   length_mm[i] ~ dnorm(l.mu[ su[i], kobs, hatch.year.f[i], age[i] ], sd = sig.l)
  #   # Optionally inform p_k by the observed k: treat k_obs as data:
  #   # k_obs[i] ~ dcat(p_k[cell[i], 1:nK])   # NIMBLE treats this as likelihood for p_k
  # }
  # 
  # ### --- Observations WITH missing k: exact mixture via zeros-trick ---
  # # logSumExp nimbleFunction must be defined in R and available
  # for(ii in 1:n_miss){
  #   i <- idx_obs_missing[ii]
  #   # build log-component vector
  #   for(k in 1:nK){
  #     logComp[i,k] <- log(p_k[cell[i], k]) +
  #       dnorm(length_mm[i],
  #             mean = l.mu[ su[i], k, hatch.year.f[i], age[i] ],
  #             sd = sig.l,
  #             log = TRUE)
  #   }
  #   logLik_mix[i] <- logSumExp(logComp[i, 1:nK])
  #   zeros[i] ~ dpois(phi[i])    # zeros trick
  #   phi[i] <- -logLik_mix[i] + C   # C is a large constant (provide as data)
  # }
  # 
  # 
  # priors for smolt age, lb.mu, g
  for(i in 1:nsu){
    g[i] ~ dnorm(g.lmean, sd = 1) # rw first year
    #smp[1:nK,i] ~ ddirch(alpha[1:nK])
    lb.mu[i] ~ dunif(12, 22) # ~15 to 20 mm, DOI: 10.1111/j.1095-8649.2009.02497.x  sd is guesstimate
  }
  
  #par[y,su,p]
  # for(i in 1:nyear) {
  #   for(j in 1:nsu) {
  #     par[i,j,1] <- log_p.year[i,j]
  #     par[i,j,2] <- log_p.year[i,(nsu+j)]
  #   }
  # }

  # LKJ prior on correlation matrix for vb-pars
  Ustar[1:npars,1:npars] ~ dlkj_corr_cholesky(1.3, npars) # eta = 1.3
  U[1:npars,1:npars] <- uppertri_mult_diag(Ustar[1:npars, 1:npars], sig.par[1:npars])

  # first year of vb-par rw, npars = nsu*number of vb-pars
  for(k in 1:npars){
     log_p.year[1,k] ~ dnorm(mu.par[k], sd = 1) # sd value on the mnorm, logscale ???
  }
  # rw for vb-par rw years > 1
  for(j in 1:nyear){
    log_p.year[(j+1), 1:npars] ~ dmnorm(log_p.year[j,1:npars], cholesky = U[1:npars, 1:npars], prec_param = 0)
  }

  for(i in 1:nsu){
    mu.par[i] ~ dnorm(l.lmean, sd = 1)
    sig.par[i] ~ dlnorm(log(l.lsd), sdlog = 0.1)
    }
  for(i in (nsu+1):npars){
    mu.par[i] ~ dnorm(mean = k.lmean, sd = 0.5)
    sig.par[i] ~ dlnorm(log(k.lsd), sdlog = 0.1)
  }

  # LKJ prior on correlation matrix for vb-pars
  # Ustar[1:npars,1:npars] ~ dlkj_corr_cholesky(1.3, npars) # eta = 1.3
  # U[1:npars,1:npars] <- uppertri_mult_diag(Ustar[1:npars, 1:npars], sig.par[1:npars])
  # 
  # for(i in 1:nsu){
  #   par[i, 1:npars] ~ dmnorm(mu.par[su[i],1:npars], cholesky = U[1:npars, 1:npars], prec_param = 0)
  # }
  # for(i in 1:nsu){
  #   mu.par[i,1] ~ dnorm(l.lmean, sd = 1)
  #   mu.par[i,2] ~ dnorm(mean = k.lmean, sd = 0.5)
  # }
  # sig.par[1] ~ dlnorm(log(l.lsd), sdlog = 0.1)
  # sig.par[2] ~ dlnorm(log(k.lsd), sdlog = 0.1)
  # # g, k & linf values (g and linf from fb)
   g.lsd <- sqrt(log(1 + (5^2) / (50^2))) # sqrt of variance (= sd) of the lognormal distr. 50 and 10 is arbitrary atm.
   g.lmean <- log(50) - 0.5 * log(g.lsd)^2  # mean of the lognorm distr.
   l.lsd <- sqrt(log(1 + (278^2) / (1378^2)))
   l.lmean <- log(1378) - 0.5 * log(l.lsd)^2 # sd of the mean based on sd:s of LKJ
   k.lsd <- sqrt(log(1 + (0.28^2) / (0.43^2)))
   k.lmean <- log(0.43) - 0.5 * log(k.lsd)^2

  sig.l ~ dexp(1/150)
  
  })
  
    # for(j in 1:nsu) {
  #     for(k in 1:nK) {
  #     for(l in 1:nyear) {
  #       
  #       # deterministic 0+ length
  #       l.mu[j,k,l,1] <- lb.mu[j] + exp(g[j])
  #       
  #       for(m in 2:(max.age[j] + 1)) {
  #         
  #         # 1 = pre-smolt, 0 = post-smolt
  #         pre_smolt[j,k,l,m] <- step(k - m)
  #         post_smolt[j,k,l,m] <- 1 - pre_smolt[j,k,l,m]
  #         
  #         # Time-varying VB parameters
  #         Linf_tmp[j,k,l,m] <- exp(par[l + (m-1), j, 1])
  #         Kval_tmp[j,k,l,m] <- exp(par[l + (m-1), j, 2])
  #         
  #         # Growth increments
  #         incr_pre[j,k,l,m]  <- pre_smolt[j,k,l,m] * exp(g[j])
  #         incr_post[j,k,l,m] <- post_smolt[j,k,l,m] *
  #           ((Linf_tmp[j,k,l,m] -
  #               l.mu[j,k,l,m-1]) *
  #              (1 - exp(-Kval_tmp[j,k,l,m])))
  #         
  #         # Final Fabens update
  #         l.mu[j,k,l,m] <- l.mu[j,k,l,m-1] +
  #           incr_pre[j,k,l,m] +
  #           incr_post[j,k,l,m]
  #       }
  #     }
  #   }
  # }
  # # # l.mu[]: to retrieve last years length ests, we need to store lengths in an array of predictors
  # #   for(j in 1:nsu){
  # #     for(k in 1:nK){
  # #       for(l in 1:nyear){
  # #         
  # #         # length prediction 0+
  # #         l.mu[j,k,l,1] <- dnorm(lb.mu[j] + exp(g[j]), sd = 10)
  # #         
  # #         for(m in 2:(max.age[j]+1)){
  # #           # length prediction 1+ and older
  # #           l.mu[j,k,l,m] <- l.mu[j,k,l,m-1] + (step(k - (m-0.999))*exp(g[j]) + 
  # #                                             (1-step(k - (m-0.999)))*(((exp(par[l+(m-1),j,1])) -
  # #                                                                         l.mu[j,k,l,m-1])*(1 - exp(-(exp(par[(l+(m-1)),j,2]))))))
  # #         }
  # #       }
  # #     }
  # #   }
  # 
  # #par[gy,su,p]
  # for(i in 1:nyear) {
  #   for(j in 1:nsu) {
  #     par[i,j,1] <- log_p.year[i,j]     
  #     par[i,j,2] <- log_p.year[i,(nsu+j)]
  #   }
  # }
  # 
  # # estimate smo.age
  # for(l in 1:nobs){
  #   smo.age[l] ~ dcat(smp[1:nK,su[l]])
  # }
  # 
  # # Priors
  
  
  # LKJ prior on correlation matrix for vb-pars
  # Ustar[1:npars,1:npars] ~ dlkj_corr_cholesky(1.3, npars) # eta = 1.3
  # U[1:npars,1:npars] <- uppertri_mult_diag(Ustar[1:npars, 1:npars], sig.par[1:npars])
  # 
  # # first year of vb-par rw, npars = nsu*number of vb-pars
  # for(k in 1:npars){
  #    log_p.year[1,k] ~ dnorm(mu.par[k],10)
  # }
  # # rw for vb-par rw years > 1
  # for(j in minyear:maxyear){
  #   log_p.year[(j+1), 1:npars] ~ dmnorm(log_p.year[j,1:npars], cholesky = U[1:npars, 1:npars], prec_param = 0)
  # }
  # 
  # # priors for smolt age, lb.mu, g
  # for(i in 1:nsu){
  #   g[i] ~ dnorm(g.par[i], sd = 1) # rw first year
  #   smp[1:nK,i] ~ ddirch(alpha[1:nK])
  #   lb.mu[i] ~ dunif(12, 22) # ~15 to 20 mm, DOI: 10.1111/j.1095-8649.2009.02497.x  sd is guesstimate
  # }
  # 
  # for(i in 1:nsu){
  #   mu.par[i] ~ dnorm(l.lmean, sd = 1)
  #   sig.par[i] ~ dlnorm(log(l.lsd), sdlog = 0.1)  
  #   }
  # for(i in (nsu+1):npars){
  #   mu.par[i] ~ dnorm(mean = k.lmean, sd = 0.5)
  #   sig.par[i] ~ dlnorm(log(k.lsd), sdlog = 0.1)
  # }
  # 
  # # g, k & linf values (g and linf from fb)
  # g.lmean <- log(50) - 0.5 * log(g.lsd)^2  # mean of the lognorm distr.
  # l.lsd <- sqrt(log(1 + (278^2) / (1378^2)))
  # l.lmean <- log(1378) - 0.5 * log(l.lsd)^2 # sd of the mean based on sd:s of LKJ
  # k.lsd <- sqrt(log(1 + (0.28^2) / (0.43^2)))
  # k.lmean <- log(0.43) - 0.5 * log(k.lsd)^2  
  # 

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


grow_one <- nimbleFunction(
  run = function(lb = double(0),g = double(0),logp = double(2),j = double(0),l = double(0),
                 max_age_l = double(0),max_age_global = double(0),k = double(0),nsu = double(0)) {
    returnType(double(1))
    
    lout <- numeric(max_age_global)
    lout[1] <- lb + exp(g)

    for(m in 2:max_age_l) {
      idx <- l + (m - 1)

      Linf_log <- logp[idx, j]
      K_log    <- logp[idx, nsu + j]

      lout[m] <- lout[m-1] +
        nimStep(k - m) * exp(g) +
        (1 - nimStep(k - m)) *
        ((exp(Linf_log) - lout[m-1]) *
           (1 - exp(-exp(K_log))))
    } 
    return(lout)
    })

logSumExp <- nimbleFunction(
  run = function(x = double(1)) {
    returnType(double(0))
    maxx <- max(x)
    return(maxx + log(sum(exp(x - maxx))))
  }
)


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

set.seed(1)
data.biph <- sallaa %>%
  mutate(yday.may1 = if_else(is.na(date), 1, yday(date) - yday("2025-05-01")+1), 
         #doy.dec = if_else(yday.may1 > 0, yday.may1, yday.may1+365)/365,
         doy.dec = round(if_else(yday.may1 > 0, yday.may1, yday.may1+365)/365,1),
         tot.age.dec = if_else(yday.may1 > 0, tot.age+doy.dec, tot.age+doy.dec-1),
         smo.age = if_else(age.type == "both", juv.age, NA),
         sea.age = if_else(is.na(sea.age), 0, sea.age),
         growth.year = trunc(tot.age.dec)+1, # + 1 to index age 0 individuals
         hatch.year = year-tot.age,
         hatch.year.f = as.numeric(factor(hatch.year, levels = sort(unique(hatch.year)), labels = seq_along(sort(unique(hatch.year))))),
         year.f = as.integer(factor(year)),
         sex = as.integer(factor(sex)),
         su = as.integer(factor(spat.unit))) %>%
  #mutate(max.age = max(tot.age), .by = su) %>%
  filter(!(is.na(smo.age))) %>%
  #slice_sample(n = 50000) %>%
  mutate(ind.id = row_number()) 

nK = max(data.biph$smo.age %>% na.omit())
nsu = length(unique(data.biph$su))
#max.age = data.biph %>% summarise(max.age = max(tot.age), .by = hatch.year.f) %>% pull(max.age) +1


#max.age = data.biph %>% summarise(max.age = max(tot.age)) %>% pull(max.age)
# maxyear = max(data.biph$hatch.year.f + data.biph$growth.year)
# minyear = min(data.biph$hatch.year.f)
ncoh <- length(unique(data.biph$hatch.year.f))
cohort_max_age <- integer(ncoh)
for(l in 1:ncoh){
  # cohort 'l' corresponds to cohort index l in your coh vector
  cohort_year <- sort(unique(data.biph$hatch.year.f))[l]   # cohort numeric value
  cohort_max_age[l] <- max( data.biph$tot.age[data.biph$hatch.year.f == cohort_year] ) + 1
}

nyear = ncoh + cohort_max_age[ncoh]
npars = 2*nsu

consts = list(nobs = nrow(data.biph),
              nK = nK,
#               alpha = c(1, 2, 3, 2, 1, 1)[1:nK],
              nsu = nsu,
              npars = 2*nsu, # n parameters in the LKJ dmnorm * number of spatial units 
              age = data.biph$tot.age,
#               #growth.year = data.biph$growth.year,
#               #doy.dec = data.biph$doy.dec,
              smo.age = data.biph$smo.age,
#               #ind.id = data.biph$ind.id,
              coh = data.biph$hatch.year.f,
#               #shfa = fec.dpars$shfa,
#               #rafa = fec.dpars$rafa,
#               #mufb = fec.dpars$mufb,
#               #sdfb = fec.dpars$sdfb,
#               #sex = data.biph$sex,
              max.age = cohort_max_age,
              nyear = nyear,
              #maxyear = maxyear,
              #minyear = minyear,
              #nyear = length(minyear:maxyear),
              ncoh = ncoh,
              max_age_global = max(cohort_max_age),
              #coh.mage = cohort_mage
              su = data.biph$su)

# initialize NA smo.ages to allow for NA in smo.age in l.mu array
# init.smo.age <- data.biph$smo.age
# na.smo.age <- which(is.na(data.biph$smo.age))
# init.smo.age[na.smo.age] <- sample(1:nK, length(na.smo.age), replace = TRUE)

inits <- function(){
  list(sig.l = rexp(1,1/250),
       lb.mu = runif(consts$nsu,13,19),
       Ustar  = diag(consts$npars),
       #Ugstar = diag(consts$nsu),
       #smo.age = init.smo.age,
       sig.par = c(rlnorm(consts$nsu,log(0.2),sdlog = 0.1),
                   rlnorm(consts$nsu,log(0.6),sdlog = 0.1)),
       mu.par = c(rnorm(consts$nsu,7,0.2),
                  rnorm(consts$nsu,-1,0.6)),
       mu.par = c(rnorm(consts$nsu,0,0.1),
                  rnorm(consts$nsu,0,0.1)),
       #g.sig.par = rlnorm(consts$nsu,0,0.1),
       g.par = rnorm(consts$nsu,4,0.2)#,
      )
  }

#grow_one <- compileNimble(grow_one) 

# build model
biph0.model <- nimbleModel(biph0.code,
                     constants = consts,
                     #inits = inits(),
                     data = data.biph %>% select(length_mm),
                     dimensions = list(
                       log_p.year = c(nyear+1, npars),
                       l.mu = c(nsu, nK, ncoh, 13)))

stochNodes <- biph0.model$getNodeNames(stochOnly = TRUE)
print(head(stochNodes, 100))

# 2) check which stochastic nodes have no data-dependent children (i.e. never affect observed nodes)
noDep <- sapply(stochNodes, function(n) {
  deps <- biph0.model$getDependencies(n, downstream = TRUE)   # all nodes dependent on n
  # if none of deps are in the observed data node set, mark
  anyObserved <- any(grepl("length_mm", deps))  # adapt to your obs node names
  !anyObserved
})
problemNodes <- stochNodes[noDep]
problemNodes

biph0.model$simulate()
biph0.model$calculate()
biph0.model$par
#Sys.time() - t
# for 10000 individs, it takes
# for model 9: 49 secs
# then removing temporal process on g: 43 secs
# then rounding plus growth: 43 secs
# then removing plus growth esr: 37 secs

# identify nodes to sample 
dataNodes <- biph0.model$getNodeNames(dataOnly = TRUE)
deterNodes <- biph0.model$getNodeNames(determOnly = TRUE)
parentNodes <- biph0.model$getParents(dataNodes, stochOnly = TRUE) #all of these should be added to monitor below to recreate other model variables...
stnodes <- biph0.model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
allvars <- biph0.model$getVarNames(nodes = stnodes)
mvars <- allvars[!(grepl("lifted",allvars))]  

# calculate vars
vs <- mvars
#des <- parentNodes[1:100]
for(i in 1:length(vs)){
  print(paste0(vs[i]," ",biph0.model$calculate(vs[i]) ))
}

exp(biph0.model$g)
biph0.model$lb.mu
biph0.model$sig.l
biph0.model$calculate("l.mu[1,2,1,1]")
biph0.model$l.mu[1,2,22,]
# for(i in 1:length(deterNodes)){
#   if(is.na(biph11.model$calculate(deterNodes[i])))
#      {print(paste0(deterNodes[i]," ",biph11.model$calculate(deterNodes[i])))}
# }
# 
# # simulate vars
# for(i in 1:length(dataNodes)){
#   biph11.model$simulate(mvars[i])
#   print(paste0(mvars[i],"",values(biph11.model, mvars[i])))
# }
# for(i in 1:length(deterNodes)){
#   biph11.model$simulate(deterNodes[i])
#   if(sum(is.na(values(biph11.model, deterNodes[i])))>0)
#   {print(deterNodes[i])}
# }

# compile model
biph0.c <- compileNimble(biph0.model)

# configure and build mcmc
#biph11.confmcmc <- configureMCMC(biph11.c, monitors = c("sig.l", "Ustar","lb.mu","smp","g","sig.par","smo.age","mu.par","par", "l.mu", "log_p.year"),
                               # useConjugacy = FALSE)
biph0.confmcmc <- configureMCMC(biph0.c, monitors = c("sig.l", "Ustar","lb.mu","g","log_p.year","l.mu"),#,"bL","bk"),
                                useConjugacy = FALSE)

biph0.mcmc <- buildMCMC(biph0.confmcmc, project = biph0.model)

# compile mcmc
biph0.mcmcc <- compileNimble(biph0.mcmc, project = biph0.model)

# MCMC samples
#t <- Sys.time()
#biph0.samples <- biph0.mcmcc$run(200,nburnin=150)
biph0.samples <- runMCMC(biph0.mcmcc, niter = 100000, nburnin = 75000, thin=10)
#biph11.samples <- runMCMC(biph11.mcmcc, niter = 200000, nburnin = 185000, thin=10)#, WAIC=TRUE) 
#Sys.time() - t
#gc()

#saveRDS(biph0.samples, file = paste0(home,"/data/biph0_samples_1025.RData"))
#biph11.samples <- as.matrix(biph11.mcmcc$mvSamples)

# Save samples
#saveRDS(biph11.samples, file = paste0(home,"/data/biph11_samples_20000_1017.RData"))


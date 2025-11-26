# Biphasic model 0.4:
# remove rw on Linf9

# Load libraries, install if needed
pkgs <- c("tidyverse", "nimble")

if(length(setdiff(pkgs,rownames(installed.packages()))) > 0){
  install.packages(setdiff(pkgs, rownames(installed.packages())), dependencies = T)
}

invisible(lapply(pkgs, library, character.only = T))

## model
biph04.code <- nimbleCode({
  
  # likelihood
  for(i in 1:nobs){
    length_mm[i] ~ dnorm(l.mu[su[i],smo.age[i],coh[i],g.year[i]], sd = sig.l)
  }
   
  for(j in 1:nsu) {
    for(k in 1:nK) {
      for(l in 1:ncoh) {
        
        l.mu[j,k,l,1] <- lb.mu[j] + exp(g[j])
        
        for(m in 2:coh.maxgage[l]) {
          
          l.mu[j,k,l,m] <- l.mu[j,k,l,m-1] + step(k - m)*exp(g[j]) +
            (1-step(k - m)) * ((exp(linf[j]) - l.mu[j,k,l,m-1])*(1 - exp(-exp(k.year[(l+(m-1)),j]))))
        }
      }
    }
  }
  
  # priors for lb.mu, g, linf and rw k.par
  for(i in 1:nsu){
    lb.mu[i] ~ dunif(12, 22) # ~15 to 20 mm, DOI: 10.1111/j.1095-8649.2009.02497.x  sd is guesstimate
    g[i] ~ dnorm(g.lmean, sd = 1) # rw first year
    linf[i] ~ dnorm(l.lmean, sd = l.lsd)
    k.par[i] ~ dnorm(mean = k.lmean, sd = 0.5)
    sig.par[i] ~ dlnorm(log(k.lsd), sdlog = 0.1)
  }
  
  # LKJ prior on correlation matrix for vb-pars
  Ustar[1:npars,1:npars] ~ dlkj_corr_cholesky(1.3, npars) # eta = 1.3
  U[1:npars,1:npars] <- uppertri_mult_diag(Ustar[1:nsu, 1:nsu], sig.par[1:nsu])

  # first year of rw for k over su:s
  for(k in 1:nsu){
     k.year[1,k] ~ dnorm(k.par[k], sd = 1) # sd value on the mnorm, logscale ???
  }
  # rw for vb-par rw years > 1
  for(j in 1:(nyear-1)){
    k.year[(j+1), 1:npars] ~ dmnorm(k.year[j,1:npars], cholesky = U[1:npars, 1:npars], prec_param = 0)
  }

  # # g, k & linf values (g and linf from fb)
   g.lsd <- sqrt(log(1 + (5^2) / (50^2))) # sqrt of variance (= sd) of the lognormal distr. 50 and 10 is arbitrary atm.
   g.lmean <- log(50) - 0.5 * log(g.lsd)^2  # mean of the lognorm distr.
   l.lsd <- sqrt(log(1 + (278^2) / (1378^2)))
   l.lmean <- log(1378) - 0.5 * log(l.lsd)^2 # sd of the mean based on sd:s of LKJ
   k.lsd <- sqrt(log(1 + (0.28^2) / (0.43^2)))
   k.lmean <- log(0.43) - 0.5 * log(k.lsd)^2

  sig.l ~ dexp(1/150)
  
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

set.seed(1)
data.biph <- sallaa %>%
  filter(!is.na(juv.age)) %>% # removes age.type = "sea.only"
  mutate(yday.may1 = if_else(is.na(date), 1, yday(date) - yday("2025-05-01")+1), 
         doy.dec = if_else(yday.may1 > 0, yday.may1, yday.may1+365)/365,
         #doy.dec = round(if_else(yday.may1 > 0, yday.may1, yday.may1+365)/365,1),
         tot.age.dec = if_else(yday.may1 > 0, tot.age+doy.dec, tot.age+doy.dec-1),
         smo.age = if_else(age.type == "both", juv.age, NA),
         sea.age = if_else(is.na(sea.age), 0, sea.age),
         g.year = trunc(tot.age.dec)+1, # + 1 to index age 0 individuals
         hatch.year = year-tot.age,
         hatch.year.f = as.numeric(factor(hatch.year, levels = sort(unique(hatch.year)), labels = seq_along(sort(unique(hatch.year))))),
         year.f = as.integer(factor(year)),
         sex = as.integer(factor(sex)),
         su = as.integer(factor(spat.unit))) %>%
  filter(!(is.na(smo.age))) %>%
  #slice_sample(n = 50000) %>%
  mutate(ind.id = row_number()) 

nK = max(data.biph$smo.age %>% na.omit())
nsu = length(unique(data.biph$su))
ncoh <- length(unique(data.biph$hatch.year.f))
coh.maxgage <- data.biph %>% summarise(max.age = max(g.year), .by = hatch.year.f) %>% arrange(hatch.year.f) %>% pull(max.age)
nyear = ncoh + coh.maxgage[ncoh]
npars = nsu

consts = list(nobs = nrow(data.biph),
              nK = nK,
#               alpha = c(1, 2, 3, 2, 1, 1)[1:nK],
              nsu = nsu,
              npars = nsu, # n parameters in the LKJ dmnorm * number of spatial units 
              #age = data.biph$tot.age,
              g.year = data.biph$g.year,
#               #doy.dec = data.biph$doy.dec,
              smo.age = data.biph$smo.age,
              coh = data.biph$hatch.year.f,
              coh.maxgage = coh.maxgage,
              nyear = nyear,
              ncoh = ncoh,
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

# build model
biph04.model <- nimbleModel(biph04.code,
                     constants = consts,
                     #inits = inits(),
                     data = data.biph %>% select(length_mm))

stochNodes <- biph04.model$getNodeNames(stochOnly = TRUE)
print(head(stochNodes, 100))

# which stochastic nodes have no data-dependent children nodes
noDep <- sapply(stochNodes, function(n) {
  deps <- biph04.model$getDependencies(n, downstream = TRUE)   # all nodes dependent on n
  # if none of deps are in the observed data node set, mark
  anyObserved <- any(grepl("length_mm", deps))  # adapt to your obs node names
  !anyObserved
})
problemNodes <- stochNodes[noDep]
problemNodes

biph04.model$simulate()
biph04.model$calculate()

# identify nodes to sample 
dataNodes <- biph04.model$getNodeNames(dataOnly = TRUE)
deterNodes <- biph04.model$getNodeNames(determOnly = TRUE)
parentNodes <- biph04.model$getParents(dataNodes, stochOnly = TRUE) #all of these should be added to monitor below to recreate other model variables...
stnodes <- biph04.model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
allvars <- biph04.model$getVarNames(nodes = stnodes)
mvars <- allvars[!(grepl("lifted",allvars))]  

# calculate vars
vs <- mvars
for(i in 1:length(vs)){
  print(paste0(vs[i]," ",biph04.model$calculate(vs[i]) ))
}

exp(biph04.model$g)
biph04.model$lb.mu
biph04.model$sig.l
biph04.model$calculate("l.mu[1,2,1,1]")
biph04.model$l.mu[1,2,22,]
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
biph04.c <- compileNimble(biph04.model)

# configure and build mcmc
#biph11.confmcmc <- configureMCMC(biph11.c, monitors = c("sig.l", "Ustar","lb.mu","smp","g","sig.par","smo.age","mu.par","par", "l.mu", "log_p.year"),
                               # useConjugacy = FALSE)
biph04.confmcmc <- configureMCMC(biph04.c, monitors = c("sig.l", "Ustar","lb.mu","g","k.year","linf","l.mu"),#,"bL","bk"),
                                useConjugacy = FALSE)

biph04.mcmc <- buildMCMC(biph04.confmcmc, project = biph04.model)

# compile mcmc
biph04.mcmcc <- compileNimble(biph04.mcmc, project = biph04.model)

# MCMC samples
#t <- Sys.time()
#biph04.samples <- biph04.mcmcc$run(200,nburnin=150)
biph04.samples <- runMCMC(biph04.mcmcc, niter = 4000000, nburnin = 3975000, nchains = 2, thin=10, samplesAsCodaMCMC = TRUE)
#biph11.samples <- runMCMC(biph11.mcmcc, niter = 200000, nburnin = 185000, thin=10)#, WAIC=TRUE) 
#Sys.time() - t
#gc()

saveRDS(biph04.samples, file = paste0(home,"/data/biph04_samples_1123.RData"))
#biph11.samples <- as.matrix(biph11.mcmcc$mvSamples)

# Save samples
#saveRDS(biph11.samples, file = paste0(home,"/data/biph11_samples_20000_1017.RData"))


# Biphasic model 0.4:
# remove rw on Linf9

# Load libraries, install if needed
pkgs <- c("tidyverse", "nimble","nimbleHMC")

if(length(setdiff(pkgs,rownames(installed.packages()))) > 0){
  install.packages(setdiff(pkgs, rownames(installed.packages())), dependencies = T)
}

invisible(lapply(pkgs, library, character.only = T))

## model
biph05.code <- nimbleCode({
  
  # likelihood both age and sea age only
  for(i in 1:nobs){
    length_mm[i] ~ dnorm(l.mu[su[i],smo.age[i],coh[i],g.year[i]], sd = sig.l)
    #length_mm[i] ~ dnorm(l.mu[su[i],smo.age[i],coh[i],g.year[i],sex[i]], sd = sig.l)
    
    # estimate smo.age
    #smo.age[i] ~ dcat(smp[1:nK,su[i]])
  }
  
  for(j in 1:nsu) {
    for(k in 1:nK) {
      for(l in 1:ncoh) {
        #for(m in 1:2) { # nsex, if I estimate sex below, do I get away with not taking into account Na sex
        
        l.mu[j,k,l,1] <- lb.mu[j] + exp(g[j])
        
        for(m in 2:coh.maxgage[l]) {
          
          l.mu[j,k,l,m] <- l.mu[j,k,l,m-1] + step(k - m)*exp(g[j]) +
            (1-step(k - m)) * ((exp(linf[j]) - l.mu[j,k,l,m-1])*(1 - exp(-exp(k.year[(l+(m-1)),j]))))
          
          # l.mu[j,k,l,m] <- l.mu[j,k,l,m-1] + step(k - m)*exp(g[j]) +
          #   (1-step(k - m)) * ((linf[j,m] - l.mu[j,k,l,m-1])*(1 - exp(-exp(k.year[(l+(m-1)),j]))))
          # linf[j,m] <- exp(linf[j,m]) + bL[m]
          # k <- exp(linf[j,m]) + bK[m]
        }
      }
    }
  }
  
  # likelihood juvenile age only
  # for(i in 1:njobs){
  #     length_mm_j[i] ~ dnorm(l.mu.j[i], sd = sig.lj)
  #     #length_mm[i] ~ dnorm(l.mu.j[su[i],coh[i],g.year[i]], sd = sig.l)
  # #}
  #      l.mu.j[i] <- lb.mu[suj[i]] + exp(g[suj[i]])*gj.year[i]
  # # for(j in 1:nsu) {
  # #    for(l in 1:ncoh) {
  # #       l.mu.j[j,l,1] <- lb.mu[j] + exp(g[j])
  # #         for(m in 2:coh.maxgage[l]) {
  # #         l.mu[j,l,m] <- l.mu[j,l,m-1] + exp(g[j])
  # #       }
  # #     }
  #   }
  
  # sex likelihood
  # for(j in 1:nobs){
  #   sex[j]~dbern(psex[1])
  # }
  # psex[1] ~ dbeta(1,1)  #prop of males (prob. to be 1)
  # psex[2] <- 1-psex[1]
  
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
    #smp[1:nK,k] ~ ddirch(alpha[1:nK])
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
  sig.lj ~ dexp(1/50)
  
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
  mutate(yday.may1 = if_else(is.na(date), 1, yday(date) - yday("2025-05-01")+1), 
         doy.dec = if_else(yday.may1 > 0, yday.may1, yday.may1+365)/365,
         tot.age.dec = if_else(yday.may1 > 0, tot.age+doy.dec, tot.age+doy.dec-1),
         smo.age = if_else(age.type == "both", juv.age, NA),
         sea.age = if_else(is.na(sea.age), 0, sea.age),
         g.year = trunc(tot.age.dec)+1, # + 1 to index age 0 individuals
         hatch.year.f = as.numeric(factor(hatch.year, levels = sort(unique(hatch.year)), labels = seq_along(sort(unique(hatch.year))))),
         year.f = as.integer(factor(year)),
         sex = as.integer(factor(sex)),
         su = as.integer(factor(spat.unit))) %>%
  #filter(!is.na(smo.age))
  mutate(ind.id = row_number())

data.b <- data.biph %>%
  filter(!age.type == "juve.only") %>%
  filter(!is.na(smo.age))

data.j <- data.biph %>%
  filter(age.type == "juve.only") %>%
  rename(length_mm_j = length_mm)

nK = max(data.b$smo.age %>% na.omit())
nsu = length(unique(data.b$su))
ncoh <- length(unique(data.b$hatch.year.f))
coh.maxgage <- data.b %>% summarise(max.age = max(g.year), .by = hatch.year.f) %>% arrange(hatch.year.f) %>% pull(max.age)
nyear = ncoh + coh.maxgage[ncoh]
npars = nsu


consts = list(nobs = nrow(data.b),
              njobs = nrow(data.j),
              nK = nK,
#               alpha = c(1, 2, 3, 2, 1, 1)[1:nK],
              nsu = nsu,
              alpha = c(1, 2, 3, 2, 1, 1)[1:nK],
              npars = nsu, # n parameters in the LKJ dmnorm * number of spatial units 
              #age = data.biph$tot.age,
              g.year = data.b$g.year,
              gj.year = data.j$g.year,
#               #doy.dec = data.biph$doy.dec,
              #smo.age = data.b$smo.age,
              coh = data.b$hatch.year.f,
              coh.maxgage = coh.maxgage,
              nyear = nyear,
              ncoh = ncoh,
              su = data.b$su,
              suj = data.j$su)

# initialize NA smo.ages to allow for NA in smo.age in l.mu array
# init.smo.age <- data.b$smo.age
# na.smo.age <- which(is.na(data.b$smo.age))
# init.smo.age[na.smo.age] <- sample(1:nK, length(na.smo.age), replace = TRUE)

inits <- function(){
  list(sig.l = rexp(1,1/250),
       lb.mu = runif(consts$nsu,13,19),
       Ustar  = diag(consts$npars),
       #Ugstar = diag(consts$nsu),
       #smo.age = init.smo.age,
       sig.par = c(rlnorm(consts$nsu,log(0.2),sdlog = 0.1),
                   rlnorm(consts$nsu,log(0.6),sdlog = 0.1))
       #linf = rnorm(1,l.lmean, sd = l.lsd),
       #k.par = rnorm(1, mean = k.lmean, sd = 0.5)
       # mu.par = c(rnorm(consts$nsu,7,0.2),
       #            rnorm(consts$nsu,-1,0.6)),
       # mu.par = c(rnorm(consts$nsu,0,0.1),
       #            rnorm(consts$nsu,0,0.1)),
       #g.sig.par = rlnorm(consts$nsu,0,0.1),
       #g.par = rnorm(consts$nsu,4,0.2)#,
      )
  }

# build model
biph05.model <- nimbleModel(biph05.code,
                            constants = consts,
                            #inits = inits(),
                            #data = data.biph %>% select(length_mm))
                            data = list(data.b %>% select(length_mm, smo.age)))
                                 #data.j %>% select(length_mm_j)
                                 #))
                                 
stochNodes <- biph05.model$getNodeNames(stochOnly = TRUE)
print(head(stochNodes, 100))

# which stochastic nodes have no data-dependent children nodes
noDep <- sapply(stochNodes, function(n) {
  deps <- biph05.model$getDependencies(n, downstream = TRUE)   # all nodes dependent on n
  # if none of deps are in the observed data node set, mark
  anyObserved <- any(grepl("length_mm", deps))  # adapt to your obs node names
  !anyObserved
})
problemNodes <- stochNodes[noDep]
problemNodes

biph05.model$simulate()
biph05.model$calculate()

# identify nodes to sample 
dataNodes <- biph05.model$getNodeNames(dataOnly = TRUE)
deterNodes <- biph05.model$getNodeNames(determOnly = TRUE)
parentNodes <- biph05.model$getParents(dataNodes, stochOnly = TRUE) #all of these should be added to monitor below to recreate other model variables...
stnodes <- biph05.model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
allvars <- biph05.model$getVarNames(nodes = stnodes)
mvars <- allvars[!(grepl("lifted",allvars))]  

# calculate vars
vs <- mvars
for(i in 1:length(vs)){
  print(paste0(vs[i]," ",biph05.model$calculate(vs[i]) ))
}

exp(biph05.model$g)
biph05.model$lb.mu
biph05.model$sig.l
biph05.model$calculate("l.mu[1,2,1,1]")
biph05.model$l.mu[1,2,22,]
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
biph05.c <- compileNimble(biph05.model)

# configure and build mcmc
#biph05.confmcmc <- configureHMC(biph05.c, monitors = c("sig.l", "Ustar","lb.mu","g","k.year","linf","l.mu","l.mu.j"),#,"bL","bk"),
biph05.confmcmc <- configureMCMC(biph05.c, monitors = c("sig.l", "Ustar","lb.mu","g","k.year","linf","l.mu","l.mu.j"),#,"bL","bk"),
                                useConjugacy = FALSE)

biph05.mcmc <- buildMCMC(biph05.confmcmc, project = biph05.model)

# compile mcmc
biph05.mcmcc <- compileNimble(biph05.mcmc, project = biph05.model)

# MCMC samples
#t <- Sys.time()
biph05.samples <- runMCMC(biph05.mcmcc, niter = 3000000, nburnin = 297500, nchains = 2, thin=10, samplesAsCodaMCMC = TRUE)
#biph05.samples <- runMCMC(biph05.mcmcc, niter = 1000, nburnin = 975, nchains = 2, thin=10, samplesAsCodaMCMC = TRUE)
#Sys.time() - t
#gc()

saveRDS(biph05.samples, file = paste0(home,"/data/biph05_samples_1123.RData"))
#biph11.samples <- as.matrix(biph11.mcmcc$mvSamples)

# Save samples
#saveRDS(biph11.samples, file = paste0(home,"/data/biph11_samples_20000_1017.RData"))


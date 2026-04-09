# DIASPARA WP2.2 LHT models - Viktor Thunell
# Salmon length at age model and mcmc
# 

slaa.code <- nimbleCode({
  
  # likelihood both age
  
  for(i in 1:nobs){
    length_mm[i] ~ dnorm(l.mu[su[i],smo.age[i],coh[i],sex[i],g.year[i]], sd = sig.l)
  }

  for(j in 1:nsu) {
    for(k in 1:nsmo) {
      for(l in 1:ncoh) {
        for(m in 1:2) {

            l.mu[j,k,l,m,1] <- lb.mu[j] + ((exp(linf.j[j]) - lb.mu[j])*(1 - exp(-exp(k.j[j]))))

            for(n in 2:coh.maxgage[l]) {

              l.mu[j,k,l,m,n] <- l.mu[j,k,l,m,n-1] +
               step(k - n)*((exp(linf.j[j]) - l.mu[j,k,l,m,n-1])*(1 - exp(-exp(k.j[j])))) + #freshwater
                (1-step(k - n)) * ((exp(linf[j,m]) - l.mu[j,k,l,m,n-1])*(1 - exp(-exp(k.year[j,m,(l+(n-1))])))) #marine
          
          # l.mu[j,k,l,m,1] <- lb.mu[j] + exp(g[j])
          # 
          # for(n in 2:coh.maxgage[l]) {
          # 
          #   l.mu[j,k,l,m,n] <- l.mu[j,k,l,m,n-1] + step(k - n)*exp(g[j]) +
          #     (1-step(k - n)) * ((exp(linf[j,m]) - l.mu[j,k,l,m,n-1])*(1 - exp(-exp(k.year[j,m,(l+(n-1))]))))

          }
        }
      }
    }
  }
  
  # likelihood juvenile age only
  for(ii in 1:njobs){
    length_mm_j[ii] ~ dnorm(l.mu.j[ii], sd = sig.lj)
    #l.mu.j[ii] <- lb.mu[suj[ii]] + exp(g[suj[ii]])*gj.year[ii]
    #l.mu.j[ii] <- lb.mu[suj[ii]] * b[suj[ii]]^gj.year[ii]
    l.mu.j[ii] <- exp(linf.j[suj[ii]])*(1 - exp(-exp(k.j[suj[ii]])*(gj.year[ii] - t0[suj[ii]])))  #freshwater 
  }
  
  # spatially varying priors
  for(i in 1:nsu){
    #lb.mu[i] ~ dunif(12, 20) # ~15 to 20 mm, DOI: 10.1111/j.1095-8649.2009.02497.x
    lb.mu[i] ~ dnorm(0, sd = 20) # ~15 to 20 mm, DOI: 10.1111/j.1095-8649.2009.02497.x
    t0[i] ~ dnorm(0, sd = 20) # ~15 to 20 mm, DOI: 10.1111/j.1095-8649.2009.02497.x
    #g[i] ~ dnorm(g.lmean, sd = 1) 
    b[i] ~ dnorm(0, sd = 1) 
    linf.j[i] ~ dnorm(lj.lmean, sd = lj.lsd)
    k.j[i] ~ dnorm(mean = kj.lmean, sd = 0.1) 
    linf[i,1] ~ dnorm(l.lmean, sd = l.lsd) # males
    linf[i,2] ~ dnorm(l.lmean, sd = l.lsd) # females
    k.par[i] ~ dnorm(mean = k.lmean, sd = 0.5)
    sig.par[i] ~ dlnorm(log(k.lsd), sdlog = 0.1)
  }
  
  #first year of rw for k over su:s (one for each sex)
  for(j in 1:nsu){
    k.year[j,1,1] ~ dnorm(k.par[j], sd = k.lsd)
    k.year[j,2,1] ~ dnorm(k.par[j], sd = k.lsd)
  }

  # transition of vB pars from year i to i+1 (1 for each sex)
  for(i in 1:(nyear-1)){
    k.year[1:nsu,1,(i+1)] ~ dmnorm(k.year[1:nsu,1,i], prec = tau_p[1:nsu,1:nsu])
    k.year[1:nsu,2,(i+1)] ~ dmnorm(k.year[1:nsu,2,i], prec = tau_p[1:nsu,1:nsu])
  }

  tau_p[1:nsu, 1:nsu] <- inverse(sigma_p[1:nsu, 1:nsu])

  for(k in 1:nsu){
    for(j in 1:nsu){
      sigma_p[k,j] <- Rnew[k,j] * sig.par[k] * sig.par[j]
    }
  }

  # Prior for correlation matrix (LKJ prior)
  phi[1]  <- eta + (nsu - 2)/2
  corY[1] ~ dbeta(phi[1], phi[1])
  r12   <- 2 * corY[1] - 1
  ##
  R[1,1]     <- 1
  R[1,2]     <- r12
  R[2,2]     <- sqrt(1 - r12^2)

  R[2:nsu,1]   <- 0

  for (m in 2:(nsu-1)) {
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
    for(jk in (m+1):nsu){
      R[jk,m] <- 0
    }
  }  #m

  Rnew[1:nsu,1:nsu] <- t(R[1:nsu,1:nsu]) %*% R[1:nsu,1:nsu]
  
  # # g, k & linf values (g and linf from fb)
  g.lsd <- sqrt(log(1 + (5^2) / (50^2))) 
  g.lmean <- log(50) - 0.5 * log(g.lsd)^2  
  l.lsd <- sqrt(log(1 + (278^2) / (1378^2)))
  l.lmean <- log(1378) - 0.5 * log(l.lsd)^2 
  k.lsd <- sqrt(log(1 + (0.28^2) / (0.43^2)))
  k.lmean <- log(0.43) - 0.5 * log(k.lsd)^2
  sig.l ~ dexp(1/150)
  sig.lj ~ dexp(1/20)
  
  # guesstimates of juv growth
  lj.lsd <- sqrt(log(1 + (50^2) / (150^2)))
  lj.lmean <- log(150) - 0.5 * log(lj.lsd)^2 
  kj.lsd <- sqrt(log(1 + (0.1^2) / (0.2^2)))
  kj.lmean <- log(0.2) - 0.5 * log(kj.lsd)^2
})

nsmo <- max(data.b$smo.age %>% na.omit())
nsu <- length(unique(data.b$su))
ncoh <- length(unique(data.b$hatch.year.f))
coh.maxgage <- data.b %>% summarise(max.age = max(g.year), .by = hatch.year.f) %>% arrange(hatch.year.f) %>% pull(max.age)
nyear <- ncoh + coh.maxgage[ncoh]

consts = list(nobs = nrow(data.b),
              njobs = nrow(data.j),
              nsmo = nsmo,
              nsu = nsu,
              g.year = data.b$g.year,
              gj.year = data.j$g.year,
              smo.age = data.b$smo.age,
              coh = data.b$hatch.year.f,
              coh.maxgage = coh.maxgage,
              nyear = nyear,
              ncoh = ncoh,
              su = data.b$su,
              sex = data.b$sex,
              suj = data.j$su,
              eta = 2)

inits <- function(){
  list(sig.l = rexp(1,1/150),
       sig.lj = rexp(1,1/20),
       lb.mu = runif(consts$nsu,13,19),
       g = runif(consts$nsu,0,1),
       #b = runif(consts$nsu,0,1),
       linf = matrix(rnorm(consts$nsu*2, log(1000), sd = 1),
                     nrow = consts$nsu, ncol = 2),
       linf.j = rnorm(consts$nsu, log(150), sd = 1),
       k.j = rnorm(consts$nsu, 0,1),
       k.par = runif(consts$nsu,0,1),
       k.year = array(runif(consts$nsu * 2 * consts$nyear,0,1), dim = c(consts$nsu, 2, consts$nyear)),
       sig.par = rlnorm(consts$nsu, log(0.2), sdlog = 0.1),
       corZ = matrix(rnorm((consts$nsu-1)*(consts$nsu-1), 0, 1),
                     nrow = consts$nsu-1, ncol = consts$nsu-1),
       corY = runif((consts$nsu-1), 0, 1)
  )}

# build model
slaa.model <- nimbleModel(slaa.code,
                          constants = consts,
                          inits = inits(),
                          data = c(data.b %>% select(length_mm),
                                   data.j %>% select(length_mm_j)),
                          buildDerivs = TRUE
                          )

slaa.model$simulate()
slaa.model$initializeInfo()
slaa.model$calculate()

dataNodes <- slaa.model$getNodeNames(dataOnly = TRUE)
deterNodes <- slaa.model$getNodeNames(determOnly = TRUE)
parentNodes <- slaa.model$getParents(dataNodes, stochOnly = TRUE) #all of these should be added to monitor below to recreate other model variables...
#biph9.model$getVarNames(nodes = parentNodes)
stnodes <- slaa.model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
allvars <- slaa.model$getVarNames(nodes = stnodes)
mvars <- allvars[!(grepl("lifted",allvars))]  

# calculate vars
for(i in 1:length(mvars)){
  print(paste0(mvars[i]," ",slaa.model$calculate(mvars[i]) ))
}


# compile model
slaa.c <- compileNimble(slaa.model)

monits = c("sig.l", "sig.lj", "R","lb.mu","k.j","k.year","linf","linf.j","l.mu","l.mu.j","corY","corZ")

# configure and build mcmc
slaa.confmcmc <- configureHMC(slaa.c, 
                              monitors = monits,
                              enableWAIC = TRUE,
                              useConjugacy = FALSE
)

slaa.mcmc <- buildMCMC(slaa.confmcmc, project = slaa.model)

# compile mcmc
slaa.mcmcc <- compileNimble(slaa.mcmc, project = slaa.model, resetFunctions = TRUE)

# MCMC samples
t <- Sys.time()
slaa.samples_b_sea <-runMCMC(slaa.mcmcc, niter = 12000, nburnin = 9000, nchains = 2, thin=1, samplesAsCodaMCMC = TRUE, WAIC = TRUE)
#slaa.samples_b <- runMCMC(slaa.mcmcc, niter = 5000, nburnin = 3500, nchains = 2, thin=1, samplesAsCodaMCMC = TRUE, WAIC = TRUE)
#slaa.samples_g <- runMCMC(slaa.mcmcc, niter = 5000, nburnin = 3500, nchains = 2, thin=2, samplesAsCodaMCMC = TRUE, WAIC = TRUE)
Sys.time() - t

# with vB juvenile growth
# slaa.samples_b_sea$WAIC #941174.1
# slaa.samples_b_sea$WAIC #807354.9 with dmnorm 00 1 lb.mu
saveRDS(slaa.samples_b_sea, file = paste0(home,"/data/slaa_samples_jvB3_",Sys.Date(),".RData"))
node.sub <- grep("linf.j\\[|lb.mu", colnames(slaa.samples_b_sea$samples[[1]]), value = TRUE) 
node.sub <- grep("linf.j\\[|k.j", colnames(slaa.samples_b_sea$samples[[1]]), value = TRUE) 
slaa.samples_b_sea$samples[, node.sub[], drop = FALSE] %>% gelman.diag() # all converged

# without adult phase
# slaa.samples_b$WAIC # 358786.1 
# node.sub <- grep("b\\[|lb.mu", colnames(slaa.samples_b$samples[[1]]), value = TRUE) 
# slaa.samples_b$samples[, node.sub[], drop = FALSE] %>% gelman.diag() # all converged
# 
# slaa.samples_g$WAIC # 320609.3
# node.sub <- grep("g\\[|lb.mu", colnames(slaa.samples_g$samples[[1]]), value = TRUE)
# slaa.samples_g$samples[, node.sub[], drop = FALSE] %>% gelman.diag() # all converged

laa.samples$samples %>%
  #gather_draws(linf.j[su], k.j[su]) %>%
  spread_draws(lb.mu[su]) %>%
  #gather_draws(linf.j[su], sep = ",") %>%
  median_qi() %>%
  mutate(v = exp(.value)) %>%
  print(n = 24)

slaa.samples_b_sea$samples %>%
  spread_draws(lb.mu[su]) %>%
  median_qi() #%>%
  expand_grid(age = 0:5) %>%
  ggplot() +
  geom_line(aes(x=age, y = lb.mu*b^age, color = factor(su)))

            
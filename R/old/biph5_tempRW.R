#Constants added to step functions
#Changes to MVN prior sigmas
#check nK
#currently using 1st 10,000 inds in data


library(dplyr)
library(nimble)
library(lubridate)


# data
sallaa <- readRDS(paste0(home,"/data/data-for-2-2/salmon-laa_2025-09-08.RData")) %>%
  rename(year = fi_year,
         site = sai_location,
         sea.age = sea_age_year,
         juv.age = juvenile_age_year,
         tot.age = tot_age_year,
         lat = fisa_y_4326,
         lon = fisa_x_4326) %>%
  #filter out non-aged individuals (~90000 individuals)
  filter(!(is.na(sea.age) & is.na(juv.age)))


# model
biph5.con.code <- nimbleCode({
  
  # likelihood
  for(i in 1:nobs){
    length_mm[i] ~ dnorm(l.mu[age.index[i], ind.id[i]], sd = sig.l)
  }
  
  # Estimate length at age
  for(i in 1:n.ind){
    #length prediction 0+
    l.mu[1,i] <- lb.mu[spat.unit[i]] +  (1-step(1-(age[i]+0.0001)))*exp(g[spat.unit[i], hatch.year.f[i]]) + step(1-(age[i]+0.0001))*exp(g[spat.unit[i], hatch.year.f[i]])*doy.dec[i] # age 0
    
    #step returns 1 if evaluates to >=0
    #0 if evaluates to <0
    
    
    #length prediction 1+ and older
    for(a in 2:(age.index[i]+1)){
      l.mu[a,i] <- l.mu[(a-1),i] + (1-step(a-(age[i]+0.0001)))*(step(smo.age[i] - (a-0.999))*exp(g[spat.unit[i], (hatch.year.f[i]+a-1)]) + (1-step(smo.age[i] - (a-0.999)))*((exp(par[spat.unit[i],(hatch.year.f[i]+a-1),1]) - l.mu[(a-1),i])*(1 - exp(-exp(par[spat.unit[i],(hatch.year.f[i]+a-1),2]))))) +
        
        step(a-(age[i]+0.0001))*(step(smo.age[i] - (a-0.999))*exp(g[spat.unit[i], (hatch.year.f[i]+a-1)])*doy.dec[i] + (1-step(smo.age[i] - (a-0.999)))*((exp(par[spat.unit[i],(hatch.year.f[i]+a-1),1]) - l.mu[(a-1),i])*(1 - exp(-exp(par[spat.unit[i],(hatch.year.f[i]+a-1),2])*doy.dec[i]))))  
      
      
      
    }
  }
  
  # estimate smo.age
  for(l in 1:n.ind){
    smo.age[l] ~ dcat(smp[1:nK,spat.unit[l]])
  }
  
  # Priors
  sig.l ~ dexp(1/100)
  
  # LKJ prior on correlation matrix, see NIMBLE manual p45.
  Ustar[1:npars,1:npars] ~ dlkj_corr_cholesky(1.3, npars) # eta = 1.3
  U[1:npars,1:npars] <- uppertri_mult_diag(Ustar[1:npars, 1:npars], sig.par[1:npars])
  
  # priors for smolt age and juvenile growth rate
  for(i in 1:nsu){
    smp[1:nK,i] ~ ddirch(alpha[1:nK])  
    lb.mu[i] ~ dnorm(17, sd = 5) # ~15 to 20 mm, https://doi.org/10.1111/j.1095-8649.2009.02497.x sd is guesstimate
    for(j in minyear:maxyear){
      par[i,j,1:npars] ~ dmnorm(mu.par[1:npars], cholesky = U[1:npars, 1:npars], prec_param = 0)
      g[i,j] ~ dnorm(mean = g.lmean, sd = g.lsd) 
    }
  }
  
  #not run!
  #  for(k in 1:npar){  #npar now nsu*number of parameters
  #    log_p.year[1,k] ~ dnorm(mu.par[k],10)
  #
  #  }
  #  
  #  # Multivariate random walk in the log scale
  #  for(y in minyear:maxyear) {
  #    log_p.year[(y+1), 1:npar] ~ dmnorm(log_p.year[y,1:npar], cholesky = U[1:npars, 1:npars], prec_param = 0)
  #  }
  #  
  
  g.lsd <- sqrt(log(1 + (10^2) / (50^2))) # sqrt of variance (= sd) of the lognormal distr. 50 and 10 is arbitrary atm.
  g.lmean <- log(50) - 0.5 * g.lsd^2  # mean of the lognorm distr.
  
  mu.par[1] ~ dnorm(mean = l.lmean, sd = 0.10)
  mu.par[2] ~ dnorm(mean = k.lmean, sd = 0.10)
  
  sig.par[1] ~ dlnorm(log(l.lsd),10)  
  sig.par[2] ~ dlnorm(log(k.lsd),10)  
  
  #sig.par[1] ~ dlnorm(0,1)  #too diffuse for Lognormal sd!
  #sig.par[2] ~ dlnorm(0,1)
  
  # K & l values from fb  
  l.lsd <- sqrt(log(1 + (278^2) / (1378^2)))
  l.lmean <- log(1378) - 0.5 * sig.par[1]^2
  k.lsd <- sqrt(log(1 + (0.28^2) / (0.43^2)))
  k.lmean <- log(0.43) - 0.5 * sig.par[2]^2  
  
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



data <- sallaa %>%
  slice_sample(n = 10000) %>%
  filter(!is.na(juv.age),
         !is.na(year),
         year > 2000) %>%
  mutate(smo.age = if_else(age.type == "both", juv.age, NA)) 

data <- data %>%
  filter(!(age.type == "both" & juv.age == 0)) %>%
  mutate(yday.may1 = if_else(is.na(date),1, yday(date) - yday("2025-05-01")+1),
         doy.dec = if_else(yday.may1 > 0, yday.may1, yday.may1+365)/365,
         #yday.may1 =  yday(date) - yday("2025-05-01")+1,
         # yday.may1.b = if_else(yday.may1 > 0, yday.may1, yday.may1+365, missing = 1),
         # yday.may1.c = if_else(age.type == "juve.only",
         #                       min(yday(date) - yday("2025-01-31")+1, yday.may1.b),
         #                       yday.may1.b),
         # doy.dec = yday.may1.c/365,
         tot.age.2 = if_else(yday.may1 > 0, tot.age+doy.dec, tot.age+doy.dec-1),
         #       ) %>% # no 0 year olds caught before May 1
         # filter(tot.age == 0, 
         #        yday.may1 <0)
         smo.age = if_else(age.type == "both", juv.age, NA),
         sea.age = if_else(is.na(sea.age), 0, sea.age),
         #age.index = as.integer(tot.age + 1), # + 1 to index age 0 individuals
         age.index = trunc(tot.age.2)+1, # + 1 to index age 0 individuals
         hatch.year = year-tot.age,
         hatch.year.f = as.numeric(factor(hatch.year, levels = sort(unique(hatch.year)), labels = seq_along(sort(unique(hatch.year))))),
         year.f = as.integer(factor(year)),
         spat.unit = as.integer(factor(spat.unit))) %>%
  select(smo.age, tot.age, length_mm, age.index, year, spat.unit, 
         hatch.year, hatch.year.f, doy.dec, tot.age.2, date, age.type, doy.dec) %>%
  mutate(ind.id = row_number()) # ind id

data %>%
  ggplot() +
  geom_point(aes(doy.dec, length_mm, color = age.type)) +
  facet_wrap(~age.type)
  
n.ind <- nrow(data) # n individuals
npars <- 2 # n parameters in the LKJ dmnorm 
nsu <- length(unique(data$spat.unit))
age.index = data$age.index
n.ages = length(unique(data$age.index))+1
nK = max(data$smo.age %>% na.omit())   #max smolt age?
# get max year by, needed as hatch.year is integer and the last year is not (simpler solution?).
maxyear <- max(data$hatch.year.f) + data %>%
  filter(year == max(year)) %>%
  filter(tot.age == max(tot.age)) %>%
  distinct(tot.age) %>% pull(tot.age)
minyear = min(data$hatch.year.f)

inits <- function(){
  list(#l.mu = matrix(abs(rnorm(n.ages * n.ind,500,25)), nrow = n.ages, ncol = n.ind), # this is not correct, only inits of =<tot.age
    sig.l = abs(rnorm(1,250,25)),
    sig.par = c(rlnorm(1,0.20,0.10),
                rlnorm(1,0.60,0.10)),
    mu.par = c(rnorm(1,7,0.1),
               rnorm(1,-1,0.1)),
    lb.mu = abs(rnorm(nsu,20,5)),
    g = matrix(rnorm(nsu * length(minyear:maxyear),4,.2), nrow = nsu, ncol = length(minyear:maxyear))
    # smo.age = round(runif(n.ind,0,6))
    #Ustar, #U, #par
  )
}


bmconsts = list(npars = npars,
                nobs = nrow(data),
                n.ind = n.ind, # change when >1 sample per ind.id
                nK = nK,
                alpha = c(1, 2, 3, 2, 1, 1, 1)[1:nK],
                nsu = nsu,
                age = data$tot.age.2,
                age.index = data$age.index,
                doy.dec = data$doy.dec,
                ind.id = data$ind.id,
                hatch.year.f = data$hatch.year.f,
                minyear = minyear,
                maxyear = maxyear,
                spat.unit = data$spat.unit)


biph5.con.model <- nimbleModel(biph5.con.code,
                               constants=bmconsts,
                               inits = inits(),
                               data = data %>% select(-ind.id,-age.index,-spat.unit,
                                                      -hatch.year, -year, -tot.age,-tot.age.2,
                                                      -hatch.year.f, -doy.dec, -date, -age.type))
#g=array(0,dim=c(12,maxyear)),
biph5.con.model$simulate() 
biph5.con.model$calculate()    

dataNodes <- biph5.con.model$getNodeNames(dataOnly = TRUE)
parentNodes <- biph5.con.model$getParents(dataNodes, stochOnly = TRUE)       #all of these should be added to monitor below to recreate other model variables...
stnodes <- biph5.con.model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
allvars<-biph5.con.model$getVarNames(nodes = stnodes)
mvars<-allvars[!(grepl("lifted",allvars))]  



for(i in 1:length(mvars)){
  print(paste0(mvars[i]," ",biph5.con.model$calculate(mvars[i]) ))
}


# compile model
biph5.con.c <- compileNimble(biph5.con.model)

# configure..
biph5.con.confmcmc <- configureMCMC(biph5.con.c, monitors = c("g","par","lb.mu","sig.par","mu.par","l.mu","sig.l","smo.age"), project = biph5.con.model,useConjugacy = FALSE)

bMCMC <- buildMCMC(biph5.con.confmcmc)

CbmMCMC <- compileNimble(bMCMC, project = biph5.con.model)  

ptm<-proc.time()
CbmMCMC$run(50000, time = TRUE, nburnin=30000, thin=10)  
print(proc.time()-ptm)  

samples <- as.matrix(CbmMCMC$mvSamples)

obs <- data %>%
  #select(length_mm,tot.age,ind.id,spat.unit) %>%
  mutate(#tot.age = round(tot.age),
         observed = length_mm)

predi_obs <- samples[,1:5000] %>%
  gather_draws(l.mu[age.i,ind.id], sep = ",") %>%
  median_qi() %>%
  rename(predicted = .value) %>%
  mutate(tot.age = age.i) %>%
  filter(tot.age %in% obs$tot.age & ind.id %in% obs$ind.id)

predi_obs %>%
  left_join(obs, by = join_by(tot.age, ind.id)) %>%
  drop_na(observed) %>%
  ggplot(aes(x = observed, y = predicted, color = factor(spat.unit))) +
  geom_point(size = 2, alpha = 0.5) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0.2, alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  theme_minimal() +
  ylim(0,1300) +
  xlim(0,1300)

predi %>%
  left_join(obs, by = join_by(tot.age, ind.id)) %>%
  #drop_na(observed) %>%
  mutate(op.rat = if_else(is.na(observed), NA, observed/predicted),
         pps = if_else(op.rat < 0.5, 1, 0)) %>%
  filter(length_mm < 250) %>%
  ggplot() +
  geom_point(aes(doy.dec,length_mm, color = factor(pps))) +
  facet_wrap(~spat.unit) 

data %>%
  filter(!is.na(smo.age)) %>%
  select(smo.age, spat.unit, ind.id) %>%
  mutate(t= "obs") %>%
  bind_rows( samples %>% 
               gather_draws(smo.age[ind.id]) %>% 
               rename(smo.age = .value) %>% 
               mutate(t="pred") %>%
               left_join(data %>% select(spat.unit, ind.id), by = "ind.id") ) %>%
  ggplot() +
  geom_density(aes(smo.age, fill = t)) +
  facet_wrap(~spat.unit)

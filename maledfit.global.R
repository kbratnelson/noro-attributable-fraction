rm(list=ls()) 

##############################
#
# Load packages
#
##############################

library(rootSolve) #ODE solver
library(here) 
library(nloptr)
library(tidyverse)
library(deSolve)
library(grid)
library(binom)
library(kableExtra)

date=format(Sys.Date(), "%y%m%d")

#setwd("/home/knbratt/novburden") # FOR CLUSTER


##############################
#
# Input model parameters
#
##############################

# Demography 

#birth.rate = (25/1000)/365   # u5 mortality rate in the UK (changing this does not appreciably affect results)

# Rate of aging from 0-11mo category to 12-24 months age category
ageout = 1/365

# Natural history
sym = 1/3                    # symptomatic period, 1/days (Lopman et al 2014)
asym = 1/15                  # asymptomatic period, 1/days (Lopman et al 2014)
wan = 1/(365*5.1)            # immunity rate, 1/days (fitted in Lopman et al. 2014)
post = 0.05                  # relative infectiousness during asymptomatic infection (Lopman et al 2014)
#ps = 1/15                   # post-symptomatic shedding period, 1/days
other = 1/5                  # duration of non-norovirus diarrheal disease episode


# Initial conditions
xstart = c(S1 = 79000, Is1 = 1000, Ips1 = 0, Ia1 = 0, 
           R1 = 0, Isr1 = 0, Ipsr1 = 0, Iar1 = 0,
           Iother1 = 10000, Icoinf1 = 10000,
           S2 = 79000, Is2 = 1000, Ips2 = 0, Ia2 = 0, 
           R2 = 0, Isr2 = 0, Ipsr2 = 0, Iar2 = 0,
           Iother2 = 10000, Icoinf2 = 10000)

RM = function(t, x, RM_p) {
  
  # 2 age groups: <1y (suffix: 1), 1-2y (suffix: 2)
  beta.noro.1 <- RM_p[1]
  beta.noro.2 <- RM_p[2]
  foi.other.1 <- RM_p[3]
  foi.other.2 <- RM_p[4]
  s.p <- RM_p[5]
  s.r <- RM_p[6]
  
  
  S1  = x[1]   # susceptible
  Is1 = x[2]   # infected-symptomatic
  Ips1 = x[3]  # infected- postsymptomatic
  Ia1 = x[4]   # infected-asymptomatic
  
  R1  = x[5]   # recovered 
  Isr1 = x[6]  # infected-symptomatic reinfected
  Ipsr1 = x[7] # infected- postsymptomatic reinfected
  Iar1 = x[8]  # infected-asymptomatic reinfected
  
  Iother1 = x[9]  # diarrhea due to another pathogen
  Icoinf1 = x[10] # diarrhea due to another pathogen, with norovirus coinfection
  
  S2  = x[11]   # susceptible
  Is2 = x[12]   # infected-symptomatic
  Ips2 = x[13]  # infected- postsymptomatic
  Ia2 = x[14]   # infected-asymptomatic
  
  R2  = x[15]   # recovered 
  Isr2 = x[16]  # infected-symptomatic reinfected
  Ipsr2 = x[17] # infected- postsymptomatic reinfected
  Iar2 = x[18]  # infected-asymptomatic reinfected
  
  Iother2 = x[19] # diarrhea due to another pathogen
  Icoinf2 = x[20] # diarrhea due to another pathogen, with norovirus coinfection
  
  
  N = S1 + Is1 + Ips1 + Ia1 + R1 + Isr1 + Ipsr1 + Iar1 + Iother1 + Icoinf1 +
    S2 + Is2 + Ips2 + Ia2 + R2 + Isr2 + Ipsr2 + Iar2 + Iother2 + Icoinf2
  
  
  # Create population equillibrium by setting the number of births equal to the number of deaths
  births.day = ageout*S2 + ageout*Is2 + ageout*Ips2 + ageout*Ia2 + 
    ageout*R2 + ageout*Isr2 + ageout*Ipsr2 + ageout*Iar2 +
    ageout*Icoinf2 + ageout*Iother2
  
  #Calculate the number of individuals in each age group that can transmit the virus at each time step
  infec.noro = Is1 + Isr1 + post*(Ips1 + Ia1 + Ipsr1 + Iar1 + Icoinf1) +
    Is2 + Isr2 + post*(Ips2 + Ia2 + Ipsr2 + Iar2 + Icoinf2)
  
  #FoI for noro
  foi.noro.1 = beta.noro.1 * (infec.noro/N)
  foi.noro.2 = beta.noro.2 * (infec.noro/N)
  
  
  dS1 = births.day + wan*R1 - foi.noro.1*S1 - foi.other.1*S1 + other*Iother1 - ageout*S1 
  dIs1 = s.p*foi.noro.1*S1 - sym*Is1 - ageout*Is1 - foi.other.1*Is1
  dIps1 = sym*Is1 - asym*Ips1 - ageout*Ips1 - foi.other.1*Ips1
  dIa1 = (1-s.p)*foi.noro.1*S1 - asym*Ia1 - ageout*Ia1 - foi.other.1*Ia1
  
  dR1 = asym*Ia1 + asym*Ips1 + asym*Iar1 + asym*Ipsr1 + other*Icoinf1 - foi.other.1*R1 - foi.noro.1*R1 - wan*R1 - ageout*R1 
  dIsr1 = s.r*foi.noro.1*R1 - sym*Isr1 - ageout*Isr1 - foi.other.1*Isr1
  dIpsr1 = sym*Isr1 - asym*Ipsr1 - ageout*Ipsr1 - foi.other.1*Ipsr1
  dIar1 = (1-s.r)*foi.noro.1*R1 - asym*Iar1 - ageout*Iar1 - foi.other.1*Iar1
  
  dIother1 = foi.other.1*S1 + foi.other.1*R1 - other*Iother1 - foi.noro.1*Iother1 - ageout*Iother1 
  dIcoinf1 = foi.noro.1*Iother1 - other*Icoinf1 - ageout*Icoinf1 + foi.other.1*(Is1 + Ia1 + Ips1 + Isr1 + Iar1 + Ipsr1)
  
  #
  
  dS2 = wan*R2 - foi.noro.2*S2 - foi.other.2*S2 + other*Iother2 + ageout*S1 - ageout*S2
  dIs2 = s.p*foi.noro.2*S2 - sym*Is2 + ageout*Is1 - ageout*Is2 - foi.other.2*Is2
  dIps2 = sym*Is2 - asym*Ips2 + ageout*Ips1 - ageout*Ips2 - foi.other.2*Ips2
  dIa2 = (1-s.p)*foi.noro.2*S2 - asym*Ia2 + ageout*Ia1 - ageout*Ia2 - foi.other.2*Ia2
  
  dR2 = asym*Ia2 + asym*Ips2 + asym*Iar2 + asym*Ipsr2 + other*Icoinf2 - foi.other.2*R2 - foi.noro.2*R2 - wan*R2 + ageout*R1 - ageout*R2
  dIsr2 = s.r*foi.noro.2*R2 - sym*Isr2 + ageout*Isr1 - ageout*Isr2 - foi.other.2*Isr2
  dIpsr2 = sym*Isr2 - asym*Ipsr2 + ageout*Ipsr1 - ageout*Ipsr2 - foi.other.2*Ipsr2
  dIar2 = (1-s.r)*foi.noro.2*R2 - asym*Iar2 + ageout*Iar1 - ageout*Iar2 - foi.other.2*Iar2
  
  dIother2 = foi.other.2*S2 + foi.other.2*R2 - other*Iother2 - foi.noro.2*Iother2 + ageout*Iother1 - ageout*Iother2
  dIcoinf2 = foi.noro.2*Iother2 - other*Icoinf2 + ageout*Icoinf1 - ageout*Icoinf2 + foi.other.2*(Is2 + Ia2 + Ips2 + Isr2 + Iar2 + Ipsr2)
  
  ##
  
  res = c(dS1, dIs1, dIps1, dIa1, dR1, dIsr1, dIpsr1, dIar1, dIother1, dIcoinf1,
          dS2, dIs2, dIps2, dIa2, dR2, dIsr2, dIpsr2, dIar2, dIother2, dIcoinf2)
  list(res) 
  
}



RM_solve = function(p) {
  
  beta.noro.1 <- p[1]
  beta.noro.2 <- p[2]
  foi.other.1 <- p[3]
  foi.other.2 <- p[4]
  s.p <- p[5]
  s.r <- p[6]
  
  out = runsteady(y = xstart, fun = RM, parms = c(beta.noro.1, beta.noro.2, foi.other.1, foi.other.2, s.p, s.r), times = c(0, 1e5))
  
  out = out$y
  
  N.1 = sum(out["S1"], out["Is1"], out["Ips1"], out["Ia1"], out["R1"], out["Isr1"], out["Ipsr1"], out["Iar1"], out["Iother1"], out["Icoinf1"])
  N.2 = sum(out["S2"], out["Is2"], out["Ips2"], out["Ia2"], out["R2"], out["Isr2"], out["Ipsr2"], out["Iar2"], out["Iother2"], out["Icoinf2"])
  N = N.1 + N.2
  
  
  # Calculate prevalence of norovirus infection among healthy (not ill) kids - CROSS-SEC (prevalent controls)
  mod.asymp.1 = sum(out["Ips1"], out["Ia1"], out["Ipsr1"], out["Iar1"]) /
    sum(out["Ips1"], out["Ia1"], out["Ipsr1"], out["Iar1"], out["S1"], out["R1"]) 
  mod.asymp.2 = sum(out["Ips2"], out["Ia2"], out["Ipsr2"], out["Iar2"]) /
    sum(out["Ips2"], out["Ia2"], out["Ipsr2"], out["Iar2"], out["S2"], out["R2"]) 
  
  # Calculate prevalence of norovirus infection among cases (sick) kids - CROSS-SEC (prevalent cases)
  mod.case.1 = sum(out["Is1"], out["Icoinf1"]) / sum(out["Is1"], out["Icoinf1"], out["Iother1"])
  mod.case.2 = sum(out["Is2"], out["Icoinf2"]) / sum(out["Is2"], out["Icoinf2"], out["Iother2"])
  
  # Calculate prevalence of norovirus infection among cases (sick) kids - COHORT (incident cases)
  #Calculate the number of individuals in each age group that can transmit the virus at each time step
  infec.noro = out["Is1"] + out["Isr1"] + post*(out["Ips1"] + out["Ia1"] + out["Ipsr1"] + out["Iar1"] + out["Icoinf1"]) +
    out["Is2"] + out["Isr2"] + post*(out["Ips2"] + out["Ia2"] + out["Ipsr2"] + out["Iar2"] + out["Icoinf2"])
  foi.noro.1 = beta.noro.1 * (infec.noro/N)
  foi.noro.2 = beta.noro.2 * (infec.noro/N)
  
  # Calculate overall AGE incidence (per person-year) for each age group
  age.inc.mod.1 = foi.noro.1*365*(out["S1"]*s.p + out["R1"]*s.r)/N.1 + 
    foi.other.1*365*out["S1"]/N.1 + 
    foi.noro.1*365*out["Iother1"]/N.1 +
    foi.other.1*365*(out["Ia1"] + out["Ips1"] +  out["Ipsr1"] + out["Iar1"])/N.1
  age.inc.mod.2 = foi.noro.2*365*(out["S2"]*s.p + out["R2"]*s.r)/N.2 + 
    foi.other.2*365*out["S2"]/N.2 + 
    foi.noro.2*365*out["Iother2"]/N.2 +
    foi.other.2*365*(out["Ia2"] + out["Ips2"] +  out["Ipsr2"] + out["Iar2"])/N.2
  
  # No AGE incidence
  #ssd = (abs(mod.asymp.1 - obs.asymp.1) + abs(mod.asymp.2 - obs.asymp.2) + abs(mod.case.1 - obs.case.1) + abs(mod.case.2 - obs.case.2))^2 # this is summary stat for one simulation - the sum squared difference between model and data
  
  # With AGE incidence
  ssd = (abs(mod.asymp.1 - obs.asymp.1) + abs(mod.asymp.2 - obs.asymp.2) + abs(mod.case.1 - obs.case.1) + abs(mod.case.2 - obs.case.2) + (abs(age.inc.mod.1-age.inc.1y)/10) +  (abs(age.inc.mod.2-age.inc.2y)/10))^2 # this is summary stat for one simulation - the sum squared difference between model and data
  
  if (is.na(ssd)) {ssd = 1e10} #if something goes wrong, we don't want to return NA, instead return a very high value for the NLL
  
  return(ssd) 
}

beta.noro.1.min <- 0
beta.noro.1.max <- 5
beta.noro.2.min <- 0
beta.noro.2.max <- 5
foi.noro.1.min <- 0
foi.noro.1.max <- 1
foi.noro.2.min <- 0
foi.noro.2.max <- 1

s.p.min <- 0.05
s.p.max <- 0.95
s.r.min <- 0.05
s.r.max <- 0.95

lb=c(beta.noro.1.min, beta.noro.2.min, foi.noro.1.min, foi.noro.2.min, s.p.min, s.r.min)
ub=c(beta.noro.1.max, beta.noro.2.max, foi.noro.1.max, foi.noro.2.max, s.p.max, s.r.max)

p0.bang = c(0.5, 0.5, 0.01, 0.01, 0.8, 0.2) # all starting values the same for global fit 
p0.ind = c(0.5, 0.5, 0.01, 0.01, 0.8, 0.2) 
p0.nep = c(0.5, 0.5, 0.01, 0.01, 0.8, 0.2)
p0.pak = c(0.5, 0.5, 0.01, 0.01, 0.8, 0.2)
p0.braz = c(0.5, 0.5, 0.01, 0.01, 0.8, 0.2)
p0.peru = c(0.5, 0.5, 0.01, 0.01, 0.8, 0.2)

#some general settings
ploton=0
tstart=proc.time() #capture current time to measure speed
plevel=0
xerr.abs=0.01
max.local.stps=250000
stps=250000
mtime=72*60*60 # runtime per scenario (in seconds, so 1*60*60 is 1 hour)

solvertype=11

#nloptr global solvers
if (solvertype==11) { alg.name="NLOPT_GN_DIRECT"   }  
if (solvertype==12) { alg.name="NLOPT_GN_CRS2_LM"  }
if (solvertype==13) {  alg.name="NLOPT_GN_ISRES"   }
if (solvertype==15) {  alg.name="NLOPT_GN_ESCH"    }  

#nloptr local gradient-free solvers
if (solvertype==21)      {        alg.name="NLOPT_LN_SBPLX"      }
if (solvertype==22)      {        alg.name="NLOPT_LN_COBYLA"     }
if (solvertype==23)      {        alg.name="NLOPT_LN_BOBYQA"     }     
if (solvertype==24)      {        alg.name="NLOPT_LN_PRAXIS"     }
if (solvertype==25)      {        alg.name="NLOPT_LN_NELDERMEAD" }
if (solvertype==26)      {        alg.name="NLOPT_GN_DIRECT_L_NOSCAL" }

#######################################################################################
# READ IN THE DATA FROM MAL-ED study to fit models to
#######################################################################################

# This data is from: Platt-Mills et al. Lancet 2011, Supplementary Table S1.
# local: maled_forfit <- read.csv('Other/maled.cacoprev.csv')

maled_forfit <- read.csv(here('Other/maled.cacoprev.csv')) 
#maled_forfit <- read.csv('maled.cacoprev.csv') # FOR CLUSTER

rownames(maled_forfit) <- maled_forfit[,1]
maled_forfit <- maled_forfit[,-1]

# Fit for all norovirus (allno) - combined g1 + g2 prevalences
# note that MAL-ED did not report prevalence in certain countries/age groups if < 0.01. In these cases, assume prevalence was 0.01.

# Also in this file is AGE prevalence, based on MAL-ED. This is needed to calculate prevaence among *cases* under each model. Briefly:
# Use Table 2 to calculate incidence (number of diarrheal eps/children*2 years). Number of children per age group calculated by taking number by site (in Table 2) and dividing by two.
# Note that this assumes all kids followed for 2 years (appears all in this table were), and equal size of 0-11mo and 12-23mo groups (true if all kids followed for full 2 year study duration). In addition, to calculate prevalence from incidence assume a duration of 5 days.

# CALCULATE 95% CIs for prevalence of shedding in each group
ci95.case.1y <- binom.confint(maled_forfit[,"case.num.1y"]*maled_forfit[,"case.allno.1y"], maled_forfit[,"case.num.1y"], conf.level = 0.95, methods = "exact")
ci95.case.1y.low <- ci95.case.1y[,"lower"]
ci95.case.1y.upp <- ci95.case.1y[,"upper"]

ci95.case.2y <- binom.confint(maled_forfit[,"case.num.2y"]*maled_forfit[,"case.allno.2y"], maled_forfit[,"case.num.2y"], conf.level = 0.95, methods = "exact")
ci95.case.2y.low <- ci95.case.2y[,"lower"]
ci95.case.2y.upp <- ci95.case.2y[,"upper"]

ci95.control.1y <- binom.confint(maled_forfit[,"control.num.1y"]*maled_forfit[,"control.allno.1y"], maled_forfit[,"control.num.1y"], conf.level = 0.95, methods = "exact")
ci95.control.1y.low <- ci95.control.1y[,"lower"]
ci95.control.1y.upp <- ci95.control.1y[,"upper"]

ci95.control.2y <- binom.confint(maled_forfit[,"control.num.2y"]*maled_forfit[,"control.allno.2y"], maled_forfit[,"control.num.2y"], conf.level = 0.95, methods = "exact")
ci95.control.2y.low <- ci95.control.2y[,"lower"]
ci95.control.2y.upp <- ci95.control.2y[,"upper"]

maled_forfit <- cbind(maled_forfit, 
                      ci95.case.1y.low, ci95.case.1y.upp,
                      ci95.control.1y.low, ci95.control.1y.upp,
                      ci95.case.2y.low, ci95.case.2y.upp,
                      ci95.control.2y.low, ci95.control.2y.upp)



#######################################################################################
# BANGLADESH # 
#######################################################################################

obs.asymp.1.pt <- maled_forfit["Bangladesh","control.allno.1y"] 
obs.asymp.1.low <- maled_forfit["Bangladesh","ci95.control.1y.low"] 
obs.asymp.1.upp <- maled_forfit["Bangladesh","ci95.control.1y.upp"] 

obs.asymp.2.pt <- maled_forfit["Bangladesh","control.allno.2y"] 
obs.asymp.2.low <- maled_forfit["Bangladesh","ci95.control.2y.low"] 
obs.asymp.2.upp <- maled_forfit["Bangladesh","ci95.control.2y.upp"] 

obs.case.1.pt <- maled_forfit["Bangladesh","case.allno.1y"] 

obs.case.2.pt <- maled_forfit["Bangladesh","case.allno.2y"] 

age.inc.1y <- maled_forfit["Bangladesh","age.inc.1y.jpm"] 
age.inc.2y <- maled_forfit["Bangladesh","age.inc.2y.jpm"]

##############################
obs.asymp.1 <- obs.asymp.1.pt  # 0-11mo target
obs.asymp.2 <- obs.asymp.2.pt  # 12-23mo target
obs.case.1 <- obs.case.1.pt    # 0-11mo target
obs.case.2 <- obs.case.2.pt    # 12-23mo target
##############################

fres <- nloptr(x0=p0.bang,eval_f=RM_solve,lb=lb,ub=ub,opts=list("algorithm"=alg.name,xtol_abs=1e-3,maxeval=stps,maxtime=mtime,print_level=1))

fpar=(fres$solution); asymout.final=fres$objective; totiterate=fres$iterations; solverstatus=fres$status
bestfit=c(fpar, asymout.final, totiterate, solverstatus)
bestfit <- as.data.frame(bestfit)
rownames(bestfit)=c("beta.noro.1", "beta.noro.2", "foi.other.1", "foi.other.2", "s.p", "s.r", "objective function value", 
                    "total iterations", "solver status")

beta.noro.1 <- bestfit[1,1]
beta.noro.2 <- bestfit[2,1]
foi.other.1 <- bestfit[3,1]
foi.other.2 <- bestfit[4,1]
s.p <- bestfit[5,1]
s.r <- bestfit[6,1]
obj.fx.val <- bestfit[7,1]

out = runsteady(y = xstart, fun = RM, parms = c(beta.noro.1, beta.noro.2, foi.other.1, foi.other.2, s.p, s.r), times = c(0, 1e5))

out = out$y


# Calculate pop sizes
N.1 = sum(out["S1"], out["Is1"], out["Ips1"], out["Ia1"], out["R1"], out["Isr1"], out["Ipsr1"], out["Iar1"], out["Iother1"], out["Icoinf1"])
N.2 = sum(out["S2"], out["Is2"], out["Ips2"], out["Ia2"], out["R2"], out["Isr2"], out["Ipsr2"], out["Iar2"], out["Iother2"], out["Icoinf2"])
N = N.1 + N.2

# Calculate prevalence of norovirus infection among healthy (not ill) kids - CROSS-SEC (prevalent controls)
mod.asymp.1 = sum(out["Ips1"], out["Ia1"], out["Ipsr1"], out["Iar1"]) /
  sum(out["Ips1"], out["Ia1"], out["Ipsr1"], out["Iar1"], out["S1"], out["R1"]) 
mod.asymp.2 = sum(out["Ips2"], out["Ia2"], out["Ipsr2"], out["Iar2"]) /
  sum(out["Ips2"], out["Ia2"], out["Ipsr2"], out["Iar2"], out["S2"], out["R2"]) 

# Calculate prevalence of norovirus infection among cases (sick) kids - CROSS-SEC (prevalent cases)
mod.case.1 = sum(out["Is1"], out["Icoinf1"]) / sum(out["Is1"], out["Icoinf1"], out["Iother1"])
mod.case.2 = sum(out["Is2"], out["Icoinf2"]) / sum(out["Is2"], out["Icoinf2"], out["Iother2"])

# Calculate prevalence of norovirus infection among cases (sick) kids - COHORT (incident cases)
infec.noro = out["Is1"] + out["Isr1"] + post*(out["Ips1"] + out["Ia1"] + out["Ipsr1"] + out["Iar1"] + out["Icoinf1"]) +
  out["Is2"] + out["Isr2"] + post*(out["Ips2"] + out["Ia2"] + out["Ipsr2"] + out["Iar2"] + out["Icoinf2"])
foi.noro.1 = beta.noro.1 * (infec.noro/N)
foi.noro.2 = beta.noro.2 * (infec.noro/N)

# Calculate overall AGE incidence (per person-year) for each age group
age.inc.mod.1 = foi.noro.1*365*(out["S1"]*s.p + out["R1"]*s.r)/N.1 + 
  foi.other.1*365*out["S1"]/N.1 + 
  foi.noro.1*365*out["Iother1"]/N.1 +
  foi.other.1*365*(out["Ia1"] + out["Ips1"] +  out["Ipsr1"] + out["Iar1"])/N.1
age.inc.mod.2 = foi.noro.2*365*(out["S2"]*s.p + out["R2"]*s.r)/N.2 + 
  foi.other.2*365*out["S2"]/N.2 + 
  foi.noro.2*365*out["Iother2"]/N.2 +
  foi.other.2*365*(out["Ia2"] + out["Ips2"] +  out["Ipsr2"] + out["Iar2"])/N.2

noro.inc.mod.1 = foi.noro.1*365*((out["S1"]*s.p + out["R1"]*s.r)/N.1)
noro.inc.mod.2 = foi.noro.2*365*((out["S2"]*s.p + out["R2"]*s.r)/N.2)

bang.pt <- c(noro.inc.mod.1, noro.inc.mod.2, 
             beta.noro.1, beta.noro.2, foi.other.1, foi.other.2, s.p, s.r, 
             age.inc.mod.1, age.inc.mod.2, age.inc.1y, age.inc.2y,
             obs.asymp.1, mod.asymp.1, obs.asymp.2, mod.asymp.2, obs.case.1, mod.case.1, obs.case.2, mod.case.2,
             obj.fx.val)

# Combine point esimate and CIs into one dataframe
#bang_all <- cbind(bang.pt, bang.low, bang.upp)

rm(noro.inc.mod.1, noro.inc.mod.2, 
   beta.noro.1, beta.noro.2, foi.other.1, foi.other.2, s.p, s.r,
   age.inc.mod.1, age.inc.mod.2, age.inc.1y, age.inc.2y,
   obs.asymp.1, mod.asymp.1, obs.asymp.2, mod.asymp.2, obs.case.1, mod.case.1, obs.case.2, mod.case.2,
   obj.fx.val)


#######################################################################################
# INDIA # 
#######################################################################################

obs.asymp.1.pt <- maled_forfit["India","control.allno.1y"] 
obs.asymp.1.low <- maled_forfit["India","ci95.control.1y.low"] 
obs.asymp.1.upp <- maled_forfit["India","ci95.control.1y.upp"] 

obs.asymp.2.pt <- maled_forfit["India","control.allno.2y"] 
obs.asymp.2.low <- maled_forfit["India","ci95.control.2y.low"] 
obs.asymp.2.upp <- maled_forfit["India","ci95.control.2y.upp"] 

obs.case.1.pt <- maled_forfit["India","case.allno.1y"] 

obs.case.2.pt <- maled_forfit["India","case.allno.2y"] 

age.inc.1y <- maled_forfit["India","age.inc.1y.jpm"] 
age.inc.2y <- maled_forfit["India","age.inc.2y.jpm"]

##############################
obs.asymp.1 <- obs.asymp.1.pt  # 0-11mo target
obs.asymp.2 <- obs.asymp.2.pt  # 12-23mo target
obs.case.1 <- obs.case.1.pt    # 0-11mo target
obs.case.2 <- obs.case.2.pt    # 12-23mo target
##############################

fres <- nloptr(x0=p0.bang,eval_f=RM_solve,lb=lb,ub=ub,opts=list("algorithm"=alg.name,xtol_abs=1e-3,maxeval=stps,maxtime=mtime,print_level=1))

fpar=(fres$solution); asymout.final=fres$objective; totiterate=fres$iterations; solverstatus=fres$status
bestfit=c(fpar, asymout.final, totiterate, solverstatus)
bestfit <- as.data.frame(bestfit)
rownames(bestfit)=c("beta.noro.1", "beta.noro.2", "foi.other.1", "foi.other.2", "s.p", "s.r", "objective function value", 
                    "total iterations", "solver status")

beta.noro.1 <- bestfit[1,1]
beta.noro.2 <- bestfit[2,1]
foi.other.1 <- bestfit[3,1]
foi.other.2 <- bestfit[4,1]
s.p <- bestfit[5,1]
s.r <- bestfit[6,1]
obj.fx.val <- bestfit[7,1]

out = runsteady(y = xstart, fun = RM, parms = c(beta.noro.1, beta.noro.2, foi.other.1, foi.other.2, s.p, s.r), times = c(0, 1e5))

out = out$y


# Calculate pop sizes
N.1 = sum(out["S1"], out["Is1"], out["Ips1"], out["Ia1"], out["R1"], out["Isr1"], out["Ipsr1"], out["Iar1"], out["Iother1"], out["Icoinf1"])
N.2 = sum(out["S2"], out["Is2"], out["Ips2"], out["Ia2"], out["R2"], out["Isr2"], out["Ipsr2"], out["Iar2"], out["Iother2"], out["Icoinf2"])
N = N.1 + N.2

# Calculate prevalence of norovirus infection among healthy (not ill) kids - CROSS-SEC (prevalent controls)
mod.asymp.1 = sum(out["Ips1"], out["Ia1"], out["Ipsr1"], out["Iar1"]) /
  sum(out["Ips1"], out["Ia1"], out["Ipsr1"], out["Iar1"], out["S1"], out["R1"]) 
mod.asymp.2 = sum(out["Ips2"], out["Ia2"], out["Ipsr2"], out["Iar2"]) /
  sum(out["Ips2"], out["Ia2"], out["Ipsr2"], out["Iar2"], out["S2"], out["R2"]) 

# Calculate prevalence of norovirus infection among cases (sick) kids - CROSS-SEC (prevalent cases)
mod.case.1 = sum(out["Is1"], out["Icoinf1"]) / sum(out["Is1"], out["Icoinf1"], out["Iother1"])
mod.case.2 = sum(out["Is2"], out["Icoinf2"]) / sum(out["Is2"], out["Icoinf2"], out["Iother2"])

# Calculate prevalence of norovirus infection among cases (sick) kids - COHORT (incident cases)
infec.noro = out["Is1"] + out["Isr1"] + post*(out["Ips1"] + out["Ia1"] + out["Ipsr1"] + out["Iar1"] + out["Icoinf1"]) +
  out["Is2"] + out["Isr2"] + post*(out["Ips2"] + out["Ia2"] + out["Ipsr2"] + out["Iar2"] + out["Icoinf2"])
foi.noro.1 = beta.noro.1 * (infec.noro/N)
foi.noro.2 = beta.noro.2 * (infec.noro/N)

# Calculate overall AGE incidence (per person-year) for each age group
age.inc.mod.1 = foi.noro.1*365*(out["S1"]*s.p + out["R1"]*s.r)/N.1 + 
  foi.other.1*365*out["S1"]/N.1 + 
  foi.noro.1*365*out["Iother1"]/N.1 +
  foi.other.1*365*(out["Ia1"] + out["Ips1"] +  out["Ipsr1"] + out["Iar1"])/N.1
age.inc.mod.2 = foi.noro.2*365*(out["S2"]*s.p + out["R2"]*s.r)/N.2 + 
  foi.other.2*365*out["S2"]/N.2 + 
  foi.noro.2*365*out["Iother2"]/N.2 +
  foi.other.2*365*(out["Ia2"] + out["Ips2"] +  out["Ipsr2"] + out["Iar2"])/N.2

noro.inc.mod.1 = foi.noro.1*365*((out["S1"]*s.p + out["R1"]*s.r)/N.1)
noro.inc.mod.2 = foi.noro.2*365*((out["S2"]*s.p + out["R2"]*s.r)/N.2)

ind.pt <- c(noro.inc.mod.1, noro.inc.mod.2, 
            beta.noro.1, beta.noro.2, foi.other.1, foi.other.2, s.p, s.r, 
            age.inc.mod.1, age.inc.mod.2, age.inc.1y, age.inc.2y,
            obs.asymp.1, mod.asymp.1, obs.asymp.2, mod.asymp.2, obs.case.1, mod.case.1, obs.case.2, mod.case.2,
            obj.fx.val)

# Combine point esimate and CIs into one dataframe
#bang_all <- cbind(bang.pt, bang.low, bang.upp)

rm(noro.inc.mod.1, noro.inc.mod.2, 
   beta.noro.1, beta.noro.2, foi.other.1, foi.other.2, s.p, s.r,
   age.inc.mod.1, age.inc.mod.2, age.inc.1y, age.inc.2y,
   obs.asymp.1, mod.asymp.1, obs.asymp.2, mod.asymp.2, obs.case.1, mod.case.1, obs.case.2, mod.case.2, 
   obj.fx.val)




#######################################################################################
# NEPAL #
#######################################################################################

obs.asymp.1.pt <- maled_forfit["Nepal","control.allno.1y"] 
obs.asymp.1.low <- maled_forfit["Nepal","ci95.control.1y.low"] 
obs.asymp.1.upp <- maled_forfit["Nepal","ci95.control.1y.upp"] 

obs.asymp.2.pt <- maled_forfit["Nepal","control.allno.2y"] 
obs.asymp.2.low <- maled_forfit["Nepal","ci95.control.2y.low"] 
obs.asymp.2.upp <- maled_forfit["Nepal","ci95.control.2y.upp"] 

obs.case.1.pt <- maled_forfit["Nepal","case.allno.1y"] 

obs.case.2.pt <- maled_forfit["Nepal","case.allno.2y"] 

age.inc.1y <- maled_forfit["Nepal","age.inc.1y.jpm"] 
age.inc.2y <- maled_forfit["Nepal","age.inc.2y.jpm"]

##############################
obs.asymp.1 <- obs.asymp.1.pt  # 0-11mo target
obs.asymp.2 <- obs.asymp.2.pt  # 12-23mo target
obs.case.1 <- obs.case.1.pt    # 0-11mo target
obs.case.2 <- obs.case.2.pt    # 12-23mo target
##############################

fres <- nloptr(x0=p0.bang,eval_f=RM_solve,lb=lb,ub=ub,opts=list("algorithm"=alg.name,xtol_abs=1e-3,maxeval=stps,maxtime=mtime,print_level=1))

fpar=(fres$solution); asymout.final=fres$objective; totiterate=fres$iterations; solverstatus=fres$status
bestfit=c(fpar, asymout.final, totiterate, solverstatus)
bestfit <- as.data.frame(bestfit)
rownames(bestfit)=c("beta.noro.1", "beta.noro.2", "foi.other.1", "foi.other.2", "s.p", "s.r", "objective function value", 
                    "total iterations", "solver status")

beta.noro.1 <- bestfit[1,1]
beta.noro.2 <- bestfit[2,1]
foi.other.1 <- bestfit[3,1]
foi.other.2 <- bestfit[4,1]
s.p <- bestfit[5,1]
s.r <- bestfit[6,1]
obj.fx.val <- bestfit[7,1]

out = runsteady(y = xstart, fun = RM, parms = c(beta.noro.1, beta.noro.2, foi.other.1, foi.other.2, s.p, s.r), times = c(0, 1e5))

out = out$y


# Calculate pop sizes
N.1 = sum(out["S1"], out["Is1"], out["Ips1"], out["Ia1"], out["R1"], out["Isr1"], out["Ipsr1"], out["Iar1"], out["Iother1"], out["Icoinf1"])
N.2 = sum(out["S2"], out["Is2"], out["Ips2"], out["Ia2"], out["R2"], out["Isr2"], out["Ipsr2"], out["Iar2"], out["Iother2"], out["Icoinf2"])
N = N.1 + N.2

# Calculate prevalence of norovirus infection among healthy (not ill) kids - CROSS-SEC (prevalent controls)
mod.asymp.1 = sum(out["Ips1"], out["Ia1"], out["Ipsr1"], out["Iar1"]) /
  sum(out["Ips1"], out["Ia1"], out["Ipsr1"], out["Iar1"], out["S1"], out["R1"]) 
mod.asymp.2 = sum(out["Ips2"], out["Ia2"], out["Ipsr2"], out["Iar2"]) /
  sum(out["Ips2"], out["Ia2"], out["Ipsr2"], out["Iar2"], out["S2"], out["R2"]) 

# Calculate prevalence of norovirus infection among cases (sick) kids - CROSS-SEC (prevalent cases)
mod.case.1 = sum(out["Is1"], out["Icoinf1"]) / sum(out["Is1"], out["Icoinf1"], out["Iother1"])
mod.case.2 = sum(out["Is2"], out["Icoinf2"]) / sum(out["Is2"], out["Icoinf2"], out["Iother2"])

# Calculate prevalence of norovirus infection among cases (sick) kids - COHORT (incident cases)
infec.noro = out["Is1"] + out["Isr1"] + post*(out["Ips1"] + out["Ia1"] + out["Ipsr1"] + out["Iar1"] + out["Icoinf1"]) +
  out["Is2"] + out["Isr2"] + post*(out["Ips2"] + out["Ia2"] + out["Ipsr2"] + out["Iar2"] + out["Icoinf2"])
foi.noro.1 = beta.noro.1 * (infec.noro/N)
foi.noro.2 = beta.noro.2 * (infec.noro/N)

# Calculate overall AGE incidence (per person-year) for each age group
age.inc.mod.1 = foi.noro.1*365*(out["S1"]*s.p + out["R1"]*s.r)/N.1 + 
  foi.other.1*365*out["S1"]/N.1 + 
  foi.noro.1*365*out["Iother1"]/N.1 +
  foi.other.1*365*(out["Ia1"] + out["Ips1"] +  out["Ipsr1"] + out["Iar1"])/N.1
age.inc.mod.2 = foi.noro.2*365*(out["S2"]*s.p + out["R2"]*s.r)/N.2 + 
  foi.other.2*365*out["S2"]/N.2 + 
  foi.noro.2*365*out["Iother2"]/N.2 +
  foi.other.2*365*(out["Ia2"] + out["Ips2"] +  out["Ipsr2"] + out["Iar2"])/N.2

noro.inc.mod.1 = foi.noro.1*365*((out["S1"]*s.p + out["R1"]*s.r)/N.1)
noro.inc.mod.2 = foi.noro.2*365*((out["S2"]*s.p + out["R2"]*s.r)/N.2)

nep.pt <- c(noro.inc.mod.1, noro.inc.mod.2, 
            beta.noro.1, beta.noro.2, foi.other.1, foi.other.2, s.p, s.r, 
            age.inc.mod.1, age.inc.mod.2, age.inc.1y, age.inc.2y,
            obs.asymp.1, mod.asymp.1, obs.asymp.2, mod.asymp.2, obs.case.1, mod.case.1, obs.case.2, mod.case.2,
            obj.fx.val)

# Combine point esimate and CIs into one dataframe
#bang_all <- cbind(bang.pt, bang.low, bang.upp)

rm(noro.inc.mod.1, noro.inc.mod.2, 
   beta.noro.1, beta.noro.2, foi.other.1, foi.other.2, s.p, s.r,
   age.inc.mod.1, age.inc.mod.2, age.inc.1y, age.inc.2y,
   obs.asymp.1, mod.asymp.1, obs.asymp.2, mod.asymp.2, obs.case.1, mod.case.1, obs.case.2, mod.case.2, 
   obj.fx.val)



#######################################################################################
# PAKISTAN #
#######################################################################################

obs.asymp.1.pt <- maled_forfit["Pakistan","control.allno.1y"] 
obs.asymp.1.low <- maled_forfit["Pakistan","ci95.control.1y.low"] 
obs.asymp.1.upp <- maled_forfit["Pakistan","ci95.control.1y.upp"] 

obs.asymp.2.pt <- maled_forfit["Pakistan","control.allno.2y"] 
obs.asymp.2.low <- maled_forfit["Pakistan","ci95.control.2y.low"] 
obs.asymp.2.upp <- maled_forfit["Pakistan","ci95.control.2y.upp"] 

obs.case.1.pt <- maled_forfit["Pakistan","case.allno.1y"] 

obs.case.2.pt <- maled_forfit["Pakistan","case.allno.2y"] 

age.inc.1y <- maled_forfit["Pakistan","age.inc.1y.jpm"] 
age.inc.2y <- maled_forfit["Pakistan","age.inc.2y.jpm"]

##############################
obs.asymp.1 <- obs.asymp.1.pt  # 0-11mo target
obs.asymp.2 <- obs.asymp.2.pt  # 12-23mo target
obs.case.1 <- obs.case.1.pt    # 0-11mo target
obs.case.2 <- obs.case.2.pt    # 12-23mo target
##############################

fres <- nloptr(x0=p0.bang,eval_f=RM_solve,lb=lb,ub=ub,opts=list("algorithm"=alg.name,xtol_abs=1e-3,maxeval=stps,maxtime=mtime,print_level=1))

fpar=(fres$solution); asymout.final=fres$objective; totiterate=fres$iterations; solverstatus=fres$status
bestfit=c(fpar, asymout.final, totiterate, solverstatus)
bestfit <- as.data.frame(bestfit)
rownames(bestfit)=c("beta.noro.1", "beta.noro.2", "foi.other.1", "foi.other.2", "s.p", "s.r", "objective function value", 
                    "total iterations", "solver status")

beta.noro.1 <- bestfit[1,1]
beta.noro.2 <- bestfit[2,1]
foi.other.1 <- bestfit[3,1]
foi.other.2 <- bestfit[4,1]
s.p <- bestfit[5,1]
s.r <- bestfit[6,1]
obj.fx.val <- bestfit[7,1]

out = runsteady(y = xstart, fun = RM, parms = c(beta.noro.1, beta.noro.2, foi.other.1, foi.other.2, s.p, s.r), times = c(0, 1e5))

out = out$y


# Calculate pop sizes
N.1 = sum(out["S1"], out["Is1"], out["Ips1"], out["Ia1"], out["R1"], out["Isr1"], out["Ipsr1"], out["Iar1"], out["Iother1"], out["Icoinf1"])
N.2 = sum(out["S2"], out["Is2"], out["Ips2"], out["Ia2"], out["R2"], out["Isr2"], out["Ipsr2"], out["Iar2"], out["Iother2"], out["Icoinf2"])
N = N.1 + N.2

# Calculate prevalence of norovirus infection among healthy (not ill) kids - CROSS-SEC (prevalent controls)
mod.asymp.1 = sum(out["Ips1"], out["Ia1"], out["Ipsr1"], out["Iar1"]) /
  sum(out["Ips1"], out["Ia1"], out["Ipsr1"], out["Iar1"], out["S1"], out["R1"]) 
mod.asymp.2 = sum(out["Ips2"], out["Ia2"], out["Ipsr2"], out["Iar2"]) /
  sum(out["Ips2"], out["Ia2"], out["Ipsr2"], out["Iar2"], out["S2"], out["R2"]) 

# Calculate prevalence of norovirus infection among cases (sick) kids - CROSS-SEC (prevalent cases)
mod.case.1 = sum(out["Is1"], out["Icoinf1"]) / sum(out["Is1"], out["Icoinf1"], out["Iother1"])
mod.case.2 = sum(out["Is2"], out["Icoinf2"]) / sum(out["Is2"], out["Icoinf2"], out["Iother2"])

# Calculate prevalence of norovirus infection among cases (sick) kids - COHORT (incident cases)
infec.noro = out["Is1"] + out["Isr1"] + post*(out["Ips1"] + out["Ia1"] + out["Ipsr1"] + out["Iar1"] + out["Icoinf1"]) +
  out["Is2"] + out["Isr2"] + post*(out["Ips2"] + out["Ia2"] + out["Ipsr2"] + out["Iar2"] + out["Icoinf2"])
foi.noro.1 = beta.noro.1 * (infec.noro/N)
foi.noro.2 = beta.noro.2 * (infec.noro/N)

# Calculate overall AGE incidence (per person-year) for each age group
age.inc.mod.1 = foi.noro.1*365*(out["S1"]*s.p + out["R1"]*s.r)/N.1 + 
  foi.other.1*365*out["S1"]/N.1 + 
  foi.noro.1*365*out["Iother1"]/N.1 +
  foi.other.1*365*(out["Ia1"] + out["Ips1"] +  out["Ipsr1"] + out["Iar1"])/N.1
age.inc.mod.2 = foi.noro.2*365*(out["S2"]*s.p + out["R2"]*s.r)/N.2 + 
  foi.other.2*365*out["S2"]/N.2 + 
  foi.noro.2*365*out["Iother2"]/N.2 +
  foi.other.2*365*(out["Ia2"] + out["Ips2"] +  out["Ipsr2"] + out["Iar2"])/N.2

noro.inc.mod.1 = foi.noro.1*365*((out["S1"]*s.p + out["R1"]*s.r)/N.1)
noro.inc.mod.2 = foi.noro.2*365*((out["S2"]*s.p + out["R2"]*s.r)/N.2)

pak.pt <- c(noro.inc.mod.1, noro.inc.mod.2, 
            beta.noro.1, beta.noro.2, foi.other.1, foi.other.2, s.p, s.r, 
            age.inc.mod.1, age.inc.mod.2, age.inc.1y, age.inc.2y,
            obs.asymp.1, mod.asymp.1, obs.asymp.2, mod.asymp.2, obs.case.1, mod.case.1, obs.case.2, mod.case.2, 
            obj.fx.val)

# Combine point esimate and CIs into one dataframe
#bang_all <- cbind(bang.pt, bang.low, bang.upp)

rm(noro.inc.mod.1, noro.inc.mod.2, 
   beta.noro.1, beta.noro.2, foi.other.1, foi.other.2, s.p, s.r,
   age.inc.mod.1, age.inc.mod.2, age.inc.1y, age.inc.2y,
   obs.asymp.1, mod.asymp.1, obs.asymp.2, mod.asymp.2, obs.case.1, mod.case.1, obs.case.2, mod.case.2, 
   obj.fx.val)




#######################################################################################
# BRAZIL #
#######################################################################################

obs.asymp.1.pt <- maled_forfit["Brazil","control.allno.1y"] 
obs.asymp.1.low <- maled_forfit["Brazil","ci95.control.1y.low"] 
obs.asymp.1.upp <- maled_forfit["Brazil","ci95.control.1y.upp"] 

obs.asymp.2.pt <- maled_forfit["Brazil","control.allno.2y"] 
obs.asymp.2.low <- maled_forfit["Brazil","ci95.control.2y.low"] 
obs.asymp.2.upp <- maled_forfit["Brazil","ci95.control.2y.upp"] 

obs.case.1.pt <- maled_forfit["Brazil","case.allno.1y"] 

obs.case.2.pt <- maled_forfit["Brazil","case.allno.2y"] 

age.inc.1y <- maled_forfit["Brazil","age.inc.1y.jpm"] 
age.inc.2y <- maled_forfit["Brazil","age.inc.2y.jpm"]

##############################
obs.asymp.1 <- obs.asymp.1.pt  # 0-11mo target
obs.asymp.2 <- obs.asymp.2.pt  # 12-23mo target
obs.case.1 <- obs.case.1.pt    # 0-11mo target
obs.case.2 <- obs.case.2.pt    # 12-23mo target
##############################

fres <- nloptr(x0=p0.bang,eval_f=RM_solve,lb=lb,ub=ub,opts=list("algorithm"=alg.name,xtol_abs=1e-3,maxeval=stps,maxtime=mtime,print_level=1))

fpar=(fres$solution); asymout.final=fres$objective; totiterate=fres$iterations; solverstatus=fres$status
bestfit=c(fpar, asymout.final, totiterate, solverstatus)
bestfit <- as.data.frame(bestfit)
rownames(bestfit)=c("beta.noro.1", "beta.noro.2", "foi.other.1", "foi.other.2", "s.p", "s.r", "objective function value", 
                    "total iterations", "solver status")

beta.noro.1 <- bestfit[1,1]
beta.noro.2 <- bestfit[2,1]
foi.other.1 <- bestfit[3,1]
foi.other.2 <- bestfit[4,1]
s.p <- bestfit[5,1]
s.r <- bestfit[6,1]
obj.fx.val <- bestfit[7,1]

out = runsteady(y = xstart, fun = RM, parms = c(beta.noro.1, beta.noro.2, foi.other.1, foi.other.2, s.p, s.r), times = c(0, 1e5))

out = out$y


# Calculate pop sizes
N.1 = sum(out["S1"], out["Is1"], out["Ips1"], out["Ia1"], out["R1"], out["Isr1"], out["Ipsr1"], out["Iar1"], out["Iother1"], out["Icoinf1"])
N.2 = sum(out["S2"], out["Is2"], out["Ips2"], out["Ia2"], out["R2"], out["Isr2"], out["Ipsr2"], out["Iar2"], out["Iother2"], out["Icoinf2"])
N = N.1 + N.2

# Calculate prevalence of norovirus infection among healthy (not ill) kids - CROSS-SEC (prevalent controls)
mod.asymp.1 = sum(out["Ips1"], out["Ia1"], out["Ipsr1"], out["Iar1"]) /
  sum(out["Ips1"], out["Ia1"], out["Ipsr1"], out["Iar1"], out["S1"], out["R1"]) 
mod.asymp.2 = sum(out["Ips2"], out["Ia2"], out["Ipsr2"], out["Iar2"]) /
  sum(out["Ips2"], out["Ia2"], out["Ipsr2"], out["Iar2"], out["S2"], out["R2"]) 

# Calculate prevalence of norovirus infection among cases (sick) kids - CROSS-SEC (prevalent cases)
mod.case.1 = sum(out["Is1"], out["Icoinf1"]) / sum(out["Is1"], out["Icoinf1"], out["Iother1"])
mod.case.2 = sum(out["Is2"], out["Icoinf2"]) / sum(out["Is2"], out["Icoinf2"], out["Iother2"])

# Calculate prevalence of norovirus infection among cases (sick) kids - COHORT (incident cases)
infec.noro = out["Is1"] + out["Isr1"] + post*(out["Ips1"] + out["Ia1"] + out["Ipsr1"] + out["Iar1"] + out["Icoinf1"]) +
  out["Is2"] + out["Isr2"] + post*(out["Ips2"] + out["Ia2"] + out["Ipsr2"] + out["Iar2"] + out["Icoinf2"])
foi.noro.1 = beta.noro.1 * (infec.noro/N)
foi.noro.2 = beta.noro.2 * (infec.noro/N)

# Calculate overall AGE incidence (per person-year) for each age group
age.inc.mod.1 = foi.noro.1*365*(out["S1"]*s.p + out["R1"]*s.r)/N.1 + 
  foi.other.1*365*out["S1"]/N.1 + 
  foi.noro.1*365*out["Iother1"]/N.1 +
  foi.other.1*365*(out["Ia1"] + out["Ips1"] +  out["Ipsr1"] + out["Iar1"])/N.1
age.inc.mod.2 = foi.noro.2*365*(out["S2"]*s.p + out["R2"]*s.r)/N.2 + 
  foi.other.2*365*out["S2"]/N.2 + 
  foi.noro.2*365*out["Iother2"]/N.2 +
  foi.other.2*365*(out["Ia2"] + out["Ips2"] +  out["Ipsr2"] + out["Iar2"])/N.2

noro.inc.mod.1 = foi.noro.1*365*((out["S1"]*s.p + out["R1"]*s.r)/N.1)
noro.inc.mod.2 = foi.noro.2*365*((out["S2"]*s.p + out["R2"]*s.r)/N.2)

braz.pt <- c(noro.inc.mod.1, noro.inc.mod.2, 
             beta.noro.1, beta.noro.2, foi.other.1, foi.other.2, s.p, s.r, 
             age.inc.mod.1, age.inc.mod.2, age.inc.1y, age.inc.2y,
             obs.asymp.1, mod.asymp.1, obs.asymp.2, mod.asymp.2, obs.case.1, mod.case.1, obs.case.2, mod.case.2, 
             obj.fx.val)

# Combine point esimate and CIs into one dataframe
#bang_all <- cbind(bang.pt, bang.low, bang.upp)

rm(noro.inc.mod.1, noro.inc.mod.2, 
   beta.noro.1, beta.noro.2, foi.other.1, foi.other.2, s.p, s.r,
   age.inc.mod.1, age.inc.mod.2, age.inc.1y, age.inc.2y,
   obs.asymp.1, mod.asymp.1, obs.asymp.2, mod.asymp.2, obs.case.1, mod.case.1, obs.case.2, mod.case.2,
   obj.fx.val)




#######################################################################################
# PERU #
#######################################################################################

obs.asymp.1.pt <- maled_forfit["Peru","control.allno.1y"] 
obs.asymp.1.low <- maled_forfit["Peru","ci95.control.1y.low"] 
obs.asymp.1.upp <- maled_forfit["Peru","ci95.control.1y.upp"] 

obs.asymp.2.pt <- maled_forfit["Peru","control.allno.2y"] 
obs.asymp.2.low <- maled_forfit["Peru","ci95.control.2y.low"] 
obs.asymp.2.upp <- maled_forfit["Peru","ci95.control.2y.upp"] 

obs.case.1.pt <- maled_forfit["Peru","case.allno.1y"] 

obs.case.2.pt <- maled_forfit["Peru","case.allno.2y"] 

age.inc.1y <- maled_forfit["Peru","age.inc.1y.jpm"] 
age.inc.2y <- maled_forfit["Peru","age.inc.2y.jpm"]

##############################
obs.asymp.1 <- obs.asymp.1.pt  # 0-11mo target
obs.asymp.2 <- obs.asymp.2.pt  # 12-23mo target
obs.case.1 <- obs.case.1.pt    # 0-11mo target
obs.case.2 <- obs.case.2.pt    # 12-23mo target
##############################

fres <- nloptr(x0=p0.bang,eval_f=RM_solve,lb=lb,ub=ub,opts=list("algorithm"=alg.name,xtol_abs=1e-3,maxeval=stps,maxtime=mtime,print_level=1))

fpar=(fres$solution); asymout.final=fres$objective; totiterate=fres$iterations; solverstatus=fres$status
bestfit=c(fpar, asymout.final, totiterate, solverstatus)
bestfit <- as.data.frame(bestfit)
rownames(bestfit)=c("beta.noro.1", "beta.noro.2", "foi.other.1", "foi.other.2", "s.p", "s.r", "objective function value", 
                    "total iterations", "solver status")

beta.noro.1 <- bestfit[1,1]
beta.noro.2 <- bestfit[2,1]
foi.other.1 <- bestfit[3,1]
foi.other.2 <- bestfit[4,1]
s.p <- bestfit[5,1]
s.r <- bestfit[6,1]
obj.fx.val <- bestfit[7,1]

out = runsteady(y = xstart, fun = RM, parms = c(beta.noro.1, beta.noro.2, foi.other.1, foi.other.2, s.p, s.r), times = c(0, 1e5))

out = out$y


# Calculate pop sizes
N.1 = sum(out["S1"], out["Is1"], out["Ips1"], out["Ia1"], out["R1"], out["Isr1"], out["Ipsr1"], out["Iar1"], out["Iother1"], out["Icoinf1"])
N.2 = sum(out["S2"], out["Is2"], out["Ips2"], out["Ia2"], out["R2"], out["Isr2"], out["Ipsr2"], out["Iar2"], out["Iother2"], out["Icoinf2"])
N = N.1 + N.2

# Calculate prevalence of norovirus infection among healthy (not ill) kids - CROSS-SEC (prevalent controls)
mod.asymp.1 = sum(out["Ips1"], out["Ia1"], out["Ipsr1"], out["Iar1"]) /
  sum(out["Ips1"], out["Ia1"], out["Ipsr1"], out["Iar1"], out["S1"], out["R1"]) 
mod.asymp.2 = sum(out["Ips2"], out["Ia2"], out["Ipsr2"], out["Iar2"]) /
  sum(out["Ips2"], out["Ia2"], out["Ipsr2"], out["Iar2"], out["S2"], out["R2"]) 

# Calculate prevalence of norovirus infection among cases (sick) kids - CROSS-SEC (prevalent cases)
mod.case.1 = sum(out["Is1"], out["Icoinf1"]) / sum(out["Is1"], out["Icoinf1"], out["Iother1"])
mod.case.2 = sum(out["Is2"], out["Icoinf2"]) / sum(out["Is2"], out["Icoinf2"], out["Iother2"])

# Calculate prevalence of norovirus infection among cases (sick) kids - COHORT (incident cases)
infec.noro = out["Is1"] + out["Isr1"] + post*(out["Ips1"] + out["Ia1"] + out["Ipsr1"] + out["Iar1"] + out["Icoinf1"]) +
  out["Is2"] + out["Isr2"] + post*(out["Ips2"] + out["Ia2"] + out["Ipsr2"] + out["Iar2"] + out["Icoinf2"])
foi.noro.1 = beta.noro.1 * (infec.noro/N)
foi.noro.2 = beta.noro.2 * (infec.noro/N)

# Calculate overall AGE incidence (per person-year) for each age group
age.inc.mod.1 = foi.noro.1*365*(out["S1"]*s.p + out["R1"]*s.r)/N.1 + 
  foi.other.1*365*out["S1"]/N.1 + 
  foi.noro.1*365*out["Iother1"]/N.1 +
  foi.other.1*365*(out["Ia1"] + out["Ips1"] +  out["Ipsr1"] + out["Iar1"])/N.1
age.inc.mod.2 = foi.noro.2*365*(out["S2"]*s.p + out["R2"]*s.r)/N.2 + 
  foi.other.2*365*out["S2"]/N.2 + 
  foi.noro.2*365*out["Iother2"]/N.2 +
  foi.other.2*365*(out["Ia2"] + out["Ips2"] +  out["Ipsr2"] + out["Iar2"])/N.2

noro.inc.mod.1 = foi.noro.1*365*((out["S1"]*s.p + out["R1"]*s.r)/N.1)
noro.inc.mod.2 = foi.noro.2*365*((out["S2"]*s.p + out["R2"]*s.r)/N.2)

peru.pt <- c(noro.inc.mod.1, noro.inc.mod.2, 
             beta.noro.1, beta.noro.2, foi.other.1, foi.other.2, s.p, s.r, 
             age.inc.mod.1, age.inc.mod.2, age.inc.1y, age.inc.2y,
             obs.asymp.1, mod.asymp.1, obs.asymp.2, mod.asymp.2, obs.case.1, mod.case.1, obs.case.2, mod.case.2, 
             obj.fx.val)

# Combine point esimate and CIs into one dataframe
#bang_all <- cbind(bang.pt, bang.low, bang.upp)

rm(noro.inc.mod.1, noro.inc.mod.2, 
   beta.noro.1, beta.noro.2, foi.other.1, foi.other.2, s.p, s.r,
   age.inc.mod.1, age.inc.mod.2, age.inc.1y, age.inc.2y,
   obs.asymp.1, mod.asymp.1, obs.asymp.2, mod.asymp.2, obs.case.1, mod.case.1, obs.case.2, mod.case.2,
   obj.fx.val)



#######################################################################################################


# RESULTS - MODEL-BASED INCIDENCE BY COUNTRY

results_byco <- rbind(bang.pt, ind.pt, nep.pt, pak.pt, braz.pt, peru.pt)

co <- c("Bangladesh", "India", "Nepal", "Pakistan", "Brazil", "Peru")

modelfits <- data.frame(cbind(co, data.frame(rbind(bang.pt, ind.pt, nep.pt, pak.pt, braz.pt, peru.pt))))

paramsumm.table <- knitr::kable(modelfits, digits = 5,
                                col.names = c("Country", "modeled incidence in C-Ys, <1y", "modeled incidence in C-Ys, 1-2y", 
                                              "beta.noro.1", "beta.noro.2", 
                                              "foi.other.1", "foi.other.2",
                                              "s.p", "s.r", 
                                              "age.inc.mod.1", "age.inc.mod.2", "age.inc.1y", "age.inc.2y",
                                              "obs.asymp.1", "mod.asymp.1", "obs.asymp.2", "mod.asymp.2", "obs.case.1", "mod.case.1", "obs.case.2", "mod.case.2", 
                                              "objective fxn value (lower is better fit)")) %>%
  kable_styling(bootstrap_options = "striped", font_size = 14) %>%
  scroll_box(width = "1000px", height = "500px")

kableExtra::save_kable(paramsumm.table, file = paste0("newstrucIotherspsronly_global_solv",solvertype,"_",date,".html", sep=""))

colnames(modelfits) <- c("Country", "mod.inc.1", "mod.inc.2", 
                         "beta.noro.1", "beta.noro.2", 
                         "foi.other.1", "foi.other.2",
                         "s.p", "s.r", 
                         "age.inc.mod.1", "age.inc.mod.2", "age.inc.1y", "age.inc.2y",
                         "obs.asymp.1", "mod.asymp.1", "obs.asymp.2", "mod.asymp.2", "obs.case.1", "mod.case.1", "obs.case.2", "mod.case.2", 
                         "objective fxn value (lower is better fit)")

outTable <- paste("newstrucIotherspsronly_global_solv",solvertype,"_",date,".csv", sep="")

write.csv(modelfits, outTable)

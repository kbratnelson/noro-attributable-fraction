rm(list=ls()) 

library(rootSolve)
library(nloptr)
library(tidyverse)
library(deSolve)
library(grid)
library(here)
library(binom)
library(kableExtra)
library(boot)

date=format(Sys.Date(), "%y%m%d")


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
  
  
  # Model equations
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
  
  adj.noro.inc.mod.1 = foi.noro.1*365*((out["S1"]*s.p + out["R1"]*s.r)/N.1)
  adj.noro.inc.mod.2 = foi.noro.2*365*((out["S2"]*s.p + out["R2"]*s.r)/N.2)
  
  # No AGE incidence
  #ssd = (abs(mod.asymp.1 - obs.asymp.1) + abs(mod.asymp.2 - obs.asymp.2) + abs(mod.case.1 - obs.case.1) + abs(mod.case.2 - obs.case.2))^2 # this is summary stat for one simulation - the sum squared difference between model and data
  
  # With AGE incidence - CROSS-SEC
  ssd = (abs(mod.asymp.1 - obs.asymp.1) + abs(mod.asymp.2 - obs.asymp.2) + abs(mod.case.1 - obs.case.1) + abs(mod.case.2 - obs.case.2) + (abs(age.inc.mod.1-age.inc.1y)/10) +  (abs(age.inc.mod.2-age.inc.2y)/10))^2 # this is summary stat for one simulation - the sum squared difference between model and data
  
  if (is.na(ssd)) {ssd = 1e10} #if something goes wrong, we don't want to return NA, instead return a very high value for the NLL
  
  return(ssd) 
}

# SPECIFY SOME UPPER AND LOWER BOUNDS FOR OPTIMIZATION (rest of bounds will be country-specific)
beta.noro.1.min <- 0
beta.noro.1.max <- 5
beta.noro.2.min <- 0
beta.noro.2.max <- 5

foi.other.1.min <- 0
foi.other.1.max <- 1
foi.other.2.min <- 0
foi.other.2.max <- 1

s.p.min <- 0
s.p.max <- 1
s.r.min <- 0
s.r.max <- 1

#some general settings for optimization
ploton=0
tstart=proc.time() #capture current time to measure speed
plevel=0
xerr.abs=0.01
max.local.stps=250000
stps=250000
mtime=72*60*60 # runtime per scenario (in seconds, so 1*60*60 is 1 hour)

# BOOTSTRAP SAMPLES
n.boot <- 100 # how many bootstrap samples?


#######################################################################################
# READ IN THE DATA FROM MAL-ED study to fit models to
#######################################################################################

# This data is from: Platt-Mills et al. Lancet 2011, Supplementary Table S1.
maled_forfit <- read.csv('Other/maled.cacoprev.csv')

rownames(maled_forfit) <- maled_forfit[,1]
maled_forfit <- maled_forfit[,-1]


maled_forfit$n.case.1 <- maled_forfit[,"case.allno.1y"]*maled_forfit[,"case.num.1y"]
maled_forfit$sd.case.1 <- sqrt(maled_forfit[,"case.allno.1y"]* maled_forfit[,"case.num.1y"] * ( 1 - maled_forfit[,"case.allno.1y"]))

maled_forfit$n.control.1 <- maled_forfit[,"control.allno.1y"]*maled_forfit[,"control.num.1y"]
maled_forfit$sd.control.1 <- sqrt(maled_forfit[,"control.allno.1y"]* maled_forfit[,"control.num.1y"] * ( 1 - maled_forfit[,"control.allno.1y"]))

maled_forfit$n.case.2 <- maled_forfit[,"case.allno.2y"]*maled_forfit[,"case.num.2y"]
maled_forfit$sd.case.2 <- sqrt(maled_forfit[,"case.allno.2y"]* maled_forfit[,"case.num.2y"] * ( 1 - maled_forfit[,"case.allno.2y"]))

maled_forfit$n.control.2 <- maled_forfit[,"control.allno.2y"]*maled_forfit[,"control.num.2y"]
maled_forfit$sd.control.2 <- sqrt(maled_forfit[,"control.allno.2y"]* maled_forfit[,"control.num.2y"] * ( 1 - maled_forfit[,"control.allno.2y"]))


#####################
# For BANGLADESH
#####################

# CHANGE THESE FOUR PARAMS FOR EACH COUNTRY
solvertype=21 #best
country = "Bangladesh"
p0 = c(1.67, 1.67, 0.016, 0.016, 0.05, 0.05)
bang.pt.all <- list() # each list element is fit from a bootstrap sample

#solver types
if (solvertype==21)      {        alg.name="NLOPT_LN_SBPLX"      }
if (solvertype==25)      {        alg.name="NLOPT_LN_NELDERMEAD" }


lb=c(beta.noro.1.min, beta.noro.2.min, foi.other.1.min, foi.other.2.min, s.p.min, s.r.min)
ub=c(beta.noro.1.max, beta.noro.2.max, foi.other.1.max, foi.other.2.max, s.p.max, s.r.max)

obs.asymp.1.pt <- maled_forfit[country,"control.allno.1y"] 
obs.asymp.1.boot <- rnorm(n.boot,	mean	=	maled_forfit[country, "n.control.1"] ,	sd	=	maled_forfit[country,"sd.control.1"])

obs.asymp.2.pt <- maled_forfit[country,"control.allno.2y"] 
obs.asymp.2.boot <- rnorm(n.boot,	mean	=	maled_forfit[country, "n.control.2"] ,	sd	=	maled_forfit[country,"sd.control.2"])

obs.case.1.pt <- maled_forfit[country,"case.allno.1y"] 
obs.case.1.boot <- rnorm(n.boot,	mean	=	maled_forfit[country, "n.case.1"] ,	sd	=	maled_forfit[country,"sd.case.1"])

obs.case.2.pt <- maled_forfit[country,"case.allno.2y"] 
obs.case.2.boot <- rnorm(n.boot,	mean	=	maled_forfit[country, "n.case.2"] ,	sd	=	maled_forfit[country,"sd.case.2"])

age.inc.1y <- maled_forfit[country,"age.inc.1y.jpm"] 
age.inc.2y <- maled_forfit[country,"age.inc.2y.jpm"]


for(i in 1:n.boot){ 
  
  ##############################
  obs.asymp.1 <- obs.asymp.1.boot[i]/maled_forfit[country,"control.num.1y"]  # 0-11mo target
  obs.asymp.2 <- obs.asymp.2.boot[i]/maled_forfit[country,"control.num.2y"]  # 12-23mo target
  obs.case.1 <- obs.case.1.boot[i]/maled_forfit[country,"case.num.1y"]    # 0-11mo target
  obs.case.2 <- obs.case.2.boot[i]/maled_forfit[country,"case.num.2y"]    # 12-23mo target
  ##############################

fres <- nloptr(x0=p0,eval_f=RM_solve,lb=lb,ub=ub,opts=list("algorithm"=alg.name,xtol_abs=1e-3,maxeval=stps,maxtime=mtime,print_level=1))

fpar=(fres$solution); asymout.final=fres$objective; totiterate=fres$iterations; solverstatus=fres$status
bestfit=c(fpar, asymout.final, totiterate, solverstatus)
bestfit <- as.data.frame(bestfit)
rownames(bestfit)=c("beta.noro.1", "beta.noro.2", "foi.other.1", "foi.other.2", "s.p", "s.r",
                    "objective function value", 
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

adj.noro.inc.mod.1 = foi.noro.1*365*((out["S1"]*s.p + out["R1"]*s.r)/N.1)
adj.noro.inc.mod.2 = foi.noro.2*365*((out["S2"]*s.p + out["R2"]*s.r)/N.2)

# Now, for validaation, use the model to calculate the AI using the OR approach - for both age groups 

#0-11mo
or.orapp.cross.1 = ((out["Is1"] + out["Isr1"] + out["Icoinf1"]) *                     
                      (out["S1"] + out["R1"])) /                          
  ((out["Ia1"] + out["Iar1"] + out["Ips1"] + out["Ipsr1"]) *                         
     out["Iother1"])   

or.orapp.nopostsym.cross.1 = ((out["Is1"] + out["Isr1"] + out["Icoinf1"]) *                     
                                (out["S1"] + out["R1"])) /                          
  ((out["Ia1"] + out["Iar1"]) *                         
     out["Iother1"])  

prev.nov.case.1 =  (out["Is1"] + out["Isr1"] + out["Icoinf1"]) /                                                                
  (out["Is1"] + out["Isr1"] + out["Icoinf1"] + out["Iother1"])

af.orapp.cross.1 = prev.nov.case.1 * (1 - (1/or.orapp.cross.1)) 
af.orapp.nopostsym.cross.1 = prev.nov.case.1 * (1 - (1/or.orapp.nopostsym.cross.1))

inc.orapp.cross.1 = af.orapp.cross.1*age.inc.1y
inc.orapp.nopostsym.cross.1 = af.orapp.cross.1*age.inc.1y


#12-23mo
or.orapp.cross.2 = ((out["Is2"] + out["Isr2"] + out["Icoinf2"]) *                     
                      (out["S2"] + out["R2"])) /                          
  ((out["Ia2"] + out["Iar2"] + out["Ips2"] + out["Ipsr2"]) *                         
     out["Iother2"])   

or.orapp.nopostsym.cross.2 = ((out["Is2"] + out["Isr2"] + out["Icoinf2"]) *                     
                                (out["S2"] + out["R2"])) /                          
  ((out["Ia2"] + out["Iar2"]) *                         
     out["Iother2"])  

prev.nov.case.2 =  (out["Is2"] + out["Isr2"] + out["Icoinf2"]) /                                                                
  (out["Is2"] + out["Isr2"] + out["Icoinf2"] + out["Iother2"])

af.orapp.cross.2 = prev.nov.case.2 * (1 - (1/or.orapp.cross.2)) 
af.orapp.nopostsym.cross.2 = prev.nov.case.2 * (1 - (1/or.orapp.nopostsym.cross.2))

inc.orapp.cross.2 = af.orapp.cross.2*age.inc.2y
inc.orapp.nopostsym.cross.2 = af.orapp.cross.2*age.inc.2y


#Calculate coinfection to report 
prop.coinf.cross.1 =  out["Icoinf1"]/(out["Icoinf1"] + out["Is1"] + out["Isr1"])
prop.coinf.cross.2 =  out["Icoinf2"]/(out["Icoinf2"] + out["Is2"] + out["Isr2"])


#Bring in file where I've already calculated incidence using the DE approach. Wasa originally using this as 'limits' for optimization procedure, but now just want to report the incidence under DE approach
limits_boundedfit <- read.csv(here("Other/limits_boundedfit.csv")) 
de.inc <- filter(limits_boundedfit, index.co==country) %>% select("norocaseprev.ai")
de.inc.1 <- de.inc[1,]
de.inc.2 <- de.inc[2,]

pt <- c(adj.noro.inc.mod.1, adj.noro.inc.mod.2, 
             de.inc.1, de.inc.2,
             inc.orapp.cross.1, inc.orapp.cross.2,
             
             prop.coinf.cross.1, prop.coinf.cross.2,
             
             beta.noro.1, beta.noro.2, foi.other.1, foi.other.2, s.p, s.r, 
             age.inc.mod.1, age.inc.mod.2, age.inc.1y, age.inc.2y,
             obs.asymp.1, mod.asymp.1, obs.asymp.2, mod.asymp.2, obs.case.1, mod.case.1, obs.case.2, mod.case.2, 
             
             or.orapp.cross.1, or.orapp.nopostsym.cross.1,
             or.orapp.cross.2, or.orapp.nopostsym.cross.2,
             af.orapp.cross.1, af.orapp.nopostsym.cross.1,
             af.orapp.cross.2, af.orapp.nopostsym.cross.2,
             inc.orapp.nopostsym.cross.1, inc.orapp.nopostsym.cross.2, 

             obj.fx.val)

bang.pt.all[[i]] <- pt

#remove objects before running model for next country

rm(adj.noro.inc.mod.1, adj.noro.inc.mod.2, 
   de.inc.1, de.inc.2,
   inc.orapp.cross.1, inc.orapp.cross.2,
   
   prop.coinf.cross.1, prop.coinf.cross.2,
   
   beta.noro.1, beta.noro.2, foi.other.1, foi.other.2, s.p, s.r, 
   age.inc.mod.1, age.inc.mod.2, 
   obs.asymp.1, mod.asymp.1, obs.asymp.2, mod.asymp.2, obs.case.1, mod.case.1, obs.case.2, mod.case.2, 
   
   or.orapp.cross.1, or.orapp.nopostsym.cross.1,
   or.orapp.cross.2, or.orapp.nopostsym.cross.2,
   af.orapp.cross.1, af.orapp.nopostsym.cross.1,
   af.orapp.cross.2, af.orapp.nopostsym.cross.2,
   inc.orapp.nopostsym.cross.1, inc.orapp.nopostsym.cross.2, 
   
   obj.fx.val)

}


#####################
# For INDIA
#####################
# CHANGE THESE FOUR PARAMS FOR EACH COUNTRY
solvertype=21 #best
country = "India"
p0 = c(1.17, 2.25, 0.01, 0.006, 0.06, 0.06) 
ind.pt.all <- list() # each list element is fit from a bootstrap sample

#solver types
if (solvertype==21)      {        alg.name="NLOPT_LN_SBPLX"      }
if (solvertype==25)      {        alg.name="NLOPT_LN_NELDERMEAD" }


lb=c(beta.noro.1.min, beta.noro.2.min, foi.other.1.min, foi.other.2.min, s.p.min, s.r.min)
ub=c(beta.noro.1.max, beta.noro.2.max, foi.other.1.max, foi.other.2.max, s.p.max, s.r.max)

obs.asymp.1.pt <- maled_forfit[country,"control.allno.1y"] 
obs.asymp.1.boot <- rnorm(n.boot,	mean	=	maled_forfit[country, "n.control.1"] ,	sd	=	maled_forfit[country,"sd.control.1"])

obs.asymp.2.pt <- maled_forfit[country,"control.allno.2y"] 
obs.asymp.2.boot <- rnorm(n.boot,	mean	=	maled_forfit[country, "n.control.2"] ,	sd	=	maled_forfit[country,"sd.control.2"])

obs.case.1.pt <- maled_forfit[country,"case.allno.1y"] 
obs.case.1.boot <- rnorm(n.boot,	mean	=	maled_forfit[country, "n.case.1"] ,	sd	=	maled_forfit[country,"sd.case.1"])

obs.case.2.pt <- maled_forfit[country,"case.allno.2y"] 
obs.case.2.boot <- rnorm(n.boot,	mean	=	maled_forfit[country, "n.case.2"] ,	sd	=	maled_forfit[country,"sd.case.2"])

age.inc.1y <- maled_forfit[country,"age.inc.1y.jpm"] 
age.inc.2y <- maled_forfit[country,"age.inc.2y.jpm"]


for(i in 1:n.boot){ 
  
  ##############################
  obs.asymp.1 <- obs.asymp.1.boot[i]/maled_forfit[country,"control.num.1y"]  # 0-11mo target
  obs.asymp.2 <- obs.asymp.2.boot[i]/maled_forfit[country,"control.num.2y"]  # 12-23mo target
  obs.case.1 <- obs.case.1.boot[i]/maled_forfit[country,"case.num.1y"]    # 0-11mo target
  obs.case.2 <- obs.case.2.boot[i]/maled_forfit[country,"case.num.2y"]    # 12-23mo target
  ##############################
  
  
  fres <- nloptr(x0=p0,eval_f=RM_solve,lb=lb,ub=ub,opts=list("algorithm"=alg.name,xtol_abs=1e-3,maxeval=stps,maxtime=mtime,print_level=1))
  
  fpar=(fres$solution); asymout.final=fres$objective; totiterate=fres$iterations; solverstatus=fres$status
  bestfit=c(fpar, asymout.final, totiterate, solverstatus)
  bestfit <- as.data.frame(bestfit)
  rownames(bestfit)=c("beta.noro.1", "beta.noro.2", "foi.other.1", "foi.other.2", "s.p", "s.r",
                      "objective function value", 
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
  
  adj.noro.inc.mod.1 = foi.noro.1*365*((out["S1"]*s.p + out["R1"]*s.r)/N.1)
  adj.noro.inc.mod.2 = foi.noro.2*365*((out["S2"]*s.p + out["R2"]*s.r)/N.2)
  
  # Now, for validaation, use the model to calculate the AI using the OR approach - for both age groups 
  
  #0-11mo
  or.orapp.cross.1 = ((out["Is1"] + out["Isr1"] + out["Icoinf1"]) *                     
                        (out["S1"] + out["R1"])) /                          
    ((out["Ia1"] + out["Iar1"] + out["Ips1"] + out["Ipsr1"]) *                         
       out["Iother1"])   
  
  or.orapp.nopostsym.cross.1 = ((out["Is1"] + out["Isr1"] + out["Icoinf1"]) *                     
                                  (out["S1"] + out["R1"])) /                          
    ((out["Ia1"] + out["Iar1"]) *                         
       out["Iother1"])  
  
  prev.nov.case.1 =  (out["Is1"] + out["Isr1"] + out["Icoinf1"]) /                                                                
    (out["Is1"] + out["Isr1"] + out["Icoinf1"] + out["Iother1"])
  
  af.orapp.cross.1 = prev.nov.case.1 * (1 - (1/or.orapp.cross.1)) 
  af.orapp.nopostsym.cross.1 = prev.nov.case.1 * (1 - (1/or.orapp.nopostsym.cross.1))
  
  inc.orapp.cross.1 = af.orapp.cross.1*age.inc.1y
  inc.orapp.nopostsym.cross.1 = af.orapp.cross.1*age.inc.1y
  
  
  #12-23mo
  or.orapp.cross.2 = ((out["Is2"] + out["Isr2"] + out["Icoinf2"]) *                     
                        (out["S2"] + out["R2"])) /                          
    ((out["Ia2"] + out["Iar2"] + out["Ips2"] + out["Ipsr2"]) *                         
       out["Iother2"])   
  
  or.orapp.nopostsym.cross.2 = ((out["Is2"] + out["Isr2"] + out["Icoinf2"]) *                     
                                  (out["S2"] + out["R2"])) /                          
    ((out["Ia2"] + out["Iar2"]) *                         
       out["Iother2"])  
  
  prev.nov.case.2 =  (out["Is2"] + out["Isr2"] + out["Icoinf2"]) /                                                                
    (out["Is2"] + out["Isr2"] + out["Icoinf2"] + out["Iother2"])
  
  af.orapp.cross.2 = prev.nov.case.2 * (1 - (1/or.orapp.cross.2)) 
  af.orapp.nopostsym.cross.2 = prev.nov.case.2 * (1 - (1/or.orapp.nopostsym.cross.2))
  
  inc.orapp.cross.2 = af.orapp.cross.2*age.inc.2y
  inc.orapp.nopostsym.cross.2 = af.orapp.cross.2*age.inc.2y
  
  
  #Calculate coinfection to report 
  prop.coinf.cross.1 =  out["Icoinf1"]/(out["Icoinf1"] + out["Is1"] + out["Isr1"])
  prop.coinf.cross.2 =  out["Icoinf2"]/(out["Icoinf2"] + out["Is2"] + out["Isr2"])
  
  
  #Bring in file where I've already calculated incidence using the DE approach. Wasa originally using this as 'limits' for optimization procedure, but now just want to report the incidence under DE approach
  limits_boundedfit <- read.csv(here("Other/limits_boundedfit.csv")) 
  de.inc <- filter(limits_boundedfit, index.co==country) %>% select("norocaseprev.ai")
  de.inc.1 <- de.inc[1,]
  de.inc.2 <- de.inc[2,]
  
  pt <- c(adj.noro.inc.mod.1, adj.noro.inc.mod.2, 
          de.inc.1, de.inc.2,
          inc.orapp.cross.1, inc.orapp.cross.2,
          
          prop.coinf.cross.1, prop.coinf.cross.2,
          
          beta.noro.1, beta.noro.2, foi.other.1, foi.other.2, s.p, s.r, 
          age.inc.mod.1, age.inc.mod.2, age.inc.1y, age.inc.2y,
          obs.asymp.1, mod.asymp.1, obs.asymp.2, mod.asymp.2, obs.case.1, mod.case.1, obs.case.2, mod.case.2, 
          
          or.orapp.cross.1, or.orapp.nopostsym.cross.1,
          or.orapp.cross.2, or.orapp.nopostsym.cross.2,
          af.orapp.cross.1, af.orapp.nopostsym.cross.1,
          af.orapp.cross.2, af.orapp.nopostsym.cross.2,
          inc.orapp.nopostsym.cross.1, inc.orapp.nopostsym.cross.2, 
          
          obj.fx.val)
  
  ind.pt.all[[i]] <- pt
  
  #remove objects before running model for next country
  
  rm(adj.noro.inc.mod.1, adj.noro.inc.mod.2, 
     de.inc.1, de.inc.2,
     inc.orapp.cross.1, inc.orapp.cross.2,
     
     prop.coinf.cross.1, prop.coinf.cross.2,
     
     beta.noro.1, beta.noro.2, foi.other.1, foi.other.2, s.p, s.r, 
     age.inc.mod.1, age.inc.mod.2,
     obs.asymp.1, mod.asymp.1, obs.asymp.2, mod.asymp.2, obs.case.1, mod.case.1, obs.case.2, mod.case.2, 
     
     or.orapp.cross.1, or.orapp.nopostsym.cross.1,
     or.orapp.cross.2, or.orapp.nopostsym.cross.2,
     af.orapp.cross.1, af.orapp.nopostsym.cross.1,
     af.orapp.cross.2, af.orapp.nopostsym.cross.2,
     inc.orapp.nopostsym.cross.1, inc.orapp.nopostsym.cross.2, 
     
     obj.fx.val)
  
}



#####################
# For NEPAL
#####################
# CHANGE THESE FOUR PARAMS FOR EACH COUNTRY
solvertype=25 #best
country = "Nepal"
p0 = c(1.67, 0.56, 0.01, 0.008, 0.15, 0.05)
nep.pt.all <- list() # each list element is fit from a bootstrap sample

#solver types
if (solvertype==21)      {        alg.name="NLOPT_LN_SBPLX"      }
if (solvertype==25)      {        alg.name="NLOPT_LN_NELDERMEAD" }


lb=c(beta.noro.1.min, beta.noro.2.min, foi.other.1.min, foi.other.2.min, s.p.min, s.r.min)
ub=c(beta.noro.1.max, beta.noro.2.max, foi.other.1.max, foi.other.2.max, s.p.max, s.r.max)

obs.asymp.1.pt <- maled_forfit[country,"control.allno.1y"] 
obs.asymp.1.boot <- rnorm(n.boot,	mean	=	maled_forfit[country, "n.control.1"] ,	sd	=	maled_forfit[country,"sd.control.1"])

obs.asymp.2.pt <- maled_forfit[country,"control.allno.2y"] 
obs.asymp.2.boot <- rnorm(n.boot,	mean	=	maled_forfit[country, "n.control.2"] ,	sd	=	maled_forfit[country,"sd.control.2"])

obs.case.1.pt <- maled_forfit[country,"case.allno.1y"] 
obs.case.1.boot <- rnorm(n.boot,	mean	=	maled_forfit[country, "n.case.1"] ,	sd	=	maled_forfit[country,"sd.case.1"])

obs.case.2.pt <- maled_forfit[country,"case.allno.2y"] 
obs.case.2.boot <- rnorm(n.boot,	mean	=	maled_forfit[country, "n.case.2"] ,	sd	=	maled_forfit[country,"sd.case.2"])

age.inc.1y <- maled_forfit[country,"age.inc.1y.jpm"] 
age.inc.2y <- maled_forfit[country,"age.inc.2y.jpm"]


for(i in 1:n.boot){ 
  
  ##############################
  obs.asymp.1 <- obs.asymp.1.boot[i]/maled_forfit[country,"control.num.1y"]  # 0-11mo target
  obs.asymp.2 <- obs.asymp.2.boot[i]/maled_forfit[country,"control.num.2y"]  # 12-23mo target
  obs.case.1 <- obs.case.1.boot[i]/maled_forfit[country,"case.num.1y"]    # 0-11mo target
  obs.case.2 <- obs.case.2.boot[i]/maled_forfit[country,"case.num.2y"]    # 12-23mo target
  ##############################
  
  fres <- nloptr(x0=p0,eval_f=RM_solve,lb=lb,ub=ub,opts=list("algorithm"=alg.name,xtol_abs=1e-3,maxeval=stps,maxtime=mtime,print_level=1))
  
  fpar=(fres$solution); asymout.final=fres$objective; totiterate=fres$iterations; solverstatus=fres$status
  bestfit=c(fpar, asymout.final, totiterate, solverstatus)
  bestfit <- as.data.frame(bestfit)
  rownames(bestfit)=c("beta.noro.1", "beta.noro.2", "foi.other.1", "foi.other.2", "s.p", "s.r",
                      "objective function value", 
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
  
  adj.noro.inc.mod.1 = foi.noro.1*365*((out["S1"]*s.p + out["R1"]*s.r)/N.1)
  adj.noro.inc.mod.2 = foi.noro.2*365*((out["S2"]*s.p + out["R2"]*s.r)/N.2)
  
  # Now, for validaation, use the model to calculate the AI using the OR approach - for both age groups 
  
  #0-11mo
  or.orapp.cross.1 = ((out["Is1"] + out["Isr1"] + out["Icoinf1"]) *                     
                        (out["S1"] + out["R1"])) /                          
    ((out["Ia1"] + out["Iar1"] + out["Ips1"] + out["Ipsr1"]) *                         
       out["Iother1"])   
  
  or.orapp.nopostsym.cross.1 = ((out["Is1"] + out["Isr1"] + out["Icoinf1"]) *                     
                                  (out["S1"] + out["R1"])) /                          
    ((out["Ia1"] + out["Iar1"]) *                         
       out["Iother1"])  
  
  prev.nov.case.1 =  (out["Is1"] + out["Isr1"] + out["Icoinf1"]) /                                                                
    (out["Is1"] + out["Isr1"] + out["Icoinf1"] + out["Iother1"])
  
  af.orapp.cross.1 = prev.nov.case.1 * (1 - (1/or.orapp.cross.1)) 
  af.orapp.nopostsym.cross.1 = prev.nov.case.1 * (1 - (1/or.orapp.nopostsym.cross.1))
  
  inc.orapp.cross.1 = af.orapp.cross.1*age.inc.1y
  inc.orapp.nopostsym.cross.1 = af.orapp.cross.1*age.inc.1y
  
  
  #12-23mo
  or.orapp.cross.2 = ((out["Is2"] + out["Isr2"] + out["Icoinf2"]) *                     
                        (out["S2"] + out["R2"])) /                          
    ((out["Ia2"] + out["Iar2"] + out["Ips2"] + out["Ipsr2"]) *                         
       out["Iother2"])   
  
  or.orapp.nopostsym.cross.2 = ((out["Is2"] + out["Isr2"] + out["Icoinf2"]) *                     
                                  (out["S2"] + out["R2"])) /                          
    ((out["Ia2"] + out["Iar2"]) *                         
       out["Iother2"])  
  
  prev.nov.case.2 =  (out["Is2"] + out["Isr2"] + out["Icoinf2"]) /                                                                
    (out["Is2"] + out["Isr2"] + out["Icoinf2"] + out["Iother2"])
  
  af.orapp.cross.2 = prev.nov.case.2 * (1 - (1/or.orapp.cross.2)) 
  af.orapp.nopostsym.cross.2 = prev.nov.case.2 * (1 - (1/or.orapp.nopostsym.cross.2))
  
  inc.orapp.cross.2 = af.orapp.cross.2*age.inc.2y
  inc.orapp.nopostsym.cross.2 = af.orapp.cross.2*age.inc.2y
  
  
  #Calculate coinfection to report 
  prop.coinf.cross.1 =  out["Icoinf1"]/(out["Icoinf1"] + out["Is1"] + out["Isr1"])
  prop.coinf.cross.2 =  out["Icoinf2"]/(out["Icoinf2"] + out["Is2"] + out["Isr2"])
  
  
  #Bring in file where I've already calculated incidence using the DE approach. Wasa originally using this as 'limits' for optimization procedure, but now just want to report the incidence under DE approach
  limits_boundedfit <- read.csv(here("Other/limits_boundedfit.csv")) 
  de.inc <- filter(limits_boundedfit, index.co==country) %>% select("norocaseprev.ai")
  de.inc.1 <- de.inc[1,]
  de.inc.2 <- de.inc[2,]
  
  pt <- c(adj.noro.inc.mod.1, adj.noro.inc.mod.2, 
          de.inc.1, de.inc.2,
          inc.orapp.cross.1, inc.orapp.cross.2,
          
          prop.coinf.cross.1, prop.coinf.cross.2,
          
          beta.noro.1, beta.noro.2, foi.other.1, foi.other.2, s.p, s.r, 
          age.inc.mod.1, age.inc.mod.2, age.inc.1y, age.inc.2y,
          obs.asymp.1, mod.asymp.1, obs.asymp.2, mod.asymp.2, obs.case.1, mod.case.1, obs.case.2, mod.case.2, 
          
          or.orapp.cross.1, or.orapp.nopostsym.cross.1,
          or.orapp.cross.2, or.orapp.nopostsym.cross.2,
          af.orapp.cross.1, af.orapp.nopostsym.cross.1,
          af.orapp.cross.2, af.orapp.nopostsym.cross.2,
          inc.orapp.nopostsym.cross.1, inc.orapp.nopostsym.cross.2, 
          
          obj.fx.val)
  
  nep.pt.all[[i]] <- pt
  
  #remove objects before running model for next country
  
  rm(adj.noro.inc.mod.1, adj.noro.inc.mod.2, 
     de.inc.1, de.inc.2,
     inc.orapp.cross.1, inc.orapp.cross.2,
     
     prop.coinf.cross.1, prop.coinf.cross.2,
     
     beta.noro.1, beta.noro.2, foi.other.1, foi.other.2, s.p, s.r, 
     age.inc.mod.1, age.inc.mod.2, 
     obs.asymp.1, mod.asymp.1, obs.asymp.2, mod.asymp.2, obs.case.1, mod.case.1, obs.case.2, mod.case.2, 
     
     or.orapp.cross.1, or.orapp.nopostsym.cross.1,
     or.orapp.cross.2, or.orapp.nopostsym.cross.2,
     af.orapp.cross.1, af.orapp.nopostsym.cross.1,
     af.orapp.cross.2, af.orapp.nopostsym.cross.2,
     inc.orapp.nopostsym.cross.1, inc.orapp.nopostsym.cross.2, 
     
     obj.fx.val)
  
}




#####################
# For PAKISTAN
#####################
# CHANGE THESE FOUR PARAMS FOR EACH COUNTRY
solvertype=21 #best
country = "Pakistan"
p0 = c(2.00, 0.71, 0.03, 0.02, 0.17, 0.06)
pak.pt.all <- list() # each list element is fit from a bootstrap sample

#solver types
if (solvertype==21)      {        alg.name="NLOPT_LN_SBPLX"      }
if (solvertype==25)      {        alg.name="NLOPT_LN_NELDERMEAD" }


lb=c(beta.noro.1.min, beta.noro.2.min, foi.other.1.min, foi.other.2.min, s.p.min, s.r.min)
ub=c(beta.noro.1.max, beta.noro.2.max, foi.other.1.max, foi.other.2.max, s.p.max, s.r.max)

obs.asymp.1.pt <- maled_forfit[country,"control.allno.1y"] 
obs.asymp.1.boot <- rnorm(n.boot,	mean	=	maled_forfit[country, "n.control.1"] ,	sd	=	maled_forfit[country,"sd.control.1"])

obs.asymp.2.pt <- maled_forfit[country,"control.allno.2y"] 
obs.asymp.2.boot <- rnorm(n.boot,	mean	=	maled_forfit[country, "n.control.2"] ,	sd	=	maled_forfit[country,"sd.control.2"])

obs.case.1.pt <- maled_forfit[country,"case.allno.1y"] 
obs.case.1.boot <- rnorm(n.boot,	mean	=	maled_forfit[country, "n.case.1"] ,	sd	=	maled_forfit[country,"sd.case.1"])

obs.case.2.pt <- maled_forfit[country,"case.allno.2y"] 
obs.case.2.boot <- rnorm(n.boot,	mean	=	maled_forfit[country, "n.case.2"] ,	sd	=	maled_forfit[country,"sd.case.2"])

age.inc.1y <- maled_forfit[country,"age.inc.1y.jpm"] 
age.inc.2y <- maled_forfit[country,"age.inc.2y.jpm"]


for(i in 1:n.boot){ 
  
  ##############################
  obs.asymp.1 <- obs.asymp.1.boot[i]/maled_forfit[country,"control.num.1y"]  # 0-11mo target
  obs.asymp.2 <- obs.asymp.2.boot[i]/maled_forfit[country,"control.num.2y"]  # 12-23mo target
  obs.case.1 <- obs.case.1.boot[i]/maled_forfit[country,"case.num.1y"]    # 0-11mo target
  obs.case.2 <- obs.case.2.boot[i]/maled_forfit[country,"case.num.2y"]    # 12-23mo target
  ##############################
  
  fres <- nloptr(x0=p0,eval_f=RM_solve,lb=lb,ub=ub,opts=list("algorithm"=alg.name,xtol_abs=1e-3,maxeval=stps,maxtime=mtime,print_level=1))
  
  fpar=(fres$solution); asymout.final=fres$objective; totiterate=fres$iterations; solverstatus=fres$status
  bestfit=c(fpar, asymout.final, totiterate, solverstatus)
  bestfit <- as.data.frame(bestfit)
  rownames(bestfit)=c("beta.noro.1", "beta.noro.2", "foi.other.1", "foi.other.2", "s.p", "s.r",
                      "objective function value", 
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
  
  adj.noro.inc.mod.1 = foi.noro.1*365*((out["S1"]*s.p + out["R1"]*s.r)/N.1)
  adj.noro.inc.mod.2 = foi.noro.2*365*((out["S2"]*s.p + out["R2"]*s.r)/N.2)
  
  # Now, for validaation, use the model to calculate the AI using the OR approach - for both age groups 
  
  #0-11mo
  or.orapp.cross.1 = ((out["Is1"] + out["Isr1"] + out["Icoinf1"]) *                     
                        (out["S1"] + out["R1"])) /                          
    ((out["Ia1"] + out["Iar1"] + out["Ips1"] + out["Ipsr1"]) *                         
       out["Iother1"])   
  
  or.orapp.nopostsym.cross.1 = ((out["Is1"] + out["Isr1"] + out["Icoinf1"]) *                     
                                  (out["S1"] + out["R1"])) /                          
    ((out["Ia1"] + out["Iar1"]) *                         
       out["Iother1"])  
  
  prev.nov.case.1 =  (out["Is1"] + out["Isr1"] + out["Icoinf1"]) /                                                                
    (out["Is1"] + out["Isr1"] + out["Icoinf1"] + out["Iother1"])
  
  af.orapp.cross.1 = prev.nov.case.1 * (1 - (1/or.orapp.cross.1)) 
  af.orapp.nopostsym.cross.1 = prev.nov.case.1 * (1 - (1/or.orapp.nopostsym.cross.1))
  
  inc.orapp.cross.1 = af.orapp.cross.1*age.inc.1y
  inc.orapp.nopostsym.cross.1 = af.orapp.cross.1*age.inc.1y
  
  
  #12-23mo
  or.orapp.cross.2 = ((out["Is2"] + out["Isr2"] + out["Icoinf2"]) *                     
                        (out["S2"] + out["R2"])) /                          
    ((out["Ia2"] + out["Iar2"] + out["Ips2"] + out["Ipsr2"]) *                         
       out["Iother2"])   
  
  or.orapp.nopostsym.cross.2 = ((out["Is2"] + out["Isr2"] + out["Icoinf2"]) *                     
                                  (out["S2"] + out["R2"])) /                          
    ((out["Ia2"] + out["Iar2"]) *                         
       out["Iother2"])  
  
  prev.nov.case.2 =  (out["Is2"] + out["Isr2"] + out["Icoinf2"]) /                                                                
    (out["Is2"] + out["Isr2"] + out["Icoinf2"] + out["Iother2"])
  
  af.orapp.cross.2 = prev.nov.case.2 * (1 - (1/or.orapp.cross.2)) 
  af.orapp.nopostsym.cross.2 = prev.nov.case.2 * (1 - (1/or.orapp.nopostsym.cross.2))
  
  inc.orapp.cross.2 = af.orapp.cross.2*age.inc.2y
  inc.orapp.nopostsym.cross.2 = af.orapp.cross.2*age.inc.2y
  
  
  #Calculate coinfection to report 
  prop.coinf.cross.1 =  out["Icoinf1"]/(out["Icoinf1"] + out["Is1"] + out["Isr1"])
  prop.coinf.cross.2 =  out["Icoinf2"]/(out["Icoinf2"] + out["Is2"] + out["Isr2"])
  
  
  #Bring in file where I've already calculated incidence using the DE approach. Wasa originally using this as 'limits' for optimization procedure, but now just want to report the incidence under DE approach
  limits_boundedfit <- read.csv(here("Other/limits_boundedfit.csv")) 
  de.inc <- filter(limits_boundedfit, index.co==country) %>% select("norocaseprev.ai")
  de.inc.1 <- de.inc[1,]
  de.inc.2 <- de.inc[2,]
  
  pt <- c(adj.noro.inc.mod.1, adj.noro.inc.mod.2, 
          de.inc.1, de.inc.2,
          inc.orapp.cross.1, inc.orapp.cross.2,
          
          prop.coinf.cross.1, prop.coinf.cross.2,
          
          beta.noro.1, beta.noro.2, foi.other.1, foi.other.2, s.p, s.r, 
          age.inc.mod.1, age.inc.mod.2, age.inc.1y, age.inc.2y,
          obs.asymp.1, mod.asymp.1, obs.asymp.2, mod.asymp.2, obs.case.1, mod.case.1, obs.case.2, mod.case.2, 
          
          or.orapp.cross.1, or.orapp.nopostsym.cross.1,
          or.orapp.cross.2, or.orapp.nopostsym.cross.2,
          af.orapp.cross.1, af.orapp.nopostsym.cross.1,
          af.orapp.cross.2, af.orapp.nopostsym.cross.2,
          inc.orapp.nopostsym.cross.1, inc.orapp.nopostsym.cross.2, 
          
          obj.fx.val)
  
  pak.pt.all[[i]] <- pt
  
  #remove objects before running model for next country
  
  rm(adj.noro.inc.mod.1, adj.noro.inc.mod.2, 
     de.inc.1, de.inc.2,
     inc.orapp.cross.1, inc.orapp.cross.2,
     
     prop.coinf.cross.1, prop.coinf.cross.2,
     
     beta.noro.1, beta.noro.2, foi.other.1, foi.other.2, s.p, s.r, 
     age.inc.mod.1, age.inc.mod.2, 
     obs.asymp.1, mod.asymp.1, obs.asymp.2, mod.asymp.2, obs.case.1, mod.case.1, obs.case.2, mod.case.2, 
     
     or.orapp.cross.1, or.orapp.nopostsym.cross.1,
     or.orapp.cross.2, or.orapp.nopostsym.cross.2,
     af.orapp.cross.1, af.orapp.nopostsym.cross.1,
     af.orapp.cross.2, af.orapp.nopostsym.cross.2,
     inc.orapp.nopostsym.cross.1, inc.orapp.nopostsym.cross.2, 
     
     obj.fx.val)
  
}



#####################
# For BRAZIL
#####################
# CHANGE THESE FOUR PARAMS FOR EACH COUNTRY
solvertype=25 #best
country = "Brazil"
p0 = c(0.56, 2.22, 0.002, 0.002, 0.15, 0.05)
braz.pt.all <- list() # each list element is fit from a bootstrap sample

#solver types
if (solvertype==21)      {        alg.name="NLOPT_LN_SBPLX"      }
if (solvertype==25)      {        alg.name="NLOPT_LN_NELDERMEAD" }


lb=c(beta.noro.1.min, beta.noro.2.min, foi.other.1.min, foi.other.2.min, s.p.min, s.r.min)
ub=c(beta.noro.1.max, beta.noro.2.max, foi.other.1.max, foi.other.2.max, s.p.max, s.r.max)

obs.asymp.1.pt <- maled_forfit[country,"control.allno.1y"] 
obs.asymp.1.boot <- rnorm(n.boot,	mean	=	maled_forfit[country, "n.control.1"] ,	sd	=	maled_forfit[country,"sd.control.1"])

obs.asymp.2.pt <- maled_forfit[country,"control.allno.2y"] 
obs.asymp.2.boot <- rnorm(n.boot,	mean	=	maled_forfit[country, "n.control.2"] ,	sd	=	maled_forfit[country,"sd.control.2"])

obs.case.1.pt <- maled_forfit[country,"case.allno.1y"] 
obs.case.1.boot <- rnorm(n.boot,	mean	=	maled_forfit[country, "n.case.1"] ,	sd	=	maled_forfit[country,"sd.case.1"])

obs.case.2.pt <- maled_forfit[country,"case.allno.2y"] 
obs.case.2.boot <- rnorm(n.boot,	mean	=	maled_forfit[country, "n.case.2"] ,	sd	=	maled_forfit[country,"sd.case.2"])

age.inc.1y <- maled_forfit[country,"age.inc.1y.jpm"] 
age.inc.2y <- maled_forfit[country,"age.inc.2y.jpm"]


for(i in 1:n.boot){ 
  
  ##############################
  obs.asymp.1 <- obs.asymp.1.boot[i]/maled_forfit[country,"control.num.1y"]  # 0-11mo target
  obs.asymp.2 <- obs.asymp.2.boot[i]/maled_forfit[country,"control.num.2y"]  # 12-23mo target
  obs.case.1 <- obs.case.1.boot[i]/maled_forfit[country,"case.num.1y"]    # 0-11mo target
  obs.case.2 <- obs.case.2.boot[i]/maled_forfit[country,"case.num.2y"]    # 12-23mo target
  ##############################
  
  fres <- nloptr(x0=p0,eval_f=RM_solve,lb=lb,ub=ub,opts=list("algorithm"=alg.name,xtol_abs=1e-3,maxeval=stps,maxtime=mtime,print_level=1))
  
  fpar=(fres$solution); asymout.final=fres$objective; totiterate=fres$iterations; solverstatus=fres$status
  bestfit=c(fpar, asymout.final, totiterate, solverstatus)
  bestfit <- as.data.frame(bestfit)
  rownames(bestfit)=c("beta.noro.1", "beta.noro.2", "foi.other.1", "foi.other.2", "s.p", "s.r",
                      "objective function value", 
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
  
  adj.noro.inc.mod.1 = foi.noro.1*365*((out["S1"]*s.p + out["R1"]*s.r)/N.1)
  adj.noro.inc.mod.2 = foi.noro.2*365*((out["S2"]*s.p + out["R2"]*s.r)/N.2)
  
  # Now, for validaation, use the model to calculate the AI using the OR approach - for both age groups 
  
  #0-11mo
  or.orapp.cross.1 = ((out["Is1"] + out["Isr1"] + out["Icoinf1"]) *                     
                        (out["S1"] + out["R1"])) /                          
    ((out["Ia1"] + out["Iar1"] + out["Ips1"] + out["Ipsr1"]) *                         
       out["Iother1"])   
  
  or.orapp.nopostsym.cross.1 = ((out["Is1"] + out["Isr1"] + out["Icoinf1"]) *                     
                                  (out["S1"] + out["R1"])) /                          
    ((out["Ia1"] + out["Iar1"]) *                         
       out["Iother1"])  
  
  prev.nov.case.1 =  (out["Is1"] + out["Isr1"] + out["Icoinf1"]) /                                                                
    (out["Is1"] + out["Isr1"] + out["Icoinf1"] + out["Iother1"])
  
  af.orapp.cross.1 = prev.nov.case.1 * (1 - (1/or.orapp.cross.1)) 
  af.orapp.nopostsym.cross.1 = prev.nov.case.1 * (1 - (1/or.orapp.nopostsym.cross.1))
  
  inc.orapp.cross.1 = af.orapp.cross.1*age.inc.1y
  inc.orapp.nopostsym.cross.1 = af.orapp.cross.1*age.inc.1y
  
  
  #12-23mo
  or.orapp.cross.2 = ((out["Is2"] + out["Isr2"] + out["Icoinf2"]) *                     
                        (out["S2"] + out["R2"])) /                          
    ((out["Ia2"] + out["Iar2"] + out["Ips2"] + out["Ipsr2"]) *                         
       out["Iother2"])   
  
  or.orapp.nopostsym.cross.2 = ((out["Is2"] + out["Isr2"] + out["Icoinf2"]) *                     
                                  (out["S2"] + out["R2"])) /                          
    ((out["Ia2"] + out["Iar2"]) *                         
       out["Iother2"])  
  
  prev.nov.case.2 =  (out["Is2"] + out["Isr2"] + out["Icoinf2"]) /                                                                
    (out["Is2"] + out["Isr2"] + out["Icoinf2"] + out["Iother2"])
  
  af.orapp.cross.2 = prev.nov.case.2 * (1 - (1/or.orapp.cross.2)) 
  af.orapp.nopostsym.cross.2 = prev.nov.case.2 * (1 - (1/or.orapp.nopostsym.cross.2))
  
  inc.orapp.cross.2 = af.orapp.cross.2*age.inc.2y
  inc.orapp.nopostsym.cross.2 = af.orapp.cross.2*age.inc.2y
  
  
  #Calculate coinfection to report 
  prop.coinf.cross.1 =  out["Icoinf1"]/(out["Icoinf1"] + out["Is1"] + out["Isr1"])
  prop.coinf.cross.2 =  out["Icoinf2"]/(out["Icoinf2"] + out["Is2"] + out["Isr2"])
  
  
  #Bring in file where I've already calculated incidence using the DE approach. Wasa originally using this as 'limits' for optimization procedure, but now just want to report the incidence under DE approach
  limits_boundedfit <- read.csv(here("Other/limits_boundedfit.csv")) 
  de.inc <- filter(limits_boundedfit, index.co==country) %>% select("norocaseprev.ai")
  de.inc.1 <- de.inc[1,]
  de.inc.2 <- de.inc[2,]
  
  pt <- c(adj.noro.inc.mod.1, adj.noro.inc.mod.2, 
          de.inc.1, de.inc.2,
          inc.orapp.cross.1, inc.orapp.cross.2,
          
          prop.coinf.cross.1, prop.coinf.cross.2,
          
          beta.noro.1, beta.noro.2, foi.other.1, foi.other.2, s.p, s.r, 
          age.inc.mod.1, age.inc.mod.2, age.inc.1y, age.inc.2y,
          obs.asymp.1, mod.asymp.1, obs.asymp.2, mod.asymp.2, obs.case.1, mod.case.1, obs.case.2, mod.case.2, 
          
          or.orapp.cross.1, or.orapp.nopostsym.cross.1,
          or.orapp.cross.2, or.orapp.nopostsym.cross.2,
          af.orapp.cross.1, af.orapp.nopostsym.cross.1,
          af.orapp.cross.2, af.orapp.nopostsym.cross.2,
          inc.orapp.nopostsym.cross.1, inc.orapp.nopostsym.cross.2, 
          
          obj.fx.val)
  
  braz.pt.all[[i]] <- pt
  
  #remove objects before running model for next country
  
  rm(adj.noro.inc.mod.1, adj.noro.inc.mod.2, 
     de.inc.1, de.inc.2,
     inc.orapp.cross.1, inc.orapp.cross.2,
     
     prop.coinf.cross.1, prop.coinf.cross.2,
     
     beta.noro.1, beta.noro.2, foi.other.1, foi.other.2, s.p, s.r, 
     age.inc.mod.1, age.inc.mod.2, 
     obs.asymp.1, mod.asymp.1, obs.asymp.2, mod.asymp.2, obs.case.1, mod.case.1, obs.case.2, mod.case.2, 
     
     or.orapp.cross.1, or.orapp.nopostsym.cross.1,
     or.orapp.cross.2, or.orapp.nopostsym.cross.2,
     af.orapp.cross.1, af.orapp.nopostsym.cross.1,
     af.orapp.cross.2, af.orapp.nopostsym.cross.2,
     inc.orapp.nopostsym.cross.1, inc.orapp.nopostsym.cross.2, 
     
     obj.fx.val)
  
}





#####################
# For PERU
#####################
# CHANGE THESE FOUR PARAMS FOR EACH COUNTRY
solvertype=21 #best
country = "Peru"
p0 = c(1.69, 1.83, 0.02, 0.02, 0.05, 0.05)
peru.pt.all <- list() # each list element is fit from a bootstrap sample

#solver types
if (solvertype==21)      {        alg.name="NLOPT_LN_SBPLX"      }
if (solvertype==25)      {        alg.name="NLOPT_LN_NELDERMEAD" }


lb=c(beta.noro.1.min, beta.noro.2.min, foi.other.1.min, foi.other.2.min, s.p.min, s.r.min)
ub=c(beta.noro.1.max, beta.noro.2.max, foi.other.1.max, foi.other.2.max, s.p.max, s.r.max)

obs.asymp.1.pt <- maled_forfit[country,"control.allno.1y"] 
obs.asymp.1.boot <- rnorm(n.boot,	mean	=	maled_forfit[country, "n.control.1"] ,	sd	=	maled_forfit[country,"sd.control.1"])

obs.asymp.2.pt <- maled_forfit[country,"control.allno.2y"] 
obs.asymp.2.boot <- rnorm(n.boot,	mean	=	maled_forfit[country, "n.control.2"] ,	sd	=	maled_forfit[country,"sd.control.2"])

obs.case.1.pt <- maled_forfit[country,"case.allno.1y"] 
obs.case.1.boot <- rnorm(n.boot,	mean	=	maled_forfit[country, "n.case.1"] ,	sd	=	maled_forfit[country,"sd.case.1"])

obs.case.2.pt <- maled_forfit[country,"case.allno.2y"] 
obs.case.2.boot <- rnorm(n.boot,	mean	=	maled_forfit[country, "n.case.2"] ,	sd	=	maled_forfit[country,"sd.case.2"])

age.inc.1y <- maled_forfit[country,"age.inc.1y.jpm"] 
age.inc.2y <- maled_forfit[country,"age.inc.2y.jpm"]


for(i in 1:n.boot){ 
  
  ##############################
  obs.asymp.1 <- obs.asymp.1.boot[i]/maled_forfit[country,"control.num.1y"]  # 0-11mo target
  obs.asymp.2 <- obs.asymp.2.boot[i]/maled_forfit[country,"control.num.2y"]  # 12-23mo target
  obs.case.1 <- obs.case.1.boot[i]/maled_forfit[country,"case.num.1y"]    # 0-11mo target
  obs.case.2 <- obs.case.2.boot[i]/maled_forfit[country,"case.num.2y"]    # 12-23mo target
  ##############################
  
  fres <- nloptr(x0=p0,eval_f=RM_solve,lb=lb,ub=ub,opts=list("algorithm"=alg.name,xtol_abs=1e-3,maxeval=stps,maxtime=mtime,print_level=1))
  
  fpar=(fres$solution); asymout.final=fres$objective; totiterate=fres$iterations; solverstatus=fres$status
  bestfit=c(fpar, asymout.final, totiterate, solverstatus)
  bestfit <- as.data.frame(bestfit)
  rownames(bestfit)=c("beta.noro.1", "beta.noro.2", "foi.other.1", "foi.other.2", "s.p", "s.r",
                      "objective function value", 
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
  
  adj.noro.inc.mod.1 = foi.noro.1*365*((out["S1"]*s.p + out["R1"]*s.r)/N.1)
  adj.noro.inc.mod.2 = foi.noro.2*365*((out["S2"]*s.p + out["R2"]*s.r)/N.2)
  
  # Now, for validaation, use the model to calculate the AI using the OR approach - for both age groups 
  
  #0-11mo
  or.orapp.cross.1 = ((out["Is1"] + out["Isr1"] + out["Icoinf1"]) *                     
                        (out["S1"] + out["R1"])) /                          
    ((out["Ia1"] + out["Iar1"] + out["Ips1"] + out["Ipsr1"]) *                         
       out["Iother1"])   
  
  or.orapp.nopostsym.cross.1 = ((out["Is1"] + out["Isr1"] + out["Icoinf1"]) *                     
                                  (out["S1"] + out["R1"])) /                          
    ((out["Ia1"] + out["Iar1"]) *                         
       out["Iother1"])  
  
  prev.nov.case.1 =  (out["Is1"] + out["Isr1"] + out["Icoinf1"]) /                                                                
    (out["Is1"] + out["Isr1"] + out["Icoinf1"] + out["Iother1"])
  
  af.orapp.cross.1 = prev.nov.case.1 * (1 - (1/or.orapp.cross.1)) 
  af.orapp.nopostsym.cross.1 = prev.nov.case.1 * (1 - (1/or.orapp.nopostsym.cross.1))
  
  inc.orapp.cross.1 = af.orapp.cross.1*age.inc.1y
  inc.orapp.nopostsym.cross.1 = af.orapp.cross.1*age.inc.1y
  
  
  #12-23mo
  or.orapp.cross.2 = ((out["Is2"] + out["Isr2"] + out["Icoinf2"]) *                     
                        (out["S2"] + out["R2"])) /                          
    ((out["Ia2"] + out["Iar2"] + out["Ips2"] + out["Ipsr2"]) *                         
       out["Iother2"])   
  
  or.orapp.nopostsym.cross.2 = ((out["Is2"] + out["Isr2"] + out["Icoinf2"]) *                     
                                  (out["S2"] + out["R2"])) /                          
    ((out["Ia2"] + out["Iar2"]) *                         
       out["Iother2"])  
  
  prev.nov.case.2 =  (out["Is2"] + out["Isr2"] + out["Icoinf2"]) /                                                                
    (out["Is2"] + out["Isr2"] + out["Icoinf2"] + out["Iother2"])
  
  af.orapp.cross.2 = prev.nov.case.2 * (1 - (1/or.orapp.cross.2)) 
  af.orapp.nopostsym.cross.2 = prev.nov.case.2 * (1 - (1/or.orapp.nopostsym.cross.2))
  
  inc.orapp.cross.2 = af.orapp.cross.2*age.inc.2y
  inc.orapp.nopostsym.cross.2 = af.orapp.cross.2*age.inc.2y
  
  
  #Calculate coinfection to report 
  prop.coinf.cross.1 =  out["Icoinf1"]/(out["Icoinf1"] + out["Is1"] + out["Isr1"])
  prop.coinf.cross.2 =  out["Icoinf2"]/(out["Icoinf2"] + out["Is2"] + out["Isr2"])
  
  
  #Bring in file where I've already calculated incidence using the DE approach. Wasa originally using this as 'limits' for optimization procedure, but now just want to report the incidence under DE approach
  limits_boundedfit <- read.csv(here("Other/limits_boundedfit.csv")) 
  de.inc <- filter(limits_boundedfit, index.co==country) %>% select("norocaseprev.ai")
  de.inc.1 <- de.inc[1,]
  de.inc.2 <- de.inc[2,]
  
  pt <- c(adj.noro.inc.mod.1, adj.noro.inc.mod.2, 
          de.inc.1, de.inc.2,
          inc.orapp.cross.1, inc.orapp.cross.2,
          
          prop.coinf.cross.1, prop.coinf.cross.2,
          
          beta.noro.1, beta.noro.2, foi.other.1, foi.other.2, s.p, s.r, 
          age.inc.mod.1, age.inc.mod.2, age.inc.1y, age.inc.2y,
          obs.asymp.1, mod.asymp.1, obs.asymp.2, mod.asymp.2, obs.case.1, mod.case.1, obs.case.2, mod.case.2, 
          
          or.orapp.cross.1, or.orapp.nopostsym.cross.1,
          or.orapp.cross.2, or.orapp.nopostsym.cross.2,
          af.orapp.cross.1, af.orapp.nopostsym.cross.1,
          af.orapp.cross.2, af.orapp.nopostsym.cross.2,
          inc.orapp.nopostsym.cross.1, inc.orapp.nopostsym.cross.2, 
          
          obj.fx.val)
  
  peru.pt.all[[i]] <- pt
  
  #remove objects before running model for next country
  
  rm(adj.noro.inc.mod.1, adj.noro.inc.mod.2, 
     de.inc.1, de.inc.2,
     inc.orapp.cross.1, inc.orapp.cross.2,
     
     prop.coinf.cross.1, prop.coinf.cross.2,
     
     beta.noro.1, beta.noro.2, foi.other.1, foi.other.2, s.p, s.r, 
     age.inc.mod.1, age.inc.mod.2, 
     obs.asymp.1, mod.asymp.1, obs.asymp.2, mod.asymp.2, obs.case.1, mod.case.1, obs.case.2, mod.case.2, 
     
     or.orapp.cross.1, or.orapp.nopostsym.cross.1,
     or.orapp.cross.2, or.orapp.nopostsym.cross.2,
     af.orapp.cross.1, af.orapp.nopostsym.cross.1,
     af.orapp.cross.2, af.orapp.nopostsym.cross.2,
     inc.orapp.nopostsym.cross.1, inc.orapp.nopostsym.cross.2, 
     
     obj.fx.val)
  
}




###############################################################

# OUTPUT MODEL RESULTS FOR POINT ESTIMATES

###############################################################


# RESULTS - MODEL-BASED INCIDENCE BY COUNTRY

# Look at model fits (obj function values), coinfection prevalence to report in manuscript to support model fit/validity

res.mat.bang <- matrix(unlist(bang.pt.all), ncol=length(bang.pt.all[[1]]), nrow=n.boot, byrow = TRUE)
res.mat.ind <- matrix(unlist(ind.pt.all), ncol=length(ind.pt.all[[1]]), nrow=n.boot, byrow = TRUE)
res.mat.nep <- matrix(unlist(nep.pt.all), ncol=length(nep.pt.all[[1]]), nrow=n.boot, byrow = TRUE)
res.mat.pak <- matrix(unlist(pak.pt.all), ncol=length(pak.pt.all[[1]]), nrow=n.boot, byrow = TRUE)
res.mat.braz <- matrix(unlist(braz.pt.all), ncol=length(braz.pt.all[[1]]), nrow=n.boot, byrow = TRUE)
res.mat.peru <- matrix(unlist(peru.pt.all), ncol=length(peru.pt.all[[1]]), nrow=n.boot, byrow = TRUE)

res <- data.frame(rbind(res.mat.bang, res.mat.ind, res.mat.nep, res.mat.pak, res.mat.braz, res.mat.peru))

res <- data.frame(cbind(res, c(rep("Bangladesh", n.boot), rep("India", n.boot), rep("Nepal", n.boot), rep("Pakistan", n.boot), rep("Brazil", n.boot), rep("Peru", n.boot))))

colnames(res) <- c("adj.noro.inc.mod.1", "adj.noro.inc.mod.2", 
                    "de.inc.1", "de.inc.2",
                    "inc.orapp.cross.1", "inc.orapp.cross.2",
                    
                    "prop.coinf.cross.1", "prop.coinf.cross.2",
                    
                    "beta.noro.1", "beta.noro.2", "foi.other.1", "foi.other.2", "s.p", "s.r", 
                    "age.inc.mod.1", "age.inc.mod.2", "age.inc.1y", "age.inc.2y",
                    "obs.asymp.1", "mod.asymp.1", "obs.asymp.2", "mod.asymp.2", "obs.case.1", "mod.case.1", "obs.case.2", "mod.case.2", 
                    
                    "or.orapp.cross.1", "or.orapp.nopostsym.cross.1",
                    "or.orapp.cross.2", "or.orapp.nopostsym.cross.2",
                    "af.orapp.cross.1", "af.orapp.nopostsym.cross.1",
                    "af.orapp.cross.2", "af.orapp.nopostsym.cross.2",
                    "inc.orapp.nopostsym.cross.1", "inc.orapp.nopostsym.cross.2", 
                    
                   "obj.fx.val",
                   "country")

outTable <- paste("newstrucIotherspsronly_withormethod_",n.boot,"bootcis_solv_",date,"full.csv", sep="")

write.csv(res, outTable)

res %>%
  group_by(country) %>%
  summarise(med = median(obj.fx.val))

res %>%
  group_by(country) %>%
  summarise(med.obj.fx.val = median(obj.fx.val),
            med.coinf.1 = median(prop.coinf.cross.1), med.coinf.2 = median(prop.coinf.cross.2),
            med.de.inc.1 = median(de.inc.1), med.de.inc.2 = median(de.inc.2), 
            med.or.inc.1 = median(inc.orapp.cross.1), med.or.inc.2 = median(inc.orapp.cross.2), 
            med.mb.inc.1 = median(adj.noro.inc.mod.1), med.mb.inc.2 = median(adj.noro.inc.mod.2),
            med.sp = median(s.p), med.sr = median(s.r),
            med.b1 = median(beta.noro.1), med.b2 = median(beta.noro.2),
            med.af.1 = median(adj.noro.inc.mod.1/age.inc.1y), med.af.2 = median(adj.noro.inc.mod.2/age.inc.2y)) -> res.summ

outTable2 <- paste("newstrucIotherspsronly_withormethod_",n.boot,"bootcis_",date,"summ.csv", sep="")

write.csv(res.summ, outTable2)


#################################################

# OUTPUT CONFIDENCE INTERVALS

#################################################


samplemedian <- function(x, d) {
  return(median(x[d]))
}

# Bangladesh
inc.1.dsn <- sapply(bang.pt.all, "[", 1)
inc.2.dsn <- sapply(bang.pt.all, "[", 2)

inc.1.dsn.med = boot(inc.1.dsn, samplemedian, R=1000) 
inc.2.dsn.med = boot(inc.2.dsn, samplemedian, R=1000) 

inc.1.dsn.ci = boot.ci(inc.1.dsn.med, type="basic")
inc.1.bang <- c(inc.1.dsn.med$t0, inc.1.dsn.ci$basic[1,4:5])

inc.2.dsn.ci = boot.ci(inc.2.dsn.med, type="basic")
inc.2.bang <- c(inc.2.dsn.med$t0, inc.2.dsn.ci$basic[1,4:5])

# India
inc.1.dsn <- sapply(ind.pt.all, "[", 1)
inc.2.dsn <- sapply(ind.pt.all, "[", 2)

inc.1.dsn.med = boot(inc.1.dsn, samplemedian, R=1000) 
inc.2.dsn.med = boot(inc.2.dsn, samplemedian, R=1000) 

inc.1.dsn.ci = boot.ci(inc.1.dsn.med, type="basic")
inc.1.ind <- c(inc.1.dsn.med$t0, inc.1.dsn.ci$basic[1,4:5])

inc.2.dsn.ci = boot.ci(inc.2.dsn.med, type="basic")
inc.2.ind <- c(inc.2.dsn.med$t0, inc.2.dsn.ci$basic[1,4:5])


# Nepal
inc.1.dsn <- sapply(nep.pt.all, "[", 1)
inc.2.dsn <- sapply(nep.pt.all, "[", 2)

inc.1.dsn.med = boot(inc.1.dsn, samplemedian, R=1000) 
inc.2.dsn.med = boot(inc.2.dsn, samplemedian, R=1000) 

inc.1.dsn.ci = boot.ci(inc.1.dsn.med, type="basic")
inc.1.nep <- c(inc.1.dsn.med$t0, inc.1.dsn.ci$basic[1,4:5])

inc.2.dsn.ci = boot.ci(inc.2.dsn.med, type="basic")
inc.2.nep <- c(inc.2.dsn.med$t0, inc.2.dsn.ci$basic[1,4:5])

# Pakistan
inc.1.dsn <- sapply(pak.pt.all, "[", 1)
inc.2.dsn <- sapply(pak.pt.all, "[", 2)

inc.1.dsn.med = boot(inc.1.dsn, samplemedian, R=1000) 
inc.2.dsn.med = boot(inc.2.dsn, samplemedian, R=1000) 

inc.1.dsn.ci = boot.ci(inc.1.dsn.med, type="basic")
inc.1.pak <- c(inc.1.dsn.med$t0, inc.1.dsn.ci$basic[1,4:5])

inc.2.dsn.ci = boot.ci(inc.2.dsn.med, type="basic")
inc.2.pak <- c(inc.2.dsn.med$t0, inc.2.dsn.ci$basic[1,4:5])

# Brazil
inc.1.dsn <- sapply(braz.pt.all, "[", 1)
inc.2.dsn <- sapply(braz.pt.all, "[", 2)

inc.1.dsn.med = boot(inc.1.dsn, samplemedian, R=1000) 
inc.2.dsn.med = boot(inc.2.dsn, samplemedian, R=1000) 

inc.1.dsn.ci = boot.ci(inc.1.dsn.med, type="basic")
inc.1.braz <- c(inc.1.dsn.med$t0, inc.1.dsn.ci$basic[1,4:5])

inc.2.dsn.ci = boot.ci(inc.2.dsn.med, type="basic")
inc.2.braz <- c(inc.2.dsn.med$t0, inc.2.dsn.ci$basic[1,4:5])


# Peru
inc.1.dsn <- sapply(peru.pt.all, "[", 1)
inc.2.dsn <- sapply(peru.pt.all, "[", 2)

inc.1.dsn.med = boot(inc.1.dsn, samplemedian, R=1000) 
inc.2.dsn.med = boot(inc.2.dsn, samplemedian, R=1000) 

inc.1.dsn.ci = boot.ci(inc.1.dsn.med, type="basic")
inc.1.per <- c(inc.1.dsn.med$t0, inc.1.dsn.ci$basic[1,4:5])

inc.2.dsn.ci = boot.ci(inc.2.dsn.med, type="basic")
inc.2.per <- c(inc.2.dsn.med$t0, inc.2.dsn.ci$basic[1,4:5])


results_byco.boot <- rbind(inc.1.bang, inc.2.bang,
                           inc.1.ind, inc.2.ind,
                           inc.1.nep, inc.2.nep,
                           inc.1.pak, inc.2.pak,
                           inc.1.braz, inc.2.braz, 
                           inc.1.per, inc.2.per)

outTable <- paste("newstrucIotherspsronly_withormethod_",n.boot,"bootcis_",date,".csv", sep="")

write.csv(results_byco.boot, outTable)

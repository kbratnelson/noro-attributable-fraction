rm(list=ls()) 

##############################
#
# Load packages
#
##############################

library(rootSolve) #ODE solver
library(here) 


##############################
#
# Input model parameters
#
##############################

# Demography 
#birth.rate = (25/1000)/365   # u5 mortality rate in the UK (changing this does not appreciably affect results)


# Natural history
sym = 1/3                    # symptomatic period, 1/days (Lopman et al 2014)
#asym = 1/15                  # asymptomatic period, 1/days (Lopman et al 2014)
#wan = 1/(365*5.1)            # immunity rate, 1/days (fitted in Lopman et al. 2014)
post = 0.05                  # relative infectiousness during asymptomatic infection (Lopman et al 2014)
#ps = 1/15                    # post-symptomatic shedding period, 1/days
other = 1/5                  # duration of non-norovirus diarrheal disease episode
coinf.postshed = 1/10        # duration of noro shedding after symptomatic infection with another pathogen


# parameters to fit/fix/vary in models
#s1                    # symptomatic fraction, primary infection - FIX
#s2                    # symptomatic fraction, secondary infection - FIX

#beta.noro           # effective contact rate - ESTIMATE
#foi.other           # effective contact rate - ESTIMATE


# CONSTRAINT 1: s1 > s2

# Initial conditions
xstart = c(S = 79000, Is = 1000, Ips = 0, Ia = 0, 
           R = 0, Isr = 0, Ipsr = 0, Iar = 0,
           Iother = 10000, Icoinf = 10000)


###########################################
#
# Function for transmission model 
#
###########################################

RM = function(t, x, parms) {
  
  S  = x[1]   # susceptible
  Is = x[2]   # infected-symptomatic
  Ips = x[3]  # infected- postsymptomatic
  Ia = x[4]   # infected-asymptomatic
  
  R  = x[5]   # recovered 
  Isr = x[6]  # infected-symptomatic reinfected
  Ipsr = x[7] # infected- postsymptomatic reinfected
  Iar = x[8]  # infected-asymptomatic reinfected
  
  Iother = x[9] # diarrhea due to another pathogen
  Icoinf = x[10] # diarrhea due to another pathogen, with norovirus coinfection
  
  N = S + Is + Ips + Ia + R + Isr + Ipsr + Iar + Iother + Icoinf
  
  # Calculate total deaths in a time step  
  #deaths = birth.rate*N
  
  # Create population equillibrium by setting the number of births equal to the number of deaths
  #births = deaths
  
  #Calculate the number of individuals in each age group that can transmit the virus at each time step
  infec.noro = Is + Isr + post*(Ips + Ia + Ipsr + Iar + Icoinf)
  
  #infec.other = Iother + Icoinf
  
  #FoI 
  FOI.noro = beta.noro * (infec.noro/N)
  
  #FoI 
  #FOI.other = beta.other * (infec.other/N) 
  
  # Model equations
  dS = wan*R - FOI.noro*S - FOI.other*S + other*Iother
  dIs = s1*FOI.noro*S - sym*Is - FOI.other*Is
  dIps = sym*Is - asym*Ips - FOI.other*Ips 
  dIa = (1-s1)*FOI.noro*S - asym*Ia - FOI.other*Ia
  
  dR = asym*Ia + asym*Ips + asym*Iar + asym*Ipsr - FOI.other*R - FOI.noro*R - wan*R + other*Icoinf 
  dIsr = s2*FOI.noro*R - sym*Isr - FOI.other*Isr
  dIpsr = sym*Isr - asym*Ipsr - FOI.other*Ipsr
  dIar = (1-s2)*FOI.noro*R - asym*Iar - FOI.other*Iar
  
  dIother = FOI.other*S + FOI.other*R - other*Iother - FOI.noro*Iother
  dIcoinf = FOI.noro*Iother - other*Icoinf + FOI.other*(Is + Ia + Ips + Isr + Iar + Ipsr)
  
  
  res = c(dS, dIs, dIps, dIa, dR, dIsr, dIpsr, dIar, dIother, dIcoinf)
  list(res)   
  
  #})
}  


# TO RUN A SWEEP


# varied parameters

s1val = c(0.3, 0.5, 0.7, 0.9)                  # symptomatic fraction, primary infection 
s2val = c(0.1, 0.3, 0.5, 0.7)                  # symptomatic fraction, secondary infection 

beta.noroval = c(0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95)                # effective contact rate
foi.otherval = c(0.001, 0.003, 0.005, 0.007, 0.01, 0.0125, 0.015, 0.0175, 0.02, 0.0225, 0.025)      # effective contact rate 

asymval = c(1/5, 1/10, 1/15, 1/20, 1/30, 1/60)
wanval = c(1/(5.1*365), 1/(3*365), 1/365, 1/182.5, 1/60)


sweep <- expand.grid(s1=s1val, 
                     s2=s2val, 
                     beta.noro = beta.noroval, 
                     FOI.other = foi.otherval, 
                     asym = asymval,
                     wan = wanval)

subset(sweep, 
       beta.noro == 0.45 & FOI.other == 0.001 |
         beta.noro == 0.5 & FOI.other ==  0.003|
         beta.noro == 0.55 & FOI.other == 0.005 |
         beta.noro == 0.60 & FOI.other == 0.007 |
         beta.noro == 0.65 & FOI.other == 0.01 |
         beta.noro == 0.70 & FOI.other == 0.0125 |
         beta.noro == 0.75 & FOI.other == 0.015 |
         beta.noro == 0.80 & FOI.other == 0.0175 |
         beta.noro == 0.85 & FOI.other == 0.02 |
         beta.noro == 0.90 & FOI.other == 0.0225 |
         beta.noro == 0.95 & FOI.other == 0.025) -> sweep

data.summary <- list(dim(sweep)[1])

for(i in 1:dim(sweep)[1]){
  data.summary[[i]] <- data.frame(nrow=1, ncol=29)
}

for(i in 1:dim(sweep)[1]){
  s1 <- sweep$s1[i]
  s2 <- sweep$s2[i]
  beta.noro <- sweep$beta.noro[i]
  FOI.other <- sweep$FOI.other[i]
  asym <- sweep$asym[i]
  wan <- sweep$wan[i]
  
  out = runsteady(y = xstart, fun = RM, times = c(0, 1e5))
  
  out = as.data.frame(out)
  
  N = sum(out["S",], out["Is",], out["Ips",], out["Ia",], out["R",], out["Isr",], out["Ipsr",], out["Iar",], out["Iother",], out["Icoinf",])
  
  # Number infectious: symptomatic cases plus asymptomatic cases (asymptomatic cases less infectious by factor of 'post')
  infec.noro = out["Is",] + out["Isr",] + 
    post*(out["Ips",] + out["Ipsr",] + out["Ia",] + out["Iar",] + out["Icoinf",])
  
  infec.other = out["Icoinf",] + out["Iother",]
  
  #FoI
  FOI.noro = beta.noro * (infec.noro/N) 
  
  #FOI.other = beta.other * (infec.other/N)
  
  # Effective reproduction number for noro
  r0 = beta.noro * (1/(sym + post*asym)) 
  
  # Modeled asymptomatic prevalence (number shedding norovirus/number without illness)
  asymp.out = sum(out["Ips",], out["Ia",], out["Ipsr",], out["Iar",]) /
    sum(out["Ips",], out["Ia",], out["Ipsr",], out["Iar",], out["S",], out["R",]) 
  
  
  # Calculate overall AGE incidence (per person-year) 
  inc.age = FOI.noro*365*(out["S",]*s1 + out["R",]*s2)/N + 
    FOI.other*365*out["S",]/N + 
    FOI.noro*365*out["Iother",]/N +
    FOI.other*365*(out["Ia",] + out["Ips",] +  out["Ipsr",] + out["Iar",])/N
  
  
  #----------------------------------------------------------
  # Calculate attributable fractions and incidence under each approach
  #----------------------------------------------------------
  
  # Detection as etiology (cross-sec - prevalent cases): Approach 1
  ## use steady state proportions in each compartment 'slice-in-time'
  af.detect.cross = (out["Is",] + out["Isr",] + out["Icoinf",]) /                                
    (out["Is",] + out["Isr",] + out["Icoinf",] + out["Iother",]) 
  
  ai.detect.cross = af.detect.cross * inc.age
  
  af.detect.unbiased.cross = (out["Is",] + out["Isr",])/                                
    (out["Is",] + out["Isr",] + out["Icoinf",] + out["Iother",]) 
  
  ai.detect.unbiased.cross = af.detect.unbiased.cross*inc.age
  
  prop.coinf.cross =  out["Icoinf",]/(out["Icoinf",] + out["Is",] + out["Isr",])
  
  # Detection as etiology (cohort - prevalent cases): Approach 1b
  ## use rate of entry into compartments as Approach 1
  #inc.detect.cohort =  FOI.noro*365*((out["S",]*s1 + out["R",]*s2 + out["Iother",])/N) 
  
  #af.detect.cohort = inc.detect.cohort / inc.age
  
  #prop.coinf.cohort =  ((FOI.noro*365*out["Iother",])/N) / (FOI.noro*365*((out["S",]*s1 + out["R",]*s2 + out["Iother",])/N)) 
  
  
  ######################
  
  
  # OR approach (cross-sec - prevalent cases): Approach 2
  ## Note: include primary asymptomatic infection (Ia), reinfected (Iar, Ipsr), and post-symptomatic shedders (Ips) among controls AND coinfected among cases (Icoinf) - FULLY BIASED
  or.orapp.cross = ((out["Is",] + out["Isr",] + out["Icoinf",]) *           # P(Nov+, cases) * 
                      (out["S",] + out["R",])) /                            # P(Nov-, ctrls) /
    ((out["Ia",] + out["Ips",] + out["Ipsr",] + out["Iar",]) *              # P(Nov+ ctrls) *
       out["Iother",])                                                      # P(Nov- cases)  
  
  # exclude reinfection (Iar) from controls 
  or.orapp.noreinf.cross = ((out["Is",] + out["Isr",] + out["Icoinf",]) *   # P(Nov+, cases) * 
                              (out["S",] + out["R",])) /                    # P(Nov-, ctrls) /
    ((out["Ia",] + out["Ips",] + out["Ipsr",]) *                            # P(Nov+ ctrls) *
       out["Iother",])                                                      # P(Nov- cases)  
  
  # exclude post-symptomatic shedders (Ips, Ipsr) from controls 
  or.orapp.nopostsym.cross = ((out["Is",] + out["Isr",] + out["Icoinf",]) *        # P(Nov+, cases) * 
                                     (out["S",] + out["R",])) /                         # P(Nov-, ctrls) /
    ((out["Ia",] + out["Iar",]) *                                                       # P(Nov+ ctrls) *
       out["Iother",])                                                                  # P(Nov- cases)
  
  # exclude asymptomatic primary infection (Ia) from controls
  or.orapp.noprimasym.cross = ((out["Is",] + out["Isr",] + out["Icoinf",]) *        # P(Nov+, cases) * 
                                             (out["S",] + out["R",])) /                         # P(Nov-, ctrls) /
    ((out["Ips",] + out["Ipsr",] + out["Iar",]) *                                               # P(Nov+ ctrls) *
       out["Iother",])                                                                          # P(Nov- cases)
  
  
  # Exclude coinfections from cases
  or.orapp.nocoinf.cross = ((out["Is",] + out["Isr",]) *        # P(Nov+, cases) * 
                                             (out["S",] + out["R",])) /        # P(Nov-, ctrls) /
    ((out["Ia",] + out["Ips",] + out["Ipsr",] + out["Iar",]) *                 # P(Nov+ ctrls) *
       out["Iother",])                                                         # P(Nov- cases)
  
  
  # Unbiased
  or.orapp.unbiased.cross = ((out["Is",] + out["Isr",]) *        # P(Nov+, cases) * 
                                             (out["S",] + out["R",])) /        # P(Nov-, ctrls) /
    (0.00001 *                                                                   # P(Nov+ ctrls) *
       out["Iother",])                                                         # P(Nov- cases)
  
  
  prev.nov.case =  (out["Is",] + out["Isr",] + out["Icoinf",]) /                                                                 # symptomatic norovirus infection /
    (out["Is",] + out["Isr",] + out["Icoinf",] + out["Iother",])
  
  af.orapp.cross = prev.nov.case * (1 - (1/or.orapp.cross)) # AF equation from GEMS methods paper (which MAL-ED also used) Blackwelder et al. CID 2012
  af.orapp.noreinf.cross = prev.nov.case * (1 - (1/or.orapp.noreinf.cross))
  af.orapp.nopostsym.cross = prev.nov.case * (1 - (1/or.orapp.nopostsym.cross))
  af.orapp.noprimasym.cross = prev.nov.case * (1 - (1/or.orapp.noprimasym.cross))
  af.orapp.nocoinf.cross = prev.nov.case * (1 - (1/or.orapp.nocoinf.cross))
  af.orapp.unbiased.cross = prev.nov.case * (1 - (1/or.orapp.unbiased.cross))
  
  ai.orapp.cross = af.orapp.cross*inc.age
  ai.orapp.noreinf.cross = af.orapp.noreinf.cross*inc.age
  ai.orapp.nopostsym.cross = af.orapp.nopostsym.cross*inc.age
  ai.orapp.noprimasym.cross = af.orapp.noprimasym.cross*inc.age
  ai.orapp.nocoinf.cross = af.orapp.nocoinf.cross*inc.age
  ai.orapp.unbiased.cross = af.orapp.unbiased.cross*inc.age
  
  
  # OR approach (cohort - incident cases): Approach 2b
  
  #or.orapp.cohort = ((FOI.noro*365*out["S",]*s1/N + FOI.noro*365*out["R",]*s2/N + FOI.noro*365*out["Iother",]/N) *     # inc. nov + cases
  #                     (out["S",]/N + out["R",]/N)) /                                                                     # prev. nov - controls
  #  ((out["Ia",]/N + out["Iar",]/N + out["Ips",]/N + out["Ipsr",]/N) *                                 # prev. nov + controls
  #     (FOI.other*365*out["S",]/N + FOI.other*365*out["R",]/N))                                           # inc. nov - cases
  
  #or.orapp.noreinf.cohort = (FOI.noro*365*out["S",]*s1/N + FOI.noro*365*out["R",]*s2/N + out["Iother",]/N) *           # inc. nov + cases
  #  (out["S",] + out["R",]) /                                                                  # prev. nov - controls
  #  (FOI.noro*out["S",]*(1-s1) + sym*out["Is",] + sym*out["Isr",]) *                           # prev. nov + controls
  #  (FOI.other*365*out["S",]/N + FOI.other*365*out["R",]/N)                                    # inc. nov - cases
  
  #or.orapp.noreinfpostsym.cohort = (FOI.noro*365*out["S",]*s1/N + FOI.noro*365*out["R",]*s2/N + out["Iother",]/N) *    # inc. nov + cases
  #  (out["S",] + out["R",]) /                                                          # prev. nov - controls
  #  (out["S",]*(1-s1)) *                                                               # prev. nov + controls
  #  (FOI.other*365*out["S",]/N + FOI.other*365*out["R",]/N)                            # inc. nov - cases
  
  #or.orapp.noreinfpostsymprimasym.cohort = (FOI.noro*365*out["S",]*s1/N + FOI.noro*365*out["R",]*s2/N + out["Iother",]/N) *    # inc. nov + cases
  #  (out["S",] + out["R",]) /                                                          # prev. nov - controls
  #  (0.001) *                                                               # prev. nov + controls
  #  (FOI.other*365*out["S",]/N + FOI.other*365*out["R",]/N)                            # inc. nov - cases
  
  #af.orapp.cohort = prev.nov.case * (1 - (1/or.orapp.cohort)) 
  #af.orapp.noreinf.cohort = prev.nov.case * (1 - (1/or.orapp.noreinf.cohort)) 
  #af.orapp.noreinfpostsym.cohort = prev.nov.case * (1 - (1/or.orapp.noreinfpostsym.cohort)) 
  #af.orapp.noreinfpostsymprimasym.cohort = prev.nov.case * (1 - (1/or.orapp.noreinfpostsymprimasym.cohort)) 
  
  #inc.orapp.cohort = af.orapp.cohort * inc.age
  #inc.orapp.noreinf.cohort = af.orapp.noreinf.cohort * inc.age
  #inc.orapp.noreinfpostsym.cohort = af.orapp.noreinfpostsym.cohort * inc.age
  #inc.orapp.noreinfpostsymprimasym.cohort = af.orapp.noreinfpostsymprimasym.cohort * inc.age
  #inc.orapp.noreinfpostsymprimasym.cohort
  
  
  ###########################
  
  # Model-based approach: Approach 3
  
  inc.mod = FOI.noro*365*((out["S",]*s1 + out["R",]*s2)/N)
  
  # Outputs
  table <- numeric(30)
  
  table[1] <- out["S",]/N
  table[2] <- out["Is",]/N
  table[3] <- out["Ips",]/N
  table[4] <- out["Ia",]/N
  table[5] <- out["R",]/N
  table[6] <- out["Isr",]/N
  table[7] <- out["Ipsr",]/N 
  table[8] <- out["Iar",]/N 
  table[9] <- out["Iother",]/N 
  table[10] <- out["Icoinf",]/N 
  
  table[11] <- af.detect.cross                        
  table[12] <- ai.detect.cross 
  table[13] <- af.detect.unbiased.cross
  table[14] <- ai.detect.unbiased.cross
  
  table[15] <- prop.coinf.cross
  
  table[16] <- prev.nov.case
  
  table[17] <- or.orapp.cross
  table[18] <- or.orapp.noreinf.cross
  table[19] <- or.orapp.nopostsym.cross
  table[20] <- or.orapp.noprimasym.cross
  table[21] <- or.orapp.nocoinf.cross
  table[22] <- or.orapp.unbiased.cross
  
  table[23] <- af.orapp.cross
  table[24] <- af.orapp.noreinf.cross
  table[25] <- af.orapp.nopostsym.cross
  table[26] <- af.orapp.noprimasym.cross
  table[27] <- af.orapp.nocoinf.cross
  table[28] <- af.orapp.unbiased.cross

  table[29] <- ai.orapp.cross
  table[30] <- ai.orapp.noreinf.cross
  table[31] <- ai.orapp.nopostsym.cross
  table[32] <- ai.orapp.noprimasym.cross
  table[33] <- ai.orapp.nocoinf.cross
  table[34] <- ai.orapp.unbiased.cross
  
  table[35] <- inc.mod
  table[36] <- inc.age
  table[37] <- inc.mod/inc.age
  
  table[38] <- r0
  
  table[39] <- asymp.out # prevalence of infection among 'controls' (AGE-free)
  table[40] <- (out["Iar",]) / (out["Ia",] + out["Ips",] + out["Ipsr",] + out["Iar",]) # proportion of all asymptomatic infection that is reinfection
  table[41] <- (out["Ips",] + out["Ipsr",]) / (out["Ia",] + out["Ips",] + out["Ipsr",] + out["Iar",]) # proportion of all asymptomatic infection that is postsymptomatic shedding
  
  table[42] <- N
  
  
  # Name it
  names(table) <- c("S", "Is", "Ips", "Ia", "R", "Isr", "Ispr", "Iar", "Iother", "Icoinf",
                    
                    "af.detect.cross",                   
                    "ai.detect.cross",
                    "af.detect.unbiased.cross",
                    "ai.detect.unbiased.cross",
                    
                    "prop.coinf.cross",
                    
                    "prev.nov.case",
                    
                    "or.orapp.cross",
                    "or.orapp.noreinf.cross",
                    "or.orapp.nopostsym.cross",
                    "or.orapp.noprimasym.cross",
                    "or.orapp.nocoinf.cross",
                    "or.orapp.unbiased.cross",
                    
                    "af.orapp.cross",
                    "af.orapp.noreinf.cross",
                    "af.orapp.nopostsym.cross",
                    "af.orapp.noprimasym.cross",
                    "af.orapp.nocoinf.cross",
                    "af.orapp.unbiased.cross",
                    
                    "ai.orapp.cross",
                    "ai.orapp.noreinf.cross",
                    "ai.orapp.nopostsym.cross",
                    "ai.orapp.noprimasym.cross",
                    "ai.orapp.nocoinf.cross",
                    "ai.orapp.unbiased.cross",
                    
                    "ai.mb",
                    "ai.age",
                    "af.mb",
                    
                    "r0",
                    
                    "asymp.out", "prop.asym.reinf", "prop.asym.postsym",
                    
                    "N")
  
  data.summary[[i]] <- table
  
}

sweep <- as.data.frame(sweep)

sweep$sweepnum<-c(1:nrow(sweep))
for(i in 1:nrow(sweep)){
  data.summary[[i]]$sweepnum<-i
}

require(dplyr)

UnlistDat = bind_rows(data.summary, .id = "sweepnum")

Merge = merge(sweep, UnlistDat[,2:ncol(UnlistDat)], by='sweepnum', all.y=TRUE)
Merge

save.image('/Users/kristinnelson/Box Sync/NoVaMod//Case-control attributable fractions/Results/simstudy/simsweep021920.RData')

















####################################################################################
# DATA EXPLORATION
####################################################################################

#load("/Users/kristinnelson/Box Sync/NoVaMod//Case-control attributable fractions/Results/simstudy/simsweep121520.RData")

# Some data exploration
library(ggplot2)

ggplot(Merge, aes(x=sweepnum, y=r0)) + geom_point()
ggplot(Merge, aes(x=beta.noro, y=r0)) + geom_point() # r0 totally dependent on beta.noro


ggplot(Merge, aes(x=sweepnum, y=asymp.out)) + geom_point()
ggplot(Merge, aes(x=beta.noro, y=asymp.out)) + geom_point() # all values of beta.noro can give asymp.out values of 0-0.4

ggplot(Merge, aes(x=sweepnum, y=inc.mod)) + geom_point()  
ggplot(Merge, aes(x=beta.noro, y=inc.mod)) + geom_point() # reasonable values of annual noro incidence (0-2 per CY) consistent with range of beta.noro 0-5

ggplot(Merge, aes(x=sweepnum, y=inc.age)) + geom_point() 

#What is the intersection of all reasonable values? - can change this to a wider or narrower range, but this is a good start
#Assume af.mod is 1-30%
subset(Merge,
       inc.age > 0.1 & inc.age < 10 &
         inc.mod > 0.001 & inc.mod < 10 &
         asymp.out > 0.001 & asymp.out < 0.5 & 
         af.mod > 0.01 & af.mod < 0.50) -> Merge_subset 
Merge_subset



# check ranges of important vars
ggplot(Merge_subset, aes(x=sweepnum, y=inc.age)) + geom_point() 
ggplot(Merge_subset, aes(x=sweepnum, y=asymp.out)) + geom_point()
ggplot(Merge_subset, aes(x=sweepnum, y=inc.mod)) + geom_point() 
ggplot(Merge_subset, aes(x=sweepnum, y=r0)) + geom_point()  
ggplot(Merge_subset, aes(x=sweepnum, y=af.mod)) + geom_point()  

# chaeck relationships between variables
ggplot(Merge_subset, aes(x=inc.age, y=inc.mod, color=FOI.other)) + geom_point()  -> ageinc.noroinc
ggplot(Merge_subset, aes(x=inc.age, y=inc.mod, color=r0)) + geom_point()  -> ageinc.noroinc.r0
ggplot(Merge_subset, aes(x=inc.age, y=inc.mod, color=beta.noro)) + geom_point()  

ggplot(Merge_subset, aes(x=inc.age, y=r0)) + geom_point()  # want r0 to increase with inc.age?
ggplot(Merge_subset, aes(x=inc.age, y=asymp.out, color=asym)) + geom_point()  -> ageinc.asympout

ggplot(Merge_subset, aes(x=inc.mod, y=asymp.out, color=interaction(s1,s2))) + geom_point()  
ggplot(Merge_subset, aes(x=inc.mod, y=asymp.out, color=asym)) + geom_point() 

ggplot(Merge_subset, aes(x=inc.age, y=prop.coinf.cohort, color = af.mod)) + geom_point() 
ggplot(Merge_subset, aes(x=inc.age, y=prop.coinf.cohort, color = asymp.out)) + geom_point() 

ggplot(Merge_subset, aes(x=inc.age, y=inc.mod, color=asymp.out)) + geom_point()  # check with different data cuts -- pairing beta.noro values with FOI.other values
# want af.mod to stay in a reasonable range and inc.mod to increase with inc.age (can level out or decline, that's fine)


  # create colors for following plots
  library(RColorBrewer)
brewer.pal(n = 3, name = "Dark2")

# Main - restrict by asym, wan
ggplot(data=subset(Merge_subset, Merge_subset$s1 == 0.7 & 
                     Merge_subset$s2 == 0.3 &
                     Merge_subset$beta.noro == 0.4 &
                     Merge_subset$wan == 1/(365*5.1) &
                     Merge_subset$asym == 1/15)) + 
  geom_line(aes(x=inc.age, y=inc.orapp.cohort, color = "Odds ratio approach (cohort)")) +
  geom_line(aes(x=inc.age, y=inc.detect.cohort, color = "Detection as etiology approach (cohort)")) +
  geom_line(aes(x=inc.age, y=inc.mod, color = "Model-based approach")) +
  scale_color_manual(values = c( "#1B9E77", "#D95F02", "#7570B3")) +
  scale_linetype_manual(values=c("dotted", "dotted", "dotted")) +
  ylab("Symptomatic norovirus disease incidence,\n per child-year") +
  xlab("Overall AGE incidence,\n per child-year") +
  theme_bw() +
  ggtitle("R0 = 1.2") +
  theme(legend.position=c(0.6,0.5),
        legend.title = element_blank()) -> main

ggplot(data=subset(Merge_subset, Merge_subset$s1 == 0.7 & 
                     Merge_subset$s2 == 0.3 &
                     Merge_subset$beta.noro == 0.4 &
                     Merge_subset$wan == 1/(365*5.1) &
                     Merge_subset$asym == 1/15)) + 
  geom_line(aes(x=inc.age, y=inc.orapp.cohort, color = "Odds ratio approach (cohort)")) +
  geom_line(aes(x=inc.age, y=inc.detect.cohort, color = "Detection as etiology approach (cohort)")) +
  geom_line(aes(x=inc.age, y=inc.mod, color = "Model-based approach")) +
  scale_color_manual(values = c( "#1B9E77", "#D95F02", "#7570B3")) +
  scale_linetype_manual(values=c("dotted", "dotted", "dotted")) +
  ylab("Symptomatic norovirus disease incidence,\n per child-year") +
  xlab("Overall AGE incidence,\n per child-year") +
  theme_bw() +
  ggtitle("asym = 15d") +
  theme(legend.position=c(0.6,0.5),
        legend.title = element_blank()) -> asym15d

# if you bump up asymptomatic period to 30d, get OR < 1 --> NEGATIVE incidence
ggplot(data=subset(Merge_subset, Merge_subset$s1 == 0.7 & 
                     Merge_subset$s2 == 0.3 &
                     Merge_subset$beta.noro == 0.4 &
                     Merge_subset$wan == 1/(365*5.1) &
                     Merge_subset$asym == 1/30)) + 
  geom_line(aes(x=inc.age, y=inc.orapp.cohort, color = "Odds ratio approach (cohort)")) +
  geom_line(aes(x=inc.age, y=inc.detect.cohort, color = "Detection as etiology approach (cohort)")) +
  geom_line(aes(x=inc.age, y=inc.mod, color = "Model-based approach")) +
  scale_color_manual(values = c( "#1B9E77", "#D95F02", "#7570B3")) +
  scale_linetype_manual(values=c("dotted", "dotted", "dotted")) +
  ylab("Symptomatic norovirus disease incidence,\n per child-year") +
  xlab("Overall AGE incidence,\n per child-year") +
  theme_bw() +
  ggtitle("asym = 30d") +
  theme(legend.position=c(0.6,0.5),
        legend.title = element_blank()) -> asym30d



library(ggpubr)
asym.change <- ggarrange(asym15d, asym30d,
                         nrow = 1, ncol = 2)
asym.change


ggplot(data=subset(Merge_subset, Merge_subset$s1 == 0.7 & 
                     Merge_subset$s2 == 0.3 &
                     Merge_subset$beta.noro == 0.4 &
                     Merge_subset$wan == 1/(365*5.1) &
                     Merge_subset$asym == 1/15)) + 
  geom_line(aes(x=inc.age, y=inc.orapp.cross, color = "Odds ratio approach (cross)")) +
  geom_line(aes(x=inc.age, y=inc.detect.cross, color = "Detection as etiology approach (cross)")) +
  geom_line(aes(x=inc.age, y=inc.mod, color = "Model-based approach")) +
  scale_color_manual(values = c( "#1B9E77", "#D95F02", "#7570B3")) +
  scale_linetype_manual(values=c("dotted", "dotted", "dotted")) +
  ylab("Symptomatic norovirus disease incidence,\n per child-year") +
  xlab("Overall AGE incidence,\n per child-year") +
  theme_bw() +
  ggtitle("R0 = 1.2") +
  theme(legend.position=c(0.6,0.3),
        legend.title = element_blank()) -> main.cross


cohort.cross <- ggarrange(main, main.cross,
                          nrow = 1, ncol = 2)
cohort.cross


# change beta
ggplot(data=subset(Merge_subset, Merge_subset$s1 == 0.7 & 
                     Merge_subset$s2 == 0.3 &
                     Merge_subset$beta.noro == 0.5 &
                     Merge_subset$wan == 1/(365*5.1) &
                     Merge_subset$asym == 1/15)) + 
  geom_line(aes(x=inc.age, y=inc.orapp.cohort, color = "Odds ratio approach (cohort)")) +
  geom_line(aes(x=inc.age, y=inc.detect.cohort, color = "Detection as etiology approach (cohort)")) +
  geom_line(aes(x=inc.age, y=inc.mod, color = "Model-based approach")) +
  scale_color_manual(values = c( "#1B9E77", "#D95F02", "#7570B3")) +
  ylab("Symptomatic norovirus disease incidence,\n per child-year") +
  xlab("Overall AGE incidence,\n per child-year") +
  ggtitle("R0 = 1.5") +  
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_blank())-> beta0.5

ggplot(data=subset(Merge_subset, Merge_subset$s1 == 0.7 & 
                     Merge_subset$s2 == 0.3 &
                     Merge_subset$beta.noro == 0.6 &
                     Merge_subset$wan == 1/(365*5.1) &
                     Merge_subset$asym == 1/15)) +   
  geom_line(aes(x=inc.age, y=inc.orapp.cohort, color = "Odds ratio approach (cohort)")) +
  geom_line(aes(x=inc.age, y=inc.detect.cohort, color = "Detection as etiology approach (cohort)")) +
  geom_line(aes(x=inc.age, y=inc.mod, color = "Model-based approach")) +
  scale_color_manual(values = c( "#1B9E77", "#D95F02", "#7570B3")) +
  ylab("Symptomatic norovirus disease incidence,\n per child-year") +
  xlab("Overall AGE incidence,\n per child-year") +
  ggtitle("R0 = 1.8") +
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_blank()) -> beta0.6

ggplot(data=subset(Merge_subset, Merge_subset$s1 == 0.7 & 
                     Merge_subset$s2 == 0.3 &
                     Merge_subset$beta.noro == 0.7 &
                     Merge_subset$wan == 1/(365*5.1) &
                     Merge_subset$asym == 1/15)) + 
  geom_line(aes(x=inc.age, y=inc.orapp.cohort, color = "Odds ratio approach (cohort)")) +
  geom_line(aes(x=inc.age, y=inc.detect.cohort, color = "Detection as etiology approach (cohort)")) +
  geom_line(aes(x=inc.age, y=inc.mod, color = "Model-based approach")) +
  scale_color_manual(values = c( "#1B9E77", "#D95F02", "#7570B3")) +
  ylab("Symptomatic norovirus disease incidence,\n per child-year") +
  xlab("Overall AGE incidence,\n per child-year") +
  ggtitle("R0 = 2.1") +
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_blank()) -> beta0.7

library(ggpubr)
changebeta <- ggarrange(main, beta0.5, beta0.6, beta0.7,
                        nrow = 2, ncol = 2)
changebeta 

#supplemental figure - changing norovirus beta (effective contact rate) and thus FOI, all else constant
ggsave(changebeta, file='/Users/kristinnelson/Box Sync/NoVaMod/Case-control attributable fractions/changebeta.pdf', 
       width = 6, height = 6, units = "in",  dpi = 300)

ggsave(main, file='/Users/kristinnelson/Box Sync/NoVaMod/Case-control attributable fractions/fig1.pdf', 
       width = 6, height = 6, units = "in",  dpi = 300)



# Compare bias in incidence, AF, OR
ggplot(data=subset(Merge_subset, Merge_subset$s1 == 0.7 & 
                     Merge_subset$s2 == 0.3 &
                     Merge_subset$beta.noro == 0.4 &
                     Merge_subset$wan == 1/(365*5.1) &
                     Merge_subset$asym == 1/15)) + 
  geom_line(aes(x=inc.age, y=inc.orapp.cohort, color = "Odds ratio approach (cohort)")) +
  geom_line(aes(x=inc.age, y=inc.detect.cohort, color = "Detection as etiology approach (cohort)")) +
  geom_line(aes(x=inc.age, y=inc.mod, color = "Model-based approach")) +
  scale_color_manual(values = c( "#1B9E77", "#D95F02", "#7570B3")) +
  ylab("Symptomatic norovirus disease incidence,\n per child-year") +
  xlab("Overall AGE incidence,\n per child-year") +
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_blank()) -> ai

ggplot(data=subset(Merge_subset, Merge_subset$s1 == 0.7 & 
                     Merge_subset$s2 == 0.3 &
                     Merge_subset$beta.noro == 0.4 &
                     Merge_subset$wan == 1/(365*5.1) &
                     Merge_subset$asym == 1/15)) + 
  geom_line(aes(x=inc.age, y=af.orapp.cohort, color = "Odds ratio approach (cohort)")) +
  geom_line(aes(x=inc.age, y=af.detect.cohort, color = "Detection as etiology approach (cohort)")) +
  geom_line(aes(x=inc.age, y=af.mod, color = "Model-based approach")) +
  scale_color_manual(values = c( "#1B9E77", "#D95F02", "#7570B3")) +
  ylab("Attributable fraction for norovirus") +
  xlab("Overall AGE incidence,\n per child-year") +
  theme_bw() +
  theme(legend.position=c(0.6,0.8),
        legend.title = element_blank()) -> af

ggplot(data=subset(Merge_subset, Merge_subset$s1 == 0.7 & 
                     Merge_subset$s2 == 0.3 &
                     Merge_subset$beta.noro == 0.4 &
                     Merge_subset$wan == 1/(365*5.1) &
                     Merge_subset$asym == 1/15)) +  
  geom_line(aes(x=inc.age, y=or.orapp.cohort), color = "#1B9E77") +
  ylab("Odds ratio") +
  xlab("Overall AGE incidence,\n per child-year") +
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_blank()) -> or

af.ai.or <- ggarrange(ai, af,
                      nrow = 1, ncol = 2)

af.ai.or

ggsave(af.ai.or, file='/Users/kristinnelson/Box Sync/NoVaMod/Case-control attributable fractions/af.ai.or.pdf', 
       width = 7, height = 4, units = "in",  dpi = 300)


# OR by compartment
library(wesanderson)
names(wes_palettes)
wes_palette("Zissou1", 3, type = "discrete") -> wes

# AF by compartment
ggplot(data=subset(Merge_subset, Merge_subset$s1 == 0.7 & 
                     Merge_subset$s2 == 0.3 &
                     Merge_subset$wan == 1/(365*5.1) &
                     Merge_subset$asym == 1/15)) +
  geom_line(aes(x=inc.age, y=af.orapp.cohort, color = "Odds ratio approach")) +
  geom_line(aes(x=inc.age, y=af.orapp.noreinf.cohort, color = "Exclude reinfection\n from controls")) +
  geom_line(aes(x=inc.age, y=af.orapp.noreinfpostsym.cohort, color = "Exclude reinfections and \n post-symptomatic shedders from controls")) +
  geom_line(aes(x=inc.age, y=af.mod, color = "Model-based approach")) +
  scale_color_manual(values = c( "#3B9AB2", "#F21A00", "#EBCC2A", "black")) +
  ylab("Attributable fraction") +
  xlab("Overall AGE incidence,\n per child-year") +
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_blank()) -> af.comp

# AI by compartment
ggplot(data=subset(Merge_subset, Merge_subset$s1 == 0.7 & 
                     Merge_subset$s2 == 0.3 &
                     Merge_subset$wan == 1/(365*5.1) &
                     Merge_subset$asym == 1/15)) +
  geom_line(aes(x=inc.age, y=inc.orapp.cohort, color = "Odds ratio approach")) +
  geom_line(aes(x=inc.age, y=inc.orapp.noreinf.cohort, color = "Exclude reinfection\n from controls")) +
  geom_line(aes(x=inc.age, y=inc.orapp.noreinfpostsym.cohort, color = "Exclude reinfections and \n post-symptomatic shedders from controls")) +
  geom_line(aes(x=inc.age, y=inc.mod, color = "Model-based approach")) +
  scale_color_manual(values = c( "#3B9AB2", "#F21A00", "#EBCC2A", "black")) +
  ylab("Symptomatic norovirus disease incidence,\n per child-year") +
  xlab("Overall AGE incidence,\n per child-year") +
  theme_bw() +
  theme(legend.position=c(0.35,0.85),
        legend.title = element_blank()) -> ai.comp

af.ai.or.comp <- ggarrange(ai.comp, af.comp,
                           nrow = 1, ncol = 2)

af.ai.or.comp

ggsave(af.ai.or.comp, file='/Users/kristinnelson/Box Sync/NoVaMod/Case-control attributable fractions/af.ai.or.comp.pdf', 
       width = 8, height = 3, units = "in",  dpi = 300)


# AF by compartment - change s1
ggplot(data=subset(Merge_subset, Merge_subset$s1 == 0.5 & 
                     Merge_subset$s2 == 0.3 &
                     Merge_subset$wan == 1/(365*5.1) &
                     Merge_subset$asym == 1/15)) +
  geom_line(aes(x=inc.age, y=af.orapp.cohort, color = "Odds ratio approach")) +
  geom_line(aes(x=inc.age, y=af.orapp.noreinf.cohort, color = "Exclude reinfection\n from controls")) +
  geom_line(aes(x=inc.age, y=af.orapp.noreinfpostsym.cohort, color = "Exclude reinfections and \n post-symptomatic shedders from controls")) +
  geom_line(aes(x=inc.age, y=af.mod, color = "Model-based approach")) +
  scale_color_manual(values = c( "#3B9AB2", "#F21A00", "#EBCC2A", "black")) +
  ylab("Attributable fraction") +
  xlab("Overall AGE incidence,\n per child-year") +
  theme_bw() +
  ggtitle("s1 = 0.5") +
  theme(legend.position="none",
        legend.title = element_blank()) -> af.comp.5

# AI by compartment
ggplot(data=subset(Merge_subset, Merge_subset$s1 == 0.7 & 
                     Merge_subset$s2 == 0.3 &
                     Merge_subset$wan == 1/(365*5.1) &
                     Merge_subset$asym == 1/15)) +
  geom_line(aes(x=inc.age, y=af.orapp.cohort, color = "Odds ratio approach")) +
  geom_line(aes(x=inc.age, y=af.orapp.noreinf.cohort, color = "Exclude reinfection\n from controls")) +
  geom_line(aes(x=inc.age, y=af.orapp.noreinfpostsym.cohort, color = "Exclude reinfections and \n post-symptomatic shedders from controls")) +
  geom_line(aes(x=inc.age, y=af.mod, color = "Model-based approach")) +
  scale_color_manual(values = c( "#3B9AB2", "#F21A00", "#EBCC2A", "black")) +
  ylab("Attributable fraction") +
  xlab("Overall AGE incidence,\n per child-year") +
  theme_bw() +
  ggtitle("s1 = 0.7") +
  theme(legend.position="none",
        legend.title = element_blank()) -> af.comp.7

ggplot(data=subset(Merge_subset, Merge_subset$s1 == 0.9 & 
                     Merge_subset$s2 == 0.3 &
                     Merge_subset$wan == 1/(365*5.1) &
                     Merge_subset$asym == 1/15)) +
  geom_line(aes(x=inc.age, y=af.orapp.cohort, color = "Odds ratio approach")) +
  geom_line(aes(x=inc.age, y=af.orapp.noreinf.cohort, color = "Exclude reinfection\n from controls")) +
  geom_line(aes(x=inc.age, y=af.orapp.noreinfpostsym.cohort, color = "Exclude reinfections and \n post-symptomatic shedders from controls")) +
  geom_line(aes(x=inc.age, y=af.mod, color = "Model-based approach")) +
  scale_color_manual(values = c( "#3B9AB2", "#F21A00", "#EBCC2A", "black")) +
  ylab("Attributable fraction") +
  xlab("Overall AGE incidence,\n per child-year") +
  theme_bw() +
  ggtitle("s1 = 0.9") +
  theme(legend.position="none",
        legend.title = element_blank()) -> af.comp.9


af.ai.or.comp <- ggarrange(af.comp.5, af.comp.7, af.comp.9,
                           nrow = 1, ncol = 3)

af.ai.or.comp


# AF by compartment - change s2
ggplot(data=subset(Merge_subset, Merge_subset$s1 == 0.7 & 
                     Merge_subset$s2 == 0.3 &
                     Merge_subset$wan == 1/(365*5.1) &
                     Merge_subset$asym == 1/15)) +
  geom_line(aes(x=inc.age, y=af.orapp.cohort, color = "Odds ratio approach")) +
  geom_line(aes(x=inc.age, y=af.orapp.noreinf.cohort, color = "Exclude reinfection\n from controls")) +
  geom_line(aes(x=inc.age, y=af.orapp.noreinfpostsym.cohort, color = "Exclude reinfections and \n post-symptomatic shedders from controls")) +
  geom_line(aes(x=inc.age, y=af.mod, color = "Model-based approach")) +
  scale_color_manual(values = c( "#3B9AB2", "#F21A00", "#EBCC2A", "black")) +
  ylab("Attributable fraction") +
  xlab("Overall AGE incidence,\n per child-year") +
  theme_bw() +
  ggtitle("s2 = 0.3") +
  theme(legend.position="none",
        legend.title = element_blank()) -> af.comp.3

ggplot(data=subset(Merge_subset, Merge_subset$s1 == 0.7 & 
                     Merge_subset$s2 == 0.5 &
                     Merge_subset$wan == 1/(365*5.1) &
                     Merge_subset$asym == 1/15)) +
  geom_line(aes(x=inc.age, y=af.orapp.cohort, color = "Odds ratio approach")) +
  geom_line(aes(x=inc.age, y=af.orapp.noreinf.cohort, color = "Exclude reinfection\n from controls")) +
  geom_line(aes(x=inc.age, y=af.orapp.noreinfpostsym.cohort, color = "Exclude reinfections and \n post-symptomatic shedders from controls")) +
  geom_line(aes(x=inc.age, y=af.mod, color = "Model-based approach")) +
  scale_color_manual(values = c( "#3B9AB2", "#F21A00", "#EBCC2A", "black")) +
  ylab("Attributable fraction") +
  xlab("Overall AGE incidence,\n per child-year") +
  theme_bw() +
  ggtitle("s2 = 0.5") +
  theme(legend.position="none",
        legend.title = element_blank()) -> af.comp.5

ggplot(data=subset(Merge_subset, Merge_subset$s1 == 0.7 & 
                     Merge_subset$s2 == 0.7 &
                     Merge_subset$wan == 1/(365*5.1) &
                     Merge_subset$asym == 1/15)) +
  geom_line(aes(x=inc.age, y=af.orapp.cohort, color = "Odds ratio approach")) +
  geom_line(aes(x=inc.age, y=af.orapp.noreinf.cohort, color = "Exclude reinfection\n from controls")) +
  geom_line(aes(x=inc.age, y=af.orapp.noreinfpostsym.cohort, color = "Exclude reinfections and \n post-symptomatic shedders from controls")) +
  geom_line(aes(x=inc.age, y=af.mod, color = "Model-based approach")) +
  scale_color_manual(values = c( "#3B9AB2", "#F21A00", "#EBCC2A", "black")) +
  ylab("Attributable fraction") +
  xlab("Overall AGE incidence,\n per child-year") +
  theme_bw() +
  ggtitle("s2 = 0.7") +
  theme(legend.position="none",
        legend.title = element_blank()) -> af.comp.7


af.ai.or.comp <- ggarrange(af.comp.3, af.comp.5, af.comp.7,
                           nrow = 1, ncol = 3)

af.ai.or.comp


ggplot(data=subset(Merge_subset, Merge_subset$s1 == 0.7 & 
                     Merge_subset$s2 == 0.3 &
                     Merge_subset$wan == 1/(365*5.1) &
                     Merge_subset$asym == 1/15)) +
  geom_line(aes(x=inc.age, y=af.orapp.cohort, color = "Odds ratio approach")) +
  geom_line(aes(x=inc.age, y=af.orapp.noreinf.cohort, color = "Exclude reinfection\n from controls")) +
  geom_line(aes(x=inc.age, y=af.orapp.noreinfpostsym.cohort, color = "Exclude reinfections and \n post-symptomatic shedders from controls")) +
  geom_line(aes(x=inc.age, y=af.mod, color = "Model-based approach")) +
  scale_color_manual(values = c( "#3B9AB2", "#F21A00", "#EBCC2A", "black")) +
  ylab("Attributable fraction") +
  xlab("Overall AGE incidence,\n per child-year") +
  theme_bw() +
  ggtitle("asym = 15d") +
  theme(legend.position="none",
        legend.title = element_blank()) -> asym.15

ggplot(data=subset(Merge_subset, Merge_subset$s1 == 0.7 & 
                     Merge_subset$s2 == 0.3 &
                     Merge_subset$wan == 1/(365*5.1) &
                     Merge_subset$asym == 1/30)) +
  geom_line(aes(x=inc.age, y=af.orapp.cohort, color = "Odds ratio approach")) +
  geom_line(aes(x=inc.age, y=af.orapp.noreinf.cohort, color = "Exclude reinfection\n from controls")) +
  geom_line(aes(x=inc.age, y=af.orapp.noreinfpostsym.cohort, color = "Exclude reinfections and \n post-symptomatic shedders from controls")) +
  geom_line(aes(x=inc.age, y=af.mod, color = "Model-based approach")) +
  scale_color_manual(values = c( "#3B9AB2", "#F21A00", "#EBCC2A", "black")) +
  ylab("Attributable fraction") +
  xlab("Overall AGE incidence,\n per child-year") +
  theme_bw() +
  ggtitle("asym = 30d") +
  theme(legend.position="none",
        legend.title = element_blank()) -> asym.30


af.ai.or.comp <- ggarrange(asym.15, asym.30,
                           nrow = 1, ncol = 2)

af.ai.or.comp
# Look at impact of s1/s2

# Calculate percentage difference between inc.mod and inc.orapp, inc.detect

mutate(Merge_subset, pc.diff.inc.OR = 100*(inc.orapp.cohort-inc.mod)/inc.mod,
       pc.diff.inc.detect = 100*(inc.detect.cohort-inc.mod)/inc.mod,
       pc.diff.af.OR = 100*(af.orapp.cohort-af.mod)/af.mod,
       pc.diff.af.detect = 100*(af.detect.cohort-af.mod)/af.mod) -> Merge_subset

library(plotly)

#OR approach bias - s1, s2, asym, wan
fig <- plot_ly(y = Merge_subset$asym, x = Merge_subset$s1, z = Merge_subset$pc.diff.inc.OR , type = "contour", 
               colorscale = list(c(0, 0.5, 1), c('yellow', 'MediumAquaMarine', 'blue'))) 
fig <- fig %>% layout(xaxis=list(title="s1"), 
                      yaxis=list(title="s2"))
contour.orbias.inc.s1asym <- fig %>%   colorbar(title='Percent difference from modeled incidence')
contour.orbias.inc.s1asym

fig <- plot_ly(y = Merge_subset$asym, x = Merge_subset$s1, z = Merge_subset$pc.diff.af.OR , type = "contour", 
               colorscale = list(c(0, 0.5, 1), c('yellow', 'MediumAquaMarine', 'blue'))) 
fig <- fig %>% layout(xaxis=list(title="s1"), 
                      yaxis=list(title="s2"))
contour.orbias.af.s1asym <- fig %>%   colorbar(title='Percent difference from modeled attributable fraction')
contour.orbias.af.s1asym


fig <- plot_ly(y = Merge_subset$s2, x = Merge_subset$s1, z = Merge_subset$pc.diff.af.OR , type = "contour", 
               colorscale = list(c(0, 0.5, 1), c('yellow', 'MediumAquaMarine', 'blue'))) 
fig <- fig %>% layout(xaxis=list(title="s1"), 
                      yaxis=list(title="s2"))
contour.orbias.inc.s1s2 <- fig %>%   colorbar(title='Percent difference from modeled attributable fraction')
contour.orbias.inc.s1s2


ggsave(contour.orbias.af.s1asym, file='/Users/kristinnelson/Box Sync/NoVaMod/Case-control attributable fractions/afbias.contour.pdf', 
       width = 4, height = 3, units = "in",  dpi = 300)

ggsave(contour.orbias.inc.s1asym, file='/Users/kristinnelson/Box Sync/NoVaMod/Case-control attributable fractions/incbias.contour.pdf', 
       width = 4, height = 3, units = "in",  dpi = 300)

# Look at impact of asymptomatic shedding duration (asym)


fig <- plot_ly(y = Merge_subset$s2, x = Merge_subset$s1, z = Merge_subset$pc.diff.af.detect , type = "contour", 
               colorscale = list(c(0, 0.5, 1), c('yellow', 'MediumAquaMarine', 'blue'))) 
fig <- fig %>% layout(xaxis=list(title="s1"), 
                      yaxis=list(title="s2"))
contour.detectbias.inc.s1s2 <- fig %>%   colorbar(title='Percent difference from modeled attributable incidence')
contour.detectbias.inc.s1s2

fig <- plot_ly(y = Merge_subset$asym, x = Merge_subset$s1, z = Merge_subset$pc.diff.af.detect , type = "contour", 
               colorscale = list(c(0, 0.5, 1), c('yellow', 'MediumAquaMarine', 'blue'))) 
fig <- fig %>% layout(xaxis=list(title="s1"), 
                      yaxis=list(title="asym"))
contour.detectbias.inc.s1asym <- fig %>%   colorbar(title='Percent difference from modeled attributable incidence')
contour.detectbias.inc.s1asym


# this is keeping beta constant while increasing FOI.other, that may not make sense. create a new dataset where beta and FOI.other increase in tandem. 
## In US, there are an estimated 179mil cases of AGE, this is an incidence of 0.55 per person-year. (179mil/328mil) Using the base parameters for models, this is roughly equivalent to a FOI.other of 0.001-0.002. (Note this is entire pop, would be higher among under-5s)
## In Pakistan, AGE incidence from MAL-ED is estimated to be > 7 per person-year, over 7 times as high. This is roughly equivalent to an FOI.other of 0.02.## Now,
## Now, how to benchmark beta.noro? Can't do it with norovirus incidence, since we're operating under the assumption these are biased. Could use R0 - use beta.noro values that give R0s ranging from >1 to 4ish. 
## This gives beta.noro range from 0.4 to 0.8. 
## So, assume that beta.noro and FOI.other increase in tandem within these ranges (linear-ish)

subset(Merge_subset, 
       beta.noro == 0.4 & FOI.other == 0.001 |
         beta.noro == 0.5 & FOI.other == 0.005 |
         beta.noro == 0.6 & FOI.other == 0.01 |
         beta.noro == 0.7 & FOI.other == 0.015 |
         beta.noro == 0.8 & FOI.other == 0.02) -> Merge_subset
Merge_subset

ggplot(data=subset(Merge_subset, Merge_subset$s1 == 0.7 & 
                     Merge_subset$s2 == 0.3 &
                     Merge_subset$wan == 1/(365*5.1) &
                     Merge_subset$asym == 1/15)) + 
  geom_line(aes(x=inc.age, y=inc.orapp.cohort, color = "Odds ratio approach (cohort)")) +
  geom_line(aes(x=inc.age, y=inc.detect.cohort, color = "Detection as etiology approach (cohort)")) +
  geom_line(aes(x=inc.age, y=inc.mod, color = "Model-based approach")) +
  scale_color_manual(values = c( "#1B9E77", "#D95F02", "#7570B3")) +
  scale_linetype_manual(values=c("dotted", "dotted", "dotted")) +
  ylab("Symptomatic norovirus disease incidence,\n per child-year") +
  xlab("Overall AGE incidence,\n per child-year") +
  theme_bw() +
  ggtitle("s1 = 0.7, s2 = 0.3") +
  theme(legend.position=c(0.4,0.9),
        legend.title = element_blank()) -> s1.high

ggplot(data=subset(Merge_subset, Merge_subset$s1 == 0.5 & 
                     Merge_subset$s2 == 0.3 &
                     Merge_subset$wan == 1/(365*5.1) &
                     Merge_subset$asym == 1/15)) + 
  geom_line(aes(x=inc.age, y=inc.orapp.cohort, color = "Odds ratio approach (cohort)")) +
  geom_line(aes(x=inc.age, y=inc.detect.cohort, color = "Detection as etiology approach (cohort)")) +
  geom_line(aes(x=inc.age, y=inc.mod, color = "Model-based approach")) +
  scale_color_manual(values = c( "#1B9E77", "#D95F02", "#7570B3")) +
  scale_linetype_manual(values=c("dotted", "dotted", "dotted")) +
  ylab("Symptomatic norovirus disease incidence,\n per child-year") +
  xlab("Overall AGE incidence,\n per child-year") +
  theme_bw() +
  ggtitle("s1 = 0.5, s2 = 0.3") +
  theme(legend.position=c(0.4,0.9),
        legend.title = element_blank()) -> s1.low

ggplot(data=subset(Merge_subset, Merge_subset$s1 == 0.7 & 
                     Merge_subset$s2 == 0.5 &
                     Merge_subset$wan == 1/(365*5.1) &
                     Merge_subset$asym == 1/15)) + 
  geom_line(aes(x=inc.age, y=inc.orapp.cohort, color = "Odds ratio approach (cohort)")) +
  geom_line(aes(x=inc.age, y=inc.detect.cohort, color = "Detection as etiology approach (cohort)")) +
  geom_line(aes(x=inc.age, y=inc.mod, color = "Model-based approach")) +
  scale_color_manual(values = c( "#1B9E77", "#D95F02", "#7570B3")) +
  scale_linetype_manual(values=c("dotted", "dotted", "dotted")) +
  ylab("Symptomatic norovirus disease incidence,\n per child-year") +
  xlab("Overall AGE incidence,\n per child-year") +
  theme_bw() +
  ggtitle("s1 = 0.7, s2 = 0.5") +
  theme(legend.position=c(0.4,0.9),
        legend.title = element_blank()) -> s2.low


s1.s2.changingbeta <- ggarrange(s1.high, s1.low, s2.low,
                                nrow = 1, ncol = 3)
s1.s2.changingbeta

subset(Merge_subset, 
       beta.noro == 0.4 & FOI.other == 0.01 |
         beta.noro == 0.5 & FOI.other == 0.02 |
         beta.noro == 0.6 & FOI.other == 0.03 |
         beta.noro == 0.7 & FOI.other == 0.04 |
         beta.noro == 0.8 & FOI.other == 0.05) -> Merge_subset



# what is impact of s1/s2 on OR approach, by setting? 
ggplot(data=(Merge_subset)) + 
  geom_point(aes(x=inc.age, y=af.orapp.cohort, color=s1)) +
  xlab("Overall AGE incidence,\n per child-year") +
  ylab("Attributable fraction OR approach") +
  theme_bw() +
  theme(legend.position=c(0.8,0.7)) -> af.s1

ggplot(data=(Merge_subset)) + 
  geom_point(aes(x=inc.age, y=af.orapp.cohort, color=s2)) +
  xlab("Overall AGE incidence,\n per child-year") +
  ylab("Attributable fraction OR approach") +
  theme_bw() +
  theme(legend.position=c(0.8,0.7)) -> af.s2

ggplot(data=(Merge_subset)) + 
  geom_point(aes(x=inc.age, y=af.orapp.cohort, color=asym)) +
  xlab("Overall AGE incidence,\n per child-year") +
  ylab("Attributable fraction OR approach") +
  theme_bw() +
  theme(legend.position=c(0.8,0.7)) -> af.asym

ggplot(data=(Merge_subset)) + 
  geom_point(aes(x=inc.age, y=af.orapp.cohort, color=wan)) +
  xlab("Overall AGE incidence,\n per child-year") +
  ylab("Attributable fraction OR approach") +
  theme_bw() +
  theme(legend.position=c(0.8,0.7)) -> af.wan

s1.s2.asym.wan <- ggarrange(af.s1, af.s2, af.asym, af.wan,
                            nrow = 2, ncol = 2)

s1.s2.asym.wan

ggsave(s1.s2.asym.wan, file='/Users/kristinnelson/Box Sync/NoVaMod/Case-control attributable fractions/af.ai.or.pdf', 
       width = 8, height = 8, units = "in",  dpi = 300)


# pc diff instead

ggplot(data=subset(Merge_subset, beta.noro == 0.5)) + 
  geom_point(aes(x=inc.age, y=pc.diff.inc.OR, color=s1)) +
  xlab("Overall AGE incidence,\n per child-year") +
  ylab("Percent difference in incidence, modeled vs. OR approach") +
  theme_bw() +
  theme(legend.position=c(0.9,0.2)) -> af.s1

ggplot(data=subset(Merge_subset, beta.noro == 0.5)) + 
  geom_point(aes(x=inc.age, y=pc.diff.inc.OR, color=s2)) +
  xlab("Overall AGE incidence,\n per child-year") +
  ylab("Percent difference in incidence, modeled vs. OR approach") +
  theme_bw() +
  theme(legend.position=c(0.9,0.2)) -> af.s2

ggplot(data=subset(Merge_subset, beta.noro == 0.5)) + 
  geom_point(aes(x=inc.age, y=pc.diff.inc.OR, color=asym)) +
  xlab("Overall AGE incidence,\n per child-year") +
  ylab("Percent difference in incidence, modeled vs. OR approach") +
  theme_bw() +
  theme(legend.position=c(0.9, 0.2)) -> af.asym

ggplot(data=subset(Merge_subset, beta.noro == 0.5)) + 
  geom_point(aes(x=inc.age, y=pc.diff.inc.OR, color=wan)) +
  xlab("Overall AGE incidence,\n per child-year") +
  ylab("Percent difference in incidence, modeled vs. OR approach") +
  theme_bw() +
  theme(legend.position=c(0.9,0.2)) -> af.wan

ggplot(data=subset(Merge_subset, beta.noro == 0.5)) + 
  geom_point(aes(x=inc.age, y=pc.diff.inc.OR, color=prop.coinf.cohort)) +
  xlab("Overall AGE incidence,\n per child-year") +
  ylab("Percent difference in incidence, modeled vs. OR approach") +
  theme_bw() +
  theme(legend.position=c(0.8,0.2)) -> af.coinf

s1.s2.asym.wan.coinf <- ggarrange(af.s1, af.s2, af.asym, af.wan, af.coinf,
                                  nrow = 2, ncol = 3)

s1.s2.asym.wan.coinf

ggsave(s1.s2.asym.wan, file='/Users/kristinnelson/Box Sync/NoVaMod/Case-control attributable fractions/af.ai.or.pdf', 
       width = 8, height = 8, units = "in",  dpi = 300)



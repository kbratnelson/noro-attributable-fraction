rm(list=ls()) #Clear workspace
library(dplyr)
library(here)
library(tidyverse)
library(cowplot)
library(ggpubr)

load("/Users/kristinnelson/Box Sync/NoVaMod//Case-control attributable fractions/Results/simstudy/simsweep021920.RData")

# create variables for further analysis/figures

mutate(Merge, 
       
       # percent difference, OR and DE methods vs MB method
       pc.diff.ai.OR.MBref = 100*(ai.orapp.cross-ai.mb)/ai.mb,
       pc.diff.ai.DE.MBref = 100*(ai.detect.cross-ai.mb)/ai.mb,
       pc.diff.af.OR.MBref = 100*(af.orapp.cross-af.mb)/af.mb,
       pc.diff.af.DE.MBref = 100*(af.detect.cross-af.mb)/af.mb,
       
       # percent difference, OR and DE methods vs unbiased OR and DE methods
       pc.diff.ai.OR.ORref = 100*(ai.orapp.cross-ai.orapp.unbiased.cross)/ai.orapp.unbiased.cross,
       pc.diff.ai.DE.DEref = 100*(ai.detect.cross-ai.detect.unbiased.cross)/ai.detect.unbiased.cross,
       pc.diff.af.OR.ORref = 100*(af.orapp.cross-ai.orapp.unbiased.cross)/ai.orapp.unbiased.cross,
       pc.diff.af.DE.DEref = 100*(af.detect.cross-ai.detect.unbiased.cross)/ai.detect.unbiased.cross,

       # percent differences OR method biased vs. OR method unbiased by compartment (FOR RESULTS TEXT)
       pc.diff.reinf.ORref = (ai.orapp.noreinf.cross-ai.orapp.unbiased.cross) / ai.orapp.unbiased.cross,
       pc.diff.postsym.ORref = (ai.orapp.nopostsym.cross-ai.orapp.unbiased.cross) / ai.orapp.unbiased.cross,
       pc.diff.primasym.ORref = (ai.orapp.noprimasym.cross-ai.orapp.unbiased.cross) / ai.orapp.unbiased.cross,
       pc.diff.coinf.ORref = (ai.orapp.nocoinf.cross-ai.orapp.unbiased.cross) / ai.orapp.unbiased.cross,
       
       # Version 2 - percent differences OR method biased vs. OR method unbiased (FOR GROUPED BAR CHARTS) - finish calculating in next 'mutate' step
       # This version take difference between biased OR and OR without given compartment (e.g. Iar), which gives change in incidence due to single bias source
       bias.reinf.ORref = ai.orapp.cross-ai.orapp.noreinf.cross,
       bias.postsym.ORref = ai.orapp.cross-ai.orapp.nopostsym.cross,
       bias.primasym.ORref = ai.orapp.cross-ai.orapp.noprimasym.cross,
       bias.coinf.ORref = ai.orapp.cross-ai.orapp.nocoinf.cross,   
       
       bias.ORMB.MBref = ai.orapp.unbiased.cross-ai.mb,   
       
       
       # percent different between OR method unbiased and MB method   
       pc.diff.unbiased.MBref = (ai.orapp.unbiased.cross-ai.mb)/ai.mb

       
       ) -> Merge

# Now, create variables for bar graphs (need to do this in a separate step from mutate command above for sum function to work properly)

mutate(Merge, 
       
       # Version 2 - percent differences OR method biased vs. OR method unbiased (FOR GROUPED BAR CHARTS) - finish calculating in next 'mutate' step
       # After getting change in incidence due to single bias source in above step, now calculate total bias (absolute value). Note that some bias is negative and some positive. 
       # This is NOT net bias, but absolute value. 
       sum = abs(bias.reinf.ORref) + abs(bias.postsym.ORref) + abs(bias.primasym.ORref) + abs(bias.coinf.ORref),
       
       # same for MB method - this time, add in MB/OR bias to get total
       sum.MB = abs(bias.reinf.ORref) + abs(bias.postsym.ORref) + abs(bias.primasym.ORref) + abs(bias.coinf.ORref) + abs(bias.ORMB.MBref),
  
       
       # Now, divide bias from each course by total absolute value of bias. This gives proportions that, when you take the absolute value of all of them, add up to 1. 
       # Use these proportions for stacked bar.
       pc.diff.reinf.ORref.2 = bias.reinf.ORref/sum,
       pc.diff.postsym.ORref.2 = bias.postsym.ORref/sum,
       pc.diff.primasym.ORref.2 = bias.primasym.ORref/sum,
       pc.diff.coinf.ORref.2 = bias.coinf.ORref/sum,

       # same for MB method, except also account for different in between OR/MB method. There are 5 sources of bias here, not 4 like above!
       pc.diff.reinf.MBref.2 = bias.reinf.ORref/sum.MB,
       pc.diff.postsym.MBref.2 = bias.postsym.ORref/sum.MB,
       pc.diff.primasym.MBref.2 = bias.primasym.ORref/sum.MB,
       pc.diff.coinf.MBref.2 = bias.coinf.ORref/sum.MB,
       pc.diff.ORMB.MBref.2 = bias.ORMB.MBref/sum.MB) -> Merge 

####################################################################################
#Subset data to sims with reasonable values
####################################################################################

# Odds ratios, under reasonable incidence values
subset(Merge,
       ai.age > 0.1 & ai.age < 7 &
         ai.mb > 0.001 & ai.mb < 10) -> Merge_subset

table(Merge_subset$beta.noro)


####################################################################################
#Results stats
####################################################################################

# Check coinfection rates

median(Merge_subset$prop.coinf.cross) # 0.41
range(Merge_subset$prop.coinf.cross) # 0.02 - 0.87


#median(Merge_subset$or.orapp.cohort) # 3.1
#range(Merge_subset$or.orapp.cohort) # 0.5 - 159
#quantile(Merge_subset$or.orapp.cohort) # IQR: 1.7 - 7.2

median(Merge_subset$or.orapp.cross) # 3.7
range(Merge_subset$or.orapp.cross) # 1.5 - 108
quantile(Merge_subset$or.orapp.cross) # IQR: 2.5 - 7.3

# Percent difference between OR/detect-as-etiology methods and model-based
# Asssume baseline values for s1 and s2 and wan and asym
subset(Merge, Merge$s1 == 0.5 & 
         Merge$s2 == 0.3 &
         Merge$wan == 1/(365*5.1) &
         Merge$asym == 1/15) -> Merge_subset2

#Create subsets for LOW, MODERATE, HIGH incidence scenarios to show in boxplots below

# LOW INCIDENCE 
subset(Merge_subset2, Merge_subset2$beta.noro==0.45) -> Merge_subset3

# MODERATE INCIDENCE 
subset(Merge_subset2, Merge_subset2$beta.noro==0.55) -> Merge_subset4

# HIGH INCIDENCE
subset(Merge_subset2, Merge_subset2$beta.noro==0.65) -> Merge_subset5

                                                         


# Reference is unbiased OR/DE method - RESULTS PARA 2
median(Merge_subset2$pc.diff.ai.DE.DEref) # -101%
range(Merge_subset2$pc.diff.ai.DE.DEref) # 7% to 192%
median(Merge_subset2$pc.diff.ai.OR.ORref) # -32%
range(Merge_subset2$pc.diff.ai.OR.ORref) # -5% to -37% 

median(Merge_subset2$pc.diff.ai.DE.MBref) # -9%
range(Merge_subset2$pc.diff.ai.DE.MBref) # -38% to 34%
median(Merge_subset2$pc.diff.ai.OR.MBref) # -38%
range(Merge_subset2$pc.diff.ai.OR.MBref) # -48% to -16%

Merge_subset2$pc.diff.ai.DE.MBref # -36% 


# Percent difference between unbiased OR approach estimates and biased OR approach

# RESULTS PARA 3
Merge_subset3$ai.age
Merge_subset3[,c("pc.diff.primasym.ORref.2", "pc.diff.postsym.ORref.2", "pc.diff.reinf.ORref.2", "pc.diff.coinf.ORref.2")]

Merge_subset4$ai.age
Merge_subset4[,c("pc.diff.primasym.ORref.2", "pc.diff.postsym.ORref.2", "pc.diff.reinf.ORref.2", "pc.diff.coinf.ORref.2")]

Merge_subset5$ai.age
Merge_subset5[,c("pc.diff.primasym.ORref.2", "pc.diff.postsym.ORref.2", "pc.diff.reinf.ORref.2", "pc.diff.coinf.ORref.2")]



####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################



# FIGURE 2: Compare attribution approaches - FULL

ggplot(data=subset(Merge, Merge$s1 == 0.7 & 
                     Merge$s2 == 0.3 &
                     Merge$wan == 1/(365*5.1) &
                     Merge$asym == 1/15)) + 
  geom_point(aes(x=ai.age, y=ai.orapp.unbiased.cross+0.75, color = "OR approach, unbiased"), shape = 15, size = 3) + #jittered so as not to fall on top of other line
  geom_point(aes(x=ai.age, y=ai.orapp.cross, color = "OR approach"), shape = 15, size = 3) +
  geom_point(aes(x=ai.age, y=ai.detect.unbiased.cross, color = "DE approach, unbiased"), shape = 16, size = 3) +
  geom_point(aes(x=ai.age, y=ai.detect.cross, color = "DE approach"), shape = 16, size = 3) +
  geom_point(aes(x=ai.age, y=ai.mb, color = "MB approach"), shape = 17, size = 3) +
  geom_line(aes(x=ai.age, y=ai.orapp.unbiased.cross+0.75, color = "OR approach, unbiased", linetype = "unbiased")) + #jittered so as not to fall on top of other line
  geom_line(aes(x=ai.age, y=ai.orapp.cross, color = "OR approach",  linetype = "biased")) +
  geom_line(aes(x=ai.age, y=ai.detect.unbiased.cross, color = "DE approach, unbiased", linetype = "unbiased")) +
  geom_line(aes(x=ai.age, y=ai.detect.cross, color = "DE approach", linetype = "biased")) +
  geom_line(aes(x=ai.age, y=ai.mb, color = "MB approach", linetype = "unbiased")) +
  scale_color_manual(values = c( "#1B9E77", "#1B9E77", "#D95F02", "#7570B3", "#7570B3")) +
  scale_linetype_manual(values=c("dotted", "solid"), guide = "none") +
  guides(color = guide_legend(override.aes = list(linetype = c("dotted", "solid", "solid", "dotted", "solid"), shape = c(16, 16, 17, 15, 15)))) +
  ylab("Norovirus-attributed AGE incidence") +
  xlab("Overall AGE incidence,\n per child-year") +
  theme_bw() +
  theme(legend.position=c(0.4,0.8),
        #legend.position="none",
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        axis.title=element_text(size=12),
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 12)) -> ai.biasedunbiased

ggsave(ai.biasedunbiased, file='or.de.mb.compare.biasedunbiased.pdf', width = 6, height = 6, units = "in",  dpi = 300)


ggplot(data=subset(Merge, Merge$s1 == 0.7 & 
                     Merge$s2 == 0.3 &
                     Merge$wan == 1/(365*5.1) &
                     Merge$asym == 1/15)) + 
  geom_point(aes(x=ai.age, y=af.orapp.unbiased.cross+0.75, color = "OR approach, unbiased"), shape = 15, size = 3) + #jittered so as not to fall on top of other line
  geom_point(aes(x=ai.age, y=af.orapp.cross, color = "OR approach"), shape = 15, size = 3) +
  geom_point(aes(x=ai.age, y=af.detect.unbiased.cross, color = "DE approach, unbiased"), shape = 16, size = 3) +
  geom_point(aes(x=ai.age, y=af.detect.cross, color = "DE approach"), shape = 16, size = 3) +
  geom_point(aes(x=ai.age, y=af.mb, color = "MB approach"), shape = 17, size = 3) +
  geom_line(aes(x=ai.age, y=af.orapp.unbiased.cross+0.75, color = "OR approach, unbiased", linetype = "unbiased")) + #jittered so as not to fall on top of other line
  geom_line(aes(x=ai.age, y=af.orapp.cross, color = "OR approach",  linetype = "biased")) +
  geom_line(aes(x=ai.age, y=af.detect.unbiased.cross, color = "DE approach, unbiased", linetype = "unbiased")) +
  geom_line(aes(x=ai.age, y=af.detect.cross, color = "DE approach", linetype = "biased")) +
  geom_line(aes(x=ai.age, y=af.mb, color = "MB approach", linetype = "unbiased")) +
  scale_color_manual(values = c( "#1B9E77", "#1B9E77", "#D95F02", "#7570B3", "#7570B3")) +
  scale_linetype_manual(values=c("dotted", "solid"), guide = "none") +
  guides(color = guide_legend(override.aes = list(linetype = c("dotted", "solid", "solid", "dotted", "solid"), shape = c(16, 16, 17, 15, 15)))) +
  ylab("Norovirus attributable fraction") +
  xlab("Overall AGE incidence,\n per child-year") +
  theme_bw() +
  theme(legend.position="none",
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        axis.title=element_text(size=12),
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 12)) -> af.biasedunbiased


fig2_unbiasedvsbias <- ggarrange(ai.biasedunbiased, af.biasedunbiased,
                        nrow = 1, ncol = 2)

fig2_unbiasedvsbias

ggsave(fig2_unbiasedvsbias, file='or.de.mb.compare.biasedunbiased.afai.pdf', width = 12, height = 6, units = "in",  dpi = 300)


# FIGURE 2: Compare attribution approaches - SIMPLIFIED

ggplot(data=subset(Merge, Merge$s1 == 0.7 & 
                     Merge$s2 == 0.3 &
                     Merge$wan == 1/(365*5.1) &
                     Merge$asym == 1/15)) + 
  geom_point(aes(x=ai.age, y=ai.orapp.cross, color = "OR approach"), shape = 15, size = 3) +
  geom_point(aes(x=ai.age, y=ai.detect.cross, color = "DE approach"), shape = 16, size = 3) +
  geom_point(aes(x=ai.age, y=ai.mb, color = "MB approach"), shape = 17, size = 3) +
  geom_line(aes(x=ai.age, y=ai.orapp.cross, color = "OR approach")) +
  geom_line(aes(x=ai.age, y=ai.detect.cross, color = "DE approach")) +
  geom_line(aes(x=ai.age, y=ai.mb, color = "MB approach")) +
  scale_color_manual(values = c( "#1B9E77", "#D95F02", "#7570B3")) +
  scale_linetype_manual(values=c("dotted", "solid"), guide = "none") +
  guides(color = guide_legend(override.aes = list(shape = c(16, 17, 15)))) +
  ylab("Norovirus-attributed AGE incidence") +
  xlab("Overall AGE incidence,\n per child-year") +
  theme_bw() +
  theme(legend.position=c(0.4,0.8),
        #legend.position="none",
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        axis.title=element_text(size=12),
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 12)) -> ai.biasedunbiased.simp

ai.biasedunbiased.simp

ggsave(ai.biasedunbiased.simp, file='ai.biasedunbiased.simp.pdf', width = 6, height = 6, units = "in",  dpi = 300)



# FIGURE 3: SHOW CONTRIBUTIONS TO BIAS IN OR APPROACH

# colors
library(wesanderson)
#names(wes_palettes)
wes_palette("FantasticFox1", 5, type = "discrete") -> wes
#"#DD8D29" "#E2D200" "#46ACC8" "#E58601" "#B40F20"

inc <- c( rep("1-low" , 4) , rep("2-medium" , 4) , rep("3-high" , 4) )
bias <- rep( c("primary asymptomatic infection", "post-infection shedding", "asymptomatic reinfection", "coinfection") , 3)
value <- c( Merge_subset3$pc.diff.primasym.ORref.2, Merge_subset3$pc.diff.postsym.ORref.2, Merge_subset3$pc.diff.reinf.ORref.2, Merge_subset3$pc.diff.coinf.ORref.2,
            Merge_subset4$pc.diff.primasym.ORref.2, Merge_subset4$pc.diff.postsym.ORref.2, Merge_subset4$pc.diff.reinf.ORref.2, Merge_subset4$pc.diff.coinf.ORref.2,
            Merge_subset5$pc.diff.primasym.ORref.2, Merge_subset5$pc.diff.postsym.ORref.2, Merge_subset5$pc.diff.reinf.ORref.2, Merge_subset5$pc.diff.coinf.ORref.2)

data <- data.frame(inc, bias, value)

#Grouped bar, relative to unbiased OR method
ggplot(data, aes(fill=bias, y=100*value, x=inc)) + 
  geom_bar(position="dodge", stat="identity", width = 0.5) +
  scale_fill_manual(values = c("#E2D200", "#46ACC8", "#E58601", "#B40F20")) +
  theme_bw() +
  ylab("Percent contribution to bias") +
  xlab("Overall AGE incidence") +
  scale_x_discrete(labels=c("1-low" = "Low\n(50 per 100 child-years)", 
                            "2-medium" = "Medium\n(200 per 100 child-years)",
                            "3-high" = "High\n(350 per 100 child-years)")) +
  geom_hline(yintercept = 0, color = "black") +
  theme(legend.position=c(0.7,0.1),
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        axis.title=element_text(size=12),
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 12)) -> fig3.groupedbar.unbiasedORref

fig3.groupedbar.unbiasedORref

ggsave(fig3.groupedbar.unbiasedORref, file='orbias.groupedbar.pdf', width = 7, height = 7, units = "in",  dpi = 300)


#Grouped bar, relative to MB method

inc <- c( rep("1-low" , 5) , rep("2-medium" , 5) , rep("3-high" , 5) )
bias <- rep( c("primary asymptomatic infection", "post-infection shedding", "asymptomatic reinfection", "coinfection", "other biases"), 3)
value <- c( Merge_subset3$pc.diff.primasym.MBref.2, Merge_subset3$pc.diff.postsym.MBref.2, Merge_subset3$pc.diff.reinf.MBref.2, Merge_subset3$pc.diff.coinf.MBref.2, Merge_subset3$pc.diff.ORMB.MBref.2, 
            Merge_subset4$pc.diff.primasym.MBref.2, Merge_subset4$pc.diff.postsym.MBref.2, Merge_subset4$pc.diff.reinf.MBref.2, Merge_subset4$pc.diff.coinf.MBref.2, Merge_subset4$pc.diff.ORMB.MBref.2,
            Merge_subset5$pc.diff.primasym.MBref.2, Merge_subset5$pc.diff.postsym.MBref.2, Merge_subset5$pc.diff.reinf.MBref.2, Merge_subset5$pc.diff.coinf.MBref.2, Merge_subset5$pc.diff.ORMB.MBref.2)

data <- data.frame(inc, bias, value)


ggplot(data, aes(fill=bias, y=100*value, x=inc)) + 
  geom_bar(position="dodge", stat="identity", width = 0.5) +
  scale_fill_manual(values = c("#E2D200", "#46ACC8", "grey40", "#E58601", "#B40F20")) +
  theme_bw() +
  ylab("Percent contribution to bias") +
  xlab("Overall AGE incidence") +
  scale_x_discrete(labels=c("1-low" = "Low\n(50 per 100 child-years)", 
                            "2-medium" = "Medium\n(200 per 100 child-years)",
                            "3-high" = "High\n(350 per 100 child-years)")) +
  geom_hline(yintercept = 0, color = "black") +
  theme(legend.position=c(0.75,0.2),
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        axis.title=element_text(size=12),
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 12)) -> fig3.groupedbar.MBref

ggsave(fig3.groupedbar.MBref, file='orbias.groupedbar.MBref.pdf', width = 7, height = 7, units = "in",  dpi = 300)

fig3.groupedbar.MBref


# OTHER FIGURE: Stacked bar, compartment sizes
inc <- c( rep("1-low" , 9) , rep("2-medium" , 9) , rep("3-high" , 9) )
comp <- rep( c("1-S", "2-Is", "3-Ia", "4-R", "5-Isr", "6-Iar", "7-Ipsr", "8-Iother", "9-Icoinf") , 3)
value <- c( Merge_subset3$S, Merge_subset3$Is, Merge_subset3$Ia, Merge_subset3$R, Merge_subset3$Isr, Merge_subset3$Iar, Merge_subset3$Ispr, Merge_subset3$Iother, Merge_subset3$Icoinf,
            Merge_subset4$S, Merge_subset4$Is, Merge_subset4$Ia, Merge_subset4$R, Merge_subset4$Isr, Merge_subset4$Iar, Merge_subset4$Ispr, Merge_subset4$Iother, Merge_subset4$Icoinf,
            Merge_subset5$S, Merge_subset5$Is, Merge_subset5$Ia, Merge_subset5$R, Merge_subset5$Isr, Merge_subset5$Iar, Merge_subset5$Ispr, Merge_subset5$Iother, Merge_subset5$Icoinf)

data <- data.frame(inc, comp, value)

ggplot(data, aes(fill=comp, y=value, x=inc)) + 
  geom_bar(position="stack", stat="identity")




############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################


# SUPPLEMENTARY FIGURE - SENSITIVITY OF BIAS TO MODEL PARAMS

ggplot(data=subset(Merge_subset, beta.noro == 0.5), aes(x=factor(s1), y=((ai.mb-ai.orapp.cross)*100))) + 
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.05), color = "#E7B800", alpha = 0.8) +
  xlab("Proportion of primary infections symptomatic") +
  ylab("Difference in incidencerelative to \nMB approach, per 100 child-years") +
  theme_bw() +
  theme(legend.position=c(0.9,0.2)) -> af.s1

ggplot(data=subset(Merge_subset, beta.noro == 0.5), aes(x=factor(s2), y=((ai.mb-ai.orapp.cross)*100))) + 
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.05), color = "#E7B800", alpha = 0.8) +
  xlab("Proportion of reinfections symptomatic") +
  ylab("Difference in incidencerelative to \nMB approach, per 100 child-years") +
  theme_bw() +
  theme(legend.position=c(0.9,0.2)) -> af.s2

ggplot(data=subset(Merge_subset, beta.noro == 0.5), aes(x=factor(1/asym), y=((ai.mb-ai.orapp.cross)*100))) + 
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.05), color = "#E7B800", alpha = 0.8) +
  xlab("Duration of asymptomatic infection, days") +
  ylab("Difference in incidencerelative to \nMB approach, per 100 child-years") +
  theme_bw() +
  theme(legend.position=c(0.9,0.2)) -> af.asym

ggplot(data=subset(Merge_subset, beta.noro == 0.5), aes(x=factor(1/wan), y=((ai.mb-ai.orapp.cross)*100))) + 
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.05), color = "#E7B800", alpha = 0.8) +
  xlab("Duration of immunity, days") +
  ylab("Difference in incidencerelative to \nMB approach, per 100 child-years") +
  theme_bw() +
  theme(legend.position=c(0.9,0.2)) -> af.wan

s1.s2.asym.wan <- ggarrange(af.s1, af.s2, af.asym, af.wan,
                            nrow = 2, ncol = 2,  labels = c("A", "B", "C", "D"))

s1.s2.asym.wan

ggsave(s1.s2.asym.wan, file='paramsens.pdf', width = 12, height = 10, units = "in",  dpi = 300)


# SUPPLEMENTARY FIGURE - CHANGING S1 AND S2

ggplot(data=subset(Merge_subset, Merge_subset$s1 == 0.7 & 
                     Merge_subset$s2 == 0.3 &
                     Merge_subset$wan == 1/(365*5.1) &
                     Merge_subset$asym == 1/15)) + 
  
  geom_point(aes(x=ai.age, y=ai.orapp.cross, color = "Odds ratio approach")) +
  geom_point(aes(x=ai.age, y=ai.detect.cross, color = "Detection as etiology approach")) +
  geom_point(aes(x=ai.age, y=ai.mb, color = "Model-based approach")) +
  geom_line(aes(x=ai.age, y=ai.orapp.cross, color = "Odds ratio approach"), linetype = "dotted") +
  geom_line(aes(x=ai.age, y=ai.detect.cross, color = "Detection as etiology approach"), linetype = "dotted") +
  geom_line(aes(x=ai.age, y=ai.mb, color = "Model-based approach"), linetype = "dotted") +
  scale_color_manual(values = c( "#1B9E77", "#D95F02", "#7570B3")) +
  scale_linetype_manual(values=c("dotted", "dotted", "dotted")) +
  ylab("Norovirus attributable incidence,\n per child-year") +
  xlab("Overall AGE incidence,\n per child-year") +
  ylim(c(0,4.5)) +
  theme_bw() +
  ggtitle("s1 = 0.7, s2 = 0.3") +
  theme(legend.position=c(0.35,0.85),
        legend.title = element_blank()) -> s1.high

ggplot(data=subset(Merge_subset, Merge_subset$s1 == 0.5 & 
                     Merge_subset$s2 == 0.3 &
                     Merge_subset$wan == 1/(365*5.1) &
                     Merge_subset$asym == 1/15)) + 
  geom_point(aes(x=ai.age, y=ai.orapp.cross, color = "Odds ratio approach")) +
  geom_point(aes(x=ai.age, y=ai.detect.cross, color = "Detection as etiology approach")) +
  geom_point(aes(x=ai.age, y=ai.mb, color = "Model-based approach")) +
  geom_line(aes(x=ai.age, y=ai.orapp.cross, color = "Odds ratio approach"), linetype = "dotted") +
  geom_line(aes(x=ai.age, y=ai.detect.cross, color = "Detection as etiology approach"), linetype = "dotted") +
  geom_line(aes(x=ai.age, y=ai.mb, color = "Model-based approach"), linetype = "dotted") +
  scale_color_manual(values = c( "#1B9E77", "#D95F02", "#7570B3")) +
  scale_linetype_manual(values=c("dotted", "dotted", "dotted")) +
  ylab("Norovirus attributable incidence,\n per child-year") +
  xlab("Overall AGE incidence,\n per child-year") +
  ylim(c(0,4.5)) +
  theme_bw() +
  ggtitle("s1 = 0.5, s2 = 0.3") +
  theme(legend.position="none",
        legend.title = element_blank()) -> s1.low

ggplot(data=subset(Merge_subset, Merge_subset$s1 == 0.7 & 
                     Merge_subset$s2 == 0.5 &
                     Merge_subset$wan == 1/(365*5.1) &
                     Merge_subset$asym == 1/15)) + 
  geom_point(aes(x=ai.age, y=ai.orapp.cross, color = "Odds ratio approach")) +
  geom_point(aes(x=ai.age, y=ai.detect.cross, color = "Detection as etiology approach")) +
  geom_point(aes(x=ai.age, y=ai.mb, color = "Model-based approach")) +
  geom_line(aes(x=ai.age, y=ai.orapp.cross, color = "Odds ratio approach"), linetype = "dotted") +
  geom_line(aes(x=ai.age, y=ai.detect.cross, color = "Detection as etiology approach"), linetype = "dotted") +
  geom_line(aes(x=ai.age, y=ai.mb, color = "Model-based approach"), linetype = "dotted") +
  scale_color_manual(values = c( "#1B9E77", "#D95F02", "#7570B3")) +
  scale_linetype_manual(values=c("dotted", "dotted", "dotted")) +
  ylab("Norovirus attributable incidence,\n per child-year") +
  xlab("Overall AGE incidence,\n per child-year") +
  ylim(c(0,4.5)) +
  theme_bw() +
  ggtitle("s1 = 0.7, s2 = 0.5") +
  theme(legend.position="none",
        legend.title = element_blank()) -> s2.low


s1.s2.changingbeta <- ggarrange(s1.high, s1.low, s2.low,
                                nrow = 2, ncol = 2, labels = c("A", "B", "C"))
s1.s2.changingbeta

ggsave(s1.s2.changingbeta, file='changes1s2.pdf', width = 10, height = 10, units = "in",  dpi = 300)

# SUPPLEMENTARY FIGURE - IMPACT OF ASYMPTOMATIC INFETION DURATION (15D, 60D)

ggplot(data=subset(Merge, Merge$s1 == 0.5 & 
                     Merge$s2 == 0.3 &
                     Merge$wan == 1/(365*5.1) &
                     Merge$asym == 1/15)) + 
  geom_point(aes(x=ai.age, y=ai.orapp.cross, color = "Odds ratio approach")) +
  geom_point(aes(x=ai.age, y=ai.detect.cross, color = "Detection as etiology approach")) +
  geom_point(aes(x=ai.age, y=ai.mb, color = "Model-based approach")) +
  geom_line(aes(x=ai.age, y=ai.orapp.cross, color = "Odds ratio approach"), linetype = "dotted") +
  geom_line(aes(x=ai.age, y=ai.detect.cross, color = "Detection as etiology approach"), linetype = "dotted") +
  geom_line(aes(x=ai.age, y=ai.mb, color = "Model-based approach"), linetype = "dotted") +
  scale_color_manual(values = c( "#1B9E77", "#D95F02", "#7570B3")) +
  scale_linetype_manual(values=c("dotted", "dotted", "dotted")) +
  ylab("Norovirus attributable fraction") +
  xlab("Overall AGE incidence,\n per child-year") +
  ylim(c(0,6)) +
  theme_bw() +
  ggtitle("15 days") +
  theme(#legend.position=c(0.4,0.9),
    legend.position="none",
    legend.title = element_blank()) -> af.15

ggplot(data=subset(Merge, Merge$s1 == 0.5 & 
                     Merge$s2 == 0.3 &
                     Merge$wan == 1/(365*5.1) &
                     Merge$asym == 1/60)) + 
  geom_point(aes(x=ai.age, y=ai.orapp.cross, color = "Odds ratio approach")) +
  geom_point(aes(x=ai.age, y=ai.detect.cross, color = "Detection as etiology approach")) +
  geom_point(aes(x=ai.age, y=ai.mb, color = "Model-based approach")) +
  geom_line(aes(x=ai.age, y=ai.orapp.cross, color = "Odds ratio approach"), linetype = "dotted") +
  geom_line(aes(x=ai.age, y=ai.detect.cross, color = "Detection as etiology approach"), linetype = "dotted") +
  geom_line(aes(x=ai.age, y=ai.mb, color = "Model-based approach"), linetype = "dotted") +
  scale_color_manual(values = c( "#1B9E77", "#D95F02", "#7570B3")) +
  scale_linetype_manual(values=c("dotted", "dotted", "dotted")) +
  ylab("Norovirus attributable fraction") +
  xlab("Overall AGE incidence,\n per child-year") +
  ylim(c(0,6)) +
  theme_bw() +
  ggtitle("60 days") +
  theme(legend.position=c(0.67,0.8),
    #legend.position="none",
    legend.title = element_blank()) -> af.60

af.asym.comp <- ggarrange(af.15, af.60,
                        nrow = 1, ncol = 2, labels = c("A", "B"))
af.asym.comp

ggsave(af.asym.comp, file='asymdurcomp.pdf', width = 8, height = 5, units = "in",  dpi = 300)
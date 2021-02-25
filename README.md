# noro-attributable-fraction

Improving norovirus burden estimates using mathematical models
================

------------------------------------------------------------------------

This repository includes source code for our analysis of three different approaches to attributing gastroenteritis to norovirus, the detection as etiology (DE) approach, the odds ratio (OR) approach, and our own model-based (MB) approach. We describe a novel method using mathematical models to calculate norovirus-attributed gastroenteritis incidence.

The models described are run in R.

![Comparison of three attribution methods](https://github.com/kbratnelson/noro-attributable-fraction/blob/main/fig2.png)

------------------------------------------------------------------------

Abstract
--------

Background: Burden estimates for norovirus vary widely depending on the method used. More generally, better methods are needed to understand the contribution of different pathogens to overall enteric disease burden.

Methods: We developed a transmission model-based framework to estimate the incidence of gastroenteritis cases attributable to norovirus. We show that two common approaches for estimating disease attributions measures – using odds ratios from multi-pathogen case-control studies (OR approach) and assuming the detection of a pathogen indicates etiology (DE approach) – result in biased norovirus incidence estimates.

Results: Compared to our model-based approach, the DE approach underestimates norovirus incidence in settings with low overall AGE burden and overestimates incidence in settings with high AGE burden. In contrast, the OR approach underestimates norovirus incidence regardless of the background force of infection. The bias in the OR approach is primarily due to coinfection and asymptomatic reinfection in settings with a high burden of AGE. When we applied our model-based approach to case-control data from the MAL-ED study, we estimate norovirus incidence in six LMIC that is higher than initially reported in the MAL-ED study

Conclusions: Existing methods for estimating the incidence of enteric disease attributable to specific pathogens are biased. Our new approach to estimating norovirus incidence improves upon existing methods and can be applied to other enteric pathogens using existing data. 

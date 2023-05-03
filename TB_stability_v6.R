################# Changes in non-brome stability through time ############


#for Mac
setwd("~/Box Sync/R work/TB_stability")
search()

#for PC
setwd("C:/Users/mdtrimas/Box Sync/R work/TB_stability")
search()


#loading packages
install.packages("vegan")
install.packages("tidyverse")
install.packages("codyn")
install.packages("ggplot2")
install.packages("reshape2")
install.packages("plotly")
install.packages("nmle")
install.packages("lme4")
install.packages("olsrr")
install.packages("car")
install.packages("patchwork")
install.packages("lmerTest")
install.packages("piecewiseSEM")
install.packages("ggrepel")
install.packages("lavaan")
install.packages("semPlot")
install.packages("MuMIn")

library(vegan)
library(tidyverse)
library(codyn)
library(ggplot2)
library(reshape2)
library(plotly)
library(nlme)
library(lme4)
library(olsrr)
library(car)
library(patchwork)
library(lmerTest)
library(piecewiseSEM)
library(ggrepel)
library(lavaan)
library(semPlot)
library(MuMIn)


##################################################
##### Read in Data ##############

#sp comp data
comp2019 <- read.csv("TBobs_speciescomp2019.csv", header = TRUE, na.strings = " ") %>%
  mutate(plot = as.factor(plot))
comp2020 <- read.csv("TB2020_speciescomp.csv", header = TRUE, na.strings = " ") %>%
  mutate(plot = as.factor(plot))
comp2021 <- read.csv("TB_2021_speciescomp.csv", header = TRUE, na.strings = " ") %>%
  mutate(plot = as.factor(plot))

#abiotic data
abiotic2019 <- read.csv("TBobs_abioticdata2019.csv", header = TRUE) %>%
  mutate(plot = as.factor(plot))
abiotic2020 <- read.csv("TB2020_abioticdata.csv", header = TRUE) %>%
  mutate(plot = as.factor(plot))
abiotic2021 <- read.csv("TB_2021_abiotics.csv", header = TRUE) %>%
  mutate(plot = as.factor(plot))

#metadata
meta <- read.csv("TB_observationalplots.csv", header = TRUE) %>%
  mutate(plot = as.factor(plot))

#species info
species <- read.csv("speciesinfo_TB.csv", header = TRUE)

#var_rat
var_rat <- read.csv("var_rat.csv", header = TRUE)

#rainfall
precip <- read.csv("precip_data.csv", header = TRUE)
precip2 <- read.csv("precip_data_v2.csv", header = TRUE)

#Set ggplot2 theme to black and white
theme_set(theme_bw())
#Update ggplot2 theme - make box around the x-axis title size 30, vertically justify x-axis title to 0.35, 
#Place a margin of 15 around the x-axis title.  
#Make the x-axis title size 30. For y-axis title, make the box size 30, put the writing at a 90 degree angle, and vertically justify the title to 0.5.  
#Add a margin of 15 and make the y-axis text size 25. Make the plot title size 30 and vertically justify it to 2.  Do not add any grid lines.  
#Do not add a legend title, and make the legend size 20
theme_update(axis.title.x=element_text(size=20, vjust=-0.35, margin=margin(t=12)),
             axis.text.x=element_text(size=20), axis.title.y=element_text(size=20, angle=90, vjust=0.5,
                                                                          margin=margin(r=15)), axis.text.y=element_text(size=20), plot.title =
               element_text(size=20, vjust=2), panel.grid.major=element_blank(),
             panel.grid.minor=element_blank(),
             legend.text=element_text(size=15))



#set colorblind friendly color palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")


#########################################################################
########### Creating functions #########################

#repeated measures anova with 2 independent variables
anova_t3 <- function(IndVars=IndVars, DepVar=DepVar, RndForm=RndForm, Data=Data){
  anova_out <- {}
  IndVarMatrix <- matrix(nrow=length(IndVars),ncol=length(IndVars))
  IndVars2x <- c(IndVars,IndVars)
  
  for(REORDER in 1:length(IndVars)){
    IndVarMatrix[REORDER,] <- IndVars2x[REORDER:(length(IndVars)+(REORDER-1))]
  }
  rm(IndVars2x)
  
  for(RUN in 1:length(IndVars)){
    model_formula_temp <- paste0(DepVar,"~", paste0(IndVarMatrix[RUN,], collapse="*"))
    model_temp <- lme(as.formula(model_formula_temp)
                      , data=Data
                      , random = as.formula(RndForm)
                      , correlation=corCompSymm(form = as.formula(RndForm))
                      , control=lmeControl(returnObject=TRUE)
                      , na.action = na.omit)
    anova_out_temp <- anova(model_temp)
    
    if(length(IndVars)==2){ ## Pulls model output for variable that is last
      if(RUN==1){ ### This currently works for 
        anova_partial_temp <- rbind(anova_out_temp[(length(IndVars)+1),],
                                    anova_out_temp[4,]
        )
      }
      if(RUN %in% 2:length(IndVars)){
        anova_partial_temp <- anova_out_temp[(length(IndVars)+1),]
        
      }
    } # End if vars == 2 statement
    anova_out <- rbind(anova_out, anova_partial_temp)
  } # End reorder loop
  return(anova_out)
} # End function



#repeated measures anova with 3 independent variables
anova_t3_3ind <- function(IndVars=IndVars, DepVar=DepVar, RndForm=RndForm, Data=Data){
  anova_out <- {}
  IndVarMatrix <- matrix(nrow=length(IndVars),ncol=length(IndVars))
  IndVars2x <- c(IndVars,IndVars)
  
  for(REORDER in 1:length(IndVars)){
    IndVarMatrix[REORDER,] <- IndVars2x[REORDER:(length(IndVars)+(REORDER-1))]
  }
  rm(IndVars2x)
  
  for(RUN in 1:length(IndVars)){
    model_formula_temp <- paste0(DepVar,"~", paste0(IndVarMatrix[RUN,], collapse="*"))
    model_temp <- lme(as.formula(model_formula_temp)
                      , data=Data
                      , random = as.formula(RndForm)
                      , correlation=corCompSymm(form = as.formula(RndForm))
                      , control=lmeControl(returnObject=TRUE)
                      , na.action = na.omit)
    anova_out_temp <- anova(model_temp)
    
    if(length(IndVars)==3){ ## Pulls model output for variable that is last
      if(RUN==1){ ### This currently works for 
        anova_partial_temp <- rbind(anova_out_temp[(length(IndVars)+1),],
                                    anova_out_temp[7:8,]
        )
      }
      if(RUN %in% 2:length(IndVars)){
        anova_partial_temp <- rbind(anova_out_temp[(length(IndVars)+1),],
                                    anova_out_temp[7,]
        )
      }
    } # End if vars==3 statement
    
    if(length(IndVars)==2){ ## Pulls model output for variable that is last
      if(RUN==1){ ### This currently works for 
        anova_partial_temp <- rbind(anova_out_temp[(length(IndVars)+1),],
                                    anova_out_temp[4,]
        )
      }
      if(RUN %in% 2:length(IndVars)){
        anova_partial_temp <- anova_out_temp[(length(IndVars)+1),]
        
      }
    } # End if vars == 2 statement
    anova_out <- rbind(anova_out, anova_partial_temp)
  } # End reorder loop
  return(anova_out)
} # End function




###################################################################################
############ Cleaning sp comp data #############
comp2019_clean <- comp2019 %>%
  subset(aerial_basal == "aerial", select = -c(location, date, aerial_basal, est_tot_cov, field_added_tot_cov, excel_added_tot_cov, field_other_tot_cov, excel_other_tot_cov, bare, dung, lichen, litter, moss, rock))

#wide to long
comp2019_long <- comp2019_clean %>%
  group_by(gradient_number, plot) %>%
  gather(symbol, cover, 4:ncol(comp2019_clean))


comp2020_clean <- comp2020 %>%
  select(-c(location, date, aerial_basal, est_total, added_total, added_total_excel, cov_other_total, cov_other_total_excel, bare, dung, lichen, litter, moss, rock))
comp2020_clean[is.na(comp2020_clean)] <- 0 

#wide to long
comp2020_long <- comp2020_clean %>%
  group_by(gradient_number, plot) %>%
  gather(symbol, cover, 4:ncol(comp2020_clean))


comp2021_clean <- comp2021 %>%
  select(-c(location, date, aerial_basal, estimated_total, added_total, rock, dung, moss, lichen, mushroom, litter, bareground, final_total))
comp2021_clean[is.na(comp2021_clean)] <- 0  

#wide to long
comp2021_long <- comp2021_clean %>%
  group_by(gradient_number, plot) %>%
  gather(symbol, cover, 4:ncol(comp2021_clean))

#full join together
comp_join <- comp2019_long %>%
  full_join(comp2020_long) %>%
  full_join(comp2021_long)

#add in relative cover for each species
totcov <- comp_join %>%
  group_by(year, plot) %>%
  summarise(totcov = sum(cover)) %>%
  ungroup()

coverclean <- full_join(comp_join, totcov) %>%
  mutate(rel_cov = (cover/totcov) * 100) %>%
  full_join(meta)


#################################################################################################
################# Changes in bromes and non-bromes across gradients over time ################
## invasion gradients worked ##

#graph year vs % cov of BRAR in BRAR gradients and BRTE in BRTE gradients
BRAR_cov <- coverclean %>%
  filter(invasive_type == "BRAR" & symbol == "BRAR") 
BRAR_cov2 <- BRAR_cov %>%
  group_by(year, invasion_percent) %>%
  summarise(avg_cov = mean(cover), se_cov = sd(cover)/sqrt(length(cover)), avg_tot = mean(totcov), se_tot = sd(totcov)/sqrt(length(totcov)), avg_rel = mean(rel_cov), se_rel = sd(rel_cov)/sqrt(length(rel_cov))) %>%
  ungroup()


BRTE_cov <- coverclean %>%
  filter(invasive_type == "BRTE" & symbol == "BRTE")
BRTE_cov2 <- BRTE_cov %>%
  group_by(year, invasion_percent) %>%
  summarise(avg_cov = mean(cover), se_cov = sd(cover)/sqrt(length(cover)), avg_tot = mean(totcov), se_tot = sd(totcov)/sqrt(length(totcov)), avg_rel = mean(rel_cov), se_rel = sd(rel_cov)/sqrt(length(rel_cov))) %>%
  ungroup()


BRARcov3 <- ggplot(data = BRAR_cov2, aes(x = invasion_percent, y = avg_rel, color = factor(year), shape = factor(year))) + 
  geom_point(size = 3) + 
  ylim(-3, 65) + 
  geom_smooth(aes(group = as.factor(year)), method = "lm", se = FALSE) + 
  scale_colour_manual(values = cbPalette) +
  geom_errorbar(aes(ymin = avg_rel - se_rel, ymax = avg_rel + se_rel), width = 7, linewidth = 0.75, position = position_dodge(.05)) +
  labs(x = "", y = "Relative % Cover", title = "BRAR", color = "Year", shape = "Year") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), legend.title = element_text(size = 16)) + 
  theme(panel.border = element_rect(linewidth = 1))

BRTEcov3 <- ggplot(data = BRTE_cov2, aes(x = invasion_percent, y = avg_rel, color = factor(year), shape = factor(year))) + 
  geom_point(size = 3) + 
  ylim(-3, 80) + 
  geom_smooth(aes(group = as.factor(year)), method = "lm", se = FALSE) + 
  scale_colour_manual(values = cbPalette) +
  geom_errorbar(aes(ymin = avg_rel - se_rel, ymax = avg_rel + se_rel), width = 7, linewidth = 0.75, position = position_dodge(.05)) +
  labs(x = "Invasion %", y = "Relative % Cover", title = "BRTE", color = "Year", shape = "Year") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), legend.title = element_text(size = 16)) + 
  theme(panel.border = element_rect(linewidth = 1))

BRARcov3 / BRTEcov3




#stats
#untransformed
mod1 <- anova_t3(IndVars = c("year", "invasion_percent"), 
         DepVar = "rel_cov", 
         RndForm = "~1 | gradient_number", 
         Data = BRAR_cov) #inv% p < 0.0001

mod2<- anova_t3(IndVars = c("year", "invasion_percent"), 
         DepVar = "rel_cov", 
         RndForm = "~1 | gradient_number", 
         Data = BRTE_cov) #inv% p < 0.0001

#check normality of residuals
resid1 <- lm(data = BRAR_cov, rel_cov ~ invasion_percent*year)
ols_plot_resid_hist(resid1) #right skew
ols_test_normality(resid1) #failed tests

resid2 <- lm(data = BRTE_cov, rel_cov ~ invasion_percent*year)
ols_plot_resid_hist(resid2) #fairly normal
ols_test_normality(resid2) #passed SW and KS tests


#try transformations - only need BRAR
BRAR_cov_trans <- BRAR_cov %>%
  mutate(ln_cov = log(rel_cov + 0.1))

resid3 <- lm(data = BRAR_cov_trans, ln_cov ~ invasion_percent*year)
ols_plot_resid_hist(resid3) #more normal
ols_test_normality(resid3) #passed tests except Cramer

#stats on transformed BRAR
mod3 <- anova_t3(IndVars = c("year", "invasion_percent"), 
         DepVar = "ln_cov", 
         RndForm = "~1 | gradient_number", 
         Data = BRAR_cov_trans) #inv% p < 0.0001; year p = 0.0219


#check linearity
plot(mod3)


## real inv % vs stability of bromes
#average across years to get real inv% 
BRARcover <- BRAR_cov %>%
  group_by(plot) %>%
  summarise(avgBRAR = mean(rel_cov)) %>%
  ungroup()

BRTEcover <- BRTE_cov %>%
  group_by(plot) %>%
  summarise(avgBRTE = mean(rel_cov)) %>%
  ungroup()

BRARcovstab <- community_stability(df = BRAR_cov, 
                                   time.var = "year", 
                                   abundance.var = "rel_cov",     
                                   replicate.var = "plot") %>%
  full_join(BRARcover) %>%
  full_join(meta) %>%
  filter(invasive_type == "BRAR")

BRTEcovstab <- community_stability(df = BRTE_cov, 
                                   time.var = "year", 
                                   abundance.var = "rel_cov",     
                                   replicate.var = "plot") %>%
  full_join(BRTEcover) %>%
  full_join(meta) %>%
  filter(invasive_type == "BRTE")





#########################################################################


####################################################################################################################
############ Stability metrics without bromes - evenness, richness, synch, turnover, tot cover, bareground #######
#no BRAR in BRAR gradients, no BRTE in BRTE gradients


####### Stability of Richness and evenness ######
BRARcomp_join_nobrome <- comp_join %>%
  filter(symbol != "BRAR")

BRTEcomp_join_nobrome <- comp_join %>%
  filter(symbol != "BRTE")

#calculate rich and even
BRARcomm_struct_nobrome <- community_structure(df = BRARcomp_join_nobrome, 
                                   time.var = "year", 
                                   abundance.var = "cover", 
                                   replicate.var = "plot", 
                                   metric = "Evar")  %>%
  full_join(meta) %>%
  filter(invasive_type == "BRAR")

BRTEcomm_struct_nobrome <- community_structure(df = BRTEcomp_join_nobrome, 
                                               time.var = "year", 
                                               abundance.var = "cover", 
                                               replicate.var = "plot", 
                                               metric = "Evar")  %>%
  full_join(meta) %>%
  filter(invasive_type == "BRTE")

#join with real brome cov
BRARcov_only <- coverclean %>% 
  filter(symbol == "BRAR") %>%
  mutate(BRARcov = rel_cov) %>%
  dplyr::select(c(plot, invasive_type, invasion_percent, year, BRARcov))

BRARcov_rich_nobrome <- BRARcov_only %>%
  full_join(BRARcomm_struct_nobrome) %>%
  filter(invasive_type == "BRAR")


BRTEcov_only <- coverclean %>% 
  filter(symbol == "BRTE") %>%
  mutate(BRTEcov = rel_cov) %>%
  dplyr::select(c(plot, invasive_type, invasion_percent, year, BRTEcov))

BRTEcov_rich_nobrome <- BRTEcov_only %>%
  full_join(BRTEcomm_struct_nobrome) %>%
  filter(invasive_type == "BRTE")


#richness stability
BRARrichstab_nobrome <- community_stability(df = BRARcov_rich_nobrome, 
                                    time.var = "year", 
                                    abundance.var = "richness",     #stability of richness
                                    replicate.var = "plot") %>%
  full_join(meta) %>%
  filter(invasive_type == "BRAR")


BRTErichstab_nobrome <- community_stability(df = BRTEcov_rich_nobrome, 
                                    time.var = "year", 
                                    abundance.var = "richness",     #stability of richness
                                    replicate.var = "plot") %>%
  full_join(meta) %>%
  filter(invasive_type == "BRTE")

#evenness stability
BRARevenstab_nobrome <- community_stability(df = BRARcov_rich_nobrome, 
                                    time.var = "year", 
                                    abundance.var = "Evar",     #stability of evenness
                                    replicate.var = "plot") %>%
  full_join(meta) %>%
  filter(invasive_type == "BRAR")


BRTEevenstab_nobrome <- community_stability(df = BRTEcov_rich_nobrome, 
                                    time.var = "year", 
                                    abundance.var = "Evar",     #stability of evenness
                                    replicate.var = "plot") %>%
  full_join(meta) %>%
  filter(invasive_type == "BRTE")

#average across years to get real inv% 
BRARcov_only2 <- BRARcov_only %>%
  group_by(plot) %>%
  summarise(avgBRAR = mean(BRARcov)) %>%
  ungroup()

BRARrichstab2_nobrome <- BRARrichstab_nobrome %>%
  full_join(BRARcov_only2) %>%
  drop_na(invasive_type)

BRARevenstab2_nobrome <- BRARevenstab_nobrome %>%
  full_join(BRARcov_only2) %>%
  drop_na(invasive_type)


BRTEcov_only2 <- BRTEcov_only %>%
  group_by(plot) %>%
  summarise(avgBRTE = mean(BRTEcov)) %>%
  ungroup()

BRTErichstab2_nobrome <- BRTErichstab_nobrome %>%
  full_join(BRTEcov_only2) %>%
  drop_na(invasive_type) %>%
  filter(stability != "Inf")

BRTEevenstab2_nobrome <- BRTEevenstab_nobrome %>%
  full_join(BRTEcov_only2) %>%
  drop_na(invasive_type)




#stats
#richness
BRARrichstab_nobrome_stats <- lmerTest::lmer(data = BRARrichstab2_nobrome, stability ~ avgBRAR +
                                       (1|gradient_number))
anova(BRARrichstab_nobrome_stats, type = 3) #p = 0.2598


BRTErichstab_nobrome_stats <- lmerTest::lmer(data = BRTErichstab2_nobrome, stability ~ avgBRTE +
                                       (1|gradient_number))
anova(BRTErichstab_nobrome_stats, type = 3) #p = 0.0153
rsquared(BRTErichstab_nobrome_stats) #r = 0.2508846

#evenness
BRARevenstab_nobrome_stats <- lmerTest::lmer(data = BRARevenstab2_nobrome, stability ~ avgBRAR +
                                       (1|gradient_number))
anova(BRARevenstab_nobrome_stats, type = 3) #p = 0.9911

BRTEevenstab_nobrome_stats <- lmerTest::lmer(data = BRTEevenstab2_nobrome, stability ~ avgBRTE +
                                       (1|gradient_number))
anova(BRTEevenstab_nobrome_stats, type = 3) #p = 0.2649


#check normality
#BRAR richness
resid1.0 <- lm(data = BRARrichstab2_nobrome, stability ~ avgBRAR)
ols_plot_resid_hist(resid1.0) #right skew sort of 
ols_test_normality(resid1.0) #passed KS


#try transformation - ln looks good
BRARrichstab2_nobrome_trans <- BRARrichstab2_nobrome %>%
  mutate(sqrt_stab = sqrt(stability)) %>%
  mutate(ln_stab = log(stability + 0.1)) %>%
  mutate(sq_stab = stability^2) %>%
  mutate(cb_stab = stability^1/3) %>%
  mutate(qd_stab = stability^1/4) %>%
  mutate(stab3 = stability^3) %>%
  mutate(stab4 = stability^4)

resid1.1 <- lm(data = BRARrichstab2_nobrome_trans, ln_stab ~ avgBRAR)
ols_plot_resid_hist(resid1.1) #better 
ols_test_normality(resid1.1) #passed SW, KS, AD

#retry stats
BRARrichstab_nobrome_stats2 <- lmerTest::lmer(data = BRARrichstab2_nobrome_trans, ln_stab ~ avgBRAR +
                                        (1|gradient_number))
anova(BRARrichstab_nobrome_stats2, type = 3) #p = 0.1562


#BRTE richness - stick with untransformed
resid1.2 <- lm(data = BRTErichstab2_nobrome, stability ~ avgBRTE)
ols_plot_resid_hist(resid1.2) #normal looking 
ols_test_normality(resid1.2) #passed KS, AD


#BRAR evenness
resid1.3 <- lm(data = BRARevenstab2_nobrome, stability ~ avgBRAR)
ols_plot_resid_hist(resid1.3) #right skew 
ols_test_normality(resid1.3) #failed

#try transformation - ln looks better
BRARevenstab2_nobrome_trans <- BRARevenstab2_nobrome %>%
  mutate(sqrt_stab = sqrt(stability)) %>%
  mutate(ln_stab = log(stability + 0.1)) %>%
  mutate(sq_stab = stability^2) %>%
  mutate(cb_stab = stability^1/3) %>%
  mutate(qd_stab = stability^1/4) %>%
  mutate(stab3 = stability^3) %>%
  mutate(stab4 = stability^4)

resid1.4 <- lm(data = BRARevenstab2_nobrome_trans, ln_stab ~ avgBRAR)
ols_plot_resid_hist(resid1.4) #right skew but better 
ols_test_normality(resid1.4) #passed KS, SW, AD

#retry stats
BRARevenstab_nobrome_stats2 <- lmerTest::lmer(data = BRARevenstab2_nobrome_trans, ln_stab ~ avgBRAR +
                                        (1|gradient_number))
anova(BRARevenstab_nobrome_stats2, type = 3) #p = 0.9335


#BRTE evenness - log trans
resid1.5 <- lm(data = BRTEevenstab2_nobrome, stability ~ avgBRTE)
ols_plot_resid_hist(resid1.5) #normalish
ols_test_normality(resid1.5) #passes SW, KS, AD

#try transformation - ln looks better
BRTEevenstab2_nobrome_trans <- BRTEevenstab2_nobrome %>%
  mutate(sqrt_stab = sqrt(stability)) %>%
  mutate(ln_stab = log(stability + 0.1)) %>%
  mutate(sq_stab = stability^2) %>%
  mutate(cb_stab = stability^1/3) %>%
  mutate(qd_stab = stability^1/4) %>%
  mutate(stab3 = stability^3) %>%
  mutate(stab4 = stability^4)

resid1.7 <- lm(data = BRTEevenstab2_nobrome_trans, ln_stab ~ avgBRTE)
ols_plot_resid_hist(resid1.7) #normal
ols_test_normality(resid1.7) #pass

#retry stats
BRTEevenstab_nobrome_stats2 <- lmerTest::lmer(data = BRTEevenstab2_nobrome_trans, ln_stab ~ avgBRTE +
                                        (1|gradient_number))
anova(BRTEevenstab_nobrome_stats2, type = 3) #p = 0.06556
rsquared(BRTEevenstab_nobrome_stats2) #r = 0.1347841


#checking other assumptions of linear models
#1)linear relationship, 2)normality, 3) no multicolinearity, 4) no auto-correlation, 5) homoscedasticity
#3 is not applicable bc we only ever have 1 predictor

#richness
#BRAR
#linearity
plot(resid(BRARrichstab_nobrome_stats2), BRARrichstab2_nobrome_trans$ln_stab) #looks linear
plot(BRARrichstab_nobrome_stats2) #no pattern so indicates linearity

#homoscedascity
BRARrichstab2_nobrome_trans$res <- residuals(BRARrichstab_nobrome_stats2)
BRARrichstab2_nobrome_trans$abs_res <- abs(BRARrichstab2_nobrome_trans$res)
BRARrichstab2_nobrome_trans$abs_res2 <- BRARrichstab2_nobrome_trans$abs_res^2
levene_brar_rich <- lm(data = BRARrichstab2_nobrome_trans, abs_res2 ~ avgBRAR)
anova(levene_brar_rich) #p = 0.4752 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(BRARrichstab_nobrome_stats2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(BRARrichstab_nobrome_stats2, retype = "normalized"))

#BRTE
#linearity
plot(resid(BRTErichstab_nobrome_stats), BRTErichstab2_nobrome$stability) #looks linear

#homoscedascity
BRTErichstab2_nobrome$res <- residuals(BRTErichstab_nobrome_stats)
BRTErichstab2_nobrome$abs_res <- abs(BRTErichstab2_nobrome$res)
BRTErichstab2_nobrome$abs_res2 <- BRTErichstab2_nobrome$abs_res^2
levene_brte_rich <- lm(data = BRTErichstab2_nobrome, abs_res2 ~ avgBRTE)
anova(levene_brte_rich) #p = 0.0559 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(BRTErichstab_nobrome_stats, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(BRTErichstab_nobrome_stats, retype = "normalized"))


#evenness
#BRAR
#linearity
plot(resid(BRARevenstab_nobrome_stats2), BRARevenstab2_nobrome_trans$ln_stab) #looks linear

#homoscedascity
BRARevenstab2_nobrome_trans$res <- residuals(BRARevenstab_nobrome_stats2)
BRARevenstab2_nobrome_trans$abs_res <- abs(BRARevenstab2_nobrome_trans$res)
BRARevenstab2_nobrome_trans$abs_res2 <- BRARevenstab2_nobrome_trans$abs_res^2
levene_brar_even <- lm(data = BRARevenstab2_nobrome_trans, abs_res2 ~ avgBRAR)
anova(levene_brar_even) #p = 0.8678 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(BRARevenstab_nobrome_stats2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(BRARevenstab_nobrome_stats2, retype = "normalized"))

#BRTE
#linearity
plot(resid(BRTEevenstab_nobrome_stats2), BRTEevenstab2_nobrome_trans$ln_stab) #looks linear

#homoscedascity
BRTEevenstab2_nobrome_trans$res <- residuals(BRTEevenstab_nobrome_stats2)
BRTEevenstab2_nobrome_trans$abs_res <- abs(BRTEevenstab2_nobrome_trans$res)
BRTEevenstab2_nobrome_trans$abs_res2 <- BRTEevenstab2_nobrome_trans$abs_res^2
levene_brte_even <- lm(data = BRTEevenstab2_nobrome_trans, abs_res2 ~ avgBRTE)
anova(levene_brte_even) #p = 0.1524 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(BRTEevenstab_nobrome_stats2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(BRTEevenstab_nobrome_stats2, retype = "normalized"))



##############################################

###############################################
####### Synchrony ######
#bromes not included
comp_met <- comp_join %>%
  full_join(meta)
BRARcompmet <- comp_met %>%
  filter(invasive_type == "BRAR")
BRTEcompmet <- comp_met %>%
  filter(invasive_type == "BRTE")
BRARcompmet2 <- BRARcompmet %>%
  filter(symbol != "BRAR")
BRTEcompmet2 <- BRTEcompmet %>%
  filter(symbol != "BRTE")


BRARsynch_nobrome <- synchrony(df = BRARcompmet2, 
                        time.var = "year", 
                        species.var = "symbol", 
                        abundance.var = "cover", 
                        replicate.var = "plot") %>%
  mutate(plot = as.factor(plot)) %>%
  full_join(meta) %>%
  filter(invasive_type == "BRAR") %>%
  rename(synch = synchrony)

BRTEsynch_nobrome <- synchrony(df = BRTEcompmet2, 
                        time.var = "year", 
                        species.var = "symbol", 
                        abundance.var = "cover", 
                        replicate.var = "plot") %>%
  mutate(plot = as.factor(plot)) %>%
  full_join(meta) %>%
  filter(invasive_type == "BRTE") %>%
  rename(synch = synchrony)

#real rel % vs synch
#bromes not included
BRARsynch_nobrome2 <- BRARsynch_nobrome %>%
  full_join(BRARcov_only2)

BRTEsynch_nobrome2 <- BRTEsynch_nobrome %>%
  full_join(BRTEcov_only2)



#stats rel cov vs synch
BRARsync_nobrome_stats <- lmerTest::lmer(data = BRARsynch_nobrome2, synch ~ avgBRAR +
                                    (1|gradient_number))
anova(BRARsync_nobrome_stats, type = 3) #p = 0.6921

BRTEsync_nobrome_stats <- lmerTest::lmer(data = BRTEsynch_nobrome2, synch ~ avgBRTE +
                                    (1|gradient_number))
anova(BRTEsync_nobrome_stats, type = 3) #p = 0.4246


#check on normality
res3 <- lm(data = BRARsynch_nobrome2, synch ~ avgBRAR)
ols_plot_resid_hist(res3) #normal
ols_test_normality(res3) #pass tests

res4 <- lm(data = BRTEsynch_nobrome2, synch ~ avgBRTE)
ols_plot_resid_hist(res4) #normal
ols_test_normality(res4) #pass some


#other assumptions
#synchrony
#BRAR
#linearity
BRARsynch_nobrome2 <- BRARsynch_nobrome2 %>%
  drop_na(invasion_percent)
plot(resid(BRARsync_nobrome_stats), BRARsynch_nobrome2$synch) #looks linear

#homoscedascity
BRARsynch_nobrome2$res <- residuals(BRARsync_nobrome_stats)
BRARsynch_nobrome2$abs_res <- abs(BRARsynch_nobrome2$res)
BRARsynch_nobrome2$abs_res2 <- BRARsynch_nobrome2$abs_res^2
levene_brar_synch <- lm(data = BRARsynch_nobrome2, abs_res2 ~ avgBRAR)
anova(levene_brar_synch) #p = 0.49 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(BRARsync_nobrome_stats, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(BRARsync_nobrome_stats, retype = "normalized"))

#BRTE
#linearity
BRTEsynch_nobrome2 <- BRTEsynch_nobrome2 %>%
  drop_na(invasion_percent)
plot(resid(BRTEsync_nobrome_stats), BRTEsynch_nobrome2$synch) #looks linear

#homoscedascity
BRTEsynch_nobrome2$res <- residuals(BRTEsync_nobrome_stats)
BRTEsynch_nobrome2$abs_res <- abs(BRTEsynch_nobrome2$res)
BRTEsynch_nobrome2$abs_res2 <- BRTEsynch_nobrome2$abs_res^2
levene_brte_synch <- lm(data = BRTEsynch_nobrome2, abs_res2 ~ avgBRTE)
anova(levene_brte_synch) #p = 0.7818 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(BRTEsync_nobrome_stats, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(BRTEsync_nobrome_stats, retype = "normalized"))






####### Turnover ######
#community turnover - start to end
noBRAR_covclean <- coverclean %>%
  filter(symbol != "BRAR" & invasive_type == "BRAR")
noBRTE_covclean <- coverclean %>%
  filter(symbol != "BRTE" & invasive_type == "BRTE")

noBRAR_turn <- turnover(df = subset(noBRAR_covclean, year != "2020"), 
                        time.var = "year",
                        species.var = "symbol",
                        abundance.var = "cover", 
                        replicate.var = "plot") %>%
  mutate(plot = as.factor(plot)) %>%
  full_join(meta) %>%
  filter(invasive_type == "BRAR")


noBRTE_turn <- turnover(df = subset(noBRTE_covclean, year != "2020"), 
                        time.var = "year",
                        species.var = "symbol",
                        abundance.var = "cover", 
                        replicate.var = "plot") %>%
  mutate(plot = as.factor(plot)) %>%
  full_join(meta) %>%
  filter(invasive_type == "BRTE")


#real rel % vs turnover
BRARcov_only2 <- BRARcov_only2 %>%
  mutate(plot = as.factor(plot))
noBRAR_turn_cov <- BRARcov_only2 %>%
  full_join(noBRAR_turn)

BRTEcov_only2 <- BRTEcov_only2 %>%
  mutate(plot = as.factor(plot))
noBRTE_turn_cov <- BRTEcov_only2 %>%
  full_join(noBRTE_turn)



#stats rel cov vs turnover
noBRAR_turn_cov_stats <- lmerTest::lmer(data = noBRAR_turn_cov, total ~ avgBRAR +
                                          (1|gradient_number))
anova(noBRAR_turn_cov_stats, type = 3) #p = 0.01998

rsquared(noBRAR_turn_cov_stats) #marg = 0.1991717; cond = 0.3561574

noBRTE_turn_cov_stats <- lmerTest::lmer(data = noBRTE_turn_cov, total ~ avgBRTE +
                                          (1|gradient_number))
anova(noBRTE_turn_cov_stats, type = 3) #p = 0.0006951

rsquared(noBRTE_turn_cov_stats) #marg = 0.3107569; cond = 0.5590348


#check on normality
resid72 <- lm(data = noBRAR_turn_cov, total ~ avgBRAR)
ols_plot_resid_hist(resid72) #normal
ols_test_normality(resid72) #pass tests

resid73 <- lm(data = noBRTE_turn_cov, total ~ avgBRTE)
ols_plot_resid_hist(resid73) #normal
ols_test_normality(resid73) #pass some tests (SW, KS)


#other assumptions
#turnover
#BRAR
#linearity
noBRAR_turn_cov <- noBRAR_turn_cov %>%
  drop_na(invasion_percent)
plot(resid(noBRAR_turn_cov_stats), noBRAR_turn_cov$total) #looks linear

#homoscedascity
noBRAR_turn_cov$res <- residuals(noBRAR_turn_cov_stats)
noBRAR_turn_cov$abs_res <- abs(noBRAR_turn_cov$res)
noBRAR_turn_cov$abs_res2 <- noBRAR_turn_cov$abs_res^2
levene_brar_turn <- lm(data = noBRAR_turn_cov, abs_res2 ~ avgBRAR)
anova(levene_brar_turn) #p = 0.106 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(noBRAR_turn_cov_stats, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(noBRAR_turn_cov_stats, retype = "normalized"))

#BRTE
#linearity
noBRTE_turn_cov <- noBRTE_turn_cov %>%
  drop_na(invasion_percent)
plot(resid(noBRTE_turn_cov_stats), noBRTE_turn_cov$total) #looks linear

#homoscedascity
noBRTE_turn_cov$res <- residuals(noBRTE_turn_cov_stats)
noBRTE_turn_cov$abs_res <- abs(noBRTE_turn_cov$res)
noBRTE_turn_cov$abs_res2 <- noBRTE_turn_cov$abs_res^2
levene_brte_turn <- lm(data = noBRTE_turn_cov, abs_res2 ~ avgBRTE)
anova(levene_brte_turn) #p = 0.4495 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(noBRTE_turn_cov_stats, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(noBRTE_turn_cov_stats, retype = "normalized"))





####### Stability of Total Cover ######
#community stability
noBRAR_covclean <- coverclean %>%
  filter(symbol != "BRAR" & invasive_type == "BRAR")

noBRAR_stab <- community_stability(df = noBRAR_covclean, 
                                   time.var = "year", 
                                   abundance.var = "cover", 
                                   replicate.var = "plot") %>%
  full_join(meta) %>%
  filter(invasive_type == "BRAR")


noBRTE_covclean <- coverclean %>%
  filter(symbol != "BRTE" & invasive_type == "BRTE")

noBRTE_stab <- community_stability(df = noBRTE_covclean, 
                                   time.var = "year", 
                                   abundance.var = "cover", 
                                   replicate.var = "plot") %>%
  full_join(meta) %>%
  filter(invasive_type == "BRTE")


#stability with real rel % of bromes
#BRAR in BRAR gradietns and BRTE in BRTE gradients only
#need to average cover across years

BRARcov_only2 <- BRARcov_only %>%
  group_by(gradient_number, plot, invasive_type, invasion_percent) %>%
  summarise(avgBRAR = mean(BRARcov)) %>%
  ungroup() %>%
  filter(invasive_type == "BRAR")

BRARcov_only3 <- BRARcov_only2 %>%
  full_join(noBRAR_stab)

BRTEcov_only2 <- BRTEcov_only %>%
  group_by(gradient_number, plot, invasive_type, invasion_percent) %>%
  summarise(avgBRTE = mean(BRTEcov)) %>%
  ungroup() %>%
  filter(invasive_type == "BRTE")

BRTEcov_only3 <- BRTEcov_only2 %>%
  full_join(noBRTE_stab)


#stats rel cov vs stability
noBRAR_stab_stats3 <- lmerTest::lmer(data = BRARcov_only3, stability ~ avgBRAR +
                                       (1|gradient_number))
anova(noBRAR_stab_stats3, type = 3) #not sig


noBRTE_stab_stats3 <- lmerTest::lmer(data = BRTEcov_only3, stability ~ avgBRTE +
                                       (1|gradient_number))
anova(noBRTE_stab_stats3, type = 3) #not sig


#check on normality
resid59 <- lm(data = BRARcov_only3, stability ~ avgBRAR)
ols_plot_resid_hist(resid59) #right skew
ols_test_normality(resid59) #fail tests

resid60 <- lm(data = BRTEcov_only3, stability ~ avgBRTE)
ols_plot_resid_hist(resid60) #right skew
ols_test_normality(resid60) #fail tests

#try transformation 
#BRAR - ln looks closest to normal
BRARcov_only3_trans <- BRARcov_only3 %>%
  mutate(sqrt_stab = sqrt(stability)) %>%
  mutate(ln_stab = log(stability + 0.1)) %>%
  mutate(sq_stab = stability^2) %>%
  mutate(cb_stab = stability^1/3) %>%
  mutate(qd_stab = stability^1/4) %>%
  mutate(stab3 = stability^3) %>%
  mutate(stab4 = stability^4)

resid61 <- lm(data = BRARcov_only3_trans, sqrt_stab ~ avgBRAR)
ols_plot_resid_hist(resid61) #right skew
ols_test_normality(resid61) #fail tests

resid62 <- lm(data = BRARcov_only3_trans, ln_stab ~ avgBRAR)
ols_plot_resid_hist(resid62) #right skew
ols_test_normality(resid62) #fail tests except KS

resid63 <- lm(data = BRARcov_only3_trans, sq_stab ~ avgBRAR)
ols_plot_resid_hist(resid63) #right skew
ols_test_normality(resid63) #fail tests

resid64 <- lm(data = BRARcov_only3_trans, cb_stab ~ avgBRAR)
ols_plot_resid_hist(resid64) #right skew
ols_test_normality(resid64) #fail tests

resid65 <- lm(data = BRARcov_only3_trans, qd_stab ~ avgBRAR)
ols_plot_resid_hist(resid65) #right skew
ols_test_normality(resid65) #fail tests

resid66 <- lm(data = BRARcov_only3_trans, stab3 ~ avgBRAR)
ols_plot_resid_hist(resid66) #right skew
ols_test_normality(resid66) #fail tests

resid67 <- lm(data = BRARcov_only3_trans, stab4 ~ avgBRAR)
ols_plot_resid_hist(resid67) #right skew
ols_test_normality(resid67) #fail tests

#BRTE - ln looks closest to normal
BRTEcov_only3_trans <- BRTEcov_only3 %>%
  mutate(sqrt_stab = sqrt(stability)) %>%
  mutate(ln_stab = log(stability + 0.1)) %>%
  mutate(sq_stab = stability^2) %>%
  mutate(cb_stab = stability^1/3) %>%
  mutate(qd_stab = stability^1/4) %>%
  mutate(stab3 = stability^3) %>%
  mutate(stab4 = stability^4)

resid68 <- lm(data = BRTEcov_only3_trans, sqrt_stab ~ avgBRTE)
ols_plot_resid_hist(resid68) #looks closer to normal
ols_test_normality(resid68) #fail tests

resid69 <- lm(data = BRTEcov_only3_trans, ln_stab ~ avgBRTE)
ols_plot_resid_hist(resid69) #looks closer to normal
ols_test_normality(resid69) #fail tests except KS and marginal SW

#try stats again with ln transformations
noBRAR_stab_stats4 <- lmerTest::lmer(data = BRARcov_only3_trans, ln_stab ~ avgBRAR +
                                       (1|gradient_number))
anova(noBRAR_stab_stats4, type = 3) #p = 0.3395

noBRTE_stab_stats4 <- lmerTest::lmer(data = BRTEcov_only3_trans, ln_stab ~ avgBRTE +
                                       (1|gradient_number))
anova(noBRTE_stab_stats4, type = 3) #p = 0.7263


#other assumptions
#cover stability
#BRAR
#linearity
plot(resid(noBRAR_stab_stats4), BRARcov_only3_trans$ln_stab) #looks linear

#homoscedascity
BRARcov_only3_trans$res <- residuals(noBRAR_stab_stats4)
BRARcov_only3_trans$abs_res <- abs(BRARcov_only3_trans$res)
BRARcov_only3_trans$abs_res2 <- BRARcov_only3_trans$abs_res^2
levene_brar_stab <- lm(data = BRARcov_only3_trans, abs_res2 ~ avgBRAR)
anova(levene_brar_stab) #p = 0.6571 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(noBRAR_stab_stats4, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(noBRAR_stab_stats4, retype = "normalized"))

#BRTE
#linearity
plot(resid(noBRTE_stab_stats4), BRTEcov_only3_trans$ln_stab) #looks linear

#homoscedascity
BRTEcov_only3_trans$res <- residuals(noBRTE_stab_stats4)
BRTEcov_only3_trans$abs_res <- abs(BRTEcov_only3_trans$res)
BRTEcov_only3_trans$abs_res2 <- BRTEcov_only3_trans$abs_res^2
levene_brte_stab <- lm(data = BRTEcov_only3_trans, abs_res2 ~ avgBRTE)
anova(levene_brte_stab) #p = 0.07563 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(noBRTE_stab_stats4, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(noBRTE_stab_stats4, retype = "normalized"))






####### Stability of Bareground Cover ######
#just want bareground cover
comp2019_clean2 <- comp2019 %>%
  subset(aerial_basal == "aerial", select = c(year, gradient_number, plot, bare, litter))

comp2020_clean2 <- comp2020 %>%
  select(c(year, gradient_number, plot, bare, litter)) 

comp2021_clean2 <- comp2021 %>%
  select(c(year, gradient_number, plot, bareground, litter)) %>%
  mutate(bare = bareground) %>%
  select(-bareground)


#full join together
comp_bare <- comp2019_clean2 %>%
  full_join(comp2020_clean2) %>%
  full_join(comp2021_clean2) %>%
  full_join(meta)

BRARcomp_bare <- comp_bare %>%
  filter(invasive_type == "BRAR")

BRTEcomp_bare <- comp_bare %>%
  filter(invasive_type == "BRTE")


##bare stability
BRAR_barestab <- community_stability(df = BRARcomp_bare, 
                                     time.var = "year", 
                                     abundance.var = "bare", 
                                     replicate.var = "plot") %>%
  full_join(BRARcov_only2) 

BRTE_barestab <- community_stability(df = BRTEcomp_bare, 
                                     time.var = "year", 
                                     abundance.var = "bare", 
                                     replicate.var = "plot") %>%
  full_join(BRTEcov_only2) 



#stats
BRAR_barestab_stats <- lmerTest::lmer(data = BRAR_barestab, stability ~ avgBRAR +
                                        (1|gradient_number))
anova(BRAR_barestab_stats, type = 3) #p = 0.4846


BRTE_barestab_stats <- lmerTest::lmer(data = BRTE_barestab, stability ~ avgBRTE +
                                        (1|gradient_number))
anova(BRTE_barestab_stats, type = 3) #p = 0.04223
rsquared(BRTE_barestab_stats) #marg = 0.1080827; cond = 0.4765484


#check normality
#BRAR - use log
resid1.18 <- lm(data = BRAR_barestab, stability ~ avgBRAR)
ols_plot_resid_hist(resid1.18) #right skew
ols_test_normality(resid1.18) #failed

#try trans
BRAR_barestab_trans <- BRAR_barestab %>%
  mutate(sqrt_stab = sqrt(stability)) %>%
  mutate(ln_stab = log10(stability)) %>%
  mutate(sq_stab = stability^2) %>%
  mutate(cb_stab = stability^1/3) %>%
  mutate(qd_stab = stability^1/4) %>%
  mutate(stab3 = stability^3) %>%
  mutate(stab4 = stability^4)

resid1.19 <- lm(data = BRAR_barestab_trans, ln_stab ~ avgBRAR)
ols_plot_resid_hist(resid1.19) #normal
ols_test_normality(resid1.19) #passed SW, KS, AD

#redo stats with ln trans
BRAR_barestab_stats2 <- lmerTest::lmer(data = BRAR_barestab_trans, ln_stab ~ avgBRAR +
                                         (1|gradient_number))
anova(BRAR_barestab_stats2, type = 3) #p = 0.4181

#BRTE - use untrans
resid1.20 <- lm(data = BRTE_barestab, stability ~ avgBRTE)
ols_plot_resid_hist(resid1.20) #normal
ols_test_normality(resid1.20) #passed SW, KS


#other assumptions
#bareground stability
#BRAR
#linearity
plot(resid(BRAR_barestab_stats2), BRAR_barestab_trans$ln_stab) #looks linear

#homoscedascity
BRAR_barestab_trans$res <- residuals(BRAR_barestab_stats2)
BRAR_barestab_trans$abs_res <- abs(BRAR_barestab_trans$res)
BRAR_barestab_trans$abs_res2 <- BRAR_barestab_trans$abs_res^2
levene_brar_bare <- lm(data = BRAR_barestab_trans, abs_res2 ~ avgBRAR)
anova(levene_brar_bare) #p = 0.5297 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(BRAR_barestab_stats2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(BRAR_barestab_stats2, retype = "normalized"))

#BRTE
#linearity
plot(resid(BRTE_barestab_stats), BRTE_barestab$stability) #looks linear

#homoscedascity
BRTE_barestab$res <- residuals(BRTE_barestab_stats)
BRTE_barestab$abs_res <- abs(BRTE_barestab$res)
BRTE_barestab$abs_res2 <- BRTE_barestab$abs_res^2
levene_brte_bare <- lm(data = BRTE_barestab, abs_res2 ~ avgBRTE)
anova(levene_brte_bare) #p = 0.1129 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(BRTE_barestab_stats, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(BRTE_barestab_stats, retype = "normalized"))




###############################################################




#####################################################################################################
############# Make combined fig of rich, even, synch, turn, tot cov, bare without bromes ##########
#richness
BRARrichstab3_nobrome <- ggplot(data = BRARrichstab2_nobrome, aes(x = avgBRAR, y = stability)) + 
  geom_point(size = 3) + 
  ylim(0, 30) +
  labs(x = "", y = "Richness Stability", title = "BRAR") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), axis.text.x = element_blank()) + 
  theme(panel.border = element_rect(linewidth = 2))

BRTErichstab3_nobrome <- ggplot(data = BRTErichstab2_nobrome, aes(x = avgBRTE, y = stability)) + 
  geom_point(size = 3) + 
  ylim(0, 30) +
  geom_smooth(method = "lm", se = FALSE) + 
  labs(x = "", y = "", title = "BRTE") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), axis.text.x = element_blank()) + 
  theme(panel.border = element_rect(linewidth = 2))

#evenness
BRARevenstab3_nobrome <- ggplot(data = BRARevenstab2_nobrome, aes(x = avgBRAR, y = stability)) + 
  geom_point(size = 3) + 
  ylim(0, 15) +
  labs(x = "", y = "Evenness Stability", title = "BRAR") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), axis.text.x = element_blank()) + 
  theme(panel.border = element_rect(linewidth = 2))

BRTEevenstab3_nobrome <- ggplot(data = BRTEevenstab2_nobrome, aes(x = avgBRTE, y = stability)) + 
  geom_point(size = 3) + 
  ylim(0, 20) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") + 
  labs(x = "", y = "", title = "BRTE") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), axis.text.x = element_blank()) + 
  theme(panel.border = element_rect(linewidth = 2))

#synchrony
BRARsynch_fig_nobrome <- ggplot(data = BRARsynch_nobrome2, aes(x = avgBRAR, y = synch)) + 
  ylim(0, 0.7) + 
  geom_point(size = 3) + 
  labs(x = "", y = "Synchrony", title = "BRAR") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), axis.text.x = element_blank()) + 
  theme(panel.border = element_rect(linewidth = 2))

BRTEsynch_fig_nobrome <- ggplot(data = BRTEsynch_nobrome2, aes(x = avgBRTE, y = synch)) + 
  ylim(0, 0.7) +
  geom_point(size = 3) + 
  labs(x = "", y = "", title = "BRTE") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), axis.text.x = element_blank()) + 
  theme(panel.border = element_rect(linewidth = 2))

#turnover
noBRAR_turn_cov2 <- ggplot(data = noBRAR_turn_cov, aes(x = avgBRAR, y = total)) + 
  ylim(0, 1.2) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE) + 
  labs(x = "", y = "Species Turnover", title = "BRAR") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), axis.text.x = element_blank()) + 
  theme(panel.border = element_rect(linewidth = 2))

noBRTE_turn_cov2 <- ggplot(data = noBRTE_turn_cov, aes(x = avgBRTE, y = total)) + 
  ylim(0, 1.2) +
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE) + 
  labs(x = "", y = "", title = "BRTE") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), axis.text.x = element_blank()) + 
  theme(panel.border = element_rect(linewidth = 2))

#total cover
noBRAR_stab3 <- ggplot(data = BRARcov_only3, aes(x = avgBRAR, y = stability)) + 
  ylim(0, 30) + 
  geom_point(size = 3) + 
  labs(x = "", y = "Cover Stability", title = "BRAR") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), axis.text.x = element_blank()) + 
  theme(panel.border = element_rect(linewidth = 2))

noBRTE_stab3 <- ggplot(data = BRTEcov_only3, aes(x = avgBRTE, y = stability)) + 
  ylim(0, 20) +
  geom_point(size = 3) + 
  labs(x = "", y = "", title = "BRTE") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), axis.text.x = element_blank()) + 
  theme(panel.border = element_rect(linewidth = 2))

#bareground
BRAR_barestab2 <- ggplot(data = BRAR_barestab, aes(x = avgBRAR, y = stability)) + 
  geom_point(size = 3) + 
  ylim(0, 25) +
  labs(x = "Relative % Cover", y = "Bareground Stability", title = "BRAR") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) + 
  theme(panel.border = element_rect(linewidth = 2))

BRTE_barestab2 <- ggplot(data = BRTE_barestab, aes(x = avgBRTE, y = stability)) + 
  geom_point(size = 3) + 
  ylim(0, 10) +
  geom_smooth(method = "lm", se = FALSE) + 
  labs(x = "Relative % Cover", y = "", title = "BRTE") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) + 
  theme(panel.border = element_rect(linewidth = 2))

BRARrichstab3_nobrome + BRTErichstab3_nobrome + BRARevenstab3_nobrome + BRTEevenstab3_nobrome + BRARsynch_fig_nobrome + BRTEsynch_fig_nobrome + 
noBRAR_turn_cov2 + noBRTE_turn_cov2 + noBRAR_stab3 + noBRTE_stab3 + BRAR_barestab2 + BRTE_barestab2 + plot_layout(ncol = 2)

#2fig   1000x2200

##########################################################################




###########################################################################
################# Functional group stability C4, C3, forbs #########
#BRAR
BRARcov_only2_grp <- BRARcov_only2 %>%
  group_by(invasion_percent, invasive_type) %>%
  summarise(avg_BRAR = mean(avgBRAR)) %>%
  ungroup()

BRARcompmet_grp <- BRARcompmet %>%
  full_join(BRARcov_only2_grp)

BRARsp <- species %>%
  full_join(BRARcompmet_grp) %>%
  filter(spec_funct_group2 != "unknown")

#BRTE
BRTEcov_only2_grp <- BRTEcov_only2 %>%
  group_by(invasion_percent, invasive_type) %>%
  summarise(avg_BRTE = mean(avgBRTE)) %>%
  ungroup()

BRTEcompmet_grp <- BRTEcompmet %>%
  full_join(BRTEcov_only2_grp)

BRTEsp <- species %>%
  full_join(BRTEcompmet_grp) %>%
  filter(spec_funct_group2 != "unknown")


#C4 grass
#BRAR
BRARsp_c4 <- BRARsp %>%
  filter(spec_funct_group4 == "C4 grass") %>%
  filter(cover != 0)

BRARsp_c4_stab <- community_stability(df = BRARsp_c4, 
                                      time.var = "year", 
                                      abundance.var = "cover", 
                                      replicate.var = "plot") %>%
  full_join(BRARcov_only2)

BRARsp_c4_stab_stats <- lmerTest::lmer(data = BRARsp_c4_stab, stability ~ avgBRAR +
                                         (1|gradient_number))

anova(BRARsp_c4_stab_stats, type = 3) #p = 0.03181

res_a <- lm(data = BRARsp_c4_stab, stability ~ avgBRAR)
ols_plot_resid_hist(res_a) #right skew
ols_test_normality(res_a) #fail

BRARsp_c4_stab_trans <- BRARsp_c4_stab %>%
  mutate(sqrt_stab = sqrt(stability)) %>%
  mutate(ln_stab = log10(stability))

res_b <- lm(data = BRARsp_c4_stab_trans, ln_stab ~ avgBRAR)
ols_plot_resid_hist(res_b) #normal
ols_test_normality(res_b) #passed KS

BRARsp_c4_stab_stats2 <- lmerTest::lmer(data = BRARsp_c4_stab_trans, ln_stab ~ avgBRAR +
                                          (1|gradient_number))
anova(BRARsp_c4_stab_stats2, type = 3) #p = 0.02457
rsquared(BRARsp_c4_stab_stats2) #marg = 0.2121843; cond = 0.274917

#BRTE
BRTEsp_c4 <- BRTEsp %>%
  filter(spec_funct_group4 == "C4 grass") %>%
  filter(cover != 0)

BRTEsp_c4_stab <- community_stability(df = BRTEsp_c4, 
                                      time.var = "year", 
                                      abundance.var = "cover", 
                                      replicate.var = "plot") %>%
  full_join(BRTEcov_only2) 
BRTEsp_c4_stab2 <- BRTEsp_c4_stab %>%
  filter(stability != "Inf") %>%
  filter(stability != "NA")

BRTEsp_c4_stab_stats <- lmerTest::lmer(data = BRTEsp_c4_stab2, stability ~ avgBRTE +
                                         (1|gradient_number))

anova(BRTEsp_c4_stab_stats, type = 3) #p = 0.4211

res_a <- lm(data = BRTEsp_c4_stab2, stability ~ avgBRTE)
ols_plot_resid_hist(res_a) #right skew
ols_test_normality(res_a) #fail

BRTEsp_c4_stab_trans <- BRTEsp_c4_stab2 %>%
  mutate(sqrt_stab = sqrt(stability)) %>%
  mutate(ln_stab = log10(stability))

res_b <- lm(data = BRTEsp_c4_stab_trans, ln_stab ~ avgBRTE)
ols_plot_resid_hist(res_b) #normal
ols_test_normality(res_b) #passed KS

BRTEsp_c4_stab_stats2 <- lmerTest::lmer(data = BRTEsp_c4_stab_trans, ln_stab ~ avgBRTE +
                                          (1|gradient_number))
anova(BRTEsp_c4_stab_stats2, type = 3) #p = 0.8744


#other assumptions
#C4 stability
#BRAR
#linearity
BRARsp_c4_stab_trans <- BRARsp_c4_stab_trans %>%
  drop_na(stability)
plot(resid(BRARsp_c4_stab_stats2), BRARsp_c4_stab_trans$ln_stab) #looks linear
plot(BRARsp_c4_stab_stats2) #random

#homoscedascity
BRARsp_c4_stab_trans$res <- residuals(BRARsp_c4_stab_stats2)
BRARsp_c4_stab_trans$abs_res <- abs(BRARsp_c4_stab_trans$res)
BRARsp_c4_stab_trans$abs_res2 <- BRARsp_c4_stab_trans$abs_res^2
levene_brar_c4 <- lm(data = BRARsp_c4_stab_trans, abs_res2 ~ avgBRAR)
anova(levene_brar_c4) #p = 0.1103 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(BRARsp_c4_stab_stats2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(BRARsp_c4_stab_stats2, retype = "normalized"))

#BRTE
#linearity
BRTEsp_c4_stab_trans <- BRTEsp_c4_stab_trans %>%
  drop_na(stability)
plot(resid(BRTEsp_c4_stab_stats2), BRTEsp_c4_stab_trans$ln_stab) #looks linear
plot(BRTEsp_c4_stab_stats2) #random

#homoscedascity
BRTEsp_c4_stab_trans$res <- residuals(BRTEsp_c4_stab_stats2)
BRTEsp_c4_stab_trans$abs_res <- abs(BRTEsp_c4_stab_trans$res)
BRTEsp_c4_stab_trans$abs_res2 <- BRTEsp_c4_stab_trans$abs_res^2
levene_brte_c4 <- lm(data = BRTEsp_c4_stab_trans, abs_res2 ~ avgBRTE)
anova(levene_brte_c4) #p = 0.08314 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(BRTEsp_c4_stab_stats2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(BRTEsp_c4_stab_stats2, retype = "normalized"))




#C3 grass - non-brome
#BRAR
BRARsp_c3_nobrome <- BRARsp %>%
  filter(symbol != "BRAR") %>%
  filter(spec_funct_group4 == "C3 grass") 

BRARsp_c3_nobrome_stab <- community_stability(df = BRARsp_c3_nobrome, 
                                      time.var = "year", 
                                      abundance.var = "cover", 
                                      replicate.var = "plot") %>%
  full_join(BRARcov_only2)

BRARsp_c3_nobrome_stab_stats <- lmerTest::lmer(data = BRARsp_c3_nobrome_stab, stability ~ avgBRAR +
                                         (1|gradient_number))

anova(BRARsp_c3_nobrome_stab_stats, type = 3) #p = 0.2318

res_a <- lm(data = BRARsp_c3_nobrome_stab, stability ~ avgBRAR)
ols_plot_resid_hist(res_a) #right skew
ols_test_normality(res_a) #fail

BRARsp_c3_stab_nobrome_trans <- BRARsp_c3_nobrome_stab %>%
  mutate(sqrt_stab = sqrt(stability)) %>%
  mutate(ln_stab = log10(stability))

res_b <- lm(data = BRARsp_c3_stab_nobrome_trans, ln_stab ~ avgBRAR)
ols_plot_resid_hist(res_b) #normal
ols_test_normality(res_b) #passed KS

BRARsp_c3_stab_nobrome_stats2 <- lmerTest::lmer(data = BRARsp_c3_stab_nobrome_trans, ln_stab ~ avgBRAR +
                                          (1|gradient_number))
anova(BRARsp_c3_stab_nobrome_stats2, type = 3) #p = 0.2657


#BRTE
BRTEsp_c3_nobrome <- BRTEsp %>%
  filter(symbol != "BRTE") %>%
  filter(spec_funct_group4 == "C3 grass") 

BRTEsp_c3_nobrome_stab <- community_stability(df = BRTEsp_c3_nobrome, 
                                              time.var = "year", 
                                              abundance.var = "cover", 
                                              replicate.var = "plot") %>%
  full_join(BRTEcov_only2)

BRTEsp_c3_nobrome_stab_stats <- lmerTest::lmer(data = BRTEsp_c3_nobrome_stab, stability ~ avgBRTE +
                                                 (1|gradient_number))

anova(BRTEsp_c3_nobrome_stab_stats, type = 3) #p = 0.001387
rsquared(BRTEsp_c3_nobrome_stab_stats) #marg = 0.2784599, cond = 0.5386053

res_a <- lm(data = BRTEsp_c3_nobrome_stab, stability ~ avgBRTE)
ols_plot_resid_hist(res_a) #pass
ols_test_normality(res_a) #none needed


#other assumptions
#C3 stability - no bromes
#BRAR
#linearity
plot(resid(BRARsp_c3_stab_nobrome_stats2), BRARsp_c3_stab_nobrome_trans$ln_stab) #looks linear
plot(BRARsp_c3_stab_nobrome_stats2) #random

#homoscedascity
BRARsp_c3_stab_nobrome_trans$res <- residuals(BRARsp_c3_stab_nobrome_stats2)
BRARsp_c3_stab_nobrome_trans$abs_res <- abs(BRARsp_c3_stab_nobrome_trans$res)
BRARsp_c3_stab_nobrome_trans$abs_res2 <- BRARsp_c3_stab_nobrome_trans$abs_res^2
levene_brar_c3 <- lm(data = BRARsp_c3_stab_nobrome_trans, abs_res2 ~ avgBRAR)
anova(levene_brar_c3) #p = 0.1936 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(BRARsp_c3_stab_nobrome_stats2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(BRARsp_c3_stab_nobrome_stats2, retype = "normalized"))

#BRTE
#linearity
plot(resid(BRTEsp_c3_nobrome_stab_stats), BRTEsp_c3_nobrome_stab$stability) #looks linear
plot(BRTEsp_c3_nobrome_stab_stats) #random

#homoscedascity
BRARsp_c3_stab_nobrome_trans$res <- residuals(BRTEsp_c3_nobrome_stab_stats)
BRARsp_c3_stab_nobrome_trans$abs_res <- abs(BRARsp_c3_stab_nobrome_trans$res)
BRARsp_c3_stab_nobrome_trans$abs_res2 <- BRARsp_c3_stab_nobrome_trans$abs_res^2
levene_brar_c3 <- lm(data = BRARsp_c3_stab_nobrome_trans, abs_res2 ~ avgBRAR)
anova(levene_brar_c3) #p = 0.2516 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(BRTEsp_c3_nobrome_stab_stats, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(BRTEsp_c3_nobrome_stab_stats, retype = "normalized"))




#Forbs
#BRAR
BRARsp_f <- BRARsp %>%
  filter(spec_funct_group4 == "forb") 

BRARsp_f_stab <- community_stability(df = BRARsp_f, 
                                      time.var = "year", 
                                      abundance.var = "cover", 
                                      replicate.var = "plot") %>%
  full_join(BRARcov_only2)

BRARsp_f_stab_stats <- lmerTest::lmer(data = BRARsp_f_stab, stability ~ avgBRAR +
                                         (1|gradient_number))

anova(BRARsp_f_stab_stats, type = 3) #p = 0.009093

res_a <- lm(data = BRARsp_f_stab, stability ~ avgBRAR)
ols_plot_resid_hist(res_a) #right skew
ols_test_normality(res_a) #fail

BRARsp_f_stab_trans <- BRARsp_f_stab %>%
  mutate(sqrt_stab = sqrt(stability)) %>%
  mutate(ln_stab = log10(stability))

res_b <- lm(data = BRARsp_f_stab_trans, ln_stab ~ avgBRAR)
ols_plot_resid_hist(res_b) #normal
ols_test_normality(res_b) #passed KS

BRARsp_f_stab_stats2 <- lmerTest::lmer(data = BRARsp_f_stab_trans, ln_stab ~ avgBRAR +
                                          (1|gradient_number))
anova(BRARsp_f_stab_stats2, type = 3) #p = 0.01687
rsquared(BRARsp_f_stab_stats2) #marg = 0.1678035; cond = 0.5377462


#BRTE
BRTEsp_f <- BRTEsp %>%
  filter(spec_funct_group4 == "forb") 

BRTEsp_f_stab <- community_stability(df = BRTEsp_f, 
                                     time.var = "year", 
                                     abundance.var = "cover", 
                                     replicate.var = "plot") %>%
  full_join(BRTEcov_only2)
BRTEsp_f_stab2 <- BRTEsp_f_stab %>%
  filter(stability != "Inf") %>%
  filter(stability != "NA")

BRTEsp_f_stab_stats <- lmerTest::lmer(data = BRTEsp_f_stab2, stability ~ avgBRTE +
                                        (1|gradient_number))

anova(BRTEsp_f_stab_stats, type = 3) #p = 0.1738

res_a <- lm(data = BRTEsp_f_stab2, stability ~ avgBRTE)
ols_plot_resid_hist(res_a) #right skew
ols_test_normality(res_a) #fail

BRTEsp_f_stab_trans <- BRTEsp_f_stab2 %>%
  mutate(sqrt_stab = sqrt(stability)) %>%
  mutate(ln_stab = log10(stability))

res_b <- lm(data = BRTEsp_f_stab_trans, ln_stab ~ avgBRTE)
ols_plot_resid_hist(res_b) #normal
ols_test_normality(res_b) #passed KS

BRTEsp_f_stab_stats2 <- lmerTest::lmer(data = BRTEsp_f_stab_trans, ln_stab ~ avgBRTE +
                                         (1|gradient_number))
anova(BRTEsp_f_stab_stats2, type = 3) #p = 0.06719
rsquared(BRTEsp_f_stab_stats2) #marg = 0.03723929; cond = 0.8007314


#other assumptions
#forb stability
#BRAR
#linearity
plot(resid(BRARsp_f_stab_stats2), BRARsp_f_stab_trans$ln_stab) #looks linear
plot(BRARsp_f_stab_stats2) #random

#homoscedascity
BRARsp_f_stab_trans$res <- residuals(BRARsp_f_stab_stats2)
BRARsp_f_stab_trans$abs_res <- abs(BRARsp_f_stab_trans$res)
BRARsp_f_stab_trans$abs_res2 <- BRARsp_f_stab_trans$abs_res^2
levene_brar_f <- lm(data = BRARsp_f_stab_trans, abs_res2 ~ avgBRAR)
anova(levene_brar_f) #p = 0.07786 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(BRARsp_f_stab_stats2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(BRARsp_f_stab_stats2, retype = "normalized"))

#BRTE
#linearity
plot(resid(BRTEsp_f_stab_stats2), BRTEsp_f_stab_trans$ln_stab) #looks linear
plot(BRTEsp_f_stab_stats2) #random

#homoscedascity
BRTEsp_f_stab_trans$res <- residuals(BRTEsp_f_stab_stats2)
BRTEsp_f_stab_trans$abs_res <- abs(BRTEsp_f_stab_trans$res)
BRTEsp_f_stab_trans$abs_res2 <- BRTEsp_f_stab_trans$abs_res^2
levene_brte_f <- lm(data = BRTEsp_f_stab_trans, abs_res2 ~ avgBRTE)
anova(levene_brte_f) #p = 0.9174 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(BRTEsp_f_stab_stats2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(BRTEsp_f_stab_stats2, retype = "normalized"))




#C3 grass - with-brome
#BRAR
BRARsp_c3 <- BRARsp %>%
  filter(spec_funct_group4 == "C3 grass") 

BRARsp_c3_stab <- community_stability(df = BRARsp_c3, 
                                              time.var = "year", 
                                              abundance.var = "cover", 
                                              replicate.var = "plot") %>%
  full_join(BRARcov_only2)

BRARsp_c3_stab_stats <- lmerTest::lmer(data = BRARsp_c3_stab, stability ~ avgBRAR +
                                                 (1|gradient_number))

anova(BRARsp_c3_stab_stats, type = 3) #p = 0.3478

res_a <- lm(data = BRARsp_c3_stab, stability ~ avgBRAR)
ols_plot_resid_hist(res_a) #right skew
ols_test_normality(res_a) #fail

BRARsp_c3_stab_trans <- BRARsp_c3_stab %>%
  mutate(sqrt_stab = sqrt(stability)) %>%
  mutate(ln_stab = log10(stability))

res_b <- lm(data = BRARsp_c3_stab_trans, ln_stab ~ avgBRAR)
ols_plot_resid_hist(res_b) #normal
ols_test_normality(res_b) #passed KS

BRARsp_c3_stab_stats2 <- lmerTest::lmer(data = BRARsp_c3_stab_trans, ln_stab ~ avgBRAR +
                                                  (1|gradient_number))
anova(BRARsp_c3_stab_stats2, type = 3) #p = 0.2829


#BRTE
BRTEsp_c3 <- BRTEsp %>%
  filter(spec_funct_group4 == "C3 grass") 

BRTEsp_c3_stab <- community_stability(df = BRTEsp_c3, 
                                      time.var = "year", 
                                      abundance.var = "cover", 
                                      replicate.var = "plot") %>%
  full_join(BRTEcov_only2)

BRTEsp_c3_stab_stats <- lmerTest::lmer(data = BRTEsp_c3_stab, stability ~ avgBRTE +
                                         (1|gradient_number))

anova(BRTEsp_c3_stab_stats, type = 3) #p = 0.2553

res_a <- lm(data = BRTEsp_c3_stab, stability ~ avgBRTE)
ols_plot_resid_hist(res_a) #right skew
ols_test_normality(res_a) #fail

BRTEsp_c3_stab_trans <- BRTEsp_c3_stab %>%
  mutate(sqrt_stab = sqrt(stability)) %>%
  mutate(ln_stab = log10(stability))

res_b <- lm(data = BRTEsp_c3_stab_trans, ln_stab ~ avgBRTE)
ols_plot_resid_hist(res_b) #normaler
ols_test_normality(res_b) #passed KS

BRTEsp_c3_stab_stats2 <- lmerTest::lmer(data = BRTEsp_c3_stab_trans, ln_stab ~ avgBRTE +
                                          (1|gradient_number))
anova(BRTEsp_c3_stab_stats2, type = 3) #p = 0.1433


#other assumptions
#c3 stability - with bromes
#BRAR
#linearity
plot(resid(BRARsp_c3_stab_stats2), BRARsp_c3_stab_trans$ln_stab) #looks linear
plot(BRARsp_c3_stab_stats2) #random

#homoscedascity
BRARsp_c3_stab_trans$res <- residuals(BRARsp_c3_stab_stats2)
BRARsp_c3_stab_trans$abs_res <- abs(BRARsp_c3_stab_trans$res)
BRARsp_c3_stab_trans$abs_res2 <- BRARsp_c3_stab_trans$abs_res^2
levene_brar_c3_2 <- lm(data = BRARsp_c3_stab_trans, abs_res2 ~ avgBRAR)
anova(levene_brar_c3_2) #p = 0.3691 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(BRARsp_c3_stab_stats2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(BRARsp_c3_stab_stats2, retype = "normalized"))

#BRTE
#linearity
plot(resid(BRTEsp_c3_stab_stats2), BRTEsp_c3_stab_trans$ln_stab) #looks linear
plot(BRTEsp_c3_stab_stats2) #random

#homoscedascity
BRTEsp_c3_stab_trans$res <- residuals(BRTEsp_c3_stab_stats2)
BRTEsp_c3_stab_trans$abs_res <- abs(BRTEsp_c3_stab_trans$res)
BRTEsp_c3_stab_trans$abs_res2 <- BRTEsp_c3_stab_trans$abs_res^2
levene_brte_c3_2 <- lm(data = BRTEsp_c3_stab_trans, abs_res2 ~ avgBRTE)
anova(levene_brte_c3_2) #p = 0.3057 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(BRTEsp_c3_stab_stats2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(BRTEsp_c3_stab_stats2, retype = "normalized"))




##############################################################3


##############################################################
############ Figure of C4, C3, forbs without bromes #########
#C4 fig
BRARc4_fig <- ggplot(data = BRARsp_c4_stab, aes(x = avgBRAR, y = stability)) +
  geom_point(size = 3) + 
  ylim(0, 6) +
  geom_smooth(method = "lm", se = FALSE) + 
  labs(x = "", y = "C4 Grass Stability", title = "BRAR") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), axis.text.x = element_blank()) + 
  theme(panel.border = element_rect(linewidth = 2))

BRTEc4_fig <- ggplot(data = BRTEsp_c4_stab2, aes(x = avgBRTE, y = stability)) +
  geom_point(size = 3) + 
  ylim(0, 6) +
  labs(x = "", y = "", title = "BRTE") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), axis.text.x = element_blank()) + 
  theme(panel.border = element_rect(linewidth = 2))

#C3 fig no brome
BRARc3_nobrome_fig <- ggplot(data = BRARsp_c3_nobrome_stab, aes(x = avgBRAR, y = stability)) +
  geom_point(size = 3) + 
  ylim(0, 8) +
  labs(x = "", y = "C3 Grass Stability", title = "BRAR") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), axis.text.x = element_blank()) + 
  theme(panel.border = element_rect(linewidth = 2))

BRTEc3_nobrome_fig <- ggplot(data = BRTEsp_c3_nobrome_stab, aes(x = avgBRTE, y = stability)) +
  geom_point(size = 3) + 
  ylim(0, 6) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "", y = "", title = "BRTE") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), axis.text.x = element_blank()) + 
  theme(panel.border = element_rect(linewidth = 2))

#forb fig
BRARf_fig <- ggplot(data = BRARsp_f_stab, aes(x = avgBRAR, y = stability)) +
  geom_point(size = 3) + 
  ylim(0, 6) +
  geom_smooth(method = "lm", se = FALSE) + 
  labs(x = "Relative % Cover", y = "Forb Stability", title = "BRAR") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) + 
  theme(panel.border = element_rect(linewidth = 2))

BRTEf_fig <- ggplot(data = BRTEsp_f_stab2, aes(x = avgBRTE, y = stability)) +
  geom_point(size = 3) + 
  ylim(0, 5) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") + 
  labs(x = "Relative % Cover", y = "", title = "BRTE") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) + 
  theme(panel.border = element_rect(linewidth = 2))

BRARc4_fig + BRTEc4_fig + BRARc3_nobrome_fig + BRTEc3_nobrome_fig + BRARf_fig + BRTEf_fig + plot_layout(ncol = 2)

#3fig 800x1000


##############################################################



##############################################################
############### Stability of light and soil moist ######
abiotic2019_2 <- abiotic2019 %>%
  dplyr::select(-c(location, date, high_clay_soil_moist, humidity, air_temp, PAR_time))
abiotic2020_2 <- abiotic2020 %>%
  dplyr::select(-c(location, date, time_start, time_end))
abiotic2021_2 <- abiotic2021 %>%
  dplyr::select(-c(location, date, time_start, time_end))

abiotics <- abiotic2019_2 %>%
  full_join(abiotic2020_2) %>%
  full_join(abiotic2021_2) %>%
  mutate(plot = as.factor(plot)) %>%
  full_join(meta) %>%
  mutate(pct_trans = (PAR_below/PAR_above) * 100)

BRARavgab <- abiotics %>%
  filter(invasive_type == "BRAR")

BRTEavgab <- abiotics %>%
  filter(invasive_type == "BRTE")


##average light and soil moist
BRARavgab_join <- BRARavgab %>%
  full_join(BRARcov_only2)

BRARavgab_group <- BRARavgab_join %>%
  group_by(gradient_number, plot, avgBRAR) %>%
  summarise(avg_light = mean(pct_trans), se_light = sd(pct_trans)/sqrt(length(pct_trans)), avg_moist = mean(soil_moist), se_moist = sd(soil_moist)/sqrt(length(soil_moist))) %>%
  ungroup()

BRTEavgab_join <- BRTEavgab %>%
  full_join(BRTEcov_only2)

BRTEavgab_group <- BRTEavgab_join %>%
  group_by(gradient_number, plot, avgBRTE) %>%
  summarise(avg_light = mean(pct_trans), se_light = sd(pct_trans)/sqrt(length(pct_trans)), avg_moist = mean(soil_moist), se_moist = sd(soil_moist)/sqrt(length(soil_moist))) %>%
  ungroup()


#graphs
#light
BRAR_avglt <- ggplot(data = BRARavgab_group, aes(x = avgBRAR, y = avg_light)) + 
  geom_point(size = 3) + 
  ylim(40, 105) +
  stat_smooth(method = "lm", se = FALSE) + 
  labs(x = "BRAR Rel Cov (%)", y = "Light Transmittance (%)", title = "BRAR") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), axis.title.x = element_blank(), axis.text.x = element_blank())+ 
  theme(panel.border = element_rect(linewidth = 2))

BRTE_avglt <- ggplot(data = BRTEavgab_group, aes(x = avgBRTE, y = avg_light)) + 
  geom_point(size = 3) + 
  ylim(40, 105) +
  stat_smooth(method = "lm", se = FALSE) + 
  labs(x = "BRTE Rel Cov (%)", y = "Light Transmittance (%)", title = "BRTE") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank()) + 
  theme(panel.border = element_rect(linewidth = 2))

#soil moist
BRAR_avgsm <- ggplot(data = BRARavgab_group, aes(x = avgBRAR, y = avg_moist)) + 
  geom_point(size = 3) + 
  ylim(0, 32) +
  labs(x = "Relative % Cover", y = "Soil Moisture (VWC%)", title = "BRAR") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) + 
  theme(panel.border = element_rect(linewidth = 2))

BRTE_avgsm <- ggplot(data = BRTEavgab_group, aes(x = avgBRTE, y = avg_moist)) + 
  geom_point(size = 3) + 
  ylim(0, 15) +
  labs(x = "Relative % Cover", y = "Soil Moisture (VWC%)", title = "BRTE") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), axis.title.y = element_blank()) + 
  theme(panel.border = element_rect(linewidth = 2))


BRAR_avglt + BRTE_avglt + BRAR_avgsm + BRTE_avgsm + plot_layout(ncol = 2)


#stats
#light
BRAR_avglt_stats <- lmerTest::lmer(data = BRARavgab_group, avg_light ~ avgBRAR +
                                       (1|gradient_number))
anova(BRAR_avglt_stats, type = 3) #p = 8.842e-7
AIC(BRAR_avglt_stats) #175.4533
rsquared(BRAR_avglt_stats)  #marg = 0.6561166, cond = 0.6748827

BRTE_avglt_stats <- lmerTest::lmer(data = BRTEavgab_group, avg_light ~ avgBRTE +
                                     (1|gradient_number))
anova(BRTE_avglt_stats, type = 3) #p = 2.28e-8
AIC(BRTE_avglt_stats) #157.2403
rsquared(BRTE_avglt_stats)  #marg = 0.7548147, cond = 0.7612477


#check normality
#BRAR - use untransformed
resid1.11 <- lm(data = BRARavgab_group, avg_light ~ avgBRAR)
ols_plot_resid_hist(resid1.11) #left skew
ols_test_normality(resid1.11) #marginal

#try trans
BRAR_avglt_trans <- BRARavgab_group %>%
  mutate(ln_lt = log10(avg_light + 0.1)) 

resid1.12 <- lm(data = BRAR_avglt_trans, ln_lt ~ avgBRAR)
ols_plot_resid_hist(resid1.12) #worse
ols_test_normality(resid1.12) #failed except KS, stick with untransformed

#BRTE - use untransformed
resid1.11 <- lm(data = BRTEavgab_group, avg_light ~ avgBRTE)
ols_plot_resid_hist(resid1.11) #normal
ols_test_normality(resid1.11) #pass

BRTE_avglt_trans <- BRTEavgab_group %>%
  mutate(ln_lt = log10(avg_light + 0.1)) 

BRTE_avglt_stats2 <- lmerTest::lmer(data = BRTE_avglt_trans, ln_lt ~ avgBRTE +
                                     (1|gradient_number))
anova(BRTE_avglt_stats2, type = 3) #p = 2.89e-8
AIC(BRTE_avglt_stats2) #-80.32344
rsquared(BRTE_avglt_stats2)  #marg = 0.7362576


#soil moist
BRAR_avgsm_stats <- lmerTest::lmer(data = BRARavgab_group, avg_moist ~ avgBRAR +
                                     (1|gradient_number))
anova(BRAR_avgsm_stats, type = 3) #p = 0.1795



BRTE_avgsm_stats <- lmerTest::lmer(data = BRTEavgab_group, avg_moist ~ avgBRTE +
                                     (1|gradient_number))
anova(BRTE_avgsm_stats, type = 3) #p = 0.222

#BRAR - use untransformed
resid1.11 <- lm(data = BRARavgab_group, avg_moist ~ avgBRAR)
ols_plot_resid_hist(resid1.11) #normal
ols_test_normality(resid1.11) #pass

#BRTE - use untransformed
resid1.11 <- lm(data = BRTEavgab_group, avg_moist ~ avgBRTE)
ols_plot_resid_hist(resid1.11) #normal
ols_test_normality(resid1.11) #pass


#other assumptions
#Avg light
#BRAR
#linearity
plot(resid(BRAR_avglt_stats), BRARavgab_group$avg_light) #looks linear
plot(BRAR_avglt_stats) #random

#homoscedascity
BRARavgab_group$res <- residuals(BRAR_avglt_stats)
BRARavgab_group$abs_res <- abs(BRARavgab_group$res)
BRARavgab_group$abs_res2 <- BRARavgab_group$abs_res^2
levene_brar_light <- lm(data = BRARavgab_group, abs_res2 ~ avgBRAR)
anova(levene_brar_light) #p = 0.2316 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(BRAR_avglt_stats, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(BRAR_avglt_stats, retype = "normalized"))

#BRTE
#linearity
plot(resid(BRTE_avglt_stats), BRTEavgab_group$avg_light) #looks linear
plot(BRTE_avglt_stats) #random

#homoscedascity
BRTEavgab_group$res <- residuals(BRTE_avglt_stats)
BRTEavgab_group$abs_res <- abs(BRTEavgab_group$res)
BRTEavgab_group$abs_res2 <- BRTEavgab_group$abs_res^2
levene_brte_light <- lm(data = BRTEavgab_group, abs_res2 ~ avgBRTE)
anova(levene_brte_light) #p = 0.2316 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(BRTE_avglt_stats, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(BRTE_avglt_stats, retype = "normalized"))



#Avg soil moisture
#BRAR
#linearity
plot(resid(BRAR_avgsm_stats), BRARavgab_group$avg_moist) #looks linear
plot(BRAR_avgsm_stats) #random

#homoscedascity
BRARavgab_group$res <- residuals(BRAR_avgsm_stats)
BRARavgab_group$abs_res <- abs(BRARavgab_group$res)
BRARavgab_group$abs_res2 <- BRARavgab_group$abs_res^2
levene_brar_sm <- lm(data = BRARavgab_group, abs_res2 ~ avgBRAR)
anova(levene_brar_sm) #p = 0.149 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(BRAR_avgsm_stats, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(BRAR_avgsm_stats, retype = "normalized"))

#BRTE
#linearity
plot(resid(BRTE_avgsm_stats), BRTEavgab_group$avg_moist) #looks linear
plot(BRTE_avgsm_stats) #random

#homoscedascity
BRTEavgab_group$res <- residuals(BRTE_avgsm_stats)
BRTEavgab_group$abs_res <- abs(BRTEavgab_group$res)
BRTEavgab_group$abs_res2 <- BRTEavgab_group$abs_res^2
levene_brte_sm <- lm(data = BRTEavgab_group, abs_res2 ~ avgBRTE)
anova(levene_brte_sm) #p = 0.054 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(BRTE_avgsm_stats, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(BRTE_avgsm_stats, retype = "normalized"))






##light stability
BRAR_pctstab <- community_stability(df = BRARavgab, 
                                    time.var = "year", 
                                    abundance.var = "pct_trans", 
                                    replicate.var = "plot") %>%
  full_join(BRARcov_only2) 


BRTE_pctstab <- community_stability(df = BRTEavgab, 
                                    time.var = "year", 
                                    abundance.var = "pct_trans", 
                                    replicate.var = "plot") %>%
  full_join(BRTEcov_only2) 



#stats
BRAR_pctstab_stats <- lmerTest::lmer(data = BRAR_pctstab, stability ~ avgBRAR +
                                       (1|gradient_number))
anova(BRAR_pctstab_stats, type = 3) #p = 0.0.0009256
AIC(BRAR_pctstab_stats) #227.8387


BRTE_pctstab_stats <- lmerTest::lmer(data = BRTE_pctstab, stability ~ avgBRTE +
                                       (1|gradient_number))
anova(BRTE_pctstab_stats, type = 3) #p = 0.0.0008559
AIC(BRTE_pctstab_stats) #222.1642


#check normality
#BRAR - use ln
resid1.11 <- lm(data = BRAR_pctstab, stability ~ avgBRAR)
ols_plot_resid_hist(resid1.11) #right skew
ols_test_normality(resid1.11) #failed

#try trans
BRAR_pctstab_trans <- BRAR_pctstab %>%
  mutate(sqrt_stab = sqrt(stability)) %>%
  mutate(ln_stab = log10(stability + 0.1)) %>%
  mutate(sq_stab = stability^2) %>%
  mutate(cb_stab = stability^1/3) %>%
  mutate(qd_stab = stability^1/4) %>%
  mutate(stab3 = stability^3) %>%
  mutate(stab4 = stability^4)

resid1.12 <- lm(data = BRAR_pctstab_trans, ln_stab ~ avgBRAR)
ols_plot_resid_hist(resid1.12) #normal
ols_test_normality(resid1.12) #passed


#redo stats with trans
BRAR_pctstab_stats2 <- lmerTest::lmer(data = BRAR_pctstab_trans, ln_stab ~ avgBRAR +
                                        (1|gradient_number))
anova(BRAR_pctstab_stats2, type = 3) #p = 1.185e-6

rsquared(BRAR_pctstab_stats2) #r sq = 0.6391579
AIC(BRAR_pctstab_stats2) #31.200043


#BRTE - use ln
resid1.13 <- lm(data = BRTE_pctstab, stability ~ avgBRTE)
ols_plot_resid_hist(resid1.13) #right skew
ols_test_normality(resid1.13) #failed

#try trans
BRTE_pctstab_trans <- BRTE_pctstab %>%
  mutate(sqrt_stab = sqrt(stability)) %>%
  mutate(ln_stab = log10(stability + 0.1)) %>%
  mutate(sq_stab = stability^2) %>%
  mutate(cb_stab = stability^1/3) %>%
  mutate(qd_stab = stability^1/4) %>%
  mutate(stab3 = stability^3) %>%
  mutate(stab4 = stability^4)

resid1.14 <- lm(data = BRTE_pctstab_trans, ln_stab ~ avgBRTE)
ols_plot_resid_hist(resid1.14) #normal
ols_test_normality(resid1.14) #passed


#redo stats with trans
BRTE_pctstab_stats2 <- lmerTest::lmer(data = BRTE_pctstab_trans, ln_stab ~ avgBRTE +
                                        (1|gradient_number))
anova(BRTE_pctstab_stats2, type = 3) #p = 3.883e-5

rsquared(BRTE_pctstab_stats2) #r sq = 0.5175495
AIC(BRTE_pctstab_stats2) #30.03156


#other assumptions
# light stability
#BRAR
#linearity
plot(resid(BRAR_pctstab_stats2), BRAR_pctstab_trans$ln_stab) #looks linear
plot(BRAR_pctstab_stats2) #random

#homoscedascity
BRAR_pctstab_trans$res <- residuals(BRAR_pctstab_stats2)
BRAR_pctstab_trans$abs_res <- abs(BRAR_pctstab_trans$res)
BRAR_pctstab_trans$abs_res2 <- BRAR_pctstab_trans$abs_res^2
levene_brar_light2 <- lm(data = BRAR_pctstab_trans, abs_res2 ~ avgBRAR)
anova(levene_brar_light2) #p = 0.5582 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(BRAR_pctstab_stats2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(BRAR_pctstab_stats2, retype = "normalized"))


#BRTE
#linearity
plot(resid(BRTE_pctstab_stats2), BRTE_pctstab_trans$ln_stab) #looks linear
plot(BRTE_pctstab_stats2) #random

#homoscedascity
BRTE_pctstab_trans$res <- residuals(BRTE_pctstab_stats2)
BRTE_pctstab_trans$abs_res <- abs(BRTE_pctstab_trans$res)
BRTE_pctstab_trans$abs_res2 <- BRTE_pctstab_trans$abs_res^2
levene_brte_light2 <- lm(data = BRTE_pctstab_trans, abs_res2 ~ avgBRTE)
anova(levene_brte_light2) #p = 0.6579 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(BRTE_pctstab_stats2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(BRTE_pctstab_stats2, retype = "normalized"))



##soil moisture stability
BRAR_smstab <- community_stability(df = BRARavgab, 
                                   time.var = "year", 
                                   abundance.var = "soil_moist", 
                                   replicate.var = "plot") %>%
  full_join(BRARcov_only2) 


BRTE_smstab <- community_stability(df = BRTEavgab, 
                                   time.var = "year", 
                                   abundance.var = "soil_moist", 
                                   replicate.var = "plot") %>%
  full_join(BRTEcov_only2) 




#stats
BRAR_smstab_stats <- lmerTest::lmer(data = BRAR_smstab, stability ~ avgBRAR +
                                      (1|gradient_number))
anova(BRAR_smstab_stats, type = 3) #p = 0.4561


BRTE_smstab_stats <- lmerTest::lmer(data = BRTE_smstab, stability ~ avgBRTE +
                                      (1|gradient_number))
anova(BRTE_smstab_stats, type = 3) #p = 0.2489


#check normality
#BRAR - use ln
resid1.15 <- lm(data = BRAR_smstab, stability ~ avgBRAR)
ols_plot_resid_hist(resid1.15) #right skew
ols_test_normality(resid1.15) #failed

#try trans
BRAR_smstab_trans <- BRAR_smstab %>%
  mutate(sqrt_stab = sqrt(stability)) %>%
  mutate(ln_stab = log(stability + 0.1)) %>%
  mutate(sq_stab = stability^2) %>%
  mutate(cb_stab = stability^1/3) %>%
  mutate(qd_stab = stability^1/4) %>%
  mutate(stab3 = stability^3) %>%
  mutate(stab4 = stability^4)

resid1.16 <- lm(data = BRAR_smstab_trans, ln_stab ~ avgBRAR)
ols_plot_resid_hist(resid1.16) #normal
ols_test_normality(resid1.16) #passed SW, KS

#redo stats with ln trans
BRAR_smstab_stats2 <- lmerTest::lmer(data = BRAR_smstab_trans, ln_stab ~ avgBRAR +
                                       (1|gradient_number))
anova(BRAR_smstab_stats2, type = 3) #p = 0.6136


#BRTE - use untransformed
resid1.17 <- lm(data = BRTE_smstab, stability ~ avgBRTE)
ols_plot_resid_hist(resid1.17) #fairly normal
ols_test_normality(resid1.17) #passed SW, KS, AD


#check other assumptions
#SM stability
#BRAR
#linearity
plot(resid(BRAR_smstab_stats2), BRAR_smstab_trans$ln_stab) #looks linear
plot(BRAR_smstab_stats2) #random

#homoscedascity
BRAR_smstab_trans$res <- residuals(BRAR_smstab_stats2)
BRAR_smstab_trans$abs_res <- abs(BRAR_smstab_trans$res)
BRAR_smstab_trans$abs_res2 <- BRAR_smstab_trans$abs_res^2
levene_brar_sm2 <- lm(data = BRAR_smstab_trans, abs_res2 ~ avgBRAR)
anova(levene_brar_sm2) #p = 0.1257 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(BRAR_smstab_stats2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(BRAR_smstab_stats2, retype = "normalized"))

#BRTE
#linearity
plot(resid(BRTE_smstab_stats), BRTE_smstab$stability) #looks linear
plot(BRTE_smstab_stats) #random

#homoscedascity
BRTE_smstab$res <- residuals(BRTE_smstab_stats)
BRTE_smstab$abs_res <- abs(BRTE_smstab$res)
BRTE_smstab$abs_res2 <- BRTE_smstab$abs_res^2
levene_brte_sm2 <- lm(data = BRTE_smstab, abs_res2 ~ avgBRTE)
anova(levene_brte_sm2) #p = 0.5975 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(BRTE_smstab_stats, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(BRTE_smstab_stats, retype = "normalized"))




###making combined figure
##light and SM
#light
BRAR_pctstab2 <- ggplot(data = BRAR_pctstab, aes(x = avgBRAR, y = stability)) + 
  geom_point(size = 3) + 
  ylim(-10, 105) +
  stat_smooth(method = "lm", se = FALSE, formula = y ~ log(x)) + 
  labs(x = "BRAR Rel Cov (%)", y = "Light Stability", title = "BRAR") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), axis.title.x = element_blank(), axis.text.x = element_blank()) + 
  theme(panel.border = element_rect(linewidth = 2))

BRTE_pctstab2 <- ggplot(data = BRTE_pctstab, aes(x = avgBRTE, y = stability)) + 
  geom_point(size = 3) + 
  ylim(-10, 100) +
  stat_smooth(method = "lm", se = FALSE, formula = y ~ log(x)) + 
  labs(x = "BRTE Rel Cov (%)", y = "Light Stability", title = "BRTE") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank()) + 
  theme(panel.border = element_rect(linewidth = 2))

#soil moist
BRAR_smstab2 <- ggplot(data = BRAR_smstab, aes(x = avgBRAR, y = stability)) + 
  geom_point(size = 3) + 
  ylim(0.8, 4) +
  #geom_smooth(method = "lm", se = FALSE) + 
  labs(x = "Relative % Cover", y = "Soil Moisture Stability", title = "BRAR") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) + 
  theme(panel.border = element_rect(linewidth = 2))

BRTE_smstab2 <- ggplot(data = BRTE_smstab, aes(x = avgBRTE, y = stability)) + 
  geom_point(size = 3) + 
  ylim(0.8, 4) +
  #geom_smooth(method = "lm", se = FALSE) + 
  labs(x = "Relative % Cover", y = "Soil Moisture Stability", title = "BRTE") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), axis.title.y = element_blank()) + 
  theme(panel.border = element_rect(linewidth = 2))


BRAR_pctstab2 + BRTE_pctstab2 + BRAR_smstab2 + BRTE_smstab2 + plot_layout(ncol = 2)

#4fig 1000x800



###############################################################




#################################################################################
########## Stability with bromes of rich, even, synch, turn, tot cov ###########

#Richness and evenness
comm_struct <- community_structure(df = comp_join, 
                                   time.var = "year", 
                                   abundance.var = "cover", 
                                   replicate.var = "plot", 
                                   metric = "Evar")  #default is Evar, can also calculate EQ and SimpsonEvenness

comm_struct2 <- comm_struct %>%
  full_join(meta)

#join with real brome cov
BRARcov_only <- coverclean %>% #recopied here
  filter(symbol == "BRAR") %>%
  mutate(BRARcov = rel_cov) %>%
  dplyr::select(c(plot, invasive_type, invasion_percent, year, BRARcov))

BRARcov_rich <- BRARcov_only %>%
  full_join(comm_struct) %>%
  filter(invasive_type == "BRAR")


BRTEcov_only <- coverclean %>% #recopied here
  filter(symbol == "BRTE") %>%
  mutate(BRTEcov = rel_cov) %>%
  dplyr::select(c(plot, invasive_type, invasion_percent, year, BRTEcov))

BRTEcov_rich <- BRTEcov_only %>%
  full_join(comm_struct) %>%
  filter(invasive_type == "BRTE")


BRARrichstab <- community_stability(df = BRARcov_rich, 
                                    time.var = "year", 
                                    abundance.var = "richness",     #stability of richness
                                    replicate.var = "plot") %>%
  full_join(meta) %>%
  filter(invasive_type == "BRAR")


BRTErichstab <- community_stability(df = BRTEcov_rich, 
                                    time.var = "year", 
                                    abundance.var = "richness",     #stability of richness
                                    replicate.var = "plot") %>%
  full_join(meta) %>%
  filter(invasive_type == "BRTE")


BRARevenstab <- community_stability(df = BRARcov_rich, 
                                    time.var = "year", 
                                    abundance.var = "Evar",     #stability of evenness
                                    replicate.var = "plot") %>%
  full_join(meta) %>%
  filter(invasive_type == "BRAR")


BRTEevenstab <- community_stability(df = BRTEcov_rich, 
                                    time.var = "year", 
                                    abundance.var = "Evar",     #stability of evenness
                                    replicate.var = "plot") %>%
  full_join(meta) %>%
  filter(invasive_type == "BRTE")




#average across years to get real inv% 
BRARcov_only2 <- BRARcov_only %>%
  group_by(plot) %>%
  summarise(avgBRAR = mean(BRARcov)) %>%
  ungroup()

BRARrichstab2 <- BRARrichstab %>%
  full_join(BRARcov_only2) %>%
  drop_na(invasive_type)

BRARevenstab2 <- BRARevenstab %>%
  full_join(BRARcov_only2) %>%
  drop_na(invasive_type)


BRTEcov_only2 <- BRTEcov_only %>%
  group_by(plot) %>%
  summarise(avgBRTE = mean(BRTEcov)) %>%
  ungroup()

BRTErichstab2 <- BRTErichstab %>%
  full_join(BRTEcov_only2) %>%
  drop_na(invasive_type) %>%
  filter(stability != "Inf")

BRTEevenstab2 <- BRTEevenstab %>%
  full_join(BRTEcov_only2) %>%
  drop_na(invasive_type)

#stats
#richness
BRARrichstab_stats <- lmerTest::lmer(data = BRARrichstab2, stability ~ avgBRAR +
                                       (1|gradient_number))
anova(BRARrichstab_stats, type = 3) #p = 0.3013


BRTErichstab_stats <- lmerTest::lmer(data = BRTErichstab2, stability ~ avgBRTE +
                                       (1|gradient_number))
anova(BRTErichstab_stats, type = 3) #p = 0.03639
rsquared(BRTErichstab_stats) #marg = 0.1932622; cond = 0.1932622

#evenness
BRARevenstab_stats <- lmerTest::lmer(data = BRARevenstab2, stability ~ avgBRAR +
                                       (1|gradient_number))
anova(BRARevenstab_stats, type = 3) #p = 0.02068

BRTEevenstab_stats <- lmerTest::lmer(data = BRTEevenstab2, stability ~ avgBRTE +
                                       (1|gradient_number))
anova(BRTEevenstab_stats, type = 3) #p = 0.2997

#check normality
#BRAR richness
resid1.0 <- lm(data = BRARrichstab2, stability ~ avgBRAR)
ols_plot_resid_hist(resid1.0) #right skew sort of 
ols_test_normality(resid1.0) #passed KS


#try transformation - ln looks good
BRARrichstab2_trans <- BRARrichstab2 %>%
  mutate(sqrt_stab = sqrt(stability)) %>%
  mutate(ln_stab = log(stability + 0.1)) %>%
  mutate(sq_stab = stability^2) %>%
  mutate(cb_stab = stability^1/3) %>%
  mutate(qd_stab = stability^1/4) %>%
  mutate(stab3 = stability^3) %>%
  mutate(stab4 = stability^4)

resid1.1 <- lm(data = BRARrichstab2_trans, ln_stab ~ avgBRAR)
ols_plot_resid_hist(resid1.1) #better 
ols_test_normality(resid1.1) #passed SW, KS, AD

#retry stats
BRARrichstab_stats2 <- lmerTest::lmer(data = BRARrichstab2_trans, ln_stab ~ avgBRAR +
                                        (1|gradient_number))
anova(BRARrichstab_stats2, type = 3) #p = 0.1922


#BRTE richness - stick with untransformed
resid1.2 <- lm(data = BRTErichstab2, stability ~ avgBRTE)
ols_plot_resid_hist(resid1.2) #normal looking 
ols_test_normality(resid1.2) #passed KS, AD


#BRAR evenness
resid1.3 <- lm(data = BRARevenstab2, stability ~ avgBRAR)
ols_plot_resid_hist(resid1.3) #right skew 
ols_test_normality(resid1.3) #failed

#try transformation - ln looks better
BRARevenstab2_trans <- BRARevenstab2 %>%
  mutate(sqrt_stab = sqrt(stability)) %>%
  mutate(ln_stab = log(stability + 0.1)) %>%
  mutate(sq_stab = stability^2) %>%
  mutate(cb_stab = stability^1/3) %>%
  mutate(qd_stab = stability^1/4) %>%
  mutate(stab3 = stability^3) %>%
  mutate(stab4 = stability^4)

resid1.4 <- lm(data = BRARevenstab2_trans, ln_stab ~ avgBRAR)
ols_plot_resid_hist(resid1.4) #right skew but better 
ols_test_normality(resid1.4) #passed KS

#retry stats
BRARevenstab_stats2 <- lmerTest::lmer(data = BRARevenstab2_trans, ln_stab ~ avgBRAR +
                                        (1|gradient_number))
anova(BRARevenstab_stats2, type = 3) #p = 0.1508


#BRTE evenness - stick with untransformed
resid1.5 <- lm(data = BRTEevenstab2, stability ~ avgBRTE)
ols_plot_resid_hist(resid1.5) #normalish
ols_test_normality(resid1.5) #passes SW, KS, AD


#check other assumptions
#richness stability - with bromes
#BRAR
#linearity
plot(resid(BRARrichstab_stats2), BRARrichstab2_trans$ln_stab) #looks linear
plot(BRARrichstab_stats2) #random

#homoscedascity
BRARrichstab2_trans$res <- residuals(BRARrichstab_stats2)
BRARrichstab2_trans$abs_res <- abs(BRARrichstab2_trans$res)
BRARrichstab2_trans$abs_res2 <- BRARrichstab2_trans$abs_res^2
levene_brar_rich2 <- lm(data = BRARrichstab2_trans, abs_res2 ~ avgBRAR)
anova(levene_brar_rich2) #p = 0.437 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(BRARrichstab_stats2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(BRARrichstab_stats2, retype = "normalized"))

#BRTE
#linearity
plot(resid(BRTErichstab_stats), BRTErichstab2$stability) #looks linear
plot(BRTErichstab_stats) #random

#homoscedascity
BRTErichstab2$res <- residuals(BRTErichstab_stats)
BRTErichstab2$abs_res <- abs(BRTErichstab2$res)
BRTErichstab2$abs_res2 <- BRTErichstab2$abs_res^2
levene_brte_rich2 <- lm(data = BRTErichstab2, abs_res2 ~ avgBRTE)
anova(levene_brte_rich2) #p = 0.06559 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(BRTErichstab_stats, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(BRTErichstab_stats, retype = "normalized"))


#evenness stability - with bromes
#BRAR
#linearity
plot(resid(BRARevenstab_stats2), BRARevenstab2_trans$ln_stab) #looks linear
plot(BRARevenstab_stats2) #random

#homoscedascity
BRARevenstab2_trans$res <- residuals(BRARevenstab_stats2)
BRARevenstab2_trans$abs_res <- abs(BRARevenstab2_trans$res)
BRARevenstab2_trans$abs_res2 <- BRARevenstab2_trans$abs_res^2
levene_brar_even2 <- lm(data = BRARevenstab2_trans, abs_res2 ~ avgBRAR)
anova(levene_brar_even2) #p = 0.2629 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(BRARevenstab_stats2, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(BRARevenstab_stats2, retype = "normalized"))

#BRTE
#linearity
plot(resid(BRTEevenstab_stats), BRTEevenstab2$stability) #looks linear
plot(BRTEevenstab_stats) #random

#homoscedascity
BRTEevenstab2$res <- residuals(BRTEevenstab_stats)
BRTEevenstab2$abs_res <- abs(BRTEevenstab2$res)
BRTEevenstab2$abs_res2 <- BRTEevenstab2$abs_res^2
levene_brte_even2 <- lm(data = BRTEevenstab2, abs_res2 ~ avgBRTE)
anova(levene_brte_even2) #p = 0.09174 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(BRTEevenstab_stats, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(BRTEevenstab_stats, retype = "normalized"))



##Synchrony
#bromes included
BRARsynch <- synchrony(df = BRARcompmet, 
                       time.var = "year", 
                       species.var = "symbol", 
                       abundance.var = "cover", 
                       replicate.var = "plot") %>%
  mutate(plot = as.factor(plot)) %>%
  full_join(meta) %>%
  filter(invasive_type == "BRAR") %>%
  rename(synch = synchrony)

BRTEsynch <- synchrony(df = BRTEcompmet, 
                       time.var = "year", 
                       species.var = "symbol", 
                       abundance.var = "cover", 
                       replicate.var = "plot") %>%
  mutate(plot = as.factor(plot)) %>%
  full_join(meta) %>%
  filter(invasive_type == "BRTE") %>%
  rename(synch = synchrony)


#real rel % vs synch
#bromes included
BRARsynch3 <- BRARsynch %>%
  full_join(BRARcov_only2)

BRTEsynch3 <- BRTEsynch %>%
  full_join(BRTEcov_only2)

#stats rel cov vs synch
BRARsync_stats <- lmerTest::lmer(data = BRARsynch3, synch ~ avgBRAR +
                                   (1|gradient_number))
anova(BRARsync_stats, type = 3) #p = 0.3476

BRTEsync_stats <- lmerTest::lmer(data = BRTEsynch3, synch ~ avgBRTE +
                                   (1|gradient_number))
anova(BRTEsync_stats, type = 3) #p = 0.22


#check on normality
res1 <- lm(data = BRARsynch3, synch ~ avgBRAR)
ols_plot_resid_hist(res1) #normal
ols_test_normality(res1) #pass tests

res2 <- lm(data = BRTEsynch3, synch ~ avgBRTE)
ols_plot_resid_hist(res2) #normal
ols_test_normality(res2) #pass 


#check other assumptions
#synchrony stability - with bromes
#BRAR
#linearity
BRARsynch3 <- BRARsynch3 %>%
  drop_na(invasion_percent)
plot(resid(BRARsync_stats), BRARsynch3$synch) #looks linear
plot(BRARsync_stats) #random

#homoscedascity
BRARsynch3$res <- residuals(BRARsync_stats)
BRARsynch3$abs_res <- abs(BRARsynch3$res)
BRARsynch3$abs_res2 <- BRARsynch3$abs_res^2
levene_brar_synch2 <- lm(data = BRARsynch3, abs_res2 ~ avgBRAR)
anova(levene_brar_synch2) #p = 0.1874 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(BRARsync_stats, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(BRARsync_stats, retype = "normalized"))

#BRTE
#linearity
BRTEsynch3 <- BRTEsynch3 %>%
  drop_na(invasion_percent)
plot(resid(BRTEsync_stats), BRTEsynch3$synch) #looks linear
plot(BRTEsync_stats) #random

#homoscedascity
BRTEsynch3$res <- residuals(BRTEsync_stats)
BRTEsynch3$abs_res <- abs(BRTEsynch3$res)
BRTEsynch3$abs_res2 <- BRTEsynch3$abs_res^2
levene_brte_synch2 <- lm(data = BRTEsynch3, abs_res2 ~ avgBRTE)
anova(levene_brte_synch2) #p = 0.1711 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(BRTEsync_stats, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(BRTEsync_stats, retype = "normalized"))




##Turnover
#community turnover - start to end
BRAR_covclean <- coverclean %>%
  filter(invasive_type == "BRAR")
BRTE_covclean <- coverclean %>%
  filter(invasive_type == "BRTE")

BRAR_turn <- turnover(df = subset(BRAR_covclean, year != "2020"), 
                        time.var = "year",
                        species.var = "symbol",
                        abundance.var = "cover", 
                        replicate.var = "plot") %>%
  mutate(plot = as.factor(plot)) %>%
  full_join(meta) %>%
  filter(invasive_type == "BRAR")


BRTE_turn <- turnover(df = subset(BRTE_covclean, year != "2020"), 
                        time.var = "year",
                        species.var = "symbol",
                        abundance.var = "cover", 
                        replicate.var = "plot") %>%
  mutate(plot = as.factor(plot)) %>%
  full_join(meta) %>%
  filter(invasive_type == "BRTE")


#real rel % vs turnover
BRARcov_only2 <- BRARcov_only2 %>%
  mutate(plot = as.factor(plot))
BRAR_turn_cov <- BRARcov_only2 %>%
  full_join(BRAR_turn)

BRTEcov_only2 <- BRTEcov_only2 %>%
  mutate(plot = as.factor(plot))
BRTE_turn_cov <- BRTEcov_only2 %>%
  full_join(BRTE_turn)


#stats rel cov vs turnover
BRAR_turn_cov_stats <- lmerTest::lmer(data = BRAR_turn_cov, total ~ avgBRAR +
                                          (1|gradient_number))
anova(BRAR_turn_cov_stats, type = 3) #p = 0.02802

rsquared(BRAR_turn_cov_stats) #marg = 0.1795923; cond = 0.3397744

BRTE_turn_cov_stats <- lmerTest::lmer(data = BRTE_turn_cov, total ~ avgBRTE +
                                          (1|gradient_number))
anova(BRTE_turn_cov_stats, type = 3) #p = 0.001718

rsquared(BRTE_turn_cov_stats) #marg = 0.2606871; cond = 0.5468884


#check on normality
resid72 <- lm(data = BRAR_turn_cov, total ~ avgBRAR)
ols_plot_resid_hist(resid72) #normal
ols_test_normality(resid72) #pass tests

resid73 <- lm(data = BRTE_turn_cov, total ~ avgBRTE)
ols_plot_resid_hist(resid73) #normal
ols_test_normality(resid73) #pass some tests (SW, KS)


#check other assumptions
#turnover stability - with bromes
#BRAR
#linearity
BRAR_turn_cov <- BRAR_turn_cov %>%
  drop_na(invasion_percent)
plot(resid(BRAR_turn_cov_stats), BRAR_turn_cov$total) #looks linear
plot(BRAR_turn_cov_stats) #random

#homoscedascity
BRAR_turn_cov$res <- residuals(BRAR_turn_cov_stats)
BRAR_turn_cov$abs_res <- abs(BRAR_turn_cov$res)
BRAR_turn_cov$abs_res2 <- BRAR_turn_cov$abs_res^2
levene_brar_turn2 <- lm(data = BRAR_turn_cov, abs_res2 ~ avgBRAR)
anova(levene_brar_turn2) #p = 0.08986 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(BRAR_turn_cov_stats, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(BRAR_turn_cov_stats, retype = "normalized"))

#BRTE
#linearity
BRTE_turn_cov <- BRTE_turn_cov %>%
  drop_na(invasion_percent)
plot(resid(BRTE_turn_cov_stats), BRTE_turn_cov$total) #looks linear
plot(BRTE_turn_cov_stats) #random

#homoscedascity
BRTE_turn_cov$res <- residuals(BRTE_turn_cov_stats)
BRTE_turn_cov$abs_res <- abs(BRTE_turn_cov$res)
BRTE_turn_cov$abs_res2 <- BRTE_turn_cov$abs_res^2
levene_brte_turn2 <- lm(data = BRTE_turn_cov, abs_res2 ~ avgBRTE)
anova(levene_brte_turn2) #p = 0.7141 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(BRTE_turn_cov_stats, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(BRTE_turn_cov_stats, retype = "normalized"))



##Total Cover
BRAR_stab <- community_stability(df = BRAR_covclean, 
                                   time.var = "year", 
                                   abundance.var = "cover", 
                                   replicate.var = "plot") %>%
  full_join(meta) %>%
  filter(invasive_type == "BRAR")


BRTE_stab <- community_stability(df = BRTE_covclean, 
                                   time.var = "year", 
                                   abundance.var = "cover", 
                                   replicate.var = "plot") %>%
  full_join(meta) %>%
  filter(invasive_type == "BRTE")


#stability with real rel % of bromes
#BRAR in BRAR gradietns and BRTE in BRTE gradients only
#need to average cover across years

BRARcov_only2 <- BRARcov_only %>%
  group_by(gradient_number, plot, invasive_type, invasion_percent) %>%
  summarise(avgBRAR = mean(BRARcov)) %>%
  ungroup() %>%
  filter(invasive_type == "BRAR")

BRARcov_only3 <- BRARcov_only2 %>%
  full_join(BRAR_stab)

BRTEcov_only2 <- BRTEcov_only %>%
  group_by(gradient_number, plot, invasive_type, invasion_percent) %>%
  summarise(avgBRTE = mean(BRTEcov)) %>%
  ungroup() %>%
  filter(invasive_type == "BRTE")

BRTEcov_only3 <- BRTEcov_only2 %>%
  full_join(BRTE_stab)

#stats rel cov vs stability
BRAR_stab_stats3 <- lmerTest::lmer(data = BRARcov_only3, stability ~ avgBRAR +
                                       (1|gradient_number))
anova(BRAR_stab_stats3, type = 3) #p = 0.4659


BRTE_stab_stats3 <- lmerTest::lmer(data = BRTEcov_only3, stability ~ avgBRTE +
                                       (1|gradient_number))
anova(BRTE_stab_stats3, type = 3) #p = 0.5546


#check on normality
resid59 <- lm(data = BRARcov_only3, stability ~ avgBRAR)
ols_plot_resid_hist(resid59) #right skew
ols_test_normality(resid59) #fail tests

resid60 <- lm(data = BRTEcov_only3, stability ~ avgBRTE)
ols_plot_resid_hist(resid60) #right skew
ols_test_normality(resid60) #fail tests

#try transformation 
#BRAR - ln looks closest to normal
BRARcov_only3_trans <- BRARcov_only3 %>%
  mutate(sqrt_stab = sqrt(stability)) %>%
  mutate(ln_stab = log(stability + 0.1)) %>%
  mutate(sq_stab = stability^2) %>%
  mutate(cb_stab = stability^1/3) %>%
  mutate(qd_stab = stability^1/4) %>%
  mutate(stab3 = stability^3) %>%
  mutate(stab4 = stability^4)

resid61 <- lm(data = BRARcov_only3_trans, sqrt_stab ~ avgBRAR)
ols_plot_resid_hist(resid61) #right skew
ols_test_normality(resid61) #fail tests

resid62 <- lm(data = BRARcov_only3_trans, ln_stab ~ avgBRAR)
ols_plot_resid_hist(resid62) #right skew
ols_test_normality(resid62) #fail tests except KS

resid63 <- lm(data = BRARcov_only3_trans, sq_stab ~ avgBRAR)
ols_plot_resid_hist(resid63) #right skew
ols_test_normality(resid63) #fail tests

resid64 <- lm(data = BRARcov_only3_trans, cb_stab ~ avgBRAR)
ols_plot_resid_hist(resid64) #right skew
ols_test_normality(resid64) #fail tests

resid65 <- lm(data = BRARcov_only3_trans, qd_stab ~ avgBRAR)
ols_plot_resid_hist(resid65) #right skew
ols_test_normality(resid65) #fail tests

resid66 <- lm(data = BRARcov_only3_trans, stab3 ~ avgBRAR)
ols_plot_resid_hist(resid66) #right skew
ols_test_normality(resid66) #fail tests

resid67 <- lm(data = BRARcov_only3_trans, stab4 ~ avgBRAR)
ols_plot_resid_hist(resid67) #right skew
ols_test_normality(resid67) #fail tests

#BRTE - ln looks closest to normal
BRTEcov_only3_trans <- BRTEcov_only3 %>%
  mutate(sqrt_stab = sqrt(stability)) %>%
  mutate(ln_stab = log(stability + 0.1)) %>%
  mutate(sq_stab = stability^2) %>%
  mutate(cb_stab = stability^1/3) %>%
  mutate(qd_stab = stability^1/4) %>%
  mutate(stab3 = stability^3) %>%
  mutate(stab4 = stability^4)

resid68 <- lm(data = BRTEcov_only3_trans, sqrt_stab ~ avgBRTE)
ols_plot_resid_hist(resid68) #looks closer to normal
ols_test_normality(resid68) #fail tests

resid69 <- lm(data = BRTEcov_only3_trans, ln_stab ~ avgBRTE)
ols_plot_resid_hist(resid69) #looks closer to normal
ols_test_normality(resid69) #fail tests except KS and marginal SW

#try stats again with ln transformations
BRAR_stab_stats4 <- lmerTest::lmer(data = BRARcov_only3_trans, ln_stab ~ avgBRAR +
                                       (1|gradient_number))
anova(BRAR_stab_stats4, type = 3) #p = 0.6794

BRTE_stab_stats4 <- lmerTest::lmer(data = BRTEcov_only3_trans, ln_stab ~ avgBRTE +
                                       (1|gradient_number))
anova(BRTE_stab_stats4, type = 3) #p = 0.7106



#check other assumptions
#cover stability - with bromes
#BRAR
#linearity
plot(resid(BRAR_stab_stats4), BRARcov_only3_trans$ln_stab) #looks linear
plot(BRAR_stab_stats4) #random

#homoscedascity
BRARcov_only3_trans$res <- residuals(BRAR_stab_stats4)
BRARcov_only3_trans$abs_res <- abs(BRARcov_only3_trans$res)
BRARcov_only3_trans$abs_res2 <- BRARcov_only3_trans$abs_res^2
levene_brar_cov2 <- lm(data = BRARcov_only3_trans, abs_res2 ~ avgBRAR)
anova(levene_brar_cov2) #p = 0.3467 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(BRAR_stab_stats4, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(BRAR_stab_stats4, retype = "normalized"))

#BRTE
#linearity
plot(resid(BRTE_stab_stats4), BRTEcov_only3_trans$ln_stab) #looks linear
plot(BRTE_stab_stats4) #random

#homoscedascity
BRTEcov_only3_trans$res <- residuals(BRTE_stab_stats4)
BRTEcov_only3_trans$abs_res <- abs(BRTEcov_only3_trans$res)
BRTEcov_only3_trans$abs_res2 <- BRTEcov_only3_trans$abs_res^2
levene_brte_cov2 <- lm(data = BRTEcov_only3_trans, abs_res2 ~ avgBRTE)
anova(levene_brte_cov2) #p = 0.9293 so > 0.05 so equal variance is met

#autocorrelation
acf(residuals(BRTE_stab_stats4, retype = "normalized"))  #lines dont go outside CI horizontal lines, so not autocorrelated
pacf(residuals(BRTE_stab_stats4, retype = "normalized"))



###############################################################





###########################################################
############# Precip data ###########
#comes from NOAA for converse county
precip$month_year <- factor(precip$month_year, levels = c("Jan_2019", "Feb_2019", "Mar_2019", "Apr_2019", "May_2019", 
                                                          "Jun_2019", "Jul_2019", "Aug_2019", "Sep_2019", "Oct_2019", 
                                                          "Nov_2019", "Dec_2019", "Oct_Jun_2019", "Total_2019", 
                                                          "Jan_2020", "Feb_2020", "Mar_2020", "Apr_2020", "May_2020", 
                                                          "Jun_2020", "Jul_2020", "Aug_2020", "Sep_2020", "Oct_2020", 
                                                          "Nov_2020", "Dec_2020", "Oct_Jun_2020", "Total_2020", 
                                                          "Jan_2021", "Feb_2021", "Mar_2021", "Apr_2021", "May_2021", 
                                                          "Jun_2021", "Jul_2021", "Aug_2021", "Sep_2021", "Oct_2021", 
                                                          "Nov_2021", "Dec_2021", "Oct_Jun_2021", "Total_2021"))


precip_fig <- ggplot(data = precip, aes(x = month_year, y = precip_mm)) +
  geom_bar(stat = "identity") +
  labs(x = "Month", y = "Precipitation (mm)") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  geom_hline(yintercept = 363.293742, linetype = "dashed", color = "red")



##oct to sept on x axis with cumulative smooth line on top
precip2$month_year <- factor(precip2$month_year, levels = c("Oct_2018", "Nov_2018", "Dec_2018", "Jan_2019", "Feb_2019", 
                                                          "Mar_2019", "Apr_2019", "May_2019", "Jun_2019", "Jul_2019",
                                                          "Aug_2019", "Sep_2019", "Oct_2019", "Nov_2019", "Dec_2019", 
                                                          "Jan_2020", "Feb_2020", "Mar_2020", "Apr_2020", "May_2020", 
                                                          "Jun_2020", "Jul_2020", "Aug_2020", "Sep_2020", "Oct_2020", 
                                                          "Nov_2020", "Dec_2020", "Jan_2021", "Feb_2021", "Mar_2021", 
                                                          "Apr_2021", "May_2021", "Jun_2021", "Jul_2021", "Aug_2021", "Sep_2021"))


precip_fig3 <- ggplot(data = precip2) +
  geom_bar(aes(x = month_year, y = precip_mm), stat = "identity") +
  geom_point(aes(x = month_year, y = cum_tot)) + 
  geom_line(aes(x = month_year, y = cum_tot), group = 1) + 
  labs(x = "Month", y = "Precipitation (mm)") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  geom_hline(yintercept = 350.1513, linetype = "dashed", color = "red") + 
  scale_x_discrete(labels=c("Oct_2018" = "Oct", "Nov_2018" = "Nov", "Dec_2018" = "Dec", "Jan_2019" = "Jan", "Feb_2019" = "Feb", 
                            "Mar_2019" = "Mar", "Apr_2019" = "Apr", "May_2019" = "May", "Jun_2019" = "Jun", "Jul_2019" = "Jul",
                            "Aug_2019" = "Aug", "Sep_2019" = "Sep", "Oct_2019" = "Oct", "Nov_2019" = "Nov", "Dec_2019" = "Dec", 
                            "Jan_2020" = "Jan", "Feb_2020" = "Feb", "Mar_2020" = "Mar", "Apr_2020" = "Apr", "May_2020" = "May", 
                            "Jun_2020" = "Jun", "Jul_2020" = "Jul", "Aug_2020" = "Aug", "Sep_2020" = "Sep", "Oct_2020" = "Oct", 
                            "Nov_2020" = "Nov", "Dec_2020" = "Dec", "Jan_2021" = "Jan", "Feb_2021" = "Feb", "Mar_2021" = "Mar", 
                            "Apr_2021" = "Apr", "May_2021" = "May", "Jun_2021" = "Jun", "Jul_2021" = "Jul", "Aug_2021" = "Aug", "Sep_2021" = "Sep")) + 
  theme(panel.border = element_rect(linewidth = 2))

#1800x700  rain



####################################################




####################################################
######### Moderator analyses ########
#are the effects of invasion on stability mediated by the impacts on light or sm stability?

#combine BRAR data
BRARrichstab_nobrome2 <- BRARrichstab_nobrome %>%
  mutate(log_rich = log10(stability)) %>%
  dplyr::select(-stability)
BRARevenstab_nobrome2 <- BRARevenstab_nobrome %>%
  mutate(log_even = log10(stability)) %>%
  dplyr::select(-stability)
BRARsynch_nobrome2 <- BRARsynch_nobrome #%>%
  #mutate(synch = synch[,1]) 
noBRAR_turn2 <- noBRAR_turn %>%
  mutate(turn = total) %>%
  dplyr::select(-total)
noBRAR_stab2 <- noBRAR_stab %>%
  mutate(log_cov = log10(stability)) %>%
  dplyr::select(-stability)
BRARsp_c4_stab2 <- BRARsp_c4_stab %>%
  mutate(log_c4 = log10(stability)) %>%
  dplyr::select(-stability)
BRARsp_c3_stab2 <- BRARsp_c3_nobrome_stab %>% 
  mutate(log_c3 = log10(stability)) %>%
  dplyr::select(-stability)
BRARsp_f_stab2 <- BRARsp_f_stab %>%
  mutate(log_f = log10(stability)) %>%
  dplyr::select(-stability)
BRAR_pctstab2 <- BRAR_pctstab %>%
  mutate(log_light = log10(stability)) %>%
  dplyr::select(-stability)
BRAR_smstab2 <- BRAR_smstab %>%
  mutate(log_sm = log10(stability)) %>%
  dplyr::select(-stability)

BRAR_med <- BRARcover %>%
  full_join(BRARrichstab_nobrome2) %>%
  full_join(BRARevenstab_nobrome2) %>%
  full_join(BRARsynch_nobrome2) %>%
  full_join(noBRAR_turn2) %>%
  full_join(noBRAR_stab2) %>%
  full_join(BRARsp_c4_stab2) %>%
  full_join(BRARsp_c3_stab2) %>%
  full_join(BRARsp_f_stab2) %>%
  full_join(BRAR_pctstab2) %>%
  full_join(BRAR_smstab2) %>%
  dplyr::select(-c(site, plot_type, plot_name, phone_lat, phone_long, garmin_lat, garmin_long, year))

BRAR_med_c4 <- BRAR_med %>%
  filter(log_c4 != "NA")


#combine BRTE data    -try logging rich and C3
BRTErichstab_nobrome2 <- BRTErichstab_nobrome %>%
  mutate(log_rich = log10(stability)) %>%
  dplyr::select(-stability)
BRTEevenstab_nobrome2 <- BRTEevenstab_nobrome %>%
  mutate(log_even = log10(stability)) %>%
  dplyr::select(-stability)
BRTEsynch_nobrome2 <- BRTEsynch_nobrome #%>%
  #mutate(synch = synch[,1]) 
noBRTE_turn2 <- noBRTE_turn %>%
  mutate(turn = total) %>%
  dplyr::select(-total)
noBRTE_stab2 <- noBRTE_stab %>%
  mutate(log_cov = log10(stability)) %>%
  dplyr::select(-stability)
BRTEsp_c4_stab2 <- BRTEsp_c4_stab %>%
  mutate(log_c4 = log10(stability)) %>%
  dplyr::select(-stability)
BRTEsp_c3_stab2 <- BRTEsp_c3_nobrome_stab %>%
  mutate(log_c3 = log10(stability)) %>%
  dplyr::select(-stability)
BRTEsp_f_stab2 <- BRTEsp_f_stab %>%
  mutate(log_f = log10(stability)) %>%
  dplyr::select(-stability)
BRTE_pctstab2 <- BRTE_pctstab %>%
  mutate(log_light = log10(stability)) %>%
  dplyr::select(-stability)
BRTE_smstab2 <- BRTE_smstab %>%
  mutate(sm = stability) %>%
  dplyr::select(-stability)

BRTE_med <- BRTEcover %>%
  full_join(BRTErichstab_nobrome2) %>%
  full_join(BRTEevenstab_nobrome2) %>%
  full_join(BRTEsynch_nobrome2) %>%
  full_join(noBRTE_turn2) %>%
  full_join(noBRTE_stab2) %>%
  full_join(BRTEsp_c4_stab2) %>%
  full_join(BRTEsp_c3_stab2) %>%
  full_join(BRTEsp_f_stab2) %>%
  full_join(BRTE_pctstab2) %>%
  full_join(BRTE_smstab2) %>%
  dplyr::select(-c(site, plot_type, plot_name, phone_lat, phone_long, garmin_lat, garmin_long, year))

BRTE_med_rich <- BRTE_med %>%
  filter(log_rich != "Inf")
BRTE_med_c4 <- BRTE_med %>%
  filter(log_c4 != "Inf" & log_c4 != "NA")
BRTE_med_f <- BRTE_med %>%
  filter(log_f != "Inf" & log_f != "NA")



#est variables
#BRAR
BRAR <- BRAR_med$avgBRAR
BRAR_rich <- BRAR_med$log_rich
BRAR_light <- BRAR_med$log_light
BRAR_sm <- BRAR_med$log_sm
BRAR_even <- BRAR_med$log_even
BRAR_synch <- BRAR_med$synch
BRAR_turn <- BRAR_med$turn
BRAR_cov <- BRAR_med$log_cov
BRAR_c3 <- BRAR_med$log_c3
BRAR_f <- BRAR_med$log_f
BRAR2 <- BRAR_med_c4$avgBRAR
BRAR_light2 <- BRAR_med_c4$log_light
BRAR_sm2 <- BRAR_med_c4$log_sm
BRAR_c4 <- BRAR_med_c4$log_c4

dat1 <- data.frame(BRAR, BRAR_rich, BRAR_light, BRAR_sm, BRAR_even, BRAR_synch, BRAR_turn, BRAR_cov, BRAR_c3, BRAR_f)
dat2 <- data.frame(BRAR2, BRAR_light2, BRAR_sm2, BRAR_c4)


#BRTE
BRTE <- BRTE_med$avgBRTE
BRTE_rich <- BRTE_med_rich$log_rich
BRTE_light <- BRTE_med$log_light
BRTE_sm <- BRTE_med$sm
BRTE_even <- BRTE_med$log_even
BRTE_synch <- BRTE_med$synch
BRTE_turn <- BRTE_med$turn
BRTE_cov <- BRTE_med$log_cov
BRTE_c3 <- BRTE_med$log_c3
BRTE_c4 <- BRTE_med_c4$log_c4
BRTE_f <- BRTE_med_f$log_f
BRTE2 <- BRTE_med_rich$avgBRTE
BRTE3 <- BRTE_med_c4$avgBRTE
BRTE4 <- BRTE_med_f$avgBRTE
BRTE_light2 <- BRTE_med_rich$log_light
BRTE_light3 <- BRTE_med_c4$log_light
BRTE_light4 <- BRTE_med_f$log_light
BRTE_sm2 <- BRTE_med_rich$sm
BRTE_sm3 <- BRTE_med_c4$sm
BRTE_sm4 <- BRTE_med_f$sm

dat3 <- data.frame(BRTE, BRTE_light, BRTE_sm, BRTE_even, BRTE_synch, BRTE_turn, BRTE_cov, BRTE_c3)
dat4 <- data.frame(BRTE2, BRTE_light2, BRTE_sm2, BRTE_rich)
dat5 <- data.frame(BRTE3, BRTE_light3, BRTE_sm3, BRTE_c4)
dat6 <- data.frame(BRTE4, BRTE_light4, BRTE_sm4, BRTE_f)


#SEMs
#BRAR
#BRAR -> abiotic -> richness

#full
BRAR_rich_full <- 'BRAR_light ~ BRAR
                        BRAR_sm ~ BRAR
                        BRAR_rich ~ BRAR_light + BRAR_sm + 0*BRAR'

BRAR_rich_full_sem <- sem(BRAR_rich_full, dat1)
modificationIndices(BRAR_rich_full_sem)

#part
BRAR_rich_part <- 'BRAR_light ~ a*BRAR
                        BRAR_sm ~ b*BRAR
                        BRAR_rich ~ c*BRAR_light + d*BRAR_sm + e*BRAR
                        light_path := c*a
                        sm_path := d*a
                        no_path := e'

BRAR_rich_part_sem <- sem(BRAR_rich_part, dat1)
semPaths(BRAR_rich_part_sem, 'std')  #plotting this gives path coefficients which is the same # as std.all in output

#no
BRAR_rich_no <- 'BRAR_rich ~ BRAR
                  BRAR_light ~ BRAR
                  BRAR_sm ~ BRAR'
BRAR_rich_no_sem <- sem(BRAR_rich_no, dat1)


AICc(BRAR_rich_full_sem, BRAR_rich_part_sem, BRAR_rich_no_sem) #AICc = -9.238633, -5.799985, -1.423548
dist(AICc(BRAR_rich_full_sem, BRAR_rich_part_sem, BRAR_rich_no_sem)[,2])

summary(BRAR_rich_full_sem, standardized = TRUE, rsq = T)
summary(BRAR_rich_part_sem, standardized = TRUE, rsq = T)
summary(BRAR_rich_no_sem, standardized = TRUE, rsq = T)

parameterEstimates(BRAR_rich_part_sem)


#BRAR -> abiotic -> evenness

#full
BRAR_even_full <- 'BRAR_light ~ BRAR
                        BRAR_sm ~ BRAR
                        BRAR_even ~ BRAR_light + BRAR_sm + 0*BRAR'

BRAR_even_full_sem <- sem(BRAR_even_full, dat1)


#part
BRAR_even_part <- 'BRAR_light ~ a*BRAR
                        BRAR_sm ~ b*BRAR
                        BRAR_even ~ c*BRAR_light + d*BRAR_sm + e*BRAR
                        light_path := c*a
                        sm_path := d*a
                        no_path := e'

BRAR_even_part_sem <- sem(BRAR_even_part, dat1)
semPaths(BRAR_even_part_sem, 'std')


#no
BRAR_even_no <- 'BRAR_even ~ BRAR
                BRAR_light ~ BRAR
                  BRAR_sm ~ BRAR'
BRAR_even_no_sem <- sem(BRAR_even_no, dat1)


AICc(BRAR_even_full_sem, BRAR_even_part_sem, BRAR_even_no_sem) #AICc = -10.614766, -6.466872, -2.090435
dist(AICc(BRAR_even_full_sem, BRAR_even_part_sem, BRAR_even_no_sem)[,2])

summary(BRAR_even_full_sem, standardized = TRUE, rsq = T)
summary(BRAR_even_part_sem, standardized = TRUE, rsq = T)
summary(BRAR_even_no_sem, standardized = TRUE, rsq = T)

parameterEstimates(BRAR_even_part_sem)


#BRAR -> abiotic -> synchrony

#full
BRAR_synch_full <- 'BRAR_light ~ BRAR
                        BRAR_sm ~ BRAR
                        BRAR_synch ~ BRAR_light + BRAR_sm + 0*BRAR'

BRAR_synch_full_sem <- sem(BRAR_synch_full, dat1)


#part
BRAR_synch_part <- 'BRAR_light ~ a*BRAR
                        BRAR_sm ~ b*BRAR
                        BRAR_synch ~ c*BRAR_light + d*BRAR_sm + e*BRAR
                        light_path := c*a
                        sm_path := d*a
                        no_path := e'

BRAR_synch_part_sem <- sem(BRAR_synch_part, dat1)


#no
BRAR_synch_no <- 'BRAR_synch ~ BRAR
                  BRAR_light ~ BRAR
                  BRAR_sm ~ BRAR'
BRAR_synch_no_sem <- sem(BRAR_synch_no, dat1)


AICc(BRAR_synch_full_sem, BRAR_synch_part_sem, BRAR_synch_no_sem) #AICc = -26.86803, -22.46963, -18.09319
dist(AICc(BRAR_synch_full_sem, BRAR_synch_part_sem, BRAR_synch_no_sem)[,2])

summary(BRAR_synch_full_sem, standardized = TRUE, rsq = T)
summary(BRAR_synch_part_sem, standardized = TRUE, rsq = T)
summary(BRAR_synch_no_sem, standardized = TRUE, rsq = T)

parameterEstimates(BRAR_synch_part_sem)


#BRAR -> abiotic -> turnover

#full
BRAR_turn_full <- 'BRAR_light ~ BRAR
                        BRAR_sm ~ BRAR
                        BRAR_turn ~ BRAR_light + BRAR_sm + 0*BRAR'

BRAR_turn_full_sem <- sem(BRAR_turn_full, dat1)
modificationIndices(BRAR_turn_full_sem)


#part
BRAR_turn_part <- 'BRAR_light ~ a*BRAR
                        BRAR_sm ~ b*BRAR
                        BRAR_turn ~ c*BRAR_light + d*BRAR_sm + e*BRAR
                        light_path := c*a
                        sm_path := d*a
                        no_path := e'

BRAR_turn_part_sem <- sem(BRAR_turn_part, dat1)


#no
BRAR_turn_no <- 'BRAR_turn ~ BRAR
                  BRAR_light ~ BRAR
                  BRAR_sm ~ BRAR'
BRAR_turn_no_sem <- sem(BRAR_turn_no, dat1)


AICc(BRAR_turn_full_sem, BRAR_turn_part_sem, BRAR_turn_no_sem) #AICc = -38.10840, -41.07802, -36.70159
dist(AICc(BRAR_turn_full_sem, BRAR_turn_part_sem, BRAR_turn_no_sem)[,2])

summary(BRAR_turn_full_sem, standardized = TRUE, rsq = T)
summary(BRAR_turn_part_sem, standardized = TRUE, rsq = T)
summary(BRAR_turn_no_sem, standardized = TRUE, rsq = T)

parameterEstimates(BRAR_turn_part_sem)


#BRAR -> abiotic -> cover

#full
BRAR_cov_full <- 'BRAR_light ~ BRAR
                        BRAR_sm ~ BRAR
                        BRAR_cov ~ BRAR_light + BRAR_sm + 0*BRAR'

BRAR_cov_full_sem <- sem(BRAR_cov_full, dat1)
modificationIndices(BRAR_cov_full_sem)


#part
BRAR_cov_part <- 'BRAR_light ~ a*BRAR
                        BRAR_sm ~ b*BRAR
                        BRAR_cov ~ c*BRAR_light + d*BRAR_sm + e*BRAR
                        light_path := c*a
                        sm_path := d*a
                        no_path := e'

BRAR_cov_part_sem <- sem(BRAR_cov_part, dat1)
semPaths(BRAR_cov_part_sem, "std")

#no
BRAR_cov_no <- 'BRAR_cov ~ BRAR
                BRAR_light ~ BRAR
                  BRAR_sm ~ BRAR'
BRAR_cov_no_sem <- sem(BRAR_cov_no, dat1)


AICc(BRAR_cov_full_sem, BRAR_cov_part_sem, BRAR_cov_no_sem) #AICc = -0.627510, 2.907321, 7.283757
dist(AICc(BRAR_cov_full_sem, BRAR_cov_part_sem, BRAR_cov_no_sem)[,2])

summary(BRAR_cov_full_sem, standardized = TRUE, rsq = T)
summary(BRAR_cov_part_sem, standardized = TRUE, rsq = T)
summary(BRAR_cov_no_sem, standardized = TRUE, rsq = T)

parameterEstimates(BRAR_cov_part_sem)


#BRAR -> abiotic -> c4

#full
BRAR_c4_full <- 'BRAR_light2 ~ BRAR2
                        BRAR_sm2 ~ BRAR2
                        BRAR_c4 ~ BRAR_light2 + BRAR_sm2 + 0*BRAR2'

BRAR_c4_full_sem <- sem(BRAR_c4_full, dat2)
modificationIndices(BRAR_c4_full_sem)


#part
BRAR_c4_part <- 'BRAR_light2 ~ a*BRAR2
                        BRAR_sm2 ~ b*BRAR2
                        BRAR_c4 ~ c*BRAR_light2 + d*BRAR_sm2 + e*BRAR2
                        light_path := c*a
                        sm_path := d*a
                        no_path := e'

BRAR_c4_part_sem <- sem(BRAR_c4_part, dat2)


#no
BRAR_c4_no <- 'BRAR_c4 ~ BRAR2
              BRAR_light2 ~ BRAR2
                  BRAR_sm2 ~ BRAR2'
BRAR_c4_no_sem <- sem(BRAR_c4_no, dat2)


AICc(BRAR_c4_full_sem, BRAR_c4_part_sem, BRAR_c4_no_sem) #AICc = -30.44380, -25.53730, -20.02804
dist(AICc(BRAR_c4_full_sem, BRAR_c4_part_sem, BRAR_c4_no_sem)[,2])

summary(BRAR_c4_full_sem, standardized = TRUE, rsq = T)
summary(BRAR_c4_part_sem, standardized = TRUE, rsq = T)
summary(BRAR_c4_no_sem, standardized = TRUE, rsq = T)

parameterEstimates(BRAR_c4_part_sem)


#BRAR -> abiotic -> c3

#full
BRAR_c3_full <- 'BRAR_light ~ BRAR
                        BRAR_sm ~ BRAR
                        BRAR_c3 ~ BRAR_light + BRAR_sm + 0*BRAR'

BRAR_c3_full_sem <- sem(BRAR_c3_full, dat1)


#part
BRAR_c3_part <- 'BRAR_light ~ a*BRAR
                        BRAR_sm ~ b*BRAR
                        BRAR_c3 ~ c*BRAR_light + d*BRAR_sm + e*BRAR
                        light_path := c*a
                        sm_path := d*a
                        no_path := e'

BRAR_c3_part_sem <- sem(BRAR_c3_part, dat1)


#no
BRAR_c3_no <- 'BRAR_c3 ~ BRAR
              BRAR_light ~ BRAR
                  BRAR_sm ~ BRAR'
BRAR_c3_no_sem <- sem(BRAR_c3_no, dat1)


AICc(BRAR_c3_full_sem, BRAR_c3_part_sem, BRAR_c3_no_sem) #AICc = -6.583123, -5.427579, -1.051142
dist(AICc(BRAR_c3_full_sem, BRAR_c3_part_sem, BRAR_c3_no_sem)[,2])

summary(BRAR_c3_full_sem, standardized = TRUE, rsq = T)
summary(BRAR_c3_part_sem, standardized = TRUE, rsq = T)
summary(BRAR_c3_no_sem, standardized = TRUE, rsq = T)

parameterEstimates(BRAR_c3_part_sem)


#BRAR -> abiotic -> forb

#full
BRAR_f_full <- 'BRAR_light ~ BRAR
                        BRAR_sm ~ BRAR
                        BRAR_f ~ BRAR_light + BRAR_sm + 0*BRAR'

BRAR_f_full_sem <- sem(BRAR_f_full, dat1)


#part
BRAR_f_part <- 'BRAR_light ~ a*BRAR
                        BRAR_sm ~ b*BRAR
                        BRAR_f ~ c*BRAR_light + d*BRAR_sm + e*BRAR
                        light_path := c*a
                        sm_path := d*a
                        no_path := e'

BRAR_f_part_sem <- sem(BRAR_f_part, dat1)


#no
BRAR_f_no <- 'BRAR_f ~ BRAR
              BRAR_light ~ BRAR
                  BRAR_sm ~ BRAR'
BRAR_f_no_sem <- sem(BRAR_f_no, dat1)


AICc(BRAR_f_full_sem, BRAR_f_part_sem, BRAR_f_no_sem) #AICc = -0.4142513, 3.6788339, 8.0552706
dist(AICc(BRAR_f_full_sem, BRAR_f_part_sem, BRAR_f_no_sem)[,2])

summary(BRAR_f_full_sem, standardized = TRUE, rsq = T)
summary(BRAR_f_part_sem, standardized = TRUE, rsq = T)
summary(BRAR_f_no_sem, standardized = TRUE, rsq = T)

parameterEstimates(BRAR_f_part_sem)



#BRTE
#BRTE -> abiotic -> richness

#full
BRTE_rich_full <- 'BRTE_light2 ~ BRTE2
                        BRTE_sm2 ~ BRTE2
                        BRTE_rich ~ BRTE_light2 + BRTE_sm2 + 0*BRTE2'

BRTE_rich_full_sem <- sem(BRTE_rich_full, dat4)

#part
BRTE_rich_part <- 'BRTE_light2 ~ a*BRTE2
                        BRTE_sm2 ~ b*BRTE2
                        BRTE_rich ~ c*BRTE_light2 + d*BRTE_sm2 + e*BRTE2
                        light_path := c*a
                        sm_path := d*a
                        no_path := e'

BRTE_rich_part_sem <- sem(BRTE_rich_part, dat4)


#no
BRTE_rich_no <- 'BRTE_rich ~ BRTE2
                 BRTE_light2 ~ BRTE2
                  BRTE_sm2 ~ BRTE2'
BRTE_rich_no_sem <- sem(BRTE_rich_no, dat4)


AICc(BRTE_rich_full_sem, BRTE_rich_part_sem, BRTE_rich_no_sem) #AICc = 48.89077, 47.07021, 52.07467
dist(AICc(BRTE_rich_full_sem, BRTE_rich_part_sem, BRTE_rich_no_sem)[,2])

summary(BRTE_rich_full_sem, standardized = TRUE, rsq = T)
summary(BRTE_rich_part_sem, standardized = TRUE, rsq = T)
summary(BRTE_rich_no_sem, standardized = TRUE, rsq = T)

parameterEstimates(BRTE_rich_part_sem)


#BRTE -> abiotic -> evenness

#full
BRTE_even_full <- 'BRTE_light ~ BRTE
                        BRTE_sm ~ BRTE
                        BRTE_even ~ BRTE_light + BRTE_sm + 0*BRTE'

BRTE_even_full_sem <- sem(BRTE_even_full, dat3)

#part
BRTE_even_part <- 'BRTE_light ~ a*BRTE
                        BRTE_sm ~ b*BRTE
                        BRTE_even ~ c*BRTE_light + d*BRTE_sm + e*BRTE
                        light_path := c*a
                        sm_path := d*a
                        no_path := e'

BRTE_even_part_sem <- sem(BRTE_even_part, dat3)


#no
BRTE_even_no <- 'BRTE_even ~ BRTE
                BRTE_light ~ BRTE
                  BRTE_sm ~ BRTE'
BRTE_even_no_sem <- sem(BRTE_even_no, dat3)


AICc(BRTE_even_full_sem, BRTE_even_part_sem, BRTE_even_no_sem) #AICc = 46.7554900, 50.2616171, 53.89370
dist(AICc(BRTE_even_full_sem, BRTE_even_part_sem, BRTE_even_no_sem)[,2])

summary(BRTE_even_full_sem, standardized = TRUE, rsq = T)
summary(BRTE_even_part_sem, standardized = TRUE, rsq = T)
summary(BRTE_even_no_sem, standardized = TRUE, rsq = T)

parameterEstimates(BRTE_even_part_sem)

#BRTE -> abiotic -> synchrony

#full
BRTE_synch_full <- 'BRTE_light ~ BRTE
                        BRTE_sm ~ BRTE
                        BRTE_synch ~ BRTE_light + BRTE_sm + 0*BRTE'

BRTE_synch_full_sem <- sem(BRTE_synch_full, dat3)

#part
BRTE_synch_part <- 'BRTE_light ~ a*BRTE
                        BRTE_sm ~ b*BRTE
                        BRTE_synch ~ c*BRTE_light + d*BRTE_sm + e*BRTE
                        light_path := c*a
                        sm_path := d*a
                        no_path := e'

BRTE_synch_part_sem <- sem(BRTE_synch_part, dat3)


#no
BRTE_synch_no <- 'BRTE_synch ~ BRTE
                  BRTE_light ~ BRTE
                  BRTE_sm ~ BRTE'
BRTE_synch_no_sem <- sem(BRTE_synch_no, dat3)


AICc(BRTE_synch_full_sem, BRTE_synch_part_sem, BRTE_synch_no_sem) #AICc = 18.62488, 22.91623, 26.54831
dist(AICc(BRTE_synch_full_sem, BRTE_synch_part_sem, BRTE_synch_no_sem)[,2])

summary(BRTE_synch_full_sem, standardized = TRUE, rsq = T)
summary(BRTE_synch_part_sem, standardized = TRUE, rsq = T)
summary(BRTE_synch_no_sem, standardized = TRUE, rsq = T)

parameterEstimates(BRTE_synch_part_sem)

#BRTE -> abiotic -> turnover

#full
BRTE_turn_full <- 'BRTE_light ~ BRTE
                        BRTE_sm ~ BRTE
                        BRTE_turn ~ BRTE_light + BRTE_sm + 0*BRTE'

BRTE_turn_full_sem <- sem(BRTE_turn_full, dat3)

#part
BRTE_turn_part <- 'BRTE_light ~ a*BRTE
                        BRTE_sm ~ b*BRTE
                        BRTE_turn ~ c*BRTE_light + d*BRTE_sm + e*BRTE
                        light_path := c*a
                        sm_path := d*a
                        no_path := e'

BRTE_turn_part_sem <- sem(BRTE_turn_part, dat3)


#no
BRTE_turn_no <- 'BRTE_turn ~ BRTE
                BRTE_light ~ BRTE
                  BRTE_sm ~ BRTE'
BRTE_turn_no_sem <- sem(BRTE_turn_no, dat3)


AICc(BRTE_turn_full_sem, BRTE_turn_part_sem, BRTE_turn_no_sem) #AICc = 24.76209, 22.40981, 26.04190
dist(AICc(BRTE_turn_full_sem, BRTE_turn_part_sem, BRTE_turn_no_sem)[,2])

summary(BRTE_turn_full_sem, standardized = TRUE, rsq = T)
summary(BRTE_turn_part_sem, standardized = TRUE, rsq = T)
summary(BRTE_turn_no_sem, standardized = TRUE, rsq = T)

parameterEstimates(BRTE_turn_part_sem)

#BRTE -> abiotic -> cover

#full
BRTE_cov_full <- 'BRTE_light ~ BRTE
                        BRTE_sm ~ BRTE
                        BRTE_cov ~ BRTE_light + BRTE_sm + 0*BRTE'

BRTE_cov_full_sem <- sem(BRTE_cov_full, dat3)

#part
BRTE_cov_part <- 'BRTE_light ~ a*BRTE
                        BRTE_sm ~ b*BRTE
                        BRTE_cov ~ c*BRTE_light + d*BRTE_sm + e*BRTE
                        light_path := c*a
                        sm_path := d*a
                        no_path := e'

BRTE_cov_part_sem <- sem(BRTE_cov_part, dat3)


#no
BRTE_cov_no <- 'BRTE_cov ~ BRTE
                BRTE_light ~ BRTE
                  BRTE_sm ~ BRTE'
BRTE_cov_no_sem <- sem(BRTE_cov_no, dat3)


AICc(BRTE_cov_full_sem, BRTE_cov_part_sem, BRTE_cov_no_sem) #AICc = 51.957537, 56.343503, 59.97559
dist(AICc(BRTE_cov_full_sem, BRTE_cov_part_sem, BRTE_cov_no_sem)[,2])

summary(BRTE_cov_full_sem, standardized = TRUE, rsq = T)
summary(BRTE_cov_part_sem, standardized = TRUE, rsq = T)
summary(BRTE_cov_no_sem, standardized = TRUE, rsq = T)

parameterEstimates(BRTE_cov_part_sem)

#BRTE -> abiotic -> c4

#full
BRTE_c4_full <- 'BRTE_light3 ~ BRTE3
                        BRTE_sm3 ~ BRTE3
                        BRTE_c4 ~ BRTE_light3 + BRTE_sm3 + 0*BRTE3'

BRTE_c4_full_sem <- sem(BRTE_c4_full, dat5)

#part
BRTE_c4_part <- 'BRTE_light3 ~ a*BRTE3
                        BRTE_sm3 ~ b*BRTE3
                        BRTE_c4 ~ c*BRTE_light3 + d*BRTE_sm3 + e*BRTE3
                        light_path := c*a
                        sm_path := d*a
                        no_path := e'

BRTE_c4_part_sem <- sem(BRTE_c4_part, dat5)


#no
BRTE_c4_no <- 'BRTE_c4 ~ BRTE3
                BRTE_light3 ~ BRTE3
                  BRTE_sm3 ~ BRTE3'
BRTE_c4_no_sem <- sem(BRTE_c4_no, dat5)


AICc(BRTE_c4_full_sem, BRTE_c4_part_sem, BRTE_c4_no_sem) #AICc = 33.92555, 38.91950, 43.68155
dist(AICc(BRTE_c4_full_sem, BRTE_c4_part_sem, BRTE_c4_no_sem)[,2])

summary(BRTE_c4_full_sem, standardized = TRUE, rsq = T)
summary(BRTE_c4_part_sem, standardized = TRUE, rsq = T)
summary(BRTE_c4_no_sem, standardized = TRUE, rsq = T)

parameterEstimates(BRTE_c4_part_sem)

#BRTE -> abiotic -> c3

#full
BRTE_c3_full <- 'BRTE_light ~ BRTE
                        BRTE_sm ~ BRTE
                        BRTE_c3 ~ BRTE_light + BRTE_sm + 0*BRTE'

BRTE_c3_full_sem <- sem(BRTE_c3_full, dat3)

#part
BRTE_c3_part <- 'BRTE_light ~ a*BRTE
                        BRTE_sm ~ b*BRTE
                        BRTE_c3 ~ c*BRTE_light + d*BRTE_sm + e*BRTE
                        light_path := c*a
                        sm_path := d*a
                        no_path := e'

BRTE_c3_part_sem <- sem(BRTE_c3_part, dat3)


#no
BRTE_c3_no <- 'BRTE_c3 ~ BRTE
              BRTE_light ~ BRTE
                  BRTE_sm ~ BRTE'
BRTE_c3_no_sem <- sem(BRTE_c3_no, dat3)


AICc(BRTE_c3_full_sem, BRTE_c3_part_sem, BRTE_c3_no_sem) #AICc = 24.00440, 23.48915, 27.12124
dist(AICc(BRTE_c3_full_sem, BRTE_c3_part_sem, BRTE_c3_no_sem)[,2])

summary(BRTE_c3_full_sem, standardized = TRUE, rsq = T)
summary(BRTE_c3_part_sem, standardized = TRUE, rsq = T)
summary(BRTE_c3_no_sem, standardized = TRUE, rsq = T)

parameterEstimates(BRTE_c3_part_sem)

#BRTE -> abiotic -> forb

#full
BRTE_f_full <- 'BRTE_light4 ~ BRTE4
                        BRTE_sm4 ~ BRTE4
                        BRTE_f ~ BRTE_light4 + BRTE_sm4 + 0*BRTE4'

BRTE_f_full_sem <- sem(BRTE_f_full, dat6)

#part
BRTE_f_part <- 'BRTE_light4 ~ a*BRTE4
                        BRTE_sm4 ~ b*BRTE4
                        BRTE_f ~ c*BRTE_light4 + d*BRTE_sm4 + e*BRTE4
                        light_path := c*a
                        sm_path := d*a
                        no_path := e'

BRTE_f_part_sem <- sem(BRTE_f_part, dat6)


#no
BRTE_f_no <- 'BRTE_f ~ BRTE4
              BRTE_light4 ~ BRTE4
                  BRTE_sm4 ~ BRTE4'
BRTE_f_no_sem <- sem(BRTE_f_no, dat6)


AICc(BRTE_f_full_sem, BRTE_f_part_sem, BRTE_f_no_sem) #AICc = 43.8772968, 48.5570161, 53.63985
dist(AICc(BRTE_f_full_sem, BRTE_f_part_sem, BRTE_f_no_sem)[,2])

summary(BRTE_f_full_sem, standardized = TRUE, rsq = T)
summary(BRTE_f_part_sem, standardized = TRUE, rsq = T)
summary(BRTE_f_no_sem, standardized = TRUE, rsq = T)

parameterEstimates(BRTE_f_part_sem)

###############################################################




##################################################################
######### Moderator analyses with avg light and sm  ##########
#combine BRAR data
BRARrichstab_nobrome2 <- BRARrichstab_nobrome %>%
  mutate(log_rich = log10(stability)) %>%
  dplyr::select(-stability)
BRARevenstab_nobrome2 <- BRARevenstab_nobrome %>%
  mutate(log_even = log10(stability)) %>%
  dplyr::select(-stability)
BRARsynch_nobrome2 <- BRARsynch_nobrome #%>%
  #mutate(synch = synch[,1]) 
noBRAR_turn2 <- noBRAR_turn %>%
  mutate(turn = total) %>%
  dplyr::select(-total)
noBRAR_stab2 <- noBRAR_stab %>%
  mutate(log_cov = log10(stability)) %>%
  dplyr::select(-stability)
BRARsp_c4_stab2 <- BRARsp_c4_stab %>%
  mutate(log_c4 = log10(stability)) %>%
  dplyr::select(-stability)
BRARsp_c3_stab2 <- BRARsp_c3_nobrome_stab %>% 
  mutate(log_c3 = log10(stability)) %>%
  dplyr::select(-stability)
BRARsp_f_stab2 <- BRARsp_f_stab %>%
  mutate(log_f = log10(stability)) %>%
  dplyr::select(-stability)
BRARavgab_group2 <- BRARavgab_group %>%
  mutate(log_light = log10(avg_light), log_sm = log10(avg_moist)) 

BRAR_med_avg <- BRARcover %>%
  full_join(BRARrichstab_nobrome2) %>%
  full_join(BRARevenstab_nobrome2) %>%
  full_join(BRARsynch_nobrome2) %>%
  full_join(noBRAR_turn2) %>%
  full_join(noBRAR_stab2) %>%
  full_join(BRARsp_c4_stab2) %>%
  full_join(BRARsp_c3_stab2) %>%
  full_join(BRARsp_f_stab2) %>%
  full_join(BRARavgab_group2) %>%
  dplyr::select(-c(site, plot_type, plot_name, phone_lat, phone_long, garmin_lat, garmin_long, year))

BRAR_med_c4_avg <- BRAR_med_avg %>%
  filter(log_c4 != "NA")


#combine BRTE data  
BRTErichstab_nobrome2 <- BRTErichstab_nobrome %>%
  mutate(log_rich = log10(stability)) %>%
  dplyr::select(-stability)
BRTEevenstab_nobrome2 <- BRTEevenstab_nobrome %>%
  mutate(log_even = log10(stability)) %>%
  dplyr::select(-stability)
BRTEsynch_nobrome2 <- BRTEsynch_nobrome #%>%
  #mutate(synch = synch[,1]) 
noBRTE_turn2 <- noBRTE_turn %>%
  mutate(log_turn = log10(total)) %>%
  dplyr::select(-total)
noBRTE_stab2 <- noBRTE_stab %>%
  mutate(log_cov = log10(stability)) %>%
  dplyr::select(-stability)
BRTEsp_c4_stab2 <- BRTEsp_c4_stab %>%
  mutate(log_c4 = log10(stability)) %>%
  dplyr::select(-stability)
BRTEsp_c3_stab2 <- BRTEsp_c3_nobrome_stab %>%
  mutate(log_c3 = log10(stability)) %>%
  dplyr::select(-stability)
BRTEsp_f_stab2 <- BRTEsp_f_stab %>%
  mutate(log_f = log10(stability)) %>%
  dplyr::select(-stability)
BRTEavgab_group2 <- BRTEavgab_group %>%
  mutate(log_light = log10(avg_light), log_sm = log10(avg_moist)) 

BRTE_med_avg <- BRTEcover %>%
  full_join(BRTErichstab_nobrome2) %>%
  full_join(BRTEevenstab_nobrome2) %>%
  full_join(BRTEsynch_nobrome2) %>%
  full_join(noBRTE_turn2) %>%
  full_join(noBRTE_stab2) %>%
  full_join(BRTEsp_c4_stab2) %>%
  full_join(BRTEsp_c3_stab2) %>%
  full_join(BRTEsp_f_stab2) %>%
  full_join(BRTEavgab_group2) %>% 
  dplyr::select(-c(site, plot_type, plot_name, phone_lat, phone_long, garmin_lat, garmin_long, year))

BRTE_med_rich_avg <- BRTE_med_avg %>%
  filter(log_rich != "Inf")
BRTE_med_c4_avg <- BRTE_med_avg %>%
  filter(log_c4 != "Inf" & log_c4 != "NA")
BRTE_med_f_avg <- BRTE_med_avg %>%
  filter(log_f != "Inf" & log_f != "NA")



#est variables
#BRAR
BRAR <- BRAR_med_avg$avgBRAR
BRAR_rich <- BRAR_med_avg$log_rich
BRAR_light <- BRAR_med_avg$log_light
BRAR_sm <- BRAR_med_avg$log_sm
BRAR_even <- BRAR_med_avg$log_even
BRAR_synch <- BRAR_med_avg$synch
BRAR_turn <- BRAR_med_avg$turn
BRAR_cov <- BRAR_med_avg$log_cov
BRAR_c3 <- BRAR_med_avg$log_c3
BRAR_f <- BRAR_med_avg$log_f
BRAR2 <- BRAR_med_c4_avg$avgBRAR
BRAR_light2 <- BRAR_med_c4_avg$log_light
BRAR_sm2 <- BRAR_med_c4_avg$log_sm
BRAR_c4 <- BRAR_med_c4_avg$log_c4

dat1_avg <- data.frame(BRAR, BRAR_rich, BRAR_light, BRAR_sm, BRAR_even, BRAR_synch, BRAR_turn, BRAR_cov, BRAR_c3, BRAR_f)
dat2_avg <- data.frame(BRAR2, BRAR_light2, BRAR_sm2, BRAR_c4)


#BRTE
BRTE <- BRTE_med_avg$avgBRTE
BRTE_rich <- BRTE_med_rich_avg$log_rich
BRTE_light <- BRTE_med_avg$log_light
BRTE_sm <- BRTE_med_avg$log_sm
BRTE_even <- BRTE_med_avg$log_even
BRTE_synch <- BRTE_med_avg$synch
BRTE_turn <- BRTE_med_avg$log_turn
BRTE_cov <- BRTE_med_avg$log_cov
BRTE_c3 <- BRTE_med_avg$log_c3
BRTE_c4 <- BRTE_med_c4_avg$log_c4
BRTE_f <- BRTE_med_f_avg$log_f
BRTE2 <- BRTE_med_rich_avg$avgBRTE
BRTE3 <- BRTE_med_c4_avg$avgBRTE
BRTE4 <- BRTE_med_f_avg$avgBRTE
BRTE_light2 <- BRTE_med_rich_avg$log_light
BRTE_light3 <- BRTE_med_c4_avg$log_light
BRTE_light4 <- BRTE_med_f_avg$log_light
BRTE_sm2 <- BRTE_med_rich_avg$log_sm
BRTE_sm3 <- BRTE_med_c4_avg$log_sm
BRTE_sm4 <- BRTE_med_f_avg$log_sm

dat3_avg <- data.frame(BRTE, BRTE_light, BRTE_sm, BRTE_even, BRTE_synch, BRTE_turn, BRTE_cov, BRTE_c3)
dat4_avg <- data.frame(BRTE2, BRTE_light2, BRTE_sm2, BRTE_rich)
dat5_avg <- data.frame(BRTE3, BRTE_light3, BRTE_sm3, BRTE_c4)
dat6_avg <- data.frame(BRTE4, BRTE_light4, BRTE_sm4, BRTE_f)


#SEMs
#BRAR
#BRAR -> abiotic -> richness

#full
BRAR_rich_full <- 'BRAR_light ~ BRAR
                        BRAR_sm ~ BRAR
                        BRAR_rich ~ BRAR_light + BRAR_sm + 0*BRAR'

BRAR_rich_full_sem <- sem(BRAR_rich_full, dat1_avg)
modificationIndices(BRAR_rich_full_sem)

#part
BRAR_rich_part <- 'BRAR_light ~ a*BRAR
                        BRAR_sm ~ b*BRAR
                        BRAR_rich ~ c*BRAR_light + d*BRAR_sm + e*BRAR
                        light_path := c*a
                        sm_path := d*a
                        no_path := e'

BRAR_rich_part_sem <- sem(BRAR_rich_part, dat1_avg)
semPaths(BRAR_rich_part_sem, 'std')  #plotting this gives path coefficients which is the same # as std.all in output

#no
BRAR_rich_no <- 'BRAR_rich ~ BRAR
                  BRAR_light ~ BRAR
                  BRAR_sm ~ BRAR'
BRAR_rich_no_sem <- sem(BRAR_rich_no, dat1_avg)


AICc(BRAR_rich_full_sem, BRAR_rich_part_sem, BRAR_rich_no_sem) #AICc = -110.4711, -106.1079, -101.3350
dist(AICc(BRAR_rich_full_sem, BRAR_rich_part_sem, BRAR_rich_no_sem)[,2])

summary(BRAR_rich_full_sem, standardized = TRUE, rsq = T)
summary(BRAR_rich_part_sem, standardized = TRUE, rsq = T)
summary(BRAR_rich_no_sem, standardized = TRUE, rsq = T)

parameterEstimates(BRAR_rich_part_sem)

#BRAR -> abiotic -> evenness

#full
BRAR_even_full <- 'BRAR_light ~ BRAR
                        BRAR_sm ~ BRAR
                        BRAR_even ~ BRAR_light + BRAR_sm + 0*BRAR'

BRAR_even_full_sem <- sem(BRAR_even_full, dat1_avg)


#part
BRAR_even_part <- 'BRAR_light ~ a*BRAR
                        BRAR_sm ~ b*BRAR
                        BRAR_even ~ c*BRAR_light + d*BRAR_sm + e*BRAR
                        light_path := c*a
                        sm_path := d*a
                        no_path := e'

BRAR_even_part_sem <- sem(BRAR_even_part, dat1_avg)
semPaths(BRAR_even_part_sem, 'std')


#no
BRAR_even_no <- 'BRAR_even ~ BRAR
                BRAR_light ~ BRAR
                  BRAR_sm ~ BRAR'
BRAR_even_no_sem <- sem(BRAR_even_no, dat1_avg)


AICc(BRAR_even_full_sem, BRAR_even_part_sem, BRAR_even_no_sem) #AICc = -98.11819, -93.73436, -88.96151
dist(AICc(BRAR_even_full_sem, BRAR_even_part_sem, BRAR_even_no_sem)[,2])

summary(BRAR_even_full_sem, standardized = TRUE, rsq = T)
summary(BRAR_even_part_sem, standardized = TRUE, rsq = T)
summary(BRAR_even_no_sem, standardized = TRUE, rsq = T)

parameterEstimates(BRAR_even_part_sem)

#BRAR -> abiotic -> synchrony

#full
BRAR_synch_full <- 'BRAR_light ~ BRAR
                        BRAR_sm ~ BRAR
                        BRAR_synch ~ BRAR_light + BRAR_sm + 0*BRAR'

BRAR_synch_full_sem <- sem(BRAR_synch_full, dat1_avg)


#part
BRAR_synch_part <- 'BRAR_light ~ a*BRAR
                        BRAR_sm ~ b*BRAR
                        BRAR_synch ~ c*BRAR_light + d*BRAR_sm + e*BRAR
                        light_path := c*a
                        sm_path := d*a
                        no_path := e'

BRAR_synch_part_sem <- sem(BRAR_synch_part, dat1_avg)


#no
BRAR_synch_no <- 'BRAR_synch ~ BRAR
                  BRAR_light ~ BRAR
                  BRAR_sm ~ BRAR'
BRAR_synch_no_sem <- sem(BRAR_synch_no, dat1_avg)


AICc(BRAR_synch_full_sem, BRAR_synch_part_sem, BRAR_synch_no_sem) #AICc = -115.1902, -111.1641, -106.3913
dist(AICc(BRAR_synch_full_sem, BRAR_synch_part_sem, BRAR_synch_no_sem)[,2])

summary(BRAR_synch_full_sem, standardized = TRUE, rsq = T)
summary(BRAR_synch_part_sem, standardized = TRUE, rsq = T)
summary(BRAR_synch_no_sem, standardized = TRUE, rsq = T)

parameterEstimates(BRAR_synch_part_sem)

#BRAR -> abiotic -> turnover

#full
BRAR_turn_full <- 'BRAR_light ~ BRAR
                        BRAR_sm ~ BRAR
                        BRAR_turn ~ BRAR_light + BRAR_sm + 0*BRAR'

BRAR_turn_full_sem <- sem(BRAR_turn_full, dat1_avg)
modificationIndices(BRAR_turn_full_sem)


#part
BRAR_turn_part <- 'BRAR_light ~ a*BRAR
                        BRAR_sm ~ b*BRAR
                        BRAR_turn ~ c*BRAR_light + d*BRAR_sm + e*BRAR
                        light_path := c*a
                        sm_path := d*a
                        no_path := e'

BRAR_turn_part_sem <- sem(BRAR_turn_part, dat1_avg)


#no
BRAR_turn_no <- 'BRAR_turn ~ BRAR
                  BRAR_light ~ BRAR
                  BRAR_sm ~ BRAR'
BRAR_turn_no_sem <- sem(BRAR_turn_no, dat1_avg)


AICc(BRAR_turn_full_sem, BRAR_turn_part_sem, BRAR_turn_no_sem) #AICc = -136.1224, -134.8026, -130.0298
dist(AICc(BRAR_turn_full_sem, BRAR_turn_part_sem, BRAR_turn_no_sem)[,2])

summary(BRAR_turn_full_sem, standardized = TRUE, rsq = T)
summary(BRAR_turn_part_sem, standardized = TRUE, rsq = T)
summary(BRAR_turn_no_sem, standardized = TRUE, rsq = T)

parameterEstimates(BRAR_turn_part_sem)

#BRAR -> abiotic -> cover

#full
BRAR_cov_full <- 'BRAR_light ~ BRAR
                        BRAR_sm ~ BRAR
                        BRAR_cov ~ BRAR_light + BRAR_sm + 0*BRAR'

BRAR_cov_full_sem <- sem(BRAR_cov_full, dat1_avg)
modificationIndices(BRAR_cov_full_sem)


#part
BRAR_cov_part <- 'BRAR_light ~ a*BRAR
                        BRAR_sm ~ b*BRAR
                        BRAR_cov ~ c*BRAR_light + d*BRAR_sm + e*BRAR
                        light_path := c*a
                        sm_path := d*a
                        no_path := e'

BRAR_cov_part_sem <- sem(BRAR_cov_part, dat1_avg)
semPaths(BRAR_cov_part_sem, "std")

#no
BRAR_cov_no <- 'BRAR_cov ~ BRAR
                BRAR_light ~ BRAR
                  BRAR_sm ~ BRAR'
BRAR_cov_no_sem <- sem(BRAR_cov_no, dat1_avg)


AICc(BRAR_cov_full_sem, BRAR_cov_part_sem, BRAR_cov_no_sem) #AICc = -90.39920, -86.00212, -81.22927
dist(AICc(BRAR_cov_full_sem, BRAR_cov_part_sem, BRAR_cov_no_sem)[,2])

summary(BRAR_cov_full_sem, standardized = TRUE, rsq = T)
summary(BRAR_cov_part_sem, standardized = TRUE, rsq = T)
summary(BRAR_cov_no_sem, standardized = TRUE, rsq = T)

parameterEstimates(BRAR_cov_part_sem)

#BRAR -> abiotic -> c4

#full
BRAR_c4_full <- 'BRAR_light2 ~ BRAR2
                        BRAR_sm2 ~ BRAR2
                        BRAR_c4 ~ BRAR_light2 + BRAR_sm2 + 0*BRAR2'

BRAR_c4_full_sem <- sem(BRAR_c4_full, dat2_avg)
modificationIndices(BRAR_c4_full_sem)


#part
BRAR_c4_part <- 'BRAR_light2 ~ a*BRAR2
                        BRAR_sm2 ~ b*BRAR2
                        BRAR_c4 ~ c*BRAR_light2 + d*BRAR_sm2 + e*BRAR2
                        light_path := c*a
                        sm_path := d*a
                        no_path := e'

BRAR_c4_part_sem <- sem(BRAR_c4_part, dat2_avg)


#no
BRAR_c4_no <- 'BRAR_c4 ~ BRAR2
              BRAR_light2 ~ BRAR2
                  BRAR_sm2 ~ BRAR2'
BRAR_c4_no_sem <- sem(BRAR_c4_no, dat2_avg)


AICc(BRAR_c4_full_sem, BRAR_c4_part_sem, BRAR_c4_no_sem) #AICc = -105.99580, -101.02872, -95.67352
dist(AICc(BRAR_c4_full_sem, BRAR_c4_part_sem, BRAR_c4_no_sem)[,2])

summary(BRAR_c4_full_sem, standardized = TRUE, rsq = T)
summary(BRAR_c4_part_sem, standardized = TRUE, rsq = T)
summary(BRAR_c4_no_sem, standardized = TRUE, rsq = T)

parameterEstimates(BRAR_c4_part_sem)

#BRAR -> abiotic -> c3

#full
BRAR_c3_full <- 'BRAR_light ~ BRAR
                        BRAR_sm ~ BRAR
                        BRAR_c3 ~ BRAR_light + BRAR_sm + 0*BRAR'

BRAR_c3_full_sem <- sem(BRAR_c3_full, dat1_avg)


#part
BRAR_c3_part <- 'BRAR_light ~ a*BRAR
                        BRAR_sm ~ b*BRAR
                        BRAR_c3 ~ c*BRAR_light + d*BRAR_sm + e*BRAR
                        light_path := c*a
                        sm_path := d*a
                        no_path := e'

BRAR_c3_part_sem <- sem(BRAR_c3_part, dat1_avg)


#no
BRAR_c3_no <- 'BRAR_c3 ~ BRAR
              BRAR_light ~ BRAR
                  BRAR_sm ~ BRAR'
BRAR_c3_no_sem <- sem(BRAR_c3_no, dat1_avg)


AICc(BRAR_c3_full_sem, BRAR_c3_part_sem, BRAR_c3_no_sem) #AICc = -94.16712, -90.98952, -86.21667
dist(AICc(BRAR_c3_full_sem, BRAR_c3_part_sem, BRAR_c3_no_sem)[,2])

summary(BRAR_c3_full_sem, standardized = TRUE, rsq = T)
summary(BRAR_c3_part_sem, standardized = TRUE, rsq = T)
summary(BRAR_c3_no_sem, standardized = TRUE, rsq = T)

parameterEstimates(BRAR_c3_part_sem)

#BRAR -> abiotic -> forb

#full
BRAR_f_full <- 'BRAR_light ~ BRAR
                        BRAR_sm ~ BRAR
                        BRAR_f ~ BRAR_light + BRAR_sm + 0*BRAR'

BRAR_f_full_sem <- sem(BRAR_f_full, dat1_avg)


#part
BRAR_f_part <- 'BRAR_light ~ a*BRAR
                        BRAR_sm ~ b*BRAR
                        BRAR_f ~ c*BRAR_light + d*BRAR_sm + e*BRAR
                        light_path := c*a
                        sm_path := d*a
                        no_path := e'

BRAR_f_part_sem <- sem(BRAR_f_part, dat1_avg)


#no
BRAR_f_no <- 'BRAR_f ~ BRAR
              BRAR_light ~ BRAR
                  BRAR_sm ~ BRAR'
BRAR_f_no_sem <- sem(BRAR_f_no, dat1_avg)


AICc(BRAR_f_full_sem, BRAR_f_part_sem, BRAR_f_no_sem) #AICc = -87.96009, -83.87868, -79.10583
dist(AICc(BRAR_f_full_sem, BRAR_f_part_sem, BRAR_f_no_sem)[,2])

summary(BRAR_f_full_sem, standardized = TRUE, rsq = T)
summary(BRAR_f_part_sem, standardized = TRUE, rsq = T)
summary(BRAR_f_no_sem, standardized = TRUE, rsq = T)

parameterEstimates(BRAR_f_part_sem)


#BRTE
#BRTE -> abiotic -> richness

#full
BRTE_rich_full <- 'BRTE_light2 ~ BRTE2
                        BRTE_sm2 ~ BRTE2
                        BRTE_rich ~ BRTE_light2 + BRTE_sm2 + 0*BRTE2'

BRTE_rich_full_sem <- sem(BRTE_rich_full, dat4_avg)

#part
BRTE_rich_part <- 'BRTE_light2 ~ a*BRTE2
                        BRTE_sm2 ~ b*BRTE2
                        BRTE_rich ~ c*BRTE_light2 + d*BRTE_sm2 + e*BRTE2
                        light_path := c*a
                        sm_path := d*a
                        no_path := e'

BRTE_rich_part_sem <- sem(BRTE_rich_part, dat4_avg)


#no
BRTE_rich_no <- 'BRTE_rich ~ BRTE2
                 BRTE_light2 ~ BRTE2
                  BRTE_sm2 ~ BRTE2'
BRTE_rich_no_sem <- sem(BRTE_rich_no, dat4_avg)


AICc(BRTE_rich_full_sem, BRTE_rich_part_sem, BRTE_rich_no_sem) #AICc = -116.7593, -114.3819, -109.0426
dist(AICc(BRTE_rich_full_sem, BRTE_rich_part_sem, BRTE_rich_no_sem)[,2])

summary(BRTE_rich_full_sem, standardized = TRUE, rsq = T)
summary(BRTE_rich_part_sem, standardized = TRUE, rsq = T)
summary(BRTE_rich_no_sem, standardized = TRUE, rsq = T)

parameterEstimates(BRTE_rich_part_sem)

#BRTE -> abiotic -> evenness

#full
BRTE_even_full <- 'BRTE_light ~ BRTE
                        BRTE_sm ~ BRTE
                        BRTE_even ~ BRTE_light + BRTE_sm + 0*BRTE'

BRTE_even_full_sem <- sem(BRTE_even_full, dat3_avg)

#part
BRTE_even_part <- 'BRTE_light ~ a*BRTE
                        BRTE_sm ~ b*BRTE
                        BRTE_even ~ c*BRTE_light + d*BRTE_sm + e*BRTE
                        light_path := c*a
                        sm_path := d*a
                        no_path := e'

BRTE_even_part_sem <- sem(BRTE_even_part, dat3_avg)


#no
BRTE_even_no <- 'BRTE_even ~ BRTE
                BRTE_light ~ BRTE
                  BRTE_sm ~ BRTE'
BRTE_even_no_sem <- sem(BRTE_even_no, dat3_avg)


AICc(BRTE_even_full_sem, BRTE_even_part_sem, BRTE_even_no_sem) #AICc = -140.8815, -138.3583, -133.9879
dist(AICc(BRTE_even_full_sem, BRTE_even_part_sem, BRTE_even_no_sem)[,2])

summary(BRTE_even_full_sem, standardized = TRUE, rsq = T)
summary(BRTE_even_part_sem, standardized = TRUE, rsq = T)
summary(BRTE_even_no_sem, standardized = TRUE, rsq = T)

parameterEstimates(BRTE_even_part_sem)

#BRTE -> abiotic -> synchrony

#full
BRTE_synch_full <- 'BRTE_light ~ BRTE
                        BRTE_sm ~ BRTE
                        BRTE_synch ~ BRTE_light + BRTE_sm + 0*BRTE'

BRTE_synch_full_sem <- sem(BRTE_synch_full, dat3_avg)

#part
BRTE_synch_part <- 'BRTE_light ~ a*BRTE
                        BRTE_sm ~ b*BRTE
                        BRTE_synch ~ c*BRTE_light + d*BRTE_sm + e*BRTE
                        light_path := c*a
                        sm_path := d*a
                        no_path := e'

BRTE_synch_part_sem <- sem(BRTE_synch_part, dat3_avg)


#no
BRTE_synch_no <- 'BRTE_synch ~ BRTE
                  BRTE_light ~ BRTE
                  BRTE_sm ~ BRTE'
BRTE_synch_no_sem <- sem(BRTE_synch_no, dat3_avg)


AICc(BRTE_synch_full_sem, BRTE_synch_part_sem, BRTE_synch_no_sem) #AICc = -168.6463, -164.3764, -160.00060
dist(AICc(BRTE_synch_full_sem, BRTE_synch_part_sem, BRTE_synch_no_sem)[,2])

summary(BRTE_synch_full_sem, standardized = TRUE, rsq = T)
summary(BRTE_synch_part_sem, standardized = TRUE, rsq = T)
summary(BRTE_synch_no_sem, standardized = TRUE, rsq = T)

parameterEstimates(BRTE_synch_part_sem)

#BRTE -> abiotic -> turnover

#full
BRTE_turn_full <- 'BRTE_light ~ BRTE
                        BRTE_sm ~ BRTE
                        BRTE_turn ~ BRTE_light + BRTE_sm + 0*BRTE'

BRTE_turn_full_sem <- sem(BRTE_turn_full, dat3_avg)

#part
BRTE_turn_part <- 'BRTE_light ~ a*BRTE
                        BRTE_sm ~ b*BRTE
                        BRTE_turn ~ c*BRTE_light + d*BRTE_sm + e*BRTE
                        light_path := c*a
                        sm_path := d*a
                        no_path := e'

BRTE_turn_part_sem <- sem(BRTE_turn_part, dat3_avg)


#no
BRTE_turn_no <- 'BRTE_turn ~ BRTE
                BRTE_light ~ BRTE
                  BRTE_sm ~ BRTE'
BRTE_turn_no_sem <- sem(BRTE_turn_no, dat3_avg)


AICc(BRTE_turn_full_sem, BRTE_turn_part_sem, BRTE_turn_no_sem) #AICc = -163.2049, -172.2777, -167.9073
dist(AICc(BRTE_turn_full_sem, BRTE_turn_part_sem, BRTE_turn_no_sem)[,2])

summary(BRTE_turn_full_sem, standardized = TRUE, rsq = T)
summary(BRTE_turn_part_sem, standardized = TRUE, rsq = T)
summary(BRTE_turn_no_sem, standardized = TRUE, rsq = T)

parameterEstimates(BRTE_turn_part_sem)

#BRTE -> abiotic -> cover

#full
BRTE_cov_full <- 'BRTE_light ~ BRTE
                        BRTE_sm ~ BRTE
                        BRTE_cov ~ BRTE_light + BRTE_sm + 0*BRTE'

BRTE_cov_full_sem <- sem(BRTE_cov_full, dat3_avg)

#part
BRTE_cov_part <- 'BRTE_light ~ a*BRTE
                        BRTE_sm ~ b*BRTE
                        BRTE_cov ~ c*BRTE_light + d*BRTE_sm + e*BRTE
                        light_path := c*a
                        sm_path := d*a
                        no_path := e'

BRTE_cov_part_sem <- sem(BRTE_cov_part, dat3_avg)


#no
BRTE_cov_no <- 'BRTE_cov ~ BRTE
                BRTE_light ~ BRTE
                  BRTE_sm ~ BRTE'
BRTE_cov_no_sem <- sem(BRTE_cov_no, dat3_avg)


AICc(BRTE_cov_full_sem, BRTE_cov_part_sem, BRTE_cov_no_sem) #AICc = -134.5811, -131.0055, -126.6351
dist(AICc(BRTE_cov_full_sem, BRTE_cov_part_sem, BRTE_cov_no_sem)[,2])

summary(BRTE_cov_full_sem, standardized = TRUE, rsq = T)
summary(BRTE_cov_part_sem, standardized = TRUE, rsq = T)
summary(BRTE_cov_no_sem, standardized = TRUE, rsq = T)

parameterEstimates(BRTE_cov_part_sem)

#BRTE -> abiotic -> c4

#full
BRTE_c4_full <- 'BRTE_light3 ~ BRTE3
                        BRTE_sm3 ~ BRTE3
                        BRTE_c4 ~ BRTE_light3 + BRTE_sm3 + 0*BRTE3'

BRTE_c4_full_sem <- sem(BRTE_c4_full, dat5_avg)

#part
BRTE_c4_part <- 'BRTE_light3 ~ a*BRTE3
                        BRTE_sm3 ~ b*BRTE3
                        BRTE_c4 ~ c*BRTE_light3 + d*BRTE_sm3 + e*BRTE3
                        light_path := c*a
                        sm_path := d*a
                        no_path := e'

BRTE_c4_part_sem <- sem(BRTE_c4_part, dat5_avg)


#no
BRTE_c4_no <- 'BRTE_c4 ~ BRTE3
                BRTE_light3 ~ BRTE3
                  BRTE_sm3 ~ BRTE3'
BRTE_c4_no_sem <- sem(BRTE_c4_no, dat5_avg)


AICc(BRTE_c4_full_sem, BRTE_c4_part_sem, BRTE_c4_no_sem) #AICc = -131.0676, -126.5454, -120.6776
dist(AICc(BRTE_c4_full_sem, BRTE_c4_part_sem, BRTE_c4_no_sem)[,2])

summary(BRTE_c4_full_sem, standardized = TRUE, rsq = T)
summary(BRTE_c4_part_sem, standardized = TRUE, rsq = T)
summary(BRTE_c4_no_sem, standardized = TRUE, rsq = T)

parameterEstimates(BRTE_c4_part_sem)

#BRTE -> abiotic -> c3

#full
BRTE_c3_full <- 'BRTE_light ~ BRTE
                        BRTE_sm ~ BRTE
                        BRTE_c3 ~ BRTE_light + BRTE_sm + 0*BRTE'

BRTE_c3_full_sem <- sem(BRTE_c3_full, dat3_avg)

#part
BRTE_c3_part <- 'BRTE_light ~ a*BRTE
                        BRTE_sm ~ b*BRTE
                        BRTE_c3 ~ c*BRTE_light + d*BRTE_sm + e*BRTE
                        light_path := c*a
                        sm_path := d*a
                        no_path := e'

BRTE_c3_part_sem <- sem(BRTE_c3_part, dat3_avg)


#no
BRTE_c3_no <- 'BRTE_c3 ~ BRTE
              BRTE_light ~ BRTE
                  BRTE_sm ~ BRTE'
BRTE_c3_no_sem <- sem(BRTE_c3_no, dat3_avg)


AICc(BRTE_c3_full_sem, BRTE_c3_part_sem, BRTE_c3_no_sem) #AICc = -164.2728, -165.8779, -161.5075
dist(AICc(BRTE_c3_full_sem, BRTE_c3_part_sem, BRTE_c3_no_sem)[,2])

summary(BRTE_c3_full_sem, standardized = TRUE, rsq = T)
summary(BRTE_c3_part_sem, standardized = TRUE, rsq = T)
summary(BRTE_c3_no_sem, standardized = TRUE, rsq = T)

parameterEstimates(BRTE_c3_part_sem)

#BRTE -> abiotic -> forb

#full
BRTE_f_full <- 'BRTE_light4 ~ BRTE4
                        BRTE_sm4 ~ BRTE4
                        BRTE_f ~ BRTE_light4 + BRTE_sm4 + 0*BRTE4'

BRTE_f_full_sem <- sem(BRTE_f_full, dat6_avg)

#part
BRTE_f_part <- 'BRTE_light4 ~ a*BRTE4
                        BRTE_sm4 ~ b*BRTE4
                        BRTE_f ~ c*BRTE_light4 + d*BRTE_sm4 + e*BRTE4
                        light_path := c*a
                        sm_path := d*a
                        no_path := e'

BRTE_f_part_sem <- sem(BRTE_f_part, dat6_avg)


#no
BRTE_f_no <- 'BRTE_f ~ BRTE4
              BRTE_light4 ~ BRTE4
                  BRTE_sm4 ~ BRTE4'
BRTE_f_no_sem <- sem(BRTE_f_no, dat6_avg)


AICc(BRTE_f_full_sem, BRTE_f_part_sem, BRTE_f_no_sem) #AICc = -131.7868, -127.1600, -121.5251
dist(AICc(BRTE_f_full_sem, BRTE_f_part_sem, BRTE_f_no_sem)[,2])

summary(BRTE_f_full_sem, standardized = TRUE, rsq = T)
summary(BRTE_f_part_sem, standardized = TRUE, rsq = T)
summary(BRTE_f_no_sem, standardized = TRUE, rsq = T)

parameterEstimates(BRTE_f_part_sem)


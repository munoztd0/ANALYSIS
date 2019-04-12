## R code for FOR REWOD_PIT
# last modified on Nov 2018 by David

# -----------------------  PRELIMINARY STUFF ----------------------------------------
#if(!require(pacman)) {
#install.packages("pacman")
#library(pacman)
#}
#pacman::p_load(car, lme4, lmerTest, pbkrtest, ggplot2, dplyr, plyr, tidyr, multcomp, mvoutlier, HH, doBy, psych, pastecs, reshape, reshape2, 
#jtools, effects, compute.es, DescTools, MBESS, afex, ez, metafor, influence.ME)

#require(lattice)
rm(list = ls())

# load library
library(lme4)
library(lmerTest)
library(ggplot2)
library(dplyr)
library(plyr)
library(pastecs)
library(plotrix)
library(Rmisc)
library(influence.ME)
#library(doBy)
#library(afex)
#library(sjstats)
#library(gamlss)
#library(multcomp)
#library(pbkrtest)
#library(mvoutlier)
#library(HH)

#SETUP

# Set working directory
analysis_path <- '~/rewod/DATABASES/'# for this to work the script needs to be sourced
setwd(analysis_path)

# open dataset
REWOD_PIT <- read.delim(file.path(analysis_path,'REWOD_PIT.txt'), header = T, sep ='') # read in dataset

## subsetting into 3 differents tasks
REWOD_PIT.all <- REWOD_PIT
REWOD_RIM <- subset (REWOD_PIT.all,task == 'Reminder') 
REWOD_PE <- subset (REWOD_PIT.all,task == 'Partial_Extinction') 
REWOD_PIT <- subset (REWOD_PIT.all,task == 'PIT') 


# define factors
REWOD_RIM$id               <- factor(REWOD_RIM$id)
REWOD_RIM$trial            <- factor(REWOD_RIM$trial)
REWOD_RIM$task              <- factor(REWOD_RIM$task)
REWOD_RIM$session          <- factor(REWOD_RIM$session)
REWOD_RIM$reward        <- factor(REWOD_RIM$reward)

REWOD_PE$id               <- factor(REWOD_PE$id)
REWOD_PE$trial            <- factor(REWOD_PE$trial)
REWOD_PE$task              <- factor(REWOD_PE$task)
REWOD_PE$session          <- factor(REWOD_PE$session)
REWOD_PE$reward        <- factor(REWOD_PE$reward)

REWOD_PIT$id               <- factor(REWOD_PIT$id)
REWOD_PIT$trial            <- factor(REWOD_PIT$trial)
REWOD_PIT$task              <- factor(REWOD_PIT$task)
REWOD_PIT$session          <- factor(REWOD_PIT$session)


# PLOTS
REWOD_PIT$Condition[REWOD_PIT$condition== 'CSplus']     <- 'CS+'
REWOD_PIT$Condition[REWOD_PIT$condition== 'CSminus']     <- 'CS-'
REWOD_PIT$Condition[REWOD_PIT$condition== 'Baseline']     <- 'Baseline'

##plot (non-averaged per participant) 
#n_grips RIM
boxplot(REWOD_RIM$n_grips ~ REWOD_RIM$trial, las = 1)
#n_grips PE
boxplot(REWOD_PE$n_grips ~ REWOD_PE$trial, las = 1)
#n_grips PIT
boxplot(REWOD_PIT$n_grips ~ REWOD_PIT$trialxcondition, las = 1)


## plot overall effect
# get means by trialxcondition
RIM.bt = ddply(REWOD_RIM, .(trial), summarise,  n_grips = mean(n_grips, na.rm = TRUE)) 
PE.bt = ddply(REWOD_PE, .(trial), summarise,  n_grips = mean(n_grips, na.rm = TRUE)) 
PIT.bt = ddply(REWOD_PIT, .(trialxcondition), summarise,  n_grips = mean(n_grips, na.rm = TRUE)) 

# get means by condition
PIT.bc = ddply(REWOD_PIT, .(condition), summarise,  n_grips = mean(n_grips, na.rm = TRUE)) 

# get means by trial & condition
PIT.bct = ddply(REWOD_PIT, .(trialxcondition, condition), summarise,  n_grips = mean(n_grips, na.rm = TRUE)) 

# get means by participant 
RIM.bs = ddply(REWOD_RIM, .(id, trial), summarise, n_grips = mean(n_grips, na.rm = TRUE)) #not condition
PE.bs = ddply(REWOD_PE, .(id, trial), summarise, n_grips = mean(n_grips, na.rm = TRUE)) #not condition
PIT.bs = ddply(REWOD_PIT, .(id, Condition), summarise, n_grips = mean(n_grips, na.rm = TRUE)) 

# ngrips average per trial
boxplot(RIM.bt$n_grips ~ RIM.bt$trial, las = 1)
boxplot(PE.bt$n_grips ~ PE.bt$trial, las = 1)
boxplot(PIT.bt$n_grips ~ PIT.bt$trialxcondition, las = 1)



##plot n_grips to see the trajectory of learning (overall average by trials)
ggplot(PIT.bt, aes(x = trialxcondition, y = n_grips, fill = I('royalblue1'), color = I('royalblue4'))) +
  geom_point() + geom_line(group=1) +
  guides(color = "none", fill = "none") +
  guides(color = "none", fill = "none") +
  theme_bw() +
  labs(
    title = "Pavlovian Instrumental Transfer",
    x = "Trial",
    y = "Number of grips"
  )

#OR different representation
ggplotRegression <- function (fit) {
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

# plot number of grips by time with regression lign
ggplotRegression(lm(n_grips ~ trialxcondition, data = PIT.bt))

##plot n_grips to see the trajectory of learning (overall average by trials) by conditions
df <- summarySE(REWOD_PIT, measurevar="n_grips", groupvars=c("trialxcondition", "Condition"))
ggplot(df, aes(x = trialxcondition, y = n_grips, color=Condition)) +
  geom_line(aes(linetype = Condition), alpha = .7, size = 1) +
  geom_point() +
  geom_errorbar(aes(ymax = n_grips + 0.5*se, ymin = n_grips - 0.5*se), width=0.1, alpha=0.7, size=0.4)+
  scale_colour_manual(values = c("CS+"="blue", "CS-"="red", "Baseline"="black")) +
  scale_linetype_manual(values = c("CS+"="dashed", "CS-"="twodash", "Baseline"="solid")) +
  scale_x_continuous(breaks=c(1:15)) + 
  scale_y_continuous(breaks=c(5,7.5,10,12.5,15)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size=20,  face="bold"), axis.title = element_text(size=14), axis.text=element_text(size=14), legend.position = c(0.9, 0.9), legend.title=element_blank()) +
  labs(
    title = "Learning Trajectory: PIT",
    x = "Trial",
    y = "Number of grips"
  )


#Bar plot (360 trial x condition)

# summarySE provides the standard deviation, standard error of the mean, and a (default 95%) confidence interval
df1 <- summarySE(REWOD_PIT, measurevar="n_grips", groupvars=c("condition"))
Condition <- c("Baseline", "CS-", "CS+")
df1 <- data.frame(df1, Condition)

  ggplot(PIT.bs, aes(x = Condition, y = n_grips)) +
    geom_jitter(width = 0.05, color="blue") +
    geom_bar(data=df1, stat="identity", fill="skyblue", alpha=0.6, width=0.35, position = position_dodge(width = 0.01)) +
    geom_errorbar(data=df1, aes(x = Condition, ymax = n_grips + se, ymin = n_grips - se), width=0.1, colour="black", alpha=1, size=0.4)+
    geom_segment(aes(x = 1.5, y = 13.5, xend = 3, yend = 13.5), size =0.4)+
    geom_segment(aes(x = 1, y = 10, xend = 2, yend = 10), size =0.4)+
    geom_segment(aes(x = 1.5, y = 10, xend = 1.5, yend = 13.5), size =0.4)+  
    geom_segment(aes(x = 1, y = 9.5, xend = 1, yend = 10), size =0.4)+  
    geom_segment(aes(x = 2, y = 9.5, xend = 2, yend = 10), size =0.4)+ 
    geom_segment(aes(x = 3, y = 13, xend = 3, yend = 13.5), size =0.4)+
    annotate("text", x = 2.25, y = 14, label = "***") +
    scale_x_discrete(limits=c("Baseline","CS-","CS+"))+
    scale_y_continuous(breaks=c(0,5, 10, 15, 20,25))+
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size=20,  face="bold"), axis.title.y = element_text(size=14), axis.text=element_text(size=14), legend.position = c(0.9, 0.85)) +
    labs(
      title = "Pavlovian Instrumental Transfer",
      x = NULL,
      y = "Number of grips"
    )
  
#text_low <- textGrob("Figure 3. Means (+/- 1 SEM) of the numbers of grips per condition during the PIT, Planned contrast: X(12) = 2.41, p = .XX, d = .21, *** = p < .001", gp=gpar(fontsize=7))
  #annotation_custom(text_low,xmin=1.5,xmax=0.9,ymin=-1.3,ymax=-1.3)+
  #gt <- ggplot_gtable(ggplot_build(p))
  #gt$layout$clip[gt$layout$name == "panel"] <- "off"
  #grid.draw(gt)
  #aspect_ratio <- 2.5
#height <- 7
#ggsave(g, height = 7 , width = 7 * aspect_ratio)
  
# ANALYSIS


## 1. number of grips: are participants gripping more on the CSplus condition? 
  #factorise trial and condition
  REWOD_PIT$trialxcondition            <- factor(REWOD_PIT$trial)
  REWOD_PIT$condition       <- factor(REWOD_PIT$condition)
  #Assumptions:
  my.model <- lmer(n_grips ~ condition + (1|id) + (1|trialxcondition), data = REWOD_PIT, REML = FALSE)  #1+cvalue\id or 1\id
  #1)Linearity (not good?)
  plot(my.model)
  #2) Absence of collinearity
  #3)Homoscedasticity AND #4)Normality of residuals
  qqnorm(residuals(my.model))
  #5) Absence of influential data points (ID =4 and ID = 22 ??)
  alt.est.id <- influence(model=my.model, group="id")
  alt.est.trialxcondition <- influence(model=my.model, group="trialxcondition")
  dfbetas(alt.est.id)
  dfbetas(alt.est.trialxcondition)
  #raising the mean by  
  #11.944444 - (7.419444+7.436111)/2

#contrasts
REWOD_PIT$cvalue[REWOD_PIT$condition== 'CSplus']     <- 2
REWOD_PIT$cvalue[REWOD_PIT$condition== 'CSminus']     <- -1
REWOD_PIT$cvalue[REWOD_PIT$condition== 'Baseline']     <- -1
REWOD_PIT$cvalue       <- factor(REWOD_PIT$cvalue)

# lmer analyis ~ condition 
main.n_grips = lmer(n_grips ~ cvalue + (1+cvalue|id) + (1|trialxcondition), data = REWOD_PIT, REML = FALSE) 
#main.n_grips = lmer(n_grips ~ cvalue + (1+cvalue|id) + (1+cvalue|trialxcondition), data = REWOD_PIT, REML = FALSE) #1+cvalue\id or 1\id
anova(main.n_grips)

  
# quick check with classical anova (! this is not reliable)
summary(aov(n_grips ~ cvalue + Error(id / (cvalue)), data = REWOD_PIT))

# model comparison
main.n_grips.0 = lmer(n_grips ~ (1+cvalue|id) + (1|trialxcondition), data = REWOD_PIT, REML = FALSE)
anova(main.n_grips.0, main.n_grips, test = 'Chisq')

#sentence => main.n_grips is signifincatly better than the null model








##############  Base cs minus ########## ??

#contrasts
REWOD_PIT$cvalue1[REWOD_PIT$condition== 'CSplus']     <- 0
REWOD_PIT$cvalue1[REWOD_PIT$condition== 'CSminus']     <- -1
REWOD_PIT$cvalue1[REWOD_PIT$condition== 'Baseline']     <- 1
REWOD_PIT$cvalue1       <- factor(REWOD_PIT$cvalue1)

# lmer analyis ~ condition 
main.n_gripsX = lmer(n_grips ~ cvalue1 + (1|id) + (1|trialxcondition), data = REWOD_PIT, REML = FALSE)
anova(main.n_gripsX)

# quick check with classical anova (! this is not reliable)
summary(aov(n_grips ~ cvalue1 + Error(id / (cvalue1)), data = REWOD_PIT))

# model comparison
main.n_grips.0X = lmer(n_grips ~ (1|id) + (1|trialxcondition), data = REWOD_PIT, REML = FALSE)
anova(main.n_grips.0X, main.n_gripsX, test = 'Chisq')

#sentence => main.n_grips is signifincatly better than the null model



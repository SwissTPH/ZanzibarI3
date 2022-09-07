# Copyright (C) 2022 Swiss Tropical and Public Health Institute
# Author: Aatreyee Mimi Das
# 
# This code is free software; you can redistribute it and/or modify it under the terms of version 2 of the GNU General Public License as published by the Free Software Foundation.
# 
# Plot probability of elimination over time for treatment of travellers and for reductions in mainland transmission

rm(list=ls())
library(ggplot2)
library(cowplot)
library(ggh4x)


setwd('~/data/')
params <- read.csv('RADZEC_data.csv')

setwd('~/interventions/')
int.tb <- read.csv('interventions.csv')

dir_figures <- '~/figures/'
save_figs <- TRUE

all_ints <- FALSE # when FALSE, we use the baseline interventions, when TRUE, we use cumulative addition of interventions (interventions are added on top of all previous interventions)

total_time <- 50 
############### Prob of elimination plot for importation treatment and beta reduction #################
setwd('~/processed_data/')

cD.df<- readRDS('D_zeroindig.rds')

scenarios <- dim(cD.df)[1]/(total_time-3)
cD.df$time <- rep(1:(total_time-3), times=scenarios)

##############################

if(all_ints){
  df.betared <- cD.df[cD.df$neighbours==100 & cD.df$followup==1 & cD.df$treatmentseeking==3
                      & cD.df$rfmda==1 & cD.df$location!='Mainland',]
} else {
  df.betared <- cD.df[cD.df$neighbours==0 & cD.df$followup==0.35 & cD.df$treatmentseeking==1
                      & cD.df$rfmda==0 & cD.df$location!='Mainland',]
}

time_int <- 7:(total_time-3) # Time for which the intervention was implemented (first 10 years is to allow to model to reach equilibrium, then we take cumulative indigenous cases over 3 years)
beta_red <- c(0.0, 0.5, 0.9) 
imptreat_red <- c(0.75, 0.9, 1)
df.betared <- df.betared[df.betared$betared %in% beta_red & df.betared$imptreat %in% imptreat_red &
                           df.betared$time %in% time_int,]
df.betared$time <- df.betared$time - 7

df.betared$imptreat <- replace(df.betared$imptreat, df.betared$imptreat==0.75, "75%" )
df.betared$imptreat <- replace(df.betared$imptreat, df.betared$imptreat==0.9, "90%" )
df.betared$imptreat <- replace(df.betared$imptreat, df.betared$imptreat==1.0, "100%" )

df.betared$imptreat <- factor(df.betared$imptreat, levels=c("75%","90%","100%"))


df.betared$betared <- replace(df.betared$betared, df.betared$betared==0.0, "No transmission\nreduction" )
df.betared$betared <- replace(df.betared$betared, df.betared$betared==0.5, "Transmission\nreduction = 50%" )
df.betared$betared <- replace(df.betared$betared, df.betared$betared==0.90, "Transmission\nreduction = 90%" )

cbp <- c("#622569", "#5b9aa0", "#b8a9c9", "#d6d4e0","#ee6b6e")

ggplot(df.betared, aes(x=time, y=prop_zero, colour=imptreat)) + geom_line() +
  labs(x="Time (years)", y="Probability of reaching elimination", colour='Proportion of infected travellers treated') +
  theme_minimal() +
  facet_nested(location~betared, nest_line = element_line(linetype = 1)) + panel_border(color = 'grey') +
  theme(panel.grid.minor = element_line(size = 0.2),panel.grid.major = element_line(size = 0.5),
        text = element_text(face="bold"),
        axis.title = element_text(size = 12), strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10), legend.position="bottom",
        legend.title=element_text(face="bold")) + 
  scale_y_continuous(limits = c(0,NA), labels = scales::percent_format()) +
  scale_color_manual(values=cbp)

if(save_figs){
  setwd(dir_figures)
  
  if(all_ints){
    filename1 <- '4imptreat_and_betared_allinterventions_prob_elim.pdf'
    filename2 <- '4imptreat_and_betared_allinterventions_prob_elim.png'
  } else {
    filename1 <- '4imptreat_and_betared_baseline_prob_elim.pdf'
    filename2 <- '4imptreat_and_betared_baseline_prob_elim.png'
  }
  ggsave(filename1, width=16, height=12, units="cm")
  ggsave(filename2, width=16, height=12, units="cm")
}


############### Prob of elimination plot for reductions in mainland transmission #################

setwd('~/interventions/')
int.tb <- read.csv('interventions_betared_mainland.csv')

setwd('~/processed_data/betared_mainland/')
df <- readRDS('D_zeroindig.rds')

if(all_ints){
  df.mainland <- df[df$neighbours==100 & df$followup==1 & df$treatmentseeking==3
                    & df$rfmda==1 & df$imptreat==0 & df$location!='Mainland',]
} else {
  df.mainland <- df[df$neighbours==0 & df$followup==0.35 & df$treatmentseeking==1
                    & df$rfmda==0 & df$imptreat==0 & df$location!='Mainland',]
}

time_int <- 7:(total_time-3)

df.mainland <- df.mainland[df.mainland$time %in% time_int,]
df.mainland$time <- df.mainland$time - 7

df.mainland$beta_red_mainland <- replace(df.mainland$beta_red_mainland, df.mainland$beta_red_mainland==0.1, "10%" )
df.mainland$beta_red_mainland <- replace(df.mainland$beta_red_mainland, df.mainland$beta_red_mainland==0.15, "15%" )
df.mainland$beta_red_mainland <- replace(df.mainland$beta_red_mainland, df.mainland$beta_red_mainland==0.2, "20%" )
df.mainland$beta_red_mainland <- replace(df.mainland$beta_red_mainland, df.mainland$beta_red_mainland==0.25, "25%" )
df.mainland$beta_red_mainland <- replace(df.mainland$beta_red_mainland, df.mainland$beta_red_mainland==0.3, "30%" )

df.mainland$beta_red_mainland <- factor(df.mainland$beta_red_mainland, levels=c("10%","15%","20%","25%","30%"))

df.mainland$beta_red_islands <- replace(df.mainland$beta_red_islands, df.mainland$beta_red_islands==0.0, 
                                        "TR on Zanzibar = 0%" )
df.mainland$beta_red_islands <- replace(df.mainland$beta_red_islands, df.mainland$beta_red_islands==0.5, 
                                        "TR on Zanzibar = 50%" )
df.mainland$beta_red_islands <- replace(df.mainland$beta_red_islands, df.mainland$beta_red_islands==0.9, 
                                        "TR on Zanzibar = 90%" )

df.mainland$beta_red_islands <- factor(df.mainland$beta_red_islands, levels=c("TR on Zanzibar = 0%",
                                                                              "TR on Zanzibar = 50%",
                                                                              "TR on Zanzibar = 90%"))


ggplot(df.mainland, aes(x=time, y=prop_zero, colour=beta_red_mainland)) + geom_line() +
  labs(x="Time (years)", y="Probability of reaching elimination", colour='Transmission reduction\non mainland Tanzania') +
  theme_minimal() + 
  facet_nested(location~beta_red_islands, nest_line = element_line(linetype = 1)) + panel_border(color = 'grey') +
  theme(panel.grid.minor = element_line(size = 0.2),panel.grid.major = element_line(size = 0.5),
        text = element_text(face="bold"),
        axis.title = element_text(size = 12), strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10), legend.position="bottom",
        legend.title=element_text(face="bold")) + 
  scale_y_continuous(limits = c(0,NA), labels = scales::percent_format()) +
  scale_color_manual(values=cbp)


setwd(dir_figures)

if(all_ints){
  filename1 <- '5betared_mainland_elimination_prob_allinterventions.pdf'
  filename2 <- '5betared_mainland_elimination_prob_allinterventions.png'
} else {
  filename1 <- '5betared_mainland_elimination_prob_baseline.pdf'
  filename2 <- '5betared_mainland_elimination_prob_baseline.png'
}
if(save_figs){
  ggsave(filename1, width=16, height=12, units="cm")
  ggsave(filename2, width=16, height=12, units="cm")
}

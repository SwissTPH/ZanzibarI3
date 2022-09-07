# Copyright (C) 2022 Swiss Tropical and Public Health Institute
# Author: Aatreyee Mimi Das
# 
# This code is free software; you can redistribute it and/or modify it under the terms of version 2 of the GNU General Public License as published by the Free Software Foundation.
# 
# Plot a barplot of baseline incidence of imported, introduced and indigenous cases

rm(list=ls()) # clear the environment
library(ggplot2)
library(ggpubr)

# Set the appropriate directory
dir_figures <- '~/figures/'
save_figs <- TRUE


setwd('~/processed_data')

P.df <- readRDS('P_incid_processed.rds')
T.df <- readRDS('T_incid_processed.rds')
D.df <- readRDS('D_incid_processed.rds')

P.df$case <- 'Imported'
T.df$case <- 'Introduced'
D.df$case <- 'Indigenous'

total_time <- 50 # Length of time for which the model was run

df <- rbind(P.df,T.df,D.df)
df$case <- factor(df$case, levels=c('Imported', 'Introduced', 'Indigenous'))
scenarios <- dim(df)[1]/total_time
df$time <- rep(1:total_time, times=scenarios)

df.baseline <- df[df$followup==0.35 & df$location!='Mainland' & df$neighbours==0 & df$treatmentseeking==1 & df$rfmda==0 & df$imptreat==0 &
                    df$betared==0 & df$time==total_time,]

# Total proportions of imported, introduced and indigenous cases on the islands and Zanzibar as a whole

props_Pemba <- df.baseline[df.baseline$location=='Pemba','median']/sum(df.baseline[df.baseline$location=='Pemba','median'])
props_Unguja <- df.baseline[df.baseline$location=='Unguja','median']/sum(df.baseline[df.baseline$location=='Unguja','median'])

total.infected <- sum(df.baseline$median)
props_ZNZ <- c(sum(df.baseline$median[1:2]), sum(df.baseline$median[3:4]), sum(df.baseline$median[5:6]))/total.infected

# Plotting
cbp1 <- c("#009E73","#0072B2", "#7800A8", "#56B4E9",
          "#F0E442","#9990999", "#D55E00", "#CC79A7")
cbp <- c("#d6d4e0", "#b8a9c9", "#622569", "#5b9aa0", "#f9d5e5")

ggplot(data=df.baseline, aes(x=location, y=median, group=case, fill=case)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8)+
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=.2,
                position=position_dodge(.8)) +
  geom_text(aes(label=sprintf("%0.0f", round(median, digits = 0))), colour='black', fontface='bold', vjust=-2.5,  position=position_dodge(.8))+
  labs(x='', y="Annual incidence of infections") +
  theme_pubr()+
  theme(axis.title = element_text(size=12, face = "bold"), strip.text.x = element_text(size = 10, face="bold"),
        strip.text.y = element_text(size = 10, face="bold"), 
        axis.text = element_text(size = 10, face="bold", colour = 'black'),
        legend.title=element_blank(), legend.text = element_text(size = 10)
  ) +
  scale_fill_manual(values=cbp) + scale_y_continuous(expand = c(0, 0), limits=c(0,13500))

if(save_figs){
  setwd(dir_figures)
  ggsave('0baseline_incidence_barplot_EB.pdf', width=16, height=9, units="cm")
}

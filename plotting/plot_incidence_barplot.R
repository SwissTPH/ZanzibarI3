# Copyright (C) 2022 Swiss Tropical and Public Health Institute
# Author: Aatreyee Mimi Das
# 
# This code is free software; you can redistribute it and/or modify it under the terms of version 2 of the GNU General Public License as published by the Free Software Foundation.
# 
# Plot a barplot of annual incidence of indigenous infections at baseline

rm(list=ls())
library(ggplot2)
library(cowplot)
library(ggh4x)
library(ggpubr)
library(gridExtra)


setwd('~/data/')
params <- read.csv('RADZEC_data.csv')

setwd('~/interventions/')
int.tb <- read.csv('interventions.csv')

# Set the appropriate directory
dir_figures <- '~/figures/'
save_figs <- TRUE

setwd('~/processed_data/')

D.df <- readRDS('D_incid_processed.rds')

df <- D.df

# Normalise to infections per 10,000 population
pop <- params$pop
df[which(df$location=='Pemba'),8:12] <- df[which(df$location=='Pemba'),8:12]*10000/pop[1]
df[which(df$location=='Unguja'),8:12] <- df[which(df$location=='Unguja'),8:12]*10000/pop[2]
df[which(df$location=='Mainland'),8:12] <- df[which(df$location=='Mainland'),8:12]*10000/pop[3]

total_time <- 50 # Relevant timepoint for calculating annual incidence
df$time <- 1:total_time 

cbp <- c("#622569", "#5b9aa0", "#b8a9c9", "#d6d4e0","#ee6b6e", "#8db39d") # Palette should be colourblind-friendly and B&W print-friendly
############# Select relevant rows #############

df <- df[df$time==total_time,]

# First figure is for each intervention implemented individually

baseline <- df[df$location!='Mainland' & df$followup==0.35 & df$neighbours==0 &
                 df$treatmentseeking==1 & df$rfmda==0 & df$imptreat==0 &
                 df$betared==0,]
baseline$scenario <- 'Baseline'

followup <- df[df$location!='Mainland' & df$followup==1 & df$neighbours==0 &
                 df$treatmentseeking==1 & df$rfmda==0 & df$imptreat==0 &
                 df$betared==0,]
followup$scenario <- '100% follow up'

neighbours <- df[df$location!='Mainland' & df$followup==0.35 & df$neighbours==100 &
                   df$treatmentseeking==1 & df$rfmda==0 & df$imptreat==0 &
                   df$betared==0,]
neighbours$scenario <- 'RCD inc. 100 neighbours'

treatseek <- df[df$location!='Mainland' & df$followup==0.35 & df$neighbours==0 &
                  df$treatmentseeking==3 & df$rfmda==0 & df$imptreat==0 &
                  df$betared==0,]
treatseek$scenario <- '3x treatment seeking rate'

rfmda <- df[df$location!='Mainland' & df$followup==0.35 & df$neighbours==0 &
              df$treatmentseeking==1 & df$rfmda==1 & df$imptreat==0 &
              df$betared==0,]
rfmda$scenario <- 'RDA'

imptreat <- df[df$location!='Mainland' & df$followup==0.35 & df$neighbours==0 &
                 df$treatmentseeking==1 & df$rfmda==0 & df$imptreat==0.9 &
                 df$betared==0,]
imptreat$scenario <- '90% treatment of travellers'

df.single <- rbind(baseline,followup,neighbours,treatseek,rfmda,imptreat)

df.single$scenario <- factor(df.single$scenario, levels=c('Baseline',
                                                          '100% follow up',
                                                          'RCD inc. 100 neighbours',
                                                          'RDA',
                                                          '3x treatment seeking rate',
                                                          '90% treatment of travellers'))


p1 <- ggplot(df.single, aes(x=location, y=median, fill=scenario)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8) +
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=.2,
                position=position_dodge(.8)) +
  labs(x='', y="Incidence of indigenous\ninfections per 10,000 population") +
  theme_classic()+
  theme(axis.title = element_text(face = "bold"), strip.text.x = element_text(size = 10, face="bold"),
        strip.text.y = element_text(size = 10, face="bold"), 
        axis.text = element_text(size = 10, face="bold", colour = 'black'),
        legend.position = "top",
        legend.title=element_blank(), legend.text = element_text(size = 10)
  ) +
  scale_fill_manual(values=cbp) + scale_y_continuous(expand = c(0, 0), limits=c(0,NA))


# Second figure is for each intervention being added cumulatively on top of last intervention
neighbours <- df[df$location!='Mainland' & df$followup==1 & df$neighbours==100 &
                   df$treatmentseeking==1 & df$rfmda==0 & df$imptreat==0 &
                   df$betared==0,]
neighbours$scenario <- '+ 100 neighbours'

rfmda <- df[df$location!='Mainland' & df$followup==1 & df$neighbours==100 &
              df$treatmentseeking==1 & df$rfmda==1 & df$imptreat==0 &
              df$betared==0,]
rfmda$scenario <- '+ RDA'

treatseek <- df[df$location!='Mainland' & df$followup==1 & df$neighbours==100 &
                  df$treatmentseeking==3 & df$rfmda==1 & df$imptreat==0 &
                  df$betared==0,]
treatseek$scenario <- '+ 3x treatment seeking rate'



imptreat <- df[df$location!='Mainland' & df$followup==1 & df$neighbours==100 &
                 df$treatmentseeking==3 & df$rfmda==1 & df$imptreat==0.9 &
                 df$betared==0,]
imptreat$scenario <- '+ 90% treatment of travellers'

df.combined <- rbind(baseline,followup,neighbours,treatseek,rfmda,imptreat)

df.combined$scenario <- factor(df.combined$scenario, levels=c('Baseline',
                                                              '100% follow up',
                                                              '+ 100 neighbours',
                                                              '+ RDA',
                                                              '+ 3x treatment seeking rate',
                                                              '+ 90% treatment of travellers'))

p2 <- ggplot(df.combined, aes(x=location, y=median, fill=scenario)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8) +
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=.2,
                position=position_dodge(.8)) +
  labs(x='', y="Incidence of indigenous\ninfections per 10,000 population") +
  theme_classic()+
  theme(axis.title = element_text(face = "bold"), strip.text.x = element_text(size = 10, face="bold"),
        strip.text.y = element_text(size = 10, face="bold"), 
        axis.text = element_text(size = 10, face="bold", colour = 'black'),
        legend.position = "top",
        legend.title=element_blank(), legend.text = element_text(size = 10)
  ) +
  scale_fill_manual(values=cbp) + scale_y_continuous(expand = c(0, 0), limits=c(0,NA))


both_plot <- plot_grid(p1, p2, labels=c("a", "b"), ncol = 1, nrow = 2)

if(save_figs){
  setwd(dir_figures)
  filename1 <- 'Incidence_barplot.pdf'
  filename2 <- 'Incidence_barplot.png'
  ggsave(filename1, width=16, height=16, units="cm")
  ggsave(filename2, width=16, height=16, units="cm")
}

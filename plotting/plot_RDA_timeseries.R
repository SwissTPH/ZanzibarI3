# Copyright (C) 2022 Swiss Tropical and Public Health Institute
# Author: Aatreyee Mimi Das
# 
# This code is free software; you can redistribute it and/or modify it under the terms of version 2 of the GNU General Public License as published by the Free Software Foundation.
# 
# Plot a timeseries plot of switching from RCD to RDA

rm(list=ls())
library(ggplot2)
library(ggh4x)
library(cowplot)

setwd('~/data/')
params <- read.csv('RADZEC_data.csv')

setwd('~/interventions/')
interventions.tb <- read.csv('interventions.csv')

dir_figures <- '~/figures/'
save_figs <- TRUE

all_ints <- FALSE # when FALSE, we use the baseline interventions, when TRUE, we use cumulative addition of interventions (interventions are added on top of all previous interventions)

total_time <- 50 

outputsdir <- '~/outputs/i3_model/'
processeddir <- '~/processed_data/'

setwd(processeddir)

P.df <- readRDS('P_incid_processed.rds')
T.df <- readRDS('T_incid_processed.rds')
D.df <- readRDS('D_incid_processed.rds')

P.df$case <- 'Imported'
T.df$case <- 'Introduced'
D.df$case <- 'Indigenous'

df <- rbind(P.df,T.df,D.df)
df$case <- factor(df$case, levels=c('Imported', 'Introduced', 'Indigenous'))
scenarios <- dim(df)[1]/total_time
df$time <- rep(1:total_time, times=scenarios)

df.zerofollowup <- df[df$followup==0 & df$neighbours==0 & df$rfmda==0 & df$imptreat==0,]
df.zerofollowup$rfmda <- 1

df <- rbind(df, df.zerofollowup) #to have a 0% follow up line for rfMDA as well


############## RDA timeseries plot ###############


interventions.tb <- cbind(interventions.tb, 1:nrow(interventions.tb))

run_row <- interventions.tb[interventions.tb$followup==0.35 &
                              interventions.tb$neighbours==0 &
                              interventions.tb$treatmentseeking==1 &
                              interventions.tb$rfmda==1 &
                              interventions.tb$imptreat==0 &
                              interventions.tb$betared==0,ncol(interventions.tb)]

run_number <- run_row[1]

setwd(outputsdir)
files_list <- c('cP_out_','cT_out_', 'cD_out_')
filename <- paste0(files_list,run_number,'.rds')
loc <- c('Pemba', 'Unguja')
case <- c('Imported', 'Introduced', 'Indigenous')

no_runs <-2

initruns <- list()

for (jj in 1:3){
  for(ii in 1:2){
    out <- readRDS(filename[jj])
    out[1:49,] <- out[2:50,]-out[1:49,]
    out[51:99,] <- out[52:100,]-out[51:99,]
    x <- data.frame(run1 = out[1:50,ii], run2 = out[51:100,ii])
    x$time <- 1:50
    x$location <- loc[ii]
    x$case <- case[jj]
    if(ii==1){
      initruns[[jj]] <- x
    }
    else{
      initruns[[jj]] <- rbind(initruns[[jj]],x)
    }
  }
}

df.runs <- rbind(initruns[[1]], initruns[[2]], initruns[[3]])

time.int <- 10:25

df.runs <- df.runs[df.runs$time %in% time.int,]
df.runs$time <- df.runs$time - 10

df.rda <- df[df$neighbours==0 & df$imptreat==0 & df$treatmentseeking==1 & 
               df$followup==0.35 & df$rfmda==1 &
               df$betared==0 & df$location!='Mainland',]

df.rda$case <- factor(df.rda$case, levels=c('Imported',
                                            'Introduced',
                                            'Indigenous'))

df.runs$case <- factor(df.runs$case, levels=c('Imported',
                                            'Introduced',
                                            'Indigenous'))

df.rda <- df.rda[df.rda$time %in% time.int,]
df.rda$time <- df.rda$time - 10
colnames(df.rda)[9:12] <- c('CI2.5', 'CI97.5', 'CI25', 'CI75')

ggplot(df.rda, aes(x=time, y=median)) + 
  facet_nested(location+case~., nest_line = element_line(linetype = 1), scales = "free_y") +
  geom_ribbon(aes(ymin=CI25, ymax=CI75, fill = '50% PI')) +
  geom_line(data=df.runs, aes(x=time, y=run1, colour='Example individual runs'), size=0.1)+
  geom_line(data=df.runs, aes(x=time, y=run2, colour='Example individual runs'), size=0.1)+
  geom_line(aes(x=time, y=median, colour='Median'), size=0.5) + 
  geom_line(aes(x=time, y=CI2.5, colour='95% PI'), size=0.5) +
  geom_line(aes(x=time, y=CI97.5, colour='95% PI'), size=0.5) +
  theme_minimal() + panel_border(color = 'grey50') +
  labs(x="Time (years)", y="Annual incidence of infections") + 
  theme(panel.grid.minor = element_blank(),
        text = element_text(face="bold"),
        axis.title = element_text(size=10), strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10), legend.position="bottom",
        legend.title=element_blank(), legend.text = element_text(size=8),
        legend.spacing = unit(-0.4, 'cm'))+
  scale_colour_manual(values=c('95% PI'='grey40', 'Median'='black', 'Example individual runs'='#c0a7c3'))+
  scale_fill_manual(values=c('50% PI'='grey80')) +
  guides(color = guide_legend(order = 2),fill = guide_legend(order = 1))

setwd(dir_figures)
filename1 <- 'RDA_uncertainty_timeseries.pdf'
filename2 <- 'RDA_uncertainty_timeseries.png'
ggsave(filename1, width=12, height=16, units="cm")
ggsave(filename2, width=12, height=16, units="cm")

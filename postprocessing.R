# Copyright (C) 2022 Swiss Tropical and Public Health Institute
# Author: Aatreyee Mimi Das
# 
# This code is free software; you can redistribute it and/or modify it under the terms of version 2 of the GNU General Public License as published by the Free Software Foundation.
# 
# Processing outputs to have the median value, 95% range and 50% range at each timepoint for each run
# and the probability of 0 indigenous cases across all 500 runs

#######FUNCTIONS########################
# For reading in files
read.in <- function(x, run_number){
  filename = paste0(x,'_out_',run_number,'.rds')
  out = readRDS(filename)
  return(out)
}

# For processing for timeseries plots
process.data <- function(x,total_time){
  no.rows <- total_time
  x_reshape <- matrix(x, nrow = no.rows, byrow=FALSE)
  print(dim(x_reshape))
  x_reshape_med <- apply(x_reshape, 1, median)
  x_reshape_025 <- apply(x_reshape, 1, function(x) quantile(x, 0.025))
  x_reshape_25 <- apply(x_reshape, 1, function(x) quantile(x, 0.25))
  x_reshape_75 <- apply(x_reshape, 1, function(x) quantile(x, 0.75))
  x_reshape_975 <- apply(x_reshape, 1, function(x) quantile(x, 0.975))
  
  x_reshape_final <- cbind(x_reshape_med, x_reshape_025, x_reshape_975,
                           x_reshape_25,x_reshape_75)
  return(x_reshape_final)
}

process.incid.data <- function(x,total_time){
  no.rows <- total_time
  x_reshape <- matrix(x, nrow = no.rows, byrow=FALSE)
  print(dim(x_reshape))
  x_reshape <- rbind(x_reshape[1,], x_reshape[2:total_time,]-x_reshape[1:(total_time-1),])
  x_reshape_med <- apply(x_reshape, 1, median)
  x_reshape_025 <- apply(x_reshape, 1, function(x) quantile(x, 0.025))
  x_reshape_25 <- apply(x_reshape, 1, function(x) quantile(x, 0.25))
  x_reshape_75 <- apply(x_reshape, 1, function(x) quantile(x, 0.75))
  x_reshape_975 <- apply(x_reshape, 1, function(x) quantile(x, 0.975))
  
  x_reshape_final <- cbind(x_reshape_med, x_reshape_025, x_reshape_975,
                           x_reshape_25,x_reshape_75)
  return(x_reshape_final)
}

# For probability of 0 indigenous cases
process.zero.cases <- function(x, total_time){
  no.rows <- total_time
  x_reshape <- matrix(x, nrow = no.rows, byrow=FALSE)
  x_diff <- x_reshape[4:total_time,]-x_reshape[1:(total_time-3),]
  zero.indig <- rowSums(x_diff==0)/dim(x_diff)[2]
  return(zero.indig)
}

# For probability of 0 indigenous cases - optomistic version where after 3 years of 0 cases, the area is certified and does not lose certification even with some indigenous cases
process.zero.cases.optimistic <- function(x,total_time){
  no.rows <- total_time
  x_reshape <- matrix(x, nrow = no.rows, byrow=FALSE)
  x_diff <- x_reshape[4:total_time,]-x_reshape[1:(total_time-3),]
  x_diff_optimistic<-x_diff
  for(kk in 1:(dim(x_diff)[2])){
    if(length(which(x_diff[,kk]==0))!=0){
      elim_ind <- min(which(x_diff[,kk]==0))
      x_diff_optimistic[elim_ind:(dim(x_diff)[1]),kk] <- 0
    }
  }
  
  zero.indig <- rowSums(x_diff_optimistic==0)/dim(x_diff_optimistic)[2]
  return(zero.indig)
}

# Directories need to be set up according to which outputs are being processed and therefore which interventions file was used
setwd('~/interventions')
interventions.tb <- read.csv('interventions.csv')
intervention_names <- colnames(interventions.tb)
loc = c('Pemba', 'Unguja', 'Mainland')

outputsdir <- '~/outputs/i3_model/'
processeddir <- '~/processed_data/'

timeseries = TRUE # postprocessing for timeseries plots
incidence = TRUE # postprocessing for incidence plots
zero_cases = TRUE # postprocessing for elimination plots

# Note, the following three Booleans should be set to false if considering any runs outside of the base 'interventions.csv'
zero_cases_optimistic = TRUE # postprocessing for optimistic elimination plots
recovery = TRUE  # postprocessing for cumulative sum of recovered cases
rdts = TRUE  # postprocessing for RDTs and ACTs used

CA <- commandArgs(TRUE)
ifelse(length(CA)==0, UseScript <- FALSE, UseScript <- TRUE)
ifelse (UseScript, jobID <- as.integer(CA)[1], jobID <- 1)

i3_list <- c('P','T','D','I')
category_inf <- i3_list[jobID]
print(list('jobID: ',jobID))

if(category_inf=='I'){
  incidence <- FALSE  # I don't have outputs for cumulative I (and this can easily be gotten from summing cP, cT and cD)
  recovery <- FALSE
}

if(category_inf!='D'){
  zero_cases <- FALSE  # to avoid running 4 times unnecessarily. zero_cases only refers to 0 indigenous cases.
  zero_cases_optimistic <- FALSE
}

if(category_inf!='P'){
  rdts <- FALSE
}

total_time = 50
#####Median values for timeseries/final prevalence plots######
if(timeseries){  
  setwd(outputsdir)
  processed_reg = list() # for processing each region separately

  
  store_processed = list() # for storing each processed data frame for concatenation at the end

  # Considering that the runs come in pairs, with each consecutive pair being 250 runs for the same scenario
  pairs = seq(from=1, to=dim(interventions.tb)[1], by=2)
  for (ii in pairs){
    print(ii)
    print(category_inf)
    P_out1 = read.in(category_inf,ii-1)
    P_out2 = read.in(category_inf,ii)
    P_out = rbind(P_out1,P_out2)
    
    for (jj in 1:3){
      print(jj)
      processed_reg[[jj]] = process.data(P_out[,jj], total_time)
      interventions = interventions.tb[ii,]
      rownames(interventions) <- NULL
      location = loc[jj]
      processed_reg[[jj]] = cbind(location, interventions, processed_reg[[jj]])
      colnames(processed_reg[[jj]]) <- c('location',intervention_names,
                                         'median', '2.5%','97.5%','25%','75%')
    }
    store_processed[[ceiling(ii/2)]] = do.call(rbind.data.frame,processed_reg) #ceiling() used to turn 'pairs' list into a list of consecutive numbers
    rm(P_out1,P_out2,P_out)
    
  }
  
  store_output = do.call(rbind.data.frame, store_processed)
  
  setwd(processeddir)
  saveRDS(store_output, file=paste0(category_inf,'_processed.rds'))
}


#####Median values for incidence plots######
if(incidence){  
  setwd(outputsdir)

  processed_reg = list() # for processing each region separately
  
  store_processed = list() # for storing each processed data frame for concatenation at the end
  
  category_inf_incid <- paste0('c',category_inf)
  
  pairs = seq(from=1, to=dim(interventions.tb)[1], by=2)
  for (ii in pairs){
    print(ii)
    print(category_inf)
    P_out1 = read.in(category_inf_incid,ii-1)
    P_out2 = read.in(category_inf_incid,ii)
    P_out = rbind(P_out1,P_out2)
    
    for (jj in 1:3){
      processed_reg[[jj]] = process.incid.data(P_out[,jj],total_time)
      interventions = interventions.tb[ii,]
      rownames(interventions) <- NULL
      location = loc[jj]
      processed_reg[[jj]] = cbind(location, interventions, processed_reg[[jj]])
      colnames(processed_reg[[jj]]) <- c('location',intervention_names,
                                         'median', '2.5%','97.5%','25%','75%')
    }
    store_processed[[ceiling(ii/2)]] = do.call(rbind.data.frame,processed_reg) #ceiling() used to turn 'pairs' list into a list of consecutive numbers
    rm(P_out1,P_out2,P_out)
  }
  
  store_output = do.call(rbind.data.frame, store_processed)

  setwd(processeddir)
  saveRDS(store_output, file=paste0(category_inf,'_incid_processed.rds'))

}


#####Probability of 0 indigenous cases###################

if(zero_cases){
  setwd(outputsdir)
  zeroindig_D = list()
  D_zeroindig = list()
  
  pairs = seq(from=1, to=dim(interventions.tb)[1], by=2)
  for (ii in pairs){
    print(ii)
    cD_out1 = read.in('cD',ii-1)
    cD_out2 = read.in('cD',ii)
    cD_out = rbind(cD_out1,cD_out2)
    
    for (jj in 1:3){
      zeroindig_D[[jj]] = process.zero.cases(cD_out[,jj], total_time)
      interventions = interventions.tb[ii,]
      rownames(interventions) <- NULL
      location = loc[jj]
      zeroindig_D[[jj]] = cbind(location, interventions, zeroindig_D[[jj]])
    }
    D_zeroindig[[ceiling(ii/2)]] = do.call(rbind.data.frame,zeroindig_D) #ceiling() used to turn 'pairs' list into a list of consecutive numbers
    colnames(D_zeroindig[[ceiling(ii/2)]]) <- c('location',intervention_names,
                                                'prop_zero')
    rm(cD_out1,cD_out2,cD_out)
  }
  D_zerocases_output = do.call(rbind.data.frame,D_zeroindig)
  D_zerocases_output$time <- 1:47
  
  setwd(processeddir)
  saveRDS(D_zerocases_output, file='D_zeroindig.rds')
}


#####Optimistic probability of 0 indigenous cases info###################

if(zero_cases_optimistic){
  setwd(outputsdir)
  zeroindig_D = list()
  D_zeroindig = list()
  
  pairs = seq(from=1, to=dim(interventions.tb)[1], by=2)
  for (ii in pairs){
    print(ii)
    cD_out1 = read.in('cD',ii-1)
    cD_out2 = read.in('cD',ii)
    cD_out = rbind(cD_out1,cD_out2)
    
    for (jj in 1:3){
      zeroindig_D[[jj]] = process.zero.cases.optimistic(cD_out[,jj], total_time)
      interventions = interventions.tb[ii,]
      rownames(interventions) <- NULL
      location = loc[jj]
      zeroindig_D[[jj]] = cbind(location, interventions, zeroindig_D[[jj]])
    }
    D_zeroindig[[ceiling(ii/2)]] = do.call(rbind.data.frame,zeroindig_D) #ceiling() used to turn 'pairs' list into a list of consecutive numbers
    colnames(D_zeroindig[[ceiling(ii/2)]]) <- c('location',intervention_names,
                                                'prop_zero')
    rm(cD_out1,cD_out2,cD_out)
  }
  D_zerocases_output = do.call(rbind.data.frame,D_zeroindig)
  D_zerocases_output$time <- 1:47
  
  setwd(processeddir)
  saveRDS(D_zerocases_output, file='D_zeroindig_optimistic.rds')
}

############# Median recovered cases and RDT and ACT use #####################3

if(recovery){  
  setwd(outputsdir)
  processed_reg = list() # for processing each region separately
  
  store_processed = list() # for storing each processed data frame for concatenation at the end

  category_inf_rec <- paste0('r',category_inf)
  
  # Considering that the runs come in pairs, with each consecutive pair being 250 runs for the same scenario
  pairs = seq(from=1, to=dim(interventions.tb)[1], by=2)
  for (ii in pairs){
    print(ii)
    print(category_inf)
    P_out1 = read.in(category_inf_rec,ii-1)
    P_out2 = read.in(category_inf_rec,ii)
    P_out = rbind(P_out1,P_out2)
    
    for (jj in 1:3){
      processed_reg[[jj]] = process.incid.data(P_out[,jj], total_time)
      interventions = interventions.tb[ii,]
      rownames(interventions) <- NULL
      location = loc[jj]
      processed_reg[[jj]] = cbind(location, interventions, processed_reg[[jj]])
      colnames(processed_reg[[jj]]) <- c('location',intervention_names,
                                         'median', '2.5%','97.5%','25%','75%')
    }
    store_processed[[ceiling(ii/2)]] = do.call(rbind.data.frame,processed_reg) #ceiling() used to turn 'pairs' list into a list of consecutive numbers
    rm(P_out1,P_out2,P_out)
  }
  
  store_output = do.call(rbind.data.frame, store_processed)
  
  setwd(processeddir)
  saveRDS(store_output, file=paste0(category_inf,'_rec_processed.rds'))
  
}

if(recovery){  
  setwd(outputsdir)
  processed_reg = list() # for processing each region separately
  
  store_processed = list() # for storing each processed data frame for concatenation at the end
  
  category_inf_rec <- paste0('r',category_inf)
  
  # considering that the runs come in pairs, with each consecutive pair being 250 runs for the same scenario
  pairs = seq(from=1, to=dim(interventions.tb)[1], by=2)
  for (ii in pairs){
    print(ii)
    print(category_inf)
    P_out1 = read.in(category_inf_rec,ii-1)
    P_out2 = read.in(category_inf_rec,ii)
    P_out = rbind(P_out1,P_out2)
    
    for (jj in 1:3){
      processed_reg[[jj]] = process.incid.data(P_out[,jj],total_time)
      interventions = interventions.tb[ii,]
      rownames(interventions) <- NULL
      location = loc[jj]
      processed_reg[[jj]] = cbind(location, interventions, processed_reg[[jj]])
      colnames(processed_reg[[jj]]) <- c('location',intervention_names,
                                         'median', '2.5%','97.5%','25%','75%')
    }
    store_processed[[ceiling(ii/2)]] = do.call(rbind.data.frame,processed_reg) # ceiling() used to turn 'pairs' list into a list of consecutive numbers
    rm(P_out1,P_out2,P_out)
  }
  
  store_output = do.call(rbind.data.frame, store_processed)
  
  setwd(processeddir)
  saveRDS(store_output, file=paste0(category_inf,'_rec_processed.rds'))
  
}

if(rdts){  
  setwd(outputsdir)
  processed_reg = list() # for processing each region separately
  
  store_processed_rdts = list() # for storing each processed data frame for concatenation at the end
  store_processed_acts = list() # for storing each processed data frame for concatenation at the end
  
  rdts <- 'rdts'
  acts <- 'acts'
  
  # Considering that the runs come in pairs, with each consecutive pair being 250 runs for the same scenario
  pairs = seq(from=1, to=dim(interventions.tb)[1], by=2)
  for (ii in pairs){
    print(ii)
    
    # RDTs
    print(category_inf)
    P_out1 = read.in(rdts,ii-1)
    P_out2 = read.in(rdts,ii)
    P_out = rbind(P_out1,P_out2)
    
    for (jj in 1:3){
      processed_reg[[jj]] = process.incid.data(P_out[,jj], total_time)
      interventions = interventions.tb[ii,]
      rownames(interventions) <- NULL
      location = loc[jj]
      processed_reg[[jj]] = cbind(location, interventions, processed_reg[[jj]])
      colnames(processed_reg[[jj]]) <- c('location',intervention_names,
                                         'median', '2.5%','97.5%','25%','75%')
    }
    store_processed_rdts[[ceiling(ii/2)]] = do.call(rbind.data.frame,processed_reg) #ceiling() used to turn 'pairs' list into a list of consecutive numbers
    rm(P_out1,P_out2,P_out)
    
    # ACTs
    print(category_inf)
    P_out1 = read.in(acts,ii-1)
    P_out2 = read.in(acts,ii)
    P_out = rbind(P_out1,P_out2)
    
    for (jj in 1:3){
      processed_reg[[jj]] = process.incid.data(P_out[,jj],total_time)
      interventions = interventions.tb[ii,]
      rownames(interventions) <- NULL
      location = loc[jj]
      processed_reg[[jj]] = cbind(location, interventions, processed_reg[[jj]])
      colnames(processed_reg[[jj]]) <- c('location',intervention_names,
                                         'median', '2.5%','97.5%','25%','75%')
    }
    store_processed_acts[[ceiling(ii/2)]] = do.call(rbind.data.frame,processed_reg) # ceiling() used to turn 'pairs' list into a list of consecutive numbers
    rm(P_out1,P_out2,P_out)
  }
  
  store_output_rdts = do.call(rbind.data.frame, store_processed_rdts)
  store_output_acts = do.call(rbind.data.frame, store_processed_acts)
  
  setwd(processeddir)
  saveRDS(store_output_rdts, file=paste0(rdts,'_rec_processed.rds'))
  saveRDS(store_output_acts, file=paste0(acts,'_rec_processed.rds'))
  
}

rm(list=ls())
gc()
source('helpers.R')
source('helpers_num_phase_diag.R')
library(dplyr)
library(parallel)
# This dataset contains the average values for the coefficients at each scenario.
load('survivor_prob_f2.5resol.RData') 
# Run in parallel for faster performance
n_co <- max(detectCores() - 8,8)
print(paste('Using',n_co,'cores!'))

# Resolution used for the plots. 
# It represents the number of values that are considered between maximum and minium values of the coefficients
Ns <- 40 
# Filter dataset to only retain highest demand, and calculate maximum and minimum values of the parameters 
# at each flooding scenario
sp <- survivor_prob %>%
  filter(demand_level == 100) %>%
  group_by(year_range,RCP,prec) %>%
  summarise(
    p = mean(p),mp = min(p),MP = max(p),
    s12 = mean(s12),ms12 = min(s12),Ms12 = max(s12),
    s13 = mean(s13),ms13 = min(s13),Ms13 = max(s13),
    s23 = mean(s23),ms23 = min(s23),Ms23 = max(s23),
    a12 = mean(a12),ma12 = min(a12),Ma12 = max(a12),
    a23 = mean(a23),ma23 = min(a23),Ma23 = max(a23))

bdir <- paste('num_phase_diag',Sys.Date(),sep='')
dir.create(bdir)

lapply(1:nrow(sp),function(row){
  print(paste('Working with scenario',sp[row,]$year_range,sp[row,]$RCP,sp[row,]$prec,'(number',row,')'))
  ##  Setup: generate a dataset with Ns values for each coefficient, between minimum and maximum values.
  ## The function automatically makes four datasets, s12 v s23, s12 v s13, s23 v 13 and a12 v a23. 
  # For the purpose of the work, only s12 v s13 is used
  coefs <- gen_coef_from_df(row,sp,Ns = Ns)
  ### s12 v s23
  print('Starting with s12 v s23')
  save_path <- sub('%','p',paste(bdir,paste(paste(sp[row,]$year_range,sp[row,]$RCP,sp[row,]$prec,'s12_s23',sep='_'),'.RData',sep=''),sep='/'))

  ## Parallel setup
  cl <- makeCluster(n_co,type = 'FORK')

  U <- parLapply(cl,1:nrow(coefs[['s12_s23']]),function(row_c){
    U <- evol_to_stable_state(row_c,coefs[['s12_s23']])
    return(U)
  }) %>%
    bind_rows
  save(U,file=save_path)
  print('Finished with s12 v s23')
  
  ### s12 v s13
  save_path <- sub('%','p',paste(bdir,paste(paste(sp[row,]$year_range,sp[row,]$RCP,sp[row,]$prec,'s12_s13',sep='_'),'.RData',sep=''),sep='/'))
  
  print('Starting with s12 v s13')
  U <- parLapply(cl,1:nrow(coefs[['s12_s13']]),function(row_c){
    U <- evol_to_stable_state(row_c,coefs[['s12_s13']])
    return(U)
  }) %>%
    bind_rows
  save(U,file=save_path)
  print('Finishing with s12 v s13')
  
  stopCluster(cl)
  print(paste('Concluded with scenario',sp[row,]$year_range,sp[row,]$RCP,sp[row,]$prec))
  
})


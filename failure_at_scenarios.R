rm(list=ls())
gc()
source('helpers.R')
source('helpers_num_phase_diag.R')
library(dplyr)
library(parallel)

load('survivor_prob_f2.5resol.RData')

# Parallelize to run faster
n_co <- 6
only_functionals <- TRUE
print(paste('Using',n_co,'cores!'))

sp <- survivor_prob %>%
  filter(demand_level == 100) 
# Dataset with initial conditions at each scenario, based on the stable values of U
load('inicon_for_scenarios.RData')

# Directory for saving trajectories generated
bdir <- 'failure_at_scenarios2023-08-04'
dir.create(bdir)


lapply(1:nrow(sp),function(row){
  print(paste('Working with scenario',sp[row,]$scenario,'(number',row,')'))
  ##  Setup
  coefs <- sp[row,] %>% ungroup %>% select(p,s12,s13,s23,a12,a23) %>% mutate(d = p /(a12*a23)) %>% as.matrix()  %>% {.[1,]}
  
  if( only_functionals & (coefs['p'] > coefs['s13'] + min(coefs['s12'],coefs['a12']*coefs['s23'])) ){
    print('Oh! This system is non-functional, skipping...')
  }else{
    save_path <- paste(bdir,paste('failure_at_',sp[row,]$scenario,'.RData',sep=''),sep='/')
    if(file.exists(save_path)){
      print(paste(save_path, 'has previously been done!'))
    }else{
      print(paste('Saving @',save_path))
      
      # Initial conditions to use:
      U0s <- init_condis[[sp[row,]$scenario]] %>%
        select(Y1,Y2,Y3) %>%
        as.matrix
      ## Parallel setup
      
      print('Starting!')
      # row_u0 <- nrow(U0s)
      time_ranges <- seq(0.25,3,.25)
      failure_start <- 0
      weight_matrix <- construct_weight_matrix_from_coefs_3x3(coefs)
      times_parms <- list(length = 1e4,tMax=15,tMin=0)         
      # i_t <- 1
      U <- lapply(1:length(time_ranges),function(i_t){
        print(paste('Working with time_lapse =',time_ranges[i_t]))
        parms_ap <- prepare_parms_approx_with_failure_period(
          weight_matrix = weight_matrix,
          dem_coef = c(0,0,coefs['d']), 
          prod_coef = c(coefs['p'],0,0),
          failure_timerange = c( failure_start, failure_start + time_ranges[i_t]),
          epsilon = 1e-6
        )
        cl <- makeCluster(n_co,type = 'FORK')
        U <- parLapply(cl,1:nrow(U0s),function(row_u0){
          U0 <- U0s[row_u0,]
          
          U <- simulate_exact_system(
            U0 = U0, dU = dU_exact, parms = parms_ap,
            times_parms = times_parms
          ) %>%
            as.data.frame %>%
            mutate(row_u0 = row_u0) 
          
          counter <- 1
          while(any(U$Y1 < 0 | U$Y1 > 1) |
                any(U$Y2 < 0 | U$Y2 > 1) |
                any(U$Y3 < 0 | U$Y3 > 1)){
            if(counter %% 10 == 0){
              print("Ups! Let's try again")
              print(paste('Trial:',counter))
            }
            counter <- counter + 1
            times_parms$length <- as.integer(1.1*times_parms$length)
            # Simulate system again
            U <- simulate_exact_system(
              U0 = U0, dU = dU_exact, parms = parms_ap,
              times_parms = times_parms
            ) %>%
              as.data.frame %>%
              mutate(row_u0 = row_u0)
          }
          ### Now calculate statistics
          U %>%
            mutate(
              D = (1+1e-6)*Y3 / (Y3 + 1e-6),
              P = ifelse(time < failure_start + time_ranges[i_t],0,(1e-6+1)*(1-Y1) / (1-Y1 + 1e-6))
            ) %>%
            summarise(
              qY1 = 0.5 * sum( (Y1 + lag(Y1)) * 
                                 ( time-lag(time) ), na.rm = TRUE ) / 
                ( max(time)-min(time) ),
              qY2 = 0.5 * sum( (Y2 + lag(Y2)) * 
                                 ( time-lag(time) ), na.rm = TRUE ) /
                ( max(time)-min(time) ),
              qY3 = 0.5 * sum( (Y3 + lag(Y3)) * 
                                 ( time-lag(time) ), na.rm = TRUE ) /
                ( max(time)-min(time) ),
              qD = 0.5 * sum( (D + lag(D)) * 
                                ( time-lag(time) ), na.rm = TRUE ) /
                ( max(time)-min(time) ),
              qP = 0.5 * sum( (P + lag(P)) * 
                                ( time-lag(time) ), na.rm = TRUE ) /
                ( max(time)-min(time) ),
              negative_fail = any(Y1 > 1 | Y1 < 0 | Y2 > 1 | Y2 < 0 | Y3 > 1 | Y3 < 0),
              tfail1e5 = first(time[which(D < (1-1e-5)*max(D))]),
              tfail99 = first(time[which(D < 0.99*max(D))]),
              tfail95 = first(time[which(D < 0.95*max(D))]),
              tfail90 = first(time[which(D < 0.90*max(D))]),
              tfail70 = first(time[which(D < 0.7*max(D))]),
              q0Y1 = first(Y1),
              q0Y2 = first(Y2),
              q0Y3 = first(Y3),
              stY1 = last(Y1),
              stY2 = last(Y2),
              stY3 = last(Y3),
              q0D = first(D),
              time_lapse = max(time)-min(time),
              .groups='drop') %>%
            mutate(qR = qY1+qY2+qY3) %>%
            mutate(qR = qY1+qY2+qY3) %>%
            mutate(U0 = sum(U0s[row_u0,]),
                   Y10 = U0s[row_u0,'Y1'],
                   Y20 = U0s[row_u0,'Y2'],
                   Y30 = U0s[row_u0,'Y3'],
                   row_u0 = row_u0) %>%
            ungroup() %>%
            mutate(
              p = coefs['p'],
              d = coefs['d'],
              a12 = coefs['a12'],
              a23 = coefs['a23'],
              s12 = coefs['s12'],
              s13 = coefs['s13'],
              s23 = coefs['s23']) %>%
            return()
          
        }) %>%
          bind_rows() %>%
          mutate(failure_duration = time_ranges[i_t])
        
        stopCluster(cl)
        print('Finished! Next time_lapse...')
        return(U)
      }) %>%
        bind_rows
      save(U,file=save_path)
      
      print(paste('Concluded with scenario',sp[row,]$scenario))
      print('################################')
    }
  }
})

rm(list=ls())
gc()
source('helpers.R')
source('helpers_num_phase_diag.R')
library(dplyr)
library(parallel)

load('survivor_prob_f2.5resol.RData')
Ns <- 5
Ns2 <- 10
n_co <- 6
print(paste('Using',n_co,'cores!'))

sp <- survivor_prob %>%
  filter(demand_level == 100) %>%
  filter(scenario == '0')

bdir <- 'final_state_u0curve_for_demo2023-08-08'
dir.create(bdir)

U0s <- expand.grid(Y1 = seq(0,1,length.out = Ns),
                   Y2 = seq(0,1,length.out = Ns),
                   Y3 = seq(0,1,length.out = Ns)) %>%
  mutate(Y1 = case_when(
    Y1 == 0 ~ Y1+1e-3,
    Y1 == 1 ~ Y1-1e-3,
    TRUE ~ Y1
  )) %>%
  mutate(Y2 = case_when(
    Y2 == 0 ~ Y2+1e-3,
    Y2 == 1 ~ Y2-1e-3,
    TRUE ~ Y2
  )) %>%
  mutate(Y3 = case_when(
    Y3 == 0 ~ Y3+1e-3,
    Y3 == 1 ~ Y3-1e-3,
    TRUE ~ Y3
  )) %>%
  as.matrix()

coefs <- sp %>% ungroup %>% select(p,s12,s13,s23,a12,a23) %>% mutate(d = p /(a12*a23)) %>% as.matrix() %>% {.[1,]} 

save_path <- paste(bdir,'varying_flow.RData',sep='/')
if(file.exists(save_path)){
  load(save_path)
}else{
  coefs_flow <- lapply(seq(.5,10,length.out = Ns2),function(frac){
    coefs[c('s12','s23','s13')] <- frac*coefs[c('s12','s23','s13')]
    return(coefs)
  }) %>% 
    do.call(rbind,.) %>%
    rbind(coefs)
  print(paste('Saving @',save_path))
  ## Parallel setup
  row_u0 <- 1
  row_c <- 1
  U <- lapply(1:nrow(coefs_flow),function(row_c){
    cl <- makeCluster(n_co,type = 'FORK')
    print(paste('Completed @',round(row_c/nrow(coefs_flow)*100,2)))
    
    U <- parLapply(cl,1:nrow(U0s),function(row_u0){
      # print(row_u0)
      U <- evol_to_stable_state(row_c = row_c,coefsi = coefs_flow,U0 = U0s[row_u0,]) %>%
        mutate(Y10 = U0s[row_u0,'Y1'],Y20 = U0s[row_u0,'Y2'], Y30 = U0s[row_u0,'Y3'])
      return(U)
    }) %>%
      bind_rows
    stopCluster(cl)
    return(U)
  }) %>% bind_rows
  
  save(U,file=save_path)
}

library(ggplot2)
U %>% 
  mutate(U0 = Y10 + a12*(Y20 + a23*Y30)) %>%
  tidyr::pivot_longer(c('Y1','Y2','Y3')) %>%
  ggplot(aes(x=U0,y=value,color=s12,group=s12)) +
  geom_line() +
  facet_wrap(~name) +
  scale_x_continuous(name = latex2exp::TeX('Initial total resource $U(t=0)$'),breaks=seq(0,10,1)) +
  scale_y_continuous(name = 'Fill level',breaks=seq(0,1,.2)) +
  scale_color_gradient(name = latex2exp::TeX('$s_{12}$')) +
  theme_minimal() + 
  theme(legend.position = c(.45,.8)) ->
  gg
print(gg)

ggsave('fills12_for_demo.svg',
       device = svglite::svglite,
       width= 120/25.4,
       height = 120/25.4,
       pointsize = 4)
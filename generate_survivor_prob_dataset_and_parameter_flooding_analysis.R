rm(list=ls())
gc()
library(dplyr)
library(stringr)
library(ggplot2)

# This section generates the dataset indicating the values of the coefficients at different flooding scenarios.
if(!file.exists('survivor_prob_f2.5resol.RData')){

  load('node_amounts_15cm.RData')
  node_amounts %>% glimpse
  df <- node_amounts
  # This fuction takes the dataset, and calculates the number of nodes of each type at each scenario.
  # Then, it calculates the number of links between nodes of each type in the multipartite network.
  # 
  survivorness <- function(sc,df){
    N1 <- df %>% filter(scenario == sc,Var1 == 'refinery') %>% select(Freq) %>% unlist %>% as.numeric
    N2 <- df %>% filter(scenario == sc,Var1 == 'terminal') %>% select(Freq) %>% unlist %>% as.numeric
    N3 <- df %>% filter(scenario == sc,Var1 == 'gas_station') %>% select(Freq) %>% unlist %>% as.numeric
    fgs_ref <- (df %>% filter(scenario == sc,Var1 == 'gas_station-refinery-road') %>% select(Freq) %>% unlist %>% as.numeric)/(df %>% filter(scenario == '0',Var1 == 'gas_station-refinery-road') %>% select(Freq) %>% unlist %>% as.numeric)
    fgs_term <- (df %>% filter(scenario == sc,Var1 == 'terminal-gas_station-road') %>% select(Freq) %>% unlist %>% as.numeric)/(df %>% filter(scenario == '0',Var1 == 'terminal-gas_station-road') %>% select(Freq) %>% unlist %>% as.numeric)
    fref_term <- (df %>% filter(scenario == sc,Var1 == 'terminal-refinery-road') %>% select(Freq) %>% unlist %>% as.numeric)/(df %>% filter(scenario == '0',Var1 == 'terminal-refinery-road') %>% select(Freq) %>% unlist %>% as.numeric)
    Mpipe <- df %>% filter(scenario == sc) %>% select(min_cut) %>% unlist %>% as.numeric %>% unique
    
    C1 <- seq(1,1.5,2.5e-2) * 38.2 # MillonGallons
    C2 <- seq(1,2,2.5e-2) * 31 # MillonGallons
    C3 <- 35 * 10**-3 # MillonGallons
    P <- 18.9 * 5 / N1 # MillionGallons Per Week
    W12M12pipe <- seq(2,4,length.out = 10) * 7 *  Mpipe # MillonGallons Per Week
    W12M12truck <- seq(3,7,length.out = 10) * 10**-3 * 7 *  (N1*N2-Mpipe) *  fref_term# MillonGallons Per Week
    W13M13 <- seq(3,11.6,length.out = 10) * 10**-3 * 7 * N1*200 * fgs_ref # MillonGallons Per Week
    W23M23 <- seq(3,7,length.out = 10) * 10**-3 * 7 * 5000 * fgs_term # MillonGallons Per Week
    W12M12 <- seq(min(W12M12pipe+W12M12truck),max(W12M12pipe+W12M12truck),length.out = 10)
    
    coef1 <- expand.grid(
      N1 = N1,
      N2 = N2,
      N3 = N3,
      C1 = C1,
      C2 = C2,
      C3 = C3,
      P = P)
    coef2 <- expand.grid(
      W12M12 = W12M12,
      W13M13 = W13M13,
      W23M23 = W23M23,
      P = P) 
    
    a_function <- function(x){
      sp <- x %>% summarise( 
        p = mean(p),
        s12 = mean(s12),
        s13=mean(s13),
        s23=mean(s23),
        a12 = mean(a12),
        a23 = mean(a23),
        N1C1 = mean(N1*C1),
        N2C2 = mean(N2*C2),
        N3C3 = mean(N3*C3),
        P = mean(N1*C1*p),
        D = mean(N3*C3*p),
        W12 = mean(s12*N1*C1),
        W13 = mean(s13*N1*C1),
        W23 = mean(s23*N2*C2)
      )
      lapply( 
        seq(0,100,2.5), 
        function(f){
          x %>% 
            summarise(
              u = mean(p*f/100 <= s13 + ifelse(s12 < a12*s23,s12,a12*s23)
              )) %>% 
            setNames(paste('survivor',f,sep=''))
        }) %>% 
        bind_cols(sp,.) %>%
        return()
    }
    
    left_join(coef1,coef2,by='P') %>%
      mutate(
        p = P/C1,
        a12 = N2*C2 / (N1 * C1),
        a23 = N3*C3 / (N2 * C2),
        s12 = W12M12/(N1*C1),
        s13 = W13M13/(N1*C1),
        s23 = W23M23/(N2*C2),
        d = p / (a12 * a23),
        D = C3*d
      ) %>%
      a_function() %>%
      mutate(scenario = sc) %>%
      return()
  }
  
  scenario_names <- unique(node_amounts$scenario)
  survivor_prob <- lapply(scenario_names,function(sc){
    print(sc)
    survivorness(sc,node_amounts)
  }) %>% bind_rows
  
  survivor_prob <- survivor_prob %>%
    tidyr::pivot_longer(
      cols = starts_with('survivor'),
      names_to = 'demand_level',
      values_to = 'survivor_prob') %>% 
    mutate(
      demand_level = as.numeric(str_extract(demand_level,'[0-9]+'))
    ) %>%
    mutate(year_range = str_extract(scenario,'[0-9]+_[0-9]+')) %>% 
    mutate(RCP = str_extract(str_remove(scenario,year_range),"[A-Z]+[0-9]+")) %>% 
    mutate(prec = str_extract(scenario,'prec[0-9]+')) %>% 
    mutate(model = str_extract(scenario,'[A-z]+[0-9]+_prec')) %>%
    mutate(year_range = str_replace(year_range,'_','-')) %>%
    mutate(prec = case_when(
      str_detect(prec,'prec500') ~ '50%',
      str_detect(prec,'prec950') ~ '95%',
      str_detect(prec,'prec999') ~ '99.9%'
    )) %>% 
    group_by(demand_level) %>%
    mutate(survivor_prob_corr = survivor_prob/survivor_prob[scenario == '0']) %>% 
    rename(survivor_prob0 = survivor_prob,
           survivor_prob = survivor_prob_corr)
  save(survivor_prob,file='survivor_prob_f2.5resol.RData')
}else{
  load('survivor_prob_f2.5resol.RData')
}

survivor_prob %>% 
  filter(scenario != '0') %>%
  group_by(year_range,prec,demand_level,RCP) %>%
  summarise( survivor_prob = mean(survivor_prob) ) %>%
  filter(year_range %in% c('2060-2080','2080-2100')) %>%
  ggplot(aes(x=demand_level,y=(1-survivor_prob)*100,color=prec)) +
  geom_line() +
  scale_x_continuous(name = 'Demand level [%]',breaks=seq(0,100,20)) +
  scale_color_discrete(
    name = 'Percentile',
    type = c('50%' = '#ffe949','95%'='orange','99.9%'='#c71f2d')
  ) +
  scale_y_continuous(name = 'Scenarios with failed demand [%]',limits=c(0,100),breaks = seq(0,100,20)) +
  facet_wrap( ~ RCP + year_range ) +
  theme_minimal() +
  theme(legend.position = c(.25,.8)) ->
  gg
print(gg)

ggsave('RCP_and_last_two_horizons.svg',
       device = svglite::svglite,
       width= 120/25.4,
       height = 90/25.4,
       pointsize = 4)


survivor_prob %>% 
  filter(scenario != '0') %>%
  filter(demand_level == 100) %>%
  group_by(year_range,prec) %>%
  summarise( survivor_mean = mean(survivor_prob),
             survivor_min = min(survivor_prob),
             survivor_max = max(survivor_prob)) %>% 
  ggplot(aes(x=year_range,y=(1-survivor_mean)*100,ymin = (1-survivor_max)*100,ymax=(1-survivor_min)*100,color=prec)) +
  geom_point(position = position_dodge(width=.5)) +
  geom_errorbar(width=0.1, position = position_dodge(width=.5)) +
  scale_x_discrete(name = 'Year range') +
  scale_color_discrete(
    name = 'Percentile',
    type = c('50%' = '#ffe949','95%'='orange','99.9%'='#c71f2d')
  ) +
  scale_y_continuous(
    name = 'Scenarios with failed demand [%]',
    limits=c(0,100)) +
  coord_flip() +
  theme_minimal() +
  theme(legend.position = c(.6,.4)) ->
  gg
print(gg)

ggsave('overview_floodfailure.svg',
       device = svglite::svglite,
       width= 90/25.4,
       height = 90/25.4,
       pointsize = 4)


survivor_prob %>% 
  filter(scenario != '0') %>%
  filter(demand_level == 100) %>%
  filter(prec == '99.9%') %>%
  group_by(year_range,RCP) %>%
  summarise( p = mean(p),
             s12 = mean(s12),
             s13 = mean(s13),
             s23 = mean(s23),
             a12 = mean(a12),
             a23 = mean(a23)) %>% 
  ungroup() %>%
  mutate(
    p = p/p[year_range == '2000-2020']*100-100,
    s12 = s12/s12[year_range == '2000-2020']*100-100,
    s13 = s13/s13[year_range == '2000-2020']*100-100,
    s23 = s23/s23[year_range == '2000-2020']*100-100,
    a12 = a12/a12[year_range == '2000-2020']*100-100,
    a23 = a23/a23[year_range == '2000-2020']*100-100
  ) %>% 
  mutate(
    RCP = RCP %>% str_replace('5','.5') %>% str_replace('RCP','RCP ')
  ) %>%
  tidyr::pivot_longer(c(p,s12,s13,s23,a12,a23),names_to = 'parameter_name',values_to = 'parameter_value') %>% 
  ggplot(aes(x = year_range, y = parameter_value, color = parameter_name,group=parameter_name)) +
  # geom_point(position = position_dodge(width=.5)) +
  geom_hline(yintercept = 100,lty='dashed') +
  geom_point() + geom_line() +
  scale_x_discrete(name = 'Year range') +
  scale_color_discrete(
    name = 'Parameter'
  ) +
  scale_y_continuous(
    name = 'Parameter change [%]',
    limits=c(-80,50),
    breaks = seq(-80,80,20)) +
  facet_wrap(~RCP,nrow=2) +
  theme_minimal() +
  theme(legend.position = c(.3,.8)) ->
  gg
print(gg)

ggsave("parameter_change.svg",
       device = svglite::svglite,
       width= 120/25.4,
       height = 90/25.4,
       pointsize = 4)



library(ggplot2)
library(dplyr)
library(stringr)
library(forcats)

##################### Now with all the scenarios

file_paths <- dir('failure_at_scenarios2023-08-04',full.names = TRUE)
fp <- file_paths[1]
load('survivor_prob_f2.5resol.RData')
df <- lapply(file_paths,function(fp){
  sc <- str_remove_all(fp,'failure_at_scenarios2023-08-04/failure_at_|.RData')
  sp <- survivor_prob %>%
    ungroup() %>%
    filter(scenario == sc,
           demand_level == 100) %>%
    select(scenario,N1C1,N2C2,N3C3,P,D,W12,W13,W23)
  envi <- new.env()
  load(fp,envi)
  envi$U %>%
    mutate(scenario = str_remove_all(fp,'failure_at_scenarios2023-08-04/failure_at_|.RData')) %>%
    mutate(year_range = str_extract(scenario,'[0-9]+_[0-9]+')) %>%
    mutate(RCP = str_extract(scenario,'RCP[0-9]+')) %>%
    mutate(prec = str_extract(scenario,'prec[0-9]+') %>% str_remove('prec')) %>%
    left_join(sp,by='scenario') %>%
    return()
}) %>%
  bind_rows

bind_rows(
  df %>% 
    filter(failure_duration == max(failure_duration)) %>%
    filter(year_range != '2000_2020'),
  df %>% 
    filter(failure_duration == max(failure_duration)) %>%
    filter(scenario == '0') %>%
    mutate(year_range = 'Original Network') %>% 
    {
      bind_rows(
        (.) %>% mutate(RCP = 'RCP45'),
        (.) %>% mutate(RCP = 'RCP85')
      )
    } %>%
    {
      bind_rows(
        (.) %>% mutate(prec = '500'),
        (.) %>% mutate(prec = '950'),
        (.) %>% mutate(prec = '999')
      )
    } 
) %>%
  filter(prec != '999') %>%
  mutate(year_range = str_replace(year_range,'_','-')) %>%
  mutate(prec = case_when(
    prec == '500' ~ 'Percentile 50%',
    prec == '950' ~ 'Percentile 95%'
  )) %>%
  mutate(RCP = str_replace(str_replace(str_replace(RCP,'45','4.5'),'85','8.5'),'RCP','RCP ')) %>%
  group_by(prec,RCP,year_range,row_u0) %>%
  summarise(U0 = mean(U0),tau = mean(tfail1e5)) %>%
  ggplot(aes(x= U0, y = tau, color=year_range,group=year_range)) +
  geom_point(size = 0.1) +
  geom_line() +
  facet_wrap(~prec + RCP, ncol = 2  ) +
  scale_color_discrete(name = NULL, breaks = c('Original Network','2020-2040','2040-2060','2060-2080','2080-2100')) +
  scale_x_continuous(name = latex2exp::TeX('Total initial resource $U(t=0)$'),breaks = seq(0.4,2,.1)) +
  scale_y_continuous(name = latex2exp::TeX('Time to demand failure $\\tau$ \\[weeks\\]'),breaks = seq(0,3,.5)) +
  theme_minimal() +
  theme(legend.position =  c(.1,.5),
        legend.background = element_rect(color='white')) ->
  gg
print(gg)
ggsave("time_to_failure_v_horizon.svg",
       device = svglite::svglite,
       width= 180/25.4,
       height = 180/25.4,
       pointsize = 6)


bind_rows(
  df %>% 
    filter(failure_duration == max(failure_duration)) %>%
    filter(year_range != '2000_2020'),
  df %>% 
    filter(failure_duration == max(failure_duration)) %>%
    filter(scenario == '0') %>%
    mutate(year_range = 'Original Network') %>% 
    {
      bind_rows(
        (.) %>% mutate(RCP = 'RCP45'),
        (.) %>% mutate(RCP = 'RCP85')
      )
    } %>%
    {
      bind_rows(
        (.) %>% mutate(prec = '500'),
        (.) %>% mutate(prec = '950'),
        (.) %>% mutate(prec = '999')
      )
    } 
) %>%
  filter(prec != '999') %>%
  mutate(year_range = str_replace(year_range,'_','-')) %>%
  mutate(prec = case_when(
    prec == '500' ~ 'Percentile 50%',
    prec == '950' ~ 'Percentile 95%'
  )) %>%
  mutate(RCP = str_replace(str_replace(str_replace(RCP,'45','4.5'),'85','8.5'),'RCP','RCP ')) %>%
  group_by(prec,RCP,year_range,row_u0) %>%
  summarise(U0 = mean(N1C1*U0),tau = mean(tfail1e5)) %>%
  ggplot(aes(x= U0, y = tau, color=year_range,group=year_range)) +
  geom_point(size = 0.5) +
  geom_line() +
  # geom_smooth(method = 'lm') +
  facet_wrap(~prec + RCP, ncol = 2  ) +
  scale_color_discrete(name = NULL, breaks = c('Original Network','2020-2040','2040-2060','2060-2080','2080-2100')) +
  scale_x_continuous(name = latex2exp::TeX('Total initial resource $U(t=0)$ \\[MG\\]'),breaks=seq(200,500,25)) +
  scale_y_continuous(name = latex2exp::TeX('Time to demand failure $\\tau$ \\[weeks\\]'),breaks = seq(0,3,.5)) +
  theme_minimal() +
  theme(legend.position =  c(.1,.5),
        legend.background = element_rect(color='white')) ->
  gg
print(gg)
ggsave("time_to_failure_v_horizon_wU.svg",
       device = svglite::svglite,
       width= 180/25.4,
       height = 180/25.4,
       pointsize = 6)


bind_rows(
  df %>% 
    # filter(failure_duration == max(failure_duration)) %>%
    filter(year_range != '2000_2020'),
  df %>% 
    # filter(failure_duration == max(failure_duration)) %>%
    filter(scenario == '0') %>%
    mutate(year_range = 'Original Network') %>% 
    {
      bind_rows(
        (.) %>% mutate(RCP = 'RCP45'),
        (.) %>% mutate(RCP = 'RCP85')
      )
    } %>%
    {
      bind_rows(
        (.) %>% mutate(prec = '500'),
        (.) %>% mutate(prec = '950'),
        (.) %>% mutate(prec = '999')
      )
    } 
) %>%
  filter(prec != '999') %>%
  mutate(year_range = str_replace(year_range,'_','-')) %>%
  mutate(prec = case_when(
    prec == '500' ~ 'Percentile 50%',
    prec == '950' ~ 'Percentile 95%'
  )) %>%
  mutate(RCP = str_replace(str_replace(str_replace(RCP,'45','4.5'),'85','8.5'),'RCP','RCP ')) %>%
  group_by(prec,RCP,year_range,failure_duration) %>%
  summarise(qD = mean(qD)) %>%
  ggplot(aes(x= failure_duration, y = qD, color=year_range,group=year_range)) +
  geom_point(size = 0.5) +
  geom_line() +
  # geom_smooth(method = 'lm') +
  facet_wrap(~prec + RCP, ncol = 2  ) +
  scale_color_discrete(name = NULL, breaks = c('Original Network','2020-2040','2040-2060','2060-2080','2080-2100')) +
  scale_x_continuous(name = latex2exp::TeX('Failure duration $\\Delta T$ \\[weeks\\]'),breaks = seq(0,3,.5)) +
  scale_y_continuous(name = latex2exp::TeX('Average demand level $Q_D$'),breaks = seq(0,1,.1)) +
  theme_minimal() +
  theme(legend.position =  c(.2,.75),
        legend.background = element_rect(color='white')) ->
  gg
print(gg)
ggsave("qD_v_horizon.svg",
       device = svglite::svglite,
       width= 180/25.4,
       height = 180/25.4,
       pointsize = 6)

bind_rows(
  df %>% 
    # filter(failure_duration == max(failure_duration)) %>%
    filter(year_range != '2000_2020'),
  df %>% 
    # filter(failure_duration == max(failure_duration)) %>%
    filter(scenario == '0') %>%
    mutate(year_range = 'Original Network') %>% 
    {
      bind_rows(
        (.) %>% mutate(RCP = 'RCP45'),
        (.) %>% mutate(RCP = 'RCP85')
      )
    } %>%
    {
      bind_rows(
        (.) %>% mutate(prec = '500'),
        (.) %>% mutate(prec = '950'),
        (.) %>% mutate(prec = '999')
      )
    } 
) %>%
  filter(prec != '999') %>%
  mutate(year_range = str_replace(year_range,'_','-')) %>%
  mutate(prec = case_when(
    prec == '500' ~ 'Percentile 50%',
    prec == '950' ~ 'Percentile 95%'
  )) %>%
  mutate(RCP = str_replace(str_replace(str_replace(RCP,'45','4.5'),'85','8.5'),'RCP','RCP ')) %>%
  group_by(prec,RCP,year_range,failure_duration) %>%
  summarise(qD = mean(N3C3*d*qD)) %>% 
  group_by(prec,RCP,failure_duration) %>%
  mutate(qD = qD - qD[year_range=='Original Network']) %>%
  filter(year_range != 'Original Network') %>%
  ggplot(aes(x= failure_duration, y = qD, color=year_range,group=year_range)) +
  geom_point(size = 0.5) +
  geom_line() +
  # geom_smooth(method = 'lm') +
  facet_wrap(~prec + RCP, ncol = 2  ) +
  scale_color_discrete(name = NULL, breaks = c('Original Network','2020-2040','2040-2060','2060-2080','2080-2100')) +
  scale_x_continuous(name = latex2exp::TeX('Failure duration $\\Delta T$ \\[weeks\\]'),breaks = seq(0,3,.5)) +
  scale_y_continuous(name = latex2exp::TeX('Average demand reduction $Q_D-Q_D^0$ \\[MG/week\\]'),limits=c(-16,0),breaks=seq(-16,0,4)) +
  theme_minimal() +
  theme(legend.position =  c(.2,.75),
        legend.background = element_rect(color='white')) ->
  gg
print(gg)
ggsave("qD_v_horizon_wU.svg",
       device = svglite::svglite,
       width= 180/25.4,
       height = 180/25.4,
       pointsize = 6)

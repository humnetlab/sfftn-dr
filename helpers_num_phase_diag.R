gen_coef_from_df <- function(row,sp,Ns){
  spr <- sp[row,]
  
  coefs_s12_s23 <- expand.grid(
    p = spr$p,
    a12 = spr$a12,
    a23 = spr$a23,
    s12 = seq(.8*min(sp$ms12),1.25*max(sp$Ms12),length.out = Ns),
    s13 = spr$s13,
    s23 = seq(.8*min(sp$ms23),1.25*max(sp$Ms23),length.out = Ns)) %>% 
    mutate(
      d = p/a12/a23
    ) %>%
    as.matrix
  
  coefs_s12_s13 <- expand.grid(
    p = spr$p,
    a12 = spr$a12,
    a23 = spr$a23,
    s12 = seq(.8*min(sp$ms12),1.25*max(sp$Ms12),length.out = Ns),
    s13 = seq(.8*min(sp$ms13),1.25*max(sp$Ms13),length.out = Ns),
    s23 = spr$s23) %>% 
    mutate(
      d = p/a12/a23
    ) %>%
    as.matrix
  
  coefs_s23_s13 <- expand.grid(
    p = spr$p,
    a12 = spr$a12,
    a23 = spr$a23,
    s12 = spr$s12,
    s13 = seq(.8*min(sp$ms13),1.25*max(sp$Ms13),length.out = Ns),
    s23 = seq(.8*min(sp$ms23),1.25*max(sp$Ms23),length.out = Ns)) %>% 
    mutate(
      d = p/a12/a23
    ) %>%
    as.matrix
  
  coefs_a12_a23 <- expand.grid(
    p = spr$p,
    a12 = seq(.8*min(sp$ma12),1.25*max(sp$Ma12),length.out = Ns),
    a23 = seq(.8*min(sp$ma23),1.25*max(sp$Ma23),length.out = Ns),
    s12 = spr$s12,
    s13 = spr$s13,
    s23 = spr$s23) %>% 
    mutate(
      d = p/a12/a23
    ) %>%
    as.matrix
  
  r <- list('s12_s23' = coefs_s12_s23,
            's12_s13' = coefs_s12_s13,
            's23_s13' = coefs_s23_s13,
            'a12_a23' = coefs_a12_a23)
  return(r)
}


evol_to_stable_state <- function(row_c,coefsi,U0 = NULL){
  coef_vec <- coefsi[row_c,]
  weight_matrix <- construct_weight_matrix_from_coefs_3x3(coef_vec)
  if(is.null(U0)){
    U0 <- rep(0.5,3)
    names(U0) <- paste('Y',1:3,sep='')
  }
  parms_ap <- prepare_parms_approx(
    weight_matrix = weight_matrix,
    dem_coef = c(0,0,coef_vec['d']), 
    prod_coef = c(coef_vec['p'],0,0),
    epsilon = 1e-6
  )
  times_parms <- list(length = 6e3,tMax=5,tMin=0) 
  
  U <- simulate_exact_system(
    U0 = U0, dU = dU_exact, parms = parms_ap,
    times_parms = times_parms
  ) %>%
    as.data.frame() 
  
  test_stability <- U %>%
    tail(n = 10) %>%
    mutate(dY1 = Y1 - lag(Y1,2),
           dY2 = Y2 - lag(Y2,2),
           dY3 = Y3 - lag(Y3,2),
           dt = time - lag(time,2)) %>%
    mutate(Y1p = dY1/dt,Y2p = dY2/dt, Y3p = dY2/dt) %>%
    tail() %>%
    summarise(modul = mean(sqrt(Y1p**2+Y2p**2+Y3p**2))) %>% 
    unlist() %>%
    {. < 1e-4}
  while(!test_stability){
    U0 <- tail(U,n=1)[c('Y1','Y2','Y3')] %>% as.matrix %>% c
    names(U0) <- c('Y1','Y2','Y3')
    U <- simulate_exact_system(
      U0 = U0, dU = dU_exact, parms = parms_ap,
      times_parms = times_parms
    ) %>%
      as.data.frame
    test_stability <- U %>%
      tail(n = 10) %>%
      mutate(dY1 = Y1 - lag(Y1,2),
             dY2 = Y2 - lag(Y2,2),
             dY3 = Y3 - lag(Y3,2),
             dt = time - lag(time,2)) %>%
      mutate(Y1p = dY1/dt,Y2p = dY2/dt, Y3p = dY2/dt) %>%
      tail() %>%
      summarise(modul = mean(sqrt(Y1p**2+Y2p**2+Y3p**2))) %>% 
      unlist() %>%
      {. < 1e-4}
  }
  
  U <- U %>%
    tail(n=1) %>%
    mutate(row_c = row_c) %>%
    bind_cols(as.data.frame(t(as.data.frame(coef_vec))))
  rownames(U) <- NULL
  return(U)
}

gen_coef_for_failures <- function(sc,df){
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
      a23 = mean(a23)
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
    as.matrix() %>%
    return()
}
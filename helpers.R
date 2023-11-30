library(deSolve)
library(rootSolve)
library(igraph)
library(Matrix)

sample_mfn <- function(layer.sizes,pref.matrix,weights_range = c(1,1)){
  n <- sum(layer.sizes)
  g <- sample_sbm(
    n = n,
    pref.matrix = pref.matrix,
    block.sizes = layer.sizes,
    directed = TRUE,
    loops = FALSE
  )
  V(g)$type <- unlist(lapply(1:length(layer.sizes),function(i) rep(i,layer.sizes[i])))
  E(g)$weight <- runif(ecount(g),weights_range[1],weights_range[2])
  layout <- lapply(1:length(layer.sizes),function(i) expand.grid(1:layer.sizes[i],i)) %>% bind_rows %>% as.matrix
  V(g)$x <- layout[,1]
  V(g)$y <- layout[,2]
  return(g)
}

correct_mfn <- function(g){
  M <- max(V(g)$type)
  m <- min(V(g)$type)
  if(any(E(g)$weight==0))
    g <- delete.edges(g,E(g)$weight==0)
  
  if(any( (degree(g,mode='out') == 0 |
           degree(g,mode='in') == 0)  &
          !is.element(
            V(g)$type, c(m,M)
          )))
    g <- delete.vertices(g,
                         (degree(g,mode='out') == 0 |
                            degree(g,mode='in') == 0)  &
                           !is.element(
                             V(g)$type, c(m,M)
                           ))
  
  if(any( degree(g,mode='in') == 0  & V(g)$type == M))
    g <- delete.vertices(g,
                         degree(g,mode='in') == 0  &
                           V(g)$type == M
    )
  
  if(any( degree(g,mode='out') == 0  & V(g)$type == m))
    g <- delete.vertices(g,
                         degree(g,mode='out') == 0  &
                           V(g)$type == m
    )

  return(g)
}

sample_mfn_corrected <- function(layer.sizes,pref.matrix,weights_range = c(1,1)){
  g <- sample_mfn(layer.sizes,pref.matrix,weights_range)
  g <- correct_mfn(g)
  while(length(unique(V(g)$type)) < length(layer_sizes)){
    g <- sample_mfn(layer.sizes,pref.matrix,weights_range)
    g <- correct_mfn(g)
  }
  return(g)
}



beta_matrix_fun <- function(parms){
  parms$M / parms$N *parms$AVW
}

prepare_parms <- function(parms){
  c(parms,list(
    prod_fun = function(u,parms){
      ifelse(
        parms$prod_coef>0,
        parms$prod_coef * (1-u) / (1 + parms$epsilon - u),
        0)
    },
    dem_fun = function(u,parms){
      ifelse(
        parms$dem_coef>0,
        parms$dem_coef * u / ( parms$epsilon + u),
        0)
    },
    flow_fun = function(u,parms){
      u <- u / parms$u_max
      M <- matrix(nrow = length(u),ncol=length(u))
      for(i in 1:(length(u)-1)){
        for(j in (i+1):length(u)){
          M[i,j] <- u[i] * (1-u[j])
          M[j,i] <- - M[i,j]
        }
      }
      diag(M) <- 0
      return(M)
    }
    
    
  ))
}

dU <- function(t,u,parms){
  #'@description Function for deSolve solver. It uses the actual state of the system to calculate the time derivative.
  #'@param t the time
  #'@param u The state of the system, as X1,...,XN (vector)
  #'@param parms parms is the parameter list. Each element has to be named. 
  prod_term <- parms$prod_fun(u,parms)
  dem_term <- parms$dem_fun(u,parms)
  flow_matrix <- parms$flow_fun(u,parms)
  result <- prod_term - dem_term - rowSums(flow_matrix * parms$beta_matrix)
  return(list(c(result)))
}

simulate_autonomous <- function(U0,dU,parms,
                                times = NULL,times_parms = list(length = 1e4,tMax=10,tMin=0)){
  if(is.null(times)) times <- seq(times_parms$tMin,times_parms$tMax,length.out = times_parms$length)
  if(is.null(names(U0))) names(U0) <- paste('Y',1:length(U0)-1,sep='')
  return(ode(y = U0,times = times,func = dU,parms = parms))
}

dU_exact <- function(t,u,parms){
  prod_term <- parms$prod_fun(t,u,parms)
  dem_term <- parms$dem_fun(t,u,parms)
  flow_term <- parms$flow_fun(t,u,parms)
  result <- prod_term - dem_term + flow_term
  for(i in 1:length(u)){
    if(u[i] > 1){
      result[i] <- -abs(result[i])
    }else if(u[i] < 0){
      result[i] <- abs(result[i])
    }
  }
  return(list(result))
}

prepare_parms_exact <- function(parms){
      parms$prod_fun <- function(t,u,parms){
        ifelse(
          parms$prod_coef > 0,
          parms$prod_coef * (1-u) / (1 + parms$epsilon - u),
          0)
      }
      parms$dem_fun <- function(t,u,parms){
        ifelse(
          parms$dem_coef > 0,
          parms$dem_coef * u / ( parms$epsilon + u),
          0)
      }
      
      parms$nei_in <- lapply(
        1:parms$N,
        function(i) neighbors(parms$g,i,'in'))
      parms$nei_out <- lapply(
        1:parms$N,
        function(i) neighbors(parms$g,i,'out'))
      parms$nei_in <- lapply(
        1:parms$N,
        function(i) neighbors(parms$g,i,'in'))
      
      parms$nei_out_weight <- lapply(
        1:parms$N,
        function(i){
          ni <- parms$nei_out[[i]]
          if(length(ni) > 0){
            r <- E(parms$g)[get.edge.ids(parms$g,c(rbind(i,ni)))]$weight
          }else{
            r <- 0
          }
          return(r)
        }
      )
      
      parms$nei_in_weight <- lapply(
        1:parms$N,
        function(i){
          ni <- parms$nei_in[[i]]
          if(length(ni) > 0){
            r <- E(parms$g)[get.edge.ids(parms$g,c(rbind(ni,i)))]$weight
          }else{
            r <- 0
          }
          return(r)
        }
      )
      
      parms$flow_fun <- function(t,u,parms){
        flow_in <- sapply(1:length(u),function(i){
          if(length(unlist(parms$nei_in[[i]]))  > 0){
            r <- sum(
              parms$flow_fun0(
                t,
                u[parms$nei_in[[i]]],
                u[i]
              ) * parms$nei_in_wei[[i]]
            )
          }else{
            r <- 0
          }
          return(r)
        })
        flow_out <- sapply(1:length(u),function(i){
          if(length(parms$nei_out[[i]])  > 0){
            r <- sum(
              parms$flow_fun0(
                t,
                u[i],
                u[parms$nei_out[[i]]]) * 
                parms$nei_out_wei[[i]]
            )
          }else{
            r <- 0
          }
          return(r)
        })
        return(flow_in - flow_out)
      }
      return(parms)
}


prepare_parms_exact_failure <- function(parms){
  parms$prod_fun <- function(t,u,parms){
    ifelse(
      parms$prod_coef > 0 & 
        ifelse(parms$failing_nodes,
               t < parms$failure_time_range[1] | t > parms$failure_time_range[2],
               TRUE),
      parms$prod_coef * (1-u) / (1 + parms$epsilon - u),
      0)
  }
  parms$dem_fun <- function(t,u,parms){
    ifelse(
      parms$dem_coef > 0,
      parms$dem_coef * u / ( parms$epsilon + u),
      0)
  }
  
  parms$nei_in <- lapply(
    1:parms$N,
    function(i) neighbors(parms$g,i,'in'))
  parms$nei_out <- lapply(
    1:parms$N,
    function(i) neighbors(parms$g,i,'out'))
  
  
  parms$nei_out_weight <- lapply(
    1:parms$N,
    function(i){
      ni <- parms$nei_out[[i]]
      if(length(ni) > 0){
        r <- E(parms$g)[get.edge.ids(parms$g,c(rbind(i,ni)))]$weight
      }else{
        r <- 0
      }
      return(r)
    }
  )
  
  parms$nei_in_weight <- lapply(
    1:parms$N,
    function(i){
      ni <- parms$nei_in[[i]]
      if(length(ni) > 0){
        r <- E(parms$g)[get.edge.ids(parms$g,c(rbind(ni,i)))]$weight
      }else{
        r <- 0
      }
      return(r)
    }
  )
  
  parms$flow_fun <- function(t,u,parms){
    flow_in <- sapply(1:length(u),function(i){
      if(length(unlist(parms$nei_in[[i]]))  > 0){
        r <- sum(
          parms$flow_fun0(
            t,
            u[parms$nei_in[[i]]],
            u[i]
          ) * parms$nei_in_wei[[i]]
        )
      }else{
        r <- 0
      }
      return(r)
    })
    flow_out <- sapply(1:length(u),function(i){
      if(length(parms$nei_out[[i]])  > 0){
        r <- sum(
          parms$flow_fun0(
            t,
            u[i],
            u[parms$nei_out[[i]]]) * 
            parms$nei_out_wei[[i]]
        )
      }else{
        r <- 0
      }
      return(r)
    })
    return(flow_in - flow_out)
  }
  return(parms)
}


prepare_parms_approx_failure <- function(parms){
  parms$prod_fun <- function(t,u,parms){
    ifelse(
      parms$prod_coef > 0, 
        ifelse(
          t < parms$failure_time_range[1] | t > parms$failure_time_range[2],
          1,
          1-parms$fraction_failed) * parms$prod_coef * (1-u) / (1 + parms$epsilon - u),
      0)
  }
  parms$dem_fun <- function(t,u,parms){
    ifelse(
      parms$dem_coef > 0,
      parms$dem_coef * u / ( parms$epsilon + u),
      0)
  }
  
  parms$nei_in <- lapply(
    1:parms$N,
    function(i) neighbors(parms$g,i,'in'))
  parms$nei_out <- lapply(
    1:parms$N,
    function(i) neighbors(parms$g,i,'out'))
  parms$nei_in <- lapply(
    1:parms$N,
    function(i) neighbors(parms$g,i,'in'))
  
  parms$nei_out_weight <- lapply(
    1:parms$N,
    function(i){
      ni <- parms$nei_out[[i]]
      if(length(ni) > 0){
        r <- E(parms$g)[get.edge.ids(parms$g,c(rbind(i,ni)))]$weight
      }else{
        r <- 0
      }
      return(r)
    }
  )
  
  parms$nei_in_weight <- lapply(
    1:parms$N,
    function(i){
      ni <- parms$nei_in[[i]]
      if(length(ni) > 0){
        r <- E(parms$g)[get.edge.ids(parms$g,c(rbind(ni,i)))]$weight
      }else{
        r <- 0
      }
      return(r)
    }
  )
  
  parms$flow_fun <- function(t,u,parms){
    flow_in <- sapply(1:length(u),function(i){
      if(length(unlist(parms$nei_in[[i]]))  > 0){
        r <- sum(
          parms$flow_fun0(
            t,
            u[parms$nei_in[[i]]],
            u[i]
          ) * parms$nei_in_wei[[i]]
        )
      }else{
        r <- 0
      }
      return(r)
    })
    flow_out <- sapply(1:length(u),function(i){
      if(length(parms$nei_out[[i]])  > 0){
        r <- sum(
          parms$flow_fun0(
            t,
            u[i],
            u[parms$nei_out[[i]]]) * 
            parms$nei_out_wei[[i]]
        )
      }else{
        r <- 0
      }
      return(r)
    })
    return(flow_in - flow_out)
  }
  return(parms)
}


simulate_autonomous <- function(U0,dU,parms,
                                times = NULL,times_parms = list(length = 1e4,tMax=10,tMin=0)){
  if(is.null(times)) times <- seq(times_parms$tMin,times_parms$tMax,length.out = times_parms$length)
  if(is.null(names(U0))) names(U0) <- paste('Y',1:length(U0)-1,sep='')
  return(ode(y = U0,times = times,func = dU,parms = parms))
}

simulate_exact_system <- function(
    U0,dU,parms,times = NULL,
    times_parms = list(length = 1e4,tMax=10,tMin=0)
){
  if(is.null(times)) times <- seq(times_parms$tMin,times_parms$tMax,length.out = times_parms$length)
  if(is.null(names(U0)) & !is.null(parms$g)) names(U0) <- paste('X',1:length(U0),'_type_',V(parms$g)$type,sep='')
  return(ode(y = U0,times = times,func = dU,parms = parms))  
}

stable_state_exact_system <- function(
    U0,dU,parms,
    times_parms = list(length = 1e4,tMax=10,tMin=0)
){
  if(is.null(names(U0)) & !is.null(parms$g)) names(U0) <- paste('X',1:length(U0),'_type_',V(parms$g)$type,sep='')
  return(stode(y = U0,func = dU,parms = parms,positive=TRUE))  
}

generate_parm_set <- function(
    layer_sizes = c(10,10,10),
    pref_matrix = matrix(
      c(0,1,1,
        0,0,1,
        0,0,0),
      byrow=TRUE,ncol=3),
    epsilon = 1e-4,
    flow_fun0 = function(t,u1,u2) u1 * (1-u2),
    weights_range = c(1,1),
    prod_level = 1, dem_level = 1,
    corrected_network = FALSE
){
  #-> INITIALIZATION OF PARAMETERS FOR EXACT SYSTEM
  parms_ex <- list(
    layer.sizes = layer_sizes,
    pref.matrix = pref_matrix,
    epsilon = epsilon,
    flow_fun0 = flow_fun0
  )
  parms_ex$N <- sum(parms_ex$layer.sizes)
  if(!corrected_network){
    parms_ex$g <- sample_mfn(parms_ex$layer.sizes,parms_ex$pref.matrix,weights_range = weights_range)  
  }else{
    parms_ex$g <- sample_mfn_corrected(parms_ex$layer.sizes,parms_ex$pref.matrix,weights_range = weights_range)  
    parms_ex$layer.sizes <- sapply( 1:length(parms_ex$layer.sizes),function(lay) sum(V(parms_ex$g)$type==lay) )
    parms_ex$N <- sum(parms_ex$layer.sizes)
  }
  
  parms_ex$prod_coef <- as.numeric(V(parms_ex$g)$type == min(V(parms_ex$g)$type)) * prod_level
  parms_ex$dem_coef <- as.numeric(V(parms_ex$g)$type == max(V(parms_ex$g)$type)) * dem_level
  parms_ex <- prepare_parms_exact(parms_ex)
  #END-> INITIALIZATION OF PARAMETERS FOR EXACT SYSTEM
  #-> INITIALIZATION OF PARAMETERS FOR APPROX SYSTEM
  parms_ap <- list(
    layer.sizes = c(1,1,1),
    pref.matrix = matrix(
      c(0,1,1,
        0,0,1,
        0,0,0),
      byrow=TRUE,ncol=3),
    epsilon = 1e-4,
    flow_fun0 = parms_ex$flow_fun0
  )
  parms_ap$N <- sum(parms_ap$layer.sizes)
  parms_ap$g <- sample_mfn(parms_ap$layer.sizes,parms_ap$pref.matrix)
  parms_ap$prod_coef <- c(1,0,0) * prod_level
  parms_ap$dem_coef <- c(0,0,1) * dem_level
  parms_ap <- prepare_parms_exact(parms_ap)
  eds <- ends(parms_ex$g,E(parms_ex$g))
  eds <- cbind(V(parms_ex$g)[eds[,1]]$type,V(parms_ex$g)[eds[,2]]$type)
  # TODO: Generalize this for an arbitrary number of layers
  parms_ap$nei_in_weight <- list(
    0,
    sum(E(parms_ex$g)$weight[eds[,1]==1 & eds[,2]==2])/parms_ex$layer.sizes[2],
    c(sum(E(parms_ex$g)$weight[eds[,1]==1 & eds[,2]==3]),sum(E(parms_ex$g)$weight[eds[,1]==2 & eds[,2]==3]))/parms_ex$layer.sizes[3])
  
  parms_ap$nei_out_weight <- list(
    c(sum(E(parms_ex$g)$weight[eds[,1]==1 & eds[,2]==2]),
      sum(E(parms_ex$g)$weight[eds[,1]==1 & eds[,2]==3]))/
      parms_ex$layer.sizes[1],
    sum(E(parms_ex$g)$weight[eds[,1]==2 & eds[,2]==3])/parms_ex$layer.sizes[2],
    0)
  #END-> INITIALIZATION OF PARAMETERS FOR APPROX SYSTEM
  return(list('exact'=parms_ex,'approx'=parms_ap))
}

generate_parm_set_failure <- function(
    layer_sizes = c(10,10,10),
    pref_matrix = matrix(
      c(0,1,1,
        0,0,1,
        0,0,0),
      byrow=TRUE,ncol=3),
    epsilon = 1e-4,
    flow_fun0 = function(t,u1,u2) u1 * (1-u2),
    weights_range = c(1,1),
    prod_level = 1, dem_level = 1,
    failure_time_range = c(1,2),
    failing_nodes_fraction = 1,
    corrected_network = FALSE
){
  #-> INITIALIZATION OF PARAMETERS FOR EXACT SYSTEM
  parms_ex <- list(
    layer.sizes = layer_sizes,
    pref.matrix = pref_matrix,
    epsilon = epsilon,
    flow_fun0 = flow_fun0,
    failure_time_range = failure_time_range
  )
  parms_ex$N <- sum(parms_ex$layer.sizes)
  if(!corrected_network){
    parms_ex$g <- sample_mfn(parms_ex$layer.sizes,parms_ex$pref.matrix,weights_range = weights_range)  
  }else{
    parms_ex$g <- sample_mfn_corrected(parms_ex$layer.sizes,parms_ex$pref.matrix,weights_range = weights_range)  
    parms_ex$layer.sizes <- sapply( 1:length(parms_ex$layer.sizes),function(lay) sum(V(parms_ex$g)$type==lay) )
    parms_ex$N <- sum(parms_ex$layer.sizes)
  }
  
  parms_ex$prod_coef <- as.numeric(V(parms_ex$g)$type == min(V(parms_ex$g)$type)) * prod_level
  parms_ex$dem_coef <- as.numeric(V(parms_ex$g)$type == max(V(parms_ex$g)$type)) * dem_level
  parms_ex$failing_nodes <- as.numeric(V(parms_ex$g)$type == min(V(parms_ex$g)$type)) * ( runif(vcount(parms_ex$g)) < failing_nodes_fraction )
  
  parms_ex <- prepare_parms_exact_failure(parms_ex)
  #END-> INITIALIZATION OF PARAMETERS FOR EXACT SYSTEM
  #-> INITIALIZATION OF PARAMETERS FOR APPROX SYSTEM
  parms_ap <- list(
    layer.sizes = c(1,1,1),
    pref.matrix = matrix(
      c(0,1,1,
        0,0,1,
        0,0,0),
      byrow=TRUE,ncol=3),
    epsilon = 1e-4,
    flow_fun0 = parms_ex$flow_fun0,
    failure_time_range = failure_time_range
  )
  parms_ap$N <- sum(parms_ap$layer.sizes)
  parms_ap$g <- sample_mfn(parms_ap$layer.sizes,parms_ap$pref.matrix)
  parms_ap$fraction_failed <- sum(parms_ex$failing_nodes) / sum(parms_ex$prod_coef > 0)
  parms_ap$prod_coef <- c(1,0,0) * prod_level
  parms_ap$dem_coef <- c(0,0,1) * dem_level
  parms_ap <- prepare_parms_approx_failure(parms_ap)
  eds <- ends(parms_ex$g,E(parms_ex$g))
  eds <- cbind(V(parms_ex$g)[eds[,1]]$type,V(parms_ex$g)[eds[,2]]$type)
  # TODO: Generalize this for an arbitrary number of layers
  parms_ap$nei_in_weight <- list(
    0,
    sum(E(parms_ex$g)$weight[eds[,1]==1 & eds[,2]==2])/parms_ex$layer.sizes[2],
    c(sum(E(parms_ex$g)$weight[eds[,1]==1 & eds[,2]==3]),sum(E(parms_ex$g)$weight[eds[,1]==2 & eds[,2]==3]))/parms_ex$layer.sizes[3])
  
  parms_ap$nei_out_weight <- list(
    c(sum(E(parms_ex$g)$weight[eds[,1]==1 & eds[,2]==2]),sum(E(parms_ex$g)$weight[eds[,1]==1 & eds[,2]==3]))/parms_ex$layer.sizes[1],
    sum(E(parms_ex$g)$weight[eds[,1]==2 & eds[,2]==3])/parms_ex$layer.sizes[2],
    0)
  #END-> INITIALIZATION OF PARAMETERS FOR APPROX SYSTEM
  return(list('exact'=parms_ex,'approx'=parms_ap))
}




run_exact_and_approx <- function(
    U0, parms_full,
    times_parms = list(tMin = 0, tMax = 10,length=5e3)
){
  
  U_exact <- simulate_exact_system(
    U0 = U0, dU = dU_exact, parms = parms_full$exact,
    times_parms = times_parms
  ) %>%
    as.data.frame() %>%
    mutate(U0_t = mean(U0)) %>%
    mutate(Y1e = select(.,contains('type_1')) %>% rowMeans(),
           Y2e = select(.,contains('type_2')) %>% rowMeans(),
           Y3e = select(.,contains('type_3')) %>% rowMeans(),
           S1e = select(.,contains('type_1')) %>% {.**2} %>% rowMeans(),
           S2e = select(.,contains('type_2')) %>% {.**2} %>% rowMeans(),
           S3e = select(.,contains('type_3')) %>% {.**2} %>% rowMeans()) %>%
    mutate(
      S1e = sqrt(abs(S1e-Y1e**2)),
      S2e = sqrt(abs(S2e-Y2e**2)),
      S3e = sqrt(abs(S3e-Y3e**2)))
  
  U0_ap <- c(mean(U0[V(parms_full$exact$g)$type == 1]),
             mean(U0[V(parms_full$exact$g)$type == 2]),
             mean(U0[V(parms_full$exact$g)$type == 3]))
  U_ap <- simulate_exact_system(
    U0 = U0_ap, dU = dU_exact, parms = parms_full$approx,
    times_parms = times_parms
  ) %>%
    as.data.frame() %>%
    rename(Y1a = X1_type_1, Y2a=X2_type_2,Y3a=X3_type_3)
  
  U <- U_exact %>% 
    left_join(U_ap,by = 'time') %>%
    mutate(
      bias1 = (Y1a-Y1e),
      bias2 = (Y2a-Y2e),
      bias3 = (Y3a-Y3e),
      S1a = sqrt(S1e**2 + bias1**2),
      S2a = sqrt(S2e**2 + bias2**2),
      S3a = sqrt(S3e**2 + bias3**2))
  return(U)
}


calculate_resilience_system_resource <- function(
    U,time_range){
  U %>%
    filter(time >= time_range[1] & time <= time_range[2]) %>%
    mutate(
      Re = select(.,starts_with('Y') & ends_with('e')) %>% rowSums(),
      Ra = select(.,starts_with('Y') & ends_with('a')) %>% rowSums()
    ) %>% 
    summarise(
      dt = first(diff(time)),
      iRe = sum(Re*dt),
      iRa = sum(Ra*dt),
      R0a = first(Ra),
      R0e = first(Re),
      DT = diff(range(time))
    ) %>%
    summarise(RRe = R0e*DT-iRe,
              RRa = R0a*DT-iRa)
}

calculate_resilience_demand_failure <- function(U,demand_level){
  U %>%
    mutate(across(contains('type_1'),~ (1-.x)/(1-.x+1e-4),.names='P_{.col}')) %>%
    mutate(Pe = select(.,contains('P_')) %>% rowMeans()) %>% 
    mutate(across(contains('type_3'),~ .x/(.x+1e-4),.names='D_{.col}')) %>%
    mutate(De = select(.,contains('D_')) %>% rowMeans()) %>% 
    mutate('Da' = Y3a/(Y3a+1e-3),'Pa'=(1-Y1a)/(1-Y1a+1e-4)) %>% summarise(
      tDa = filter(.,Da < demand_level) %>% summarise(first(time)) %>% unlist,
      tDe = filter(.,De < demand_level) %>% summarise(first(time)) %>% unlist
    )
}


prepare_parms_approx <- function(
    epsilon = 1e-4,
    prod_coef = c(1,0,0), dem_coef = c(0,0,1),
    flow_fun0 = function(x,parms) {
      x * c((parms$weight_matrix  * (parms$weight_matrix < 0)) %*% (1-x)) +
              (1-x) * c((parms$weight_matrix  * (parms$weight_matrix > 0)) %*% x)
      },
    weight_matrix = matrix(c(0,-1,-1,1,0,-1,1,1,0),ncol=3,byrow=TRUE)
){
  parms <- list()
  parms$prod_coef <- prod_coef
  parms$dem_coef <- dem_coef
  parms$epsilon <- epsilon
  parms$weight_matrix <- weight_matrix
  parms$flow_fun <- function(t,x,parms) flow_fun0(x,parms)
  parms$prod_fun <- function(t,x,parms) {
    parms$prod_coef * ifelse(
      parms$prod_coef > 0,
      (1-x)/(1-x+parms$epsilon),
      0) 
  }
  parms$dem_fun <- function(t,x,parms) {
    parms$dem_coef * ifelse(
      parms$dem_coef > 0,
      x/(x+parms$epsilon),
      0) 
  }
  
  return(parms)
}

construct_weight_matrix_from_coefs_3x3 <- function(coef_vec){
  matrix(
    nrow = 3, ncol = 3, byrow=TRUE,
    data = c(
      0, - coef_vec['s12'], - coef_vec['s13'],
      coef_vec['s12']/coef_vec['a12'], 0, -coef_vec['s23'],
      coef_vec['s13']/(coef_vec['a12']*coef_vec['a23']), coef_vec['s23']/coef_vec['a23'], 0))
}


prepare_parms_approx_with_failure_period <- function(
    epsilon = 1e-4,
    prod_coef = c(1,0,0), dem_coef = c(0,0,1),
    flow_fun0 = function(x,parms) {
      x * c((parms$weight_matrix  * (parms$weight_matrix < 0)) %*% (1-x)) +
        (1-x) * c((parms$weight_matrix  * (parms$weight_matrix > 0)) %*% x)
    },
    weight_matrix = matrix(c(0,-1,-1,1,0,-1,1,1,0),ncol=3,byrow=TRUE),
    failure_timerange = c(1,2)
){
  parms <- list()
  parms$prod_coef <- prod_coef
  parms$dem_coef <- dem_coef
  parms$epsilon <- epsilon
  parms$weight_matrix <- weight_matrix
  parms$flow_fun <- function(t,x,parms) flow_fun0(x,parms)
  parms$failure_timerange <- failure_timerange
  parms$prod_fun <- function(t,x,parms) {
    if(t < parms$failure_timerange[1] | t > parms$failure_timerange[2]){
    r <- parms$prod_coef * ifelse(
      parms$prod_coef > 0,
      (1-x)/(1-x+parms$epsilon),
      0) 
    }else{
      r <- rep(0,length(parms$prod_coef))
    }
  }
  parms$dem_fun <- function(t,x,parms) {
    parms$dem_coef * ifelse(
      parms$dem_coef > 0,
      x/(x+parms$epsilon),
      0) 
  }
  
  return(parms)
}

construct_weight_matrix_from_coefs_3x3 <- function(coef_vec){
  matrix(
    nrow = 3, ncol = 3, byrow=TRUE,
    data = c(
      0, - coef_vec['s12'], - coef_vec['s13'],
      coef_vec['s12']/coef_vec['a12'], 0, -coef_vec['s23'],
      coef_vec['s13']/(coef_vec['a12']*coef_vec['a23']), coef_vec['s23']/coef_vec['a23'], 0))
}


generate_parm_set_for_demo <- function(
    layer_sizes = c(5,29,100),
    pref_matrix = matrix(
      c(0,1,1,
        0,0,1,
        0,0,0),
      byrow=TRUE,ncol=3),
    epsilon = 1e-6,
    flow_fun0 = function(t,u1,u2) u1 * (1-u2),
    prod_level = 18.9,
    capac = c(47.75,46.5,1.19),
    weights_range = c(1,1),
    weights_center = c(0.77,0.1,0.06),
    corrected_network = FALSE
){
  dem_level <- prod_level*layer_sizes[1]/layer_sizes[3]
  #-> INITIALIZATION OF PARAMETERS FOR EXACT SYSTEM
  parms_ex <- list(
    layer.sizes = layer_sizes,
    pref.matrix = pref_matrix,
    epsilon = epsilon,
    flow_fun0 = flow_fun0,
    capac = capac
  )
  parms_ex$N <- sum(parms_ex$layer.sizes)
  parms_ex$g <- sample_mfn_for_demo(parms_ex$layer.sizes,parms_ex$pref.matrix,weights_range = weights_range, weights_center = weights_center)  
  
  
  parms_ex$prod_coef <- as.numeric(V(parms_ex$g)$type == min(V(parms_ex$g)$type)) * prod_level
  parms_ex$dem_coef <- as.numeric(V(parms_ex$g)$type == max(V(parms_ex$g)$type)) * dem_level
  parms_ex <- prepare_parms_exact_for_demo(parms_ex)
  
  #END-> INITIALIZATION OF PARAMETERS FOR EXACT SYSTEM
  #-> INITIALIZATION OF PARAMETERS FOR APPROX SYSTEM
  parms_ap <- list(
    layer.sizes = c(1,1,1),
    pref.matrix = matrix(
      c(0,1,1,
        0,0,1,
        0,0,0),
      byrow=TRUE,ncol=3),
    epsilon = 1e-4,
    flow_fun0 = parms_ex$flow_fun0
  )
  parms_ap$N <- sum(parms_ap$layer.sizes)
  parms_ap$g <- sample_mfn(parms_ap$layer.sizes,parms_ap$pref.matrix)
  parms_ap$prod_coef <- c(1,0,0) * prod_level/capac[1]
  parms_ap$dem_coef <- c(0,0,1) * dem_level/capac[3]
  parms_ap <- prepare_parms_exact(parms_ap)
  eds <- ends(parms_ex$g,E(parms_ex$g))
  eds <- cbind(V(parms_ex$g)[eds[,1]]$type,V(parms_ex$g)[eds[,2]]$type)
  # TODO: Generalize this for an arbitrary number of layers
  parms_ap$nei_in_weight <- list(
    0,
    sum(E(parms_ex$g)$weight[eds[,1]==1 & eds[,2]==2])/parms_ex$layer.sizes[2],
    c(sum(E(parms_ex$g)$weight[eds[,1]==1 & eds[,2]==3]),sum(E(parms_ex$g)$weight[eds[,1]==2 & eds[,2]==3]))/parms_ex$layer.sizes[3])
  
  parms_ap$nei_out_weight <- list(
    c(sum(E(parms_ex$g)$weight[eds[,1]==1 & eds[,2]==2]),
      sum(E(parms_ex$g)$weight[eds[,1]==1 & eds[,2]==3]))/
      parms_ex$layer.sizes[1],
    sum(E(parms_ex$g)$weight[eds[,1]==2 & eds[,2]==3])/parms_ex$layer.sizes[2],
    0)
  #END-> INITIALIZATION OF PARAMETERS FOR APPROX SYSTEM
  return(list('exact'=parms_ex,'approx'=parms_ap))
}

sample_mfn_for_demo <- function(layer.sizes,pref.matrix,weights_range = c(1,1),weights_center = c(0.77,0.1,0.06)){
  n <- sum(layer.sizes)
  g <- sample_sbm(
    n = n,
    pref.matrix = pref.matrix,
    block.sizes = layer.sizes,
    directed = TRUE,
    loops = FALSE
  )
  V(g)$type <- unlist(lapply(1:length(layer.sizes),function(i) rep(i,layer.sizes[i])))
  eds <- ends(g,E(g))
  eds_type <- cbind(V(g)$type[eds[,1]],V(g)$type[eds[,2]]) 
  E(g)$type <- apply(eds_type,1,function(x){
    ifelse(all(x %in% c(1,2)),
           '12',
           ifelse(all(x %in% c(2,3)),
                  '23','13'))
  })
  E(g)$weight <- runif(ecount(g),weights_range[1],weights_range[2]) * 
    ifelse(
      E(g)$type == '12', 
      weights_center[1],
      ifelse(E(g)$type == '13',
             weights_center[2],
             weights_center[3]))
  layout <- lapply(1:length(layer.sizes),function(i) expand.grid(1:layer.sizes[i],i)) %>% bind_rows %>% as.matrix
  V(g)$x <- layout[,1]
  V(g)$y <- layout[,2]
  return(g)
}

prepare_parms_exact_for_demo <- function(parms){
  parms$prod_fun <- function(t,u,parms){
    ifelse(
      parms$prod_coef > 0,
      parms$prod_coef/parms$capac[1] * (1-u) / (1 + parms$epsilon - u),
      0)
  }
  parms$dem_fun <- function(t,u,parms){
    ifelse(
      parms$dem_coef > 0,
      parms$dem_coef/parms$capac[3] * u / ( parms$epsilon + u),
      0)
  }
  capac_per_node <- c(
    rep(parms$capac[1],parms$layer.sizes[1]),
    rep(parms$capac[2],parms$layer.sizes[2]),
    rep(parms$capac[3],parms$layer.sizes[3]))
  
  parms$nei_in <- lapply(
    1:parms$N,
    function(i) neighbors(parms$g,i,'in'))
  parms$nei_out <- lapply(
    1:parms$N,
    function(i) neighbors(parms$g,i,'out'))
  parms$nei_in <- lapply(
    1:parms$N,
    function(i) neighbors(parms$g,i,'in'))
  
  parms$nei_out_weight <- lapply(
    1:parms$N,
    function(i){
      ni <- parms$nei_out[[i]]
      if(length(ni) > 0){
        r <- E(parms$g)[get.edge.ids(parms$g,c(rbind(i,ni)))]$weight/capac_per_node[i]
      }else{
        r <- 0
      }
      return(r)
    }
  )
  
  parms$nei_in_weight <- lapply(
    1:parms$N,
    function(i){
      ni <- parms$nei_in[[i]]
      if(length(ni) > 0){
        r <- E(parms$g)[get.edge.ids(parms$g,c(rbind(ni,i)))]$weight/capac_per_node[i]
      }else{
        r <- 0
      }
      return(r)
    }
  )
  
  parms$flow_fun <- function(t,u,parms){
    flow_in <- sapply(1:length(u),function(i){
      if(length(unlist(parms$nei_in[[i]]))  > 0){
        r <- sum(
          parms$flow_fun0(
            t,
            u[parms$nei_in[[i]]],
            u[i]
          ) * parms$nei_in_wei[[i]]
        )
      }else{
        r <- 0
      }
      return(r)
    })
    flow_out <- sapply(1:length(u),function(i){
      if(length(parms$nei_out[[i]])  > 0){
        r <- sum(
          parms$flow_fun0(
            t,
            u[i],
            u[parms$nei_out[[i]]]) * 
            parms$nei_out_wei[[i]]
        )
      }else{
        r <- 0
      }
      return(r)
    })
    return(flow_in - flow_out)
  }
  return(parms)
}

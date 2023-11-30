rm(list=ls())
gc()
source('attack_network_helpers.R')

g <- read.graph('G_revised.gml',format='gml')
load('sea_levels_nodes.RData')

scenarios_names <- setdiff(colnames(df),c('id','longitude','latitude','type'))

g3N0 <- get_multiplex(g,FALSE)
eds <- ends(g3N0,E(g3N0))

node_amounts <- as.data.frame(
  table(
    c(
      V(g3N0)$type,
      paste(
        V(g3N0)$type[eds[,1]],
        V(g3N0)$type[eds[,2]],
        E(g3N0)$etype,
        sep='-'
      )
    )
  )
)

g_prod <- delete.vertices(
  g, 
  !V(g)$type %in% c('refinery','terminal','product')
)

fn <- add.vertices(g_prod,2) %>%
  add_edges(c(
    c(t(cbind(vcount(g_prod)+1,which(V(g_prod)$type=='refinery')))),
    c(t(cbind(vcount(g_prod)+2,which(V(g_prod)$type=='terminal')))))) %>% as.undirected()

node_amounts$scenario <- '0'
node_amounts$min_cut  <- min_cut(fn,vcount(g_prod)+1,vcount(g_prod)+2)
k <<- 1

node_amounts <- lapply(scenarios_names,function(sc){
  print(round(k/length(scenarios_names)*100,2))
  k <<- k+1
  scenario <- !is.na(df[,sc]) & df[,sc] >= 0.15
  g3N <- get_multiplex(g,scenario)
  eds <- ends(g3N,E(g3N))
  nm <- as.data.frame(
    table(
      c(
        V(g3N)$type,
        paste(
          V(g3N)$type[eds[,1]],
          V(g3N)$type[eds[,2]],
          E(g3N)$etype,
          sep='-'
        )
      )
    )
  )
  g_prod <- delete.vertices(
    g, 
    !V(g)$type %in% c('refinery','terminal','product') | scenario
  )
  
  fn <- add.vertices(g_prod,2) %>%
    add_edges(c(
      c(t(cbind(vcount(g_prod)+1,which(V(g_prod)$type=='refinery')))),
      c(t(cbind(vcount(g_prod)+2,which(V(g_prod)$type=='terminal')))))) %>% as.undirected()
  nm$scenario <- sc
  nm$min_cut <- min_cut(fn,vcount(g_prod)+1,vcount(g_prod)+2)
  return(nm)
}) %>% 
  bind_rows(node_amounts)

save(node_amounts,file='node_amounts_15cm.RData')

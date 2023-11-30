library(igraph)
library(dplyr)
library(tidygraph)

get_multiplex <- function(g,scenario){
  
  g <- delete.vertices(g,v = V(g)[scenario])
  
  LINK_ASSETS <- c('product','road')
  NODE_ASSETS <- c('gas_station','refinery','terminal')
  
  gNA <- induced_subgraph(g, V(g)$type %in% NODE_ASSETS)
  gNA <- as.undirected(gNA)
  
  for(la in LINK_ASSETS){
    
    eids <- grep(la,E(g)$type)
    types <- ends(g,eids,FALSE) %>% c %>% unique %>% {V(g)$type[.]} %>% unique
    
    
    g_ <- induced_subgraph(g,V(g)$type %in% types)
    g_ <- delete.vertices(g_, degree(g_) == 0)
    cl_ <- clusters(g_)
    rm_clu <- c()
    for(clu in 1:cl_$no){
      if(!la %in% V(g_)$type[cl_$membership == clu]){
        rm_clu <- c(rm_clu,clu)
      }else{
        to_comb <- V(g_)$id[cl_$membership == clu & (V(g_)$type %in% NODE_ASSETS)]
        if(length(to_comb) > 1){
          eseq <- combn(to_comb,2) %>% c
          gNA <- add_edges(gNA, match(eseq,V(gNA)$id),etype = la)
        }
      }
    }
  }
  
  g3N <- induced_subgraph(
    gNA, 
    V(gNA)$type %in% c('refinery','terminal','gas_station')) %>% 
    as_tbl_graph()
  return(g3N)
}

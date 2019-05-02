require(here)
source("pckgs_and_useful_wrappers.R")

#S. Cerevisiae graph
SC_net <- readxl::read_xlsx(
  path = here("SCerevisiae_data/PLOSONE_Gold_standard_networks.xlsx"),
  sheet = 3) %>% 
  set_names(c("from","to"))

SC_net <- graph_from_data_frame(SC_net)




#E. Coli graph
EC_edges <- readxl::read_xlsx(
  path = here("EColi_data/PNAS_Gold_standard_networks.xlsx"),
  sheet = 2) %>% 
  set_names(c("from","from_id","to","to_id",
              "sign","evidence","evidence_lvl")) %>%
  mutate(sign = ifelse(sign=='+',1,-1)) %>%
  separate_rows(from, from_id, sep = ";") %>%
  separate_rows(to, to_id, sep = ";") 
  

EC_nodes <-   EC_edges %>% 
  select(from,from_id) %>% 
  rename(gene = from, gene_id = from_id) %>%
  union_all(
    EC_edges %>% 
      select(to,to_id) %>% 
      rename(gene = to, gene_id = to_id)
  ) %>% 
  distinct()

EC_edges <- EC_edges %>% select(-from_id,-to_id)

EC_net <- graph_from_data_frame(
  d = EC_edges,
  directed = T,
  vertices = EC_nodes) %>% simplify()

rm(EC_edges,EC_nodes)




##make smaller graphs
source("marbach_algo.R")


#largest connected_component of SCerevisiae net
largest_connected_component_SC <- 
  which(
    components(SC_net)$membership == 
      which.max( components(SC_net)$csize )
  )
#SCerevisiae 10 nodes
set.seed(1204)
ini_node <- sample(1, x = largest_connected_component_SC)
SC_20 <- marbach_algo(n = 20, G = as.undirected(SC_net), initial_node = ini_node)
SC_20 <- induced_subgraph(SC_net, vids = SC_20)

#SCerevisiae 50 nodes
set.seed(1213)
ini_node <- sample(1, x = largest_connected_component_SC)
SC_50 <- marbach_algo(n = 50, G = as.undirected(SC_net), initial_node = ini_node)
SC_50 <- induced_subgraph(SC_net, vids = SC_50)

#SCerevisiae 200 nodes
set.seed(0909)
ini_node <- sample(1, x = largest_connected_component_SC)
SC_200 <- marbach_algo(n = 200, G = as.undirected(SC_net), initial_node = ini_node)
SC_200 <- induced_subgraph(SC_net, vids = SC_200)


#largest connected_component of EColi net
largest_connected_component_EC <- 
  which(
    components(EC_net)$membership == 
      which.max( components(EC_net)$csize )
  )
#EColi 10 nodes
set.seed(0720)
ini_node <- sample(1, x = largest_connected_component_EC)
EC_20 <- marbach_algo(n = 20, G = as.undirected(EC_net), initial_node = ini_node)
EC_20 <- induced_subgraph(EC_net, vids = EC_20)


#EColi 50 nodes
set.seed(1004)
ini_node <- sample(1, x = largest_connected_component_EC)
EC_50 <- marbach_algo(n = 50, G = as.undirected(EC_net), initial_node = ini_node)
EC_50 <- induced_subgraph(EC_net, vids = EC_50)


#EColi 200 nodes
set.seed(1104)
ini_node <- sample(1, x = largest_connected_component_EC)
EC_200 <- marbach_algo(n = 200, G = as.undirected(EC_net), initial_node = ini_node)
EC_200 <- induced_subgraph(EC_net, vids = EC_200)


rm(add_node_marbach,marbach_algo,modularity,i,ini_node,largest_connected_component_EC,largest_connected_component_SC)


for(i in ls()){
  saveRDS(get(i),here(paste0("Robjects/graph_objects/",i,".RDS")))
  #print(i)
  #print(components(get(i))$no)
}

#str_saveobj(ls(), subdir = "Robjects/graph_objects")
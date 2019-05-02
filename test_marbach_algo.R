
G<-grg.game(10000,0.5)
A<-get.adjacency(G)

subset_genes <- sample(1:nrow(A),size = 3)
membership_v <- rep(1,nrow(A))
membership_v[subset_genes] <- 2

system.time(modularity(subset_genes,A))
system.time(a<-igraph::modularity(G,membership_v))

add_edge_marbach(subset_genes, A,3)

marbach_algo(5,A)

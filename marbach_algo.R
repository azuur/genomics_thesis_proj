#VERSIONES CON MATR DE ADYACENCIA
# modularity <- function(subset_genes, A){
#   s_vec <- rep(-1,nrow(A))
#   s_vec[subset_genes] <- 1
#   p_vec <- rowSums(as.matrix(A))
#   B <- A - p_vec%*%t(p_vec)/(sum(A))
#   as.numeric(unlist(t(s_vec)%*%B%*%s_vec/(2*sum(A))))
# }
# 
# add_node_marbach <- function(subset_genes, adj_mat, n){
#   if(n<=0){ return(subset_genes) }
# 
#   candidates <- setdiff(
#     c(which(rowSums(adj_mat[subset_genes,])>0),which(colSums(adj_mat[,subset_genes])>0)), 
#       subset_genes)
#   modularity_vec <- sapply(X = candidates,
#                            FUN = function(x) modularity(c(subset_genes,x),adj_mat)
#                            )
#   add_node_marbach(c(subset_genes,candidates[which.max(modularity_vec)]), adj_mat, n-1)
# }
# 
# marbach_algo<- function(n, adj_mat, seed=NULL, initial_node = NULL){
#   if(!is.null(seed)){ set.seed(seed) }
#   subset_genes <- ifelse(is.null(initial_node),
#                          sample(1:nrow(adj_mat),size = 1),
#                          initial_node)
#   add_node_marbach(subset_genes, adj_mat, n-length(subset_genes))
# }








#VERSIONES CON OBJETO TIPO IGRAPH
modularity <- function(subset_genes, G){
  membership_v <- rep(1,vcount(G))
  membership_v[subset_genes] <- 2
  igraph::modularity(G, membership_v)
}

add_node_marbach <- function(subset_genes, G, n){
  if(n<=0){ return(subset_genes) }
  
  candidates <- setdiff(
    unique(unlist(
      ego(graph = G, order = 1, nodes = subset_genes, mode = "all")
    )),
    subset_genes
  )
  if(is.null(candidates)){ candidates <- set.diff(1:vcount(G),subset_genes)}
  modularity_vec <- sapply(X = candidates,
                           FUN = function(x) modularity(c(subset_genes,x), G)
  )
  add_node_marbach( c(subset_genes, candidates[which.max(modularity_vec)]), G, n-1)
}

marbach_algo<- function(n, G, seed=NULL, initial_node = NULL){
  if(!is.null(seed)){ set.seed(seed) }
  subset_genes <- ifelse(is.null(initial_node), 
                         sample(1:igraph::vcount(G),size = 1), 
                         initial_node)
  add_node_marbach(subset_genes, G, n-length(subset_genes))
}


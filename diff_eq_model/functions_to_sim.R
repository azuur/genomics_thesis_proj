# graph <- igraph::graph_from_isomorphism_class(3,number = 7,directed = T)
# idx <- 1
# r <- tao <- lambda.RNA <- lambda.protein <- 1
# b.values <- runif(8)
# c.values <- b.values/b.values


grn_diff_eq_clos <- 
  function(
    graph,idx,
    b.values = NULL, c.values=NULL, 
    tao = NULL, r = NULL,
    lambda.RNA = NULL, lambda.protein = NULL
  ){
    
    #differential equation for protein synth + degradation
    diff_eq_protein <- function(
      RNA,protein
    ){
      res <- rep(0, length(RNA))
      res[idx] <- r*RNA[idx] - lambda.protein*protein[idx]
      res
    }
    
    #regulators (parents, subsets of parents)
    parents_idx <- igraph::ego(graph,order = 1,nodes = idx, mode = "in")
    parents_idx <- as.numeric(parents_idx[[1]])
    reg_sets <- do.call(
      c, 
      lapply(0:length(parents_idx), function(x) combn(parents_idx,x,simplify = F))
    )
    
    #probability of RNA polymerase binding to promoter
    prob_bind <- function(protein){
      sum(b.values*sapply(reg_sets, function(idx) prod(protein[idx])) )/
        sum(c.values*sapply(reg_sets, function(idx) prod(protein[idx])) )
    }
    
    #differential equation for RNA transcription + degradation
    diff_eq_RNA <- function(
      RNA,protein
    ){
      res <- rep(0, length(RNA))
      res[idx] <- tao*prob_bind(protein) - lambda.RNA*RNA[idx]
      res
    }
    
    return(list(diff_eq_RNA = diff_eq_RNA, diff_eq_protein = diff_eq_protein))
    
  }


all_grn_diff_eq_clos <- 
  function(
    graph,
    list.of.b.values = NULL, list.of.c.values=NULL, 
    list.of.tao = NULL, list.of.r = NULL,
    list.of.lambda.RNA = NULL, list.of.lambda.protein = NULL
  ){
    
    result <-
      lapply(1:igraph::vcount(graph),
             function(i) grn_diff_eq_clos(
               graph,i,
               b.values = list.of.b.values[[i]], c.values=list.of.c.values[[i]], 
               tao = list.of.tao[[i]], r = list.of.r[[i]],
               lambda.RNA = list.of.lambda.RNA[[i]], lambda.protein = list.of.lambda.protein[[i]])
      )
    
    return(result)
    
  }


simulate_diff_eq <- function(
  all_diff_eq, 
  RNA, protein, delta.t = 0.1, timesteps = 1L, 
  noise.sd = 0,
  knockout.function = NULL, path = F, history = list(RNA = NULL, protein = NULL)){
  
  if(timesteps < 1){
    if(path == F){
      return(list(RNA = RNA, protein = protein))
    } else {
      return(history)
    }
  }
  
  if(!is.null(knockout.function)){
    controlled <- knockout.function(RNA,protein)
    RNA <- controlled$RNA
    protein <- controlled$protein
  }
  
  RNA_diff_eq <- lapply(all_diff_eq, function(l) l[[1]])
  protein_diff_eq <- lapply(all_diff_eq, function(l) l[[2]])
  
  delta.RNA <- lapply(seq_along(RNA_diff_eq),function(i) RNA_diff_eq[[i]](RNA,protein)  )
  delta.protein <- lapply(seq_along(protein_diff_eq),function(i) protein_diff_eq[[i]](RNA,protein)  )
  
  delta.RNA <- colSums(do.call(rbind,delta.RNA))*delta.t
  delta.protein <- colSums(do.call(rbind,delta.protein))*delta.t
  
  RNA1 <- RNA + delta.RNA + rnorm(n = length(RNA),mean = 0,sd = noise.sd*delta.t)
  protein1 <- protein + delta.protein + rnorm(n = length(protein),mean = 0,sd = noise.sd*delta.t)
  
  if(!is.null(knockout.function)){
    controlled <- knockout.function(RNA1,protein1)
    RNA1 <- controlled$RNA
    protein1 <- controlled$protein
  }
  
  if(path == T){
    history <- list(RNA = rbind(history$RNA,RNA1),protein = rbind(history$protein,protein1))
  }
  
  simulate_diff_eq(all_diff_eq, RNA1, protein1, delta.t, timesteps - 1, noise.sd, knockout.function, path, history)
  
}









genes <- c("O","S","N","C","Gc","G")

ADJ <- matrix(0,6,6)

ADJ[c(1,2,3,4,5),1] <- 1
ADJ[c(1,2,3),2] <- 1
ADJ[c(1,2,3,6),3] <- 1
ADJ[c(1,4),4] <- 1
ADJ[c(4,6),5] <- 1
ADJ[c(1,3,6),6] <- 1

colnames(ADJ) <- rownames(ADJ) <-
  genes

graph <- igraph::graph_from_adjacency_matrix(ADJ)

plot(graph)

reg.values.fun <- function(graph,idx,regulators = NULL,values = NULL){
  parents <- igraph::ego(graph,order = 1,nodes = idx, mode = "in")
  parents <- names(parents[[1]])
  reg_sets <- do.call(
    c, 
    lapply(0:length(parents), function(x) combn(parents,x,simplify = F))
  )
  
  vals <- rep(0, length(reg_sets))
  reg_sets_idx <- sapply(regulators, function(x) which(sapply(reg_sets, function(y) setequal(x,y))))
  
  vals[reg_sets_idx] <- values
  vals
}

A.value <- 10

b.regulators<-list(
  list(NULL,c("O","S"),c("O","S","N")),
  list(NULL,c("O","S"),c("O","S","N")),
  list(NULL,c("O","S"),c("O","S","N")),
  list(NULL,c("C")),
  list(NULL,c("C"),c("G")),
  list(NULL,c("O"),c("G"))
)

b.values <- list(
  c(0.001+A.value,0.005,0.025),
  c(0.001,0.005,0.025),
  c(0.001,0.1,0.1),
  c(0.001,2),
  c(0.001,0.1,0.1),
  c(0.1,1,0.00025)
)  

c.regulators<-list(
  list(NULL,c("O"),c("O","S"),c("O","S","N"),c("O","C"),c("Gc")),
  list(NULL,c("O"),c("O","S"),c("O","S","N")),
  list(NULL,c("O"),c("O","S"),c("O","S","N"),c("O","G")),
  list(NULL,c("C"),c("O","C")),
  list(NULL,c("C"),c("G")),
  list(NULL,c("O"),c("G"),c("N"))
)

c.values <- list(
  c(1+A.value,0.001,0.005,0.025,10,10),
  c(1,0.001,0.005,0.025),
  c(1,0.001,0.1,0.1,10),
  c(1,2,5),
  c(1,0.1,0.1),
  c(1,1,0.00025,15)
)  


b.val.list <-
  lapply(seq_along(b.regulators), function(i) reg.values.fun(graph,i,b.regulators[[i]],b.values[[i]]))

c.val.list <-
  lapply(seq_along(c.regulators), function(i) reg.values.fun(graph,i,c.regulators[[i]],c.values[[i]]))




modelo <- all_grn_diff_eq_clos(
  graph,
  b.val.list,c.val.list,
  rep(list(1),6),rep(list(1),6),
  rep(list(0.1),6),rep(list(0.1),6)
)



knockout.function <- function(RNA,protein){
  list(RNA = RNA, protein = RNA)
}


#ini_RNA <- c(0.4,7,7,1.9,1.5,0.1)

#ini_RNA <- runif(6,3,7)

#ini_RNA <- rep(0,6)

ini_RNA <- c(5,5,5,0.1,0.1,0.1)

ini_RNA <- c(2.1181439,7.3786152,7.0095017,1.5095752,1.4720853,0.2048936)


# simulate_diff_eq(modelo, ini_RNA, ini_RNA, delta.t = 0.01,timesteps = 500,
#                  knockout.function = knockout.function, path = T)

res <- list(RNA = matrix(ini_RNA,nrow = 1), protein = matrix(ini_RNA,nrow = 1))
for(i in 1:20){
  res <- simulate_diff_eq(modelo, 
                          pmax(
                            res$RNA[nrow(res$RNA),]+
                              #(which.max(runif(6))==c(1:6))*
                              #c(0,0,0,1,0,0)*
                              #rnorm(6,0,res$RNA[nrow(res$RNA),]),
                              rnorm(6,0,1),
                            0),
                          res$RNA[nrow(res$RNA),],
                          delta.t = 0.1,timesteps = 1000,noise.sd = 0,#0.4/2^i,
                          knockout.function = knockout.function, path = T,history = res)
}
res1 <- list(RNA = matrix(ini_RNA,nrow = 1), protein = matrix(ini_RNA,nrow = 1))
for(i in 1:5){
  res1 <- simulate_diff_eq(modelo, 
                           res1$RNA[nrow(res1$RNA),], 
                           res1$RNA[nrow(res1$RNA),],
                           delta.t = 0.1,timesteps = 1000,noise.sd = 0,
                           knockout.function = knockout.function, path = T,history = res1)
}


res$RNA[nrow(res$RNA),]
res1$RNA[nrow(res1$RNA),]

{
  plot(1:nrow(res$RNA), res$RNA[,1],type = "l",ylim = c(0,10))
  lines(1:nrow(res$RNA), res$RNA[,2])
  lines(1:nrow(res$RNA), res$RNA[,3])
  lines(1:nrow(res$RNA), res$RNA[,4])
  lines(1:nrow(res$RNA), res$RNA[,5])
  lines(1:nrow(res$RNA), res$RNA[,6])
  # lines(1:nrow(res1$RNA), res1$RNA[,1])
  # lines(1:nrow(res1$RNA), res1$RNA[,2])
  # lines(1:nrow(res1$RNA), res1$RNA[,3])
  # lines(1:nrow(res1$RNA), res1$RNA[,4])
  # lines(1:nrow(res1$RNA), res1$RNA[,5])
  # lines(1:nrow(res1$RNA), res1$RNA[,6])
}  


plot(1:1000, res$RNA[4002:5001,1],type = "l",ylim = c(0,10))
lines(1:1000, res$RNA[4002:5001,2])
lines(1:1000, res$RNA[4002:5001,3])
lines(1:1000, res$RNA[4002:5001,4])
lines(1:1000, res$RNA[4002:5001,5])
lines(1:1000, res$RNA[4002:5001,6])

res$RNA[3998:4010,]





elast <- 10000

(0.5*1^(1-elast) + 0.5*1.5^(1-elast))^(1/(1-elast))

require(here)
source("pckgs_and_useful_wrappers.R")


#poner opción isnull
readobj(subdir="Robjects/graph_objects")


source("BN_simulator.R")
source("types_of_BN_functions.R")
# 
# set.seed(10)
# sim_info <- lapply(X = 1:200, FUN = create_additive_T1_params)
# functional_forms_list <- do.call(c, lapply(X = sim_info, function(l) l[[1]]))
# param_list = do.call(c, lapply(X = sim_info, function(l) list(l[[2]])))
# all_g <- all_gene_closures(SC_200, 
#                               functional_forms_list = functional_forms_list,
#                               param_list = param_list)


# {simulator1 <- sim_closure(SC_200, all_g)
# errores <- runif(200*5000,-1,1) %>% matrix(nrow=200)
# sims <- apply(errores, MARGIN = 2, simulator1$f)
# 
# system.time(bnlearn::pc.stable(x = as.data.frame(t(errores)[,1:70])))
# }
#EN REDES PEQUEÑAS
#Decisión 0. Número de réplicas
#10mil?
#Decision 1. Tamaño de muestra (5 casos)
#n=20,50,100,500,1000
#Decisión 2. Errores (6 casos)
#normales, uniformes. 1-r2 = 0.3, 0.5, 0.8
#Decisión 3. Tipos de relación (3 casos)
#solo lineales, 
#solo relU/sigmoide/otra-monótona, 
#por tercios al azar (lineal, monótona, raraza)


#REDES GRANDES
#Decisión 0. Número de réplicas
#10mil?
#Decision 1. Tamaño de muestra (4 casos)
#n=50,100,500,1000
#Decisión 2. Errores (2 casos)
#normales, uniformes. varianza chica, varianza media, varianza grande
#Decisión 3. Tipos de relación (2 casos)
#solo lineales/relU/sigmoide/otra-monótona, 
#por tercios al azar (monótonas, rarazas)

#GUARDAR
#puntajes
#redes resultantes
## ~ ~ ~ idea... en PC, calcular tasas de falsos positivos bajo indep, y ajustar alpha con eso
#tiempo de


#LINEALES solo
netw <- SC_200
n <- vcount(netw)
functional_forms <- rep(list("linear"), n)
parameters <- rep(list(list(coef_distribution = function(n) runif(n, -1, 1))),n)

set.seed(10)
all_g <- all_gene_closures(netw, 
                           functional_forms_list = functional_forms,
                           param_list = parameters)

simulator <- sim_closure(netw, all_g)

# A<-sapply(simulator$metadata, 
#        function(l){
#          rw <- rep(0,n)
#          rw[l$parent_nodes] <- l$coefficients
#          rw
#          }
#        ) %>% t()
# B_inv <- solve(diag(n)-A)
# multiplicadores_diag <-sapply(1:nrow(B_inv), function(x) sum(B_inv[x,]^2))
# summary(1/multiplicadores_diag)
# c(1/multiplicadores_diag) %>%  split(rowSums(A!=0)) %>% lapply(summary)

# 
# errores <- runif(200*5000,-1,1) %>% matrix(nrow=200)
# sims <- apply(errores, MARGIN = 2, simulator$f)
# hist(apply(sims,1,var) - multiplicadores_diag*(1-(-1))^2/12)




#https://math.stackexchange.com/questions/47543/getting-the-inverse-of-a-lower-upper-triangular-matrix
A<-sapply(simulator$metadata, 
          function(l){
            rw <- rep(0,n)
            rw[l$parent_nodes] <- l$coefficients
            rw
          }
) %>% t()

ord <- as.numeric(topo_sort(netw,mode = "out"))
# A <- t(as.matrix(get.adjacency(netw)))[ord,ord]
# sum(t(A)%*%A)

r<-0.3
edges <- igraph::get.edgelist(netw, names=F)
vs<-rep(1/3,length(ord))
for(x in ord){
  parents_x <- edges[edges[,2]==x,1]
  if(length(parents_x)>0){
    B_inv <- solve(diag(n)-A)
    vs[x] <- (r/(1-r))*(A[x,parents_x]%*%
      (B_inv%*%diag(vs)%*%t(B_inv))[parents_x,parents_x]%*%
      A[x,parents_x])
  }
}

b <- sqrt(3*vs)
set.seed(10)
errores <- sapply(1:1000, function(x)runif(n,min = -b,max = b)) %>% t()
vs/apply(errores, MARGIN = 2,FUN = var)

y <- apply(errores, MARGIN = 1, simulator$f) %>% t()
vs/apply(y, MARGIN = 2,FUN = var)
summary(vs/apply(y, MARGIN = 2,FUN = var))

y1<-y
saveobj(y1)




























#Solo sigmoides
netw <- SC_200
n <- vcount(netw)
functional_forms <- rep(list("sigmoid"), n)
parameters <- rep(list(list(coef_distribution = function(n) runif(n, -1, 1))),n)

set.seed(10)
all_g <- all_gene_closures(netw, 
                           functional_forms_list = functional_forms,
                           param_list = parameters)

simulator <- sim_closure(netw, all_g)

# A<-sapply(simulator$metadata, 
#        function(l){
#          rw <- rep(0,n)
#          rw[l$parent_nodes] <- l$coefficients
#          rw
#          }
#        ) %>% t()
# B_inv <- solve(diag(n)-A)
# multiplicadores_diag <-sapply(1:nrow(B_inv), function(x) sum(B_inv[x,]^2))
# summary(1/multiplicadores_diag)
# c(1/multiplicadores_diag) %>%  split(rowSums(A!=0)) %>% lapply(summary)

# 
# errores <- runif(200*5000,-1,1) %>% matrix(nrow=200)
# sims <- apply(errores, MARGIN = 2, simulator$f)
# hist(apply(sims,1,var) - multiplicadores_diag*(1-(-1))^2/12)




#https://math.stackexchange.com/questions/47543/getting-the-inverse-of-a-lower-upper-triangular-matrix
A<-sapply(simulator$metadata, 
          function(l){
            rw <- rep(0,n)
            rw[l$parent_nodes] <- l$coefficients
            rw
          }
) %>% t()

ord <- as.numeric(topo_sort(netw,mode = "out"))
# A <- t(as.matrix(get.adjacency(netw)))[ord,ord]
# sum(t(A)%*%A)

r<-0.3
edges <- igraph::get.edgelist(netw, names=F)
vs<-rep(1/3,length(ord))
for(x in ord){
  parents_x <- edges[edges[,2]==x,1]
  if(length(parents_x)>0){
    B_inv <- solve(diag(n)-A/2)
    vs[x] <- (1/4)*(r/(1-r))*(A[x,parents_x]%*%
                          (B_inv%*%diag(vs)%*%t(B_inv))[parents_x,parents_x]%*%
                          A[x,parents_x])
  }
}

b <- sqrt(3*vs)
set.seed(10)
errores <- sapply(1:1000, function(x)runif(n,min = -b,max = b)) %>% t()
vs/apply(errores, MARGIN = 2,FUN = var)

y <- apply(errores, MARGIN = 1, simulator$f) %>% t()
vs/apply(y, MARGIN = 2,FUN = var)
summary(vs/apply(y, MARGIN = 2,FUN = var))
hist(vs/apply(y, MARGIN = 2,FUN = var))

y2<-y
saveobj(y2)

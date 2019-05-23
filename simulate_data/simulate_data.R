require(here)
source("pckgs_and_useful_wrappers.R")
require(foreach)

#poner opción isnull
readobj(subdir="Robjects/graph_objects")

source("BN_simulator.R")
source("types_of_BN_functions.R")




graph_options <- c(#"EC_10","EC_20","EC_50","EC_200",
  "SC_10","SC_20","SC_50","SC_200")
type_options <- c("linear","sigmoid")
noise_options <- c("uniform","gaussian")
r_options <- c(0.2,0.5,0.8)
sample_size_options <- c(20, 50, 100, 500, 1000)
n_sim <- 1e3



# dir.create(here::here("Robjects","simulated_data"),showWarnings = T)
# for(netw_str in graph_options){
#   dir.create(here::here("Robjects","simulated_data",netw_str), showWarnings = T)
#   for(type in type_options){
#     dir.create(here::here("Robjects","simulated_data",netw_str,type), showWarnings = T)
#     for(noise in noise_options){
#       dir.create(here::here("Robjects","simulated_data",netw_str,type,noise), showWarnings = T)
#       for(r in r_options){
#         dir.create(here::here("Robjects","simulated_data",netw_str,type,noise,paste0("noise_to_sig_",r)), showWarnings = T)
#         for(sample_size in sample_size_options){
#           dir.create(here::here("Robjects","simulated_data",netw_str,type,noise,paste0("noise_to_sig_",r),paste0("sample_size_",sample_size)), showWarnings = T)
#         }
#       }
#     }
#   }
# }


# netw_str <- graph_options[1]
# type <- type_options[1]
# noise <- noise_options[1]
# r <- r_options[1]
# sample_size <- sample_size_options[1]



iter <- 0

for(sample_size in sample_size_options){
  print(paste0("sample_size=",sample_size))
  for(netw_str in graph_options){
    print(paste0("netw_str=",netw_str))
    for(type in type_options){
      print(paste0("type=",type))
      for(noise in noise_options){
        print(paste0("noise=",noise))
        for(r in r_options){
          print(paste0("r=",r))
          
          {
            iter <- iter + 1
            str_readobj(netw_str,"Robjects/graph_objects")
            netw <- get(netw_str)
            
            #simulador
            n <- vcount(netw)
            functional_forms <- rep(list(type), n)
            parameters <- rep(list(list(coef_distribution = function(n) runif(n, -1, 1))), n)
            
            set.seed(1308+iter)
            all_g <- all_gene_closures(netw, 
                                       functional_forms_list = functional_forms,
                                       param_list = parameters)
            
            simulator <- sim_closure(netw, all_g)
            
            
            #calibrar varianza
            W<-sapply(simulator$metadata, 
                      function(l){
                        rw <- rep(0,n)
                        rw[l$parent_nodes] <- l$coefficients
                        rw
                      }
            ) %>% t()
            
            vs <- calibrate_variances(netw, W, r, ini_var = 1, type = type)
            
            set.seed(1204+iter)
            random_order <- sample(n)
            col_names <- paste0("X",random_order)
            
            
            
            
            cl <- parallel::makeForkCluster(5)
            doParallel::registerDoParallel(cl)
            
            
            semilla <- 0.5
            path <- here::here(
              "Robjects","simulated_data",
              netw_str,
              type,
              paste0(noise,"_noise"),
              paste0(r,"_noise_to_sig"),
              paste0("sample_size_",sample_size)
            )
            dir.create(path, showWarnings = T, recursive = T)
            
            foreach(i = 1:n_sim) %dopar% {
              semilla <- 3.7*semilla*(1-semilla)
              set.seed(semilla)
              
              if(noise == "uniform"){
                b <- sqrt(3*vs)
                errores <- sapply(1:sample_size, 
                                  function(x)runif(n,min = -b,max = b)
                ) %>% t()
              }
              
              if(noise == "gaussian"){
                standdev <- sqrt(vs)
                errores <- sapply(1:sample_size, 
                                  function(x)rnorm(n,mean = 0, sd =vs)) %>% t()
              }
              
              assign(paste0("sim",i),
                     apply(errores, MARGIN = 1, simulator$f) %>% 
                       t() %>% 
                       `[`(, random_order) %>%
                       `colnames<-`(col_names),
                     envir = environment()
              )
              #vs/apply(errores, MARGIN = 2,FUN = var)
              #median(vs[random_order]/apply(get(paste0("sim",i)), MARGIN = 2,FUN = var))
              # if(abs(r-median(vs[random_order]/apply(get(paste0("sim",i)), MARGIN = 2,FUN = var))) >0.1 |
              #    abs(1-median(vs/apply(errores, MARGIN = 2,FUN = var)))>0.1){
              #   warning("Something fishy with variances of errors/sims.")
              # }
              
              saveRDS(get(paste0("sim",i),envir = environment()),
                      paste0(path,"/sim",i,".RDS")
                      
              )
              
              # 
              # str_saveobj(paste0("sim",i),
              #             subdir = 
              #               paste0(
              #                 "Robjects/simulated_data/",
              #                 netw_str,"/",
              #                 type,"/",
              #                 noise,"_noise/",
              #                 r,"_noise_to_sig/",
              #                 "sample_size",sample_size
              #               )
              # )
              # rm(list=paste0("sim",i))
              # 
              
              #gc()
              
            }
          }
        }
      }
    }
  }
}


# 
# 
# 
# 
# 
# #simulador
# netw <- SC_200
# n <- vcount(netw)
# functional_forms <- rep(list("linear"), n)
# parameters <- rep(list(list(coef_distribution = function(n) runif(n, -1, 1))),n)
# 
# set.seed(10)
# all_g <- all_gene_closures(netw, 
#                            functional_forms_list = functional_forms,
#                            param_list = parameters)
# 
# simulator <- sim_closure(netw, all_g)
# 
# 
# #calibrar varianza
# A<-sapply(simulator$metadata, 
#           function(l){
#             rw <- rep(0,n)
#             rw[l$parent_nodes] <- l$coefficients
#             rw
#           }
# ) %>% t()
# 
# ord <- as.numeric(topo_sort(netw,mode = "out"))
# 
# r<-0.3
# edges <- igraph::get.edgelist(netw, names=F)
# vs<-rep(1,length(ord))
# for(x in ord){
#   parents_x <- edges[edges[,2]==x,1]
#   if(length(parents_x)>0){
#     B_inv <- solve(diag(n)-A)
#     vs[x] <- (r/(1-r))*(A[x,parents_x]%*%
#       (B_inv%*%diag(vs)%*%t(B_inv))[parents_x,parents_x]%*%
#       A[x,parents_x])
#   }
# }
# 
# b <- sqrt(3*vs)
# set.seed(10)
# errores <- sapply(1:1000, function(x)runif(n,min = -b,max = b)) %>% t()
# vs/apply(errores, MARGIN = 2,FUN = var)
# 
# y <- apply(errores, MARGIN = 1, simulator$f) %>% t()
# vs/apply(y, MARGIN = 2,FUN = var)
# 
# vs/apply(errores, MARGIN = 2,FUN = var)
# summary(vs/apply(y, MARGIN = 2,FUN = var))
# summary(vs/apply(errores, MARGIN = 2,FUN = var))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #Solo sigmoides
# netw <- SC_200
# n <- vcount(netw)
# functional_forms <- rep(list("sigmoid"), n)
# parameters <- rep(list(list(coef_distribution = function(n) runif(n, -1, 1))),n)
# 
# set.seed(10)
# all_g <- all_gene_closures(netw, 
#                            functional_forms_list = functional_forms,
#                            param_list = parameters)
# 
# simulator <- sim_closure(netw, all_g)
# 
# # A<-sapply(simulator$metadata, 
# #        function(l){
# #          rw <- rep(0,n)
# #          rw[l$parent_nodes] <- l$coefficients
# #          rw
# #          }
# #        ) %>% t()
# # B_inv <- solve(diag(n)-A)
# # multiplicadores_diag <-sapply(1:nrow(B_inv), function(x) sum(B_inv[x,]^2))
# # summary(1/multiplicadores_diag)
# # c(1/multiplicadores_diag) %>%  split(rowSums(A!=0)) %>% lapply(summary)
# 
# # 
# # errores <- runif(200*5000,-1,1) %>% matrix(nrow=200)
# # sims <- apply(errores, MARGIN = 2, simulator$f)
# # hist(apply(sims,1,var) - multiplicadores_diag*(1-(-1))^2/12)
# 
# 
# 
# 
# #https://math.stackexchange.com/questions/47543/getting-the-inverse-of-a-lower-upper-triangular-matrix
# A<-sapply(simulator$metadata, 
#           function(l){
#             rw <- rep(0,n)
#             rw[l$parent_nodes] <- l$coefficients
#             rw
#           }
# ) %>% t()
# 
# ord <- as.numeric(topo_sort(netw,mode = "out"))
# # A <- t(as.matrix(get.adjacency(netw)))[ord,ord]
# # sum(t(A)%*%A)
# 
# r<-0.3
# edges <- igraph::get.edgelist(netw, names=F)
# vs<-rep(1,length(ord))
# for(x in ord){
#   parents_x <- edges[edges[,2]==x,1]
#   if(length(parents_x)>0){
#     B_inv <- solve(diag(n)-A/2)
#     vs[x] <- (1/4)*(r/(1-r))*(A[x,parents_x]%*%
#                           (B_inv%*%diag(vs)%*%t(B_inv))[parents_x,parents_x]%*%
#                           A[x,parents_x])
#   }
# }
# 
# b <- sqrt(3*vs)
# set.seed(10)
# errores <- sapply(1:1000, function(x)runif(n,min = -b,max = b)) %>% t()
# vs/apply(errores, MARGIN = 2,FUN = var)
# 
# y <- apply(errores, MARGIN = 1, simulator$f) %>% t()
# vs/apply(y, MARGIN = 2,FUN = var)
# summary(vs/apply(y, MARGIN = 2,FUN = var))
# hist(vs/apply(y, MARGIN = 2,FUN = var))
# 
# y2<-y
# saveobj(y2)

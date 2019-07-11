require(here)
source("pckgs_and_useful_wrappers.R")
require(foreach)

#poner opción isnull
readobj(subdir="Robjects/graph_objects")

source("BN_simulator.R")
source("types_of_BN_functions.R")




graph_options <- c(#"EC_10","EC_20","EC_50","EC_200",
  #"SC_10",
  #"SC_20",
  "SC_50"#,
  #"SC_200"
  )
type_options <- c("linear"#,
                  #"sigmoid"
                  )
noise_options <- c(#"uniform",
                   "gaussian")
r_options <- c(0.2,
               #0.5,
               0.8
               )
sample_size_options <- c(20
                         #, 50, 100, 500
                         )#,c(200)#
n_sim <- 1e3


# 
netw_str <- graph_options[1]
type <- type_options[1]
noise <- noise_options[1]
r <- r_options[2]
sample_size <- sample_size_options[1]


for(netw_str in graph_options){

  seed <- make.seed(netw_str)
  
  str_readobj(netw_str,"Robjects/graph_objects")
  netw <- get(netw_str)
  
  for(sample_size in sample_size_options){
    
    for(type in type_options){
      #simulador
      n <- vcount(netw)
      functional_forms <- rep(list(type), n)
      parameters <- rep(list(list(coef_distribution = function(n) runif(n, -1, 1))), n)
      
      set.seed(1308+seed)
      all_g <- all_gene_closures(netw, 
                                 functional_forms_list = functional_forms,
                                 param_list = parameters)
      
      simulator <- sim_closure(netw, all_g)
      
      for(noise in noise_options){
        
        for(r in r_options){
          print(paste0("sample_size=",sample_size))
          print(paste0("netw_str=",netw_str))
          print(paste0("type=",type))
          print(paste0("noise=",noise))
          print(paste0("r=",r))
          
          {
            seed <- make.seed(paste0(sample_size,netw_str,type,noise,r))
            
            
            #calibrar varianza
            W<-sapply(simulator$metadata, 
                      function(l){
                        rw <- rep(0,n)
                        rw[l$parent_nodes] <- l$coefficients
                        rw
                      }
            ) %>% t()
            
            vs <- calibrate_variances(netw, W, r, ini_var = 1, type = type)
            
            set.seed(1204+seed)
            random_order <- sample(n)
            col_names <- paste0("X",random_order)
            
            semilla <- seed/10^ceiling(log(seed,10))
            path <- file.path("/media","adrian","bodega","thesis",
              "Robjects","simulated_data",
              netw_str,
              type,
              paste0(noise,"_noise"),
              paste0(r,"_noise_to_sig"),
              paste0("sample_size_",sample_size)
            )
            dir.create(path, showWarnings = T, recursive = T)
            
            
            cl <- parallel::makeForkCluster(3)
            doParallel::registerDoParallel(cl)
            parallel::clusterSetRNGStream(cl = cl, semilla)
            
            foreach(i = 1:n_sim) %dopar% {
              # semilla <- 3.7*semilla*(1-semilla)
              # set.seed(semilla)
              
              if(noise == "uniform"){
                b <- sqrt(3*vs)
                errores <- sapply(1:sample_size, 
                                  function(x)runif(n,min = -b,max = b)
                ) %>% t()
              }
              
              if(noise == "gaussian"){
                standdev <- sqrt(vs)
                errores <- sapply(1:sample_size, 
                                  function(x)rnorm(n,mean = 0, sd = standdev)) %>% t()
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
                      paste0(file.path(path,"sim"),i,".RDS")
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
              rm(list=paste0("sim",i))
              # 
              
              #gc()
              
            }
            
            parallel::stopCluster(cl)
            rm(cl)
            gc()
          }
        }
      }
    }
  }
}


#file.path()

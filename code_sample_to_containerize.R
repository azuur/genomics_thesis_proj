require(here)
source("pckgs_and_useful_wrappers.R")
require(foreach)
require(parallel)
require(doParallel)
require(bnlearn)
require(pcalg)
require(ParallelPC)
source("narromi.R")

#poner opción isnull
readobj(subdir="Robjects/graph_objects")

source("BN_simulator.R")
source("types_of_BN_functions.R")


source("other_dependencies.R")




graph_options <- c("SC_20")
type_options <- c("sigmoid")
noise_options <- c("gaussian")
r_options <- c(0.8)
sample_size_options <- c(100)
n_sim <- 5



for(sample_size in sample_size_options){
  
  for(netw_str in graph_options){
    
    for(type in type_options){
      
      for(noise in noise_options){
        
        for(r in r_options){
          print(paste0("sample_size=",sample_size))
          print(paste0("netw_str=",netw_str))
          print(paste0("type=",type))
          print(paste0("noise=",noise))
          print(paste0("r=",r))
          
          {
            seed <- make.seed(paste0(sample_size,netw_str,type,noise,r))
            str_readobj(netw_str,"Robjects/graph_objects")
            netw <- get(netw_str)
            
            #simulador
            n <- vcount(netw)
            functional_forms <- rep(list(type), n)
            parameters <- rep(list(list(coef_distribution = function(n) runif(n, -1, 1))), n)
            
            set.seed(1308+seed)
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
            
            set.seed(1204+seed)
            random_order <- sample(n)
            col_names <- paste0("X",random_order)
            
            
            semilla <- seed/10^ceiling(log(seed,10))
            path <- here::here("/results","data",
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
          
          
          algoritmos <- list(
            mi_splines = function(x){
              ord <- min(as.integer(nrow(x)^(1/3)),4)
              fastGeneMI::get.mim.bspline(x, 
                                          order = ord, 
                                          n.cores = 5)
            },
            mi_mm = function(x)fastGeneMI::get.mim.MM(x,discretisation = "equalwidth", n.cores = 5),
            MRNET_splines = minet::mrnet,
            CLR_splines = parmigene::clr,
            ARACNE_splines = parmigene::aracne.a,
            MRNET_mm = minet::mrnet,
            CLR_mm = parmigene::clr,
            ARACNE_mm = parmigene::aracne.a,
            NARROMI = function(x){
              cl <- makeCluster(5)
              res <- NARROMI(x, cl = cl)
              stopCluster(cl)
              rm(cl)
              res},
            GENIE3 = function(x)GENIE3::GENIE3(exprMatrix = t(x),nCores = 5),
            TIGRESS = function(x)tigress::tigress(expdata = x, usemulticore = T)
          ) 
          
          data_path <- path
          
          
          estimates_path <- here::here("/results","estimates",
                                      netw_str,
                                      type,
                                      paste0(noise,"_noise"),
                                      paste0(r,"_noise_to_sig"),
                                      paste0("sample_size_",sample_size)
          )
          
          for(alg in seq_along(algoritmos)){
            dir.create(file.path(estimates_path,names(algoritmos)[alg]), showWarnings = T, recursive = T)
          }
          algo_times <- data.frame(
            times_mi_splines = rep(as.numeric(NA), n_sim),
            times_mi_mm = rep(as.numeric(NA), n_sim),
            times_MRNET_splines = rep(as.numeric(NA), n_sim),
            times_CLR_splines = rep(as.numeric(NA), n_sim),
            times_ARACNE_splines = rep(as.numeric(NA), n_sim),
            times_MRNET_mm = rep(as.numeric(NA), n_sim),
            times_CLR_mm = rep(as.numeric(NA), n_sim),
            times_ARACNE_mm = rep(as.numeric(NA), n_sim),
            times_NARROMI = rep(as.numeric(NA), n_sim),
            times_GENIE3 = rep(as.numeric(NA), n_sim),
            times_TIGRESS = rep(as.numeric(NA), n_sim))
          
          for(i in 1:n_sim){
            dat <- readRDS(file = file.path(data_path,paste0("sim",i,".RDS")))
            
            
            algo_times$times_mi_splines <- system.time(
              assign(paste0("mi_splines",i),value = algoritmos$mi_splines(dat))
            )
            algo_times$times_MRNET_splines <- system.time(
              assign(paste0("MRNET_splines",i),value = algoritmos$MRNET_splines(get(paste0("mi_splines",i))))
            )
            algo_times$times_CLR_splines <- system.time(
              assign(paste0("CLR_splines",i),value = algoritmos$CLR_splines(get(paste0("mi_splines",i))))
            )
            algo_times$times_ARACNE_splines <- system.time(
              assign(paste0("ARACNE_splines",i),value = algoritmos$ARACNE_splines(get(paste0("mi_splines",i))))
            )
            algo_times$times_mi_mm <- system.time(
              assign(paste0("mi_mm",i),value = algoritmos$mi_mm(dat))
            )
            algo_times$times_MRNET_mm <- system.time(
              assign(paste0("MRNET_mm",i),value = algoritmos$MRNET_mm(get(paste0("mi_mm",i))))
            )
            algo_times$times_CLR_mm <- system.time(
              assign(paste0("CLR_mm",i),value = algoritmos$CLR_mm(get(paste0("mi_mm",i))))
            )
            algo_times$times_ARACNE_mm <- system.time(
              assign(paste0("ARACNE_mm",i),value = algoritmos$ARACNE_mm(get(paste0("mi_mm",i))))
            )
            algo_times$times_NARROMI <- system.time(
              assign(paste0("NARROMI",i),value = algoritmos$NARROMI(dat))
            )
            algo_times$times_GENIE3 <- system.time(
              assign(paste0("GENIE3",i),value = algoritmos$GENIE3(dat))
            )
            algo_times$times_TIGRESS <- system.time(
              assign(paste0("TIGRESS",i),value = algoritmos$TIGRESS(dat))
            )
            
            
            
            for(alg in seq_along(algoritmos)){
              saveRDS(get(paste0(names(algoritmos)[alg],i)),
                      paste0(file.path(estimates_path,
                                       names(algoritmos)[alg],
                                       names(algoritmos)[alg]
                      ),
                      i,".RDS")
              )            
            }
            
            rm(list=paste0(c("mi_splines","MRNET_splines",
                             "CLR_splines","ARACNE_splines",
                             "mi_mm","MRNET_mm",
                             "CLR_mm","ARACNE_mm",
                             "NARROMI","GENIE3",
                             "TIGRESS"),i))
          }
          saveRDS(algo_times,file.path(estimates_path,"times.RDS"))
          
          
          
        }
      }
    }
  }
}




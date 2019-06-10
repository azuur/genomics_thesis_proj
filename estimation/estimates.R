require(here)
source("pckgs_and_useful_wrappers.R")
require(foreach)
require(parallel)
require(doParallel)
source("narromi.R")



# graph_options <- c(#"EC_10","EC_20","EC_50","EC_200",
#   "SC_10","SC_20",
#   "SC_50")#,
#   #"SC_200")
# type_options <- c("linear",
#                   "sigmoid")
# noise_options <- c("uniform",
#                    "gaussian")
# r_options <- c(0.2,0.5,
#   0.8)
# sample_size_options <- c(
#   #20),
#   50, 100, 1000)
n_sim <- 1e3

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


# netw_str <- graph_options[1]
# type <- type_options[1]
# noise <- noise_options[1]
# r <- r_options[1]
# sample_size <- sample_size_options[1]

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
          
          data_path <- file.path("/media","adrian","bodega","thesis",
                            "Robjects","simulated_data",
                            netw_str,
                            type,
                            paste0(noise,"_noise"),
                            paste0(r,"_noise_to_sig"),
                            paste0("sample_size_",sample_size)
          )
          
          estimates_path <- file.path("/media","adrian","bodega","thesis",
                            "Robjects","estimates",
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
          
          for(i in 725:n_sim){
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

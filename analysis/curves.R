require(here)
source("pckgs_and_useful_wrappers.R")
require(precrec)

graph_options <- c(#"EC_10","EC_20","EC_50","EC_200",
  "SC_10","SC_20",
  "SC_50")#,
#"SC_200")
type_options <- c("linear",
                  "sigmoid")
noise_options <- c("uniform",
                   "gaussian"
)
r_options <- c(0.2,0.5,
               0.8
)
sample_size_options <- c(20, 
  50, 
  100, 
  500
)
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



netw_str <- graph_options[1]
type <- type_options[1]
noise <- noise_options[1]
r <- r_options[1]
sample_size <- sample_size_options[1]
alg <- 1
i <- 1


for(sample_size in sample_size_options){
  
  for(netw_str in graph_options){
    str_readobj(netw_str,"Robjects/graph_objects")
    netw <- get(netw_str)
    
    A <- get.adjacency(netw, type = "both")
    
    A_vec <- as.logical(A)
    
    
    for(type in type_options){
      
      for(noise in noise_options){
        
        for(r in r_options){
          
          for(alg in seq_along(algoritmos)){
            
            list_of_curves <- list()
            
            
            for(i in 1:n_sim){
              
              estimates_path <- file.path("/media","adrian","bodega","thesis",
                                          "Robjects","estimates",
                                          netw_str,
                                          type,
                                          paste0(noise,"_noise"),
                                          paste0(r,"_noise_to_sig"),
                                          paste0("sample_size_",sample_size)
              )
              
              cur <- readRDS(file = paste0(file.path(estimates_path,
                                                     names(algoritmos)[alg],
                                                     names(algoritmos)[alg]),
                                           i,".RDS"))
              if(names(algoritmos)[alg] == "NARROMI"){
                cur <- cur$net
              }
              if(names(algoritmos)[alg] == "TIGRESS"){
                cur <- cur[[length(cur)]]
              }
              
              alg_m_vec <- abs(c(cur))
              alg_m_vec <- alg_m_vec/(max(alg_m_vec)+1e-10)
              
              list_of_curves[[i]] <- precrec::evalmod(scores = alg_m_vec, labels = as.numeric(A_vec))
              
              
            }
            
            roc_path <- file.path("/media","adrian","bodega","thesis",
                                            "Robjects","analysis",
                                            netw_str,
                                            type,
                                            paste0(noise,"_noise"),
                                            paste0(r,"_noise_to_sig"),
                                            paste0("sample_size_",sample_size)
            )
            
            dir.create(file.path(roc_path,names(algoritmos)[alg],"curves"), showWarnings = T, recursive = T)
            saveRDS(list_of_curves,file.path(roc_path,names(algoritmos)[alg],"curves","list_of_curves.RDS"))
            
            
            
            
            
            probabilities_path <- file.path("/media","adrian","bodega","thesis",
                                            "Robjects","analysis",
                                            netw_str,
                                            type,
                                            paste0(noise,"_noise"),
                                            paste0(r,"_noise_to_sig"),
                                            paste0("sample_size_",sample_size)
            )
            
            alg_avg <- readRDS(file.path(probabilities_path,names(algoritmos)[alg],"probabilities.RDS"))
            alg_avg_vec <- abs(c(alg_avg))
            alg_avg_vec <- alg_avg_vec/(max(alg_avg_vec)+1e-10)
            
            curve_avg <- precrec::evalmod(scores = alg_avg_vec, labels = as.numeric(A_vec))
            
            saveRDS(curve_avg,file.path(roc_path,names(algoritmos)[alg],"curves","curve_avg.RDS"))
            
            
          }
        }
      }
    }
  }
}
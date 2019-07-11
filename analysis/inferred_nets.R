require(here)
source("pckgs_and_useful_wrappers.R")
require(precrec)


graph_options <- c(#"EC_10","EC_20","EC_50","EC_200",
  #"SC_10","SC_20",
  "SC_50")#,
#"SC_200")
type_options <- c("linear"#,
                  #"sigmoid"
)
noise_options <- c("uniform",
  "gaussian"
  )

r_options <- c(0.2,
               #0.5,
               0.8)
sample_size_options <- c(20,
                         50,
                         100,
                         500
)

n_sim <- 1e3

algoritmos <- list(
  # mi_splines = function(x){
  #   ord <- min(as.integer(nrow(x)^(1/3)),4)
  #   fastGeneMI::get.mim.bspline(x, 
  #                               order = ord, 
  #                               n.cores = 5)
  # },
  mi_mm = function(x)fastGeneMI::get.mim.MM(x,discretisation = "equalwidth", n.cores = 5),
  # MRNET_splines = minet::mrnet,
  # CLR_splines = parmigene::clr,
  # ARACNE_splines = parmigene::aracne.a,
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



# netw_str <- graph_options[3]
# type <- type_options[2]
# noise <- noise_options[1]
# r <- r_options[2]
# sample_size <- sample_size_options[1]
# alg <- 11

for(sample_size in sample_size_options){
  
  for(netw_str in graph_options){
    str_readobj(netw_str,"Robjects/graph_objects")
    netw <- get(netw_str)
    
    A <- get.adjacency(netw, type = "both")
    
    A_vec <- as.logical(A)
    
    GS_density <- sum(A_vec)/length(A_vec)
    
    
    for(type in type_options){
      
      for(noise in noise_options){
        
        for(r in r_options){
          
          data_path <- file.path("/media","adrian","bodega","thesis",
                                 "Robjects","simulated_data",
                                 netw_str,
                                 type,
                                 paste0(noise,"_noise"),
                                 paste0(r,"_noise_to_sig"),
                                 paste0("sample_size_",sample_size)
          )
          
          dat <- readRDS(file = file.path(data_path,paste0("sim1.RDS")))
          
          ord <- order(as.numeric(gsub("X","",colnames(dat))))
          
          for(alg in seq_along(algoritmos)){
            #alg <- 7

          
            q90_nets <- list()
            GS_density_nets <- list()
            
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
                cur <- 1/
                  (1e-10 + 
                     (ifelse(is.na(cur$sig), max(cur$sig,na.rm = T), cur$sig) %>%
                        `diag<-`(Inf))
                   )
              }
              if(names(algoritmos)[alg] == "TIGRESS"){
                cur <- cur[[length(cur)]]
              }
              
              if(names(algoritmos)[alg] == "GENIE3"){
                ord <- order(as.numeric(gsub("X","",colnames(cur))))
              }
              
              cur <- cur[ord,ord]
              
              #cur <- cur + t(cur)
              
              diag(cur) <- 0
              
              # thr90 <- quantile(abs(cur), 0.9)
              # thrdens <- quantile(abs(cur), 1-GS_density)

              thr90 <- abs(cur)[
                order(abs(cur),decreasing = T,na.last = T)[round(0.1*length(A_vec))]
                ] 
              thrdens <- abs(cur)[
                order(abs(cur),decreasing = T,na.last = T)[sum(A_vec)+1]
                ]              
                 
              thr90 <- min(max(thr90,1e-10),1-1e-10)
              thrdens <- min(max(thrdens,1e-10),1-1e-10)   
              
              q90_nets[[i]] <- 2*(abs(cur) >= thr90) + t(abs(cur) >= thr90)
              GS_density_nets[[i]] <- 2*(abs(cur) >= thrdens) + t(abs(cur) >= thrdens)
              
            }
            
            nets_path <- file.path("/media","adrian","bodega","thesis",
                                  "Robjects","analysis","thresholded_nets",
                                  netw_str,
                                  type,
                                  paste0(noise,"_noise"),
                                  paste0(r,"_noise_to_sig"),
                                  paste0("sample_size_",sample_size)
            )
            
            dir.create(file.path(nets_path,names(algoritmos)[alg]), showWarnings = T, recursive = T)
            saveRDS(q90_nets,file.path(nets_path,names(algoritmos)[alg],"q90_nets.RDS"))
            saveRDS(GS_density_nets,file.path(nets_path,names(algoritmos)[alg],"GS_density_nets.RDS"))
            
            
          }
        }
      }
    }
  }
}
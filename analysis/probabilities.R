require(here)
source("pckgs_and_useful_wrappers.R")

graph_options <- c(#"EC_10","EC_20","EC_50","EC_200",
  #"SC_10","SC_20",
  "SC_50")#,
#"SC_200")
type_options <- c(#"linear",
  "sigmoid")
noise_options <- c("uniform"#,
                   #"gaussian"
                   )
r_options <- c(0.2#,0.5,
               #0.8
               )
sample_size_options <- c(20#, 
                         #50, 100, 500
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


# netw_str <- graph_options[1]
# type <- type_options[1]
# noise <- noise_options[1]
# r <- r_options[1]
# sample_size <- sample_size_options[1]

sample_size <- sample_size_options[1]
netw_str <- graph_options[1]
type <- type_options[1]
noise <- noise_options[1]
r <- r_options[1]
alg <- seq_along(algoritmos)[6]

i <- 1



estimates_path <- file.path("/media","adrian","bodega","thesis",
                            "Robjects","estimates",
                            netw_str,
                            type,
                            paste0(noise,"_noise"),
                            paste0(r,"_noise_to_sig"),
                            paste0("sample_size_",sample_size)
)
readRDS(file = paste0(file.path(estimates_path,
                                names(algoritmos)[alg],
                                names(algoritmos)[alg]),
                      i,".RDS")
)
# 
# for(sample_size in sample_size_options){
#   
#   for(netw_str in graph_options){
#     
#     for(type in type_options){
#       
#       for(noise in noise_options){
#         
#         for(r in r_options){
#           
#           for(alg in seq_along(algoritmos)){
#             
#                         
#           }
#         }
#       }
#     }
#   }
# }
# 
#           
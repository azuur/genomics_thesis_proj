#source("simulate_data/simulate_data_SC.R")


# first_iteration <- 725
# graph_options <- c(#"EC_10","EC_20","EC_50","EC_200",
#   #"SC_10","SC_20",
#   "SC_50")#,
# #"SC_200")
# type_options <- c("linear")
# noise_options <- c("uniform")
# r_options <- c(0.8)
# sample_size_options <- c(500)
# #n_sim <- 1e3
# source("estimation/estimates.R")




rm(list=ls())

# 
# graph_options <- c(#"EC_10","EC_20","EC_50","EC_200",
#   #"SC_10","SC_20",
#   "SC_50")#,
# #"SC_200")
# type_options <- c("linear",
#                   "sigmoid")
# noise_options <- c("uniform",
#                    "gaussian")
# r_options <- c(0.2,0.5,
#                0.8)
# sample_size_options <- c(100)
# #n_sim <- 1e3
# source("estimation/estimates.R")



graph_options <- c(#"EC_10","EC_20","EC_50","EC_200",
  "SC_10","SC_20",
  "SC_50")#,
#"SC_200")
type_options <- c("linear",
                  "sigmoid")
noise_options <- c("uniform",
                   "gaussian")

r_options <- c(0.2,
               0.5,
               0.8)
sample_size_options <- c(20, 
                         50, 
                         100, 
                         500)

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
source("analysis/curves.R")














graph_options <- c(#"EC_10","EC_20","EC_50","EC_200",
  #"SC_10","SC_20",
  #"SC_50")#,
  "SC_200")
type_options <- c("linear",
                  "sigmoid")
noise_options <- c("uniform",
                   "gaussian")

r_options <- c(0.2,
               0.5,
               0.8)
sample_size_options <- c(20)#, 
                         #50, 
                         #100, 
                         #500)

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
source("analysis/curves.R")

rm(list=ls())






















graph_options <- c(#"EC_10","EC_20","EC_50","EC_200",
  "SC_10","SC_20",
  "SC_50")#,
#"SC_200")
type_options <- c("linear",
                  "sigmoid")
noise_options <- c("uniform",
                   "gaussian")

r_options <- c(0.2,
               0.5,
               0.8)
sample_size_options <- c(20, 
                         50, 
                         100, 
                         500)

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
source("analysis/inferred_nets.R")














graph_options <- c(#"EC_10","EC_20","EC_50","EC_200",
  #"SC_10","SC_20",
  #"SC_50")#,
  "SC_200")
type_options <- c("linear",
                  "sigmoid")
noise_options <- c("uniform",
                   "gaussian")

r_options <- c(0.2,
               0.5,
               0.8)
sample_size_options <- c(20)#, 
#50, 
#100, 
#500)

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
source("analysis/inferred_nets.R")
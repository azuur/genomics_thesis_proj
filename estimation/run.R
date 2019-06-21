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

first_iteration<-1
graph_options <- c(#"EC_10","EC_20","EC_50","EC_200",
  #"SC_10","SC_20",
  "SC_50")#,
#"SC_200")
type_options <- c("linear",
                  "sigmoid")
noise_options <- c("uniform",
                   "gaussian")
r_options <- c(0.2,0.5,
               0.8)
sample_size_options <- c(100)
#n_sim <- 1e3
source("estimation/estimates.R")



# X=dat; lambda = 1; alpha = 0.05; beta = 0.05; t = 0.6; cl =NULL

# 
# j<-9
#  res <- narromi(X[,j,drop=F], X[,-j,drop=F], lambda, alpha, beta, t)
# # net <- split(res$net, 1:length(res$net)<j)
# # net <- c(net[["TRUE"]],0,net[["FALSE"]])
# # net_value <- split(res$net_value, 1:length(res$net_value)<j)
# # net_value <- c(net_value[["TRUE"]],0,net_value[["FALSE"]])
# # sig <- split(res$sig, 1:length(res$sig)<j)
# # sig <- c(sig[["TRUE"]],0,sig[["FALSE"]])
# 
# 
# y = X[,j]; X = X[,-j,drop=F]
# net_value <- 
#   G1 <- 
#   apply(X = X, MARGIN = 2, FUN = function(x) narromi_cmi(x,y))
# #    apply(X = X, MARGIN = 2, FUN = function(x) infotheo::mutinformation(x,y))
# 
# idx <- which(abs(G1) >= alpha)
# X<-X[,idx, drop=F]
# 
# ini_B <- lm(y~ -1+X)$coefficients
# ini_B[is.na(ini_B)] <- 0.004
# ini_B <- rep(1, ncol(X))
# MTE::LADlasso(y,X,
#          ini_B,
#          lambda = lambda/nrow(X))
# LP_TGN(y,X,lambda)
# 
# 
# 
# 
# if(T){
#   if(T){
#     if(F){
#       print("pipi")
#     } else{ next }
#   }
#   print("popo")
# }

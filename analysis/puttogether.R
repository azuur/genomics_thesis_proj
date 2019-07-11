require(here)
source("pckgs_and_useful_wrappers.R")
require(ggthemes)
require(precrec)


graph_options <- c(#"EC_10","EC_20","EC_50","EC_200",
  #"SC_10","SC_20",
  "SC_50")#,
#"SC_200")
type_options <- c("linear"#,
                  #"sigmoid"
)
noise_options <- c("uniform",
  "gaussian")

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



# netw_str <- graph_options[1]
# type <- type_options[1]
# noise <- noise_options[1]
# r <- r_options[2]
# sample_size <- sample_size_options[1]
# alg <- 5





eval_measures_cutoff <- list()
eval_measures_global <- list()

for(sample_size in sample_size_options){
  
  for(netw_str in graph_options){
    str_readobj(netw_str,"Robjects/graph_objects")
    netw <- get(netw_str)
    
    A <- get.adjacency(netw, type = "both")
    
    A_vec <- as.logical(A)
    
    for(type in type_options){
      
      for(noise in noise_options){
        
        for(r in r_options){
          #alg <- 7
          for(alg in seq_along(algoritmos)){
            
            
            nets_path <- file.path("/media","adrian","bodega","thesis",
                                   "Robjects","analysis","thresholded_nets",
                                   netw_str,
                                   type,
                                   paste0(noise,"_noise"),
                                   paste0(r,"_noise_to_sig"),
                                   paste0("sample_size_",sample_size)
            )
            list_of_nets <- 
              readRDS(file.path(nets_path,names(algoritmos)[alg],"q90_nets.RDS"))
            
            list_of_nets <- lapply(list_of_nets, function(x) x>0)
            mat <- list_of_nets[[1]]
            #,ncol= sqrt(length(list_of_nets[[1]])))
            
            if(isSymmetric(mat)){
              A <- (as.matrix(A) + t(as.matrix(A)))/2
              A_vec <- as.logical(A)
            }
            
            #PPV/precision
            PPV <- sapply(list_of_nets,
                          function(x){
                            sum(x==T & A_vec==1)/sum(x)
                          })
            mean_PPV <- mean(PPV)
            sd_PPV <- sd(PPV)
            
            
            #NPV
            NPV <- sapply(list_of_nets,
                          function(x){
                            sum(x==F & A_vec==0)/sum(1-x)
                          })
            mean_NPV <- mean(NPV)
            sd_NPV <- sd(NPV)
            
            
            
            #TPR/sensitivity/recall
            sensitivity <- sapply(list_of_nets,
                                  function(x){
                                    sum(x==T & A_vec==1)/sum(A_vec)
                                  })
            mean_sensitivity <- mean(sensitivity)
            sd_sensitivity <- sd(sensitivity)
            
            
            
            
            #TNR,specificity, 1-FPR
            specificity <- sapply(list_of_nets,
                                  function(x){
                                    sum(x==F & A_vec==0)/sum(1-A_vec)
                                  })
            mean_specificity <- mean(specificity)
            sd_specificity <- sd(specificity)
            
            
            #F1Score
            F1 <- 2*PPV*sensitivity/(PPV+sensitivity)
            mean_F1 <- mean(F1)
            sd_F1 <- sd(F1)
            
            eval_measures_cutoff[[length(eval_measures_cutoff)+1]] <-
              data.frame(netw = netw_str,
                         sample_size= sample_size,
                         r=r,
                         noise = noise,
                         type=type,
                         algorithm=names(algoritmos)[alg],
                         estimator_cutoff="q90",
                         mean_PPV = mean_PPV,
                         sd_PPV = sd_PPV,
                         mean_NPV = mean_NPV,
                         sd_NPV = sd_NPV,
                         mean_sensitivity = mean_sensitivity,
                         sd_sensitivity = sd_sensitivity,
                         mean_specificity = mean_specificity,
                         sd_specificity = sd_specificity,
                         mean_F1=mean_F1,
                         sd_F1= sd_F1
              )
            
            
            list_of_nets <-
              readRDS(file.path(nets_path,names(algoritmos)[alg],"GS_density_nets.RDS"))
            
            list_of_nets <- lapply(list_of_nets, function(x) x>0)
            mat <- list_of_nets[[1]]
            
            #PPV/precision
            PPV <- sapply(list_of_nets,
                          function(x){
                            sum(x==T & A_vec==1)/sum(x)
                          })
            mean_PPV <- mean(PPV)
            sd_PPV <- sd(PPV)
            
            
            #NPV
            NPV <- sapply(list_of_nets,
                          function(x){
                            sum(x==F & A_vec==0)/sum(1-x)
                          })
            mean_NPV <- mean(NPV)
            sd_NPV <- sd(NPV)
            
            
            
            #TPR/sensitivity/recall
            sensitivity <- sapply(list_of_nets,
                                  function(x){
                                    sum(x==T & A_vec==1)/sum(A_vec)
                                  })
            mean_sensitivity <- mean(sensitivity)
            sd_sensitivity <- sd(sensitivity)
            
            
            
            
            #TNR,specificity, 1-FPR
            specificity <- sapply(list_of_nets,
                                  function(x){
                                    sum(x==F & A_vec==0)/sum(1-A_vec)
                                  })
            mean_specificity <- mean(specificity)
            sd_specificity <- sd(specificity)
            
            
            #F1Score
            F1 <- 2*PPV*sensitivity/(PPV+sensitivity)
            mean_F1 <- mean(F1)
            sd_F1 <- sd(F1)
            
            eval_measures_cutoff[[length(eval_measures_cutoff)+1]] <-
              data.frame(netw = netw_str,
                         sample_size= sample_size,
                         r=r,
                         noise = noise,
                         type=type,
                         algorithm=names(algoritmos)[alg],
                         estimator_cutoff="GS_density",
                         mean_PPV = mean_PPV,
                         sd_PPV = sd_PPV,
                         mean_NPV = mean_NPV,
                         sd_NPV = sd_NPV,
                         mean_sensitivity = mean_sensitivity,
                         sd_sensitivity = sd_sensitivity,
                         mean_specificity = mean_specificity,
                         sd_specificity = sd_specificity,
                         mean_F1=mean_F1,
                         sd_F1= sd_F1
              )
            
            roc_path <- file.path("/media","adrian","bodega","thesis",
                                  "Robjects","analysis","ROC_curves",
                                  netw_str,
                                  type,
                                  paste0(noise,"_noise"),
                                  paste0(r,"_noise_to_sig"),
                                  paste0("sample_size_",sample_size)
            )

            list_of_curves <- readRDS(file.path(roc_path,names(algoritmos)[alg],"list_of_curves.RDS"))

            #AUROC
            AUROC <- sapply(list_of_curves,
                            function(x){
                              precrec::auc(x) %>%
                                filter(curvetypes == "ROC") %>%
                                select(aucs) %>% unlist()
                            })
            mean_AUROC <- mean(AUROC)
            p90_AUROC <- quantile(AUROC,0.9)
            p10_AUROC <- quantile(AUROC,0.1)
            sd_AUROC <- sd(AUROC)

            #AUPRC
            AUPRC <- sapply(list_of_curves,
                            function(x){
                              precrec::auc(x) %>%
                                filter(curvetypes == "PRC") %>%
                                select(aucs) %>% unlist()
                            })
            mean_AUPRC <- mean(AUPRC)
            sd_AUPRC <- sd(AUPRC)

            eval_measures_global[[length(eval_measures_global)+1]] <-
              data.frame(netw = netw_str,
                         sample_size= sample_size,
                         r=r,
                         noise = noise,
                         type=type,
                         algorithm=names(algoritmos)[alg],
                         mean_AUROC = mean_AUROC,
                         p90_AUROC = p90_AUROC,
                         p10_AUROC = p10_AUROC,
                         sd_AUROC = sd_AUROC,
                         mean_AUPRC = mean_AUPRC,
                         sd_AUPRC = sd_AUPRC
              )
            print(length(eval_measures_cutoff))
          }
        }
      }
    }
  }
}

df_eval_measures_cutoff <- do.call(rbind, eval_measures_cutoff)
df_eval_measures_global <- do.call(rbind, eval_measures_global)





overall_results_path <- file.path("/media","adrian","bodega","thesis",
                                  "Robjects","analysis","summaries")
dir.create(overall_results_path,showWarnings = F,recursive = T)
write.table(df_eval_measures_cutoff, 
            file.path(overall_results_path,"df_eval_measures_cutoff.csv"),
            row.names = F, sep = ";", dec = ".")
write.table(df_eval_measures_global,
            file.path(overall_results_path,"df_eval_measures_global.csv"),
            row.names = F, sep = ";", dec = ".")
saveRDS(df_eval_measures_cutoff, 
        file.path(overall_results_path,"df_eval_measures_cutoff.RDS"))
saveRDS(df_eval_measures_global,
        file.path(overall_results_path,"df_eval_measures_global.RDS"))



AUROCs <- list()

for(sample_size in sample_size_options){
  
  for(netw_str in graph_options){
    str_readobj(netw_str,"Robjects/graph_objects")
    netw <- get(netw_str)
    
    A <- get.adjacency(netw, type = "both")
    
    A_vec <- as.logical(A)
    
    for(type in type_options){
      
      for(noise in noise_options){
        
        for(r in r_options){
          #alg <- 7
          for(alg in seq_along(algoritmos)){
            roc_path <- file.path("/media","adrian","bodega","thesis",
                                  "Robjects","analysis","ROC_curves",
                                  netw_str,
                                  type,
                                  paste0(noise,"_noise"),
                                  paste0(r,"_noise_to_sig"),
                                  paste0("sample_size_",sample_size)
            )
            
            list_of_curves <- readRDS(file.path(roc_path,names(algoritmos)[alg],"list_of_curves.RDS"))
            
            #AUROC
            AUROC <- sapply(list_of_curves,
                            function(x){
                              precrec::auc(x) %>%
                                filter(curvetypes == "ROC") %>%
                                select(aucs) %>% unlist()
                            })
            AUROCs[[length(AUROCs)+1]] <-
              data.frame(netw = netw_str,
                         sample_size= sample_size,
                         r=r,
                         noise = noise,
                         type=type,
                         algorithm=names(algoritmos)[alg],
                         AUROC = AUROC
              )
            print(length(AUROCs))
          }
        }
      }
    }
  }
} 
df_AUROCs <- do.call(rbind, AUROCs)





overall_results_path <- file.path("/media","adrian","bodega","thesis",
                                  "Robjects","analysis","summaries")
dir.create(overall_results_path,showWarnings = F,recursive = T)
write.table(df_AUROCs, 
            file.path(overall_results_path,"df_AUROCs.csv"),
            row.names = F, sep = ";", dec = ".")
saveRDS(df_AUROCs, 
        file.path(overall_results_path,"df_AUROCs.RDS"))


# 
# df_eval_measures_cutoff <-readRDS(file.path(overall_results_path,"df_eval_measures_cutoff.RDS"))
# 
# 
# 
# df_eval_measures_cutoff %>% 
#   filter(type == "linear", 
#          estimator_cutoff=="GS_density",
#          netw == "SC_50") %>% 
#   select(-type, -estimator_cutoff,-netw) %>% 
#   View()





































# 
# 
# 
# netw_str <- graph_options[3]
# type <- type_options[1]
# noise <- noise_options[2]
# r <- r_options[2]
# sample_size <- sample_size_options[1]
# alg <- 11
# 
# 
# 
# 
# 
# for(sample_size in sample_size_options){
#   
#   for(netw_str in graph_options){
#     str_readobj(netw_str,"Robjects/graph_objects")
#     netw <- get(netw_str)
#     
#     A <- get.adjacency(netw, type = "both")
#     
#     A_vec <- as.logical(A)
#     
#     for(type in type_options){
#       
#       for(noise in noise_options){
#         
#         for(r in r_options){
#           
#           for(alg in seq_along(algoritmos)){
#             
#             roc_path <- file.path("/media","adrian","bodega","thesis",
#                                   "Robjects","analysis",
#                                   netw_str,
#                                   type,
#                                   paste0(noise,"_noise"),
#                                   paste0(r,"_noise_to_sig"),
#                                   paste0("sample_size_",sample_size)
#             )
#             
#             list_of_curves <- readRDS(file.path(roc_path,names(algoritmos)[alg],"curves","list_of_curves.RDS"))
#             
#             
#             curve_avg <- readRDS(file.path(roc_path,names(algoritmos)[alg],"curves","curve_avg.RDS"))
#             
#             df_roc_curves <- 
#               lapply(1:length(list_of_curves),
#                      function(i){
#                        data.frame(
#                          X=list_of_curves[[i]]$rocs[[1]]$x,
#                          Y=list_of_curves[[i]]$rocs[[1]]$y,
#                          I=as.character(i))
#                      }) %>%
#               do.call(what = rbind.data.frame, args = .) %>%
#               rbind(data.frame(
#                 X=curve_avg$rocs[[1]]$x,
#                 Y=curve_avg$rocs[[1]]$y,
#                 I="avg"))
#             
#             df_roc_curves %>% 
#               ggplot() +
#               geom_line(aes(x = X, y = Y, group = I), 
#                         colour = "orange", 
#                         alpha = 0.1) +
#               # geom_line(data = df_roc_curves %>% filter(I == "avg"),
#               #           aes(x = X, y = Y), colour = "black") +
#               ggthemes::theme_economist() #+
#             # geom_density(data = data_frame(outcome), aes(x = outcome), colour = "black", adjust = 0.8) +
#             # ggtitle("Actual outcomes and posterior predictive replications") +
#             # annotate("text", x = 0.2, y = 5, label = "Density of actual outcomes", hjust = 0) +
#             # annotate("text", x = 0.2, y = 3.5, label = "Posterior replications", colour = "orange", hjust = 0) 
#             
#           }
#         }
#       }
#     }
#   }
# }
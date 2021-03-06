# overall_results_path <- file.path("/media","adrian","bodega","thesis",
#                                   "Robjects","analysis","summaries")
# 
# df_eval_measures_cutoff <- readRDS(
#   file.path(overall_results_path,"df_eval_measures_cutoff.RDS")
#   )
# df_eval_measures_global <- readRDS( 
#         file.path(overall_results_path,"df_eval_measures_global.RDS")
#         )
# 
# df_eval_measures_global %>%
#   filter(noise == "gaussian", type == "linear", netw == "SC_50") %>%
#   arrange(desc(mean_AUROC)) %>% View()



require(here)
source("pckgs_and_useful_wrappers.R")
require(ggthemes)
require(RColorBrewer)
require(gganimate)

graph_options <- c(#"EC_10","EC_20","EC_50","EC_200",
  #"SC_10","SC_20",
  "SC_50")#,
#"SC_200")
type_options <- c("linear"#,
                  #"sigmoid"
)
noise_options <- c(#"uniform",
  "gaussian",
  "laplacian")

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
# alg <- 10





# str_readobj(netw_str,"Robjects/graph_objects")
# netw <- get(netw_str)
# 
# A <- get.adjacency(netw, type = "both")
# 
# A_vec <- as.logical(A)



for(alg in seq_along(algoritmos)){
  
  plot_path <- file.path("/media","adrian","bodega","thesis",
                         "Robjects","analysis",
                         "plots","ss_plots","samples",names(algoritmos)[alg])
  
  dir.create(plot_path, showWarnings = T, recursive = T)
  
  for(netw_str in graph_options){

    for(type in type_options){
      
      for(noise in noise_options){
        
        for(r in r_options){
          df_roc_curves <- data.frame(X=numeric(0),
                                      Y=numeric(0),
                                      I=character(0),
                                      sample_size=numeric(0))
          
          for(sample_size in sample_size_options){
            roc_path <- file.path("/media","adrian","bodega","thesis",
                                  "Robjects","analysis","ROC_curves",
                                  netw_str,
                                  type,
                                  paste0(noise,"_noise"),
                                  paste0(r,"_noise_to_sig"),
                                  paste0("sample_size_",sample_size)
            )
            
            list_of_curves <- readRDS(file.path(roc_path,names(algoritmos)[alg],"list_of_curves.RDS"))
            
            df_roc_curves <- 
              df_roc_curves %>% union_all( 
                lapply(1:length(list_of_curves),
                       function(i){
                         data.frame(
                           X=list_of_curves[[i]]$rocs[[1]]$x,
                           Y=list_of_curves[[i]]$rocs[[1]]$y,
                           I=as.character(i))
                       }) %>%
                  do.call(what = rbind.data.frame, args = .) %>% 
                  mutate(sample_size=sample_size))
          }
          
          
          capt <- paste0("netw_str = ",netw_str,"; ",
                         "type = ",type,"; ",
                         "noise = ",noise,"; ",
                         "r = ",r)
          set.seed(1308)
          g1 <- (df_roc_curves %>% 
                   filter(I %in% sample(size = 100, 1:1000)) %>%
                   mutate(I = as.factor(paste(I,sample_size))) %>%
                   mutate(sample_size = as.factor(sample_size))) %>%
            ggplot() +
            geom_line(aes(x = X, y = Y, group = I,color = sample_size), 
                      alpha = 0.1) +
            theme_minimal() + 
            #scale_color_fivethirtyeight() +
            #scale_color_brewer(palette = "Reds") +
            scale_color_stata() +
            xlab("1-specificity (1-TN/N)") + ylab("sensitivity (TP/P)") +
            labs(title = paste0("ROCs for ",names(algoritmos)[alg]),
                 caption = capt) + 
            guides(colour = guide_legend(override.aes = list(alpha = 1)))
          g1
          
          ggsave(file.path(plot_path,
                           paste(c(netw_str,
                                   type,
                                   paste0(noise,"_noise"),
                                   paste0(r,"_noise_to_sig")),collapse ="_")
          ), 
          g1,
          device = "png")
          
          rm(df_roc_curves)
          
        }
      }
    }
  }
}
  






for(alg in seq_along(algoritmos)){
  
  plot_path <- file.path("/media","adrian","bodega","thesis",
                         "Robjects","analysis",
                         "plots","ss_plots","avg",names(algoritmos)[alg])
  
  dir.create(plot_path, showWarnings = T, recursive = T)
  
  for(netw_str in graph_options){
    
    for(type in type_options){
      
      for(noise in noise_options){
        
        for(r in r_options){
          df_roc_curves <- data.frame(X=numeric(0),
                                      Y=numeric(0),
                                      sample_size=numeric(0))
          
          for(sample_size in sample_size_options){
            roc_path <- file.path("/media","adrian","bodega","thesis",
                                  "Robjects","analysis","ROC_curves",
                                  netw_str,
                                  type,
                                  paste0(noise,"_noise"),
                                  paste0(r,"_noise_to_sig"),
                                  paste0("sample_size_",sample_size)
            )
            
            curve_avg <- readRDS(file.path(roc_path,names(algoritmos)[alg],"curve_avg.RDS"))          
            
            df_roc_curves <- 
              df_roc_curves %>% union_all( 
                data.frame(
                  X=curve_avg$rocs[[1]]$x,
                  Y=curve_avg$rocs[[1]]$y
                ) %>% 
                  mutate(sample_size=sample_size))
          }
          
          
          capt <- paste0("netw_str = ",netw_str,"; ",
                         "type = ",type,"; ",
                         "noise = ",noise,"; ",
                         "r = ",r)
          set.seed(1308)
          g1 <- (df_roc_curves %>% 
                   #filter(I %in% sample(size = 100, 1:1000)) %>%
                  # mutate(I = as.factor(paste(I,sample_size))) %>%
                   mutate(sample_size = as.factor(sample_size))) %>%
            ggplot() +
            geom_line(aes(x = X, y = Y,color = sample_size), 
                      alpha = 1) +
            theme_minimal() + 
            #scale_color_fivethirtyeight() +
            #scale_color_brewer(palette = "Reds") +
            scale_color_stata() +
            xlab("1-specificity (1-TN/N)") + ylab("sensitivity (TP/P)") +
            labs(title = paste0("ROCs for ",names(algoritmos)[alg]),
                 caption = capt) + 
            guides(colour = guide_legend(override.aes = list(alpha = 1)))
          g1
          
          ggsave(file.path(plot_path,
                           paste(c(netw_str,
                                   type,
                                   paste0(noise,"_noise"),
                                   paste0(r,"_noise_to_sig")),collapse ="_")
          ), 
          g1,
          device = "png")
          
          rm(df_roc_curves)
          
        }
      }
    }
  }
}











overall_results_path <- file.path("/media","adrian","bodega","thesis",
                                  "Robjects","analysis","summaries")

# C:\Users\ASUS PC\Documents\summaries

# overall_results_path <- file.path("C:","Users","ASUS PC","Documents","summaries")

df_eval_measures_global <- readRDS(file.path(overall_results_path,"df_eval_measures_global.RDS"))

df_list_plot <- list()
for(k in seq_along(sample_size_options)){
  df_list_plot[[k]] <- 
    df_eval_measures_global %>% 
    filter(noise == "laplacian", sample_size == sample_size_options[k]) %>%
    filter(r==0.2) %>%
    select(algorithm,sample_size,mean_AUROC,p10_AUROC,p90_AUROC) %>%
    rename(mean_AUROC02 = mean_AUROC,
           p10_AUROC02 = p10_AUROC,
           p90_AUROC02 = p90_AUROC) %>%
    left_join(
      df_eval_measures_global %>%
        filter(noise == "laplacian", sample_size == sample_size_options[k]) %>%
        filter(r==0.8) %>%
        select(algorithm,sample_size,mean_AUROC,p10_AUROC,p90_AUROC) %>%
        rename(mean_AUROC08 = mean_AUROC,
               p10_AUROC08 = p10_AUROC,
               p90_AUROC08 = p90_AUROC)
    )
}

df_plot <- do.call(rbind,df_list_plot)

dfp<-df_plot %>% 
  filter(sample_size < 500) %>%
  left_join(
    df_plot %>% 
      filter(sample_size < 500) %>%
      dplyr::group_by(algorithm) %>% 
    dplyr::summarize(pos = max(pmax(mean_AUROC02,mean_AUROC08)))
  ) %>%
  arrange(pos) %>%
  mutate(algorithm = factor(algorithm,levels = unique(algorithm)))

plot_AUROC <-
  ggplot(dfp) +
  geom_segment( aes(x=algorithm, xend=algorithm, y=mean_AUROC02, yend=mean_AUROC08), color="grey") +
  geom_rect( aes(xmin=algorithm, xmax=algorithm, ymin=p10_AUROC02, ymax=p90_AUROC02), color="green", position = "dodge") +
  geom_rect( aes(xmin=algorithm, xmax=algorithm, ymin=p10_AUROC08, ymax=p90_AUROC08), color="red", position = "dodge") +
  # geom_segment( aes(x=algorithm, xend=algorithm, y=p10_AUROC02, yend=p90_AUROC02), color="green") +
  # geom_segment( aes(x=algorithm, xend=algorithm, y=p10_AUROC08, yend=p90_AUROC08), color="red") +
  geom_point( aes(x=algorithm, y=mean_AUROC02), color=rgb(0.2,0.7,0.1,0.5), size=2 ) +
  geom_point( aes(x=algorithm, y=mean_AUROC08), color=rgb(0.7,0.2,0.1,0.5), size=2 ) +
  coord_flip()







overall_results_path <- file.path("/media","adrian","bodega","thesis",
                                  "Robjects","analysis","summaries")

# C:\Users\ASUS PC\Documents\summaries

# overall_results_path <- file.path("C:","Users","ASUS PC","Documents","summaries")


df_eval_measures_global <- readRDS(file.path(overall_results_path,"df_eval_measures_global.RDS"))


df_list_plot <- list()
for(k in seq_along(sample_size_options)){
  df_list_plot[[k]] <- 
    df_eval_measures_global %>% 
    filter(noise == "laplacian", sample_size == sample_size_options[k]) %>%
    select(algorithm,sample_size,r,mean_AUROC,p10_AUROC,p90_AUROC) %>%
    mutate(r = as.factor(r), sample_size = as.factor(sample_size))
}

df_plot <- do.call(rbind,df_list_plot)

levels(df_plot$algorithm) <- 
  toupper(gsub("_mm", "", levels(df_plot$algorithm)))

dfp<-df_plot %>%
  # mutate(algorithm = ifelse(algorithm == "ARACNE_mm","ARACNE",algorithm)) %>%
  # mutate(algorithm = ifelse(algorithm == "MRNET_mm","MRNET",algorithm)) %>%
  # mutate(algorithm = ifelse(algorithm == "mI_mm","MutInf",algorithm)) %>%
  # mutate(algorithm = ifelse(algorithm == "CLR_mm","CLR",algorithm)) %>% 
  #  filter(sample_size != 500) %>%
  left_join(
    df_plot %>% 
      #      filter(sample_size != 500) %>%
      dplyr::group_by(algorithm) %>% 
      dplyr::summarize(pos = max(mean_AUROC))
  ) %>%
  arrange(pos) %>%
  mutate(algorithm = factor(algorithm,levels = unique(algorithm)))

dfp <- dfp %>%
  left_join(
    df_plot %>%
      group_by(algorithm,sample_size) %>%
      summarize(y = min(mean_AUROC), yend = max(mean_AUROC)) %>%
      ungroup()
  )

dfp <- dfp %>% 
  rename(`Sample Size` = sample_size, FVU = r) %>% 
  mutate(algorithm = ifelse(algorithm=="mi_mm", "MUT INF",algorithm))%>% 
  mutate(algorithm = ifelse(algorithm=="CLR_mm", "CLR",algorithm))%>% 
  mutate(algorithm = ifelse(algorithm=="MRNET_mm", "MRNET",algorithm))%>% 
  mutate(algorithm = ifelse(algorithm=="ARACNE_mm", "ARACNE",algorithm))


mutinf <- c("ARACNE","CLR","MRNET", "MUT INF")

plot_AUROC_MI <-
  ggplot(dfp %>% 
           filter(algorithm %in% mutinf) %>% arrange(pos), 
         aes(x=algorithm, 
                  y = mean_AUROC, 
                  color = `Sample Size`, shape = FVU
                  #color = r, shape = sample_size
                  )
         ) +
  theme_bw() +
  ylim(c(0.5,1)) +
  theme(plot.title = element_text(size = rel(1.4), face = "bold"),
        plot.subtitle = element_text(size = rel(1.4)),
        legend.title=element_text(size=rel(1.2)), 
        legend.text=element_text(size=rel(1.2)),
        axis.text.y = element_text(size = rel(1.8)),
        axis.text.x = element_text(size = rel(1.5))) +
  geom_point(size = 4, position = position_dodge(width = 1)) +
  geom_errorbar(
    aes(ymin = p10_AUROC, ymax = p90_AUROC),
    alpha = 1,
    width = 0.5,
    size = 1.5, position = position_dodge(width = 1)) + 
  scale_color_discrete(name = "Sample Size") +
  scale_shape_discrete(name = "Frac. of Var. Unexplained") +
  labs(x= "", y = "AUROC", title ="AUROC for Mutual Information-based Algorithms - Laplacian noise",
       subtitle = "Sample averages & bars between quantiles 0.1 and 0.9")+
  # geom_segment(
  #   aes(x=algorithm,
  #       xend=algorithm,
  #       y=y,
  #       yend=yend), 
  #   color="black") +
  coord_flip() 

plot_AUROC_MI



plot_AUROC_REG <-
  ggplot(dfp %>% 
           filter(!algorithm %in% mutinf) %>% arrange(pos), 
         aes(x=algorithm, 
             y = mean_AUROC, 
             color = `Sample Size`, shape = FVU
             #color = r, shape = sample_size
         )
  ) +
  theme_bw() +
  ylim(c(0.5,1)) +
  theme(plot.title = element_text(size = rel(1.4), face = "bold"),
        plot.subtitle = element_text(size = rel(1.4)),
        legend.title=element_text(size=rel(1.2)), 
        legend.text=element_text(size=rel(1.2)),
        axis.text.y = element_text(size = rel(1.8)),
        axis.text.x = element_text(size = rel(1.5))) +
  geom_point(size = 4, position = position_dodge(width = 1)) +
  geom_errorbar(
    aes(ymin = p10_AUROC, ymax = p90_AUROC),
    alpha = 1,
    width = 0.5,
    size = 1.5, position = position_dodge(width = 1)) + 
  scale_color_discrete(name = "Sample Size") +
  scale_shape_discrete(name = "Frac. of Var. Unexplained") +
  labs(x= "", y = "AUROC", title ="AUROC for Regression-based Algorithms - Laplacian noise",
       subtitle = "Sample averages & bars between quantiles 0.1 and 0.9")+
  # geom_segment(
  #   aes(x=algorithm,
  #       xend=algorithm,
  #       y=y,
  #       yend=yend), 
  #   color="black") +
  coord_flip() 

plot_AUROC_REG


ggpubr::ggarrange(plot_AUROC_REG,
                  plot_AUROC_MI, align = "v", nrow= 2, ncol=1
)











dfp %>%  rename(AUROC = mean_AUROC, 
         Algorithm = algorithm,
         `Sample Size` = sample_size,
         FVU = r) 

plot_AUROC <-
  ggplot(dfp, aes(x=Algorithm, 
                  y = AUROC, 
                  color = `Sample Size`, shape = FVU
                  )
         ) +
  theme_bw() + 
  scale_color_brewer(palette="Dark2") +
  geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
  geom_errorbar(
    aes(ymin = p10_AUROC, ymax = p90_AUROC),
    alpha = 0.7,
    width = 0.3,
    size = 0.5, position = position_dodge(width = 0.5)) +
  coord_flip() + 
  labs(color = "Sample Size",
       shape = "Frac. of Var. Unexplained",
       caption = "Dots represent means. Bars range between quantiles 0.1 and 0.9. All are estimated with 1000 simulated datasets.") + 
#  theme(legend.position="bottom") + 
  ggtitle("Distributions of AUROC",
          "Multivariate normal model over 50 node subgraph") + 
  theme(
    plot.title = element_text(size = 16),
    plot.subtitle = element_text(size = 14),
    plot.caption = element_text(size = 12)
  )


animated_AUROC <-
  plot_AUROC +
  transition_states(`Sample Size`, 
                    state_length = 50,
                    wrap = TRUE) +
  shadow_mark(alpha=alpha/4) 

animated_AUROC <- 
  animated_AUROC %>% 
  animate(fps = 1, width = 1103*0.7,height = 881*0.7, units = "px")


>>>>>>> refs/remotes/origin/master

plot_path <- file.path("/media","adrian","bodega","thesis",
                       "Robjects","analysis",
                       "plots")

ggsave(file.path(plot_path,"AUROCs"), plot_AUROC, device = "png",
       width = 1103/37,height = 881/37, units = "cm")
anim_save(file.path(plot_path,"anim_AUROCs"), animated_AUROC)



3520*881/1103




















df_eval_measures_cutoff <- 
  readRDS(file.path(overall_results_path,"df_eval_measures_cutoff.RDS"))


df_list_plot <- list()

for(k in seq_along(sample_size_options)){
  df_list_plot[[k]] <- 
    df_eval_measures_cutoff %>% 
    filter(estimator_cutoff == "GS_density", 
           sample_size == sample_size_options[k]) %>%
    filter(r==0.2) %>%
    select(algorithm,sample_size,mean_PPV,sd_PPV) %>%
    rename(mean_PPV02 = mean_PPV, sd_PPV02 = sd_PPV) %>%
    left_join(
      df_eval_measures_cutoff  %>% 
        filter(estimator_cutoff == "GS_density", 
               sample_size == sample_size_options[k]) %>%
        filter(r==0.8) %>%
        select(algorithm,sample_size,mean_PPV,sd_PPV) %>%
        rename(mean_PPV08 = mean_PPV, sd_PPV08 = sd_PPV)
    )
}

df_plot <- do.call(rbind,df_list_plot)

dfp<-df_plot %>% 
  left_join(
    df_plot %>% 
      dplyr::group_by(algorithm) %>% 
      dplyr::summarize(pos = max(pmax(mean_PPV02,mean_PPV08)))
  ) %>%
  arrange(pos) %>%
  mutate(algorithm = factor(algorithm,levels = unique(algorithm)))

plot_PPV <-
  ggplot(dfp) +
  geom_segment( aes(x=algorithm, xend=algorithm, y=mean_PPV02, yend=mean_PPV08), color="grey") +
  geom_point( aes(x=algorithm, y=mean_PPV02), color=rgb(0.2,0.7,0.1,0.5), size=3 ) +
  geom_point( aes(x=algorithm, y=mean_PPV08), color=rgb(0.7,0.2,0.1,0.5), size=3 ) +
  coord_flip() 

plot_path <- file.path("/media","adrian","bodega","thesis",
                       "Robjects","analysis",
                       "plots")

ggsave(file.path(plot_path,"PPVs"), plot_PPV, device = "png")



# 
# 
# roc_path <- file.path("/media","adrian","bodega","thesis",
#                       "Robjects","analysis",
#                       netw_str,
#                       type,
#                       paste0(noise,"_noise"),
#                       paste0(r,"_noise_to_sig"),
#                       paste0("sample_size_",sample_size)
# )
# 
# list_of_curves <- readRDS(file.path(roc_path,names(algoritmos)[alg],"curves","list_of_curves.RDS"))
# #curve_avg <- readRDS(file.path(roc_path,names(algoritmos)[alg],"curves","curve_avg.RDS"))
# 
# df_roc_curves <- 
#   lapply(1:length(list_of_curves),
#          function(i){
#            data.frame(
#              X=list_of_curves[[i]]$rocs[[1]]$x,
#              Y=list_of_curves[[i]]$rocs[[1]]$y,
#              I=as.character(i))
#          }) %>%
#   do.call(what = rbind.data.frame, args = .) %>% 
#   mutate(sample_size=sample_size)
# 
# 
# sample_size <- sample_size_options[2]
# roc_path <- file.path("/media","adrian","bodega","thesis",
#                       "Robjects","analysis",
#                       netw_str,
#                       type,
#                       paste0(noise,"_noise"),
#                       paste0(r,"_noise_to_sig"),
#                       paste0("sample_size_",sample_size)
# )
# 
# list_of_curves <- readRDS(file.path(roc_path,names(algoritmos)[alg],"curves","list_of_curves.RDS"))
# #curve_avg <- readRDS(file.path(roc_path,names(algoritmos)[alg],"curves","curve_avg.RDS"))
# 
# df_roc_curves <- df_roc_curves %>% union_all( 
#   lapply(1:length(list_of_curves),
#          function(i){
#            data.frame(
#              X=list_of_curves[[i]]$rocs[[1]]$x,
#              Y=list_of_curves[[i]]$rocs[[1]]$y,
#              I=as.character(i))
#          }) %>%
#   do.call(what = rbind.data.frame, args = .) %>% 
#     mutate(sample_size=sample_size))
# 
# 
# sample_size <- sample_size_options[3]
# roc_path <- file.path("/media","adrian","bodega","thesis",
#                       "Robjects","analysis",
#                       netw_str,
#                       type,
#                       paste0(noise,"_noise"),
#                       paste0(r,"_noise_to_sig"),
#                       paste0("sample_size_",sample_size)
# )
# 
# list_of_curves <- readRDS(file.path(roc_path,names(algoritmos)[alg],"curves","list_of_curves.RDS"))
# #curve_avg <- readRDS(file.path(roc_path,names(algoritmos)[alg],"curves","curve_avg.RDS"))
# 
# 

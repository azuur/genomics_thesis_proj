require(here)
source("functions.R")



net <- fread("DREAM5_NetworkInference_GoldStandard_Network1.tsv")
net_s <- fread("net1_gold_standard_signed.tsv",header = F)

geneids <- fread("net1_gene_ids.tsv")

tfs <- fread("net1_transcription_factors.tsv")
chip_features <- fread("net1_chip_features.tsv")
datos <- fread("net1_expression_data.tsv")


summary(net)




##next step: calculate net
#transitity/clustering
#all graphs from network science book chap 9

net_s
sqrt(dim(net)[1])




#datos_SC <- readxl::read_excel(here("SCerevisiae_data/PLOSONE_Their_database.xlsx"),)
nombres_replicas_SC <- unlist(datos_SC[,1])

datos_SC <- datos_SC %>% 
  select(-1) %>%
  t() %>% 
  as.data.frame()

datos

datos_SC %>% sapply(FUN = var) %>% summary()
datos_SC %>% sapply(FUN = var) %>% hist()
datos_SC %>% sapply(FUN = var) %>% boxplot()
datos_SC %>% sapply(FUN = mean) %>% summary()
datos_SC %>% sapply(FUN = mean) %>% hist()
datos_SC %>% sapply(FUN = mean) %>% boxplot()
datos_SC %>% sapply(FUN = sd) %>% summary()
datos_SC %>% sapply(FUN = sd) %>% hist()
datos_SC %>% sapply(FUN = sd) %>% boxplot()
datos_SC %>% sapply(FUN = function(x)sd(x)/mean(x)) %>% summary()
datos_SC %>% sapply(FUN = function(x)sd(x)/mean(x)) %>% hist()
datos_SC %>% sapply(FUN = function(x)sd(x)/mean(x)) %>% boxplot()
datos_SC %>% sapply(FUN = function(x)var(x/max(x))) %>% summary()
datos_SC %>% sapply(FUN = function(x)var(x/max(x))) %>% hist()
datos_SC %>% sapply(FUN = function(x)var(x/max(x))) %>% boxplot()
datos_SC %>% sapply(FUN = function(x)sd(x/max(x))/mean(x/max(x))) %>% summary()
datos_SC %>% sapply(FUN = function(x)sd(x/max(x))/mean(x/max(x))) %>% hist()
datos_SC %>% sapply(FUN = function(x)sd(x/max(x))/mean(x/max(x))) %>% boxplot()


datos_EC %>% apply(MARGIN = 2, FUN = var) %>% summary()
datos_EC %>% apply(MARGIN = 2, FUN = var) %>% hist()
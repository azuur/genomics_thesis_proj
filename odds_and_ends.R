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

install.packages("graph")

#GENIE3
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GENIE3")
#BiocManager::install("GENIE3", version = "3.8")


#NARROMI
#sudo apt-get install libglpk-dev
install.packages("Rglpk")
source("narromi.R")

#TIGRESS
devtools::install_github("jpvert/tigress")

#MRNET
#BiocManager::install("minet", version = "3.8")
BiocManager::install("minet")
install.packages("parmigene")
devtools::install_bitbucket('Jonathan-Ish-Horowicz/fastGeneMI')

#CLR
#BiocManager::install("minet", version = "3.8")
BiocManager::install("minet")
install.packages("parmigene")
devtools::install_bitbucket('Jonathan-Ish-Horowicz/fastGeneMI')

#ARACNE
#BiocManager::install("minet", version = "3.8")
BiocManager::install("minet")
install.packages("parmigene")
devtools::install_bitbucket('Jonathan-Ish-Horowicz/fastGeneMI')

#PC ALG
install.packages("bnlearn")
#sudo apt-get install libv8-dev
BiocManager::install("graph")
BiocManager::install("RBGL")
install.packages("pcalg")
install.packages("ParallelPC")

#CORR/MI
#minet
devtools::install_bitbucket('Jonathan-Ish-Horowicz/fastGeneMI')


#BDEU
install.packages("bnlearn")
install.packages("pcalg")


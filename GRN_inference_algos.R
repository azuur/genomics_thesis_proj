install.packages("graph")

#GENIE3
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GENIE3", version = "3.8")

#NARROMI


#TIGRESS
devtools::install_github("jpvert/tigress")

#MRNET
BiocManager::install("minet", version = "3.8")
install.packages("parmigene")
devtools::install_bitbucket('Jonathan-Ish-Horowicz/fastGeneMI')

#CLR
BiocManager::install("minet", version = "3.8")
install.packages("parmigene")
devtools::install_bitbucket('Jonathan-Ish-Horowicz/fastGeneMI')

#ARACNE
BiocManager::install("minet", version = "3.8")
install.packages("parmigene")
devtools::install_bitbucket('Jonathan-Ish-Horowicz/fastGeneMI')

#PC ALG
install.packages("bnlearn")
install.packages("pcalg")

#CORR/MI
#minet
devtools::install_bitbucket('Jonathan-Ish-Horowicz/fastGeneMI')


#BDEU
install.packages("bnlearn")
install.packages("pcalg")


#test
require(here)
source("pckgs_and_useful_wrappers.R")
source("narromi.R")
require(parallel)
require(doParallel)
require(doRNG)
readobj(y1)
colnames(y1) <- paste0("Y",1:ncol(y1))






results_y1<-list()
time_y1 <- list()


#GENIE3
time_y1[["GENIE3"]] <- system.time({
  GENIE3 <- GENIE3::GENIE3(exprMatrix = t(y1),nCores = 5)
saveobj(GENIE3)
rm(GENIE3)
})
print(time_y1)
#NARROMI
cl <- makeCluster(5)
time_y1[["NARROMI"]] <- system.time({
  NARROMI <- NARROMI(X = y1, cl = cl)
saveobj(NARROMI)
rm(NARROMI)
})
stopCluster(cl)
rm(cl)
print(time_y1)

#TIGRESS
time_y1[["TIGRESS"]] <- system.time({
  TIGRESS <- tigress::tigress(expdata = y1, usemulticore = T)
saveobj(TIGRESS)
rm(TIGRESS)
})
print(time_y1)

#MI matrix
#BiocManager::install("minet", version = "3.8")
time_y1[["mi_minet"]] <- system.time({
  mi_minet <- minet::build.mim(dataset = y1, estimator = "mi.mm", disc = "equalwidth", nbins = nrow(y1)^(1/3))
saveobj(mi_minet)
rm(mi_minet)
})
print(time_y1)
time_y1[["mi_parmigene"]] <- system.time({
  mi_parmigene <- parmigene::knnmi.all(t(y1))
saveobj(mi_parmigene)
rm(mi_parmigene)
})
print(time_y1)
time_y1[["mi_fastGeneMI_bsp"]] <- system.time({
  mi_fastGeneMI_bsp <- fastGeneMI::get.mim.bspline(y1, order = 4, n.cores = 5)
saveobj(mi_fastGeneMI_bsp)
rm(mi_fastGeneMI_bsp)
})
print(time_y1)
time_y1[["mi_fastGeneMI_mm"]] <- system.time({
  mi_fastGeneMI_mm <- fastGeneMI::get.mim.MM(y1, discretisation = "equalwidth", n.cores = 5)
saveobj(mi_fastGeneMI_mm)
rm(mi_fastGeneMI_mm)
})
print(time_y1)


#MRNET
time_y1[["MRNET_minet"]] <- system.time({
  MRNET_minet <- minet::mrnet(results_y1[["mi_minet"]])
saveobj(MRNET_minet)
rm(MRNET_minet)
  })
print(time_y1)
time_y1[["MRNET_parmigene"]] <- system.time({
  MRNET_parmigene <- parmigene::mrnet(results_y1[["mi_parmigene"]])
saveobj(MRNET_parmigene)
rm(MRNET_parmigene)
})
print(time_y1)
time_y1[["MRNET_parmigene_fgm_bsp"]] <- system.time({
  MRNET_parmigene_fgm_bsp <- parmigene::mrnet(results_y1[["mi_fastGeneMI_bsp"]]) 
saveobj(MRNET_parmigene_fgm_bsp)
rm(MRNET_parmigene_fgm_bsp)
  })
print(time_y1)
time_y1[["MRNET_parmigene_fgm_mm"]] <- system.time({
  MRNET_parmigene_fgm_mm <- parmigene::mrnet(results_y1[["mi_fastGeneMI_mm"]])
saveobj(MRNET_parmigene_fgm_mm)
rm(MRNET_parmigene_fgm_mm)
  })
print(time_y1)
#CLR
#BiocManager::install("minet", version = "3.8")
time_y1[["CLR_minet"]] <- system.time({
  CLR_minet <- minet::clr(results_y1[["mi_minet"]])
saveobj(CLR_minet)
rm(CLR_minet)
  })
print(time_y1)
time_y1[["CLR_parmigene"]] <- system.time({
  CLR_parmigene <- parmigene::clr(results_y1[["mi_parmigene"]])
saveobj(CLR_parmigene)
rm(CLR_parmigene)
  })
print(time_y1)
time_y1[["CLR_parmigene_fgm_bsp"]] <- system.time({
  CLR_parmigene_fgm_bsp <- parmigene::clr(results_y1[["mi_fastGeneMI_bsp"]]) 
saveobj(CLR_parmigene_fgm_bsp)
rm(CLR_parmigene_fgm_bsp)
  })
print(time_y1)
time_y1[["CLR_parmigene_fgm_mm"]] <- system.time({
  CLR_parmigene_fgm_mm <- parmigene::clr(results_y1[["mi_fastGeneMI_mm"]])
saveobj(CLR_parmigene_fgm_mm)
rm(CLR_parmigene_fgm_mm)
  })
print(time_y1)

#ARACNE
#BiocManager::install("minet", version = "3.8")
time_y1[["ARACNE_minet"]] <- system.time({
  ARACNE_minet <-  minet::aracne(results_y1[["mi_minet"]], eps = 0.05)
saveobj(ARACNE_minet)
rm(ARACNE_minet)
  })
print(time_y1)
time_y1[["ARACNE_parmigene"]] <- system.time({
  ARACNE_parmigene <- parmigene::aracne.a(results_y1[["mi_parmigene"]])
saveobj(ARACNE_parmigene)
rm(ARACNE_parmigene)
  })
print(time_y1)
time_y1[["ARACNE_parmigene_fgm_bsp"]] <- system.time({
  ARACNE_parmigene_fgm_bsp <- parmigene::aracne.a(results_y1[["mi_fastGeneMI_bsp"]]) 
saveobj(ARACNE_parmigene_fgm_bsp)
rm(ARACNE_parmigene_fgm_bsp)
  })
print(time_y1)
time_y1[["ARACNE_parmigene_fgm_mm"]] <- system.time({
  ARACNE_parmigene_fgm_mm <- parmigene::aracne.a(results_y1[["mi_fastGeneMI_mm"]])
saveobj(ARACNE_parmigene_fgm_mm)
rm(ARACNE_parmigene_fgm_mm)
  })
print(time_y1)

#PC ALG
cl <- makeCluster(5)
time_y1[["PCstable_0"]] <- system.time({
  PCstable_0 <- bnlearn::pc.stable(as.data.frame(y1[,173:193]), cluster = cl)
saveobj(PCstable_0)
rm(PCstable_0)
  })
stopCluster(cl)
rm(cl)
print(time_y1)

cl <- makeCluster(5)
time_y1[["PCstable_1"]] <- system.time({
  PCstable_1 <- bnlearn::pc.stable(as.data.frame(y1), cluster = cl,test = "zf")
saveobj(PCstable_1)
rm(PCstable_1)
  })
stopCluster(cl)
rm(cl)
print(time_y1)

cl <- makeCluster(5)
time_y1[["PCstable_2"]] <- system.time({
  PCstable_2 <- bnlearn::pc.stable(as.data.frame(y1), cluster = cl,test =  "mi-g")
saveobj(PCstable_2)
rm(PCstable_2)
  })
stopCluster(cl)
rm(cl)
print(time_y1)



#BDEU
time_y1[["BDEU"]] <- system.time({
  supl_data <- bnlearn::discretize(as.data.frame(y1))
  BDEU <- bnlearn::hc(supl_data, score = "bde", optimized = T)
  saveobj(BDEU)
  rm(BDEU)
})
print(time_y1)




saveobj(time_y1)




















results_y2<-list()
time_y2 <- list()

#GENIE3
require(parallel)
require(doParallel)
require(doRNG)
time_y2[["GENIE3"]] <- system.time({
  GENIE3 <- GENIE3::GENIE3(exprMatrix = t(y2),nCores = 5)
  saveobj(GENIE3)
  rm(GENIE3)
})
print(time_y2)
#NARROMI
cl <- makeCluster(5)
time_y2[["NARROMI"]] <- system.time({
  NARROMI <- NARROMI(X = y2, cl = cl)
  saveobj(NARROMI)
  rm(NARROMI)
})
stopCluster(cl)
rm(cl)
print(time_y2)

#TIGRESS
time_y2[["TIGRESS"]] <- system.time({
  TIGRESS <- tigress::tigress(expdata = y2, usemulticore = T)
  saveobj(TIGRESS)
  rm(TIGRESS)
})
print(time_y2)

#MI matrix
#BiocManager::install("minet", version = "3.8")
time_y2[["mi_minet"]] <- system.time({
  mi_minet <- minet::build.mim(dataset = y2, estimator = "mi.mm", disc = "equalwidth", nbins = nrow(y2)^(1/3))
  saveobj(mi_minet)
  rm(mi_minet)
})
print(time_y2)
time_y2[["mi_parmigene"]] <- system.time({
  mi_parmigene <- parmigene::knnmi.all(t(y2))
  saveobj(mi_parmigene)
  rm(mi_parmigene)
})
print(time_y2)
time_y2[["mi_fastGeneMI_bsp"]] <- system.time({
  mi_fastGeneMI_bsp <- fastGeneMI::get.mim.bspline(y2, order = 4, n.cores = 5)
  saveobj(mi_fastGeneMI_bsp)
  rm(mi_fastGeneMI_bsp)
})
print(time_y2)
time_y2[["mi_fastGeneMI_mm"]] <- system.time({
  mi_fastGeneMI_mm <- fastGeneMI::get.mim.MM(y2, discretisation = "equalwidth", n.cores = 5)
  saveobj(mi_fastGeneMI_mm)
  rm(mi_fastGeneMI_mm)
})
print(time_y2)


#MRNET
time_y2[["MRNET_minet"]] <- system.time({
  MRNET_minet <- minet::mrnet(results_y2[["mi_minet"]])
  saveobj(MRNET_minet)
  rm(MRNET_minet)
})
print(time_y2)
time_y2[["MRNET_parmigene"]] <- system.time({
  MRNET_parmigene <- parmigene::mrnet(results_y2[["mi_parmigene"]])
  saveobj(MRNET_parmigene)
  rm(MRNET_parmigene)
})
print(time_y2)
time_y2[["MRNET_parmigene_fgm_bsp"]] <- system.time({
  MRNET_parmigene_fgm_bsp <- parmigene::mrnet(results_y2[["mi_fastGeneMI_bsp"]]) 
  saveobj(MRNET_parmigene_fgm_bsp)
  rm(MRNET_parmigene_fgm_bsp)
})
print(time_y2)
time_y2[["MRNET_parmigene_fgm_mm"]] <- system.time({
  MRNET_parmigene_fgm_mm <- parmigene::mrnet(results_y2[["mi_fastGeneMI_mm"]])
  saveobj(MRNET_parmigene_fgm_mm)
  rm(MRNET_parmigene_fgm_mm)
})
print(time_y2)
#CLR
#BiocManager::install("minet", version = "3.8")
time_y2[["CLR_minet"]] <- system.time({
  CLR_minet <- minet::clr(results_y2[["mi_minet"]])
  saveobj(CLR_minet)
  rm(CLR_minet)
})
print(time_y2)
time_y2[["CLR_parmigene"]] <- system.time({
  CLR_parmigene <- parmigene::clr(results_y2[["mi_parmigene"]])
  saveobj(CLR_parmigene)
  rm(CLR_parmigene)
})
print(time_y2)
time_y2[["CLR_parmigene_fgm_bsp"]] <- system.time({
  CLR_parmigene_fgm_bsp <- parmigene::clr(results_y2[["mi_fastGeneMI_bsp"]]) 
  saveobj(CLR_parmigene_fgm_bsp)
  rm(CLR_parmigene_fgm_bsp)
})
print(time_y2)
time_y2[["CLR_parmigene_fgm_mm"]] <- system.time({
  CLR_parmigene_fgm_mm <- parmigene::clr(results_y2[["mi_fastGeneMI_mm"]])
  saveobj(CLR_parmigene_fgm_mm)
  rm(CLR_parmigene_fgm_mm)
})
print(time_y2)

#ARACNE
#BiocManager::install("minet", version = "3.8")
time_y2[["ARACNE_minet"]] <- system.time({
  ARACNE_minet <-  minet::aracne(results_y2[["mi_minet"]], eps = 0.05)
  saveobj(ARACNE_minet)
  rm(ARACNE_minet)
})
print(time_y2)
time_y2[["ARACNE_parmigene"]] <- system.time({
  ARACNE_parmigene <- parmigene::aracne.a(results_y2[["mi_parmigene"]])
  saveobj(ARACNE_parmigene)
  rm(ARACNE_parmigene)
})
print(time_y2)
time_y2[["ARACNE_parmigene_fgm_bsp"]] <- system.time({
  ARACNE_parmigene_fgm_bsp <- parmigene::aracne.a(results_y2[["mi_fastGeneMI_bsp"]]) 
  saveobj(ARACNE_parmigene_fgm_bsp)
  rm(ARACNE_parmigene_fgm_bsp)
})
print(time_y2)
time_y2[["ARACNE_parmigene_fgm_mm"]] <- system.time({
  ARACNE_parmigene_fgm_mm <- parmigene::aracne.a(results_y2[["mi_fastGeneMI_mm"]])
  saveobj(ARACNE_parmigene_fgm_mm)
  rm(ARACNE_parmigene_fgm_mm)
})
print(time_y2)

#PC ALG
cl <- makeCluster(5)
time_y2[["PCstable_0"]] <- system.time({
  PCstable_0 <- bnlearn::pc.stable(as.data.frame(y2[,1:100]), cluster = cl)
  saveobj(PCstable_0)
  rm(PCstable_0)
})
stopCluster(cl)
rm(cl)
print(time_y2)

cl <- makeCluster(5)
time_y2[["PCstable_1"]] <- system.time({
  PCstable_1 <- bnlearn::pc.stable(as.data.frame(y2), cluster = cl,test = "zf")
  saveobj(PCstable_1)
  rm(PCstable_1)
})
stopCluster(cl)
rm(cl)
print(time_y2)

cl <- makeCluster(5)
time_y2[["PCstable_2"]] <- system.time({
  PCstable_2 <- bnlearn::pc.stable(as.data.frame(y2), cluster = cl,test =  "mi-g")
  saveobj(PCstable_2)
  rm(PCstable_2)
})
stopCluster(cl)
rm(cl)
print(time_y2)



#BDEU
time_y2[["BDEU"]] <- system.time({
  supl_data <- bnlearn::discretize(as.data.frame(y2))
  BDEU <- bnlearn::hc(supl_data, score = "bde", optimized = T)
  saveobj(BDEU)
  rm(BDEU)
})
print(time_y2)




saveobj(time_y2)



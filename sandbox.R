i <- 1

path <- file.path("/media","adrian","bodega","thesis",
                  "Robjects","simulated_data",
                  netw_str,
                  type,
                  paste0(noise,"_noise"),
                  paste0(r,"_noise_to_sig"),
                  paste0("sample_size_",sample_size)
)
sim1 <- 
  readRDS(
        paste0(file.path(path,"sim"),i,".RDS")
)


npath <- file.path("/media","adrian","bodega","thesis",
                  "Robjects","simulated_data",
                  netw_str,
                  type,
                  paste0(noise,"_noise"),
                  paste0(r,"_noise_to_sig"),
                  paste0("sample_size_",sample_size,"_new_test")
                  )
sim1n <- 
  readRDS(
    paste0(file.path(npath,"sim"),i,".RDS")
  )



rm(list=ls())
Sys.time()
source("estimation/estimates.R")
rm(list=ls())
Sys.time()
source("analysis/probabilities.R")
rm(list=ls())
Sys.time()
source("analysis/curves.R")
rm(list=ls())
Sys.time()
source("analysis/puttogether.R")
rm(list=ls())
Sys.time()
source("analysis/inferred_nets.R")



data_path <- file.path("/media","adrian","bodega","thesis",
                       "Robjects","simulated_data",
                       netw_str,
                       type,
                       paste0(noise,"_noise"),
                       paste0(r,"_noise_to_sig"),
                       paste0("sample_size_",sample_size)
)

dat1 <- readRDS(file = file.path(data_path,paste0("sim183.RDS")))

data_path <- file.path("/media","adrian","bodega","thesis",
                       "Robjects","simulated_data",
                       netw_str,
                       type,
                       paste0(noise,"_noise"),
                       paste0(r,"_noise_to_sig"),
                       paste0("sample_size_",sample_size,"_old")
)

dat2 <- readRDS(file = file.path(data_path,paste0("sim183.RDS")))

dat1[1,1]
dat2[1,1]
sum(dat1==dat2)/length(dat1)

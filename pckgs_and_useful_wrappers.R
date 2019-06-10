require(here)
require(tidyverse)
require(magrittr)
require(igraph)
require(digest)




filelist <- function(subdir = "Robjects"){
  if(require(here)){
    str <- list.files(path = here::here(subdir), pattern = ".RDS", ignore.case = T)
  } else{
    str <- list.files(path = paste0(getwd(),"/",subdir), pattern = ".RDS", ignore.case = T)
  }
  gsub(".RDS","",str)
}

str_saveobj <-function(x, subdir = "Robjects"){
  if(require(here)){
    dir.create(here::here(subdir), showWarnings = FALSE, recursive = T)
    assign("obj",get(x,envir = .GlobalEnv),envir = environment())
    saveRDS(obj,here::here(paste0(subdir,"/",x,".RDS")))
  } else{
    dir.create(paste0(getwd(),"/",subdir), showWarnings = FALSE, recursive = T)
    assign("obj",get(x,envir = .GlobalEnv),envir = environment())
    saveRDS(obj,paste0(getwd(),"/",subdir,"/",x,".RDS"))
  }
} 

str_readobj <-function(x, subdir = "Robjects"){
  if(require(here)){
    assign(x,readRDS(here::here(paste0(subdir,"/",x,".RDS"))),envir = .GlobalEnv)
  } else {
    assign(x,readRDS(paste0(getwd(),"/",subdir,"/",x,".RDS")),envir = .GlobalEnv)
  }
} 

str_savelist <- function(..., subdir = "Robjects"){
  l <- c(...)
  for(i in l) str_saveobj(i, subdir = subdir)
}

str_readlist <- function(..., subdir = "Robjects"){
  l <- c(...)
  if( is.null(l) ){ l <- filelist(subdir)}
  for(i in l) str_readobj(i, subdir = subdir)
}

saveobj <- function(..., subdir = "Robjects"){
  strs <- sapply(substitute(list(...))[-1], deparse)
  str_savelist(strs, subdir = subdir)
}

readobj <- function(..., subdir = "Robjects"){
  strs <- sapply(substitute(list(...))[-1], deparse)
  if( length(strs) == 0 ){ strs <- NULL}
  str_readlist(strs, subdir = subdir)
}

calibrate_variances <- function(graph, W, r = 0.5, ini_var = 1, type = "linear"){
  
  edges <- igraph::get.edgelist(graph, names=F)
  ord <- as.numeric(topo_sort(graph,mode = "out"))
  vs<-rep(ini_var,length(ord))
  
  if(type == "linear"){
    for(x in ord){
      parents_x <- edges[edges[,2]==x,1]
      if(length(parents_x)>0){
        B_inv <- solve(diag(n)-W)
        vs[x] <- (r/(1-r))*(W[x,parents_x]%*%
                              (B_inv%*%diag(vs)%*%t(B_inv))[parents_x,parents_x]%*%
                              W[x,parents_x])
      }
    }
  } else if(type == "sigmoid"){
    for(x in ord){
      parents_x <- edges[edges[,2]==x,1]
      if(length(parents_x)>0){
        B_inv <- solve(diag(n)-W/2)
        vs[x] <- (1/4)*(r/(1-r))*(W[x,parents_x]%*%
                                    (B_inv%*%diag(vs)%*%t(B_inv))[parents_x,parents_x]%*%
                                    W[x,parents_x])
      }
    }
  } else{ stop("Currently only 'linear' and 'sigmoid' types supported.")}
  vs[is.na(vs)] <- mean(vs,na.rm = T)/r
  vs
}

#modified from Ben Bolver @
#https://stackoverflow.com/questions/10910698/questions-about-set-seed-in-r
make.seed <- function(x) {
  require("digest")
  hexval <- paste0("0x",digest(x,"crc32"))
  intval <- type.convert(hexval) %% 2147483647#.Machine$integer.max
  intval
}

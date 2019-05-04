require(here)
require(tidyverse)
require(magrittr)
require(igraph)


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
    dir.create(here::here(subdir), showWarnings = FALSE)
    assign("obj",get(x,envir = .GlobalEnv),envir = environment())
    saveRDS(obj,here::here(paste0(subdir,"/",x,".RDS")))
  } else{
    dir.create(paste0(getwd(),"/",subdir), showWarnings = FALSE)
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
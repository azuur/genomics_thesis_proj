NARROMI <- function(X, lambda = 1, alpha = 0.05, beta = 0.05, t = 0.6, cl =NULL){

  if(is.null(cl)){
    l <- lapply(1:ncol(X),
                function(i){
                  res <- narromi(X[,i,drop=F], X[,-i, drop=F], lambda, alpha, beta, t)
                  net <- split(res$net, 1:length(res$net)<i)
                  net <- c(net[["TRUE"]],0,net[["FALSE"]])
                  net_value <- split(res$net_value, 1:length(res$net_value)<i)
                  net_value <- c(net_value[["TRUE"]],0,net_value[["FALSE"]])
                  sig <- split(res$sig, 1:length(res$sig)<i)
                  sig <- c(sig[["TRUE"]],0,sig[["FALSE"]])
                  return(list(net = net, net_value = net_value, sig = sig))
                })
  } else{
    clusterExport(cl, 
                  c("X", "lambda", "alpha", "beta", "t"), 
                  envir = environment()
                  )
    
    clusterExport(cl, 
                  c("narromi", "narromi_cmi","reoptim","LP_TGN"), 
                  envir = .GlobalEnv
    )
    l <- parLapply(cl,
                   1:ncol(X),
                   function(i){
                     res <- narromi(X[,i,drop=F], X[,-i,drop=F], lambda, alpha, beta, t)
                     net <- split(res$net, 1:length(res$net)<i)
                     net <- c(net[["TRUE"]],0,net[["FALSE"]])
                     net_value <- split(res$net_value, 1:length(res$net_value)<i)
                     net_value <- c(net_value[["TRUE"]],0,net_value[["FALSE"]])
                     sig <- split(res$sig, 1:length(res$sig)<i)
                     sig <- c(sig[["TRUE"]],0,sig[["FALSE"]])
                     return(list(net = net, net_value = net_value, sig = sig))
                   })
  }
  
  net <- do.call(rbind, lapply(l, function(x) x$net))
  net_value <- do.call(rbind, lapply(l, function(x) x$net_value))
  sig <- do.call(rbind, lapply(l, function(x) x$sig))
  return(list(net = net, net_value = net_value, sig = sig))
}

narromi <- function(y, X, lambda = 1, alpha = 0.05, beta = 0.05, t = 0.6){
  net <- matrix(0, nrow = 1, ncol = ncol(X))
  
  net_value <- 
    G1 <- 
    apply(X = X, MARGIN = 2, FUN = function(x) narromi_cmi(x,y))
#    apply(X = X, MARGIN = 2, FUN = function(x) infotheo::mutinformation(x,y))
    
  idx <- which(abs(G1) >= alpha)
  if( length(idx) == 0){ return(
    list(net = rep(0,ncol(X)), net_value = rep(0,ncol(X)), sig = rep(NA,ncol(X)))
    )}
  X1 <- X[,idx, drop=F]
  
  J <- reoptim(y, X1, lambda, beta)
  
  net[idx] <- J$sparse
  net_value[idx] <- J$value
  
  z <- abs(net_value)
  z <- (z - min(z))/(max(z) - min(z))
  z <- 0.5*log2((1 + z)/(1 - z))
  z[z == max(z)] <- max(z[z != max(z)])
  sig1 <- 1 - pnorm(z, mean(z), sd(z))
  
  z <- abs(G1)
  z <- (z - min(z))/(max(z) - min(z))
  z <- 0.5*log2((1 + z)/(1 - z))
  z[z == max(z)] <- max(z[z != max(z)])
  sig2 <- 1 - pnorm(z, mean(z), sd(z))
  
  sig <- sqrt(sig1^2 + sig2^2)
  
  net_value <- sign(net_value)*(abs(net_value)*t + G1*(1-t))
  
  
  return(list(net = net, net_value = net_value, sig = sig))
}


narromi_cmi <- function(...){
  args <- list(...)
  
  if(length(args) == 2){
#    cmiv <- -0.5*log(1-cor(args[[1]],args[[2]])^2)
    c1 <- det(as.matrix(cov(as.matrix(args[[1]]))))
    c2 <- det(as.matrix(cov(as.matrix(args[[2]]))))
    c3 <- det(cov(do.call(cbind,args)))
    cmiv <- 0.5*log(c1*c2/c3)
  } else if(length(args) == 3){
    c1 <- det(cov(cbind(args[[1]],args[[3]])))
    c2 <- det(cov(cbind(args[[2]],args[[3]])))
    c3 <- det(as.matrix(cov(as.matrix(args[[3]]))))
    c4 <- det(cov(do.call(cbind,args)))
    cmiv <- 0.5*log((c1*c2)/(c3*c4))
  }
  
  if( is.infinite(cmiv) ){cmiv <- 1}
  cmiv <- abs(cmiv)
  
  return(cmiv)
}


reoptim <- function(y, X, lambda, alpha){
  
  J <- LP_TGN(y, X, lambda)
  sparse <- 
    value <-
    matrix(0, ncol = ncol(J), nrow = nrow(J))
  
  
  idx <- which(abs(J) >= alpha)
  idx_c <- which(abs(J) < alpha)
  if(length(idx)>0){
    
    J1 <- LP_TGN(y, X[,idx, drop=F], lambda)
    
    idx1 <- which(abs(J1) >= alpha)
    idx_c1 <- which(abs(J1) < alpha)
    
    #if(length(idx1)==0){break}
    
    old_idx <- idx
    
    idx <- old_idx[idx1]
    idx_c <- old_idx[idx_c1]
    
    sparse[idx] <- J1[idx1]
    sparse[idx_c] <- 0
    value[idx_c] <- J1[idx_c1]
    
    while(length(idx_c) > 0 & length(idx1) > 0){
      
      J1 <- LP_TGN(y, X[,idx, drop=F], lambda)
      
      idx1 <- which(abs(J1) >= alpha)
      idx_c1 <- which(abs(J1) < alpha)
      
      old_idx <- idx
      
      idx <- old_idx[idx1]
      idx_c <- old_idx[idx_c1]
      
      sparse[idx] <- J1[idx1]
      sparse[idx_c] <- 0
      value[idx_c] <- J1[idx_c1]
    }
    
    value <- value + sparse
  }
  
  
  return(list(sparse = sparse, value = value))
}


#this works the same as MTE::LADLasso, with lambdaLP_TGN = n*lambda_LADLasso 
#THAT WAS THE MISTAKE WITH OUR PREVIOUS VERSION OF NARROMI...
#ALSO... LP_TGN is 4-6 times faster than MTE::LADLasso
LP_TGN <- function(Y, X, lambda){
  
  if(ncol(X)==0){stop("Empty X matrix provided.")}
  
  if( is.matrix(Y) ){ reps = ncol(Y) } else{ reps = 1 }

  A <- cbind(Matrix::sparseMatrix(i = 1:length(Y), j = 1:length(Y), x = 1),
             Matrix::sparseMatrix(i = 1:length(Y), j = 1:length(Y), x = -1),
             Matrix::bdiag( rep(list(X), reps) ),
             Matrix::bdiag( rep(list(-X), reps) ))
  
  n_aux_vars <- ncol(A)
  
  A <- rbind(A,
             Matrix::sparseMatrix(i = 1:n_aux_vars, j = 1:n_aux_vars, x = 1)
             )
  
  b <- c(Y, rep(0,n_aux_vars))
  
  f <- c(rep(1, 2*length(Y)),
         rep(lambda, 2*ncol(X)*reps)
  )
  
  symb <- c(rep("==", length(Y)),rep(">=", n_aux_vars))
  
  x <- Rglpk::Rglpk_solve_LP(obj = f, mat = A, dir = symb, rhs = b,
                             control = list("tm_limit" = 60*1e3))
  
  if(x$status > 0){
    if(is.matrix(Y)){
      if(ncol(Y) > 1){
        warning("Optimal solution not found with GLPK solver. Returning 0 vector.")
        x$solution[-c( 1:(2*length(Y)) )] <- 0
      } else{ Y <- c(Y) }
    }
    warning("Optimal solution not found with GLPK solver. Used MTE::LADlasso.")
    x$solution[-c( 1:(2*length(Y)) )] <- MTE::LADlasso(
      Y, 
      X, 
      rep(1,ncol(X)), 
      lambda = lambda/nrow(X),
      adaptive = F)$beta
  }

  res <- x$solution[-c( 1:(2*length(Y)) )]
  res <- res[1:(length(res)/2)] + res[(length(res)/2 + 1):length(res)]
  res <- matrix(res, ncol = ncol(X), byrow = T)
  res
}


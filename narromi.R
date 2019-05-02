narromi <- function(y, X, lambda, alpha, beta, t){
  net <- matrix(0, nrow = 1, ncol = ncol(X))
  
  net_value <- 
    G1 <- 
    apply(X = X, MARGIN = 2, FUN = function(x) narromi_cmi(x,y))
#    apply(X = X, MARGIN = 2, FUN = function(x) infotheo::mutinformation(x,y))
    
  idx <- which(abs(G1) >= alpha)
  if( length(idx) == 0){ return(
    list(net = rep(0,length(y)), net_value = rep(0,length(y)), sig = rep(NA,length(y)))
    )}
  X1 <- X[,idx]
  
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
  
  J1 <- LP_TGN(y, X[,idx], lambda)
  
  idx1 <- which(abs(J1) >= alpha)
  idx_c1 <- which(abs(J1) < alpha)
  
  old_idx <- idx
  
  idx <- old_idx[idx1]
  idx_c <- old_idx[idx_c1]
  
  sparse[idx] <- J1[idx1]
  sparse[idx_c] <- 0
  value[idx_c] <- J1[idx_c1]
  
  while(length(idx_c) > 0){
    
    J1 <- LP_TGN(y, X[,idx], lambda)
    
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
  
  return(list(sparse = sparse, value = value))
}


#this works the same as MTE::LADLasso, with lambdaLP_TGN = n*lambda_LADLasso 
#THAT WAS THE MISTAKE WITH OUR PREVIOUS VERSION OF NARROMI...
#ALSO... LP_TGN is 4-6 times faster than MTE::LADLasso
LP_TGN <- function(Y, X, lambda){
  
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
  
  x <- Rglpk::Rglpk_solve_LP(obj = f, mat = A, dir = symb, rhs = b)

  res <- x$solution[-c( 1:(2*length(Y)) )]
  res <- res[1:(length(res)/2)] + res[(length(res)/2 + 1):length(res)]
  res <- matrix(res, ncol = ncol(X), byrow = T)
  res
}


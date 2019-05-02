#linear_regulation
#Use inv(I- coef*ADJ) matrix



library(igraph)



#build closure for individual node
one_gene_closure <- function(x, edges, functional_form, params = NULL){
  if(!functional_form %in% c("linear","additive","flattened_line","sigmoid","quadratic","cos_sigmoid","sin_sigmoid")){
    stop("Unknown functional form.")
  }

  if(!is.null(params$coefs) & !is.null(params$coef_distribution)){ stop("Only one of coefs/coef_distribution must be supplied.")}
  
  
  parents_x <- edges[edges[,2]==x,1]
  
  if(functional_form =="linear"){
    
    if(is.null(params$coefs)){ 
      coefficients <- params$coef_distribution(length(parents_x)) 
      coefficient_origin = params$coef_distribution
    } else{
      coefficients <- rep_len(params$coefs, length(parents_x))
      coefficient_origin = params$coefs
    }
    
    reg_x <- function(vec,errors){ 
      sum(coefficients*vec[parents_x]) + errors[x]
    }
    
    metadata <- list(functional_form = functional_form, 
                     parent_nodes = parents_x,
                     coefficients = coefficients,
                     coeficient_origin = coefficient_origin
    )
    
  } else if(functional_form =="additive"){

    if(is.null(params$coefs)){ 
      coefficients <- params$coef_distribution(length(parents_x)) 
      coefficient_origin = params$coef_distribution
    } else{
      coefficients <- rep_len(params$coefs, length(parents_x))
      coefficient_origin = params$coefs
    }
    
    fs <- rep_len(params$functions_add, length(parents_x))
    
    reg_x <- function(vec,errors){ 
      sum(coefficients*
            sapply(seq_along(parents_x),FUN = function(i) fs[[i]](vec[parents_x][i]))
          ) + errors[x]
    }
    
    metadata <- list(functional_form = functional_form, 
                     parent_nodes = parents_x,
                     coefficients = coefficients,
                     functions_to_add = fs,
                     coeficient_origin = coefficient_origin
    )
    
  } else if(functional_form =="flattened_line"){

    if(is.null(params$coefs)){ 
      coefficients <- params$coef_distribution(length(parents_x)) 
      coefficient_origin = params$coef_distribution
    } else{
      coefficients <- rep_len(params$coefs, length(parents_x))
      coefficient_origin = params$coefs
    }

    reg_x <- function(vec,errors){ 
      v <- sum(coefficients*vec[parents_x]) 
      v*(v <= params$threshold) + params$threshold*(v > params$threshold) + errors[x]
    }
    
    metadata <- list(functional_form = functional_form, 
                     parent_nodes = parents_x,
                     coefficients = coefficients,
                     threshold = params$threshold,
                     coeficient_origin = coefficient_origin
    )
    
  } else if(functional_form =="sigmoid"){

    if(is.null(params$coefs)){ 
      coefficients <- params$coef_distribution(length(parents_x)) 
      coefficient_origin = params$coef_distribution
    } else{
      coefficients <- rep_len(params$coefs, length(parents_x))
      coefficient_origin = params$coefs
    }
    
    reg_x <- function(vec,errors){ 
      v <- sum(coefficients*vec[parents_x]) + errors[x]
      1/(1+ exp(-v))
    }
    
    metadata <- list(functional_form = functional_form, 
                     parent_nodes = parents_x,
                     coefficients = coefficients,
                     coeficient_origin = coefficient_origin
    )
    
  } else if(functional_form =="quadratic"){

    if(is.null(params$coefs)){ 
      coefficients <- params$coef_distribution(length(parents_x)) 
      coefficient_origin = params$coef_distribution
    } else{
      coefficients <- rep_len(params$coefs, length(parents_x))
      coefficient_origin = params$coefs
    }
    
    coefficients <- matrix(coefficients,nrow=length(parents_x))
    
    reg_x <- function(vec,errors){ 
      t(vec[parents_x])%*%coefficients%*%vec[parents_x] + errors[x]
    }
    
    metadata <- list(functional_form = functional_form, 
                     parent_nodes = parents_x,
                     coefficients = coefficients,
                     coeficient_origin = coefficient_origin
    )
    
  } else if(functional_form =="cos_sigmoid"){
    
    if(is.null(params$coefs)){ 
      coefficients <- params$coef_distribution(length(parents_x)) 
      coefficient_origin = params$coef_distribution
    } else{
      coefficients <- rep_len(params$coefs, length(parents_x))
      coefficient_origin = params$coefs
    }
    
    reg_x <- function(vec,errors){ 
      v <- sum(coefficients*vec[parents_x]) + errors[x]
      cos(params$domain_max/(1+ exp(-v)))
    }
    
    metadata <- list(functional_form = functional_form, 
                     parent_nodes = parents_x,
                     domain_max = params$domain_max, 
                     coefficients = coefficients,
                     coeficient_origin = coefficient_origin
    )
    
  } else if(functional_form =="sin_sigmoid"){
    
    if(is.null(params$coefs)){ 
      coefficients <- params$coef_distribution(length(parents_x)) 
      coefficient_origin = params$coef_distribution
    } else{
      coefficients <- rep_len(params$coefs, length(parents_x))
      coefficient_origin = params$coefs
    }
    
    reg_x <- function(vec,errors){ 
      v <- sum(coefficients*vec[parents_x]) + errors[x]
      sin(params$domain_max/(1+ exp(-v)))
    }
    
    metadata <- list(functional_form = functional_form, 
                     parent_nodes = parents_x,
                     domain_max = params$domain_max, 
                     coefficients = coefficients,
                     coeficient_origin = coefficient_origin
    )
    
  }
  
  if(length(parents_x)==0){
    reg_x <- function(vec,errors){ errors[x] }
  }
  
  list(closure = reg_x, metadata = metadata)
}



#build closure for whole graph
all_gene_closures <- function(graph, functional_forms_list, param_list){
  edges <- igraph::get.edgelist(graph, names=F)

  l <- list()
  for(x in 1:igraph::vcount(graph)){
    l[[x]] <- one_gene_closure(x, edges, functional_forms_list[[x]], param_list[[x]])
  }

  list(
    closures = lapply(l,FUN = function(x) x$closure),
    metadata = lapply(l,FUN = function(x) x$metadata)
    )
}



sim_closure <- function(graph, all_gene_closures){
  sorted_nodes <- as.numeric(topo_sort(graph,mode = "out"))
  
  sim_function <- function(errors){
    if(length(errors) != length(sorted_nodes)){ stop("Need error vector of appropriate length.") }
    vec <- rep(0, length(sorted_nodes))
    for(x in sorted_nodes){
      vec[x] <- all_gene_cl$closures[[x]](vec, errors)
    }
    return(vec)
  }
  
  list(f = sim_function, metadata = all_gene_closures$metadata)
}








# # set.seed(144)
# g <- sample_pa(20, power=1, directed=T)
# 
# 
# params = list(coef_distribution = function(x) runif(x,-1,1))
# params = list(coefs = c(0.5,0.2))
# one_gene_closure(1,get.edgelist(g),functional_form = "linear",params)
# 
# params <- rep(list(params),20)
# a<-all_gene_closures(g,rep(list("linear"),20),param_list = params)
# a
# 
# s<-sim_closure(g,a)
# 
# 
# 
# s$f(rep(1,20))
# 
# 
# 
# 

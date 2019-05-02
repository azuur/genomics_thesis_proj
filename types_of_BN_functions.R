#TYPE1
#additive, with random priority of quadratic, right-truncated log(x+1), or reflected sqrt
create_additive_T1_params <- function(x){
  functional_form <- list("additive")
  params <- list(
    coef_distribution = function(n) runif(n, -1, 1),
  functions_add = sample(
    list(function(x) ifelse(x<0, -x^2, x^2), 
         function(x) ifelse(x<(-0.5), log(0.5), log(x+1)),
         function(x) ifelse(x<0, -(-x)^(1/2), x^(1/2))),
    6, replace = T)
  )
  
  list(functional_form, params)
}

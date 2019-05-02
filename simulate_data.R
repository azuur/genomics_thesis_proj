require(here)
source("pckgs_and_useful_wrappers.R")


#poner opción isnull
#readobj(subdir="Robjects/graph_objects")
SC_20 <- readRDS("~/Documents/Code_thesis/Robjects/graph_objects/SC_20.RDS")


source("BN_simulator.R")
source("types_of_BN_functions.R")

set.seed(10)
sim_info <- lapply(X = 1:20, FUN = create_additive_T1_params)
functional_forms_list <- do.call(c, lapply(X = sim_info, function(l) l[[1]]))
param_list = do.call(c, lapply(X = sim_info, function(l) list(l[[2]])))
all_g <- all_gene_closures(SC_20, 
                              functional_forms_list = functional_forms_list,
                              param_list = param_list)


simulator1 <- sim_closure(SC_20, sim_info)

errores <- runif(20,-1,1)
simulator1$f(errores)



require(here)
source("pckgs_and_useful_wrappers.R")


#poner opción isnull
readobj(subdir="Robjects/graph_objects")


source("BN_simulator.R")
source("types_of_BN_functions.R")

set.seed(10)
sim_info <- lapply(X = 1:200, FUN = create_additive_T1_params)
functional_forms_list <- do.call(c, lapply(X = sim_info, function(l) l[[1]]))
param_list = do.call(c, lapply(X = sim_info, function(l) list(l[[2]])))
all_g <- all_gene_closures(SC_200, 
                              functional_forms_list = functional_forms_list,
                              param_list = param_list)


{simulator1 <- sim_closure(SC_200, all_g)
errores <- runif(200*5000,-1,1) %>% matrix(nrow=200)
sims <- apply(errores, MARGIN = 2, simulator1$f)

system.time(bnlearn::pc.stable(x = as.data.frame(t(errores)[,1:70])))
}
#EN REDES PEQUEÑAS
#Decisión 0. Número de réplicas
#10mil?
#Decision 1. Tamaño de muestra (5 casos)
#n=20,50,100,500,1000
#Decisión 2. Errores (6 casos)
#normales, uniformes. varianza chica, varianza media, varianza grande
datos_SC <- readxl::read_excel(here("SCerevisiae_data/PLOSONE_Their_database.xlsx"))


datos_SC %>% apply(MARGIN = 2, FUN = var) %>% summary()
datos_SC %>% apply(MARGIN = 2, FUN = var) %>% hist()
datos_EC %>% apply(MARGIN = 2, FUN = var) %>% summary()
datos_EC %>% apply(MARGIN = 2, FUN = var) %>% hist()
#Decisión 3. Tipos de relación (3 casos)
#solo lineales, 
#solo relU/sigmoide/otra-monótona, 
#por tercios al azar (lineal, monótona, raraza)


#REDES GRANDES
#Decisión 0. Número de réplicas
#10mil?
#Decision 1. Tamaño de muestra (4 casos)
#n=50,100,500,1000
#Decisión 2. Errores (2 casos)
#normales, uniformes. varianza chica, varianza media, varianza grande
datos_SC %>% apply(MARGIN = 2, FUN = var) %>% summary()
datos_SC %>% apply(MARGIN = 2, FUN = var) %>% hist()
datos_EC %>% apply(MARGIN = 2, FUN = var) %>% summary()
datos_EC %>% apply(MARGIN = 2, FUN = var) %>% hist()
#Decisión 3. Tipos de relación (2 casos)
#solo lineales/relU/sigmoide/otra-monótona, 
#por tercios al azar (monótonas, rarazas)

#GUARDAR
#puntajes
#redes resultantes
## ~ ~ ~ idea... en PC, calcular tasas de falsos positivos bajo indep, y ajustar alpha con eso
#tiempo de


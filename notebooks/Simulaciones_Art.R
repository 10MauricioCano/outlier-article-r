options(warn = -1)
set.seed(50)
## Required Packages
my_packages = c("expm","gridExtra","tidyverse","knitr","kableExtra",
               "IRdisplay","mrfDepth","e1071","pracma","MASS","mixtools",
                "doParallel","foreach","extraDistr","rrcov","Rfast",
                "dobin","FNN","parallel","depthTools","OutliersO3", "mvoutlier",
                "fds", "caret", "tictoc", "REPPlab", "rJava", "PPcovMcd", "ICSOutlier", "doRNG")
not_installed = my_packages[!(my_packages %in% installed.packages()[ , "Package"])]
if (length(not_installed)) install.packages(not_installed, dependencies = TRUE)
for (q in 1:length(my_packages)) {
  library(my_packages[q], character.only = TRUE)
}
## Working Directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
## Load outlier detection codes and auxiliary routines
source("OL_skew_self_V2.R")

############### Originalmente son 50 iteraciones
iter = 50
alpha = c(0.05, 0.1, 0.2, 0.3, 0.4)
dimension = c(10, 30, 50)

############### Originalmente es 10*dimensión
cant_datos = 10*dimension
outl_dist = c(3,4)
outl_concent = c(0.5)
#Funciones High Dimensional Outlier detection----
sd2 <- function(x){ 
  n = length(x)
  sd(x) / sqrt(n) * sqrt(n-1)
}
quantil_id_robust <-function(j, n, d, t){
  muestra = matrix(rnorm(d * n), n, d)
  outlier = rnorm(d)
  pto = t / sqrt(sum(outlier^2)) * outlier
  v = rnorm(d)
  proy = muestra %*% v
  abs(pto %*% v - median(proy)) / max(mad(proy), 1e-6)
}
calculate_ab_robust <- function(d, n, alpha = 0.05) {
  EK = 50 
  p = 1 - 1 / EK
  u = (1 - p) * (1 - alpha)
  v = 1 - alpha / (1 - alpha) * u
  C_nd = sqrt(qchisq(0.95 ^ (1 / n), d)) 
  
  # Simulamos los cuantiles
  num_simulations = 100  
  a0_b0 = sort(unlist(lapply(1:num_simulations, function(j) quantil_id_robust(j, n = n, d = d, t = C_nd))))
  a0 = a0_b0[ floor(length(a0_b0) * u)]
  b0 = a0_b0[ floor(length(a0_b0) * v)]
  
  return(list(a0 = a0, b0 = b0, C_nd = C_nd))
}

final_robust <- function(j, n, d, t = C_nd, rsigma, a, b){
  
  tryCatch(
  
    {VP <- 0
    outlier <-  rnorm(d)
    pto <-  t / sqrt(sum(outlier^2)) * outlier %*% rsigma
    muestra <-  matrix(rnorm(d * n), n, d) %*% rsigma 
    v <-  rnorm(d)
    proy <-  muestra %*% v
    pto_proy <-  abs(pto %*% v - median(proy)) / max(mad(proy), 1e-6)
    
    while(pto_proy > a && pto_proy < b){
      v <-  rnorm(d)
      proy <-  muestra %*% v
      pto_proy <-  abs(pto %*% v - median(proy)) / max(mad(proy), 1e-6)
    }
    if(pto_proy > b){VP <-  1}
    VP},
    error = function(e) {
      message("Ocurrió un error: ", e$message)
      NA  # Retorna NA en caso de error
    },
    warning = function(w) {
      message("Advertencia: ", w$message)
      Inf  # Retorna infinito en caso de advertencia
    },
    finally = {
      message("Finalizando la operación.")
    }
  )
  
}

#########################################################################################################
####### Multivariate normal distribution N(0,I) contaminated with (alpha/2)*N(0 + delta,lambda*I) #######
#########################################################################################################

n_cores <- max(1, parallel::detectCores() - 1)   # deja 1 libre
registerDoParallel(n_cores)
#registerDoParallel(detectCores())

#res4 = foreach (i2 = 1:length(dimension), .combine = rbind) %do% {
  #res3 = foreach (i1 = 1:length(alpha), .combine = rbind) %do% {
    #res2 = foreach (i3 = 1:length(outl_dist), .combine = rbind) %do% {
      #res1 = foreach (i4 = 1:length(outl_concent), .combine = rbind) %do% {
        #res0 = foreach (i5 = 1:iter, .combine = rbind) %dopar% {

n_metrics <- 14            # 7 métodos × (TPR,FPR)
n_meta    <- 5             # modelo, mode, p, alpha, delta

metric_names <- paste0(rep(c("W1","W2","W3","W4","W5","W6","W7"), each = 2),
                       c(".TPR",".FPR"))

colnames_all <- c("modelo","mode","p","alpha","delta", metric_names)
resul.t <- matrix(NA, 0, length(colnames_all))
colnames(resul.t) <- colnames_all

modelos <- c("FDCM", "FICM")      # 2 × 3 = 6 escenarios
sim.modes <- 1:3

#Niveles de clúster ----
safeCM <- function(pred, ref, lev = c("0","1")) {
  pred <- factor(pred, levels = lev)
  ref  <- factor(ref,  levels = lev)
  if (length(pred) != length(ref))
    stop("pred y ref tienen longitudes distintas")
  confusionMatrix(pred, ref)$table
}

#r_i-----
calculate_r_i <- function(x) {
  num_components <- ncol(x)
  r_i <- matrix(NA, nrow = nrow(x), ncol = num_components)
  for (i in 1:num_components) {
    xi <- x[, i]  # Usar una variable temporal para la columna
    mu <- abs(xi - median(xi))
    s <- median(mu)
    r_i[, i] <- mu / s
  }
  return(r_i)
}

#High Dimensional Outlier Detection ----

robust_outlier_detection <- function(X, label.X, a0, b0, C_nd) {
  n <- nrow(X)
  d <- ncol(X)
  
  num_projections <- dimension[i2]  
  VP <- rep(0, n)
  # 
  for (i in 1:n) {
    pto <- X[i, ]
    detected <- FALSE
    for (j in 1:num_projections) {
      v <- rnorm(d)
      proy <- X %*% v
      pto_proy <- abs(pto %*% v - median(proy)) / max(mad(proy), 1e-6)
      if (pto_proy > b0) {
        VP[i] <- 1
        detected <- TRUE
        break
      }
    }
    if (!detected) {
      VP[i] <- 0
    }
  }
  
 
  confusion_mat <- safeCM(VP, label.X)
  return(confusion_mat)
}

# Peña y Prieto Transformation (2001) para identificación de clusters----

pena_prieto_transformation <- function(X, dob_base) {
  dob_base <- calculated_dobin$rotation
  s1 <- cov(X)
  inv1 <- solve(t(dob_base) %*% s1 %*% dob_base)
  I <- diag(nrow(dob_base))
  Q2 <- I - (solve(t(dob_base) %*% s1 %*% dob_base) * dob_base %*% t(dob_base) %*% s1)
  T3 <- X %*% Q2
  return(T3)
}

# Cálculo de Outlyingness Tradicional-----
calculate_outlyingness <- function(selected_data) {
  num_components <- ncol(selected_data)
  at <- matrix(NA, nrow = nrow(selected_data), ncol = num_components)
  for (i in 1:num_components) {
    x <- selected_data[, i]
    mu <- abs(x - median(x))
    s <- median(mu)
    at[,i] <- mu / s
  }
  return(apply(at, 1, max))
}

#Outlyingness ajustado por componente SDC------

SDC <- function(x, r_i){
  
  mediana <- apply(x, 2, median)
  p <- ncol(x)
  n <- nrow(x)
  ind_1 <- floor((n + p - 1) / 2)
  ind_2 <- ceiling((n + p - 1) / 2) + 1
  beta <- qnorm(0.5 * ((n + p - 1) / (2 * n) + 1))
  
  mad_modificado <- numeric(length = p)
  
  # Calcular la MAD modificada para cada componente
  for (i in 1:p) {
    desv_abs <- abs(x[, i] - mediana[i])
    desv_abs_ordenado <- sort(desv_abs)
    mad_modificado[i] <- (desv_abs_ordenado[ind_1] + desv_abs_ordenado[ind_2]) / (2 * beta)
  }
  
  # Outlyingness de observación i en dirección j c_ij
  c_ij <- matrix(NA, nrow = n, ncol = p)
  for (i in 1:p) {
    c_ij[, i] <- abs(x[, i] - mediana[i]) / mad_modificado[i]
  }
  
  # Calcular el máximo c_ij para cada observación
  max_c_ij <- apply(c_ij, 1, max)
  
  # Evitar divisiones por cero
  max_c_ij[max_c_ij == 0] <- 1e-6
  
  # Calcular alpha_SDC
  
  alpha_SDC <- sweep(c_ij, 1, max_c_ij, "/")
  
  
  # Calcular r_ij
  r_ij <- alpha_SDC * r_i + (1 - alpha_SDC) * c_ij
  
  # Calcular el outlyingness ajustado
  outlyingness_SDC <- apply(r_ij, 1, max)
  
  return(outlyingness_SDC)
}

#Cálculo del outlyingness por componentes maximizados----
SDM <- function(x2){
  mediana <- apply(x2, 2, median)
  p <- ncol(x2)
  n <- nrow(x2)
  ind_1 <- floor((n + p - 1) / 2)
  ind_2 <- ceiling((n + p - 1) / 2) + 1
  beta <- qnorm(0.5 * ((n + p - 1) / (2 * n) + 1))
  mad_componentes_seleccionados <- apply(x2, 2, mad)
  reescalado <- sweep(x2,1,mad_componentes_seleccionados, FUN = '/')
  #Se calcula el respectivo r_i
  mu_reescalado <- abs(reescalado-median(reescalado))
  s_reescalado <- median(mu_reescalado)
  r_i_SDM <- mu_reescalado/s_reescalado
  
  
  #Cálculo del alpha para el método SDM----
  u_i <- apply(r_i_SDM, 1, max)
  alpha_SDM <- matrix(NA, nrow = nrow(x2), ncol = ncol(x2))
  alpha_SDM <- r_i_SDM / u_i
  
  # Calcular la MAD modificada para cada componente
  mad_modificado <- numeric(length = p)
  for (i in 1:ncol(x2)) {
    desv_abs <- abs(x2[, i] - mediana[i])
    desv_abs_ordenado <- sort(desv_abs)
    mad_modificado[i] <- (desv_abs_ordenado[ind_1] + desv_abs_ordenado[ind_2]) / (2 * beta)
  }
  
  # Outlyingness de observación i en dirección j c_ij
  c_ij <- matrix(NA,nrow = nrow(x2), ncol = ncol(x2))
  
  for (i in 1:ncol(x2)) {
    c_ij[,] <- abs(x2[,]-mediana[i])/mad_modificado[i]
  }
  
  # Outlyingness adaptado de observación i en dirección que maximiza outlyngness
  
  r_ij_SDM <- matrix(NA, nrow = nrow(x2), ncol= ncol(x2))
  
  r_ij_SDM[,] <- alpha_SDM[,] * r_i[,] + (1 - alpha_SDM) * c_ij[,]
  
  outlyingness_SDM <- apply((r_ij_SDM), 1, max)
  
  return(outlyingness_SDM)
}
#Ejecución de escenarios ----
for (modelo in modelos) {
  for (mode in sim.modes) {
    for (i2 in seq_along(dimension)) {
      for (i1 in seq_along(alpha)) {
        for (i3 in seq_along(outl_dist)) {
          for (i4 in seq_along(outl_concent)) {
            
            resul_mat <- foreach(i5 = 1:iter,
                             .combine   = rbind,      # apila filas
                             .packages  = my_packages, # ← reutiliza vector
                             .export      = c("safeCM","calculate_r_i",
                                              "SDC","SDM",
                                              "pena_prieto_transformation",
                                              "calculate_outlyingness"),
                             .options.RNG = 50) %dorng% {
                               
              
              if (modelo == "FDCM") {
                data.sim <- GenAtip(cant_datos[i2], dimension[i2],
                                    c(alpha[i1], outl_dist[i3],
                                      outl_concent[i4]),
                                    sim.mode = mode)
              } else {           
                data.sim <- GenAtip_FICM(cant_datos[i2], dimension[i2],
                                         c(alpha[i1], outl_dist[i3],
                                           outl_concent[i4]),
                                         sim.mode = mode)
              }
              
              X       <- data.sim$x
              label.X <- factor(data.sim$lbl)
          #-----------------------#-----------------------#-----------------------#-----------------------#-----------------------
         
          
          #-----------------------#-----------------------#-----------------------#---------------------------------------------
          ##
          #pairs(X[,1:5], pch = 20)
          ##
          calculated_dobin<-dobin(X)
          dob_knn <- X%*%calculated_dobin$rotation
          
          
          #Exploratory Projection Pursuit REPPLAB----
          
          X_EPP <- EPPlab(X, PPalg="PSO", PPindex="FriedmanTukey", n.simu=1, maxiter=1000, sphere=TRUE)
          
          out <- EPPlabOutlier(X_EPP, which = 1:ncol(X), k = 3, location = mean,
                               scale = sd)
          aux_repp <- unlist(out$outlier)
          threshold_repp <- quantile(aux_repp, 0.95)
          
          # Generar etiquetas predichas basadas en el umbral
          predicted_labels_repp <- ifelse(aux_repp > threshold_repp, "1", "0")
          
          
          # Alinear niveles de factores
          #levels_to_use <- union(levels(as.factor(predicted_labels_repp)), levels(as.factor(label.X)))
          #predicted_labels_repp <- factor(predicted_labels_repp, levels = levels_to_use)
          #label.X <- factor(label.X, levels = levels_to_use)
          
          
          
          
          
          # Calcular la matriz de confusión
          
          W1.cm <- safeCM(predicted_labels_repp, label.X)
          
          
          
          
          datos_limp <- sum(data.sim$lbl == 0)   # 0 = fila sin ninguna celda contaminada
          
          # Evita problemas si no queda ninguna fila limpia (puede ocurrir en FICM)
          if (datos_limp < 2) {
            datos_limp <- 2         # valor mínimo para que mad/quantiles no fallen
          }
          
          ab_values <- calculate_ab_robust(d = dimension[i2], n = datos_limp, alpha = 0.05)
          a0 <- ab_values$a0
          b0 <- ab_values$b0
          C_nd <- ab_values$C_nd
          W2.cm <- robust_outlier_detection(X, label.X, a0, b0, C_nd)
          
          
          # Projection Pursuit Minimum Covariance Determinant----
          out_values <- PPcovMcd(X)
          Projection_pursuit_MCD <- mahalanobis(X, center = out_values$center, cov = out_values$cov)
          threshold_MCD <- quantile(Projection_pursuit_MCD, 0.95)
          predicted_labels_MCD <- ifelse(Projection_pursuit_MCD > threshold_MCD, "1", "0")

          levels_to_use <- union(levels(as.factor(predicted_labels_MCD)), levels(as.factor(label.X)))
          predicted_labels_MCD <- factor(predicted_labels_MCD, levels = levels_to_use)
          label.X <- factor(label.X, levels = levels_to_use)

          W3.cm <- safeCM(predicted_labels_MCD, label.X)

          # Invariant Coordinate Selection ICS----
          data_ics <- ics2(X)
          
          # Ajustar el nivel de significancia para aumentar la sensibilidad
          resultado_ICS <- ics.outlier(data_ics, level.test = 0.1)
          
          # Acceder al slot 'outliers'
          predicted_labels <- resultado_ICS@outliers
          
          # Verificar la longitud de predicted_labels
          
          # Convertir predicted_labels y label.X a factores con los mismos niveles
          #predicted_labels <- factor(predicted_labels, levels = c(0, 1))
          #label.X <- factor(as.numeric(as.character(label.X)), levels = c(0, 1))
          
          # Calcular la matriz de confusión
          W4.cm <- safeCM(predicted_labels, label.X)
          
# DOBIN with S-Orthogonal Projection----
          
          
         # Selección de componentes basados en Kurtosis y Skewness
            component_selection <- function(transformation) {
            kurtosis_per_component <- apply(transformation, 2, kurtosis)
            skewness_per_component <- apply(transformation, 2, skewness)
            sum_coef_per_component <- kurtosis_per_component^2 + skewness_per_component^2
            
            ordered_indexes <- order(-sum_coef_per_component)
            cutoff <- median(sum_coef_per_component[ordered_indexes])
            selected_components <- which(sum_coef_per_component[ordered_indexes] > cutoff)
            
            return(transformation[, selected_components])
          }
          
          #-----------------------#-----------------------#-----------------------#-----------------------#----------------------------------------
          
          
          T3 <- pena_prieto_transformation(X, dob_base)
          selected_data <- component_selection(T3)
          at_m <- calculate_outlyingness(selected_data)
          
# Huberizar los datos ----
          
          cH <- apply(abs(selected_data), 2, quantile, probs = 0.975, na.rm = TRUE)
          
          mediana <- apply(selected_data, 2, median, na.rm= TRUE)
          mad_orig <- apply(selected_data,2,mad, na.rm=TRUE)
          c_ij_huber <- matrix(NA, nrow = nrow(selected_data), ncol = ncol(selected_data))
          
          for (i in 1:ncol(selected_data)){
            c_ij_huber[,i] <- (selected_data[,i] - mediana[i]) / mad_orig[i]
          }
          
          sel_Hub <- matrix(0,nrow = nrow(selected_data), ncol = ncol(selected_data))
          
          for (i in 1:ncol(selected_data)) {
            sel_Hub[, i] <- ifelse(c_ij_huber[, i] < -cH[i], mediana[i] - cH[i] * mad_orig[i], 
                                 ifelse(c_ij_huber[, i] > cH[i], mediana[i] + cH[i] * mad_orig[i], X[, i]))
          }
          
          selected_data <- sel_Hub
          
          # Cálculo de r_i-----
          r_i <- calculate_r_i(selected_data)

          
          #-----------------------#-----------------------#-----------------------#-----------------------#-------------------------------------
          
          # Prueba del método
          
          # Con cálculo de outlyingness tradicional
          T3 <- pena_prieto_transformation(X, dob_base)
          selected_data <- component_selection(T3)
          at_m <- calculate_outlyingness(selected_data)
          
          cl <- kmeans(at_m, centers=2)
          if(sum(cl$cluster == 1) < sum(cl$cluster == 2)) {
            cl$cluster <- ifelse(cl$cluster == 1, 2, 1)
          }
          vec3 <- as.factor(data.sim$lbl)
          
          cl$cluster <- cl$cluster - 1
          W5.cm <- safeCM(cl$cluster, vec3)
          #W5.cm <- confusionMatrix(as.factor(cl$cluster), vec3)$table
          #print(W5.cm)
          
          # Recalcula r_i con el nuevo selected_data
          r_i <- calculate_r_i(selected_data)
          
          # Con Cálculo de Outlyingness ajustado por componente
          outlyingness_ajustado_SDC <- SDC(selected_data, r_i)
          
          # Clustering usando el outlyingness ajustado
          cl <- kmeans(outlyingness_ajustado_SDC, centers = 2)
          if (sum(cl$cluster == 1) < sum(cl$cluster == 2)) {
            cl$cluster <- ifelse(cl$cluster == 1, 2, 1)
          }
          vec4 <- as.factor(data.sim$lbl)
          
          cl$cluster <- cl$cluster - 1
          W6.cm <- safeCM(cl$cluster, vec4)
          
          #W6.cm <- confusionMatrix(as.factor(cl$cluster), vec4)$table
          #print(W6.cm)
          
          # Cálculo del Outlyingness por maximización SDM
          
          outlyingness_ajustado_SDM <- SDM(selected_data)
          
          # Clustering usando el outlyingness maximizado
          cl <- kmeans(outlyingness_ajustado_SDM, centers = 2)
          if (sum(cl$cluster == 1) < sum(cl$cluster == 2)) {
            cl$cluster <- ifelse(cl$cluster == 1, 2, 1)
          }
          vec5 <- as.factor(data.sim$lbl)
          
          cl$cluster <- cl$cluster - 1
          W7.cm <- safeCM(cl$cluster, vec5)
          
          #W7.cm <- confusionMatrix(as.factor(cl$cluster), vec5)$table
          #print(W7.cm)
          
          
          #-----------------------#-----------------------#-----------------------#-------------------------------
          
          # RESULTS
                      c(W1.cm[2,2]/sum(W1.cm[,2]), W1.cm[2,1]/sum(W1.cm[,1]),
                          W2.cm[2,2]/sum(W2.cm[,2]), W2.cm[2,1]/sum(W2.cm[,1]),
                          W3.cm[2,2]/sum(W3.cm[,2]), W3.cm[2,1]/sum(W3.cm[,1]),
                          W4.cm[2,2]/sum(W4.cm[,2]), W4.cm[2,1]/sum(W4.cm[,1]),
                          W5.cm[2,2]/sum(W5.cm[,2]), W5.cm[2,1]/sum(W5.cm[,1]),
                          W6.cm[2,2]/sum(W6.cm[,2]), W6.cm[2,1]/sum(W6.cm[,1]),
                          W7.cm[2,2]/sum(W7.cm[,2]), W7.cm[2,1]/sum(W7.cm[,1])
                        )
          
          
          
            }
            metric_means <- colMeans(resul_mat, na.rm = TRUE)           # numérico
            resul1 <- c(modelo, mode, dimension[i2], alpha[i1],
                        outl_dist[i3], metric_means)                # 5 + 14 = 19
            resul.t <- rbind(resul.t, resul1) 
      }}}}}}

#print("Experiment Symmetric TypeA has finished")
save(resul.t, file = "Symm.RData")
#save(resul.t,file = 'datos.xls')

write.csv(resul.t, "resul_t.csv", row.names=FALSE)
write.csv(resul, "resul.csv", row.names=FALSE)


stopImplicitCluster()

###########
### FIN ###
###########


library(here)

# Always resolve paths from the repo root (RStudio Project)

# Outputs directory inside the repo

data_dir <- here::here("data", "raw")
stopifnot(dir.exists(data_dir))


run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- here::here("outputs", paste0("run_", run_id))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


## Required Packages----


my_packages = c("expm","gridExtra","tidyverse","knitr","kableExtra",
                "IRdisplay","mrfDepth","e1071","pracma","MASS","mixtools",
                "doParallel","foreach","extraDistr","rrcov","Rfast",
                "dobin","FNN","parallel","depthTools", "mvoutlier",
                "fds", "caret", "tictoc", "REPPlab", "PPcovMcd", "ICSOutlier", "data.table"
                , "R.matlab")

missing <- my_packages[!vapply(my_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing packages: ", paste(missing, collapse = ", "),
       ". Install them (in renv) with install.packages(...).")
}
invisible(lapply(my_packages, library, character.only = TRUE))

# Función auxiliar para convertir datos a matriz----
convert_to_matrix <- function(data) {
  if (is.data.frame(data)) {
    return(as.matrix(data))
  } else if (is.matrix(data)) {
    return(data)
  } else if (is.list(data)) {
    return(do.call(cbind, data))
  } else {
    stop("Unsupported data type")
  }
}


# Función para crear matrices de confusión por cada simulación----

cm_to_row <- function(cm, sim_id, dataset, method, positive = "2", negative = "1") {
  if (is.null(dim(cm)) || any(dim(cm) != c(2,2))) {
    cm <- matrix(NA_integer_, 2, 2, dimnames = list(c(negative, positive), c(negative, positive)))
  }
  if (is.null(rownames(cm))) rownames(cm) <- c(negative, positive)
  if (is.null(colnames(cm))) colnames(cm) <- c(negative, positive)
  
  data.frame(
    simulation = sim_id, dataset = dataset, method = method,
    TN = cm[negative, negative],
    FP = cm[positive, negative],
    FN = cm[negative, positive],
    TP = cm[positive, positive],
    stringsAsFactors = FALSE
  )
}




# Cálculo de métricas de desempeño----

calculate_metrics <- function(conf_matrix) {
  TN <- conf_matrix["1","1"]; FP <- conf_matrix["2","1"]
  FN <- conf_matrix["1","2"]; TP <- conf_matrix["2","2"]
  
  den_sens <- TP + FN
  den_spec <- TN + FP
  den_acc  <- TN + FP + FN + TP
  
  sensitivity <- ifelse(den_sens > 0, TP/den_sens, NA_real_)
  specificity <- ifelse(den_spec > 0, TN/den_spec, NA_real_)
  accuracy    <- ifelse(den_acc  > 0, (TP+TN)/den_acc, NA_real_)
  balanced_accuracy <- ifelse(is.na(sensitivity) | is.na(specificity),
                              NA_real_, (sensitivity + specificity)/2)
  
  c(Sensitivity = sensitivity,
    Specificity = specificity,
    Accuracy = accuracy,
    Balanced_Accuracy = balanced_accuracy)
}

rm(list = ls(pattern = "^predicted_labels_"), inherits = TRUE)

# Simulaciones----

all_cm <- list()
all_metrics <- list()
all_raw <- list()
set.seed(123)

simulations_to_run <- 2:4
#Borrar resultados de simulaciones anteriores----


for(simulation_proposed in simulations_to_run) {

rm(list = ls(pattern = "^predicted_labels_"), inherits = TRUE)
  
  # Simulation 2 ----
  if (simulation_proposed == 2) {
    file_path <- file.path(data_dir, "wine.mat")
    mat_data <- readMat(file_path)
    x3 <- convert_to_matrix(mat_data$X)
    fil <- nrow(x3) / 2
    vec <- as.factor(mat_data$y + 1)
    true_labels <- factor(vec, levels = c(1, 2))
    dataset_name <- "Wine"
    
    # Simulation 3 ----
  } else if (simulation_proposed == 3) {
    file_path <- file.path(data_dir, "satimage-2.mat")
    mat_data <- readMat(file_path)
    x3 <- mat_data$X
    fil <- nrow(x3) / 2
    vec <- as.factor(mat_data$y + 1)
    true_labels <- factor(vec, levels = c(1, 2))
    dataset_name <- "Satimage-2"
    
    # Simulation 4 ----
  } else if (simulation_proposed == 4) {
    file_path <- file.path(data_dir, "breastw.mat")
    mat_data <- readMat(file_path)
    x3 <- mat_data$X
    fil <- nrow(x3) / 2
    vec <- as.factor(mat_data$y + 1)
    true_labels <- factor(vec, levels = c(1, 2))
    dataset_name <- "Breastw"
    
  } else {
    stop("Simulations go from 2 to 4.")
  }
  

x3 <- convert_to_matrix(x3)


# Conjunto de datos no huberizados ----
X <- x3  
d <- ncol(X)
n <- nrow(X)

# Huberización de datos: se calcula valor de corte cH (percentil 97.5)-----

huber <- function(x){
  x <- as.matrix(x)
  
  mediana_huber <- apply(x, 2, median, na.rm = TRUE)
  mad_huber <- pmax(1e-6,apply(x, 2, mad))
  
  cH <- numeric(ncol(x))
  c_ij_huber <- matrix(0, nrow = nrow(x), ncol = ncol(x))
  x_hub <- matrix(0, ncol = ncol(x), nrow = nrow(x))
  
  for (i in seq_len(ncol(x))){
    c_ij_huber[, i] <- (x[, i] - mediana_huber[i]) / mad_huber[i]
    cH[i] <- quantile(abs(c_ij_huber[, i]), probs = 0.975, na.rm = TRUE)
    c_aux <- pmax(-cH[i], pmin(cH[i], c_ij_huber[, i]))
    x_hub[, i] <- mediana_huber[i] + c_aux * mad_huber[i]
   
  }
  return(x_hub)
}

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
  abs(pto %*% v - median(proy)) / mad(proy)
}




calculate_ab_robust <- function(d, n, alpha = 0.05) {
  EK = 50  # Puedes ajustar este valor
  p = 1 - 1 / EK
  u = (1 - p) * (1 - alpha)
  v = 1 - alpha / (1 - alpha) * u
  C_nd = sqrt(qchisq(0.95 ^ (1 / n), d)) 
  
  # Simulamos los cuantiles
  num_simulations = 100  # Ajusta este número según tus recursos computacionales
  a0_b0 = sort(unlist(lapply(1:num_simulations, function(j) quantil_id_robust(j, n = n, d = d, t = C_nd))))
  a0 = a0_b0[ floor(length(a0_b0) * u)]
  b0 = a0_b0[ floor(length(a0_b0) * v)]
  
  return(list(a0 = a0, b0 = b0, C_nd = C_nd))
}

final_robust <- function(j, n, d, t = C_nd, rsigma, a, b){
  VP = 0
  outlier = rnorm(d)
  pto = t / sqrt(sum(outlier^2)) * outlier %*% rsigma
  muestra = matrix(rnorm(d * n), n, d) %*% rsigma 
  v = rnorm(d)
  proy = muestra %*% v
  pto_proy = abs(pto %*% v - median(proy)) / mad(proy)
  while(pto_proy > a && pto_proy < b){
    v = rnorm(d)
    proy = muestra %*% v
    pto_proy = abs(pto %*% v - median(proy)) / mad(proy)
  }
  if(pto_proy > b){VP = 1}
  VP
}


## ------------------------------------------------------------------------------
## 2. APLICACIÓN DE MÉTODOS (epp, hdod, pp-mcd, ics)----
## ------------------------------------------------------------------------------



# if (simulation_proposed == 1){
#   vec<- c(rep(1,sample1),rep(2,sample2))
#   true_labels <- factor(vec,         levels = c(1,2)) 
# } else if (simulation_proposed == 2){
#   fil<-nrow(x3)/2
#   vec <- as.factor(mat_data$y + 1)
#   true_labels <- factor(vec,         levels = c(1,2))  
# } else if(simulation_proposed == 3){
#   vec <- as.factor(mat_data$y + 1)
#   true_labels <- factor(vec,         levels = c(1,2)) 
# } else if (simulation_proposed == 4) {
#   fil <- nrow(x3)/2
#   vec <- as.factor(mat_data$y + 1)
#   true_labels <- factor(vec,         levels = c(1,2)) 
# } 

#### 2.1. epp (Exploratory Projection Pursuit)----

X_EPP <- EPPlab(X, PPalg="PSO", PPindex="FriedmanTukey", n.simu=1, maxiter=1000, sphere=TRUE)

out <- EPPlabOutlier(X_EPP, which = 1:ncol(X), k = 3, location = mean,
                     scale = sd)
aux_repp <- unlist(out$outlier)
threshold_repp <- quantile(aux_repp, 0.95)

# Generar etiquetas predichas basadas en el umbral
predicted_labels_repp <- ifelse(aux_repp > threshold_repp, "2", "1")
# Alinear niveles de factores
#levels_to_use <- union(levels(as.factor(predicted_labels_repp)), levels(as.factor(true_labels)))
#predicted_labels_repp <- factor(predicted_labels_repp, levels = levels_to_use)
#true_labels <- factor(true_labels, levels = levels_to_use)

predicted_labels_repp <- factor(predicted_labels_repp, levels = c("1","2"))
true_labels <- factor(true_labels, levels = c("1","2"))

stopifnot(length(predicted_labels_repp) == length(true_labels))
# Calcular la matriz de confusión
W1.cm <- confusionMatrix(predicted_labels_repp, true_labels)$table

#### 2.2 hdod (High Dimensional Outlier Detection)----

robust_outlier_detection <- function(X, true_labels, a0, b0, C_nd, num_projections = 50) {
  d <- ncol(X)
  n <- nrow(X)
  
  
  # Creamos un vector (0=inlier / 1=outlier)
  VP <- rep(0, n)
  
  for (i in seq_len(n)) {
    pto <- X[i, ]
    for (j in seq_len(num_projections)) {
      v <- rnorm(d)
      proy <- X %*% v
      pto_proy <- abs(pto %*% v - median(proy)) / mad(proy)
      if (pto_proy > b0) {
        VP[i] <- 1
        break
      }
    }
  }
  
  predicted_labels_hdod <- ifelse(VP == 1, "2", "1")
  
  # Mismos niveles que true_labels (1=inlier,2=outlier)
  predicted_labels_hdod <- factor(predicted_labels_hdod, levels = c("1", "2"))
  true_labels_factor     <- factor(true_labels,         levels = c("1", "2"))
  
  confusion_mat <- confusionMatrix(predicted_labels_hdod, true_labels_factor)$table
  return(confusion_mat)
}



ab_values <- calculate_ab_robust(
  d = ncol(X),
  n = nrow(X),
  alpha = 0.05
)
a0   <- ab_values$a0
b0   <- ab_values$b0
C_nd <- ab_values$C_nd

W2.cm <- robust_outlier_detection(X, true_labels, a0, b0, C_nd, num_projections = 50)


#### 2.3. pp-mcd (Projection Pursuit + Minimum Covariance Determinant)-----

# tryCatch para que, si PPcovMcd da error, no detenga el flujo
maybe_W3 <- tryCatch(
  {
    out_values <- PPcovMcd(X)
    d2 <- mahalanobis(X, center = out_values$center, cov = out_values$cov)
    thr <- quantile(d2, 0.95, na.rm = TRUE)
    pred <- factor(ifelse(d2 > thr, "2", "1"), levels = c("1","2"))
    truth <- factor(true_labels, levels = c("1","2"))
    list(cm = confusionMatrix(pred, truth)$table, status = "OK")
  },
  error = function(e) {
    cat("PP-MCD falló:", conditionMessage(e), "\n")
    dummy <- matrix(NA_integer_, 2,2, dimnames=list(c("1","2"), c("1","2")))
    list(cm = dummy, status = "FAIL")
  }
)

W3.cm <- maybe_W3$cm
W3.status <- maybe_W3$status


#### 2.4. ics (Invariant Coordinate Selection)----
data_ics <- ics2(X)

# Ajustar el nivel de significancia para aumentar la sensibilidad
resultado_ICS <- ics.outlier(data_ics, level.test = 0.3)


idx <- resultado_ICS@outliers
predicted_labels_ics <- rep("1", nrow(X))
if (length(idx)) predicted_labels_ics[idx] <- "2"
predicted_labels_ics <- factor(predicted_labels_ics, levels=c("1","2"))
W4.cm <- confusionMatrix(predicted_labels_ics, true_labels)$table

## ------------------------------------------------------------------------------
## 3. APLICACIÓN DE DOBIN + SDC / SDM----
## ------------------------------------------------------------------------------


## 3.1. DOBIN base (S-Ort)
dob <- dobin(X)
dob_base <- if (!is.null(dob$basis)) dob$basis else dob$rotation

# DOBIN with S-Orthogonal Projection
# Peña y Prieto Transformation (2001) para identificación de clusters
pena_prieto_transformation <- function(X, dob_base) {
  #dob_base <- calculated_dobin$rotation
  s1 <- cov(X)
  inv1 <- solve(t(dob_base) %*% s1 %*% dob_base)
  I <- diag(nrow(dob_base))
  Q2 <- I - (inv1 %*% dob_base %*% t(dob_base) %*% s1)
  T3 <- X %*% Q2
  return(T3)
}

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

# Cálculo de Outlyingness Tradicional
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
#-----------------------#-----------------------#-----------------------#-----------------------#----------------------------------------
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

T3 <- pena_prieto_transformation(X, dob_base)
selected_data <- component_selection(T3)
at_m <- calculate_outlyingness(selected_data)

# Cálculo de r_i-----
r_i <- calculate_r_i(component_selection(T3))

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
    mad_modificado[i] <- pmax(1e-6,(desv_abs_ordenado[ind_1] + desv_abs_ordenado[ind_2]) / (2 * beta))
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

SDM <- function(x2){
  mediana <- apply(x2, 2, median)
  p <- ncol(x2); n <- nrow(x2)
  ind_1 <- floor((n + p - 1) / 2)
  ind_2 <- ceiling((n + p - 1) / 2) + 1
  beta <- qnorm(0.5 * ((n + p - 1) / (2 * n) + 1))
  mad_componentes_seleccionados <- apply(x2, 2, mad)
  reescalado <- sweep(x2,2,mad_componentes_seleccionados, FUN = '/')
  #Se calcula el respectivo r_i
  mu_reescalado <- abs(sweep(reescalado, 2, apply(reescalado, 2, median), "-"))
  s_reescalado <- pmax(1e-6, apply(mu_reescalado, 2, median))
  r_i_SDM <- sweep(mu_reescalado, 2, s_reescalado, "/")
  
  
  #Cálculo del alpha para el método SDM----
  u_i <- pmax(1e-6, apply(r_i_SDM, 1, max))
  alpha_SDM <- matrix(NA, nrow = nrow(x2), ncol = ncol(x2))
  alpha_SDM <- sweep(r_i_SDM, 1, u_i, "/")
  
  # Calcular la MAD modificada para cada componente
  mad_modificado <- numeric(length = p)
  for (i in 1:ncol(x2)) {
    desv_abs <- abs(x2[, i] - mediana[i])
    desv_abs_ordenado <- sort(desv_abs)
    mad_modificado[i] <- (desv_abs_ordenado[ind_1] + desv_abs_ordenado[ind_2]) / (2 * beta)
  }
  
  # Outlyingness de observación i en dirección j c_ij
  c_ij <- matrix(NA_real_, n, p)
  
  for (i in 1:ncol(x2)) {
    c_ij[, i] <- abs(x2[, i] - mediana[i]) / mad_modificado[i]
  }
  
  
  # Outlyingness adaptado de observación i en dirección que maximiza outlyngness
  
  r_ij_SDM <- matrix(NA, nrow = nrow(x2), ncol= ncol(x2))
  
  r_ij_SDM[,] <- alpha_SDM[,] * r_i_SDM[,] + (1 - alpha_SDM) * c_ij[,]
  
  outlyingness_SDM <- apply((r_ij_SDM), 1, max)
  
  return(outlyingness_SDM)
}

# Cálculo del outlyingness ajustado por componente SDC
# outlyingness_ajustado_SDC <- SDC(selected_data, r_i)


#-----------------------#-----------------------#-----------------------#-----------------------#-------------------------------------

# Prueba del método

# Con cálculo de outlyingness tradicional

T3 <- pena_prieto_transformation(X, dob_base)
selected_data <- component_selection(T3)
at_m <- calculate_outlyingness(selected_data)


#Etiquetado de outliers----

cluster_to_labels <- function(score_vec, alpha = 0.02, alpha_max = 0.05) {
  alpha <- min(alpha, alpha_max)
  n <- length(score_vec)
  k <- max(1L, ceiling(alpha * n))
  ord <- order(score_vec, decreasing = TRUE)
  pred <- rep("1", n)
  pred[ord[seq_len(k)]] <- "2"
  factor(pred, levels = c("1","2"))
}


pred_dso    <- cluster_to_labels(at_m)
W5.cm <- confusionMatrix(pred_dso, true_labels)$table
#W5.cm <- confusionMatrix(as.factor(cl$cluster), as.factor(vec))$table
#print(W5.cm)


# Recalcula r_i con el nuevo selected_data
selected_data_h <- huber(selected_data)
r_i <- calculate_r_i(selected_data_h)

# Con Cálculo de Outlyingness ajustado por componente
outlyingness_ajustado_SDC <- SDC(selected_data_h, r_i)

#-------------------------------------------------------#

# Clustering usando el outlyingness ajustado---
# cl <- kmeans(outlyingness_ajustado_SDC, centers = 2)
# if (sum(cl$cluster == 1) < sum(cl$cluster == 2)) {
#   cl$cluster <- ifelse(cl$cluster == 1, 2, 1)
# }
# 
# W6.cm <- confusionMatrix(as.factor(cl$cluster), as.factor(vec))$table
# print(W6.cm)

pred_sdc    <- cluster_to_labels(outlyingness_ajustado_SDC)
W6.cm <- confusionMatrix(pred_sdc, true_labels)$table

# Cálculo del Outlyingness por maximización SDM----

outlyingness_ajustado_SDM <- SDM(selected_data_h)

# Clustering usando el outlyingness maximizado
# cl <- kmeans(outlyingness_ajustado_SDM, centers = 2)
# if (sum(cl$cluster == 1) < sum(cl$cluster == 2)) {
#   cl$cluster <- ifelse(cl$cluster == 1, 2, 1)
# }
# 
# W7.cm <- confusionMatrix(as.factor(cl$cluster), as.factor(vec))$table
# print(W7.cm)

pred_sdm    <- cluster_to_labels(outlyingness_ajustado_SDM)
W7.cm <- confusionMatrix(pred_sdm, true_labels)$table

#Guardar cada simulación sin sobreescribir ----

sim_id <- simulation_proposed
key <- paste0("sim_", sim_id, "_", dataset_name)

#Consolidar matrices de confusión en un dataframe

cm_block <- rbind(
  cm_to_row(W1.cm, sim_id, dataset_name, "REPPLAB"),
  cm_to_row(W2.cm, sim_id, dataset_name, "HDOD"),
  cm_to_row(W3.cm, sim_id, dataset_name, "PP-MCD"),
  cm_to_row(W4.cm, sim_id, dataset_name, "ICS"),
  cm_to_row(W5.cm, sim_id, dataset_name, "DSO"),
  cm_to_row(W6.cm, sim_id, dataset_name, "DSO-SDC"),
  cm_to_row(W7.cm, sim_id, dataset_name, "DSO-SDM")
)


cms_list <- list(W1.cm, W2.cm, W3.cm, W4.cm, W5.cm, W6.cm, W7.cm)
m <- do.call(rbind, lapply(cms_list, calculate_metrics))

metrics_block <- cbind(
  cm_block[, c("simulation","dataset","method")],
  as.data.frame(m, stringsAsFactors = FALSE)
)

metrics_block[] <- lapply(metrics_block, function(x) {
  if (is.factor(x)) as.character(x) else x
})


#Guardar resultados en listas para que no se pierda nada

all_cm[[key]] <- cm_block
all_metrics[[key]] <- metrics_block

# Escribe archivo por simulación, por trazabilidad
write.csv(cm_block, file.path(out_dir, paste0("cm_", key, ".csv")), row.names = FALSE)
#write.csv(metrics_block, file.path(out_dir, paste0("metrics_", key, ".csv")), row.names = FALSE)


# Calcula métricas a partir de cada matriz de confusión
    metrics_W1 <- calculate_metrics(W1.cm)
    
    metrics_W2 <- calculate_metrics(W2.cm)
    
    if(any(is.na(W3.cm))) {
      
      metrics_W3 <- list(Sensitivity=NA, Specificity=NA, Accuracy=NA)
    } else {
      metrics_W3 <- calculate_metrics(W3.cm)
    }
    
    metrics_W4 <- calculate_metrics(W4.cm)
    
    metrics_W5 <- calculate_metrics(W5.cm)
    
    metrics_W6 <- calculate_metrics(W6.cm)
    
    metrics_W7 <- calculate_metrics(W7.cm)

# Imprimir resultados----
   #cat('Confusion Matrix REPPlab'); print(W1.cm)
   #cat("Metrics for REPPlab:\n"); print(metrics_W1)
   #
   #cat('Confusion Matrix HDOD'); print(W2.cm)
   #cat("Metrics for HDOD:\n"); print(metrics_W2)
   #
   #cat('Confusion Matrix PP-MCD'); print(W3.cm)
   #cat("Metrics for PP-MCD:\n"); print(metrics_W3)
   #
   #cat('Confusion Matrix ICS'); print(W4.cm)
   #cat("Metrics for ICS:\n"); print(metrics_W4)
   #
   #cat('Confusion Matrix DSO'); print(W5.cm)
   #cat("Metrics for DSO:\n"); print(metrics_W5)
   #
   #cat('Confusion Matrix DSO-SDC'); print(W6.cm)
   #cat("Metrics for DSO-SDC:\n"); print(metrics_W6)
   #
   #cat('Confusion Matrix DSO-SDM'); print(W7.cm)
   #cat("Metrics for DSO-SDM:\n"); print(metrics_W7)

}

# Consolidar y exportar resultados----
cm_all <- do.call(rbind, all_cm)
metrics_all <- do.call(rbind, all_metrics)

write.csv(cm_all, file.path(out_dir, "confusion_ALL.csv"), row.names = FALSE)
write.csv(metrics_all, file.path(out_dir, "metrics_ALL.csv"), row.names = FALSE)
#write.csv(metrics_block, file.path(out_dir, paste0("metrics_", key, ".csv")),
          #row.names = FALSE)

saveRDS(list(cm = all_cm, metrics = all_metrics),
        file.path(out_dir, "all_results.rds"))


cat("\nOK. Resultados en: ", out_dir, "\n")

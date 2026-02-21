options(warn = -1)
set.seed(50)

library(here)
source(here::here("R", "utils_analysis.R"))

## Required packages (NO auto-install) ----
my_packages <- c(
  "expm","gridExtra","tidyverse","knitr","kableExtra","IRdisplay","mrfDepth","e1071","pracma",
  "MASS","mixtools","doParallel","foreach","extraDistr","rrcov","Rfast","dobin","FNN","parallel",
  "depthTools","OutliersO3","mvoutlier","fds","caret","tictoc","REPPlab","rJava","PPcovMcd",
  "ICSOutlier","doRNG"
)

missing <- my_packages[!vapply(my_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing packages: ", paste(missing, collapse = ", "),
       ". Install them (in renv) with install.packages(...).")
}
invisible(lapply(my_packages, library, character.only = TRUE))

## Paths ----
data_dir <- here::here("data", "raw")
sim_path <- file.path(data_dir, "simulations_data.rds")
stopifnot(file.exists(sim_path))

## Load pre-generated simulations (NO re-simulation) ----
sim_data <- readRDS(sim_path)

## Output directory ----
out_dir <- create_run_directory()


## Helpers (aligned with Real_Data.R)----

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

# calculate_metrics <- function(conf_matrix, positive = "2", negative = "1") {
#   cm <- as.matrix(conf_matrix)
#   if (is.null(dim(cm)) || any(dim(cm) != c(2,2))) {
#     return(c(Sensitivity = NA_real_,
#              Specificity = NA_real_,
#              Accuracy = NA_real_,
#              Balanced_Accuracy = NA_real_))
#   }
#   rn <- rownames(cm); cn <- colnames(cm)
#   has_names <- !is.null(rn) && !is.null(cn) && all(c(negative, positive) %in% rn) && all(c(negative, positive) %in% cn)
#   if (has_names) {
#     TN <- cm[negative, negative]; FP <- cm[positive, negative]
#     FN <- cm[negative, positive]; TP <- cm[positive, positive]
#   } else {
#     TN <- cm[1,1]; FP <- cm[2,1]
#     FN <- cm[1,2]; TP <- cm[2,2]
#   }
#   den_sens <- TP + FN
#   den_spec <- TN + FP
#   den_acc  <- TN + FP + FN + TP
#   sensitivity <- ifelse(den_sens > 0, TP/den_sens, NA_real_)
#   specificity <- ifelse(den_spec > 0, TN/den_spec, NA_real_)
#   accuracy    <- ifelse(den_acc  > 0, (TP+TN)/den_acc, NA_real_)
#   balanced_accuracy <- ifelse(is.na(sensitivity) | is.na(specificity), NA_real_, (sensitivity + specificity)/2)
#   c(Sensitivity = sensitivity,
#     Specificity = specificity,
#     Accuracy = accuracy,
#     Balanced_Accuracy = balanced_accuracy)
# }

# cluster_to_labels <- function(score_vec, centers = 2, nstart = 50, seed = 123) {
#   if (any(!is.finite(score_vec))) {
#     stop("cluster_to_labels(): score_vec contains non-finite values (NA/NaN/Inf).")
#   }
#   set.seed(seed)
#   km <- stats::kmeans(as.matrix(score_vec), centers = centers, nstart = nstart)
#   tab <- table(km$cluster)
#   out_cluster <- as.integer(names(tab)[which.min(tab)])
#   pred <- ifelse(km$cluster == out_cluster, "2", "1")
#   factor(pred, levels = c("1","2"))
# }

## DSO / Peña-Prieto / cellwise adjustments (kept from the original Simulaciones_Art.R)----

sd2 <- function(x){ 
  n = length(x)
  sd(x) / sqrt(n) * sqrt(n-1)
}

pena_prieto_transformation <- function(X, dob_base, ridge = 1e-8) {
  X <- as.matrix(X)
  s1 <- cov(X)
  A  <- t(dob_base) %*% s1 %*% dob_base
  I  <- diag(nrow(dob_base))
  inv1 <- tryCatch(
    solve(A),
    error = function(e) { qr.solve(A + ridge * diag(nrow(A))) }
  )
  Q2 <- I - (dob_base %*% inv1 %*% t(dob_base) %*% s1)
  T3 <- X %*% Q2
  return(T3)
}

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

SDC <- function(x, r_i){

  mediana <- apply(x, 2, median)
  p <- ncol(x)
  n <- nrow(x)
  ind_1 <- max(1L, min(n, floor((n + p - 1) / 2)))
  ind_2 <- max(1L, min(n, ceiling((n + p - 1) / 2) + 1))
  beta <- qnorm(0.5 * ((n + p - 1) / (2 * n) + 1))

  mad_modificado <- numeric(length = p)

  # Calcular la MAD modificada para cada componente
  for (i in 1:p) {
    desv_abs <- abs(x[, i] - mediana[i])
    desv_abs_ordenado <- sort(desv_abs)
    mad_modificado[i] <- (desv_abs_ordenado[ind_1] + desv_abs_ordenado[ind_2]) / (2 * beta)
  }
  mad_modificado <- pmax(1e-6, mad_modificado)
  
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
  p <- ncol(x2)
  n <- nrow(x2)
  ind_1 <- max(1L, min(n, floor((n + p - 1) / 2)))
  ind_2 <- max(1L, min(n, ceiling((n + p - 1) / 2) + 1))
  beta <- qnorm(0.5 * ((n + p - 1) / (2 * n) + 1))
  mad_componentes_seleccionados <- apply(x2, 2, mad)
  reescalado <- sweep(x2, 2, mad_componentes_seleccionados, FUN = '/')
  #Se calcula el respectivo r_i
  mu_reescalado <- abs(sweep(reescalado, 2, apply(reescalado, 2, median), "-"))
  s_reescalado  <- pmax(1e-6, apply(mu_reescalado, 2, median))
  r_i_SDM <- sweep(mu_reescalado, 2, s_reescalado, "/")


  #Cálculo del alpha para el método SDM----
  u_i <- pmax(1e-6, apply(r_i_SDM, 1, max))
  alpha_SDM <- sweep(r_i_SDM, 1, u_i, "/")

  # Calcular la MAD modificada para cada componente
  mad_modificado <- numeric(length = p)
  for (i in 1:ncol(x2)) {
    desv_abs <- abs(x2[, i] - mediana[i])
    desv_abs_ordenado <- sort(desv_abs)
    mad_modificado[i] <- (desv_abs_ordenado[ind_1] + desv_abs_ordenado[ind_2]) / (2 * beta)
  }
  mad_modificado <- pmax(1e-6, mad_modificado)

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

## Experiment design (analysis-only; data are pre-generated)----

set.seed(50)
iter <- 50
alpha <- c(0.05, 0.1, 0.2, 0.3, 0.4)
dimension <- c(10, 30, 50)
cant_datos <- 10*dimension
outl_dist <- c(3,4)
outl_concent <- c(0.5)
modelos <- c("FDCM", "FICM")
sim.modes <- 1:3

n_cores <- max(1, parallel::detectCores() - 1)
doParallel::registerDoParallel(n_cores)

all_cm <- list()
all_metrics <- list()

## Run analysis----
for (modelo in modelos) {
  for (mode in sim.modes) {
    for (i2 in seq_along(dimension)) {
      for (i1 in seq_along(alpha)) {
        for (i3 in seq_along(outl_dist)) {
          for (i4 in seq_along(outl_concent)) {
            
            dataset_name <- paste0(modelo, "_mode", mode,
                                   "_p", dimension[i2],
                                   "_alpha", alpha[i1],
                                   "_delta", outl_dist[i3],
                                   "_c", outl_concent[i4])
            key <- dataset_name
            
            if (!key %in% names(sim_data)) {
              stop("Scenario not found in sim_data: ", key, "\n",
                   "Regenerate data with notebooks/Generate_Simulated_Data.R")
            }
            
            scenario_list <- sim_data[[key]]
            if (length(scenario_list) < iter) {
              stop("Scenario ", key, " has only ", length(scenario_list),
                   " replications; expected at least ", iter)
            }
            
            cms_rows <- vector("list", iter)
            
            for (sim_id in seq_len(iter)) {
              
              data.sim <- scenario_list[[sim_id]]
              X <- data.sim$x
              
              # Ground truth: lbl=0 clean (inlier), lbl=1 outlier
              true_labels <- ifelse(as.integer(data.sim$lbl) == 1L, "2", "1")
              true_labels <- factor(true_labels, levels = c("1","2"))
              
              # -------------------- Methods --------------------
              
              # DOBIN transform (used by your DSO pipeline)
              calculated_dobin <- dobin::dobin(X)

              # REPPlab (EPP)
              X_EPP <- REPPlab::EPPlab(X, PPalg="PSO", PPindex="FriedmanTukey",
                                       n.simu=1, maxiter=1000, sphere=TRUE)
              out <- REPPlab::EPPlabOutlier(X_EPP, which = 1:ncol(X), k = 3,
                                            location = mean, scale = sd)
              aux_repp <- unlist(out$outlier)
              threshold_repp <- stats::quantile(aux_repp, 0.95)
              predicted_labels_repp <- ifelse(aux_repp > threshold_repp, "2", "1")
              predicted_labels_repp <- factor(predicted_labels_repp, levels = c("1","2"))
              W1.cm <- caret::confusionMatrix(predicted_labels_repp, true_labels)$table
              
              # HDO (depth-based)
              out2 <- depthTools::hDepth(X)
              threshold_hdod <- stats::quantile(out2, 0.05)
              predicted_labels_hdod <- ifelse(out2 < threshold_hdod, "2", "1")
              predicted_labels_hdod <- factor(predicted_labels_hdod, levels = c("1","2"))
              W2.cm <- caret::confusionMatrix(predicted_labels_hdod, true_labels)$table
              
              # PP-MCD
              W3.safe <- tryCatch({
                out_values <- PPcovMcd::PPcovMcd(X)
                Projection_pursuit_MCD <- stats::mahalanobis(X, center = out_values$center, cov = out_values$cov)
                threshold_MCD <- stats::quantile(Projection_pursuit_MCD, 0.95)
                pred <- ifelse(Projection_pursuit_MCD > threshold_MCD, "2", "1")
                list(cm = caret::confusionMatrix(factor(pred, levels=c("1","2")), true_labels)$table, ok = TRUE)
              }, error = function(e) {
                list(cm = matrix(NA_integer_, 2, 2, dimnames = list(c("1","2"), c("1","2"))), ok = FALSE)
              })
              W3.cm <- W3.safe$cm
              
              # ICS
              data_ics <- ICSOutlier::ics2(X)
              resultado_ICS <- ICSOutlier::ics.outlier(data_ics, level.test = 0.3)
              idx_ics <- resultado_ICS@outliers
              predicted_labels_ics <- rep("1", nrow(X))
              if (length(idx_ics)) predicted_labels_ics[idx_ics] <- "2"
              predicted_labels_ics <- factor(predicted_labels_ics, levels = c("1","2"))
              W4.cm <- caret::confusionMatrix(predicted_labels_ics, true_labels)$table
              
              
              # DSO 
              X_pp <- pena_prieto_transformation(X, calculated_dobin$rotation)
              selected_data <- X_pp
              
              at_m <- calculate_outlyingness(selected_data)
              pred_dso <- cluster_to_labels(at_m)
              W5.cm <- caret::confusionMatrix(pred_dso, true_labels)$table
              
              r_i <- calculate_r_i(selected_data)
              outlyingness_ajustado_SDC <- SDC(selected_data, r_i)
              pred_sdc <- cluster_to_labels(outlyingness_ajustado_SDC)
              W6.cm <- caret::confusionMatrix(pred_sdc, true_labels)$table
              
              outlyingness_ajustado_SDM <- SDM(selected_data)
              pred_sdm <- cluster_to_labels(outlyingness_ajustado_SDM)
              W7.cm <- caret::confusionMatrix(pred_sdm, true_labels)$table
              
              # -------------------- Store --------------------
              cms_rows[[sim_id]] <- rbind(
                cm_to_row(W1.cm, sim_id, dataset_name, "REPPLAB"),
                cm_to_row(W2.cm, sim_id, dataset_name, "HDO"),
                cm_to_row(W3.cm, sim_id, dataset_name, "PP-MCD"),
                cm_to_row(W4.cm, sim_id, dataset_name, "ICS"),
                cm_to_row(W5.cm, sim_id, dataset_name, "DSO"),
                cm_to_row(W6.cm, sim_id, dataset_name, "DSO-SDC"),
                cm_to_row(W7.cm, sim_id, dataset_name, "DSO-SDM")
              )
            }
            
            cm_block <- do.call(rbind, cms_rows)
            
            metrics_rows <- vector("list", nrow(cm_block))
            idx <- 0L
            for (sim_id in sort(unique(cm_block$simulation))) {
              sub_sim <- cm_block[cm_block$simulation == sim_id, ]
              for (mth in unique(sub_sim$method)) {
                row <- sub_sim[sub_sim$method == mth, ]
                cm <- matrix(c(row$TN, row$FN, row$FP, row$TP), nrow=2, byrow=TRUE,
                             dimnames=list(c("1","2"), c("1","2")))
                met <- calculate_metrics(cm)
                idx <- idx + 1L
                metrics_rows[[idx]] <- data.frame(
                  simulation = sim_id, dataset = dataset_name, method = mth,
                  t(met), stringsAsFactors = FALSE
                )
              }
            }
            metrics_rows <- metrics_rows[seq_len(idx)]
            metrics_block <- do.call(rbind, metrics_rows)
            
            all_cm[[key]] <- cm_block
            all_metrics[[key]] <- metrics_block
            
            write.csv(cm_block, file.path(out_dir, paste0("cm_", key, ".csv")), row.names = FALSE)
            write.csv(metrics_block, file.path(out_dir, paste0("metrics_", key, ".csv")), row.names = FALSE)
          }
        }
      }
    }
  }
}

## Consolidate and export----

cm_all <- do.call(rbind, all_cm)
metrics_all <- do.call(rbind, all_metrics)

write.csv(cm_all, file.path(out_dir, "confusion_ALL.csv"), row.names = FALSE)
write.csv(metrics_all, file.path(out_dir, "metrics_ALL.csv"), row.names = FALSE)

saveRDS(
  list(cm = all_cm, metrics = all_metrics),
  file = file.path(out_dir, "all_results.rds")
)

cat("\nOK. Results in: ", out_dir, "\n")


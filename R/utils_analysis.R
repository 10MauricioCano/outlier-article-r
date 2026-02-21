# ============================================================
# utils_analysis.R
# Shared analysis utilities for Real_Data.R and Simulaciones_Art.R
# ============================================================

#-------------------------------------------------------------
# Create timestamped output directory
#-------------------------------------------------------------
create_run_directory <- function(base_path = here::here("outputs")) {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  run_dir <- file.path(base_path, paste0("run_", timestamp))
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  return(run_dir)
}

#-------------------------------------------------------------
# Cluster top alpha proportion as outliers
#-------------------------------------------------------------

cluster_to_labels <- function(score_vec, centers = 2, nstart = 50, seed = 123) {
  # Convert a numeric outlyingness score into labels via 1D k-means.
  # Outliers are defined as the cluster with the highest mean score.
  stopifnot(is.numeric(score_vec))
  if (length(score_vec) < 2L) {
    return(factor(rep("1", length(score_vec)), levels = c("1","2")))
  }
  if (any(!is.finite(score_vec))) {
    stop("cluster_to_labels(): score_vec contains non-finite values (NA/NaN/Inf).")
  }
  if (length(unique(score_vec)) < 2L) {
    # No separability in the score vector; default to all inliers.
    return(factor(rep("1", length(score_vec)), levels = c("1","2")))
  }
  
  # Preserve RNG state for reproducibility without contaminating the caller.
  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (old_seed_exists) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (old_seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  
  set.seed(seed)
  km <- stats::kmeans(score_vec, centers = centers, nstart = nstart)
  
  cl <- km$cluster
  cl_means <- tapply(score_vec, cl, mean)
  out_cluster <- as.integer(names(which.max(cl_means)))
  
  pred <- ifelse(cl == out_cluster, "2", "1")
  factor(pred, levels = c("1","2"))
}

# cluster_to_labels <- function(score_vec, alpha = 0.02, alpha_max = 0.05) {
#   alpha <- min(alpha, alpha_max)
#   n <- length(score_vec)
#   k <- max(1L, ceiling(alpha * n))
#   ord <- order(score_vec, decreasing = TRUE)
#   pred <- rep("1", n)
#   pred[ord[seq_len(k)]] <- "2"
#   factor(pred, levels = c("1", "2"))
# }

#-------------------------------------------------------------
# Compute metrics from a 2x2 confusion matrix
#-------------------------------------------------------------
calculate_metrics <- function(conf_matrix, positive = "2", negative = "1") {
  cm <- as.matrix(conf_matrix)
  if (is.null(dim(cm)) || any(dim(cm) != c(2, 2))) {
    return(c(Sensitivity = NA_real_, Specificity = NA_real_,
             Accuracy = NA_real_, Balanced_Accuracy = NA_real_))
  }
  rn <- rownames(cm); cn <- colnames(cm)
  has_names <- !is.null(rn) && !is.null(cn) &&
    all(c(negative, positive) %in% rn) && all(c(negative, positive) %in% cn)
  if (has_names) {
    TN <- cm[negative, negative]; FP <- cm[positive, negative]
    FN <- cm[negative, positive]; TP <- cm[positive, positive]
  } else {
    TN <- cm[1,1]; FP <- cm[2,1]; FN <- cm[1,2]; TP <- cm[2,2]
  }
  den_sens <- TP + FN; den_spec <- TN + FP; den_acc <- TN + FP + FN + TP
  sensitivity       <- ifelse(den_sens > 0, TP / den_sens, NA_real_)
  specificity       <- ifelse(den_spec > 0, TN / den_spec, NA_real_)
  accuracy          <- ifelse(den_acc  > 0, (TP + TN) / den_acc, NA_real_)
  balanced_accuracy <- ifelse(is.na(sensitivity) | is.na(specificity),
                              NA_real_, (sensitivity + specificity) / 2)
  c(Sensitivity = sensitivity, Specificity = specificity,
    Accuracy = accuracy, Balanced_Accuracy = balanced_accuracy)
}

#-------------------------------------------------------------
# Compute confusion matrix and Balanced Accuracy
#-------------------------------------------------------------
compute_confusion_metrics <- function(pred_labels, true_labels) {
  
  cm <- caret::confusionMatrix(pred_labels, true_labels)$table
  
  TP <- cm["2","2"]
  TN <- cm["1","1"]
  FP <- cm["2","1"]
  FN <- cm["1","2"]
  
  TPR <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
  TNR <- ifelse((TN + FP) > 0, TN / (TN + FP), 0)
  
  BA <- (TPR + TNR) / 2
  
  list(
    confusion = cm,
    TPR = TPR,
    TNR = TNR,
    BA  = BA
  )
}

#-------------------------------------------------------------
# Save consolidated outputs
#-------------------------------------------------------------
save_consolidated_results <- function(conf_list, metrics_df, run_dir) {
  
  write.csv(
    do.call(rbind, conf_list),
    file = file.path(run_dir, "confusion_ALL.csv"),
    row.names = TRUE
  )
  
  write.csv(
    metrics_df,
    file = file.path(run_dir, "metrics_ALL.csv"),
    row.names = FALSE
  )
  
  saveRDS(
    list(confusion = conf_list, metrics = metrics_df),
    file = file.path(run_dir, "all_results.rds")
  )
}

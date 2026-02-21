library(here)

# One-time generator for simulated scenarios used by notebooks/Simulaciones_Art.R
# It sources R/OL_skew_self_V2.R (kept intact) and saves all replications to data/raw/simulations_data.rds

data_dir <- here::here("data", "raw")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

## Required Packages----
my_packages = c("extraDistr")
missing <- my_packages[!vapply(my_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing packages: ", paste(missing, collapse = ", "),
       ". Install them (in renv) with install.packages(...).")
}
invisible(lapply(my_packages, library, character.only = TRUE))

## Source simulation generator (DO NOT EDIT)----
source(here::here("R", "OL_skew_self_V2.R"))

## Design (must match Simulaciones_Art.R)----
set.seed(50)
iter <- 50
alpha <- c(0.05, 0.1, 0.2, 0.3, 0.4)
dimension <- c(10, 30, 50)
cant_datos <- 10*dimension
outl_dist <- c(3,4)
outl_concent <- c(0.5)
modelos <- c("FDCM", "FICM")
sim.modes <- 1:3

sim_data <- list()
meta <- list()

# Use deterministic seeds per replication for full reproducibility
base_seed <- 1000L
rep_id <- 0L

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
            scenario_list <- vector("list", iter)
            for (sim_id in seq_len(iter)) {
              rep_id <- rep_id + 1L
              set.seed(base_seed + rep_id)
              if (modelo == "FDCM") {
                scenario_list[[sim_id]] <- GenAtip(cant_datos[i2], dimension[i2],
                                                  c(alpha[i1], outl_dist[i3], outl_concent[i4]),
                                                  sim.mode = mode)
              } else {
                scenario_list[[sim_id]] <- GenAtip_FICM(cant_datos[i2], dimension[i2],
                                                      c(alpha[i1], outl_dist[i3], outl_concent[i4]),
                                                      sim.mode = mode)
              }
            }
            sim_data[[dataset_name]] <- scenario_list
            meta[[dataset_name]] <- data.frame(
              dataset = dataset_name,
              modelo = modelo,
              mode = mode,
              p = dimension[i2],
              n = cant_datos[i2],
              alpha = alpha[i1],
              delta = outl_dist[i3],
              concentration = outl_concent[i4],
              iter = iter,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }
}

meta_df <- do.call(rbind, meta)
saveRDS(sim_data, file.path(data_dir, "simulations_data.rds"))
write.csv(meta_df, file.path(data_dir, "simulations_data_meta.csv"), row.names = FALSE)
cat("OK. Saved:", file.path(data_dir, "simulations_data.rds"), "\n")

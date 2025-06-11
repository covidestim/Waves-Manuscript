rm(list = ls())
gc()

library(tidyverse)
library(sf)
library(spdep)
library(spatialreg)
library(Matrix)

setwd("~/Desktop/repos/Waves/Waves-Manuscript/")

## Hexgrid pop
## New hexgrid with Meta 30m population
# hexgrid_pop <- st_read("Data/data-sources/hex_pop_meta_new.shp") |> 
#   rename(population = meta_pop) |>
#   filter(population != 0) ## Filter out the 0 population hexes to never estimate infectionsPC there

## Pre-Omicron
## Remember to rebuild the full dataset from the four parts in code 01-hexAllocation_Infections_Meta30m.R
# hexgrid_preomicron <- st_read("Data/data-products/hexid-observations_preomicron06-05-25.geojson")

## Pre-Omicron expanded dataset
# hex_spacetime <- expand.grid(hexid = as.character(unique(hexgrid_pop$hexid)),
#                              date = seq.Date(from = min(hexgrid_preomicron$date),
#                                              to = max(hexgrid_preomicron$date), 
#                                              by = "day")) |> 
#   left_join(hexgrid_preomicron |> 
#               st_drop_geometry() |> 
#               select(hexid, date, infections, infectionsPC)) |>  ## just join the desired columns, we have the geometry
#   # mutate(Time = as.numeric(date - min(date)) + 1,
#   #        # ID = as.numeric(ID),
#   #        logpopulation = log10(population+1)) |>
#   full_join(hexgrid_pop |> 
#               st_drop_geometry() |> 
#               mutate(hexid = as.character(hexid), 
#                      ID = 1:n()), ## Creating a numeric index for the hexgrid, 
#                                   ## watch to this when back joining anything from the model
#             by = "hexid") 

## Saving the grid to facilitate re-runs
# vroom::vroom_write(x = hex_spacetime, 
#                    file = "Data/data-products/hex_spacetime.csv")

hex_spacetime <- vroom::vroom("Data/data-products/hex_spacetime.csv") |> 
  mutate(hexid = as.character(hexid))

## Flag if it is a rerun
is.rerun <- FALSE

## This should be run in batches, as it may take long time and not be suit for memory size.
## The fitting strategy is to run the model for all dates and rerun over the problematic dates, when the model didn't fit well
## One very common issue is returning a bad run from INLA, 
## which only the iid part fits to the model and not the Besag model.
## The Newton-Raphson method could not converged, the method is used to some numerical integration on the hyperparameters. 
## Another alternative is to log-transform the data before imputing to the INLA model, 
## Hävard Rue suggested this on a help desk response at INLA email list. 
## Yet another alternative is to set the prior parameterization to the model that can afford huge variance on the data ## Omicron dataset case. We can set the shape and rate parameters of log-Gamma prior to have mean = 1 and variance = 1.
# 
# hexgrid <- hex_spacetime |>
#   left_join(hexgrid_pop |> select(hexid, geometry)) |>
#   st_as_sf() |>
#   group_by(ID) |>
#   collapse::fsummarise(geometry = st_union(geometry))
# 
# ## Diagonostic plot to check if the hex ID's are correct, should be a gradient of one color
# ggplot()+
#   geom_sf(data = hexgrid, 
#           aes(fill = ID))

# sf::st_write(obj = hexgrid,
#          dsn = "Data/data-products/hexgrid.geojson",
#          delete_dsn = T,
#          delete_layer = T)

hexgrid <- sf::st_read("Data/data-products/hexgrid.geojson")

## Neighbors
hexes_nb <- spdep::poly2nb(hexgrid, queen = TRUE, row.names = hexgrid$hexid)
# hexes_nb2 <- st_touches(st_geometry(hexes))

nb2INLA("Data/data-products/hexes_adjmat.graph", nb = hexes_nb)
hexes_graph <- INLA::inla.read.graph("Data/data-products/hexes_adjmat.graph")

# wt_B <- nb2mat(neighbours = hexes_nb, style = "B", zero.policy = T)
# wt_W <- nb2mat(neighbours = hexes_nb, style = "W", zero.policy = T)
# 
# listB <- nb2listw(neighbours = hexes_nb, style = "B", zero.policy = T)
# listW <- nb2listw(neighbours = hexes_nb, style = "W", zero.policy = T)
# 
# ## Weight list as a Sparse C Matrix
# B <- as(listB, "CsparseMatrix")
# W <- as(listW, "CsparseMatrix")

## CAR INLA equivalent model
library(INLA)
compute_list <- control.compute(hyperpar = T, 
                                return.marginals = T, 
                                return.marginals.predictor = T,
                                config = T,
                                openmp.strategy = "huge",
                                dic = T, 
                                cpo = T, 
                                waic = T)

predictor_list <- control.predictor(compute = T)

# library(foreach)
# library(doParallel)
library(dplyr)
library(tibble)

hyper_smooth <- list(
  prec = list(
    prior = "pc.prec",
    param = c(0.2, 0.01),  # P(SD > 0.2) = 1% (STRONG smoothing)
    initial = 5            # Start with high precision (1/exp(4) ≈ 0.018 SD)
  )
)

hyper_smooth_bym2 <- list(
  phi = list(prior = "pc", param = c(0.95, 0.7)),  # 50% prob ϕ > 0.95
  prec = list(prior = "pc.prec", param = c(0.2, 0.01))
)

# diag.eps = 1e-3

# # Set up parallel backend
# cl <- makeCluster(detectCores() - 2)
# registerDoParallel(cl)
# 
# # # Preprocess adjacency matrix on all workers
# clusterEvalQ(cl, {
#   library(INLA)
#   library(dplyr)
#   library(tibble)
#   library(base)
#   library(stats)
#   hex_graph <- inla.read.graph("Data/data-products/hexes_adjmat.graph")
# })
# 
# ## Export required objects
# clusterExport(cl, c("hex_spacetime", "hyper_smooth", "hyper_smooth_bym2", "compute_list", "predictor_list", "W"))

weeks <- sort(unique(na.omit(hex_spacetime$date)))

## Uncomment when is not a rerun
CAR_list <- list()

# alpha_peak <- as.Date("2020-11-19")
# delta_peak <- as.Date("2021-09-04")
# # 
# test_dates <- c(c(alpha_peak-63,
#                   alpha_peak-45,
#                   alpha_peak-24,
#                   alpha_peak),
#                 c(delta_peak-63,
#                   delta_peak-45,
#                   delta_peak-24,
#                   delta_peak))
# 
# which(weeks == alpha_peak)
# which(weeks == delta_peak)

## This fitting takes very long, it is process of prune fitting. 
## Where it starts fitting something and narrows down to a best fit, which parameteres are defined at the while loop

# CAR_list <- foreach(i = 1:length(weeks),
#                     .combine = "list",
#                     .multicombine = T) %dopar% {

for (i in 20:length(weeks)) {
  
  # if(is.na(dates_to_rerun[i]))next
  
  current_date <- weeks[i]
  cat("Starting model for date: ", as.character(current_date),"! \n")
  hex_week <- hex_spacetime %>% 
    filter(date == current_date)
  
  best_model <- NULL
  # counter <- 0
  # sd_values <- c(0.0001, 0.03)  # Start with high value
  
  # while ((sd(sd_values)<0.0025 || sd(sd_values) > 0.010) & counter <= 50) {  # Flipped condition
  tryCatch({
    set.seed(.Random.seed)
    
    best_model <- inla(
      as.formula(infectionsPC ~ 1 + 
                   f(ID, 
                     model = "bym2",
                     graph = hexes_graph,
                     scale.model = TRUE,
                     # diagonal = diag.eps,
                     constr = TRUE,
                     hyper = hyper_smooth_bym2)),
      data = as.data.frame(hex_week),
      family = "gaussian",
      control.inla = control.inla(strategy = "gaussian", int.strategy = "eb"),
      # control.inla = control.inla(strategy = "gaussian", h = h.value, int.strategy = "eb"),
      # control.mode = control.mode(restart = TRUE),
      control.compute = compute_list,
      control.predictor = predictor_list,
      # control.fixed = list(prec.intercept = 0.1),
      # num.threads = 1,  # Prevent internal threading conflicts
      # num.threads = 2,
      # verbose = T
    )
    
    # sd_values <- best_model$summary.fitted.values$sd
    
  }, error = function(e) {
    set.seed(.Random.seed)
    
    best_model <- inla(
      as.formula(infectionsPC ~ 1 + 
                   f(ID, 
                     model = "bym2",
                     graph = hexes_graph,
                     # diagonal = diag.eps,
                     scale.model = TRUE,
                     constr = TRUE,
                     hyper = hyper_smooth_bym2)),
      data = as.data.frame(hex_week),
      family = "gaussian",
      control.inla = control.inla(strategy = "gaussian", int.strategy = "eb"),
      # control.inla = control.inla(strategy = "gaussian", h = h.value, int.strategy = "eb"),
      # control.mode = control.mode(restart = TRUE),
      control.compute = compute_list,
      control.predictor = predictor_list,
      # control.fixed = list(prec.intercept = 0.1),
      # num.threads = 1,  # Prevent internal threading conflicts
      # num.threads = 2,
      # verbose = T
    )
    # sd_values <- best_model$summary.fitted.values$sd
  })
  
  # counter <- counter + 1
  # hess.start <- which(best_model$logfile == 'Eigenvalues of the Hessian')
  # hess.min <- min(as.numeric(best_model$logfile[(hess.start+1):(hess.start+3)]))
  # h.value <- h.trials + 0.001
  # h.trials <- h.value
  # trials <- trials + 1
  # }
  
  cat("Finished CAR model for week ", as.character(current_date),"! \n")
  
  CAR_list[[i]] <-
    # result <- 
    cbind(hex_week, 
          best_model$summary.fitted.values |> 
            rownames_to_column(var = "INLApred"))
  # )
  

  rm(best_model, hex_week)
  gc()
  # result
  # )
}

# stopCluster(cl)

# ## Saving as list object
# save(list = CAR_list, 
#      file = "Data/data-products/CAR_list_meta30m.RDS", 
#      compress = "xz", 
#      compression_level = 9)

## Turning into a df
CAR_df <- bind_rows(CAR_list)

## Saving the df, remember to change the name if the dataset is "preomicron" or "omicronera". The pattern nomenclature to files are tsa_preomicron.csv or tsa_omicronera.csv
dataset <- "all_dates_meta30m_run_preomicron_daily"

vroom::vroom_write(x = CAR_df, 
                   file = paste0("Data/data-products/tsa_", 
                                 dataset,
                                 ".csv"))

# sf::st_write(obj = CAR_df,
#              dsn = "Data/data-products/tsa_meta30m_run_preomicron_daily.geojson",
#              delete_dsn = T,
#              delete_layer = T)

#

## Breakdowns of each peaks
# breaks_plt1 <- seq(130,310, 20)
# labels_plt1 <- c("<150",seq(150, 290, 20), '>300')
# limits_plt1 <- c(0,350)
# color_option <- "magma"

# CAR_df <- bind_rows(CAR_list) |> 
#   left_join(hexgrid) |> 
#   st_as_sf()

# ggplot() +
#   geom_sf(data = CAR_df |> 
#             filter(date %in% test_dates),
#           aes(fill = mean))+
#   # scale_fill_viridis_c()+
#   scale_fill_viridis_b(option = color_option,
#                        name = "Estimated Infections/100k/day",
#                        direction = -1,
#                        breaks = breaks_plt1,
#                        labels = labels_plt1,
#                        na.value = "steelblue4",
#                        limits = limits_plt1
#   )+
#   theme_minimal()+
#   theme(legend.position = "bottom",
#         legend.title.position = "top",
#         legend.key.width = grid::unit(1, "in"))+
#   facet_wrap(.~date, nrow = 2)


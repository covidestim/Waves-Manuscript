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
# is.rerun <- FALSE

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
  # openmp.strategy = "huge",
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
    
  ## Uncomment when is not a rerun
  # CAR_list <- list()
    
  # weeks <- sort(unique(na.omit(hex_spacetime$date)))
    
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
    
## This fitting takes very long, it is done by running the model to all the dates, after checking which dates the model converge to a smooth solution or when it flips to iid part. 
## Where it starts fitting something and narrows down to a best fit, which parameteres are defined at the while loop
    
# CAR_list <- foreach(i = 1:length(weeks),
#                     .combine = "list",
#                     .multicombine = T) %dopar% {

## When rerunning some dates, uncomment below
CAR_df <- vroom::vroom("Data/data-products/tsa_meta30m_run_preomicron_daily.csv") 
  
# ## Rebuilding as list to run the models
CAR_list <- CAR_df |> 
mutate(hexid = as.character(hexid)) |> 
group_split(date)
    
weeks <- sort(unique(na.omit(CAR_df$date)))
length(weeks)

detect_sd_aberration_with_recovery <- function(upper_sd,
  median_sd,
  median_threshold = 5,
  range_threshold = 5,
  recovery_margin = 0.75) {
n <- length(median_sd)
if (!is.numeric(upper_sd) || !is.numeric(median_sd) || length(upper_sd) != n || n < 3) {
stop("upper_sd and median_sd must be same-length numeric vectors ≥ 3")
}

result <- rep(FALSE, n)
in_recovery <- FALSE
baseline_median <- NA

for (i in 2:(n - 1)) {
# Primary aberration detection
drop_prev  <- (median_sd[i - 1] - median_sd[i]) > median_threshold
drop_next  <- (median_sd[i + 1] - median_sd[i]) > median_threshold
median_drop <- drop_prev || drop_next

spike_prev <- (upper_sd[i] - upper_sd[i - 1]) > range_threshold
spike_next <- (upper_sd[i] - upper_sd[i + 1]) > range_threshold
upper_spike <- spike_prev && spike_next

# Flag primary aberration
if (median_drop || upper_spike) {
result[i] <- TRUE
in_recovery <- TRUE
# Set baseline to last good median
baseline_median <- median_sd[i - 1]
} else if (in_recovery) {
# Continue flagging until recovery
if (!is.na(baseline_median) && median_sd[i] < baseline_median * recovery_margin) {
result[i] <- TRUE
} else {
in_recovery <- FALSE  # Median recovered
}
}
}

return(result)
}

## List to rerun if the model is not well fitted to the data
sd_values <- data.frame(weeks = sort(unique(CAR_df$date)),
upper_sd = sapply(CAR_list, function(x){x <- range(x$sd)[2]}),
median_sd = sapply(CAR_list, function(x){x <- median(x$sd)}),
lower_sd = sapply(CAR_list, function(x){x <- range(x$sd)[1]}))

sd_values$aberration <- detect_sd_aberration_with_recovery(
  upper_sd = sd_values$upper_sd,
  median_sd = sd_values$median_sd,
  median_threshold = 1.5,
  range_threshold = 5,
  recovery_margin = 0.85  # you can tune this!
)

dates_to_rerun <- ifelse(sd_values$aberration, sd_values$weeks, NA)
length_dates_to_rerun <- length(na.omit(dates_to_rerun))

cat("Will rerun for :", length_dates_to_rerun, "dates! \n")
    
for (i in 1:length(weeks)) {
      
  if(is.na(dates_to_rerun[i]))next
      
  current_date <- weeks[i]
  cat("Starting model for date: ", as.character(current_date),"! \n")
  hex_week <- hex_spacetime %>% 
                filter(date == current_date)
      
  best_model <- NULL
  ## counter for while
  counter <- 0
  while_condition <- TRUE
      
    while ((while_condition) & counter <= 10) {  # Flipped condition

      cat("Attempt to fit the model number: ", counter, "!\n")
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
                control.mode = control.mode(restart = TRUE),
                control.compute = compute_list,
                control.predictor = predictor_list,
                # control.fixed = list(prec.intercept = 0.1),
                num.threads = 2,  # Prevent internal threading conflicts
                # verbose = T
              )
              
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
                    control.mode = control.mode(restart = TRUE),
                    control.compute = compute_list,
                    control.predictor = predictor_list,
                    # control.fixed = list(prec.intercept = 0.1),
                    num.threads = 2,  # Prevent internal threading conflicts
                    # verbose = T
              )
            })
                
      counter <- counter + 1
      # Update while condition
        if(i>=2){
          ## Create a vector of size 2 that will check for big drop in the median value for the sd
          vec_1 <- data.frame(
            upper_sd  = range(CAR_list[[i-1]]$sd)[2],
            median_sd = median(CAR_list[[i-1]]$sd))
            
          vec_2 <- data.frame(
            upper_sd  = range(best_model$summary.fitted.values$sd)[2],
            median_sd = median(best_model$summary.fitted.values$sd))
              
          vec_3 <- data.frame(
              upper_sd  = range(CAR_list[[i+1]]$sd)[2],
              median_sd = median(CAR_list[[i+1]]$sd))
                
          vec <- rbind(vec_1, vec_2, vec_3)
                
          # Position 2 is the rerun position
          # 1) big drop in median vs. either neighbor
          drop_prev <- (vec$median_sd[1] - vec$median_sd[2]) > 3
          drop_next <- (vec$median_sd[3] - vec$median_sd[2]) > 3
          median_drop <- drop_prev || drop_next
              
          # 2) big spike in upper vs. both neighbors
          spike_prev <- (vec$upper_sd[2] - vec$upper_sd[1]) > 3
          spike_next <- (vec$upper_sd[2] - vec$upper_sd[3]) > 3
          upper_spike <- spike_prev && spike_next
                      
          while_condition <- median_drop || upper_spike
        }else{
          while_condition <- TRUE}
    }
              
  cat("Finished CAR model for week ", as.character(current_date),"! \n")
              
  CAR_list[[i]] <-
    cbind(hex_week, 
      best_model$summary.fitted.values |> 
      rownames_to_column(var = "INLApred"))
                
  rm(best_model, hex_week)
  gc()
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
dataset <- "meta30m_run_preomicron_daily"
              
vroom::vroom_write(x = CAR_df, 
  file = paste0("Data/data-products/tsa_", 
  dataset,
  ".csv"))
                
# sf::st_write(obj = CAR_df,
#              dsn = "Data/data-products/tsa_meta30m_run_preomicron_daily.geojson",
#              delete_dsn = T,
#              delete_layer = T)

#

# sd_values <- data.frame(weeks = sort(unique(CAR_df$date)),
# upper_sd = sapply(CAR_list, function(x){x <- range(x$sd)[2]}),
# median_sd = sapply(CAR_list, function(x){x <- median(x$sd)}),
# lower_sd = sapply(CAR_list, function(x){x <- range(x$sd)[1]}))

# check_large_drop_recursive <- function(vec, threshold = 5) {
#   # Internal recursive function
#   helper <- function(v) {
#     if (length(v) < 2) {
#       return(logical(0))  # No comparison possible
#     } else {
#       drop_forward <- (v[2] - v[1]) > threshold
#       drop_backward <- (v[1] - v[2]) < threshold
#       drop <- drop_forward == drop_backward
#       return(c(drop, helper(v[-1])))
#     }
#   }
  
#   # Input check
#   if (!is.numeric(vec) || length(vec) < 2) {
#     stop("Input must be a numeric vector with at least two elements.")
#   }
  
#   # Run helper function and pad result to match input length
#   result <- helper(vec)
#   return(c(result, FALSE))  # Last element has no next value
# }

# sd_values <- sd_values |> 
# mutate(aberration_upper = check_large_drop_recursive(upper_sd, threshold = 2),
# aberration_median = check_large_drop_recursive(median_sd, threshold = 2))

# alpha_peak <- as.Date("2020-11-19")
# delta_peak <- as.Date("2021-09-04")

# ggplot(data = sd_values, 
#   aes(x = weeks, y = median_sd,
#     ymax = upper_sd, ymin = lower_sd,
#     color = if_else(aberration, "firebrick", "steelblue")
#   ))+
#   geom_vline(xintercept = alpha_peak, color = "firebrick")+
#   annotate(geom = "text", 
#   label = "1st Wave Peak", 
#   x = alpha_peak + 5, y = 80, 
#   color = "firebrick", 
#   angle = 90)+
#   geom_vline(xintercept = delta_peak, color = "steelblue")+
#   annotate(geom = "text", 
#   label = "2nd Wave Peak", 
#   x = delta_peak + 5, y = 80, color = "steelblue", 
#   angle = 90)+
#   # geom_label()+
#   # geom_point()+
#   geom_pointrange()+
#   theme_minimal()+
#   theme(legend.title = element_blank())
  
# ## Weeks vector
# weeks <- sort(unique(na.omit(CAR_df$date)))
  
# ## List to rerun if the model is not well fitted to the data
# dates_to_rerun <- ifelse(check_large_drop_recursive(sd_values$median_sd), weeks, NA)
# dates_to_rerun
# ## How many reruns
# length(na.omit(dates_to_rerun))
                
# ## Breakdowns of each peaks
# breaks_plt1 <- seq(130,310, 20)
# labels_plt1 <- c("<150",seq(150, 290, 20), '>300')
# limits_plt1 <- c(0,350)
# color_option <- "magma"
                
# CAR_df <- bind_rows(CAR_list) |> 
#                 left_join(hexgrid) |> 
#                 st_as_sf()
                
# ggplot() +
#   geom_sf(data =  CAR_list[[119]]|> 
#                   left_join(hexgrid) |> 
#                   st_as_sf(),
#           aes(fill = mean))+
#           # scale_fill_viridis_c()+
#           scale_fill_viridis_b(option = color_option,
#             name = "Estimated Infections/100k/day",
#             direction = -1,
#             breaks = breaks_plt1,
#             labels = labels_plt1,
#             na.value = "steelblue4",
#             limits = limits_plt1
#           )+
#           theme_minimal()+
#           theme(legend.position = "bottom",
#           legend.title.position = "top",
#           legend.key.width = grid::unit(1, "in"))+
#           facet_wrap(.~date, nrow = 2)
                  

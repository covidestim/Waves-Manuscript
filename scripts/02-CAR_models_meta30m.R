rm(list = ls())
gc()

library(tidyverse)
library(sf)
library(spdep)
library(spatialreg)
library(Matrix)

## Hexgrid pop
## New hexgrid with Meta 30m population
hexgrid_pop <- st_read("Data/data-sources/hexgrid_meta30m_population.geojson") |> 
  rename(population = metapop_30m)

## Only keep hexes with a less than 10 cumulative infections per capita
hexgrid_preomicron_cum <- vroom::vroom("Data/data-products/hexid-observations_preomicron05-31-25NOGEOM.csv") |> 
  mutate(hexid = as.character(hexid)) |> 
  group_by(hexid) |> 
  summarise(cum_infections = sum(infections, na.rm = T)) |> 
  right_join(hexgrid_pop |>  
               mutate(hexid = as.character(hexid)) |> 
               select(hexid, population))  |> 
  mutate(cum_infectionsPC = cum_infections/population) |> 
  sf::st_as_sf() |> 
  st_transform(crs = 26915)|>
  dplyr::mutate(logpopulation = log10(population),
                cum_incidence = exp(log10(cum_infections) - logpopulation),
                log_incidence = log10(cum_incidence+1))

# ### Cumulative plot test
ggplot() +
  geom_sf(hexgrid_preomicron_cum |>
            # filter(population >0 & cum_infectionsPC <=10)|> ## uncomment if you wanna check the flatten figures
            st_transform(crs = 'ESRI:102009'),
          mapping=aes(fill = log10(population+1))) +
  geom_sf(us_states,
          mapping=aes(),
          color = "black",
          fill = "transparent")+
  scico::scale_fill_scico(name = "Cumulative Infections (log10 scale) \n (March 2020 - December 2021)",
                       palette = "lapaz",
                       direction = 1,
                       na.value = "grey70",
                       breaks = seq(0,7,1),
                       labels = scales::label_math(),
                       limits = c(0,7)
  )+
  theme_minimal()+
  theme(legend.position = "bottom",
        legend.title.position = "top",
        legend.title = element_text(hjust = 0.5),
        legend.key.width = grid::unit(3, "cm"),
        axis.text = element_text(size = 6))

## Neighbors
hexes_nb <- spdep::poly2nb(hexgrid_pop, queen = TRUE, row.names = hexgrid_pop$hexid)
# hexes_nb2 <- st_touches(st_geometry(hexes))

nb2INLA("Data/data-products/hexes_adjmat.graph", nb = hexes_nb)
hexes_graph <- INLA::inla.read.graph("Data/data-products/hexes_adjmat.graph")

# wt_B <- nb2mat(neighbours = hexes_nb, style = "B", zero.policy = T)
wt_W <- nb2mat(neighbours = hexes_nb, style = "W", zero.policy = T)
# 
# listB <- nb2listw(neighbours = hexes_nb, style = "B", zero.policy = T)
listW <- nb2listw(neighbours = hexes_nb, style = "W", zero.policy = T)
# 
# ## Weight list as a Sparse C Matrix
# B <- as(listB, "CsparseMatrix")
W <- as(listW, "CsparseMatrix")

## Pre-Omicron
## Remember to rebuild the full dataset from the four parts in code 01-hexAllocation_Infections_Meta30m.R
hexgrid_preomicron <- vroom::vroom("Data/data-products/hexid-observations_preomicron05-31-25NOGEOM.csv") |> 
  mutate(hexid = as.character(hexid)) 

# ### Cumulative plot test
ggplot() +
  geom_sf(dplyr::full_join(hexgrid_preomicron |> 
                             filter(date == "2021-09-08"), 
                           metapop  |>  
                             mutate(hexid = as.character(hexid)) |> 
                             dplyr::select(!metapop_30m), 
                           by=c("hexid")) %>% 
            sf::st_as_sf() |> 
            st_transform(crs = 'ESRI:102009'),
          mapping=aes(fill = log10(infections+1))) +
  geom_sf(us_states,
          mapping=aes(),
          color = "black",
          fill = "transparent")+
  scico::scale_fill_scico(name = "Cumulative Infections (log10 scale) \n (March 2020 - December 2021)",
                          palette = "lapaz",
                          direction = 1,
                          na.value = "grey70",
                          # breaks = seq(0,7,1),
                          # labels = scales::label_math(),
                          # limits = c(0,7)
  )+
  theme_minimal()+
  theme(legend.position = "bottom",
        legend.title.position = "top",
        legend.title = element_text(hjust = 0.5),
        legend.key.width = grid::unit(3, "cm"),
        axis.text = element_text(size = 6))

## Pre-Omicron expanded dataset
hex_spacetime <- expand.grid(hexid = as.character(unique(hexgrid_pop$hexid)),
                             date = seq.Date(from = min(hexgrid_preomicron$date),
                                             to = max(hexgrid_preomicron$date), 
                                             by = "day")) |> 
  left_join(hexgrid_preomicron |> 
              st_drop_geometry() |> 
              select(hexid, date, infections, population, infectionsPC)) |>  ## just join the desired columns, we have the geometry
  mutate(Time = as.numeric(date - min(date)) + 1,
         ID = as.numeric(hexid),
         logpopulation = log10(population+1))

## Certifying the correct number of unique hex; 7352(?)
length(unique(na.omit(hex_spacetime$hexid)))

## Returning the data.frame into a list format
# CAR_list_rerun <- CAR_df_rerun |> 
#   group_split(date) |> 
#   as.list()

## Reruns accouting
## 499:628 ok!
## 398:498 ok!
## 297:397 ok!
## 196:296 ok!
## 95:195 ok!
## 1:94 ok!

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

library(foreach)
library(doParallel)
library(dplyr)
library(tibble)

# # Set up parallel backend
cl <- makeCluster(detectCores() - 2)
registerDoParallel(cl)

hyper_smooth <- list(
  prec = list(
    prior = "pc.prec",
    param = c(0.2, 0.01),  # P(SD > 0.2) = 1% (STRONG smoothing)
    initial = 5            # Start with high precision (1/exp(4) ≈ 0.018 SD)
  )
)

hyper_smooth_bym2 <- list(
  phi = list(prior = "pc", param = c(0.95, 0.5)),  # 50% prob ϕ > 0.95
  prec = list(prior = "pc.prec", param = c(0.2, 0.01))
)

# diag.eps = 1e-3

# # Preprocess adjacency matrix on all workers
clusterEvalQ(cl, {
  library(INLA)
  library(dplyr)
  library(tibble)
  library(base)
  library(stats)
  hex_graph <- inla.read.graph("Data/data-products/hexes_adjmat.graph")
})
# 
# # Export required objects
clusterExport(cl, c("hex_spacetime", "hyper_smooth", "hyper_smooth_bym2", "compute_list", "predictor_list"))

weeks <- sort(unique(na.omit(hex_spacetime$date)))

hex_graph <- INLA::inla.read.graph("Data/data-products/hexes_adjmat.graph")

CAR_df <- vroom::vroom("Data/data-products/tsa_meta30m_run_preomicron_daily.csv")

CAR_list <- CAR_df |> 
  group_split(date)

## Loading from the RDS
CAR_list <- load("Data/data-products/CAR_list_meta30m.RDS")

## Rerun flag, set this to TRUE is the dates to rerun are a lot to get better model estimates
is.rerun <- TRUE

## setting the hexid column on the same class
CAR_list <- lapply(CAR_list, function(x){x <- x |> mutate(hexid = as.integer(hexid))})

## List to rerun if the model is not well fitted to the data
dates_to_rerun <- sapply(CAR_list, function(x){ifelse(sd(x$sd)<0.0025 || sd(x$sd)>0.010, x$date, NA)})
dates_to_rerun
length(na.omit(dates_to_rerun))

sd_values <- data.frame(weeks = test_dates, 
                        sd = sapply(CAR_list, function(x){sd(x$sd)}),
                        upper = sapply(CAR_list, function(x){max(x$sd)}),
                        lower = sapply(CAR_list, function(x){sd(x$sd)})) |> 
  mutate(id = row_number())

ggplot(data = sd_values, 
       aes(x = weeks, y = sd, 
           label = id,
           # ymin = lower, ymax = upper,
           color = if_else(sd <= 0.01 & sd >= 0.0025, "firebrick", "steelblue"),
       ))+
  geom_label()+
  # geom_pointrange()+
  theme_minimal()

## Peak dates
alpha_peak <- as.Date("2020-11-19")
delta_peak <- as.Date("2021-09-04")

## Alpha wave convergence analysis
## Alpha Peak Movie
j <- which(weeks == alpha_peak)

sd_values_alpha <- sd_values |> 
  filter(weeks %in% weeks[seq(j-90,j+60, 2)])

ggplot(data = sd_values_alpha, 
       aes(x = weeks, y = sd, 
           label = id,
           # ymin = lower, ymax = upper,
           color = if_else(sd <= 0.01 & sd >= 0.0025, "blue", "red"),
       ))+
  geom_label()+
  # geom_pointrange()+
  theme_minimal()

## Alpha wave convergence analysis
## Alpha Peak Movie
j <- which(weeks == delta_peak)

sd_values_delta <- sd_values |> 
  filter(weeks %in% weeks[seq(j-90,j+60, 2)])

ggplot(data = sd_values_delta, 
       aes(x = weeks, y = sd, 
           label = id,
           # ymin = lower, ymax = upper,
           color = if_else(sd <= 0.01 & sd >= 0.0025, "blue", "red"),
       ))+
  geom_label()+
  # geom_pointrange()+
  theme_minimal()

## Uncomment when is not a rerun
CAR_list <- list()

alpha_peak <- as.Date("2020-11-19")
delta_peak <- as.Date("2021-09-04")

test_dates <- c(c(alpha_peak-63,
                  alpha_peak-45,
                  alpha_peak-24,
                  alpha_peak),
                c(delta_peak-63,
                  delta_peak-45,
                  delta_peak-24,
                  delta_peak))

## This fitting takes very long, it is process of prune fitting. Where it starts fitting something and narrows down to a best fit, which parameteres are defined at the while loop

CAR_list <- foreach(i = 1:length(weeks),
                    .combine = c,
                    .multicombine = TRUE) %dopar% {

# for (i in 1:length(weeks)) {
  
  ## To clean unused variables here
  # gc()
  
  # if(is.na(dates_to_rerun[i]))next
  
  current_date <- weeks[i]
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
                       model = "besag",
                       graph = "Data/data-products/hexes_adjmat.graph",
                       scale.model = TRUE,
                       # diagonal = diag.eps,
                       constr = TRUE,
                       hyper = hyper_smooth)),
        data = as.data.frame(hex_week),
        family = "gaussian",
        control.inla = control.inla(strategy = "adaptive"),
        # control.inla = control.inla(strategy = "gaussian", h = h.value, int.strategy = "eb"),
        # control.mode = control.mode(restart = TRUE),
        control.compute = compute_list,
        control.predictor = predictor_list,
        # control.fixed = list(prec.intercept = 0.1),
        num.threads = 6,  # Prevent internal threading conflicts
        verbose = T
      )
      
      # sd_values <- best_model$summary.fitted.values$sd
      
    }, error = function(e) {
      set.seed(.Random.seed)
      
      best_model <- inla(
        as.formula(infectionsPC ~ 1 + 
                     f(ID, 
                       model = "besag",
                       graph = hex_graph,
                       # diagonal = diag.eps,
                       scale.model = TRUE,
                       constr = TRUE,
                       hyper = hyper_smooth)),
        data = as.data.frame(hex_week),
        family = "gaussian",
        control.inla = control.inla(strategy = "gaussian", int.strategy = "eb"),
        # control.inla = control.inla(strategy = "gaussian", h = h.value, int.strategy = "eb"),
        # control.mode = control.mode(restart = TRUE),
        control.compute = compute_list,
        control.predictor = predictor_list,
        # control.fixed = list(prec.intercept = 0.1),
        num.threads = 6,  # Prevent internal threading conflicts
        verbose = T
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
  
  # CAR_list[[i]] <-
    list(
    cbind(hex_week, 
          best_model$summary.fitted.values |> 
            rownames_to_column(var = "INLApred"))
    )
  
  # gc()
  rm(best_model, hex_week)
  # )
}

# stopCluster(cl)

## Saving as list object
save(list = CAR_list, 
     file = "Data/data-products/CAR_list_meta30m.RDS", 
     compress = "xz", 
     compression_level = 9)

load("Data/data-products/CAR_list_meta30m.RDS")

## Turning into a df
CAR_df <- bind_rows(CAR_list)

## Saving the df, remember to change the name if the dataset is "preomicron" or "omicronera". The pattern nomenclature to files are tsa_preomicron.csv or tsa_omicronera.csv
dataset <- "meta30m_run_preomicron_daily"

vroom::vroom_write(x = CAR_df, 
                   file = paste0("data-products/tsa_", 
                                 dataset,
                                 ".csv"))

#

# fitted_values <- inla.tmarginal(exp, CAR_model$marginals.fitted.values)
color_option <- "magma"
na_color <- "grey70"
## Breakdowns of each peaks
breaks_plt <- c(0,seq(250,2500, 250))
labels_plt <- c("250< ",seq(250,2250, 250), ' >2,500')
limits_plt <- c(0,2500)

# CAR_model <- best_model
family <- "gaussian"
model <- "besag2"

## Figure2
# Source: https://en.wikipedia.org/wiki/Federal_Information_Processing_Standard_state_code
excludes = c(
  "02", "60", "03", "81", "07", "64",
  "14", "66", "84", "86", "67", "89",
  "68", "71", "76", "69", "70", "95",
  "43", "72", "74", "78", "79", "15", "11"
)

us_states <- tigris::states(cb = T) |> 
  dplyr::filter(!STATEFP %in% excludes) |> 
  tigris::shift_geometry()|> 
  st_transform(crs = 26915)

# i=325

alpha_peak <- as.Date("2020-11-19")
delta_peak <- as.Date("2021-09-04")

test_dates <- c(c(alpha_peak-63,
                  alpha_peak-45,
                  alpha_peak-24,
                  alpha_peak),
                c(delta_peak-63,
                  delta_peak-45,
                  delta_peak-24,
                  delta_peak))

index <- which(weeks %in% test_dates)

CAR_dates <- bind_rows(CAR_list[index]) |> 
  left_join(hexes) |> 
  st_as_sf()

hex_spacetime_date <- hex_spacetime |> 
  filter(date %in% test_dates)|> 
  left_join(hexes) |> 
  st_as_sf()

hexgrid_preomicron_date <- hexgrid_preomicron |> 
  st_drop_geometry() |> 
  filter(date %in% test_dates)|> 
  left_join(hexes) |> 
  st_as_sf()

hexObservations <- hexObservationsAllNoMissingGeom|> 
  select(-geometry) |> 
  filter(date %in% test_dates)|> 
  mutate(hexid = as.character(hexid)) |> 
  left_join(hexes) |> 
  st_as_sf()

## Breakdowns of each peaks
breaks_plt1 <- seq(0,500, 50)
labels_plt1 <- c("<50",seq(50,450, 50), '>500')
limits_plt1 <- c(0,500)
color_option <- "magma"

ggplot() +
  geom_sf(data = hexgrid_preomicron_date |> 
            st_transform(crs = 26915),
          aes(fill = infectionsPC))+
  scale_fill_viridis_b(option = color_option,
                       name = "Estimated Infections/100k/day",
                       direction = -1,
                       breaks = breaks_plt1,
                       labels = labels_plt1,
                       na.value = "steelblue4",
                       limits = limits_plt1
  )+
  theme_minimal()+
  theme(legend.position = "bottom",
        legend.title.position = "top",
        legend.key.width = grid::unit(1, "in"))+
  facet_wrap(.~date, nrow = 2)
# +
  # guides(fill = guide_bins(title = "Infections per capita/100k",
  #                          # labels = scales::label_math(expr = 10^., format = "force"),
  #                          title.position = "top",
  #                          title.vjust = 0.5))+
  # labs(title = "InfectionsPC", subtitle = delta_peak)+
  # labs(subtitle = paste0("CAR model with ", model," implementation and ", family, " likelihood"),
       # caption = "*Proportional weights means weights ranging from 1 to 0")

CAR_df <- bind_rows(CAR_list) |> 
  left_join(metapop |> 
              mutate(hexid = as.character(hexid))) |> 
  st_as_sf() |> 
  st_transform(crs = 26915)

# test <- cbind(metapop  |>  
#                 mutate(hexid = as.character(hexid)) |> 
#                 dplyr::select(!metapop_30m),
#               mean = best_model$summary.fitted.values$mean) |> 
#   sf::st_as_sf()

ggplot() +
  geom_sf(data = CAR_df,
          aes(fill = mean))+
  scale_fill_viridis_c()+
  # scale_fill_viridis_b(option = color_option,
  #                      name = "Estimated Infections/100k/day",
  #                      direction = -1,
  #                      breaks = scales::breaks_pretty(n = 10),
  #                      # labels = labels_plt1,
  #                      na.value = "steelblue4",
  #                      # limits = limits_plt1
  # )+
  theme_minimal()+
  theme(legend.position = "bottom",
        legend.title.position = "top",
        legend.key.width = grid::unit(1, "in"))+
  facet_wrap(.~date, nrow = 2)


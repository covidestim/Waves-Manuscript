################################################################################################
#### Example code to Håvard Rue to help me with the Besag2 and BYM2 model
#### Code created by Rafael Lopes, doubts: rafael.lopes@yale.edu
################################################################################################

## Remember to change the file paths
rm(list = ls())
gc()

## Hexes to be kept
hexid_to_keep <- sf::st_read("Data/data-sources/hexid_to_keep.geojson")

## New hexgrid and new hexgrid with Meta 30m population
hexes <- sf::st_read("Data/data-sources/hexes.geojson") |> 
  filter(as.integer(hexid) < 7662,
         ## Taking out the isolated hex at Keywest
         as.integer(hexid) != 6545,
         ## Now filtering out the high infections per capita hex
         as.character(hexid) %in% as.character(hexid_to_keep$hexid)) |> 
  st_transform(crs = 26915)

# Neighbors
hexes_nb <- spdep::poly2nb(hexes, queen = TRUE, row.names = hexes$hexid)
# hexes_nb2 <- st_touches(st_geometry(hexes))

nb2INLA("Data/data-products/hexes_adjmat.graph", nb = hexes_nb)
hexes_graph <- INLA::inla.read.graph("Data/data-products/hexes_adjmat.graph")

# # wt_B <- nb2mat(neighbours = hexes_nb, style = "B", zero.policy = T)
# wt_W <- nb2mat(neighbours = hexes_nb, style = "W", zero.policy = T)
# # 
# # listB <- nb2listw(neighbours = hexes_nb, style = "B", zero.policy = T)
# listW <- nb2listw(neighbours = hexes_nb, style = "W", zero.policy = T)
# # 
# # ## Weight list as a Sparse C Matrix
# # B <- as(listB, "CsparseMatrix")
# W <- as(listW, "CsparseMatrix")

## Hexgrid data
hex_preomicron <- vroom::vroom("Data/data-products/minimal_model/hexgrid_preomicron_testDates.csv")
dates_unique <- unique(hex_preomicron$date)

## Pre-Omicron expanded dataset
hex_spacetime <- expand.grid(hexid = as.character(unique(hexes$hexid)),
                             date = dates_unique) |> 
  left_join(hex_preomicron |> 
              mutate(hexid = as.character(hexid)) |> 
              select(hexid, date, infections, population)) |>  ## just join the desired columns, we have the geometry
  mutate(Time = as.numeric(date - min(date)) + 1,
         ID = as.numeric(hexid),
         ## Ceiling infections and having NA's where date was missing, mostly over Nebraska
         infections = case_when(infections >= population ~ population,
                                infections < population ~ infections,
                                infections > 1 ~ infections,
                                infections <= 1 ~ NA),
         infectionsPC = (infections/population)*1e5,
         logpopulation = log10(population+1)) |> 
  left_join(hexes) |> 
  st_as_sf() |> 
  st_transform(crs=26915)

## Uncomment when is not a rerun
CAR_list <- list()

test_dates <- unique(hex_spacetime$date)

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

## Hyperparameters
hyper_smooth <- list(
  prec = list(
    prior = "pc.prec",
    param = c(0.2, 0.01),  # P(SD > 0.2) = 1% (STRONG smoothing)
    initial = 4            # Start with high precision (1/exp(4) ≈ 0.018 SD)
  )
)

hyper_smooth_bym2 <- list(
  phi = list(prior = "pc", param = c(0.95, 0.5)),  # 50% prob ϕ > 0.95
  prec = list(prior = "pc.prec", param = c(0.2, 0.01))
)

## Diagonal minimal adding to avoid numerical singularity
diag.eps = 1e-3

## graph
hex_graph <- INLA::inla.read.graph("Data/data-products/hexes_adjmat.graph")

## List to hold model output
CAR_list <- list()

for (i in 1:length(test_dates)) {
  
  current_date <- test_dates[i]
  hex_week <- hex_spacetime %>% 
    filter(date == current_date)
  
  set.seed(.Random.seed)
  
  best_model <- inla(
    as.formula(infectionsPC ~ 1 + 
                 f(ID, 
                   model = "besag2",
                   graph = "Data/data-products/hexes_adjmat.graph",
                   scale.model = TRUE,
                   diagonal = diag.eps,
                   constr = TRUE,
                   hyper = hyper_smooth
                   )),
    data = as.data.frame(hex_week),
    family = "gaussian",
    control.inla = control.inla(strategy = "gaussian"),
    # control.inla = control.inla(strategy = "gaussian", h = h.value, int.strategy = "eb"),
    control.mode = control.mode(restart = TRUE),
    control.compute = compute_list,
    control.predictor = predictor_list,
    control.fixed = list(prec.intercept = 0.1),
    num.threads = 6,  # Prevent internal threading conflicts
    # verbose = T
  )
  
  cat("Finished CAR model for week ", as.character(current_date),"! \n")
  
  CAR_list[[i]] <-
    # list(
    cbind(hex_week, 
          best_model$summary.fitted.values |> 
            rownames_to_column(var = "INLApred"))
  
  # gc()
  rm(best_model, hex_week)
  # )
}

## Binding everything
CAR_df <- bind_rows(CAR_list)

# |> 
#   left_join(hexes) |> 
#   st_as_sf() |> 
#   st_transform(crs = 26915)

## Breakdowns of each peaks
# breaks_plt1 <- seq(0,500, 50)
# labels_plt1 <- c("<50",seq(50,450, 50), '>500')
# limits_plt1 <- c(0,500)
color_option <- "magma"

## Test plot
ggplot() +
  geom_sf(data = CAR_df,
          aes(fill = mean))+
  scale_fill_viridis_c(option = color_option, direction = -1)+
  theme_minimal()+
  theme(legend.position = "bottom",
        legend.title.position = "top",
        legend.key.width = grid::unit(1, "in"))+
  facet_wrap(.~date, nrow = 2)

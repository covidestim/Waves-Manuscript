rm(list = ls())
gc()

library(tidyverse)
library(sf)
library(spdep)
library(spatialreg)
library(Matrix)
library(INLA)
library(vroom)

setwd("~/Desktop/repos/Waves/Waves-Manuscript/")

# Load spacetime and hexgrid
hex_spacetime <- vroom::vroom("Data/data-products/hex_spacetime.csv") %>%
  mutate(hexid = as.character(hexid))
hexgrid     <- sf::st_read("Data/data-products/hexgrid.geojson")

# Build neighbor graph
hexes_nb      <- spdep::poly2nb(hexgrid, queen = TRUE, row.names = hexgrid$hexid)
nb2INLA("Data/data-products/hexes_adjmat.graph", nb = hexes_nb)
hexes_graph   <- INLA::inla.read.graph("Data/data-products/hexes_adjmat.graph")

# INLA control lists
compute_list <- control.compute(hyperpar = TRUE, return.marginals = TRUE,
                                 return.marginals.predictor = TRUE, config = TRUE,
                                 dic = TRUE, cpo = TRUE, waic = TRUE)
predictor_list <- control.predictor(compute = TRUE)
hyper_smooth_bym2 <- list(
  phi  = list(prior = "pc", param = c(0.95, 0.7)),
  prec = list(prior = "pc.prec", param = c(0.2, 0.01))
)

# Define aberration detection with recovery
detect_sd_aberration_with_recovery <- function(upper_sd, median_sd,
                                               median_threshold = 1.5,
                                               range_threshold  = 5,
                                               recovery_margin  = 0.85) {
  n <- length(median_sd)
  result <- rep(FALSE, n)
  in_recovery <- FALSE
  baseline_median <- NA

  for (i in 2:(n-1)) {
    # median drop
    drop_prev  <- (median_sd[i-1] - median_sd[i]) > median_threshold
    drop_next  <- (median_sd[i+1] - median_sd[i]) > median_threshold
    median_drop <- drop_prev || drop_next
    # upper spike
    spike_prev <- (upper_sd[i] - upper_sd[i-1]) > range_threshold
    spike_next <- (upper_sd[i] - upper_sd[i+1]) > range_threshold
    upper_spike <- spike_prev && spike_next

    if (median_drop || upper_spike) {
      result[i] <- TRUE
      in_recovery <- TRUE
      baseline_median <- median_sd[i-1]
    } else if (in_recovery) {
      if (median_sd[i] < baseline_median * recovery_margin) {
        result[i] <- TRUE
      } else {
        in_recovery <- FALSE
      }
    }
  }
  return(result)
}

# Load results of prior runs
CAR_df   <- vroom::vroom("Data/data-products/tsa_meta30m_run_preomicron_daily.csv")
CAR_list <- CAR_df %>% mutate(hexid = as.character(hexid)) %>% group_split(date)
weeks    <- sort(unique(CAR_df$date))

# Compute SD summary per week
sd_values <- tibble(
  weeks     = weeks,
  upper_sd  = map_dbl(CAR_list, ~ max(.x$sd, na.rm=TRUE)),
  median_sd = map_dbl(CAR_list, ~ median(.x$sd, na.rm=TRUE)),
  lower_sd  = map_dbl(CAR_list, ~ min(.x$sd, na.rm=TRUE))
)

# Detect dates to rerun
sd_values <- sd_values %>%
  mutate(aberration = detect_sd_aberration_with_recovery(upper_sd, median_sd))

dates_to_rerun <- sd_values %>% filter(aberration) %>% pull(weeks)
cat("Will rerun for:", length(dates_to_rerun), "dates!\n")

# Loop over all weeks, rerunning only flagged dates
for (current_date in weeks) {
  if (!(current_date %in% dates_to_rerun)) next
  cat("Starting model for date:", as.character(current_date), "!\n")
  hex_week <- hex_spacetime %>% filter(date == current_date)

  best_model <- NULL
  counter    <- 0
  keep_fitting <- TRUE
  while (keep_fitting && counter <= 10) {
    cat("Attempt to fit model number:", counter, "for", as.character(current_date), "\n")
    counter <- counter + 1
    res <- try(
      inla(
        infectionsPC ~ 1 + f(ID, model="bym2", graph=hexes_graph,
                              scale.model=TRUE, constr=TRUE,
                              hyper=hyper_smooth_bym2),
        data = as.data.frame(hex_week), family = "gaussian",
        control.inla = control.inla(strategy="gaussian", int.strategy = "eb"),
        control.mode = control.mode(restart = TRUE),
        control.compute = compute_list,
        control.predictor = predictor_list,
        num.threads = 2
      ), silent = TRUE)

    if (inherits(res, "try-error")) next
    best_model <- res

    # Check fit SD stats
    prev_sd   <- sd_values$upper_sd[sd_values$weeks == current_date]
    new_sd    <- range(best_model$summary.fitted.values$sd)
    new_med   <- median(best_model$summary.fitted.values$sd)
    # check primary aberration on this single date
    # (for stopping loop if model converged)
    diff_drop <- (prev_sd - new_med) > 1.5
    diff_spike <- (new_sd[2] - prev_sd) > 5
    keep_fitting <- (diff_drop || diff_spike)
  }

  cat("Finished CAR model for week", as.character(current_date), "!\n")
  CAR_list[[which(weeks == current_date)]] <-
    bind_cols(hex_week,
              best_model$summary.fitted.values %>%
                rownames_to_column("INLApred"))
  rm(best_model, hex_week)
  gc()
}

# Save outputs
CAR_df <- bind_rows(CAR_list)
vroom::vroom_write(CAR_df,
  file = "Data/data-products/tsa_meta30m_run_preomicron_daily_updated.csv")

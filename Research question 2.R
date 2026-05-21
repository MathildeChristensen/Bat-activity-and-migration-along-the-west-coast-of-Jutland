# Question 2:
# How does bat activity vary among monitoring stations, and how can these 
# differences be interpreted in relation to distance to the coast and geographic position?

  
#               PREPARING DATA

# Making sure the datasets are data.table
setDT(dat_model_total)
setDT(stationLoc)

# Keeping only the three focal species/groups used in the analysis
dat_model_subset3 = dat_model_total[species=="ENV"|species=="Pnat"|species=="Ppyg",]

# Rename day to julian_day for clarity.
# julian_day represents day of year, where 1 = January 1.
setnames(dat_model_subset3, "day", "julian_day")
setnames(dat_model_subset3, "day", "julian_day")

# The original data are hourly observations.
# Here, activity is summed to daily activity per station, species and year.
data_daily2 <- dat_model_subset3[
  , .(daily_activity_min = sum(n_act_mins, na.rm = TRUE)),
  by = .(bat_year, species, julian_day, station_id, distance_coast)]

# Checking the structure of the model dataset
str(data_daily2)

# Converting grouping variables to factors
data_daily2[, station_id := as.factor(station_id)]
data_daily2[, species := as.factor(species)]
data_daily2[, bat_year := factor(bat_year, levels = c("2023", "2024", "2025"))]


#              RESPONSE VARIABLE CHECKS

# The response variable is daily_activity_min:
# the total number of minutes with bat activity per night/day combination.

hist(data_daily2$daily_activity_min)

# The histogram is strongly right-skewed.
# This means that most observations have low activity,
# while a few observations have very high activity.

# Counting how many observations have zero activity
data_daily2[daily_activity_min == 0, .N]

# Many zeros are expected because bats are not detected every night at every station.

# Checking whether observations reach the maximum possible number of minutes.
# If a night has approximately 14 dark hours, the maximum would be 60 * 14 minutes.
data_daily2[daily_activity_min == 60 * 14, .N]

# If very few or no observations reach this maximum,
# the upper bound is unlikely to strongly affect the analysis.

# Comparing mean and standard deviation
data_daily2[, mean(daily_activity_min)]
data_daily2[, sd(daily_activity_min)]

#              MODEL FORMULA

# Modelling daily bat activity as a function of station, species,
# and seasonal timing throughout the year.

formula_dist <- daily_activity_min ~
  
  # Difference in overall activity level among stations and species
  station_id * species +
  
  # Species- and year-specific seasonal smooth.
  # The cyclic cubic spline connects day 366 and day 1,
  # allowing seasonal activity patterns to vary smoothly across the year.
  s(julian_day,
    bat_year,
    by = species,
    bs = "fs",
    xt = list(bs = "cc"))

#              NEGATIVE BINOMIAL MODEL

# bam() is used instead of gam() because it is faster for large datasets.
# The negative binomial family is used because the response variable is overdispersed count data.

model_nb1a <- bam(
  formula_dist,
  family = nb(),
  data = data_daily2,
  method = "fREML",
  discrete = TRUE,
  control = gam.control(trace = TRUE),
  knots = list(julian_day = c(1, 366)))

#              MODEL DIAGNOSTICS

# Simulating residuals to check whether the model fits the data 
sim_res2 <- simulateResiduals(model_nb1a)
plot(sim_res2)

# Testing whether overdispersion is still present after fitting the model
testDispersion(sim_res2)

# If overdispersion is no longer problematic, the negative binomial model
# is a better choice than a Poisson model.

#              OBSERVED VS PREDICTED ACTIVITY

# Compare observed daily activity with model-predicted activity.
# Points close to the 1:1 line indicate good agreement between observations and predictions.
# Points above the line indicate that the model underpredicts observed activity.
# Points below the line indicate that the model overpredicts observed activity.

data_daily2 %>% 
  as_tibble() %>% 
  mutate(
    pred = predict(model_nb1a, ., type = "response") %>% as.vector()
  ) %>% 
  ggplot(aes(pred, daily_activity_min)) +
  geom_bin2d(binwidth = 1) +
  stat_summary_bin(fun.data = mean_se) +
  geom_abline(slope = 1, intercept = 0, color = "firebrick") +
  scale_fill_viridis_c(trans = "log10")

# If high observed values lie above the line, the model underpredicts activity peaks.
# This may suggest that some high-activity nights are driven by factors not included
# in the model, such as weather conditions or migration pulses.

# Testing if there are many zeros 
testZeroInflation(model_nb1a)
#if value close to 1, then zero-inflation is not major problem

# Checking the k-selection
# If the diagnostic "k' index" is near 1 or p-value < 0.05, they k may be too small
# If the effective degrees of freedom (edf) is close to k-1, the k could be raised
k.check(model_nb1a)
#however: allowing more freedom for the day smooths will may make the lines very wiggly

# Plotting model terms
plot(model_nb1a, all.terms = T, pages = 1)

# Model summary
summary(model_nb1a)


#              PREDICTIONS BY STATION

# Get the observed combinations of species, year and station.
# This avoids predicting combinations that do not exist in the original data.

observed_combos2 <- data_daily2 |>
  distinct(species, bat_year, station_id)

# Creating a prediction grid.
# For each observed species-year-station combination,
# predictions are made across all observed julian days.

smooth_grid2 <- observed_combos2 |>
  crossing(julian_day = sort(unique(data_daily2$julian_day))) |>
  mutate(.row = row_number())

# Drawing 500 fitted samples from the model.
# These are not 500 new datasets.
# They are 500 plausible model predictions that represent uncertainty
# in the fitted model.

fit2 <- fitted_samples(
  model_nb1a,
  data = smooth_grid2,
  n = 500,
  scale = "response"
) %>%
  left_join(
    smooth_grid2 %>% select(station_id, julian_day, bat_year, species, .row),
    by = ".row")

#              SUMMARISE PREDICTED ACTIVITY BY STATION

# Summarizing predicted activity across the year and across years.
# This gives one overall predicted activity estimate per station and species.

species_act_station <- fit2 %>%
  group_by(species, station_id) %>%
  summarize(
    pred_median = median(.fitted, na.rm = TRUE),
    pred_lower  = quantile(.fitted, 0.025, na.rm = TRUE),
    pred_upper  = quantile(.fitted, 0.975, na.rm = TRUE),
    .groups = "drop")

# Pred_median = central model estimate
# Pred_lower and pred_upper = uncertainty interval around the model estimate

#              ORDER STATIONS NORTH TO SOUTH

setDT(stationLoc)
setDT(species_act_station)

# Renaming station column if needed
if ("station" %in% names(stationLoc) && !"station_id" %in% names(stationLoc)) {
  setnames(stationLoc, "station", "station_id")}

# Ordering stations by latitude from north to south
station_order <- stationLoc[
  order(-lat),
  station_id]

# Applying this order to the plotting dataset
species_act_station[
  , station := factor(station_id, levels = station_order)]

# Removing unused species levels 
species_act_station <- species_act_station[!is.na(species)]
species_act_station[, species := droplevels(species)]

#              PLOT: AVERAGE PREDICTED ACTIVITY BY STATION

# This plot compares overall predicted activity among stations.
# Stations are ordered from north to south.
# This helps visualize whether activity changes geographically along the coast.

ggplot(species_act_station,
       aes(x = station, y = pred_median)) +
  geom_pointrange(
    aes(
      ymin = pred_lower,
      ymax = pmin(pred_upper, 50)
    ),
    alpha = 0.35
  ) +
  geom_point(size = 1.8) +
  coord_cartesian(ylim = c(0, 10)) +
  facet_grid(species ~ .) +
  labs(
    x = "Stations (north → south)",
    y = "Predicted activity (min per night)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Interpretation:
# Higher predicted values indicate stations where the model estimates higher
# average bat activity.
# Because stations are ordered north to south, spatial patterns can be compared
# visually across the coastline.



# Creating seasonal prediction curves for each station.
observed_combos20 <- data_daily2 %>%
  distinct(species, bat_year, station_id)

# Creating prediction grid for each observed station-species-year combination
smooth_grid20 <- observed_combos20 %>%
  crossing(julian_day = sort(unique(data_daily2$julian_day))) %>%
  mutate(.row20 = row_number())

# Drawing 500 fitted samples to represent uncertainty in the predicted seasonal curves
fit20 <- fitted_samples(
  model_nb1a,
  data = smooth_grid20,
  n = 500,
  scale = "response"
) %>%
  left_join(
    smooth_grid20 %>%
      select(station_id, julian_day, bat_year, species, .row20),
    by = c(".row" = ".row20"))

# Summarising predictions for each station, species and day of year
plot_dat_station20 <- fit20 %>%
  group_by(station_id, species, julian_day) %>%
  summarise(
    pred_median20 = median(.fitted, na.rm = TRUE),
    pred_lower20  = quantile(.fitted, 0.025, na.rm = TRUE),
    pred_upper20  = quantile(.fitted, 0.975, na.rm = TRUE),
    .groups = "drop")

# Adding station information
plot_dat_station20 <- plot_dat_station20 %>%
  left_join(stationLoc, by = "station_id")

# Creating station order from north to south using latitude
station_order20 <- plot_dat_station20 %>%
  group_by(station_id) %>%
  summarise(lat = mean(lat, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(lat)) %>%
  pull(station_id)

# Applying north to south order
plot_dat_station20 <- plot_dat_station20 %>%
  mutate(
    station_id = factor(station_id, levels = unique(station_order20)))

#             FINAL PLOT

ggplot(plot_dat_station20,
       aes(x = julian_day,
           y = pred_median20,
           color = species,
           fill = species,
           group = species)) +
  
  geom_ribbon(
    aes(ymin = pred_lower20, ymax = pred_upper20),
    alpha = 0.2,
    color = NA
  ) +
  
  geom_line(linewidth = 0.8) +
  
  facet_wrap(~ station_id, ncol = 4) +
  
  scale_x_continuous(
    breaks = c(15, 74, 135, 196, 258, 319),
    labels = c("Jan","Mar","May","Jul","Sep","Nov")
  ) +
  
  labs(
    x = "Month",
    y = "Predicted activity (min per night)",
    color = "Species",
    fill = "Species"
  ) +
  
  theme_bw() +
  
  theme(
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  
  coord_cartesian(ylim = c(0, 300))

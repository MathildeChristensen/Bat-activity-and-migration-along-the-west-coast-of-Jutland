# Question 1: 
# What are the seasonal activity patterns of migrating bat species on the west 
# coast of Jutland?

#             REPRODUCIBILITY AND PACKAGES

# Reproducibility
set.seed(42)

# Installed packages
install.packages(c(
  "glmmTMB","lme4","mgcViz","gratia","data.table",
  "tidyverse","corrgram","DHARMa","suncalc","ggrepel",
  "patchwork","viridis"))

# Loading libraries
library(data.table)
library(tidyverse)
library(ggrepel)
library(corrgram)
library(lme4)
library(glmmTMB)
library(mgcv)
library(gratia)
library(DHARMa)
library(ggplot2)
library(patchwork)
library(viridis)
library(suncalc)
library(lubridate)

#               UPLOADING DATA

# Upload data

# Combining data
dat_model_total = rbind(dat_model_land, dat_model_land_new)

#              PREPARING VARIABLES

# Taking the data and making it into a data table
dat_model_total <- as.data.table(dat_model_total)

# Removing objects not used in the statistical model
rm(dat, dat_model_land, dat_model, summary, locations_all_buoys, 
   removed_dates, species_offshore, species, sun, colours)

# Making species as a category
dat_model_total[, species := as.factor(species)]

# Making dates 
dat_model_total[, date := as.IDate(date)]

# Making dates as a yearly number e.g. 1. jan as number 1 and +1 so days are moved 1 day forward
dat_model_total[, julian_day := yday(date)+1]


# Defining monitoring years based on calendar year
dat_model_total[, bat_year := factor(year(date), levels = c(2023, 2024, 2025), ordered = TRUE)] 

# Months are ordered to follow the biological bat season (Apr-Mar),
# while years are still defined as calendar years.

# Making month labels from the dates
dat_model_total[, month := month(date, label = TRUE)]

# Reordering months to follow the bat season (Apr-Mar)
dat_model_total[, bat_month := factor(month,levels = c("Apr","May","Jun","Jul","Aug","Sep",
                                                       "Oct","Nov","Dec","Jan","Feb","Mar"),
                                                        ordered = TRUE)]

# Making the date for the next day, because a night goes over midnight
dat_model_total[, date_end := date + 1]

# Centering time around midnight
dat_model_total[, hour_cent := ifelse(hour <= 8, hour, hour - 24)]

# Upload station data

# Matching stationdata to the big dataset
dat_model_total <- stationLoc[, .(station_id, lat, long, distance_coast)][
  dat_model_total,
  on = "station_id"]

# Renaming a station
dat_model_total <- dat_model_total %>%
  mutate(station_id = recode(
    station_id,
    "Skjern" = "Skjern-Traekfaerge"))

#              SUNRISE, SUNSET AND NIGHTDATA

# Adding sunset to each station - only getting night activity
dat_model_total[, sunset := getSunlightTimes(
  date = date,
  lat = lat[1],
  lon = long[1],
  keep = "sunset",
  tz = "UTC"
)$sunset, by = station_id]

# Adding sunrise to each station - only getting night activity
dat_model_total[, sunrise := getSunlightTimes(
  date = date_end,
  lat = lat[1],
  lon = long[1],
  keep = "sunrise",
  tz = "UTC"
)$sunrise, by = station_id]

# Only used as exploratory data
dat_model_total[, lunar_phase := getMoonIllumination(
  date = date_end,
  keep = "phase"
)$phase, by = station_id]

# Rounding sunrise and sunset hours based on minutes
dat_model_total[, rounded_sunrise := hour(sunrise)]
dat_model_total[, rounded_sunset  := hour(sunset)]

# Creating darkness column: TRUE if time is before sunrise or after sunset
dat_model_total[, darkness := hour < rounded_sunrise | hour > rounded_sunset]

# Only keeping activity in dark hours
dat_model_total= dat_model_total[darkness == T,]

# Night length in seconds
dat_model_total[, night_length :=
                  as.numeric(difftime(sunrise, sunset), units = "secs")]

# If time is after sunset the same date is used, if time is after midnight the next date is used
dat_model_total[, hour_tform := ymd_hms(
  paste0(as.Date(ifelse(hour >= hour(sunset), date, date_end)),
         sprintf(" %02d:00:00", hour)),
  tz = "UTC")]

# Calculating seconds from sunset of the observation
dat_model_total[, seconds_from_sunset := as.numeric(difftime(hour_tform, sunset, tz = "UTC"), 
                                                    units = "secs")]

# Bringing hour to a scaled level by time from sunset: 0 is sunset, 0.5 middle 
#of night, 1 is sunrise
dat_model_total[, scaled_hour := seconds_from_sunset/night_length]

#              CHECKING PLOT

# Quick check of scaling outcome
ggplot(dat_model_total, aes(y = scaled_hour, x = date)) +
  geom_point()+
  facet_grid(.~bat_year, scales = "free_x")

# Overview of all factors
summary(dat_model_total)


#             DESCRIPTIVE STATISTICS

# Total activity minutes by station
dat_model_total[, sum(n_act_mins)/60, by = station_id][order(V1)]

# Total activity minutes by per month and species (across several years)
dat_model_total[, sum(n_act_mins), by = .(bat_month, species)]


#                OVERVIEW PLOTS

# Plotting activity in the night for the months
ggplot(dat_model_total, aes(x = hour_cent, y = n_act_mins)) +
  geom_point() +
  scale_x_continuous(breaks = seq(-9, 8, 3),
                     labels = function(x) (x + 24) %% 24) +
  facet_grid(.~bat_month)+
  labs(x = "Hour of night (midnight centered)")


# Activity for year and species
ggplot(dat_model_total, aes(x = julian_day, y = n_act_mins)) +
  geom_point() +
  facet_grid(bat_year~species)+
  labs(x = "julian_day")


#                 CORRELATION ANALYSIS

# Selecting numeric variables that may be relevant as predictors
num_predictors = dat_model_total[, c("julian_day", "scaled_hour", 
                                     "mean_temp", "wind_speed","wind_direction", "precip",  
                                     "cloud_coverage", "atm_pressure", "lunar_phase")]#,

# Selecting numeric variables that may be relevant as predictors
corr_res = corrgram(num_predictors)

# Showing only correlations stronger than +/- 0.15
# Strong correlations between predictors may indicate collinearity
filtered_matrix <- ifelse(corr_res > 0.15 | corr_res < -0.15, corr_res, NA) 
filtered_matrix 

# Cleaning up temporary objects 
rm(num_predictors, corr_res, filtered_matrix)


#            DATASET FOR THE MODEL

# Keeping only the three focal species/groups used in the model
dat_model_subset = dat_model_total[species=="ENV"|species=="Pnat"|species=="Ppyg",]


# Aggregating data to sum activity minutes
data_daily <- dat_model_subset[
  , .(daily_activity_min = sum(n_act_mins, na.rm = TRUE)),
  by = .(bat_year, species, day, station_id)]

# Checking the structure of the model dataset
str(dat_model_subset)

# Making grouping variables categorical
data_daily$station_id = as.factor(data_daily$station_id)
data_daily$bat_year = factor(data_daily$bat_year, ordered = TRUE)

dat_model_subset[, bat_year := factor(bat_year, levels = c("2023", "2024", "2025"))]

# Checking the distribution of daily activity
hist(data_daily$daily_activity_min) 
# very right skewed, meaning lots of data around 0

# Number of observations with zero activity
data_daily[daily_activity_min == 0, .N] #many true zeros

# Checking whether any observations reach the maximum possible number of minutes
data_daily[daily_activity_min == 60*14, .N] 

# Comparing mean and standard deviation
# For a Poisson distribution, the mean and variance are expected to be similar.
# If the standard deviation is much larger than the mean, this indicates overdispersion.
data_daily[, mean(daily_activity_min)]
data_daily[, sd(daily_activity_min)]

# Because the response is overdispersed, a negative binomial model is expected
# to fit better than a Poisson model.

#              MODEL FORMULA

# Model daily bat activity as a function of species, seasonal timing,
# year-specific seasonal patterns, and station-level differences.

formula <- daily_activity_min ~
  
  # Difference in overall activity level between species
  species +
  
  # Species-specific seasonal smooth.
  # The cyclic cubic spline connects day 366 and day 1.
  s(day, by = species, bs = "cc") +
  
  # Year-specific deviations in the seasonal pattern
  s(day, bat_year, bs = "fs", xt = list(bs = "cc")) +
  
  # Species-specific random effect of station
  s(station_id, by = species, bs = "re")


#              POISSON MODEL

# Fitting an initial Poisson GAMM.
# bam() is used because it is faster for larger datasets.
model_poisson <- bam(
  formula,
  family = poisson(),
  data = data_daily,
  method = "fREML",
  discrete = TRUE,
  control = gam.control(trace = TRUE),
  knots = list(day = c(1, 366))
)

# Model diagnostics
sim_res <- simulateResiduals(model_poisson)
plot(sim_res)

# Testing for overdispersion
testDispersion(sim_res)

# The Poisson model showed overdispersion, so a negative binomial model was fitted.


#              NEGATIVE BINOMIAL MODEL

model_nb <- bam(
  formula,
  family = nb(),
  data = data_daily,
  method = "fREML",
  discrete = TRUE,
  control = gam.control(trace = TRUE),
  knots = list(day = c(1, 366))
)

# Model diagnostics
sim_res <- simulateResiduals(model_nb)
plot(sim_res)

# Testing whether overdispersion is still problematic
testDispersion(sim_res)

# The negative binomial model fits better and reduces the overdispersion problem.

#              OBSERVED VS PREDICTED ACTIVITY

# Comparing observed daily activity with model-predicted activity.
# If points lie close to the 1:1 line, observed and predicted values are similar.
# Points above the line indicate that the model underpredicts observed activity.
# Points below the line indicate that the model overpredicts observed activity.

data_daily %>% 
  as_tibble() %>% 
  mutate(
    pred = predict(model_nb, ., type = "response") %>% as.vector()
  ) %>% 
  ggplot(aes(pred, daily_activity_min)) +
  geom_bin2d(binwidth = 1) +
  stat_summary_bin(fun.data = mean_se) +
  geom_abline(slope = 1, intercept = 0, color = "firebrick") +
  scale_fill_viridis_c(trans = "log10")

# The model appears to underpredict some of the highest activity values.

# Testing whether there are more zeros than expected by the model
testZeroInflation(model_nb)

# A value close to 1 indicates that zero-inflation is not a major problem.


# Checking whether the basis dimension k is large enough
# If the k-index is low and the p-value is < 0.05, k may be too small.
# If the effective degrees of freedom are close to the maximum possible,
# the smooth may need more flexibility.
k.check(model_nb)

# Increasing k can make the seasonal curves more flexible,
# but may also make them biologically unrealistically wiggly.


# Plotting model smooths
plot(model_nb, all.terms = TRUE, pages = 1)

#               ALTERNATIVE MODEL FORMULA

# Adding year and species together
formula2 = daily_activity_min ~
  #main effect of species
  species +
  
  #interaction of day year and species
  s(day, bat_year, by = species, bs = "fs", xt = list(bs = "cc"))+
  
  #random effect for station by species
  s(station_id, by = species, bs = "re") 

#              ALTERNATIVE MODEL
model_nb2 = bam(formula2,
                family = nb(), 
                data = data_daily, 
                method = "fREML",
                discrete = T, # Uses a fast approximation method (good for large data)
                control = gam.control(trace = T), #for console output
                #to let out model know what the max and min value of our day cycle is
                knots = list(day = c(1, 366)))

# Comparing fits
AIC(model_nb, model_nb2)
#the second model is actually a lot better as the AIC is much smaller


# Model output
summary(model_nb2)

# Plotting the effects
plot(model_nb2, all.terms = T, pages = 1)

#                   PLOTS

# Get observed combinations
observed_combos <- data_daily |>
  distinct(species, bat_year, station_id)

# Expand only observed combinations across sequences
smooth_grid <- observed_combos |>
  crossing(day = evenly(data_daily$day)) |>
  mutate(.row = row_number())

# Drawing 500 samples from the fitted model to represent uncertainty in the
# model-predicted mean activity.
fit = fitted_samples(model_nb2, data = smooth_grid, n = 500, scale = "response") %>%
  left_join(smooth_grid %>% 
              select(station_id, day, bat_year, species, .row), by = ".row")

# Predictions are first made for all observed combinations of species, year and station.
# They are then summarized across stations to show the overall seasonal pattern
# for each species and year.
species_act_per_day_year <- fit %>%
  group_by(day, species, bat_year) %>%
  summarize(
    pred_median = median(.fitted),
    pred_lower = quantile(.fitted, 0.025),
    pred_upper = quantile(.fitted, 0.975),
    .groups = "drop")

#.           FINAL PLOT
ggplot(species_act_per_day_year,
         aes(x = day, y = pred_median, group = bat_year)) +
  
geom_line(aes(colour = bat_year), linewidth = 1.5) +
  
geom_ribbon(aes(ymin = pred_lower,
                  ymax = pred_upper,
                  fill = bat_year),
              alpha = 0.2) +
  
  facet_grid(. ~ species) +
  
  scale_x_continuous(
    breaks = c(15, 75, 135, 196, 258, 319),
    labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov")
  ) +
  
  coord_cartesian(ylim = c(0, 100)) +
  
  labs(
    x = "Month",
    y = "Predicted activity (min per night)",
    colour = "Year",
    fill = "Year"
  ) +
  
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
  
  
  

# Question 3:
# Do seasonal peaks in bat activity differ between migratory,
# foraging, and mixed habitats along the west coast of Jutland?


# Adding habitat type to the modelling dataset using station_id.
dat_model_total <- stationLoc[, .(station_id, lat, long, habitat_type)][
  dat_model_total,
  on = .(station_id)]

# Keeping only the three focal species/groups used in the analysis
dat_model_subset10 <- dat_model_total[
  species %in% c("ENV", "Pnat", "Ppyg")]

# Aggregating data information
data_daily10 <- dat_model_subset10[
  ,
  .(daily_activity_min = sum(n_act_mins, na.rm = TRUE)),
  by = .(bat_year, species, julian_day, station_id, habitat_type)]


#                PREPARING VARIABLES

# Converting group variables to factors
data_daily10[, species := factor(species, levels = c("ENV", "Pnat", "Ppyg"))]

data_daily10[, station_id := factor(station_id)]

# Making habitat variables
data_daily10[, habitat_type := factor(
  habitat_type,
  levels = c("F", "F+M", "M"))]

data_daily10[, bat_year := factor(bat_year, ordered = TRUE)]


# Defining biological seasons
# Winter is excluded due to low activity
data_daily10[, season := fifelse(
  julian_day >= 60 & julian_day < 152, "Spring",
  fifelse(
    julian_day >= 152 & julian_day < 244, "Summer",
    fifelse(
      julian_day >= 244 & julian_day <= 335, "Autumn",
      NA_character_)))]

data_daily10[, season := factor(
  season,
  levels = c("Spring", "Summer", "Autumn"))]

# Removing unused factors
data_daily10 <- droplevels(data_daily10)


# Making a combined species-habitat factor
data_daily10[, sp_hab := interaction(species, habitat_type, drop = TRUE)]

#                   MODEL FORMULA
formula_habitat10 <- daily_activity_min ~
  
  # Overall differences among species/groups
  species +
  
  # Overall differences among habitat types
  habitat_type +
  
  # Species- and habitat-specific seasonal smooths
  s(julian_day,
    by = sp_hab,
    bs = "tp",
    k = 20) +
  
  # Random effect of station
  s(station_id, bs = "re") +
  
  # Random effect of year
  s(bat_year, bs = "re")

#             NEGATIVE BINOMIAL

# Negative binomial is used because the response variable
# is overdispersed count-like data with many zeros.
model_nb10 <- bam(
  formula = formula_habitat10,
  family = nb(),
  data = data_daily10,
  method = "fREML",
  discrete = TRUE)

#              MODEL DIAGNOSTICS
# Model summary
summary(model_nb10)

# Simulating residuals to evaluate model assumptions
sim_res10 <- simulateResiduals(model_nb10)
plot(sim_res10)

# Testing overdispersion
testDispersion(sim_res10)

# Testing for zero inflation
testZeroInflation(sim_res10)

# Checking whether k is large enough
k.check(model_nb10)

# Plot to see smooth terms
plot(model_nb10, all.terms = TRUE, pages = 1)


#              CREATING PREDICTION GRID
# The prediction grid is not raw observed data.
# It is an artificial dataset used to ask the model:
# "What is the expected activity for each Julian day,
# species and habitat type?"

pred_grid10 <- expand.grid(
  julian_day = seq(
    min(data_daily10$julian_day, na.rm = TRUE),
    max(data_daily10$julian_day, na.rm = TRUE),
    by = 1
  ),
  species = levels(data_daily10$species),
  habitat_type = levels(data_daily10$habitat_type),
  station_id = levels(data_daily10$station_id)[1],
  bat_year = levels(data_daily10$bat_year)[1])


# Matching factor levels to the model data
pred_grid10$species <- factor(
  pred_grid10$species,
  levels = levels(data_daily10$species))

pred_grid10$habitat_type <- factor(
  pred_grid10$habitat_type,
  levels = levels(data_daily10$habitat_type))

pred_grid10$station_id <- factor(
  pred_grid10$station_id,
  levels = levels(data_daily10$station_id))

pred_grid10$bat_year <- factor(
  pred_grid10$bat_year,
  levels = levels(data_daily10$bat_year))

# Recreating the species × habitat interaction used in the model
pred_grid10$sp_hab <- interaction(
  pred_grid10$species,
  pred_grid10$habitat_type,
  drop = TRUE)

pred_grid10$sp_hab <- factor(
  pred_grid10$sp_hab,
  levels = levels(data_daily10$sp_hab))


#               PREDICTING ACTIVITY
# The random effects for station and year are excluded, so pred represent average
# population-level activity rather than one specific station or year
pred10 <- predict(
  model_nb10,
  newdata = pred_grid10,
  type = "link",
  se.fit = TRUE,
  exclude = c("s(station_id)", "s(bat_year)"))

# Storing fitted values and standard errors
pred_grid10$fit_link <- pred10$fit
pred_grid10$se_link  <- pred10$se.fit

# Back-transform from link scale to response scale.
# After exp(), pred is interpreted as expected activity
# in minutes per night.
pred_grid10$pred <- exp(pred_grid10$fit_link)

# lower and upper are approximate 95% confidence intervals
# for the expected activity.
pred_grid10$lower <- exp(
  pred_grid10$fit_link - 1.96 * pred_grid10$se_link)

pred_grid10$upper <- exp(
  pred_grid10$fit_link + 1.96 * pred_grid10$se_link)

# Seasons are added to the prediction data, since the model uses Julian day as a 
# continuous seasonal variable
pred_season <- pred_grid10 %>%
  mutate(
    season = case_when(
      julian_day >= 60  & julian_day < 152 ~ "Spring",
      julian_day >= 152 & julian_day < 244 ~ "Summer",
      julian_day >= 244 & julian_day <= 335 ~ "Autumn",
      TRUE ~ NA_character_
    ),
    season = factor(
      season,
      levels = c("Spring", "Summer", "Autumn")
    )
  ) %>%
  filter(!is.na(season))

#               EXTRACTING SEASONAL PEAK ACTIVITY
# For each species × habitat × season, this selects the Julian day
# with the highest predicted activity.
# These peaks are model-based estimates and not raw observed peaks

# Each row in season_peaks therefore represents the highest expected
# nightly activity within a season for one species and one habitat type.
season_peaks <- pred_season %>%
  group_by(species, habitat_type, season) %>%
  slice_max(
    order_by = pred,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  mutate(
    habitat_label = recode(
      habitat_type,
      "F"   = "Foraging",
      "F+M" = "Mixed",
      "M"   = "Migratory"
    ),
    habitat_label = factor(
      habitat_label,
      levels = c("Foraging", "Mixed", "Migratory")))


#            FINAL PLOT
# Plotting predicted peak activity by season and habitat type.
# Each facet represents one species/group.

p_peaks <- ggplot(
  season_peaks,
  aes(x = season, y = pred, fill = habitat_label)
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  facet_wrap(~ species, scales = "free_y") +
  labs(
    x = "Season",
    y = "Predicted peak activity (min/night)",
    fill = "Habitat type"
  ) +
  scale_fill_manual(
    values = c(
      "Foraging" = "bisque3",
      "Mixed" = "darkkhaki",
      "Migratory" = "darkcyan"
    )
  ) +
  theme_bw()

p_peaks



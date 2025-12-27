library(lubridate)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(readr)
library(patchwork)
library(tidyr)

urls <- c(
  "https://ourworldindata.org/grapher/temperature-anomaly.csv?v=1&csvType=full&useColumnShortNames=true",
  "https://scrippsco2.ucsd.edu/assets/data/atmospheric/stations/in_situ_co2/monthly/monthly_in_situ_co2_mlo.csv"
)

destinations <- c(
  "data/temperature-anomaly.csv",
  "data/monthly_in_situ_co2_mlo.csv"
)

# Download files
mapply(download.file, urls, destinations)

#-------------------------------------------------------------------------------
# Select required columns and rename column names
temp_anomaly <- read.csv("data/temperature-anomaly.csv")

temp_anomaly <- temp_anomaly |>
  filter(Entity == "World") |>
  select(
    Year, near_surface_temperature_anomaly,
    near_surface_temperature_anomaly_lower,
    near_surface_temperature_anomaly_upper
  )

# Show original column names
colnames(temp_anomaly)

# Replace column names
colnames(temp_anomaly) <- c(
  "year",
  "global_average_rise",
  "lower_limit",
  "upper_limit"
)

temp_anomaly
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Adjust monthly_in_situ_co2_mlo.csv dataset:
# remove header etc. and select only year and co2_seas_adj_filled

# 1. Read file as raw text
raw <- readLines("data/monthly_in_situ_co2_mlo.csv")

# 2. Find first row containing a year (4 digits)
start <- grep("^\\s*[0-9]{4},", raw)[1]

# 3. Read data starting from that row
df <- read_csv(
  "data/monthly_in_situ_co2_mlo.csv",
  skip = start - 1,
  col_names = FALSE,
  show_col_types = FALSE
)

# 4. Create new column names
names(df) <- c(
  "year", "month", "date_excel", "date_fraction",
  "co2", "co2_seas_adj", "co2_fit", "co2_fit_seas_adj",
  "co2_filled", "co2_seas_adj_filled", "station"
)

# Replace missing values (-99.99) with NA
df <- df |>
  mutate(across(where(is.numeric), ~ ifelse(.x == -99.99, NA, .x)))

# Remove NA values and select only year and co2_seas_adj_filled
co2_concentration <- df |>
  filter(!is.na(co2_seas_adj_filled)) |>
  select(year, co2_seas_adj_filled)

co2_concentration

# Normalized CO2 and temperature trend
trend <- co2_concentration |>
  select(year, co2_seas_adj_filled) |>
  left_join(temp_anomaly |> select(year, global_average_rise), by = "year") |>
  mutate(
    co2_norm = scales::rescale(co2_seas_adj_filled),
    temp_norm = scales::rescale(global_average_rise)
  )
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# save as rdata in rdas.
save(temp_anomaly, file = "rdas/temp_anomaly.rda")
save(co2_concentration, file = "rdas/co2_concentration.rda")
save(trend, file = "rdas/trend.rda")
#-------------------------------------------------------------------------------

today <- Sys.Date()
present_year <- year(today)
previous_year <- present_year - 1

airborne_farction <- 0.47

temp_present <- trend |> filter(year == present_year) |> select(global_average_rise)
avg_temp_present <- mean(temp_present$global_average_rise, na.rm = TRUE)

ppm_last_year <- trend |> filter(year == previous_year) |> select(co2_seas_adj_filled)
avg_ppm_last_year <- mean(ppm_last_year$co2_seas_adj_filled, na.rm = TRUE)

ppm_now <- trend |> filter(year == present_year) |> select(co2_seas_adj_filled)
avg_ppm_now <- mean(ppm_now$co2_seas_adj_filled, na.rm = TRUE)

net_emis <- avg_ppm_now - avg_ppm_last_year

gross_emis <- net_emis / airborne_farction

current_absorbtion <- gross_emis - net_emis
#-------------------------------------------------------------------------------
# Sigmoid model

# Sigmoid form: R(T) = R_max / (1 + exp(k * (T - T_mid)))

# Choose initial shape parameters (can be tuned later)
T_mid <- 1.7    # temperature where the curve is roughly halfway down
k     <- 10     # steepness of the transition

# Derive R_max so that R(T_now) = R_now
R_max <- current_absorbtion * (1 + exp(k * (avg_temp_present - T_mid)))

R_log <- function(T) {
  R_max / (1 + exp(k * (T - T_mid)))
}

# Temperature range
T_seq <- seq(0, 2.5, by = 0.01)
R_seq <- R_log(T_seq)

# Put your data in a dataframe (ggplot always works with dataframes)
temp_reduct <- data.frame(
  T = T_seq,
  R = R_seq
)
# save as rdata in rdas.
save(temp_reduct, file = "rdas/temp_reduct.rda")
#-------------------------------------------------------------------------------

summary <- trend |>
  group_by(year) |>
  summarise(
    mean_value = mean(co2_seas_adj_filled, na.rm = TRUE),
    mean_rise  = mean(global_average_rise, na.rm = TRUE)
  )
co2 <- summary$mean_value

# create a function that follows the global average temp in relation to co2 emission.
model_log <- lm(global_average_rise ~ log(co2_seas_adj_filled), data = trend)

b0 <- coef(model_log)[1]
b1 <- coef(model_log)[2]

global_average <- function(co2) {
  b0 + b1 * log(co2)
}
calc_temp_rise <- global_average(co2)

calc_temp_con <- data.frame(
  calc_temp_rise,
  co2
)
# save as rdata in rdas.
save(calc_temp_con, file = "rdas/calc_temp_con.rda")

tau_fast <- 4
tau_slow <- 400
C <- avg_ppm_now
w <- 0.7
years <- temp_anomaly$year
T_fast <- temp_anomaly$global_average_rise
T_slow <- temp_anomaly$global_average_rise

# 3. Equilibrium temperature from log(C)
Y <- b0 + b1 * log(C)

# 4. Update fast and slow components
Temp_fast <- T_fast + (Y - T_fast) / tau_fast
Temp_slow <- T_slow + (Y - T_slow) / tau_slow

# 5. Weighted total temperature
temp_total <- w  * Temp_fast + (1 - w) * Temp_slow

temp_relaxation <- data.frame(
  year = years,
  Total_temp = temp_total,
  Land_surface_temp  = Temp_fast,
  Deep_ocean_temp  = Temp_slow
) |>
  pivot_longer(cols = c(Total_temp, Land_surface_temp, Deep_ocean_temp),
               names_to = "component",
               values_to = "temperature")
# save as rdata in rdas.
save(temp_relaxation, file = "rdas/temp_relaxation.rda")

p1 <- ggplot() +
  # model curve
  geom_line(
    data = calc_temp_con,
    aes(x = co2, y = calc_temp_rise),
    linewidth = 0.7,
    color = "steelblue"
  ) +
  
  # real-world smooth
  geom_smooth(
    data = summary,
    aes(x = mean_value, y = mean_rise),
    linewidth = 0.7,
    color = "darkred",
    se = TRUE
  ) +
  geom_hline(yintercept = 1.2,
             color = "darkgreen", linetype = "dashed") +
  
  labs(
    x = "CO₂ concentration",
    y = "Temperature rise (°C)",
    title = "Calculated vs Real Temperature Rise"
  ) +
  theme_economist(base_size = 14)

p2 <- ggplot(temp_reduct, aes(x = T, y = R)) +
  geom_line(linewidth = 0.7) +
  geom_vline(xintercept = 1.2,
             color = "darkgreen", linetype = "dashed") +
  annotate("point",
           x = avg_temp_present,
           y = current_absorbtion,
           color = "darkgreen",
           size = 3) +
  geom_vline(xintercept = 2,
             color = "red", linetype = "dashed") +
  labs(
    x = "Temperature increase (°C above pre-industrial)",
    y = "Net CO₂ uptake (ppm/year)",
    title = "Logistic (S-shaped) model for CO₂ uptake"
  ) +
  theme_economist(base_size = 14)
p3 <- ggplot(temp_relaxation, aes(x = year, y = temperature, color = component)) +
  geom_line(size = 0.7) +
  scale_color_manual(values = c(
    Total_temp = "black",
    Land_surface_temp  = "red",
    Deep_ocean_temp  = "blue"
  )) +
  geom_hline(yintercept = 1.2,
             color = "darkgreen", linetype = "dashed") +
  labs(
    title = "Fast and slow components of total temperature",
    x = "Year",
    y = "Temperature (°C)",
    color = ""
  ) +
  theme_economist(base_size = 14)

cat("Model: global_average =", round(b0, 3), "+", round(b1, 3), "* log(CO2)\n")
cat("Temperature rise now :", avg_temp_present, "(°C)", "\n")
cat("Net CO2 emission :", net_emis, "(ppm)", "\n")
cat("Gross CO2 emission :", gross_emis, "(ppm)", "\n")
cat("Current CO2 absorption :", current_absorbtion, "(ppm)", "\n")

p1 / p2 + p3

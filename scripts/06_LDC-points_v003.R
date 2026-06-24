# Created: 2026-06-23
# Updated: 2026-06-23

# Purpose: Retain only points that include AERO data to create LDC points v003.

# horizontal_flux_total_MD is the column of interest; it is the median.
#   MN is mean, LPI is lower prediction interval, UPI is upper prediction interval,
#   STD is standard deviation.

# This AERO data is the latest available (for real this time), which is more up to date
#   than what is currently loaded into LDC and part of the downloaded geoindicators data.
#   It was generated in September 2025.


library(tidyverse)

# Load data ---------------------------------------------------------------

aero.raw <- read_csv("data/raw/all_aero_AIM_LMF_20260610.csv")

ldc.001.sj <- read_csv("data/GIS-exports/003_LDC001-Eco3-MLRA-SpatialJoin_export.csv")


# Data wrangling ----------------------------------------------------------

# Create version for joining
aero.join <- aero.raw |>
  select(PrimaryKey, horizontal_flux_total_MD)

# Join with LDC points to create v004
ldc.004 <- ldc.001.sj |>
  left_join(aero.join) |>
  filter(!is.na(horizontal_flux_total_MD))

# Arrange columns
ldc.004 <- ldc.004 |>
  select(-NA_L3CODE, -NA_L2CODE, -NA_L1CODE, -MLRA_ID, -LRRSYM, -LRR_NAME)


# Write to CSV ------------------------------------------------------------

write_csv(ldc.004,
  file = "data/versions-from-R/06_LDC-points_v003.csv",
  na = ""
)


# Created: 2026-06-18
# Updated: 2026-06-18

# Purpose: Write the LDC points to shapefile.

# From 2025-06-01 download of geoindicators.csv.


library(tidyverse)
library(sf)

# Load data ---------------------------------------------------------------

geoindicators.raw <- read_csv("data/raw/ldc-data-2026-06-01/all/geoindicators.csv")

# Data wrangling ----------------------------------------------------------

# Retain Primary Key and coordinate columns only
geoindicators <- geoindicators.raw |>
  rename(
    PrimaryKey = `Primary Key`,
    Latitude_NAD1983 = `Latitude (decimal degrees, NAD83)`,
    Longitude_NAD1983 = `Longitude (decimal degrees, NAD83)`
  ) |>
  select(PrimaryKey, Latitude_NAD1983, Longitude_NAD1983)


# Convert to shapefile ----------------------------------------------------

ldc.points <- st_as_sf(geoindicators,
  coords = c("Longitude_NAD1983", "Latitude_NAD1983"),
  crs = 4269
)

st_write(ldc.points, "data/LDC-points/02_LDC-points.shp",
  delete_layer = TRUE
)

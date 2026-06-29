# Created: 2026-06-29
# Updated: 2026-06-29

# Purpose: Drop rows/points with no SOLUS data.

# Same as 12.2.R from Project v03.


library(tidyverse)

# Load data ---------------------------------------------------------------

ldc.cetwi.solus.raw <- read_csv("data/versions-from-R/14.1_LDC-CETWI_v011.csv")


# Data wrangling ----------------------------------------------------------

ldc.cetwi.solus <- ldc.cetwi.solus.raw |> 
  filter(!is.na(sandtotal_0_cm))


# Write to CSV ------------------------------------------------------------

write_csv(ldc.cetwi.solus,
          file = "data/versions-from-R/14.2_LDC-CETWI-SOLUS_v011.csv")

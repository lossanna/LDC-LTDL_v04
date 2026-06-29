# Created: 2026-06-29
# Updated: 2026-06-29

# Purpose: Condense CETWI data into a single column. In the case of multiple CETWI values for
#   a single point, average the CETWI values. Drop rows/points with no CETWI data.

# Same as 12.1.R from Project v03.

library(tidyverse)

# Load data ---------------------------------------------------------------

ldc.cetwi.solus.raw <- read_csv("data/GIS-exports/011_LDC010-CETWI-SOLUS_export.csv")


# Data wrangling ----------------------------------------------------------

# Convert all CETWI columns to numeric
ldc.cetwi.solus <- ldc.cetwi.solus.raw |> 
  mutate(across(starts_with("CETWI_"), as.numeric))

# Count the number of CETWI values per row
ldc.cetwi.solus <- ldc.cetwi.solus |> 
  mutate(n_CETWI = rowSums(across(starts_with("CETWI_"), ~ !is.na(.)))) 

# Examine how many CETWI values per row there are
count(ldc.cetwi.solus, n_CETWI)


# Multiple CETWI values ---------------------------------------------------

# Extract rows with more than one CETWI value
multi.cetwi <- ldc.cetwi.solus |>  
  filter(n_CETWI > 1)

# Combine CETWI values into just 3 columns
multi.cetwi <- multi.cetwi |> 
  mutate(row_id = row_number()) |>  
  pivot_longer(
    cols = starts_with("CETWI_"),
    names_to = "CETWI_col",
    values_to = "value"
  ) |>
  filter(!is.na(value)) |>
  group_by(row_id) |>
  mutate(CETWI_index = row_number()) |>
  ungroup() |>
  pivot_wider(
    id_cols = row_id,
    names_from = CETWI_index,
    values_from = value,
    names_prefix = "CETWI_"
  ) |>
  right_join(multi.cetwi |> mutate(row_id = row_number()) |> select(-starts_with("CETWI_")),
             by = "row_id") |>
  select(-row_id)

# Find average for cases of multiple CETWI values
multi.cetwi <- multi.cetwi |> 
  mutate(CETWI_avg = rowMeans(across(c(CETWI_1, CETWI_2, CETWI_3)), na.rm = TRUE))

# Create version for joining
multi.cetwi.join <- multi.cetwi |> 
  select(-CETWI_1, -CETWI_2, -CETWI_3) |> 
  rename(CETWI = CETWI_avg)


# Single CETWI value ------------------------------------------------------

single.cetwi <- ldc.cetwi.solus |> 
  filter(n_CETWI == 1)

single.cetwi <- single.cetwi |> 
  mutate(CETWI = rowSums(across(starts_with("CETWI_")), na.rm = TRUE)) |> 
  select(-starts_with("CETWI_"))


# Combine -----------------------------------------------------------------

# Join the multi and single CETWI dataframes
ldc.cetwi <- bind_rows(single.cetwi, multi.cetwi.join) |> 
  arrange(LDCpointID)


# Write to CSV ------------------------------------------------------------

write_csv(ldc.cetwi,
          file = "data/versions-from-R/14.1_LDC-CETWI_v011.csv")

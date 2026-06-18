# Created: 2026-06-18
# Updated: 2026-06-18

# Purpose: Create LDC points v001, which is the most recent monitoring instance for each
#   unique point in space.

# There are still the same issues with Project v02 and v03, where there are 22 points 
#   without a date visited.

# There are 26 points with multiple DateVisited values (54 total, including multiples).
#   I took the last instance this time (rather than the first) because the dates 
#   embedded in the Primary Keys refer to when it was loaded in the database, and there's
#   a chance that things loaded later are more accurate or had corrections.
# (Project v03 had 28 points with this issue, or 62 total including multiples).


library(tidyverse)

# Load data ---------------------------------------------------------------

ldc.pts <- read_csv("data/GIS-exports/001_LDCpts_export.csv")
ldc.countoverlapping <- read_csv("data/GIS-exports/001_LDCpts-CountOverlapping_export.csv")
ldc.overlaptable <- read_csv("data/GIS-exports/001_LDCpts-OverlapTable_export.csv")
geoindicators.raw <- read_csv("data/raw/ldc-data-2026-06-01/all/geoindicators.csv")

# Join tables -------------------------------------------------------------

# Join with OverlapTable
ldc.join <- ldc.pts |>
  left_join(ldc.overlaptable)

#   Look for NAs
apply(ldc.join, 2, anyNA)

# Join with CountOverlapping
ldc.join <- ldc.join |>
  left_join(ldc.countoverlapping)

#   Look for NAs
apply(ldc.join, 2, anyNA)

# Remove unnecessary cols and rename COUNT_
ldc.join <- ldc.join |>
  select(-COUNT_FC, -ORIG_NAME) |>
  rename(ldc_count = COUNT_)


# Add DateVisited column --------------------------------------------------

# Cols from geoindicators
geoindicators.join <- geoindicators.raw |>
  rename(
    ProjectKey = `Project Key`,
    PrimaryKey = `Primary Key`,
    DateVisited = `Date Visited`
  ) |>
  select(ProjectKey, PrimaryKey, DateVisited)

# Join with LDC points
ldc.join <- ldc.join |>
  left_join(geoindicators.join)


# Extract rows with no DateVisited ----------------------------------------

# NA for DateVisited
ldc.date.na <- ldc.join |>
  filter(is.na(DateVisited))
# For now, these rows should just be deleted. (still 22 points)


# Extract rows of most recent monitoring ----------------------------------

# Extract the most recent point for LDC plots that were monitored multiple times
most.recent <- ldc.join |>
  group_by(OVERLAP_OID) |>
  filter(DateVisited == max(DateVisited)) |>
  ungroup()
length(unique(most.recent$OVERLAP_OID)) == nrow(most.recent) # FALSE
#   this means that there are some points that have the same DateVisited, so multiple rows
#     for those cases are created

# Separate out points where there is only one most recent date for DateVisited
#   these ones are fine and don't need fixing
most.recent.single <- most.recent |>
  group_by(OVERLAP_OID) |>
  filter(n() == 1) |>
  ungroup()


## Multiple points/rows for DateVisited ------------------------------------

# Separate out points where DateVisited is the same for multiple rows
most.recent.multiple <- most.recent |>
  group_by(OVERLAP_OID) |>
  filter(n() > 1) |>
  ungroup() |>
  arrange(OVERLAP_OID)

# Retain only the last instance of duplicate rows
#   The date in the Primary Key is the date when the point was loaded into the database,
#     so I am going to assume that later entries might be more accurate.
most.recent.multiple.fixed <- most.recent.multiple |>
  group_by(OVERLAP_OID) |>
  slice_tail(n = 1) |>
  ungroup()


## Combine all with corrections -------------------------------------------

# Remove rows with NA for DateVisited and use correction when there are multiple most recent rows
most.recent.combined <- most.recent.single |>
  bind_rows(most.recent.multiple.fixed) |>
  filter(!is.na(DateVisited)) |>
  select(ORIG_OID, OVERLAP_OID, ProjectKey, PrimaryKey, DateVisited)


# Write LDC001 to CSV -----------------------------------------------------

# All columns
write_csv(most.recent.combined,
  file = "data/versions-from-R/03_LDC-points_v001.csv"
)


save.image("RData/03_LDC-points_v001.RData")

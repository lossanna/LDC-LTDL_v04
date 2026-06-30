# Created: 2026-06-30
# Updated: 2026-06-30

# Purpose: Collate Landscape Data Commons geospecies data from the   
#   four batches of downloads into single table, and write new CSV.


library(tidyverse)

# Load data ---------------------------------------------------------------

# LDC v011
ldc.011 <- read_csv("data/versions-from-R/14.4_LDC-points_v011.csv")

# All
geoindicators.all <- read_csv("data/versions-from-R/14.3_geoindicators.csv")

# Batch 1
geoindicators1 <- read_csv("data/raw/ldc-data-2026-06-01/ldc-lossanna-dot-nmsu-at-gmail-dot-com-20260601-180018/geoindicators.csv")
geospecies1 <- read_csv("data/raw/ldc-data-2026-06-01/ldc-lossanna-dot-nmsu-at-gmail-dot-com-20260601-180018/geospecies.csv")

# Batch 2
geoindicators2 <- read_csv("data/raw/ldc-data-2026-06-01/ldc-lossanna-dot-nmsu-at-gmail-dot-com-20260601-180224/geoindicators.csv")
geospecies2 <- read_csv("data/raw/ldc-data-2026-06-01/ldc-lossanna-dot-nmsu-at-gmail-dot-com-20260601-180224/geospecies.csv")

# Batch 3
geoindicators3 <- read_csv("data/raw/ldc-data-2026-06-01/ldc-lossanna-dot-nmsu-at-gmail-dot-com-20260601-180331/geoindicators.csv")
geospecies3 <- read_csv("data/raw/ldc-data-2026-06-01/ldc-lossanna-dot-nmsu-at-gmail-dot-com-20260601-180331/geospecies.csv")

# Batch 4
geoindicators4 <- read_csv("data/raw/ldc-data-2026-06-01/ldc-lossanna-dot-nmsu-at-gmail-dot-com-20260601-180523/geoindicators.csv")
geospecies4 <- read_csv("data/raw/ldc-data-2026-06-01/ldc-lossanna-dot-nmsu-at-gmail-dot-com-20260601-180523/geospecies.csv")


# Data wrangling ----------------------------------------------------------

# Combine geoindicators to check for row count (1 row per plot)
geoindicators <- bind_rows(geoindicators1, geoindicators2, geoindicators3, geoindicators4) %>% 
  distinct(.keep_all = TRUE)
nrow(geoindicators) == nrow(geoindicators.all) # all 62,441 plots included

# Combine geospecies
geospecies <- bind_rows(geospecies1, geospecies2, geospecies3, geospecies4) %>% 
  distinct(.keep_all = TRUE)


# # Check for missing primary keys from LDC points v011
setdiff(ldc.011$PrimaryKey, geospecies$`Primary Key`)



# Rename geospecies cols --------------------------------------------------

geospecies <- geospecies |> 
  rename(
    ProjectKey = `Project Key`,
    PrimaryKey = `Primary Key`,
    DateVisited = `Date Visited`,
    ScientificName = `Scientific Name`,
    SpeciesCover_AH = `AH Species Cover`,
    SpeciesCover_AH_n = `AH Species Cover Count`,
    MeanSpeciesHgt = `Mean Species Height (cm)`,
    MeanSpeciesHgt_n = `Mean Species Height Count (n)`,
    Woody = `Growth Habitat`,
    Lifeform = `Growth Habitat Subcategory`,
    SpeciesKey = `Species Key`,
    DatabaseKey = `Database Key`,
    DateLoaded = `Date Loaded in Database`
  )



# Write to CSV ------------------------------------------------------------

write_csv(geospecies,
          file = "data/versions-from-R/16_geospecies.csv",
          na = "")

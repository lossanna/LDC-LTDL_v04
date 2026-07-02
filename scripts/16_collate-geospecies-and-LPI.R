# Created: 2026-06-30
# Updated: 2026-07-02

# Purpose: Collate Landscape Data Commons geospecies and LPI data from the   
#   four batches of downloads into single table, and write new CSV.


library(tidyverse)

# Load data ---------------------------------------------------------------

# All
geoindicators.all <- read_csv("data/versions-from-R/14.3_geoindicators.csv")

# Batch 1
geoindicators1 <- read_csv("data/raw/ldc-data-2026-06-01/ldc-lossanna-dot-nmsu-at-gmail-dot-com-20260601-180018/geoindicators.csv")
geospecies1 <- read_csv("data/raw/ldc-data-2026-06-01/ldc-lossanna-dot-nmsu-at-gmail-dot-com-20260601-180018/geospecies.csv")
lpi1 <- read_csv("data/raw/ldc-data-2026-06-01/ldc-lossanna-dot-nmsu-at-gmail-dot-com-20260601-180018/data-lpi.csv")

# Batch 2
geoindicators2 <- read_csv("data/raw/ldc-data-2026-06-01/ldc-lossanna-dot-nmsu-at-gmail-dot-com-20260601-180224/geoindicators.csv")
geospecies2 <- read_csv("data/raw/ldc-data-2026-06-01/ldc-lossanna-dot-nmsu-at-gmail-dot-com-20260601-180224/geospecies.csv")
lpi2 <- read_csv("data/raw/ldc-data-2026-06-01/ldc-lossanna-dot-nmsu-at-gmail-dot-com-20260601-180224/data-lpi.csv")

# Batch 3
geoindicators3 <- read_csv("data/raw/ldc-data-2026-06-01/ldc-lossanna-dot-nmsu-at-gmail-dot-com-20260601-180331/geoindicators.csv")
geospecies3 <- read_csv("data/raw/ldc-data-2026-06-01/ldc-lossanna-dot-nmsu-at-gmail-dot-com-20260601-180331/geospecies.csv")
lpi3 <- read_csv("data/raw/ldc-data-2026-06-01/ldc-lossanna-dot-nmsu-at-gmail-dot-com-20260601-180331/data-lpi.csv")

# Batch 4
geoindicators4 <- read_csv("data/raw/ldc-data-2026-06-01/ldc-lossanna-dot-nmsu-at-gmail-dot-com-20260601-180523/geoindicators.csv")
geospecies4 <- read_csv("data/raw/ldc-data-2026-06-01/ldc-lossanna-dot-nmsu-at-gmail-dot-com-20260601-180523/geospecies.csv")
lpi4 <- read_csv("data/raw/ldc-data-2026-06-01/ldc-lossanna-dot-nmsu-at-gmail-dot-com-20260601-180523/data-lpi.csv")


# Data wrangling ----------------------------------------------------------

# Combine geoindicators to check for row count (1 row per plot)
geoindicators <- bind_rows(geoindicators1, geoindicators2, geoindicators3, geoindicators4) |> 
  distinct(.keep_all = TRUE)
nrow(geoindicators) == nrow(geoindicators.all) # all 62,441 plots included

# Combine geospecies
geospecies <- bind_rows(geospecies1, geospecies2, geospecies3, geospecies4) |>  
  distinct(.keep_all = TRUE)

# Combine LPI
lpi <- bind_rows(lpi1, lpi2, lpi3, lpi4) |> 
  distinct(.keep_all = TRUE)

# Check for missing primary keys
setdiff(geospecies$`Primary Key`, geoindicators.all$PrimaryKey)
setdiff(lpi$`Primary Key`, geoindicators.all$PrimaryKey)



# Rename geospecies cols --------------------------------------------------

colnames(geospecies)
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



# Rename LPI columns ------------------------------------------------------

colnames(lpi)
lpi <- lpi |> 
  rename(
    ProjectKey = `Project Key`,
    PrimaryKey = `Primary Key`,
    LineKey = `Line Key`,
    RecordKey = `Record Key`,
    ShrubShape = `Shrub Shape`,
    FormType = `Form Type`,
    FormDate = `Form Date`,
    DateVisited = `Date Visited`,
    LineLength = `Line Length`,
    SpacingIntervalAmt = `Spacing Interval Amount`,
    SpacingType = `Spacing Type`,
    ShowCheckbox = `Show Checkbox`,
    CheckboxLabel = `Checkbox Label`,
    PointLocation = `Point Location`,
    PointNumber = `Point Number`,
    DatabaseKey = `Database Key`,
    DateLoaded = `Date Loaded in Database`
  )

# Write to CSV ------------------------------------------------------------

# geospecies
write_csv(geospecies,
          file = "data/versions-from-R/16_geospecies.csv",
          na = "")

# LPI
write_csv(lpi,
          file = "data/versions-from-R/16_data-lpi.csv",
          na = "")

# Created: 2026-06-29
# Updated: 2026-06-30

# Purpose: Gather veg cover columns needed for propensity score matching (PSM)
#   from LDC geoindicators data.

# Same as 10.R script from Project v03, except additional geoindicators.csv table is 
#   written out.

# To install terradactyl: remotes::install_github(repo = 'Landscape-Data-Commons/terradactyl')


library(tidyverse)
library(terradactyl)

# Load data ---------------------------------------------------------------

ldc.010.raw <- read_csv("data/versions-from-R/13_LDC-points_v010.csv")
geoindicators.raw <- read_csv("data/raw/ldc-data-2026-06-01/all/geoindicators.csv")


# Data wrangling ----------------------------------------------------------

# Adjust column names of geoindicators
col_rename_map <- c(
  "Project Key" = "ProjectKey",
  "Primary Key" = "PrimaryKey",
  "Date Visited" = "DateVisited",
  "Ecological Site ID" = "EcoSiteID",
  "Latitude (decimal degrees, NAD83)" = "Latitude",
  "Longitude (decimal degrees, NAD83)" = "Longitude",
  "Location Status" = "LocationStatus",
  "Location Type" = "LocationType",
  "Latitude, Actual (decimal degrees, NAD83)" = "LatActual",
  "Longitude, Actual (decimal degrees, NAD83)" = "LonActual",
  "Bare Soil (% First Hit)" = "BareSoil_FH",
  "Total Foliar Cover (%)" = "TotalFoliarCover",
  "Annual Forb Cover (% Any Hit)" = "AnnForbCover_AH",
  "Annual Graminoid Cover (% Any Hit)" = "AnnGramCover_AH",
  "Forb Cover (% Any Hit)" = "ForbCover_AH",
  "Annual Forb and Graminoid Cover (% Any Hit)" = "AnnForbGramCover_AH",
  "Graminoid Cover (% Any Hit)" = "GramCover_AH",
  "Perennial Forb Cover (% Any Hit)" = "PerForbCover_AH",
  "Perennial Forb and Graminoid Cover (% Any Hit)" = "PerForbGramCover_AH",
  "Perennial Graminoid Cover (% Any Hit)" = "PerGramCover_AH",
  "Shrub Cover (% Any Hit)" = "ShrubCover_AH",
  "FH Cyanobacteria Cover (% First Hit)" = "CyanobacteriaCover_FH",
  "Deposited Soil Cover (% First Hit)" = "DepositedSoilCover_FH",
  "Duff Cover (% First Hit)" = "DuffCover_FH",
  "Embedded Litter Cover (% First Hit)" = "EmbeddedLitterCover_FH",
  "Herbaceous Litter Cover (% First Hit)" = "HerbLitterCover_FH",
  "Lichen Cover (% First Hit)" = "LichenCover_FH",
  "Moss Cover (% First Hit)" = "MossCover_FH",
  "Rock Cover (% First Hit)" = "RockCover_FH",
  "Total Litter Cover (% First Hit)" = "TotalLitterCover_FH",
  "Vagrant Lichen Cover (% First Hit)" = "VagrantLichenCover_FH",
  "Water Cover (% First Hit)" = "WaterCover_FH",
  "Woody Litter Cover (% First Hit)" = "WoodyLitterCover_FH",
  "Canopy Gaps 25 - 50 cm (%)" = "Gap25_50",
  "Canopy Gaps 51-100 cm (%)" = "Gap51_100",
  "Canopy Gaps 101 - 200 cm (%)" = "Gap101_200",
  "Canopy Gaps > 200 cm (%)" = "Gap200plus",
  "Canopy Gaps > 25 cm (%)" = "Gap25plus",
  "Mean Forb Height (cm)" = "MeanForbHgt",
  "Mean Graminoid Height (cm)" = "MeanGramHgt",
  "Mean Herbaceous Plant Height (cm)" = "MeanHerbHgt",
  "Mean Perennial Forb Height (cm)" = "MeanPerForbHgt",
  "Mean Perennial Forb Graminoid Height (cm)" = "MeanPFbGrHgt",
  "Mean Perennial Graminoid Height (cm)" = "MeanPerGramHgt",
  "Mean Woody Plant Height (cm)" = "MeanWoodyHgt",
  "Total Annual Production (Rangeland Health)" = "TotAnnualProduction_RH",
  "Bare Ground (Rangeland Health)" = "BareGround_RH",
  "Biotic Integrity (Rangeland Health)" = "BioticIntegrity_RH",
  "Comments: Biotic Integrity (Rangeland Health)" = "BioticIntegrity_comments",
  "Comments: Hydrologic Function (Rangeland Health)" = "HydrologicFunction_comments",
  "Comments: Soil and Site Stability (Rangeland Health)" = "SoilSiteStability_comments",
  "Compaction (Rangeland Health)" = "Compaction_RH",
  "Proportion of Dead or Dying Plant Parts (Rangeland Health)" = "PropDeadDyingPlants_RH",
  "Functional/Sructural Groups (Rangeland Health)" = "FunctionalStructuralGroups_RH",
  "Gullies (Rangeland Health)" = "Gullies_RH",
  "Hydrologic Function (Rangeland Health)" = "HydrologicFunction_RH",
  "Invasive Plants (Rangeland Health)" = "InvasivePlants_RH",
  "Litter Amount (Rangeland Health)" = "LitterAmount_RH",
  "Litter Movement (Rangeland Health)" = "LitterMovement_RH",
  "Pedestals/Terracettes (Rangeland Health)" = "Pedestals_RH",
  "Plant Community Composition (Rangeland Health)" = "PlantCommunityComposition_RH",
  "Perennial Reproductive Capability (Rangeland Health)" = "PerReproCapactiy_RH",
  "Rills (Rangeland Health)" = "Rills_RH",
  "Soil Site Stability (Rangeland Health)" = "SoilSiteStability_RH",
  "Soil Surface Loss/Degradation (Rangeland Health)" = "SoilSurfaceLoss_RH",
  "Soil Surface Erosion Resistance (Rangeland Health)" = "SoilErosionResistance_RH",
  "Water Flow Patterns (Rangeland Health)" = "WaterFlowPatterns_RH",
  "Wind Scoured Areas (Rangeland Health)" = "WindScouredAreas_RH",
  "Mean Soil Stability: Surface" = "MeanSoilStability_Surface",
  "Mean Soil Stability: Protected Samples" = "MeanSoilStability_Protected",
  "Mean Soil Stability: Unprotected Samples" = "MeanSoilStabilityUnprotected",
  "MLRA Description" = "MLRADesc_LDC",
  "MLRA Symbol" = "MLRASym_LDC",
  "Ecoregion Level I" = "EcoLvl1_LDC",
  "Ecoregion Level II" = "EcoLvl2_LDC",
  "Ecoregion Level III" = "EcoLvl3_LDC",
  "Ecoregion Level IV" = "EcoLvl4_LDC",
  "State" = "State",
  "MODIS IGBP Name" = "MODISName",
  "Database Key" = "DBKey",
  "Date Loaded in Database" = "DateLoad",
  "Total Horizontal Flux" = "TotalHorizontalFlux",
  "Total Vertical Flux" = "TotalVerticalFlux",
  "PM 2.5 Vertical Flux" = "PM25Flux",
  "PM 10 Vertical Flux" = "PM10Flux",
  "Long-Term Mean Precipitation" = "LongTermMeanPrecip",
  "Long-Term Mean Runoff" = "LongTermMeanRunoff",
  "Long-Term Mean Sediment Yield" = "LongTermMeanSedimentYield",
  "Long-Term Mean Soil Loss" = "LongTermMeanSoilLoss"
)

geoindicators <- geoindicators.raw |>
  rename(!!!setNames(names(col_rename_map), col_rename_map))

# Identify completely empty columns
empty_cols <- geoindicators |>
  summarise(across(everything(), ~ all(is.na(.)))) |>
  pivot_longer(everything(), names_to = "column", values_to = "is_empty") |>
  filter(is_empty) |>
  pull(column)

geoindicators <- geoindicators |>
  select(-all_of(empty_cols))


# Columns of interest -----------------------------------------------------

# Gather columns of interest
psm.col <- geoindicators |>
  select(
    PrimaryKey, BareSoil_FH, TotalFoliarCover,
    ForbCover_AH, GramCover_AH, ShrubCover_AH,
    AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH,
    Gap101_200, Gap200plus
  )


# Create column of gaps >100
psm.col <- psm.col |>
  mutate(Gap100plus = Gap101_200 + Gap200plus)
summary(psm.col$Gap100plus)


# Combine -----------------------------------------------------------------

ldc.010.psm <- ldc.010.raw |>
  left_join(psm.col)


# Write to CSV ------------------------------------------------------------

# LDC points
write_csv(ldc.010.psm,
  file = "data/versions-from-R/14.3_LDC-with-PSM-cols_v010.csv",
  na = ""
)

# geoindicators with corrected col names
write_csv(geoindicators,
          file = "data/versions-from-R/14.3_geoindicators.csv",
          na = "")


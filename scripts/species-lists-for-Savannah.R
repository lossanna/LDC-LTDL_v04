# Created: 2026-08-31
# Updated: 2026-08-31

# Purpose: Compile species lists for 5 plots that Savannah requested:
#   - 13052016230692112013-09-01
#   - 19051816540882232019-09-01
#   - ID_BUFO_LUP_2021_NAT-038_V12023-09-01
#   - ID_BUFO_LUP_2021_NAT_018_V12022-09-01
#   - ID_USFO_LUP_2023_USFO2023-009_V12023-09-01

library(tidyverse)

# Load data ---------------------------------------------------------------

all.matched.cover <- read_csv("data/versions-from-R/17_all-matched-data-with-species-cover_v012.csv")
crosswalk <- read_csv("data/versions-from-R/17_species-crosswalk.csv")


# Compile species lists ---------------------------------------------------

# Primary keys
primary.key <- c("13052016230692112013-09-01", "19051816540882232019-09-01",
                 "ID_BUFO_LUP_2021_NAT-038_V12023-09-01", "ID_BUFO_LUP_2021_NAT_018_V12022-09-01",
                 "ID_USFO_LUP_2023_USFO2023-009_V12023-09-01")

# Narrow matched data to primary keys
dat <- all.matched.cover |> 
  filter(PrimaryKey %in% primary.key) |> 
  select(PrimaryKey, trt_control, CurrentPLANTSCode)

# Create list of species info
species.info <- crosswalk |> 
  select(CurrentPLANTSCode, Scientific, GrowthHabit, GrowthHabitSub, Duration, Nonnative) |> 
  distinct(.keep_all = TRUE) |> 
  arrange(CurrentPLANTSCode)

# Narrow down to species in plots of interest
species.select <- species.info |> 
  filter(CurrentPLANTSCode %in% dat$CurrentPLANTSCode)

#   Remove duplicate rows
species.select <- species.select |> 
  filter(!Scientific %in% c("Artemisia tridentata subsp. tridentata",
                            "Artemisia tridentata subsp. wyomingensis",
                            "Anisantha tectorum",
                            "Bromus tectorum var. glabratus",
                            "Ranunculus testiculatus",
                            "Delphinium sonnei",
                            "Chrysothamnus nauseosus",
                            "Elymus cinereus",
                            "Leptodactylon pungens",
                            "Agropyron smithii",
                            "Elymus smithii",
                            "Poa canbyi",
                            "Poa juncifolia",
                            "Agropyron spicatum",
                            "Festuca octoflora"))

# Add species info to dat
dat <- dat |> 
  left_join(species.select) |> 
  arrange(CurrentPLANTSCode) |> 
  arrange(trt_control) |> 
  arrange(PrimaryKey)



# Write to CSV ------------------------------------------------------------

write_csv(dat,
          file = "data/data-wrangling-intermediate/species-list-for-Savannah.csv")

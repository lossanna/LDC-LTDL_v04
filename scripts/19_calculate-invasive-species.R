# Created: 2026-07-08
# Updated: 2026-07-09

# Purpose: Identify prominent invasive species across ecoregions.


library(tidyverse)

# Load data ---------------------------------------------------------------

all.matched.cover <- read_csv("data/versions-from-R/17_all-matched-data-with-species-cover_v012.csv")
all.matched <- read_csv("data/versions-from-R/17_all-matched-data_v012.csv")
ldc.011.raw <- read_csv("data/versions-from-R/14.4_LDC-points_v011.csv")
geospecies <- read_csv("data/versions-from-R/16_geospecies.csv")
crosswalk <- read_csv("data/versions-from-R/17_species-crosswalk.csv")


# Prepare LDC and geospecies data -----------------------------------------

# LDc version for joining
ldc.011 <- ldc.011.raw |> 
  select(EcoLvl3, Category, PrimaryKey) 

# geospecies version for joining
geospecies.join <- geospecies |> 
  filter(!is.na(SpeciesCover_AH)) |> 
  select(PrimaryKey, Species, ScientificName, SpeciesCover_AH, SpeciesCover_AH_n)


# NW Forested Mts / Western Cordillera ------------------------------------

## Blue Mountains ---------------------------------------------------------

# Species found in ecoregion
blue.mts.species <- ldc.011 |> 
  filter(EcoLvl3 == "Blue Mountains") |> 
  left_join(geospecies.join)

# Most abundant species
blue.mts.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) # BRTE


## Middle Rockies ---------------------------------------------------------

# Species found in ecoregion
m.rockies.species <- ldc.011 |> 
  filter(EcoLvl3 == "Middle Rockies") |> 
  left_join(geospecies.join)

# Most abundant species
m.rockies.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) |> 
  print(n = 20) # BRTE



## Southern Rockies -------------------------------------------------------

# Species found in ecoregion
s.rockies.species <- ldc.011 |> 
  filter(EcoLvl3 == "Southern Rockies") |> 
  left_join(geospecies.join)

# Most abundant species
s.rockies.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) |> 
  print(n = 30) # BRTE





# Great Plains / West-Central Semiarid Prairies ---------------------------

## Northwestern Great Plains ----------------------------------------------

# Species found in ecoregion
ngp.species <- ldc.011 |> 
  filter(EcoLvl3 == "Northwestern Great Plains") |> 
  left_join(geospecies.join)

# Most abundant species
ngp.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) # BRTE



# Cold Deserts ------------------------------------------------------------

## Snake River Plain ------------------------------------------------------

# Species found in ecoregion
srp.species <- ldc.011 |> 
  filter(EcoLvl3 == "Snake River Plain") |> 
  left_join(geospecies.join)

# Most abundant species
srp.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) # BRTE


## Northern Basin and Range -----------------------------------------------

# Species found in ecoregion
nbr.species <- ldc.011 |> 
  filter(EcoLvl3 == "Northern Basin and Range") |> 
  left_join(geospecies.join)

# Most abundant species
nbr.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) # BRTE


## Central Basin and Range ------------------------------------------------

# Species found in ecoregion
cbr.species <- ldc.011 |> 
  filter(EcoLvl3 == "Central Basin and Range") |> 
  left_join(geospecies.join)

# Most abundant species
cbr.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) # BRTE


## Wyoming Basin ----------------------------------------------------------

# Species found in ecoregion
wyb.species <- ldc.011 |> 
  filter(EcoLvl3 == "Wyoming Basin") |> 
  left_join(geospecies.join)

# Most abundant species
wyb.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) # BRTE


## Colorado Plateaus ------------------------------------------------------

# Species found in ecoregion
cop.species <- ldc.011 |> 
  filter(EcoLvl3 == "Colorado Plateaus") |> 
  left_join(geospecies.join)

# Most abundant species
cop.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) # BRTE


## Arizona/New Mexico Plateau ---------------------------------------------

# Species found in ecoregion
az.nm.plat.species <- ldc.011 |> 
  filter(EcoLvl3 == "Arizona/New Mexico Plateau") |> 
  left_join(geospecies.join)

# Most abundant species
az.nm.plat.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) |> 
  print(n = 20) # BRTE



# Warm Deserts ------------------------------------------------------------

## Mojave Basin and Range -------------------------------------------------

# Species found in ecoregion
mojave.species <- ldc.011 |> 
  filter(EcoLvl3 == "Mojave Basin and Range") |> 
  left_join(geospecies.join)

# Most abundant species
mojave.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) # BRRU2


## Chihuahuan Desert ------------------------------------------------------

# Species found in ecoregion
chihuahuan.species <- ldc.011 |> 
  filter(EcoLvl3 == "Chihuahuan Desert") |> 
  left_join(geospecies.join)

# Most abundant species
chihuahuan.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) # PRGL2 (technically native, but probably what the herbicide is for)



# Temperate Sierras / Upper Gila ------------------------------------------

## Arizona/New Mexico Mountains -------------------------------------------

# Species found in ecoregion
az.nm.mts.species <- ldc.011 |> 
  filter(EcoLvl3 == "Arizona/New Mexico Mountains") |> 
  left_join(geospecies.join)

# Most abundant species
az.nm.mts.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) # BRRU2





# List of primary keys ----------------------------------------------------

# Primary keys only
primarykeys <- all.matched |> 
  select(PrimaryKey) |> 
  distinct(.keep_all = TRUE)
primarykeys <- primarykeys$PrimaryKey  



# BRTE --------------------------------------------------------------------

# BRTE rows only
brte <- all.matched.cover |> 
  filter(CurrentPLANTSCode == "BRTE") |> 
  select(PrimaryKey, CurrentPLANTSCode, Cover_AH) |> 
  distinct(.keep_all = TRUE)

# Plots without BRTE
brte0 <- data.frame(
  PrimaryKey = primarykeys,
  CurrentPLANTSCode = "BRTE",
  Cover_AH = 0
) |> 
  filter(!PrimaryKey %in% brte$PrimaryKey)
  
# Combine
all.matched.brte <- brte |> 
  bind_rows(brte0)

# Rename column
all.matched.brte <- all.matched.brte |> 
  select(-CurrentPLANTSCode) |> 
  rename(BRTE_cover = Cover_AH)

# Join the rest of all.matched cols
all.matched.brte <- all.matched |> 
  left_join(all.matched.brte)

# Check for NAs
all.matched.brte |> 
  filter(is.na(BRTE_cover))



# BRRU2 -------------------------------------------------------------------

# BRRU2 rows only
brru2 <- all.matched.cover |> 
  filter(CurrentPLANTSCode == "BRRU2") |> 
  select(PrimaryKey, CurrentPLANTSCode, Cover_AH) |> 
  distinct(.keep_all = TRUE)

# Plots without BRRU2
brru20 <- data.frame(
  PrimaryKey = primarykeys,
  CurrentPLANTSCode = "BRRU2",
  Cover_AH = 0
) |> 
  filter(!PrimaryKey %in% brru2$PrimaryKey)

# Combine
all.matched.brru2 <- brru2 |> 
  bind_rows(brru20)

# Rename column
all.matched.brru2 <- all.matched.brru2 |> 
  select(-CurrentPLANTSCode) |> 
  rename(BRRU2_cover = Cover_AH)

# Join the rest of all.matched cols
all.matched.brru2 <- all.matched |> 
  left_join(all.matched.brru2)

# Check for NAs
all.matched.brru2 |> 
  filter(is.na(BRRU2_cover))



# Prosopis ----------------------------------------------------------------

# Prosopis codes
prosopis.codes <- crosswalk |> 
  filter(str_detect(Scientific, "Prosopis"))

# Prosopis rows only
prosopis <- all.matched.cover |> 
  filter(CurrentPLANTSCode %in% prosopis.codes$CurrentPLANTSCode) |> 
  select(PrimaryKey, CurrentPLANTSCode, Cover_AH) |> 
  distinct(.keep_all = TRUE)

# Sum to get total
prosopis <- prosopis |> 
  group_by(PrimaryKey) |> 
  summarise(Prosopis_cover = sum(Cover_AH))

# Plots without prosopis
prosopis0 <- data.frame(
  PrimaryKey = primarykeys,
  Prosopis_cover = 0
) |> 
  filter(!PrimaryKey %in% prosopis$PrimaryKey)

# Combine
all.matched.prosopis <- prosopis |> 
  bind_rows(prosopis0)

# Join the rest of all.matched cols
all.matched.prosopis <- all.matched |> 
  left_join(all.matched.prosopis)

# Check for NAs
all.matched.prosopis |> 
  filter(is.na(Prosopis_cover))




# Artemisia ---------------------------------------------------------------

# Artemisia codes
artemisia.codes <- crosswalk |> 
  filter(str_detect(Scientific, "Artemisia"))

# Artemisia rows only
artemisia <- all.matched.cover |> 
  filter(CurrentPLANTSCode %in% artemisia.codes$CurrentPLANTSCode) |> 
  select(PrimaryKey, CurrentPLANTSCode, Cover_AH) |> 
  distinct(.keep_all = TRUE)

# Sum to get total
artemisia <- artemisia |> 
  group_by(PrimaryKey) |> 
  summarise(Artemisia_cover = sum(Cover_AH))

# Plots without artemisia
artemisia0 <- data.frame(
  PrimaryKey = primarykeys,
  Artemisia_cover = 0
) |> 
  filter(!PrimaryKey %in% artemisia$PrimaryKey)

# Combine
all.matched.artemisia <- artemisia |> 
  bind_rows(artemisia0)

# Join the rest of all.matched cols
all.matched.artemisia <- all.matched |> 
  left_join(all.matched.artemisia)

# Check for NAs
all.matched.artemisia |> 
  filter(is.na(Artemisia_cover))



# Write to CSV ------------------------------------------------------------

# BRTE
write_csv(all.matched.brte,
          file = "data/versions-from-R/19_BRTE-cover_all-models_v012.csv",
          na = "")


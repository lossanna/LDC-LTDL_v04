# Created: 2026-07-08
# Updated: 2026-07-13

# Purpose: Identify prominent species across ecoregions (invasive species or targets of
#   land treatments). Find cover of these species for each LDC point/primary key.
#   Also, recalculate functional group cover for each LDC point.

# If a species/group of interest shows up in the top 30 of both raw summed cover and
#   number of plots, then it will be included in permutation tests.

# Results:
#   Blue Mountains: BRTE, Artemisia, PJ
#   Middle Rockies: BRTE, Artemisia
#   Southern Rockies: Artemisia, PJ
#   NW Great Plains: BRTE, Artemisia
#   Snake River Plain: BRTE, Artemisia
#   Northern Basin and Range: BRTE, Artemisia, PJ
#   Central Basin and Range: BRTE, Artemisia, PJ
#   Wyoming Basin: BRTE, Artemisia
#   CO Plateaus: BRTE, Artemisia, PJ
#   AZ/NM Plateau: BRTE, Artemisia, PJ
#   Mojave: BRRU2, LATR2
#   Chihuahuan: LATR2, Prosopis
#   AZ/NM Mountains: BRRU2, PJ


library(tidyverse)

# Load data ---------------------------------------------------------------

all.matched.cover <- read_csv("data/versions-from-R/17_all-matched-data-with-species-cover_v012.csv")
all.matched <- read_csv("data/versions-from-R/17_all-matched-data_v012.csv")
ldc.011.raw <- read_csv("data/versions-from-R/14.4_LDC-points_v011.csv")
geospecies <- read_csv("data/versions-from-R/16_geospecies.csv")
crosswalk <- read_csv("data/versions-from-R/17_species-crosswalk.csv")


# Prepare LDC and geospecies data -----------------------------------------

# LDC version for joining
ldc.011 <- ldc.011.raw |> 
  select(EcoLvl3, Category, PrimaryKey) 

# geospecies version for joining
geospecies.join <- geospecies |> 
  filter(!is.na(SpeciesCover_AH)) |> 
  select(PrimaryKey, Species, ScientificName, SpeciesCover_AH, SpeciesCover_AH_n)



# Identify common species -------------------------------------------------

## NW Forested Mts / Western Cordillera -----------------------------------

### Blue Mountains --------------------------------------------------------

# Species found in ecoregion
blue.mts.species <- ldc.011 |> 
  filter(EcoLvl3 == "Blue Mountains") |> 
  left_join(geospecies.join)

# Most abundant species
blue.mts.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) |> 
  print(n = 30) # BRTE, Artemisia, Juniperus

blue.mts.species |> 
  group_by(Species, ScientificName) |> 
  summarise(n = n()) |> 
  arrange(desc(n)) |> 
  print(n = 30) # BRTE, Artemisia, Juniperus


### Middle Rockies --------------------------------------------------------

# Species found in ecoregion
m.rockies.species <- ldc.011 |> 
  filter(EcoLvl3 == "Middle Rockies") |> 
  left_join(geospecies.join)

# Most abundant species
m.rockies.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) |> 
  print(n = 30) # Artemisia, BRTE, Pinus

m.rockies.species |> 
  group_by(Species, ScientificName) |> 
  summarise(n = n()) |> 
  arrange(desc(n)) |> 
  print(n = 30) # Artemisia, BRTE


### Southern Rockies ------------------------------------------------------

# Species found in ecoregion
s.rockies.species <- ldc.011 |> 
  filter(EcoLvl3 == "Southern Rockies") |> 
  left_join(geospecies.join)

# Most abundant species
s.rockies.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) |> 
  print(n = 30) # Artemisia, Pinus

s.rockies.species |> 
  group_by(Species, ScientificName) |> 
  summarise(n = n()) |> 
  arrange(desc(n)) |> 
  print(n = 30) # Artemisia, Pinus



## Great Plains / West-Central Semiarid Prairies --------------------------

### Northwestern Great Plains ---------------------------------------------

# Species found in ecoregion
ngp.species <- ldc.011 |> 
  filter(EcoLvl3 == "Northwestern Great Plains") |> 
  left_join(geospecies.join)

# Most abundant species
ngp.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) |> 
  print(n = 30) # BRTE, Artemisia, Pinus, Juniperus

ngp.species |> 
  group_by(Species, ScientificName) |> 
  summarise(n = n()) |> 
  arrange(desc(n)) |> 
  print(n = 30) # Artemisia, BRTE



## Cold Deserts -----------------------------------------------------------

### Snake River Plain -----------------------------------------------------

# Species found in ecoregion
srp.species <- ldc.011 |> 
  filter(EcoLvl3 == "Snake River Plain") |> 
  left_join(geospecies.join)

# Most abundant species
srp.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) |> 
  print(n = 30) # BRTE, Artemisia

srp.species |> 
  group_by(Species, ScientificName) |> 
  summarise(n = n()) |> 
  arrange(desc(n)) |> 
  print(n = 30) # BRTE, Artemisia


### Northern Basin and Range ----------------------------------------------

# Species found in ecoregion
nbr.species <- ldc.011 |> 
  filter(EcoLvl3 == "Northern Basin and Range") |> 
  left_join(geospecies.join)

# Most abundant species
nbr.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) |> 
  print(n = 30) # BRTE, Artemisia, Juniperus

nbr.species |> 
  group_by(Species, ScientificName) |> 
  summarise(n = n()) |> 
  arrange(desc(n)) |> 
  print(n = 30) # BRTE, Artemisia, Juniperus


### Central Basin and Range -----------------------------------------------

# Species found in ecoregion
cbr.species <- ldc.011 |> 
  filter(EcoLvl3 == "Central Basin and Range") |> 
  left_join(geospecies.join)

# Most abundant species
cbr.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) |> 
  print(n = 30) # BRTE, Artemisia, Juniperus

cbr.species |> 
  group_by(Species, ScientificName) |> 
  summarise(n = n()) |> 
  arrange(desc(n)) |> 
  print(n = 30) # BRTE, Artemisia, Juniperus
        

### Wyoming Basin ---------------------------------------------------------

# Species found in ecoregion
wyb.species <- ldc.011 |> 
  filter(EcoLvl3 == "Wyoming Basin") |> 
  left_join(geospecies.join)

# Most abundant species
wyb.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) |> 
  print(n = 30) # Artemisia, BRTE, Juniperus

wyb.species |> 
  group_by(Species, ScientificName) |> 
  summarise(n = n()) |> 
  arrange(desc(n)) |> 
  print(n = 30) # Artemisia, BRTE


### Colorado Plateaus -----------------------------------------------------

# Species found in ecoregion
cop.species <- ldc.011 |> 
  filter(EcoLvl3 == "Colorado Plateaus") |> 
  left_join(geospecies.join)

# Most abundant species
cop.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) |> 
  print(n = 30) # Juniperus, BRTE, Pinus, Artemisia

cop.species |> 
  group_by(Species, ScientificName) |> 
  summarise(n = n()) |> 
  arrange(desc(n)) |> 
  print(n = 30) # Juniperus, BRTE, Artemisia, Pinus


### Arizona/New Mexico Plateau --------------------------------------------

# Species found in ecoregion
az.nm.plat.species <- ldc.011 |> 
  filter(EcoLvl3 == "Arizona/New Mexico Plateau") |> 
  left_join(geospecies.join)

# Most abundant species
az.nm.plat.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) |> 
  print(n = 30) # Artemisia, Pinus, Juniperus, BRTE

az.nm.plat.species |> 
  group_by(Species, ScientificName) |> 
  summarise(n = n()) |> 
  arrange(desc(n)) |> 
  print(n = 30) # Artemisia, Pinus, Juniperus, BRTE



## Warm Deserts -----------------------------------------------------------

### Mojave Basin and Range ------------------------------------------------

# Species found in ecoregion
mojave.species <- ldc.011 |> 
  filter(EcoLvl3 == "Mojave Basin and Range") |> 
  left_join(geospecies.join)

# Most abundant species
mojave.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) |> 
  print(n = 30) # BRRU2, LATR12

mojave.species |> 
  group_by(Species, ScientificName) |> 
  summarise(n = n()) |> 
  arrange(desc(n)) |> 
  print(n = 30) # LATR12, BRRU2


### Chihuahuan Desert -----------------------------------------------------

# Species found in ecoregion
chihuahuan.species <- ldc.011 |> 
  filter(EcoLvl3 == "Chihuahuan Desert") |> 
  left_join(geospecies.join)

# Most abundant species
chihuahuan.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) |> 
  print(n = 30) # LATR2, Prosopis

chihuahuan.species |> 
  group_by(Species, ScientificName) |> 
  summarise(n = n()) |> 
  arrange(desc(n)) |> 
  print(n = 30) # Prosopis, LATR2



## Temperate Sierras / Upper Gila -----------------------------------------

### Arizona/New Mexico Mountains ------------------------------------------

# Species found in ecoregion
az.nm.mts.species <- ldc.011 |> 
  filter(EcoLvl3 == "Arizona/New Mexico Mountains") |> 
  left_join(geospecies.join)

# Most abundant species
az.nm.mts.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) |> 
  print(n = 30) # Pinus, Juniperus, BRRU2

az.nm.mts.species |> 
  group_by(Species, ScientificName) |> 
  summarise(n = n()) |> 
  arrange(desc(n)) |> 
  print(n = 30) # Juniperus, Pinus, BRRU2


# List of primary keys ----------------------------------------------------

# Primary keys only
primarykeys <- all.matched |> 
  select(PrimaryKey) |> 
  distinct(.keep_all = TRUE)
primarykeys <- primarykeys$PrimaryKey  



# Species of interest -----------------------------------------------------

## BRTE -------------------------------------------------------------------

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



## BRRU2 ------------------------------------------------------------------

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



## Prosopis ---------------------------------------------------------------

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

# Plots without Prosopis
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



## Artemisia --------------------------------------------------------------

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

# Plots without Artemisia
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


## Pinus & Juniperus ------------------------------------------------------

# PJ codes
pj.codes <- crosswalk |> 
  filter(str_detect(Scientific, "Pinus|Juniperus"))

# pj rows only
pj <- all.matched.cover |> 
  filter(CurrentPLANTSCode %in% pj.codes$CurrentPLANTSCode) |> 
  select(PrimaryKey, CurrentPLANTSCode, Cover_AH) |> 
  distinct(.keep_all = TRUE)

# Sum to get total
pj <- pj |> 
  group_by(PrimaryKey) |> 
  summarise(PJ_cover = sum(Cover_AH))

# Plots without pj
pj0 <- data.frame(
  PrimaryKey = primarykeys,
  PJ_cover = 0
) |> 
  filter(!PrimaryKey %in% pj$PrimaryKey)

# Combine
all.matched.pj <- pj |> 
  bind_rows(pj0)

# Join the rest of all.matched cols
all.matched.pj <- all.matched |> 
  left_join(all.matched.pj)

# Check for NAs
all.matched.pj |> 
  filter(is.na(PJ_cover))


## LATR2 ------------------------------------------------------------------

# LATR2 rows only
latr2 <- all.matched.cover |> 
  filter(CurrentPLANTSCode == "LATR2") |> 
  select(PrimaryKey, CurrentPLANTSCode, Cover_AH) |> 
  distinct(.keep_all = TRUE)

# Plots without LATR2
latr20 <- data.frame(
  PrimaryKey = primarykeys,
  CurrentPLANTSCode = "LATR2",
  Cover_AH = 0
) |> 
  filter(!PrimaryKey %in% latr2$PrimaryKey)

# Combine
all.matched.latr2 <- latr2 |> 
  bind_rows(latr20)

# Rename column
all.matched.latr2 <- all.matched.latr2 |> 
  select(-CurrentPLANTSCode) |> 
  rename(LATR2_cover = Cover_AH)

# Join the rest of all.matched cols
all.matched.latr2 <- all.matched |> 
  left_join(all.matched.latr2)

# Check for NAs
all.matched.latr2 |> 
  filter(is.na(LATR2_cover))


# Recalculate functional group cover --------------------------------------

# Note that GrowthHabitSub2 collapses Tree, Succulent, and Subshrub all into Shrub
#   category.

# Create grouping column
funct.cover <- all.matched.cover |> 
  mutate(FunctionalGroup = paste(Duration, GrowthHabitSub2)) |> 
  select(PrimaryKey, CurrentPLANTSCode, Cover_AH, FunctionalGroup) |> 
  distinct(.keep_all = TRUE)

unique(funct.cover$FunctionalGroup)


## Annual forb ------------------------------------------------------------

# Annual forb
annual.forb <- funct.cover |> 
  filter(FunctionalGroup == "Annual Forb") |> 
  group_by(PrimaryKey) |> 
  summarise(AnnForbCover_AH = sum(Cover_AH))

# Plots with 0
annual.forb0 <- data.frame(
  PrimaryKey = primarykeys,
  AnnForbCover_AH = 0
) |> 
  filter(!PrimaryKey %in% annual.forb$PrimaryKey)

# Combine
annual.forb <- bind_rows(annual.forb, annual.forb0)


## Annual grass -----------------------------------------------------------

# Annual grass
annual.gram <- funct.cover |> 
  filter(FunctionalGroup == "Annual Graminoid") |> 
  group_by(PrimaryKey) |> 
  summarise(AnnGramCover_AH = sum(Cover_AH))

# Plots with 0
annual.gram0 <- data.frame(
  PrimaryKey = primarykeys,
  AnnGramCover_AH = 0
) |> 
  filter(!PrimaryKey %in% annual.gram$PrimaryKey)

# Combine
annual.gram <- bind_rows(annual.gram, annual.gram0)


## Perennial forb ---------------------------------------------------------

# Perennial forb
perennial.forb <- funct.cover |> 
  filter(FunctionalGroup == "Perennial Forb") |> 
  group_by(PrimaryKey) |> 
  summarise(PerForbCover_AH = sum(Cover_AH))

# Plots with 0
perennial.forb0 <- data.frame(
  PrimaryKey = primarykeys,
  PerForbCover_AH = 0
) |> 
  filter(!PrimaryKey %in% perennial.forb$PrimaryKey)

# Combine
perennial.forb <- bind_rows(perennial.forb, perennial.forb0)


## Perennial grass --------------------------------------------------------

# Perennial grass
perennial.gram <- funct.cover |> 
  filter(FunctionalGroup == "Perennial Graminoid") |> 
  group_by(PrimaryKey) |> 
  summarise(PerGramCover_AH = sum(Cover_AH))

# Plots with 0
perennial.gram0 <- data.frame(
  PrimaryKey = primarykeys,
  PerGramCover_AH = 0
) |> 
  filter(!PrimaryKey %in% perennial.gram$PrimaryKey)

# Combine
perennial.gram <- bind_rows(perennial.gram, perennial.gram0)


## Shrub ------------------------------------------------------------------

# Shrub
shrub <- funct.cover |> 
  filter(FunctionalGroup == "Perennial Shrub") |> 
  group_by(PrimaryKey) |> 
  summarise(PerShrubCover_AH = sum(Cover_AH))

# Plots with 0
shrub0 <- data.frame(
  PrimaryKey = primarykeys,
  PerShrubCover_AH = 0
) |> 
  filter(!PrimaryKey %in% shrub$PrimaryKey)

# Combine
shrub <- bind_rows(shrub, shrub0)



## Combine groups ---------------------------------------------------------

# Combine groups
all.groups <- annual.forb |> 
  left_join(annual.gram) |> 
  left_join(perennial.forb) |> 
  left_join(perennial.gram) |> 
  left_join(shrub)

#   Check for NAs
apply(all.groups, 2, anyNA)

# Add in all.matched cols
all.matched.funct.cover <- all.matched |> 
  left_join(all.groups)


# Combine all -------------------------------------------------------------

# Combine
all.matched.species.funct <- all.matched.funct.cover |> 
  left_join(all.matched.brte) |> 
  left_join(all.matched.brru2) |> 
  left_join(all.matched.prosopis) |> 
  left_join(all.matched.artemisia) |> 
  left_join(all.matched.pj) |> 
  left_join(all.matched.latr2)

# Check for NAs
apply(all.matched.species.funct, 2, anyNA)

# Check for matching length with all.matched
nrow(all.matched.species.funct) == nrow(all.matched)



# Write to CSV ------------------------------------------------------------

write_csv(all.matched.species.funct,
          file = "data/versions-from-R/19_species-of-interest-and-functional-group-cover_all-models_v012.csv",
          na = "")


save.image("RData/19_calculate-species-of-interest-and-functional-group-cover.RData")

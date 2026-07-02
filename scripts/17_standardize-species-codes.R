# Created: 2026-06-30
# Updated: 2026-07-02

# Purpose: Standardize plant species names/codes.

library(tidyverse)
library(foreign)

# Load data ---------------------------------------------------------------

load("RData/15.1_matched-data.RData")
plants.db.raw <- read.dbf("data/raw/tblNationalPlants/tblNationalPlants.dbf")
geospecies <- read_csv("data/versions-from-R/16_geospecies.csv")


# Combine all matched data ------------------------------------------------

# Bind rows
all.matched <- bind_rows(mget(ls(pattern = "\\.matched$")),
                         .id = "Model") |> 
  mutate(Model = as.integer(gsub("\\D+", "", Model)))

# List of relevant species
matched.species <- geospecies |> 
  filter(PrimaryKey %in% all.matched$PrimaryKey) |> 
  select(Species, ScientificName, Duration, Woody, Lifeform, SpeciesKey, DatabaseKey) |> 
  distinct(.keep_all = TRUE)

#   Without SpeciesKey and DatabaseKey cols
ms.sans.key <- matched.species |> 
  select(-SpeciesKey, -DatabaseKey) |> 
  distinct(.keep_all = TRUE)

#   Code and scientific name only
ms.code.sci <- matched.species |> 
  select(Species, ScientificName) |> 
  distinct(.keep_all = TRUE)


# Data wrangling for National Plants DB -----------------------------------

## Compile list of all possible codes (from any column) -------------------

# NameCode is unique. Other columns (CurrentPLA, TaxonCode, and OriginalPu have
#   code-like values).

# Check if all NameCodes are included in other code columns
setdiff(plants.db.raw$CurrentPLA, plants.db.raw$NameCode)
setdiff(plants.db.raw$TaxonCode, plants.db.raw$NameCode)
setdiff(plants.db.raw$OriginalPu, plants.db.raw$NameCode)

codes.namecodes.missing1 <- plants.db.raw |> 
  filter(CurrentPLA %in% setdiff(plants.db.raw$CurrentPLA, plants.db.raw$NameCode)) |> 
  filter(!is.na(CurrentPLA))

codes.namecodes.missing2 <- plants.db.raw |> 
  filter(TaxonCode %in% setdiff(plants.db.raw$TaxonCode, plants.db.raw$NameCode)) |> 
  filter(!is.na(TaxonCode))

codes.namecodes.missing3 <- plants.db.raw |> 
  filter(OriginalPu %in% setdiff(plants.db.raw$OriginalPu, plants.db.raw$NameCode)) |> 
  filter(!is.na(OriginalPu))

# Look for instances of these missing codes in matched.species
matched.species |> 
  filter(Species %in% c(codes.namecodes.missing1$CurrentPLA,
                        codes.namecodes.missing2$TaxonCode,
                        codes.namecodes.missing3$OriginalPu)) # all matched.species codes
#                                       are included as NameCodes in plants database
  

## Narrow down to relevant codes ------------------------------------------

# Narrow down to NameCodes that appear in matched.species
plants.db <- plants.db.raw |> 
  filter(NameCode %in% matched.species$Species)


# Identify instances of multiple NameCode for the same scientific name
plants.multiple.scientific <- count(plants.db, Scientific) |> 
  arrange(desc(n)) |> 
  filter(n > 1)

# Narrow down to relevant names in matched data
multiple.scientific <- plants.db |> 
  filter(Scientific %in% intersect(matched.species$ScientificName, 
                                   plants.multiple.scientific$Scientific)) |> 
  filter(NameCode %in% matched.species$Species)



# Explore relevant species ------------------------------------------------

# Multiple rows for same Species code
ms.multiple.code.codes <- count(ms.code.sci, Species) |> 
  arrange(desc(n)) |> 
  filter(n > 1) 

ms.multiple.code <- ms.code.sci |> 
  filter(Species %in% ms.multiple.code.codes$Species) |> 
  arrange(Species) 


# Multiple codes for the same scientific name
ms.multiple.sci.codes <- count(ms.code.sci, ScientificName) |> 
  arrange(desc(n)) |> 
  filter(n > 1) |> 
  filter(!is.na(ScientificName))

ms.multiple.sci <- ms.code.sci |> 
  filter(ScientificName %in% ms.multiple.sci.codes$ScientificName) |> 
  arrange(ScientificName)



# By model ----------------------------------------------------------------

## Northwest Forested Mountains / Western Cordillera ----------------------

### 1. Blue Mountains: Herbicide ------------------------------------------

# Join with geospecies
model01.geospecies <- model01.matched |> 
  left_join(geospecies)

# Create species list
model01.species <- model01.geospecies |> 
  select(Species, ScientificName) |> 
  distinct(.keep_all = TRUE)

# Look for codes with multiple scientific names
model01.multi.code.codes <- count(model01.species, Species) |> 
  arrange(desc(n)) |> 
  filter(n > 1)

model01.multi.code <- model01.species |> 
  filter(Species %in% model01.multi.code.codes$Species) |> 
  arrange(Species) # these refer to the same species; no fix needed

# Look for scientific names with multiple codes
count(model01.species, ScientificName) |> 
  arrange(desc(n)) |> 
  filter(n > 1) # none

# Save as modelxx.geospecies.fixed
model01.geospecies.fixed <- model01.geospecies


### 2. Blue Mountains: Post-burn herbicide --------------------------------

# Join with geospecies
model02.geospecies <- model02.matched |> 
  left_join(geospecies)

# Create species list
model02.species <- model02.geospecies |> 
  select(Species, ScientificName) |> 
  distinct(.keep_all = TRUE)

# Look for codes with multiple scientific names
count(model02.species, Species) |> 
  arrange(desc(n)) |> 
  filter(n > 1) # none

# Look for scientific names with multiple codes
count(model02.species, ScientificName) |> 
  arrange(desc(n)) |> 
  filter(n > 1) # none

# Save as modelxx.geospecies.fixed
model02.geospecies.fixed <- model02.geospecies


### 3. Middle Rockies: Herbicide ------------------------------------------

# Join with geospecies
model03.geospecies <- model03.matched |> 
  left_join(geospecies)

# Create species list
model03.species <- model03.geospecies |> 
  select(Species, ScientificName) |> 
  distinct(.keep_all = TRUE)

# Look for codes with multiple scientific names
model03.multi.code.codes <- count(model03.species, Species) |> 
  arrange(desc(n)) |> 
  filter(n > 1)

model03.multi.code <- model03.species |> 
  filter(Species %in% model03.multi.code.codes$Species) |> 
  arrange(Species) # "Aster Species" vs. "Asteraceae"? Same for Astragalus, Carex,
#         Castilleja, Erigeron, Phlox
#   For all of these instances, the genus alone is a valid species, but there are also
#     lots of other species in the genus, so the "_ species" entries could be refering
#     to individuals not identified to the species level, and should therefore
#     be considered different.

#   Change code for "Aster species" to "ASTERspp" so it is different than ASTER,
#     as well as other 5 relevant cases
model03.geospecies.fixed <- model03.geospecies |> 
  mutate(Species = case_when(
    ScientificName == "Aster Species" ~ "ASTERspp",
    ScientificName == "Astragalus species" ~ "ASTRAspp",
    ScientificName == "Carex species" ~ "CAREXspp",
    ScientificName == "Erigeron species" ~ "ERIGE2spp",
    ScientificName == "Phlox species" ~ "PHLOXspp",
    TRUE ~ Species
  ))

# Look for scientific names with multiple codes
count(model03.species, ScientificName) |> 
  arrange(desc(n)) |> 
  filter(n > 1) # none



### 4. Southern Rockies: Herbicide ----------------------------------------

# Join with geospecies
model04.geospecies <- model04.matched |> 
  left_join(geospecies)

# Create species list
model04.species <- model04.geospecies |> 
  select(Species, ScientificName) |> 
  distinct(.keep_all = TRUE)

# Look for codes with multiple scientific names
model04.multi.code.codes <- count(model04.species, Species) |> 
  arrange(desc(n)) |> 
  filter(n > 1)

model04.multi.code <- model04.species |> 
  filter(Species %in% model04.multi.code.codes$Species) |> 
  arrange(Species) # these refer to the same species; no fix needed

# Look for scientific names with multiple codes
count(model04.species, ScientificName) |> 
  arrange(desc(n)) |> 
  filter(n > 1) # none

# Save as modelxx.geospecies.fixed
model04.geospecies.fixed <- model04.geospecies


### 5. Southern Rockies: Prescribed burn ----------------------------------

# Join with geospecies
model05.geospecies <- model05.matched |> 
  left_join(geospecies)

# Create species list
model05.species <- model05.geospecies |> 
  select(Species, ScientificName) |> 
  distinct(.keep_all = TRUE)

# Look for codes with multiple scientific names
model05.multi.code.codes <- count(model05.species, Species) |> 
  arrange(desc(n)) |> 
  filter(n > 1)

model05.multi.code <- model05.species |> 
  filter(Species %in% model05.multi.code.codes$Species) |> 
  arrange(Species) # these refer to the same species; no fix needed

# Look for scientific names with multiple codes
model05.multi.sci.name <- count(model05.species, ScientificName) |> 
  arrange(desc(n)) |> 
  filter(n > 1) 

model05.multi.sci <- model05.species |> 
  filter(ScientificName %in% model05.multi.sci.name$ScientificName) |> 
  filter(!is.na(ScientificName)) # these are different varieties




## Great Plains / West-Central Semiarid Prairies --------------------------

### 6. Northwestern Great Plains: Prescribed burn -------------------------

# Join with geospecies
model06.geospecies <- model06.matched |> 
  left_join(geospecies)

# Create species list
model06.species <- model06.geospecies |> 
  select(Species, ScientificName) |> 
  distinct(.keep_all = TRUE)

# Look for codes with multiple scientific names
model06.multi.code.codes <- count(model06.species, Species) |> 
  arrange(desc(n)) |> 
  filter(n > 1)

model06.multi.code <- model06.species |> 
  filter(Species %in% model06.multi.code.codes$Species) |> 
  arrange(Species) # these refer to the same species; no fix needed

# Look for scientific names with multiple codes
count(model06.species, ScientificName) |> 
  arrange(desc(n)) |> 
  filter(n > 1) # none


## Cold Deserts -----------------------------------------------------------

### 7. Snake River Plain: Post-burn aerial seeding ------------------------

# Join with geospecies
model07.geospecies <- model07.matched |> 
  left_join(geospecies)

# Create species list
model07.species <- model07.geospecies |> 
  select(Species, ScientificName) |> 
  distinct(.keep_all = TRUE)

# Look for codes with multiple scientific names
model07.multi.code.codes <- count(model07.species, Species) |> 
  arrange(desc(n)) |> 
  filter(n > 1)

model07.multi.code <- model07.species |> 
  filter(Species %in% model07.multi.code.codes$Species) |> 
  arrange(Species) # these refer to the same species; no fix needed

# Look for scientific names with multiple codes
model07.multi.sci.name <- count(model07.species, ScientificName) |> 
  arrange(desc(n)) |> 
  filter(n > 1) 

model07.multi.sci <- model07.species |> 
  filter(ScientificName %in% model07.multi.sci.name$ScientificName) |> 
  filter(!is.na(ScientificName)) |> 
  arrange(ScientificName) # half have an extra letter to denote they are perennial;
#       idk I guess I will just leave this? Possibly they could be referring to different
#       species across the ecoregion.


### 8. Snake River Plain: Post-burn aerial & drill seeding ----------------

# Join with geospecies
model08.geospecies <- model08.matched |> 
  left_join(geospecies)

# Create species list
model08.species <- model08.geospecies |> 
  select(Species, ScientificName) |> 
  distinct(.keep_all = TRUE)

# Look for codes with multiple scientific names
count(model08.species, Species) |> 
  arrange(desc(n)) |> 
  filter(n > 1) # none
 
# Look for scientific names with multiple codes
model08.multi.sci.name <- count(model08.species, ScientificName) |> 
  arrange(desc(n)) |> 
  filter(n > 1) 
 
model08.multi.sci <- model08.species |> 
  filter(ScientificName %in% model08.multi.sci.name$ScientificName) |> 
  filter(!is.na(ScientificName)) |> 
  arrange(ScientificName) # this is the same situation as model 7, and I am just leaving it
 

### 9. Snake River Plain: Post-burn closure -------------------------------
 
# Join with geospecies
model09.geospecies <- model09.matched |> 
  left_join(geospecies)
 
# Create species list
model09.species <- model09.geospecies |> 
  select(Species, ScientificName) |> 
  distinct(.keep_all = TRUE)
 
# Look for codes with multiple scientific names
count(model09.species, Species) |> 
   arrange(desc(n)) |> 
   filter(n > 1)
 
# Look for scientific names with multiple codes
model09.multi.sci.name <- count(model09.species, ScientificName) |> 
  arrange(desc(n)) |> 
  filter(n > 1) 

model09.multi.sci <- model09.species |> 
  filter(ScientificName %in% model09.multi.sci.name$ScientificName) |> 
  filter(!is.na(ScientificName)) |> 
  arrange(ScientificName) # this is the same situation as model 7, and I am just leaving it


### 10. Snake River Plain: Post-burn drill seeding ------------------------

# Join with geospecies
model10.geospecies <- model10.matched |> 
  left_join(geospecies)

# Create species list
model10.species <- model10.geospecies |> 
  select(Species, ScientificName) |> 
  distinct(.keep_all = TRUE)

# Look for codes with multiple scientific names
count(model10.species, Species) |> 
  arrange(desc(n)) |> 
  filter(n > 1)

# Look for scientific names with multiple codes
model10.multi.sci.name <- count(model10.species, ScientificName) |> 
  arrange(desc(n)) |> 
  filter(n > 1) 

model10.multi.sci <- model10.species |> 
  filter(ScientificName %in% model10.multi.sci.name$ScientificName) |> 
  filter(!is.na(ScientificName)) |> 
  arrange(ScientificName) # this is the same situation as model 7 (sans ASTRA), and I am just leaving it


### 11. Snake River Plain: Post-burn herbicide ----------------------------

# Join with geospecies
model11.geospecies <- model11.matched |> 
  left_join(geospecies)

# Create species list
model11.species <- model11.geospecies |> 
  select(Species, ScientificName) |> 
  distinct(.keep_all = TRUE)

# Look for codes with multiple scientific names
count(model11.species, Species) |> 
  arrange(desc(n)) |> 
  filter(n > 1)

# Look for scientific names with multiple codes
count(model11.species, ScientificName) |> 
  arrange(desc(n)) |> 
  filter(n > 1) # none


### 12. Northern Basin & Range: Drill seeding -----------------------------

# Join with geospecies
model12.geospecies <- model12.matched |> 
  left_join(geospecies)

# Create species list
model12.species <- model12.geospecies |> 
  select(Species, ScientificName) |> 
  distinct(.keep_all = TRUE)

# Look for codes with multiple scientific names
model12.multi.code.codes <- count(model12.species, Species) |> 
  arrange(desc(n)) |> 
  filter(n > 1)

model12.multi.code <- model12.species |> 
  filter(Species %in% model12.multi.code.codes$Species) |> 
  arrange(Species) # these refer to the same species; no fix needed

# Look for scientific names with multiple codes
model12.multi.sci.name <- count(model12.species, ScientificName) |> 
  arrange(desc(n)) |> 
  filter(n > 1) 

model12.multi.sci <- model12.species |> 
  filter(ScientificName %in% model12.multi.sci.name$ScientificName) |> 
  filter(!is.na(ScientificName)) |> 
  arrange(ScientificName) # same deal where half are marked perennial with P;
#       just leaving this for now



### 13. Northern Basin & Range: Drill seeding & soil disturbance ----------

# Join with geospecies
model13.geospecies <- model13.matched |> 
  left_join(geospecies)

# Create species list
model13.species <- model13.geospecies |> 
  select(Species, ScientificName) |> 
  distinct(.keep_all = TRUE)

# Look for codes with multiple scientific names
model13.multi.code.codes <- count(model13.species, Species) |> 
  arrange(desc(n)) |> 
  filter(n > 1)

model13.multi.code <- model13.species |> 
  filter(Species %in% model13.multi.code.codes$Species) |> 
  arrange(Species) # these refer to the same species; no fix needed

# Look for scientific names with multiple codes
model13.multi.sci.name <- count(model13.species, ScientificName) |> 
  arrange(desc(n)) |> 
  filter(n > 1) 

model13.multi.sci <- model13.species |> 
  filter(ScientificName %in% model13.multi.sci.name$ScientificName) |> 
  filter(!is.na(ScientificName)) |> 
  arrange(ScientificName) # one is marked with PF (for perennial forb, probably);
#       just leaving this for now


### 14. Northern Basin & Range: Herbicide ---------------------------------

# Join with geospecies
model14.geospecies <- model14.matched |> 
  left_join(geospecies)

# Create species list
model14.species <- model14.geospecies |> 
  select(Species, ScientificName) |> 
  distinct(.keep_all = TRUE)

# Look for codes with multiple scientific names
model14.multi.code.codes <- count(model14.species, Species) |> 
  arrange(desc(n)) |> 
  filter(n > 1)

model14.multi.code <- model14.species |> 
  filter(Species %in% model14.multi.code.codes$Species) |> 
  arrange(Species) # these refer to the same species; no fix needed

# Look for scientific names with multiple codes
model14.multi.sci.name <- count(model14.species, ScientificName) |> 
  arrange(desc(n)) |> 
  filter(n > 1) 

model14.multi.sci <- model14.species |> 
  filter(ScientificName %in% model14.multi.sci.name$ScientificName) |> 
  filter(!is.na(ScientificName)) |> 
  arrange(ScientificName) # this actually needs a fix


### 15. Northern Basin & Range: Prescribed burn ---------------------------

# Join with geospecies
model15.geospecies <- model15.matched |> 
  left_join(geospecies)

# Create species list
model15.species <- model15.geospecies |> 
  select(Species, ScientificName) |> 
  distinct(.keep_all = TRUE)

# Look for codes with multiple scientific names
model15.multi.code.codes <- count(model15.species, Species) |> 
  arrange(desc(n)) |> 
  filter(n > 1)

model15.multi.code <- model15.species |> 
  filter(Species %in% model15.multi.code.codes$Species) |> 
  arrange(Species) # LUPINPF is included here

# Look for scientific names with multiple codes
model15.multi.sci.name <- count(model15.species, ScientificName) |> 
  arrange(desc(n)) |> 
  filter(n > 1) 

model15.multi.sci <- model15.species |> 
  filter(ScientificName %in% model15.multi.sci.name$ScientificName) |> 
  filter(!is.na(ScientificName)) |> 
  arrange(ScientificName)


### 16. Northern Basin & Range: Vegetation disturbance --------------------

# Join with geospecies
model16.geospecies <- model16.matched |> 
  left_join(geospecies)

# Create species list
model16.species <- model16.geospecies |> 
  select(Species, ScientificName) |> 
  distinct(.keep_all = TRUE)

# Look for codes with multiple scientific names
model16.multi.code.codes <- count(model16.species, Species) |> 
  arrange(desc(n)) |> 
  filter(n > 1)

model16.multi.code <- model16.species |> 
  filter(Species %in% model16.multi.code.codes$Species) |> 
  arrange(Species) # LUPINPF included

# Look for scientific names with multiple codes
model16.multi.sci.name <- count(model16.species, ScientificName) |> 
  arrange(desc(n)) |> 
  filter(n > 1) 

model16.multi.sci <- model16.species |> 
  filter(ScientificName %in% model16.multi.sci.name$ScientificName) |> 
  filter(!is.na(ScientificName)) |> 
  arrange(ScientificName) 








save.image("RData/17_standardize-species-codes.RData")

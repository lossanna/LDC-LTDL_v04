# Created: 2026-06-30
# Updated: 2026-07-06

# Purpose: Standardize plant species names/codes, and assign cover values for all matched data.
#   Include codes in geospecies with no LPI data (these species were found during the
#   species inventory but not measured on the transect for LPI). Assign these species a 
#   cover value of 0.5% to be able to calculate diversity later on.

# devtools::install_github("landscape-data-commons/trex", build_vignettes = TRUE)

library(tidyverse)
library(foreign)
library(terradactyl)

# Load data ---------------------------------------------------------------

load("RData/15.1_matched-data.RData")
plants.db.raw <- read.dbf("data/raw/tblNationalPlants/tblNationalPlants.dbf")
geospecies <- read_csv("data/versions-from-R/16_geospecies.csv")
lpi <- read_csv("data/versions-from-R/16_data-lpi.csv")


# Combine all matched data ------------------------------------------------

# Bind rows
all.matched <- bind_rows(mget(ls(pattern = "\\.matched$")),
                         .id = "Model") |> 
  mutate(Model = as.integer(gsub("\\D+", "", Model)))

# List of relevant species
matched.species <- geospecies |> 
  filter(PrimaryKey %in% all.matched$PrimaryKey) |> 
  select(Species, ScientificName, Duration, Woody, Lifeform) |> 
  distinct(.keep_all = TRUE)



# Data wrangling for National Plants DB -----------------------------------

# Narrow down to NameCodes that appear in matched.species
plants.db <- plants.db.raw |> 
  filter(NameCode %in% matched.species$Species) |> 
  rename(CurrentPLANTSCode = CurrentPLA,
         GrowthHabitat = GrowthHabi,
         GrowthHabitat_sub = GrowthHa_1)


## Create crosswalk -------------------------------------------------------

# Cols for grouping
crosswalk1 <- plants.db |> 
  select(NameCode, Scientific, CurrentPLANTSCode, GrowthHabitat, GrowthHabitat_sub, Duration,
         Nonnative, Invasive) |> 
  arrange(NameCode) |> 
  arrange(CurrentPLANTSCode) 

# Inspect rows with multiple same CurrentPLANTSCode
multiple.codes.codes <- count(crosswalk1, CurrentPLANTSCode) |> 
  arrange(desc(n)) |> 
  filter(n > 1)

multiple.codes <- crosswalk1 |> 
  filter(CurrentPLANTSCode %in% multiple.codes.codes$CurrentPLANTSCode)


### Example of multiple CurrentPLANTSCode: Lupin --------------------------

# This section is to demonstrate the situation of multiple CurrentPLANTSCode values
#   that have different NameCodes, and how to handle them.

# Example: LUPIN
lupin <- multiple.codes |> 
  filter(CurrentPLANTSCode == "LUPIN")
lupin

# Lupinus is both a genus and species (frustratingly, the species has no species name). 
# When the scientific name is "Lupinus", we assume this refers to the species. However, 
#   when there are qualifiers in parentheses (i.e., "annual forb") we assume that this 
#   refers to individuals not identified to the species level.
# Therefore, "Lupinus (annual forb)" should not be grouped with "Lupinus", but the
#   CurrentPLANTSCode has all of them listed as "LUPIN". Normally, we would group
#   by CurrentPLANTSCode to collapse species with multiple/former names into a single 
#   group, but we do not want to do that in this instance.
# Hence, NameCodes with parentheses in the scientific name usually refer to individuals
#   not indentified to the species level, and should remain as a separate group.


### Handle scientific names with parentheses ------------------------------

# Inspect scientific names with parentheses
sci.parenth <- crosswalk1 |> 
  filter(str_detect(Scientific, fixed("("))) # all of these are not identified to the species
#       level, so their original NameCode should be retained

#   Change CurrentPLANTSCode to NameCode for scientific names with parentheses
sci.parenth.fix <- sci.parenth |> 
  mutate(CurrentPLANTSCode = NameCode)


# Remove codes of those not identified to species level
crosswalk2 <- crosswalk1 |> 
  filter(!NameCode %in% sci.parenth$NameCode)

#   Add back fixed version where CurrentPLANTSCode is NameCode
crosswalk2 <- crosswalk2 |> 
  bind_rows(sci.parenth.fix)


# Ensure other cols are standardized for each CurrentPLANTSCode
species.info <- crosswalk2 |> 
  select(-NameCode, -Scientific) |> 
  distinct(.keep_all = TRUE)

#   Extract problem codes
problem.codes.codes <- count(species.info, CurrentPLANTSCode) |> 
  filter(n > 1)

problem.codes <- crosswalk2 |> 
  filter(CurrentPLANTSCode %in% problem.codes.codes$CurrentPLANTSCode)

#   Apply fixes to problem codes
# Erigeron tracyi is also called Erigeron colomexicanus; code is ERCO28
# idk what Lepidium cf. montanum is, so we can just retain its original NameCode
# POAG is the same as POPRP2, which is invasive in the lower 48
problem.codes.fixed <- problem.codes |> 
  filter(NameCode %in% c("ERITRA", "LEPCFMON", "POAG")) |> 
  mutate(CurrentPLANTSCode = c("ERCO28", "LEPCFMON", "POPRP2"),
         Nonnative = c("NATIVE", "NATIVE", "EXOTIC"),
         Invasive = c(NA, NA, "INVASIVE"))

# Remove problem codes from crosswalk and add back fix
crosswalk3 <- crosswalk2 |> 
  filter(!NameCode %in% problem.codes.fixed$NameCode) |> 
  bind_rows(problem.codes.fixed)


# Save final version
crosswalk <- crosswalk3



# Apply crosswalk fix -----------------------------------------------------

# Filter LPI data for primary keys in all.matched
all.matched.lpi <- lpi |> 
  filter(PrimaryKey %in% all.matched$PrimaryKey)

# Rename Code col as NameCode
all.matched.lpi <- all.matched.lpi |> 
  rename(NameCode = Code)

# Join crosswalk
all.matched.lpi.joined <- all.matched.lpi |> 
  left_join(crosswalk) 

# Change names so pct_cover() function will work
#   Example LPI data
lpi_tall <- terradactyl::tall_lpi_sample 

#   Inspect column names
colnames(lpi_tall)
colnames(all.matched.lpi.joined)

#   Change differently named columns to match example LPI data
all.matched.lpi.joined <- all.matched.lpi.joined |> 
  rename(code = NameCode,
         PointNbr = PointNumber,
         RecKey = RecordKey,
         LineLengthAmount = LineLength,
         SpacingIntervalAmount = SpacingIntervalAmt,
         PointLoc = PointLocation,
         source = Source,
         DBKey = DatabaseKey,
         layer = Layer)



# Recalculate LPI data ----------------------------------------------------

# Recalculate cover
lpi.recalculate1 <- pct_cover(lpi_tall = all.matched.lpi.joined,
                              tall = TRUE, 
                              hit = "any",
                              indicator_variables = "CurrentPLANTSCode")

# Drop 0s and rename cols
lpi.recalculate2 <- lpi.recalculate1 |> 
  filter(percent > 0) |> 
  rename(CurrentPLANTSCode = indicator,
         Cover_AH = percent) 

# Save final version
lpi.recalculate <- lpi.recalculate2



# Create cover table for non-LPI species ----------------------------------

# Apply crosswalk fix
#   Filter geospecies for primary keys in all.matched
all.matched.geospecies <- geospecies |> 
  filter(PrimaryKey %in% all.matched$PrimaryKey)

#   Rename Species col to NameCode
all.matched.geospecies <- all.matched.geospecies |> 
  rename(NameCode = Species) 

#   Join crosswalk
all.matched.geospecies.joined <- all.matched.geospecies |> 
  left_join(crosswalk) 

#   Collapse into groups via CurrentPLANTSCode
all.matched.geospecies.fixed <- all.matched.geospecies.joined |> 
  select(PrimaryKey, SpeciesCover_AH, CurrentPLANTSCode) |> 
  distinct(.keep_all = TRUE)


# Narrow down to only plants with NA species cover (not included in LPI)
all.matched.nonlpi <- all.matched.geospecies.fixed |> 
  filter(is.na(SpeciesCover_AH))

# Assign 0.5% cover and rename col
all.matched.nonlpi <- all.matched.nonlpi |> 
  mutate(SpeciesCover_AH = 0.5) |> 
  rename(Cover_AH = SpeciesCover_AH)



# Combine cover data ------------------------------------------------------

# Combine and add species info
cover <- lpi.recalculate |> 
  bind_rows(all.matched.nonlpi) 

# Create df of species info with fixed codes
species.info.fixed <- crosswalk |> 
  select(-NameCode, -Scientific) |> 
  distinct(.keep_all = TRUE)

# Add species info to cover
cover <- cover |> 
  left_join(species.info.fixed)


# Add to species cover data to all.matched --------------------------------

all.matched.cover <- all.matched |> 
  left_join(cover,
            relationship = "many-to-many") # relationship is many-to-many because all primary
#   keys in all.matched will have multiple rows in cover (many species), and rows in
#   cover might match multiple rows in all.matched because there are multiple instances of
#   the same primary key because the same point might be used in multiple models.



# Write to CSV ------------------------------------------------------------

write_csv(all.matched.cover,
          file = "data/versions-from-R/17_all-matched-data-with-species-cover_v012.csv",
          na = "")

save.image("RData/17_standardize-species-codes.RData")

# Created: 2026-07-08
# Updated: 2026-07-08

# Purpose: Sum invasive species cover for each LDC point. (Species categorized as invasive
#   via the National Plants database from Cecelia.)


library(tidyverse)

# Load data ---------------------------------------------------------------

all.matched.cover <- read_csv("data/versions-from-R/17_all-matched-data-with-species-cover_v012.csv")


# Data wrangling ----------------------------------------------------------

# Invasive species only, and remove duplicate instances of primary keys
invasive <- all.matched.cover |> 
  filter(Invasive == "INVASIVE") |> 
  select(PrimaryKey, CurrentPLANTSCode, Cover_AH, Invasive) |> 
  distinct(.keep_all = TRUE)
  
# Sum invasive cover by PrimaryKey
invasive.summed <- invasive |> 
  group_by(PrimaryKey) |> 
  summarise(InvasiveCover = sum(Cover_AH))


# Assign missing PrimaryKeys invasive cover value of 0
invasive0 <- data.frame(
  PrimaryKey = unique(all.matched.cover$PrimaryKey),
  InvasiveCover = 0
)

invasive0 <- invasive0 |> 
  filter(!PrimaryKey %in% invasive$PrimaryKey)


# Combine
invasive.all <- invasive.summed |> 
  bind_rows(invasive0)

#   Check for no duplicate PrimaryKeys
length(unique(invasive.all$PrimaryKey)) == nrow(invasive.all)


# Add back other cols
all.matched.join <- all.matched.cover |> 
  select(-CurrentPLANTSCode, -Cover_AH, -GrowthHabitat, -GrowthHabitat_sub,
         -Duration, -Nonnative, -Invasive) |> 
  distinct(.keep_all = TRUE)

invasive.all <- all.matched.join |> 
  left_join(invasive.all)



# Write to CSV ------------------------------------------------------------

write_csv(invasive.all,
          file = "data/versions-from-R/19_invasive-cover_all-models.csv",
          na = "")


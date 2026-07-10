# Created: 2026-07-08
# Updated: 2026-07-10

# Purpose: Calculate Shannon diversity for each LDC point.

library(tidyverse)
library(vegan)

# Load data ---------------------------------------------------------------

all.matched.cover <- read_csv("data/versions-from-R/17_all-matched-data-with-species-cover_v012.csv")


# Calculate Shannon diversity ---------------------------------------------

# Remove duplicate instances of primary keys
cover <- all.matched.cover |> 
  select(PrimaryKey, CurrentPLANTSCode, Cover_AH) |> 
  distinct(.keep_all = TRUE)

# Write function
calculate_shannon <- function(df) {
  
  comm <- df |>
    group_by(PrimaryKey, CurrentPLANTSCode) |>
    summarise(Cover = sum(Cover_AH), .groups = "drop") |>
    pivot_wider(
      names_from = CurrentPLANTSCode,
      values_from = Cover,
      values_fill = 0
    )
  
  shannon <- diversity(
    as.matrix(comm[, -1]),
    index = "shannon"
  )
  
  tibble(
    PrimaryKey = comm$PrimaryKey,
    Shannon = shannon
  )
}
  

# Run function
all.diversity <- calculate_shannon(cover)

# Add back other cols
all.matched.join <- all.matched.cover |> 
  select(-CurrentPLANTSCode, -Cover_AH, -GrowthHabit, -GrowthHabitSub, -GrowthHabitSub2,
         -Duration, -Nonnative, -Invasive) |> 
  distinct(.keep_all = TRUE)

all.diversity <- all.matched.join |> 
  left_join(all.diversity)


# Write to CSV ------------------------------------------------------------

write_csv(all.diversity,
          file = "data/versions-from-R/18_shannon-diversity_all-models_v012.csv",
          na = "")

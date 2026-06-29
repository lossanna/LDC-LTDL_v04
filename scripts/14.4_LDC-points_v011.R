# Created: 2026-06-29
# Updated: 2026-06-29

# Purpose: Recalculate groups to find treatments per ecoregion with 30+ points, now that
#   CETWI, SOLUS, and geoindicators data have been added (and some LDC points therefore removed).

# Same as 12.3.R script from Project v03.

# Result: Nothing with the treatment groups actually changed; only 18 points without
#   gap data were dropped (same as Project v03 result).


library(tidyverse)

# Load data ---------------------------------------------------------------

ldc.cetwi.solus.raw <- read_csv("data/versions-from-R/14.2_LDC-CETWI-SOLUS_v011.csv")
ldc.010.psm <- read_csv("data/versions-from-R/14.3_LDC-with-PSM-cols_v010.csv")
  
# Convert date columns ----------------------------------------------------

# Convert relevant columns to date format
ldc.cetwi.solus <- ldc.cetwi.solus.raw |>
  mutate(
    DateVisited = as.Date(DateVisited, format = "%m/%d/%Y"),
    MR_trt_comp = as.Date(MR_trt_comp, format = "%m/%d/%Y"),
    MR_wildfire = as.Date(MR_wildfire, format = "%m/%d/%Y")
  )


# Add cols from geoindicators.csv -----------------------------------------

# Add via ldc.010.psm (see 14.3.R script)
ldc.010.psm.join <- ldc.010.psm |> 
  select(LDCpointID, BareSoil_FH, TotalFoliarCover, ForbCover_AH, GramCover_AH,
         ShrubCover_AH, Gap100plus)
  
ldc.cetwi.solus <- ldc.cetwi.solus |> 
  left_join(ldc.010.psm.join)


# Check for NAs
apply(ldc.cetwi.solus, 2, anyNA)

#   Check which is missing gap data
ldc.cetwi.solus |>
  filter(is.na(Gap100plus)) # only one treated missing gap data; Chihuahua herbicide

# Drop points missing gap data
ldc.cetwi.solus <- ldc.cetwi.solus |> 
  filter(!is.na(Gap100plus))


# LDC points by Ecoregion 3 -----------------------------------------------

# Count table of at least 30 points per treatment group
level3.trt30 <- ldc.cetwi.solus |>
  group_by(EcoLvl1, EcoLvl2, EcoLvl3, Category) |>
  count(Trt_Type_Sub) |>
  filter(
    n >= 30,
    !is.na(Trt_Type_Sub)
  ) |>
  ungroup()
length(unique(level3.trt30$EcoLvl3)) # 13
unique(level3.trt30$Trt_Type_Sub)


# LDC points with at least 30 points per treatment group
eco3.trt30 <- level3.trt30 |>
  left_join(ldc.cetwi.solus)

#   Control equivalents
eco3.trt30.ctrl.count <- ldc.cetwi.solus |>
  filter(EcoLvl3 %in% eco3.trt30$EcoLvl3 & str_detect(Category, "Control")) |>
  group_by(EcoLvl3, Category) |>
  summarise(
    n = n(),
    .groups = "keep"
  ) |>
  ungroup()

eco3.trt30.ctrl <- ldc.cetwi.solus |>
  filter(EcoLvl3 %in% eco3.trt30$EcoLvl3 & str_detect(Category, "Control")) |>
  left_join(eco3.trt30.ctrl.count)

#   Combine & order cols
eco3.trt30.all <- eco3.trt30 |>
  bind_rows(eco3.trt30.ctrl) |>
  select(
    EcoLvl3, Category, Trt_Type_Sub, n, MR_trt_comp, LDCpointID, DateVisited,
    ProjectKey, PrimaryKey,
    Trt_Type_Major, recent_trt_count, FirePolyID, USGS_Assigned_ID, MR_wildfire,
    Fire_freq, Fire_freq_post_trt,
    EcoLvl1, NA_L1KEY, EcoLvl2, NA_L2KEY, NA_L3KEY, MLRARSYM, MLRA_name,
    horizontal_flux_total_MD, n_CETWI, CETWI, sandtotal_0_cm,
    BareSoil_FH, TotalFoliarCover, ForbCover_AH, GramCover_AH, ShrubCover_AH, Gap100plus
  ) |>
  rename(Treatment_count = n) |>
  arrange(LDCpointID)
nrow(eco3.trt30.all) / nrow(ldc.cetwi.solus.raw) # 99.96% used
nrow(eco3.trt30.all) == nrow(ldc.cetwi.solus) # only points that were dropped are ones without gap data



# Write to CSV ------------------------------------------------------------

write_csv(eco3.trt30.all,
          file = "data/versions-from-R/14.4_LDC-points_v011.csv",
          na = ""
)

write_csv(level3.trt30,
          file = "data/versions-from-R/14.4_eco3-trt30_count-table.csv")

save.image("RData/14.4_LDC-points_v011.RData")

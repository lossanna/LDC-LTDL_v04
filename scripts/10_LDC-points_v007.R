# Created: 2026-06-25
# Updated: 2026-06-25

# Purpose: Create LDC points v007, which is groups of 30+ points per treatment,
#   based on LDC points v006 (no treatment combinations).

# Eliminating treatment combinations does not increase the number of models (it
#   actually reduces it), so LDC points v007 will be discontinued/abandoned.


library(tidyverse)

# Load data ---------------------------------------------------------------

ldc.006.raw <- read_csv("data/versions-from-R/09.2_LDC-points_v006.csv")
ldc.004.raw <- read_csv("data/GIS-exports/004_LDCpts004_export.csv")


# Convert date columns ----------------------------------------------------

# Convert relevant columns to date format
ldc.006 <- ldc.006.raw |>
  rename(MR_wildfire = Wildfire_Date,
         MR_trt_comp = comp_date_est)


# Correct treatment info for TRT -> Fire -> LDC ---------------------------

# Examine post-burn control points
post.burn.ctrl <- ldc.006 |>
  filter(Category == "Control, post-burn")

#   See that some have treatment info
unique(post.burn.ctrl$Trt_Type_Major)
post.burn.ctrl.trt <- post.burn.ctrl |>
  filter(!is.na(Trt_Type_Major))

#   Examine dates
post.burn.ctrl.trt |>
  select(DateVisited, MR_trt_comp, MR_wildfire)
post.burn.ctrl.trt$MR_wildfire >= post.burn.ctrl.trt$MR_trt_comp # all have wildfire occur after
#                                                       or on the same date as the most recent treatment
#   In the case of TRT -> Fire -> LDC, we assume that the effects of the most recent wildfire
#     wiped out any treatment effects, which is why they are categorized as controls.


# Remove treatment info for post-burn control situations of TRT -> Fire -> LDC
ldc.007 <- ldc.006 |>
  mutate(
    Trt_Type_Major = case_when(
      Category == "Control, post-burn" ~ NA,
      TRUE ~ Trt_Type_Major
    ),
    Trt_Type_Sub = case_when(
      Category == "Control, post-burn" ~ NA,
      TRUE ~ Trt_Type_Sub
    ),
    MR_trt_comp = case_when(
      Category == "Control, post-burn" ~ NA,
      TRUE ~ MR_trt_comp
    )
  )


# LDC points by Ecoregion 3 -----------------------------------------------

# Join with LDC v004 to get Ecoregion cols
ldc.004.join <- ldc.004.raw |> 
  select(PrimaryKey, NA_L1NAME, NA_L2NAME, NA_L3NAME, MLRA_NAME)

ldc.007 <- ldc.007 |> 
  left_join(ldc.004.join)

# Rename ecoregion & MLRA columns
ldc.007 <- ldc.007 |>
  rename(
    EcoLvl1 = NA_L1NAME,
    EcoLvl2 = NA_L2NAME,
    EcoLvl3 = NA_L3NAME,
    MLRA_name = MLRA_NAME
  )

# Count table of at least 30 points per treatment group
level3.trt30 <- ldc.007 |>
  group_by(EcoLvl1, EcoLvl2, EcoLvl3, Category) |>
  count(Trt_Type_Sub) |>
  filter(
    n >= 30,
    !is.na(Trt_Type_Sub)
  ) |>
  ungroup()
length(unique(level3.trt30$EcoLvl3)) # 13
unique(level3.trt30$Trt_Type_Sub)




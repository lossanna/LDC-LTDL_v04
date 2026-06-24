# Created: 2026-06-24
# Updated: 2026-06-24

# Purpose: Create LDC points v005, which is groups of 30+ points per treatment.


library(tidyverse)

# Load data ---------------------------------------------------------------

ldc.004.raw <- read_csv("data/GIS-exports/004_LDCpts004_export.csv")


# Convert date columns ----------------------------------------------------

# Convert relevant columns to date format
ldc.004 <- ldc.004.raw |>
  mutate(
    DateVisited = as.Date(DateVisited, format = "%m/%d/%Y"),
    MR_trt_comp = as.Date(MR_trt_comp, format = "%m/%d/%Y"),
    MR_wildfire = as.Date(MR_wildfire, format = "%m/%d/%Y")
  )


# Correct treatment info for TRT -> Fire -> LDC ---------------------------

# Examine post-burn control points
post.burn.ctrl <- ldc.004 |>
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
ldc.005 <- ldc.004 |>
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
    ),
    recent_trt_count = case_when(
      Category == "Control, post-burn" ~ NA,
      TRUE ~ recent_trt_count
    )
  )


# LDC points by Ecoregion 3 -----------------------------------------------

# Rename ecoregion & MLRA columns
ldc.005 <- ldc.005 |>
  rename(
    EcoLvl1 = NA_L1NAME,
    EcoLvl2 = NA_L2NAME,
    EcoLvl3 = NA_L3NAME,
    MLRA_name = MLRA_NAME
  )

# Count table of at least 30 points per treatment group
level3.trt30 <- ldc.005 |>
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
  left_join(ldc.005)

#   Control equivalents
eco3.trt30.ctrl.count <- ldc.005 |>
  filter(EcoLvl3 %in% eco3.trt30$EcoLvl3 & str_detect(Category, "Control")) |>
  group_by(EcoLvl3, Category) |>
  summarise(
    n = n(),
    .groups = "keep"
  ) |>
  ungroup()

eco3.trt30.ctrl <- ldc.005 |>
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
    horizontal_flux_total_MD
  ) |>
  rename(Treatment_count = n) |>
  arrange(LDCpointID)
nrow(eco3.trt30.all) / nrow(ldc.004) # 89% of points used


# Write to CSV ------------------------------------------------------------

write_csv(eco3.trt30.all,
  file = "data/versions-from-R/08_LDC-points_v005.csv",
  na = ""
)


save.image("RData/08_LDC-points_v005.RData")

# Created: 2026-06-24
# Updated: 2026-06-24

# Purpose: Assign land history (treatment + wildfire info) to each LDC point,
#   based on treatment info v006 (only single most recent treatment included).


library(tidyverse)

# Load data ---------------------------------------------------------------

ldc.wildfire.sjoin.raw <- read_csv("data/GIS-exports/004_LDC003-Wildfires004-SpatialJoin_export.csv")
wildfires.assigned <- read_csv("data/versions-from-R/07.1_wildfires-with-assigned-dates_v004.csv")
treatment.info.006 <- read_csv("data/versions-from-R/09.1_Treatment-info_v006.csv")


# LDC cols alone ----------------------------------------------------------

ldc <- treatment.info.006 |>
  select(LDCpointID, ProjectKey, PrimaryKey, DateVisited)


# LDC points and land treatments ------------------------------------------

# LDC points that overlap with land treatments
ldc.trt <- treatment.info.006 |>
  filter(!is.na(Trt_Type_Sub))

# LDC points without land treatments
ldc.no.trt <- treatment.info.006 |>
  filter(is.na(Trt_Type_Sub))


# Join wildfire dates with LDC points -------------------------------------

# Join with wildfire dates
ldc.wildfire <- ldc.wildfire.sjoin.raw |>
  select(
    LDCpointID, ProjectKey, PrimaryKey, DateVisited,
    FirePolyID, USGS_Assigned_ID
  ) |>
  left_join(wildfires.assigned)

# Filter out fires that occurred after LDC date
ldc.wildfire <- ldc.wildfire |>
  mutate(DateVisited = as.Date(DateVisited, format = "%m/%d/%Y")) |>
  filter(DateVisited >= Wildfire_Date)


# Calculate number of times burned
ldc.wildfire.freq <- ldc.wildfire |>
  group_by(LDCpointID, PrimaryKey) |>
  summarise(
    Fire_freq = n(),
    .groups = "keep"
  )


# LDC points not burned
ldc.not.burned <- ldc |>
  filter(!LDCpointID %in% ldc.wildfire$LDCpointID)

#   Check for duplicate LDC points
count(ldc.not.burned, LDCpointID) |>
  arrange(desc(n))


# Join wildfires with treatments ------------------------------------------

# Join LDC + treatment with LDC + wildfire
ldc.trt.wildfire <- ldc.wildfire |>
  left_join(treatment.info.006)

# Filter out fires that occurred after most recent treatment date
ldc.trt.wildfire.post.trt <- ldc.trt.wildfire |>
  filter(comp_date_est > Wildfire_Date)

# Calculate number of times burned pre-treatment
ldc.wildfire.freq.post.trt <- ldc.trt.wildfire.post.trt |>
  group_by(LDCpointID, PrimaryKey) |>
  summarise(
    Fire_freq_post_trt = n(),
    .groups = "keep"
  )


# Categories --------------------------------------------------------------

## 1. Control, post-burn (Fire -> LDC) ------------------------------------

# Does not overlap with land treatment, previously burned
ldc.control.burned1 <- ldc.trt.wildfire |>
  filter(is.na(Trt_Type_Sub)) |>
  left_join(ldc.wildfire.freq) |>
  mutate(Fire_freq_post_trt = NA) |>
  select(-Fire_Year, -Date_type, -Date_est) |>
  group_by(LDCpointID, PrimaryKey) |>
  slice_max(Wildfire_Date, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(Category = "Control, post-burn")


## 2. Control, never burned -----------------------------------------------

# Does not overlap with land treatment, never burned
ldc.control.not.burned2 <- ldc.not.burned |>
  filter(LDCpointID %in% ldc.no.trt$LDCpointID) |>
  mutate(
    FirePolyID = NA,
    USGS_Assigned_ID = NA,
    Wildfire_Date = NA,
    Fire_freq = NA,
    Fire_freq_post_trt = NA,
    Trt_Type_Major = NA,
    Trt_Type_Sub = NA,
    comp_date_est = NA
  ) |>
  mutate(Category = "Control, never burned")


## 3. Treated, never burned (TRT -> LDC) ----------------------------------

# Overlaps with land treatment, never burned
ldc.treated.not.burned3 <- ldc.not.burned |>
  filter(LDCpointID %in% ldc.trt$LDCpointID) |>
  left_join(treatment.info.006) |>
  mutate(
    FirePolyID = NA,
    USGS_Assigned_ID = NA,
    Wildfire_Date = NA,
    Fire_freq = NA,
    Fire_freq_post_trt = NA
  ) |>
  mutate(Category = "Treated, never burned")


## 5/6. Control, post-burn (TRT -> Fire -> LDC) ---------------------------

# Overlaps with land treatment, but fire happened in between TRT and LDC
#   (can also be [6] Fire -> TRT -> Fire -> LDC)
ldc.control.burned56 <- ldc.trt.wildfire |>
  filter(!is.na(Trt_Type_Sub)) |>
  filter(comp_date_est <= Wildfire_Date) |>
  left_join(ldc.wildfire.freq) |>
  left_join(ldc.wildfire.freq.post.trt) |>
  select(-Fire_Year, -Date_type, -Date_est) |>
  group_by(LDCpointID, PrimaryKey) |>
  slice_max(Wildfire_Date, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(Category = "Control, post-burn")


## 4. Treated, post-burn (MR Fire -> TRT -> LDC) ---------------------------

# Overlaps with land treatment, previously burned
ldc.treated.burned4 <- ldc.trt.wildfire |>
  filter(!is.na(Trt_Type_Sub)) |>
  filter(comp_date_est > Wildfire_Date) |>
  filter(!LDCpointID %in% ldc.control.burned56$LDCpointID) |>
  left_join(ldc.wildfire.freq) |>
  left_join(ldc.wildfire.freq.post.trt) |>
  select(-Fire_Year, -Date_type, -Date_est) |>
  group_by(LDCpointID, PrimaryKey) |>
  slice_max(Wildfire_Date, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(Category = "Treated, post-burn")


# Combine for v006 --------------------------------------------------------

# Combine categories
ldc.006 <- bind_rows(
  ldc.control.burned1, ldc.control.not.burned2,
  ldc.treated.not.burned3, ldc.treated.burned4,
  ldc.control.burned56
) |>
  arrange(LDCpointID)




# Write to CSV ------------------------------------------------------------

# LDC points 
write_csv(ldc.006,
          file = "data/versions-from-R/09.2_LDC-points_v006.csv",
          na = ""
)

save.image("RData/09.2_assign-land-history_v006.RData")

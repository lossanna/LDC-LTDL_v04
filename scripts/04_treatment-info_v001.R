# Created: 2026-06-19
# Updated: 2026-06-19

# Purpose: Complete initial data cleaning to create Treatment_Info v001.

# Data cleaning:
#   Filter for polygons, implemented plans, confirmed features only.
#   Clean dates (Init_Date, Comp_Date).


library(tidyverse)

# Load data ---------------------------------------------------------------

treatment.info.raw <- read_csv("data/versions-from-R/01_Treatment-info_v000.csv")


# Initial filtering -------------------------------------------------------

# Polygons only
treatment.info <- treatment.info.raw |>
  filter(Trt_Feature_Type == "Polygon")

# Implemented plans only
treatment.info <- treatment.info |>
  filter(Plan_Imp == "Implemented")

# Confirmed features only
count(treatment.info, Feature_Status)
treatment.info$Feature_Status[treatment.info$Feature_Status == "confirmed"] <- "Confirmed"
treatment.info <- treatment.info |>
  filter(Feature_Status == "Confirmed")


# Check for treatment hierarchy -------------------------------------------

# Hierarchy test
hierarchy.violations <- treatment.info |>
  distinct(Trt_Type_Major, Trt_Type_Sub, Treatment_Type) |>
  summarise(
    n_sub = n_distinct(Trt_Type_Sub),
    n_major = n_distinct(Trt_Type_Major),
    .by = Treatment_Type
  ) |>
  filter(n_sub > 1 | n_major > 1)
hierarchy.violations

# Investigate problem Treatment_Types
treatment.info |>
  filter(Treatment_Type %in% hierarchy.violations$Treatment_Type) |>
  count(Treatment_Type, Trt_Type_Major, Trt_Type_Sub) # clearly some were just miscategorized

# Remove problem Treatment_Types to be able to make hierarchy table
treatment.info.hierarchy <- treatment.info |>
  filter(
    !(Treatment_Type == "Ground Seeding: Drill" & Trt_Type_Major == "Vegetation/Soil Manipulation") &
      !(Treatment_Type == "Timber: Timber Harvest" & Trt_Type_Major == "Other") &
      !(Treatment_Type == "Aerial Seeding" & Trt_Type_Major == "Vegetation/Soil Manipulation") &
      !(Treatment_Type == "General Treatment: Vegetation Removal" & Trt_Type_Major == "Vegetation/Soil Manipulation") &
      !(Treatment_Type == "Vegetation Disturbance: Thinning" & Trt_Type_Sub == "Timber")
  )

# Build hierarchy table
major.to.treatment <- treatment.info.hierarchy |>
  distinct(Trt_Type_Major, Trt_Type_Sub, Treatment_Type) |>
  arrange(Treatment_Type) |>
  arrange(Trt_Type_Sub) |>
  arrange(Trt_Type_Major)

#   Check for matching lengths to verify hierarchy (should have only one row per Treatment_Type)
length(unique(treatment.info$Treatment_Type)) == nrow(major.to.treatment)


# Write treatment hierarchy table to CSV ----------------------------------

write_csv(major.to.treatment,
  file = "data/data-wrangling-intermediate/04_treatment-hierarchy-table_v001.csv"
)


# Init_Date ---------------------------------------------------------------

# Split date into three columns
init.date <- treatment.info |>
  select(Trt_ID, Init_Date) |>
  separate(Init_Date, into = c("init_m", "init_d", "init_y"), sep = "/", remove = FALSE)

# Check for NA
apply(init.date, 2, anyNA) # no NAs

# Check for occurrences of 00
count(init.date, init_m)
count(init.date, init_d) |>
  print(n = 33)

# Check 0s
init.date0 <- init.date |>
  filter(init_d == "0" | init_m == "0") # 0 is not mixed with 00

# Create estimate month and day for 00 (unknown)
init.date <- init.date |>
  mutate(init_m_est = case_when(
    init_m == "00" ~ "6",
    init_m == "0" ~ "6",
    TRUE ~ init_m
  )) |>
  mutate(init_d_est = case_when(
    init_d == "00" ~ "15",
    init_d == "0" ~ "15",
    TRUE ~ init_d
  )) |>
  mutate(init_est = case_when(
    init_m == "00" & init_d == "00" ~ "estimated month and day",
    init_m == "0" & init_d == "0" ~ "estimated month and day",
    init_m == "00" ~ "estimated month",
    init_m == "0" ~ "estimated month",
    init_d == "00" ~ "estimated day",
    init_d == "0" ~ "estimated day",
    TRUE ~ "actual"
  )) |>
  mutate(init_m_est = str_pad(init_m_est, width = 2, pad = "0")) |>
  mutate(init_d_est = str_pad(init_d_est, width = 2, pad = "0")) |>
  mutate(init_date_est = paste0(init_y, "-", init_m_est, "-", init_d_est)) |>
  mutate(init_date_est = as.Date(init_date_est))

# Look for NA dates
init.date.na <- init.date |>
  filter(is.na(init_date_est)) # NAs created because of invalid dates for the month

# Adjust date to fit month for NAs
init.date.na.fix <- init.date.na |>
  mutate(
    init_d_est = c("30", "30", "28", "28", "30", "20"),
    init_est = rep("estimated day", nrow(init.date.na))
  ) |>
  mutate(init_date_est = paste0(init_y, "-", init_m_est, "-", init_d_est)) |>
  mutate(init_date_est = as.Date(init_date_est))

# Replace fixed NA estimated dates
init.date <- init.date |>
  filter(!is.na(init_date_est)) |>
  bind_rows(init.date.na.fix) |>
  arrange(Trt_ID)

# Join with treatment.info
treatment.info <- treatment.info |>
  left_join(init.date)


# Comp_Date ---------------------------------------------------------------

# Split date into three columns
comp.date <- treatment.info |>
  select(Trt_ID, Comp_Date) |>
  separate(Comp_Date, into = c("comp_m", "comp_d", "comp_y"), sep = "/", remove = FALSE)

# Check for NA
apply(comp.date, 2, anyNA) # no NAs

# Check for occurrences of 00
count(comp.date, comp_m)
count(comp.date, comp_d) |>
  print(n = 33)

# Check 0s
comp.date0 <- comp.date |>
  filter(comp_d == "0" | comp_m == "0") # 0 is not mixed with 00

# Create estimate month and day for 00 (unknown)
comp.date <- comp.date |>
  mutate(comp_m_est = case_when(
    comp_m == "00" ~ "6",
    comp_m == "0" ~ "6",
    TRUE ~ comp_m
  )) |>
  mutate(comp_d_est = case_when(
    comp_d == "00" ~ "15",
    comp_d == "0" ~ "15",
    TRUE ~ comp_d
  )) |>
  mutate(comp_est = case_when(
    comp_m == "00" & comp_d == "00" ~ "estimated month and day",
    comp_m == "0" & comp_d == "0" ~ "estimated month and day",
    comp_m == "00" ~ "estimated month",
    comp_m == "0" ~ "estimated month",
    comp_d == "00" ~ "estimated day",
    comp_d == "0" ~ "estimated day",
    TRUE ~ "actual"
  )) |>
  mutate(comp_m_est = str_pad(comp_m_est, width = 2, pad = "0")) |>
  mutate(comp_d_est = str_pad(comp_d_est, width = 2, pad = "0")) |>
  mutate(comp_date_est = paste0(comp_y, "-", comp_m_est, "-", comp_d_est)) |>
  mutate(comp_date_est = as.Date(comp_date_est))

# Look for NA dates
comp.date.na <- comp.date |>
  filter(is.na(comp_date_est))
count(comp.date.na, comp_m_est)
count(comp.date.na, comp_d_est) # NAs created because of invalid dates for the month,
#                                   except in one instance, where there is no year

# Adjust date to fit month for NAs
comp.date.na.fix <- comp.date.na |>
  filter(Trt_ID != 73162) |>
  mutate(comp_d_est = case_when(
    comp_m_est == "02" ~ "28",
    TRUE ~ "30"
  )) |>
  mutate(comp_est = rep("estimated day", nrow(comp.date.na) - 1)) |>
  mutate(comp_date_est = paste0(comp_y, "-", comp_m_est, "-", comp_d_est)) |>
  mutate(comp_date_est = as.Date(comp_date_est))

# Investigate Trt_ID 73162 (missing year completed)
trtid.73162 <- treatment.info |>
  filter(Trt_ID == 73162)

#   Assign year to Trt_ID 73162
trtid.73162.fix <- comp.date.na |>
  filter(Trt_ID == 73162) |>
  mutate(
    comp_y = "2020",
    comp_est = "actual",
    comp_date_est = as.Date("2020-02-22")
  ) # confirmed and GIS version of table says 2020

# Replace fixed NA estimated dates
comp.date <- comp.date |>
  filter(!is.na(comp_date_est)) |>
  bind_rows(comp.date.na.fix, trtid.73162.fix) |>
  arrange(Trt_ID)

# Join with treatment.info
treatment.info <- treatment.info |>
  left_join(comp.date)


# Look for negative values of time elapsed --------------------------------

# Add column of days elapsed between treatment initiation and completion
treatment.info <- treatment.info |>
  mutate(init_comp_elapsed = comp_date_est - init_date_est)

# Examine elapse of negative time
neg.elapse <- treatment.info |>
  filter(init_comp_elapsed < 0)


## No dates est -----------------------------------------------------------

# Negative time elapsed with no date estimates (original must be wrong)
neg.elapse.actual <- neg.elapse |>
  filter(init_est == "actual" & comp_est == "actual")
#   Will just go with comp date given since it is the actual


## Init date est ----------------------------------------------------------

# Negative time elapsed with only initiation date estimated
neg.elapse.init.est <- neg.elapse |>
  filter(init_est != "actual" & comp_est == "actual")
#   Will just go with comp date given since it is the actual


# Comp date est -----------------------------------------------------------

# Negative time elapsed with only completion date estimated
neg.elapse.comp.est <- neg.elapse |>
  filter(init_est == "actual" & comp_est != "actual") |>
  select(
    Trt_ID, Init_Date, Comp_Date, Trt_Type_Sub, Treatment_Type, init_est, init_date_est, init_y,
    comp_est, comp_date_est, comp_y
  )

# Check for matching years
neg.elapse.comp.est$init_y == neg.elapse.comp.est$comp_y # all have same init & comp year

# Replace est. comp date with actual init date
neg.elapse.comp.est.fix <- neg.elapse |>
  filter(init_est == "actual" & comp_est != "actual") |>
  mutate(comp_date_est = init_date_est) |>
  mutate(init_comp_elapsed = comp_date_est - init_date_est)


## Both dates est ---------------------------------------------------------

# Negative time elapsed with both dates estimated
neg.elapse.both.est <- neg.elapse |>
  filter(init_est != "actual" & comp_est != "actual") |>
  select(
    Trt_ID, Init_Date, Comp_Date, Trt_Type_Sub, Treatment_Type, init_est,
    init_date_est, comp_est, comp_date_est, init_comp_elapsed
  )

#   When comp day only is estimated (not month), retain estimated date
#   When both comp and init day only (not month) are estimated, retain date

# When comp month & day are estimated and only init day is estimated, replace est. comp date with init date
neg.elapse.both.est.fix <- neg.elapse |>
  filter(init_est == "estimated day" & comp_est == "estimated month and day") |>
  mutate(comp_date_est = init_date_est) |>
  mutate(init_comp_elapsed = comp_date_est - init_date_est)


# Check for unique Trt_ID -------------------------------------------------

treatment.info |>
  summarise(all_unique = n() == n_distinct(Trt_ID))


## Add in fixes -----------------------------------------------------------

treatment.info <- treatment.info |>
  filter(!Trt_ID %in% c(neg.elapse.both.est.fix$Trt_ID, neg.elapse.comp.est.fix$Trt_ID)) |>
  bind_rows(neg.elapse.both.est.fix, neg.elapse.comp.est.fix) |>
  arrange(Trt_ID)


# Separate out columns for GIS join ---------------------------------------

treatment.info.gisjoin <- treatment.info |>
  select(Prj_ID, Trt_ID, init_date_est, comp_date_est, Trt_Type_Major, Trt_Type_Sub, Treatment_Type)


# Write out as v001 -------------------------------------------------------

write_csv(treatment.info,
  file = "data/versions-from-R/04_Treatment-info_v001.csv"
)

write_csv(treatment.info.gisjoin,
  file = "data/versions-from-R/04_Treatment-info_v001-gisjoin.csv"
)


save.image("RData/04_treatment-info_v001.RData")

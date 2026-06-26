# Created: 2026-06-25
# Updated: 2026-06-26

# Purpose: Examine treatment polygons that overlap with LDC points and find most recent
#   treatment and recent treatment combos (within one year of most recent), and make
#   appropriate categories. Only consider treatments implemented before LDC monitoring.

# Filter out rows where comp_est_date > DateVisited.
# Filter out "Cultural Protection" and "Other" treatments (based on Trt_Type_Major).

# Apply the following treatment combinations (most recent treatment first):
#   1. Aerial Seeding, Drill Seeding
#   2. Drill Seeding, Soil Disturbance
#   3. Aerial Seeding, Soil Disturbance


# Differences with 07.1.R:
#   1. Only three treatment combos are included, since Project v03 showed that 
#       these are the only combos with enough points to create a model for PSM.
#   2. Different in that no arbitrary application of treatment preference is applied
#       for points with multiple of the same comp_est_date of the most recent treatments.
#   3. Cultural Protection/Other and rows where comp_est_date > DateVisited are
#       filtered out completely and not added back in as control points with treatment
#       info removed.


library(tidyverse)

# Load data ---------------------------------------------------------------

ldc.trt.sjoin.raw <- read_csv("data/GIS-exports/004_LDC003-TrtPoly002-SpatialJoin_export.csv")


# Preliminary data wrangling ----------------------------------------------

# Convert date cols
ldc.trt.sjoin <- ldc.trt.sjoin.raw |>
  mutate(
    DateVisited = as.Date(DateVisited, format = "%m/%d/%Y"),
    init_date_est = as.Date(init_date_est, format = "%m/%d/%Y"),
    comp_date_est = as.Date(comp_date_est, format = "%m/%d/%Y")
  )

# Separate out points that overlap with treatment polygons
ldc.trt.overlap <- ldc.trt.sjoin |>
  filter(!is.na(TrtPolyID))


# Find instances where comp_est_date > DateVisited ------------------------

# Find instances where comp_date_est <= Date_Visited
trt.pre <- ldc.trt.overlap |>
  filter(comp_date_est <= DateVisited)


# Filter out "Other" and "Cultural Protection" ----------------------------

# Remove because these treatments are unlikely to affect veg, or categories are too broad
trt.pre.filtered <- trt.pre |>
  filter(!Trt_Type_Major %in% c("Cultural Protection", "Other"))



# Extract rows of recent treatments and treatment combos ------------------

# Extract the most recent treatment polygon(s) for each LDC point, and polygons within 1 year
#   of most recent
recent.combo <- trt.pre.filtered |>
  group_by(PrimaryKey, LDCpointID) |>
  mutate(most_recent_comp = max(comp_date_est)) |>
  filter(comp_date_est >= most_recent_comp - years(1)) |>
  mutate(recent_trt_count = n()) |>
  ungroup() |>
  arrange(Trt_Type_Sub) |>
  arrange(desc(comp_date_est)) |>
  arrange(LDCpointID)

# Separate out polygons with overlaps (multiple recent treatments)
recent.overlaps <- recent.combo |>
  filter(recent_trt_count > 1)



# Create treatment combinations -------------------------------------------

# Build hierarchy table (all most recent)
recent.trt.table.all <- recent.combo |>
  distinct(Trt_Type_Major, Trt_Type_Sub, Treatment_Type) |>
  arrange(Treatment_Type) |>
  arrange(Trt_Type_Sub) |>
  arrange(Trt_Type_Major)

#   Check for matching lengths to verify hierarchy (should have only one row per Treatment_Type)
length(unique(recent.combo$Treatment_Type)) == nrow(recent.trt.table.all)

# Build hierarchy table (overlapping treatments only)
recent.trt.table.overlap <- recent.overlaps |>
  distinct(Trt_Type_Major, Trt_Type_Sub, Treatment_Type) |>
  arrange(Treatment_Type) |>
  arrange(Trt_Type_Sub) |>
  arrange(Trt_Type_Major)


# Begin assembling treatment combos based on Trt_Type_Sub
#   (most recent treatment listed first; multiple most recent treatments are alphabetized)
recent.combo.sub <- recent.overlaps |>
  select(LDCpointID, PrimaryKey, Trt_Type_Sub, comp_date_est, most_recent_comp, recent_trt_count) |>
  distinct(.keep_all = TRUE) |>
  arrange(Trt_Type_Sub) |>
  arrange(desc(comp_date_est)) |>
  group_by(LDCpointID, PrimaryKey) |>
  summarise(
    combo_sub = paste(unique(Trt_Type_Sub), collapse = ", "),
    .groups = "drop"
  ) |>
  ungroup() |>
  mutate(sub_count = str_count(combo_sub, ",") + 1) |>
  filter(sub_count > 1) |>
  arrange(LDCpointID)

# Examine combo options
combo.sub.types <- recent.combo.sub |>
  select(-PrimaryKey, -LDCpointID) |>
  distinct(.keep_all = TRUE) |>
  arrange(combo_sub) |>
  arrange(sub_count)

count(recent.combo.sub, combo_sub) |>
  arrange(desc(n)) |>
  print(n = 20)


# Apply combo selection ---------------------------------------------------

# Based on Project v03 results (and 08.R results from this version), there are 
#   three treatment combinations with enough points per ecoregion to graduate to
#   the propensity score matching step (n => 30), so I have selected those only.

# Will consider the following combos:
#   (most recent treatment listed first)
#   1. Aerial Seeding, Drill Seeding
#   2. Drill Seeding, Soil Disturbance
#   3. Aerial Seeding, Soil Disturbance

combo.selection <- c(
  "Aerial Seeding, Drill Seeding", "Drill Seeding, Soil Disturbance",
  "Aerial Seeding, Soil Disturbance"
)


## Separate out points with selected combos -------------------------------

# Separate out points with selected combos (of points with any combo)
combos.select <- recent.combo.sub |>
  filter(combo_sub %in% combo.selection)

# Assign combo_major
combos.select <- combos.select |>
  mutate(combo_major = case_when(
    combo_sub %in% c(
      "Aerial Seeding, Drill Seeding"
    ) ~ "Seeding",
    combo_sub %in% c(
      "Aerial Seeding, Soil Disturbance",
      "Drill Seeding, Soil Disturbance"
    ) ~ " Seeding, Vegetation/Soil Manipulation",
  ))

# Join to get other cols
combos.select <- combos.select |>
  left_join(recent.overlaps)

# Save version that retains all rows
combos.select.allrows <- combos.select

# Save TrtPolyIDs & Trt_IDs
combos.select.trtpolyid <- combos.select |>
  select(LDCpointID, PrimaryKey, TrtPolyID, Trt_ID)

# Retain only combo cols and remove instances of multiple rows
#   (now ready for eventual bind_rows())
combos.select <- combos.select |>
  select(
    -TrtPolyID, -Prj_ID, -Trt_ID, -init_date_est, -comp_date_est,
    -Trt_Type_Major, -Trt_Type_Sub, -Treatment_Type, -recent_trt_count
  ) |>
  distinct(.keep_all = TRUE) |>
  rename(
    Trt_Type_Sub = combo_sub,
    Trt_Type_Major = combo_major,
    recent_trt_count = sub_count
  ) |>
  select(
    LDCpointID, ProjectKey, PrimaryKey, DateVisited,
    Trt_Type_Major, Trt_Type_Sub, most_recent_comp, recent_trt_count
  )


## Combine all other points without selected combos -----------------------

# Get list of all other points without selected combos
non.combo <- trt.pre.filtered |>
  filter(!PrimaryKey %in% combos.select$PrimaryKey)

# Extract rows with most recent treatment polygon(s)
non.combo.mr <- non.combo |>
  group_by(PrimaryKey) |>
  filter(comp_date_est == max(comp_date_est)) |>
  mutate(most_recent_trt_count = n()) |>
  ungroup()


# Only one most recent treatment ------------------------------------------

# Separate out points with only one most recent treatment
#   (treatment info is fine as is)
non.combo.mr.single <- non.combo.mr |>
  filter(most_recent_trt_count == 1)

# Save TrtPolyIDs & Trt_IDs
non.combo.mr.single.trtpolyid <- non.combo.mr.single |>
  select(LDCpointID, PrimaryKey, TrtPolyID, Trt_ID)

# Standardize cols for eventual bind_rows()
#   (now ready for eventual bind_rows())
non.combo.mr.single <- non.combo.mr.single |>
  mutate(
    most_recent_comp = comp_date_est,
    recent_trt_count = 1
  ) |>
  select(
    LDCpointID, ProjectKey, PrimaryKey, DateVisited,
    Trt_Type_Major, Trt_Type_Sub, most_recent_comp, recent_trt_count
  )


# Combinations from single most recent date -------------------------------

# Separate out points with multiple most recent treatments
non.combo.mr.multiple <- non.combo.mr |>
  filter(most_recent_trt_count > 1)


## Same Trt_Type_Sub ------------------------------------------------------

# Separate points that have overlapping polygons with multiple/same comp_date_est
#   and the same Trt_Type_Sub
non.combo.same <- non.combo.mr.multiple |>
  arrange(Trt_Type_Sub) |>
  group_by(PrimaryKey, LDCpointID) |>
  summarise(
    treatments_sub = paste(unique(Trt_Type_Sub), collapse = ", "),
    .groups = "drop"
  ) |>
  ungroup() |>
  mutate(sub_count = str_count(treatments_sub, ",") + 1) |>
  filter(sub_count == 1)


# Join to get other cols
non.combo.same <- non.combo.same |>
  left_join(non.combo.mr.multiple)

# Save TrtPolyIDs & Trt_IDs
non.combo.same.trtpolyid <- non.combo.same |>
  select(LDCpointID, PrimaryKey, TrtPolyID, Trt_ID)

# Remove some cols and instances of multiple rows
#   (now ready for eventual bind_rows())
non.combo.same <- non.combo.same |>
  mutate(
    most_recent_comp = comp_date_est,
    recent_trt_count = 1
  ) |>
  select(
    LDCpointID, ProjectKey, PrimaryKey, DateVisited,
    Trt_Type_Major, Trt_Type_Sub, most_recent_comp, recent_trt_count
  ) |>
  distinct(.keep_all = TRUE)


# I am  going to drop all of the points with different Trt_Type_Sub because it's 
#   impossible to tell which treatment might be having the effect.


# Combine all with treatments ---------------------------------------------

# Combine all
with.combos <- combos.select |>
  bind_rows(non.combo.mr.single) |>
  bind_rows(non.combo.same) |>
  arrange(LDCpointID)

# Check for single instance of primary keys
length(unique(with.combos$PrimaryKey)) == nrow(with.combos)


# TrtPolyIDs --------------------------------------------------------------

# Combine
trtpolyid <- combos.select.trtpolyid |> 
  bind_rows(non.combo.mr.single.trtpolyid, non.combo.same.trtpolyid) |> 
  arrange(TrtPolyID)


# Control points ----------------------------------------------------------

# Make df of LDC points that don't overlap with treatment polygons with
#   standardized cols for bind_rows()
ldc.control <- ldc.trt.sjoin |>
  filter(is.na(TrtPolyID)) |>
  mutate(
    most_recent_comp = NA,
    recent_trt_count = NA
  ) |>
  select(
    LDCpointID, ProjectKey, PrimaryKey, DateVisited,
    Trt_Type_Major, Trt_Type_Sub, most_recent_comp, recent_trt_count
  )

#   Check there are no duplicate primary keys
count(ldc.control, PrimaryKey) |>
  arrange(desc(n))



# bind_rows() to make treatment info v008 ---------------------------------

# Combine
treatment.info.008 <- with.combos |>
  bind_rows(ldc.control) |>
  arrange(LDCpointID)


# Write v008 to CSV -------------------------------------------------------

# Treatment info v008
write_csv(treatment.info.008,
          file = "data/versions-from-R/11.1_Treatment-info_v008.csv",
          na = ""
)

# TrtPolyID list
write_csv(trtpolyid,
          file = "data/versions-from-R/11.1_TrtPolyID-for-treatment-info-v008.csv"
)


save.image("RData/11.1_treatment-info_v008.RData")

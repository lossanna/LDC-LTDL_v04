# Created: 2026-06-24
# Updated: 2026-06-24

# Purpose: Examine treatment polygons that overlap with LDC points and find most recent
#   treatment and make appropriate categories. Only consider treatments implemented 
#   before LDC monitoring.

# Unlike 07.2.R script, this version does NOT include treatment combinations.

# Filter out rows where comp_est_date > DateVisited.
# Filter out "Cultural Protection" and "Other" treatments (based on Trt_Type_Major).

library(tidyverse)
library(readxl)

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


# Treatment based on single most recent date ------------------------------

# Separate out the rows with most recent treatment polygon(s)
most.recent <- trt.pre.filtered |>
  group_by(PrimaryKey) |>
  filter(comp_date_est == max(comp_date_est)) |>
  ungroup() 

nrow(most.recent) / nrow(trt.pre.filtered) # 62%


# Separate points that have overlapping polygons with multiple of the same comp_date_est
most.recent.multiple.id <- most.recent |> 
  count(LDCpointID) |> 
  filter(n > 1)

most.recent.multiple <- most.recent |> 
  filter(LDCpointID %in% most.recent.multiple.id$LDCpointID)

nrow(most.recent.multiple) / nrow(trt.pre.filtered) # 13%


# Retain only points with a single most recent treatment
most.recent.single <- most.recent |> 
  filter(!PrimaryKey %in% most.recent.multiple$PrimaryKey)

nrow(most.recent.single) / nrow(trt.pre.filtered) # 49%


# Organize cols
most.recent.single <- most.recent.single |> 
  select(
    LDCpointID, ProjectKey, PrimaryKey, DateVisited,
    TrtPolyID, Trt_Type_Major, Trt_Type_Sub, comp_date_est
  ) 



# Create df of control points ---------------------------------------------

# Make df of LDC points that don't overlap with treatment polygons with
#   standardized cols for bind_rows()
ldc.control <- ldc.trt.sjoin |>
  filter(is.na(TrtPolyID)) |>
  select(
    LDCpointID, ProjectKey, PrimaryKey, DateVisited,
    TrtPolyID, Trt_Type_Major, Trt_Type_Sub, comp_date_est
  )


# bind_rows() to make treatment info v006 ---------------------------------

# Combine
treatment.info.006 <- most.recent.single |>
  bind_rows(ldc.control) |>
  arrange(LDCpointID) 

# Check for duplicate primary keys
count(treatment.info.006, PrimaryKey) |> 
  arrange(desc(n))


# Write v005 to CSV -------------------------------------------------------

# Treatment info v005
write_csv(treatment.info.006,
          file = "data/versions-from-R/09.1_Treatment-info_v006.csv",
          na = ""
)


save.image("RData/09.1_treatment-info_v006.RData")

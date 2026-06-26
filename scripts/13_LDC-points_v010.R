# Created: 2026-06-26
# Updated: 2026-06-26

# Purpose: Create LDC points v010, which only includes control points in MLRAs that
#   also contain treated points.

# Script is basically the same as 09.R script from Project v03.


library(tidyverse)

# Load data ---------------------------------------------------------------

trt.mlra.raw <- read_csv("data/GIS-exports/010_TrtPoly009-MLRA-PairwiseIntersect_export.csv")
ldc009.trtpolyid.raw <- read_csv("data/versions-from-R/12_TrtPolyID_v009.csv")
ldc.009.raw <- read_csv("data/versions-from-R/12_LDC-points_v009.csv")


# Data wrangling ----------------------------------------------------------

# LDCpointID & TrtPolyID cols for joining
trtid.join <- ldc009.trtpolyid.raw |>
  select(LDCpointID, TrtPolyID)

# Join with LDC009 and filter for treated only
ldc.009.trt <- ldc.009.raw |>
  left_join(trtid.join) |>
  filter(!Category %in% c("Control, never burned", "Control, post-burn"))


# Create list of MLRAs with treatment polygons based on LDC009
mlra.009 <- trt.mlra.raw |>
  filter(TrtPolyID %in% ldc.009.trt$TrtPolyID) |>
  select(MLRA_NAME) |>
  distinct(.keep_all = TRUE)

#   Inspect to ensure there are only control points
setdiff(ldc.009.raw$MLRA_name, mlra.009$MLRA_NAME)
mlra.ctrl.inspect <- ldc.009.raw |>
  filter(MLRA_name %in% setdiff(ldc.009.raw$MLRA_name, mlra.009$MLRA_NAME))
unique(mlra.ctrl.inspect$Category)


# Compare MLRA list vs. one based just on LDC location (not the pairwise intersect)
mlra.ldc.009 <- ldc.009.raw |>
  filter(!Category %in% c("Control, never burned", "Control, post-burn")) |>
  select(MLRA_name) |>
  distinct(.keep_all = TRUE)

setdiff(mlra.009$MLRA_NAME, mlra.ldc.009$MLRA_name) # Cascade Mountains, Eastern Slope
#   I checked on the map and there are only six points that fall in this MLRA; they are
#     all very close to the border with another MLRA that has a lot of points (Malheur High
#     Plateau. There is just a tiny bit of overlap with an aerial seeding in the
#     Cascade Mountains E Slope MLRA.
#   It is probably easiest just to get rid of these points, because I'm not sure how
#     they will be matched with propensity score matching (I think I would have to manually
#     change the MLRA classification, and I don't to do that).


# Filter to include LDC points only in relevant MLRAs
ldc.010 <- ldc.009.raw |>
  filter(MLRA_name %in% mlra.ldc.009$MLRA_name)


# Write to CSV ------------------------------------------------------------

write_csv(ldc.010,
  file = "data/versions-from-R/13_LDC-points_v010.csv",
  na = ""
)

# Created: 2026-09-22
# Updated: 2026-09-22

# Purpose: Gather a list of the relevant TrtPolyID/Trt_ID values for LDC points v011,
#   which are the points used for PSM.

# Further filtering is not needed based on TrtPolyID, but I want a count of the number of
#   treatment polygons that were included as PSM input.


library(tidyverse)

# Load data ---------------------------------------------------------------

ldc.011 <- read_csv("data/versions-from-R/14.4_LDC-points_v011.csv")
trtpolyid.009 <- read_csv("data/versions-from-R/12_TrtPolyID_v009.csv")

# Data wrangling ----------------------------------------------------------

trtpolyid.011 <- trtpolyid.009 |> 
  filter(PrimaryKey %in% ldc.011$PrimaryKey)

length(unique(trtpolyid.011$TrtPolyID))
length(unique(trtpolyid.011$Trt_ID))

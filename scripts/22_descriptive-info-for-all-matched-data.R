# Created: 2026-10-27
# Updated: 2026-10-27

# Purpose: Calculate some summary statistics and information for all matched data
#   (data used for the 42 models of final analysis).


library(tidyverse)

# Load data ---------------------------------------------------------------

trtpolyid <- read_csv("data/versions-from-R/11.1_TrtPolyID-for-treatment-info-v008.csv")
all.matched <- read_csv("data/versions-from-R/17_all-matched-data_v012.csv")


# Count of LDC points -----------------------------------------------------

length(unique(all.matched$PrimaryKey)) # 10,337


# Origin of LDC points ----------------------------------------------------

count(all.matched, ProjectKey) # all came from BLM AIM


# Count of treatment polygons across all models ---------------------------

all.matched.trt <- trtpolyid |> 
  filter(PrimaryKey %in% all.matched$PrimaryKey)

length(unique(all.matched.trt$TrtPolyID)) # 1,819 treatment polygons

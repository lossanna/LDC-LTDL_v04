# Created: 2026-06-18
# Updated: 2026-06-18

# Purpose: Write out Treatment Info v000, which has fixed the issues created by 
#   overflow text of a single cell (a problem of Excel).

library(tidyverse)

# Load data ---------------------------------------------------------------

treatment.info.raw <- read_csv("data/raw/Treatment_Info_R.csv")
trtid.fix <- read_csv("data/data-wrangling-intermediate/01_treatment-info_fix-rows.csv")

# Data wrangling ----------------------------------------------------------

# Fix rows that need replacement in Treatment_Info_R.csv
treatment.info <- treatment.info.raw |> 
  filter(!Trt_ID %in% trtid.fix$Trt_ID)
treatment.info <- treatment.info |> 
  bind_rows(trtid.fix) |> 
  arrange(Trt_ID)


# Write out fixed CSVs as version 000 -------------------------------------

write_csv(treatment.info,
          file = "data/versions-from-R/01_Treatment-info_v000.csv",
          na = "")

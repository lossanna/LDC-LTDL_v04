# Created: 2026-06-19
# Updated: 2026-06-19

# Purpose: Create the equivalent of a Treatment_info table (GIS join cols only) for the
#   prescribed fires identified in the USGS Combined Wildland Fire Dataset that are
#   missing from the LTDL dataset (same as Project v02).


library(tidyverse)

# Load data ---------------------------------------------------------------

pf.missing.raw <- read_csv("data/GIS-exports/002_Prescribed-fires-missing_export.csv")
treatment.info.001.gisjoin <- read_csv("data/versions-from-R/04_Treatment-info_v001-gisjoin.csv")
treatment.info.001 <- read_csv("data/versions-from-R/04_Treatment-info_v001.csv")


# Remove duplicate rows ---------------------------------------------------

# There are sometimes multiple rows for the same fire because that's how the polygons
#   were split up by the Symmetrical Difference tool during geoprocessing.

# Remove Shape_Length and Shape_Area cols to remove duplicate rows
pf.missing <- pf.missing.raw |>
  select(-Shape_Length, -Shape_Area) |>
  distinct(.keep_all = TRUE)

# Check for matching lengths
length(unique(pf.missing$USGS_Assigned_ID)) == nrow(pf.missing)


# Assign Trt_ID and TrtPolyID ---------------------------------------------

# Figure out range of current Trt_ID from LTDL
range(treatment.info.001$Trt_ID)

# Assign Trt_ID beyond current range
#   Use ID of PrescribedFires002, because that is the ID of original polygons
pf.missing <- pf.missing |>
  mutate(Trt_ID = FID_PrescribedFires002 + 73350)


# Based on TrtPoly001 (checked in ArcGIS), the max TrtPolyID value is 30397

# Assign Trt_ID beyond current range
#   Use ID of PrescribedFires002, because that is the ID of original polygons
pf.missing <- pf.missing |>
  mutate(TrtPolyID = FID_PrescribedFires002 + 30397)


# Assign treatment --------------------------------------------------------

# Major
treatment.info.001 |>
  filter(str_detect(Trt_Type_Major, "Burn")) |>
  select(Trt_Type_Major) |>
  unique()

# Sub
treatment.info.001 |>
  filter(str_detect(Trt_Type_Sub, "Burn")) |>
  select(Trt_Type_Sub) |>
  unique()

# Treatment_Type
treatment.info.001 |>
  filter(str_detect(Treatment_Type, "Burn")) |>
  select(Treatment_Type) |>
  unique()

# Assign Major, Sub, and Treatment_Type
pf.missing <- pf.missing |>
  mutate(
    Trt_Type_Major = "Prescribed Burn",
    Trt_Type_Sub = "Prescribed Burn",
    Treatment_Type = "Prescribed Burn"
  )

# Save intermediate
pf.missing1 <- pf.missing


# Inspect date columns ----------------------------------------------------

# Inspect Listed_Fire_Dates
unique(pf.missing$Listed_Fire_Dates)

# Split into multiple columns
pf.missing <- pf.missing1 |>
  separate_wider_delim(Listed_Fire_Dates,
    delim = " | ",
    names = c("DateA", "DateB", "DateC"),
    too_few = "align_start",
    cols_remove = FALSE
  )

# Split again to separate date vs type of date
pf.missing <- pf.missing |>
  separate_wider_delim(
    cols = c(DateA, DateB, DateC),
    delim = ": ",
    names_sep = ""
  )

# Inspect A, B, and C
unique(pf.missing$DateA1)
unique(pf.missing$DateB1)
unique(pf.missing$DateC1) # do not need listed upload date


# Based on end date -------------------------------------------------------

# Fires with end date in col A
end.date.A <- pf.missing |>
  filter(DateA1 == "Listed Prescribed Fire End Date(s)")

#   Remove other date cols and rename A
end.date.A <- end.date.A |>
  select(-DateA1, -DateB1, -DateB2, -DateC1, -DateC2) |>
  rename(Rx_End_Dates = DateA2)

# Fires with end date in col B
end.date.B <- pf.missing |>
  filter(DateB1 == "Listed Prescribed Fire End Date(s)")

#   Remove other date cols and rename B
end.date.B <- end.date.B |>
  select(-DateA1, -DateA2, -DateB1, -DateC1, -DateC2) |>
  rename(Rx_End_Dates = DateB2)

# Combine A & B
end.date <- end.date.A |>
  bind_rows(end.date.B)

# Fires with a single end date
end.date.single <- end.date |>
  filter(!str_detect(Rx_End_Dates, ",")) |>
  mutate(Rx_End_Date = str_sub(Rx_End_Dates, 1, 10)) |>
  mutate(Rx_End_Date = as.Date(Rx_End_Dates)) |>
  select(-Rx_End_Dates)

#   Add comp_est column
end.date.single <- end.date.single |>
  mutate(comp_est = "actual")


## Fires with multiple end dates ------------------------------------------

# Fires with multiple end dates
end.date.multiple <- end.date |>
  filter(str_detect(Rx_End_Dates, ","))

# Count the number of dates included in the col
end.date.multiple <- end.date.multiple |>
  mutate(Rx_End_Date_count = str_count(Rx_End_Dates, ","))
max(end.date.multiple$Rx_End_Date_count)

# Split dates into multiple columns
end.date.multiple <- end.date.multiple |>
  separate_wider_delim(Rx_End_Dates,
    delim = ", ",
    names = c(
      "RxEndDate1", "RxEndDate2", "RxEndDate3", "RxEndDate4",
      "RxEndDate5", "RxEndDate6", "RxEndDate7"
    ),
    too_few = "align_start",
    cols_remove = TRUE
  )

# Split again to separate count in parentheses
end.date.multiple <- end.date.multiple |>
  separate_wider_delim(RxEndDate1,
    delim = " (",
    names = c("RxEndDate1", "RxEndDate1_count")
  ) |>
  separate_wider_delim(RxEndDate2,
    delim = " (",
    names = c("RxEndDate2", "RxEndDate2_count")
  ) |>
  separate_wider_delim(RxEndDate3,
    delim = " (",
    names = c("RxEndDate3", "RxEndDate3_count")
  ) |>
  separate_wider_delim(RxEndDate4,
    delim = " (",
    names = c("RxEndDate4", "RxEndDate4_count")
  ) |>
  separate_wider_delim(RxEndDate5,
    delim = " (",
    names = c("RxEndDate5", "RxEndDate5_count")
  ) |>
  separate_wider_delim(RxEndDate6,
    delim = " (",
    names = c("RxEndDate6", "RxEndDate6_count")
  ) |>
  separate_wider_delim(RxEndDate7,
    delim = " (",
    names = c("RxEndDate7", "RxEndDate7_count")
  )

# Remove close parentheses
end.date.multiple <- end.date.multiple |>
  mutate(
    RxEndDate1_count = str_sub(RxEndDate1_count, 1, -2),
    RxEndDate2_count = str_sub(RxEndDate2_count, 1, -2),
    RxEndDate3_count = str_sub(RxEndDate3_count, 1, -2),
    RxEndDate4_count = str_sub(RxEndDate4_count, 1, -2),
    RxEndDate5_count = str_sub(RxEndDate5_count, 1, -2),
    RxEndDate6_count = str_sub(RxEndDate6_count, 1, -2),
    RxEndDate7_count = str_sub(RxEndDate7_count, 1, -2)
  )

# pivot_longer()
#   Pivot to name/value
end.date.multiple.pivot <- end.date.multiple |>
  pivot_longer(
    cols = starts_with("RxEndDate"),
    names_to = "name",
    values_to = "value"
  )

#   Parse index and type
end.date.multiple.pivot <- end.date.multiple.pivot |>
  mutate(
    index = str_extract(name, "\\d+"),
    type  = if_else(str_detect(name, "_count$"), "count", "Rx_End_Date")
  )

#   Widen back to paired cols
end.date.multiple.pivot <- end.date.multiple.pivot |>
  select(-name) |>
  pivot_wider(
    names_from = type,
    values_from = value
  ) |>
  filter(!is.na(Rx_End_Date))

# Select date with highest count
end.date.multiple.pivot <- end.date.multiple.pivot |>
  group_by(TrtPolyID) |>
  slice_max(count, n = 1, with_ties = TRUE) |>
  ungroup() |>
  select(TrtPolyID, Rx_End_Date, count) |>
  distinct(.keep_all = TRUE)

# In the case of ties (same count), select most recent date
end.date.multiple.pivot <- end.date.multiple.pivot |>
  mutate(Rx_End_Date = as.Date(Rx_End_Date)) |>
  group_by(TrtPolyID) |>
  filter(count == max(count)) |>
  slice_max(Rx_End_Date, n = 1, with_ties = FALSE) |>
  ungroup()

# Check that all TrtPolyID were retained
setdiff(end.date.multiple$TrtPolyID, end.date.multiple.pivot$TrtPolyID)

# Join with end.date to get all cols
end.date.multiple.filtered <- end.date.multiple.pivot |>
  left_join(end.date) |>
  select(-Rx_End_Dates, -count) |>
  mutate(Rx_End_Date = as.Date(Rx_End_Date))

# Add comp_est column
end.date.multiple.filtered <- end.date.multiple.filtered |>
  mutate(comp_est = "actual (most common & most recent of multiple dates)")


## Combine single and multiple --------------------------------------------

# Combine
end.date.assigned <- end.date.single |>
  bind_rows(end.date.multiple.filtered)

# Check for NAs
apply(end.date.assigned, 2, anyNA)

# Check for matching lengths
nrow(end.date) == nrow(end.date.assigned)

# Check for missing IDs
setdiff(end.date$TrtPolyID, end.date.assigned$TrtPolyID)


## Based on start date ----------------------------------------------------

# Fires with start date (only in col A)
start.date <- pf.missing |>
  filter(DateA1 == "Listed Prescribed Fire Start Date(s)")

#   Remove other date cols and rename start date col
start.date <- start.date |>
  select(-DateA1, -DateB1, -DateB2, -DateC1, -DateC2) |>
  rename(Rx_Start_Dates = DateA2)

# Fires with a single start date
start.date.single <- start.date |>
  filter(!str_detect(Rx_Start_Dates, ",")) |>
  mutate(Rx_Start_Date = str_sub(Rx_Start_Dates, 1, 10)) |>
  mutate(Rx_Start_Date = as.Date(Rx_Start_Dates)) |>
  select(-Rx_Start_Dates)

#   Add init_est col
start.date.single <- start.date.single |>
  mutate(init_est = "actual")


## Fires with multiple start dates ----------------------------------------

# Fires with multiple start dates
start.date.multiple <- start.date |>
  filter(str_detect(Rx_Start_Dates, ","))

# Count the number of dates included in the col
start.date.multiple <- start.date.multiple |>
  mutate(Rx_Start_Date_count = str_count(Rx_Start_Dates, ","))
max(start.date.multiple$Rx_Start_Date_count)

# Split dates into multiple columns
start.date.multiple <- start.date.multiple |>
  separate_wider_delim(Rx_Start_Dates,
    delim = ", ",
    names = c(
      "RxStartDate1", "RxStartDate2", "RxStartDate3", "RxStartDate4",
      "RxStartDate5"
    ),
    too_few = "align_start",
    cols_remove = TRUE
  )

# Split again to separate count in parentheses
start.date.multiple <- start.date.multiple |>
  separate_wider_delim(RxStartDate1,
    delim = " (",
    names = c("RxStartDate1", "RxStartDate1_count")
  ) |>
  separate_wider_delim(RxStartDate2,
    delim = " (",
    names = c("RxStartDate2", "RxStartDate2_count")
  ) |>
  separate_wider_delim(RxStartDate3,
    delim = " (",
    names = c("RxStartDate3", "RxStartDate3_count")
  ) |>
  separate_wider_delim(RxStartDate4,
    delim = " (",
    names = c("RxStartDate4", "RxStartDate4_count")
  ) |>
  separate_wider_delim(RxStartDate5,
    delim = " (",
    names = c("RxStartDate5", "RxStartDate5_count")
  )

# Remove close parentheses
start.date.multiple <- start.date.multiple |>
  mutate(
    RxStartDate1_count = str_sub(RxStartDate1_count, 1, -2),
    RxStartDate2_count = str_sub(RxStartDate2_count, 1, -2),
    RxStartDate3_count = str_sub(RxStartDate3_count, 1, -2),
    RxStartDate4_count = str_sub(RxStartDate4_count, 1, -2),
    RxStartDate5_count = str_sub(RxStartDate5_count, 1, -2)
  )

# pivot_longer()
#   Pivot to name/value
start.date.multiple.pivot <- start.date.multiple |>
  pivot_longer(
    cols = starts_with("RxStartDate"),
    names_to = "name",
    values_to = "value"
  )

#   Parse index and type
start.date.multiple.pivot <- start.date.multiple.pivot |>
  mutate(
    index = str_extract(name, "\\d+"),
    type  = if_else(str_detect(name, "_count$"), "count", "Rx_Start_Date")
  )

#   Widen back to paired cols
start.date.multiple.pivot <- start.date.multiple.pivot |>
  select(-name) |>
  pivot_wider(
    names_from = type,
    values_from = value
  ) |>
  filter(!is.na(Rx_Start_Date))

# Select date with highest count
start.date.multiple.pivot <- start.date.multiple.pivot |>
  group_by(TrtPolyID) |>
  filter(!is.na(count)) |>
  slice_max(count, n = 1, with_ties = TRUE) |>
  ungroup() |>
  select(TrtPolyID, Rx_Start_Date, count) |>
  distinct(.keep_all = TRUE)

# In the case of ties (same count), select most recent date
start.date.multiple.pivot <- start.date.multiple.pivot |>
  mutate(Rx_Start_Date = as.Date(Rx_Start_Date)) |>
  group_by(TrtPolyID) |>
  filter(count == max(count)) |>
  slice_max(Rx_Start_Date, n = 1, with_ties = FALSE) |>
  ungroup()

# Check that all TrtPolyID were retained
setdiff(start.date.multiple$TrtPolyID, start.date.multiple.pivot$TrtPolyID)

# Join with start.date to get all cols
start.date.multiple.filtered <- start.date.multiple.pivot |>
  left_join(start.date) |>
  select(-Rx_Start_Dates, -count) |>
  mutate(Rx_Start_Date = as.Date(Rx_Start_Date))

# Add init_est column
start.date.multiple.filtered <- start.date.multiple.filtered |>
  mutate(init_est = "actual (most common & most recent of multiple dates)")


## Combine single and multiple --------------------------------------------

# Combine
start.date.assigned <- start.date.single |>
  bind_rows(start.date.multiple.filtered)

# Check for NAs
apply(start.date.assigned, 2, anyNA)

# Check for matching lengths
nrow(start.date) == nrow(start.date.assigned)

# Check for missing IDs
setdiff(start.date$TrtPolyID, start.date.assigned$TrtPolyID)


# Assign end date based on start date -------------------------------------

# Calculate fire duration for fires with start & end dates to get estimate
#   Filter and join with end.date.assigned
fire.duration <- start.date.multiple.filtered |>
  filter(TrtPolyID %in% end.date.assigned$TrtPolyID) |>
  left_join(end.date.assigned) |>
  mutate(Fire_duration = as.numeric(Rx_End_Date - Rx_Start_Date))

summary(fire.duration$Fire_duration)
count(fire.duration, Fire_duration) |>
  arrange(desc(n))


# Fires with a start date but no end date
end.date.missing <- start.date.assigned |>
  filter(!TrtPolyID %in% end.date.assigned$TrtPolyID)

# Assign end date as 1 week after start date
end.date.missing.assigned <- end.date.missing |>
  mutate(Rx_End_Date = Rx_Start_Date + weeks(1))

# Add comp_est column
end.date.missing.assigned <- end.date.missing.assigned |>
  mutate(comp_est = "estimated from start date")

# Join with start.date.assigned to get init_est col
end.date.missing.assigned <- end.date.missing.assigned |>
  left_join(start.date.assigned)


# Assign start date based on end date -------------------------------------

# Fires with an end date but no start date
start.date.missing <- end.date.assigned |>
  filter(!TrtPolyID %in% start.date.assigned$TrtPolyID)

# Assign start date as 1 week before end date
start.date.missing.assigned <- start.date.missing |>
  mutate(Rx_Start_Date = Rx_End_Date - weeks(1))

#   Add init_est column
start.date.missing.assigned <- start.date.missing.assigned |>
  mutate(init_est = "estimated from end date")

# Join with end.date.assigned to get comp_est col
start.date.missing.assigned <- start.date.missing.assigned |>
  left_join(end.date.assigned)


# Combine those with dates ------------------------------------------------

# Known start and end dates
start.date.end.date1 <- start.date.assigned |>
  left_join(end.date.assigned) |>
  filter(!is.na(Rx_End_Date))

start.date.end.date2 <- end.date.assigned |>
  left_join(start.date.assigned) |>
  filter(!is.na(Rx_Start_Date))

start.date.end.date <- start.date.end.date1 |>
  bind_rows(start.date.end.date2) |>
  distinct(.keep_all = TRUE)

# Combine with assigned dates
by.date <- start.date.end.date |>
  bind_rows(end.date.missing.assigned) |>
  bind_rows(start.date.missing.assigned)

# Check for NAs
apply(by.date, 2, anyNA)

# Check for duplicate IDs
length(unique(by.date$TrtPolyID)) == nrow(by.date)


# Fires without start or end date -----------------------------------------

# Other fires
missing.date <- pf.missing |>
  filter(!TrtPolyID %in% by.date$TrtPolyID) |>
  select(-DateA1, -DateA2, -DateB1, -DateB2, -DateC1, -DateC2)

# Figure out when most fires are started
start.month <- by.date |>
  mutate(start_month = as.numeric(str_sub(Rx_Start_Date, 6, 7)))
summary(start.month$start_month)
hist(start.month$start_month)

# Assign start date as Sept 1 and add init_est col
missing.date <- missing.date |>
  mutate(
    Rx_Start_Date = as.Date(paste0(Fire_Year, "-09-01")),
    init_est = "estimated from fire year"
  )

# Assign end date as Dec 31 and add comp_est col
missing.date <- missing.date |>
  mutate(
    Rx_End_Date = as.Date(paste0(Fire_Year, "-12-31")),
    comp_est = "estimated from fire year"
  )


# Combine all -------------------------------------------------------------

# Combine
pf.missing.assigned <- by.date |>
  bind_rows(missing.date) |>
  arrange(TrtPolyID)

# Check for matching lengths
nrow(pf.missing) == nrow(pf.missing.assigned)

# Check for missing IDs
setdiff(pf.missing.raw$USGS_Assigned_ID, pf.missing.assigned$USGS_Assigned_ID)


# Create treatment info table equivalent ----------------------------------

# GIS-join cols only
treatment.info.pf <- pf.missing.assigned |>
  mutate(Prj_ID = FID_PrescribedFires002 + max(treatment.info.001$Prj_ID)) |>
  rename(
    init_date_est = Rx_Start_Date,
    comp_date_est = Rx_End_Date
  ) |>
  select(
    Prj_ID, Trt_ID, init_date_est, comp_date_est, Trt_Type_Major, Trt_Type_Sub,
    Treatment_Type, TrtPolyID
  )


# Write to CSV ------------------------------------------------------------

# Prescribed fires as treatment info table
write_csv(treatment.info.pf,
  file = "data/versions-from-R/05_prescribed-fires-added-treatment-info_v002.csv"
)


save.image("RData/05_prescribed-fires-added_v002.RData")

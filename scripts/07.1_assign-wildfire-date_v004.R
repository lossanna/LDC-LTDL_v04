# Created: 2026-06-24
# Updated: 2026-06-24

# Purpose: Assign a date to each wildfire polygon. Same as 06.1.R from Project v03.


library(tidyverse)

# Load data ---------------------------------------------------------------

ldc.wildfire.sj.raw <- read_csv("data/GIS-exports/004_LDC003-Wildfires004-SpatialJoin_export.csv")
overlap.table.raw <- read_csv("data/GIS-exports/004_Overlapping-fires-table_export.csv")


# Separate wildfire polygons ----------------------------------------------

# Fire polygons only
wildfires.all <- ldc.wildfire.sj.raw |>
  filter(!is.na(FirePolyID)) |>
  select(-LDCpointID, -ProjectKey, -PrimaryKey, -DateVisited) |>
  distinct(.keep_all = TRUE)

# Remove fires excluded from summary rasters
wildfires <- wildfires.all |>
  select(
    FirePolyID, USGS_Assigned_ID, Fire_Year, Listed_Fire_Dates, Overlap_Within_1_or_2_Flag,
    Circle_Flag, Exclude_From_Summary_Rasters
  ) |>
  filter(Exclude_From_Summary_Rasters == "No") |> 
  distinct(.keep_all = TRUE)

wildfires1 <- wildfires


# Inspect Overlap_Within_1_or_2_Flag col
unique(wildfires$Overlap_Within_1_or_2_Flag)

# Filter overlap table for overlaps with other wildfires only
unique(overlap.table.raw$Overlapping_Assigned_Fire_Type)
overlap.table <- overlap.table.raw |>
  filter(!Overlapping_Assigned_Fire_Type %in% c("Prescribed Fire", "Unknown - Likely Prescribed Fire")) |>
  filter(USGS_Assigned_ID %in% wildfires$USGS_Assigned_ID)

# Examine overlaps of over 90% and within 2 years of each other
overlap90 <- overlap.table |>
  filter(Percent_Overlap >= 90) |>
  filter(Years_Before_or_After_Fire <= 2 & Years_Before_or_After_Fire >= -2)
count(overlap90, USGS_Assigned_ID) |>
  arrange(desc(n))

# There aren't that many fires with potentially incorrect overlaps (duplicates),
#   so I will just leave them.


# Examine fire dates ------------------------------------------------------

# Calculate the max number of dates in a single cell
wildfires <- wildfires1 |>
  mutate(Listed_Fire_Dates_count = str_count(Listed_Fire_Dates, "Listed"))
max(wildfires$Listed_Fire_Dates_count, na.rm = TRUE)

# Split into multiple columns
wildfires <- wildfires |>
  separate_wider_delim(Listed_Fire_Dates,
    delim = " | ",
    names = c(
      "DateA", "DateB", "DateC", "DateD", "DateE",
      "DateF", "DateG", "DateH", "DateI"
    ),
    too_few = "align_start",
    cols_remove = FALSE
  )

# Split again to separate date vs type of date
wildfires <- wildfires |>
  separate_wider_delim(
    cols = c(
      DateA, DateB, DateC, DateD, DateE, DateF,
      DateG, DateH, DateI
    ),
    delim = ": ",
    names_sep = ""
  )

# Inspect DateA through DateI
unique(wildfires$DateA1)
unique(wildfires$DateB1)
unique(wildfires$DateC1)
unique(wildfires$DateD1)
unique(wildfires$DateE1)
unique(wildfires$DateF1)
unique(wildfires$DateG1)
unique(wildfires$DateH1)
unique(wildfires$DateI1) # do not need upload date


# Count number of date types
unique(c(
  wildfires$DateA1, wildfires$DateB1, wildfires$DateC1, wildfires$DateD1,
  wildfires$DateE1, wildfires$DateF1, wildfires$DateG1, wildfires$DateH1
))

date.count <- wildfires1 |>
  mutate(
    Ignition_count = str_count(Listed_Fire_Dates, "Ignition"),
    Discovery_count = str_count(Listed_Fire_Dates, "Discovery"),
    Control_count = str_count(Listed_Fire_Dates, "Controlled"),
    Containment_count = str_count(Listed_Fire_Dates, "Containment"),
    Other_count = str_count(Listed_Fire_Dates, "Other"),
    Out_count = str_count(Listed_Fire_Dates, "Out")
  ) 

sum(date.count$Ignition_count, na.rm = TRUE)
sum(date.count$Ignition_count, na.rm = TRUE) / nrow(wildfires1) # 60%

sum(date.count$Discovery_count, na.rm = TRUE)
sum(date.count$Discovery_count, na.rm = TRUE) / nrow(wildfires1) # 86%

sum(date.count$Control_count, na.rm = TRUE)
sum(date.count$Control_count, na.rm = TRUE) / nrow(wildfires1) # 68%

sum(date.count$Containment_count, na.rm = TRUE)
sum(date.count$Containment_count, na.rm = TRUE) / nrow(wildfires1) # 49%

sum(date.count$Other_count, na.rm = TRUE)
sum(date.count$Other_count, na.rm = TRUE) / nrow(wildfires1) # 46%

sum(date.count$Out_count, na.rm = TRUE)
sum(date.count$Out_count, na.rm = TRUE) / nrow(wildfires1) # 1%


# Based on Wildfire Out date ----------------------------------------------

# Fires with out date in col B
out.date.B <- wildfires |>
  filter(DateB1 == "Listed Wildfire Out Date(s)")

#   Remove other date cols and rename B
out.date.B <- out.date.B |>
  select(FirePolyID, USGS_Assigned_ID, Fire_Year, DateB2) |>
  rename(Wildfire_Dates = DateB2)


# Fires with end date in col C
out.date.C <- wildfires |>
  filter(DateC1 == "Listed Wildfire Out Date(s)")

#   Remove other date cols and rename C
out.date.C <- out.date.C |>
  select(FirePolyID, USGS_Assigned_ID, Fire_Year, DateC2) |>
  rename(Wildfire_Dates = DateC2)


# Fires with end date in col D
out.date.D <- wildfires |>
  filter(DateD1 == "Listed Wildfire Out Date(s)")

#   Remove other date cols and rename D
out.date.D <- out.date.D |>
  select(FirePolyID, USGS_Assigned_ID, Fire_Year, DateD2) |>
  rename(Wildfire_Dates = DateD2)


# Fires with end date in col E
out.date.E <- wildfires |>
  filter(DateE1 == "Listed Wildfire Out Date(s)")

#   Remove other date cols and rename E
out.date.E <- out.date.E |>
  select(FirePolyID, USGS_Assigned_ID, Fire_Year, DateE2) |>
  rename(Wildfire_Dates = DateE2)


# Combine all
out.date <- bind_rows(out.date.B, out.date.C, out.date.D, out.date.E)

# Add Date_type col
out.date <- out.date |>
  mutate(Date_type = "Wildfire Out")


## Single date ------------------------------------------------------------

# Fires with single date
out.date.single <- out.date |>
  filter(!str_detect(Wildfire_Dates, ",")) |>
  mutate(Wildfire_Date = str_sub(Wildfire_Dates, 1, 10)) |>
  mutate(Wildfire_Date = as.Date(Wildfire_Date)) |>
  select(-Wildfire_Dates)

# Add Date_est column
out.date.single <- out.date.single |>
  mutate(Date_est = "actual")


## Multiple dates ---------------------------------------------------------

# Fires with multiple end dates
out.date.multiple <- out.date |>
  filter(str_detect(Wildfire_Dates, ","))

# Count the number of dates included in the col
out.date.multiple <- out.date.multiple |>
  mutate(Wildfire_Date_count = str_count(Wildfire_Dates, ","))
max(out.date.multiple$Wildfire_Date_count)

# Split dates into multiple columns
out.date.multiple <- out.date.multiple |>
  separate_wider_delim(Wildfire_Dates,
    delim = ", ",
    names = c("OutDate1", "OutDate2", "OutDate3", "OutDate4"),
    too_few = "align_start",
    cols_remove = TRUE
  )

# Split again to separate count in parentheses
out.date.multiple <- out.date.multiple |>
  separate_wider_delim(OutDate1,
    delim = " (",
    names = c("OutDate1", "OutDate1_count")
  ) |>
  separate_wider_delim(OutDate2,
    delim = " (",
    names = c("OutDate2", "OutDate2_count")
  ) |>
  separate_wider_delim(OutDate3,
    delim = " (",
    names = c("OutDate3", "OutDate3_count")
  ) |>
  separate_wider_delim(OutDate4,
    delim = " (",
    names = c("OutDate4", "OutDate4_count")
  )

# Remove close parentheses
out.date.multiple <- out.date.multiple |>
  mutate(
    OutDate1_count = str_sub(OutDate1_count, 1, -2),
    OutDate2_count = str_sub(OutDate2_count, 1, -2),
    OutDate3_count = str_sub(OutDate3_count, 1, -2),
    OutDate4_count = str_sub(OutDate4_count, 1, -2)
  )

# pivot_longer()
#   Pivot to name/value
out.date.multiple.pivot <- out.date.multiple |>
  pivot_longer(
    cols = starts_with("OutDate"),
    names_to = "name",
    values_to = "value"
  ) |> 
  distinct(.keep_all = TRUE)

#   Parse index and type
out.date.multiple.pivot <- out.date.multiple.pivot |>
  mutate(
    index = str_extract(name, "\\d+"),
    type  = if_else(str_detect(name, "_count$"), "count", "Wildfire_Date")
  )

#   Widen back to paired cols
out.date.multiple.pivot <- out.date.multiple.pivot |>
  select(-name) |>
  pivot_wider(
    names_from = type,
    values_from = value
  ) |>
  filter(!is.na(Wildfire_Date))

# Select date with highest count
out.date.multiple.pivot <- out.date.multiple.pivot |>
  group_by(FirePolyID) |>
  slice_max(count, n = 1, with_ties = TRUE) |>
  ungroup() |>
  select(FirePolyID, Wildfire_Date, count) |>
  distinct(.keep_all = TRUE)

# In the case of ties (same count), select most recent date
out.date.multiple.pivot <- out.date.multiple.pivot |>
  mutate(Wildfire_Date = as.Date(Wildfire_Date)) |>
  group_by(FirePolyID) |>
  filter(count == max(count)) |>
  slice_max(Wildfire_Date, n = 1, with_ties = FALSE) |>
  ungroup()

# Check that all FirePolyID were retained
setdiff(out.date.multiple$FirePolyID, out.date.multiple.pivot$FirePolyID)

# Join with out.date to get all cols
out.date.multiple.filtered <- out.date.multiple.pivot |>
  left_join(out.date) |>
  select(-Wildfire_Dates, -count) |>
  mutate(Wildfire_Date = as.Date(Wildfire_Date))

# Add Date_est column
out.date.multiple.filtered <- out.date.multiple.filtered |>
  mutate(Date_est = "actual (most common & most recent of multiple dates)")


## Combine single and multiple --------------------------------------------

# Combine
out.date.assigned <- out.date.single |>
  bind_rows(out.date.multiple.filtered)

# Check for NAs
apply(out.date.assigned, 2, anyNA)

# Check for matching lengths
nrow(out.date) == nrow(out.date.assigned)

# Check for missing IDs
setdiff(out.date$FirePolyID, out.date.assigned$FirePolyID)


# Based on Wildfire Containment date --------------------------------------

# Fires with containment date in col A
contain.date.A <- wildfires |>
  filter(DateA1 == "Listed Wildfire Containment Date(s)")

#   Remove other date cols and rename A
contain.date.A <- contain.date.A |>
  select(FirePolyID, USGS_Assigned_ID, Fire_Year, DateA2) |>
  rename(Wildfire_Dates = DateA2)


# Fires with containment date in col B
contain.date.B <- wildfires |>
  filter(DateB1 == "Listed Wildfire Containment Date(s)")

#   Remove other date cols and rename B
contain.date.B <- contain.date.B |>
  select(FirePolyID, USGS_Assigned_ID, Fire_Year, DateB2) |>
  rename(Wildfire_Dates = DateB2)


# Fires with containment date in col C
contain.date.C <- wildfires |>
  filter(DateC1 == "Listed Wildfire Containment Date(s)")

#   Remove other date cols and rename C
contain.date.C <- contain.date.C |>
  select(FirePolyID, USGS_Assigned_ID, Fire_Year, DateC2) |>
  rename(Wildfire_Dates = DateC2)


# Fires with containment date in col D
contain.date.D <- wildfires |>
  filter(DateD1 == "Listed Wildfire Containment Date(s)")

#   Remove other date cols and rename D
contain.date.D <- contain.date.D |>
  select(FirePolyID, USGS_Assigned_ID, Fire_Year, DateD2) |>
  rename(Wildfire_Dates = DateD2)


# Combine all
contain.date <- bind_rows(contain.date.A, contain.date.B, contain.date.C, contain.date.D)

# Add Date_type col
contain.date <- contain.date |>
  mutate(Date_type = "Wildfire Containment")


## Single date ------------------------------------------------------------

# Fires with single date
contain.date.single <- contain.date |>
  filter(!str_detect(Wildfire_Dates, ",")) |>
  mutate(Wildfire_Date = str_sub(Wildfire_Dates, 1, 10)) |>
  mutate(Wildfire_Date = as.Date(Wildfire_Date)) |>
  select(-Wildfire_Dates)

# Add Date_est column
contain.date.single <- contain.date.single |>
  mutate(Date_est = "actual")


## Multiple dates ---------------------------------------------------------

# Fires with multiple end dates
contain.date.multiple <- contain.date |>
  filter(str_detect(Wildfire_Dates, ","))

# Count the number of dates included in the col
contain.date.multiple <- contain.date.multiple |>
  mutate(Wildfire_Date_count = str_count(Wildfire_Dates, ","))
max(contain.date.multiple$Wildfire_Date_count)

# Split dates into multiple columns
contain.date.multiple <- contain.date.multiple |>
  separate_wider_delim(Wildfire_Dates,
    delim = ", ",
    names = c(
      "ContainDate1", "ContainDate2", "ContainDate3", "ContainDate4",
      "ContainDate5", "ContainDate6", "ContainDate7", "ContainDate8"
    ),
    too_few = "align_start",
    cols_remove = TRUE
  )

# Split again to separate count in parentheses
contain.date.multiple <- contain.date.multiple |>
  separate_wider_delim(ContainDate1,
    delim = " (",
    names = c("ContainDate1", "ContainDate1_count")
  ) |>
  separate_wider_delim(ContainDate2,
    delim = " (",
    names = c("ContainDate2", "ContainDate2_count")
  ) |>
  separate_wider_delim(ContainDate3,
    delim = " (",
    names = c("ContainDate3", "ContainDate3_count")
  ) |>
  separate_wider_delim(ContainDate4,
    delim = " (",
    names = c("ContainDate4", "ContainDate4_count")
  ) |>
  separate_wider_delim(ContainDate5,
    delim = " (",
    names = c("ContainDate5", "ContainDate5_count")
  ) |>
  separate_wider_delim(ContainDate6,
    delim = " (",
    names = c("ContainDate6", "ContainDate6_count")
  ) |>
  separate_wider_delim(ContainDate7,
    delim = " (",
    names = c("ContainDate7", "ContainDate7_count")
  ) |>
  separate_wider_delim(ContainDate8,
    delim = " (",
    names = c("ContainDate8", "ContainDate8_count")
  )

# Remove close parentheses
contain.date.multiple <- contain.date.multiple |>
  mutate(
    ContainDate1_count = str_sub(ContainDate1_count, 1, -2),
    ContainDate2_count = str_sub(ContainDate2_count, 1, -2),
    ContainDate3_count = str_sub(ContainDate3_count, 1, -2),
    ContainDate4_count = str_sub(ContainDate4_count, 1, -2),
    ContainDate5_count = str_sub(ContainDate5_count, 1, -2),
    ContainDate6_count = str_sub(ContainDate6_count, 1, -2),
    ContainDate7_count = str_sub(ContainDate7_count, 1, -2),
    ContainDate8_count = str_sub(ContainDate8_count, 1, -2)
  )

# pivot_longer()
#   Pivot to name/value
contain.date.multiple.pivot <- contain.date.multiple |>
  pivot_longer(
    cols = starts_with("ContainDate"),
    names_to = "name",
    values_to = "value"
  )

#   Parse index and type
contain.date.multiple.pivot <- contain.date.multiple.pivot |>
  mutate(
    index = str_extract(name, "\\d+"),
    type  = if_else(str_detect(name, "_count$"), "count", "Wildfire_Date")
  )

#   Widen back to paired cols
contain.date.multiple.pivot <- contain.date.multiple.pivot |>
  select(-name) |>
  pivot_wider(
    names_from = type,
    values_from = value
  ) |>
  filter(!is.na(Wildfire_Date))

# Select date with highest count
contain.date.multiple.pivot <- contain.date.multiple.pivot |>
  group_by(FirePolyID) |>
  slice_max(count, n = 1, with_ties = TRUE) |>
  ungroup() |>
  select(FirePolyID, Wildfire_Date, count) |>
  distinct(.keep_all = TRUE)

# In the case of ties (same count), select most recent date
contain.date.multiple.pivot <- contain.date.multiple.pivot |>
  mutate(Wildfire_Date = as.Date(Wildfire_Date)) |>
  group_by(FirePolyID) |>
  filter(count == max(count)) |>
  slice_max(Wildfire_Date, n = 1, with_ties = FALSE) |>
  ungroup()

# Check that all FirePolyID were retained
setdiff(contain.date.multiple$FirePolyID, contain.date.multiple.pivot$FirePolyID)

# Join with contain.date to get all cols
contain.date.multiple.filtered <- contain.date.multiple.pivot |>
  left_join(contain.date) |>
  select(-Wildfire_Dates, -count) |>
  mutate(Wildfire_Date = as.Date(Wildfire_Date))

# Add Date_est column
contain.date.multiple.filtered <- contain.date.multiple.filtered |>
  mutate(Date_est = "actual (most common & most recent of multiple dates)")


## Combine single and multiple --------------------------------------------

# Combine
contain.date.assigned <- contain.date.single |>
  bind_rows(contain.date.multiple.filtered)

# Check for NAs
apply(contain.date.assigned, 2, anyNA)

# Check for matching lengths
nrow(contain.date) == nrow(contain.date.assigned)

# Check for missing IDs
setdiff(contain.date$FirePolyID, contain.date.assigned$FirePolyID)

# Remove fires already assigned
contain.date.filtered <- contain.date.assigned |>
  filter(!FirePolyID %in% out.date.assigned$FirePolyID)


# Based on Wildfire Controlled date ---------------------------------------

# Fires with controlled date in col A
control.date.A <- wildfires |>
  filter(DateA1 == "Listed Wildfire Controlled Date(s)")

#   Remove other date cols and rename A
control.date.A <- control.date.A |>
  select(FirePolyID, USGS_Assigned_ID, Fire_Year, DateA2) |>
  rename(Wildfire_Dates = DateA2)


# Fires with controlled date in col B
control.date.B <- wildfires |>
  filter(DateB1 == "Listed Wildfire Controlled Date(s)")

#   Remove other date cols and rename B
control.date.B <- control.date.B |>
  select(FirePolyID, USGS_Assigned_ID, Fire_Year, DateB2) |>
  rename(Wildfire_Dates = DateB2)


# Fires with controlled date in col C
control.date.C <- wildfires |>
  filter(DateC1 == "Listed Wildfire Controlled Date(s)")

#   Remove other date cols and rename C
control.date.C <- control.date.C |>
  select(FirePolyID, USGS_Assigned_ID, Fire_Year, DateC2) |>
  rename(Wildfire_Dates = DateC2)


# Combine all
control.date <- bind_rows(control.date.A, control.date.B, control.date.C)

# Add Date_type col
control.date <- control.date |>
  mutate(Date_type = "Wildfire Controlled")


## Single date ------------------------------------------------------------

# Fires with single date
control.date.single <- control.date |>
  filter(!str_detect(Wildfire_Dates, ",")) |>
  mutate(Wildfire_Date = str_sub(Wildfire_Dates, 1, 10)) |>
  mutate(Wildfire_Date = as.Date(Wildfire_Date)) |>
  select(-Wildfire_Dates)

# Add Date_est column
control.date.single <- control.date.single |>
  mutate(Date_est = "actual")


## Multiple dates ---------------------------------------------------------

# Fires with multiple end dates
control.date.multiple <- control.date |>
  filter(str_detect(Wildfire_Dates, ","))

# Count the number of dates included in the col
control.date.multiple <- control.date.multiple |>
  mutate(Wildfire_Date_count = str_count(Wildfire_Dates, ","))
max(control.date.multiple$Wildfire_Date_count)

# Split dates into multiple columns
control.date.multiple <- control.date.multiple |>
  separate_wider_delim(Wildfire_Dates,
    delim = ", ",
    names = c(
      "ControlDate1", "ControlDate2", "ControlDate3", "ControlDate4",
      "ControlDate5", "ControlDate6", "ControlDate7", "ControlDate8",
      "ControlDate9", "ControlDate10", "ControlDate11"
    ),
    too_few = "align_start",
    cols_remove = TRUE
  )

# Split again to separate count in parentheses
control.date.multiple <- control.date.multiple |>
  separate_wider_delim(ControlDate1,
    delim = " (",
    names = c("ControlDate1", "ControlDate1_count")
  ) |>
  separate_wider_delim(ControlDate2,
    delim = " (",
    names = c("ControlDate2", "ControlDate2_count")
  ) |>
  separate_wider_delim(ControlDate3,
    delim = " (",
    names = c("ControlDate3", "ControlDate3_count")
  ) |>
  separate_wider_delim(ControlDate4,
    delim = " (",
    names = c("ControlDate4", "ControlDate4_count")
  ) |>
  separate_wider_delim(ControlDate5,
    delim = " (",
    names = c("ControlDate5", "ControlDate5_count")
  ) |>
  separate_wider_delim(ControlDate6,
    delim = " (",
    names = c("ControlDate6", "ControlDate6_count")
  ) |>
  separate_wider_delim(ControlDate7,
    delim = " (",
    names = c("ControlDate7", "ControlDate7_count")
  ) |>
  separate_wider_delim(ControlDate8,
    delim = " (",
    names = c("ControlDate8", "ControlDate8_count")
  ) |>
  separate_wider_delim(ControlDate9,
    delim = " (",
    names = c("ControlDate9", "ControlDate9_count")
  ) |>
  separate_wider_delim(ControlDate10,
    delim = " (",
    names = c("ControlDate10", "ControlDate10_count")
  ) |>
  separate_wider_delim(ControlDate11,
    delim = " (",
    names = c("ControlDate11", "ControlDate11_count")
  )

# Remove close parentheses
control.date.multiple <- control.date.multiple |>
  mutate(
    ControlDate1_count = str_sub(ControlDate1_count, 1, -2),
    ControlDate2_count = str_sub(ControlDate2_count, 1, -2),
    ControlDate3_count = str_sub(ControlDate3_count, 1, -2),
    ControlDate4_count = str_sub(ControlDate4_count, 1, -2),
    ControlDate5_count = str_sub(ControlDate5_count, 1, -2),
    ControlDate6_count = str_sub(ControlDate6_count, 1, -2),
    ControlDate7_count = str_sub(ControlDate7_count, 1, -2),
    ControlDate8_count = str_sub(ControlDate8_count, 1, -2),
    ControlDate9_count = str_sub(ControlDate9_count, 1, -2),
    ControlDate10_count = str_sub(ControlDate10_count, 1, -2),
    ControlDate11_count = str_sub(ControlDate11_count, 1, -2)
  )

# pivot_longer()
#   Pivot to name/value
control.date.multiple.pivot <- control.date.multiple |>
  pivot_longer(
    cols = starts_with("ControlDate"),
    names_to = "name",
    values_to = "value"
  )

#   Parse index and type
control.date.multiple.pivot <- control.date.multiple.pivot |>
  mutate(
    index = str_extract(name, "\\d+"),
    type  = if_else(str_detect(name, "_count$"), "count", "Wildfire_Date")
  )

#   Widen back to paired cols
control.date.multiple.pivot <- control.date.multiple.pivot |>
  select(-name) |>
  pivot_wider(
    names_from = type,
    values_from = value
  ) |>
  filter(!is.na(Wildfire_Date))

# Select date with highest count
control.date.multiple.pivot <- control.date.multiple.pivot |>
  group_by(FirePolyID) |>
  slice_max(count, n = 1, with_ties = TRUE) |>
  ungroup() |>
  select(FirePolyID, Wildfire_Date, count) |>
  distinct(.keep_all = TRUE)

# In the case of ties (same count), select most recent date
control.date.multiple.pivot <- control.date.multiple.pivot |>
  mutate(Wildfire_Date = as.Date(Wildfire_Date)) |>
  group_by(FirePolyID) |>
  filter(count == max(count)) |>
  slice_max(Wildfire_Date, n = 1, with_ties = FALSE) |>
  ungroup()

# Check that all FirePolyID were retained
setdiff(control.date.multiple$FirePolyID, control.date.multiple.pivot$FirePolyID)

# Join with control.date to get all cols
control.date.multiple.filtered <- control.date.multiple.pivot |>
  left_join(control.date) |>
  select(-Wildfire_Dates, -count) |>
  mutate(Wildfire_Date = as.Date(Wildfire_Date))

# Add Date_est column
control.date.multiple.filtered <- control.date.multiple.filtered |>
  mutate(Date_est = "actual (most common & most recent of multiple dates)")


## Combine single and multiple --------------------------------------------

# Combine
control.date.assigned <- control.date.single |>
  bind_rows(control.date.multiple.filtered)

# Check for NAs
apply(control.date.assigned, 2, anyNA)

# Check for matching lengths
nrow(control.date) == nrow(control.date.assigned)

# Check for missing IDs
setdiff(control.date$FirePolyID, control.date.assigned$FirePolyID)

# Remove fires already assigned
control.date.filtered <- control.date.assigned |>
  filter(!FirePolyID %in% out.date.assigned$FirePolyID) |>
  filter(!FirePolyID %in% contain.date.filtered$FirePolyID)


# Based on Wildfire Discovery date ----------------------------------------

# Fires with discovery date in col A
discovery.date.A <- wildfires |>
  filter(DateA1 == "Listed Wildfire Discovery Date(s)")

#   Remove other date cols and rename A
discovery.date.A <- discovery.date.A |>
  select(FirePolyID, USGS_Assigned_ID, Fire_Year, DateA2) |>
  rename(Wildfire_Dates = DateA2)


# Fires with discovery date in col B
discovery.date.B <- wildfires |>
  filter(DateB1 == "Listed Wildfire Discovery Date(s)")

#   Remove other date cols and rename B
discovery.date.B <- discovery.date.B |>
  select(FirePolyID, USGS_Assigned_ID, Fire_Year, DateB2) |>
  rename(Wildfire_Dates = DateB2)


# Combine all
discovery.date <- bind_rows(discovery.date.A, discovery.date.B)

# Add Date_type col
discovery.date <- discovery.date |>
  mutate(Date_type = "Wildfire Discovery")


## Single date ------------------------------------------------------------

# Fires with single date
discovery.date.single <- discovery.date |>
  filter(!str_detect(Wildfire_Dates, ",")) |>
  mutate(Wildfire_Date = str_sub(Wildfire_Dates, 1, 10)) |>
  mutate(Wildfire_Date = as.Date(Wildfire_Date)) |>
  select(-Wildfire_Dates)

# Add Date_est column
discovery.date.single <- discovery.date.single |>
  mutate(Date_est = "actual")


## Multiple dates ---------------------------------------------------------

# Fires with multiple end dates
discovery.date.multiple <- discovery.date |>
  filter(str_detect(Wildfire_Dates, ","))

# Count the number of dates included in the col
discovery.date.multiple <- discovery.date.multiple |>
  mutate(Wildfire_Date_count = str_count(Wildfire_Dates, ","))
max(discovery.date.multiple$Wildfire_Date_count)

# Split dates into multiple columns
discovery.date.multiple <- discovery.date.multiple |>
  separate_wider_delim(Wildfire_Dates,
    delim = ", ",
    names = c(
      "DiscoveryDate1", "DiscoveryDate2", "DiscoveryDate3", "DiscoveryDate4",
      "DiscoveryDate5", "DiscoveryDate6", "DiscoveryDate7", "DiscoveryDate8",
      "DiscoveryDate9"
    ),
    too_few = "align_start",
    cols_remove = TRUE
  )

# Split again to separate count in parentheses
discovery.date.multiple <- discovery.date.multiple |>
  separate_wider_delim(DiscoveryDate1,
    delim = " (",
    names = c("DiscoveryDate1", "DiscoveryDate1_count")
  ) |>
  separate_wider_delim(DiscoveryDate2,
    delim = " (",
    names = c("DiscoveryDate2", "DiscoveryDate2_count")
  ) |>
  separate_wider_delim(DiscoveryDate3,
    delim = " (",
    names = c("DiscoveryDate3", "DiscoveryDate3_count")
  ) |>
  separate_wider_delim(DiscoveryDate4,
    delim = " (",
    names = c("DiscoveryDate4", "DiscoveryDate4_count")
  ) |>
  separate_wider_delim(DiscoveryDate5,
    delim = " (",
    names = c("DiscoveryDate5", "DiscoveryDate5_count")
  ) |>
  separate_wider_delim(DiscoveryDate6,
    delim = " (",
    names = c("DiscoveryDate6", "DiscoveryDate6_count")
  ) |>
  separate_wider_delim(DiscoveryDate7,
    delim = " (",
    names = c("DiscoveryDate7", "DiscoveryDate7_count")
  ) |>
  separate_wider_delim(DiscoveryDate8,
    delim = " (",
    names = c("DiscoveryDate8", "DiscoveryDate8_count")
  ) |>
  separate_wider_delim(DiscoveryDate9,
    delim = " (",
    names = c("DiscoveryDate9", "DiscoveryDate9_count")
  )

# Remove close parentheses
discovery.date.multiple <- discovery.date.multiple |>
  mutate(
    DiscoveryDate1_count = str_sub(DiscoveryDate1_count, 1, -2),
    DiscoveryDate2_count = str_sub(DiscoveryDate2_count, 1, -2),
    DiscoveryDate3_count = str_sub(DiscoveryDate3_count, 1, -2),
    DiscoveryDate4_count = str_sub(DiscoveryDate4_count, 1, -2),
    DiscoveryDate5_count = str_sub(DiscoveryDate5_count, 1, -2),
    DiscoveryDate6_count = str_sub(DiscoveryDate6_count, 1, -2),
    DiscoveryDate7_count = str_sub(DiscoveryDate7_count, 1, -2),
    DiscoveryDate8_count = str_sub(DiscoveryDate8_count, 1, -2),
    DiscoveryDate9_count = str_sub(DiscoveryDate9_count, 1, -2)
  )

# pivot_longer()
#   Pivot to name/value
discovery.date.multiple.pivot <- discovery.date.multiple |>
  pivot_longer(
    cols = starts_with("DiscoveryDate"),
    names_to = "name",
    values_to = "value"
  )

#   Parse index and type
discovery.date.multiple.pivot <- discovery.date.multiple.pivot |>
  mutate(
    index = str_extract(name, "\\d+"),
    type  = if_else(str_detect(name, "_count$"), "count", "Wildfire_Date")
  )

#   Widen back to paired cols
discovery.date.multiple.pivot <- discovery.date.multiple.pivot |>
  select(-name) |>
  pivot_wider(
    names_from = type,
    values_from = value
  ) |>
  filter(!is.na(Wildfire_Date))

# Select date with highest count
discovery.date.multiple.pivot <- discovery.date.multiple.pivot |>
  group_by(FirePolyID) |>
  slice_max(count, n = 1, with_ties = TRUE) |>
  ungroup() |>
  select(FirePolyID, Wildfire_Date, count) |>
  distinct(.keep_all = TRUE)

# In the case of ties (same count), select most recent date
discovery.date.multiple.pivot <- discovery.date.multiple.pivot |>
  mutate(Wildfire_Date = as.Date(Wildfire_Date)) |>
  group_by(FirePolyID) |>
  filter(count == max(count)) |>
  slice_max(Wildfire_Date, n = 1, with_ties = FALSE) |>
  ungroup()

# Check that all FirePolyID were retained
setdiff(discovery.date.multiple$FirePolyID, discovery.date.multiple.pivot$FirePolyID)

# Join with discovery.date to get all cols
discovery.date.multiple.filtered <- discovery.date.multiple.pivot |>
  left_join(discovery.date) |>
  select(-Wildfire_Dates, -count) |>
  mutate(Wildfire_Date = as.Date(Wildfire_Date))

# Add Date_est column
discovery.date.multiple.filtered <- discovery.date.multiple.filtered |>
  mutate(Date_est = "actual (most common & most recent of multiple dates)")


## Combine single and multiple --------------------------------------------

# Combine
discovery.date.assigned <- discovery.date.single |>
  bind_rows(discovery.date.multiple.filtered)

# Check for NAs
apply(discovery.date.assigned, 2, anyNA)

# Check for matching lengths
nrow(discovery.date) == nrow(discovery.date.assigned)

# Check for missing IDs
setdiff(discovery.date$FirePolyID, discovery.date.assigned$FirePolyID)

# Remove fires already assigned
discovery.date.filtered <- discovery.date.assigned |>
  filter(!FirePolyID %in% out.date.assigned$FirePolyID) |>
  filter(!FirePolyID %in% contain.date.filtered$FirePolyID) |>
  filter(!FirePolyID %in% control.date.filtered$FirePolyID)


# Based on Ignition date --------------------------------------------------

# Fires with ignition date in col A
ignition.date <- wildfires |>
  filter(DateA1 == "Listed Ignition Date(s)")

#   Remove other date cols and rename A
ignition.date <- ignition.date |>
  select(FirePolyID, USGS_Assigned_ID, Fire_Year, DateA2) |>
  rename(Wildfire_Dates = DateA2)

# Add Date_type col
ignition.date <- ignition.date |>
  mutate(Date_type = "Ignition")


## Single date ------------------------------------------------------------

# Fires with single date
ignition.date.single <- ignition.date |>
  filter(!str_detect(Wildfire_Dates, ",")) |>
  mutate(Wildfire_Date = str_sub(Wildfire_Dates, 1, 10)) |>
  mutate(Wildfire_Date = as.Date(Wildfire_Date)) |>
  select(-Wildfire_Dates)

# Add Date_est column
ignition.date.single <- ignition.date.single |>
  mutate(Date_est = "actual")


## Multiple dates ---------------------------------------------------------

# Fires with multiple end dates
ignition.date.multiple <- ignition.date |>
  filter(str_detect(Wildfire_Dates, ","))

# Count the number of dates included in the col
ignition.date.multiple <- ignition.date.multiple |>
  mutate(Wildfire_Date_count = str_count(Wildfire_Dates, ","))
max(ignition.date.multiple$Wildfire_Date_count)

# Split dates into multiple columns
ignition.date.multiple <- ignition.date.multiple |>
  separate_wider_delim(Wildfire_Dates,
    delim = ", ",
    names = c(
      "IgnitionDate1", "IgnitionDate2", "IgnitionDate3", "IgnitionDate4",
      "IgnitionDate5", "IgnitionDate6"
    ),
    too_few = "align_start",
    cols_remove = TRUE
  )

# Split again to separate count in parentheses
ignition.date.multiple <- ignition.date.multiple |>
  separate_wider_delim(IgnitionDate1,
    delim = " (",
    names = c("IgnitionDate1", "IgnitionDate1_count")
  ) |>
  separate_wider_delim(IgnitionDate2,
    delim = " (",
    names = c("IgnitionDate2", "IgnitionDate2_count")
  ) |>
  separate_wider_delim(IgnitionDate3,
    delim = " (",
    names = c("IgnitionDate3", "IgnitionDate3_count")
  ) |>
  separate_wider_delim(IgnitionDate4,
    delim = " (",
    names = c("IgnitionDate4", "IgnitionDate4_count")
  ) |>
  separate_wider_delim(IgnitionDate5,
    delim = " (",
    names = c("IgnitionDate5", "IgnitionDate5_count")
  ) |>
  separate_wider_delim(IgnitionDate6,
    delim = " (",
    names = c("IgnitionDate6", "IgnitionDate6_count")
  )

# Remove close parentheses
ignition.date.multiple <- ignition.date.multiple |>
  mutate(
    IgnitionDate1_count = str_sub(IgnitionDate1_count, 1, -2),
    IgnitionDate2_count = str_sub(IgnitionDate2_count, 1, -2),
    IgnitionDate3_count = str_sub(IgnitionDate3_count, 1, -2),
    IgnitionDate4_count = str_sub(IgnitionDate4_count, 1, -2),
    IgnitionDate5_count = str_sub(IgnitionDate5_count, 1, -2),
    IgnitionDate6_count = str_sub(IgnitionDate6_count, 1, -2)
  )

# pivot_longer()
#   Pivot to name/value
ignition.date.multiple.pivot <- ignition.date.multiple |>
  pivot_longer(
    cols = starts_with("IgnitionDate"),
    names_to = "name",
    values_to = "value"
  )

#   Parse index and type
ignition.date.multiple.pivot <- ignition.date.multiple.pivot |>
  mutate(
    index = str_extract(name, "\\d+"),
    type  = if_else(str_detect(name, "_count$"), "count", "Wildfire_Date")
  )

#   Widen back to paired cols
ignition.date.multiple.pivot <- ignition.date.multiple.pivot |>
  select(-name) |>
  pivot_wider(
    names_from = type,
    values_from = value
  ) |>
  filter(!is.na(Wildfire_Date))

# Select date with highest count
ignition.date.multiple.pivot <- ignition.date.multiple.pivot |>
  group_by(FirePolyID) |>
  slice_max(count, n = 1, with_ties = TRUE) |>
  ungroup() |>
  select(FirePolyID, Wildfire_Date, count) |>
  distinct(.keep_all = TRUE)

# In the case of ties (same count), select most recent date
ignition.date.multiple.pivot <- ignition.date.multiple.pivot |>
  mutate(Wildfire_Date = as.Date(Wildfire_Date)) |>
  group_by(FirePolyID) |>
  filter(count == max(count)) |>
  slice_max(Wildfire_Date, n = 1, with_ties = FALSE) |>
  ungroup()

# Check that all FirePolyID were retained
setdiff(ignition.date.multiple$FirePolyID, ignition.date.multiple.pivot$FirePolyID)

# Join with ignition.date to get all cols
ignition.date.multiple.filtered <- ignition.date.multiple.pivot |>
  left_join(ignition.date) |>
  select(-Wildfire_Dates, -count) |>
  mutate(Wildfire_Date = as.Date(Wildfire_Date))

# Add Date_est column
ignition.date.multiple.filtered <- ignition.date.multiple.filtered |>
  mutate(Date_est = "actual (most common & most recent of multiple dates)")


## Combine single and multiple --------------------------------------------

# Combine
ignition.date.assigned <- ignition.date.single |>
  bind_rows(ignition.date.multiple.filtered)

# Check for NAs
apply(ignition.date.assigned, 2, anyNA)

# Check for matching lengths
nrow(ignition.date) == nrow(ignition.date.assigned)

# Check for missing IDs
setdiff(ignition.date$FirePolyID, ignition.date.assigned$FirePolyID)

# Remove fires already assigned
ignition.date.filtered <- ignition.date.assigned |>
  filter(!FirePolyID %in% out.date.assigned$FirePolyID) |>
  filter(!FirePolyID %in% contain.date.filtered$FirePolyID) |>
  filter(!FirePolyID %in% control.date.filtered$FirePolyID) |>
  filter(!FirePolyID %in% discovery.date.filtered$FirePolyID)


# Based on Other date -----------------------------------------------------

# Fires with other date in col A
other.date.A <- wildfires |>
  filter(DateA1 == "Listed Other Fire Date(s)")

#   Remove non-other date cols and rename A
other.date.A <- other.date.A |>
  select(FirePolyID, USGS_Assigned_ID, Fire_Year, DateA2) |>
  rename(Wildfire_Dates = DateA2)


# Fires with other date in col B
other.date.B <- wildfires |>
  filter(DateB1 == "Listed Other Fire Date(s)")

#   Remove non-other date cols and rename B
other.date.B <- other.date.B |>
  select(FirePolyID, USGS_Assigned_ID, Fire_Year, DateB2) |>
  rename(Wildfire_Dates = DateB2)


# Fires with other date in col C
other.date.C <- wildfires |>
  filter(DateC1 == "Listed Other Fire Date(s)")

#   Remove non-other date cols and rename C
other.date.C <- other.date.C |>
  select(FirePolyID, USGS_Assigned_ID, Fire_Year, DateC2) |>
  rename(Wildfire_Dates = DateC2)


# Fires with other date in col D
other.date.D <- wildfires |>
  filter(DateD1 == "Listed Other Fire Date(s)")

#   Remove non-other date cols and rename D
other.date.D <- other.date.D |>
  select(FirePolyID, USGS_Assigned_ID, Fire_Year, DateD2) |>
  rename(Wildfire_Dates = DateD2)


# Fires with other date in col E
other.date.E <- wildfires |>
  filter(DateE1 == "Listed Other Fire Date(s)")

#   Remove non-other date cols and rename E
other.date.E <- other.date.E |>
  select(FirePolyID, USGS_Assigned_ID, Fire_Year, DateE2) |>
  rename(Wildfire_Dates = DateE2)


# Fires with other date in col F
other.date.F <- wildfires |>
  filter(DateF1 == "Listed Other Fire Date(s)")

#   Remove non-other date cols and rename F
other.date.F <- other.date.F |>
  select(FirePolyID, USGS_Assigned_ID, Fire_Year, DateF2) |>
  rename(Wildfire_Dates = DateF2)


# Fires with other date in col G
other.date.G <- wildfires |>
  filter(DateG1 == "Listed Other Fire Date(s)")

#   Remove non-other date cols and rename G
other.date.G <- other.date.G |>
  select(FirePolyID, USGS_Assigned_ID, Fire_Year, DateG2) |>
  rename(Wildfire_Dates = DateG2)


# Fires with other date in col H
other.date.H <- wildfires |>
  filter(DateH1 == "Listed Other Fire Date(s)")

#   Remove non-other date cols and rename H
other.date.H <- other.date.H |>
  select(FirePolyID, USGS_Assigned_ID, Fire_Year, DateH2) |>
  rename(Wildfire_Dates = DateH2)


# Combine all
other.date <- bind_rows(
  other.date.A, other.date.B, other.date.C, other.date.D,
  other.date.E, other.date.F, other.date.G, other.date.H
)

# Remove fires already assigned
other.date <- other.date |>
  filter(!FirePolyID %in% out.date.assigned$FirePolyID) |>
  filter(!FirePolyID %in% contain.date.filtered$FirePolyID) |>
  filter(!FirePolyID %in% control.date.filtered$FirePolyID) |>
  filter(!FirePolyID %in% discovery.date.filtered$FirePolyID) |>
  filter(!FirePolyID %in% ignition.date.filtered$FirePolyID)

# Add Date_type col
other.date <- other.date |>
  mutate(Date_type = "Other Fire")

# These dates don't really seem remotely close to the Fire_Year, so I'll not use any
#   of them for assignment.


# Combine all with assignment so far --------------------------------------

# Assigned a Wildfire_Date
fire.assigned <- bind_rows(
  out.date.assigned, contain.date.filtered, control.date.filtered,
  discovery.date.filtered, ignition.date.filtered
)

# Check for duplicate IDs
length(unique(fire.assigned$FirePolyID)) == nrow(fire.assigned)


# Calculate fire duration -------------------------------------------------

# From Discovery Date to Control Date
discovery.control <- discovery.date.assigned |>
  filter(FirePolyID %in% control.date.assigned$FirePolyID) |>
  rename(Discovery_Date = Wildfire_Date) |>
  select(-Date_type, -Date_est) |>
  left_join(control.date.assigned) |>
  rename(Control_Date = Wildfire_Date) |>
  select(-Date_type, -Date_est) |>
  mutate(Fire_duration = as.numeric(Control_Date - Discovery_Date))
summary(discovery.control$Fire_duration)


# From Discovery Date to Containment Date
discovery.contain <- discovery.date.assigned |>
  filter(FirePolyID %in% contain.date.assigned$FirePolyID) |>
  rename(Discovery_Date = Wildfire_Date) |>
  select(-Date_type, -Date_est) |>
  left_join(contain.date.assigned) |>
  rename(Contain_Date = Wildfire_Date) |>
  select(-Date_type, -Date_est) |>
  mutate(Fire_duration = as.numeric(Contain_Date - Discovery_Date))
summary(discovery.contain$Fire_duration)


# From Ignition Date to Control Date
ignition.control <- ignition.date.assigned |>
  filter(FirePolyID %in% control.date.assigned$FirePolyID) |>
  rename(Ignition_Date = Wildfire_Date) |>
  select(-Date_type, -Date_est) |>
  left_join(control.date.assigned) |>
  rename(Control_Date = Wildfire_Date) |>
  select(-Date_type, -Date_est) |>
  mutate(Fire_duration = as.numeric(Control_Date - Ignition_Date))
summary(ignition.control$Fire_duration)


# From Ignition Date to Containment Date
ignition.contain <- ignition.date.assigned |>
  filter(FirePolyID %in% contain.date.assigned$FirePolyID) |>
  rename(Ignition_Date = Wildfire_Date) |>
  select(-Date_type, -Date_est) |>
  left_join(contain.date.assigned) |>
  rename(Contain_Date = Wildfire_Date) |>
  select(-Date_type, -Date_est) |>
  mutate(Fire_duration = as.numeric(Contain_Date - Ignition_Date))
summary(ignition.contain$Fire_duration)


# Most fires last only about a week or less, so I won't adjust dates for those
#   assigned via ignition or discovery dates (for simplicity).


# Assign date based on Fire Year ------------------------------------------

# No date assigned yet
fire.no.date <- wildfires1 |>
  filter(!FirePolyID %in% fire.assigned$FirePolyID)

# Assign date as Dec 31
year.date.assigned <- fire.no.date |>
  select(FirePolyID, USGS_Assigned_ID, Fire_Year) |>
  mutate(
    Date_type = "Fire Year",
    Wildfire_Date = as.Date(paste0(Fire_Year, "-12-31")),
    Date_est = "estimated from Fire Year"
  )


# Combine all -------------------------------------------------------------

# Combine
wildfires.assigned <- fire.assigned |>
  bind_rows(year.date.assigned) |>
  arrange(FirePolyID)

# With all fire info columns
wildfires.assigned.allcol <- wildfires.assigned |>
  left_join(wildfires.all)

# Check for duplicate IDs
length(unique(wildfires.assigned$FirePolyID)) == nrow(wildfires.assigned)

# Check for missing IDs
setdiff(wildfires$FirePolyID, wildfires.assigned$FirePolyID)


# Write to CSV ------------------------------------------------------------

write_csv(wildfires.assigned,
  file = "data/versions-from-R/07.1_wildfires-with-assigned-dates_v004.csv"
)

write_csv(wildfires.assigned.allcol,
  file = "data/versions-from-R/07.1_wildfires-with-assigned-dates_all-fire-info_v004.csv"
)


save.image("RData/07.1_assign-wildfire-date_v004.RData")

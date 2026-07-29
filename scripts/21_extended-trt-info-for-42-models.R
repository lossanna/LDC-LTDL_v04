# Created: 2026-07-16
# Updated: 2026-07-28

# Purpose: Examine treatment information cols for each model.

library(tidyverse)

# Load data ---------------------------------------------------------------

trt.info.raw <- read_csv("data/raw/Treatment_Info_R.csv")
load("RData/15.1_matched-data.RData")
trtpolyid <- read_csv("data/versions-from-R/11.1_TrtPolyID-for-treatment-info-v008.csv")


# NW Forested Mts / Western Cordillera ------------------------------------

## Blue Mountains ---------------------------------------------------------

### 1. Herbicide ----------------------------------------------------------

# Identify relevant Trt_IDs
model01.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model01.matched$PrimaryKey)

# Attach treatment info
model01.trt <- model01.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw)

# Remove duplicate rows and add model number col
model01.trt <- model01.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 1", .before = Trt_ID) |> 
  mutate(Ecoregion = "Blue Mountains", .after = Model) |> 
  mutate(Treatment = "Herbicide", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model01.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model01-treatment-info.csv",
          na = "")

# Dates
summary(model01.matched$MR_trt_comp)
summary(model01.matched$DateVisited)


### 2. Post-burn herbicide ------------------------------------------------

# Identify relevant Trt_IDs
model02.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model02.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model02.trt <- model02.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Herbicide")

# Remove duplicate rows and add model number col
model02.trt <- model02.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 2", .before = Trt_ID) |> 
  mutate(Ecoregion = "Blue Mountains", .after = Model) |> 
  mutate(Treatment = "Post-burn herbicide", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model02.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model02-treatment-info.csv",
          na = "")

# Dates
summary(model02.matched$MR_trt_comp)
summary(model02.matched$DateVisited)



## Middle Rockies ---------------------------------------------------------

### 3. Herbicide ----------------------------------------------------------

# Identify relevant Trt_IDs
model03.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model03.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model03.trt <- model03.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Herbicide")

# Remove duplicate rows and add model number col
model03.trt <- model03.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 3", .before = Trt_ID) |> 
  mutate(Ecoregion = "Middle Rockies", .after = Model) |> 
  mutate(Treatment = "Herbicide", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model03.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model03-treatment-info.csv",
          na = "")

# Dates
summary(model03.matched$MR_trt_comp)
summary(model03.matched$DateVisited)



## Southern Rockies -------------------------------------------------------

### 4. Herbicide ----------------------------------------------------------

# Identify relevant Trt_IDs
model04.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model04.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model04.trt <- model04.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Herbicide")

# Remove duplicate rows and add model number col
model04.trt <- model04.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 4", .before = Trt_ID) |> 
  mutate(Ecoregion = "Southern Rockies", .after = Model) |> 
  mutate(Treatment = "Herbicide", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model04.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model04-treatment-info.csv",
          na = "")

# Dates
summary(model04.matched$MR_trt_comp)
summary(model04.matched$DateVisited)


### 5. Prescribed burn ----------------------------------------------------

# Identify relevant Trt_IDs
model05.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model05.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model05.trt <- model05.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Prescribed Burn")

# Remove duplicate rows and add model number col
model05.trt <- model05.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 5", .before = Trt_ID) |> 
  mutate(Ecoregion = "Southern Rockies", .after = Model) |> 
  mutate(Treatment = "Prescribed burn", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model05.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model05-treatment-info.csv",
          na = "")

# Dates
summary(model05.matched$MR_trt_comp)
summary(model05.matched$DateVisited)




# Great Plains / West-Central Semiarid Prairies ---------------------------

## Northwestern Great Plains ----------------------------------------------

### 6. Prescribed burn ----------------------------------------------------

# Identify relevant Trt_IDs
model06.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model06.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model06.trt <- model06.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Prescribed Burn")

# Remove duplicate rows and add model number col
model06.trt <- model06.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 6", .before = Trt_ID) |> 
  mutate(Ecoregion = "Northwestern Great Plains", .after = Model) |> 
  mutate(Treatment = "Prescribed burn", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model06.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model06-treatment-info.csv",
          na = "")

# Dates
summary(model06.matched$MR_trt_comp)
summary(model06.matched$DateVisited)




# Cold Deserts ------------------------------------------------------------

## Snake River Plain ------------------------------------------------------

### 7. Post-burn aerial seeding -------------------------------------------

# Identify relevant Trt_IDs
model07.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model07.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model07.trt <- model07.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Aerial Seeding")

# Remove duplicate rows and add model number col
model07.trt <- model07.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 7", .before = Trt_ID) |> 
  mutate(Ecoregion = "Snake River Plain", .after = Model) |> 
  mutate(Treatment = "Post-burn aerial seeding", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model07.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model07-treatment-info.csv",
          na = "")

# Dates
summary(model07.matched$MR_trt_comp)
summary(model07.matched$DateVisited)


### 8. Post-burn aerial & drill seeding -----------------------------------

# Identify relevant Trt_IDs
model08.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model08.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model08.trt <- model08.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub %in% c("Aerial Seeding", "Drill Seeding"))

# Remove duplicate rows and add model number col
model08.trt <- model08.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 8", .before = Trt_ID) |> 
  mutate(Ecoregion = "Snake River Plain", .after = Model) |> 
  mutate(Treatment = "Post-burn aerial & drill seeding", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model08.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model08-treatment-info.csv",
          na = "")

# Dates
summary(model08.matched$MR_trt_comp)
summary(model08.matched$DateVisited)


### 9. Post-burn closure --------------------------------------------------

# Identify relevant Trt_IDs
model09.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model09.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model09.trt <- model09.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Closure")

# Remove duplicate rows and add model number col
model09.trt <- model09.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 9", .before = Trt_ID) |> 
  mutate(Ecoregion = "Snake River Plain", .after = Model) |> 
  mutate(Treatment = "Post-burn closure", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model09.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model09-treatment-info.csv",
          na = "")

# Dates
summary(model09.matched$MR_trt_comp)
summary(model09.matched$DateVisited)


### 10. Post-burn drill seeding -------------------------------------------------

# Identify relevant Trt_IDs
model10.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model10.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model10.trt <- model10.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Drill Seeding")

# Remove duplicate rows and add model number col
model10.trt <- model10.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 10", .before = Trt_ID) |> 
  mutate(Ecoregion = "Snake River Plain", .after = Model) |> 
  mutate(Treatment = "Post-burn drill seeding", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model10.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model10-treatment-info.csv",
          na = "")

# Dates
summary(model10.matched$MR_trt_comp)
summary(model10.matched$DateVisited)


### 11. Post-burn herbicide -----------------------------------------------

# Identify relevant Trt_IDs
model11.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model11.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model11.trt <- model11.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Herbicide")

# Remove duplicate rows and add model number col
model11.trt <- model11.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 11", .before = Trt_ID) |> 
  mutate(Ecoregion = "Snake River Plain", .after = Model) |> 
  mutate(Treatment = "Post-burn herbicide", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model11.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model11-treatment-info.csv",
          na = "")

# Dates
summary(model11.matched$MR_trt_comp)
summary(model11.matched$DateVisited)



## Northern Basin and Range -----------------------------------------------

### 12. Drill seeding -----------------------------------------------------

# Identify relevant Trt_IDs
model12.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model12.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model12.trt <- model12.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Drill Seeding")

# Remove duplicate rows and add model number col
model12.trt <- model12.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 12", .before = Trt_ID) |> 
  mutate(Ecoregion = "Nothern Basin and Range", .after = Model) |> 
  mutate(Treatment = "Drill seeding", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model12.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model12-treatment-info.csv",
          na = "")

# Dates
summary(model12.matched$MR_trt_comp)
summary(model12.matched$DateVisited)


### 13. Drill seeding & soil disturbance ----------------------------------

# Identify relevant Trt_IDs
model13.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model13.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model13.trt <- model13.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub %in% c("Drill Seeding", "Soil Disturbance"))

# Remove duplicate rows and add model number col
model13.trt <- model13.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 13", .before = Trt_ID) |> 
  mutate(Ecoregion = "Nothern Basin and Range", .after = Model) |> 
  mutate(Treatment = "Drill seeding & soil disturbance", .after = Ecoregion) |> 
  arrange(Trt_ID) |> 
  arrange(Trt_Type_Sub)

# Write to CSV
write_csv(model13.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model13-treatment-info.csv",
          na = "")

# Dates
summary(model13.matched$MR_trt_comp)
summary(model13.matched$DateVisited)


### 14. Herbicide ---------------------------------------------------------

# Identify relevant Trt_IDs
model14.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model14.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model14.trt <- model14.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Herbicide")

# Remove duplicate rows and add model number col
model14.trt <- model14.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 14", .before = Trt_ID) |> 
  mutate(Ecoregion = "Nothern Basin and Range", .after = Model) |> 
  mutate(Treatment = "Herbicide", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model14.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model14-treatment-info.csv",
          na = "")

# Dates
summary(model14.matched$MR_trt_comp)
summary(model14.matched$DateVisited)






summary(model13.matched$MR_trt_comp)

save.image("RData/21_extended-trt-info-for-42-models.RData")

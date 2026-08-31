# Created: 2026-07-16
# Updated: 2026-08-06

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


### 15. Prescribed burn ---------------------------------------------------

# Identify relevant Trt_IDs
model15.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model15.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model15.trt <- model15.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Prescribed Burn")

# Remove duplicate rows and add model number col
model15.trt <- model15.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 15", .before = Trt_ID) |> 
  mutate(Ecoregion = "Nothern Basin and Range", .after = Model) |> 
  mutate(Treatment = "Prescribed burn", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model15.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model15-treatment-info.csv",
          na = "")

# Dates
summary(model15.matched$MR_trt_comp)
summary(model15.matched$DateVisited)


### 16. Vegetation disturbance --------------------------------------------

# Identify relevant Trt_IDs
model16.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model16.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model16.trt <- model16.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Vegetation Disturbance")

# Remove duplicate rows and add model number col
model16.trt <- model16.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 16", .before = Trt_ID) |> 
  mutate(Ecoregion = "Nothern Basin and Range", .after = Model) |> 
  mutate(Treatment = "Vegetation disturbance", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model16.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model16-treatment-info.csv",
          na = "")

# Dates
summary(model16.matched$MR_trt_comp)
summary(model16.matched$DateVisited)


### 17. Post-burn aerial seeding ------------------------------------------

# Identify relevant Trt_IDs
model17.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model17.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model17.trt <- model17.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Aerial Seeding")

# Remove duplicate rows and add model number col
model17.trt <- model17.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 17", .before = Trt_ID) |> 
  mutate(Ecoregion = "Nothern Basin and Range", .after = Model) |> 
  mutate(Treatment = "Post-burn aerial seeding", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model17.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model17-treatment-info.csv",
          na = "")

# Dates
summary(model17.matched$MR_trt_comp)
summary(model17.matched$DateVisited)


### 18. Post-burn aerial & drill seeding ----------------------------------

# Identify relevant Trt_IDs
model18.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model18.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model18.trt <- model18.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub %in% c("Aerial Seeding", "Drill Seeding"))

# Remove duplicate rows and add model number col
model18.trt <- model18.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 18", .before = Trt_ID) |> 
  mutate(Ecoregion = "Nothern Basin and Range", .after = Model) |> 
  mutate(Treatment = "Post-burn aerial & drill seeding", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model18.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model18-treatment-info.csv",
          na = "")

# Dates
summary(model18.matched$MR_trt_comp)
summary(model18.matched$DateVisited)


### 19. Post-burn closure -------------------------------------------------

# Identify relevant Trt_IDs
model19.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model19.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model19.trt <- model19.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Closure")

# Remove duplicate rows and add model number col
model19.trt <- model19.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 19", .before = Trt_ID) |> 
  mutate(Ecoregion = "Nothern Basin and Range", .after = Model) |> 
  mutate(Treatment = "Post-burn closure", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model19.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model19-treatment-info.csv",
          na = "")

# Dates
summary(model19.matched$MR_trt_comp)
summary(model19.matched$DateVisited)


### 20. Post-burn drill seeding -------------------------------------------------

# Identify relevant Trt_IDs
model20.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model20.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model20.trt <- model20.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Drill Seeding")

# Remove duplicate rows and add model number col
model20.trt <- model20.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 20", .before = Trt_ID) |> 
  mutate(Ecoregion = "Nothern Basin and Range", .after = Model) |> 
  mutate(Treatment = "Post-burn drill seeding", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model20.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model20-treatment-info.csv",
          na = "")

# Dates
summary(model20.matched$MR_trt_comp)
summary(model20.matched$DateVisited)


### 21. Post-burn herbicide -----------------------------------------------

# Identify relevant Trt_IDs
model21.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model21.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model21.trt <- model21.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Herbicide")

# Remove duplicate rows and add model number col
model21.trt <- model21.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 21", .before = Trt_ID) |> 
  mutate(Ecoregion = "Nothern Basin and Range", .after = Model) |> 
  mutate(Treatment = "Post-burn herbicide", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model21.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model21-treatment-info.csv",
          na = "")

# Dates
summary(model21.matched$MR_trt_comp)
summary(model21.matched$DateVisited)


### 22. Post-burn seedling planting ---------------------------------------

# Identify relevant Trt_IDs
model22.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model22.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model22.trt <- model22.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Seedling Planting")

# Remove duplicate rows and add model number col
model22.trt <- model22.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 22", .before = Trt_ID) |> 
  mutate(Ecoregion = "Nothern Basin and Range", .after = Model) |> 
  mutate(Treatment = "Post-burn seedling planting", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model22.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model22-treatment-info.csv",
          na = "")

# Dates
summary(model22.matched$MR_trt_comp)
summary(model22.matched$DateVisited)



## Central Basin and Range ------------------------------------------------

### 23. Drill seeding & soil disturbance ----------------------------------

# Identify relevant Trt_IDs
model23.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model23.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model23.trt <- model23.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub %in% c("Drill Seeding", "Soil Disturbance"))

# Remove duplicate rows and add model number col
model23.trt <- model23.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 23", .before = Trt_ID) |> 
  mutate(Ecoregion = "Central Basin and Range", .after = Model) |> 
  mutate(Treatment = "Drill seeding & soil disturbance", .after = Ecoregion) |> 
  arrange(Trt_ID) |> 
  arrange(Trt_Type_Sub)

# Write to CSV
write_csv(model23.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model23-treatment-info.csv",
          na = "")

# Dates
summary(model23.matched$MR_trt_comp)
summary(model23.matched$DateVisited)


### 24. Prescribed burn ---------------------------------------------------

# Identify relevant Trt_IDs
model24.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model24.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model24.trt <- model24.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Prescribed Burn")

# Remove duplicate rows and add model number col
model24.trt <- model24.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 24", .before = Trt_ID) |> 
  mutate(Ecoregion = "Central Basin and Range", .after = Model) |> 
  mutate(Treatment = "Prescribed burn", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model24.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model24-treatment-info.csv",
          na = "")

# Dates
summary(model24.matched$MR_trt_comp)
summary(model24.matched$DateVisited)


### 25. Vegetation disturbance --------------------------------------------

# Identify relevant Trt_IDs
model25.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model25.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model25.trt <- model25.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Vegetation Disturbance")

# Remove duplicate rows and add model number col
model25.trt <- model25.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 25", .before = Trt_ID) |> 
  mutate(Ecoregion = "Central Basin and Range", .after = Model) |> 
  mutate(Treatment = "Vegetation disturbance", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model25.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model25-treatment-info.csv",
          na = "")

# Dates
summary(model25.matched$MR_trt_comp)
summary(model25.matched$DateVisited)


### 26. Post-burn aerial seeding ------------------------------------------

# Identify relevant Trt_IDs
model26.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model26.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model26.trt <- model26.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Aerial Seeding")

# Remove duplicate rows and add model number col
model26.trt <- model26.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 26", .before = Trt_ID) |> 
  mutate(Ecoregion = "Central Basin and Range", .after = Model) |> 
  mutate(Treatment = "Post-burn aerial seeding", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model26.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model26-treatment-info.csv",
          na = "")

# Dates
summary(model26.matched$MR_trt_comp)
summary(model26.matched$DateVisited)



### 27. Post-burn drill seeding -------------------------------------------

# Identify relevant Trt_IDs
model27.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model27.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model27.trt <- model27.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Drill Seeding")

# Remove duplicate rows and add model number col
model27.trt <- model27.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 27", .before = Trt_ID) |> 
  mutate(Ecoregion = "Central Basin and Range", .after = Model) |> 
  mutate(Treatment = "Post-burn drill seeding", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model27.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model27-treatment-info.csv",
          na = "")

# Dates
summary(model27.matched$MR_trt_comp)
summary(model27.matched$DateVisited)


### 28. Post-burn ground seeding ------------------------------------------

# Identify relevant Trt_IDs
model28.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model28.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model28.trt <- model28.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Ground Seeding")

# Remove duplicate rows and add model number col
model28.trt <- model28.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 28", .before = Trt_ID) |> 
  mutate(Ecoregion = "Central Basin and Range", .after = Model) |> 
  mutate(Treatment = "Post-burn ground seeding", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model28.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model28-treatment-info.csv",
          na = "")

# Dates
summary(model28.matched$MR_trt_comp)
summary(model28.matched$DateVisited)


### 29. Post-burn herbicide -----------------------------------------------

# Identify relevant Trt_IDs
model29.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model29.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model29.trt <- model29.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Herbicide")

# Remove duplicate rows and add model number col
model29.trt <- model29.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 29", .before = Trt_ID) |> 
  mutate(Ecoregion = "Central Basin and Range", .after = Model) |> 
  mutate(Treatment = "Post-burn herbicide", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model29.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model29-treatment-info.csv",
          na = "")

# Dates
summary(model29.matched$MR_trt_comp)
summary(model29.matched$DateVisited)



## Wyoming Basin ----------------------------------------------------------

### 30. Prescribed burn ---------------------------------------------------

# Identify relevant Trt_IDs
model30.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model30.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model30.trt <- model30.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Prescribed Burn")

# Remove duplicate rows and add model number col
model30.trt <- model30.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 30", .before = Trt_ID) |> 
  mutate(Ecoregion = "Wyoming Basin", .after = Model) |> 
  mutate(Treatment = "Prescribed burn", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model30.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model30-treatment-info.csv",
          na = "")

# Dates
summary(model30.matched$MR_trt_comp)
summary(model30.matched$DateVisited)



## Colorado Plateaus ------------------------------------------------------

### 31. Aerial seeding & soil disturbance ---------------------------------

# Identify relevant Trt_IDs
model31.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model31.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model31.trt <- model31.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub %in% c("Aerial Seeding", "Soil Disturbance"))

# Remove duplicate rows and add model number col
model31.trt <- model31.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 31", .before = Trt_ID) |> 
  mutate(Ecoregion = "Colorado Plateaus", .after = Model) |> 
  mutate(Treatment = "Aerial seeding & soil disturbance", .after = Ecoregion) |> 
  arrange(Trt_ID) |> 
  arrange(Trt_Type_Sub)

# Write to CSV
write_csv(model31.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model31-treatment-info.csv",
          na = "")

# Dates
summary(model31.matched$MR_trt_comp)
summary(model31.matched$DateVisited)


### 32. Herbicide ---------------------------------------------------------

# Identify relevant Trt_IDs
model32.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model32.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model32.trt <- model32.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Herbicide")

# Remove duplicate rows and add model number col
model32.trt <- model32.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 32", .before = Trt_ID) |> 
  mutate(Ecoregion = "Colorado Plateaus", .after = Model) |> 
  mutate(Treatment = "Herbicide", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model32.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model32-treatment-info.csv",
          na = "")

# Dates
summary(model32.matched$MR_trt_comp)
summary(model32.matched$DateVisited)


### 33. Prescribed burn ---------------------------------------------------

# Identify relevant Trt_IDs
model33.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model33.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model33.trt <- model33.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Prescribed Burn")

# Remove duplicate rows and add model number col
model33.trt <- model33.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 33", .before = Trt_ID) |> 
  mutate(Ecoregion = "Colorado Plateaus", .after = Model) |> 
  mutate(Treatment = "Prescribed Burn", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model33.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model33-treatment-info.csv",
          na = "")

# Dates
summary(model33.matched$MR_trt_comp)
summary(model33.matched$DateVisited)


### 34. Soil disturbance --------------------------------------------------

# Identify relevant Trt_IDs
model34.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model34.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model34.trt <- model34.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Soil Disturbance")

# Remove duplicate rows and add model number col
model34.trt <- model34.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 34", .before = Trt_ID) |> 
  mutate(Ecoregion = "Colorado Plateaus", .after = Model) |> 
  mutate(Treatment = "Soil disturbance", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model34.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model34-treatment-info.csv",
          na = "")

# Dates
summary(model34.matched$MR_trt_comp)
summary(model34.matched$DateVisited)


### 35. Vegetation disturbance --------------------------------------------

# Identify relevant Trt_IDs
model35.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model35.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model35.trt <- model35.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Vegetation Disturbance")

# Remove duplicate rows and add model number col
model35.trt <- model35.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 35", .before = Trt_ID) |> 
  mutate(Ecoregion = "Colorado Plateaus", .after = Model) |> 
  mutate(Treatment = "Vegetation disturbance", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model35.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model35-treatment-info.csv",
          na = "")

# Dates
summary(model35.matched$MR_trt_comp)
summary(model35.matched$DateVisited)


### 36. Post-burn aerial seeding ------------------------------------------

# Identify relevant Trt_IDs
model36.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model36.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model36.trt <- model36.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Aerial Seeding")

# Remove duplicate rows and add model number col
model36.trt <- model36.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 36", .before = Trt_ID) |> 
  mutate(Ecoregion = "Colorado Plateaus", .after = Model) |> 
  mutate(Treatment = "Post-burn aerial seeding", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model36.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model36-treatment-info.csv",
          na = "")

# Dates
summary(model36.matched$MR_trt_comp)
summary(model36.matched$DateVisited)



## Arizona/New Mexico Plateau ---------------------------------------------

### 37. Herbicide ---------------------------------------------------------

# Identify relevant Trt_IDs
model37.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model37.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model37.trt <- model37.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Herbicide")

# Remove duplicate rows and add model number col
model37.trt <- model37.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 37", .before = Trt_ID) |> 
  mutate(Ecoregion = "Arizona/New Mexico Plateau", .after = Model) |> 
  mutate(Treatment = "Herbicide", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model37.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model37-treatment-info.csv",
          na = "")

# Dates
summary(model37.matched$MR_trt_comp)
summary(model37.matched$DateVisited)


### 38. Prescribed burn ---------------------------------------------------

# Identify relevant Trt_IDs
model38.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model38.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model38.trt <- model38.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Prescribed Burn")

# Remove duplicate rows and add model number col
model38.trt <- model38.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 38", .before = Trt_ID) |> 
  mutate(Ecoregion = "Arizona/New Mexico Plateau", .after = Model) |> 
  mutate(Treatment = "Prescribed burn", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model38.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model38-treatment-info.csv",
          na = "")

# Dates
summary(model38.matched$MR_trt_comp)
summary(model38.matched$DateVisited)


### 39. Soil disturbance --------------------------------------------------

# Identify relevant Trt_IDs
model39.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model39.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model39.trt <- model39.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Soil Disturbance")

# Remove duplicate rows and add model number col
model39.trt <- model39.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 39", .before = Trt_ID) |> 
  mutate(Ecoregion = "Arizona/New Mexico Plateau", .after = Model) |> 
  mutate(Treatment = "Soil disturbance", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model39.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model39-treatment-info.csv",
          na = "")

# Dates
summary(model39.matched$MR_trt_comp)
summary(model39.matched$DateVisited)




# Warm Deserts ------------------------------------------------------------

## Mojave Basin and Range -------------------------------------------------

### 40. Post-burn aerial seeding ------------------------------------------

# Identify relevant Trt_IDs
model40.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model40.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model40.trt <- model40.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Aerial Seeding")

# Remove duplicate rows and add model number col
model40.trt <- model40.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 40", .before = Trt_ID) |> 
  mutate(Ecoregion = "Mojave Basin and Range", .after = Model) |> 
  mutate(Treatment = "Post-burn aerial seeding", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model40.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model40-treatment-info.csv",
          na = "")

# Dates
summary(model40.matched$MR_trt_comp)
summary(model40.matched$DateVisited)



## Chihuahuan Desert ------------------------------------------------------

### 41. Herbicide ---------------------------------------------------------

# Identify relevant Trt_IDs
model41.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model41.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model41.trt <- model41.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Herbicide")

# Remove duplicate rows and add model number col
model41.trt <- model41.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 41", .before = Trt_ID) |> 
  mutate(Ecoregion = "Chihuahuan Desert", .after = Model) |> 
  mutate(Treatment = "Herbicide", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model41.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model41-treatment-info.csv",
          na = "")

# Dates
summary(model41.matched$MR_trt_comp)
summary(model41.matched$DateVisited)




# Temperate Sierras / Upper Gila ------------------------------------------

## Arizona/New Mexico Mountains -------------------------------------------

### 42. Prescribed Burn ---------------------------------------------------

# Identify relevant Trt_IDs
model42.trtid <- trtpolyid |> 
  filter(PrimaryKey %in% model42.matched$PrimaryKey)

# Attach treatment info and filter for specific Trt_Type_Sub
model42.trt <- model42.trtid |> 
  select(PrimaryKey, Trt_ID) |> 
  left_join(trt.info.raw) |> 
  filter(Trt_Type_Sub == "Prescribed Burn")

# Remove duplicate rows and add model number col
model42.trt <- model42.trt |> 
  select(-PrimaryKey) |> 
  distinct(.keep_all = TRUE) |> 
  mutate(Model = "Model 42", .before = Trt_ID) |> 
  mutate(Ecoregion = "Arizona/New Mexico Mountains", .after = Model) |> 
  mutate(Treatment = "Prescribed burn", .after = Ecoregion) |> 
  arrange(Trt_ID)

# Write to CSV
write_csv(model42.trt,
          file = "data/data-wrangling-intermediate/21_treatment-info/model42-treatment-info.csv",
          na = "")

# Dates
summary(model42.matched$MR_trt_comp)
summary(model42.matched$DateVisited)



save.image("RData/21_extended-trt-info-for-42-models.RData")

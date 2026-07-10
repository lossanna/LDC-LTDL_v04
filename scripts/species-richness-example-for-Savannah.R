# Created: 2026-07-10
# Updated: 2026-07-10

# Purpose: Gather LDC points from Snake River Plain closure treatment and equivalent controls
#   to then calculate richness and demonstrate significant difference between
#   treated and control.

library(tidyverse)

# Load data ---------------------------------------------------------------

all.matched.cover <- read_csv("data/versions-from-R/17_all-matched-data-with-species-cover_v012.csv")
all.matched.diversity <- read_csv("data/versions-from-R/18_shannon-diversity_all-models_v012.csv")
geoindicators <- read_csv("data/versions-from-R/14.3_geoindicators.csv")



# Data wrangling ----------------------------------------------------------

# Ecological Site ID col
ecositeid <- geoindicators |> 
  select(PrimaryKey, EcoSiteID)


# SRP closure (richness) --------------------------------------------------

# Filter to just SRP closure points
srp.closure <- all.matched.cover |> 
  filter(Model == 9) |> 
  select(PrimaryKey, trt_control, MLRA_name, CurrentPLANTSCode, Cover_AH)

srp.closure.points <- srp.closure |> 
  select(PrimaryKey, trt_control, MLRA_name) |> 
  distinct(.keep_all = TRUE)

# Calculate richness for each point
srp.closure.richness <- srp.closure |> 
  group_by(PrimaryKey) |> 
  summarise(richness = n())

# Join to assign Ecological Site ID and trt_control
srp.closure.richness <- srp.closure.richness |> 
  left_join(ecositeid) |> 
  left_join(srp.closure.points)


# Count number of points per Ecological Site ID
count(srp.closure.richness, EcoSiteID) |> 
  arrange(desc(n))

srp.closure.richness |> 
  group_by(EcoSiteID, trt_control) |> 
  summarise(n = n()) |> 
  arrange(desc(n)) |> 
  print(n = 20)

# Use Ecological Site R011BY010ID 
dat <- srp.closure.richness |> 
  filter(str_detect(EcoSiteID, "011AY014ID"))

count(dat, trt_control)

# t-test
t.test(dat$richness ~ dat$trt_control) # NS

t.test(srp.closure.richness$richness ~ srp.closure.richness$trt_control)


# Count number of points per MLRA
srp.closure.richness |> 
  group_by(MLRA_name, trt_control) |> 
  summarise(n = n()) |> 
  arrange(desc(n))

dat <- srp.closure.richness |> 
  filter(MLRA_name == "Snake River Plains")

# t-test
t.test(dat$richness ~ dat$trt_control)



# SRP closure (Shannon) ---------------------------------------------------

srp.closure.shannon <- all.matched.diversity |> 
  filter(Model == 9) |> 
  select(PrimaryKey, trt_control, Shannon)


# SRP aerial seeding (richness) -------------------------------------------

# Filter to just SRP aerial seeding points
srp.aerial <- all.matched.cover |> 
  filter(Model == 7) |> 
  select(PrimaryKey, trt_control, CurrentPLANTSCode, Cover_AH)

srp.aerial.points <- srp.aerial |> 
  select(PrimaryKey, trt_control) |> 
  distinct(.keep_all = TRUE)

# Ecological Site ID col
ecositeid <- geoindicators |> 
  select(PrimaryKey, EcoSiteID)


# Calculate richness for each point
srp.aerial.richness <- srp.aerial |> 
  group_by(PrimaryKey) |> 
  summarise(richness = n())

# Join to assign Ecological Site ID and trt_control
srp.aerial.richness <- srp.aerial.richness |> 
  left_join(ecositeid) |> 
  left_join(srp.aerial.points)


# Count number of points per Ecological Site ID
count(srp.aerial.richness, EcoSiteID) |> 
  arrange(desc(n))

srp.aerial.richness |> 
  group_by(EcoSiteID, trt_control) |> 
  summarise(n = n()) |> 
  arrange(desc(n))

# Use Ecological Site R011BY010ID 
dat <- srp.aerial.richness |> 
  filter(EcoSiteID == "R011XY001ID")

count(dat, trt_control)

# t-test
t.test(dat$richness ~ dat$trt_control)

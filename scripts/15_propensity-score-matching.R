# Created: 2026-06-29
# Updated: 2026-06-29

# Purpose: Propensity score matching incorporating CETWI and SOLUS data with updated
#   numbering from PSM 1 (results are the same as 13_propensity-score-matching-1.R).

# https://kosukeimai.github.io/MatchIt/articles/MatchIt.html

library(tidyverse)
library(MatchIt)
library(cobalt)
library(marginaleffects)
library(ggsignif)

# Load data ---------------------------------------------------------------

ldc.011.raw <- read_csv("data/versions-from-R/14.4_LDC-points_v011.csv")
level3.trt30 <- read_csv("data/versions-from-R/14.4_eco3-trt30_count-table.csv")


# Data wrangling ----------------------------------------------------------

# Create binary col for treatment and fire
ldc.011 <- ldc.011.raw |>
  mutate(trt_binary = if_else(is.na(Trt_Type_Sub), 0, 1))

# Natural log transformation of horizontal flux (q) when greater than 0
ldc.011 <- ldc.011 |>
  mutate(
    ln_q = if_else(horizontal_flux_total_MD == 0, 0, log(horizontal_flux_total_MD)),
    .after = horizontal_flux_total_MD
  )


# Check for NAs
apply(ldc.011, 2, anyNA)



# NW Forested Mts / Western Cordillera ------------------------------------

## Blue Mountains ---------------------------------------------------------

### 1. Herbicide ----------------------------------------------------------

# Filter data
model01.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Herbicide" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Blue Mountains")

# PSM
model01.psm <- matchit(
  data = model01.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model01.psm
summary(model01.psm) # 96 treated matched

# Diagnostic love plot
model01.loveplot <- love.plot(model01.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "1. Blue Mountains: Herbicide") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model01.loveplot

# eCDF plots
bal.plot(model01.psm, which = "both", type = "ecdf")
bal.plot(model01.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model01.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model01.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model01.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model01.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model01.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model01.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model01.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model01.psm, which = "both")
bal.plot(model01.psm, var = "BareSoil_FH", which = "both")
bal.plot(model01.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model01.psm, var = "ForbCover_AH", which = "both")
bal.plot(model01.psm, var = "GramCover_AH", which = "both")
bal.plot(model01.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model01.psm, var = "Gap100plus", which = "both")
bal.plot(model01.psm, var = "CETWI", which = "both")
bal.plot(model01.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model01.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model01.psm, type = "qq")

# Matched data
model01.matched <- match_data(model01.psm)

# Create trt_control variable
model01.matched <- model01.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Herbicide", "Control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Control", "Herbicide")))

# Center and scale numeric variables
model01.matched <- model01.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model01.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model01.matched,
  weights = weights
)

# G computation to estimate marginal effects
model01.pred <- avg_predictions(
  model = model01.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 1, .before = trt_control)
model01.pred

# Estimation of average treatment effect
model01.comp <- avg_comparisons(
  model = model01.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 1, .before = term)
model01.comp

# Plot
model01.plot <- model01.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "1. Blue Mountains: Herbicide",
       x = NULL)
model01.plot


### 2. Post-burn herbicide ------------------------------------------------

# Filter data
model02.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Herbicide" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post")) |>
  filter(EcoLvl3 == "Blue Mountains")

# PSM
model02.psm <- matchit(
  data = model02.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model02.psm
summary(model02.psm) # 42 treated matched

# Diagnostic love plot
model02.loveplot <- love.plot(model02.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "2. Blue Mountains: Post-burn herbicide") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model02.loveplot

# eCDF plots
bal.plot(model02.psm, which = "both", type = "ecdf")
bal.plot(model02.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model02.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model02.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model02.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model02.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model02.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model02.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model02.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(model02.psm, type = "density")
bal.plot(model02.psm, which = "both")
bal.plot(model02.psm, var = "BareSoil_FH", which = "both")
bal.plot(model02.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model02.psm, var = "ForbCover_AH", which = "both")
bal.plot(model02.psm, var = "GramCover_AH", which = "both")
bal.plot(model02.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model02.psm, var = "Gap100plus", which = "both")
bal.plot(model02.psm, var = "CETWI", which = "both")
bal.plot(model02.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model02.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model02.psm, type = "qq")

# Matched data
model02.matched <- match_data(model02.psm)

# Create trt_control variable
model02.matched <- model02.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Post-burn herbicide", "Post-burn control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn herbicide")))

# Center and scale numeric variables
model02.matched <- model02.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model02.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model02.matched,
  weights = weights
)

# G computation to estimate marginal effects
model02.pred <- avg_predictions(
  model = model02.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 2, .before = trt_control)
model02.pred

# Estimation of average treatment effect
model02.comp <- avg_comparisons(
  model = model02.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 2, .before = term)
model02.comp

# Plot
model02.plot <- model02.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "2. Blue Mountains: Post-burn herbicide",
       x = NULL)
model02.plot



## Middle Rockies ---------------------------------------------------------

### 3. Herbicide ----------------------------------------------------------

# Filter data
model03.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Herbicide" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Middle Rockies")

# PSM
model03.psm <- matchit(
  data = model03.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model03.psm
summary(model03.psm) # 33 treated matched

# Diagnostic love plot
model03.loveplot <- love.plot(model03.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "3. Middle Rockies: Herbicide") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model03.loveplot

# eCDF plots
bal.plot(model03.psm, which = "both", type = "ecdf")
bal.plot(model03.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model03.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model03.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model03.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model03.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model03.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model03.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model03.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model03.psm, which = "both")
bal.plot(model03.psm, var = "BareSoil_FH", which = "both")
bal.plot(model03.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model03.psm, var = "ForbCover_AH", which = "both")
bal.plot(model03.psm, var = "GramCover_AH", which = "both")
bal.plot(model03.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model03.psm, var = "Gap100plus", which = "both")
bal.plot(model03.psm, var = "CETWI", which = "both")
bal.plot(model03.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model03.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model03.psm, type = "qq")

# Matched data
model03.matched <- match_data(model03.psm)

# Create trt_control variable
model03.matched <- model03.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Herbicide", "Control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Control", "Herbicide")))

# Center and scale numeric variables
model03.matched <- model03.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model03.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model03.matched,
  weights = weights
)

# G computation to estimate marginal effects
model03.pred <- avg_predictions(
  model = model03.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 3, .before = trt_control)
model03.pred

# Estimation of average treatment effect
model03.comp <- avg_comparisons(
  model = model03.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 3, .before = term)
model03.comp

# Plot
model03.plot <- model03.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "3. Middle Rockies: Herbicide",
       x = NULL)
model03.plot



## Southern Rockies -------------------------------------------------------

### 4. Herbicide ----------------------------------------------------------

# Filter data
model04.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Herbicide" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Southern Rockies")

# PSM
model04.psm <- matchit(
  data = model04.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model04.psm
summary(model04.psm) # 76 treated matched

# Diagnostic love plot
model04.loveplot <- love.plot(model04.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "4. Southern Rockies: Herbicide") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model04.loveplot

# eCDF plots
bal.plot(model04.psm, which = "both", type = "ecdf")
bal.plot(model04.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model04.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model04.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model04.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model04.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model04.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model04.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model04.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model04.psm, which = "both")
bal.plot(model04.psm, var = "BareSoil_FH", which = "both")
bal.plot(model04.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model04.psm, var = "ForbCover_AH", which = "both")
bal.plot(model04.psm, var = "GramCover_AH", which = "both")
bal.plot(model04.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model04.psm, var = "Gap100plus", which = "both")
bal.plot(model04.psm, var = "CETWI", which = "both")
bal.plot(model04.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model04.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model04.psm, type = "qq")

# Matched data
model04.matched <- match_data(model04.psm)

# Create trt_control variable
model04.matched <- model04.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Herbicide", "Control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Control", "Herbicide")))

# Center and scale numeric variables
model04.matched <- model04.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model04.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model04.matched,
  weights = weights
)

# G computation to estimate marginal effects
model04.pred <- avg_predictions(
  model = model04.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 4, .before = trt_control)
model04.pred

# Estimation of average treatment effect
model04.comp <- avg_comparisons(
  model = model04.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 4, .before = term)
model04.comp

# Plot
model04.plot <- model04.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "4. Southern Rockies: Herbicide",
       x = NULL)
model04.plot


### 5. Prescribed burn ----------------------------------------------------

# Filter data
model05.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Prescribed Burn" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Southern Rockies")

# PSM
model05.psm <- matchit(
  data = model05.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model05.psm
summary(model05.psm) # 37 treated matched

# Diagnostic love plot
model05.loveplot <- love.plot(model05.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "5. Southern Rockies: Prescribed burn") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model05.loveplot

# eCDF plots
bal.plot(model05.psm, which = "both", type = "ecdf")
bal.plot(model05.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model05.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model05.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model05.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model05.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model05.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model05.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model05.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model05.psm, which = "both")
bal.plot(model05.psm, var = "BareSoil_FH", which = "both")
bal.plot(model05.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model05.psm, var = "ForbCover_AH", which = "both")
bal.plot(model05.psm, var = "GramCover_AH", which = "both")
bal.plot(model05.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model05.psm, var = "Gap100plus", which = "both")
bal.plot(model05.psm, var = "CETWI", which = "both")
bal.plot(model05.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model05.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model05.psm, type = "qq")

# Matched data
model05.matched <- match_data(model05.psm)

# Create trt_control variable
model05.matched <- model05.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Prescribed burn", "Control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Control", "Prescribed burn")))

# Center and scale numeric variables
model05.matched <- model05.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model05.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model05.matched,
  weights = weights
)

# G computation to estimate marginal effects
model05.pred <- avg_predictions(
  model = model05.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 5, .before = trt_control)
model05.pred

# Estimation of average treatment effect
model05.comp <- avg_comparisons(
  model = model05.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 5, .before = term)
model05.comp # p = 0.03

# Plot
model05.plot <- model05.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "5. Southern Rockies: Prescribed burn",
       x = NULL) + 
  geom_signif(
    comparisons = list(c("Control", "Prescribed burn")),
    annotations = c("*")
  )
model05.plot




# Great Plains / West-Central Semiarid Prairies ---------------------------

## Northwestern Great Plains ----------------------------------------------

### 6. Prescribed burn ----------------------------------------------------

# Filter data
model06.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Prescribed Burn" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Northwestern Great Plains")

# PSM
model06.psm <- matchit(
  data = model06.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model06.psm
summary(model06.psm) # 51 treated matched

# Diagnostic love plot
model06.loveplot <- love.plot(model06.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "6. Northwestern Great Plains: Prescribed burn") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model06.loveplot

# eCDF plots
bal.plot(model06.psm, which = "both", type = "ecdf")
bal.plot(model06.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model06.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model06.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model06.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model06.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model06.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model06.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model06.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model06.psm, which = "both")
bal.plot(model06.psm, var = "BareSoil_FH", which = "both")
bal.plot(model06.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model06.psm, var = "ForbCover_AH", which = "both")
bal.plot(model06.psm, var = "GramCover_AH", which = "both")
bal.plot(model06.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model06.psm, var = "Gap100plus", which = "both")
bal.plot(model06.psm, var = "CETWI", which = "both")
bal.plot(model06.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model06.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model06.psm, type = "qq")

# Matched data
model06.matched <- match_data(model06.psm)

# Create trt_control variable
model06.matched <- model06.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Prescribed burn", "Control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Control", "Prescribed burn")))

# Center and scale numeric variables
model06.matched <- model06.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model06.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model06.matched,
  weights = weights
)

# G computation to estimate marginal effects
model06.pred <- avg_predictions(
  model = model06.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 6, .before = trt_control)
model06.pred

# Estimation of average treatment effect
model06.comp <- avg_comparisons(
  model = model06.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 6, .before = term)
model06.comp # p < 0.001

# Plot
model06.plot <- model06.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "6. Northwestern Great Plains: Prescribed burn",
       x = NULL) +
  geom_signif(
    comparisons = list(c("Control", "Prescribed burn")),
    annotations = c("***")
  )
model06.plot




# Cold Deserts ------------------------------------------------------------

## Snake River Plain ------------------------------------------------------

### 7. Post-burn aerial seeding -------------------------------------------

# Filter data
model07.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Aerial Seeding" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post")) |>
  filter(EcoLvl3 == "Snake River Plain")

# PSM
model07.psm <- matchit(
  data = model07.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model07.psm
summary(model07.psm) # 153 treated matched

# Diagnostic love plot
model07.loveplot <- love.plot(model07.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "7. Snake River Plain: Post-burn aerial seeding") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model07.loveplot

# eCDF plots
bal.plot(model07.psm, which = "both", type = "ecdf")
bal.plot(model07.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model07.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model07.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model07.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model07.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model07.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model07.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model07.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model07.psm, which = "both")
bal.plot(model07.psm, var = "BareSoil_FH", which = "both")
bal.plot(model07.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model07.psm, var = "ForbCover_AH", which = "both")
bal.plot(model07.psm, var = "GramCover_AH", which = "both")
bal.plot(model07.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model07.psm, var = "Gap100plus", which = "both")
bal.plot(model07.psm, var = "CETWI", which = "both")
bal.plot(model07.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model07.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model07.psm, type = "qq")

# Matched data
model07.matched <- match_data(model07.psm)

# Create trt_control variable
model07.matched <- model07.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Post-burn aerial seeding", "Post-burn control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn aerial seeding")))

# Center and scale numeric variables
model07.matched <- model07.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model07.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model07.matched,
  weights = weights
)

# G computation to estimate marginal effects
model07.pred <- avg_predictions(
  model = model07.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 7, .before = trt_control)
model07.pred

# Estimation of average treatment effect
model07.comp <- avg_comparisons(
  model = model07.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 7, .before = term)
model07.comp # p = 0.011

# Plot
model07.plot <- model07.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "7. Snake River Plain: Post-burn aerial seeding",
       x = NULL) +
  geom_signif(
    comparisons = list(c("Post-burn control", "Post-burn aerial seeding")),
    annotations = c("**")
  )
model07.plot


### 8. Post burn aerial & drill seeding -----------------------------------

# Filter data
model08.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Aerial Seeding, Drill Seeding" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post")) |>
  filter(EcoLvl3 == "Snake River Plain")

# PSM
model08.psm <- matchit(
  data = model08.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model08.psm
summary(model08.psm) # 77 treated matched

# Diagnostic love plot
model08.loveplot <- love.plot(model08.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "8. Snake River Plain: Post-burn aerial & drill seeding") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model08.loveplot

# eCDF plots
bal.plot(model08.psm, which = "both", type = "ecdf")
bal.plot(model08.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model08.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model08.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model08.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model08.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model08.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model08.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model08.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model08.psm, which = "both")
bal.plot(model08.psm, var = "BareSoil_FH", which = "both")
bal.plot(model08.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model08.psm, var = "ForbCover_AH", which = "both")
bal.plot(model08.psm, var = "GramCover_AH", which = "both")
bal.plot(model08.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model08.psm, var = "Gap100plus", which = "both")
bal.plot(model08.psm, var = "CETWI", which = "both")
bal.plot(model08.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model08.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model08.psm, type = "qq")

# Matched data
model08.matched <- match_data(model08.psm)

# Create trt_control variable
model08.matched <- model08.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Post-burn aerial & drill seeding", "Post-burn control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn aerial & drill seeding")))

# Center and scale numeric variables
model08.matched <- model08.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model08.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model08.matched,
  weights = weights
)

# G computation to estimate marginal effects
model08.pred <- avg_predictions(
  model = model08.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 8, .before = trt_control)
model08.pred

# Estimation of average treatment effect
model08.comp <- avg_comparisons(
  model = model08.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 8, .before = term)
model08.comp

# Plot
model08.plot <- model08.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "8. Snake River Plain: Post-burn aerial & drill seeding",
       x = NULL)
model08.plot


### 9. Post-burn closure --------------------------------------------------

# Filter data
model09.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Closure" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post")) |>
  filter(EcoLvl3 == "Snake River Plain")

# PSM
model09.psm <- matchit(
  data = model09.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model09.psm
summary(model09.psm) # 87 treated matched

# Diagnostic love plot
model09.loveplot <- love.plot(model09.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "9. Snake River Plain: Post-burn closure") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model09.loveplot

# eCDF plots
bal.plot(model09.psm, which = "both", type = "ecdf")
bal.plot(model09.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model09.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model09.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model09.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model09.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model09.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model09.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model09.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model09.psm, which = "both")
bal.plot(model09.psm, var = "BareSoil_FH", which = "both")
bal.plot(model09.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model09.psm, var = "ForbCover_AH", which = "both")
bal.plot(model09.psm, var = "GramCover_AH", which = "both")
bal.plot(model09.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model09.psm, var = "Gap100plus", which = "both")
bal.plot(model09.psm, var = "CETWI", which = "both")
bal.plot(model09.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model09.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model09.psm, type = "qq")

# Matched data
model09.matched <- match_data(model09.psm)

# Create trt_control variable
model09.matched <- model09.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Post-burn closure", "Post-burn control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn closure")))

# Center and scale numeric variables
model09.matched <- model09.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model09.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model09.matched,
  weights = weights
)

# G computation to estimate marginal effects
model09.pred <- avg_predictions(
  model = model09.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 9, .before = trt_control)
model09.pred

# Estimation of average treatment effect
model09.comp <- avg_comparisons(
  model = model09.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 9, .before = term)
model09.comp # p < 0.001

# Plot
model09.plot <- model09.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "9. Snake River Plain: Post-burn closure",
       x = NULL) +
  geom_signif(
    comparisons = list(c("Post-burn control", "Post-burn closure")),
    annotations = c("***")
  )
model09.plot


### 10. Post-burn drill seeding -------------------------------------------------

# Filter data
model10.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Drill Seeding" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post")) |>
  filter(EcoLvl3 == "Snake River Plain")

# PSM
model10.psm <- matchit(
  data = model10.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model10.psm
summary(model10.psm) # 36 treated matched

# Diagnostic love plot
model10.loveplot <- love.plot(model10.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "10. Snake River Plain: Post-burn drill seeding") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model10.loveplot

# eCDF plots
bal.plot(model10.psm, which = "both", type = "ecdf")
bal.plot(model10.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model10.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model10.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model10.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model10.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model10.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model10.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model10.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model10.psm, which = "both")
bal.plot(model10.psm, var = "BareSoil_FH", which = "both")
bal.plot(model10.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model10.psm, var = "ForbCover_AH", which = "both")
bal.plot(model10.psm, var = "GramCover_AH", which = "both")
bal.plot(model10.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model10.psm, var = "Gap100plus", which = "both")
bal.plot(model10.psm, var = "CETWI", which = "both")
bal.plot(model10.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model10.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model10.psm, type = "qq")

# Matched data
model10.matched <- match_data(model10.psm)

# Create trt_control variable
model10.matched <- model10.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Post-burn drill seeding", "Post-burn control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn drill seeding")))

# Center and scale numeric variables
model10.matched <- model10.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model10.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model10.matched,
  weights = weights
)

# G computation to estimate marginal effects
model10.pred <- avg_predictions(
  model = model10.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 10, .before = trt_control)
model10.pred

# Estimation of average treatment effect
model10.comp <- avg_comparisons(
  model = model10.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 10, .before = term)
model10.comp

# Plot
model10.plot <- model10.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "10. Snake River Plain: Post-burn drill seeding",
       x = NULL)
model10.plot


### 11. Post-burn herbicide -----------------------------------------------

# Filter data
model11.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Herbicide" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post")) |>
  filter(EcoLvl3 == "Snake River Plain")

# PSM
model11.psm <- matchit(
  data = model11.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model11.psm
summary(model11.psm) # 28 treated matched

# Diagnostic love plot
model11.loveplot <- love.plot(model11.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "11. Snake River Plain: Post-burn herbicide") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model11.loveplot

# eCDF plots
bal.plot(model11.psm, which = "both", type = "ecdf")
bal.plot(model11.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model11.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model11.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model11.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model11.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model11.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model11.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model11.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model11.psm, which = "both")
bal.plot(model11.psm, var = "BareSoil_FH", which = "both")
bal.plot(model11.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model11.psm, var = "ForbCover_AH", which = "both")
bal.plot(model11.psm, var = "GramCover_AH", which = "both")
bal.plot(model11.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model11.psm, var = "Gap100plus", which = "both")
bal.plot(model11.psm, var = "CETWI", which = "both")
bal.plot(model11.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model11.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model11.psm, type = "qq")

# Matched data
model11.matched <- match_data(model11.psm)

# Create trt_control variable
model11.matched <- model11.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Post-burn herbicide", "Post-burn control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn herbicide")))

# Center and scale numeric variables
model11.matched <- model11.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model11.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model11.matched,
  weights = weights
)

# G computation to estimate marginal effects
model11.pred <- avg_predictions(
  model = model11.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 11, .before = trt_control)
model11.pred

# Estimation of average treatment effect
model11.comp <- avg_comparisons(
  model = model11.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 11, .before = term)
model11.comp

# Plot
model11.plot <- model11.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "11. Snake River Plain: Post-burn herbicide",
       x = NULL)
model11.plot



## Northern Basin and Range -----------------------------------------------

### 12. Drill seeding -----------------------------------------------------

# Filter data
model12.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Drill Seeding" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Northern Basin and Range")

# PSM
model12.psm <- matchit(
  data = model12.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model12.psm
summary(model12.psm) # 82 treated matched

# Diagnostic love plot
model12.loveplot <- love.plot(model12.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "12. Northern Basin and Range: Drill seeding") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model12.loveplot

# eCDF plots
bal.plot(model12.psm, which = "both", type = "ecdf")
bal.plot(model12.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model12.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model12.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model12.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model12.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model12.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model12.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model12.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model12.psm, which = "both")
bal.plot(model12.psm, var = "BareSoil_FH", which = "both")
bal.plot(model12.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model12.psm, var = "ForbCover_AH", which = "both")
bal.plot(model12.psm, var = "GramCover_AH", which = "both")
bal.plot(model12.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model12.psm, var = "Gap100plus", which = "both")
bal.plot(model12.psm, var = "CETWI", which = "both")
bal.plot(model12.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model12.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model12.psm, type = "qq")

# Matched data
model12.matched <- match_data(model12.psm)

# Create trt_control variable
model12.matched <- model12.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Drill seeding", "Control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Control", "Drill seeding")))

# Center and scale numeric variables
model12.matched <- model12.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model12.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model12.matched,
  weights = weights
)

# G computation to estimate marginal effects
model12.pred <- avg_predictions(
  model = model12.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 12, .before = trt_control)
model12.pred

# Estimation of average treatment effect
model12.comp <- avg_comparisons(
  model = model12.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 12, .before = term)
model12.comp

# Plot
model12.plot <- model12.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "12. Northern Basin and Range: Drill seeding",
       x = NULL)
model12.plot


### 13. Drill seeding & soil disturbance ----------------------------------

# Filter data
model13.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Drill Seeding, Soil Disturbance" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Northern Basin and Range")

# PSM
model13.psm <- matchit(
  data = model13.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model13.psm
summary(model13.psm) # 61 treated matched

# Diagnostic love plot
model13.loveplot <- love.plot(model13.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "13. Northern Basin and Range: Drill seeding & soil disturbance") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model13.loveplot

# eCDF plots
bal.plot(model13.psm, which = "both", type = "ecdf")
bal.plot(model13.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model13.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model13.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model13.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model13.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model13.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model13.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model13.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model13.psm, which = "both")
bal.plot(model13.psm, var = "BareSoil_FH", which = "both")
bal.plot(model13.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model13.psm, var = "ForbCover_AH", which = "both")
bal.plot(model13.psm, var = "GramCover_AH", which = "both")
bal.plot(model13.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model13.psm, var = "Gap100plus", which = "both")
bal.plot(model13.psm, var = "CETWI", which = "both")
bal.plot(model13.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model13.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model13.psm, type = "qq")

# Matched data
model13.matched <- match_data(model13.psm)

# Create trt_control variable
model13.matched <- model13.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Drill seeding & soil disturbance", "Control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Control", "Drill seeding & soil disturbance")))

# Center and scale numeric variables
model13.matched <- model13.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model13.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model13.matched,
  weights = weights
)

# G computation to estimate marginal effects
model13.pred <- avg_predictions(
  model = model13.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 13, .before = trt_control)
model13.pred

# Estimation of average treatment effect
model13.comp <- avg_comparisons(
  model = model13.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 13, .before = term)
model13.comp

# Plot
model13.plot <- model13.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "13. Northern Basin and Range: Drill seeding & soil disturbance",
       x = NULL)
model13.plot


### 14. Herbicide ---------------------------------------------------------

# Filter data
model14.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Herbicide" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Northern Basin and Range")

# PSM
model14.psm <- matchit(
  data = model14.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model14.psm
summary(model14.psm) # 259 treated matched

# Diagnostic love plot
model14.loveplot <- love.plot(model14.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "14. Northern Basin and Range: Herbicide") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model14.loveplot

# eCDF plots
bal.plot(model14.psm, which = "both", type = "ecdf")
bal.plot(model14.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model14.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model14.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model14.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model14.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model14.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model14.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model14.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model14.psm, which = "both")
bal.plot(model14.psm, var = "BareSoil_FH", which = "both")
bal.plot(model14.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model14.psm, var = "ForbCover_AH", which = "both")
bal.plot(model14.psm, var = "GramCover_AH", which = "both")
bal.plot(model14.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model14.psm, var = "Gap100plus", which = "both")
bal.plot(model14.psm, var = "CETWI", which = "both")
bal.plot(model14.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model14.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model14.psm, type = "qq")

# Matched data
model14.matched <- match_data(model14.psm)

# Create trt_control variable
model14.matched <- model14.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Herbicide", "Control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Control", "Herbicide")))

# Center and scale numeric variables
model14.matched <- model14.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model14.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model14.matched,
  weights = weights
)

# G computation to estimate marginal effects
model14.pred <- avg_predictions(
  model = model14.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 14, .before = trt_control)
model14.pred

# Estimation of average treatment effect
model14.comp <- avg_comparisons(
  model = model14.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 14, .before = term)
model14.comp

# Plot
model14.plot <- model14.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "14. Northern Basin and Range: Herbicide",
       x = NULL)
model14.plot


### 15. Prescribed burn ---------------------------------------------------

# Filter data
model15.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Prescribed Burn" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Northern Basin and Range")

# PSM
model15.psm <- matchit(
  data = model15.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model15.psm
summary(model15.psm) # 158 treated matched

# Diagnostic love plot
model15.loveplot <- love.plot(model15.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "15. Northern Basin and Range: Prescribed burn") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model15.loveplot

# eCDF plots
bal.plot(model15.psm, which = "both", type = "ecdf")
bal.plot(model15.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model15.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model15.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model15.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model15.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model15.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model15.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model15.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model15.psm, which = "both")
bal.plot(model15.psm, var = "BareSoil_FH", which = "both")
bal.plot(model15.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model15.psm, var = "ForbCover_AH", which = "both")
bal.plot(model15.psm, var = "GramCover_AH", which = "both")
bal.plot(model15.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model15.psm, var = "Gap100plus", which = "both")
bal.plot(model15.psm, var = "CETWI", which = "both")
bal.plot(model15.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model15.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model15.psm, type = "qq")

# Matched data
model15.matched <- match_data(model15.psm)

# Create trt_control variable
model15.matched <- model15.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Prescribed burn", "Control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Control", "Prescribed burn")))

# Center and scale numeric variables
model15.matched <- model15.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model15.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model15.matched,
  weights = weights
)

# G computation to estimate marginal effects
model15.pred <- avg_predictions(
  model = model15.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 15, .before = trt_control)
model15.pred

# Estimation of average treatment effect
model15.comp <- avg_comparisons(
  model = model15.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 15, .before = term)
model15.comp

# Plot
model15.plot <- model15.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "15. Northern Basin and Range: Prescribed burn",
       x = NULL)
model15.plot


### 16. Vegetation disturbance --------------------------------------------

# Filter data
model16.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Vegetation Disturbance" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Northern Basin and Range")

# PSM
model16.psm <- matchit(
  data = model16.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model16.psm
summary(model16.psm) # 105 treated matched

# Diagnostic love plot
model16.loveplot <- love.plot(model16.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "16. Northern Basin and Range: Vegetation disturbance") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model16.loveplot

# eCDF plots
bal.plot(model16.psm, which = "both", type = "ecdf")
bal.plot(model16.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model16.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model16.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model16.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model16.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model16.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model16.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model16.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model16.psm, which = "both")
bal.plot(model16.psm, var = "BareSoil_FH", which = "both")
bal.plot(model16.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model16.psm, var = "ForbCover_AH", which = "both")
bal.plot(model16.psm, var = "GramCover_AH", which = "both")
bal.plot(model16.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model16.psm, var = "Gap100plus", which = "both")
bal.plot(model16.psm, var = "CETWI", which = "both")
bal.plot(model16.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model16.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model16.psm, type = "qq")

# Matched data
model16.matched <- match_data(model16.psm)

# Create trt_control variable
model16.matched <- model16.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Vegetation disturbance", "Control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Control", "Vegetation disturbance")))

# Center and scale numeric variables
model16.matched <- model16.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model16.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model16.matched,
  weights = weights
)

# G computation to estimate marginal effects
model16.pred <- avg_predictions(
  model = model16.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 16, .before = trt_control)
model16.pred

# Estimation of average treatment effect
model16.comp <- avg_comparisons(
  model = model16.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 16, .before = term)
model16.comp # p = 0.003

# Plot
model16.plot <- model16.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 16,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "16. Northern Basin and Range: Vegetation disturbance",
       x = NULL) +
  geom_signif(
    comparisons = list(c("Control", "Vegetation disturbance")),
    annotations = c("**")
  )
model16.plot


### 17. Post-burn aerial seeding ------------------------------------------

# Filter data
model17.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Aerial Seeding" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post")) |>
  filter(EcoLvl3 == "Northern Basin and Range")

# PSM
model17.psm <- matchit(
  data = model17.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model17.psm
summary(model17.psm) # 542 treated matched

# Diagnostic love plot
model17.loveplot <- love.plot(model17.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "17. Northern Basin and Range: Post-burn aerial seeding") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model17.loveplot

# eCDF plots
bal.plot(model17.psm, which = "both", type = "ecdf")
bal.plot(model17.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model17.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model17.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model17.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model17.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model17.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model17.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model17.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model17.psm, which = "both")
bal.plot(model17.psm, var = "BareSoil_FH", which = "both")
bal.plot(model17.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model17.psm, var = "ForbCover_AH", which = "both")
bal.plot(model17.psm, var = "GramCover_AH", which = "both")
bal.plot(model17.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model17.psm, var = "Gap100plus", which = "both")
bal.plot(model17.psm, var = "CETWI", which = "both")
bal.plot(model17.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model17.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model17.psm, type = "qq")

# Matched data
model17.matched <- match_data(model17.psm)

# Create trt_control variable
model17.matched <- model17.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Post-burn aerial seeding", "Post-burn control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn aerial seeding")))

# Center and scale numeric variables
model17.matched <- model17.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model17.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model17.matched,
  weights = weights
)

# G computation to estimate marginal effects
model17.pred <- avg_predictions(
  model = model17.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 17, .before = trt_control)
model17.pred

# Estimation of average treatment effect
model17.comp <- avg_comparisons(
  model = model17.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 17, .before = term)
model17.comp

# Plot
model17.plot <- model17.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "17. Northern Basin and Range: Post-burn aerial seeding",
       x = NULL)
model17.plot


### 18. Post-burn aerial & drill seeding ----------------------------------

# Filter data
model18.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Aerial Seeding, Drill Seeding" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post")) |>
  filter(EcoLvl3 == "Northern Basin and Range")

# PSM
model18.psm <- matchit(
  data = model18.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model18.psm
summary(model18.psm) # 88 treated matched

# Diagnostic love plot
model18.loveplot <- love.plot(model18.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "18. Northern Basin and Range: Post-burn aerial & drill seeding") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model18.loveplot

# eCDF plots
bal.plot(model18.psm, which = "both", type = "ecdf")
bal.plot(model18.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model18.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model18.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model18.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model18.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model18.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model18.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model18.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model18.psm, which = "both")
bal.plot(model18.psm, var = "BareSoil_FH", which = "both")
bal.plot(model18.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model18.psm, var = "ForbCover_AH", which = "both")
bal.plot(model18.psm, var = "GramCover_AH", which = "both")
bal.plot(model18.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model18.psm, var = "Gap100plus", which = "both")
bal.plot(model18.psm, var = "CETWI", which = "both")
bal.plot(model18.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model18.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model18.psm, type = "qq")

# Matched data
model18.matched <- match_data(model18.psm)

# Create trt_control variable
model18.matched <- model18.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Post-burn aerial & drill seeding", "Post-burn control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn aerial & drill seeding")))

# Center and scale numeric variables
model18.matched <- model18.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model18.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model18.matched,
  weights = weights
)

# G computation to estimate marginal effects
model18.pred <- avg_predictions(
  model = model18.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 18, .before = trt_control)
model18.pred

# Estimation of average treatment effect
model18.comp <- avg_comparisons(
  model = model18.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 18, .before = term)
model18.comp

# Plot
model18.plot <- model18.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "18. Northern Basin and Range: Post-burn aerial & drill seeding",
       x = NULL)
model18.plot


### 19. Post-burn closure -------------------------------------------------

# Filter data
model19.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Closure" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post")) |>
  filter(EcoLvl3 == "Northern Basin and Range")

# PSM
model19.psm <- matchit(
  data = model19.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model19.psm
summary(model19.psm) # 151 treated matched

# Diagnostic love plot
model19.loveplot <- love.plot(model19.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "19. Northern Basin and Range: Post-burn closure") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model19.loveplot

# eCDF plots
bal.plot(model19.psm, which = "both", type = "ecdf")
bal.plot(model19.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model19.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model19.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model19.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model19.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model19.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model19.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model19.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model19.psm, which = "both")
bal.plot(model19.psm, var = "BareSoil_FH", which = "both")
bal.plot(model19.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model19.psm, var = "ForbCover_AH", which = "both")
bal.plot(model19.psm, var = "GramCover_AH", which = "both")
bal.plot(model19.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model19.psm, var = "Gap100plus", which = "both")
bal.plot(model19.psm, var = "CETWI", which = "both")
bal.plot(model19.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model19.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model19.psm, type = "qq")

# Matched data
model19.matched <- match_data(model19.psm)

# Create trt_control variable
model19.matched <- model19.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Post-burn closure", "Post-burn control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn closure")))

# Center and scale numeric variables
model19.matched <- model19.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model19.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model19.matched,
  weights = weights
)

# G computation to estimate marginal effects
model19.pred <- avg_predictions(
  model = model19.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 19, .before = trt_control)
model19.pred

# Estimation of average treatment effect
model19.comp <- avg_comparisons(
  model = model19.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 19, .before = term)
model19.comp # p < 0.001

# Plot
model19.plot <- model19.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "19. Northern Basin and Range: Post-burn closure",
       x = NULL) +
  geom_signif(
    comparisons = list(c("Post-burn control", "Post-burn closure")),
    annotations = c("***")
  )
model19.plot


### 20. Post-burn drill seeding -------------------------------------------------

# Filter data
model20.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Drill Seeding" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post")) |>
  filter(EcoLvl3 == "Northern Basin and Range")

# PSM
model20.psm <- matchit(
  data = model20.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model20.psm
summary(model20.psm) # 217 treated matched

# Diagnostic love plot
model20.loveplot <- love.plot(model20.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "20. Northern Basin and Range: Post-burn drill seeding") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model20.loveplot

# eCDF plots
bal.plot(model20.psm, which = "both", type = "ecdf")
bal.plot(model20.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model20.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model20.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model20.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model20.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model20.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model20.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model20.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model20.psm, which = "both")
bal.plot(model20.psm, var = "BareSoil_FH", which = "both")
bal.plot(model20.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model20.psm, var = "ForbCover_AH", which = "both")
bal.plot(model20.psm, var = "GramCover_AH", which = "both")
bal.plot(model20.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model20.psm, var = "Gap100plus", which = "both")
bal.plot(model20.psm, var = "CETWI", which = "both")
bal.plot(model20.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model20.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model20.psm, type = "qq")

# Matched data
model20.matched <- match_data(model20.psm)

# Create trt_control variable
model20.matched <- model20.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Post-burn drill seeding", "Post-burn control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn drill seeding")))

# Center and scale numeric variables
model20.matched <- model20.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model20.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model20.matched,
  weights = weights
)

# G computation to estimate marginal effects
model20.pred <- avg_predictions(
  model = model20.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 20, .before = trt_control)
model20.pred

# Estimation of average treatment effect
model20.comp <- avg_comparisons(
  model = model20.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 20, .before = term)
model20.comp

# Plot
model20.plot <- model20.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "20. Northern Basin and Range: Post-burn drill seeding",
       x = NULL)
model20.plot


### 21. Post-burn herbicide -----------------------------------------------

# Filter data
model21.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Herbicide" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post")) |>
  filter(EcoLvl3 == "Northern Basin and Range")

# PSM
model21.psm <- matchit(
  data = model21.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model21.psm
summary(model21.psm) # 370 treated matched

# Diagnostic love plot
model21.loveplot <- love.plot(model21.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "21. Northern Basin and Range: Post-burn herbicide") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model21.loveplot

# eCDF plots
bal.plot(model21.psm, which = "both", type = "ecdf")
bal.plot(model21.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model21.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model21.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model21.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model21.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model21.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model21.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model21.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model21.psm, which = "both")
bal.plot(model21.psm, var = "BareSoil_FH", which = "both")
bal.plot(model21.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model21.psm, var = "ForbCover_AH", which = "both")
bal.plot(model21.psm, var = "GramCover_AH", which = "both")
bal.plot(model21.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model21.psm, var = "Gap100plus", which = "both")
bal.plot(model21.psm, var = "CETWI", which = "both")
bal.plot(model21.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model21.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model21.psm, type = "qq")

# Matched data
model21.matched <- match_data(model21.psm)

# Create trt_control variable
model21.matched <- model21.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Post-burn herbicide", "Post-burn control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn herbicide")))

# Center and scale numeric variables
model21.matched <- model21.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model21.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model21.matched,
  weights = weights
)

# G computation to estimate marginal effects
model21.pred <- avg_predictions(
  model = model21.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 21, .before = trt_control)
model21.pred

# Estimation of average treatment effect
model21.comp <- avg_comparisons(
  model = model21.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 21, .before = term)
model21.comp

# Plot
model21.plot <- model21.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "21. Northern Basin and Range: Post-burn herbicide",
       x = NULL)
model21.plot


### 22. Post-burn seedling planting ---------------------------------------

# Filter data
model22.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Seedling Planting" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post")) |>
  filter(EcoLvl3 == "Northern Basin and Range")

# PSM
model22.psm <- matchit(
  data = model22.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model22.psm
summary(model22.psm) # 56 treated matched

# Diagnostic love plot
model22.loveplot <- love.plot(model22.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "22. Northern Basin and Range: Post-burn seedling planting") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model22.loveplot

# eCDF plots
bal.plot(model22.psm, which = "both", type = "ecdf")
bal.plot(model22.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model22.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model22.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model22.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model22.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model22.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model22.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model22.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model22.psm, which = "both")
bal.plot(model22.psm, var = "BareSoil_FH", which = "both")
bal.plot(model22.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model22.psm, var = "ForbCover_AH", which = "both")
bal.plot(model22.psm, var = "GramCover_AH", which = "both")
bal.plot(model22.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model22.psm, var = "Gap100plus", which = "both")
bal.plot(model22.psm, var = "CETWI", which = "both")
bal.plot(model22.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model22.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model22.psm, type = "qq")

# Matched data
model22.matched <- match_data(model22.psm)

# Create trt_control variable
model22.matched <- model22.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Post-burn seedling planting", "Post-burn control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn seedling planting")))

# Center and scale numeric variables
model22.matched <- model22.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model22.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model22.matched,
  weights = weights
)

# G computation to estimate marginal effects
model22.pred <- avg_predictions(
  model = model22.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 22, .before = trt_control)
model22.pred

# Estimation of average treatment effect
model22.comp <- avg_comparisons(
  model = model22.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 22, .before = term)
model22.comp # p = 0.03

# Plot
model22.plot <- model22.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "22. Northern Basin and Range: Post-burn seedling planting",
       x = NULL) + 
  geom_signif(
    comparisons = list(c("Post-burn control", "Post-burn seedling planting")),
    annotations = c("*")
  )
model22.plot



## Central Basin and Range ------------------------------------------------

### 23. Drill seeding & soil disturbance ----------------------------------

# Filter data
model23.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Drill Seeding, Soil Disturbance" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never burned")) |>
  filter(EcoLvl3 == "Central Basin and Range")

# PSM
model23.psm <- matchit(
  data = model23.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
) # warning: fitted probabilities numerically 0 or 1 occurred 
model23.psm
summary(model23.psm) # 36 treated matched

# Diagnostic love plot
model23.loveplot <- love.plot(model23.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "23. Central Basin and Range: Drill seeding & soil disturbance") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model23.loveplot

# eCDF plots
bal.plot(model23.psm, which = "both", type = "ecdf")
bal.plot(model23.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model23.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model23.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model23.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model23.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model23.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model23.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model23.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model23.psm, which = "both")
bal.plot(model23.psm, var = "BareSoil_FH", which = "both")
bal.plot(model23.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model23.psm, var = "ForbCover_AH", which = "both")
bal.plot(model23.psm, var = "GramCover_AH", which = "both")
bal.plot(model23.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model23.psm, var = "Gap100plus", which = "both")
bal.plot(model23.psm, var = "CETWI", which = "both")
bal.plot(model23.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model23.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model23.psm, type = "qq")

# Matched data
model23.matched <- match_data(model23.psm)

# Create trt_control variable
model23.matched <- model23.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Drill seeding & soil disturbance", "Control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Control", "Drill seeding & soil disturbance")))

# Center and scale numeric variables
model23.matched <- model23.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model23.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model23.matched,
  weights = weights
)

# G computation to estimate marginal effects
model23.pred <- avg_predictions(
  model = model23.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 23, .before = trt_control)
model23.pred # NS for both control & treated (not different from 0)

# Estimation of average treatment effect
model23.comp <- avg_comparisons(
  model = model23.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 23, .before = term)
model23.comp

# Plot
model23.plot <- model23.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "23. Central Basin and Range: Drill seeding & soil disturbance",
       x = NULL)
model23.plot


### 24. Prescribed burn ---------------------------------------------------

# Filter data
model24.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Prescribed Burn" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never burned")) |>
  filter(EcoLvl3 == "Central Basin and Range")

# PSM
model24.psm <- matchit(
  data = model24.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model24.psm
summary(model24.psm) # 31 treated matched

# Diagnostic love plot
model24.loveplot <- love.plot(model24.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "24. Central Basin and Range: Prescribed burn") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model24.loveplot

# eCDF plots
bal.plot(model24.psm, which = "both", type = "ecdf")
bal.plot(model24.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model24.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model24.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model24.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model24.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model24.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model24.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model24.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model24.psm, which = "both")
bal.plot(model24.psm, var = "BareSoil_FH", which = "both")
bal.plot(model24.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model24.psm, var = "ForbCover_AH", which = "both")
bal.plot(model24.psm, var = "GramCover_AH", which = "both")
bal.plot(model24.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model24.psm, var = "Gap100plus", which = "both")
bal.plot(model24.psm, var = "CETWI", which = "both")
bal.plot(model24.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model24.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model24.psm, type = "qq")

# Matched data
model24.matched <- match_data(model24.psm)

# Create trt_control variable
model24.matched <- model24.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Prescribed burn", "Control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Control", "Prescribed burn")))

# Center and scale numeric variables
model24.matched <- model24.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model24.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model24.matched,
  weights = weights
)

# G computation to estimate marginal effects
model24.pred <- avg_predictions(
  model = model24.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 24, .before = trt_control)
model24.pred

# Estimation of average treatment effect
model24.comp <- avg_comparisons(
  model = model24.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 24, .before = term)
model24.comp

# Plot
model24.plot <- model24.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "24. Central Basin and Range: Prescribed burn",
       x = NULL)
model24.plot


### 25. Vegetation disturbance --------------------------------------------

# Filter data
model25.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Vegetation Disturbance" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never burned")) |>
  filter(EcoLvl3 == "Central Basin and Range")

# PSM
model25.psm <- matchit(
  data = model25.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model25.psm
summary(model25.psm) # 99 treated matched

# Diagnostic love plot
model25.loveplot <- love.plot(model25.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "25. Central Basin and Range: Vegetation disturbance") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model25.loveplot

# eCDF plots
bal.plot(model25.psm, which = "both", type = "ecdf")
bal.plot(model25.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model25.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model25.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model25.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model25.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model25.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model25.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model25.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model25.psm, which = "both")
bal.plot(model25.psm, var = "BareSoil_FH", which = "both")
bal.plot(model25.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model25.psm, var = "ForbCover_AH", which = "both")
bal.plot(model25.psm, var = "GramCover_AH", which = "both")
bal.plot(model25.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model25.psm, var = "Gap100plus", which = "both")
bal.plot(model25.psm, var = "CETWI", which = "both")
bal.plot(model25.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model25.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model25.psm, type = "qq")

# Matched data
model25.matched <- match_data(model25.psm)

# Create trt_control variable
model25.matched <- model25.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Vegetation disturbance", "Control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Control", "Vegetation disturbance")))

# Center and scale numeric variables
model25.matched <- model25.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model25.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model25.matched,
  weights = weights
)

# G computation to estimate marginal effects
model25.pred <- avg_predictions(
  model = model25.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 25, .before = trt_control)
model25.pred

# Estimation of average treatment effect
model25.comp <- avg_comparisons(
  model = model25.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 25, .before = term)
model25.comp # p = 0.001

# Plot
model25.plot <- model25.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "25. Central Basin and Range: Vegetation disturbance",
       x = NULL) +
  geom_signif(
    comparisons = list(c("Control", "Vegetation disturbance")),
    annotations = c("**")
  )
model25.plot


### 26. Post-burn aerial seeding ------------------------------------------

# Filter data
model26.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Aerial Seeding" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post-burn")) |>
  filter(EcoLvl3 == "Central Basin and Range")

# PSM
model26.psm <- matchit(
  data = model26.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model26.psm
summary(model26.psm) # 345 treated matched

# Diagnostic love plot
model26.loveplot <- love.plot(model26.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "26. Central Basin and Range: Post-burn aerial seeding") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model26.loveplot

# eCDF plots
bal.plot(model26.psm, which = "both", type = "ecdf")
bal.plot(model26.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model26.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model26.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model26.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model26.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model26.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model26.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model26.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model26.psm, which = "both")
bal.plot(model26.psm, var = "BareSoil_FH", which = "both")
bal.plot(model26.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model26.psm, var = "ForbCover_AH", which = "both")
bal.plot(model26.psm, var = "GramCover_AH", which = "both")
bal.plot(model26.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model26.psm, var = "Gap100plus", which = "both")
bal.plot(model26.psm, var = "CETWI", which = "both")
bal.plot(model26.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model26.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model26.psm, type = "qq")

# Matched data
model26.matched <- match_data(model26.psm)

# Create trt_control variable
model26.matched <- model26.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Post-burn aerial seeding", "Post-burn control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn aerial seeding")))

# Center and scale numeric variables
model26.matched <- model26.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model26.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model26.matched,
  weights = weights
)

# G computation to estimate marginal effects
model26.pred <- avg_predictions(
  model = model26.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 26, .before = trt_control)
model26.pred

# Estimation of average treatment effect
model26.comp <- avg_comparisons(
  model = model26.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 26, .before = term)
model26.comp

# Plot
model26.plot <- model26.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "26. Central Basin and Range: Post-burn aerial seeding",
       x = NULL)
model26.plot


### 27. Post-burn drill seeding -------------------------------------------

# Filter data
model27.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Drill Seeding" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post-burn")) |>
  filter(EcoLvl3 == "Central Basin and Range")

# PSM
model27.psm <- matchit(
  data = model27.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model27.psm
summary(model27.psm) # 88 treated matched

# Diagnostic love plot
model27.loveplot <- love.plot(model27.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "27. Central Basin and Range: Post-burn drill seeding") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model27.loveplot

# eCDF plots
bal.plot(model27.psm, which = "both", type = "ecdf")
bal.plot(model27.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model27.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model27.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model27.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model27.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model27.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model27.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model27.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model27.psm, which = "both")
bal.plot(model27.psm, var = "BareSoil_FH", which = "both")
bal.plot(model27.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model27.psm, var = "ForbCover_AH", which = "both")
bal.plot(model27.psm, var = "GramCover_AH", which = "both")
bal.plot(model27.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model27.psm, var = "Gap100plus", which = "both")
bal.plot(model27.psm, var = "CETWI", which = "both")
bal.plot(model27.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model27.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model27.psm, type = "qq")

# Matched data
model27.matched <- match_data(model27.psm)

# Create trt_control variable
model27.matched <- model27.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Post-burn drill seeding", "Post-burn control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn drill seeding")))

# Center and scale numeric variables
model27.matched <- model27.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model27.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model27.matched,
  weights = weights
)

# G computation to estimate marginal effects
model27.pred <- avg_predictions(
  model = model27.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 27, .before = trt_control)
model27.pred # Control NS (not different from 0)

# Estimation of average treatment effect
model27.comp <- avg_comparisons(
  model = model27.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 27, .before = term)
model27.comp # p = 0.01

# Plot
model27.plot <- model27.pred |>
  
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "27. Central Basin and Range: Post-burn drill seeding",
       x = NULL) +
  geom_signif(
    comparisons = list(c("Post-burn control", "Post-burn drill seeding")),
    annotations = c("**")
  )
model27.plot


### 28. Post-burn ground seeding ------------------------------------------

# Filter data
model28.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Ground Seeding" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post-burn")) |>
  filter(EcoLvl3 == "Central Basin and Range")

# PSM
model28.psm <- matchit(
  data = model28.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model28.psm
summary(model28.psm) # 38 treated matched

# Diagnostic love plot
model28.loveplot <- love.plot(model28.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "28. Central Basin and Range: Post-burn ground seeding") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model28.loveplot

# eCDF plots
bal.plot(model28.psm, which = "both", type = "ecdf")
bal.plot(model28.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model28.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model28.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model28.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model28.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model28.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model28.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model28.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model28.psm, which = "both")
bal.plot(model28.psm, var = "BareSoil_FH", which = "both")
bal.plot(model28.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model28.psm, var = "ForbCover_AH", which = "both")
bal.plot(model28.psm, var = "GramCover_AH", which = "both")
bal.plot(model28.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model28.psm, var = "Gap100plus", which = "both")
bal.plot(model28.psm, var = "CETWI", which = "both")
bal.plot(model28.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model28.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model28.psm, type = "qq")

# Matched data
model28.matched <- match_data(model28.psm)

# Create trt_control variable
model28.matched <- model28.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Post-burn ground seeding", "Post-burn control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn ground seeding")))

# Center and scale numeric variables
model28.matched <- model28.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model28.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model28.matched,
  weights = weights
)

# G computation to estimate marginal effects
model28.pred <- avg_predictions(
  model = model28.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 28, .before = trt_control)
model28.pred # NS for both control & treated (not different from 0)

# Estimation of average treatment effect
model28.comp <- avg_comparisons(
  model = model28.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 28, .before = term)
model28.comp

# Plot
model28.plot <- model28.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "28. Central Basin and Range: Post-burn ground seeding",
       x = NULL)
model28.plot


### 29. Post-burn herbicide -----------------------------------------------

# Filter data
model29.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Herbicide" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post-burn")) |>
  filter(EcoLvl3 == "Central Basin and Range")

# PSM
model29.psm <- matchit(
  data = model29.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model29.psm
summary(model29.psm) # 80 treated matched

# Diagnostic love plot
model29.loveplot <- love.plot(model29.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "29. Central Basin and Range: Post-burn herbicide") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model29.loveplot

# eCDF plots
bal.plot(model29.psm, which = "both", type = "ecdf")
bal.plot(model29.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model29.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model29.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model29.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model29.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model29.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model29.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model29.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model29.psm, which = "both")
bal.plot(model29.psm, var = "BareSoil_FH", which = "both")
bal.plot(model29.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model29.psm, var = "ForbCover_AH", which = "both")
bal.plot(model29.psm, var = "GramCover_AH", which = "both")
bal.plot(model29.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model29.psm, var = "Gap100plus", which = "both")
bal.plot(model29.psm, var = "CETWI", which = "both")
bal.plot(model29.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model29.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model29.psm, type = "qq")

# Matched data
model29.matched <- match_data(model29.psm)

# Create trt_control variable
model29.matched <- model29.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Post-burn herbicide", "Post-burn control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn herbicide")))

# Center and scale numeric variables
model29.matched <- model29.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model29.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model29.matched,
  weights = weights
)

# G computation to estimate marginal effects
model29.pred <- avg_predictions(
  model = model29.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 29, .before = trt_control)
model29.pred # NS for both control & treated (not different from 0)

# Estimation of average treatment effect
model29.comp <- avg_comparisons(
  model = model29.lm,
  variables = "trt_control",
  vcov = ~subclass
) |>
  mutate(Model = 29, .before = term)
model29.comp

# Plot
model29.plot <- model29.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "29. Central Basin and Range: Post-burn herbicide",
       x = NULL)
model29.plot



## Wyoming Basin ----------------------------------------------------------

### 30. Prescribed burn ---------------------------------------------------

# Filter data
model30.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Prescribed Burn" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Wyoming Basin")

# PSM
model30.psm <- matchit(
  data = model30.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model30.psm
summary(model30.psm) # 37 treated matched

# Diagnostic love plot
model30.loveplot <- love.plot(model30.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "30. Wyoming Basin: Prescribed burn") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model30.loveplot

# eCDF plots
bal.plot(model30.psm, which = "both", type = "ecdf")
bal.plot(model30.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model30.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model30.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model30.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model30.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model30.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model30.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model30.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model30.psm, which = "both")
bal.plot(model30.psm, var = "BareSoil_FH", which = "both")
bal.plot(model30.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model30.psm, var = "ForbCover_AH", which = "both")
bal.plot(model30.psm, var = "GramCover_AH", which = "both")
bal.plot(model30.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model30.psm, var = "Gap100plus", which = "both")
bal.plot(model30.psm, var = "CETWI", which = "both")
bal.plot(model30.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model30.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model30.psm, type = "qq")

# Matched data
model30.matched <- match_data(model30.psm)

# Create trt_control variable
model30.matched <- model30.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Prescribed burn", "Control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Control", "Prescribed burn")))

# Center and scale numeric variables
model30.matched <- model30.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model30.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model30.matched,
  weights = weights
)

# G computation to estimate marginal effects
model30.pred <- avg_predictions(
  model = model30.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 30, .before = trt_control)
model30.pred

# Estimation of average treatment effect
model30.comp <- avg_comparisons(
  model = model30.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 30, .before = term)
model30.comp

# Plot
model30.plot <- model30.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "30. Wyoming Basin: Prescribed burn",
       x = NULL)
model30.plot



## Colorado Plateaus ------------------------------------------------------

### 31. Aerial seeding & soil disturbance ---------------------------------

# Filter data
model31.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Aerial Seeding, Soil Disturbance" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never burned")) |>
  filter(EcoLvl3 == "Colorado Plateaus")

# PSM
model31.psm <- matchit(
  data = model31.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model31.psm
summary(model31.psm) # 60 treated matched

# Diagnostic love plot
model31.loveplot <- love.plot(model31.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "31. Colorado Plateaus: Aerial seeding & soil disturbance") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model31.loveplot

# eCDF plots
bal.plot(model31.psm, which = "both", type = "ecdf")
bal.plot(model31.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model31.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model31.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model31.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model31.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model31.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model31.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model31.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model31.psm, which = "both")
bal.plot(model31.psm, var = "BareSoil_FH", which = "both")
bal.plot(model31.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model31.psm, var = "ForbCover_AH", which = "both")
bal.plot(model31.psm, var = "GramCover_AH", which = "both")
bal.plot(model31.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model31.psm, var = "Gap100plus", which = "both")
bal.plot(model31.psm, var = "CETWI", which = "both")
bal.plot(model31.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model31.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model31.psm, type = "qq")

# Matched data
model31.matched <- match_data(model31.psm)

# Create trt_control variable
model31.matched <- model31.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Aerial seeding & soil disturbance", "Control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Control", "Aerial seeding & soil disturbance")))

# Center and scale numeric variables
model31.matched <- model31.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model31.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model31.matched,
  weights = weights
)

# G computation to estimate marginal effects
model31.pred <- avg_predictions(
  model = model31.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 31, .before = trt_control)
model31.pred

# Estimation of average treatment effect
model31.comp <- avg_comparisons(
  model = model31.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 31, .before = term)
model31.comp

# Plot
model31.plot <- model31.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "31. Colorado Plateaus: Aerial seeding & soil disturbance",
       x = NULL)
model31.plot


### 32. Herbicide ---------------------------------------------------------

# Filter data
model32.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Herbicide" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Colorado Plateaus")

# PSM
model32.psm <- matchit(
  data = model32.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model32.psm
summary(model32.psm) # 47 treated matched

# Diagnostic love plot
model32.loveplot <- love.plot(model32.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "32. Colorado Plateaus: Herbicide") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model32.loveplot

# eCDF plots
bal.plot(model32.psm, which = "both", type = "ecdf")
bal.plot(model32.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model32.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model32.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model32.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model32.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model32.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model32.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model32.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model32.psm, which = "both")
bal.plot(model32.psm, var = "BareSoil_FH", which = "both")
bal.plot(model32.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model32.psm, var = "ForbCover_AH", which = "both")
bal.plot(model32.psm, var = "GramCover_AH", which = "both")
bal.plot(model32.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model32.psm, var = "Gap100plus", which = "both")
bal.plot(model32.psm, var = "CETWI", which = "both")
bal.plot(model32.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model32.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model32.psm, type = "qq")

# Matched data
model32.matched <- match_data(model32.psm)

# Create trt_control variable
model32.matched <- model32.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Herbicide", "Control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Control", "Herbicide")))

# Center and scale numeric variables
model32.matched <- model32.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model32.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model32.matched,
  weights = weights
)

# G computation to estimate marginal effects
model32.pred <- avg_predictions(
  model = model32.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 32, .before = trt_control)
model32.pred

# Estimation of average treatment effect
model32.comp <- avg_comparisons(
  model = model32.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 32, .before = term)
model32.comp

# Plot
model32.plot <- model32.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "32. Colorado Plateaus: Herbicide",
       x = NULL)
model32.plot


### 33. Prescribed burn ---------------------------------------------------

# Filter data
model33.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Prescribed Burn" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Colorado Plateaus")

# PSM
model33.psm <- matchit(
  data = model33.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model33.psm
summary(model33.psm) # 76 treated matched

# Diagnostic love plot
model33.loveplot <- love.plot(model33.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "33. Colorado Plateaus: Prescribed burn") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model33.loveplot

# eCDF plots
bal.plot(model33.psm, which = "both", type = "ecdf")
bal.plot(model33.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model33.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model33.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model33.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model33.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model33.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model33.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model33.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model33.psm, which = "both")
bal.plot(model33.psm, var = "BareSoil_FH", which = "both")
bal.plot(model33.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model33.psm, var = "ForbCover_AH", which = "both")
bal.plot(model33.psm, var = "GramCover_AH", which = "both")
bal.plot(model33.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model33.psm, var = "Gap100plus", which = "both")
bal.plot(model33.psm, var = "CETWI", which = "both")
bal.plot(model33.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model33.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model33.psm, type = "qq")

# Matched data
model33.matched <- match_data(model33.psm)

# Create trt_control variable
model33.matched <- model33.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Prescribed burn", "Control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Control", "Prescribed burn")))

# Center and scale numeric variables
model33.matched <- model33.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model33.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model33.matched,
  weights = weights
)

# G computation to estimate marginal effects
model33.pred <- avg_predictions(
  model = model33.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 33, .before = trt_control)
model33.pred

# Estimation of average treatment effect
model33.comp <- avg_comparisons(
  model = model33.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 33, .before = term)
model33.comp # p = 0.014

# Plot
model33.plot <- model33.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "33. Colorado Plateaus: Prescribed burn",
       x = NULL) + 
  geom_signif(
    comparisons = list(c("Control", "Prescribed burn")),
    annotations = c("*")
  )
model33.plot


### 34. Soil disturbance --------------------------------------------------

# Filter data
model34.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Soil Disturbance" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Colorado Plateaus")

# PSM
model34.psm <- matchit(
  data = model34.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model34.psm
summary(model34.psm) # 40 treated matched

# Diagnostic love plot
model34.loveplot <- love.plot(model34.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "34. Colorado Plateaus: Soil disturbance") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model34.loveplot

# eCDF plots
bal.plot(model34.psm, which = "both", type = "ecdf")
bal.plot(model34.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model34.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model34.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model34.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model34.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model34.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model34.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model34.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model34.psm, which = "both")
bal.plot(model34.psm, var = "BareSoil_FH", which = "both")
bal.plot(model34.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model34.psm, var = "ForbCover_AH", which = "both")
bal.plot(model34.psm, var = "GramCover_AH", which = "both")
bal.plot(model34.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model34.psm, var = "Gap100plus", which = "both")
bal.plot(model34.psm, var = "CETWI", which = "both")
bal.plot(model34.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model34.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model34.psm, type = "qq")

# Matched data
model34.matched <- match_data(model34.psm)

# Create trt_control variable
model34.matched <- model34.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Soil disturbance", "Control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Control", "Soil disturbance")))

# Center and scale numeric variables
model34.matched <- model34.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model34.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model34.matched,
  weights = weights
)

# G computation to estimate marginal effects
model34.pred <- avg_predictions(
  model = model34.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 34, .before = trt_control)
model34.pred

# Estimation of average treatment effect
model34.comp <- avg_comparisons(
  model = model34.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 34, .before = term)
model34.comp

# Plot
model34.plot <- model34.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "34. Colorado Plateaus: Soil disturbance",
       x = NULL)
model34.plot


### 35. Vegetation disturbance --------------------------------------------

# Filter data
model35.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Vegetation Disturbance" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Colorado Plateaus")

# PSM
model35.psm <- matchit(
  data = model35.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model35.psm
summary(model35.psm) # 34 treated matched

# Diagnostic love plot
model35.loveplot <- love.plot(model35.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "35. Colorado Plateaus: Vegetation disturbance") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model35.loveplot

# eCDF plots
bal.plot(model35.psm, which = "both", type = "ecdf")
bal.plot(model35.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model35.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model35.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model35.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model35.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model35.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model35.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model35.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model35.psm, which = "both")
bal.plot(model35.psm, var = "BareSoil_FH", which = "both")
bal.plot(model35.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model35.psm, var = "ForbCover_AH", which = "both")
bal.plot(model35.psm, var = "GramCover_AH", which = "both")
bal.plot(model35.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model35.psm, var = "Gap100plus", which = "both")
bal.plot(model35.psm, var = "CETWI", which = "both")
bal.plot(model35.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model35.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model35.psm, type = "qq")

# Matched data
model35.matched <- match_data(model35.psm)

# Create trt_control variable
model35.matched <- model35.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Vegetation disturbance", "Control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Control", "Vegetation disturbance")))

# Center and scale numeric variables
model35.matched <- model35.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model35.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model35.matched,
  weights = weights
)

# G computation to estimate marginal effects
model35.pred <- avg_predictions(
  model = model35.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 35, .before = trt_control)
model35.pred

# Estimation of average treatment effect
model35.comp <- avg_comparisons(
  model = model35.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 35, .before = term)
model35.comp

# Plot
model35.plot <- model35.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "35. Colorado Plateaus: Vegetation disturbance",
       x = NULL)
model35.plot


### 36. Post-burn aerial seeding ------------------------------------------

# Filter data
model36.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Aerial Seeding" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post")) |>
  filter(EcoLvl3 == "Colorado Plateaus")

# PSM
model36.psm <- matchit(
  data = model36.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model36.psm
summary(model36.psm) # 58 treated matched

# Diagnostic love plot
model36.loveplot <- love.plot(model36.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "36. Colorado Plateaus: Post-burn aerial seeding") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model36.loveplot

# eCDF plots
bal.plot(model36.psm, which = "both", type = "ecdf")
bal.plot(model36.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model36.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model36.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model36.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model36.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model36.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model36.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model36.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model36.psm, which = "both")
bal.plot(model36.psm, var = "BareSoil_FH", which = "both")
bal.plot(model36.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model36.psm, var = "ForbCover_AH", which = "both")
bal.plot(model36.psm, var = "GramCover_AH", which = "both")
bal.plot(model36.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model36.psm, var = "Gap100plus", which = "both")
bal.plot(model36.psm, var = "CETWI", which = "both")
bal.plot(model36.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model36.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model36.psm, type = "qq")

# Matched data
model36.matched <- match_data(model36.psm)

# Create trt_control variable
model36.matched <- model36.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Post-burn aerial seeding", "Post-burn control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn aerial seeding")))

# Center and scale numeric variables
model36.matched <- model36.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model36.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model36.matched,
  weights = weights
)

# G computation to estimate marginal effects
model36.pred <- avg_predictions(
  model = model36.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 36, .before = trt_control)
model36.pred

# Estimation of average treatment effect
model36.comp <- avg_comparisons(
  model = model36.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 36, .before = term)
model36.comp

# Plot
model36.plot <- model36.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "36. Colorado Plateaus: Post-burn aerial seeding",
       x = NULL)
model36.plot



## Arizona/New Mexico Plateau ---------------------------------------------

### 37. Herbicide ---------------------------------------------------------

# Filter data
model37 <- ldc.011 |>
  filter(Trt_Type_Sub == "Herbicide" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Arizona/New Mexico Plateau")

# PSM
model37.psm <- matchit(
  data = model37,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model37.psm
summary(model37.psm) # 49 treated matched

# Diagnostic love plot
model37.loveplot <- love.plot(model37.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "37. AZ/NM Plateau: Herbicide") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model37.loveplot

# eCDF plots
bal.plot(model37.psm, which = "both", type = "ecdf")
bal.plot(model37.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model37.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model37.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model37.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model37.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model37.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model37.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model37.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model37.psm, which = "both")
bal.plot(model37.psm, var = "BareSoil_FH", which = "both")
bal.plot(model37.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model37.psm, var = "ForbCover_AH", which = "both")
bal.plot(model37.psm, var = "GramCover_AH", which = "both")
bal.plot(model37.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model37.psm, var = "Gap100plus", which = "both")
bal.plot(model37.psm, var = "CETWI", which = "both")
bal.plot(model37.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model37.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model37.psm, type = "qq")

# Matched data
model37.matched <- match_data(model37.psm)

# Create trt_control variable
model37.matched <- model37.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Herbicide", "Control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Control", "Herbicide")))

# Center and scale numeric variables
model37.matched <- model37.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model37.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model37.matched,
  weights = weights
)

# G computation to estimate marginal effects
model37.pred <- avg_predictions(
  model = model37.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 37, .before = trt_control)
model37.pred

# Estimation of average treatment effect
model37.comp <- avg_comparisons(
  model = model37.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 37, .before = term)
model37.comp

# Plot
model37.plot <- model37.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(x = NULL,
       title = "37. AZ/NM Plateau: Herbicide")
model37.plot


### 38. Prescribed burn ---------------------------------------------------

# Filter data
model38.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Prescribed Burn" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Arizona/New Mexico Plateau")

# PSM
model38.psm <- matchit(
  data = model38.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model38.psm
summary(model38.psm) # 28 treated matched

# Diagnostic love plot
model38.loveplot <- love.plot(model38.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "38. AZ/NM Plateau: Prescribed burn") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model38.loveplot

# eCDF plots
bal.plot(model38.psm, which = "both", type = "ecdf")
bal.plot(model38.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model38.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model38.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model38.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model38.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model38.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model38.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model38.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model38.psm, which = "both")
bal.plot(model38.psm, var = "BareSoil_FH", which = "both")
bal.plot(model38.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model38.psm, var = "ForbCover_AH", which = "both")
bal.plot(model38.psm, var = "GramCover_AH", which = "both")
bal.plot(model38.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model38.psm, var = "Gap100plus", which = "both")
bal.plot(model38.psm, var = "CETWI", which = "both")
bal.plot(model38.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model38.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model38.psm, type = "qq")

# Matched data
model38.matched <- match_data(model38.psm)

# Create trt_control variable
model38.matched <- model38.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Prescribed burn", "Control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Control", "Prescribed burn")))

# Center and scale numeric variables
model38.matched <- model38.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model38.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model38.matched,
  weights = weights
)

# G computation to estimate marginal effects
model38.pred <- avg_predictions(
  model = model38.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 38, .before = trt_control)
model38.pred

# Estimation of average treatment effect
model38.comp <- avg_comparisons(
  model = model38.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 38, .before = term)
model38.comp

# Plot
model38.plot <- model38.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "38. AZ/NM Plateau: Prescribed burn",
       x = NULL)
model38.plot


### 39. Soil Disturbance --------------------------------------------------

# Filter data
model39.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Soil Disturbance" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Arizona/New Mexico Plateau")

# PSM
model39.psm <- matchit(
  data = model39.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model39.psm
summary(model39.psm) # 41 treated matched

# Diagnostic love plot
model39.loveplot <- love.plot(model39.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "39. AZ/NM Plateau: Soil disturbance") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model39.loveplot

# eCDF plots
bal.plot(model39.psm, which = "both", type = "ecdf")
bal.plot(model39.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model39.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model39.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model39.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model39.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model39.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model39.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model39.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model39.psm, which = "both")
bal.plot(model39.psm, var = "BareSoil_FH", which = "both")
bal.plot(model39.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model39.psm, var = "ForbCover_AH", which = "both")
bal.plot(model39.psm, var = "GramCover_AH", which = "both")
bal.plot(model39.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model39.psm, var = "Gap100plus", which = "both")
bal.plot(model39.psm, var = "CETWI", which = "both")
bal.plot(model39.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model39.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model39.psm, type = "qq")

# Matched data
model39.matched <- match_data(model39.psm)

# Create trt_control variable
model39.matched <- model39.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Soil disturbance", "Control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Control", "Soil disturbance")))

# Center and scale numeric variables
model39.matched <- model39.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model39.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model39.matched,
  weights = weights
)

# G computation to estimate marginal effects
model39.pred <- avg_predictions(
  model = model39.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 39, .before = trt_control)
model39.pred

# Estimation of average treatment effect
model39.comp <- avg_comparisons(
  model = model39.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 39, .before = term)
model39.comp # p = 0.009

# Plot
model39.plot <- model39.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "39. AZ/NM Plateau: Soil disturbance",
       x = NULL) +
  geom_signif(
    comparisons = list(c("Control", "Soil disturbance")),
    annotations = c("**")
  )
model39.plot




# Warm Deserts ------------------------------------------------------------

## Mojave Basin and Range -------------------------------------------------

### 40. Post-burn Aerial Seeding ------------------------------------------

# Filter data
model40.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Aerial Seeding" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post")) |>
  filter(EcoLvl3 == "Mojave Basin and Range")

# PSM
model40.psm <- matchit(
  data = model40.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model40.psm
summary(model40.psm) # 65 treated matched

# Diagnostic love plot
model40.loveplot <- love.plot(model40.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "40. Mojave Basin and Range: Post-burn aerial seeding") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model40.loveplot

# eCDF plots
bal.plot(model40.psm, which = "both", type = "ecdf")
bal.plot(model40.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model40.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model40.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model40.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model40.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model40.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model40.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model40.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model40.psm, which = "both")
bal.plot(model40.psm, var = "BareSoil_FH", which = "both")
bal.plot(model40.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model40.psm, var = "ForbCover_AH", which = "both")
bal.plot(model40.psm, var = "GramCover_AH", which = "both")
bal.plot(model40.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model40.psm, var = "Gap100plus", which = "both")
bal.plot(model40.psm, var = "CETWI", which = "both")
bal.plot(model40.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model40.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model40.psm, type = "qq")

# Matched data
model40.matched <- match_data(model40.psm)

# Create trt_control variable
model40.matched <- model40.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Post-burn aerial seeding", "Post-burn control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn aerial seeding")))

# Center and scale numeric variables
model40.matched <- model40.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model40.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model40.matched,
  weights = weights
)

# G computation to estimate marginal effects
model40.pred <- avg_predictions(
  model = model40.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 40, .before = trt_control)
model40.pred

# Estimation of average treatment effect
model40.comp <- avg_comparisons(
  model = model40.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 40, .before = term)
model40.comp # p = 0.009

# Plot
model40.plot <- model40.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "40. Mojave Basin and Range: Post-burn aerial seeding",
       x = NULL) +
  geom_signif(
    comparisons = list(c("Post-burn control", "Post-burn aerial seeding")),
    annotations = c("**")
  )
model40.plot



## Chihuahuan Desert ------------------------------------------------------

### 41. Herbicide ---------------------------------------------------------

# Filter data
model41.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Herbicide" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Chihuahuan Desert")

# PSM
model41.psm <- matchit(
  data = model41.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model41.psm
summary(model41.psm) # 64 treated matched

# Diagnostic love plot
model41.loveplot <- love.plot(model41.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "41. Chihuahuan Desert: Herbicide") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model41.loveplot

# eCDF plots
bal.plot(model41.psm, which = "both", type = "ecdf")
bal.plot(model41.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model41.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model41.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model41.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model41.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model41.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model41.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model41.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model41.psm, which = "both")
bal.plot(model41.psm, var = "BareSoil_FH", which = "both")
bal.plot(model41.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model41.psm, var = "ForbCover_AH", which = "both")
bal.plot(model41.psm, var = "GramCover_AH", which = "both")
bal.plot(model41.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model41.psm, var = "Gap100plus", which = "both")
bal.plot(model41.psm, var = "CETWI", which = "both")
bal.plot(model41.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model41.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model41.psm, type = "qq")

# Matched data
model41.matched <- match_data(model41.psm)

# Create trt_control variable
model41.matched <- model41.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Herbicide", "Control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Control", "Herbicide")))

# Center and scale numeric variables
model41.matched <- model41.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model41.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model41.matched,
  weights = weights
)

# G computation to estimate marginal effects
model41.pred <- avg_predictions(
  model = model41.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 41, .before = trt_control)
model41.pred # Herbicide NS (not different from 0)

# Estimation of average treatment effect
model41.comp <- avg_comparisons(
  model = model41.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 41, .before = term)
model41.comp

# Plot
model41.plot <- model41.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "41. Chihuahuan Desert: Herbicide",
       x = NULL)
model41.plot




# Temperate Sierras / Upper Gila ------------------------------------------

## Arizona/New Mexico Mountains -------------------------------------------

### 42. Prescribed Burn ---------------------------------------------------

# Filter data
model42.dat <- ldc.011 |>
  filter(Trt_Type_Sub == "Prescribed Burn" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Arizona/New Mexico Mountains")

# PSM
model42.psm <- matchit(
  data = model42.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model42.psm
summary(model42.psm) # 30 treated matched

# Diagnostic love plot
model42.loveplot <- love.plot(model42.psm, stars = "std",
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "42. AZ/NM Mountains: Prescribed burn") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model42.loveplot

# eCDF plots
bal.plot(model42.psm, which = "both", type = "ecdf")
bal.plot(model42.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model42.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model42.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model42.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model42.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model42.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model42.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model42.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model42.psm, which = "both")
bal.plot(model42.psm, var = "BareSoil_FH", which = "both")
bal.plot(model42.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model42.psm, var = "ForbCover_AH", which = "both")
bal.plot(model42.psm, var = "GramCover_AH", which = "both")
bal.plot(model42.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model42.psm, var = "Gap100plus", which = "both")
bal.plot(model42.psm, var = "CETWI", which = "both")
bal.plot(model42.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model42.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model42.psm, type = "qq")

# Matched data
model42.matched <- match_data(model42.psm)

# Create trt_control variable
model42.matched <- model42.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Prescribed burn", "Control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Control", "Prescribed burn")))

# Center and scale numeric variables
model42.matched <- model42.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model42.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model42.matched,
  weights = weights
)

# G computation to estimate marginal effects
model42.pred <- avg_predictions(
  model = model42.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 42, .before = trt_control) 
model42.pred

# Estimation of average treatment effect
model42.comp <- avg_comparisons(
  model = model42.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 42, .before = term) 
model42.comp

# Plot
model42.plot <- model42.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(x = NULL,
       title = "42. AZ/NM Mountains: Prescribed burn")
model42.plot




# Save matched data -------------------------------------------------------

# Matched data
save(list = ls(pattern = "\\.matched$"), 
     file = "RData/20.1_matched-data-2.RData")

# PSM
save(list = ls(pattern = "\\.psm$"), 
     file = "RData/20.1_PSM-2.RData")

# G computation
save(list = ls(pattern = "\\.pred$"), 
     file = "RData/20.1_g-computation-2.RData")

# Average treatment effect
save(list = ls(pattern = "\\.comp$"), 
     file = "RData/20.1_average-treatment-effect-2.RData")




# Write out figures -------------------------------------------------------

## Blue Mountains ---------------------------------------------------------

# 1. Blue Mountains: Herbicide
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model01_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model01.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model01_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model01.plot
dev.off()


# 2. Blue Mountains: Post-burn herbicide
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model02_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model02.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model02_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model02.plot
dev.off()



## Middle Rockies ---------------------------------------------------------

# 3. Middle Rockies: Herbicide
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model03_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model03.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model03_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model03.plot
dev.off()



## Southern Rockies -------------------------------------------------------

# 4. Southern Rockies: Herbicide
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model04_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model04.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model04_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model04.plot
dev.off()


# 5. Southern Rockies: Prescribed burn
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model05_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model05.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model05_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model05.plot
dev.off()



## Northwestern Great Plains ----------------------------------------------

# 6. NW Great Plains: Prescribed burn
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model06_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model06.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model06_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model06.plot
dev.off()



## Snake River Plain ------------------------------------------------------

# 7. Snake River Plain: Post-burn aerial seeding
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model07_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model07.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model07_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model07.plot
dev.off()


# 8. Snake River Plain: Post-burn aerial & drill seeding
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model08_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model08.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model08_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model08.plot
dev.off()


# 9. Snake River Plain: Post-burn closure
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model09_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model09.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model09_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model09.plot
dev.off()


# 10. Snake River Plain: Post-burn drill seeding
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model10_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model10.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model10_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model10.plot
dev.off()


# 11. Snake River Plain: Post-burn herbicide
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model11_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model11.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model11_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model11.plot
dev.off()



## Northern Basin and Range -----------------------------------------------

# 12. Northern BR: Drill seeding
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model12_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model12.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model12_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model12.plot
dev.off()


# 13. Northern BR: Drill seeding & soil disturbance
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model13_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model13.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model13_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model13.plot
dev.off()


# 14. Northern BR: Herbicide
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model14_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model14.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model14_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model14.plot
dev.off()


# 15. Northern BR: Prescribed burn
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model15_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model15.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model15_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model15.plot
dev.off()


# 16. Northern BR: Vegetation disturbance
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model16_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model16.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model16_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model16.plot
dev.off()


# 17. Northern BR: Post-burn aerial seeding
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model17_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model17.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model17_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model17.plot
dev.off()


# 18. Northern BR: Post-burn aerial and drill seeding
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model18_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model18.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model18_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model18.plot
dev.off()


# 19. Northern BR: Post-burn closure
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model19_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model19.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model19_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model19.plot
dev.off()


# 20. Northern BR: Post-burn drill seeding
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model20_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model20.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model20_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model20.plot
dev.off()


# 21. Northern BR: Post-burn herbicide
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model21_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model21.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model21_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model21.plot
dev.off()


# 22. Northern BR: Post-burn seedling planting
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model22_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model22.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model22_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model22.plot
dev.off()



## Central Basin and Range ------------------------------------------------

# 23. Central BR: Drill seeding & soil disturbance
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model23_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model23.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model23_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model23.plot
dev.off()


# 24. Central BR: Prescribed burn
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model24_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model24.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model24_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model24.plot
dev.off()


# 25. Central BR: Vegetation disturbance
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model25_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model25.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model25_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model25.plot
dev.off()


# 26. Central BR: Post-burn aerial seeding
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model26_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model26.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model26_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model26.plot
dev.off()


# 27. Central BR: Post-burn drill seeding
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model27_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model27.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model27_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model27.plot
dev.off()


# 28. Central BR: Post-burn ground seeding
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model28_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model28.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model28_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model28.plot
dev.off()


# 29. Central BR: Post-burn herbicide
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model29_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model29.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model29_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model29.plot
dev.off()



## Wyoming Basin ----------------------------------------------------------

# 30. Wyoming Basin: Prescribed burn
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model30_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model30.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model30_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model30.plot
dev.off()



## Colorado Plateaus ------------------------------------------------------

# 31. CO Plateaus: Aerial seeding & soil disturbance
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model31_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model31.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model31_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model31.plot
dev.off()


# 32. CO Plateaus: Herbicide
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model32_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model32.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model32_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model32.plot
dev.off()


# 33. CO Plateaus: Prescribed burn
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model33_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model33.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model33_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model33.plot
dev.off()


# 34. CO Plateaus: Soil disturbance
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model34_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model34.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model34_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model34.plot
dev.off()


# 35. CO Plateaus: Vegetation disturbance
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model35_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model35.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model35_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model35.plot
dev.off()


# 36. CO Plateaus: Post-burn aerial seeding
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model36_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model36.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model36_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model36.plot
dev.off()



## AZ/NM Plateau ----------------------------------------------------------

# 37. AZ/NM Plateau: Herbicide
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model37_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model37.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model37_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model37.plot
dev.off()


# 38. AZ/NM Plateau: Prescribed burn
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model38_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model38.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model38_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model38.plot
dev.off()


# 39. AZ/NM Plateau: Soil disturbance
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model39_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model39.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model39_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model39.plot
dev.off()



## Mojave Basin and Range -------------------------------------------------

# 40. Mojave BR: Post-burn aerial seeding
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model40_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model40.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model40_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model40.plot
dev.off()



## Chihuahuan Desert ------------------------------------------------------

# 41. Chihuahuan Desert: Herbicide
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model41_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model41.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model41_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model41.plot
dev.off()



## AZ/NM Mountains --------------------------------------------------------

# 42. AZ/NM Mountains: Prescribed burn
#   Love plot
tiff("figures/2026-06_PSM-and-permutation-tests/model42_loveplot.tiff",
     units = "in", width = 6, height = 4, res = 150)
model42.loveplot
dev.off()

#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests/model42_average-treatment-effect.tiff",
     units = "in", width = 6, height = 4, res = 150)
model42.plot
dev.off()



save.image("RData/15_propensity-score-matching.RData")

# Created: 2026-06-30
# Updated: 2026-06-30

# Purpose: Create summary tables of PSM results for the 42 models.

# Basically the same as 20.3.R script from Project v03.


library(tidyverse)
library(MatchIt)
library(marginaleffects)

# Load data ---------------------------------------------------------------

load("RData/15.1_matched-data.RData")
load("RData/15.1_PSM.RData")
load("RData/15.1_g-computation.RData")
load("RData/15.1_average-treatment-effect.RData")


# Count table -------------------------------------------------------------

# Create table of grouping cols
count.table0 <- mget(ls(pattern = "\\.matched$")) |>
  bind_rows(.id = "Model") |> 
  select(Model, EcoLvl3, trt_control) |> 
  distinct(.keep_all = TRUE) |> 
  filter(!trt_control %in% c("Control", "Post-burn control"))
count.table0$Model <- str_extract(count.table0$Model, "\\d+") |> as.integer()

# Add sample sizes from PSM
count.table.psm <- count.table0 |> 
  mutate(
    Control_n_before = c(
      summary(model01.psm)[[2]][2, 1], summary(model02.psm)[[2]][2, 1],
      summary(model03.psm)[[2]][2, 1], summary(model04.psm)[[2]][2, 1],
      summary(model05.psm)[[2]][2, 1], summary(model06.psm)[[2]][2, 1],
      summary(model07.psm)[[2]][2, 1], summary(model08.psm)[[2]][2, 1],
      summary(model09.psm)[[2]][2, 1], summary(model10.psm)[[2]][2, 1],
      summary(model11.psm)[[2]][2, 1], summary(model12.psm)[[2]][2, 1],
      summary(model13.psm)[[2]][2, 1], summary(model14.psm)[[2]][2, 1],
      summary(model15.psm)[[2]][2, 1], summary(model16.psm)[[2]][2, 1],
      summary(model17.psm)[[2]][2, 1], summary(model18.psm)[[2]][2, 1],
      summary(model19.psm)[[2]][2, 1], summary(model20.psm)[[2]][2, 1],
      summary(model21.psm)[[2]][2, 1], summary(model22.psm)[[2]][2, 1],
      summary(model23.psm)[[2]][2, 1], summary(model24.psm)[[2]][2, 1],
      summary(model25.psm)[[2]][2, 1], summary(model26.psm)[[2]][2, 1],
      summary(model27.psm)[[2]][2, 1], summary(model28.psm)[[2]][2, 1],
      summary(model29.psm)[[2]][2, 1], summary(model30.psm)[[2]][2, 1],
      summary(model31.psm)[[2]][2, 1], summary(model32.psm)[[2]][2, 1],
      summary(model33.psm)[[2]][2, 1], summary(model34.psm)[[2]][2, 1],
      summary(model35.psm)[[2]][2, 1], summary(model36.psm)[[2]][2, 1],
      summary(model37.psm)[[2]][2, 1], summary(model38.psm)[[2]][2, 1],
      summary(model39.psm)[[2]][2, 1], summary(model11.psm)[[2]][2, 1],
      summary(model41.psm)[[2]][2, 1], summary(model42.psm)[[2]][2, 1]
    )
  ) |> 
  mutate(
    Treated_n_before = c(
      summary(model01.psm)[[2]][2, 2], summary(model02.psm)[[2]][2, 2],
      summary(model03.psm)[[2]][2, 2], summary(model04.psm)[[2]][2, 2],
      summary(model05.psm)[[2]][2, 2], summary(model06.psm)[[2]][2, 2],
      summary(model07.psm)[[2]][2, 2], summary(model08.psm)[[2]][2, 2],
      summary(model09.psm)[[2]][2, 2], summary(model10.psm)[[2]][2, 2],
      summary(model11.psm)[[2]][2, 2], summary(model12.psm)[[2]][2, 2],
      summary(model13.psm)[[2]][2, 2], summary(model14.psm)[[2]][2, 2],
      summary(model15.psm)[[2]][2, 2], summary(model16.psm)[[2]][2, 2],
      summary(model17.psm)[[2]][2, 2], summary(model18.psm)[[2]][2, 2],
      summary(model19.psm)[[2]][2, 2], summary(model20.psm)[[2]][2, 2],
      summary(model21.psm)[[2]][2, 2], summary(model22.psm)[[2]][2, 2],
      summary(model23.psm)[[2]][2, 2], summary(model24.psm)[[2]][2, 2],
      summary(model25.psm)[[2]][2, 2], summary(model26.psm)[[2]][2, 2],
      summary(model27.psm)[[2]][2, 2], summary(model28.psm)[[2]][2, 2],
      summary(model29.psm)[[2]][2, 2], summary(model30.psm)[[2]][2, 2],
      summary(model31.psm)[[2]][2, 2], summary(model32.psm)[[2]][2, 2],
      summary(model33.psm)[[2]][2, 2], summary(model34.psm)[[2]][2, 2],
      summary(model35.psm)[[2]][2, 2], summary(model36.psm)[[2]][2, 2],
      summary(model37.psm)[[2]][2, 2], summary(model38.psm)[[2]][2, 2],
      summary(model39.psm)[[2]][2, 2], summary(model11.psm)[[2]][2, 2],
      summary(model41.psm)[[2]][2, 2], summary(model42.psm)[[2]][2, 2]
    )
  ) |> 
  mutate(
    Control_n_after = c(
      summary(model01.psm)[[2]][4, 1], summary(model02.psm)[[2]][4, 1],
      summary(model03.psm)[[2]][4, 1], summary(model04.psm)[[2]][4, 1],
      summary(model05.psm)[[2]][4, 1], summary(model06.psm)[[2]][4, 1],
      summary(model07.psm)[[2]][4, 1], summary(model08.psm)[[2]][4, 1],
      summary(model09.psm)[[2]][4, 1], summary(model10.psm)[[2]][4, 1],
      summary(model11.psm)[[2]][4, 1], summary(model12.psm)[[2]][4, 1],
      summary(model13.psm)[[2]][4, 1], summary(model14.psm)[[2]][4, 1],
      summary(model15.psm)[[2]][4, 1], summary(model16.psm)[[2]][4, 1],
      summary(model17.psm)[[2]][4, 1], summary(model18.psm)[[2]][4, 1],
      summary(model19.psm)[[2]][4, 1], summary(model20.psm)[[2]][4, 1],
      summary(model21.psm)[[2]][4, 1], summary(model22.psm)[[2]][4, 1],
      summary(model23.psm)[[2]][4, 1], summary(model24.psm)[[2]][4, 1],
      summary(model25.psm)[[2]][4, 1], summary(model26.psm)[[2]][4, 1],
      summary(model27.psm)[[2]][4, 1], summary(model28.psm)[[2]][4, 1],
      summary(model29.psm)[[2]][4, 1], summary(model30.psm)[[2]][4, 1],
      summary(model31.psm)[[2]][4, 1], summary(model32.psm)[[2]][4, 1],
      summary(model33.psm)[[2]][4, 1], summary(model34.psm)[[2]][4, 1],
      summary(model35.psm)[[2]][4, 1], summary(model36.psm)[[2]][4, 1],
      summary(model37.psm)[[2]][4, 1], summary(model38.psm)[[2]][4, 1],
      summary(model39.psm)[[2]][4, 1], summary(model11.psm)[[2]][4, 1],
      summary(model41.psm)[[2]][4, 1], summary(model42.psm)[[2]][4, 1]
    )
  ) |> 
  mutate(
    Treated_n_after = c(
      summary(model01.psm)[[2]][4, 2], summary(model02.psm)[[2]][4, 2],
      summary(model03.psm)[[2]][4, 2], summary(model04.psm)[[2]][4, 2],
      summary(model05.psm)[[2]][4, 2], summary(model06.psm)[[2]][4, 2],
      summary(model07.psm)[[2]][4, 2], summary(model08.psm)[[2]][4, 2],
      summary(model09.psm)[[2]][4, 2], summary(model10.psm)[[2]][4, 2],
      summary(model11.psm)[[2]][4, 2], summary(model12.psm)[[2]][4, 2],
      summary(model13.psm)[[2]][4, 2], summary(model14.psm)[[2]][4, 2],
      summary(model15.psm)[[2]][4, 2], summary(model16.psm)[[2]][4, 2],
      summary(model17.psm)[[2]][4, 2], summary(model18.psm)[[2]][4, 2],
      summary(model19.psm)[[2]][4, 2], summary(model20.psm)[[2]][4, 2],
      summary(model21.psm)[[2]][4, 2], summary(model22.psm)[[2]][4, 2],
      summary(model23.psm)[[2]][4, 2], summary(model24.psm)[[2]][4, 2],
      summary(model25.psm)[[2]][4, 2], summary(model26.psm)[[2]][4, 2],
      summary(model27.psm)[[2]][4, 2], summary(model28.psm)[[2]][4, 2],
      summary(model29.psm)[[2]][4, 2], summary(model30.psm)[[2]][4, 2],
      summary(model31.psm)[[2]][4, 2], summary(model32.psm)[[2]][4, 2],
      summary(model33.psm)[[2]][4, 2], summary(model34.psm)[[2]][4, 2],
      summary(model35.psm)[[2]][4, 2], summary(model36.psm)[[2]][4, 2],
      summary(model37.psm)[[2]][4, 2], summary(model38.psm)[[2]][4, 2],
      summary(model39.psm)[[2]][4, 2], summary(model11.psm)[[2]][4, 2],
      summary(model41.psm)[[2]][4, 2], summary(model42.psm)[[2]][4, 2]
    )
  )



## G computation ----------------------------------------------------------

# Create table of grouping cols
gcomp.pred.0 <- mget(ls(pattern = "\\.matched$")) |>
  bind_rows(.id = "Model") |> 
  select(Model, EcoLvl3, trt_control) |> 
  distinct(.keep_all = TRUE)
gcomp.pred.0$Model <- str_extract(gcomp.pred.0$Model, "\\d+") |> as.integer()

# Bind rows with values
gcomp.pred <- bind_rows(mget(ls(pattern = "\\.pred$")))

# Add other cols
gcomp.pred <- gcomp.pred |> 
  left_join(gcomp.pred.0)
gcomp.pred <- gcomp.pred[, c(1, 11, 2:10)]

# Add significance stars
gcomp.pred <- gcomp.pred |> 
  mutate(sig = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE            ~ ""
  ), .before = s.value)

# Add back-transformation for estimate
gcomp.pred <- gcomp.pred |> 
  mutate(estimate_bt = exp(estimate),
         .after = estimate)



## Average treatment effect -----------------------------------------------

# Create table of grouping cols
avg.comp0 <- mget(ls(pattern = "\\.matched$")) |>
  bind_rows(.id = "Model") |> 
  select(Model, EcoLvl3) |> 
  distinct(.keep_all = TRUE)
avg.comp0$Model <- str_extract(avg.comp0$Model, "\\d+") |> as.integer()

# Bind rows with values
avg.comp <- bind_rows(mget(ls(pattern = "\\.comp$")))

# Add other cols
avg.comp <- avg.comp0 |> 
  left_join(avg.comp) |> 
  select(-term)

# Add significance stars
avg.comp <- avg.comp |> 
  mutate(sig = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE            ~ ""
  ), .before = s.value)



# Created: 2026-07-08
# Updated: 2026-07-14

# Purpose: Run permutation tests for functional group cover, Shannon diversity,
#   and species/groups of interest (defined in 19.R).


library(tidyverse)
library(ggsignif)
library(RColorBrewer)
library(ggpubr)
library(gridExtra)

# Load data ---------------------------------------------------------------

all.matched.diversity <- read_csv("data/versions-from-R/18_shannon-diversity_all-models_v012.csv")
all.matched.species.funct <- read_csv("data/versions-from-R/19_species-of-interest-and-functional-group-cover_all-models_v012.csv")


# Data wrangling ----------------------------------------------------------

# Combine diversity & species/function group dfs
all.matched <- all.matched.diversity |> 
  left_join(all.matched.species.funct)



# NW Forested Mts / Western Cordillera ------------------------------------

## Blue Mountains ---------------------------------------------------------

### 1. Herbicide ----------------------------------------------------------

# Filter for model
model01.matched <- all.matched |> 
  filter(Model == 1) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover, PJ_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Control", "Herbicide")))

#   pivot_longer() for cover & shannon cols
model01.matched <- model01.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model01.matched <- model01.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia",
             indicators == "PJ_cover" ~ "Pinus and Juniperus"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia",
                                        "Pinus and Juniperus")))

# Calculate observed mean difference
model01.diff <- model01.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model01.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model01.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model01.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values01 <- model01.perm |>
  inner_join(model01.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values01 # p = 0.002 for shannon

# Boxplot
model01.bp <- model01.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "1. Blue Mountains: Herbicide") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")),
               "Pinus and Juniperus" = expression(italic("Pinus") ~ "&" ~ italic("Juniperus")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model01.bp

# Plot frequency distribution
#   Annual forb
model01.annforb <- model01.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model01.diff$obs_diff[model01.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model01.annforb

#   Annual grass
model01.anngrass <- model01.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model01.diff$obs_diff[model01.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model01.anngrass

#   Perennial forb
model01.perforb <- model01.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model01.diff$obs_diff[model01.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model01.perforb

#   Perennial grass
model01.pergrass <- model01.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model01.diff$obs_diff[model01.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model01.pergrass

#   Shrub
model01.shrub <- model01.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model01.diff$obs_diff[model01.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model01.shrub

#   Shannon diversity
model01.shannon <- model01.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model01.diff$obs_diff[model01.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity (**)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model01.shannon

#   Bromus tectorum
model01.brte <- model01.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model01.diff$obs_diff[model01.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model01.brte

#   Artemisia
model01.artemisia <- model01.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model01.diff$obs_diff[model01.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model01.artemisia

#   Pinus and Juniperus
model01.pj <- model01.perm |> 
  filter(indicators == "Pinus and Juniperus") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model01.diff$obs_diff[model01.diff$indicators == "Pinus and Juniperus"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Pinus") ~ "&" ~ italic("Juniperus"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model01.pj

# Combine plots
grid.arrange(
  model01.bp, model01.annforb, model01.anngrass,
  model01.perforb, model01.pergrass, model01.shrub,
  model01.shannon, model01.brte, model01.artemisia, 
  model01.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)


### 2. Post-burn herbicide ------------------------------------------------

# Filter for model
model02.matched <- all.matched |> 
  filter(Model == 2) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover, PJ_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn herbicide")))

#   pivot_longer() for cover & shannon cols
model02.matched <- model02.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model02.matched <- model02.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia",
             indicators == "PJ_cover" ~ "Pinus and Juniperus"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia",
                                        "Pinus and Juniperus")))

# Calculate observed mean difference
model02.diff <- model02.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model02.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model02.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model02.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values02 <- model02.perm |>
  inner_join(model02.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values02

# Boxplot
model02.bp <- model02.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "2. Blue Mountains: Post-burn herbicide") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")),
               "Pinus and Juniperus" = expression(italic("Pinus") ~ "&" ~ italic("Juniperus")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model02.bp

# Plot frequency distribution
#   Annual forb
model02.annforb <- model02.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model02.diff$obs_diff[model02.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model02.annforb

#   Annual grass
model02.anngrass <- model02.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model02.diff$obs_diff[model02.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model02.anngrass

#   Perennial forb
model02.perforb <- model02.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model02.diff$obs_diff[model02.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model02.perforb

#   Perennial grass
model02.pergrass <- model02.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model02.diff$obs_diff[model02.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model02.pergrass

#   Shrub
model02.shrub <- model02.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model02.diff$obs_diff[model02.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model02.shrub

#   Shannon diversity
model02.shannon <- model02.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model02.diff$obs_diff[model02.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model02.shannon

#   Bromus tectorum
model02.brte <- model02.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model02.diff$obs_diff[model02.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model02.brte

#   Artemisia
model02.artemisia <- model02.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model02.diff$obs_diff[model02.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model02.artemisia

#   Pinus and Juniperus
model02.pj <- model02.perm |> 
  filter(indicators == "Pinus and Juniperus") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model02.diff$obs_diff[model02.diff$indicators == "Pinus and Juniperus"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Pinus") ~ "&" ~ italic("Juniperus"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model02.pj

# Combine plots
grid.arrange(
  model02.bp, model02.annforb, model02.anngrass,
  model02.perforb, model02.pergrass, model02.shrub,
  model02.shannon, model02.brte, model02.artemisia, 
  model02.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)



## Middle Rockies ---------------------------------------------------------

### 3. Herbicide ----------------------------------------------------------

# Filter for model
model03.matched <- all.matched |> 
  filter(Model == 3) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Control", "Herbicide")))

#   pivot_longer() for cover & shannon cols
model03.matched <- model03.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model03.matched <- model03.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia")))

# Calculate observed mean difference
model03.diff <- model03.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model03.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model03.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model03.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values03 <- model03.perm |>
  inner_join(model03.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values03

# Boxplot
model03.bp <- model03.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "3. Middle Rockies: Herbicide") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model03.bp

# Plot frequency distribution
#   Annual forb
model03.annforb <- model03.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model03.diff$obs_diff[model03.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model03.annforb

#   Annual grass
model03.anngrass <- model03.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model03.diff$obs_diff[model03.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model03.anngrass

#   Perennial forb
model03.perforb <- model03.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model03.diff$obs_diff[model03.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model03.perforb

#   Perennial grass
model03.pergrass <- model03.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model03.diff$obs_diff[model03.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model03.pergrass

#   Shrub
model03.shrub <- model03.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model03.diff$obs_diff[model03.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model03.shrub

#   Shannon diversity
model03.shannon <- model03.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model03.diff$obs_diff[model03.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model03.shannon

#   Bromus tectorum
model03.brte <- model03.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model03.diff$obs_diff[model03.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model03.brte

#   Artemisia
model03.artemisia <- model03.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model03.diff$obs_diff[model03.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model03.artemisia

# Combine plots
grid.arrange(
  model03.bp, model03.annforb, model03.anngrass,
  model03.perforb, model03.pergrass, model03.shrub,
  model03.shannon, model03.brte, model03.artemisia, 
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(NA, 8, 8, 9, 9, NA)
  )
)



## Southern Rockies -------------------------------------------------------

### 4. Herbicide ----------------------------------------------------------

# Filter for model
model04.matched <- all.matched |> 
  filter(Model == 4) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, Artemisia_cover, PJ_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Control", "Herbicide")))

#   pivot_longer() for cover & shannon cols
model04.matched <- model04.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model04.matched <- model04.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "Artemisia_cover" ~ "Artemisia",
             indicators == "PJ_cover" ~ "Pinus and Juniperus"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Artemisia",
                                        "Pinus and Juniperus")))

# Calculate observed mean difference
model04.diff <- model04.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model04.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model04.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model04.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values04 <- model04.perm |>
  inner_join(model04.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values04 # p = 0.014 for Artemisia

# Boxplot
model04.bp <- model04.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "4. Southern Rockies: Herbicide") +
   geom_signif(
    y_position = 65,
    xmin = 5.8,
    xmax = 6.2, 
    annotations = c("*")
  ) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Artemisia" = expression(italic("Artemisia")),
               "Pinus and Juniperus" = expression(italic("Pinus") ~ "&" ~ italic("Juniperus")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model04.bp

# Plot frequency distribution
#   Annual forb
model04.annforb <- model04.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model04.diff$obs_diff[model04.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model04.annforb

#   Annual grass
model04.anngrass <- model04.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model04.diff$obs_diff[model04.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model04.anngrass

#   Perennial forb
model04.perforb <- model04.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model04.diff$obs_diff[model04.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model04.perforb

#   Perennial grass
model04.pergrass <- model04.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model04.diff$obs_diff[model04.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model04.pergrass

#   Shrub
model04.shrub <- model04.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model04.diff$obs_diff[model04.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model04.shrub

#   Shannon diversity
model04.shannon <- model04.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model04.diff$obs_diff[model04.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model04.shannon

#   Artemisia
model04.artemisia <- model04.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model04.diff$obs_diff[model04.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia") ~ "(*)")) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model04.artemisia

#   Pinus and Juniperus
model04.pj <- model04.perm |> 
  filter(indicators == "Pinus and Juniperus") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model04.diff$obs_diff[model04.diff$indicators == "Pinus and Juniperus"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Pinus") ~ "&" ~ italic("Juniperus"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model04.pj

# Combine plots
grid.arrange(
  model04.bp, model04.annforb, model04.anngrass,
  model04.perforb, model04.pergrass, model04.shrub,
  model04.shannon, model04.artemisia, model04.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(NA, 8, 8, 9, 9, NA)
  )
)


### 5. Prescribed burn ----------------------------------------------------

# Filter for model
model05.matched <- all.matched |> 
  filter(Model == 5) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, Artemisia_cover, PJ_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Control", "Prescribed burn")))

#   pivot_longer() for cover & shannon cols
model05.matched <- model05.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model05.matched <- model05.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "Artemisia_cover" ~ "Artemisia",
             indicators == "PJ_cover" ~ "Pinus and Juniperus"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Artemisia",
                                        "Pinus and Juniperus")))

# Calculate observed mean difference
model05.diff <- model05.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model05.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model05.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model05.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values05 <- model05.perm |>
  inner_join(model05.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values05

# Boxplot
model05.bp <- model05.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "5. Southern Rockies: Prescribed burn") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Artemisia" = expression(italic("Artemisia")),
               "Pinus and Juniperus" = expression(italic("Pinus") ~ "&" ~ italic("Juniperus")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model05.bp

# Plot frequency distribution
#   Annual forb
model05.annforb <- model05.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model05.diff$obs_diff[model05.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model05.annforb

#   Annual grass
model05.anngrass <- model05.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model05.diff$obs_diff[model05.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model05.anngrass

#   Perennial forb
model05.perforb <- model05.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model05.diff$obs_diff[model05.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model05.perforb

#   Perennial grass
model05.pergrass <- model05.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model05.diff$obs_diff[model05.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model05.pergrass

#   Shrub
model05.shrub <- model05.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model05.diff$obs_diff[model05.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model05.shrub

#   Shannon diversity
model05.shannon <- model05.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model05.diff$obs_diff[model05.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model05.shannon

#   Artemisia
model05.artemisia <- model05.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model05.diff$obs_diff[model05.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model05.artemisia

#   Pinus and Juniperus
model05.pj <- model05.perm |> 
  filter(indicators == "Pinus and Juniperus") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model05.diff$obs_diff[model05.diff$indicators == "Pinus and Juniperus"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Pinus") ~ "&" ~ italic("Juniperus"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model05.pj

# Combine plots
grid.arrange(
  model05.bp, model05.annforb, model05.anngrass,
  model05.perforb, model05.pergrass, model05.shrub,
  model05.shannon, model05.artemisia, model05.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(NA, 8, 8, 9, 9, NA)
  )
)



# Great Plains / West-Central Semiarid Prairies ---------------------------

## Northwestern Great Plains ----------------------------------------------

### 6. Prescribed burn ----------------------------------------------------

# Filter for model
model06.matched <- all.matched |> 
  filter(Model == 6) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover, 
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Control", "Prescribed burn")))

#   pivot_longer() for cover & shannon cols
model06.matched <- model06.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model06.matched <- model06.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia")))

# Calculate observed mean difference
model06.diff <- model06.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model06.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model06.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model06.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values06 <- model06.perm |>
  inner_join(model06.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values06 # p = 0.02 for shannon; p = 0.004 for Artemisia

# Boxplot
model06.bp <- model06.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "6. Northwestern Great Plains: Prescribed burn") +
  geom_signif(
    y_position = 40,
    xmin = 6.8,
    xmax = 7.2, 
    annotations = c("**")
  ) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model06.bp

# Plot frequency distribution
#   Annual forb
model06.annforb <- model06.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model06.diff$obs_diff[model06.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model06.annforb

#   Annual grass
model06.anngrass <- model06.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model06.diff$obs_diff[model06.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model06.anngrass

#   Perennial forb
model06.perforb <- model06.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model06.diff$obs_diff[model06.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model06.perforb

#   Perennial grass
model06.pergrass <- model06.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model06.diff$obs_diff[model06.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model06.pergrass

#   Shrub
model06.shrub <- model06.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model06.diff$obs_diff[model06.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model06.shrub

#   Shannon diversity
model06.shannon <- model06.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model06.diff$obs_diff[model06.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity (*)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model06.shannon

#   Bromus tectorum
model06.brte <- model06.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model06.diff$obs_diff[model06.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model06.brte

#   Artemisia
model06.artemisia <- model06.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model06.diff$obs_diff[model06.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia") ~ "(**)")) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model06.artemisia

# Combine plots
grid.arrange(
  model06.bp, model06.annforb, model06.anngrass,
  model06.perforb, model06.pergrass, model06.shrub,
  model06.shannon, model06.brte, model06.artemisia, 
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(NA, 8, 8, 9, 9, NA)
  )
)




# Cold Deserts ------------------------------------------------------------

## Snake River Plain ------------------------------------------------------

### 7. Post-burn aerial seeding -------------------------------------------

# Filter for model
model07.matched <- all.matched |> 
  filter(Model == 7) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn aerial seeding")))

#   pivot_longer() for cover & shannon cols
model07.matched <- model07.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model07.matched <- model07.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia")))

# Calculate observed mean difference
model07.diff <- model07.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model07.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model07.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model07.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values07 <- model07.perm |>
  inner_join(model07.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values07 # p = 0.013 for perennial forb; p = 0.005 for perennial grass
#             p < 0.001 for shannon

# Boxplot
model07.bp <- model07.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "7. Snake River Plain: Post-burn aerial seeding") +
  geom_signif(
    y_position = 55,
    xmin = 2.8,
    xmax = 3.2, 
    annotations = c("*")
  ) +
  geom_signif(
    y_position = 135,
    xmin = 3.8,
    xmax = 4.2, 
    annotations = c("**")
  ) +
  ylim(0, 145) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model07.bp

# Plot frequency distribution
#   Annual forb
model07.annforb <- model07.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model07.diff$obs_diff[model07.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model07.annforb

#   Annual grass
model07.anngrass <- model07.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model07.diff$obs_diff[model07.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model07.anngrass

#   Perennial forb
model07.perforb <- model07.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model07.diff$obs_diff[model07.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb (*)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model07.perforb

#   Perennial grass
model07.pergrass <- model07.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model07.diff$obs_diff[model07.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass (**)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model07.pergrass

#   Shrub
model07.shrub <- model07.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model07.diff$obs_diff[model07.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model07.shrub

#   Shannon diversity
model07.shannon <- model07.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model07.diff$obs_diff[model07.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity (***)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model07.shannon

#   Bromus tectorum
model07.brte <- model07.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model07.diff$obs_diff[model07.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model07.brte

#   Artemisia
model07.artemisia <- model07.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model07.diff$obs_diff[model07.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model07.artemisia

# Combine plots
grid.arrange(
  model07.bp, model07.annforb, model07.anngrass,
  model07.perforb, model07.pergrass, model07.shrub,
  model07.shannon, model07.brte, model07.artemisia,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(NA, 8, 8, 9, 9, NA)
  )
)


### 8. Post-burn aerial & drill seeding -----------------------------------

# Filter for model
model08.matched <- all.matched |> 
  filter(Model == 8) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", 
                                                      "Post-burn aerial & drill seeding")))

#   pivot_longer() for cover & shannon cols
model08.matched <- model08.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model08.matched <- model08.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia")))

# Calculate observed mean difference
model08.diff <- model08.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model08.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model08.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model08.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values08 <- model08.perm |>
  inner_join(model08.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values08 # p = 0.04 for BRTE

# Boxplot
model08.bp <- model08.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "8. Snake River Plain: Post-burn aerial & drill seeding") +
  geom_signif(
    y_position = 100,
    xmin = 5.8,
    xmax = 6.2, 
    annotations = c("*")
  ) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model08.bp

# Plot frequency distribution
#   Annual forb
model08.annforb <- model08.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model08.diff$obs_diff[model08.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model08.annforb

#   Annual grass
model08.anngrass <- model08.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model08.diff$obs_diff[model08.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model08.anngrass

#   Perennial forb
model08.perforb <- model08.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model08.diff$obs_diff[model08.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model08.perforb

#   Perennial grass
model08.pergrass <- model08.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model08.diff$obs_diff[model08.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model08.pergrass

#   Shrub
model08.shrub <- model08.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model08.diff$obs_diff[model08.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model08.shrub

#   Shannon diversity
model08.shannon <- model08.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model08.diff$obs_diff[model08.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model08.shannon

#   Bromus tectorum
model08.brte <- model08.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model08.diff$obs_diff[model08.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum") ~ "(*)")) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model08.brte

#   Artemisia
model08.artemisia <- model08.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model08.diff$obs_diff[model08.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model08.artemisia

# Combine plots
grid.arrange(
  model08.bp, model08.annforb, model08.anngrass,
  model08.perforb, model08.pergrass, model08.shrub,
  model08.shannon, model08.brte, model08.artemisia,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(NA, 8, 8, 9, 9, NA)
  )
)


### 9. Post-burn closure --------------------------------------------------

# Filter for model
model09.matched <- all.matched |> 
  filter(Model == 9) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn closure")))

#   pivot_longer() for cover & shannon cols
model09.matched <- model09.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model09.matched <- model09.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia")))

# Calculate observed mean difference
model09.diff <- model09.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model09.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model09.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model09.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values09 <- model09.perm |>
  inner_join(model09.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values09 # p = 0.003 for annual grass; p = 0.012 for perennial grass;
#             p < 0.001 for shannon; p = 0.015 for BRTE

# Boxplot
model09.bp <- model09.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "9. Snake River Plain: Post-burn closure") +
  geom_signif(
    y_position = 140,
    xmin = 1.8,
    xmax = 2.2, 
    annotations = c("**")
  ) +
  geom_signif(
    y_position = 105,
    xmin = 3.8,
    xmax = 4.2, 
    annotations = c("*")
  ) +
  geom_signif(
    y_position = 105,
    xmin = 5.8,
    xmax = 6.2, 
    annotations = c("*")
  ) +
  ylim(0, 150) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model09.bp

# Plot frequency distribution
#   Annual forb
model09.annforb <- model09.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model09.diff$obs_diff[model09.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model09.annforb

#   Annual grass
model09.anngrass <- model09.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model09.diff$obs_diff[model09.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass (**)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model09.anngrass

#   Perennial forb
model09.perforb <- model09.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model09.diff$obs_diff[model09.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model09.perforb

#   Perennial grass
model09.pergrass <- model09.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model09.diff$obs_diff[model09.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass (*)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model09.pergrass

#   Shrub
model09.shrub <- model09.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model09.diff$obs_diff[model09.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model09.shrub

#   Shannon diversity
model09.shannon <- model09.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model09.diff$obs_diff[model09.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity (***)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model09.shannon

#   Bromus tectorum
model09.brte <- model09.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model09.diff$obs_diff[model09.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum") ~ "(*)")) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model09.brte

#   Artemisia
model09.artemisia <- model09.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model09.diff$obs_diff[model09.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model09.artemisia

# Combine plots
grid.arrange(
  model09.bp, model09.annforb, model09.anngrass,
  model09.perforb, model09.pergrass, model09.shrub,
  model09.shannon, model09.brte, model09.artemisia,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(NA, 8, 8, 9, 9, NA)
  )
)


### 10. Post-burn drill seeding -------------------------------------------------

# Filter for model
model10.matched <- all.matched |> 
  filter(Model == 10) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn drill seeding")))

#   pivot_longer() for cover & shannon cols
model10.matched <- model10.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model10.matched <- model10.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia")))

# Calculate observed mean difference
model10.diff <- model10.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model10.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model10.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model10.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values10 <- model10.perm |>
  inner_join(model10.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values10 # p = 0.04 for perennial grass

# Boxplot
model10.bp <- model10.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "10. Snake River Plain: Post-burn drill seeding") +
  geom_signif(
    y_position = 100,
    xmin = 3.8,
    xmax = 4.2, 
    annotations = c("*")
  ) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model10.bp

# Plot frequency distribution
#   Annual forb
model10.annforb <- model10.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model10.diff$obs_diff[model10.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model10.annforb

#   Annual grass
model10.anngrass <- model10.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model10.diff$obs_diff[model10.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model10.anngrass

#   Perennial forb
model10.perforb <- model10.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model10.diff$obs_diff[model10.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model10.perforb

#   Perennial grass
model10.pergrass <- model10.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model10.diff$obs_diff[model10.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass (*)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model10.pergrass

#   Shrub
model10.shrub <- model10.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model10.diff$obs_diff[model10.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model10.shrub

#   Shannon diversity
model10.shannon <- model10.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model10.diff$obs_diff[model10.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model10.shannon

#   Bromus tectorum
model10.brte <- model10.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model10.diff$obs_diff[model10.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model10.brte

#   Artemisia
model10.artemisia <- model10.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model10.diff$obs_diff[model10.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model10.artemisia

# Combine plots
grid.arrange(
  model10.bp, model10.annforb, model10.anngrass,
  model10.perforb, model10.pergrass, model10.shrub,
  model10.shannon, model10.brte, model10.artemisia,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(NA, 8, 8, 9, 9, NA)
  )
)


### 11. Post-burn herbicide -----------------------------------------------

# Filter for model
model11.matched <- all.matched |> 
  filter(Model == 11) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn herbicide")))

#   pivot_longer() for cover & shannon cols
model11.matched <- model11.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model11.matched <- model11.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia")))

# Calculate observed mean difference
model11.diff <- model11.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model11.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model11.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model11.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values11 <- model11.perm |>
  inner_join(model11.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values11

# Boxplot
model11.bp <- model11.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "11. Snake River Plain: Post-burn herbicide") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model11.bp

# Plot frequency distribution
#   Annual forb
model11.annforb <- model11.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model11.diff$obs_diff[model11.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model11.annforb

#   Annual grass
model11.anngrass <- model11.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model11.diff$obs_diff[model11.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model11.anngrass

#   Perennial forb
model11.perforb <- model11.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model11.diff$obs_diff[model11.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model11.perforb

#   Perennial grass
model11.pergrass <- model11.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model11.diff$obs_diff[model11.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model11.pergrass

#   Shrub
model11.shrub <- model11.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model11.diff$obs_diff[model11.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model11.shrub

#   Shannon diversity
model11.shannon <- model11.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model11.diff$obs_diff[model11.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model11.shannon

#   Bromus tectorum
model11.brte <- model11.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model11.diff$obs_diff[model11.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model11.brte

#   Artemisia
model11.artemisia <- model11.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model11.diff$obs_diff[model11.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model11.artemisia

# Combine plots
grid.arrange(
  model11.bp, model11.annforb, model11.anngrass,
  model11.perforb, model11.pergrass, model11.shrub,
  model11.shannon, model11.brte, model11.artemisia,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(NA, 8, 8, 9, 9, NA)
  )
)



## Northern Basin and Range -----------------------------------------------

### 12. Drill seeding -----------------------------------------------------

# Filter for model
model12.matched <- all.matched |> 
  filter(Model == 12) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover, PJ_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Control", "Drill seeding")))

#   pivot_longer() for cover & shannon cols
model12.matched <- model12.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model12.matched <- model12.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia",
             indicators == "PJ_cover" ~ "Pinus and Juniperus"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia",
                                        "Pinus and Juniperus")))

# Calculate observed mean difference
model12.diff <- model12.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model12.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model12.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model12.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values12 <- model12.perm |>
  inner_join(model12.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values12 # p = 0.007 for annual grass; p < 0.001 for perennial grass, 
#     p = 0.004 for BRTE

# Boxplot
model12.bp <- model12.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "12. Northern Basin and Range: Drill seeding") +
  geom_signif(
    y_position = 120,
    xmin = 1.8,
    xmax = 2.2, 
    annotations = c("**")
  ) +
  geom_signif(
    y_position = 110,
    xmin = 3.8,
    xmax = 4.2, 
    annotations = c("***")
  ) +
  geom_signif(
    y_position = 95,
    xmin = 5.8,
    xmax = 6.2, 
    annotations = c("**")
  ) +
  ylim(0, 130) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")),
               "Pinus and Juniperus" = expression(italic("Pinus") ~ "&" ~ italic("Juniperus")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model12.bp

# Plot frequency distribution
#   Annual forb
model12.annforb <- model12.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model12.diff$obs_diff[model12.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model12.annforb

#   Annual grass
model12.anngrass <- model12.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model12.diff$obs_diff[model12.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass (**)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model12.anngrass

#   Perennial forb
model12.perforb <- model12.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model12.diff$obs_diff[model12.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model12.perforb

#   Perennial grass
model12.pergrass <- model12.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model12.diff$obs_diff[model12.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass (***)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model12.pergrass

#   Shrub
model12.shrub <- model12.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model12.diff$obs_diff[model12.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model12.shrub

#   Shannon diversity
model12.shannon <- model12.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model12.diff$obs_diff[model12.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model12.shannon

#   Bromus tectorum
model12.brte <- model12.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model12.diff$obs_diff[model12.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum") ~ "(**)")) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model12.brte

#   Artemisia
model12.artemisia <- model12.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model12.diff$obs_diff[model12.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model12.artemisia

#   Pinus and Juniperus
model12.pj <- model12.perm |> 
  filter(indicators == "Pinus and Juniperus") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model12.diff$obs_diff[model12.diff$indicators == "Pinus and Juniperus"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Pinus") ~ "&" ~ italic("Juniperus"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model12.pj

# Combine plots
grid.arrange(
  model12.bp, model12.annforb, model12.anngrass,
  model12.perforb, model12.pergrass, model12.shrub,
  model12.shannon, model12.brte, model12.artemisia, 
  model12.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)


### 13. Drill seeding & soil disturbance ----------------------------------

# Filter for model
model13.matched <- all.matched |> 
  filter(Model == 13) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover, PJ_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Control", "Drill seeding & soil disturbance")))

#   pivot_longer() for cover & shannon cols
model13.matched <- model13.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model13.matched <- model13.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia",
             indicators == "PJ_cover" ~ "Pinus and Juniperus"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia",
                                        "Pinus and Juniperus")))

# Calculate observed mean difference
model13.diff <- model13.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model13.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model13.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model13.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values13 <- model13.perm |>
  inner_join(model13.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values13 # p < 0.001 for perennial forb; p = 0.03 for PJ

# Boxplot
model13.bp <- model13.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "13. Northern Basin and Range: Drill seeding & soil disturbance") +
  geom_signif(
    y_position = 40,
    xmin = 2.8,
    xmax = 3.2, 
    annotations = c("***")
  ) +
  geom_signif(
    y_position = 30,
    xmin = 7.8,
    xmax = 8.2, 
    annotations = c("*")
  ) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")),
               "Pinus and Juniperus" = expression(italic("Pinus") ~ "&" ~ italic("Juniperus")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model13.bp

# Plot frequency distribution
#   Annual forb
model13.annforb <- model13.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model13.diff$obs_diff[model13.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model13.annforb

#   Annual grass
model13.anngrass <- model13.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model13.diff$obs_diff[model13.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model13.anngrass

#   Perennial forb
model13.perforb <- model13.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model13.diff$obs_diff[model13.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb (***)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model13.perforb

#   Perennial grass
model13.pergrass <- model13.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model13.diff$obs_diff[model13.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model13.pergrass

#   Shrub
model13.shrub <- model13.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model13.diff$obs_diff[model13.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model13.shrub

#   Shannon diversity
model13.shannon <- model13.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model13.diff$obs_diff[model13.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model13.shannon

#   Bromus tectorum
model13.brte <- model13.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model13.diff$obs_diff[model13.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model13.brte

#   Artemisia
model13.artemisia <- model13.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model13.diff$obs_diff[model13.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model13.artemisia

#   Pinus and Juniperus
model13.pj <- model13.perm |> 
  filter(indicators == "Pinus and Juniperus") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model13.diff$obs_diff[model13.diff$indicators == "Pinus and Juniperus"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Pinus") ~ "&" ~ italic("Juniperus") ~ "(*)")) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model13.pj

# Combine plots
grid.arrange(
  model13.bp, model13.annforb, model13.anngrass,
  model13.perforb, model13.pergrass, model13.shrub,
  model13.shannon, model13.brte, model13.artemisia, 
  model13.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)


### 14. Herbicide ---------------------------------------------------------

# Filter for model
model14.matched <- all.matched |> 
  filter(Model == 14) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover, PJ_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Control", "Herbicide")))

#   pivot_longer() for cover & shannon cols
model14.matched <- model14.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model14.matched <- model14.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia",
             indicators == "PJ_cover" ~ "Pinus and Juniperus"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia",
                                        "Pinus and Juniperus")))
# Calculate observed mean difference
model14.diff <- model14.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model14.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model14.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model14.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values14 <- model14.perm |>
  inner_join(model14.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values14 # p = 0.04 for perennial forb; p = 0.03 for PJ

# Boxplot
model14.bp <- model14.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "14. Northern Basin and Range: Herbicide") +
  geom_signif(
    y_position = 80,
    xmin = 2.8,
    xmax = 3.2, 
    annotations = c("*")
  ) +
  geom_signif(
    y_position = 60,
    xmin = 7.8,
    xmax = 8.2, 
    annotations = c("*")
  ) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")),
               "Pinus and Juniperus" = expression(italic("Pinus") ~ "&" ~ italic("Juniperus")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model14.bp

# Plot frequency distribution
#   Annual forb
model14.annforb <- model14.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model14.diff$obs_diff[model14.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model14.annforb

#   Annual grass
model14.anngrass <- model14.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model14.diff$obs_diff[model14.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model14.anngrass

#   Perennial forb
model14.perforb <- model14.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model14.diff$obs_diff[model14.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb (*)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model14.perforb

#   Perennial grass
model14.pergrass <- model14.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model14.diff$obs_diff[model14.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model14.pergrass

#   Shrub
model14.shrub <- model14.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model14.diff$obs_diff[model14.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model14.shrub

#   Shannon diversity
model14.shannon <- model14.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model14.diff$obs_diff[model14.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model14.shannon

#   Bromus tectorum
model14.brte <- model14.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model14.diff$obs_diff[model14.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model14.brte

#   Artemisia
model14.artemisia <- model14.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model14.diff$obs_diff[model14.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model14.artemisia

#   Pinus and Juniperus
model14.pj <- model14.perm |> 
  filter(indicators == "Pinus and Juniperus") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model14.diff$obs_diff[model14.diff$indicators == "Pinus and Juniperus"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Pinus") ~ "&" ~ italic("Juniperus") ~ "(*)")) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model14.pj

# Combine plots
grid.arrange(
  model14.bp, model14.annforb, model14.anngrass,
  model14.perforb, model14.pergrass, model14.shrub,
  model14.shannon, model14.brte, model14.artemisia, 
  model14.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)


### 15. Prescribed burn ---------------------------------------------------

# Filter for model
model15.matched <- all.matched |> 
  filter(Model == 15) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover, PJ_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Control", "Prescribed burn")))

#   pivot_longer() for cover & shannon cols
model15.matched <- model15.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model15.matched <- model15.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia",
             indicators == "PJ_cover" ~ "Pinus and Juniperus"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia",
                                        "Pinus and Juniperus")))

# Calculate observed mean difference
model15.diff <- model15.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model15.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model15.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model15.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values15 <- model15.perm |>
  inner_join(model15.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values15 # p = 0.0011 for perennial grass; p = 0.04 for shannon;
# p = 0.02 for Artemisia

# Boxplot
model15.bp <- model15.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "15. Northern Basin and Range: Prescribed burn") +
  geom_signif(
    y_position = 107,
    xmin = 3.8,
    xmax = 4.2, 
    annotations = c("**")
  ) +
  geom_signif(
    y_position = 60,
    xmin = 6.8,
    xmax = 7.2, 
    annotations = c("*")
  ) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")),
               "Pinus and Juniperus" = expression(italic("Pinus") ~ "&" ~ italic("Juniperus")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model15.bp

# Plot frequency distribution
#   Annual forb
model15.annforb <- model15.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model15.diff$obs_diff[model15.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model15.annforb

#   Annual grass
model15.anngrass <- model15.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model15.diff$obs_diff[model15.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model15.anngrass

#   Perennial forb
model15.perforb <- model15.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model15.diff$obs_diff[model15.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model15.perforb

#   Perennial grass
model15.pergrass <- model15.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model15.diff$obs_diff[model15.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass (**)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model15.pergrass

#   Shrub
model15.shrub <- model15.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model15.diff$obs_diff[model15.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model15.shrub

#   Shannon diversity
model15.shannon <- model15.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model15.diff$obs_diff[model15.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity (*)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model15.shannon

#   Bromus tectorum
model15.brte <- model15.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model15.diff$obs_diff[model15.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model15.brte

#   Artemisia
model15.artemisia <- model15.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model15.diff$obs_diff[model15.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia") ~ "(*)")) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model15.artemisia

#   Pinus and Juniperus
model15.pj <- model15.perm |> 
  filter(indicators == "Pinus and Juniperus") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model15.diff$obs_diff[model15.diff$indicators == "Pinus and Juniperus"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Pinus") ~ "&" ~ italic("Juniperus"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model15.pj

# Combine plots
grid.arrange(
  model15.bp, model15.annforb, model15.anngrass,
  model15.perforb, model15.pergrass, model15.shrub,
  model15.shannon, model15.brte, model15.artemisia, 
  model15.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)


### 16. Vegetation disturbance --------------------------------------------

# Filter for model
model16.matched <- all.matched |> 
  filter(Model == 16) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover, PJ_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Control", "Vegetation disturbance")))

#   pivot_longer() for cover & shannon cols
model16.matched <- model16.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model16.matched <- model16.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia",
             indicators == "PJ_cover" ~ "Pinus and Juniperus"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia",
                                        "Pinus and Juniperus")))

# Calculate observed mean difference
model16.diff <- model16.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model16.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model16.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model16.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values16 <- model16.perm |>
  inner_join(model16.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values16

# Boxplot
model16.bp <- model16.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "16. Northern Basin and Range: Vegetation disturbance") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")),
               "Pinus and Juniperus" = expression(italic("Pinus") ~ "&" ~ italic("Juniperus")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model16.bp

# Plot frequency distribution
#   Annual forb
model16.annforb <- model16.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model16.diff$obs_diff[model16.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model16.annforb

#   Annual grass
model16.anngrass <- model16.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model16.diff$obs_diff[model16.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model16.anngrass

#   Perennial forb
model16.perforb <- model16.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model16.diff$obs_diff[model16.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model16.perforb

#   Perennial grass
model16.pergrass <- model16.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model16.diff$obs_diff[model16.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model16.pergrass

#   Shrub
model16.shrub <- model16.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model16.diff$obs_diff[model16.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model16.shrub

#   Shannon diversity
model16.shannon <- model16.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model16.diff$obs_diff[model16.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model16.shannon

#   Bromus tectorum
model16.brte <- model16.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model16.diff$obs_diff[model16.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model16.brte

#   Artemisia
model16.artemisia <- model16.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model16.diff$obs_diff[model16.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model16.artemisia

#   Pinus and Juniperus
model16.pj <- model16.perm |> 
  filter(indicators == "Pinus and Juniperus") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model16.diff$obs_diff[model16.diff$indicators == "Pinus and Juniperus"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Pinus") ~ "&" ~ italic("Juniperus"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model16.pj

# Combine plots
grid.arrange(
  model16.bp, model16.annforb, model16.anngrass,
  model16.perforb, model16.pergrass, model16.shrub,
  model16.shannon, model16.brte, model16.artemisia, 
  model16.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)


### 17. Post-burn aerial seeding ------------------------------------------

# Filter for model
model17.matched <- all.matched |> 
  filter(Model == 17) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover, PJ_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn aerial seeding")))

#   pivot_longer() for cover & shannon cols
model17.matched <- model17.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model17.matched <- model17.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia",
             indicators == "PJ_cover" ~ "Pinus and Juniperus"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia",
                                        "Pinus and Juniperus")))

# Calculate observed mean difference
model17.diff <- model17.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model17.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model17.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model17.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values17 <- model17.perm |>
  inner_join(model17.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values17 # p = 0.001 for shannon

# Boxplot
model17.bp <- model17.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "17. Northern Basin and Range: Post-burn aerial seeding") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")),
               "Pinus and Juniperus" = expression(italic("Pinus") ~ "&" ~ italic("Juniperus")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model17.bp

# Plot frequency distribution
#   Annual forb
model17.annforb <- model17.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model17.diff$obs_diff[model17.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model17.annforb

#   Annual grass
model17.anngrass <- model17.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model17.diff$obs_diff[model17.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model17.anngrass

#   Perennial forb
model17.perforb <- model17.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model17.diff$obs_diff[model17.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model17.perforb

#   Perennial grass
model17.pergrass <- model17.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model17.diff$obs_diff[model17.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model17.pergrass

#   Shrub
model17.shrub <- model17.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model17.diff$obs_diff[model17.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model17.shrub

#   Shannon diversity
model17.shannon <- model17.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model17.diff$obs_diff[model17.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity (**)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model17.shannon

#   Bromus tectorum
model17.brte <- model17.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model17.diff$obs_diff[model17.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model17.brte

#   Artemisia
model17.artemisia <- model17.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model17.diff$obs_diff[model17.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model17.artemisia

#   Pinus and Juniperus
model17.pj <- model17.perm |> 
  filter(indicators == "Pinus and Juniperus") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model17.diff$obs_diff[model17.diff$indicators == "Pinus and Juniperus"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Pinus") ~ "&" ~ italic("Juniperus"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model17.pj

# Combine plots
grid.arrange(
  model17.bp, model17.annforb, model17.anngrass,
  model17.perforb, model17.pergrass, model17.shrub,
  model17.shannon, model17.brte, model17.artemisia, 
  model17.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)


### 18. Post-burn aerial & drill seeding ----------------------------------

# Filter for model
model18.matched <- all.matched |> 
  filter(Model == 18) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover, PJ_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn aerial & drill seeding")))

#   pivot_longer() for cover & shannon cols
model18.matched <- model18.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model18.matched <- model18.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia",
             indicators == "PJ_cover" ~ "Pinus and Juniperus"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia",
                                        "Pinus and Juniperus")))

# Calculate observed mean difference
model18.diff <- model18.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model18.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model18.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model18.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values18 <- model18.perm |>
  inner_join(model18.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values18 # p = 0.04 for shannon

# Boxplot
model18.bp <- model18.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "18. Northern Basin and Range: Post-burn aerial & drill seeding") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")),
               "Pinus and Juniperus" = expression(italic("Pinus") ~ "&" ~ italic("Juniperus")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model18.bp

# Plot frequency distribution
#   Annual forb
model18.annforb <- model18.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model18.diff$obs_diff[model18.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model18.annforb

#   Annual grass
model18.anngrass <- model18.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model18.diff$obs_diff[model18.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model18.anngrass

#   Perennial forb
model18.perforb <- model18.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model18.diff$obs_diff[model18.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model18.perforb

#   Perennial grass
model18.pergrass <- model18.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model18.diff$obs_diff[model18.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model18.pergrass

#   Shrub
model18.shrub <- model18.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model18.diff$obs_diff[model18.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model18.shrub

#   Shannon diversity
model18.shannon <- model18.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model18.diff$obs_diff[model18.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity (*)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model18.shannon

#   Bromus tectorum
model18.brte <- model18.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model18.diff$obs_diff[model18.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model18.brte

#   Artemisia
model18.artemisia <- model18.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model18.diff$obs_diff[model18.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model18.artemisia

#   Pinus and Juniperus
model18.pj <- model18.perm |> 
  filter(indicators == "Pinus and Juniperus") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model18.diff$obs_diff[model18.diff$indicators == "Pinus and Juniperus"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Pinus") ~ "&" ~ italic("Juniperus"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model18.pj

# Combine plots
grid.arrange(
  model18.bp, model18.annforb, model18.anngrass,
  model18.perforb, model18.pergrass, model18.shrub,
  model18.shannon, model18.brte, model18.artemisia, 
  model18.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)


### 19. Post-burn closure -------------------------------------------------

# Filter for model
model19.matched <- all.matched |> 
  filter(Model == 19) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover, PJ_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn closure")))

#   pivot_longer() for cover & shannon cols
model19.matched <- model19.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model19.matched <- model19.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia",
             indicators == "PJ_cover" ~ "Pinus and Juniperus"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia",
                                        "Pinus and Juniperus")))

# Calculate observed mean difference
model19.diff <- model19.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model19.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model19.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model19.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values19 <- model19.perm |>
  inner_join(model19.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values19 # p < 0.001 for annual & perennial grass; 
#   p = 0.002 for BRTE; p = 0.048 for Artemisia; p = 0.03 for PJ

# Boxplot
model19.bp <- model19.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "19. Northern Basin and Range: Post-burn closure") +
  geom_signif(
    y_position = 165,
    xmin = 1.8,
    xmax = 2.2, 
    annotations = c("***")
  ) +
  geom_signif(
    y_position = 135,
    xmin = 3.8,
    xmax = 4.2, 
    annotations = c("***")
  ) +
  geom_signif(
    y_position = 100,
    xmin = 5.8,
    xmax = 6.2, 
    annotations = c("***")
  ) +
  geom_signif(
    y_position = 65,
    xmin = 6.8,
    xmax = 7.2, 
    annotations = c("*")
  ) +
  geom_signif(
    y_position = 35,
    xmin = 7.8,
    xmax = 8.2, 
    annotations = c("*")
  ) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")),
               "Pinus and Juniperus" = expression(italic("Pinus") ~ "&" ~ italic("Juniperus")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model19.bp 

# Plot frequency distribution
#   Annual forb
model19.annforb <- model19.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model19.diff$obs_diff[model19.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model19.annforb

#   Annual grass
model19.anngrass <- model19.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model19.diff$obs_diff[model19.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass (***)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model19.anngrass

#   Perennial forb
model19.perforb <- model19.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model19.diff$obs_diff[model19.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model19.perforb

#   Perennial grass
model19.pergrass <- model19.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model19.diff$obs_diff[model19.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass (***)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model19.pergrass

#   Shrub
model19.shrub <- model19.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model19.diff$obs_diff[model19.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model19.shrub

#   Shannon diversity
model19.shannon <- model19.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model19.diff$obs_diff[model19.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model19.shannon

#   Bromus tectorum
model19.brte <- model19.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model19.diff$obs_diff[model19.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum") ~ "(***)")) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model19.brte

#   Artemisia
model19.artemisia <- model19.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model19.diff$obs_diff[model19.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia") ~ "(*)")) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model19.artemisia

#   Pinus and Juniperus
model19.pj <- model19.perm |> 
  filter(indicators == "Pinus and Juniperus") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model19.diff$obs_diff[model19.diff$indicators == "Pinus and Juniperus"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Pinus") ~ "&" ~ italic("Juniperus") ~ "(*)")) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model19.pj

# Combine plots
grid.arrange(
  model19.bp, model19.annforb, model19.anngrass,
  model19.perforb, model19.pergrass, model19.shrub,
  model19.shannon, model19.brte, model19.artemisia, 
  model19.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)


### 20. Post-burn drill seeding -------------------------------------------------

# Filter for model
model20.matched <- all.matched |> 
  filter(Model == 20) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover, PJ_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn drill seeding")))

#   pivot_longer() for cover & shannon cols
model20.matched <- model20.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model20.matched <- model20.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia",
             indicators == "PJ_cover" ~ "Pinus and Juniperus"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia",
                                        "Pinus and Juniperus")))

# Calculate observed mean difference
model20.diff <- model20.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model20.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model20.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model20.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values20 <- model20.perm |>
  inner_join(model20.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values20 # p = 0.009 for perennial forb; p = 0.02 for shannon;
#   p = 0.0006 for PJ

# Boxplot
model20.bp <- model20.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "20. Northern Basin and Range: Post-burn drill seeding") +
  geom_signif(
    y_position = 65,
    xmin = 2.8,
    xmax = 3.2, 
    annotations = c("**")
  ) +
  geom_signif(
    y_position = 30,
    xmin = 7.8,
    xmax = 8.2, 
    annotations = c("***")
  ) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")),
               "Pinus and Juniperus" = expression(italic("Pinus") ~ "&" ~ italic("Juniperus")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model20.bp

# Plot frequency distribution
#   Annual forb
model20.annforb <- model20.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model20.diff$obs_diff[model20.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model20.annforb

#   Annual grass
model20.anngrass <- model20.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model20.diff$obs_diff[model20.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model20.anngrass

#   Perennial forb
model20.perforb <- model20.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model20.diff$obs_diff[model20.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb (**)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model20.perforb

#   Perennial grass
model20.pergrass <- model20.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model20.diff$obs_diff[model20.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model20.pergrass

#   Shrub
model20.shrub <- model20.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model20.diff$obs_diff[model20.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model20.shrub

#   Shannon diversity
model20.shannon <- model20.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model20.diff$obs_diff[model20.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity (*)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model20.shannon

#   Bromus tectorum
model20.brte <- model20.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model20.diff$obs_diff[model20.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model20.brte

#   Artemisia
model20.artemisia <- model20.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model20.diff$obs_diff[model20.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model20.artemisia

#   Pinus and Juniperus
model20.pj <- model20.perm |> 
  filter(indicators == "Pinus and Juniperus") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model20.diff$obs_diff[model20.diff$indicators == "Pinus and Juniperus"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Pinus") ~ "&" ~ italic("Juniperus") ~ "(***)")) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model20.pj

# Combine plots
grid.arrange(
  model20.bp, model20.annforb, model20.anngrass,
  model20.perforb, model20.pergrass, model20.shrub,
  model20.shannon, model20.brte, model20.artemisia, 
  model20.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)


### 21. Post-burn herbicide -----------------------------------------------

# Filter for model
model21.matched <- all.matched |> 
  filter(Model == 21) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover, PJ_cover, 
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn herbicide")))

#   pivot_longer() for cover & shannon cols
model21.matched <- model21.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model21.matched <- model21.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia",
             indicators == "PJ_cover" ~ "Pinus and Juniperus"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia",
                                        "Pinus and Juniperus")))

# Calculate observed mean difference
model21.diff <- model21.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model21.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model21.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model21.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values21 <- model21.perm |>
  inner_join(model21.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values21 # p < 0.0001 for shannon; p = 0.04 for BRTE

# Boxplot
model21.bp <- model21.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "21. Northern Basin and Range: Post-burn herbicide") +
  geom_signif(
    y_position = 110,
    xmin = 5.8,
    xmax = 6.2, 
    annotations = c("*")
  ) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")),
               "Pinus and Juniperus" = expression(italic("Pinus") ~ "&" ~ italic("Juniperus")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model21.bp

# Plot frequency distribution
#   Annual forb
model21.annforb <- model21.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model21.diff$obs_diff[model21.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model21.annforb

#   Annual grass
model21.anngrass <- model21.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model21.diff$obs_diff[model21.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model21.anngrass

#   Perennial forb
model21.perforb <- model21.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model21.diff$obs_diff[model21.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model21.perforb

#   Perennial grass
model21.pergrass <- model21.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model21.diff$obs_diff[model21.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model21.pergrass

#   Shrub
model21.shrub <- model21.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model21.diff$obs_diff[model21.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model21.shrub

#   Shannon diversity
model21.shannon <- model21.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model21.diff$obs_diff[model21.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity (***)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model21.shannon

#   Bromus tectorum
model21.brte <- model21.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model21.diff$obs_diff[model21.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum") ~ "(*)")) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model21.brte

#   Artemisia
model21.artemisia <- model21.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model21.diff$obs_diff[model21.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model21.artemisia

#   Pinus and Juniperus
model21.pj <- model21.perm |> 
  filter(indicators == "Pinus and Juniperus") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model21.diff$obs_diff[model21.diff$indicators == "Pinus and Juniperus"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Pinus") ~ "&" ~ italic("Juniperus"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model21.pj

# Combine plots
grid.arrange(
  model21.bp, model21.annforb, model21.anngrass,
  model21.perforb, model21.pergrass, model21.shrub,
  model21.shannon, model21.brte, model21.artemisia, 
  model21.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)


### 22. Post-burn seedling planting ---------------------------------------

# Filter for model
model22.matched <- all.matched |> 
  filter(Model == 22) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover, PJ_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn seedling planting")))

#   pivot_longer() for cover & shannon cols
model22.matched <- model22.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model22.matched <- model22.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia",
             indicators == "PJ_cover" ~ "Pinus and Juniperus"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia",
                                        "Pinus and Juniperus")))

# Calculate observed mean difference
model22.diff <- model22.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model22.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model22.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model22.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values22 <- model22.perm |>
  inner_join(model22.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values22 # p = 0.03 for perennial grass; p = 0.04 for BRTE

# Boxplot
model22.bp <- model22.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "22. Northern Basin and Range: Post-burn seedling planting") +
  geom_signif(
    y_position = 115,
    xmin = 3.8,
    xmax = 4.2, 
    annotations = c("*")
  ) +
  geom_signif(
    y_position = 90,
    xmin = 5.8,
    xmax = 6.2, 
    annotations = c("*")
  ) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")),
               "Pinus and Juniperus" = expression(italic("Pinus") ~ "&" ~ italic("Juniperus")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model22.bp

# Plot frequency distribution
#   Annual forb
model22.annforb <- model22.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model22.diff$obs_diff[model22.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model22.annforb

#   Annual grass
model22.anngrass <- model22.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model22.diff$obs_diff[model22.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model22.anngrass

#   Perennial forb
model22.perforb <- model22.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model22.diff$obs_diff[model22.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model22.perforb

#   Perennial grass
model22.pergrass <- model22.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model22.diff$obs_diff[model22.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass (*)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model22.pergrass

#   Shrub
model22.shrub <- model22.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model22.diff$obs_diff[model22.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model22.shrub

#   Shannon diversity
model22.shannon <- model22.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model22.diff$obs_diff[model22.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity (*)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model22.shannon

#   Bromus tectorum
model22.brte <- model22.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model22.diff$obs_diff[model22.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum") ~ "(*)")) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model22.brte

# Combine plots
#   Artemisia
model22.artemisia <- model22.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model22.diff$obs_diff[model22.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model22.artemisia

#   Pinus and Juniperus
model22.pj <- model22.perm |> 
  filter(indicators == "Pinus and Juniperus") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model22.diff$obs_diff[model22.diff$indicators == "Pinus and Juniperus"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Pinus") ~ "&" ~ italic("Juniperus"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model22.pj

# Combine plots
grid.arrange(
  model22.bp, model22.annforb, model22.anngrass,
  model22.perforb, model22.pergrass, model22.shrub,
  model22.shannon, model22.brte, model22.artemisia, 
  model22.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)



## Central Basin and Range ------------------------------------------------

### 23. Drill seeding & soil disturbance ----------------------------------

# Filter for model
model23.matched <- all.matched |> 
  filter(Model == 23) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover, PJ_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Control", "Drill seeding & soil disturbance")))

#   pivot_longer() for cover & shannon cols
model23.matched <- model23.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model23.matched <- model23.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia",
             indicators == "PJ_cover" ~ "Pinus and Juniperus"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia",
                                        "Pinus and Juniperus")))

# Calculate observed mean difference
model23.diff <- model23.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model23.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model23.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model23.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values23 <- model23.perm |>
  inner_join(model23.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values23 # p = 0.015 for perennial grass

# Boxplot
model23.bp <- model23.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "23. Central Basin and Range: Drill seeding & soil disturbance") +
  geom_signif(
    y_position = 95,
    xmin = 3.8,
    xmax = 4.2, 
    annotations = c("*")
  ) +
  ylim(0, 105) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")),
               "Pinus and Juniperus" = expression(italic("Pinus") ~ "&" ~ italic("Juniperus")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model23.bp

# Plot frequency distribution
#   Annual forb
model23.annforb <- model23.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model23.diff$obs_diff[model23.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model23.annforb

#   Annual grass
model23.anngrass <- model23.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model23.diff$obs_diff[model23.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model23.anngrass

#   Perennial forb
model23.perforb <- model23.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model23.diff$obs_diff[model23.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model23.perforb

#   Perennial grass
model23.pergrass <- model23.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model23.diff$obs_diff[model23.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass (*)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model23.pergrass

#   Shrub
model23.shrub <- model23.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model23.diff$obs_diff[model23.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model23.shrub

#   Shannon diversity
model23.shannon <- model23.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model23.diff$obs_diff[model23.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model23.shannon

#   Bromus tectorum
model23.brte <- model23.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model23.diff$obs_diff[model23.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model23.brte

#   Artemisia
model23.artemisia <- model23.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model23.diff$obs_diff[model23.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model23.artemisia

#   Pinus and Juniperus
model23.pj <- model23.perm |> 
  filter(indicators == "Pinus and Juniperus") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model23.diff$obs_diff[model23.diff$indicators == "Pinus and Juniperus"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Pinus") ~ "&" ~ italic("Juniperus"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model23.pj

# Combine plots
grid.arrange(
  model23.bp, model23.annforb, model23.anngrass,
  model23.perforb, model23.pergrass, model23.shrub,
  model23.shannon, model23.brte, model23.artemisia, 
  model23.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)


### 24. Prescribed burn ---------------------------------------------------

# Filter for model
model24.matched <- all.matched |> 
  filter(Model == 24) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover, PJ_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Control", "Prescribed burn")))

#   pivot_longer() for cover & shannon cols
model24.matched <- model24.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model24.matched <- model24.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia",
             indicators == "PJ_cover" ~ "Pinus and Juniperus"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia",
                                        "Pinus and Juniperus")))

# Calculate observed mean difference
model24.diff <- model24.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model24.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model24.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model24.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values24 <- model24.perm |>
  inner_join(model24.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values24 # p = 0.04 for perennial grass; p = 0.047 for shannon

# Boxplot
model24.bp <- model24.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "24. Central Basin and Range: Prescribed burn") +
  geom_signif(
    y_position = 55,
    xmin = 3.8,
    xmax = 4.2, 
    annotations = c("*")
  ) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")),
               "Pinus and Juniperus" = expression(italic("Pinus") ~ "&" ~ italic("Juniperus")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model24.bp

# Plot frequency distribution
#   Annual forb
model24.annforb <- model24.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model24.diff$obs_diff[model24.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model24.annforb

#   Annual grass
model24.anngrass <- model24.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model24.diff$obs_diff[model24.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model24.anngrass

#   Perennial forb
model24.perforb <- model24.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model24.diff$obs_diff[model24.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model24.perforb

#   Perennial grass
model24.pergrass <- model24.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model24.diff$obs_diff[model24.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass (*)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model24.pergrass

#   Shrub
model24.shrub <- model24.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model24.diff$obs_diff[model24.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model24.shrub

#   Shannon diversity
model24.shannon <- model24.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model24.diff$obs_diff[model24.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity (*)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model24.shannon

#   Bromus tectorum
model24.brte <- model24.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model24.diff$obs_diff[model24.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model24.brte

#   Artemisia
model24.artemisia <- model24.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model24.diff$obs_diff[model24.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model24.artemisia

#   Pinus and Juniperus
model24.pj <- model24.perm |> 
  filter(indicators == "Pinus and Juniperus") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model24.diff$obs_diff[model24.diff$indicators == "Pinus and Juniperus"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Pinus") ~ "&" ~ italic("Juniperus"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model24.pj

# Combine plots
grid.arrange(
  model24.bp, model24.annforb, model24.anngrass,
  model24.perforb, model24.pergrass, model24.shrub,
  model24.shannon, model24.brte, model24.artemisia, 
  model24.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)


### 25. Vegetation disturbance --------------------------------------------

# Filter for model
model25.matched <- all.matched |> 
  filter(Model == 25) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover, PJ_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Control", "Vegetation disturbance")))

#   pivot_longer() for cover & shannon cols
model25.matched <- model25.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model25.matched <- model25.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia",
             indicators == "PJ_cover" ~ "Pinus and Juniperus"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia",
                                        "Pinus and Juniperus")))

# Calculate observed mean difference
model25.diff <- model25.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model25.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model25.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model25.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values25 <- model25.perm |>
  inner_join(model25.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values25 # p = 0.001 for shannon; p = 0.03 for Artemisia;
#   p = 0.002 for PJ

# Boxplot
model25.bp <- model25.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "25. Central Basin and Range: Vegetation disturbance") +
  geom_signif(
    y_position = 60,
    xmin = 6.8,
    xmax = 7.2, 
    annotations = c("*")
  ) +
  geom_signif(
    y_position = 90,
    xmin = 7.8,
    xmax = 8.2, 
    annotations = c("**")
  ) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")),
               "Pinus and Juniperus" = expression(italic("Pinus") ~ "&" ~ italic("Juniperus")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model25.bp

# Plot frequency distribution
#   Annual forb
model25.annforb <- model25.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model25.diff$obs_diff[model25.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model25.annforb

#   Annual grass
model25.anngrass <- model25.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model25.diff$obs_diff[model25.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model25.anngrass

#   Perennial forb
model25.perforb <- model25.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model25.diff$obs_diff[model25.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model25.perforb

#   Perennial grass
model25.pergrass <- model25.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model25.diff$obs_diff[model25.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model25.pergrass

#   Shrub
model25.shrub <- model25.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model25.diff$obs_diff[model25.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model25.shrub

#   Shannon diversity
model25.shannon <- model25.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model25.diff$obs_diff[model25.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity (***)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model25.shannon

#   Bromus tectorum
model25.brte <- model25.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model25.diff$obs_diff[model25.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model25.brte

#   Artemisia
model25.artemisia <- model25.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model25.diff$obs_diff[model25.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia") ~ "(*)")) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model25.artemisia

#   Pinus and Juniperus
model25.pj <- model25.perm |> 
  filter(indicators == "Pinus and Juniperus") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model25.diff$obs_diff[model25.diff$indicators == "Pinus and Juniperus"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Pinus") ~ "&" ~ italic("Juniperus") ~ "(**)")) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model25.pj

# Combine plots
grid.arrange(
  model25.bp, model25.annforb, model25.anngrass,
  model25.perforb, model25.pergrass, model25.shrub,
  model25.shannon, model25.brte, model25.artemisia, 
  model25.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)


### 26. Post-burn aerial seeding ------------------------------------------

# Filter for model
model26.matched <- all.matched |> 
  filter(Model == 26) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover, PJ_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn aerial seeding")))

#   pivot_longer() for cover & shannon cols
model26.matched <- model26.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model26.matched <- model26.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia",
             indicators == "PJ_cover" ~ "Pinus and Juniperus"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia",
                                        "Pinus and Juniperus")))

# Calculate observed mean difference
model26.diff <- model26.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model26.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model26.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model26.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values26 <- model26.perm |>
  inner_join(model26.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values26 # p < 0.0001 for perennial forb; p = 0.0003 for shannon;
#   p = 0.006 for Artemisia

# Boxplot
model26.bp <- model26.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "26. Central Basin and Range: Post-burn aerial seeding") +
  geom_signif(
    y_position = 60,
    xmin = 2.8,
    xmax = 3.2, 
    annotations = c("***")
  ) +
  geom_signif(
    y_position = 55,
    xmin = 6.8,
    xmax = 7.2, 
    annotations = c("**")
  ) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")),
               "Pinus and Juniperus" = expression(italic("Pinus") ~ "&" ~ italic("Juniperus")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model26.bp

# Plot frequency distribution
#   Annual forb
model26.annforb <- model26.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model26.diff$obs_diff[model26.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model26.annforb

#   Annual grass
model26.anngrass <- model26.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model26.diff$obs_diff[model26.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model26.anngrass

#   Perennial forb
model26.perforb <- model26.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model26.diff$obs_diff[model26.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb (***)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model26.perforb

#   Perennial grass
model26.pergrass <- model26.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model26.diff$obs_diff[model26.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model26.pergrass

#   Shrub
model26.shrub <- model26.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model26.diff$obs_diff[model26.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model26.shrub

#   Shannon diversity
model26.shannon <- model26.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model26.diff$obs_diff[model26.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity (***)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model26.shannon

#   Bromus tectorum
model26.brte <- model26.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model26.diff$obs_diff[model26.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model26.brte

#   Artemisia
model26.artemisia <- model26.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model26.diff$obs_diff[model26.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia") ~ "(**)")) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model26.artemisia

#   Pinus and Juniperus
model26.pj <- model26.perm |> 
  filter(indicators == "Pinus and Juniperus") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model26.diff$obs_diff[model26.diff$indicators == "Pinus and Juniperus"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Pinus") ~ "&" ~ italic("Juniperus"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model26.pj

# Combine plots
grid.arrange(
  model26.bp, model26.annforb, model26.anngrass,
  model26.perforb, model26.pergrass, model26.shrub,
  model26.shannon, model26.brte, model26.artemisia, 
  model26.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)


### 27. Post-burn drill seeding -------------------------------------------

# Filter for model
model27.matched <- all.matched |> 
  filter(Model == 27) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover, PJ_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn drill seeding")))

#   pivot_longer() for cover & shannon cols
model27.matched <- model27.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model27.matched <- model27.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia",
             indicators == "PJ_cover" ~ "Pinus and Juniperus"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia",
                                        "Pinus and Juniperus")))

# Calculate observed mean difference
model27.diff <- model27.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model27.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model27.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model27.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values27 <- model27.perm |>
  inner_join(model27.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values27 # p = 0.0098 for PJ

# Boxplot
model27.bp <- model27.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "27. Central Basin and Range: Post-burn drill seeding") +
  geom_signif(
    y_position = 63,
    xmin = 7.8,
    xmax = 8.2, 
    annotations = c("**")
  ) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")),
               "Pinus and Juniperus" = expression(italic("Pinus") ~ "&" ~ italic("Juniperus")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model27.bp

# Plot frequency distribution
#   Annual forb
model27.annforb <- model27.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model27.diff$obs_diff[model27.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model27.annforb

#   Annual grass
model27.anngrass <- model27.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model27.diff$obs_diff[model27.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model27.anngrass

#   Perennial forb
model27.perforb <- model27.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model27.diff$obs_diff[model27.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model27.perforb

#   Perennial grass
model27.pergrass <- model27.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model27.diff$obs_diff[model27.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model27.pergrass

#   Shrub
model27.shrub <- model27.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model27.diff$obs_diff[model27.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model27.shrub

#   Shannon diversity
model27.shannon <- model27.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model27.diff$obs_diff[model27.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model27.shannon

#   Bromus tectorum
model27.brte <- model27.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model27.diff$obs_diff[model27.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model27.brte

#   Artemisia
model27.artemisia <- model27.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model27.diff$obs_diff[model27.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model27.artemisia

#   Pinus and Juniperus
model27.pj <- model27.perm |> 
  filter(indicators == "Pinus and Juniperus") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model27.diff$obs_diff[model27.diff$indicators == "Pinus and Juniperus"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Pinus") ~ "&" ~ italic("Juniperus") ~ "(**)")) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model27.pj

# Combine plots
grid.arrange(
  model27.bp, model27.annforb, model27.anngrass,
  model27.perforb, model27.pergrass, model27.shrub,
  model27.shannon, model27.brte, model27.artemisia, 
  model27.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)


### 28. Post-burn ground seeding ------------------------------------------

# Filter for model
model28.matched <- all.matched |> 
  filter(Model == 28) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover, PJ_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn ground seeding")))

#   pivot_longer() for cover & shannon cols
model28.matched <- model28.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model28.matched <- model28.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia",
             indicators == "PJ_cover" ~ "Pinus and Juniperus"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia",
                                        "Pinus and Juniperus")))

# Calculate observed mean difference
model28.diff <- model28.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model28.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model28.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model28.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values28 <- model28.perm |>
  inner_join(model28.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values28

# Boxplot
model28.bp <- model28.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "28. Central Basin and Range: Post-burn ground seeding") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")),
               "Pinus and Juniperus" = expression(italic("Pinus") ~ "&" ~ italic("Juniperus")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model28.bp

# Plot frequency distribution
#   Annual forb
model28.annforb <- model28.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model28.diff$obs_diff[model28.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model28.annforb

#   Annual grass
model28.anngrass <- model28.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model28.diff$obs_diff[model28.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model28.anngrass

#   Perennial forb
model28.perforb <- model28.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model28.diff$obs_diff[model28.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model28.perforb

#   Perennial grass
model28.pergrass <- model28.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model28.diff$obs_diff[model28.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model28.pergrass

#   Shrub
model28.shrub <- model28.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model28.diff$obs_diff[model28.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model28.shrub

#   Shannon diversity
model28.shannon <- model28.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model28.diff$obs_diff[model28.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model28.shannon

#   Bromus tectorum
model28.brte <- model28.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model28.diff$obs_diff[model28.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model28.brte

#   Artemisia
model28.artemisia <- model28.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model28.diff$obs_diff[model28.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model28.artemisia

#   Pinus and Juniperus
model28.pj <- model28.perm |> 
  filter(indicators == "Pinus and Juniperus") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model28.diff$obs_diff[model28.diff$indicators == "Pinus and Juniperus"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Pinus") ~ "&" ~ italic("Juniperus"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model28.pj

# Combine plots
grid.arrange(
  model28.bp, model28.annforb, model28.anngrass,
  model28.perforb, model28.pergrass, model28.shrub,
  model28.shannon, model28.brte, model28.artemisia, 
  model28.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)


### 29. Post-burn herbicide -----------------------------------------------

# Filter for model
model29.matched <- all.matched |> 
  filter(Model == 29) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover, PJ_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn herbicide")))

#   pivot_longer() for cover & shannon cols
model29.matched <- model29.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model29.matched <- model29.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia",
             indicators == "PJ_cover" ~ "Pinus and Juniperus"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia",
                                        "Pinus and Juniperus")))

# Calculate observed mean difference
model29.diff <- model29.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model29.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model29.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model29.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values29 <- model29.perm |>
  inner_join(model29.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values29

# Boxplot
model29.bp <- model29.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "29. Central Basin and Range: Post-burn herbicide") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")),
               "Pinus and Juniperus" = expression(italic("Pinus") ~ "&" ~ italic("Juniperus")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model29.bp

# Plot frequency distribution
#   Annual forb
model29.annforb <- model29.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model29.diff$obs_diff[model29.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model29.annforb

#   Annual grass
model29.anngrass <- model29.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model29.diff$obs_diff[model29.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model29.anngrass

#   Perennial forb
model29.perforb <- model29.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model29.diff$obs_diff[model29.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model29.perforb

#   Perennial grass
model29.pergrass <- model29.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model29.diff$obs_diff[model29.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model29.pergrass

#   Shrub
model29.shrub <- model29.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model29.diff$obs_diff[model29.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model29.shrub

#   Shannon diversity
model29.shannon <- model29.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model29.diff$obs_diff[model29.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model29.shannon

#   Bromus tectorum
model29.brte <- model29.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model29.diff$obs_diff[model29.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model29.brte

#   Artemisia
model29.artemisia <- model29.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model29.diff$obs_diff[model29.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model29.artemisia

#   Pinus and Juniperus
model29.pj <- model29.perm |> 
  filter(indicators == "Pinus and Juniperus") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model29.diff$obs_diff[model29.diff$indicators == "Pinus and Juniperus"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Pinus") ~ "&" ~ italic("Juniperus"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model29.pj

# Combine plots
grid.arrange(
  model29.bp, model29.annforb, model29.anngrass,
  model29.perforb, model29.pergrass, model29.shrub,
  model29.shannon, model29.brte, model29.artemisia, 
  model29.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)


## Wyoming Basin ----------------------------------------------------------

### 30. Prescribed burn ---------------------------------------------------

# Filter for model
model30.matched <- all.matched |> 
  filter(Model == 30) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Control", "Prescribed burn")))

#   pivot_longer() for cover & shannon cols
model30.matched <- model30.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model30.matched <- model30.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia")))

# Calculate observed mean difference
model30.diff <- model30.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model30.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model30.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model30.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values30 <- model30.perm |>
  inner_join(model30.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values30

# Boxplot
model30.bp <- model30.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "30. Wyoming Basin: Prescribed burn") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model30.bp

# Plot frequency distribution
#   Annual forb
model30.annforb <- model30.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model30.diff$obs_diff[model30.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model30.annforb

#   Annual grass
model30.anngrass <- model30.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model30.diff$obs_diff[model30.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model30.anngrass

#   Perennial forb
model30.perforb <- model30.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model30.diff$obs_diff[model30.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model30.perforb

#   Perennial grass
model30.pergrass <- model30.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model30.diff$obs_diff[model30.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model30.pergrass

#   Shrub
model30.shrub <- model30.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model30.diff$obs_diff[model30.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model30.shrub

#   Shannon diversity
model30.shannon <- model30.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model30.diff$obs_diff[model30.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model30.shannon

#   Bromus tectorum
model30.brte <- model30.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model30.diff$obs_diff[model30.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model30.brte

#   Artemisia
model30.artemisia <- model30.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model30.diff$obs_diff[model30.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model30.artemisia

# Combine plots
grid.arrange(
  model30.bp, model30.annforb, model30.anngrass,
  model30.perforb, model30.pergrass, model30.shrub,
  model30.shannon, model30.brte, model30.artemisia,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(NA, 8, 8, 9, 9, NA)
  )
)



## Colorado Plateaus ------------------------------------------------------

### 31. Aerial seeding & soil disturbance ---------------------------------

# Filter for model
model31.matched <- all.matched |> 
  filter(Model == 31) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover, PJ_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Control", "Aerial seeding & soil disturbance")))

#   pivot_longer() for cover & shannon cols
model31.matched <- model31.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model31.matched <- model31.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia",
             indicators == "PJ_cover" ~ "Pinus and Juniperus"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia",
                                        "Pinus and Juniperus")))

# Calculate observed mean difference
model31.diff <- model31.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model31.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model31.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model31.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values31 <- model31.perm |>
  inner_join(model31.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values31 # p = 0.04 for perennial forb; p = 0.007 for shannon

# Boxplot
model31.bp <- model31.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "31. Colorado Plateaus: Aerial seeding & soil disturbance") +
  geom_signif(
    y_position = 52,
    xmin = 2.8,
    xmax = 3.2, 
    annotations = c("*")
  ) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")),
               "Pinus and Juniperus" = expression(italic("Pinus") ~ "&" ~ italic("Juniperus")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model31.bp

# Plot frequency distribution
#   Annual forb
model31.annforb <- model31.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model31.diff$obs_diff[model31.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model31.annforb

#   Annual grass
model31.anngrass <- model31.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model31.diff$obs_diff[model31.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model31.anngrass

#   Perennial forb
model31.perforb <- model31.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model31.diff$obs_diff[model31.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb (*)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model31.perforb

#   Perennial grass
model31.pergrass <- model31.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model31.diff$obs_diff[model31.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model31.pergrass

#   Shrub
model31.shrub <- model31.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model31.diff$obs_diff[model31.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model31.shrub

#   Shannon diversity
model31.shannon <- model31.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model31.diff$obs_diff[model31.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity (**)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model31.shannon

#   Bromus tectorum
model31.brte <- model31.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model31.diff$obs_diff[model31.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model31.brte

#   Artemisia
model31.artemisia <- model31.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model31.diff$obs_diff[model31.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model31.artemisia

#   Pinus and Juniperus
model31.pj <- model31.perm |> 
  filter(indicators == "Pinus and Juniperus") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model31.diff$obs_diff[model31.diff$indicators == "Pinus and Juniperus"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Pinus") ~ "&" ~ italic("Juniperus"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model31.pj

# Combine plots
grid.arrange(
  model31.bp, model31.annforb, model31.anngrass,
  model31.perforb, model31.pergrass, model31.shrub,
  model31.shannon, model31.brte, model31.artemisia, 
  model31.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)


### 32. Herbicide ---------------------------------------------------------

# Filter for model
model32.matched <- all.matched |> 
  filter(Model == 32) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover, PJ_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Control", "Herbicide")))

#   pivot_longer() for cover & shannon cols
model32.matched <- model32.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model32.matched <- model32.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia",
             indicators == "PJ_cover" ~ "Pinus and Juniperus"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia",
                                        "Pinus and Juniperus")))

# Calculate observed mean difference
model32.diff <- model32.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model32.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model32.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model32.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values32 <- model32.perm |>
  inner_join(model32.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values32 # p < 0.0001 for Artemisia

# Boxplot
model32.bp <- model32.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "32. Colorado Plateaus: Herbicide") +
  geom_signif(
    y_position = 75,
    xmin = 6.8,
    xmax = 7.2, 
    annotations = c("***")
  ) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")),
               "Pinus and Juniperus" = expression(italic("Pinus") ~ "&" ~ italic("Juniperus")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model32.bp

# Plot frequency distribution
#   Annual forb
model32.annforb <- model32.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model32.diff$obs_diff[model32.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model32.annforb

#   Annual grass
model32.anngrass <- model32.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model32.diff$obs_diff[model32.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model32.anngrass

#   Perennial forb
model32.perforb <- model32.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model32.diff$obs_diff[model32.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model32.perforb

#   Perennial grass
model32.pergrass <- model32.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model32.diff$obs_diff[model32.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model32.pergrass

#   Shrub
model32.shrub <- model32.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model32.diff$obs_diff[model32.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model32.shrub

#   Shannon diversity
model32.shannon <- model32.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model32.diff$obs_diff[model32.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model32.shannon

#   Bromus tectorum
model32.brte <- model32.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model32.diff$obs_diff[model32.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model32.brte

#   Artemisia
model32.artemisia <- model32.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model32.diff$obs_diff[model32.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia") ~ "(***)")) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model32.artemisia

#   Pinus and Juniperus
model32.pj <- model32.perm |> 
  filter(indicators == "Pinus and Juniperus") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model32.diff$obs_diff[model32.diff$indicators == "Pinus and Juniperus"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Pinus") ~ "&" ~ italic("Juniperus"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model32.pj

# Combine plots
grid.arrange(
  model32.bp, model32.annforb, model32.anngrass,
  model32.perforb, model32.pergrass, model32.shrub,
  model32.shannon, model32.brte, model32.artemisia, 
  model32.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)


### 33. Prescribed burn ---------------------------------------------------

# Filter for model
model33.matched <- all.matched |> 
  filter(Model == 33) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover, PJ_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Control", "Prescribed burn")))

#   pivot_longer() for cover & shannon cols
model33.matched <- model33.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model33.matched <- model33.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia",
             indicators == "PJ_cover" ~ "Pinus and Juniperus"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia",
                                        "Pinus and Juniperus")))

# Calculate observed mean difference
model33.diff <- model33.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model33.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model33.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model33.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values33 <- model33.perm |>
  inner_join(model33.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values33

# Boxplot
model33.bp <- model33.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "33. Colorado Plateaus: Prescribed burn") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")),
               "Pinus and Juniperus" = expression(italic("Pinus") ~ "&" ~ italic("Juniperus")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model33.bp

# Plot frequency distribution
#   Annual forb
model33.annforb <- model33.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model33.diff$obs_diff[model33.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model33.annforb

#   Annual grass
model33.anngrass <- model33.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model33.diff$obs_diff[model33.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model33.anngrass

#   Perennial forb
model33.perforb <- model33.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model33.diff$obs_diff[model33.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model33.perforb

#   Perennial grass
model33.pergrass <- model33.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model33.diff$obs_diff[model33.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model33.pergrass

#   Shrub
model33.shrub <- model33.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model33.diff$obs_diff[model33.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model33.shrub

#   Shannon diversity
model33.shannon <- model33.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model33.diff$obs_diff[model33.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model33.shannon

#   Bromus tectorum
model33.brte <- model33.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model33.diff$obs_diff[model33.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model33.brte

#   Artemisia
model33.artemisia <- model33.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model33.diff$obs_diff[model33.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model33.artemisia

#   Pinus and Juniperus
model33.pj <- model33.perm |> 
  filter(indicators == "Pinus and Juniperus") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model33.diff$obs_diff[model33.diff$indicators == "Pinus and Juniperus"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Pinus") ~ "&" ~ italic("Juniperus"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model33.pj

# Combine plots
grid.arrange(
  model33.bp, model33.annforb, model33.anngrass,
  model33.perforb, model33.pergrass, model33.shrub,
  model33.shannon, model33.brte, model33.artemisia, 
  model33.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)


### 34. Soil disturbance --------------------------------------------------

# Filter for model
model34.matched <- all.matched |> 
  filter(Model == 34) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover, PJ_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Control", "Soil disturbance")))

#   pivot_longer() for cover & shannon cols
model34.matched <- model34.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model34.matched <- model34.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia",
             indicators == "PJ_cover" ~ "Pinus and Juniperus"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia",
                                        "Pinus and Juniperus")))

# Calculate observed mean difference
model34.diff <- model34.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model34.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model34.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model34.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values34 <- model34.perm |>
  inner_join(model34.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values34 # p = 0.047 for Artemisia

# Boxplot
model34.bp <- model34.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "34. Colorado Plateaus: Soil disturbance") +
  geom_signif(
    y_position = 75,
    xmin = 6.8,
    xmax = 7.2, 
    annotations = c("*")
  ) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")),
               "Pinus and Juniperus" = expression(italic("Pinus") ~ "&" ~ italic("Juniperus")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model34.bp

# Plot frequency distribution
#   Annual forb
model34.annforb <- model34.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model34.diff$obs_diff[model34.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model34.annforb

#   Annual grass
model34.anngrass <- model34.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model34.diff$obs_diff[model34.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model34.anngrass

#   Perennial forb
model34.perforb <- model34.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model34.diff$obs_diff[model34.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model34.perforb

#   Perennial grass
model34.pergrass <- model34.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model34.diff$obs_diff[model34.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model34.pergrass

#   Shrub
model34.shrub <- model34.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model34.diff$obs_diff[model34.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model34.shrub

#   Shannon diversity
model34.shannon <- model34.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model34.diff$obs_diff[model34.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model34.shannon

#   Bromus tectorum
model34.brte <- model34.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model34.diff$obs_diff[model34.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model34.brte

#   Artemisia
model34.artemisia <- model34.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model34.diff$obs_diff[model34.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia") ~ "(*)")) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model34.artemisia

#   Pinus and Juniperus
model34.pj <- model34.perm |> 
  filter(indicators == "Pinus and Juniperus") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model34.diff$obs_diff[model34.diff$indicators == "Pinus and Juniperus"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Pinus") ~ "&" ~ italic("Juniperus"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model34.pj

# Combine plots
grid.arrange(
  model34.bp, model34.annforb, model34.anngrass,
  model34.perforb, model34.pergrass, model34.shrub,
  model34.shannon, model34.brte, model34.artemisia, 
  model34.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)


### 35. Vegetation disturbance --------------------------------------------

# Filter for model
model35.matched <- all.matched |> 
  filter(Model == 35) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover, PJ_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Control", "Vegetation disturbance")))

#   pivot_longer() for cover & shannon cols
model35.matched <- model35.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model35.matched <- model35.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia",
             indicators == "PJ_cover" ~ "Pinus and Juniperus"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia",
                                        "Pinus and Juniperus")))

# Calculate observed mean difference
model35.diff <- model35.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model35.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model35.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model35.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values35 <- model35.perm |>
  inner_join(model35.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values35 # p = 0.002 for annual grass; p = 0.006 for BRTE

# Boxplot
model35.bp <- model35.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "35. Colorado Plateaus: Vegetation disturbance") +
  geom_signif(
    y_position = 100,
    xmin = 1.8,
    xmax = 2.2, 
    annotations = c("**")
  ) +
  geom_signif(
    y_position = 65,
    xmin = 5.8,
    xmax = 6.2, 
    annotations = c("**")
  ) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")),
               "Pinus and Juniperus" = expression(italic("Pinus") ~ "&" ~ italic("Juniperus")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model35.bp

# Plot frequency distribution
#   Annual forb
model35.annforb <- model35.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model35.diff$obs_diff[model35.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model35.annforb

#   Annual grass
model35.anngrass <- model35.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model35.diff$obs_diff[model35.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass (**)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model35.anngrass

#   Perennial forb
model35.perforb <- model35.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model35.diff$obs_diff[model35.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model35.perforb

#   Perennial grass
model35.pergrass <- model35.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model35.diff$obs_diff[model35.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model35.pergrass

#   Shrub
model35.shrub <- model35.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model35.diff$obs_diff[model35.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model35.shrub

#   Shannon diversity
model35.shannon <- model35.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model35.diff$obs_diff[model35.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model35.shannon

#   Bromus tectorum
model35.brte <- model35.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model35.diff$obs_diff[model35.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum") ~ "(**)")) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model35.brte

#   Artemisia
model35.artemisia <- model35.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model35.diff$obs_diff[model35.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model35.artemisia

#   Pinus and Juniperus
model35.pj <- model35.perm |> 
  filter(indicators == "Pinus and Juniperus") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model35.diff$obs_diff[model35.diff$indicators == "Pinus and Juniperus"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Pinus") ~ "&" ~ italic("Juniperus"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model35.pj

# Combine plots
grid.arrange(
  model35.bp, model35.annforb, model35.anngrass,
  model35.perforb, model35.pergrass, model35.shrub,
  model35.shannon, model35.brte, model35.artemisia, 
  model35.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)


### 36. Post-burn aerial seeding ------------------------------------------

# Filter for model
model36.matched <- all.matched |> 
  filter(Model == 36) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover, PJ_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn aerial seeding")))

#   pivot_longer() for cover & shannon cols
model36.matched <- model36.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model36.matched <- model36.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia",
             indicators == "PJ_cover" ~ "Pinus and Juniperus"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia",
                                        "Pinus and Juniperus")))

# Calculate observed mean difference
model36.diff <- model36.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model36.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model36.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model36.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values36 <- model36.perm |>
  inner_join(model36.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values36 # p < 0.001 for shannon; p = 0.04 for Artemisia

# Boxplot
model36.bp <- model36.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "36. Colorado Plateaus: Post-burn aerial seeding") +
  geom_signif(
    y_position = 50,
    xmin = 6.8,
    xmax = 7.2, 
    annotations = c("*")
  ) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")),
               "Pinus and Juniperus" = expression(italic("Pinus") ~ "&" ~ italic("Juniperus")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model36.bp

# Plot frequency distribution
#   Annual forb
model36.annforb <- model36.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model36.diff$obs_diff[model36.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model36.annforb

#   Annual grass
model36.anngrass <- model36.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model36.diff$obs_diff[model36.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model36.anngrass

#   Perennial forb
model36.perforb <- model36.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model36.diff$obs_diff[model36.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model36.perforb

#   Perennial grass
model36.pergrass <- model36.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model36.diff$obs_diff[model36.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model36.pergrass

#   Shrub
model36.shrub <- model36.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model36.diff$obs_diff[model36.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model36.shrub

#   Shannon diversity
model36.shannon <- model36.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model36.diff$obs_diff[model36.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity (***)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model36.shannon

#   Bromus tectorum
model36.brte <- model36.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model36.diff$obs_diff[model36.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model36.brte

#   Artemisia
model36.artemisia <- model36.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model36.diff$obs_diff[model36.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia") ~ "(*)")) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model36.artemisia

#   Pinus and Juniperus
model36.pj <- model36.perm |> 
  filter(indicators == "Pinus and Juniperus") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model36.diff$obs_diff[model36.diff$indicators == "Pinus and Juniperus"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Pinus") ~ "&" ~ italic("Juniperus"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model36.pj

# Combine plots
grid.arrange(
  model36.bp, model36.annforb, model36.anngrass,
  model36.perforb, model36.pergrass, model36.shrub,
  model36.shannon, model36.brte, model36.artemisia, 
  model36.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)



## Arizona/New Mexico Plateau ---------------------------------------------

### 37. Herbicide ---------------------------------------------------------

# Filter for model
model37.matched <- all.matched |> 
  filter(Model == 37) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover, PJ_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Control", "Herbicide")))

#   pivot_longer() for cover & shannon cols
model37.matched <- model37.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model37.matched <- model37.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia",
             indicators == "PJ_cover" ~ "Pinus and Juniperus"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia",
                                        "Pinus and Juniperus")))

# Calculate observed mean difference
model37.diff <- model37.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model37.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model37.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model37.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values37 <- model37.perm |>
  inner_join(model37.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values37

# Boxplot
model37.bp <- model37.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "37. AZ/NM Plateau: Herbicide") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")),
               "Pinus and Juniperus" = expression(italic("Pinus") ~ "&" ~ italic("Juniperus")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model37.bp

# Plot frequency distribution
#   Annual forb
model37.annforb <- model37.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model37.diff$obs_diff[model37.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model37.annforb

#   Annual grass
model37.anngrass <- model37.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model37.diff$obs_diff[model37.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model37.anngrass

#   Perennial forb
model37.perforb <- model37.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model37.diff$obs_diff[model37.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model37.perforb

#   Perennial grass
model37.pergrass <- model37.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model37.diff$obs_diff[model37.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model37.pergrass

#   Shrub
model37.shrub <- model37.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model37.diff$obs_diff[model37.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model37.shrub

#   Shannon diversity
model37.shannon <- model37.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model37.diff$obs_diff[model37.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model37.shannon

#   Bromus tectorum
model37.brte <- model37.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model37.diff$obs_diff[model37.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model37.brte

#   Artemisia
model37.artemisia <- model37.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model37.diff$obs_diff[model37.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model37.artemisia

#   Pinus and Juniperus
model37.pj <- model37.perm |> 
  filter(indicators == "Pinus and Juniperus") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model37.diff$obs_diff[model37.diff$indicators == "Pinus and Juniperus"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Pinus") ~ "&" ~ italic("Juniperus"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model37.pj

# Combine plots
grid.arrange(
  model37.bp, model37.annforb, model37.anngrass,
  model37.perforb, model37.pergrass, model37.shrub,
  model37.shannon, model37.brte, model37.artemisia, 
  model37.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)


### 38. Prescribed burn ---------------------------------------------------

# Filter for model
model38.matched <- all.matched |> 
  filter(Model == 38) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover, PJ_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Control", "Prescribed burn")))

#   pivot_longer() for cover & shannon cols
model38.matched <- model38.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model38.matched <- model38.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia",
             indicators == "PJ_cover" ~ "Pinus and Juniperus"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia",
                                        "Pinus and Juniperus")))

# Calculate observed mean difference
model38.diff <- model38.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model38.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model38.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model38.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values38 <- model38.perm |>
  inner_join(model38.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values38 # p = 0.02 for annual grass; p = 0.02 for BRTE

# Boxplot
model38.bp <- model38.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "38. AZ/NM Plateau: Prescribed burn") +
  geom_signif(
    y_position = 48,
    xmin = 1.8,
    xmax = 2.2, 
    annotations = c("*")
  ) +
  geom_signif(
    y_position = 45,
    xmin = 5.8,
    xmax = 6.2, 
    annotations = c("*")
  ) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")),
               "Pinus and Juniperus" = expression(italic("Pinus") ~ "&" ~ italic("Juniperus")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model38.bp

# Plot frequency distribution
#   Annual forb
model38.annforb <- model38.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model38.diff$obs_diff[model38.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model38.annforb

#   Annual grass
model38.anngrass <- model38.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model38.diff$obs_diff[model38.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass (*)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model38.anngrass

#   Perennial forb
model38.perforb <- model38.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model38.diff$obs_diff[model38.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model38.perforb

#   Perennial grass
model38.pergrass <- model38.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model38.diff$obs_diff[model38.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model38.pergrass

#   Shrub
model38.shrub <- model38.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model38.diff$obs_diff[model38.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model38.shrub

#   Shannon diversity
model38.shannon <- model38.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model38.diff$obs_diff[model38.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model38.shannon

#   Bromus tectorum
model38.brte <- model38.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model38.diff$obs_diff[model38.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum") ~ "(*)")) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model38.brte

#   Artemisia
model38.artemisia <- model38.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model38.diff$obs_diff[model38.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model38.artemisia

#   Pinus and Juniperus
model38.pj <- model38.perm |> 
  filter(indicators == "Pinus and Juniperus") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model38.diff$obs_diff[model38.diff$indicators == "Pinus and Juniperus"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Pinus") ~ "&" ~ italic("Juniperus"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model38.pj

# Combine plots
grid.arrange(
  model38.bp, model38.annforb, model38.anngrass,
  model38.perforb, model38.pergrass, model38.shrub,
  model38.shannon, model38.brte, model38.artemisia, 
  model38.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)


### 39. Soil disturbance --------------------------------------------------

# Filter for model
model39.matched <- all.matched |> 
  filter(Model == 39) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRTE_cover, Artemisia_cover, PJ_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Control", "Soil disturbance")))

#   pivot_longer() for cover & shannon cols
model39.matched <- model39.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model39.matched <- model39.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRTE_cover" ~ "Bromus tectorum",
             indicators == "Artemisia_cover" ~ "Artemisia",
             indicators == "PJ_cover" ~ "Pinus and Juniperus"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum",
                                        "Artemisia",
                                        "Pinus and Juniperus")))

# Calculate observed mean difference
model39.diff <- model39.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model39.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model39.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model39.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values39 <- model39.perm |>
  inner_join(model39.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values39 # p = 0.40 for shannon

# Boxplot
model39.bp <- model39.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "39. AZ/NM Mountains: Soil disturbance") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")),
               "Artemisia" = expression(italic("Artemisia")),
               "Pinus and Juniperus" = expression(italic("Pinus") ~ "&" ~ italic("Juniperus")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model39.bp

# Plot frequency distribution
#   Annual forb
model39.annforb <- model39.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model39.diff$obs_diff[model39.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model39.annforb

#   Annual grass
model39.anngrass <- model39.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model39.diff$obs_diff[model39.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model39.anngrass

#   Perennial forb
model39.perforb <- model39.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model39.diff$obs_diff[model39.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model39.perforb

#   Perennial grass
model39.pergrass <- model39.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model39.diff$obs_diff[model39.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model39.pergrass

#   Shrub
model39.shrub <- model39.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model39.diff$obs_diff[model39.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model39.shrub

#   Shannon diversity
model39.shannon <- model39.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model39.diff$obs_diff[model39.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity (*)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model39.shannon

#   Bromus tectorum
model39.brte <- model39.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model39.diff$obs_diff[model39.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model39.brte

#   Artemisia
model39.artemisia <- model39.perm |> 
  filter(indicators == "Artemisia") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model39.diff$obs_diff[model39.diff$indicators == "Artemisia"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Artemisia"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model39.artemisia

#   Pinus and Juniperus
model39.pj <- model39.perm |> 
  filter(indicators == "Pinus and Juniperus") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model39.diff$obs_diff[model39.diff$indicators == "Pinus and Juniperus"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Pinus") ~ "&" ~ italic("Juniperus"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model39.pj

# Combine plots
grid.arrange(
  model39.bp, model39.annforb, model39.anngrass,
  model39.perforb, model39.pergrass, model39.shrub,
  model39.shannon, model39.brte, model39.artemisia, 
  model39.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)




# Warm Deserts ------------------------------------------------------------

## Mojave Basin and Range -------------------------------------------------

### 40. Post-burn aerial seeding ------------------------------------------

# Filter for model
model40.matched <- all.matched |> 
  filter(Model == 40) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRRU2_cover, LATR2_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn aerial seeding")))

#   pivot_longer() for cover & shannon cols
model40.matched <- model40.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model40.matched <- model40.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRRU2_cover" ~ "Bromus rubens",
             indicators == "LATR2_cover" ~ "Larrea tridentata"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus rubens",
                                        "Larrea tridentata")))

# Calculate observed mean difference
model40.diff <- model40.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model40.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model40.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model40.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values40 <- model40.perm |>
  inner_join(model40.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values40 # p = 0.007 for BRRU2; p = 0.0012 for LATR2

# Boxplot
model40.bp <- model40.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "40. Mojave Basin and Range: Post-burn aerial seeding") +
  geom_signif(
    y_position = 80,
    xmin = 5.8,
    xmax = 6.2, 
    annotations = c("**")
  ) +
  geom_signif(
    y_position = 18,
    xmin = 6.8,
    xmax = 7.2, 
    annotations = c("**")
  ) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus rubens" = expression(italic("Bromus rubens")),
               "Larrea tridentata" = expression(italic("Larrea tridentata")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model40.bp

# Plot frequency distribution
#   Annual forb
model40.annforb <- model40.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model40.diff$obs_diff[model40.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model40.annforb

#   Annual grass
model40.anngrass <- model40.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model40.diff$obs_diff[model40.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model40.anngrass

#   Perennial forb
model40.perforb <- model40.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model40.diff$obs_diff[model40.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model40.perforb

#   Perennial grass
model40.pergrass <- model40.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model40.diff$obs_diff[model40.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model40.pergrass

#   Shrub
model40.shrub <- model40.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model40.diff$obs_diff[model40.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model40.shrub

#   Shannon diversity
model40.shannon <- model40.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model40.diff$obs_diff[model40.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model40.shannon

#   Bromus rubens
model40.brru2 <- model40.perm |> 
  filter(indicators == "Bromus rubens") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model40.diff$obs_diff[model40.diff$indicators == "Bromus rubens"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus rubens") ~ "(**)")) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model40.brru2

#   Larrea tridentata
model40.latr2 <- model40.perm |> 
  filter(indicators == "Larrea tridentata") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model40.diff$obs_diff[model40.diff$indicators == "Larrea tridentata"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Larrea tridentata") ~ "(**)")) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model40.latr2

# Combine plots
grid.arrange(
  model40.bp, model40.annforb, model40.anngrass,
  model40.perforb, model40.pergrass, model40.shrub,
  model40.shannon, model40.brru2, model40.latr2,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(NA, 8, 8, 9, 9, NA)
  )
)



## Chihuahuan Desert ------------------------------------------------------

### 41. Herbicide ---------------------------------------------------------

# Filter for model
model41.matched <- all.matched |> 
  filter(Model == 41) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, LATR2_cover, Prosopis_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Control", "Herbicide")))

#   pivot_longer() for cover & shannon cols
model41.matched <- model41.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model41.matched <- model41.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "Prosopis_cover" ~ "Prosopis",
             indicators == "LATR2_cover" ~ "Larrea tridentata"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Prosopis",
                                        "Larrea tridentata")))

# Calculate observed mean difference
model41.diff <- model41.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model41.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model41.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model41.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values41 <- model41.perm |>
  inner_join(model41.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values41

# Boxplot
model41.bp <- model41.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "41. Chihuahuan Desert: Herbicide") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Prosopis" = expression(italic("Prosopis")),
               "Larrea tridentata" = expression(italic("Larrea tridentata")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model41.bp

# Plot frequency distribution
#   Annual forb
model41.annforb <- model41.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model41.diff$obs_diff[model41.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model41.annforb

#   Annual grass
model41.anngrass <- model41.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model41.diff$obs_diff[model41.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model41.anngrass

#   Perennial forb
model41.perforb <- model41.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model41.diff$obs_diff[model41.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model41.perforb

#   Perennial grass
model41.pergrass <- model41.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model41.diff$obs_diff[model41.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model41.pergrass

#   Shrub
model41.shrub <- model41.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model41.diff$obs_diff[model41.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model41.shrub

#   Shannon diversity
model41.shannon <- model41.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model41.diff$obs_diff[model41.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model41.shannon

#   Prosopis 
model41.prosopis <- model41.perm |> 
  filter(indicators == "Prosopis") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model41.diff$obs_diff[model41.diff$indicators == "Prosopis"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Prosopis"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model41.prosopis

#   Larrea tridentata
model41.latr2 <- model41.perm |> 
  filter(indicators == "Larrea tridentata") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model41.diff$obs_diff[model41.diff$indicators == "Larrea tridentata"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Larrea tridentata"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model41.latr2

# Combine plots
grid.arrange(
  model41.bp, model41.annforb, model41.anngrass,
  model41.perforb, model41.pergrass, model41.shrub,
  model41.shannon, model41.prosopis, model41.latr2,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(NA, 8, 8, 9, 9, NA)
  )
)




# Temperate Sierras / Upper Gila ------------------------------------------

## Arizona/New Mexico Mountains -------------------------------------------

### 42. Prescribed Burn ---------------------------------------------------

# Filter for model
model42.matched <- all.matched |> 
  filter(Model == 42) |> 
  select(LDCpointID, PrimaryKey, trt_control, Shannon, BRRU2_cover, PJ_cover,
         AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH, PerShrubCover_AH) |> 
  mutate(trt_control = factor(trt_control, levels = c("Control", "Prescribed burn")))

#   pivot_longer() for cover & shannon cols
model42.matched <- model42.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model42.matched <- model42.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "PerShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "BRRU2_cover" ~ "Bromus rubens",
             indicators == "PJ_cover" ~ "Pinus and Juniperus"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus rubens",
                                        "Pinus and Juniperus")))

# Calculate observed mean difference
model42.diff <- model42.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model42.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model42.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model42.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values42 <- model42.perm |>
  inner_join(model42.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values42 # p = 0.006 for shannon; p = 0.02 for BRRU2

# Boxplot
model42.bp <- model42.matched |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "42. AZ/NM Mountains: Prescribed burn") +
  theme(legend.title = element_blank()) +
  geom_signif(
    y_position = 30,
    xmin = 5.8,
    xmax = 6.2,
    annotations = c("*")
  ) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus rubens" = expression(italic("Bromus rubens")),
               "Pinus and Juniperus" = expression(italic("Pinus") ~ "&" ~ italic("Juniperus")))) +
  theme(plot.margin = margin(10, 10, 10, 10))
model42.bp

# Plot frequency distribution
#   Annual forb
model42.annforb <- model42.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model42.diff$obs_diff[model42.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model42.annforb

#   Annual grass
model42.anngrass <- model42.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model42.diff$obs_diff[model42.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model42.anngrass

#   Perennial forb
model42.perforb <- model42.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model42.diff$obs_diff[model42.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model42.perforb

#   Perennial grass
model42.pergrass <- model42.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model42.diff$obs_diff[model42.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model42.pergrass

#   Shrub
model42.shrub <- model42.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model42.diff$obs_diff[model42.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model42.shrub

#   Shannon diversity
model42.shannon <- model42.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model42.diff$obs_diff[model42.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity (**)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model42.shannon

#   Bromus rubens
model42.brru2 <- model42.perm |> 
  filter(indicators == "Bromus rubens") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model42.diff$obs_diff[model42.diff$indicators == "Bromus rubens"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus rubens") ~ "(*)")) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model42.brru2

#   Pinus and Juniperus
model42.pj <- model42.perm |> 
  filter(indicators == "Pinus and Juniperus") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model42.diff$obs_diff[model42.diff$indicators == "Pinus and Juniperus"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Pinus") ~ "&" ~ italic("Juniperus"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model42.pj

# Combine plots
grid.arrange(
  model42.bp, model42.annforb, model42.anngrass,
  model42.perforb, model42.pergrass, model42.shrub,
  model42.shannon, model42.brru2, model42.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(NA, 8, 8, 9, 9, NA)
  )
)




# Write out figures -------------------------------------------------------

## Blue Mountains ---------------------------------------------------------

# 1. Blue Mountains: Herbicide
tiff("figures/2026-06_PSM-and-permutation-tests/model01_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model01.bp, model01.annforb, model01.anngrass,
  model01.perforb, model01.pergrass, model01.shrub,
  model01.shannon, model01.brte, model01.artemisia, 
  model01.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)
dev.off()


# 2. Blue Mountains: Post-burn herbicide
tiff("figures/2026-06_PSM-and-permutation-tests/model02_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model02.bp, model02.annforb, model02.anngrass,
  model02.perforb, model02.pergrass, model02.shrub,
  model02.shannon, model02.brte, model02.artemisia, 
  model02.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)
dev.off()


## Middle Rockies ---------------------------------------------------------

# 3. Middle Rockies: Herbicide
tiff("figures/2026-06_PSM-and-permutation-tests/model03_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model03.bp, model03.annforb, model03.anngrass,
  model03.perforb, model03.pergrass, model03.shrub,
  model03.shannon, model03.brte, model03.artemisia, 
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(NA, 8, 8, 9, 9, NA)
  )
)
dev.off()


## Southern Rockies -------------------------------------------------------

# 4. Southern Rockies: Herbicide
tiff("figures/2026-06_PSM-and-permutation-tests/model04_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model04.bp, model04.annforb, model04.anngrass,
  model04.perforb, model04.pergrass, model04.shrub,
  model04.shannon, model04.artemisia, model04.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(NA, 8, 8, 9, 9, NA)
  )
)
dev.off()

# 5. Southern Rockies: Prescribed burn
tiff("figures/2026-06_PSM-and-permutation-tests/model05_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model05.bp, model05.annforb, model05.anngrass,
  model05.perforb, model05.pergrass, model05.shrub,
  model05.shannon, model05.artemisia, model05.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(NA, 8, 8, 9, 9, NA)
  )
)
dev.off()



## Northwestern Great Plains ----------------------------------------------

# 6. NW Great Plains: Prescribed burn
tiff("figures/2026-06_PSM-and-permutation-tests/model06_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model06.bp, model06.annforb, model06.anngrass,
  model06.perforb, model06.pergrass, model06.shrub,
  model06.shannon, model06.brte, model06.artemisia, 
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(NA, 8, 8, 9, 9, NA)
  )
)
dev.off()


## Snake River Plain ------------------------------------------------------

# 7. Snake River Plain: Post-burn aerial seeding
tiff("figures/2026-06_PSM-and-permutation-tests/model07_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model07.bp, model07.annforb, model07.anngrass,
  model07.perforb, model07.pergrass, model07.shrub,
  model07.shannon, model07.brte, model07.artemisia,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(NA, 8, 8, 9, 9, NA)
  )
)
dev.off()

# 8. Snake River Plain: Post-burn aerial & drill seeding
tiff("figures/2026-06_PSM-and-permutation-tests/model08_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model08.bp, model08.annforb, model08.anngrass,
  model08.perforb, model08.pergrass, model08.shrub,
  model08.shannon, model08.brte, model08.artemisia,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(NA, 8, 8, 9, 9, NA)
  )
)
dev.off()

# 9. Snake River Plain: Post-burn closure
tiff("figures/2026-06_PSM-and-permutation-tests/model09_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model09.bp, model09.annforb, model09.anngrass,
  model09.perforb, model09.pergrass, model09.shrub,
  model09.shannon, model09.brte, model09.artemisia,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(NA, 8, 8, 9, 9, NA)
  )
)
dev.off()

# 10. Snake River Plain: Post-burn drill seeding
tiff("figures/2026-06_PSM-and-permutation-tests/model10_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model10.bp, model10.annforb, model10.anngrass,
  model10.perforb, model10.pergrass, model10.shrub,
  model10.shannon, model10.brte, model10.artemisia,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(NA, 8, 8, 9, 9, NA)
  )
)
dev.off()

# 11. Snake River Plain: Post-burn herbicide
tiff("figures/2026-06_PSM-and-permutation-tests/model11_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model11.bp, model11.annforb, model11.anngrass,
  model11.perforb, model11.pergrass, model11.shrub,
  model11.shannon, model11.brte, model11.artemisia,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(NA, 8, 8, 9, 9, NA)
  )
)
dev.off()


## Northern Basin and Range -----------------------------------------------

# 12. Northern BR: Drill seeding
tiff("figures/2026-06_PSM-and-permutation-tests/model12_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model12.bp, model12.annforb, model12.anngrass,
  model12.perforb, model12.pergrass, model12.shrub,
  model12.shannon, model12.brte, model12.artemisia, 
  model12.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)
dev.off()

# 13. Northern BR: Drill seeding & soil disturbance
tiff("figures/2026-06_PSM-and-permutation-tests/model13_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model13.bp, model13.annforb, model13.anngrass,
  model13.perforb, model13.pergrass, model13.shrub,
  model13.shannon, model13.brte, model13.artemisia, 
  model13.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)
dev.off()

# 14. Northern BR: Herbicide
tiff("figures/2026-06_PSM-and-permutation-tests/model14_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model14.bp, model14.annforb, model14.anngrass,
  model14.perforb, model14.pergrass, model14.shrub,
  model14.shannon, model14.brte, model14.artemisia, 
  model14.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)
dev.off()

# 15. Northern BR: Prescribed burn
tiff("figures/2026-06_PSM-and-permutation-tests/model15_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model15.bp, model15.annforb, model15.anngrass,
  model15.perforb, model15.pergrass, model15.shrub,
  model15.shannon, model15.brte, model15.artemisia, 
  model15.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)
dev.off()

# 16. Northern BR: Vegetation disturbance
tiff("figures/2026-06_PSM-and-permutation-tests/model16_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model16.bp, model16.annforb, model16.anngrass,
  model16.perforb, model16.pergrass, model16.shrub,
  model16.shannon, model16.brte, model16.artemisia, 
  model16.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)
dev.off()

# 17. Northern BR: Post-burn aerial seeding
tiff("figures/2026-06_PSM-and-permutation-tests/model17_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model17.bp, model17.annforb, model17.anngrass,
  model17.perforb, model17.pergrass, model17.shrub,
  model17.shannon, model17.brte, model17.artemisia, 
  model17.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)
dev.off()

# 18. Northern BR: Post-burn aerial and drill seeding
tiff("figures/2026-06_PSM-and-permutation-tests/model18_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model18.bp, model18.annforb, model18.anngrass,
  model18.perforb, model18.pergrass, model18.shrub,
  model18.shannon, model18.brte, model18.artemisia, 
  model18.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)
dev.off()

# 19. Northern BR: Post-burn closure
tiff("figures/2026-06_PSM-and-permutation-tests/model19_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model19.bp, model19.annforb, model19.anngrass,
  model19.perforb, model19.pergrass, model19.shrub,
  model19.shannon, model19.brte, model19.artemisia, 
  model19.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)
dev.off()

# 20. Northern BR: Post-burn drill seeding
tiff("figures/2026-06_PSM-and-permutation-tests/model20_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model20.bp, model20.annforb, model20.anngrass,
  model20.perforb, model20.pergrass, model20.shrub,
  model20.shannon, model20.brte, model20.artemisia, 
  model20.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)
dev.off()

# 21. Northern BR: Post-burn herbicide
tiff("figures/2026-06_PSM-and-permutation-tests/model21_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model21.bp, model21.annforb, model21.anngrass,
  model21.perforb, model21.pergrass, model21.shrub,
  model21.shannon, model21.brte, model21.artemisia, 
  model21.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)
dev.off()

# 22. Northern BR: Post-burn seedling planting
tiff("figures/2026-06_PSM-and-permutation-tests/model22_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model22.bp, model22.annforb, model22.anngrass,
  model22.perforb, model22.pergrass, model22.shrub,
  model22.shannon, model22.brte, model22.artemisia, 
  model22.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)
dev.off()


## Central Basin and Range ------------------------------------------------

# 23. Central BR: Drill seeding & soil disturbance
tiff("figures/2026-06_PSM-and-permutation-tests/model23_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model23.bp, model23.annforb, model23.anngrass,
  model23.perforb, model23.pergrass, model23.shrub,
  model23.shannon, model23.brte, model23.artemisia, 
  model23.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)
dev.off()

# 24. Central BR: Prescribed burn
tiff("figures/2026-06_PSM-and-permutation-tests/model24_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model24.bp, model24.annforb, model24.anngrass,
  model24.perforb, model24.pergrass, model24.shrub,
  model24.shannon, model24.brte, model24.artemisia, 
  model24.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)
dev.off()

# 25. Central BR: Vegetation disturbance
tiff("figures/2026-06_PSM-and-permutation-tests/model25_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model25.bp, model25.annforb, model25.anngrass,
  model25.perforb, model25.pergrass, model25.shrub,
  model25.shannon, model25.brte, model25.artemisia, 
  model25.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)
dev.off()

# 26. Central BR: Post-burn aerial seeding
tiff("figures/2026-06_PSM-and-permutation-tests/model26_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model26.bp, model26.annforb, model26.anngrass,
  model26.perforb, model26.pergrass, model26.shrub,
  model26.shannon, model26.brte, model26.artemisia, 
  model26.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)
dev.off()

# 27. Central BR: Post-burn drill seeding
tiff("figures/2026-06_PSM-and-permutation-tests/model27_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model27.bp, model27.annforb, model27.anngrass,
  model27.perforb, model27.pergrass, model27.shrub,
  model27.shannon, model27.brte, model27.artemisia, 
  model27.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)
dev.off()

# 28. Central BR: Post-burn ground seeding
tiff("figures/2026-06_PSM-and-permutation-tests/model28_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model28.bp, model28.annforb, model28.anngrass,
  model28.perforb, model28.pergrass, model28.shrub,
  model28.shannon, model28.brte, model28.artemisia, 
  model28.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)
dev.off()

# 29. Central BR: Post-burn herbicide
tiff("figures/2026-06_PSM-and-permutation-tests/model29_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model29.bp, model29.annforb, model29.anngrass,
  model29.perforb, model29.pergrass, model29.shrub,
  model29.shannon, model29.brte, model29.artemisia, 
  model29.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)
dev.off()


## Wyoming Basin ----------------------------------------------------------

# 30. Wyoming Basin: Prescribed burn
tiff("figures/2026-06_PSM-and-permutation-tests/model30_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model30.bp, model30.annforb, model30.anngrass,
  model30.perforb, model30.pergrass, model30.shrub,
  model30.shannon, model30.brte, model30.artemisia,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(NA, 8, 8, 9, 9, NA)
  )
)
dev.off()


## Colorado Plateaus ------------------------------------------------------

# 31. CO Plateaus: Aerial seeding & soil disturbance
tiff("figures/2026-06_PSM-and-permutation-tests/model31_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model31.bp, model31.annforb, model31.anngrass,
  model31.perforb, model31.pergrass, model31.shrub,
  model31.shannon, model31.brte, model31.artemisia, 
  model31.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)
dev.off()

# 32. CO Plateaus: Herbicide
tiff("figures/2026-06_PSM-and-permutation-tests/model32_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model32.bp, model32.annforb, model32.anngrass,
  model32.perforb, model32.pergrass, model32.shrub,
  model32.shannon, model32.brte, model32.artemisia, 
  model32.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)
dev.off()

# 33. CO Plateaus: Prescribed burn
tiff("figures/2026-06_PSM-and-permutation-tests/model33_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model33.bp, model33.annforb, model33.anngrass,
  model33.perforb, model33.pergrass, model33.shrub,
  model33.shannon, model33.brte, model33.artemisia, 
  model33.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)
dev.off()

# 34. CO Plateaus: Soil disturbance
tiff("figures/2026-06_PSM-and-permutation-tests/model34_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model34.bp, model34.annforb, model34.anngrass,
  model34.perforb, model34.pergrass, model34.shrub,
  model34.shannon, model34.brte, model34.artemisia, 
  model34.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)
dev.off()

# 35. CO Plateaus: Vegetation disturbance
tiff("figures/2026-06_PSM-and-permutation-tests/model35_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model35.bp, model35.annforb, model35.anngrass,
  model35.perforb, model35.pergrass, model35.shrub,
  model35.shannon, model35.brte, model35.artemisia, 
  model35.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)
dev.off()

# 36. CO Plateaus: Post-burn aerial seeding
tiff("figures/2026-06_PSM-and-permutation-tests/model36_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model36.bp, model36.annforb, model36.anngrass,
  model36.perforb, model36.pergrass, model36.shrub,
  model36.shannon, model36.brte, model36.artemisia, 
  model36.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)
dev.off()


## AZ/NM Plateau ----------------------------------------------------------

# 37. AZ/NM Plateau: Herbicide
tiff("figures/2026-06_PSM-and-permutation-tests/model37_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model37.bp, model37.annforb, model37.anngrass,
  model37.perforb, model37.pergrass, model37.shrub,
  model37.shannon, model37.brte, model37.artemisia, 
  model37.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)
dev.off()

# 38. AZ/NM Plateau: Prescribed burn
tiff("figures/2026-06_PSM-and-permutation-tests/model38_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model41.bp, model41.annforb, model41.anngrass,
  model41.perforb, model41.pergrass, model41.shrub,
  model41.shannon, model41.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()

# 39. AZ/NM Plateau: Soil disturbance
tiff("figures/2026-06_PSM-and-permutation-tests/model39_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model38.bp, model38.annforb, model38.anngrass,
  model38.perforb, model38.pergrass, model38.shrub,
  model38.shannon, model38.brte, model38.artemisia, 
  model38.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(8, 8, 9, 9, 10, 10)
  )
)
dev.off()


## Mojave Basin and Range -------------------------------------------------

# 40. Mojave BR: Post-burn aerial seeding
tiff("figures/2026-06_PSM-and-permutation-tests/model40_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model40.bp, model40.annforb, model40.anngrass,
  model40.perforb, model40.pergrass, model40.shrub,
  model40.shannon, model40.brru2, model40.latr2,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(NA, 8, 8, 9, 9, NA)
  )
)
dev.off()


## Chihuahuan Desert ------------------------------------------------------

# 41. Chihuahuan Desert: Herbicide
tiff("figures/2026-06_PSM-and-permutation-tests/model41_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model41.bp, model41.annforb, model41.anngrass,
  model41.perforb, model41.pergrass, model41.shrub,
  model41.shannon, model41.prosopis, model41.latr2,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(NA, 8, 8, 9, 9, NA)
  )
)
dev.off()


## AZ/NM Mountains --------------------------------------------------------

# 42. AZ/NM Mountains: Prescribed burn
tiff("figures/2026-06_PSM-and-permutation-tests/model42_permutation.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model42.bp, model42.annforb, model42.anngrass,
  model42.perforb, model42.pergrass, model42.shrub,
  model42.shannon, model42.brru2, model42.pj,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(2, 2, 3, 3, 4, 4),
    c(5, 5, 6, 6, 7, 7),
    c(NA, 8, 8, 9, 9, NA)
  )
)
dev.off()


save.image("RData/20_permutation-tests.RData")

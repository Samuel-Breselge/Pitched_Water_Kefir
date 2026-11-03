#### processing of biolector results #### 
rm(list=ls())
library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(here)

df_merged = read.csv(here("data", "Data_BioLector.csv"), check.names = FALSE)

# remove outliers 
# E06 # outlier in biomass 
df_merged %>%
  filter(combination == "B.acola S.cer L.hil G.oxy") %>%
  slice_max(Biomass, n = 1)

# A07 # outlier in biomass
df_merged %>%
  filter(combination == "B.aqui S.cer L.sat Z.mob G.jap") %>%
  slice_max(Biomass, n = 1)

df_merged <- df_merged %>%
  filter(!Sample %in% c("A07", "D02"))


#### dot plots AAB #### 
#### plots pH #### 
df_subset <- df_merged %>%
  filter(combination %in% c("B.fer S.cer L.hil", "B.fer S.cer L.hil G.oxy"))

# Calculate mean ± SD pH per Time[h] and combination
df_summary <- df_subset %>%
  group_by(`Time[h]`, combination) %>%
  summarise(
    mean_pH = mean(pH, na.rm = TRUE),
    sd_pH   = sd(pH, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(sd_pH = ifelse(is.na(sd_pH), 0, sd_pH))

df_pH_below_4.4 <- df_summary %>%
  filter(mean_pH < 4.4) %>%
  group_by(combination) %>%
  slice_min(`Time[h]`, with_ties = FALSE)

df_pH_below_4.4

# Ribbon + line
p1_AAB_pH <- ggplot() +
  geom_ribbon(
    data = df_summary,
    aes(x = `Time[h]`,
        ymin = mean_pH - sd_pH,
        ymax = mean_pH + sd_pH,
        fill = combination),
    alpha = 0.30,
    colour = NA
  ) +
  geom_line(
    data = df_summary,
    aes(x = `Time[h]`, y = mean_pH, colour = combination),
    linewidth = 1
  ) +
  labs(x = "Time [h]", y = "pH", title = "a") +
  theme_bw()

p1_AAB_pH

df_subset <- df_merged %>%
  filter(combination %in% c("B.acola S.cer L.sat Z.mob", "B.acola S.cer L.sat Z.mob G.jap"))

# Calculate mean ± SD pH per Time[h] and combination
df_summary <- df_subset %>%
  group_by(`Time[h]`, combination) %>%
  summarise(
    mean_pH = mean(pH, na.rm = TRUE),
    sd_pH   = sd(pH, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(sd_pH = ifelse(is.na(sd_pH), 0, sd_pH))

df_pH_below_4.4 <- df_summary %>%
  filter(mean_pH < 4.4) %>%
  group_by(combination) %>%
  slice_min(`Time[h]`, with_ties = FALSE)

df_pH_below_4.4

p2_AAB_pH <- ggplot() +
  geom_ribbon(
    data = df_summary,
    aes(x = `Time[h]`,
        ymin = mean_pH - sd_pH,
        ymax = mean_pH + sd_pH,
        fill = combination),
    alpha = 0.30,
    colour = NA
  ) +
  geom_line(
    data = df_summary,
    aes(x = `Time[h]`, y = mean_pH, colour = combination),
    linewidth = 1
  ) +
  labs(x = "Time [h]", y = "pH", title = "b") +
  theme_bw()

p2_AAB_pH


#### plots Biomass #### 
df_subset <- df_merged %>%
  filter(combination %in% c("B.fer S.cer L.hil", "B.fer S.cer L.hil G.oxy"))

df_summary <- df_subset %>%
  group_by(`Time[h]`, combination) %>%
  summarise(
    mean_Biomass = mean(Biomass, na.rm = TRUE),
    sd_Biomass   = sd(Biomass, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(sd_Biomass = ifelse(is.na(sd_Biomass), 0, sd_Biomass))

p1_AAB_biomass <- ggplot() +
  geom_ribbon(
    data = df_summary,
    aes(x = `Time[h]`,
        ymin = mean_Biomass - sd_Biomass,
        ymax = mean_Biomass + sd_Biomass,
        fill = combination),
    alpha = 0.30,
    colour = NA
  ) +
  geom_line(
    data = df_summary,
    aes(x = `Time[h]`, y = mean_Biomass, colour = combination),
    linewidth = 1
  ) +
  labs(x = "Time [h]", y = "Biomass (RU)", title = "c") +
  ylim(200, 1400) +
  theme_bw()

p1_AAB_biomass

df_subset <- df_merged %>%
  filter(combination %in% c("B.acola S.cer L.sat Z.mob", "B.acola S.cer L.sat Z.mob G.jap"))

df_summary <- df_subset %>%
  group_by(`Time[h]`, combination) %>%
  summarise(
    mean_Biomass = mean(Biomass, na.rm = TRUE),
    sd_Biomass   = sd(Biomass, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(sd_Biomass = ifelse(is.na(sd_Biomass), 0, sd_Biomass))

p2_AAB_biomass <- ggplot() +
  geom_ribbon(
    data = df_summary,
    aes(x = `Time[h]`,
        ymin = mean_Biomass - sd_Biomass,
        ymax = mean_Biomass + sd_Biomass,
        fill = combination),
    alpha = 0.30,
    colour = NA
  ) +
  geom_line(
    data = df_summary,
    aes(x = `Time[h]`, y = mean_Biomass, colour = combination),
    linewidth = 1
  ) +
  labs(x = "Time [h]", y = "Biomass (RU)", title = "d") +
  ylim(200, 1400) +
  theme_bw()

p2_AAB_biomass


#### plots DO #### 
df_subset <- df_merged %>%
  filter(combination %in% c("B.fer S.cer L.hil", "B.fer S.cer L.hil G.oxy"))

df_summary <- df_subset %>%
  group_by(`Time[h]`, combination) %>%
  summarise(
    mean_DO = mean(DO, na.rm = TRUE),
    sd_DO   = sd(DO, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(sd_DO = ifelse(is.na(sd_DO), 0, sd_DO))

p1_AAB_DO <- ggplot() +
  geom_ribbon(
    data = df_summary,
    aes(x = `Time[h]`,
        ymin = mean_DO - sd_DO,
        ymax = mean_DO + sd_DO,
        fill = combination),
    alpha = 0.30,
    colour = NA
  ) +
  geom_line(
    data = df_summary,
    aes(x = `Time[h]`, y = mean_DO, colour = combination),
    linewidth = 1
  ) +
  labs(x = "Time [h]", y = "Dissolved oxygen (%)", title = "e") +
  theme_bw()

p1_AAB_DO

df_subset <- df_merged %>%
  filter(combination %in% c("B.acola S.cer L.sat Z.mob", "B.acola S.cer L.sat Z.mob G.jap"))

df_summary <- df_subset %>%
  group_by(`Time[h]`, combination) %>%
  summarise(
    mean_DO = mean(DO, na.rm = TRUE),
    sd_DO   = sd(DO, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(sd_DO = ifelse(is.na(sd_DO), 0, sd_DO))

p2_AAB_DO <- ggplot() +
  geom_ribbon(
    data = df_summary,
    aes(x = `Time[h]`,
        ymin = mean_DO - sd_DO,
        ymax = mean_DO + sd_DO,
        fill = combination),
    alpha = 0.30,
    colour = NA
  ) +
  geom_line(
    data = df_summary,
    aes(x = `Time[h]`, y = mean_DO, colour = combination),
    linewidth = 1
  ) +
  labs(x = "Time [h]", y = "Dissolved oxygen (%)", title = "f") +
  theme_bw()

p2_AAB_DO

# Arrange the plots
combined_plot <- (p1_AAB_pH | p2_AAB_pH) / 
  (p1_AAB_biomass | p2_AAB_biomass) / 
  (p1_AAB_DO | p2_AAB_DO)

# Save to PDF
ggsave(
  filename = here("output", "Figure9.pdf"), 
  combined_plot, width = 12, height = 9)



#### dot plots BIF #### 
#### plots pH #### 
df_subset <- df_merged %>%
  filter(combination %in% c(
    "B.fer S.cer L.hil G.oxy",
    "B.acola S.cer L.hil G.oxy",
    "B.aqui S.cer L.hil G.oxy",
    "S.cer L.hil G.oxy"
  ))

df_summary <- df_subset %>%
  group_by(`Time[h]`, combination) %>%
  summarise(
    mean_pH = mean(pH, na.rm = TRUE),
    sd_pH   = sd(pH, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(sd_pH = ifelse(is.na(sd_pH), 0, sd_pH))

p1_Bif_pH <- ggplot() +
  geom_ribbon(
    data = df_summary,
    aes(x = `Time[h]`,
        ymin = mean_pH - sd_pH,
        ymax = mean_pH + sd_pH,
        fill = combination),
    alpha = 0.30,
    colour = NA
  ) +
  geom_line(
    data = df_summary,
    aes(x = `Time[h]`, y = mean_pH, colour = combination),
    linewidth = 1
  ) +
  labs(x = "Time [h]", y = "pH", title = "a") +
  xlim(0, 48) +
  theme_bw()

p1_Bif_pH

df_subset <- df_merged %>%
  filter(combination %in% c(
    "B.fer S.cer L.sat Z.mob G.jap",
    "B.acola S.cer L.sat Z.mob G.jap",
    "B.aqui S.cer L.sat Z.mob G.jap",
    "S.cer L.sat Z.mob G.jap"
  ))

df_summary <- df_subset %>%
  group_by(`Time[h]`, combination) %>%
  summarise(
    mean_pH = mean(pH, na.rm = TRUE),
    sd_pH   = sd(pH, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(sd_pH = ifelse(is.na(sd_pH), 0, sd_pH))

p2_Bif_pH <- ggplot() +
  geom_ribbon(
    data = df_summary,
    aes(x = `Time[h]`,
        ymin = mean_pH - sd_pH,
        ymax = mean_pH + sd_pH,
        fill = combination),
    alpha = 0.30,
    colour = NA
  ) +
  geom_line(
    data = df_summary,
    aes(x = `Time[h]`, y = mean_pH, colour = combination),
    linewidth = 1
  ) +
  labs(x = "Time [h]", y = "pH", title = "b") +
  xlim(0, 48) +
  theme_bw()

p2_Bif_pH


#### plots Biomass #### 
df_subset <- df_merged %>%
  filter(combination %in% c(
    "B.fer S.cer L.hil G.oxy",
    "B.acola S.cer L.hil G.oxy",
    "B.aqui S.cer L.hil G.oxy",
    "S.cer L.hil G.oxy"
  ))

df_summary <- df_subset %>%
  group_by(`Time[h]`, combination) %>%
  summarise(
    mean_Biomass = mean(Biomass, na.rm = TRUE),
    sd_Biomass   = sd(Biomass, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(sd_Biomass = ifelse(is.na(sd_Biomass), 0, sd_Biomass))

p1_Bif_biomass <- ggplot() +
  geom_ribbon(
    data = df_summary,
    aes(x = `Time[h]`,
        ymin = mean_Biomass - sd_Biomass,
        ymax = mean_Biomass + sd_Biomass,
        fill = combination),
    alpha = 0.30,
    colour = NA
  ) +
  geom_line(
    data = df_summary,
    aes(x = `Time[h]`, y = mean_Biomass, colour = combination),
    linewidth = 1
  ) +
  labs(x = "Time [h]", y = "Biomass (RU)", title = "c") +
  ylim(200, 1400) +
  theme_bw()

p1_Bif_biomass

df_subset <- df_merged %>%
  filter(combination %in% c(
    "B.fer S.cer L.sat Z.mob G.jap",
    "B.acola S.cer L.sat Z.mob G.jap",
    "B.aqui S.cer L.sat Z.mob G.jap",
    "S.cer L.sat Z.mob G.jap"
  ))

df_summary <- df_subset %>%
  group_by(`Time[h]`, combination) %>%
  summarise(
    mean_Biomass = mean(Biomass, na.rm = TRUE),
    sd_Biomass   = sd(Biomass, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(sd_Biomass = ifelse(is.na(sd_Biomass), 0, sd_Biomass))

p2_Bif_biomass <- ggplot() +
  geom_ribbon(
    data = df_summary,
    aes(x = `Time[h]`,
        ymin = mean_Biomass - sd_Biomass,
        ymax = mean_Biomass + sd_Biomass,
        fill = combination),
    alpha = 0.30,
    colour = NA
  ) +
  geom_line(
    data = df_summary,
    aes(x = `Time[h]`, y = mean_Biomass, colour = combination),
    linewidth = 1
  ) +
  labs(x = "Time [h]", y = "Biomass (RU)", title = "d") +
  ylim(200, 1400) +
  theme_bw()

p2_Bif_biomass


#### plots DO #### 
df_subset <- df_merged %>%
  filter(combination %in% c(
    "B.fer S.cer L.hil G.oxy",
    "B.acola S.cer L.hil G.oxy",
    "B.aqui S.cer L.hil G.oxy",
    "S.cer L.hil G.oxy"
  ))

df_summary <- df_subset %>%
  group_by(`Time[h]`, combination) %>%
  summarise(
    mean_DO = mean(DO, na.rm = TRUE),
    sd_DO   = sd(DO, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(sd_DO = ifelse(is.na(sd_DO), 0, sd_DO))

p1_Bif_DO <- ggplot() +
  geom_ribbon(
    data = df_summary,
    aes(x = `Time[h]`,
        ymin = mean_DO - sd_DO,
        ymax = mean_DO + sd_DO,
        fill = combination),
    alpha = 0.30,
    colour = NA
  ) +
  geom_line(
    data = df_summary,
    aes(x = `Time[h]`, y = mean_DO, colour = combination),
    linewidth = 1
  ) +
  labs(x = "Time [h]", y = "Dissolved oxygen (%)", title = "e") +
  theme_bw()

p1_Bif_DO

df_subset <- df_merged %>%
  filter(combination %in% c(
    "B.fer S.cer L.sat Z.mob G.jap",
    "B.acola S.cer L.sat Z.mob G.jap",
    "B.aqui S.cer L.sat Z.mob G.jap",
    "S.cer L.sat Z.mob G.jap"
  ))

df_summary <- df_subset %>%
  group_by(`Time[h]`, combination) %>%
  summarise(
    mean_DO = mean(DO, na.rm = TRUE),
    sd_DO   = sd(DO, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(sd_DO = ifelse(is.na(sd_DO), 0, sd_DO))

p2_Bif_DO <- ggplot() +
  geom_ribbon(
    data = df_summary,
    aes(x = `Time[h]`,
        ymin = mean_DO - sd_DO,
        ymax = mean_DO + sd_DO,
        fill = combination),
    alpha = 0.30,
    colour = NA
  ) +
  geom_line(
    data = df_summary,
    aes(x = `Time[h]`, y = mean_DO, colour = combination),
    linewidth = 1
  ) +
  labs(x = "Time [h]", y = "Dissolved oxygen (%)", title = "f") +
  theme_bw()

p2_Bif_DO

# Arrange the plots
combined_plot <- (p1_Bif_pH | p2_Bif_pH) / 
  (p1_Bif_biomass | p2_Bif_biomass) / 
  (p1_Bif_DO | p2_Bif_DO)

# Save to PDF
ggsave(
  filename = here("output", "FigureS9 - Bifs BioLector.pdf"), 
  combined_plot, width = 12, height = 9)

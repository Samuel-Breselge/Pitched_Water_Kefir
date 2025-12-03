#### processing of biolector results #### 
rm(list=ls())
library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(here)

df_merged = read.csv(here("data", "Data_BioLector.csv"), check.names = F)

#remove outliers 
#E06 #outlyier in biomass 
df_merged %>%
  filter(combination == "B.acola S.cer L.hil G.oxy") %>%
  slice_max(Biomass, n = 1)
#A07 #outlyier in biomass
df_merged %>%
  filter(combination == "B.aqui S.cer L.sat Z.mob G.jap") %>%
  slice_max(Biomass, n = 1)

df_merged <- df_merged %>%
  filter(!Sample %in% c("A07", "D02"))


#### dot plots AAB #### 
#### plots pH #### 
df_subset <- df_merged %>%
  filter(combination %in% c("B.fer S.cer L.hil", "B.fer S.cer L.hil G.oxy"))

# Calculate mean pH per Time[h] and combination
df_summary <- df_subset %>%
  group_by(`Time[h]`, combination) %>%
  summarise(mean_pH = mean(pH, na.rm = TRUE), .groups = "drop")

df_pH_below_4.4 <- df_summary %>%
  filter(mean_pH < 4.4) %>%
  group_by(combination) %>%
  slice_min(`Time[h]`, with_ties = FALSE)

df_pH_below_4.4

# Line plot with individual dots
p1_AAB_pH <- ggplot() +
  # Individual points
  geom_jitter(data = df_subset, aes(x = `Time[h]`, y = pH, color = combination), width = 0.4, alpha = 0.3) +
  # Mean line
  geom_line(data = df_summary, aes(x = `Time[h]`, y = mean_pH, group = combination), 
            color = "black", linewidth = 1) +
  labs(x = "Time [h]", y = "pH", title = "a") +
  #ylim(4, 6) + #limits of BioLector plate
  theme_bw()

p1_AAB_pH

df_subset <- df_merged %>%
  filter(combination %in% c("B.acola S.cer L.sat Z.mob", "B.acola S.cer L.sat Z.mob G.jap"))

# Calculate mean pH per Time[h] and combination
df_summary <- df_subset %>%
  group_by(`Time[h]`, combination) %>%
  summarise(mean_pH = mean(pH, na.rm = TRUE), .groups = "drop")

df_pH_below_4.4 <- df_summary %>%
  filter(mean_pH < 4.4) %>%
  group_by(combination) %>%
  slice_min(`Time[h]`, with_ties = FALSE)

df_pH_below_4.4

# Line plot with individual dots
p2_AAB_pH <- ggplot() +
  # Individual points
  geom_jitter(data = df_subset, aes(x = `Time[h]`, y = pH, color = combination), width = 0.4, alpha = 0.3) +
  # Mean line
  geom_line(data = df_summary, aes(x = `Time[h]`, y = mean_pH, group = combination), 
            color = "black", linewidth = 1) +
  labs(x = "Time [h]", y = "pH", title = "b") +
  #ylim(4, 6) + #limits of BioLector plate
  theme_bw()

p2_AAB_pH

#### plots Biomass #### 
df_subset <- df_merged %>%
  filter(combination %in% c("B.fer S.cer L.hil", "B.fer S.cer L.hil G.oxy"))

# Calculate mean biomass per Time[h] and combination
df_summary <- df_subset %>%
  group_by(`Time[h]`, combination) %>%
  summarise(mean_Biomass = mean(Biomass, na.rm = TRUE), .groups = "drop")

# Line plot with individual dots
p1_AAB_biomass <- ggplot() +
  # Individual points
  geom_jitter(data = df_subset, aes(x = `Time[h]`, y = Biomass, color = combination), width = 0.4, alpha = 0.3) +
  # Mean line
  geom_line(data = df_summary, aes(x = `Time[h]`, y = mean_Biomass, group = combination), 
            color = "black", linewidth = 1) +
  labs(x = "Time [h]", y = "Biomass in RU", title = "c") +
  ylim(200, 1400) + 
  theme_bw()

p1_AAB_biomass

df_subset <- df_merged %>%
  filter(combination %in% c("B.acola S.cer L.sat Z.mob", "B.acola S.cer L.sat Z.mob G.jap"))

# Calculate mean biomass per Time[h] and combination
df_summary <- df_subset %>%
  group_by(`Time[h]`, combination) %>%
  summarise(mean_Biomass = mean(Biomass, na.rm = TRUE), .groups = "drop")

# Line plot with individual dots
p2_AAB_biomass <- ggplot() +
  # Individual points
  geom_jitter(data = df_subset, aes(x = `Time[h]`, y = Biomass, color = combination), width = 0.4, alpha = 0.3) +
  # Mean line
  geom_line(data = df_summary, aes(x = `Time[h]`, y = mean_Biomass, group = combination), 
            color = "black", linewidth = 1) +
  labs(x = "Time [h]", y = "Biomass in RU", title = "d") +
  ylim(200, 1400) + 
  theme_bw()

p2_AAB_biomass


#### plots DO #### 
df_subset <- df_merged %>%
  filter(combination %in% c("B.fer S.cer L.hil", "B.fer S.cer L.hil G.oxy"))

# Calculate mean biomass per Time[h] and combination
df_summary <- df_subset %>%
  group_by(`Time[h]`, combination) %>%
  summarise(mean_DO = mean(DO, na.rm = TRUE), .groups = "drop")

# Line plot with individual dots
p1_AAB_DO <- ggplot() +
  # Individual points
  geom_jitter(data = df_subset, aes(x = `Time[h]`, y = DO, color = combination), width = 0.4, alpha = 0.3) +
  # Mean line
  geom_line(data = df_summary, aes(x = `Time[h]`, y = mean_DO, group = combination), 
            color = "black", linewidth = 1) +
  labs(x = "Time [h]", y = "Dissolved Oxygen", title = "e") +
  #ylim(200, 1400) + 
  theme_bw()

p1_AAB_DO

df_subset <- df_merged %>%
  filter(combination %in% c("B.acola S.cer L.sat Z.mob", "B.acola S.cer L.sat Z.mob G.jap"))

# Calculate mean biomass per Time[h] and combination
df_summary <- df_subset %>%
  group_by(`Time[h]`, combination) %>%
  summarise(mean_DO = mean(DO, na.rm = TRUE), .groups = "drop")

# Line plot with individual dots
p2_AAB_DO <- ggplot() +
  # Individual points
  geom_jitter(data = df_subset, aes(x = `Time[h]`, y = DO, color = combination), width = 0.4, alpha = 0.3) +
  # Mean line
  geom_line(data = df_summary, aes(x = `Time[h]`, y = mean_DO, group = combination), 
            color = "black", linewidth = 1) +
  labs(x = "Time [h]", y = "Dissolved Oxygen", title = "f") +
  #ylim(200, 1400) + 
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
  filter(combination %in% c("B.fer S.cer L.hil G.oxy", "B.acola S.cer L.hil G.oxy", "B.aqui S.cer L.hil G.oxy", "S.cer L.hil G.oxy"))

# Calculate mean pH per Time[h] and combination
df_summary <- df_subset %>%
  group_by(`Time[h]`, combination) %>%
  summarise(mean_pH = mean(pH, na.rm = TRUE), .groups = "drop")

# Line plot with individual dots
p1_Bif_pH <- ggplot() +
  # Individual points
  geom_jitter(data = df_subset, aes(x = `Time[h]`, y = pH, color = combination), width = 0.4, alpha = 0.2, size = 0.2) +
  # Mean line
  geom_line(data = df_summary, aes(x = `Time[h]`, y = mean_pH, color = combination), linewidth = 1) +
  labs(x = "Time [h]", y = "pH", title = "a") +
  #ylim(4, 6) + #limits of BioLector plate
  xlim(0, 48) + 
  theme_bw()

p1_Bif_pH

df_subset <- df_merged %>%
  filter(combination %in% c("B.fer S.cer L.sat Z.mob G.jap", "B.acola S.cer L.sat Z.mob G.jap", "B.aqui S.cer L.sat Z.mob G.jap", "S.cer L.sat Z.mob G.jap"))

# Calculate mean pH per Time[h] and combination
df_summary <- df_subset %>%
  group_by(`Time[h]`, combination) %>%
  summarise(mean_pH = mean(pH, na.rm = TRUE), .groups = "drop")

# Line plot with individual dots
p2_Bif_pH <- ggplot() +
  # Individual points
  geom_jitter(data = df_subset, aes(x = `Time[h]`, y = pH, color = combination), width = 0.4, alpha = 0.2, size = 0.2) +
  # Mean line
  geom_line(data = df_summary, aes(x = `Time[h]`, y = mean_pH, color = combination), linewidth = 1) +
  labs(x = "Time [h]", y = "pH", title = "b") +
  #ylim(4, 6) + #limits of BioLector plate
  xlim(0, 48) + 
  theme_bw()

p2_Bif_pH

#### plots Biomass #### 
df_subset <- df_merged %>%
  filter(combination %in% c("B.fer S.cer L.hil G.oxy", "B.acola S.cer L.hil G.oxy", "B.aqui S.cer L.hil G.oxy", "S.cer L.hil G.oxy"))

# Calculate mean biomass per Time[h] and combination
df_summary <- df_subset %>%
  group_by(`Time[h]`, combination) %>%
  summarise(mean_Biomass = mean(Biomass, na.rm = TRUE), .groups = "drop")

# Line plot with individual dots
p1_Bif_biomass <- ggplot() +
  # Individual points
  geom_jitter(data = df_subset, aes(x = `Time[h]`, y = Biomass, color = combination), width = 0.4, alpha = 0.2, size = 0.2) +
  # Mean line
  geom_line(data = df_summary, aes(x = `Time[h]`, y = mean_Biomass, color = combination), linewidth = 1) +
  labs(x = "Time [h]", y = "Biomass in RU", title = "c") +
  ylim(200, 1400) + 
  theme_bw()

p1_Bif_biomass

df_subset <- df_merged %>%
  filter(combination %in% c("B.fer S.cer L.sat Z.mob G.jap", "B.acola S.cer L.sat Z.mob G.jap", "B.aqui S.cer L.sat Z.mob G.jap", "S.cer L.sat Z.mob G.jap"))

# Calculate mean biomass per Time[h] and combination
df_summary <- df_subset %>%
  group_by(`Time[h]`, combination) %>%
  summarise(mean_Biomass = mean(Biomass, na.rm = TRUE), .groups = "drop")

# Line plot with individual dots
p2_Bif_biomass <- ggplot() +
  # Individual points
  geom_jitter(data = df_subset, aes(x = `Time[h]`, y = Biomass, color = combination), width = 0.4, alpha = 0.2, size = 0.2) +
  # Mean line
  geom_line(data = df_summary, aes(x = `Time[h]`, y = mean_Biomass, color = combination), linewidth = 1) +
  labs(x = "Time [h]", y = "Biomass in RU", title = "d") +
  ylim(200, 1400) + 
  theme_bw()

p2_Bif_biomass


#### plots DO #### 
df_subset <- df_merged %>%
  filter(combination %in% c("B.fer S.cer L.hil G.oxy", "B.acola S.cer L.hil G.oxy", "B.aqui S.cer L.hil G.oxy", "S.cer L.hil G.oxy"))

# Calculate mean biomass per Time[h] and combination
df_summary <- df_subset %>%
  group_by(`Time[h]`, combination) %>%
  summarise(mean_DO = mean(DO, na.rm = TRUE), .groups = "drop")

# Line plot with individual dots
p1_Bif_DO <- ggplot() +
  # Individual points
  geom_jitter(data = df_subset, aes(x = `Time[h]`, y = DO, color = combination), width = 0.4, alpha = 0.2, size = 0.2) +
  # Mean line
 # geom_line(data = df_summary, aes(x = `Time[h]`, y = mean_DO, color = combination), linewidth = 1) +
  labs(x = "Time [h]", y = "Dissolved Oxygen", title = "e") +
  #ylim(200, 1400) + 
  theme_bw()

p1_Bif_DO

df_subset <- df_merged %>%
  filter(combination %in% c("B.fer S.cer L.sat Z.mob G.jap", "B.acola S.cer L.sat Z.mob G.jap", "B.aqui S.cer L.sat Z.mob G.jap", "S.cer L.sat Z.mob G.jap"))

# Calculate mean biomass per Time[h] and combination
df_summary <- df_subset %>%
  group_by(`Time[h]`, combination) %>%
  summarise(mean_DO = mean(DO, na.rm = TRUE), .groups = "drop")

# Line plot with individual dots
p2_Bif_DO <- ggplot() +
  # Individual points
  geom_jitter(data = df_subset, aes(x = `Time[h]`, y = DO, color = combination), width = 0.4, alpha = 0.2, size = 0.2) +
  # Mean line
  geom_line(data = df_summary, aes(x = `Time[h]`, y = mean_DO, color = combination), linewidth = 1) +
  labs(x = "Time [h]", y = "Dissolved Oxygen", title = "f") +
  #ylim(200, 1400) + 
  theme_bw()

p2_Bif_DO

# Arrange the plots
combined_plot <- (p1_Bif_pH | p2_Bif_pH) / 
  (p1_Bif_biomass | p2_Bif_biomass) / 
  (p1_Bif_DO | p2_Bif_DO)

# Save to PDF
ggsave(
  filename = here("output", "FigureS18 - Bifs BioLector.pdf"), 
  combined_plot, width = 12, height = 9)



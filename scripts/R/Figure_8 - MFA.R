rm(list=ls())

library(dplyr)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(here)

#####################################################################
df_VOC_01 <- read.csv(
  here("data", "Data_VOC_in_RU.csv"),
  check.names = FALSE, 
  row.names = 1
)

#assign levels
df_VOC_01$sample_name <- factor(df_VOC_01$sample_name, levels = c(
  "media",
  "B.fer S.cer L.hil G.oxy",
  "B.fer H.val L.par G.jap Z.mob",
  "B.fer P.fer L.hil A.lov Z.mob",
  "B.acola S.cer L.sat G.jap Z.mob",
  "B.acola S.cer L.sat G.jap",
  "WK grains (30g/L)"
))

# Calculate mean for each chemical per sample_name
df_VOC_mean <- df_VOC_01 %>%
  group_by(sample_name) %>%
  summarise(across(1:83, mean, na.rm = TRUE))

names(df_VOC_mean)[names(df_VOC_mean) == "sample_name"] <- "sample"

#remove media control
df_VOC_mean <- df_VOC_mean[df_VOC_mean$sample != "media", ]

#remove columns that are all zero
df_VOC_mean <- df_VOC_mean %>%
  select(sample, where(~ is.numeric(.x) && !all(.x == 0)))

# importing sample data # 
df_sample = read.csv(here("data", "Taste_panel_Matrix.csv"))

# importing metabolites data # 
df_metabolites <- read.csv(here("data", "Data_metabolites_taste_panel.csv"))

# preparing the metabolites data # 
df_filtered <- df_metabolites %>%
  filter(fermenation == "A")

df_filtered <- df_filtered %>%
  select(ID, Sucrose, Glucose, Ethanol, Fructose, Glycerol, Lactic_Acid, Acetic_Acid)

df_metabolites_mean <- df_filtered %>%
  group_by(ID) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

df_metabolites_mean <- df_metabolites_mean %>%
  filter(ID != "fSB370")

df_metabolites_mean <- df_metabolites_mean %>%
  mutate(
    sample = case_when(
      ID == "fSB363" ~ "B.fer S.cer L.hil G.oxy",
      ID == "fSB365" ~ "B.fer H.val L.par G.jap Z.mob",
      ID == "fSB366" ~ "B.acola S.cer L.sat G.jap",
      ID == "fSB367" ~ "B.acola S.cer L.sat G.jap Z.mob",
      ID == "fSB368" ~ "WK grains (30g/L)",
      ID == "fSB369" ~ "B.fer P.fer L.hil A.lov Z.mob"
    )
  )
df_metabolites_mean$ID = NULL

# importing panel data # 
df_panel = read.csv(here("data", "Data_tastepanel.csv"))
df_panel$sample_name = NULL
df_panel$sample_codes_as_entered_by_participants = NULL

df_panel_mean <- df_panel %>%
  group_by(sample) %>%  
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  ungroup()

#### creating MFA data set #### 
# Make sure all dataframes have a clean `sample` column
df_MFA <- df_sample %>%
  left_join(df_VOC_mean, by = "sample") %>%
  left_join(df_panel_mean, by = "sample") %>%
  left_join(df_metabolites_mean, by = "sample")

# Col. 1 = sample names to be used
# Col. 2 to 14 presence/absence table of microbes
# Col. 15 to 95 is VOC data 
# Col. 96 to 100 is hedonic data 
# Col. 101 to 116 is ODP data 
# Col. 117 to 123 is metabolite data 

#running MFA
# Assume Composition block is columns 2 to 13
df_MFA[, 2:14] <- lapply(df_MFA[, 2:14], as.factor)

# Remove the sample name column for MFA but keep row names for plotting
df_MFA_data <- df_MFA[, -1]  # drops column 1
rownames(df_MFA_data) <- df_MFA$sample

# Define block structure
group_sizes <- c(
  13,   # Composition
  81,   # VOCs
  5,    # Hedonic
  16,   # ODP
  7     # Metabolites
)

# Variable types:
# Composition is categorical, the rest are numeric
# So we'll use type = c("n", "s", "s", "s", "s")

# Run MFA
res_mfa <- MFA(
  df_MFA_data,
  group = group_sizes,
  type = c("n", "s", "s", "s", "s"),
  name.group = c("Composition", "VOCs", "Hedonic", "ODP", "Metabolites"),
  graph = TRUE
)

# Scree plot (how many dimensions to keep)
fviz_screeplot(res_mfa)

# Plot samples
fviz_mfa_ind(res_mfa, habillage = "Z.mob", addEllipses = TRUE)


# Plot variables
fviz_mfa_var(res_mfa, "group")  # Groups
fviz_mfa_var(res_mfa, "quanti.var")  # Quantitative vars
fviz_mfa_var(res_mfa, "quanti.var", repel = TRUE)
fviz_mfa_var(res_mfa, "quali.var")   # Qualitative vars (your composition)

# Contribution of groups to dimensions
res_mfa$group$contrib

# Contribution of individual variables
res_mfa$quanti.var$contrib

# Load libraries
library(ggplot2)
library(patchwork)

# Quantitative vars plot (with repel)
p_quanti <- fviz_mfa_var(
  res_mfa,
  "quanti.var",
  repel = TRUE
) +
  ggtitle("(a) Quantitative Variables")

# Qualitative vars: custom filter to show only _1 levels (presence)

# Extract coordinates
quali_var <- as.data.frame(res_mfa$quali.var$coord)
quali_var$Variable <- rownames(quali_var)

# Keep only presence levels
quali_var_1 <- subset(quali_var, grepl("_1$", Variable))

# Clean labels: remove _1 suffix
quali_var_1$Microbe <- sub("_1$", "", quali_var_1$Variable)

# Load ggrepel if not already loaded
library(ggrepel)

# Custom ggplot for filtered qualitative vars with non-overlapping labels
p_quali <- ggplot(quali_var_1, aes(x = Dim.1, y = Dim.2, label = Microbe)) +
  geom_point(color = "steelblue") +
  geom_text_repel(size = 4) +
  theme_minimal() +
  ggtitle("(b) Composition (Presence Only)") +
  xlab("Dim1") +
  ylab("Dim2")
# Combine vertically with patchwork
combined_plot <- p_quanti / p_quali + plot_annotation(
  title = "MFA: Variables"
)

# Show combined plot
combined_plot

# Save if needed
ggsave(
  filename = here("output", "Figure8.pdf"),
  combined_plot,
  width = 15, height = 20, dpi = 300
)
res_mfa$group$contrib

rm(list = ls())

# Libraries ----
# Install ggrepel once if needed:
# install.packages("ggrepel")

library(tidyverse)
library(patchwork)
library(here)
library(readxl)
library(ggrepel)

# Data import ----
df_metabolites <- read_csv(
  here("data", "Data_metabolites.csv")
)

df_01 <- read_excel(
  here("data", "Data_fermentations.xlsx")
)

# Clean fermentation metadata ----
df_01 <- df_01 %>%
  rename(sample = `#`)

# Combine data frames ----
df_combined <- df_metabolites %>%
  left_join(
    df_01 %>% select(sample, average_pannel, `pH 72h`),
    by = "sample",
    relationship = "many-to-one"
  )

# Convert metabolite concentrations from µg/mL to g/L ----
metabolite_columns <- c(
  "Ethanol",
  "Fructose",
  "Glucose",
  "Glycerol",
  "Mannitol",
  "Sucrose",
  "Lactic_Acid",
  "Acetic_Acid"
)

df_combined <- df_combined %>%
  mutate(
    across(all_of(metabolite_columns), ~ .x / 1000)
  )

# Global factor level definitions ----
yeast_levels <- c(
  "media",
  "Br. bruxellensis",
  "H. valbyensis",
  "P. fermentans",
  "S. cerevisiae",
  "WK grains"
)

bif_levels <- c(
  "media",
  "B. aquikefiricola",
  "B. fermentum",
  "WK grains"
)

lab_levels <- c(
  "media",
  "Lac. paracasei",
  "Len. hilgardii",
  "Leu. pseudomesenteroides",
  "Liq. satsumensis",
  "WK grains"
)

aab_levels <- c(
  "media",
  "A. lovaniensis",
  "A. orientalis",
  "G. japonicus",
  "G. oxydans",
  "WK grains"
)

other_levels <- c(
  "media",
  "Z. mobilis",
  "none",
  "WK grains"
)

# Assign factors and readable measurement names ----
df_combined <- df_combined %>%
  mutate(
    Yeast = factor(Yeast, levels = yeast_levels),
    Bif   = factor(Bif, levels = bif_levels),
    LAB   = factor(LAB, levels = lab_levels),
    AAB   = factor(AAB, levels = aab_levels),
    other = factor(other, levels = other_levels),
    
    pH_day3 = `pH 72h`,
    Aroma_score = average_pannel
  )

# Replace missing values representing ND with zero ----
# Keep this step only where NA means "not detected".
df_combined <- df_combined %>%
  mutate(
    across(
      c(Mannitol, Lactic_Acid, Acetic_Acid),
      ~ replace_na(.x, 0)
    )
  )

# Calculate totals: g/L ----
df_combined <- df_combined %>%
  mutate(
    Total_sugars = Fructose + Glucose + Sucrose,
    Total_organic_acids = Lactic_Acid + Acetic_Acid
  )

# Average starting sugars for each experimental set: g/L ----
media_sugars_by_set <- df_combined %>%
  filter(short_name == "media") %>%
  group_by(`experimental set`) %>%
  summarise(
    Initial_sugars = mean(Total_sugars, na.rm = TRUE),
    .groups = "drop"
  )

print(media_sugars_by_set)

# Calculate apparent sugar consumption and yields ----
df_combined <- df_combined %>%
  left_join(
    media_sugars_by_set,
    by = "experimental set",
    relationship = "many-to-one"
  ) %>%
  mutate(
    # Apparent sugar consumption: g/L
    Sugar_consumed = if_else(
      short_name == "media",
      NA_real_,
      Initial_sugars - Total_sugars
    ),
    
    Sugar_consumed_percent =
      100 * Sugar_consumed / Initial_sugars,
    
    # Apparent yields: g/g sugar consumed
    Ethanol_yield = if_else(
      Sugar_consumed > 0,
      Ethanol / Sugar_consumed,
      NA_real_
    ),
    
    Organic_acid_yield = if_else(
      Sugar_consumed > 0,
      Total_organic_acids / Sugar_consumed,
      NA_real_
    ),
    
    # Measured organic acids relative to ethanol: g/g
    Acid_to_ethanol_ratio = if_else(
      short_name != "media" & Ethanol > 0,
      Total_organic_acids / Ethanol,
      NA_real_
    )
  )

# Export combined data ----
write_csv(
  df_combined,
  here("data", "Data_combined.csv")
)

# PCA variables ----
# Use individual measurements, excluding derived totals and yields.
pca_variables <- c(
  metabolite_columns,
  "pH_day3",
  "Aroma_score"
)

# Select defined and grain-based fermentations ----
# Exclude unfermented media.
df_pca_candidates <- df_combined %>%
  filter(
    !is.na(Yeast),
    Yeast != "media"
  )

# Check that PCA variables are numeric ----
numeric_check <- vapply(
  df_pca_candidates[pca_variables],
  is.numeric,
  logical(1)
)

if (!all(numeric_check)) {
  stop(
    "These PCA variables are not numeric: ",
    paste(names(numeric_check)[!numeric_check], collapse = ", ")
  )
}

# Identify incomplete samples ----
df_pca_excluded <- df_pca_candidates %>%
  filter(
    !if_all(all_of(pca_variables), ~ is.finite(.x))
  ) %>%
  select(sample, short_name, all_of(pca_variables))

if (nrow(df_pca_excluded) > 0) {
  message(
    nrow(df_pca_excluded),
    " samples excluded because PCA measurements are missing or non-finite."
  )
  
  print(df_pca_excluded)
}

# Retain complete samples ----
df_pca <- df_pca_candidates %>%
  filter(
    if_all(all_of(pca_variables), ~ is.finite(.x))
  ) %>%
  droplevels()

if (nrow(df_pca) < 3) {
  stop("At least three complete samples are needed for this PCA.")
}

# Inspect sample numbers ----
df_pca %>%
  count(Yeast, name = "Samples_in_PCA") %>%
  print()

# Check variable variance ----
pca_sd <- vapply(
  df_pca[pca_variables],
  sd,
  numeric(1)
)

if (any(!is.finite(pca_sd) | pca_sd == 0)) {
  stop(
    "PCA variables with zero or undefined variance: ",
    paste(
      names(pca_sd)[!is.finite(pca_sd) | pca_sd == 0],
      collapse = ", "
    )
  )
}

# Run PCA ----
# Standardisation gives each measurement equal initial variance.
pca_result <- prcomp(
  df_pca %>% select(all_of(pca_variables)),
  center = TRUE,
  scale. = TRUE
)

# Variance explained ----
variance_explained <- 100 * pca_result$sdev^2 /
  sum(pca_result$sdev^2)

pca_variance <- tibble(
  PC = paste0("PC", seq_along(variance_explained)),
  Variance_percent = variance_explained,
  Cumulative_percent = cumsum(variance_explained)
)

print(pca_variance)

# Combine PCA scores with sample metadata ----
# Different symbols for grain inoculum concentrations ----
df_pca_scores <- bind_cols(
  df_pca,
  as_tibble(pca_result$x)
) %>%
  mutate(
    Grain_symbol = if_else(
      Yeast == "WK grains",
      paste0("grains_", readr::parse_number(as.character(grains))),
      "Defined community"
    )
  )

# Colour palettes ----
# Grain-based fermentations are black in every panel.
make_group_colours <- function(group_column) {
  
  group_levels <- levels(df_pca_scores[[group_column]])
  defined_levels <- setdiff(group_levels, "WK grains")
  
  colours <- setNames(
    scales::hue_pal()(length(defined_levels)),
    defined_levels
  )
  
  if ("WK grains" %in% group_levels) {
    colours <- c(colours, "WK grains" = "black")
  }
  
  colours
}

# PCA score plotting function ----
plot_pca <- function(group_column, legend_title) {
  
  ggplot(
    df_pca_scores,
    aes(
      x = PC1,
      y = PC2,
      colour = .data[[group_column]],
      shape = Grain_symbol
    )
  ) +
    geom_hline(yintercept = 0, colour = "grey85") +
    geom_vline(xintercept = 0, colour = "grey85") +
    stat_ellipse(
      aes(
        x = PC1,
        y = PC2,
        colour = .data[[group_column]],
        group = .data[[group_column]]
      ),
      inherit.aes = FALSE,
      type = "norm",
      level = 0.80,
      linewidth = 0.6,
      alpha = 0.7,
      show.legend = FALSE
    ) +
    geom_point(size = 2.8, alpha = 0.85) +
    scale_colour_manual(
      values = make_group_colours(group_column),
      labels = function(x) stringr::str_wrap(x, width = 25),
      drop = TRUE
    ) +
    scale_shape_manual(
      values = c(
        "Defined community" = 16,  # Filled circle
        "grains_15"  = 17,        # Filled triangle
        "grains_30"  = 15,        # Filled square
        "grains_60"  = 18,        # Filled diamond
        "grains_120" = 8          # Star
      ),
      guide = "none"
    ) +
    coord_equal() +
    labs(
      x = sprintf("PC1 (%.1f%%)", variance_explained[1]),
      y = sprintf("PC2 (%.1f%%)", variance_explained[2]),
      colour = legend_title,
      title = paste("Coloured by", legend_title)
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      legend.key.width = grid::unit(0.45, "cm"),
      legend.spacing.x = grid::unit(0.15, "cm"),
      legend.box.margin = margin(t = 6),
      plot.title = element_text(face = "bold"),
      plot.margin = margin(10, 18, 10, 10)
    ) +
    guides(
      colour = guide_legend(
        ncol = 2,
        byrow = TRUE,
        override.aes = list(size = 3, shape = 16)
      )
    )
}

# Create five community panels ----
p_pca_yeast <- plot_pca("Yeast", "Yeast")
p_pca_bif   <- plot_pca("Bif", "Bifidobacterium")
p_pca_lab   <- plot_pca("LAB", "LAB")
p_pca_aab   <- plot_pca("AAB", "AAB")
p_pca_other <- plot_pca("other", "Z. mobilis")

# Readable variable labels ----
variable_labels <- c(
  Ethanol = "Ethanol",
  Fructose = "Fructose",
  Glucose = "Glucose",
  Glycerol = "Glycerol",
  Mannitol = "Mannitol",
  Sucrose = "Sucrose",
  Lactic_Acid = "Lactic acid",
  Acetic_Acid = "Acetic acid",
  pH_day3 = "pH after 3 days",
  Aroma_score = "Aroma score"
)

# Extract PCA loadings ----
df_pca_loadings <- as_tibble(
  pca_result$rotation,
  rownames = "Variable"
) %>%
  mutate(
    Label = unname(variable_labels[Variable])
  )

# Symmetric plotting range with space for labels ----
loading_limit <- 1.5 * max(
  abs(c(df_pca_loadings$PC1, df_pca_loadings$PC2))
)

# Loadings plot with transparent, repelled labels ----
p_pca_loadings <- ggplot(df_pca_loadings) +
  geom_hline(yintercept = 0, colour = "grey85") +
  geom_vline(xintercept = 0, colour = "grey85") +
  geom_segment(
    aes(
      x = 0,
      y = 0,
      xend = PC1,
      yend = PC2
    ),
    arrow = grid::arrow(
      length = grid::unit(0.18, "cm")
    ),
    colour = "steelblue4",
    linewidth = 0.6
  ) +
  geom_text_repel(
    aes(
      x = PC1,
      y = PC2,
      label = Label
    ),
    seed = 123,
    size = 3.8,
    box.padding = 0.6,
    point.padding = 0.3,
    force = 2,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.colour = "grey55",
    segment.size = 0.3
  ) +
  coord_equal(
    xlim = c(-loading_limit, loading_limit),
    ylim = c(-loading_limit, loading_limit),
    clip = "off"
  ) +
  labs(
    title = "Variable loadings",
    x = "PC1 loading",
    y = "PC2 loading"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.margin = margin(15, 25, 15, 25)
  )

# Combine score plots and loadings ----

p_pca_overview <- wrap_plots(
  p_pca_yeast,
  p_pca_bif,
  p_pca_lab,
  p_pca_aab,
  p_pca_other,
  p_pca_loadings,
  ncol = 3
) +
  plot_annotation(tag_levels = "a")


# Scree plot ----
p_pca_scree <- ggplot(
  pca_variance,
  aes(
    x = factor(PC, levels = PC),
    y = Variance_percent
  )
) +
  geom_col(fill = "steelblue4", width = 0.7) +
  labs(
    title = "Variance explained by each principal component",
    x = NULL,
    y = "Variance explained (%)"
  ) +
  theme_classic(base_size = 12)

print(p_pca_scree)

# Save plots and results ----
dir.create(here("results"), showWarnings = FALSE)

ggsave(
  here("output", "Figure4.pdf"),
  p_pca_overview,
  width = 12,
  height = 11,
  units = "in",
  limitsize = FALSE
)

ggsave(
  here("results", "PCA_fermentation_loadings.pdf"),
  p_pca_loadings,
  width = 8,
  height = 8
)

ggsave(
  here("results", "PCA_fermentation_scree.pdf"),
  p_pca_scree,
  width = 7,
  height = 5
)

write_csv(
  df_pca_scores,
  here("results", "PCA_fermentation_scores.csv")
)

write_csv(
  df_pca_loadings,
  here("results", "PCA_fermentation_loadings.csv")
)

write_csv(
  pca_variance,
  here("results", "PCA_fermentation_variance.csv")
)

write_csv(
  df_pca_excluded,
  here("results", "PCA_fermentation_excluded_samples.csv")
)

#notes: 
# 15 g/L triangle, 30 g/L square, 60 g/L diamond, 120 g/L star.


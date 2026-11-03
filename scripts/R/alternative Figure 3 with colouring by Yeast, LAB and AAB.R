rm(list = ls())

# Libraries ----
library(tidyverse)     # dplyr, ggplot2, readr, etc.
library(rstatix)
library(multcompView)  # multcompLetters()
library(patchwork)     # wrap_plots(), plot_annotation()
library(here)

# Data import ----
df_metabolites <- read_csv(here("data", "Data_metabolites.csv"))

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

# Assign factors ----
df_metabolites <- df_metabolites %>%
  mutate(
    Yeast = factor(Yeast, levels = yeast_levels),
    Bif   = factor(Bif,   levels = bif_levels),
    LAB   = factor(LAB,   levels = lab_levels),
    AAB   = factor(AAB,   levels = aab_levels),
    other = factor(other, levels = other_levels)
  )

# Yeast colours matching the PCA
yeast_colours <- c(
  "media"            = "grey60",
  "Br. bruxellensis" = "#F8766D",  # Salmon
  "H. valbyensis"    = "#7CAE00",  # Green
  "P. fermentans"    = "#00BFC4",  # Blue/turquoise
  "S. cerevisiae"    = "#C77CFF",  # Purple
  "WK grains"       = "black"
)

aab_colours <- c(
  "media"          = "grey60",
  "A. lovaniensis" = "#F8766D",
  "A. orientalis"  = "#7CAE00",
  "G. japonicus"   = "#00BFC4",
  "G. oxydans"     = "#C77CFF",
  "WK grains"     = "black"
)

lab_colours <- c(
  "media"                    = "grey60",
  "Lac. paracasei"           = "#F8766D",
  "Len. hilgardii"           = "#7CAE00",
  "Leu. pseudomesenteroides" = "#00BFC4",
  "Liq. satsumensis"         = "#C77CFF",
  "WK grains"               = "black"
)

# Grain concentration symbols matching the PCA
df_metabolites <- df_metabolites %>%
  mutate(
    Grain_symbol = if_else(
      Yeast == "WK grains",
      paste0("grains_", readr::parse_number(as.character(grains))),
      "Other"
    )
  )

grain_shapes <- c(
  "Other"      = 16,  # Circle
  "grains_15"  = 17,  # Triangle
  "grains_30"  = 15,  # Square
  "grains_60"  = 18,  # Diamond
  "grains_120" = 8    # Star
)

# Replace ND with 0 for key metabolites ----
df_metabolites$Mannitol[is.na(df_metabolites$Mannitol)]       <- 0
df_metabolites$Lactic_Acid[is.na(df_metabolites$Lactic_Acid)] <- 0
df_metabolites$Acetic_Acid[is.na(df_metabolites$Acetic_Acid)] <- 0

# Notes:
# Normality of residuals:
# shapiro.test(residuals(anova_res))
# Homogeneity of variances:
# car::leveneTest(formula, data = df_metabolites)


#### Combined plots: ethanol and organic acids ----

chemical_columns   <- c("Ethanol", "Acetic_Acid", "Lactic_Acid")
grouping_variables <- c("Yeast", "Bif", "LAB", "AAB", "other")

plot_list  <- list()
plot_index <- 1

for (chem in chemical_columns) {
  colour_factor <- switch(
    chem,
    "Ethanol"     = "Yeast",
    "Acetic_Acid" = "AAB",
    "Lactic_Acid" = "LAB"
  )
  
  point_colours <- switch(
    colour_factor,
    "Yeast" = yeast_colours,
    "AAB"   = aab_colours,
    "LAB"   = lab_colours
  )
  for (group in grouping_variables) {
    
    if (length(unique(df_metabolites[[group]])) <= 1) next
    
    # Kruskal-Wallis
    res.kruskal <- df_metabolites %>%
      kruskal_test(reformulate(group, response = chem))
    
    # Pairwise Wilcoxon
    pwc <- df_metabolites %>%
      pairwise_wilcox_test(
        formula = reformulate(group, response = chem),
        p.adjust.method = "BH"
      ) %>%
      add_xy_position(x = group)
    
    # Prepare p-value matrix for CLD
    pwc_clean <- pwc %>%
      select(group1, group2, p.adj)
    
    sorted_for_letters <- df_metabolites %>%
      group_by(.data[[group]]) %>%
      summarise(med = median(.data[[chem]], na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(med)) %>%
      pull(!!sym(group)) %>%
      as.character()
    
    pval_matrix <- matrix(
      1,
      nrow = length(sorted_for_letters),
      ncol = length(sorted_for_letters),
      dimnames = list(sorted_for_letters, sorted_for_letters)
    )
    
    for (i in seq_len(nrow(pwc_clean))) {
      g1   <- as.character(pwc_clean$group1[i])
      g2   <- as.character(pwc_clean$group2[i])
      pval <- pwc_clean$p.adj[i]
      pval_matrix[g1, g2] <- pval
      pval_matrix[g2, g1] <- pval
    }
    
    # CLD letters
    cld <- multcompLetters(pval_matrix)
    
    cld_df <- data.frame(
      group = names(cld$Letters),
      label = cld$Letters,
      y     = max(df_metabolites[[chem]], na.rm = TRUE) + 1
    )
    colnames(cld_df)[1] <- group
    
    subtitle_text <- get_test_label(res.kruskal, detailed = TRUE)
    
    p <- ggplot(df_metabolites, aes_string(x = group, y = chem)) +
      geom_boxplot(outlier.shape = NA, fill = "lightblue") +
      geom_jitter(
        aes(colour = !!rlang::sym(colour_factor), shape = Grain_symbol),
        width = 0.2,
        height = 0,
        size = 2,
        alpha = 0.30
      ) +
      scale_colour_manual(values = point_colours) +
      scale_shape_manual(values = grain_shapes) +
      guides(colour = "none", shape = "none") +
      geom_text(
        data = cld_df,
        aes_string(x = group, y = "y", label = "label"),
        vjust = 0,
        fontface = "bold"
      ) +
      theme_bw() +
      xlab("") +
      ylab(paste(gsub("_", " ", chem), "(µg/mL)")) +
      labs(subtitle = subtitle_text) +
      theme(
        axis.text.x    = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.subtitle  = element_text(size = 9),
        plot.caption   = element_text(size = 8, face = "italic", hjust = 0)
      )
    
    plot_list[[plot_index]] <- p
    plot_index <- plot_index + 1
  }
}

combined_plot <- wrap_plots(plot_list, ncol = 5) +
  plot_annotation(tag_levels = "a")

ggsave(
  filename = here("output", "Figure3_multi_col.pdf"),
  combined_plot,
  width  = 18,
  height = 18
)

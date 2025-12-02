rm(list = ls())

# Libraries ----
library(tidyverse)     # dplyr, ggplot2, readr, etc.
library(rstatix)
library(multcompView)  # multcompLetters()
library(patchwork)     # wrap_plots(), plot_annotation()
library(here)

# Data import ----
df_metabolites <- read_csv(here("data", "Data_Metabolites.csv"))

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
      geom_jitter(aes(color = grains), width = 0.2, alpha = 0.4) +
      guides(color = "none") +
      geom_text(
        data = cld_df,
        aes_string(x = group, y = "y", label = "label"),
        vjust = 0,
        fontface = "bold"
      ) +
      theme_bw() +
      xlab("") +
      ylab(paste(chem, "in µg/ml")) +
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
  filename = here("output", "Figure3.pdf"),
  combined_plot,
  width  = 18,
  height = 18
)


#### Combined plots: sugars and sugar alcohols ----

chemical_columns   <- c("Sucrose", "Fructose", "Glucose", "Glycerol", "Mannitol")
grouping_variables <- c("Yeast", "Bif", "LAB", "AAB", "other")

plot_list  <- list()
plot_index <- 1

for (chem in chemical_columns) {
  for (group in grouping_variables) {
    
    if (length(unique(df_metabolites[[group]])) <= 1) next
    
    res.kruskal <- df_metabolites %>%
      kruskal_test(reformulate(group, response = chem))
    
    pwc <- df_metabolites %>%
      pairwise_wilcox_test(
        formula = reformulate(group, response = chem),
        p.adjust.method = "BH"
      ) %>%
      add_xy_position(x = group)
    
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
      geom_jitter(aes(color = grains), width = 0.2, alpha = 0.4) +
      guides(color = "none") +
      geom_text(
        data = cld_df,
        aes_string(x = group, y = "y", label = "label"),
        vjust = 0,
        fontface = "bold"
      ) +
      theme_bw() +
      xlab("") +
      ylab(paste(chem, "in µg/ml")) +
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
  filename = here("output", "Figure2.pdf"),
  combined_plot,
  width  = 18,
  height = 25
)

rm(list = ls())

# Libraries ----
library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(rstatix)
library(patchwork)
library(multcompView)
library(here)

# Import data ----
df_01 <- read_excel(here("data", "Fermentation_data.xlsx"))

# Harmonise names ----
df_01$Bif <- gsub("Bifidobacterium fermentum \\(B2\\)", "B. fermentum", df_01$Bif)
df_01$Bif <- gsub("Bifidobacterium aquikefiricola \\(B1\\)", "B. aquikefiricola", df_01$Bif)

df_01$Yeast <- gsub("Brettanomyces bruxellensis \\(Y1\\)", "Br. bruxellensis", df_01$Yeast)
df_01$Yeast <- gsub("Hanseniaspora valbyensis \\(Y2\\)", "H. valbyensis", df_01$Yeast)
df_01$Yeast <- gsub("Pichia fermentans \\(Y3\\)", "P. fermentans", df_01$Yeast)
df_01$Yeast <- gsub("Saccharomyces cerevisiae \\(Y4\\)", "S. cerevisiae", df_01$Yeast)

df_01$LAB <- gsub("Lacticaseibacillus paracasei \\(L1\\)", "Lac. paracasei", df_01$LAB)
df_01$LAB <- gsub("Lentilactobacillus hilgardii \\(L2\\)", "Len. hilgardii", df_01$LAB)
df_01$LAB <- gsub("Leuconostoc pseudomesenteroides \\(L3\\)", "Leu. pseudomesenteroides", df_01$LAB)
df_01$LAB <- gsub("Liquorilactobacillus satsumensis \\(L4\\)", "Liq. satsumensis", df_01$LAB)

df_01$AAB <- gsub("Acetobacter lovaniensis \\(A1\\)", "A. lovaniensis", df_01$AAB)
df_01$AAB <- gsub("Acetobacter orientalis \\(A2\\)", "A. orientalis", df_01$AAB)
df_01$AAB <- gsub("Gluconobacter japonicus \\(A3\\)", "G. japonicus", df_01$AAB)
df_01$AAB <- gsub("Gluconobacter oxydans \\(A4\\)", "G. oxydans", df_01$AAB)

df_01$other <- gsub("Zymomonas mobilis \\(Z1\\)", "Z. mobilis", df_01$other)

# Experimental set selection ----
df_selected_set <- df_01 %>%
  filter(`experimental set` %in% c("set1", "set3", "set4", "set5"))

df_selected_set$Yeast <- factor(
  df_selected_set$Yeast,
  levels = c("Br. bruxellensis",
             "H. valbyensis",
             "P. fermentans",
             "S. cerevisiae",
             "WK grains")
)

df_selected_set$Bif <- factor(
  df_selected_set$Bif,
  levels = c("B. aquikefiricola",
             "B. fermentum",
             "WK grains")
)

df_selected_set$LAB <- factor(
  df_selected_set$LAB,
  levels = c("Lac. paracasei",
             "Len. hilgardii",
             "Leu. pseudomesenteroides",
             "Liq. satsumensis",
             "WK grains")
)

df_selected_set$AAB <- factor(
  df_selected_set$AAB,
  levels = c("A. lovaniensis",
             "A. orientalis",
             "G. japonicus",
             "G. oxydans",
             "WK grains")
)

df_selected_set$other <- factor(
  df_selected_set$other,
  levels = c("Z. mobilis",
             "none",
             "WK grains")
)

# Remove rows with "media" in combined name short ----
df_selected_set_wo_media <- df_selected_set %>%
  filter(!str_detect(`combined name short`, "media"))

# pH summaries (0h and 72h) ----
df_summary_pH_0h <- df_selected_set %>%
  group_by(`combined name short`) %>%
  summarise(mean_pH_0h = round(mean(`pH 0h`, na.rm = TRUE), 1))

df_summary_pH_72h <- df_selected_set %>%
  group_by(`combined name short`) %>%
  summarise(mean_pH_72h = round(mean(`pH 72h`, na.rm = TRUE), 1))

df_summary_pH_72h_AAB <- df_selected_set %>%
  group_by(AAB) %>%
  summarise(mean_pH_72h = round(mean(`pH 72h`, na.rm = TRUE), 1))

df_summary_pH_72h_yeast <- df_selected_set %>%
  group_by(Yeast) %>%
  summarise(mean_pH_72h = round(mean(`pH 72h`, na.rm = TRUE), 1))


#### Figure 1: average panel scores ----

df_02 <- df_selected_set %>%
  filter(`combined name short` != "media")

df_04 <- df_02 %>%
  rename(smell_z_score = `smell z-score`)

df_04$grains <- ifelse(
  grepl("^WK grains", df_04$`combined name short`),
  df_04$`combined name short`,
  NA
)

measurement_column <- c("average_pannel")
grouping_variables <- c("Yeast", "Bif", "LAB", "AAB", "other")

plot_list <- list()
plot_index <- 1

for (chem in measurement_column) {
  for (group in grouping_variables) {
    
    if (length(unique(df_04[[group]])) <= 1) next
    
    # Kruskal-Wallis
    res.kruskal <- df_04 %>%
      kruskal_test(reformulate(group, response = chem))
    
    # Pairwise Wilcoxon
    pwc <- df_04 %>%
      pairwise_wilcox_test(
        formula = reformulate(group, response = chem),
        p.adjust.method = "BH"
      ) %>%
      add_xy_position(x = group)
    
    # Prepare p-value matrix for CLD
    pwc_clean <- pwc %>%
      dplyr::select(group1, group2, p.adj)
    
    sorted_for_letters <- df_04 %>%
      group_by(.data[[group]]) %>%
      summarise(med = median(.data[[chem]], na.rm = TRUE)) %>%
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
      g1 <- as.character(pwc_clean$group1[i])
      g2 <- as.character(pwc_clean$group2[i])
      pval <- pwc_clean$p.adj[i]
      pval_matrix[g1, g2] <- pval
      pval_matrix[g2, g1] <- pval
    }
    
    # CLD letters
    cld <- multcompLetters(pval_matrix)
    
    cld_df <- data.frame(
      group = names(cld$Letters),
      label = cld$Letters,
      y = max(df_04[[chem]], na.rm = TRUE) + 1
    )
    colnames(cld_df)[1] <- group
    
    subtitle_text <- get_test_label(res.kruskal, detailed = TRUE)
    
    # Plot
    p <- ggplot(df_04, aes_string(x = group, y = chem)) +
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
      ylab(chem) +
      labs(subtitle = subtitle_text) +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.subtitle = element_text(size = 9),
        plot.caption = element_text(size = 8, face = "italic", hjust = 0)
      )
    
    plot_list[[plot_index]] <- p
    plot_index <- plot_index + 1
  }
}

combined_plot <- wrap_plots(plot_list, ncol = 5) +
  plot_annotation(tag_levels = "a")

ggsave(
  filename = here("output", "Figure1.pdf"),
  combined_plot,
  width  = 45,
  height = 12,
  units  = "cm"
)


#### Figure S12: cells/ml ----

df_04 <- df_selected_set %>%
  filter(`combined name short` != "media") %>%
  rename(
    cells_per_ml_day3 = `events/ml (day3)`,
    cells_per_ml_day0 = `events/ml (day0)`
  )

df_04$grains <- ifelse(
  grepl("^WK grains", df_04$`combined name short`),
  df_04$`combined name short`,
  NA
)

df_summary_cells_0h <- df_04 %>%
  group_by(`combined name short`) %>%
  summarise(mean_cells_per_ml_day0 = round(mean(cells_per_ml_day0, na.rm = TRUE) / 1e6, 1))

df_summary_cells_per_ml_day3 <- df_04 %>%
  group_by(`combined name short`) %>%
  summarise(mean_cells_per_ml_day3 = round(mean(cells_per_ml_day3, na.rm = TRUE) / 1e6, 1))

measurement_column <- c("cells_per_ml_day0", "cells_per_ml_day3")
grouping_variables <- c("Yeast", "Bif", "LAB", "AAB", "other")

plot_list <- list()
plot_index <- 1

for (chem in measurement_column) {
  for (group in grouping_variables) {
    
    if (length(unique(df_04[[group]])) <= 1) next
    
    res.kruskal <- df_04 %>%
      kruskal_test(reformulate(group, response = chem))
    
    pwc <- df_04 %>%
      pairwise_wilcox_test(
        formula = reformulate(group, response = chem),
        p.adjust.method = "BH"
      ) %>%
      add_xy_position(x = group)
    
    pwc_clean <- pwc %>%
      dplyr::select(group1, group2, p.adj)
    
    sorted_for_letters <- df_04 %>%
      group_by(.data[[group]]) %>%
      summarise(med = median(.data[[chem]], na.rm = TRUE)) %>%
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
      g1 <- as.character(pwc_clean$group1[i])
      g2 <- as.character(pwc_clean$group2[i])
      pval <- pwc_clean$p.adj[i]
      pval_matrix[g1, g2] <- pval
      pval_matrix[g2, g1] <- pval
    }
    
    cld <- multcompLetters(pval_matrix)
    
    cld_df <- data.frame(
      group = names(cld$Letters),
      label = cld$Letters,
      y = 9e8
    )
    colnames(cld_df)[1] <- group
    
    subtitle_text <- get_test_label(res.kruskal, detailed = TRUE)
    
    p <- ggplot(df_04, aes_string(x = group, y = chem)) +
      geom_boxplot(outlier.shape = NA, fill = "lightblue") +
      geom_jitter(aes(color = grains), width = 0.2, alpha = 0.4) +
      guides(color = "none") +
      geom_hline(yintercept = 2e6, linetype = "dashed", color = "red") +
      geom_text(
        data = cld_df,
        aes_string(x = group, y = "y", label = "label"),
        vjust = 0,
        fontface = "bold"
      ) +
      theme_bw() +
      xlab("") +
      ylab(chem) +
      labs(subtitle = subtitle_text) +
      scale_y_log10(limits = c(1e6, 1e9)) +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.subtitle = element_text(size = 9),
        plot.caption = element_text(size = 8, face = "italic", hjust = 0)
      )
    
    plot_list[[plot_index]] <- p
    plot_index <- plot_index + 1
  }
}

combined_plot <- wrap_plots(plot_list, ncol = 5) +
  plot_annotation(tag_levels = "a")

ggsave(
  filename = here("output", "FigureS12 cells per ml.pdf"),
  combined_plot,
  width  = 18,
  height = 10
)


#### Figure S11: pH ----

df_selected_set <- df_01 %>%
  filter(`experimental set` %in% c("set1", "set3", "set4", "set5"))

df_selected_set$Yeast <- factor(
  df_selected_set$Yeast,
  levels = c("media",
             "Br. bruxellensis",
             "H. valbyensis",
             "P. fermentans",
             "S. cerevisiae",
             "WK grains")
)

df_selected_set$Bif <- factor(
  df_selected_set$Bif,
  levels = c("media",
             "B. aquikefiricola",
             "B. fermentum",
             "WK grains")
)

df_selected_set$LAB <- factor(
  df_selected_set$LAB,
  levels = c("media",
             "Lac. paracasei",
             "Len. hilgardii",
             "Leu. pseudomesenteroides",
             "Liq. satsumensis",
             "WK grains")
)

df_selected_set$AAB <- factor(
  df_selected_set$AAB,
  levels = c("media",
             "A. lovaniensis",
             "A. orientalis",
             "G. japonicus",
             "G. oxydans",
             "WK grains")
)

df_selected_set$other <- factor(
  df_selected_set$other,
  levels = c("media",
             "Z. mobilis",
             "none",
             "WK grains")
)

df_04 <- df_selected_set

df_04$grains <- ifelse(
  grepl("^WK grains", df_04$`combined name short`),
  df_04$`combined name short`,
  NA
)

df_04 <- df_04 %>%
  rename(
    pH_day0 = `pH 0h`,
    pH_day3 = `pH 72h`
  )

measurement_column <- c("pH_day0", "pH_day3")
grouping_variables <- c("Yeast", "Bif", "LAB", "AAB", "other")

plot_list <- list()
plot_index <- 1

for (chem in measurement_column) {
  for (group in grouping_variables) {
    
    if (length(unique(df_04[[group]])) <= 1) next
    
    res.kruskal <- df_04 %>%
      kruskal_test(reformulate(group, response = chem))
    
    pwc <- df_04 %>%
      pairwise_wilcox_test(
        formula = reformulate(group, response = chem),
        p.adjust.method = "BH"
      ) %>%
      add_xy_position(x = group)
    
    pwc_clean <- pwc %>%
      dplyr::select(group1, group2, p.adj)
    
    sorted_for_letters <- df_04 %>%
      group_by(.data[[group]]) %>%
      summarise(med = median(.data[[chem]], na.rm = TRUE)) %>%
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
      g1 <- as.character(pwc_clean$group1[i])
      g2 <- as.character(pwc_clean$group2[i])
      pval <- pwc_clean$p.adj[i]
      pval_matrix[g1, g2] <- pval
      pval_matrix[g2, g1] <- pval
    }
    
    cld <- multcompLetters(pval_matrix)
    
    cld_df <- data.frame(
      group = names(cld$Letters),
      label = cld$Letters,
      y = max(df_04[[chem]], na.rm = TRUE) + 1
    )
    colnames(cld_df)[1] <- group
    
    subtitle_text <- get_test_label(res.kruskal, detailed = TRUE)
    
    p <- ggplot(df_04, aes_string(x = group, y = chem)) +
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
      ylab(chem) +
      labs(subtitle = subtitle_text) +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.subtitle = element_text(size = 9),
        plot.caption = element_text(size = 8, face = "italic", hjust = 0)
      )
    
    plot_list[[plot_index]] <- p
    plot_index <- plot_index + 1
  }
}

combined_plot <- wrap_plots(plot_list, ncol = 5) +
  plot_annotation(tag_levels = "a")

ggsave(
  filename = here("output", "FigureS11 pH.pdf"),
  combined_plot,
  width  = 18,
  height = 10
)

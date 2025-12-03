#### analysis of taste panel results #### 
rm(list = ls())

#### loading pkgs ####
library(readxl)
library(ggplot2)
library(dplyr)
library(rstatix)
library(patchwork)
library(multcompView)
library(here)

df_01 = read.csv(here("data", "Data_tastepanel.csv"))

#assign levels
df_01$sample <- factor(df_01$sample, levels = c(
  "B.fer S.cer L.hil G.oxy",
  "B.fer H.val L.par G.jap Z.mob",
  "B.fer P.fer L.hil A.lov Z.mob",
  "B.acola S.cer L.sat G.jap Z.mob",
  "B.acola S.cer L.sat G.jap",
  "WK grains (30g/L)"
))


#### Sensory Acceptance Test plots ####
#### bar plot  #### 
#### selecting the set #### 
descriptor_column <- c("appearance", "liking_of_aroma", "liking_of_flavour", "liking_of_viscosity", "overall_acceptability")
grouping_variables <- c("sample")

plot_list <- list()
plot_index <- 1

for (descriptor in descriptor_column) {
  for (group in grouping_variables) {
    
    if (length(unique(df_01[[group]])) <= 1) next
    
    # Perform Kruskal-Wallis and pairwise Wilcoxon test
    res.kruskal <- df_01 %>%
      kruskal_test(reformulate(group, response = descriptor))
    
    pwc <- df_01 %>%
      pairwise_wilcox_test(
        formula = reformulate(group, response = descriptor),
        p.adjust.method = "BH"
      ) %>%
      add_xy_position(x = group)
    
    # Create CLD using multcompLetters
    pwc_clean <- pwc %>%
      dplyr::select(group1, group2, p.adj)

    # Use original plotting order
    plotting_levels <- levels(df_01[[group]])
    
    # Use descending median order ONLY for assigning letters
    sorted_for_letters <- df_01 %>%
      group_by(.data[[group]]) %>%
      summarise(med = median(.data[[descriptor]], na.rm = TRUE)) %>%
      arrange(desc(med)) %>%
      pull(!!sym(group)) %>%
      as.character()
    
    # Build p-value matrix using sorted order for letters
    pval_matrix <- matrix(1, nrow = length(sorted_for_letters), ncol = length(sorted_for_letters),
                          dimnames = list(sorted_for_letters, sorted_for_letters))
    
    for (i in 1:nrow(pwc_clean)) {
      g1 <- as.character(pwc_clean$group1[i])
      g2 <- as.character(pwc_clean$group2[i])
      pval <- pwc_clean$p.adj[i]
      pval_matrix[g1, g2] <- pval
      pval_matrix[g2, g1] <- pval
    }
    
    # Assign letters
    cld <- multcompLetters(pval_matrix)
    
    # Create CLD dataframe using original group levels (for plot)
    cld_df <- data.frame(
      group = names(cld$Letters),
      label = cld$Letters,
      y = max(df_01[[descriptor]], na.rm = TRUE) + 1
    )
    colnames(cld_df)[1] <- group
    
    
    cld_df <- data.frame(
      group = names(cld$Letters),
      label = cld$Letters
    )
    colnames(cld_df)[1] <- group
    
    # Create summary data (mean + sem)
    summary_df <- df_01 %>%
      group_by(across(all_of(group))) %>%
      summarise(
        mean = mean(.data[[descriptor]], na.rm = TRUE),
        sem = sd(.data[[descriptor]], na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
      )
    
    # Merge for label positioning
    cld_df <- left_join(cld_df, summary_df, by = group)
    cld_df$y <- cld_df$mean + cld_df$sem + 0.5  # Offset for label position
    
    # Generate subtitle and caption
    subtitle_text <- get_test_label(res.kruskal, detailed = TRUE)
    caption_text <- get_pwc_label(pwc)
    
    # Build barplot with SEM
    p <- ggplot(summary_df, aes_string(x = group, y = "mean")) +
      geom_bar(stat = "identity", fill = "skyblue") +
      geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.2) +
      geom_text(data = cld_df, aes_string(x = group, y = "y", label = "label"),
                vjust = 0, fontface = "bold") +
      theme_bw() +
      xlab("") +
      ylab(descriptor) +
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

# Combine all plots
combined_barplot_Sensory_Acceptance_Test <- wrap_plots(plot_list, ncol = 5) +
  plot_annotation(tag_levels = 'a') 

# Save
ggsave(
  filename = here("output", "Figure4.pdf")  , 
  combined_barplot_Sensory_Acceptance_Test, 
  width = 20, height = 8)

#### Optimized Descriptive Analysis plots ####
#### bar plot  #### 
#### selecting the set #### 
#removing spaces from column names
descriptor_column <- c("off_flavour", 
                       "flavour_intensity", 
                       "kefir_like_flavour", 
                       "medicinal_flavour", 
                       "yeasty_flavour", 
                       "fruity_estery_flavour", 
                       "berry_flavour", 
                       "apple_juice_flavour", 
                       "apple_cider_flavour", 
                       "pear_flavour", 
                       "vinegar_flavour", 
                       "acetone_flavour", 
                       "sweet_taste", 
                       "sour_taste", 
                       "effervescence", 
                       "alcohol_after_taste")

grouping_variables <- c("sample")

plot_list <- list()
plot_index <- 1

for (descriptor in descriptor_column) {
  for (group in grouping_variables) {
    
    if (length(unique(df_01[[group]])) <= 1) next
    
    # Perform Kruskal-Wallis and pairwise Wilcoxon test
    res.kruskal <- df_01 %>%
      kruskal_test(reformulate(group, response = descriptor))
    
    pwc <- df_01 %>%
      pairwise_wilcox_test(
        formula = reformulate(group, response = descriptor),
        p.adjust.method = "BH"
      ) %>%
      add_xy_position(x = group)
    
    # Create CLD using multcompLetters
    pwc_clean <- pwc %>%
      dplyr::select(group1, group2, p.adj)
    
    # Use original plotting order
    plotting_levels <- levels(df_01[[group]])
    
    # Use descending median order ONLY for assigning letters
    sorted_for_letters <- df_01 %>%
      group_by(.data[[group]]) %>%
      summarise(med = median(.data[[descriptor]], na.rm = TRUE)) %>%
      arrange(desc(med)) %>%
      pull(!!sym(group)) %>%
      as.character()
    
    # Build p-value matrix using sorted order for letters
    pval_matrix <- matrix(1, nrow = length(sorted_for_letters), ncol = length(sorted_for_letters),
                          dimnames = list(sorted_for_letters, sorted_for_letters))
    
    for (i in 1:nrow(pwc_clean)) {
      g1 <- as.character(pwc_clean$group1[i])
      g2 <- as.character(pwc_clean$group2[i])
      pval <- pwc_clean$p.adj[i]
      pval_matrix[g1, g2] <- pval
      pval_matrix[g2, g1] <- pval
    }
    
    # Assign letters
    cld <- multcompLetters(pval_matrix)
    
    # Create CLD dataframe using original group levels (for plot)
    cld_df <- data.frame(
      group = names(cld$Letters),
      label = cld$Letters,
      y = max(df_01[[descriptor]], na.rm = TRUE) + 1
    )
    colnames(cld_df)[1] <- group
    
    cld_df <- data.frame(
      group = names(cld$Letters),
      label = cld$Letters
    )
    colnames(cld_df)[1] <- group
    
    # Create summary data (mean + sem)
    summary_df <- df_01 %>%
      group_by(across(all_of(group))) %>%
      summarise(
        mean = mean(.data[[descriptor]], na.rm = TRUE),
        sem = sd(.data[[descriptor]], na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
      )
    
    # Merge for label positioning
    cld_df <- left_join(cld_df, summary_df, by = group)
    cld_df$y <- cld_df$mean + cld_df$sem + 0.5  # Offset for label position
    
    # Generate subtitle and caption
    subtitle_text <- get_test_label(res.kruskal, detailed = TRUE)
    caption_text <- get_pwc_label(pwc)
    
    # Build barplot with SEM
    p <- ggplot(summary_df, aes_string(x = group, y = "mean")) +
      geom_bar(stat = "identity", fill = "skyblue") +
      geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.2) +
      geom_text(data = cld_df, aes_string(x = group, y = "y", label = "label"),
                vjust = 0, fontface = "bold") +
      theme_bw() +
      xlab("") +
      ylab(descriptor) +
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

# Combine all plots
combined_barplot_Optimized_Descriptive_Test <- wrap_plots(plot_list, ncol = 4) +
  plot_annotation(tag_levels = 'a') 

# Save
ggsave(
  filename = here("output", "Figure5.pdf"), 
       combined_barplot_Optimized_Descriptive_Test, 
       width = 16, height = 20)


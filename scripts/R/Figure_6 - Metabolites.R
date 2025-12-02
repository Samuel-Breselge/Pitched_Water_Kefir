#### taste panel plots #### 
rm(list=ls())

library(dplyr)
library(ggplot2)
library(rstatix)
library(patchwork)
library(multcompView)
library(here)

df_01 = read.csv(here("data", "Data_metabolites_taste_panel.csv"))

df_01 <- df_01 %>%
  mutate(across(c(Sucrose, Glucose, Fructose, Ethanol, Glycerol, Lactic_Acid, Acetic_Acid), as.numeric))

#assign levels
df_01$short_name <- factor(df_01$short_name, levels = c(
  "media",
  "B.fer S.cer L.hil G.oxy",
  "B.fer P.fer L.hil A.lov Z.mob",
  "B.fer H.val L.par G.jap Z.mob",
  "B.acola S.cer L.sat G.jap Z.mob",
  "B.acola S.cer L.sat G.jap", 
  "WK grains (30g/L)"
))


#### 
chemical_columns <- c("Sucrose", "Glucose", "Fructose", "Ethanol", "Glycerol", "Lactic_Acid", "Acetic_Acid")
grouping_variables <- c("short_name")

plot_list <- list()
plot_index <- 1

for (chem in chemical_columns) {
  for (group in grouping_variables) {
    
    if (length(unique(df_01[[group]])) <= 1) next
    
    # Perform Kruskal-Wallis and pairwise Wilcoxon test
    res.kruskal <- df_01 %>%
      kruskal_test(reformulate(group, response = chem))
    
    pwc <- df_01 %>%
      pairwise_wilcox_test(
        formula = reformulate(group, response = chem),
        p.adjust.method = "BH"
      ) %>%
      add_xy_position(x = group)
    
    # Create CLD using multcompLetters (from pairwise p.adj)
    pwc_clean <- pwc %>%
      dplyr::select(group1, group2, p.adj)
    
    # Use original plotting order
    plotting_levels <- levels(df_01[[group]])
    
    # Use descending median order ONLY for assigning letters
    sorted_for_letters <- df_01 %>%
      group_by(.data[[group]]) %>%
      summarise(med = median(.data[[chem]], na.rm = TRUE)) %>%
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
      y = max(df_01[[chem]], na.rm = TRUE) + 1
    )
    colnames(cld_df)[1] <- group
    
    
    cld_df <- data.frame(
      group = names(cld$Letters),
      label = cld$Letters,
      y = max(df_01[[chem]], na.rm = TRUE) + 1
    )
    colnames(cld_df)[1] <- group
    
    # Generate subtitle and caption using rstatix helpers
    subtitle_text <- get_test_label(res.kruskal, detailed = TRUE)
    caption_text <- get_pwc_label(pwc)
    
    # Build the plot
    p <- ggplot(df_01, aes_string(x = group, y = chem)) +
      geom_boxplot(outlier.shape = NA,  fill = "lightblue") +
      geom_jitter(width = 0.2, alpha = 0.4) +
      geom_text(data = cld_df, aes_string(x = group, y = "y", label = "label"),
                vjust = 0, fontface = "bold") +
      theme_bw() +
      xlab("") +
      ylab(paste(chem, "in µg/ml")) +
      labs(#caption = caption_text, 
        subtitle = subtitle_text
      ) +
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
combined_plot_taste_panel <- wrap_plots(plot_list, ncol = 4)
# Stack and label
final_plot_taste_panel <- (combined_plot_taste_panel) +
  plot_annotation(tag_levels = 'a')  # This will automatically tag with a), b)


# Save
ggsave(
  filename = here("output", "Figure6.pdf"),
  final_plot_taste_panel, 
  width = 15, height = 12)


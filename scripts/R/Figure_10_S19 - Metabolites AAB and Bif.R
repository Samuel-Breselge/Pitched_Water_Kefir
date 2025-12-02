#### Bif and AAB plots #### 
rm(list = ls())

library(dplyr)
library(ggplot2)
library(rstatix)
library(multcompView)
library(patchwork)
library(here)

df_01 = read.csv(here("data", "Data_metabolites_Bif_and_AAB.csv"))

#remove AAB samples 
df_01 <- df_01[df_01$`set.name` != "influence of AAB", ]
df_01 <- df_01[df_01$other != "Zymomonas mobilis (Z1)", ]


df_01 <- df_01 %>%
  mutate(across(c(Sucrose, Glucose, Fructose, Ethanol, Glycerol, Lactic_Acid, Acetic_Acid
  ), as.numeric))

unique(df_01$short_name)

#assign levels
df_01$short_name <- factor(df_01$short_name, levels = c(
  "media",
  "B.fer S.cer L.hil G.oxy",
  "B.acola S.cer L.hil G.oxy",
  "B.aqui S.cer L.hil G.oxy",
  "S.cer L.hil G.oxy"
))


#### 
chemical_columns <- c("Sucrose", "Glucose", "Fructose", "Ethanol", "Glycerol", "Lactic_Acid", "Acetic_Acid"
)
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
      geom_boxplot(outlier.shape = NA, fill = "lightblue") +
      geom_jitter(, width = 0.2, alpha = 0.4) +
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
combined_plot_Bif1 <- wrap_plots(plot_list, ncol = 1)


#### selecting second Bif set #### 
df_01 = read.csv(here("data", "Data_metabolites_Bif_and_AAB.csv"))
#remove AAB samples 
df_01 <- df_01[df_01$`set.name` != "influence of AAB", ]
df_01 <- df_01[df_01$other != "none", ]

df_01 <- df_01 %>%
  mutate(across(c(Sucrose, Glucose, Fructose, Ethanol, Glycerol
                  , Lactic_Acid, Acetic_Acid
  ), as.numeric))

unique(df_01$short_name)

#assign levels
df_01$short_name <- factor(df_01$short_name, levels = c(
  "media",
  "B.fer S.cer L.sat G.jap Z.mob",
  "B.acola S.cer L.sat G.jap Z.mob",
  "B.aqui S.cer L.sat G.jap Z.mob",
  "S.cer L.sat G.jap Z.mob"
))

#### 
chemical_columns <- c("Sucrose", "Glucose", "Fructose", "Ethanol", "Glycerol"
                      , "Lactic_Acid", "Acetic_Acid"
)
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
      geom_boxplot(outlier.shape = NA, fill = "lightblue") +
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
combined_plot_Bif2 <- wrap_plots(plot_list, ncol = 1)


library(patchwork)

# Stack and label
final_plot_Bif <- (combined_plot_Bif1 | combined_plot_Bif2) +
  plot_annotation(tag_levels = 'a')  # This will automatically tag with a), b)

# Display
print(final_plot_Bif)

ggsave(
  filename = here("output", "FigureS19 - Bif Metabolites.pdf"),
  plot = final_plot_Bif,
  width = 8,   # adjust depending on your plot size
  height = 35,
  units = "in"
)




#### selecting first AAB set #### 
rm(list = ls())

library(dplyr)
library(ggplot2)
library(rstatix)
library(multcompView)
library(patchwork)
library(here)


df_01 = read.csv(here("data", "Data_metabolites_Bif_and_AAB.csv"))
# Define the names you want to keep
keep_names <- c("media", "B.fer S.cer L.hil", "B.fer S.cer L.hil G.oxy")

# Filter the dataframe
df_01 <- df_01 %>%
  filter(`short_name` %in% keep_names)

df_01 <- df_01 %>%
  mutate(across(c(Sucrose, Glucose, Fructose, Ethanol, Glycerol
                  , Lactic_Acid, Acetic_Acid
  ), as.numeric))


#assign levels
df_01$short_name <- factor(df_01$short_name, levels = c(
  "media",
  "B.fer S.cer L.hil G.oxy",
  "B.fer S.cer L.hil"
))


#### 
chemical_columns <- c("Sucrose", "Glucose", "Fructose", "Ethanol", "Glycerol"
                      , "Lactic_Acid", "Acetic_Acid"
)
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
      geom_boxplot(outlier.shape = NA, fill = "lightblue") +
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
combined_plot_AAB1 <- wrap_plots(plot_list, ncol = 1)


#### selecting second AAB set #### 
df_01 = read.csv(here("data", "Data_metabolites_Bif_and_AAB.csv"))

# Define the names you want to keep
keep_names <- c("media", "B.acola S.cer L.sat G.jap Z.mob", "B.acola S.cer L.sat Z.mob")

# Filter the dataframe
df_01 <- df_01 %>%
  filter(`short_name` %in% keep_names)

df_01 <- df_01 %>%
  mutate(across(c(Sucrose, Glucose, Fructose, Ethanol, Glycerol
                  , Lactic_Acid, Acetic_Acid
  ), as.numeric))

#assign levels
df_01$short_name <- factor(df_01$short_name, levels = c(
  "media",
  "B.acola S.cer L.sat G.jap Z.mob",
  "B.acola S.cer L.sat Z.mob"
))

#### 
chemical_columns <- c("Sucrose", "Glucose", "Fructose", "Ethanol", "Glycerol"
                      , "Lactic_Acid", "Acetic_Acid"
)
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
      geom_boxplot(outlier.shape = NA, fill = "lightblue") +
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
combined_plot_AAB2 <- wrap_plots(plot_list, ncol = 1)

library(patchwork)

# Stack and label
final_plot_AAB <- (combined_plot_AAB1 | combined_plot_AAB2) +
  plot_annotation(tag_levels = 'a')  # This will automatically tag with a), b)

ggsave(
  filename = here("output", "Figure10 - AAB Metabolites.pdf"),
  plot = final_plot_AAB,
  width = 8,   # adjust depending on your plot size
  height = 35,
  units = "in"
)



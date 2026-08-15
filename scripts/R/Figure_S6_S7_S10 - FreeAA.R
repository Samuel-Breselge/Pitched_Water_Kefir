#### Bif and AAB plots #### 
rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(rstatix)
library(patchwork)
library(multcompView)
library(here)

chemical_columns <- c(
  "Alanine",
  "Arginine",
  "Aspartic_Acid",
  "Cysteic_Acid",
  "Cysteine",
  "GABA",
  "Glutamic_Acid",
  "Glycine",
  "Histidine",
  "Leucine",
  "Lysine",
  "Proline",
  "Serine",
  "Threonine",
  "Valine"
)

grouping_variables <- c("short_name")

df_00 = read.csv(here("data", "Data_FreeAA.csv"), check.names = F)
  
  #### selecting first Bif set #### 
  #remove AAB samples 
  df_01 <- df_00[df_00$`set name` != "influence of AAB", ]
  df_01 <- df_01[df_01$other != "Zymomonas mobilis (Z1)", ]
  
  #assign levels
  df_01$short_name <- factor(df_01$short_name, levels = c(
    "media",
    "B.fer S.cer L.hil G.oxy",
    "B.acola S.cer L.hil G.oxy",
    "B.aqui S.cer L.hil G.oxy",
    "S.cer L.hil G.oxy"
  ))

#### bar plot Bif 1 #### 
  df_01$Norleucine = NULL
  df_01$Ammonia = NULL
  # Step 1: Summarize amino acid columns by short_name (mean or sum, depending on your need)
  df_summary <- df_01 %>%
    group_by(short_name) %>%
    summarise(across(2:21, ~ mean(.x, na.rm = TRUE)))
  
  # Step 2: Convert to long format for plotting
  df_long <- df_summary %>%
    pivot_longer(cols = -short_name, names_to = "amino_acid", values_to = "amount")
  
  # Step 3: Plot as stacked bar chart
  p1 <- ggplot(df_long, aes(x = short_name, y = amount, fill = amino_acid)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(
      #title = "Free Amino Acid Composition per Condition",
      x = "",
      y = "Free amino acid content (nmol/ml)",
      fill = "Amino Acid"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # improve readability
  
  p1

  
####
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
        ylab(paste(gsub("_", " ", chem), "(nmol/ml)")) +
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
  combined_plot_Bif1 <- wrap_plots(plot_list, ncol = 5)



#### selecting second Bif set #### 
#remove AAB samples 
df_01 <- df_00[df_00$`set name` != "influence of AAB", ]
df_01 <- df_01[df_01$other != "none", ]



#assign levels
df_01$short_name <- factor(df_01$short_name, levels = c(
  "media",
  "B.fer S.cer L.sat G.jap Z.mob",
  "B.acola S.cer L.sat G.jap Z.mob",
  "B.aqui S.cer L.sat G.jap Z.mob",
  "S.cer L.sat G.jap Z.mob"
))

#### bar plot Bif 2 #### 
df_01$Norleucine = NULL
df_01$Ammonia = NULL
# Step 1: Summarize amino acid columns by short_name (mean or sum, depending on your need)
df_summary <- df_01 %>%
  group_by(short_name) %>%
  summarise(across(2:21, ~ mean(.x, na.rm = TRUE)))  

# Step 2: Convert to long format for plotting
df_long <- df_summary %>%
  pivot_longer(cols = -short_name, names_to = "amino_acid", values_to = "amount")

# Step 3: Plot as stacked bar chart
p2 <- ggplot(df_long, aes(x = short_name, y = amount, fill = amino_acid)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    #title = "Free Amino Acid Composition per Condition",
    x = "",
    y = "Free amino acid content (nmol/ml)",
    fill = "Amino Acid"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # improve readability

p2

#### 
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
      ylab(paste(gsub("_", " ", chem), "(nmol/ml)")) +
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
combined_plot_Bif2 <- wrap_plots(plot_list, ncol = 5)

library(patchwork)

# Stack and label
final_plot_Bif <- (combined_plot_Bif1 / combined_plot_Bif2) +
  plot_annotation(tag_levels = '1')  # This will automatically tag with 1, 2, 3... 

# Display
print(final_plot_Bif)

ggsave(
  filename = here("output","FigureS10 - Bif FreeAA.pdf"),
  plot = final_plot_Bif,
  width = 18,   # adjust depending on your plot size
  height = 30,
  units = "in"
)




#### selecting first AAB set #### 
# Define the names you want to keep
keep_names <- c("media", "B.fer S.cer L.hil", "B.fer S.cer L.hil G.oxy")

# Filter the dataframe
df_01 <- df_00 %>%
  filter(`short_name` %in% keep_names)

#assign levels
df_01$short_name <- factor(df_01$short_name, levels = c(
  "media",
  "B.fer S.cer L.hil G.oxy",
  "B.fer S.cer L.hil"
))


#### bar plot AAB1 #### 
df_01$Norleucine = NULL
df_01$Ammonia = NULL
# Step 1: Summarize amino acid columns by short_name (mean or sum, depending on your need)
df_summary <- df_01 %>%
  group_by(short_name) %>%
  summarise(across(2:21, ~ mean(.x, na.rm = TRUE)))  

# Step 2: Convert to long format for plotting
df_long <- df_summary %>%
  pivot_longer(cols = -short_name, names_to = "amino_acid", values_to = "amount")

# Step 3: Plot as stacked bar chart
p3 <- ggplot(df_long, aes(x = short_name, y = amount, fill = amino_acid)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    #title = "Free Amino Acid Composition per Condition",
    x = "",
    y = "Free amino acid content (nmol/ml)",
    fill = "Amino Acid"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # improve readability

p3

#### 
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
      ylab(paste(gsub("_", " ", chem), "(nmol/ml)")) +
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
combined_plot_AAB1 <- wrap_plots(plot_list, ncol = 5)


#### selecting second AAB set #### 
# Define the names you want to keep
keep_names <- c("media", "B.acola S.cer L.sat G.jap Z.mob", "B.acola S.cer L.sat Z.mob")

# Filter the dataframe
df_01 <- df_00 %>%
  filter(`short_name` %in% keep_names)

#assign levels
df_01$short_name <- factor(df_01$short_name, levels = c(
  "media",
  "B.acola S.cer L.sat G.jap Z.mob",
  "B.acola S.cer L.sat Z.mob"
))

#### bar plot AAB2 #### 
df_01$Norleucine = NULL
df_01$Ammonia = NULL
# Step 1: Summarize amino acid columns by short_name (mean or sum, depending on your need)
df_summary <- df_01 %>%
  group_by(short_name) %>%
  summarise(across(2:21, ~ mean(.x, na.rm = TRUE)))  

# Step 2: Convert to long format for plotting
df_long <- df_summary %>%
  pivot_longer(cols = -short_name, names_to = "amino_acid", values_to = "amount")

# Step 3: Plot as stacked bar chart
p4 <- ggplot(df_long, aes(x = short_name, y = amount, fill = amino_acid)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    #title = "Free Amino Acid Composition per Condition",
    x = "",
    y = "Free amino acid cContent (nmol/ml)",
    fill = "Amino Acid"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # improve readability

p4

####
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
      ylab(paste(gsub("_", " ", chem), "(nmol/ml)")) +
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
combined_plot_AAB2 <- wrap_plots(plot_list, ncol = 5)

# Stack and label
final_plot_AAB <- (combined_plot_AAB1 / combined_plot_AAB2) +
  plot_annotation(tag_levels = '1')  

ggsave(
  filename = here("output", "FigureS7 - AAB FreeAA.pdf"),
  plot = final_plot_AAB,
  width = 18,   # adjust depending on your plot size
  height = 30,
  units = "in"
)


#### combined bar plots #### 
# Remove legends from all but p2
p1_legend <- p1  + theme(legend.position = "bottom") + labs(fill = NULL) 
p2_nolegend <- p2 + theme(legend.position = "none")
p3_nolegend <- p3 + theme(legend.position = "none")
p4_nolegend <- p4 + theme(legend.position = "none")

# Combine plots in a 2x2 grid
combined_plot <- (p3_nolegend | p4_nolegend) / (p1_legend | p2_nolegend)

# Add subplot tags
combined_plot <- combined_plot +
  plot_annotation(tag_levels = 'a')

# Print and save
print(combined_plot)
ggsave(
  filename = here("output", "FigureS6 - FreeAA.pdf"), 
  combined_plot, width = 10, height = 12, dpi = 600)


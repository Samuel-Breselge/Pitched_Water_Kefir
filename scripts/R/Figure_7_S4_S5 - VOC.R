rm(list=ls())

library(pheatmap)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
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

# Then make the matrix (excluding the combined_name column if needed)
m_VOC_01 <- as.matrix(df_VOC_01[, 1:83])

m_VOC_01 = t(m_VOC_01)

#getting the meta data
df_VOC_info = data.frame(row.names = colnames(m_VOC_01), 
                         sample_name = df_VOC_01$sample_name, 
                         tech_rep = df_VOC_01$VOC_replicate, 
                         exp_rep = df_VOC_01$Exp_replicate)

#import VOC class info
df_VOC_class_info <- read.csv(
  here("data", "VOC_classes.csv"),
  check.names = FALSE, 
  row.names = 1
)

# Calculate mean for each chemical per sample_name
df_VOC_mean <- df_VOC_01 %>%
  group_by(sample_name) %>%
  summarise(across(1:83, mean, na.rm = TRUE))

df_VOC_long <- df_VOC_mean %>%
  pivot_longer(
    cols = 2:84,  # 1st col is sample_name, 2:84 are chemicals
    names_to = "chemical",
    values_to = "mean_abundance"
  )

df_VOC_class_info <- df_VOC_class_info %>%
  rownames_to_column(var = "chemical")


df_VOC_long <- df_VOC_long %>%
  left_join(df_VOC_class_info, by = "chemical")
#####################################################################
#### bar plots ####
VOC_barplot <- ggplot(df_VOC_long, aes(x = sample_name, y = mean_abundance, fill = class)) +
  geom_bar(stat = "identity") +
  #facet_wrap(~ class, scales = "free_y") +  # free_y allows each facet to scale independently
  labs(
    x = "",
    y = "Mean abundance (RU)",
    fill = "Chemical Class",
    title = ""
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    strip.text = element_text(face = "bold")
  )

ggsave(
  here("output", "FigureS4 - VOC barplot.pdf"), 
  plot = VOC_barplot, 
  scale = 1,  width =17, height = 14, units = 'cm', dpi = 900)

###
# -------------------------------------------------------------
# Get top 7 most abundant chemicals
top7_chemicals <- df_VOC_long %>%
  group_by(chemical) %>%
  summarise(total_abundance = sum(mean_abundance, na.rm = TRUE)) %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 7) %>%
  pull(chemical)

# Filter to top 7
df_VOC_long_top7 <- df_VOC_long %>%
  filter(chemical %in% top7_chemicals)

# Make the plot
VOC_barplot_top7 <- ggplot(df_VOC_long_top7, aes(x = sample_name, y = mean_abundance, fill = chemical)) +
  geom_bar(stat = "identity") +
  labs(
    x = "",
    y = "Mean abundance (RU)",
    fill = "Chemical",
    title = ""
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(
  filename = here("output", "FigureS5 - VOC top 7.pdf"), 
  plot = VOC_barplot_top7, 
  scale = 1,  width =17, height = 14, units = 'cm', dpi = 900)

#####################################################################
#### heat maps with average data #### 
# ) Create initial annotation
df_VOC_class_info <- data.frame(
  row.names = df_VOC_class_info$chemical,
  class = df_VOC_class_info$class
)

m_VOC_mean <- df_VOC_mean %>%
  tibble::column_to_rownames(var = "sample_name") %>%
  as.matrix()


# 2) Make initial heatmap to get dendrogram
plot04 <- pheatmap(
  t(m_VOC_mean),  
  annotation_row = df_VOC_class_info,
  cluster_cols = TRUE,
  cutree_rows = 10,
  scale = "row"
)

# 3) Cut tree and get membership
row_clusters <- cutree(plot04$tree_row, k = 10)

# 4) Get order of rows (top-to-bottom)
row_order <- plot04$tree_row$order

# 5) Map clusters to tree order: assign ordered cluster IDs
cluster_df <- data.frame(
  name = rownames(t(m_VOC_mean))[row_order],
  cluster = row_clusters[row_order]
)

cluster_levels <- unique(cluster_df$cluster)
remap <- setNames(seq_along(cluster_levels), cluster_levels)
row_clusters_ordered <- remap[as.character(row_clusters)]

# 6) Add ordered cluster info with lowercase name
df_VOC_class_info$cluster <- as.factor(row_clusters_ordered)

# 7) Redraw heatmap with new annotation
plot04 <- pheatmap(
  t(m_VOC_mean),
  annotation_row = df_VOC_class_info,
  cluster_cols = TRUE,
  cutree_rows = 10,
  scale = "row"
)

# 8) Save
ggsave(
  filename = here("output", "Figure7.pdf"),
  plot = plot04,
  scale = 1,
  width = 28, height = 35, units = 'cm',
  dpi = 900
)
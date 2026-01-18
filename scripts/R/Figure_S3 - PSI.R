rm(list=ls())

df_01 = read.csv(here("data", "Data_PSI.csv"), check.names = F)

# Assign factor levels in the same order as replacements
df_01$yeast <- factor(df_01$yeast, levels = c(
  "B.bru B.acola L.sat A.lov",
  "H.val B.acola L.sat A.lov",
  "P.fer B.acola L.sat A.lov",
  "S.cer B.acola L.sat A.lov",
  "Z.mob B.acola L.sat A.lov",
  "Grains (30g/L)"
))

# Define the shape codes (0–25 are available in base R)
shape_vals <- c(
  "B.bru B.acola L.sat A.lov" = 0,  
  "H.val B.acola L.sat A.lov" = 1,  
  "P.fer B.acola L.sat A.lov" = 2, 
  "S.cer B.acola L.sat A.lov" = 3, 
  "Z.mob B.acola L.sat A.lov" = 5,   
  "Grains (30g/L)" = 6              
)

# Load ggplot2
library(ggplot2)

# Create the ggplot
ggplot(df_01, aes(x = time_in_days, y = `psi normalized at 23C`, group = yeast, shape = yeast)) +
  geom_line(aes(color = yeast), size = 1) +  # Connect points with lines, color by Sample
  geom_point(size = 2) +       # Add points with different shapes
  labs(
    x = "Time in Days",
    y = "PSI at 23°C",
    title = ""
  ) +
  scale_x_continuous(breaks = seq(0, 28, by = 2)) +
  scale_y_continuous(breaks = seq(0, 60, by = 5), limits = c(0, 60)) +
  scale_shape_manual(values = shape_vals) +  # Assign shapes
  theme_bw() +  # Use a clean theme
  theme(
    legend.title = element_blank(),  # Remove legend title
    legend.text = element_text(size = 9),
    plot.title = element_text(hjust = 0.5, size = 14)
  ) +
  facet_wrap(~ condition)  # Split the plot into separate panels by 'set'

# Save the plot as a PDF
ggsave(
  here("output", "FigureS3 - PSI.pdf"), 
  width = 12, height = 6, device = "pdf")



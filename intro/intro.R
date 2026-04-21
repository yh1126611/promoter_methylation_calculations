library(ggplot2)
library(scales)
library(data.table)
library(dplyr)

# Nature figure theme (linewidth=0.25 per Nature guidelines)
nature_theme <- theme(
  panel.background = element_blank(),
  plot.background = element_blank(),
  panel.grid = element_blank(),
  axis.line.y.left = element_line(linewidth = 0.25, color = "black"),
  axis.line.x.bottom = element_line(linewidth = 0.25, color = "black"),
  axis.ticks = element_blank(),
  axis.title = element_text(size = 6, color = "black"),
  axis.text = element_text(size = 5, color = "black"),
  legend.position = "none"
)

# Color schemes
gc_colors <- c("GC" = "#E7298A", "None" = "#BBBBBB")
mp_gradient_colors <- c(low = "#9EC9E2", high = "#0D4A70")

# Plot dimensions
axis_line_width <- 0.25
font_title_size <- 6
font_text_size <- 5
line_dotted_width <- 0.25

# Fig. 1b
ggplot(df, aes(x = genome_id, y = length, fill = type, group = type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("#BBBBBB", "transparent")) +
  labs(y = "Proportion", x = "Genome") +
  nature_theme +
  theme(axis.line = element_line(linewidth = axis_line_width))

# Fig. 1c
ggplot(df, aes(x = genome_id, y = proportion, fill = feature_color, group = feature_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_identity() +
  labs(y = "Proportion", x = "Genome") +
  nature_theme +
  theme(axis.line = element_line(linewidth = axis_line_width))

# Fig. 1d
ggplot(df, aes(x = dinucleotide, y = count, fill = "black", group = dinucleotide)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_identity() +
  labs(y = "Length (Mbp)", x = "Dinucleotide") +
  scale_y_continuous(labels = label_number(scale = 1e-6)) +
  nature_theme

# Fig. 1e
ggplot(df, aes(x = MP)) + 
  geom_histogram(aes(fill = after_stat(x)), binwidth = 0.1) +
  scale_fill_gradient(low = mp_gradient_colors["low"], 
                      high = mp_gradient_colors["high"], 
                      guide = "none") +
  labs(y = "Frequency", x = "Modification Probability (%)") +
  scale_y_continuous(labels = comma) +
  nature_theme

# Fig. 1f
ggplot(stack_dataframe_sorted[stack_dataframe_sorted$element != "MT", ], 
       aes(x = element, y = count, fill = modification_probability)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_gradient(low = mp_gradient_colors["low"], high = mp_gradient_colors["high"]) + 
  labs(y = "Number of CpGs", x = "Chromosome", fill = "Modification Probability") +
  scale_y_continuous(labels = comma) +
  nature_theme +
  theme(legend.position = "right",
        legend.key.size = unit(0.25, "cm"),
        legend.title = element_text(size = 8),
        legend.text = element_blank(),
        legend.ticks = element_blank())

# Fig. 1g
ggplot(
  stack_dataframe_sorted[stack_dataframe_sorted$chromosome != "MT", ], 
  aes(x = element, y = count, fill = modification_probability)
) + 
  geom_bar(stat = "identity", position = "fill") + 
  scale_fill_gradient(low = mp_gradient_colors["low"], high = mp_gradient_colors["high"]) +  
  labs(y = "Proportion of CpGs", x = "Chromosome", fill = "Modification Probability") +
  scale_y_continuous(labels = comma) +
  nature_theme +
  theme(legend.position = "right",
        legend.key.size = unit(0.25, "cm"),
        legend.title = element_text(size = 8),
        legend.text = element_blank(),
        legend.ticks = element_blank())

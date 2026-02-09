library(ggplot2)
library(data.table)
library(cowplot)
library(dplyr)
library(scales)
library(ggpubr)
setwd("D:/R/WD/MP_TSS/intro")

### Bar 1. GC ###
df=fread("feature_bar_GC_corrected.txt")
dfh = fread("feature_bar_GC_corrected_hypo.txt")
colnames(df)=c("mp", "count", "chromosome", "type", "id")
colnames(dfh)=c("mp", "count", "chromosome", "type", "id")
df=rbind(df, dfh)

df$id <- factor(df$id, levels = c("Hypothetical", "hg38"))

ggplot() +
  geom_bar(data=df, aes(x=id, y=count, fill=type), stat="identity", position="fill") +
  scale_fill_manual(values = c("GC"="#E7298A", "None"="#BBBBBB"), name="Type") + 
  labs(y="Proportion", x="Genome", size=24) +
  theme(
    #plot.title = element_text(hjust=0.5, size=24),
    panel.background=element_blank(),
    plot.background=element_blank(),
    panel.grid = element_blank(),
    axis.line=element_line(linewidth=.1),
    axis.ticks=element_blank(),
    axis.title=element_text(size=6),
    axis.text=element_text(size=5, color="black"),
    legend.position="none"
  )->gc_genome

### Bar 2. CpG 및 any GC-dimer proportion ###
df=fread("feature_bar_GCdimers_corrected.txt")
dfh = fread("feature_bar_GCdimers_corrected_hypo.txt")
colnames(df)=c("mp", "count", "chromosome", "type", "id")
colnames(dfh)=c("mp", "count", "chromosome", "type", "id")
df=rbind(df, dfh)

# Normalize counts within each id for 100% stacking
df <- df %>%
  group_by(id) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

df$id <- factor(df$id, levels = c("Hypothetical", "hg38"))
df$type <- factor(df$type, levels = c("CpG", "GC", "None"))
df <- df[order(df$type, df$mp), ]

# Create a new column for fill color
df <- df %>%
  mutate(fill_color = case_when(
    type == "None" ~ "#BBBBBB",
    type == "GC"   ~ "#004444",
    type == "CpG"  ~ scales::col_numeric(c("#FFFFFF", "#008888"), domain = range(mp, na.rm = TRUE))(mp),
    TRUE ~ "#000000"  # Fallback color
  ))

# Plot
ggplot(df, aes(x = id, y = prop, fill = fill_color, group = type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_identity() +
  labs(y="Proportion", x="Genome") +
  theme(
    panel.background=element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.y.left=element_line(),
    axis.line.x.bottom=element_line(),
    axis.ticks=element_blank(),
    axis.title=element_text(size=20),
    axis.text=element_text(size=16),
    legend.position="none"
  )->dimer_cpg_genome

### Bar 2. Redo w/o gradation for CpG ###
df=fread("feature_bar_GCdimers_corrected.txt")
dfh = fread("feature_bar_GCdimers_corrected_hypo.txt")
colnames(df)=c("mp", "count", "chromosome", "type", "id")
colnames(dfh)=c("mp", "count", "chromosome", "type", "id")
df=rbind(df, dfh)

# Normalize counts within each id for 100% stacking
df <- df %>%
  group_by(id) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

df$id <- factor(df$id, levels = c("Hypothetical", "hg38"))
df$type <- factor(df$type, levels = c("CpG", "GC", "None"))
df <- df[order(df$type, df$mp), ]

# Create a new column for fill color
df <- df %>%
  mutate(fill_color = case_when(
    type == "None" ~ "#BBBBBB",
    type == "GC"   ~ "#004444",
    type == "CpG"  ~ "#004444",
    TRUE ~ "#000000"  # Fallback color
  ))

# Plot
ggplot(df, aes(x = id, y = prop, fill = fill_color, group = type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_identity() +
  labs(y="Proportion", x="Genome") +
  theme(
    #plot.title = element_text(hjust=0.5, size=18),
    panel.background=element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    axis.line=element_line(linewidth=.1),
    axis.ticks=element_blank(),
    axis.title=element_text(size=6, color="black"),
    axis.text=element_text(size=5, color="black"),
    legend.position="none"
  )->dimer_cpg_genome

### Bar 3. Comparison of 4 dimers ###
df=fread("feature_bar_fourDimers.txt")
colnames(df)=c("mp", "count", "chromosome", "type", "id")

# Create a new column for fill color
df <- df %>%
  mutate(fill_color = case_when(
    type == "None" ~ "#BBBBBB",
    type == "CpC"   ~ "#004444",
    type == "GpC"   ~ "#004444",
    type == "GpG"   ~ "#004444",
    type == "CpG"  ~ scales::col_numeric(c("#FFFFFF", "#008888"), domain = range(mp, na.rm = TRUE))(mp),
    TRUE ~ "#000000"  # Fallback color
  ))

df$type <- factor(df$type, levels = c("CpC", "CpG", "GpC", "GpG"))

# Plot
ggplot(df, aes(x = type, y = count, fill = fill_color, group = type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_identity() +
  labs(y="Length (bp)", x="Dinucleotide") +
  scale_y_continuous(labels = label_comma()) +
  theme(
    panel.background=element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.y.left=element_line(),
    axis.line.x.bottom=element_line(),
    axis.ticks=element_blank(),
    axis.title=element_text(size=5, color="black"),
    axis.text=element_text(size=5, color="black"),
  ) + 
  geom_hline(yintercept = 373972181, linetype = "dotted", linewidth=.1, color = "red")->dimer

### Bar 3. Redo w/o gradation for CpG ###

df=fread("feature_bar_fourDimers.txt")
colnames(df)=c("mp", "count", "chromosome", "type", "id")

# Create a new column for fill color
df <- df %>%
  mutate(fill_color = case_when(
    type == "None" ~ "#BBBBBB",
    type == "CpC"   ~ "#004444",
    type == "GpC"   ~ "#004444",
    type == "GpG"   ~ "#004444",
    type == "CpG"  ~ "#004444",
    TRUE ~ "#000000"  # Fallback color
  ))

df$type <- factor(df$type, levels = c("CpC", "CpG", "GpC", "GpG"))

# Plot
ggplot(df, aes(x = type, y = count, fill = fill_color, group = type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_identity() +
  labs(y="Length (Mbp)", x="Dinucleotide") +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6)) +
  theme(
    panel.background=element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    axis.line=element_line(linewidth=.1),
    axis.ticks=element_blank(),
    axis.title=element_text(size=6, color="black"),
    axis.text=element_text(size=5, color="black")
  ) + 
  geom_hline(yintercept = 373972181, linetype = "dotted", linewidth=.25, color = "red")->dimer

### 4. genome-wide MP distribution ###
df = fread("CG_placedOnly.mp")
colnames(df)=c("Chromosome", "Start", "End", "MP")

genome_distribution = ggplot(df, aes(x = MP)) + 
  geom_histogram(aes(fill = after_stat(x)),binwidth = 0.1) +
  scale_fill_gradient(low = "#FFFFFF", high = "#009999", guide = "none") +
  theme(
    plot.background = element_blank(),
    panel.background = element_rect(fill = "#BBBBBB"),
    panel.grid = element_blank(),
    axis.line.y.left = element_blank(),
    axis.line.x.bottom = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_text(size=6, color="black"),
    axis.text=element_text(size=5, color="black"),
    legend.position = "none"  # Fixed syntax (added comma above)
  ) +
  labs(y = "Frequency", x = "Modification probability (%)") +
  scale_y_continuous(labels = comma)
#ggsave("genome_distribution.png", genome_distribution, units="mm", height=90, width=90)

### 5, 6. Chromosomal modification distribution ###
stack_dataframe_sorted = fread("stack_dataframe_sorted.txt") # A data with columns denoting: mod. prob. (in order 100~0), count (of sites with that mod. prob.), element (category)
colnames(stack_dataframe_sorted) = c("modification_probability", "count", "element")

# Set factor levels
order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM")
stack_dataframe_sorted$element <- factor(stack_dataframe_sorted$element, levels = order)
levels(stack_dataframe_sorted$element) <- paste(c(1:22, "X", "Y", "MT"))
stack_dataframe_sorted <- stack_dataframe_sorted[order(stack_dataframe_sorted$modification_probability), ]

# Proportional height
chromosomal_gradation=ggplot(stack_dataframe_sorted[element!="MT",], aes(x=element, y=count, fill=modification_probability)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_gradient(low="#FFFFFF", high="#009999") + 
  theme(axis.ticks = element_blank(),
        plot.background=element_blank(),
        panel.background=element_rect(fill="#BBBBBB"),
        panel.grid=element_blank(),
        legend.position = "none",
        legend.background = element_blank(),
        legend.key.size = unit(.25, 'cm'),
        legend.title=element_text(size=8),
        legend.text=element_blank(),
        legend.ticks=element_blank(),
        axis.title=element_text(size=6, color="black"),
        axis.text=element_text(size=5, color="black")) + 
  labs(y="Number of CpGs", x="Chromosome", fill="Modification probability") + 
  scale_y_continuous(labels = comma)

# Equal height
chromosomal_gradation_100=ggplot(stack_dataframe_sorted[element!="MT",], aes(x=element, y=count, fill=modification_probability)) + 
  geom_bar(stat = "identity", position = "fill") + 
  scale_fill_gradient(low="#FFFFFF", high="#009999") + 
  theme(axis.ticks = element_blank(),
        plot.background=element_blank(),
        panel.background=element_rect(fill="#BBBBBB"),
        panel.grid=element_blank(),
        legend.position = "none",
        legend.background = element_blank(),
        legend.key.size = unit(.25, 'cm'),
        legend.title=element_text(size=5),
        legend.text=element_blank(),
        legend.ticks=element_blank(),
        axis.title=element_text(size=6, color="black"),
        axis.text=element_text(size=5, color="black")) + 
  labs(y="Proportion of CpGs", x="Chromosome", fill="Modification probability") + 
  scale_y_continuous(labels = comma)

# ggsave("chromosomal_gradation.png", chromosomal_gradation, units="mm", height=90, width=180)

### Legend ###
df=rbind(fread("feature_bar_GC_corrected.txt"), fread("feature_bar_GCdimers_corrected.txt"), fread("feature_bar_fourDimers.txt"))
colnames(df)=c("mp", "count", "chromosome", "type", "id")

for_legend <- ggplot() +
  geom_bar(data = df[is.na(mp)], aes(x = id, y = count, fill = type), stat = "identity") +
  scale_fill_manual(values=c("GC" = "#E7298A", "CpC" = "#004444"), name=" ", breaks=c("GC", "CpC"), labels=c("GC"="GC content", "CpC"="GC-comprised dimer")) +
  theme(legend.text = element_text(size=6), legend.key.size=unit(.5, "lines")) +
  # ↑ above line manipulates legend. values→color, name→title, breaks→order, labels→item names
  guides(fill = guide_legend(order = 1)) +  # ← Discrete legend first
  ggnewscale::new_scale_fill() +
  geom_bar(data = df[!is.na(mp)], aes(x = id, y = count, fill = mp), stat = "identity") +
  scale_fill_gradientn(colors = c("#FFFFFF", "#008888"),name="MP (0-100%)") + 
  guides(fill = guide_colourbar(label=FALSE, order=2)) +
  theme(legend.position = "top",
        legend.ticks = element_blank(),
        legend.background = element_blank(),
        legend.title=element_text(size=6))

legend <- as_ggplot(get_legend(for_legend));legend

### Combine ###
#combined = plot_grid(gc_genome, legend, dimer_cpg_genome, dimer, ncol=2, labels = c("a", "", "b", "c"), label_size=12, rel_widths=c(1, 1, 1, 1), align = "hv", axis = "bt")
combined_top = plot_grid(gc_genome, dimer_cpg_genome, dimer, genome_distribution, ncol=2, labels=c("", "", "", ""), label_size=20, rel_widths=c(1,1,1,1), align="hv", axis="bt")
combined_bot = plot_grid(chromosomal_gradation, chromosomal_gradation_100, legend, ncol=1,labels=c("","", ""), label_size=20, rel_heights=c(10,10,1), align="v", axis="lr")
combined_all = plot_grid(combined_top, combined_bot, ncol=1, rel_heights=c(20,21))

### Recombine ###
combined_top = plot_grid(gc_genome, dimer_cpg_genome, dimer, ncol=3, labels=c("a","b","c"), label_size=8)
combined_mid_left = plot_grid(genome_distribution, ncol=1, labels=c("d"), label_size=8)
combined_mid_right = plot_grid(chromosomal_gradation, chromosomal_gradation_100, ncol=1, align="v", axis="lr", labels=c("e","f"), label_size=8)
combined_mid = plot_grid(combined_mid_left, combined_mid_right, ncol=2)
combined_bot = plot_grid(legend, ncol=1)
combined_all = plot_grid(labels=c(""), label_size=24, combined_top, combined_mid, combined_bot, ncol=1, rel_heights=c(10,15,1))

### Re-recombine ###
fig3=align_plots(dimer + theme(plot.margin = margin(5.5,0,5.5,5.5)))
combined_top = plot_grid(gc_genome, dimer_cpg_genome, fig3[[1]], ncol=3, labels=c("a","b","c"), label_size=8, rel_widths=c(1,1,1))
fig4=align_plots(genome_distribution + theme(plot.margin = margin(5.5,0,5.5,5.5)))
combined_mid = plot_grid(fig4[[1]], ncol=1, labels=c("d"), label_size=8)
chromosomal_gradations = align_plots(chromosomal_gradation + theme(plot.margin = margin(5.5,0,5.5,5.5)), chromosomal_gradation_100 + theme(plot.margin = margin(5.5,0,5.5,5.5)), align="v", axis="lr")
combined_bot = plot_grid(chromosomal_gradations[[1]], chromosomal_gradations[[2]], ncol=1, labels=c("e", "f"), label_size=8)

combined_all = plot_grid(combined_top, combined_mid, combined_bot, legend, ncol=1, rel_heights=c(12.5,12.5,25,1))

height_total=12.5+17.5+25+1
width_total=30

ggsave("combined_intro_shaded_poster.png", combined_all, height=130, width=85, units = "mm")
ggsave("combined_intro_shaded_poster.pdf", combined_all, height=130, width=85, units = "mm")

### 
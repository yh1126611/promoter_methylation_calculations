library(ggplot2)
library(data.table)
library(cowplot)
library(dplyr)
library(scales)
library(grid)
library(ggpubr)
library(ggnewscale)
library(svgparser)
library(png)

setwd("D:/R/WD/MP_TSS/tvds")

human_svg=svgparser::read_svg("D:/R/WD/MP_TSS/fig7/svgs/human_sized.svg")

wholeY=158.9
wholeX=48.2
justY=66.3
justX=27.9

ggplot(data.frame(x=c(0.5,0.5),y=c(0.5,0.5)))+
  #geom_rect(xmin = (wholeX-justX)/2, xmax = (wholeX-justX)/2+justX, ymin = (wholeY-justY)/2, ymax = (wholeY-justY)/2+justY)+
  annotation_custom(human_svg, xmin=0, xmax=1, ymin=0, ymax=1)+
  geom_point(mapping=aes(x=x,y=y), color="transparent")+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))+
  scale_x_continuous(limits=c(0,0.2), expand=c(0,0))+
  theme(panel.background=element_blank(),
        plot.background=element_blank(),
        panel.grid = element_blank(),
        axis.line=element_blank(),
        axis.line.x=element_blank(),
        axis.line.y = element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank()
        )->just_human;just_human


compute_tvd_pvalues <- function(wg_file, element_file, n_perm = 1000) {
  library(data.table)
  compute_tvd <- function(x, y, breaks = 99) { # Nested function
    if (length(y) == 0) { warning("Empty y vector - returning NA"); return(NA) }
    if (anyNA(y)) { warning("NA values in y - removing them"); y <- y[!is.na(y)] }
    y_range <- range(y); x_range <- range(x)
    combined_range <- if (y_range[1] < x_range[1] || y_range[2] > x_range[2]) {
      warning("Y values outside X range - adjusting breaks"); range(c(x, y))
    } else x_range
    breaks_vec <- seq(combined_range[1], combined_range[2], length.out = breaks + 1)
    hist_x <- hist(x, breaks = breaks_vec, plot = FALSE)
    hist_y <- hist(y, breaks = breaks_vec, plot = FALSE)
    p_x <- hist_x$counts / length(x); p_y <- hist_y$counts / length(y)
    0.5 * sum(abs(p_x - p_y))
  }
  wg_distrbtn <- fread(wg_file)[, V4]; mp <- fread(element_file) # Load population & sample
  p_values <- data.frame(spcdist = integer(), p_value = numeric()) # Create empty data frame “p_values”
  for (spcdist in -10000:10000) {
    spcdist_distrbtn <- mp[V2 == spcdist, V3]
    if (length(spcdist_distrbtn) == 0) { warning(paste("No data for spcdist", spcdist)); next }
    obs_tvd <- compute_tvd(wg_distrbtn, spcdist_distrbtn); n <- length(spcdist_distrbtn)
    perm_tvd <- replicate(n_perm, {
      shuffled <- sample(wg_distrbtn, size = 2 * n)
      compute_tvd(shuffled[1:n], shuffled[(n + 1):(2 * n)])
    })
    p_values <- rbind(p_values, data.frame(spcdist = spcdist, p_value = mean(perm_tvd >= obs_tvd)))
  }
  return(p_values)
}

plot_pvalue_heatmap <- function(p_values_df, title) {
  ggplot_heatmap_data <- p_values_df; ggplot_heatmap_data$one <- "1"
  ggplot(ggplot_heatmap_data, aes(spcdist, one, fill = p_value, height = 1)) +
    geom_tile() +
    scale_x_continuous(expand=c(0,0))+
    scale_fill_gradientn(colours=c("#FF4B33", "white", "#BBBBBB"), values=rescale(c(0, 0.05, 1)), limits=c(0, 1), breaks=c(0, 0.05, 1), labels=c("0", "5", "100"), name="p-value (%)") + labs(title=title) +
    theme(plot.background=element_blank(), panel.background=element_blank(), panel.grid=element_blank(), axis.ticks=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), plot.title = element_text(size = 5, vjust=-2), legend.position="none")
  pvalue_heatmap <- ggplot(ggplot_heatmap_data, aes(spcdist, one, fill = p_value, height = 1)) +
    geom_tile() +
    scale_x_continuous(expand=c(0,0))+
    scale_fill_gradientn(colours=c("#FF4B33", "white", "#BBBBBB"), values=rescale(c(0, 0.05, 1)), limits=c(0, 1), breaks=c(0, 0.05, 1), labels=c("0", "5", "100"), name="p-value (%)") + labs(title=title) +
    theme(plot.background=element_blank(), panel.background=element_blank(), panel.grid=element_blank(), axis.ticks=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), plot.title = element_text(size = 5, vjust=-2), legend.position="none")
  return(pvalue_heatmap)
}

data_preproc_median = function(df){
  df[V4 == "-", V2 := V2 * -1] # If reverse (-) strand, flip the coordinates
  df = data.frame(summarise(group_by(df, V2),med=median(V3),SD=sd(V3)))
  colnames(df)=c("Distance","Median","SD")
  #df = summarise(group_by(mutate(df, x_bin = round(Distance / 100) * 100),x_bin),Average = mean(Average,na.rm=TRUE))
  # ↑ Downsamples the data with BIN size of 100 bp. Change the value at: “x_bin = round(Distance / n) * n” to change the BIN size.
  colnames(df)=c("Distance","Median")
  return(df)
}

data_preproc_mean = function(df){
  df[V4 == "-", V2 := V2 * -1] # If reverse (-) strand, flip the coordinates
  df = data.frame(summarise(group_by(df, V2),avg=mean(V3),SD=sd(V3)))
  colnames(df)=c("Distance","Mean","SD")
  #df = summarise(group_by(mutate(df, x_bin = round(Distance / 100) * 100),x_bin),Average = mean(Average,na.rm=TRUE))
  # ↑ Downsamples the data with BIN size of 100 bp. Change the value at: “x_bin = round(Distance / n) * n” to change the BIN size.
  colnames(df)=c("Distance","Mean")
  return(df)
}

pvalues_TSS = compute_tvd_pvalues("CG_placedOnly.mp", "MP_TSS_hg38.tsv") #⌛
save.image(file = "pvalues6.RData")
pvalues_TES = compute_tvd_pvalues("CG_placedOnly.mp", "MP_TES_hg38.tsv") #⌛
save.image(file = "pvalues5.RData")
pvalues_promoter_ENCODE = compute_tvd_pvalues("CG_placedOnly.mp", "MP_promoter_ENCODE.tsv") #⌛
save.image(file = "pvalues4.RData")
pvalues_promoter_NCBI = compute_tvd_pvalues("CG_placedOnly.mp", "MP_promoter_hg38.tsv") #⌛
save.image(file = "pvalues3.RData")
pvalues_enhancer_NCBI = compute_tvd_pvalues("CG_placedOnly.mp", "MP_enhancer_hg38.tsv") #⌛
save.image(file = "pvalues2.RData")
pvalues_silencer_NCBI = compute_tvd_pvalues("CG_placedOnly.mp", "MP_silencer_hg38.tsv") #⌛
save.image(file = "pvalues1.RData")
pvalues_enhancer_ENCODE = compute_tvd_pvalues("CG_placedOnly.mp", "MP_enhancer_ENCODE.tsv") #⌛
save.image(file = "pvalues0.RData")
pvalues_transcriptTSS = compute_tvd_pvalues("CG_placedOnly.mp", "MP_transcriptTSS_hg38.tsv") #⌛
save.image(file = "pvalues-1.RData")
pvalues_transcriptTES = compute_tvd_pvalues("CG_placedOnly.mp", "MP_transcriptTES_hg38.tsv") #⌛
save.image(file = "pvalues-2.RData")
pvalues_nonGeneTranscriptTSS = compute_tvd_pvalues("CG_placedOnly.mp", "MP_nonGeneTranscriptTSS_hg38.tsv") #⌛
save.image(file = "pvalues-3.RData")
pvalues_nonGeneTranscriptTES = compute_tvd_pvalues("CG_placedOnly.mp", "MP_nonGeneTranscriptTES_hg38.tsv") #⌛
save.image(file = "pvalues-4.RData")

heatmap_TSS=plot_pvalue_heatmap(pvalues_TSS, "GSS")
heatmap_TTS=plot_pvalue_heatmap(pvalues_TES, "GTS")
heatmap_promoter_ENCODE=plot_pvalue_heatmap(pvalues_promoter_ENCODE, "Promoter")
heatmap_promoter_NCBI=plot_pvalue_heatmap(pvalues_promoter_NCBI, "Promoter")
heatmap_enhancer_ENCODE=plot_pvalue_heatmap(pvalues_enhancer_ENCODE, "Enhancer")
heatmap_enhancer_NCBI=plot_pvalue_heatmap(pvalues_enhancer_NCBI, "Enhancer")
heatmap_silencer_NCBI=plot_pvalue_heatmap(pvalues_silencer_NCBI, "Silencer")

heatmap_transcriptTSS=plot_pvalue_heatmap(pvalues_transcriptTSS, "TSS (GSS + ATSS)")
heatmap_transcriptTTS=plot_pvalue_heatmap(pvalues_transcriptTES, "TTS (GTS + ATTS)")
heatmap_nonGeneTranscriptTSS=plot_pvalue_heatmap(pvalues_nonGeneTranscriptTSS, "ATSS")
heatmap_nonGeneTranscriptTTS=plot_pvalue_heatmap(pvalues_nonGeneTranscriptTES, "ATTS")



# Transcript TSS

df = data_preproc_median(fread("MP_transcriptTSS_hg38.tsv"))
ggplot(df, aes(x = Distance, y = Median)) + 
  geom_line(linewidth=.1, alpha=1) +
  # geom_rect(aes(ymax=100,ymin=0,xmin=-2500,xmax=2500),fill="transparent", color="red", linewidth=.1, linetype="dashed") +
  #geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,25,50,75,100), expand=c(0,0)) +
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Median MP (%)", x = "Distance from TTS (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = c("-10", "0", "10"), limits=c(-10000,10000), breaks=c(-10000,0,10000), expand=c(0,0)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 5, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_transcriptTSS_median

df = data_preproc_mean(fread("MP_transcriptTSS_hg38.tsv"))
ggplot(df, aes(x = Distance, y = Mean)) + 
  geom_line(linewidth=.1, alpha=1) +
  #geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,25,50,75,100), expand=c(0,0)) +
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Average MP (%)", x = "Distance from TTS (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = c("-10", "0", "10"), limits=c(-10000,10000), breaks=c(-10000,0,10000), expand=c(0,0)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 5, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_transcriptTSS_mean

# Transcript TTS

df = data_preproc_median(fread("MP_transcriptTES_hg38.tsv"))
ggplot(df, aes(x = Distance, y = Median)) + 
  geom_line(linewidth=.1, alpha=1) +
  #geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,25,50,75,100), expand=c(0,0)) +
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Median MP (%)", x = "Distance from TTS (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = c("-10", "0", "10"), limits=c(-10000,10000), breaks=c(-10000,0,10000), expand=c(0,0)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 5, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_transcriptTTS_median

df = data_preproc_mean(fread("MP_transcriptTES_hg38.tsv"))
ggplot(df, aes(x = Distance, y = Mean)) + 
  geom_line(linewidth=.1, alpha=1) +
  #geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,25,50,75,100), expand=c(0,0)) +
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Average MP (%)", x = "Distance from TTS (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = c("-10", "0", "10"), limits=c(-10000,10000), breaks=c(-10000,0,10000), expand=c(0,0)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 5, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_transcriptTTS_mean

# Non-gene transcript TSS

df = data_preproc_median(fread("MP_nonGeneTranscriptTSS_hg38.tsv"))
ggplot(df, aes(x = Distance, y = Median)) + 
  geom_line(linewidth=.1, alpha=1) +
  #geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,25,50,75,100), expand=c(0,0)) +
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Median MP (%)", x = "Distance from TTS (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = c("-10", "0", "10"), limits=c(-10000,10000), breaks=c(-10000,0,10000), expand=c(0,0)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 5, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_nonGeneTranscriptTSS_median

df = data_preproc_mean(fread("MP_nonGeneTranscriptTSS_hg38.tsv"))
ggplot(df, aes(x = Distance, y = Mean)) + 
  geom_line(linewidth=.1, alpha=1) +
  #geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,25,50,75,100), expand=c(0,0)) +
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Average MP (%)", x = "Distance from TTS (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = c("-10", "0", "10"), limits=c(-10000,10000), breaks=c(-10000,0,10000), expand=c(0,0)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 5, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_nonGeneTranscriptTSS_mean

# Non-gene transcript TTS

df = data_preproc_median(fread("MP_nonGeneTranscriptTES_hg38.tsv"))
ggplot(df, aes(x = Distance, y = Median)) + 
  geom_line(linewidth=.1, alpha=1) +
  #geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,25,50,75,100), expand=c(0,0)) +
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Median MP (%)", x = "Distance from TTS (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = c("-10", "0", "10"), limits=c(-10000,10000), breaks=c(-10000,0,10000), expand=c(0,0)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 5, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_nonGeneTranscriptTTS_median

df = data_preproc_mean(fread("MP_nonGeneTranscriptTES_hg38.tsv"))
ggplot(df, aes(x = Distance, y = Mean)) + 
  geom_line(linewidth=.1, alpha=1) +
  #geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,25,50,75,100), expand=c(0,0)) +
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Average MP (%)", x = "Distance from TTS (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = c("-10", "0", "10"), limits=c(-10000,10000), breaks=c(-10000,0,10000), expand=c(0,0)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 5, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_nonGeneTranscriptTTS_mean

# TSS

df = data_preproc_median(fread("MP_TSS_hg38.tsv"))
ggplot(df, aes(x = Distance, y = Median)) + 
  geom_line(linewidth=.1, alpha=1) +
  #geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,25,50,75,100), expand=c(0,0)) +
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Median MP (%)", x = "Distance from TSS (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = c("-10", "0", "10"), limits=c(-10000,10000), breaks=c(-10000,0,10000), expand=c(0,0)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 5, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_TSS_median

df = data_preproc_mean(fread("MP_TSS_hg38.tsv"))
ggplot(df, aes(x = Distance, y = Mean)) + 
  geom_line(linewidth=.1, alpha=1) +
  #geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,25,50,75,100), expand=c(0,0)) + 
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Average MP (%)", x = "Distance from TSS (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = c("-10", "0", "10"), limits=c(-10000,10000), breaks=c(-10000,0,10000), expand=c(0,0)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 5, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_TSS_mean

# TTS

df = data_preproc_median(fread("MP_TES_hg38.tsv"))
ggplot(df, aes(x = Distance, y = Median)) + 
  geom_line(linewidth=.1, alpha=1) +
  #geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,25,50,75,100), expand=c(0,0)) +
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Median MP (%)", x = "Distance from TTS (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = c("-10", "0", "10"), limits=c(-10000,10000), breaks=c(-10000,0,10000), expand=c(0,0)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 5, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_TTS_median

df = data_preproc_mean(fread("MP_TES_hg38.tsv"))
ggplot(df, aes(x = Distance, y = Mean)) + 
  geom_line(linewidth=.1, alpha=1) +
  #geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,25,50,75,100), expand=c(0,0)) +
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Average MP (%)", x = "Distance from TTS (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = c("-10", "0", "10"), limits=c(-10000,10000), breaks=c(-10000,0,10000), expand=c(0,0)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 5, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_TTS_mean

# Promoter NCBI

df = data_preproc_median(fread("MP_promoter_hg38.tsv"))
ggplot(df, aes(x = Distance, y = Median)) + 
  geom_line(linewidth=.1, alpha=1) +
  geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,25,50,75,100), expand=c(0,0)) +
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Median MP (%)", x = "Distance from center (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = c("-10", "0", "10"), limits=c(-10000,10000), breaks=c(-10000,0,10000), expand=c(0,0)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 5, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_promoterNCBI_median

df = data_preproc_mean(fread("MP_promoter_hg38.tsv"))
ggplot(df, aes(x = Distance, y = Mean)) + 
  geom_line(linewidth=.1, alpha=1) +
  geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,25,50,75,100), expand=c(0,0)) +
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Average MP (%)", x = "Distance from center (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = c("-10", "0", "10"), limits=c(-10000,10000), breaks=c(-10000,0,10000), expand=c(0,0)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 5, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_promoterNCBI_mean

# Promoter ENCODE

df = data_preproc_median(fread("MP_promoter_ENCODE.tsv"))
ggplot(df, aes(x = Distance, y = Median)) + 
  geom_line(linewidth=.1, alpha=1) +
  #geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,25,50,75,100), expand=c(0,0)) + 
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Median MP (%)", x = "Distance from center (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = c("-10", "0", "10"), limits=c(-10000,10000), breaks=c(-10000,0,10000), expand=c(0,0)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 5, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_promoterENCODE_median

df = data_preproc_mean(fread("MP_promoter_ENCODE.tsv"))
ggplot(df, aes(x = Distance, y = Mean)) + 
  geom_line(linewidth=.1, alpha=1) +
  #geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,25,50,75,100), expand=c(0,0)) +
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Average MP (%)", x = "Distance from center (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = c("-10", "0", "10"), limits=c(-10000,10000), breaks=c(-10000,0,10000), expand=c(0,0)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 5, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_promoterENCODE_mean

# Enhancer NCBI

df = data_preproc_median(fread("MP_enhancer_hg38.tsv"))
ggplot(df, aes(x = Distance, y = Median)) + 
  geom_line(linewidth=.1, alpha=1) +
  #geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,25,50,75,100), expand=c(0,0)) +
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Median MP (%)", x = "Distance from center (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = c("-10", "0", "10"), limits=c(-10000,10000), breaks=c(-10000,0,10000), expand=c(0,0)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 5, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_enhancerNCBI_median

df = data_preproc_mean(fread("MP_enhancer_hg38.tsv"))
ggplot(df, aes(x = Distance, y = Mean)) + 
  geom_line(linewidth=.1, alpha=1) +
  #geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,25,50,75,100), expand=c(0,0)) +
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Average MP (%)", x = "Distance from center (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = c("-10", "0", "10"), limits=c(-10000,10000), breaks=c(-10000,0,10000), expand=c(0,0)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 5, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_enhancerNCBI_mean


# Enhancer ENCODE

df = data_preproc_median(fread("MP_enhancer_ENCODE.tsv"))
ggplot(df, aes(x = Distance, y = Median)) + 
  geom_line(linewidth=.1, alpha=1) +
  #geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,25,50,75,100), expand=c(0,0)) +
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Median MP (%)", x = "Distance from center (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = c("-10", "0", "10"), limits=c(-10000,10000), breaks=c(-10000,0,10000), expand=c(0,0)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 5, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_enhancerENCODE_median

df = data_preproc_mean(fread("MP_enhancer_ENCODE.tsv"))
ggplot(df, aes(x = Distance, y = Mean)) + 
  geom_line(linewidth=.1, alpha=1) +
  #geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,25,50,75,100), expand=c(0,0)) +
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Average MP (%)", x = "Distance from center (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = c("-10", "0", "10"), limits=c(-10000,10000), breaks=c(-10000,0,10000), expand=c(0,0)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 5, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_enhancerENCODE_mean


# Silencer NCBI

df = data_preproc_median(fread("MP_silencer_hg38.tsv"))
ggplot(df, aes(x = Distance, y = Median)) + 
  geom_line(linewidth=.1, alpha=1) +
  #geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,25,50,75,100), expand=c(0,0)) +
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Median MP (%)", x = "Distance from center (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = c("-10", "0", "10"), limits=c(-10000,10000), breaks=c(-10000,0,10000), expand=c(0,0)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 5, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_silencerNCBI_median

df = data_preproc_mean(fread("MP_silencer_hg38.tsv"))
ggplot(df, aes(x = Distance, y = Mean)) + 
  geom_line(linewidth=.1, alpha=1) +
  #geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,25,50,75,100), expand=c(0,0)) +
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Average MP (%)", x = "Distance from center (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = c("-10", "0", "10"), limits=c(-10000,10000), breaks=c(-10000,0,10000), expand=c(0,0)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 5, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_silencerNCBI_mean


# Blank ★ 여기부터 시작하면 됨.
blank_plot <- ggplot() + theme_void()
blank <- ggplot() + theme_void()

## Now let's align & combine all those plots into a multi-panel figure ##
# 1. Align plots:

aligned_transcriptTSS <- align_plots(
  heatmap_transcriptTSS + theme(plot.margin = margin(0,0,0,10)),
  line_transcriptTSS_median + theme(plot.margin = margin(0,0,5,5)),
  line_transcriptTSS_mean + theme(plot.margin = margin(0,0,5,5)),
  align = "v",
  axis = "lr"
)

aligned_transcriptTTS <- align_plots(
  heatmap_transcriptTTS + theme(plot.margin = margin(0,0,0,10)),
  line_transcriptTTS_median + theme(plot.margin = margin(0,0,5,5)),
  line_transcriptTTS_mean + theme(plot.margin = margin(0,0,5,5)),
  align = "v",
  axis = "lr"
)

aligned_nonGeneTranscriptTSS <- align_plots(
  heatmap_nonGeneTranscriptTSS + theme(plot.margin = margin(0,0,0,10)),
  line_nonGeneTranscriptTSS_median + theme(plot.margin = margin(0,0,5,5)),
  line_nonGeneTranscriptTSS_mean + theme(plot.margin = margin(0,0,5,5)),
  align = "v",
  axis = "lr"
)

aligned_nonGeneTranscriptTTS <- align_plots(
  heatmap_nonGeneTranscriptTTS + theme(plot.margin = margin(0,0,0,10)),
  line_nonGeneTranscriptTTS_median + theme(plot.margin = margin(0,0,5,5)),
  line_nonGeneTranscriptTTS_mean + theme(plot.margin = margin(0,0,5,5)),
  align = "v",
  axis = "lr"
)

aligned_TSS <- align_plots(
  heatmap_TSS + theme(plot.margin = margin(0,0,0,10)),
  line_TSS_median + theme(plot.margin = margin(0,0,5,5)),
  line_TSS_mean + theme(plot.margin = margin(0,0,5,5)),
  align = "v",
  axis = "lr"
)

aligned_TTS <- align_plots(
  heatmap_TTS + theme(plot.margin = margin(0,0,0,10)),
  line_TTS_median + theme(plot.margin = margin(0,0,5,5)),
  line_TTS_mean + theme(plot.margin = margin(0,0,5,5)),
  align = "v",
  axis = "lr"
)

aligned_promoterNCBI <- align_plots(
  heatmap_promoter_NCBI + theme(plot.margin = margin(0,0,0,10)),
  line_promoterNCBI_median + theme(plot.margin = margin(0,0,5,5)),
  line_promoterNCBI_mean + theme(plot.margin = margin(0,0,5,5)),
  align = "v",
  axis = "lr"
)

aligned_promoterENCODE <- align_plots(
  heatmap_promoter_ENCODE + theme(plot.margin = margin(0,0,0,10)),
  line_promoterENCODE_median + theme(plot.margin = margin(0,0,5,5)),
  line_promoterENCODE_mean + theme(plot.margin = margin(0,0,5,5)),
  align = "v",
  axis = "lr"
)

aligned_enhancerNCBI <- align_plots(
  heatmap_enhancer_NCBI + theme(plot.margin = margin(0,0,0,10)),
  line_enhancerNCBI_median + theme(plot.margin = margin(0,0,5,5)),
  line_enhancerNCBI_mean + theme(plot.margin = margin(0,0,5,5)),
  align = "v",
  axis = "lr"
)

aligned_enhancerENCODE <- align_plots(
  heatmap_enhancer_ENCODE + theme(plot.margin = margin(0,0,0,10)),
  line_enhancerENCODE_median + theme(plot.margin = margin(0,0,5,5)),
  line_enhancerENCODE_mean + theme(plot.margin = margin(0,0,5,5)),
  align = "v",
  axis = "lr"
)

aligned_silencerNCBI <- align_plots(
  heatmap_silencer_NCBI + theme(plot.margin = margin(0,0,0,10)),
  line_silencerNCBI_median + theme(plot.margin = margin(0,0,5,5)),
  line_silencerNCBI_mean + theme(plot.margin = margin(0,0,5,5)),
  align = "v",
  axis = "lr"
)

save.image("D:/R/WD/MP_TSS/tvds/pvalues_rightbeforecombine.RData")
load("D:/R/WD/MP_TSS/tvds/pvalues_rightbeforecombine.RData")

# 2. Combine plots:
median_label = ggdraw() + draw_grob(textGrob("          Median MP (%)", rot = 90, y=.5, x=1, gp = gpar(fontsize = 5, fontface = "plain")))
mean_label = ggdraw() + draw_grob(textGrob("          Mean MP (%)", rot = 90, y=.5, x=1, gp = gpar(fontsize = 5, fontface = "plain")))

NCBI=plot_grid(blank_plot, aligned_promoterNCBI[[1]], aligned_enhancerNCBI[[1]], aligned_silencerNCBI[[1]],
               mean_label, aligned_promoterNCBI[[3]], aligned_enhancerNCBI[[3]], aligned_silencerNCBI[[3]],
               median_label, aligned_promoterNCBI[[2]], aligned_enhancerNCBI[[2]], aligned_silencerNCBI[[2]],
               ncol=4, rel_heights=c(6,15,15), rel_widths=c(1,15,15,15))
NCBI_halfTitled=plot_grid(blank_plot, NCBI,
                          labels=c("RefSeq regulator annotations"), label_size=6, label_fontface="bold", label_x=0.5, hjust=0.5, label_y=0.5, vjust=0.5,
                          ncol=1, rel_heights=c(2, 37.5))
NCBI_titled=plot_grid(NCBI_halfTitled, blank_plot, blank_plot,
                      labels=c("", "Distance from center (kb)"), label_size=5, label_fontface="plain", label_x=0.5, hjust=0.5, label_y=0.5, vjust=0.5,
                      ncol=1, rel_heights=c(2+37.5 ,1,.5))

ENCODE=plot_grid(blank_plot, aligned_promoterENCODE[[1]], aligned_enhancerENCODE[[1]],
               mean_label, aligned_promoterENCODE[[3]], aligned_enhancerENCODE[[3]],
               median_label, aligned_promoterENCODE[[2]], aligned_enhancerENCODE[[2]],
               ncol=3, rel_heights=c(6,15,15), rel_widths=c(1,15,15))
ENCODE_halfTitled=plot_grid(blank_plot, ENCODE,
                          labels=c("ENCODE regulator annotations"), label_size=6, label_fontface="bold", label_x=0.5, hjust=0.5, label_y=0.5, vjust=0.5,
                          ncol=1, rel_heights=c(2, 37.5))
ENCODE_titled=plot_grid(ENCODE_halfTitled, blank_plot, blank_plot,
                      labels=c("", "Distance from center (kb)"), label_size=5, label_fontface="plain", label_x=0.5, hjust=0.5, label_y=0.5, vjust=0.5,
                      ncol=1, rel_heights=c(2+37.5 ,1,.5))

RefSeq1=plot_grid(blank_plot, aligned_TSS[[1]], aligned_TTS[[1]],
                 mean_label, aligned_TSS[[3]], aligned_TTS[[3]],
                 median_label, aligned_TSS[[2]], aligned_TTS[[2]],
                 ncol=3, rel_heights=c(6,15,15), rel_widths=c(1,15,15))
RefSeq1_titled=plot_grid(RefSeq1, blank_plot, blank_plot,
                      labels=c("", "Distance from site (kb)"), label_size=5, label_fontface="plain", label_x=0.5, hjust=0.5, label_y=0.5, vjust=0.5,
                      ncol=1, rel_heights=c(37.5 ,1,.5))

RefSeq2=plot_grid(blank_plot, aligned_nonGeneTranscriptTSS[[1]], aligned_nonGeneTranscriptTTS[[1]],
                  mean_label, aligned_nonGeneTranscriptTSS[[3]], aligned_nonGeneTranscriptTTS[[3]],
                  median_label, aligned_nonGeneTranscriptTSS[[2]], aligned_nonGeneTranscriptTTS[[2]],
                  ncol=3, rel_heights=c(6,15,15), rel_widths=c(1,15,15))
RefSeq2_titled=plot_grid(RefSeq2, blank_plot, blank_plot,
                         labels=c("", "Distance from site (kb)"), label_size=5, label_fontface="plain", label_x=0.5, hjust=0.5, label_y=0.5, vjust=0.5,
                         ncol=1, rel_heights=c(37.5 ,1,.5))

RefSeq3=plot_grid(blank_plot, aligned_transcriptTSS[[1]], aligned_transcriptTTS[[1]],
                 mean_label, aligned_transcriptTSS[[3]], aligned_transcriptTTS[[3]],
                 median_label, aligned_transcriptTSS[[2]], aligned_transcriptTTS[[2]],
                 ncol=3, rel_heights=c(6,15,15), rel_widths=c(1,15,15))
RefSeq3_titled=plot_grid(RefSeq3, blank_plot, blank_plot,
                        labels=c("", "Distance from site (kb)"), label_size=5, label_fontface="plain", label_x=0.5, hjust=0.5, label_y=0.5, vjust=0.5,
                        ncol=1, rel_heights=c(37.5 ,1,.5))

plot_legendOnly_heatmap <- function(p_values_df){
  ggplot_heatmap_data <- p_values_df; ggplot_heatmap_data$one <- "1"
  ggplot(data.frame(spcdist=1,p_value=1, one="1"), aes(spcdist, one, fill=p_value, height=1)) +
    geom_tile() + 
    scale_fill_gradientn(colours = c("#FF4B33", "white", "#BBBBBB"), values = scales::rescale(c(0, 0.05, 1)), limits = c(0, 1), breaks = c(1, 0.05), labels = c("1", "0.05"), name = "P-value") +
    guides(fill = guide_colourbar(label=TRUE)) +
    theme(legend.position = "top",
      legend.ticks = element_blank(),
      legend.background = element_blank(),
      legend.title.position = "top",
      #legend.text.position = "bottom",
      legend.title=element_text(size=5),
      legend.text=element_text(size=5),
      legend.key.width = unit(4, "mm"), legend.key.height=unit(4, "mm")
    )}
# Compute only legend
heatmap_fakeForLegend = plot_legendOnly_heatmap(data.frame(spcdist=1, p_value=1))
heatmap_legend = as_ggplot(get_legend(heatmap_fakeForLegend)); heatmap_legend

#human_heatmap_legend = align_plots(just_human + theme(plot.margin = margin(0,0,0,0)), heatmap_legend + theme(plot.margin = margin(0,0,0,0)), align="v", axis="c")
#plot_grid(human_heatmap_legend[[1]], human_heatmap_legend[[2]], ncol=2)

plot_grid(blank, just_human, blank, ncol=3, rel_widths=c(1,1,1), align="h", axis="tb")->kindof_legend
plot_grid(blank, kindof_legend, heatmap_legend, blank, ncol=1, rel_heights=c(1,1,1,1), align="v", axis="lr")->real_legend

row1=plot_grid(NCBI_titled, blank_plot, ENCODE_titled, real_legend,
               align="v",axis="r",
               labels=c("a", "", "b"), label_size=8,
               ncol=4, rel_widths=c(1+15+15+15, 2, 1+15+15, 20))
# The sum of rel_width of row 2 (↓) is 99. So adjust the proportion of the two blank_plots enclosing heatmap_legend to make the sum of rel_width of row 1 99 as well.

row2=plot_grid(RefSeq1_titled, blank_plot, RefSeq2_titled, blank_plot, RefSeq3_titled, blank_plot,
               labels=c("c", "", "d", "", "e"), label_size=8,
               ncol=6, rel_widths=c(1+15+15, 2, 1+15+15, 2, 1+15+15, 2))

row2_sized=plot_grid(blank_plot, row2,
                     ncol=1, rel_heights=c(2, 37.5+1+.5),
                     labels=c("RefSeq gene/transcript annotations"), label_size=6, label_fontface="bold", label_x=.5, hjust=.5, label_y=.5, vjust=.5)

tvd=plot_grid(row1, blank_plot, row2_sized, ncol=1, rel_heights=c(2+37.5+1+.5,  .5,  2+37.5+1+.5))

# Width was 180 when total rel_widths were: 1+15+15+15 + 2 + 1+15+15 + 2 + 1+15+15 + 2 + 5
# Now the toal rel_widths are 1+15+15+15 + 2 + 1+15+15 + 8.5 + 5 + 6.5 / Next line: 1+15+15 + 2 + 1+15+15 + 2 + 1+15+15 + 2 (Note: both are same)
# So new width should be: 180 / (1+15+15+15 + 2 + 1+15+15 + 2 + 1+15+15 + 2 + 5) * (1+15+15 + 2 + 1+15+15 + 2 + 1+15+15 + 2) = 149.7479
# But should be 180 so multiply (180/149.7479) to it.

# Height was 50.65 when total rel_heights were: 2+37.5+1+.5
# Now the total rel_heights are 2+37.5+1+.5  +  .5  +  2+37.5+1+.5
# So new height should be 50.65 / (2+37.5+1+.5) * (2+37.5+1+.5  +  .5  +  2+37.5+1+.5) = 101.9177
# But WIDTH should be 180 so multiply the same factor (180/149.7479) to it i.e. 101.9177 * (180/149.7479) = 122.5071

ggsave("tvd.png", tvd, width=180, height=122.5071, units="mm")
ggsave("tvd.pdf", tvd, width=180, height=122.5071, units="mm")

save(tvd, file="tvd.RData")





####


# ★ line w/ Violin + Box
library(data.table)
library(dplyr)
library(stringr)
library(ggfortify)
library(tidyr)
library(cowplot)
library(ggplot2)
library(ggpubr)
library(scales)
library(grid)
library(png)
library(gridGraphics)
library(introdataviz)
library(ggnewscale)
library(distrEx)
setwd("D:/R/WD/MP_TSS/tvds")

blank = ggplot() + theme_void()
human_thin<-readPNG("human_thin.png")
zf_thin<-readPNG("zebra_finch_thin.png")

## Human gene

#

T2T_line_data = fread("MP_TSS_T2T-CHM13v2.0.tsv")
T2T_line_data[V4 == "-", V2 := V2 * -1] # If reverse (-) strand, flip the coordinates
T2T_line_data = data.frame(summarise(group_by(T2T_line_data, V2),"Mean"=mean(V3),"Median"=median(V3),"SD"=sd(V3)))
T2T_line_data=cbind(T2T_line_data, "T2T-CHM13v2.0")
colnames(T2T_line_data)=c("Distance","Mean","Median","SD","Assembly")

hg38_line_data = fread("MP_TSS_GRCh38.p14.tsv")
hg38_line_data[V4 == "-", V2 := V2 * -1] # If reverse (-) strand, flip the coordinates
hg38_line_data = data.frame(summarise(group_by(hg38_line_data, V2),"Mean"=mean(V3),"Median"=median(V3),"SD"=sd(V3)))
hg38_line_data=cbind(hg38_line_data, "GRCh38.p14")
colnames(hg38_line_data)=c("Distance","Mean","Median","SD","Assembly")

plot_line_mean_simple<-function(df, xmin, xmax, image, TEXT){
  ggplot(df, aes(x=Distance, y=Mean, group=Assembly)) +
    annotation_raster(image, xmin=-Inf, ymin=-Inf, xmax=Inf, ymax=Inf) +
    geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=-1000, xmax=-250), color="transparent", fill="#8FB4C6", alpha=.1) +
    geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=250, xmax=1000), color="transparent", fill="#CFE1EC", alpha=.1) +
    geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=-250, xmax=-150), color="transparent", fill="#7CCBA2", alpha=.1) +
    geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=150, xmax=250), color="transparent", fill="#B7E6A5", alpha=.1) +
    geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=-150, xmax=35), color="transparent", fill="#F07462", alpha=.1) +
    geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=35, xmax=150), color="transparent", fill="#FABF7B", alpha=.1) +
    geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=-35, xmax=0), color="transparent", fill="#AB1866", alpha=.1) +
    geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=0, xmax=35), color="transparent", fill="#D12959", alpha=.1) +
    geom_line(linewidth=.1, alpha=1) +
    #geom_vline(xintercept=-1000, color="#003147", linetype="dashed", linewidth=.1) +
    #geom_vline(xintercept=1000, color="#003147", linetype="dashed", linewidth=.1) +
    #geom_vline(xintercept=-250, color="#7CCBA2", linetype="dashed", linewidth=.1) +
    #geom_vline(xintercept=250, color="#7CCBA2", linetype="dashed", linewidth=.1) +
    #geom_vline(xintercept=-150, color="#F07462", linetype="dashed", linewidth=.1) +
    #geom_vline(xintercept=150, color="#F07462", linetype="dashed", linewidth=.1) +
    #geom_vline(xintercept=-35, color="#6E005F", linetype="dashed", linewidth=.1) +
    #geom_vline(xintercept=35, color="#6E005F", linetype="dashed", linewidth=.1) +
    #annotate("text", y=100, x=2500, vjust=1, hjust="inward", label = TEXT, color = "black", size = 1.94) + 
    scale_color_identity() +
    ylim(0,100) + scale_x_continuous(labels = label_comma(), breaks = c(xmin, 0, xmax), limits = c(xmin, xmax)) +
    labs(y="Mean MP (%)", x="Distance from GSS (bp)") + 
    theme(legend.position="none",
          panel.background = element_blank(),
          panel.grid = element_blank(),
          plot.background = element_blank(),
          axis.line.y.left = element_line(linewidth=0.1, color="black"),
          axis.line.x.bottom = element_line(linewidth=0.1, color="black"),
          axis.ticks = element_blank(),
          axis.text = element_text(size=5, color="black"),
          axis.title = element_text(size=5, color="black"))
}

#

get_box_data=function(filename){
  df=fread(filename)
  assembly=sub("MP_T[SE]S_(.*?)\\.tsv", "\\1", filename)
  df=cbind(df, assembly)
  colnames(df)=c("TSS", "Distance", "MP", "Strand", "Assembly")
  df$TrueDist = ifelse(df$Strand=="+", df$Distance,
                       ifelse(df$Strand=="-", df$Distance * -1, df$Distance))
  boxed_df=data.frame()
  # Condition
  boxed_df = rbind(boxed_df, cbind(df[-35<=df$TrueDist & df$TrueDist<0,], "5'", "35"))
  boxed_df = rbind(boxed_df, cbind(df[0<df$TrueDist & df$TrueDist<=35], "3'", "35"))
  boxed_df = rbind(boxed_df, cbind(df[-150<=df$TrueDist & df$TrueDist<0,], "5'", "150"))
  boxed_df = rbind(boxed_df, cbind(df[0<df$TrueDist & df$TrueDist<=150,], "3'", "150"))
  boxed_df = rbind(boxed_df, cbind(df[-250<=df$TrueDist & df$TrueDist<0,], "5'", "250"))
  boxed_df = rbind(boxed_df, cbind(df[0<df$TrueDist & df$TrueDist<=250], "3'", "250"))
  boxed_df = rbind(boxed_df, cbind(df[-1000<=df$TrueDist & df$TrueDist<0,], "5'", "1,000"))
  boxed_df = rbind(boxed_df, cbind(df[0<df$TrueDist & df$TrueDist<=1000,], "3'", "1,000"))
  colnames(boxed_df) = c("TSS", "Distance", "MP", "Strand", "Assembly", "TrueDist", "Condition", "Group")
  boxed_df=boxed_df[,c("Group", "Condition", "MP")]
  return(boxed_df)
}

T2T_box_data=get_box_data("MP_TSS_T2T-CHM13v2.0.tsv")
T2T_box_data_multicolorized = data.frame()
T2T_box_data_multicolorized = rbind(T2T_box_data_multicolorized, cbind(T2T_box_data[T2T_box_data$Group=="1,000" & T2T_box_data$Condition=="5'",], "#8FB4C6"))
T2T_box_data_multicolorized = rbind(T2T_box_data_multicolorized, cbind(T2T_box_data[T2T_box_data$Group=="1,000" & T2T_box_data$Condition=="3'",], "#CFE1EC"))
T2T_box_data_multicolorized = rbind(T2T_box_data_multicolorized, cbind(T2T_box_data[T2T_box_data$Group=="250" & T2T_box_data$Condition=="5'",], "#7CCBA2"))
T2T_box_data_multicolorized = rbind(T2T_box_data_multicolorized, cbind(T2T_box_data[T2T_box_data$Group=="250" & T2T_box_data$Condition=="3'",], "#B7E6A5"))
T2T_box_data_multicolorized = rbind(T2T_box_data_multicolorized, cbind(T2T_box_data[T2T_box_data$Group=="150" & T2T_box_data$Condition=="5'",], "#F07462"))
T2T_box_data_multicolorized = rbind(T2T_box_data_multicolorized, cbind(T2T_box_data[T2T_box_data$Group=="150" & T2T_box_data$Condition=="3'",], "#FABF7B"))
T2T_box_data_multicolorized = rbind(T2T_box_data_multicolorized, cbind(T2T_box_data[T2T_box_data$Group=="35" & T2T_box_data$Condition=="5'",], "#AB1866"))
T2T_box_data_multicolorized = rbind(T2T_box_data_multicolorized, cbind(T2T_box_data[T2T_box_data$Group=="35" & T2T_box_data$Condition=="3'",], "#D12959"))
colnames(T2T_box_data_multicolorized)=c("Group", "Condition", "MP", "Color")
temp=cbind("Genome / Gene", "5'", fread("CG_placedOnly_T2T-CHM13v2.0.mp")[,V4], "#AAAAAA"); colnames(temp)=c("Group", "Condition", "MP", "Color"); T2T_box_data_multicolorized = rbind(T2T_box_data_multicolorized, temp)
temp=cbind("Genome / Gene", "3'", fread("CG_geneOnly_T2T-CHM13v2.0.mp")[,V4], "#dddddd"); colnames(temp)=c("Group", "Condition", "MP", "Color"); T2T_box_data_multicolorized = rbind(T2T_box_data_multicolorized, temp)
T2T_box_data_multicolorized$MP = as.numeric(T2T_box_data_multicolorized$MP)
T2T_box_data_multicolorized$Group = factor(T2T_box_data_multicolorized$Group, levels=c("Genome / Gene", "35", "150", "250", "1,000"))
T2T_box_data_multicolorized$Condition = factor(T2T_box_data_multicolorized$Condition, levels=c("5'", "3'"))

hg38_box_data=get_box_data("MP_TSS_GRCh38.p14.tsv")
hg38_box_data_multicolorized = data.frame()
hg38_box_data_multicolorized = rbind(hg38_box_data_multicolorized, cbind(hg38_box_data[hg38_box_data$Group=="1,000" & hg38_box_data$Condition=="5'",], "#8FB4C6"))
hg38_box_data_multicolorized = rbind(hg38_box_data_multicolorized, cbind(hg38_box_data[hg38_box_data$Group=="1,000" & hg38_box_data$Condition=="3'",], "#CFE1EC"))
hg38_box_data_multicolorized = rbind(hg38_box_data_multicolorized, cbind(hg38_box_data[hg38_box_data$Group=="250" & hg38_box_data$Condition=="5'",], "#7CCBA2"))
hg38_box_data_multicolorized = rbind(hg38_box_data_multicolorized, cbind(hg38_box_data[hg38_box_data$Group=="250" & hg38_box_data$Condition=="3'",], "#B7E6A5"))
hg38_box_data_multicolorized = rbind(hg38_box_data_multicolorized, cbind(hg38_box_data[hg38_box_data$Group=="150" & hg38_box_data$Condition=="5'",], "#F07462"))
hg38_box_data_multicolorized = rbind(hg38_box_data_multicolorized, cbind(hg38_box_data[hg38_box_data$Group=="150" & hg38_box_data$Condition=="3'",], "#FABF7B"))
hg38_box_data_multicolorized = rbind(hg38_box_data_multicolorized, cbind(hg38_box_data[hg38_box_data$Group=="35" & hg38_box_data$Condition=="5'",], "#AB1866"))
hg38_box_data_multicolorized = rbind(hg38_box_data_multicolorized, cbind(hg38_box_data[hg38_box_data$Group=="35" & hg38_box_data$Condition=="3'",], "#D12959"))
colnames(hg38_box_data_multicolorized)=c("Group", "Condition", "MP", "Color")
temp=cbind("Genome / Gene", "5'", fread("CG_placedOnly_GRCh38.p14.mp")[,V4], "#AAAAAA"); colnames(temp)=c("Group", "Condition", "MP", "Color"); hg38_box_data_multicolorized = rbind(hg38_box_data_multicolorized, temp)
temp=cbind("Genome / Gene", "3'", fread("CG_geneOnly_GRCh38.p14.mp")[,V4], "#dddddd"); colnames(temp)=c("Group", "Condition", "MP", "Color"); hg38_box_data_multicolorized = rbind(hg38_box_data_multicolorized, temp)
hg38_box_data_multicolorized$MP = as.numeric(hg38_box_data_multicolorized$MP)
hg38_box_data_multicolorized$Group = factor(hg38_box_data_multicolorized$Group, levels=c("Genome / Gene", "35", "150", "250", "1,000"))
hg38_box_data_multicolorized$Condition = factor(hg38_box_data_multicolorized$Condition, levels=c("5'", "3'"))

plot_viol_simple<-function(df){
  ggplot(df) +
    geom_split_violin(mapping=aes(x=Group,y=MP,fill=Color), color="transparent", adjust=2, scale="width") +
    scale_fill_manual(values=c("#8FB4C6"="#8FB4C6", "#CFE1EC"="#CFE1EC", "#7CCBA2"="#7CCBA2", "#B7E6A5"="#B7E6A5", "#F07462"="#F07462", "#FABF7B"="#FABF7B", "#AB1866"="#AB1866", "#D12959"="#D12959", "#AAAAAA"="#AAAAAA", "#dddddd"="#dddddd"), na.value="red") +
    new_scale_fill() +
    geom_boxplot(mapping=aes(x=Group,y=MP,fill=Condition), color="#000000", width=0.25, position=position_dodge(), outlier.shape=NA, linewidth=0.1) + 
    scale_fill_manual(values=c("5'"="transparent"), na.value="transparent") +
    labs(y="MP (%)", x="Interval from GSS (bp)") +
    scale_y_continuous(limits=c(0,100), breaks=c(0,25,50,75,100)) +
    theme(legend.position="none", plot.background=element_blank(), panel.background=element_blank(), panel.grid=element_blank(), axis.title=element_text(size=5, color="black"), axis.ticks=element_blank(), axis.text=element_text(size=5, color="black"), axis.line=element_line(linewidth=0.1, color="black"))
}

#

bin_size=50; bound=2000

T2T_tvd_data=fread("MP_TSS_T2T-CHM13v2.0.tsv")[V2<=bound&V2>=(bound*-1),]
T2T_tvd_data$V2=ifelse(T2T_tvd_data$V4=="-", T2T_tvd_data$V2*-1, T2T_tvd_data$V2)
T2T_tvd_data$Dir = ifelse(T2T_tvd_data$V2<0, "5'", "3'")
T2T_tvd_data$Bin = abs(T2T_tvd_data$V2) %/% bin_size * bin_size # ← binning
T2T_tvd_data$Bin = abs(T2T_tvd_data$Bin)
colnames(T2T_tvd_data)=c("ID", "Dist", "MP", "Strand", "Dir", "Bin")
T2T_tvd_data=T2T_tvd_data[order(T2T_tvd_data$Bin),]
tvd_test <- function(df, value_col = "MP", dir_col = "Dir", bin_col = "Bin", n_perm = 1000) {
  calc_tvd <- function(x, y, breaks = NULL) {
    if (is.null(breaks)) {
      min_val <- min(c(x, y))
      max_val <- max(c(x, y))
      breaks <- seq(min_val, max_val, length.out = 30)
    }
    px <- hist(x, breaks = breaks, plot = FALSE)$counts
    py <- hist(y, breaks = breaks, plot = FALSE)$counts
    px <- px / sum(px)
    py <- py / sum(py)
    0.5 * sum(abs(px - py))
  }
  df[[dir_col]] <- as.character(df[[dir_col]])
  unique_bins <- unique(df[[bin_col]])
  results_list <- lapply(unique_bins, function(bin_val) {
    subdf <- subset(df, df[[bin_col]] == bin_val)
    subdf[[dir_col]] <- as.character(subdf[[dir_col]])
    
    group_vals <- unique(subdf[[dir_col]])
    if(length(group_vals) != 2) {
      cat(sprintf("Bin: %s - Skipped (does not have exactly two groups)\n", bin_val))
      return(data.frame(Bin = bin_val, TVD = NA, p_value = NA, Larger_Dir = NA))
    }
    
    x1 <- subdf[[value_col]][subdf[[dir_col]] == group_vals[1]]
    x2 <- subdf[[value_col]][subdf[[dir_col]] == group_vals[2]]
    
    if(length(x1) == 0 || length(x2) == 0) {
      cat(sprintf("Bin: %s - Skipped (empty group)\n", bin_val))
      return(data.frame(Bin = bin_val, TVD = NA, p_value = NA, Larger_Dir = NA))
    }
    
    obs_tvd <- calc_tvd(x1, x2)
    pooled <- c(x1, x2)
    n1 <- length(x1)
    
    perm_tvd <- replicate(n_perm, {
      perm_labels <- sample(rep(group_vals, c(n1, length(x2))))
      calc_tvd(pooled[perm_labels == group_vals[1]], pooled[perm_labels == group_vals[2]])
    })
    
    p_value <- mean(perm_tvd >= obs_tvd)
    
    mean1 <- mean(x1)
    mean2 <- mean(x2)
    
    larger_dir <- ifelse(mean1 > mean2, group_vals[1], group_vals[2])
    
    cat(sprintf("Bin: %s, TVD: %.4f, p-value: %.4f, Larger Dir: %s\n", bin_val, obs_tvd, p_value, larger_dir))
    
    data.frame(Bin = bin_val, TVD = obs_tvd, p_value = p_value, Larger_Dir = larger_dir)
  })
  do.call(rbind, results_list)
}
T2T_pvalues=tvd_test(T2T_tvd_data)

hg38_tvd_data=fread("MP_TSS_GRCh38.p14.tsv")[V2<=bound&V2>=(bound*-1),]
hg38_tvd_data$V2=ifelse(hg38_tvd_data$V4=="-", hg38_tvd_data$V2*-1, hg38_tvd_data$V2)
hg38_tvd_data$Dir = ifelse(hg38_tvd_data$V2<0, "5'", "3'")
hg38_tvd_data$Bin = abs(hg38_tvd_data$V2) %/% bin_size * bin_size # ← binning
hg38_tvd_data$Bin = abs(hg38_tvd_data$Bin)
colnames(hg38_tvd_data)=c("ID", "Dist", "MP", "Strand", "Dir", "Bin")
hg38_tvd_data=hg38_tvd_data[order(hg38_tvd_data$Bin),]
tvd_test <- function(df, value_col = "MP", dir_col = "Dir", bin_col = "Bin", n_perm = 1000) {
  calc_tvd <- function(x, y, breaks = NULL) {
    if (is.null(breaks)) {
      min_val <- min(c(x, y))
      max_val <- max(c(x, y))
      breaks <- seq(min_val, max_val, length.out = 30)
    }
    px <- hist(x, breaks = breaks, plot = FALSE)$counts
    py <- hist(y, breaks = breaks, plot = FALSE)$counts
    px <- px / sum(px)
    py <- py / sum(py)
    0.5 * sum(abs(px - py))
  }
  df[[dir_col]] <- as.character(df[[dir_col]])
  unique_bins <- unique(df[[bin_col]])
  results_list <- lapply(unique_bins, function(bin_val) {
    subdf <- subset(df, df[[bin_col]] == bin_val)
    subdf[[dir_col]] <- as.character(subdf[[dir_col]])
    
    group_vals <- unique(subdf[[dir_col]])
    if(length(group_vals) != 2) {
      cat(sprintf("Bin: %s - Skipped (does not have exactly two groups)\n", bin_val))
      return(data.frame(Bin = bin_val, TVD = NA, p_value = NA, Larger_Dir = NA))
    }
    
    x1 <- subdf[[value_col]][subdf[[dir_col]] == group_vals[1]]
    x2 <- subdf[[value_col]][subdf[[dir_col]] == group_vals[2]]
    
    if(length(x1) == 0 || length(x2) == 0) {
      cat(sprintf("Bin: %s - Skipped (empty group)\n", bin_val))
      return(data.frame(Bin = bin_val, TVD = NA, p_value = NA, Larger_Dir = NA))
    }
    
    obs_tvd <- calc_tvd(x1, x2)
    pooled <- c(x1, x2)
    n1 <- length(x1)
    
    perm_tvd <- replicate(n_perm, {
      perm_labels <- sample(rep(group_vals, c(n1, length(x2))))
      calc_tvd(pooled[perm_labels == group_vals[1]], pooled[perm_labels == group_vals[2]])
    })
    
    p_value <- mean(perm_tvd >= obs_tvd)
    
    mean1 <- mean(x1)
    mean2 <- mean(x2)
    
    larger_dir <- ifelse(mean1 > mean2, group_vals[1], group_vals[2])
    
    cat(sprintf("Bin: %s, TVD: %.4f, p-value: %.4f, Larger Dir: %s\n", bin_val, obs_tvd, p_value, larger_dir))
    
    data.frame(Bin = bin_val, TVD = obs_tvd, p_value = p_value, Larger_Dir = larger_dir)
  })
  do.call(rbind, results_list)
}
hg38_pvalues=tvd_test(hg38_tvd_data)

plot_pvalue_heatmap <- function(results_df, bin_col = "Bin", pvalue_col = "p_value", larger_dir_col = "Larger_Dir") {
  df <- results_df %>%
    mutate(
      fill_color = case_when(
        results_df$Larger_Dir == "5'" ~ 1 - .data[[pvalue_col]],   # invert p-value for red gradient
        results_df$Larger_Dir == "3'" ~ -(1 - .data[[pvalue_col]]), # invert p-value for blue gradient as negative
        TRUE ~ NA_real_
      )
    )
  ggplot(df, aes(x = Bin, y = 1, fill = fill_color)) +
    geom_tile(color="black") +
    scale_fill_gradientn(
      colors = c("blue", "white", "grey", "white", "red"),
      values = scales::rescale(c(-1, -0.95, 0, 0.95, 1)),
      limits = c(-1, 1),
      breaks = c(-1, -0.95, 0, 0.95, 1),
      labels = c("", "0.05 (3'>5')", "0", "0.05 5'>3'", ""),
      name = "p-value"
    )+
    theme_minimal() +
    theme(
      plot.background = element_blank(), panel.background = element_blank(), panel.grid = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.x=element_text(color="black", size=5),
      axis.title.x=element_text(color="black", size=5),
      legend.position="none"
    ) +
    labs(x = "Distance from GSS (bp)") +
    scale_x_continuous(labels=label_comma())
}

#

plot_VSline <- function(df){
  left=df[df$Dir=="5'"]; right=df[df$Dir=="3'"]
  left=data.frame(summarise(group_by(left, Bin),"Mean"=mean(MP))); right=data.frame(summarise(group_by(right, Bin),"Mean"=mean(MP)))
  both=data.frame(); both=rbind(both, cbind(left, "Dir"="5'")); both=rbind(both, cbind(right, "Dir"="3'"))
  ggplot(both)+
    geom_line(mapping=aes(x=Bin, y=Mean, group=Dir, color=Dir),linewidth=.25, alpha=.5)+
    scale_x_continuous(limits=c(0, bound), labels=label_comma())+
    scale_color_manual(name="Location respect to TSS", values=c("5'"="red","3'"="blue"),na.value="transparent") +
    scale_y_continuous(limits=c(trunc(min(both$Mean)), trunc(max(both$Mean))+1), breaks=c(trunc(min(both$Mean)), trunc(max(both$Mean))+1), name="\nMean MP (%)")+
    theme(plot.background=element_blank(), panel.background = element_blank(), panel.grid = element_blank(), axis.ticks=element_blank(),
          axis.title.y=element_text(size=5, color="black"), axis.text.y=element_text(size=5, color="black"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.line.y.left=element_line(linewidth=.1, color="black"), axis.line.x.bottom=element_blank(),
          legend.position="none"
    )
}

tile_legend=as_ggplot(get_legend(ggplot(data.frame("x"=c(-1,0,1),"y"=c(-1,0,1),"Location"=c("3’", "5’", "3’"))) + geom_line(mapping=aes(x=x,y=y,color=Location,group=Location)) +
  scale_color_manual(values=c("5’"="red", "3’"="blue"), name="Position respect to G(T)SS:", guide = guide_legend(order = 1)) +
  new_scale_fill() +
  geom_tile(mapping=aes(x=x,y=y,fill=x)) +
  scale_fill_gradientn(colors=c("blue","white","grey","white","red"), values = scales::rescale(c(-1, -0.95, 0, 0.95, 1)), breaks=c(-.95, 0, .95), labels=c("0.05 (3’ > 5’)", "1", "0.05 (5’ > 3’)"), name="P-value\n\n\n"
                       ) +
  theme(legend.ticks = element_blank(), legend.position="top",
    legend.title=element_text(size=5, color="black"),
    legend.text=element_text(size=5, color="black"),
    legend.key=element_rect(fill="transparent", color=NA),
    legend.key.width = unit(6, "mm"), legend.key.height=unit(2,"mm"),
    legend.background=element_rect(fill="transparent")
  ))); tile_legend

#

#save.image("D:/R/WD/MP_TSS/tvds/pvalues_human_gene.RData")
load("D:/R/WD/MP_TSS/tvds/pvalues_human_gene.RData")

T2T_line=plot_line_mean_simple(T2T_line_data,-2500,2500,human_thin,"T2T-CHM13v2.0")
T2T_viol=plot_viol_simple(T2T_box_data_multicolorized)
T2T_vsline=plot_VSline(T2T_tvd_data)
T2T_tile=plot_pvalue_heatmap(T2T_pvalues)

T2T_lineNviol = align_plots(T2T_line + theme(plot.margin = margin(t=0,b=5.5,l=5.5,r=5.5)), T2T_viol + theme(plot.margin = margin(t=0,b=0,l=5.5,r=5.5)), align="h", axis="tb")
T2T_vsNtile = align_plots(T2T_vsline + theme(plot.margin = margin(b=0)), T2T_tile + theme(plot.margin = margin(t=0)), align="v", axis="lr")

hg38_line=plot_line_mean_simple(hg38_line_data,-2500,2500,human_thin, "GRCh38.p14")
hg38_viol=plot_viol_simple(hg38_box_data_multicolorized)
hg38_vsline=plot_VSline(hg38_tvd_data)
hg38_tile=plot_pvalue_heatmap(hg38_pvalues)

hg38_lineNviol = align_plots(hg38_line + theme(plot.margin = margin(t=0,b=5.5,l=5.5,r=5.5)), hg38_viol + theme(plot.margin = margin(t=0,b=0,l=5.5,r=5.5)), align="h", axis="tb")
hg38_vsNtile = align_plots(hg38_vsline + theme(plot.margin = margin(b=0)), hg38_tile + theme(plot.margin = margin(t=0)), align="v", axis="lr")

row1=plot_grid(T2T_lineNviol[[1]], T2T_lineNviol[[2]], ncol=2, rel_widths=c(3.5,6.5))
human_gene_T2T=plot_grid(blank, row1, T2T_vsNtile[[1]], T2T_vsNtile[[2]], ncol=1, rel_heights=c(1,5,2.5,1.5),
                         labels=c("Gene start site (T2T-CHM13v2.0)"), label_fontface="plain", label_size=6, label_x=0.5, hjust=0.5)
row1=plot_grid(hg38_lineNviol[[1]], hg38_lineNviol[[2]], ncol=2, rel_widths=c(3.5,6.5))
human_gene_hg38=plot_grid(blank, row1, hg38_vsNtile[[1]], hg38_vsNtile[[2]], ncol=1, rel_heights=c(1,5,2.5,1.5),
                          labels=c("Gene start site (GRCh38.p14)"), label_fontface="plain", label_size=6, label_x=0.5, hjust=0.5)
human_gene=plot_grid(human_gene_T2T, blank, human_gene_hg38, ncol=3, rel_heights=c(1+5+2.5+1.5), rel_widths=c(20,1,20), # 0.5를 지웠음
                     labels=c("a", "", "b"), label_size=8)
#ggsave("human_gene.png", human_gene, width=180, height=55, units="mm")

tvd_human_gene = plot_grid(human_gene, ncol=1, rel_heights=c(55)) # 50을 지웠음
ggsave("tvd_human_gene.png", tvd_human_gene, height=105/105*55/10.5*10, width=180, units="mm")

save.image("D:/R/WD/MP_TSS/tvds/pvalues_human_gene.RData")
save(human_gene, file="human_gene.RData")

## Human transcript

#

T2T_line_data = fread("MP_transcriptTSS_T2T-CHM13v2.0.tsv")
T2T_line_data[V4 == "-", V2 := V2 * -1] # If reverse (-) strand, flip the coordinates
T2T_line_data = data.frame(summarise(group_by(T2T_line_data, V2),"Mean"=mean(V3),"Median"=median(V3),"SD"=sd(V3)))
T2T_line_data=cbind(T2T_line_data, "T2T-CHM13v2.0")
colnames(T2T_line_data)=c("Distance","Mean","Median","SD","Assembly")

hg38_line_data = fread("MP_transcriptTSS_GRCh38.p14.tsv")
hg38_line_data[V4 == "-", V2 := V2 * -1] # If reverse (-) strand, flip the coordinates
hg38_line_data = data.frame(summarise(group_by(hg38_line_data, V2),"Mean"=mean(V3),"Median"=median(V3),"SD"=sd(V3)))
hg38_line_data=cbind(hg38_line_data, "GRCh38.p14")
colnames(hg38_line_data)=c("Distance","Mean","Median","SD","Assembly")

plot_line_mean_simple<-function(df, xmin, xmax, image, TEXT){
  ggplot(df, aes(x=Distance, y=Mean, group=Assembly)) +
    annotation_raster(image, xmin=-Inf, ymin=-Inf, xmax=Inf, ymax=Inf) +
    geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=-1000, xmax=-250), color="transparent", fill="#8FB4C6", alpha=.1) +
    geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=250, xmax=1000), color="transparent", fill="#CFE1EC", alpha=.1) +
    geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=-250, xmax=-150), color="transparent", fill="#7CCBA2", alpha=.1) +
    geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=150, xmax=250), color="transparent", fill="#B7E6A5", alpha=.1) +
    geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=-150, xmax=35), color="transparent", fill="#F07462", alpha=.1) +
    geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=35, xmax=150), color="transparent", fill="#FABF7B", alpha=.1) +
    geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=-35, xmax=0), color="transparent", fill="#AB1866", alpha=.1) +
    geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=0, xmax=35), color="transparent", fill="#D12959", alpha=.1) +
    geom_line(linewidth=.1, alpha=1) +
    #geom_vline(xintercept=-1000, color="#003147", linetype="dashed", linewidth=.1) +
    #geom_vline(xintercept=1000, color="#003147", linetype="dashed", linewidth=.1) +
    #geom_vline(xintercept=-250, color="#7CCBA2", linetype="dashed", linewidth=.1) +
    #geom_vline(xintercept=250, color="#7CCBA2", linetype="dashed", linewidth=.1) +
    #geom_vline(xintercept=-150, color="#F07462", linetype="dashed", linewidth=.1) +
    #geom_vline(xintercept=150, color="#F07462", linetype="dashed", linewidth=.1) +
    #geom_vline(xintercept=-35, color="#6E005F", linetype="dashed", linewidth=.1) +
    #geom_vline(xintercept=35, color="#6E005F", linetype="dashed", linewidth=.1) +
    #annotate("text", y=100, x=2500, vjust=1, hjust="inward", label = TEXT, color = "black", size = 1.94) + 
    scale_color_identity() +
    ylim(0,100) + scale_x_continuous(labels = label_comma(), breaks = c(xmin, 0, xmax), limits = c(xmin, xmax)) +
    labs(y="Mean MP (%)", x="Distance from TSS (bp)") + 
    theme(legend.position="none",
          panel.background = element_blank(),
          panel.grid = element_blank(),
          plot.background = element_blank(),
          axis.line.y.left = element_line(linewidth=0.1, color="black"),
          axis.line.x.bottom = element_line(linewidth=0.1, color="black"),
          axis.ticks = element_blank(),
          axis.text = element_text(size=5, color="black"),
          axis.title = element_text(size=5, color="black"))
}

#

get_box_data=function(filename){
  df=fread(filename)
  assembly=sub("MP_transcriptT[SE]S_(.*?)\\.tsv", "\\1", filename)
  df=cbind(df, assembly)
  colnames(df)=c("TSS", "Distance", "MP", "Strand", "Assembly")
  df$TrueDist = ifelse(df$Strand=="+", df$Distance,
                       ifelse(df$Strand=="-", df$Distance * -1, df$Distance))
  boxed_df=data.frame()
  # Condition
  boxed_df = rbind(boxed_df, cbind(df[-35<=df$TrueDist & df$TrueDist<0,], "5'", "35"))
  boxed_df = rbind(boxed_df, cbind(df[0<df$TrueDist & df$TrueDist<=35], "3'", "35"))
  boxed_df = rbind(boxed_df, cbind(df[-150<=df$TrueDist & df$TrueDist<0,], "5'", "150"))
  boxed_df = rbind(boxed_df, cbind(df[0<df$TrueDist & df$TrueDist<=150,], "3'", "150"))
  boxed_df = rbind(boxed_df, cbind(df[-250<=df$TrueDist & df$TrueDist<0,], "5'", "250"))
  boxed_df = rbind(boxed_df, cbind(df[0<df$TrueDist & df$TrueDist<=250], "3'", "250"))
  boxed_df = rbind(boxed_df, cbind(df[-1000<=df$TrueDist & df$TrueDist<0,], "5'", "1,000"))
  boxed_df = rbind(boxed_df, cbind(df[0<df$TrueDist & df$TrueDist<=1000,], "3'", "1,000"))
  colnames(boxed_df) = c("TSS", "Distance", "MP", "Strand", "Assembly", "TrueDist", "Condition", "Group")
  boxed_df=boxed_df[,c("Group", "Condition", "MP")]
  return(boxed_df)
}

T2T_box_data=get_box_data("MP_transcriptTSS_T2T-CHM13v2.0.tsv")
T2T_box_data_multicolorized = data.frame()
T2T_box_data_multicolorized = rbind(T2T_box_data_multicolorized, cbind(T2T_box_data[T2T_box_data$Group=="1,000" & T2T_box_data$Condition=="5'",], "#8FB4C6"))
T2T_box_data_multicolorized = rbind(T2T_box_data_multicolorized, cbind(T2T_box_data[T2T_box_data$Group=="1,000" & T2T_box_data$Condition=="3'",], "#CFE1EC"))
T2T_box_data_multicolorized = rbind(T2T_box_data_multicolorized, cbind(T2T_box_data[T2T_box_data$Group=="250" & T2T_box_data$Condition=="5'",], "#7CCBA2"))
T2T_box_data_multicolorized = rbind(T2T_box_data_multicolorized, cbind(T2T_box_data[T2T_box_data$Group=="250" & T2T_box_data$Condition=="3'",], "#B7E6A5"))
T2T_box_data_multicolorized = rbind(T2T_box_data_multicolorized, cbind(T2T_box_data[T2T_box_data$Group=="150" & T2T_box_data$Condition=="5'",], "#F07462"))
T2T_box_data_multicolorized = rbind(T2T_box_data_multicolorized, cbind(T2T_box_data[T2T_box_data$Group=="150" & T2T_box_data$Condition=="3'",], "#FABF7B"))
T2T_box_data_multicolorized = rbind(T2T_box_data_multicolorized, cbind(T2T_box_data[T2T_box_data$Group=="35" & T2T_box_data$Condition=="5'",], "#AB1866"))
T2T_box_data_multicolorized = rbind(T2T_box_data_multicolorized, cbind(T2T_box_data[T2T_box_data$Group=="35" & T2T_box_data$Condition=="3'",], "#D12959"))
colnames(T2T_box_data_multicolorized)=c("Group", "Condition", "MP", "Color")
temp=cbind("Genome / Gene", "5'", fread("CG_placedOnly_T2T-CHM13v2.0.mp")[,V4], "#AAAAAA"); colnames(temp)=c("Group", "Condition", "MP", "Color"); T2T_box_data_multicolorized = rbind(T2T_box_data_multicolorized, temp)
temp=cbind("Genome / Gene", "3'", fread("CG_geneOnly_T2T-CHM13v2.0.mp")[,V4], "#dddddd"); colnames(temp)=c("Group", "Condition", "MP", "Color"); T2T_box_data_multicolorized = rbind(T2T_box_data_multicolorized, temp)
T2T_box_data_multicolorized$MP = as.numeric(T2T_box_data_multicolorized$MP)
T2T_box_data_multicolorized$Group = factor(T2T_box_data_multicolorized$Group, levels=c("Genome / Gene", "35", "150", "250", "1,000"))
T2T_box_data_multicolorized$Condition = factor(T2T_box_data_multicolorized$Condition, levels=c("5'", "3'"))

hg38_box_data=get_box_data("MP_transcriptTSS_GRCh38.p14.tsv")
hg38_box_data_multicolorized = data.frame()
hg38_box_data_multicolorized = rbind(hg38_box_data_multicolorized, cbind(hg38_box_data[hg38_box_data$Group=="1,000" & hg38_box_data$Condition=="5'",], "#8FB4C6"))
hg38_box_data_multicolorized = rbind(hg38_box_data_multicolorized, cbind(hg38_box_data[hg38_box_data$Group=="1,000" & hg38_box_data$Condition=="3'",], "#CFE1EC"))
hg38_box_data_multicolorized = rbind(hg38_box_data_multicolorized, cbind(hg38_box_data[hg38_box_data$Group=="250" & hg38_box_data$Condition=="5'",], "#7CCBA2"))
hg38_box_data_multicolorized = rbind(hg38_box_data_multicolorized, cbind(hg38_box_data[hg38_box_data$Group=="250" & hg38_box_data$Condition=="3'",], "#B7E6A5"))
hg38_box_data_multicolorized = rbind(hg38_box_data_multicolorized, cbind(hg38_box_data[hg38_box_data$Group=="150" & hg38_box_data$Condition=="5'",], "#F07462"))
hg38_box_data_multicolorized = rbind(hg38_box_data_multicolorized, cbind(hg38_box_data[hg38_box_data$Group=="150" & hg38_box_data$Condition=="3'",], "#FABF7B"))
hg38_box_data_multicolorized = rbind(hg38_box_data_multicolorized, cbind(hg38_box_data[hg38_box_data$Group=="35" & hg38_box_data$Condition=="5'",], "#AB1866"))
hg38_box_data_multicolorized = rbind(hg38_box_data_multicolorized, cbind(hg38_box_data[hg38_box_data$Group=="35" & hg38_box_data$Condition=="3'",], "#D12959"))
colnames(hg38_box_data_multicolorized)=c("Group", "Condition", "MP", "Color")
temp=cbind("Genome / Gene", "5'", fread("CG_placedOnly_GRCh38.p14.mp")[,V4], "#AAAAAA"); colnames(temp)=c("Group", "Condition", "MP", "Color"); hg38_box_data_multicolorized = rbind(hg38_box_data_multicolorized, temp)
temp=cbind("Genome / Gene", "3'", fread("CG_geneOnly_GRCh38.p14.mp")[,V4], "#dddddd"); colnames(temp)=c("Group", "Condition", "MP", "Color"); hg38_box_data_multicolorized = rbind(hg38_box_data_multicolorized, temp)
hg38_box_data_multicolorized$MP = as.numeric(hg38_box_data_multicolorized$MP)
hg38_box_data_multicolorized$Group = factor(hg38_box_data_multicolorized$Group, levels=c("Genome / Gene", "35", "150", "250", "1,000"))
hg38_box_data_multicolorized$Condition = factor(hg38_box_data_multicolorized$Condition, levels=c("5'", "3'"))

plot_viol_simple<-function(df){
  ggplot(df) +
    geom_split_violin(mapping=aes(x=Group,y=MP,fill=Color), color="transparent", adjust=2, scale="width") +
    scale_fill_manual(values=c("#8FB4C6"="#8FB4C6", "#CFE1EC"="#CFE1EC", "#7CCBA2"="#7CCBA2", "#B7E6A5"="#B7E6A5", "#F07462"="#F07462", "#FABF7B"="#FABF7B", "#AB1866"="#AB1866", "#D12959"="#D12959", "#AAAAAA"="#AAAAAA", "#dddddd"="#dddddd"), na.value="red") +
    new_scale_fill() +
    geom_boxplot(mapping=aes(x=Group,y=MP,fill=Condition), color="#000000", width=0.25, position=position_dodge(), outlier.shape=NA, linewidth=0.1) + 
    scale_fill_manual(values=c("5'"="transparent"), na.value="transparent") +
    labs(y="MP (%)", x="Interval from TSS (bp)") +
    scale_y_continuous(limits=c(0,100), breaks=c(0,25,50,75,100)) +
    theme(legend.position="none", plot.background=element_blank(), panel.background=element_blank(), panel.grid=element_blank(), axis.title=element_text(size=5, color="black"), axis.ticks=element_blank(), axis.text=element_text(size=5, color="black"), axis.line=element_line(linewidth=0.1, color="black"))
}

#

bin_size=50; bound=2000

T2T_tvd_data=fread("MP_transcriptTSS_T2T-CHM13v2.0.tsv")[V2<=bound&V2>=(bound*-1),]
T2T_tvd_data$V2=ifelse(T2T_tvd_data$V4=="-", T2T_tvd_data$V2*-1, T2T_tvd_data$V2)
T2T_tvd_data$Dir = ifelse(T2T_tvd_data$V2<0, "5'", "3'")
T2T_tvd_data$Bin = abs(T2T_tvd_data$V2) %/% bin_size * bin_size # ← binning
T2T_tvd_data$Bin = abs(T2T_tvd_data$Bin)
colnames(T2T_tvd_data)=c("ID", "Dist", "MP", "Strand", "Dir", "Bin")
T2T_tvd_data=T2T_tvd_data[order(T2T_tvd_data$Bin),]
tvd_test <- function(df, value_col = "MP", dir_col = "Dir", bin_col = "Bin", n_perm = 1000) {
  calc_tvd <- function(x, y, breaks = NULL) {
    if (is.null(breaks)) {
      min_val <- min(c(x, y))
      max_val <- max(c(x, y))
      breaks <- seq(min_val, max_val, length.out = 30)
    }
    px <- hist(x, breaks = breaks, plot = FALSE)$counts
    py <- hist(y, breaks = breaks, plot = FALSE)$counts
    px <- px / sum(px)
    py <- py / sum(py)
    0.5 * sum(abs(px - py))
  }
  df[[dir_col]] <- as.character(df[[dir_col]])
  unique_bins <- unique(df[[bin_col]])
  results_list <- lapply(unique_bins, function(bin_val) {
    subdf <- subset(df, df[[bin_col]] == bin_val)
    subdf[[dir_col]] <- as.character(subdf[[dir_col]])
    
    group_vals <- unique(subdf[[dir_col]])
    if(length(group_vals) != 2) {
      cat(sprintf("Bin: %s - Skipped (does not have exactly two groups)\n", bin_val))
      return(data.frame(Bin = bin_val, TVD = NA, p_value = NA, Larger_Dir = NA))
    }
    
    x1 <- subdf[[value_col]][subdf[[dir_col]] == group_vals[1]]
    x2 <- subdf[[value_col]][subdf[[dir_col]] == group_vals[2]]
    
    if(length(x1) == 0 || length(x2) == 0) {
      cat(sprintf("Bin: %s - Skipped (empty group)\n", bin_val))
      return(data.frame(Bin = bin_val, TVD = NA, p_value = NA, Larger_Dir = NA))
    }
    
    obs_tvd <- calc_tvd(x1, x2)
    pooled <- c(x1, x2)
    n1 <- length(x1)
    
    perm_tvd <- replicate(n_perm, {
      perm_labels <- sample(rep(group_vals, c(n1, length(x2))))
      calc_tvd(pooled[perm_labels == group_vals[1]], pooled[perm_labels == group_vals[2]])
    })
    
    p_value <- mean(perm_tvd >= obs_tvd)
    
    mean1 <- mean(x1)
    mean2 <- mean(x2)
    
    larger_dir <- ifelse(mean1 > mean2, group_vals[1], group_vals[2])
    
    cat(sprintf("Bin: %s, TVD: %.4f, p-value: %.4f, Larger Dir: %s\n", bin_val, obs_tvd, p_value, larger_dir))
    
    data.frame(Bin = bin_val, TVD = obs_tvd, p_value = p_value, Larger_Dir = larger_dir)
  })
  do.call(rbind, results_list)
}
T2T_pvalues=tvd_test(T2T_tvd_data)

hg38_tvd_data=fread("MP_transcriptTSS_GRCh38.p14.tsv")[V2<=bound&V2>=(bound*-1),]
hg38_tvd_data$V2=ifelse(hg38_tvd_data$V4=="-", hg38_tvd_data$V2*-1, hg38_tvd_data$V2)
hg38_tvd_data$Dir = ifelse(hg38_tvd_data$V2<0, "5'", "3'")
hg38_tvd_data$Bin = abs(hg38_tvd_data$V2) %/% bin_size * bin_size # ← binning
hg38_tvd_data$Bin = abs(hg38_tvd_data$Bin)
colnames(hg38_tvd_data)=c("ID", "Dist", "MP", "Strand", "Dir", "Bin")
hg38_tvd_data=hg38_tvd_data[order(hg38_tvd_data$Bin),]
tvd_test <- function(df, value_col = "MP", dir_col = "Dir", bin_col = "Bin", n_perm = 1000) {
  calc_tvd <- function(x, y, breaks = NULL) {
    if (is.null(breaks)) {
      min_val <- min(c(x, y))
      max_val <- max(c(x, y))
      breaks <- seq(min_val, max_val, length.out = 30)
    }
    px <- hist(x, breaks = breaks, plot = FALSE)$counts
    py <- hist(y, breaks = breaks, plot = FALSE)$counts
    px <- px / sum(px)
    py <- py / sum(py)
    0.5 * sum(abs(px - py))
  }
  df[[dir_col]] <- as.character(df[[dir_col]])
  unique_bins <- unique(df[[bin_col]])
  results_list <- lapply(unique_bins, function(bin_val) {
    subdf <- subset(df, df[[bin_col]] == bin_val)
    subdf[[dir_col]] <- as.character(subdf[[dir_col]])
    
    group_vals <- unique(subdf[[dir_col]])
    if(length(group_vals) != 2) {
      cat(sprintf("Bin: %s - Skipped (does not have exactly two groups)\n", bin_val))
      return(data.frame(Bin = bin_val, TVD = NA, p_value = NA, Larger_Dir = NA))
    }
    
    x1 <- subdf[[value_col]][subdf[[dir_col]] == group_vals[1]]
    x2 <- subdf[[value_col]][subdf[[dir_col]] == group_vals[2]]
    
    if(length(x1) == 0 || length(x2) == 0) {
      cat(sprintf("Bin: %s - Skipped (empty group)\n", bin_val))
      return(data.frame(Bin = bin_val, TVD = NA, p_value = NA, Larger_Dir = NA))
    }
    
    obs_tvd <- calc_tvd(x1, x2)
    pooled <- c(x1, x2)
    n1 <- length(x1)
    
    perm_tvd <- replicate(n_perm, {
      perm_labels <- sample(rep(group_vals, c(n1, length(x2))))
      calc_tvd(pooled[perm_labels == group_vals[1]], pooled[perm_labels == group_vals[2]])
    })
    
    p_value <- mean(perm_tvd >= obs_tvd)
    
    mean1 <- mean(x1)
    mean2 <- mean(x2)
    
    larger_dir <- ifelse(mean1 > mean2, group_vals[1], group_vals[2])
    
    cat(sprintf("Bin: %s, TVD: %.4f, p-value: %.4f, Larger Dir: %s\n", bin_val, obs_tvd, p_value, larger_dir))
    
    data.frame(Bin = bin_val, TVD = obs_tvd, p_value = p_value, Larger_Dir = larger_dir)
  })
  do.call(rbind, results_list)
}
hg38_pvalues=tvd_test(hg38_tvd_data)

plot_pvalue_heatmap <- function(results_df, bin_col = "Bin", pvalue_col = "p_value", larger_dir_col = "Larger_Dir") {
  df <- results_df %>%
    mutate(
      fill_color = case_when(
        results_df$Larger_Dir == "5'" ~ 1 - .data[[pvalue_col]],   # invert p-value for red gradient
        results_df$Larger_Dir == "3'" ~ -(1 - .data[[pvalue_col]]), # invert p-value for blue gradient as negative
        TRUE ~ NA_real_
      )
    )
  ggplot(df, aes(x = Bin, y = 1, fill = fill_color)) +
    geom_tile(color="black") +
    scale_fill_gradientn(
      colors = c("blue", "white", "grey", "white", "red"),
      values = scales::rescale(c(-1, -0.95, 0, 0.95, 1)),
      limits = c(-1, 1),
      breaks = c(-1, -0.95, 0, 0.95, 1),
      labels = c("", "0.05 (3'>5')", "0", "0.05 5'>3'", ""),
      name = "p-value"
    )+
    theme_minimal() +
    theme(
      plot.background = element_blank(), panel.background = element_blank(), panel.grid = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.x=element_text(color="black", size=5),
      axis.title.x=element_text(color="black", size=5),
      legend.position="none"
    ) +
    labs(x = "Distance from TSS (bp)") +
    scale_x_continuous(labels=label_comma())
}

#

plot_VSline <- function(df){
  left=df[df$Dir=="5'"]; right=df[df$Dir=="3'"]
  left=data.frame(summarise(group_by(left, Bin),"Mean"=mean(MP))); right=data.frame(summarise(group_by(right, Bin),"Mean"=mean(MP)))
  both=data.frame(); both=rbind(both, cbind(left, "Dir"="5'")); both=rbind(both, cbind(right, "Dir"="3'"))
  ggplot(both)+
    geom_line(mapping=aes(x=Bin, y=Mean, group=Dir, color=Dir),linewidth=.25, alpha=.5)+
    scale_x_continuous(limits=c(0, bound), labels=label_comma())+
    scale_color_manual(name="Location respect to TSS", values=c("5'"="red","3'"="blue"),na.value="transparent") +
    scale_y_continuous(limits=c(trunc(min(both$Mean)), trunc(max(both$Mean))+1), breaks=c(trunc(min(both$Mean)), trunc(max(both$Mean))+1), name="\nMean MP (%)")+
    theme(plot.background=element_blank(), panel.background = element_blank(), panel.grid = element_blank(), axis.ticks=element_blank(),
          axis.title.y=element_text(size=5, color="black"), axis.text.y=element_text(size=5, color="black"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.line.y.left=element_line(linewidth=.1, color="black"), axis.line.x.bottom=element_blank(),
          legend.position="none"
    )
}

tile_legend=as_ggplot(get_legend(ggplot(data.frame("x"=c(-1,0,1),"y"=c(-1,0,1),"Location"=c("3’", "5’", "3’"))) + geom_line(mapping=aes(x=x,y=y,color=Location,group=Location)) +
                                   scale_color_manual(values=c("5’"="red", "3’"="blue"), name="Position respect to G(T)SS:", guide = guide_legend(order = 1)) +
                                   new_scale_fill() +
                                   geom_tile(mapping=aes(x=x,y=y,fill=x)) +
                                   scale_fill_gradientn(colors=c("blue","white","grey","white","red"), values = scales::rescale(c(-1, -0.95, 0, 0.95, 1)), breaks=c(-.95, 0, .95), labels=c("0.05 (3’ > 5’)", "1", "0.05 (5’ > 3’)"), name="P-value\n\n\n"
                                   ) +
                                   theme(legend.ticks = element_blank(), legend.position="top",
                                         legend.title=element_text(size=5, color="black"),
                                         legend.text=element_text(size=5, color="black"),
                                         legend.key=element_rect(fill="transparent", color=NA),
                                         legend.key.width = unit(6, "mm"), legend.key.height=unit(2,"mm"),
                                         legend.background=element_rect(fill="transparent")
                                   ))); tile_legend

#
load("pvalues_human_transcript.RData")
load("human_gene.RData")

T2T_line=plot_line_mean_simple(T2T_line_data,-2500,2500,human_thin,"T2T-CHM13v2.0")
T2T_viol=plot_viol_simple(T2T_box_data_multicolorized)
T2T_vsline=plot_VSline(T2T_tvd_data)
T2T_tile=plot_pvalue_heatmap(T2T_pvalues)

T2T_lineNviol = align_plots(T2T_line + theme(plot.margin = margin(t=0,b=5.5,l=5.5,r=5.5)), T2T_viol + theme(plot.margin = margin(t=0,b=0,l=5.5,r=5.5)), align="h", axis="tb")
T2T_vsNtile = align_plots(T2T_vsline + theme(plot.margin = margin(b=0)), T2T_tile + theme(plot.margin = margin(t=0)), align="v", axis="lr")

hg38_line=plot_line_mean_simple(hg38_line_data,-2500,2500,human_thin, "GRCh38.p14")
hg38_viol=plot_viol_simple(hg38_box_data_multicolorized)
hg38_vsline=plot_VSline(hg38_tvd_data)
hg38_tile=plot_pvalue_heatmap(hg38_pvalues)

hg38_lineNviol = align_plots(hg38_line + theme(plot.margin = margin(t=0,b=5.5,l=5.5,r=5.5)), hg38_viol + theme(plot.margin = margin(t=0,b=0,l=5.5,r=5.5)), align="h", axis="tb")
hg38_vsNtile = align_plots(hg38_vsline + theme(plot.margin = margin(b=0)), hg38_tile + theme(plot.margin = margin(t=0)), align="v", axis="lr")

row1=plot_grid(T2T_lineNviol[[1]], T2T_lineNviol[[2]], ncol=2, rel_widths=c(3.5,6.5))
human_transcript_T2T=plot_grid(blank, row1, T2T_vsNtile[[1]], T2T_vsNtile[[2]], ncol=1, rel_heights=c(1,5,2.5,1.5),
                         labels=c("Transcript start site (T2T-CHM13v2.0)"), label_fontface="plain", label_size=6, label_x=0.5, hjust=0.5)
row1=plot_grid(hg38_lineNviol[[1]], hg38_lineNviol[[2]], ncol=2, rel_widths=c(3.5,6.5))
human_transcript_hg38=plot_grid(blank, row1, hg38_vsNtile[[1]], hg38_vsNtile[[2]], ncol=1, rel_heights=c(1,5,2.5,1.5),
                          labels=c("Transcript start site (GRCh38.p14)"), label_fontface="plain", label_size=6, label_x=0.5, hjust=0.5)
human_transcript=plot_grid(blank, blank, blank, human_transcript_T2T, blank, human_transcript_hg38, ncol=3, rel_heights=c(0.5, 1+5+2.5+1.5), rel_widths=c(20,1,20),
                     labels=c("", "", "", "c", "", "d"), label_size=8)
#ggsave("human_gene.png", human_gene, width=180, height=55, units="mm")

tvd_human_transcript = plot_grid(human_gene, human_transcript, tile_legend, ncol=1, rel_heights=c(52.38095,55,15)) # 50을 지우고 55→52.38095로 고치고 10→15로 고침.
ggsave("tvd_human_transcript.png", tvd_human_transcript, height=170 / (50+55+55+10) * (52.38095+55+15), width=180, units="mm")
ggsave("tvd_human_transcript.svg", tvd_human_transcript, height=170 / (50+55+55+10) * (52.38095+55+15), width=180, units="mm")
ggsave("tvd_human_transcript.pdf", tvd_human_transcript, height=170 / (50+55+55+10) * (52.38095+55+15), width=180, units="mm")

save.image("D:/R/WD/MP_TSS/tvds/pvalues_human_transcript.RData")
save(human_transcript, file="human_transcript.RData")










### Extra for VGP meeting only ###
load("D:/R/WD/MP_TSS/tvds/pvalues_rightbeforecombine.RData")

load("pvalues0.RData")

plot_pvalue_heatmap <- function(p_values_df, title) {
  ggplot_heatmap_data <- p_values_df; ggplot_heatmap_data$one <- "1"
  ggplot(ggplot_heatmap_data, aes(spcdist, one, fill = p_value, height = 1)) +
    geom_tile() +
    scale_fill_gradientn(colours=c("#FF4B33", "white", "#BBBBBB"), values=rescale(c(0, 0.05, 1)), limits=c(0, 1), breaks=c(0, 0.05, 1), labels=c("0", "5", "100"), name="p-value (%)") + labs(title=title) +
    theme(plot.background=element_blank(), panel.background=element_blank(), panel.grid=element_blank(), axis.ticks=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), plot.title = element_text(size = 12, vjust=-2), legend.position="none")
  pvalue_heatmap <- ggplot(ggplot_heatmap_data, aes(spcdist, one, fill = p_value, height = 1)) +
    geom_tile() +
    scale_fill_gradientn(colours=c("#FF4B33", "white", "#BBBBBB"), values=rescale(c(0, 0.05, 1)), limits=c(0, 1), breaks=c(0, 0.05, 1), labels=c("0", "5", "100"), name="p-value (%)") + labs(title=title) +
    theme(plot.background=element_blank(), panel.background=element_blank(), panel.grid=element_blank(), axis.ticks=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), plot.title = element_text(size = 12, vjust=-2), legend.position="none")
  return(pvalue_heatmap)
}

plot_legendOnly_heatmap <- function(p_values_df){
  ggplot_heatmap_data <- p_values_df; ggplot_heatmap_data$one <- "1"
  ggplot(data.frame(spcdist=1,p_value=1, one="1"), aes(spcdist, one, fill=p_value, height=1)) +
    geom_tile() + 
    scale_fill_gradientn(colours = c("#FF4B33", "white", "#BBBBBB"), values = scales::rescale(c(0, 0.05, 1)), limits = c(0, 1), breaks = c(1, 0.05), labels = c("1", "0.05"), name = "P-value") +
    guides(fill = guide_colourbar(label=TRUE)) +
    theme(#legend.position = "top",
      #legend.ticks = element_blank(),
      legend.background = element_blank(),
      #legend.title.position = "left",
      #legend.text.position = "bottom",
      legend.title=element_text(size=10),
      legend.text=element_text(size=8),
      legend.key.width = unit(2, "mm"), legend.key.height=unit(2.5, "mm")
    )}
# Compute only legend
heatmap_fakeForLegend = plot_legendOnly_heatmap(data.frame(spcdist=1, p_value=1))
heatmap_legend = as_ggplot(get_legend(heatmap_fakeForLegend)); heatmap_legend

data_preproc_median = function(df){
  df[V4 == "-", V2 := V2 * -1] # If reverse (-) strand, flip the coordinates
  df = data.frame(summarise(group_by(df, V2),med=median(V3),SD=sd(V3)))
  colnames(df)=c("Distance","Median","SD")
  #df = summarise(group_by(mutate(df, x_bin = round(Distance / 100) * 100),x_bin),Average = mean(Average,na.rm=TRUE))
  # ↑ Downsamples the data with BIN size of 100 bp. Change the value at: “x_bin = round(Distance / n) * n” to change the BIN size.
  colnames(df)=c("Distance","Median")
  return(df)
}

data_preproc_mean = function(df){
  df[V4 == "-", V2 := V2 * -1] # If reverse (-) strand, flip the coordinates
  df = data.frame(summarise(group_by(df, V2),avg=mean(V3),SD=sd(V3)))
  colnames(df)=c("Distance","Mean","SD")
  #df = summarise(group_by(mutate(df, x_bin = round(Distance / 100) * 100),x_bin),Average = mean(Average,na.rm=TRUE))
  # ↑ Downsamples the data with BIN size of 100 bp. Change the value at: “x_bin = round(Distance / n) * n” to change the BIN size.
  colnames(df)=c("Distance","Mean")
  return(df)
}

heatmap_TSS=plot_pvalue_heatmap(pvalues_TSS, "Start")
heatmap_TTS=plot_pvalue_heatmap(pvalues_TES, "End")
heatmap_promoter_ENCODE=plot_pvalue_heatmap(pvalues_promoter_ENCODE, "Promoter")
heatmap_promoter_NCBI=plot_pvalue_heatmap(pvalues_promoter_NCBI, "Promoter")
heatmap_enhancer_ENCODE=plot_pvalue_heatmap(pvalues_enhancer_ENCODE, "Enhancer")
heatmap_enhancer_NCBI=plot_pvalue_heatmap(pvalues_enhancer_NCBI, "Enhancer")
heatmap_silencer_NCBI=plot_pvalue_heatmap(pvalues_silencer_NCBI, "Silencer")

df = data_preproc_median(fread("MP_TSS_hg38.tsv"))
ggplot(df, aes(x = Distance, y = Median)) + 
  geom_line(linewidth=.1, alpha=1) +
  #geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,50,100)) +
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Median MP (%)", x = "Distance from TSS (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = label_comma(), limits=c(-10000,10000), breaks=c(-10000,0,10000)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 10, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_TSS_median

df = data_preproc_mean(fread("MP_TSS_hg38.tsv"))
ggplot(df, aes(x = Distance, y = Mean)) + 
  geom_line(linewidth=.1, alpha=1) +
  #geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  geom_rect(aes(ymax=100,ymin=0,xmin=-2500,xmax=2500),fill="transparent", color="red", linewidth=.1, linetype="dashed") +
  scale_y_continuous(limits=c(0,100), breaks=c(0,50,100)) + 
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Average MP (%)", x = "Distance from TSS (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = label_comma(), limits=c(-10000,10000), breaks=c(-10000,0,10000)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 10, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_TSS_mean

# TTS

df = data_preproc_median(fread("MP_TES_hg38.tsv"))
ggplot(df, aes(x = Distance, y = Median)) + 
  geom_line(linewidth=.1, alpha=1) +
  #geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,50,100)) +
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Median MP (%)", x = "Distance from TTS (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = label_comma(), limits=c(-10000,10000), breaks=c(-10000,0,10000)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 10, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_TTS_median

df = data_preproc_mean(fread("MP_TES_hg38.tsv"))
ggplot(df, aes(x = Distance, y = Mean)) + 
  geom_line(linewidth=.1, alpha=1) +
  #geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,50,100)) +
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Average MP (%)", x = "Distance from TTS (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = label_comma(), limits=c(-10000,10000), breaks=c(-10000,0,10000)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 10, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_TTS_mean

# Promoter NCBI

df = data_preproc_median(fread("MP_promoter_hg38.tsv"))
ggplot(df, aes(x = Distance, y = Median)) + 
  geom_line(linewidth=.1, alpha=1) +
  geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,50,100)) +
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Median MP (%)", x = "Distance from center (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = label_comma(), limits=c(-10000,10000), breaks=c(-10000,0,10000)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 10, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_promoterNCBI_median

df = data_preproc_mean(fread("MP_promoter_hg38.tsv"))
ggplot(df, aes(x = Distance, y = Mean)) + 
  geom_line(linewidth=.1, alpha=1) +
  geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,50,100)) +
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Average MP (%)", x = "Distance from center (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = label_comma(), limits=c(-10000,10000), breaks=c(-10000,0,10000)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 10, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_promoterNCBI_mean

# Promoter ENCODE

df = data_preproc_median(fread("MP_promoter_ENCODE.tsv"))
ggplot(df, aes(x = Distance, y = Median)) + 
  geom_line(linewidth=.1, alpha=1) +
  #geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,50,100)) + 
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Median MP (%)", x = "Distance from center (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = label_comma(), limits=c(-10000,10000), breaks=c(-10000,0,10000)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 10, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_promoterENCODE_median

df = data_preproc_mean(fread("MP_promoter_ENCODE.tsv"))
ggplot(df, aes(x = Distance, y = Mean)) + 
  geom_line(linewidth=.1, alpha=1) +
  #geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,50,100)) +
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Average MP (%)", x = "Distance from center (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = label_comma(), limits=c(-10000,10000), breaks=c(-10000,0,10000)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 10, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_promoterENCODE_mean

# Enhancer NCBI

df = data_preproc_median(fread("MP_enhancer_hg38.tsv"))
ggplot(df, aes(x = Distance, y = Median)) + 
  geom_line(linewidth=.1, alpha=1) +
  #geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,50,100)) +
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Median MP (%)", x = "Distance from center (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = label_comma(), limits=c(-10000,10000), breaks=c(-10000,0,10000)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 10, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_enhancerNCBI_median

df = data_preproc_mean(fread("MP_enhancer_hg38.tsv"))
ggplot(df, aes(x = Distance, y = Mean)) + 
  geom_line(linewidth=.1, alpha=1) +
  #geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,50,100)) +
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Average MP (%)", x = "Distance from center (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = label_comma(), limits=c(-10000,10000), breaks=c(-10000,0,10000)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 10, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_enhancerNCBI_mean


# Enhancer ENCODE

df = data_preproc_median(fread("MP_enhancer_ENCODE.tsv"))
ggplot(df, aes(x = Distance, y = Median)) + 
  geom_line(linewidth=.1, alpha=1) +
  #geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,50,100)) +
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Median MP (%)", x = "Distance from center (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = label_comma(), limits=c(-10000,10000), breaks=c(-10000,0,10000)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 10, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_enhancerENCODE_median

df = data_preproc_mean(fread("MP_enhancer_ENCODE.tsv"))
ggplot(df, aes(x = Distance, y = Mean)) + 
  geom_line(linewidth=.1, alpha=1) +
  #geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,50,100)) +
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Average MP (%)", x = "Distance from center (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = label_comma(), limits=c(-10000,10000), breaks=c(-10000,0,10000)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 10, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_enhancerENCODE_mean


# Silencer NCBI

df = data_preproc_median(fread("MP_silencer_hg38.tsv"))
ggplot(df, aes(x = Distance, y = Median)) + 
  geom_line(linewidth=.1, alpha=1) +
  #geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,50,100)) +
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Median MP (%)", x = "Distance from center (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = label_comma(), limits=c(-10000,10000), breaks=c(-10000,0,10000)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 10, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_silencerNCBI_median

df = data_preproc_mean(fread("MP_silencer_hg38.tsv"))
ggplot(df, aes(x = Distance, y = Mean)) + 
  geom_line(linewidth=.1, alpha=1) +
  #geom_smooth(method="loess", color="#DDDDDD", linewidth=.2) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,50,100)) +
  # coord_fixed(ratio = 200) + 
  # ↑ x-axis length ÷ y-axis length (e.g. 20,000 ÷ 100 = 200)
  labs(y = "Average MP (%)", x = "Distance from center (bp)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(labels = label_comma(), limits=c(-10000,10000), breaks=c(-10000,0,10000)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(linewidth=.1),
    axis.line.x.bottom = element_line(linewidth=.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    axis.text = element_text(size = 10, color="black"),
    legend.position = "none" # Comment out this line if you DO want to include a legend
  )->line_silencerNCBI_mean

blank_plot <- ggplot() + theme_void()

## Now let's align & combine all those plots into a multi-panel figure ##
# 1. Align plots:

aligned_TSS <- align_plots(
  heatmap_TSS + theme(plot.margin = margin(0,0,0,10)),
  line_TSS_median + theme(plot.margin = margin(0,0,5,5)),
  line_TSS_mean + theme(plot.margin = margin(0,0,5,5)),
  align = "v",
  axis = "lr"
)

aligned_TTS <- align_plots(
  heatmap_TTS + theme(plot.margin = margin(0,0,0,10)),
  line_TTS_median + theme(plot.margin = margin(0,0,5,5)),
  line_TTS_mean + theme(plot.margin = margin(0,0,5,5)),
  align = "v",
  axis = "lr"
)

aligned_promoterNCBI <- align_plots(
  heatmap_promoter_NCBI + theme(plot.margin = margin(0,0,0,10)),
  line_promoterNCBI_median + theme(plot.margin = margin(0,0,5,5)),
  line_promoterNCBI_mean + theme(plot.margin = margin(0,0,5,5)),
  align = "v",
  axis = "lr"
)

aligned_promoterENCODE <- align_plots(
  heatmap_promoter_ENCODE + theme(plot.margin = margin(0,0,0,10)),
  line_promoterENCODE_median + theme(plot.margin = margin(0,0,5,5)),
  line_promoterENCODE_mean + theme(plot.margin = margin(0,0,5,5)),
  align = "v",
  axis = "lr"
)

aligned_enhancerNCBI <- align_plots(
  heatmap_enhancer_NCBI + theme(plot.margin = margin(0,0,0,10)),
  line_enhancerNCBI_median + theme(plot.margin = margin(0,0,5,5)),
  line_enhancerNCBI_mean + theme(plot.margin = margin(0,0,5,5)),
  align = "v",
  axis = "lr"
)

aligned_enhancerENCODE <- align_plots(
  heatmap_enhancer_ENCODE + theme(plot.margin = margin(0,0,0,10)),
  line_enhancerENCODE_median + theme(plot.margin = margin(0,0,5,5)),
  line_enhancerENCODE_mean + theme(plot.margin = margin(0,0,5,5)),
  align = "v",
  axis = "lr"
)

aligned_silencerNCBI <- align_plots(
  heatmap_silencer_NCBI + theme(plot.margin = margin(0,0,0,10)),
  line_silencerNCBI_median + theme(plot.margin = margin(0,0,5,5)),
  line_silencerNCBI_mean + theme(plot.margin = margin(0,0,5,5)),
  align = "v",
  axis = "lr"
)

# 2. Combine plots:
median_label = ggdraw() + draw_grob(textGrob("          Median MP (%)", rot = 90, y=.5, x=1, gp = gpar(fontsize = 10, fontface = "plain")))
mean_label = ggdraw() + draw_grob(textGrob("          Mean MP (%)", rot = 90, y=.5, x=1, gp = gpar(fontsize = 10, fontface = "plain")))

NCBI=plot_grid(blank_plot, aligned_promoterNCBI[[1]], aligned_enhancerNCBI[[1]], aligned_silencerNCBI[[1]],
               mean_label, aligned_promoterNCBI[[3]], aligned_enhancerNCBI[[3]], aligned_silencerNCBI[[3]],
               median_label, aligned_promoterNCBI[[2]], aligned_enhancerNCBI[[2]], aligned_silencerNCBI[[2]],
               ncol=4, rel_heights=c(6,15,15), rel_widths=c(1,15,15,15))
NCBI_halfTitled=plot_grid(blank_plot, NCBI,
                          labels=c("RefSeq regulators"), label_size=16, label_fontface="bold", label_x=0.5, hjust=0.5, label_y=0.5, vjust=0.5,
                          ncol=1, rel_heights=c(4, 36))
NCBI_titled=plot_grid(NCBI_halfTitled, blank_plot,
                      labels=c("", "Distance from center (bp)"), label_size=10, label_fontface="plain", label_x=0.5, hjust=0.5, label_y=0.5, vjust=0.5,
                      ncol=1, rel_heights=c(4+36 ,2))

ENCODE=plot_grid(blank_plot, aligned_promoterENCODE[[1]], aligned_enhancerENCODE[[1]],
                 mean_label, aligned_promoterENCODE[[3]], aligned_enhancerENCODE[[3]],
                 median_label, aligned_promoterENCODE[[2]], aligned_enhancerENCODE[[2]],
                 ncol=3, rel_heights=c(6,15,15), rel_widths=c(1,15,15))
ENCODE_halfTitled=plot_grid(blank_plot, ENCODE,
                            labels=c("ENCODE regulators"), label_size=16, label_fontface="bold", label_x=0.5, hjust=0.5, label_y=0.5, vjust=0.5,
                            ncol=1, rel_heights=c(4, 36))
ENCODE_titled=plot_grid(ENCODE_halfTitled, blank_plot,
                        labels=c("", "Distance from center (bp)"), label_size=10, label_fontface="plain", label_x=0.5, hjust=0.5, label_y=0.5, vjust=0.5,
                        ncol=1, rel_heights=c(4+36 ,2))

RefSeq=plot_grid(blank_plot, aligned_TSS[[1]], aligned_TTS[[1]],
                 mean_label, aligned_TSS[[3]], aligned_TTS[[3]],
                 median_label, aligned_TSS[[2]], aligned_TTS[[2]],
                 ncol=3, rel_heights=c(6,15,15), rel_widths=c(1,15,15))
RefSeq_halfTitled=plot_grid(blank_plot, RefSeq,
                            labels=c("RefSeq genes"), label_size=16, label_fontface="bold", label_x=0.5, hjust=0.5, label_y=0.5, vjust=0.5,
                            ncol=1, rel_heights=c(4, 36))
RefSeq_titled=plot_grid(RefSeq_halfTitled, blank_plot,
                        labels=c("", "Distance from site (bp)"), label_size=10, label_fontface="plain", label_x=0.5, hjust=0.5, label_y=0.5, vjust=0.5,
                        ncol=1, rel_heights=c(4+36 ,2))

tvd=plot_grid(NCBI_titled, blank_plot, ENCODE_titled, blank_plot, RefSeq_titled, blank_plot, heatmap_legend,
              labels=c("", "", "", "", ""), label_size=8,
              ncol=7, rel_widths=c(1+15+15+15, 2, 1+15+15, 2, 1+15+15, 2, 5))

ggsave("VGP_TVD.png", tvd, units="mm", width=330, height=330/180*(50*((4+36+2)/39.5)))

tvd_cut=plot_grid(NCBI_titled, blank_plot, ENCODE_titled, blank_plot, heatmap_legend,
                  labels=c("", "", "", "", ""), label_size=8,
                  ncol=5, rel_widths=c(1+15+15+15, 2, 1+15+15, 2, 5))

ggsave("VGP_TVD_cut.png", tvd_cut, units="mm", width=330/(1+15+15+15+ 2+ 1+15+15+ 2+ 1+15+15+ 2+ 5)*(1+15+15+15+ 2+ 1+15+15+ 2+ 5), height=330/180*(50*((4+36+2)/39.5)))


# 위에서 human gene - line section 돌리고 와 ○
library(png)
blank = ggplot() + theme_void()
human_thin<-readPNG("human_thin.png")
zf_thin<-readPNG("zebra_finch_thin.png")

T2T_line_data = fread("MP_TSS_T2T-CHM13v2.0.tsv")
T2T_line_data[V4 == "-", V2 := V2 * -1] # If reverse (-) strand, flip the coordinates
T2T_line_data = data.frame(summarise(group_by(T2T_line_data, V2),"Mean"=mean(V3),"Median"=median(V3),"SD"=sd(V3)))
T2T_line_data=cbind(T2T_line_data, "T2T-CHM13v2.0")
colnames(T2T_line_data)=c("Distance","Mean","Median","SD","Assembly")

hg38_line_data = fread("MP_TSS_GRCh38.p14.tsv")
hg38_line_data[V4 == "-", V2 := V2 * -1] # If reverse (-) strand, flip the coordinates
hg38_line_data = data.frame(summarise(group_by(hg38_line_data, V2),"Mean"=mean(V3),"Median"=median(V3),"SD"=sd(V3)))
hg38_line_data=cbind(hg38_line_data, "GRCh38.p14")
colnames(hg38_line_data)=c("Distance","Mean","Median","SD","Assembly")

plot_line_mean_simple<-function(df, xmin, xmax, image, TEXT){
  ggplot(df, aes(x=Distance, y=Mean, group=Assembly)) +
    annotation_raster(image, xmin=-Inf, ymin=-Inf, xmax=Inf, ymax=Inf) +
    geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=-1000, xmax=-250), color="transparent", fill="#8FB4C6", alpha=.1) +
    geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=250, xmax=1000), color="transparent", fill="#CFE1EC", alpha=.1) +
    geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=-250, xmax=-150), color="transparent", fill="#7CCBA2", alpha=.1) +
    geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=150, xmax=250), color="transparent", fill="#B7E6A5", alpha=.1) +
    geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=-150, xmax=35), color="transparent", fill="#F07462", alpha=.1) +
    geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=35, xmax=150), color="transparent", fill="#FABF7B", alpha=.1) +
    geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=-35, xmax=0), color="transparent", fill="#AB1866", alpha=.1) +
    geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=0, xmax=35), color="transparent", fill="#D12959", alpha=.1) +
    geom_line(linewidth=.1, alpha=1) +
    #geom_vline(xintercept=-1000, color="#003147", linetype="dashed", linewidth=.1) +
    #geom_vline(xintercept=1000, color="#003147", linetype="dashed", linewidth=.1) +
    #geom_vline(xintercept=-250, color="#7CCBA2", linetype="dashed", linewidth=.1) +
    #geom_vline(xintercept=250, color="#7CCBA2", linetype="dashed", linewidth=.1) +
    #geom_vline(xintercept=-150, color="#F07462", linetype="dashed", linewidth=.1) +
    #geom_vline(xintercept=150, color="#F07462", linetype="dashed", linewidth=.1) +
    #geom_vline(xintercept=-35, color="#6E005F", linetype="dashed", linewidth=.1) +
    #geom_vline(xintercept=35, color="#6E005F", linetype="dashed", linewidth=.1) +
    #annotate("text", y=100, x=2500, vjust=1, hjust="inward", label = TEXT, color = "black", size = 1.94) + 
    scale_color_identity() +
    ylim(0,100) + scale_x_continuous(labels = label_comma(), breaks = c(xmin, 0, xmax), limits = c(xmin, xmax)) +
    labs(y="Mean MP (%)", x="Distance from TSS (bp)") + 
    theme(legend.position="none",
          panel.background = element_blank(),
          panel.grid = element_blank(),
          plot.background = element_blank(),
          axis.line.y.left = element_line(linewidth=0.1, color="black"),
          axis.line.x.bottom = element_line(linewidth=0.1, color="black"),
          axis.ticks = element_blank(),
          axis.text = element_text(size=10, color="black"),
          axis.title = element_text(size=10, color="black"))
}

T2T_line=plot_line_mean_simple(T2T_line_data,-2500,2500,human_thin,"T2T-CHM13v2.0")
T2T_lineOnly = align_plots(T2T_line + theme(plot.margin = margin(t=0,b=5.5,l=5.5,r=5.5)))
hg38_line=plot_line_mean_simple(hg38_line_data,-2500,2500,human_thin, "GRCh38.p14")
hg38_lineOnly = align_plots(hg38_line + theme(plot.margin = margin(t=0,b=5.5,l=5.5,r=5.5)))
gene_lines=plot_grid(blank, blank, T2T_lineOnly[[1]], hg38_lineOnly[[1]], ncol=2, rel_heights=c(0.5,5),
                     labels=c("hs1 (T2T v1)", "hg38 (Traditional reference)"), label_fontface="plain", label_size=12, label_x=0.5, hjust=0.5, label_y=0.5, vjust=0.5)
gene_lines_moreTitled=plot_grid(blank, gene_lines, ncol=1, rel_heights=c(0.25, 5.5),
                                labels=c("Gene body start site"), label_fontface="bold", label_size=16, label_x=0, hjust=0, label_y=0, vjust=0)

ggsave("VGP_gene_lines.png", gene_lines_moreTitled, units="mm", width=300, height=(300/30.73*(55*(5.75/10))/2) )

#

load("pvalues_human_gene.RData")
# 위에서 human gene - line section 및 box & viol section 돌리고 와 

T2T_line=plot_line_mean_simple(T2T_line_data,-2500,2500,human_thin, "T2T-CHM13v2.0")
T2T_viol=plot_viol_simple(T2T_box_data_multicolorized)
T2T_vsline=plot_VSline(T2T_tvd_data)
T2T_tile=plot_pvalue_heatmap(T2T_pvalues)

T2T_lineNviol = align_plots(T2T_line + theme(plot.margin = margin(t=0,b=5.5,l=5.5,r=5.5)), T2T_viol + theme(plot.margin = margin(t=0,b=0,l=5.5,r=5.5)), align="h", axis="tb")
T2T_vsNtile = align_plots(T2T_vsline + theme(plot.margin = margin(b=0)), T2T_tile + theme(plot.margin = margin(t=0)), align="v", axis="lr")

hg38_line=plot_line_mean_simple(hg38_line_data,-2500,2500,human_thin, "GRCh38.p14")
hg38_viol=plot_viol_simple(hg38_box_data_multicolorized)
hg38_vsline=plot_VSline(hg38_tvd_data)
hg38_tile=plot_pvalue_heatmap(hg38_pvalues)

hg38_lineNviol = align_plots(hg38_line + theme(plot.margin = margin(t=0,b=5.5,l=5.5,r=5.5)), hg38_viol + theme(plot.margin = margin(t=0,b=0,l=5.5,r=5.5)), align="h", axis="tb")
hg38_vsNtile = align_plots(hg38_vsline + theme(plot.margin = margin(b=0)), hg38_tile + theme(plot.margin = margin(t=0)), align="v", axis="lr")

row1=plot_grid(T2T_lineNviol[[1]], T2T_lineNviol[[2]], ncol=2, rel_widths=c(3.5,6.5))
human_gene_T2T=plot_grid(blank, row1, T2T_vsNtile[[1]], T2T_vsNtile[[2]], ncol=1, rel_heights=c(1,5,3,1),
                         labels=c("Gene start site (T2T-CHM13v2.0)"), label_fontface="plain", label_size=6, label_x=0.5, hjust=0.5)
row1=plot_grid(hg38_lineNviol[[1]], hg38_lineNviol[[2]], ncol=2, rel_widths=c(3.5,6.5))
human_gene_hg38=plot_grid(blank, row1, hg38_vsNtile[[1]], hg38_vsNtile[[2]], ncol=1, rel_heights=c(1,5,3,1),
                          labels=c("Gene start site (GRCh38.p14)"), label_fontface="plain", label_size=6, label_x=0.5, hjust=0.5)
human_gene_comp=plot_grid(blank, blank, blank, human_gene_T2T, blank, human_gene_hg38, ncol=3, rel_heights=c(0.5, 1+5+3+1), rel_widths=c(20,1,20),
                     labels=c("", "", "", "d", "", "e"), label_size=8)
human_gene_comp_wLegend=plot_grid(human_gene_comp, tile_legend, ncol=1, rel_heights=c(55, 10))

ggsave("VGP_gene_comp.png", human_gene_comp_wLegend, units="mm", width=300, height=300/180*(55+10))


# 위에서 human transcript - line section 돌리고 와. ○

library(png)
blank = ggplot() + theme_void()
human_thin<-readPNG("human_thin.png")
zf_thin<-readPNG("zebra_finch_thin.png")

T2T_line_data = fread("MP_transcriptTSS_T2T-CHM13v2.0.tsv")
T2T_line_data[V4 == "-", V2 := V2 * -1] # If reverse (-) strand, flip the coordinates
T2T_line_data = data.frame(summarise(group_by(T2T_line_data, V2),"Mean"=mean(V3),"Median"=median(V3),"SD"=sd(V3)))
T2T_line_data=cbind(T2T_line_data, "T2T-CHM13v2.0")
colnames(T2T_line_data)=c("Distance","Mean","Median","SD","Assembly")

hg38_line_data = fread("MP_transcriptTSS_GRCh38.p14.tsv")
hg38_line_data[V4 == "-", V2 := V2 * -1] # If reverse (-) strand, flip the coordinates
hg38_line_data = data.frame(summarise(group_by(hg38_line_data, V2),"Mean"=mean(V3),"Median"=median(V3),"SD"=sd(V3)))
hg38_line_data=cbind(hg38_line_data, "GRCh38.p14")
colnames(hg38_line_data)=c("Distance","Mean","Median","SD","Assembly")

plot_line_mean_simple<-function(df, xmin, xmax, image, TEXT){
  ggplot(df, aes(x=Distance, y=Mean, group=Assembly)) +
    annotation_raster(image, xmin=-Inf, ymin=-Inf, xmax=Inf, ymax=Inf) +
    geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=-1000, xmax=-250), color="transparent", fill="#8FB4C6", alpha=.1) +
    geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=250, xmax=1000), color="transparent", fill="#CFE1EC", alpha=.1) +
    geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=-250, xmax=-150), color="transparent", fill="#7CCBA2", alpha=.1) +
    geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=150, xmax=250), color="transparent", fill="#B7E6A5", alpha=.1) +
    geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=-150, xmax=35), color="transparent", fill="#F07462", alpha=.1) +
    geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=35, xmax=150), color="transparent", fill="#FABF7B", alpha=.1) +
    geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=-35, xmax=0), color="transparent", fill="#AB1866", alpha=.1) +
    geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=0, xmax=35), color="transparent", fill="#D12959", alpha=.1) +
    geom_line(linewidth=.1, alpha=1) +
    #geom_vline(xintercept=-1000, color="#003147", linetype="dashed", linewidth=.1) +
    #geom_vline(xintercept=1000, color="#003147", linetype="dashed", linewidth=.1) +
    #geom_vline(xintercept=-250, color="#7CCBA2", linetype="dashed", linewidth=.1) +
    #geom_vline(xintercept=250, color="#7CCBA2", linetype="dashed", linewidth=.1) +
    #geom_vline(xintercept=-150, color="#F07462", linetype="dashed", linewidth=.1) +
    #geom_vline(xintercept=150, color="#F07462", linetype="dashed", linewidth=.1) +
    #geom_vline(xintercept=-35, color="#6E005F", linetype="dashed", linewidth=.1) +
    #geom_vline(xintercept=35, color="#6E005F", linetype="dashed", linewidth=.1) +
    #annotate("text", y=100, x=2500, vjust=1, hjust="inward", label = TEXT, color = "black", size = 1.94) + 
    scale_color_identity() +
    ylim(0,100) + scale_x_continuous(labels = label_comma(), breaks = c(xmin, 0, xmax), limits = c(xmin, xmax)) +
    labs(y="Mean MP (%)", x="Distance from TSS (bp)") + 
    theme(legend.position="none",
          panel.background = element_blank(),
          panel.grid = element_blank(),
          plot.background = element_blank(),
          axis.line.y.left = element_line(linewidth=0.1, color="black"),
          axis.line.x.bottom = element_line(linewidth=0.1, color="black"),
          axis.ticks = element_blank(),
          axis.text = element_text(size=10, color="black"),
          axis.title = element_text(size=10, color="black"))
}


blank = ggplot() + theme_void()
human_thin<-readPNG("human_thin.png")
zf_thin<-readPNG("zebra_finch_thin.png")

T2T_line=plot_line_mean_simple(T2T_line_data,-2500,2500,human_thin,"T2T-CHM13v2.0")
T2T_lineOnly = align_plots(T2T_line + theme(plot.margin = margin(t=0,b=5.5,l=5.5,r=5.5)))
hg38_line=plot_line_mean_simple(hg38_line_data,-2500,2500,human_thin, "GRCh38.p14")
hg38_lineOnly = align_plots(hg38_line + theme(plot.margin = margin(t=0,b=5.5,l=5.5,r=5.5)))
gene_lines=plot_grid(blank, blank, T2T_lineOnly[[1]], hg38_lineOnly[[1]], ncol=2, rel_heights=c(0.5,5),
                     labels=c("hs1 (T2T v1)", "hg38 (Traditional reference)"), label_fontface="plain", label_size=12, label_x=0.5, hjust=0.5, label_y=0.25, vjust=0)
gene_lines_moreTitled=plot_grid(blank, gene_lines, ncol=1, rel_heights=c(0.25, 5.5),
                                labels=c("Transcript start site"), label_fontface="bold", label_size=16, label_x=0, hjust=0, label_y=0, vjust=0)

ggsave("VGP_transcript_lines.png", gene_lines_moreTitled, units="mm", width=300, height=(300/30.73*(55*(5.75/10))/2) )

#

load("pvalues_human_transcript.RData")
# 위에서 human transcript - line section 및 box & viol section 돌리고 와

T2T_line=plot_line_mean_simple(T2T_line_data,-2500,2500,human_thin, "T2T-CHM13v2.0")
T2T_viol=plot_viol_simple(T2T_box_data_multicolorized)
T2T_vsline=plot_VSline(T2T_tvd_data)
T2T_tile=plot_pvalue_heatmap(T2T_pvalues)

T2T_lineNviol = align_plots(T2T_line + theme(plot.margin = margin(t=0,b=5.5,l=5.5,r=5.5)), T2T_viol + theme(plot.margin = margin(t=0,b=0,l=5.5,r=5.5)), align="h", axis="tb")
T2T_vsNtile = align_plots(T2T_vsline + theme(plot.margin = margin(b=0)), T2T_tile + theme(plot.margin = margin(t=0)), align="v", axis="lr")

hg38_line=plot_line_mean_simple(hg38_line_data,-2500,2500,human_thin, "GRCh38.p14")
hg38_viol=plot_viol_simple(hg38_box_data_multicolorized)
hg38_vsline=plot_VSline(hg38_tvd_data)
hg38_tile=plot_pvalue_heatmap(hg38_pvalues)

hg38_lineNviol = align_plots(hg38_line + theme(plot.margin = margin(t=0,b=5.5,l=5.5,r=5.5)), hg38_viol + theme(plot.margin = margin(t=0,b=0,l=5.5,r=5.5)), align="h", axis="tb")
hg38_vsNtile = align_plots(hg38_vsline + theme(plot.margin = margin(b=0)), hg38_tile + theme(plot.margin = margin(t=0)), align="v", axis="lr")

row1=plot_grid(T2T_lineNviol[[1]], T2T_lineNviol[[2]], ncol=2, rel_widths=c(3.5,6.5))
human_gene_T2T=plot_grid(blank, row1, T2T_vsNtile[[1]], T2T_vsNtile[[2]], ncol=1, rel_heights=c(1,5,3,1),
                         labels=c("Transcript start site (T2T-CHM13v2.0)"), label_fontface="plain", label_size=6, label_x=0.5, hjust=0.5)
row1=plot_grid(hg38_lineNviol[[1]], hg38_lineNviol[[2]], ncol=2, rel_widths=c(3.5,6.5))
human_gene_hg38=plot_grid(blank, row1, hg38_vsNtile[[1]], hg38_vsNtile[[2]], ncol=1, rel_heights=c(1,5,3,1),
                          labels=c("Transcript start site (GRCh38.p14)"), label_fontface="plain", label_size=6, label_x=0.5, hjust=0.5)
human_gene_comp=plot_grid(blank, blank, blank, human_gene_T2T, blank, human_gene_hg38, ncol=3, rel_heights=c(0.5, 1+5+3+1), rel_widths=c(20,1,20),
                          labels=c("", "", "", "", "", ""), label_size=8)
human_gene_comp_wLegend=plot_grid(human_gene_comp, tile_legend, ncol=1, rel_heights=c(55, 10))

ggsave("VGP_transcript_comp.png", human_gene_comp_wLegend, units="mm", width=300, height=300/180*(55+10))


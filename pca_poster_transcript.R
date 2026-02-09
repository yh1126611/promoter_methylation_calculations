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
setwd("D:/R/WD/MP_TSS/pca/transcript")

files_TSS <- list.files(pattern = "^MP_transcriptTSS_.*\\.tsv$")
files_TTS <- list.files(pattern = "^MP_transcriptTES_.*\\.tsv$")

pca_preprocessor = function(file){
  df=fread(file)
  df[V4 == "-", V2 := V2 * -1] # If reverse (-) strand, flip the coordinates
  df = data.frame(summarise(group_by(df, V2),"Mean"=mean(V3),"Median"=median(V3),"SD"=sd(V3)))
  colnames(df)=c("Distance","Mean","Median","SD")
  # 1. Create the full sequence of distances
  all_distances <- seq(-10000, 10000)
  # 2. Initialize an empty matrix for results
  result_mat <- matrix(NA, nrow=3, ncol=length(all_distances))
  rownames(result_mat) <- c("Mean", "Median", "SD")
  colnames(result_mat) <- as.character(all_distances)
  # 3. Match your distances to the columns of the matrix
  matched_cols <- match(df$Distance, all_distances)
  # 4. Fill in the matrix with your values
  result_mat["Mean", matched_cols]   <- df$Mean
  result_mat["Median", matched_cols] <- df$Median
  result_mat["SD", matched_cols]     <- df$SD
  # 5. Convert to a data frame if you prefer
  species_summary <- as.data.frame(result_mat)
  # 6. Add species name and class
  species_name = str_match(file, "^MP_transcriptT[ES]S_(.*)\\.tsv$")[,2]
  classification = substr(species_name, 1, 1)
  species_summary$"Assembly" = species_name
  species_summary$"Class" = classification
  return(species_summary)
} # 이 부분에서 종마다 class/order를 assign해서 df에 넣어줌.

normalize_numbers <- function(df){
  # Identify numeric columns
  numeric_cols <- sapply(df, is.numeric)
  # If there are no numeric columns, return the original dataframe
  if (all(!numeric_cols)) return(df)
  # Row-wise min-max normalization for numeric columns, scaling to [0, 100]
  normalized_numeric <- t(apply(df[, numeric_cols, drop = FALSE], 1, function(x) {
    if (max(x) == min(x)) {
      rep(0.5, length(x))  # Avoid division by zero; assign 0.5 if all values are equal
    } else {
      (x - min(x)) / (max(x) - min(x))
    }
  }))
  # Ensure the result is a data frame with proper column names
  normalized_numeric <- as.data.frame(normalized_numeric)
  colnames(normalized_numeric) <- colnames(df)[numeric_cols]
  # Recombine with non-numeric columns
  df_normalized <- cbind(
    normalized_numeric,
    df[, !numeric_cols, drop = FALSE]
  )
  # Reorder columns to match the original dataframe
  df_normalized <- df_normalized[, colnames(df)]
  return(df_normalized)
}

data_matrix_TSS = bind_rows(lapply(files_TSS, pca_preprocessor)) # ⌛
data_matrix_TTS = bind_rows(lapply(files_TTS, pca_preprocessor)) # ⌛
# Now, in data_matrix, columns 1-20,001 are MP, 20,002 is Assembly name and 20,003 is class.
#save.image("D:/R/WD/MP_TSS/pca/transcript/pca_transcript_data.RData")
load("D:/R/WD/MP_TSS/pca/transcript/pca_transcript_data.RData")

# Add EGAP
data_matrix_TTS=bind_rows(data_matrix_TTS, pca_preprocessor("MP_transcriptTES_bTaeGut7v0.4_MT_rDNA.EGAPx.v0.1_mat.tsv"))
data_matrix_TSS=bind_rows(data_matrix_TSS, pca_preprocessor("MP_transcriptTSS_bTaeGut7v0.4_MT_rDNA.EGAPx.v0.1_mat.tsv"))

# Rename columns
colnames(data_matrix_TSS)[1:20001] <- paste0(colnames(data_matrix_TSS)[1:20001], "TSS")
colnames(data_matrix_TTS)[1:20001] <- paste0(colnames(data_matrix_TTS)[1:20001], "TTS")
# Add Statistic column if not present (Mean, Median, SD)
data_matrix_TSS$Statistic <- rep(c("Mean", "Median", "SD"), length(unique(data_matrix_TSS$Assembly)))
data_matrix_TTS$Statistic <- rep(c("Mean", "Median", "SD"), length(unique(data_matrix_TTS$Assembly)))
# Merge
data_matrix <- merge(
  data_matrix_TSS, data_matrix_TTS,
  by = c("Assembly", "Class", "Statistic"),
  all = TRUE
)
# Send columns "Class" and "Statistic" to the back
data_matrix <- data_matrix[, c(setdiff(colnames(data_matrix), c("Class", "Statistic")), "Assembly", "Statistic")]

# ↑ 한 번만

data_matrix_normalized = normalize_numbers(data_matrix)
data_matrix_TSS_normalized = normalize_numbers(data_matrix_TSS)
data_matrix_TTS_normalized = normalize_numbers(data_matrix_TTS)
#data_matrix_normalized=data_matrix # If you DON'T want to normalize
#data_matrix_TSS_normalized=data_matrix_TSS # If you DON'T want to normalize
#data_matrix_TTS_normalized=data_matrix_TTS # If you DON'T want to normalize

comprehensive <- fread("C:/Users/biopop/Documents/Downloads/MP analysis for all VGP Freeze1.0 species - Sheet1.tsv")
data_matrix_normalized <- merge(data_matrix_normalized, comprehensive, by="Assembly", all.x=TRUE)
data_matrix_TSS_normalized <- merge(data_matrix_TSS_normalized, comprehensive, by="Assembly", all.x=TRUE)
data_matrix_TTS_normalized <- merge(data_matrix_TTS_normalized, comprehensive, by="Assembly", all.x=TRUE)
# Now, in data_matrix, "Assembly" becomes column 1, columns 2-20,002 are MP, 20,011 is color.

# Remove outliers
data_matrix_normalized <- filter(data_matrix_normalized, !`Common name` %in% c("Wood mouse", "African elephant", "Yellowfin tuna", "Crucian carp", "Northern goshawk", "European flounder") & `Assembly` != "bTaeGut7")
data_matrix_TSS_normalized <- filter(data_matrix_TSS_normalized, !`Common name` %in% c("Wood mouse", "African elephant", "Yellowfin tuna", "Crucian carp", "Northern goshawk", "European flounder") & `Assembly` != "bTaeGut7")
data_matrix_TTS_normalized <- filter(data_matrix_TTS_normalized, !`Common name` %in% c("Wood mouse", "African elephant", "Yellowfin tuna", "Crucian carp", "Northern goshawk", "European flounder") & `Assembly` != "bTaeGut7")
data_matrix_normalized$`Class (강)`[data_matrix_normalized$Species=="Homo sapiens"]="Human"
data_matrix_TSS_normalized$`Class (강)`[data_matrix_TSS_normalized$Species=="Homo sapiens"]="Human"
data_matrix_TTS_normalized$`Class (강)`[data_matrix_TTS_normalized$Species=="Homo sapiens"]="Human"

data_mean <- data_matrix_normalized[seq(1, nrow(data_matrix_normalized), by = 3), ]
data_median <- data_matrix_normalized[seq(2, nrow(data_matrix_normalized), by = 3), ]
data_SD <- data_matrix_normalized[seq(3, nrow(data_matrix_normalized), by = 3), ]
data_mean_TSS <- data_matrix_TSS_normalized[seq(1, nrow(data_matrix_TSS_normalized), by = 3), ]
data_median_TSS <- data_matrix_TSS_normalized[seq(2, nrow(data_matrix_TSS_normalized), by = 3), ]
data_SD_TSS <- data_matrix_TSS_normalized[seq(3, nrow(data_matrix_TSS_normalized), by = 3), ]
data_mean_TTS <- data_matrix_TTS_normalized[seq(1, nrow(data_matrix_TTS_normalized), by = 3), ]
data_median_TTS <- data_matrix_TTS_normalized[seq(2, nrow(data_matrix_TTS_normalized), by = 3), ]
data_SD_TTS <- data_matrix_TTS_normalized[seq(3, nrow(data_matrix_TTS_normalized), by = 3), ]


# Function for conducting & plotting PCA
plot_PCA <- function(pca_result, whole_data, title){
  explained_var <- (pca_result$sdev^2) / sum(pca_result$sdev^2)
  percent_var <- round(100 * explained_var, 1)
  
  assembly_size <- c("GRCh38.p14"=1,"T2T-CHM13v2.0"=1, "bTaeGut7"=1, "bTaeGut1.4"=1)
  
  scores <- as.data.frame(pca_result$x)
  scores$Assembly   <- whole_data$Assembly
  scores$Color      <- whole_data$Color
  scores$Class      <- whole_data$`Class (강)`
  scores$common_name <- whole_data$`Common name`
  
  assembly_color <- c(
    "Human"="#000000", "Mammalia"="#E69F00","Aves"="#00796B","Reptilia"="#A6D609",
    "Amphibia"="#984EA3","Actinopterygii"="#56B4E9","Sarcopterygii"="#A6761D","Chondrichthyes"="#0072B2"
  )
  
  scores_human <- subset(scores, Class == "Human")
  
  plotted <- ggplot(scores, aes(x = PC1, y = PC2)) +
    geom_text(aes(label = common_name), vjust = -1, size = 3, color = "transparent") +
    geom_point(aes(color = Class), size = 0.5) +
    geom_point(data = scores_human, aes(color = Class), size = 0.5) +
    scale_color_manual(values = assembly_color) +
    scale_size_manual(values = assembly_size, na.value = 0.4) +
    scale_shape_identity() +
    labs(
      x = paste0("PC1 (", percent_var[1], "%)"),
      y = paste0("PC2 (", percent_var[2], "%)"),
      title = title
    ) +
    coord_cartesian(clip = "off") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 5.5),
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.line.y.left = element_line(),
      axis.line.x.bottom = element_line(),
      panel.grid = element_blank(),
      axis.line = element_line(linewidth = 0.1),
      axis.ticks = element_blank(),
      axis.title = element_text(size = 5.5),
      axis.text = element_text(size = 5, color = "black")
    )
  
  return(plotted)
}

plot_PCA_byTissue <- function(pca_result,whole_data,title){
  explained_var <- (pca_result$sdev^2) / sum(pca_result$sdev^2)
  percent_var <- round(100 * explained_var, 1)
  
  scores <- as.data.frame(pca_result$x)
  scores$Assembly <- whole_data$Assembly
  scores$Color <- whole_data$Color
  scores$Tissue <- whole_data$"Tissue"
  scores$common_name <- whole_data$"Common name"
  tissue_color=c("Multiple"="#000000", "Blood"="#A50026", "Muscle"="#EEDD88", "Organ"="#364B9A", "Blood, muscle"="#F67E4B", "Blood, organ"="#762A83", "Muscle, organ"="#6EA6CD","Brain"="#1B7837","Fibroblast"="#009988","Fat"="#AAAA00","Other"="#888888","N/A"="#CCCCCC")
  # Plot
  ggplot(scores, aes(x=PC1, y=PC2, label = common_name)) +
    geom_text(vjust = -1, size = 3, color="transparent") + # Change vjust to 0 and color="transparent" to aes(color=Color) to display names instead of dots
    geom_point(aes(color=Tissue), size=0.4) + # Change aes(color=Color) to color="transparent" to display names instead of dots
    # ↑ If you want to fix size instead of setting it by Size column of scores data frame, remove "size=Size" inside the aes() argument and add "size=<your size>" inside the geom_point() argument
    scale_color_manual(values=tissue_color) +
    scale_shape_identity() +
    labs(x = paste0("PC1 (", percent_var[1], "%)"), y = paste0("PC2 (", percent_var[2], "%)"), title=title) +
    coord_cartesian(clip = "off") + #To prevent clipping-off of geom_texts that go out of boundaries
    theme(
      legend.position="none",
      plot.title = element_text(size=5.5),
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.line.y.left = element_line(),
      axis.line.x.bottom = element_line(),
      panel.grid = element_blank(),
      axis.line=element_line(linewidth=0.1),
      axis.ticks = element_blank(),
      axis.title=element_text(size=5.5),
      axis.text=element_text(size=5, color="black")
    )->plotted
  return(plotted)
}


# Create legend
class_order <- c("Human", "Non-human mammal", "Bird","Reptile", "Amphibian", "Lobe-finned fish", "Ray-finned fish", "Shark")
df <- data.frame(X = 1,Y = 1,Class = c("Human", "Non-human mammal", "Bird","Reptile", "Amphibian", "Lobe-finned fish", "Ray-finned fish", "Shark"))
df$Class <- factor(df$Class, levels = class_order)
fake <- ggplot(df, aes(x = X, y = Y, color = Class)) +
  geom_point(size = 3) +
  scale_color_manual(name="Class:", values = c(
    "Human" = "#000000",
    "Non-human mammal" = "#E69F00",
    "Bird" = "#00796B",
    "Reptile" = "#A6D609",
    "Amphibian" = "#984EA3",
    "Lobe-finned fish" = "#A6761D",
    "Ray-finned fish" = "#56B4E9",
    "Shark" = "#0072B2"
  )) +
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 2))) +
  theme(
    #legend.position = "top",
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill="transparent", color=NA),
    legend.title = element_text(size=6, face="bold"),
    legend.text = element_text(size=6),
    legend.key.size=unit(.75,"lines")
  )
legend = as_ggplot(get_legend(fake));legend

# Create legend for tissue
class_order <- c("Multiple","Blood","Muscle","Organ","Blood, muscle","Blood, organ","Muscle, organ","Brain","Fibroblast","Fat","Other","N/A")
df <- data.frame(X = 1,Y = 1,Class = c("Multiple","Blood","Muscle","Organ","Blood, muscle","Blood, organ","Muscle, organ","Brain","Fibroblast","Fat","Other","N/A"))
df$Class <- factor(df$Class, levels = class_order)
fake <- ggplot(df, aes(x = X, y = Y, color = Class)) +
  geom_point(size = 3) +
  scale_color_manual(name="Tissue:", values = c("Multiple"="#000000", "Blood"="#A50026", "Muscle"="#EEDD88", "Organ"="#364B9A", "Blood, muscle"="#F67E4B", "Blood, organ"="#762A83", "Muscle, organ"="#6EA6CD","Brain"="#1B7837","Fibroblast"="#009988","Fat"="#AAAA00","Other"="#888888","N/A"="#CCCCCC")) +
  guides(color = guide_legend(ncol = 3, override.aes = list(size = 2))) +
  theme(
    #legend.position = "top",
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill="transparent", color=NA),
    legend.title = element_text(size=6, face="bold"),
    legend.text = element_text(size=6),
    legend.key.size=unit(.75, "lines")
  )
legend_tissue = as_ggplot(get_legend(fake));legend_tissue


## PCA on tissue ##

# PCA on mean
pca_result_mean = prcomp(data_mean[,2:20002],scale.=TRUE)
pca_plot_mean = plot_PCA_byTissue(pca_result_mean,data_mean, "Mean")
pca_result_mean_TSS = prcomp(data_mean_TSS[,2:20002],scale.=TRUE)
pca_plot_mean_TSS = plot_PCA_byTissue(pca_result_mean_TSS, data_mean_TSS, "Mean"); #pca_plot_mean_TSS
pca_result_mean_TTS = prcomp(data_mean_TTS[,2:20002],scale.=TRUE)
pca_plot_mean_TTS = plot_PCA_byTissue(pca_result_mean_TTS, data_mean_TTS, "Mean"); #pca_plot_mean_TTS

#PCA on median
pca_result_median = prcomp(data_median[,2:20002],scale.=TRUE)
pca_plot_median = plot_PCA_byTissue(pca_result_median, data_median, "Median"); #pca_plot_median
pca_result_median_TSS = prcomp(data_median_TSS[,2:20002],scale.=TRUE)
pca_plot_median_TSS = plot_PCA_byTissue(pca_result_median_TSS, data_median_TSS, "Median"); #pca_plot_median_TSS
pca_result_median_TTS = prcomp(data_median_TTS[,2:20002],scale.=TRUE)
pca_plot_median_TTS = plot_PCA_byTissue(pca_result_median_TTS, data_median_TTS, "Median"); #pca_plot_median_TTS

#PCA on SD
pca_result_SD = prcomp(data_SD[,2:20002],scale.=TRUE)
pca_plot_SD = plot_PCA_byTissue(pca_result_SD, data_SD, "SD"); #pca_plot_SD
pca_result_SD_TSS = prcomp(data_SD_TSS[,2:20002],scale.=TRUE)
pca_plot_SD_TSS = plot_PCA_byTissue(pca_result_SD_TSS, data_SD_TSS, "SD"); #pca_plot_SD_TSS
pca_result_SD_TTS = prcomp(data_SD_TTS[,2:20002],scale.=TRUE)
pca_plot_SD_TTS = plot_PCA_byTissue(pca_result_SD_TTS, data_SD_TTS, "SD"); #pca_plot_SD_TTS

# Align & Combine
aligned_TSS = align_plots(
  pca_plot_mean_TSS + theme(plot.margin = margin(0,0,0,0)),
  pca_plot_median_TSS + theme(plot.margin = margin(0,0,0,0)),
  pca_plot_SD_TSS + theme(plot.margin = margin(0,0,0,0)),
  align = "h", axis = "tb")
aligned_TTS = align_plots(
  pca_plot_mean_TTS + theme(plot.margin = margin(0,0,0,0)),
  pca_plot_median_TTS + theme(plot.margin = margin(0,0,0,0)),
  pca_plot_SD_TTS + theme(plot.margin = margin(0,0,0,0)),
  align = "h", axis = "tb")
aligned = align_plots(
  pca_plot_mean + theme(plot.margin = margin(0,0,0,0)),
  pca_plot_median + theme(plot.margin = margin(0,0,0,0)),
  pca_plot_SD + theme(plot.margin = margin(0,0,0,0)),
  align = "h", axis = "tb")


blank_plot <- ggplot() + theme_void()

tss_part = plot_grid(
  labels=c("", "Transcript start site (TSS)", ""),
  label_size=6, label_x=0.5, hjust=0.5, label_y=1, vjust=0.5,
  blank_plot, blank_plot, blank_plot,
  aligned_TSS[[1]],  aligned_TSS[[2]], aligned_TSS[[3]],
  ncol=3, rel_heights = c(1,11))
tts_part = plot_grid(
  labels=c("", "Transcript termination site (TTS)", ""),
  label_size=6, label_x=0.5, hjust=0.5, label_y=1, vjust=0.5,
  blank_plot, blank_plot, blank_plot,
  aligned_TTS[[1]],  aligned_TTS[[2]], aligned_TTS[[3]],
  ncol=3, rel_heights = c(1,11))
tsstts_part = plot_grid(
  labels=c("", "Transcript start and termination site (TSS + TTS)", ""),
  label_size=6, label_x=0.5, hjust=0.5, label_y=1, vjust=0.5,
  blank_plot, blank_plot, blank_plot,
  aligned[[1]],  aligned[[2]], aligned[[3]],
  ncol=3, rel_heights = c(1,11))

tissue = plot_grid(labels=c("", "", "", "", "d", "", "e", "", "", "", "", "", "f"), label_size=8, label_y=1, vjust=0.5, blank_plot, blank_plot, blank_plot, blank_plot, tss_part, blank_plot, tts_part, blank_plot, blank_plot, blank_plot, blank_plot, blank_plot, tsstts_part, blank_plot, legend_tissue, blank_plot, ncol=4, rel_widths=c(30,1,30,1), rel_heights=c(1,12,1,12), align="h", axis="t"); all
#ggsave("pca_tissue_pres.png", all, width=183, height=73.79, units="mm")



## PCA on real ##

# PCA on mean
pca_result_mean = prcomp(data_mean[,2:20002],scale.=TRUE)
pca_plot_mean = plot_PCA(pca_result_mean,data_mean, "Mean")
pca_result_mean_TSS = prcomp(data_mean_TSS[,2:20002],scale.=TRUE)
pca_plot_mean_TSS = plot_PCA(pca_result_mean_TSS, data_mean_TSS, "Mean"); #pca_plot_mean_TSS
pca_result_mean_TTS = prcomp(data_mean_TTS[,2:20002],scale.=TRUE)
pca_plot_mean_TTS = plot_PCA(pca_result_mean_TTS, data_mean_TTS, "Mean"); #pca_plot_mean_TTS

#PCA on median
pca_result_median = prcomp(data_median[,2:20002],scale.=TRUE)
pca_plot_median = plot_PCA(pca_result_median, data_median, "Median"); #pca_plot_median
pca_result_median_TSS = prcomp(data_median_TSS[,2:20002],scale.=TRUE)
pca_plot_median_TSS = plot_PCA(pca_result_median_TSS, data_median_TSS, "Median"); #pca_plot_median_TSS
pca_result_median_TTS = prcomp(data_median_TTS[,2:20002],scale.=TRUE)
pca_plot_median_TTS = plot_PCA(pca_result_median_TTS, data_median_TTS, "Median"); #pca_plot_median_TTS

#PCA on SD
pca_result_SD = prcomp(data_SD[,2:20002],scale.=TRUE)
pca_plot_SD = plot_PCA(pca_result_SD, data_SD, "SD"); #pca_plot_SD
pca_result_SD_TSS = prcomp(data_SD_TSS[,2:20002],scale.=TRUE)
pca_plot_SD_TSS = plot_PCA(pca_result_SD_TSS, data_SD_TSS, "SD"); #pca_plot_SD_TSS
pca_result_SD_TTS = prcomp(data_SD_TTS[,2:20002],scale.=TRUE)
pca_plot_SD_TTS = plot_PCA(pca_result_SD_TTS, data_SD_TTS, "SD"); #pca_plot_SD_TTS

# Align & Combine
aligned_TSS = align_plots(
  pca_plot_mean_TSS + theme(plot.margin = margin(0,0,0,0)),
  pca_plot_median_TSS + theme(plot.margin = margin(0,0,0,0)),
  pca_plot_SD_TSS + theme(plot.margin = margin(0,0,0,0)),
  align = "h", axis = "tb")
aligned_TTS = align_plots(
  pca_plot_mean_TTS + theme(plot.margin = margin(0,0,0,0)),
  pca_plot_median_TTS + theme(plot.margin = margin(0,0,0,0)),
  pca_plot_SD_TTS + theme(plot.margin = margin(0,0,0,0)),
  align = "h", axis = "tb")
aligned = align_plots(
  pca_plot_mean + theme(plot.margin = margin(0,0,0,0)),
  pca_plot_median + theme(plot.margin = margin(0,0,0,0)),
  pca_plot_SD + theme(plot.margin = margin(0,0,0,0)),
  align = "h", axis = "tb")


blank_plot <- ggplot() + theme_void()

tss_part = plot_grid(
  labels=c("", "Transcript start site (TSS)", ""),
  label_size=6, label_x=0.5, hjust=0.5, label_y=1, vjust=0.5,
  blank_plot, blank_plot, blank_plot,
  aligned_TSS[[1]],  aligned_TSS[[2]], aligned_TSS[[3]],
  ncol=3, rel_heights = c(1,11))
tts_part = plot_grid(
  labels=c("", "Transcript termination site (TTS)", ""),
  label_size=6, label_x=0.5, hjust=0.5, label_y=1, vjust=0.5,
  blank_plot, blank_plot, blank_plot,
  aligned_TTS[[1]],  aligned_TTS[[2]], aligned_TTS[[3]],
  ncol=3, rel_heights = c(1,11))
tsstts_part = plot_grid(
  labels=c("", "Transcript start and termination site (TSS + TTS)", ""),
  label_size=6, label_x=0.5, hjust=0.5, label_y=1, vjust=0.5,
  blank_plot, blank_plot, blank_plot,
  aligned[[1]],  aligned[[2]], aligned[[3]],
  ncol=3, rel_heights = c(1,11))

real = plot_grid(labels=c("", "", "", "", "a", "", "b", "", "", "", "", "", "c"), label_size=8, label_y=1, vjust=0.5, blank_plot, blank_plot, blank_plot, blank_plot, tss_part, blank_plot, tts_part, blank_plot, blank_plot, blank_plot, blank_plot, blank_plot, tsstts_part, blank_plot, legend, blank_plot, ncol=4, rel_widths=c(30,1,30,1), rel_heights=c(1,12,1,12))
#ggsave("pca_pres.png", all, width=183, height=73.79, units="mm")


all = plot_grid(real, tissue, ncol=1)
ggsave("pca_pres_transcript.png", all, width=183, height=147.58, units="mm")
ggsave("pca_pres_transcript.pdf", all, width=183, height=147.58, units="mm")



## PCA on mammals only ##
plot_PCA <- function(pca_result,whole_data,title){
  assembly_colors = c("mCynVol1"="#BB7C23","mGorGor1"="#FF5C00","mMacNem1"="#BB7C23","mNycCou1"="#BB7C23","mPanPan1"="#FF5C00","mPanTro3"="#FF5C00","mPonAbe1"="#FF5C00","mPonPyg2"="#FF5C00","mSymSyn1"="#FF5C00","T2T-CHM13v2.0"="#E6194B","GRCh38.p14"="#E6194B",
                      "mChiNiv1"="#5E270B", "mMarFla1"="#5E270B", "mOchPri1"="#5E270B", "mEriEur2"="#5E270B", "mSorAra1"="#5E270B")
  assembly_size = c("T2T-CHM13v2.0"=1,"GRCh38.p14"=1)
  explained_var <- (pca_result$sdev^2) / sum(pca_result$sdev^2)
  percent_var <- round(100 * explained_var, 1)
  
  scores <- as.data.frame(pca_result$x)
  scores$Assembly <- whole_data$Assembly
  scores$Color <- whole_data$Color
  scores$common_name <- whole_data$"Common name"
  # Plot
  ggplot(scores, aes(x=PC1, y=PC2, label = common_name)) +
    geom_text(vjust = 0, size = 2, color="transparent") + # Change vjust to 0 and color="transparent" to aes(color=Color) to display names instead of dots
    geom_point(aes(color=Assembly, size=Assembly)) + # Change aes(color=Color) to color="transparent" to display names instead of dots
    # ↑ If you want to fix size instead of setting it by Size column of scores data frame, remove "size=Size" inside the aes() argument and add "size=<your size>" inside the geom_point() argument
    # scale_color_identity() + # Use the actual hex codes from the Color column
    scale_color_manual(values=assembly_colors, na.value="#AAAAAA") +
    scale_size_manual(values=assembly_size, na.value=0.4) +
    scale_shape_identity() +
    labs(x = paste0("PC1 (", percent_var[1], "%)"), y = paste0("PC2 (", percent_var[2], "%)"), title=title) +
    coord_cartesian(clip = "off") + #To prevent clipping-off of geom_texts that go out of boundaries
    #geom_text(data=subset(scores, Assembly %in% c("GRCh38.p14", "T2T-CHM13v2.0")), aes(label="★", color=Color), size=5)  +
    theme(
      legend.position="none",
      plot.title = element_text(size=5.5),
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.line.y.left = element_line(),
      axis.line.x.bottom = element_line(),
      panel.grid = element_blank(),
      axis.line=element_line(linewidth=0.1),
      axis.ticks = element_blank(),
      axis.title=element_text(size=5.5),
      axis.text=element_text(size=5, color="black")
    )->plotted
  return(plotted)
}
df = data.frame("Col"=c("Human","Non-human ape","Primatomorpha","Rodent"))
df$Col <- factor(df$Col, levels=c("Human","Non-human ape","Primatomorpha","Rodent"))
fake <- ggplot(df, aes(x = Col, y = Col, color = Col)) +
  geom_point(size = 3) +
  scale_color_manual(name="Lineage:", values = c(
    "Human" = "#E6194B",
    "Non-human ape" = "#FF5C00",
    "Primatomorpha" = "#BB7C23",
    "Rodent"="#5E270B"
  )) +
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 2))) +
  theme(
    #legend.position = "top",
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill="transparent", color=NA),
    legend.title = element_text(size=6, face="bold"),
    legend.text = element_text(size=6),
    legend.key.size=unit(.75,"lines")
  ); legend = as_ggplot(get_legend(fake));legend

data_mean = data_mean[!is.na(data_mean$`Class (강)`) & data_mean$"Class (강)"=="Mammalia",]
data_mean_TSS = data_mean_TSS[!is.na(data_mean_TSS$`Class (강)`) & data_mean_TSS$"Class (강)"=="Mammalia",]
data_mean_TTS = data_mean_TTS[!is.na(data_mean_TTS$`Class (강)`) & data_mean_TTS$"Class (강)"=="Mammalia",]
data_median = data_median[!is.na(data_median$`Class (강)`) & data_median$"Class (강)"=="Mammalia",]
data_median_TSS = data_median_TSS[!is.na(data_median_TSS$`Class (강)`) & data_median_TSS$"Class (강)"=="Mammalia",]
data_median_TTS = data_median_TTS[!is.na(data_median_TTS$`Class (강)`) & data_median_TTS$"Class (강)"=="Mammalia",]
data_SD = data_SD[!is.na(data_SD$`Class (강)`) & data_SD$"Class (강)"=="Mammalia",]
data_SD_TSS = data_SD_TSS[!is.na(data_SD_TSS$`Class (강)`) & data_SD_TSS$"Class (강)"=="Mammalia",]
data_SD_TTS = data_SD_TTS[!is.na(data_SD_TTS$`Class (강)`) & data_SD_TTS$"Class (강)"=="Mammalia",]

pca_result_mean = prcomp(data_mean[,2:20002],scale.=TRUE)
pca_plot_mean = plot_PCA(pca_result_mean,data_mean, "Mean")
pca_result_mean_TSS = prcomp(data_mean_TSS[,2:20002],scale.=TRUE)
pca_plot_mean_TSS = plot_PCA(pca_result_mean_TSS, data_mean_TSS, "Mean"); #pca_plot_mean_TSS
pca_result_mean_TTS = prcomp(data_mean_TTS[,2:20002],scale.=TRUE)
pca_plot_mean_TTS = plot_PCA(pca_result_mean_TTS, data_mean_TTS, "Mean"); #pca_plot_mean_TTS
pca_result_median = prcomp(data_median[,2:20002],scale.=TRUE)
pca_plot_median = plot_PCA(pca_result_median, data_median, "Median"); #pca_plot_median
pca_result_median_TSS = prcomp(data_median_TSS[,2:20002],scale.=TRUE)
pca_plot_median_TSS = plot_PCA(pca_result_median_TSS, data_median_TSS, "Median"); #pca_plot_median_TSS
pca_result_median_TTS = prcomp(data_median_TTS[,2:20002],scale.=TRUE)
pca_plot_median_TTS = plot_PCA(pca_result_median_TTS, data_median_TTS, "Median"); #pca_plot_median_TTS
pca_result_SD = prcomp(data_SD[,2:20002],scale.=TRUE)
pca_plot_SD = plot_PCA(pca_result_SD, data_SD, "SD"); #pca_plot_SD
pca_result_SD_TSS = prcomp(data_SD_TSS[,2:20002],scale.=TRUE)
pca_plot_SD_TSS = plot_PCA(pca_result_SD_TSS, data_SD_TSS, "SD"); #pca_plot_SD_TSS
pca_result_SD_TTS = prcomp(data_SD_TTS[,2:20002],scale.=TRUE)
pca_plot_SD_TTS = plot_PCA(pca_result_SD_TTS, data_SD_TTS, "SD"); #pca_plot_SD_TTS
aligned_TSS = align_plots(
  pca_plot_mean_TSS + theme(plot.margin = margin(0,0,0,0)),
  pca_plot_median_TSS + theme(plot.margin = margin(0,0,0,0)),
  pca_plot_SD_TSS + theme(plot.margin = margin(0,0,0,0)),
  align = "h", axis = "tb")
aligned_TTS = align_plots(
  pca_plot_mean_TTS + theme(plot.margin = margin(0,0,0,0)),
  pca_plot_median_TTS + theme(plot.margin = margin(0,0,0,0)),
  pca_plot_SD_TTS + theme(plot.margin = margin(0,0,0,0)),
  align = "h", axis = "tb")
aligned = align_plots(
  pca_plot_mean + theme(plot.margin = margin(0,0,0,0)),
  pca_plot_median + theme(plot.margin = margin(0,0,0,0)),
  pca_plot_SD + theme(plot.margin = margin(0,0,0,0)),
  align = "h", axis = "tb")
blank_plot <- ggplot() + theme_void()
tss_part = plot_grid(
  labels=c("", "Transcript start site (TSS)", ""),
  label_size=6, label_x=0.5, hjust=0.5, label_y=1, vjust=0.5,
  blank_plot, blank_plot, blank_plot,
  aligned_TSS[[1]],  aligned_TSS[[2]], aligned_TSS[[3]],
  ncol=3, rel_heights = c(1,11))
tts_part = plot_grid(
  labels=c("", "Transcript termination site (TTS)", ""),
  label_size=6, label_x=0.5, hjust=0.5, label_y=1, vjust=0.5,
  blank_plot, blank_plot, blank_plot,
  aligned_TTS[[1]],  aligned_TTS[[2]], aligned_TTS[[3]],
  ncol=3, rel_heights = c(1,11))
tsstts_part = plot_grid(
  labels=c("", "Transcript start and termination site (TSS + TTS)", ""),
  label_size=6, label_x=0.5, hjust=0.5, label_y=1, vjust=0.5,
  blank_plot, blank_plot, blank_plot,
  aligned[[1]],  aligned[[2]], aligned[[3]],
  ncol=3, rel_heights = c(1,11))
all = plot_grid(labels=c("", "", "", "", "a", "", "b", "", "", "", "", "", "c"), label_size=8, label_y=1, vjust=0.5, blank_plot, blank_plot, blank_plot, blank_plot, tss_part, blank_plot, tts_part, blank_plot, blank_plot, blank_plot, blank_plot, blank_plot, tsstts_part, blank_plot, legend, blank_plot, ncol=4, rel_widths=c(30,1,30,1), rel_heights=c(1,12,1,12)); all
ggsave("pca_pres_mammals_transcript.png", all, width=183, height=73.79, units="mm")
ggsave("pca_pres_mammals_transcript.pdf", all, width=183, height=73.79, units="mm")




#PCA on bird only
rm(list=ls())

load("D:/R/WD/MP_TSS/pca/transcript/pca_transcript_data.RData")


# Rename columns
colnames(data_matrix_TSS)[1:20001] <- paste0(colnames(data_matrix_TSS)[1:20001], "TSS")
colnames(data_matrix_TTS)[1:20001] <- paste0(colnames(data_matrix_TTS)[1:20001], "TTS")
# Add Statistic column if not present (Mean, Median, SD)
data_matrix_TSS$Statistic <- rep(c("Mean", "Median", "SD"), length(unique(data_matrix_TSS$Assembly)))
data_matrix_TTS$Statistic <- rep(c("Mean", "Median", "SD"), length(unique(data_matrix_TTS$Assembly)))
# Merge
data_matrix <- merge(
  data_matrix_TSS, data_matrix_TTS,
  by = c("Assembly", "Class", "Statistic"),
  all = TRUE
)
# Send columns "Class" and "Statistic" to the back
data_matrix <- data_matrix[, c(setdiff(colnames(data_matrix), c("Class", "Statistic")), "Assembly", "Statistic")]

# ↑ 한 번만

data_matrix_normalized = normalize_numbers(data_matrix)
data_matrix_TSS_normalized = normalize_numbers(data_matrix_TSS)
data_matrix_TTS_normalized = normalize_numbers(data_matrix_TTS)
#data_matrix_normalized=data_matrix # If you DON'T want to normalize
#data_matrix_TSS_normalized=data_matrix_TSS # If you DON'T want to normalize
#data_matrix_TTS_normalized=data_matrix_TTS # If you DON'T want to normalize

comprehensive <- fread("C:/Users/biopop/Documents/Downloads/MP analysis for all VGP Freeze1.0 species - Sheet1.tsv")
data_matrix_normalized <- merge(data_matrix_normalized, comprehensive, by="Assembly", all.x=TRUE)
data_matrix_TSS_normalized <- merge(data_matrix_TSS_normalized, comprehensive, by="Assembly", all.x=TRUE)
data_matrix_TTS_normalized <- merge(data_matrix_TTS_normalized, comprehensive, by="Assembly", all.x=TRUE)
# Now, in data_matrix, "Assembly" becomes column 1, columns 2-20,002 are MP, 20,011 is color.

# Remove outliers
data_matrix_normalized <- filter(data_matrix_normalized, !`Common name` %in% c("Wood mouse", "African elephant", "Yellowfin tuna", "Crucian carp", "Northern goshawk", "European flounder"))
data_matrix_TSS_normalized <- filter(data_matrix_TSS_normalized, !`Common name` %in% c("Wood mouse", "African elephant", "Yellowfin tuna", "Crucian carp", "Northern goshawk", "European flounder"))
data_matrix_TTS_normalized <- filter(data_matrix_TTS_normalized, !`Common name` %in% c("Wood mouse", "African elephant", "Yellowfin tuna", "Crucian carp", "Northern goshawk", "European flounder"))

data_mean <- data_matrix_normalized[seq(1, nrow(data_matrix_normalized), by = 3), ]
data_median <- data_matrix_normalized[seq(2, nrow(data_matrix_normalized), by = 3), ]
data_SD <- data_matrix_normalized[seq(3, nrow(data_matrix_normalized), by = 3), ]
data_mean_TSS <- data_matrix_TSS_normalized[seq(1, nrow(data_matrix_TSS_normalized), by = 3), ]
data_median_TSS <- data_matrix_TSS_normalized[seq(2, nrow(data_matrix_TSS_normalized), by = 3), ]
data_SD_TSS <- data_matrix_TSS_normalized[seq(3, nrow(data_matrix_TSS_normalized), by = 3), ]
data_mean_TTS <- data_matrix_TTS_normalized[seq(1, nrow(data_matrix_TTS_normalized), by = 3), ]
data_median_TTS <- data_matrix_TTS_normalized[seq(2, nrow(data_matrix_TTS_normalized), by = 3), ]
data_SD_TTS <- data_matrix_TTS_normalized[seq(3, nrow(data_matrix_TTS_normalized), by = 3), ]

plot_PCA <- function(pca_result,whole_data,title){
  assembly_colors = c("bTaeGut7"="#4363D8","bTaeGut1.4"="#4363DB","bAgePho1"="#4363D8","bAmmCau1"="#4363D8","bAmmNel1"="#4363D8","bMelGeo1"="#4363D8","bHaeMex1"="#4363D8","bHelExo1"="#48AAAD", "bDixPip1"="#48AAAD")
  assembly_size = c("bTaeGut7"=1, "bTaeGut1.4"=1)
  explained_var <- (pca_result$sdev^2) / sum(pca_result$sdev^2)
  percent_var <- round(100 * explained_var, 1)
  
  scores <- as.data.frame(pca_result$x)
  scores$Assembly <- whole_data$Assembly
  scores$Color <- whole_data$Color
  scores$common_name <- whole_data$"Common name"
  # Plot
  ggplot(scores, aes(x=PC1, y=PC2, label = common_name)) +
    geom_point(aes(color=Assembly, size=Assembly)) + # Change aes(color=Color) to color="transparent" to display names instead of dots
    geom_text(scores[scores$Assembly %in% c("bTaeGut7", "bTaeGut1.4"),], mapping=aes(label=Assembly), vjust = -1, size = 1.6, color="black") + # Change vjust to 0 and color="transparent" to aes(color=Color) to display names instead of dots
    # ↑ If you want to fix size instead of setting it by Size column of scores data frame, remove "size=Size" inside the aes() argument and add "size=<your size>" inside the geom_point() argument
    # scale_color_identity() + # Use the actual hex codes from the Color column
    scale_color_manual(values=assembly_colors, na.value="#AAAAAA") +
    scale_size_manual(values=assembly_size, na.value=0.4) +
    scale_shape_identity() +
    labs(x = paste0("PC1 (", percent_var[1], "%)"), y = paste0("PC2 (", percent_var[2], "%)"), title=title) +
    coord_cartesian(clip = "off") + #To prevent clipping-off of geom_texts that go out of boundaries
    #geom_text(data=subset(scores, Assembly %in% c("GRCh38.p14", "T2T-CHM13v2.0")), aes(label="★", color=Color), size=5)  +
    theme(
      legend.position="none",
      plot.title = element_text(size=5.5),
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.line.y.left = element_line(),
      axis.line.x.bottom = element_line(),
      panel.grid = element_blank(),
      axis.line=element_line(linewidth=0.1),
      axis.ticks = element_blank(),
      axis.title=element_text(size=5.5),
      axis.text=element_text(size=5, color="black")
    )->plotted
  return(plotted)
};

df = data.frame("Col"=c("Songbird", "Non-songbird passerine"))
df$Col <- factor(df$Col, levels=c("Songbird","Non-songbird passerine"))
fake <- ggplot(df, aes(x = Col, y = Col, color = Col)) +
  geom_point(size = 3) +
  scale_color_manual(name="Lineage:", values = c(
    "Songbird" = "#4363D8",
    "Non-songbird passerine"="#48AAAD"
  )) +
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 2))) +
  theme(
    #legend.position = "top",
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill="transparent", color=NA),
    legend.title = element_text(size=6, face="bold"),
    legend.text = element_text(size=6),
    legend.key.size=unit(.75,"lines"),
  ); legend = as_ggplot(get_legend(fake));legend

data_mean = data_mean[!is.na(data_mean$`Class (강)`) & data_mean$"Class (강)"=="Aves",]
data_mean_TSS = data_mean_TSS[!is.na(data_mean_TSS$`Class (강)`) & data_mean_TSS$"Class (강)"=="Aves",]
data_mean_TTS = data_mean_TTS[!is.na(data_mean_TTS$`Class (강)`) & data_mean_TTS$"Class (강)"=="Aves",]
data_median = data_median[!is.na(data_median$`Class (강)`) & data_median$"Class (강)"=="Aves",]
data_median_TSS = data_median_TSS[!is.na(data_median_TSS$`Class (강)`) & data_median_TSS$"Class (강)"=="Aves",]
data_median_TTS = data_median_TTS[!is.na(data_median_TTS$`Class (강)`) & data_median_TTS$"Class (강)"=="Aves",]
data_SD = data_SD[!is.na(data_SD$`Class (강)`) & data_SD$"Class (강)"=="Aves",]
data_SD_TSS = data_SD_TSS[!is.na(data_SD_TSS$`Class (강)`) & data_SD_TSS$"Class (강)"=="Aves",]
data_SD_TTS = data_SD_TTS[!is.na(data_SD_TTS$`Class (강)`) & data_SD_TTS$"Class (강)"=="Aves",]

pca_result_mean = prcomp(data_mean[,2:20002],scale.=TRUE)
pca_plot_mean = plot_PCA(pca_result_mean,data_mean, "Mean")
pca_result_mean_TSS = prcomp(data_mean_TSS[,2:20002],scale.=TRUE)
pca_plot_mean_TSS = plot_PCA(pca_result_mean_TSS, data_mean_TSS, "Mean"); #pca_plot_mean_TSS
pca_result_mean_TTS = prcomp(data_mean_TTS[,2:20002],scale.=TRUE)
pca_plot_mean_TTS = plot_PCA(pca_result_mean_TTS, data_mean_TTS, "Mean"); #pca_plot_mean_TTS
pca_result_median = prcomp(data_median[,2:20002],scale.=TRUE)
pca_plot_median = plot_PCA(pca_result_median, data_median, "Median"); #pca_plot_median
pca_result_median_TSS = prcomp(data_median_TSS[,2:20002],scale.=TRUE)
pca_plot_median_TSS = plot_PCA(pca_result_median_TSS, data_median_TSS, "Median"); #pca_plot_median_TSS
pca_result_median_TTS = prcomp(data_median_TTS[,2:20002],scale.=TRUE)
pca_plot_median_TTS = plot_PCA(pca_result_median_TTS, data_median_TTS, "Median"); #pca_plot_median_TTS
pca_result_SD = prcomp(data_SD[,2:20002],scale.=TRUE)
pca_plot_SD = plot_PCA(pca_result_SD, data_SD, "SD"); #pca_plot_SD
pca_result_SD_TSS = prcomp(data_SD_TSS[,2:20002],scale.=TRUE)
pca_plot_SD_TSS = plot_PCA(pca_result_SD_TSS, data_SD_TSS, "SD"); #pca_plot_SD_TSS
pca_result_SD_TTS = prcomp(data_SD_TTS[,2:20002],scale.=TRUE)
pca_plot_SD_TTS = plot_PCA(pca_result_SD_TTS, data_SD_TTS, "SD"); #pca_plot_SD_TTS
aligned_TSS = align_plots(
  pca_plot_mean_TSS + theme(plot.margin = margin(0,0,0,0)),
  pca_plot_median_TSS + theme(plot.margin = margin(0,0,0,0)),
  pca_plot_SD_TSS + theme(plot.margin = margin(0,0,0,0)),
  align = "h", axis = "tb")
aligned_TTS = align_plots(
  pca_plot_mean_TTS + theme(plot.margin = margin(0,0,0,0)),
  pca_plot_median_TTS + theme(plot.margin = margin(0,0,0,0)),
  pca_plot_SD_TTS + theme(plot.margin = margin(0,0,0,0)),
  align = "h", axis = "tb")
aligned = align_plots(
  pca_plot_mean + theme(plot.margin = margin(0,0,0,0)),
  pca_plot_median + theme(plot.margin = margin(0,0,0,0)),
  pca_plot_SD + theme(plot.margin = margin(0,0,0,0)),
  align = "h", axis = "tb")
blank_plot <- ggplot() + theme_void()
tss_part = plot_grid(
  labels=c("", "Transcript start sites (TSS)", ""),
  label_size=6, label_x=0.5, hjust=0.5, label_y=1, vjust=0.5,
  blank_plot, blank_plot, blank_plot,
  aligned_TSS[[1]],  aligned_TSS[[2]], aligned_TSS[[3]],
  ncol=3, rel_heights = c(1,11))
tts_part = plot_grid(
  labels=c("", "Transcript termination sites (TTS)", ""),
  label_size=6, label_x=0.5, hjust=0.5, label_y=1, vjust=0.5,
  blank_plot, blank_plot, blank_plot,
  aligned_TTS[[1]],  aligned_TTS[[2]], aligned_TTS[[3]],
  ncol=3, rel_heights = c(1,11))
tsstts_part = plot_grid(
  labels=c("", "Transcript start and termination sites (TSS + TTS)", ""),
  label_size=6, label_x=0.5, hjust=0.5, label_y=1, vjust=0.5,
  blank_plot, blank_plot, blank_plot,
  aligned[[1]],  aligned[[2]], aligned[[3]],
  ncol=3, rel_heights = c(1,11))
all = plot_grid(labels=c("", "", "", "", "a", "", "b", "", "", "", "", "", "c"), label_size=8, label_y=1, vjust=0.5, blank_plot, blank_plot, blank_plot, blank_plot, tss_part, blank_plot, tts_part, blank_plot, blank_plot, blank_plot, blank_plot, blank_plot, tsstts_part, blank_plot, legend, blank_plot, ncol=4, rel_widths=c(30,1,30,1), rel_heights=c(1,12,1,12)); all
ggsave("pca_pres_birds_transcript.png", all, width=183, height=73.79, units="mm")
ggsave("pca_pres_birds_transcript.pdf", all, width=183, height=73.79, units="mm")








###









### LINE ###
library(png)
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
library(svgparser)
setwd("D:/R/WD/MP_TSS/pca/transcript")

load("line_transcript_data.RData")

TSS_files <- list.files(pattern = "^MP_transcriptTSS_.*\\.tsv$")
TTS_files <- list.files(pattern = "^MP_transcriptTES_.*\\.tsv$")
line_preprocessor = function(file){
  df=fread(file)
  df[V4 == "-", V2 := V2 * -1] # If reverse (-) strand, flip the coordinates
  df = data.frame(summarise(group_by(df, V2),"Mean"=mean(V3),"Median"=median(V3),"SD"=sd(V3)))
  colnames(df)=c("Distance","Mean","Median","SD")
  df = summarise(group_by(mutate(df, x_bin = round(Distance / 100) * 100),x_bin),Mean = mean(Mean,na.rm=TRUE), Median=median(Median,na.rm=TRUE))
  df <- df %>% rename(Distance = x_bin)
  # ↑ Downsamples the data with BIN size of 100 bp. Change the value at: “x_bin = round(Distance / n) * n” to change the BIN size.
  return(df)
}


plot_line_mean<-function(df, ALPHA, SVG){
  class_colors=c("Amphibia"="#984EA3", "Aves"="#00796B", "Reptilia"="#A6D609", "Actinopterygii"="#56B4E9", "Sarcopterygii"="#A6761D", "Chondrichthyes"="#0072B2", "Mammalia"="#E69F00")
  ggplot(df, aes(x=Distance, y=Mean, group=Assembly)) +
    #annotation_raster(img, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + 
    annotation_custom(SVG, xmin = 4000, xmax = 9000, ymin = 0, ymax = 33) +
    geom_line(linewidth=.1, alpha=ALPHA, mapping=aes(color=`Class (강)`)) +
    scale_color_manual(values=class_colors, na.value="#111111") +
    # geom_line(data=subset(df, Assembly %in% c("GRCh38.p14", "T2T-CHM13v2.0", "bTaeGut7", "mPanTro3", "mPanPan1", "mGorGor1", "mPonAbe1", "mPonPyg2", "mSymSyn1", "mMacNem1", "mNycCou1", "mCynVol1", "mApoSyl1", "mChiNiv1", "mMarFla1", "mOchPri1", "mEriEur2", "mSorAra1", "bAgePho1", "bAmmCau1", "bAmmNel1", "bMelGeo1", "bHaeMex1")), aes(color=Color), linewidth=.1, alpha=1) +
    #geom_line(data=subset(df, Assembly %in% c("GRCh38.p14", "T2T-CHM13v2.0", "bTaeGut7", "bTaeGut1.4")), color="black", linewidth=.15, alpha=1) +
    scale_y_continuous(limits=c(0,100), expand=c(0,0)) +
    scale_x_continuous(labels = label_comma(), breaks = c(-10000, 0, 10000), limits = c(-10000, 10000), expand=c(0,0)) +
    #geom_text(data=label_data, aes(label=`Common name`), hjust=0, nudge_x=-5, size=3) +
    labs(y="Average MP (%)", x="Distance (bp)") + 
    theme(legend.position="none",
          panel.background = element_blank(),
          panel.grid = element_blank(),
          plot.background = element_blank(),
          axis.line.y.left = element_line(linewidth=0.1),
          axis.line.x.bottom = element_line(linewidth=0.1),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank())
}

plot_line_median<-function(df, ALPHA, SVG){
  class_colors=c("Amphibia"="#984EA3", "Aves"="#00796B", "Reptilia"="#A6D609", "Actinopterygii"="#56B4E9", "Sarcopterygii"="#A6761D", "Chondrichthyes"="#0072B2", "Mammalia"="#E69F00")
  ggplot(df, aes(x=Distance, y=Median, group=Assembly)) +
    #annotation_raster(img, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + 
    annotation_custom(SVG, xmin = 4000, xmax = 9000, ymin = 0, ymax = 33) +
    geom_line(linewidth=.1, alpha=ALPHA, mapping=aes(color=`Class (강)`)) +
    scale_color_manual(values=class_colors, na.value="#111111") +
    # geom_line(data=subset(df, Assembly %in% c("GRCh38.p14", "T2T-CHM13v2.0", "bTaeGut7", "mPanTro3", "mPanPan1", "mGorGor1", "mPonAbe1", "mPonPyg2", "mSymSyn1", "mMacNem1", "mNycCou1", "mCynVol1", "mApoSyl1", "mChiNiv1", "mMarFla1", "mOchPri1", "mEriEur2", "mSorAra1", "bAgePho1", "bAmmCau1", "bAmmNel1", "bMelGeo1", "bHaeMex1")), aes(color=Color), linewidth=.1, alpha=1) +
    #geom_line(data=subset(df, Assembly %in% c("GRCh38.p14", "T2T-CHM13v2.0", "bTaeGut7", "bTaeGut1.4")), color="black", linewidth=.15, alpha=1) +
    scale_y_continuous(limits=c(0,100), expand=c(0,0)) +
    scale_x_continuous(labels = label_comma(), breaks = c(-10000, 0, 10000), limits = c(-10000, 10000), expand=c(0,0)) +
    #geom_text(data=label_data, aes(label=`Common name`), hjust=0, nudge_x=-5, size=3) +
    labs(y="Average MP (%)", x="Distance (bp)") + 
    theme(legend.position="none",
          panel.background = element_blank(),
          panel.grid = element_blank(),
          plot.background = element_blank(),
          axis.line.y.left = element_line(linewidth=0.1),
          axis.line.x.bottom = element_line(linewidth=0.1),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank())
}




combined_data_TSS <- bind_rows(lapply(TSS_files, function(x) {
  df <- line_preprocessor(x)
  df$Assembly <- gsub("^MP_transcriptTSS_|\\.tsv$", "", basename(x))
  df$Class <- substr(gsub("^MP_transcriptTSS_|\\.tsv$", "", basename(x)), 1, 1)
  df
})) # 이 부분에서 종마다 class/order를 assign해서 df에 넣어줌. ⌛

combined_data_TTS <- bind_rows(lapply(TTS_files, function(x) {
  df <- line_preprocessor(x)
  df$Assembly <- gsub("^MP_transcriptTES_|\\.tsv$", "", basename(x))
  df$Class <- substr(gsub("^MP_transcriptTES_|\\.tsv$", "", basename(x)), 1, 1)
  df
})) # 이 부분에서 종마다 class/order를 assign해서 df에 넣어줌. ⌛

#save.image("D:/R/WD/MP_TSS/pca/transcript/line_transcript_data.RData")

# Add EGAP
create_line_data=function(x) {
  df <- line_preprocessor(x)
  df$Assembly <- gsub("^MP_transcriptTSS_|\\.tsv$", "", basename(x))
  df$Class <- substr(gsub("^MP_transcriptTSS_|\\.tsv$", "", basename(x)), 1, 1)
  df
}
combined_data_TSS=bind_rows(combined_data_TSS, create_line_data("MP_transcriptTSS_bTaeGut7v0.4_MT_rDNA.EGAPx.v0.1_mat.tsv"))
create_line_data=function(x) {
  df <- line_preprocessor(x)
  df$Assembly <- gsub("^MP_transcriptTES_|\\.tsv$", "", basename(x))
  df$Class <- substr(gsub("^MP_transcriptTES_|\\.tsv$", "", basename(x)), 1, 1)
  df
}
combined_data_TTS=bind_rows(combined_data_TTS, create_line_data("MP_transcriptTES_bTaeGut7v0.4_MT_rDNA.EGAPx.v0.1_mat.tsv"))

normalize_by_assembly <- function(df) {
  # Ensure required columns exist
  if (!all(c("Mean", "Median", "Assembly") %in% colnames(df))) {
    stop("The input df must contain 'Mean', 'Median', and 'Assembly' columns.")
  }
  
  # Function for min-max normalization to [0, 100] within a group
  min_max_100 <- function(x) {
    if (max(x) == min(x)) {
      rep(50, length(x))  # All identical: assign 50 (midpoint of 0 and 100)
    } else {
      (x - min(x)) / (max(x) - min(x)) * 100
    }
  }
  # Apply normalization by Assembly group
  df$Mean_norm <- ave(df$Mean, df$Assembly, FUN = min_max_100)
  df$Median_norm <- ave(df$Median, df$Assembly, FUN = min_max_100)
  return(df)
}

normalize_by_assembly(combined_data_TSS)->combined_data_TSS_normalized # Skip this function if you DON'T want to normalize ea. assemblies
normalize_by_assembly(combined_data_TTS)->combined_data_TTS_normalized # Skip this function if you DON'T want to normalize ea. assemblies

comprehensive <- fread("C:/Users/biopop/Documents/Downloads/MP analysis for all VGP Freeze1.0 species - Sheet1.tsv")
combined_data_TSS_normalized <- merge(combined_data_TSS_normalized, comprehensive, by="Assembly", all.x=TRUE)
combined_data_TTS_normalized <- merge(combined_data_TTS_normalized, comprehensive, by="Assembly", all.x=TRUE)
combined_data_TSS_filtered <- filter(combined_data_TSS_normalized, !`Common name` %in% c("Wood mouse", "African elephant", "Yellowfin tuna", "Crucian carp", "Northern goshawk", "European flounder") & `Assembly` != "bTaeGut7")
combined_data_TTS_filtered <- filter(combined_data_TTS_normalized, !`Common name` %in% c("Wood mouse", "African elephant", "Yellowfin tuna", "Crucian carp", "Northern goshawk", "European flounder") & `Assembly` != "bTaeGut7")

label_data = combined_data_TSS_filtered %>% group_by(Species) %>% filter(Distance == max(Distance)) %>% ungroup()
label_data_m = label_data[label_data$Class=="m",]
label_data_b = label_data[label_data$Class=="b",]
label_data_r = label_data[label_data$Class=="r",]
label_data_a = label_data[label_data$Class=="a",]
label_data_f = label_data[label_data$Class=="f",]
label_data_s = label_data[label_data$Class=="s",]
label_data_fs = label_data[label_data$Class %in% c("f", "s"), ]

label_data_TTS = combined_data_TTS_filtered %>% group_by(Species) %>% filter(Distance == max(Distance)) %>% ungroup()
label_data_m_TTS = label_data_TTS[label_data_TTS$Class=="m",]
label_data_b_TTS = label_data_TTS[label_data_TTS$Class=="b",]
label_data_r_TTS = label_data_TTS[label_data_TTS$Class=="r",]
label_data_a_TTS = label_data_TTS[label_data_TTS$Class=="a",]
label_data_f_TTS = label_data_TTS[label_data_TTS$Class=="f",]
label_data_s_TTS = label_data_TTS[label_data_TTS$Class=="s",]
label_data_fs_TTS = label_data_TTS[label_data_TTS$Class %in% c("f", "s"), ]



ggplot()+
  scale_x_continuous(labels = c("-10","0","10"), breaks = c(-10000, 0, 10000), limits = c(-10000, 10000), name="Distance from TSS (kbp)", expand=c(0,0)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(size=5, color="black"),
    axis.title.x = element_blank())->Xaxis_line;Xaxis_line

ggplot()+
  scale_y_continuous(labels = c("0","25","50","75","100"), breaks = c(0,25,50,75,100), limits = c(0,100), name="Median MP (%)", expand=c(0,0)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size=5, color="black"),
    axis.title.y = element_blank()
  )->Yaxis_line;Yaxis_line

blank <- ggplot() + theme_void()

empty = svgparser::read_svg("https://upload.wikimedia.org/wikipedia/commons/1/1d/No_image.svg")
mammal = svgparser::read_svg("svgs/01_mammals_silhouette_01192026.svg")
bird=svgparser::read_svg("svgs/02_birds_silhouette_01192026.svg")
reptile=svgparser::read_svg("svgs/03_reptiles_silhouette_01192026.svg")
amphibian=svgparser::read_svg("svgs/04_amphibians_silhouette_01192026.svg")
coel=svgparser::read_svg("svgs/05_coelacanth_silhouette_01192026.svg")
fish=svgparser::read_svg("svgs/06_ray-finned_Fishes_silhouette_01192026.svg")
shark=svgparser::read_svg("svgs/07_Cartilaginous_fish_silhouette_01192026.svg")


# TSS median
plot_line_median(combined_data_TSS_filtered, 0.2, empty) -> all_line
plot_line_median(combined_data_TSS_filtered[combined_data_TSS_filtered$Class %in% c("m", "h", "G", "T"),], 1, empty) -> mammal_line
plot_line_median(combined_data_TSS_filtered[combined_data_TSS_filtered$Class=="b",], 1,empty) -> bird_line
plot_line_median(combined_data_TSS_filtered[combined_data_TSS_filtered$Class=="r",], 1,empty) -> reptile_line
plot_line_median(combined_data_TSS_filtered[combined_data_TSS_filtered$Class=="a",], 1,empty) -> amphibian_line
plot_line_median(combined_data_TSS_filtered[combined_data_TSS_filtered$Assembly=="fLatCha1",], 1,empty) -> coel_line
plot_line_median(combined_data_TSS_filtered[combined_data_TSS_filtered$Class=="f"&combined_data_TSS_filtered$Assembly!="fLatCha1",], 1,empty) -> fish_line
plot_line_median(combined_data_TSS_filtered[combined_data_TSS_filtered$Class=="s",], 1,empty) -> shark_line
median_TSS_row = align_plots(Yaxis_line + theme(plot.margin = margin(8,0,0,0)), all_line + theme(plot.margin = margin(8,8,0,0)), mammal_line + theme(plot.margin = margin(8,8,0,0)), bird_line + theme(plot.margin = margin(8,8,0,0)), reptile_line + theme(plot.margin = margin(8,8,0,0)), amphibian_line + theme(plot.margin = margin(8,8,0,0)), coel_line + theme(plot.margin = margin(8,8,0,0)), fish_line + theme(plot.margin = margin(8,8,0,0)), shark_line + theme(plot.margin = margin(8,8,0,0)), align="h", axis="tb")

# TSS mean
plot_line_mean(combined_data_TSS_filtered, 0.2,empty) -> all_line
plot_line_mean(combined_data_TSS_filtered[combined_data_TSS_filtered$Class %in% c("m", "h", "G", "T"),], 1,mammal) -> mammal_line
plot_line_mean(combined_data_TSS_filtered[combined_data_TSS_filtered$Class=="b",], 1,bird) -> bird_line
plot_line_mean(combined_data_TSS_filtered[combined_data_TSS_filtered$Class=="r",], 1,reptile) -> reptile_line
plot_line_mean(combined_data_TSS_filtered[combined_data_TSS_filtered$Class=="a",], 1,amphibian) -> amphibian_line
plot_line_mean(combined_data_TSS_filtered[combined_data_TSS_filtered$Assembly=="fLatCha1",], 1,coel) -> coel_line
plot_line_mean(combined_data_TSS_filtered[combined_data_TSS_filtered$Class=="f"&combined_data_TSS_filtered$Assembly!="fLatCha1",], 1,fish) -> fish_line
plot_line_mean(combined_data_TSS_filtered[combined_data_TSS_filtered$Class=="s",], 1,shark) -> shark_line
mean_TSS_row = align_plots(Yaxis_line + theme(plot.margin = margin(8,0,0,0)), all_line + theme(plot.margin = margin(8,8,0,0)), mammal_line + theme(plot.margin = margin(8,8,0,0)), bird_line + theme(plot.margin = margin(8,8,0,0)), reptile_line + theme(plot.margin = margin(8,8,0,0)), amphibian_line + theme(plot.margin = margin(8,8,0,0)), coel_line + theme(plot.margin = margin(8,8,0,0)), fish_line + theme(plot.margin = margin(8,8,0,0)), shark_line + theme(plot.margin = margin(8,8,0,0)), align="h", axis="tb")

# Axis
axis_row = align_plots(mammal_line + theme(plot.margin=margin(8,8,0,0)), Xaxis_line + theme(plot.margin=margin(0,8,0,0)), align="v", axis="lr")
axis_blank = align_plots(Yaxis_line + theme(plot.margin = margin(8,0,0,0)), blank + theme(plot.margin = margin(0,0,0,0)), align="v", axis="lr")

TSS=plot_grid(mean_TSS_row[[1]], mean_TSS_row[[2]], mean_TSS_row[[3]], mean_TSS_row[[4]], mean_TSS_row[[5]], mean_TSS_row[[6]], mean_TSS_row[[7]], mean_TSS_row[[8]], mean_TSS_row[[9]],
              median_TSS_row[[1]], median_TSS_row[[2]], median_TSS_row[[3]], median_TSS_row[[4]], median_TSS_row[[5]], median_TSS_row[[6]], median_TSS_row[[7]], median_TSS_row[[8]], median_TSS_row[[9]],
              axis_blank[[2]], axis_row[[2]], axis_row[[2]], axis_row[[2]], axis_row[[2]], axis_row[[2]], axis_row[[2]], axis_row[[2]], axis_row[[2]],
              ncol=9, rel_heights=c(10,10,0.5), rel_widths=c(1,10,10,10,10,10,10,10,10))
rotated_label_median <- ggdraw() + draw_grob(textGrob("Median MP (%)", rot = 90, gp = gpar(fontsize = 5.5, fontface = "plain")))
rotated_label_mean <- ggdraw() + draw_grob(textGrob("Mean MP (%)", rot = 90, gp = gpar(fontsize = 5.5, fontface = "plain")))
rotated_label = plot_grid(rotated_label_mean, rotated_label_median, blank, ncol=1, rel_heights=c(10,10,0.5))
TSS_quarterTitled = plot_grid(rotated_label, TSS, ncol=2, rel_widths=c(3,81))
TSS_halfTitled = plot_grid(TSS_quarterTitled, blank, ncol=1, rel_heights=c(20.5,3), labels=c("", "Distance from TSS (kbp)"), label_size=5.5, label_x=0.5, hjust=0.5, label_y=0.5, vjust=0.5, label_fontface="plain")
TSS_titled = plot_grid(blank, TSS_halfTitled, ncol=1, rel_heights=c(1,23.5), labels=c("Transcript start site (TSS) in all classes"), label_size=7, label_x=0.5, hjust=0.5, label_y=0.5, vjust=0.5)



# TTS median
plot_line_median(combined_data_TTS_filtered, 0.2,empty) -> all_line
plot_line_median(combined_data_TTS_filtered[combined_data_TTS_filtered$Class %in% c("m", "h", "G", "T"),], 1,empty) -> mammal_line
plot_line_median(combined_data_TTS_filtered[combined_data_TTS_filtered$Class=="b",], 1,empty) -> bird_line
plot_line_median(combined_data_TTS_filtered[combined_data_TTS_filtered$Class=="r",], 1,empty) -> reptile_line
plot_line_median(combined_data_TTS_filtered[combined_data_TTS_filtered$Class=="a",], 1,empty) -> amphibian_line
plot_line_median(combined_data_TTS_filtered[combined_data_TTS_filtered$Assembly=="fLatCha1",], 1,empty) -> coel_line
plot_line_median(combined_data_TTS_filtered[combined_data_TTS_filtered$Class=="f"&combined_data_TSS_filtered$Assembly!="fLatCha1",], 1,empty) -> fish_line
plot_line_median(combined_data_TTS_filtered[combined_data_TTS_filtered$Class=="s",], 1,empty) -> shark_line
median_TTS_row = align_plots(Yaxis_line + theme(plot.margin = margin(8,0,0,0)), all_line + theme(plot.margin = margin(8,8,0,0)), mammal_line + theme(plot.margin = margin(8,8,0,0)), bird_line + theme(plot.margin = margin(8,8,0,0)), reptile_line + theme(plot.margin = margin(8,8,0,0)), amphibian_line + theme(plot.margin = margin(8,8,0,0)), coel_line + theme(plot.margin = margin(8,8,0,0)), fish_line + theme(plot.margin = margin(8,8,0,0)), shark_line + theme(plot.margin = margin(8,8,0,0)), align="h", axis="tb")

# TTS mean
plot_line_mean(combined_data_TTS_filtered, 0.2,empty) -> all_line
plot_line_mean(combined_data_TTS_filtered[combined_data_TTS_filtered$Class %in% c("m", "h", "G", "T"),], 1,empty) -> mammal_line
plot_line_mean(combined_data_TTS_filtered[combined_data_TTS_filtered$Class=="b",], 1,empty) -> bird_line
plot_line_mean(combined_data_TTS_filtered[combined_data_TTS_filtered$Class=="r",], 1,empty) -> reptile_line
plot_line_mean(combined_data_TTS_filtered[combined_data_TTS_filtered$Class=="a",], 1,empty) -> amphibian_line
plot_line_mean(combined_data_TTS_filtered[combined_data_TTS_filtered$Assembly=="fLatCha1",], 1,empty) -> coel_line
plot_line_mean(combined_data_TTS_filtered[combined_data_TTS_filtered$Class=="f"&combined_data_TSS_filtered$Assembly!="fLatCha1",], 1,empty) -> fish_line
plot_line_mean(combined_data_TTS_filtered[combined_data_TTS_filtered$Class=="s",], 1,empty) -> shark_line
mean_TTS_row = align_plots(Yaxis_line + theme(plot.margin = margin(8,0,0,0)), all_line + theme(plot.margin = margin(8,8,0,0)), mammal_line + theme(plot.margin = margin(8,8,0,0)), bird_line + theme(plot.margin = margin(8,8,0,0)), reptile_line + theme(plot.margin = margin(8,8,0,0)), amphibian_line + theme(plot.margin = margin(8,8,0,0)), coel_line + theme(plot.margin = margin(8,8,0,0)), fish_line + theme(plot.margin = margin(8,8,0,0)), shark_line + theme(plot.margin = margin(8,8,0,0)), align="h", axis="tb")

# Axis
axis_row = align_plots(mammal_line + theme(plot.margin=margin(8,8,0,0)), Xaxis_line + theme(plot.margin=margin(0,8,0,0)), align="v", axis="lr")
axis_blank = align_plots(Yaxis_line + theme(plot.margin = margin(8,0,0,0)), blank + theme(plot.margin = margin(0,0,0,0)), align="v", axis="lr")

TTS=plot_grid(mean_TTS_row[[1]], mean_TTS_row[[2]], mean_TTS_row[[3]], mean_TTS_row[[4]], mean_TTS_row[[5]], mean_TTS_row[[6]], mean_TTS_row[[7]], mean_TTS_row[[8]], mean_TTS_row[[9]],
              median_TTS_row[[1]], median_TTS_row[[2]], median_TTS_row[[3]], median_TTS_row[[4]], median_TTS_row[[5]], median_TTS_row[[6]], median_TTS_row[[7]], median_TTS_row[[8]], median_TTS_row[[9]],
              axis_blank[[2]], axis_row[[2]], axis_row[[2]], axis_row[[2]], axis_row[[2]], axis_row[[2]], axis_row[[2]], axis_row[[2]], axis_row[[2]],
              ncol=9, rel_heights=c(10,10,0.5), rel_widths=c(1,10,10,10,10,10,10,10,10))
rotated_label_median <- ggdraw() + draw_grob(textGrob("Median MP (%)", rot = 90, gp = gpar(fontsize = 5.5, fontface = "plain")))
rotated_label_mean <- ggdraw() + draw_grob(textGrob("Mean MP (%)", rot = 90, gp = gpar(fontsize = 5.5, fontface = "plain")))
rotated_label = plot_grid(rotated_label_mean, rotated_label_median, blank, ncol=1, rel_heights=c(10,10,0.5))
TTS_quarterTitled = plot_grid(rotated_label, TTS, ncol=2, rel_widths=c(3,81))
TTS_halfTitled = plot_grid(TTS_quarterTitled, blank, ncol=1, rel_heights=c(20.5,3), labels=c("", "Distance from TTS (kbp)"), label_size=5.5, label_x=0.5, hjust=0.5, label_y=0.5, vjust=0.5, label_fontface="plain")
TTS_titled = plot_grid(blank, TTS_halfTitled, ncol=1, rel_heights=c(1,23.5), labels=c("Transcript termination site (TTS) in all classes"), label_size=7, label_x=0.5, hjust=0.5, label_y=0.5, vjust=0.5)


complete = plot_grid(labels=c("a", "", "b"), label_size=8, TSS_titled, blank, TTS_titled, ncol=1, rel_heights=c(24.5,1,24.5))

widths = 3+81
heights=24.5 + 1 + 24.5

final_width = 183
final_height = 183 / widths * heights
ggsave("line_transcript.png", complete, units="mm", height=final_height, width=final_width)
ggsave("line_transcript.pdf", complete, units="mm", height=final_height, width=final_width)







### Both (Supplementary Fig. 4) ###

# Mammal - TSS

plot_line_mean<-function(df, ALPHA, SVG){
  ggplot(df, aes(x=Distance, y=Mean, group=Assembly)) +
    #annotation_raster(img, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + 
    annotation_custom(SVG, xmin = 4000, xmax = 9000, ymin = 0, ymax = 33) +
    geom_line(data=df[!df$Assembly %in% c("GRCh38.p14", "T2T-CHM13v2.0", "mGorGor1", "mMacNem1", "mNycCou1", "mPanPan1", "mPanTro3", "mPonAbe1", "mPonPyg2", "mSymSyn1"),],linewidth=.1, alpha=ALPHA, color="#F0E442") +
    geom_line(data=df[df$Assembly %in% c("GRCh38.p14", "T2T-CHM13v2.0", "mGorGor1", "mMacNem1", "mNycCou1", "mPanPan1", "mPanTro3", "mPonAbe1", "mPonPyg2", "mSymSyn1"),],linewidth=.1, alpha=ALPHA, color="#D55E00") +
    # geom_line(data=subset(df, Assembly %in% c("GRCh38.p14", "T2T-CHM13v2.0", "bTaeGut7", "mPanTro3", "mPanPan1", "mGorGor1", "mPonAbe1", "mPonPyg2", "mSymSyn1", "mMacNem1", "mNycCou1", "mCynVol1", "mApoSyl1", "mChiNiv1", "mMarFla1", "mOchPri1", "mEriEur2", "mSorAra1", "bAgePho1", "bAmmCau1", "bAmmNel1", "bMelGeo1", "bHaeMex1")), aes(color=Color), linewidth=.1, alpha=1) +
    #geom_line(data=subset(df, Assembly %in% c("GRCh38.p14", "T2T-CHM13v2.0", "bTaeGut7", "bTaeGut1.4")), color="black", linewidth=.15, alpha=1) +
    scale_y_continuous(limits=c(0,100), expand=c(0,0)) +
    scale_x_continuous(labels = label_comma(), breaks = c(-10000, 0, 10000), limits = c(-10000, 10000), expand=c(0,0)) +
    #geom_text(data=label_data, aes(label=`Common name`), hjust=0, nudge_x=-5, size=3) +
    labs(y="Average MP (%)", x="Distance (bp)") + 
    theme(legend.position="none",
          panel.background = element_blank(),
          panel.grid = element_blank(),
          plot.background = element_blank(),
          axis.line.y.left = element_line(linewidth=0.1),
          axis.line.x.bottom = element_line(linewidth=0.1),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank())
}

plot_line_median<-function(df, ALPHA, SVG){
  ggplot(df, aes(x=Distance, y=Median, group=Assembly)) +
    #annotation_raster(img, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + 
    annotation_custom(SVG, xmin = 4000, xmax = 9000, ymin = 0, ymax = 33) +
    geom_line(data=df[!df$Assembly %in% c("GRCh38.p14", "T2T-CHM13v2.0", "mGorGor1", "mMacNem1", "mNycCou1", "mPanPan1", "mPanTro3", "mPonAbe1", "mPonPyg2", "mSymSyn1"),],linewidth=.1, alpha=ALPHA, color="#F0E442") +
    geom_line(data=df[df$Assembly %in% c("GRCh38.p14", "T2T-CHM13v2.0", "mGorGor1", "mMacNem1", "mNycCou1", "mPanPan1", "mPanTro3", "mPonAbe1", "mPonPyg2", "mSymSyn1"),],linewidth=.1, alpha=ALPHA, color="#D55E00") +
    # geom_line(data=subset(df, Assembly %in% c("GRCh38.p14", "T2T-CHM13v2.0", "bTaeGut7", "mPanTro3", "mPanPan1", "mGorGor1", "mPonAbe1", "mPonPyg2", "mSymSyn1", "mMacNem1", "mNycCou1", "mCynVol1", "mApoSyl1", "mChiNiv1", "mMarFla1", "mOchPri1", "mEriEur2", "mSorAra1", "bAgePho1", "bAmmCau1", "bAmmNel1", "bMelGeo1", "bHaeMex1")), aes(color=Color), linewidth=.1, alpha=1) +
    #geom_line(data=subset(df, Assembly %in% c("GRCh38.p14", "T2T-CHM13v2.0", "bTaeGut7", "bTaeGut1.4")), color="black", linewidth=.15, alpha=1) +
    scale_y_continuous(limits=c(0,100), expand=c(0,0)) +
    scale_x_continuous(labels = label_comma(), breaks = c(-10000, 0, 10000), limits = c(-10000, 10000), expand=c(0,0)) +
    #geom_text(data=label_data, aes(label=`Common name`), hjust=0, nudge_x=-5, size=3) +
    labs(y="Average MP (%)", x="Distance (bp)") + 
    theme(legend.position="none",
          panel.background = element_blank(),
          panel.grid = element_blank(),
          plot.background = element_blank(),
          axis.line.y.left = element_line(linewidth=0.1),
          axis.line.x.bottom = element_line(linewidth=0.1),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank())
}

# Create legend
class_order <- c("Primate", "Non-primate mammal")
df <- data.frame(X = 1,Y = 1,Class = c("Primate", "Non-primate mammal"))
df$Class <- factor(df$Class, levels = class_order)
fake <- ggplot(df, aes(x = X, y = Y, color = Class)) +
  geom_point(size = 3) +
  scale_color_manual(name="Order:", values = c(
    "Primate"="#D55E00",
    "Non-primate mammal"="#F0E442"
  )) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 2))) +
  theme(
    #legend.position = "top",
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill="transparent", color=NA),
    legend.title = element_text(size=6, face="bold"),
    legend.text = element_text(size=6),
    legend.key.size=unit(.75,"lines")
  )
legend = as_ggplot(get_legend(fake));legend

plot_PCA <- function(pca_result,whole_data,title){
  assembly_colors = c("GRCh38.p14"="#D55E00", "T2T-CHM13v2.0"="#D55E00", "mGorGor1"="#D55E00", "mMacNem1"="#D55E00", "mNycCou1"="#D55E00", "mPanPan1"="#D55E00", "mPanTro3"="#D55E00", "mPonAbe1"="#D55E00", "mPonPyg2"="#D55E00", "mSymSyn1"="#D55E00")
  explained_var <- (pca_result$sdev^2) / sum(pca_result$sdev^2)
  percent_var <- round(100 * explained_var, 1)
  
  scores <- as.data.frame(pca_result$x)
  scores$Assembly <- whole_data$Assembly
  scores$Color <- whole_data$Color
  scores$common_name <- whole_data$"Common name"
  # Plot
  ggplot(scores, aes(x=PC1, y=PC2, label = common_name)) +
    geom_text(vjust = 0, size = 2, color="transparent") + # Change vjust to 0 and color="transparent" to aes(color=Color) to display names instead of dots
    geom_point(aes(color=Assembly), size=0.3) + # Change aes(color=Color) to color="transparent" to display names instead of dots
    # ↑ If you want to fix size instead of setting it by Size column of scores data frame, remove "size=Size" inside the aes() argument and add "size=<your size>" inside the geom_point() argument
    # scale_color_identity() + # Use the actual hex codes from the Color column
    scale_color_manual(values=assembly_colors, na.value="#F0E442") +
    scale_shape_identity() +
    labs(x = paste0("PC1 (", percent_var[1], "%)"), y = paste0("PC2 (", percent_var[2], "%)"), title=title) +
    coord_cartesian(clip = "off") + #To prevent clipping-off of geom_texts that go out of boundaries
    #geom_text(data=subset(scores, Assembly %in% c("GRCh38.p14", "T2T-CHM13v2.0")), aes(label="★", color=Color), size=5)  +
    theme(
      legend.position="none",
      plot.title = element_text(size=5.5),
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.line.y.left = element_line(),
      axis.line.x.bottom = element_line(),
      panel.grid = element_blank(),
      axis.line=element_line(linewidth=0.1),
      axis.ticks = element_blank(),
      axis.title=element_text(size=5.5),
      axis.text=element_blank(),
      plot.margin = unit(c(8,8,0,0), "pt")
    )->plotted
  return(plotted)
}


# Rename columns
colnames(data_matrix_TSS)[1:20001] <- paste0(colnames(data_matrix_TSS)[1:20001], "TSS")
colnames(data_matrix_TTS)[1:20001] <- paste0(colnames(data_matrix_TTS)[1:20001], "TTS")
# Add Statistic column if not present (Mean, Median, SD)
data_matrix_TSS$Statistic <- rep(c("Mean", "Median", "SD"), length(unique(data_matrix_TSS$Assembly)))
data_matrix_TTS$Statistic <- rep(c("Mean", "Median", "SD"), length(unique(data_matrix_TTS$Assembly)))
# Merge
data_matrix <- merge(
  data_matrix_TSS, data_matrix_TTS,
  by = c("Assembly", "Class", "Statistic"),
  all = TRUE
)
# Send columns "Class" and "Statistic" to the back
data_matrix <- data_matrix[, c(setdiff(colnames(data_matrix), c("Class", "Statistic")), "Assembly", "Statistic")]


data_matrix_normalized = normalize_numbers(data_matrix)
data_matrix_TSS_normalized = normalize_numbers(data_matrix_TSS)
data_matrix_TTS_normalized = normalize_numbers(data_matrix_TTS)
#data_matrix_normalized=data_matrix # If you DON'T want to normalize
#data_matrix_TSS_normalized=data_matrix_TSS # If you DON'T want to normalize
#data_matrix_TTS_normalized=data_matrix_TTS # If you DON'T want to normalize

comprehensive <- fread("C:/Users/biopop/Documents/Downloads/MP analysis for all VGP Freeze1.0 species - Sheet1.tsv")
data_matrix_normalized <- merge(data_matrix_normalized, comprehensive, by="Assembly", all.x=TRUE)
data_matrix_TSS_normalized <- merge(data_matrix_TSS_normalized, comprehensive, by="Assembly", all.x=TRUE)
data_matrix_TTS_normalized <- merge(data_matrix_TTS_normalized, comprehensive, by="Assembly", all.x=TRUE)
# Now, in data_matrix, "Assembly" becomes column 1, columns 2-20,002 are MP, 20,011 is color.

# Remove outliers
data_matrix_normalized <- filter(data_matrix_normalized, !`Common name` %in% c("Wood mouse", "African elephant", "Yellowfin tuna", "Crucian carp", "Northern goshawk", "European flounder") & `Assembly` != "bTaeGut7")
data_matrix_TSS_normalized <- filter(data_matrix_TSS_normalized, !`Common name` %in% c("Wood mouse", "African elephant", "Yellowfin tuna", "Crucian carp", "Northern goshawk", "European flounder") & `Assembly` != "bTaeGut7")
data_matrix_TTS_normalized <- filter(data_matrix_TTS_normalized, !`Common name` %in% c("Wood mouse", "African elephant", "Yellowfin tuna", "Crucian carp", "Northern goshawk", "European flounder") & `Assembly` != "bTaeGut7")
data_matrix_normalized$`Class (강)`[data_matrix_normalized$Species=="Homo sapiens"]="Human"
data_matrix_TSS_normalized$`Class (강)`[data_matrix_TSS_normalized$Species=="Homo sapiens"]="Human"
data_matrix_TTS_normalized$`Class (강)`[data_matrix_TTS_normalized$Species=="Homo sapiens"]="Human"

data_mean <- data_matrix_normalized[seq(1, nrow(data_matrix_normalized), by = 3), ]
data_median <- data_matrix_normalized[seq(2, nrow(data_matrix_normalized), by = 3), ]
data_SD <- data_matrix_normalized[seq(3, nrow(data_matrix_normalized), by = 3), ]
data_mean_TSS <- data_matrix_TSS_normalized[seq(1, nrow(data_matrix_TSS_normalized), by = 3), ]
data_median_TSS <- data_matrix_TSS_normalized[seq(2, nrow(data_matrix_TSS_normalized), by = 3), ]
data_SD_TSS <- data_matrix_TSS_normalized[seq(3, nrow(data_matrix_TSS_normalized), by = 3), ]
data_mean_TTS <- data_matrix_TTS_normalized[seq(1, nrow(data_matrix_TTS_normalized), by = 3), ]
data_median_TTS <- data_matrix_TTS_normalized[seq(2, nrow(data_matrix_TTS_normalized), by = 3), ]
data_SD_TTS <- data_matrix_TTS_normalized[seq(3, nrow(data_matrix_TTS_normalized), by = 3), ]


data_mean = data_mean[!is.na(data_mean$`Class (강)`) & data_mean$"Class (강)"=="Mammalia",]
data_mean_TSS = data_mean_TSS[!is.na(data_mean_TSS$`Class (강)`) & data_mean_TSS$"Class (강)"=="Mammalia",]
data_mean_TTS = data_mean_TTS[!is.na(data_mean_TTS$`Class (강)`) & data_mean_TTS$"Class (강)"=="Mammalia",]
data_median = data_median[!is.na(data_median$`Class (강)`) & data_median$"Class (강)"=="Mammalia",]
data_median_TSS = data_median_TSS[!is.na(data_median_TSS$`Class (강)`) & data_median_TSS$"Class (강)"=="Mammalia",]
data_median_TTS = data_median_TTS[!is.na(data_median_TTS$`Class (강)`) & data_median_TTS$"Class (강)"=="Mammalia",]
data_SD = data_SD[!is.na(data_SD$`Class (강)`) & data_SD$"Class (강)"=="Mammalia",]
data_SD_TSS = data_SD_TSS[!is.na(data_SD_TSS$`Class (강)`) & data_SD_TSS$"Class (강)"=="Mammalia",]
data_SD_TTS = data_SD_TTS[!is.na(data_SD_TTS$`Class (강)`) & data_SD_TTS$"Class (강)"=="Mammalia",]

pca_result_mean = prcomp(data_mean[,2:20002],scale.=TRUE)
pca_plot_mean = plot_PCA(pca_result_mean,data_mean, "Mean") #pca_plot_mean
pca_result_mean_TSS = prcomp(data_mean_TSS[,2:20002],scale.=TRUE)
pca_plot_mean_TSS = plot_PCA(pca_result_mean_TSS, data_mean_TSS, "Mean"); #pca_plot_mean_TSS
pca_result_mean_TTS = prcomp(data_mean_TTS[,2:20002],scale.=TRUE)
pca_plot_mean_TTS = plot_PCA(pca_result_mean_TTS, data_mean_TTS, "Mean"); #pca_plot_mean_TTS
pca_result_median = prcomp(data_median[,2:20002],scale.=TRUE)
pca_plot_median = plot_PCA(pca_result_median, data_median, "Median"); #pca_plot_median
pca_result_median_TSS = prcomp(data_median_TSS[,2:20002],scale.=TRUE)
pca_plot_median_TSS = plot_PCA(pca_result_median_TSS, data_median_TSS, "Median"); #pca_plot_median_TSS
pca_result_median_TTS = prcomp(data_median_TTS[,2:20002],scale.=TRUE)
pca_plot_median_TTS = plot_PCA(pca_result_median_TTS, data_median_TTS, "Median"); #pca_plot_median_TTS
pca_result_SD = prcomp(data_SD[,2:20002],scale.=TRUE)
pca_plot_SD = plot_PCA(pca_result_SD, data_SD, "SD"); #pca_plot_SD
pca_result_SD_TSS = prcomp(data_SD_TSS[,2:20002],scale.=TRUE)
pca_plot_SD_TSS = plot_PCA(pca_result_SD_TSS, data_SD_TSS, "SD"); #pca_plot_SD_TSS
pca_result_SD_TTS = prcomp(data_SD_TTS[,2:20002],scale.=TRUE)
pca_plot_SD_TTS = plot_PCA(pca_result_SD_TTS, data_SD_TTS, "SD"); #pca_plot_SD_TTS


plot_line_median(combined_data_TSS_filtered[combined_data_TSS_filtered$Class %in% c("m", "h", "G", "T"),], 1, empty) -> mammal_line
plot_line_median(combined_data_TSS_filtered[combined_data_TSS_filtered$Assembly %in% c("GRCh38.p14", "T2T-CHM13v2.0", "mGorGor1", "mMacNem1", "mNycCou1", "mPanPan1", "mPanTro3", "mPonAbe1", "mPonPyg2", "mSymSyn1"),], 1, empty) -> primate_line
plot_line_median(combined_data_TSS_filtered[combined_data_TSS_filtered$Class %in% c("m") & !combined_data_TSS_filtered$Assembly %in% c("GRCh38.p14", "T2T-CHM13v2.0", "mGorGor1", "mMacNem1", "mNycCou1", "mPanPan1", "mPanTro3", "mPonAbe1", "mPonPyg2", "mSymSyn1"),], 1, empty) -> non_primate_line
median_TSS_row = align_plots(Yaxis_line + theme(plot.margin = margin(8,0,0,0)), mammal_line + theme(plot.margin = margin(8,8,0,0)), primate_line + theme(plot.margin = margin(8,8,0,0)), non_primate_line + theme(plot.margin = margin(8,8,0,0)), align="h", axis="tb")

plot_line_mean(combined_data_TSS_filtered[combined_data_TSS_filtered$Class %in% c("m", "h", "G", "T"),], 1, empty) -> mammal_line
plot_line_mean(combined_data_TSS_filtered[combined_data_TSS_filtered$Assembly %in% c("GRCh38.p14", "T2T-CHM13v2.0", "mGorGor1", "mMacNem1", "mNycCou1", "mPanPan1", "mPanTro3", "mPonAbe1", "mPonPyg2", "mSymSyn1"),], 1, empty) -> primate_line
plot_line_mean(combined_data_TSS_filtered[combined_data_TSS_filtered$Class %in% c("m") & !combined_data_TSS_filtered$Assembly %in% c("GRCh38.p14", "T2T-CHM13v2.0", "mGorGor1", "mMacNem1", "mNycCou1", "mPanPan1", "mPanTro3", "mPonAbe1", "mPonPyg2", "mSymSyn1"),], 1, empty) -> non_primate_line
mean_TSS_row = align_plots(Yaxis_line + theme(plot.margin = margin(8,0,0,0)), mammal_line + theme(plot.margin = margin(8,8,0,0)), primate_line + theme(plot.margin = margin(8,8,0,0)), non_primate_line + theme(plot.margin = margin(8,8,0,0)), align="h", axis="tb")


# Axis + plottin together w/ PCA
axis_row = align_plots(mammal_line + theme(plot.margin=margin(8,8,0,0)), Xaxis_line + theme(plot.margin=margin(0,8,0,0)), align="v", axis="lr")
axis_blank = align_plots(Yaxis_line + theme(plot.margin = margin(8,0,0,0)), blank + theme(plot.margin = margin(0,0,0,0)), align="v", axis="lr")

swinded_label=ggdraw() + draw_grob(textGrob("Distance from TSS (kbp)", gp = gpar(fontsize = 5.5, fontface = "plain")))

TSS=plot_grid(mean_TSS_row[[1]], mean_TSS_row[[2]], mean_TSS_row[[3]], mean_TSS_row[[4]], blank, pca_plot_mean_TSS, pca_plot_SD_TSS,
              median_TSS_row[[1]], median_TSS_row[[2]], median_TSS_row[[3]], median_TSS_row[[4]], blank, pca_plot_median_TSS, legend,
              axis_blank[[2]], axis_row[[2]], axis_row[[2]], axis_row[[2]], blank, blank, blank,
              blank, blank, swinded_label, blank, blank, blank, blank,
              ncol=7, rel_heights=c(10,10,0.5, 2), rel_widths=c(1,10,10,10,1,10,10),
              labels=c("","","","","","b"), label_size=8)
rotated_label_median <- ggdraw() + draw_grob(textGrob("Median MP (%)", rot = 90, gp = gpar(fontsize = 5.5, fontface = "plain")))
rotated_label_mean <- ggdraw() + draw_grob(textGrob("Mean MP (%)", rot = 90, gp = gpar(fontsize = 5.5, fontface = "plain")))
rotated_label = plot_grid(rotated_label_mean, rotated_label_median, blank, ncol=1, rel_heights=c(10,10,3.5), labels=c("a"), label_size=8)
TSS_quarterTitled = plot_grid(rotated_label, TSS, ncol=2, rel_widths=c(2,52))
TSS_titled = plot_grid(blank, TSS_quarterTitled, ncol=1, rel_heights=c(1,22.5), labels=c("Transcript start site (TSS) in mammals"), label_size=7, label_x=0.5, hjust=0.5, label_y=0.5, vjust=0.5)

# Mammal - TTS

plot_line_median(combined_data_TTS_filtered[combined_data_TTS_filtered$Class %in% c("m", "h", "G", "T"),], 1, empty) -> mammal_line
plot_line_median(combined_data_TTS_filtered[combined_data_TTS_filtered$Assembly %in% c("GRCh38.p14", "T2T-CHM13v2.0", "mGorGor1", "mMacNem1", "mNycCou1", "mPanPan1", "mPanTro3", "mPonAbe1", "mPonPyg2", "mSymSyn1"),], 1, empty) -> primate_line
plot_line_median(combined_data_TTS_filtered[combined_data_TTS_filtered$Class %in% c("m") & !combined_data_TTS_filtered$Assembly %in% c("GRCh38.p14", "T2T-CHM13v2.0", "mGorGor1", "mMacNem1", "mNycCou1", "mPanPan1", "mPanTro3", "mPonAbe1", "mPonPyg2", "mSymSyn1"),], 1, empty) -> non_primate_line
median_TTS_row = align_plots(Yaxis_line + theme(plot.margin = margin(8,0,0,0)), mammal_line + theme(plot.margin = margin(8,8,0,0)), primate_line + theme(plot.margin = margin(8,8,0,0)), non_primate_line + theme(plot.margin = margin(8,8,0,0)), align="h", axis="tb")

plot_line_mean(combined_data_TTS_filtered[combined_data_TTS_filtered$Class %in% c("m", "h", "G", "T"),], 1, empty) -> mammal_line
plot_line_mean(combined_data_TTS_filtered[combined_data_TTS_filtered$Assembly %in% c("GRCh38.p14", "T2T-CHM13v2.0", "mGorGor1", "mMacNem1", "mNycCou1", "mPanPan1", "mPanTro3", "mPonAbe1", "mPonPyg2", "mSymSyn1"),], 1, empty) -> primate_line
plot_line_mean(combined_data_TTS_filtered[combined_data_TTS_filtered$Class %in% c("m") & !combined_data_TTS_filtered$Assembly %in% c("GRCh38.p14", "T2T-CHM13v2.0", "mGorGor1", "mMacNem1", "mNycCou1", "mPanPan1", "mPanTro3", "mPonAbe1", "mPonPyg2", "mSymSyn1"),], 1, empty) -> non_primate_line
mean_TTS_row = align_plots(Yaxis_line + theme(plot.margin = margin(8,0,0,0)), mammal_line + theme(plot.margin = margin(8,8,0,0)), primate_line + theme(plot.margin = margin(8,8,0,0)), non_primate_line + theme(plot.margin = margin(8,8,0,0)), align="h", axis="tb")


# Axis + plottin together w/ PCA
axis_row = align_plots(mammal_line + theme(plot.margin=margin(8,8,0,0)), Xaxis_line + theme(plot.margin=margin(0,8,0,0)), align="v", axis="lr")
axis_blank = align_plots(Yaxis_line + theme(plot.margin = margin(8,0,0,0)), blank + theme(plot.margin = margin(0,0,0,0)), align="v", axis="lr")

swinded_label=ggdraw() + draw_grob(textGrob("Distance from TTS (kbp)", gp = gpar(fontsize = 5.5, fontface = "plain")))

TTS=plot_grid(mean_TTS_row[[1]], mean_TTS_row[[2]], mean_TTS_row[[3]], mean_TTS_row[[4]], blank, pca_plot_mean_TTS, pca_plot_SD_TTS,
              median_TTS_row[[1]], median_TTS_row[[2]], median_TTS_row[[3]], median_TTS_row[[4]], blank, pca_plot_median_TTS, legend,
              axis_blank[[2]], axis_row[[2]], axis_row[[2]], axis_row[[2]], blank, blank, blank,
              blank, blank, swinded_label, blank, blank, blank, blank,
              ncol=7, rel_heights=c(10,10,0.5, 2), rel_widths=c(1,10,10,10,1,10,10),
              labels=c("","","","","","d"), label_size=8)
rotated_label_median <- ggdraw() + draw_grob(textGrob("Median MP (%)", rot = 90, gp = gpar(fontsize = 5.5, fontface = "plain")))
rotated_label_mean <- ggdraw() + draw_grob(textGrob("Mean MP (%)", rot = 90, gp = gpar(fontsize = 5.5, fontface = "plain")))
rotated_label = plot_grid(rotated_label_mean, rotated_label_median, blank, ncol=1, rel_heights=c(10,10,3.5), labels=c("c"), label_size=8)
TTS_quarterTitled = plot_grid(rotated_label, TTS, ncol=2, rel_widths=c(2,52))
TTS_titled = plot_grid(blank, TTS_quarterTitled, ncol=1, rel_heights=c(1,22.5), labels=c("Transcript termination site (TTS) in mammals"), label_size=7, label_x=0.5, hjust=0.5, label_y=0.5, vjust=0.5)


# Bird - TSS

plot_line_mean<-function(df, ALPHA, SVG){
  ggplot(df, aes(x=Distance, y=Mean, group=Assembly)) +
    #annotation_raster(img, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + 
    annotation_custom(SVG, xmin = 4000, xmax = 9000, ymin = 0, ymax = 33) +
    geom_line(data=df[df$Assembly %in% c("bTaeGut1.4", "bTaeGut7", "bAgePho1", "bAmmCau1", "bAmmNel1", "bDixPip1", "bHaeMex1", "bMelGeo1"),],linewidth=.1, alpha=ALPHA, color="#3e6939") +
    geom_line(data=df[!df$Assembly %in% c("bTaeGut1.4", "bTaeGut7", "bAgePho1", "bAmmCau1", "bAmmNel1", "bDixPip1", "bHaeMex1", "bMelGeo1"),],linewidth=.1, alpha=ALPHA, color="#bbdb95") +
    # geom_line(data=subset(df, Assembly %in% c("GRCh38.p14", "T2T-CHM13v2.0", "bTaeGut7", "mPanTro3", "mPanPan1", "mGorGor1", "mPonAbe1", "mPonPyg2", "mSymSyn1", "mMacNem1", "mNycCou1", "mCynVol1", "mApoSyl1", "mChiNiv1", "mMarFla1", "mOchPri1", "mEriEur2", "mSorAra1", "bAgePho1", "bAmmCau1", "bAmmNel1", "bMelGeo1", "bHaeMex1")), aes(color=Color), linewidth=.1, alpha=1) +
    #geom_line(data=subset(df, Assembly %in% c("GRCh38.p14", "T2T-CHM13v2.0", "bTaeGut7", "bTaeGut1.4")), color="black", linewidth=.15, alpha=1) +
    scale_y_continuous(limits=c(0,100), expand=c(0,0)) +
    scale_x_continuous(labels = label_comma(), breaks = c(-10000, 0, 10000), limits = c(-10000, 10000), expand=c(0,0)) +
    #geom_text(data=label_data, aes(label=`Common name`), hjust=0, nudge_x=-5, size=3) +
    labs(y="Average MP (%)", x="Distance (bp)") + 
    theme(legend.position="none",
          panel.background = element_blank(),
          panel.grid = element_blank(),
          plot.background = element_blank(),
          axis.line.y.left = element_line(linewidth=0.1),
          axis.line.x.bottom = element_line(linewidth=0.1),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank())
}

plot_line_median<-function(df, ALPHA, SVG){
  ggplot(df, aes(x=Distance, y=Median, group=Assembly)) +
    #annotation_raster(img, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + 
    annotation_custom(SVG, xmin = 4000, xmax = 9000, ymin = 0, ymax = 33) +
    geom_line(data=df[df$Assembly %in% c("bTaeGut1.4", "bTaeGut7", "bAgePho1", "bAmmCau1", "bAmmNel1", "bDixPip1", "bHaeMex1", "bMelGeo1"),],linewidth=.1, alpha=ALPHA, color="#3e6939") +
    geom_line(data=df[!df$Assembly %in% c("bTaeGut1.4", "bTaeGut7", "bAgePho1", "bAmmCau1", "bAmmNel1", "bDixPip1", "bHaeMex1", "bMelGeo1"),],linewidth=.1, alpha=ALPHA, color="#bbdb95") +
    # geom_line(data=subset(df, Assembly %in% c("GRCh38.p14", "T2T-CHM13v2.0", "bTaeGut7", "mPanTro3", "mPanPan1", "mGorGor1", "mPonAbe1", "mPonPyg2", "mSymSyn1", "mMacNem1", "mNycCou1", "mCynVol1", "mApoSyl1", "mChiNiv1", "mMarFla1", "mOchPri1", "mEriEur2", "mSorAra1", "bAgePho1", "bAmmCau1", "bAmmNel1", "bMelGeo1", "bHaeMex1")), aes(color=Color), linewidth=.1, alpha=1) +
    #geom_line(data=subset(df, Assembly %in% c("GRCh38.p14", "T2T-CHM13v2.0", "bTaeGut7", "bTaeGut1.4")), color="black", linewidth=.15, alpha=1) +
    scale_y_continuous(limits=c(0,100), expand=c(0,0)) +
    scale_x_continuous(labels = label_comma(), breaks = c(-10000, 0, 10000), limits = c(-10000, 10000), expand=c(0,0)) +
    #geom_text(data=label_data, aes(label=`Common name`), hjust=0, nudge_x=-5, size=3) +
    labs(y="Average MP (%)", x="Distance (bp)") + 
    theme(legend.position="none",
          panel.background = element_blank(),
          panel.grid = element_blank(),
          plot.background = element_blank(),
          axis.line.y.left = element_line(linewidth=0.1),
          axis.line.x.bottom = element_line(linewidth=0.1),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank())
}

# Create legend
class_order <- c("Passerine", "Non-passerine bird")
df <- data.frame(X = 1,Y = 1,Class = c("Passerine", "Non-passerine bird"))
df$Class <- factor(df$Class, levels = class_order)
fake <- ggplot(df, aes(x = X, y = Y, color = Class)) +
  geom_point(size = 3) +
  scale_color_manual(name="Order:", values = c(
    "Passerine"="#3e6939",
    "Non-passerine bird"="#bbdb95"
  )) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 2))) +
  theme(
    #legend.position = "top",
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill="transparent", color=NA),
    legend.title = element_text(size=6, face="bold"),
    legend.text = element_text(size=6),
    legend.key.size=unit(.75,"lines")
  )
legend = as_ggplot(get_legend(fake));legend

plot_PCA <- function(pca_result,whole_data,title){
  assembly_colors = c("bTaeGut1.4"="#3e6939", "bTaeGut7"="#3e6939", "bAgePho1"="#3e6939", "bAmmCau1"="#3e6939", "bAmmNel1"="#3e6939", "bDixPip1"="#3e6939", "bHaeMex1"="#3e6939", "bMelGeo1"="#3e6939")
  explained_var <- (pca_result$sdev^2) / sum(pca_result$sdev^2)
  percent_var <- round(100 * explained_var, 1)
  
  scores <- as.data.frame(pca_result$x)
  scores$Assembly <- whole_data$Assembly
  scores$Color <- whole_data$Color
  scores$common_name <- whole_data$"Common name"
  # Plot
  ggplot(scores, aes(x=PC1, y=PC2, label = common_name)) +
    geom_text(vjust = 0, size = 2, color="transparent") + # Change vjust to 0 and color="transparent" to aes(color=Color) to display names instead of dots
    geom_point(aes(color=Assembly), size=0.3) + # Change aes(color=Color) to color="transparent" to display names instead of dots
    # ↑ If you want to fix size instead of setting it by Size column of scores data frame, remove "size=Size" inside the aes() argument and add "size=<your size>" inside the geom_point() argument
    # scale_color_identity() + # Use the actual hex codes from the Color column
    scale_color_manual(values=assembly_colors, na.value="#bbdb95") +
    scale_shape_identity() +
    labs(x = paste0("PC1 (", percent_var[1], "%)"), y = paste0("PC2 (", percent_var[2], "%)"), title=title) +
    coord_cartesian(clip = "off") + #To prevent clipping-off of geom_texts that go out of boundaries
    #geom_text(data=subset(scores, Assembly %in% c("GRCh38.p14", "T2T-CHM13v2.0")), aes(label="★", color=Color), size=5)  +
    theme(
      legend.position="none",
      plot.title = element_text(size=5.5),
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.line.y.left = element_line(),
      axis.line.x.bottom = element_line(),
      panel.grid = element_blank(),
      axis.line=element_line(linewidth=0.1),
      axis.ticks = element_blank(),
      axis.title=element_text(size=5.5),
      axis.text=element_blank(),
      plot.margin = unit(c(8,8,0,0), "pt")
    )->plotted
  return(plotted)
}

data_matrix_normalized = normalize_numbers(data_matrix)
data_matrix_TSS_normalized = normalize_numbers(data_matrix_TSS)
data_matrix_TTS_normalized = normalize_numbers(data_matrix_TTS)
#data_matrix_normalized=data_matrix # If you DON'T want to normalize
#data_matrix_TSS_normalized=data_matrix_TSS # If you DON'T want to normalize
#data_matrix_TTS_normalized=data_matrix_TTS # If you DON'T want to normalize

comprehensive <- fread("C:/Users/biopop/Documents/Downloads/MP analysis for all VGP Freeze1.0 species - Sheet1.tsv")
data_matrix_normalized <- merge(data_matrix_normalized, comprehensive, by="Assembly", all.x=TRUE)
data_matrix_TSS_normalized <- merge(data_matrix_TSS_normalized, comprehensive, by="Assembly", all.x=TRUE)
data_matrix_TTS_normalized <- merge(data_matrix_TTS_normalized, comprehensive, by="Assembly", all.x=TRUE)
# Now, in data_matrix, "Assembly" becomes column 1, columns 2-20,002 are MP, 20,011 is color.

# Remove outliers
data_matrix_normalized <- filter(data_matrix_normalized, !`Common name` %in% c("Wood mouse", "African elephant", "Yellowfin tuna", "Crucian carp", "Northern goshawk", "European flounder") & `Assembly` != "bTaeGut7")
data_matrix_TSS_normalized <- filter(data_matrix_TSS_normalized, !`Common name` %in% c("Wood mouse", "African elephant", "Yellowfin tuna", "Crucian carp", "Northern goshawk", "European flounder") & `Assembly` != "bTaeGut7")
data_matrix_TTS_normalized <- filter(data_matrix_TTS_normalized, !`Common name` %in% c("Wood mouse", "African elephant", "Yellowfin tuna", "Crucian carp", "Northern goshawk", "European flounder") & `Assembly` != "bTaeGut7")
data_matrix_normalized$`Class (강)`[data_matrix_normalized$Species=="Homo sapiens"]="Human"
data_matrix_TSS_normalized$`Class (강)`[data_matrix_TSS_normalized$Species=="Homo sapiens"]="Human"
data_matrix_TTS_normalized$`Class (강)`[data_matrix_TTS_normalized$Species=="Homo sapiens"]="Human"

data_mean <- data_matrix_normalized[seq(1, nrow(data_matrix_normalized), by = 3), ]
data_median <- data_matrix_normalized[seq(2, nrow(data_matrix_normalized), by = 3), ]
data_SD <- data_matrix_normalized[seq(3, nrow(data_matrix_normalized), by = 3), ]
data_mean_TSS <- data_matrix_TSS_normalized[seq(1, nrow(data_matrix_TSS_normalized), by = 3), ]
data_median_TSS <- data_matrix_TSS_normalized[seq(2, nrow(data_matrix_TSS_normalized), by = 3), ]
data_SD_TSS <- data_matrix_TSS_normalized[seq(3, nrow(data_matrix_TSS_normalized), by = 3), ]
data_mean_TTS <- data_matrix_TTS_normalized[seq(1, nrow(data_matrix_TTS_normalized), by = 3), ]
data_median_TTS <- data_matrix_TTS_normalized[seq(2, nrow(data_matrix_TTS_normalized), by = 3), ]
data_SD_TTS <- data_matrix_TTS_normalized[seq(3, nrow(data_matrix_TTS_normalized), by = 3), ]

data_mean = data_mean[!is.na(data_mean$`Class (강)`) & data_mean$"Class (강)"=="Aves",]
data_mean_TSS = data_mean_TSS[!is.na(data_mean_TSS$`Class (강)`) & data_mean_TSS$"Class (강)"=="Aves",]
data_mean_TTS = data_mean_TTS[!is.na(data_mean_TTS$`Class (강)`) & data_mean_TTS$"Class (강)"=="Aves",]
data_median = data_median[!is.na(data_median$`Class (강)`) & data_median$"Class (강)"=="Aves",]
data_median_TSS = data_median_TSS[!is.na(data_median_TSS$`Class (강)`) & data_median_TSS$"Class (강)"=="Aves",]
data_median_TTS = data_median_TTS[!is.na(data_median_TTS$`Class (강)`) & data_median_TTS$"Class (강)"=="Aves",]
data_SD = data_SD[!is.na(data_SD$`Class (강)`) & data_SD$"Class (강)"=="Aves",]
data_SD_TSS = data_SD_TSS[!is.na(data_SD_TSS$`Class (강)`) & data_SD_TSS$"Class (강)"=="Aves",]
data_SD_TTS = data_SD_TTS[!is.na(data_SD_TTS$`Class (강)`) & data_SD_TTS$"Class (강)"=="Aves",]

pca_result_mean = prcomp(data_mean[,2:20002],scale.=TRUE)
pca_plot_mean = plot_PCA(pca_result_mean,data_mean, "Mean") #pca_plot_mean
pca_result_mean_TSS = prcomp(data_mean_TSS[,2:20002],scale.=TRUE)
pca_plot_mean_TSS = plot_PCA(pca_result_mean_TSS, data_mean_TSS, "Mean"); #pca_plot_mean_TSS
pca_result_mean_TTS = prcomp(data_mean_TTS[,2:20002],scale.=TRUE)
pca_plot_mean_TTS = plot_PCA(pca_result_mean_TTS, data_mean_TTS, "Mean"); #pca_plot_mean_TTS
pca_result_median = prcomp(data_median[,2:20002],scale.=TRUE)
pca_plot_median = plot_PCA(pca_result_median, data_median, "Median"); #pca_plot_median
pca_result_median_TSS = prcomp(data_median_TSS[,2:20002],scale.=TRUE)
pca_plot_median_TSS = plot_PCA(pca_result_median_TSS, data_median_TSS, "Median"); #pca_plot_median_TSS
pca_result_median_TTS = prcomp(data_median_TTS[,2:20002],scale.=TRUE)
pca_plot_median_TTS = plot_PCA(pca_result_median_TTS, data_median_TTS, "Median"); #pca_plot_median_TTS
pca_result_SD = prcomp(data_SD[,2:20002],scale.=TRUE)
pca_plot_SD = plot_PCA(pca_result_SD, data_SD, "SD"); #pca_plot_SD
pca_result_SD_TSS = prcomp(data_SD_TSS[,2:20002],scale.=TRUE)
pca_plot_SD_TSS = plot_PCA(pca_result_SD_TSS, data_SD_TSS, "SD"); #pca_plot_SD_TSS
pca_result_SD_TTS = prcomp(data_SD_TTS[,2:20002],scale.=TRUE)
pca_plot_SD_TTS = plot_PCA(pca_result_SD_TTS, data_SD_TTS, "SD"); #pca_plot_SD_TTS


plot_line_median(combined_data_TSS_filtered[combined_data_TSS_filtered$Class %in% c("b"),], 1, empty) -> bird_line
plot_line_median(combined_data_TSS_filtered[combined_data_TSS_filtered$Assembly %in% c("bTaeGut1.4", "bTaeGut7", "bAgePho1", "bAmmCau1", "bAmmNel1", "bDixPip1", "bHaeMex1", "bMelGeo1"),], 1, empty) -> pass_line
plot_line_median(combined_data_TSS_filtered[combined_data_TSS_filtered$Class %in% c("b") & !combined_data_TSS_filtered$Assembly %in% c("bTaeGut1.4", "bTaeGut7", "bAgePho1", "bAmmCau1", "bAmmNel1", "bDixPip1", "bHaeMex1", "bMelGeo1"),], 1, empty) -> non_pass_line
median_TSS_row = align_plots(Yaxis_line + theme(plot.margin = margin(8,0,0,0)), bird_line + theme(plot.margin = margin(8,8,0,0)), pass_line + theme(plot.margin = margin(8,8,0,0)), non_pass_line + theme(plot.margin = margin(8,8,0,0)), align="h", axis="tb")

plot_line_mean(combined_data_TSS_filtered[combined_data_TSS_filtered$Class %in% c("b"),], 1, empty) -> bird_line
plot_line_mean(combined_data_TSS_filtered[combined_data_TSS_filtered$Assembly %in% c("bTaeGut1.4", "bTaeGut7", "bAgePho1", "bAmmCau1", "bAmmNel1", "bDixPip1", "bHaeMex1", "bMelGeo1"),], 1, empty) -> pass_line
plot_line_mean(combined_data_TSS_filtered[combined_data_TSS_filtered$Class %in% c("b") & !combined_data_TSS_filtered$Assembly %in% c("bTaeGut1.4", "bTaeGut7", "bAgePho1", "bAmmCau1", "bAmmNel1", "bDixPip1", "bHaeMex1", "bMelGeo1"),], 1, empty) -> non_pass_line
mean_TSS_row = align_plots(Yaxis_line + theme(plot.margin = margin(8,0,0,0)), bird_line + theme(plot.margin = margin(8,8,0,0)), pass_line + theme(plot.margin = margin(8,8,0,0)), non_pass_line + theme(plot.margin = margin(8,8,0,0)), align="h", axis="tb")


# Axis + plottin together w/ PCA
axis_row = align_plots(bird_line + theme(plot.margin=margin(8,8,0,0)), Xaxis_line + theme(plot.margin=margin(0,8,0,0)), align="v", axis="lr")
axis_blank = align_plots(Yaxis_line + theme(plot.margin = margin(8,0,0,0)), blank + theme(plot.margin = margin(0,0,0,0)), align="v", axis="lr")

swinded_label=ggdraw() + draw_grob(textGrob("Distance from TSS (kbp)", gp = gpar(fontsize = 5.5, fontface = "plain")))

TSS=plot_grid(mean_TSS_row[[1]], mean_TSS_row[[2]], mean_TSS_row[[3]], mean_TSS_row[[4]], blank, pca_plot_mean_TSS, pca_plot_SD_TSS,
              median_TSS_row[[1]], median_TSS_row[[2]], median_TSS_row[[3]], median_TSS_row[[4]], blank, pca_plot_median_TSS, legend,
              axis_blank[[2]], axis_row[[2]], axis_row[[2]], axis_row[[2]], blank, blank, blank,
              blank, blank, swinded_label, blank, blank, blank, blank,
              ncol=7, rel_heights=c(10,10,0.5, 2), rel_widths=c(1,10,10,10,1,10,10),
              labels=c("","","","","","b"), label_size=8)
rotated_label_median <- ggdraw() + draw_grob(textGrob("Median MP (%)", rot = 90, gp = gpar(fontsize = 5.5, fontface = "plain")))
rotated_label_mean <- ggdraw() + draw_grob(textGrob("Mean MP (%)", rot = 90, gp = gpar(fontsize = 5.5, fontface = "plain")))
rotated_label = plot_grid(rotated_label_mean, rotated_label_median, blank, ncol=1, rel_heights=c(10,10,3.5), labels=c("a"), label_size=8)
TSS_quarterTitled = plot_grid(rotated_label, TSS, ncol=2, rel_widths=c(2,52))
TSS_titled_bird = plot_grid(blank, TSS_quarterTitled, ncol=1, rel_heights=c(1,22.5), labels=c("Transcript start site (TSS) in birds"), label_size=7, label_x=0.5, hjust=0.5, label_y=0.5, vjust=0.5)


# Bird - TTS

plot_line_median(combined_data_TTS_filtered[combined_data_TTS_filtered$Class %in% c("b"),], 1, empty) -> bird_line
plot_line_median(combined_data_TTS_filtered[combined_data_TTS_filtered$Assembly %in% c("bTaeGut1.4", "bTaeGut7", "bAgePho1", "bAmmCau1", "bAmmNel1", "bDixPip1", "bHaeMex1", "bMelGeo1"),], 1, empty) -> pass_line
plot_line_median(combined_data_TTS_filtered[combined_data_TTS_filtered$Class %in% c("b") & !combined_data_TTS_filtered$Assembly %in% c("bTaeGut1.4", "bTaeGut7", "bAgePho1", "bAmmCau1", "bAmmNel1", "bDixPip1", "bHaeMex1", "bMelGeo1"),], 1, empty) -> non_pass_line
median_TTS_row = align_plots(Yaxis_line + theme(plot.margin = margin(8,0,0,0)), bird_line + theme(plot.margin = margin(8,8,0,0)), pass_line + theme(plot.margin = margin(8,8,0,0)), non_pass_line + theme(plot.margin = margin(8,8,0,0)), align="h", axis="tb")

plot_line_mean(combined_data_TTS_filtered[combined_data_TTS_filtered$Class %in% c("b"),], 1, empty) -> bird_line
plot_line_mean(combined_data_TTS_filtered[combined_data_TTS_filtered$Assembly %in% c("bTaeGut1.4", "bTaeGut7", "bAgePho1", "bAmmCau1", "bAmmNel1", "bDixPip1", "bHaeMex1", "bMelGeo1"),], 1, empty) -> pass_line
plot_line_mean(combined_data_TTS_filtered[combined_data_TTS_filtered$Class %in% c("b") & !combined_data_TTS_filtered$Assembly %in% c("bTaeGut1.4", "bTaeGut7", "bAgePho1", "bAmmCau1", "bAmmNel1", "bDixPip1", "bHaeMex1", "bMelGeo1"),], 1, empty) -> non_pass_line
mean_TTS_row = align_plots(Yaxis_line + theme(plot.margin = margin(8,0,0,0)), bird_line + theme(plot.margin = margin(8,8,0,0)), pass_line + theme(plot.margin = margin(8,8,0,0)), non_pass_line + theme(plot.margin = margin(8,8,0,0)), align="h", axis="tb")


# Axis + plottin together w/ PCA
axis_row = align_plots(bird_line + theme(plot.margin=margin(8,8,0,0)), Xaxis_line + theme(plot.margin=margin(0,8,0,0)), align="v", axis="lr")
axis_blank = align_plots(Yaxis_line + theme(plot.margin = margin(8,0,0,0)), blank + theme(plot.margin = margin(0,0,0,0)), align="v", axis="lr")

swinded_label=ggdraw() + draw_grob(textGrob("Distance from TSS (kbp)", gp = gpar(fontsize = 5.5, fontface = "plain")))

TTS=plot_grid(mean_TTS_row[[1]], mean_TTS_row[[2]], mean_TTS_row[[3]], mean_TTS_row[[4]], blank, pca_plot_mean_TTS, pca_plot_SD_TTS,
              median_TTS_row[[1]], median_TTS_row[[2]], median_TTS_row[[3]], median_TTS_row[[4]], blank, pca_plot_median_TTS, legend,
              axis_blank[[2]], axis_row[[2]], axis_row[[2]], axis_row[[2]], blank, blank, blank,
              blank, blank, swinded_label, blank, blank, blank, blank,
              ncol=7, rel_heights=c(10,10,0.5, 2), rel_widths=c(1,10,10,10,1,10,10),
              labels=c("","","","","","d"), label_size=8)
rotated_label_median <- ggdraw() + draw_grob(textGrob("Median MP (%)", rot = 90, gp = gpar(fontsize = 5.5, fontface = "plain")))
rotated_label_mean <- ggdraw() + draw_grob(textGrob("Mean MP (%)", rot = 90, gp = gpar(fontsize = 5.5, fontface = "plain")))
rotated_label = plot_grid(rotated_label_mean, rotated_label_median, blank, ncol=1, rel_heights=c(10,10,3.5), labels=c("c"), label_size=8)
TTS_quarterTitled = plot_grid(rotated_label, TTS, ncol=2, rel_widths=c(2,52))
TTS_titled_bird = plot_grid(blank, TTS_quarterTitled, ncol=1, rel_heights=c(1,22.5), labels=c("Transcript termination site (TTS) in birds"), label_size=7, label_x=0.5, hjust=0.5, label_y=0.5, vjust=0.5)








supp4 = plot_grid(TSS_titled, blank, TTS_titled,
                     ncol=1, rel_heights=c(23.5, 1, 23.5))

supp5 = plot_grid(TSS_titled_bird, blank, TTS_titled_bird,
                  ncol=1, rel_heights=c(23.5, 1, 23.5))

widths = 54
heights=sum(23.5, 1, 23.5)

final_width=135
final_height=135/widths*heights

ggsave("supp4.png", supp4, units="mm", height=final_height, width=final_width)
ggsave("supp5.png", supp5, units="mm", height=final_height, width=final_width)
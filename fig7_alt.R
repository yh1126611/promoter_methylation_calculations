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
library(gmodels)
library(svgparser)

setwd("D:/R/WD/MP_TSS/fig7")

##

# Intro

df = fread("MP_transcriptTSS_T2T-CHM13v2.0.tsv")
colnames(df)=c("ID", "Dist", "MP", "Strand")
df$Dist = ifelse(df$Strand=="-", df$Dist * -1, ifelse(df$Strand=="+", df$Dist, "Wrong"))
df = data.frame(summarise(group_by(df, Dist), "Mean"=mean(MP), "Assembly"="T2T-CHM13v2.0"))
df$Dist=as.numeric(df$Dist)
ggplot(df) +
  #geom_rect(ymin=0,ymax=87.5,xmin=-874,xmax=948,fill="#C7E8FB")+
  #geom_rect(ymin=0,ymax=50,xmin=-140,xmax=0,fill="#FFD1DC")+
  #geom_segment(y=0,yend=87.5,x=-874,xend=-874,color="#009FE9",linewidth=0.2,linetype="dashed") + geom_segment(y=0,yend=87.5,x=948,xend=948,color="#009FE9",linewidth=0.2,linetype="dashed") +
  #geom_segment(y=0,yend=50,x=-140,xend=-140,color="#DC143C",linewidth=0.2,linetype="dashed") + geom_segment(y=0,yend=50,x=0,xend=0,color="#DC143C",linewidth=0.2,linetype="dashed") +
  #annotate("text", y=100, x=-2000, vjust=1, hjust=0, label = "Human (T2T-CHM13v2.0)", color = "#000000", size = 1.94, fontface="bold") +
  geom_line(data=df[df$Dist<0,], mapping=aes(x=Dist, y=Mean, group=Assembly), linewidth=0.1, color="#FF02EB") +
  geom_line(data=df[df$Dist>=0,], mapping=aes(x=Dist, y=Mean, group=Assembly), linewidth=0.1, color="black") +
  #annotate("text", y=df[df$Dist==0,]$Mean + 5, x=0, vjust=0, hjust=.5, label = "▼", color = "black", size = 2.25) +
  annotate("text", y=df[df$Dist==-2000,]$Mean+10, x=-2000+100, vjust=0, hjust=0, label = "5’", color = "#FF02EB", size = 1.94) +
  annotate("text", y=df[df$Dist==2000,]$Mean+10, x=2000-100, vjust=0, hjust=1, label = "3’", color = "black", size = 1.94) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,25,50,75,100), name="\nMean MP (%)", expand=c(0,0)) +
  scale_x_continuous(limits=c(-2000, 2000), breaks=c(-2000,0,2000), labels=label_comma(), name="Distance from TSS (bp)", expand=c(0,0)) +
  theme(plot.background=element_blank(), panel.background = element_blank(), panel.grid=element_blank(),
        axis.line=element_line(color="black", linewidth=0.1), axis.ticks = element_blank(),
        axis.title=element_text(color="black", size=5), axis.text=element_text(color="black", size=5))->human_basic

df = fread("MP_transcriptTSS_bTaeGut7.tsv")
colnames(df)=c("ID", "Dist", "MP", "Strand")
df$Dist = ifelse(df$Strand=="-", df$Dist * -1, ifelse(df$Strand=="+", df$Dist, "Wrong"))
df = data.frame(summarise(group_by(df, Dist), "Mean"=mean(MP), "Assembly"="bTaeGut7"))
df$Dist=as.numeric(df$Dist)
ggplot(df) +
  #geom_rect(ymin=0,ymax=87.5,xmin=-874,xmax=948,fill="#C7E8FB")+
  #geom_rect(ymin=0,ymax=50,xmin=-140,xmax=0,fill="#FFD1DC")+
  #geom_segment(y=0,yend=87.5,x=-874,xend=-874,color="#009FE9",linewidth=0.2,linetype="dashed") + geom_segment(y=0,yend=87.5,x=948,xend=948,color="#009FE9",linewidth=0.2,linetype="dashed") +
  #geom_segment(y=0,yend=50,x=-140,xend=-140,color="#DC143C",linewidth=0.2,linetype="dashed") + geom_segment(y=0,yend=50,x=0,xend=0,color="#DC143C",linewidth=0.2,linetype="dashed") +
  #annotate("text", y=100, x=-2000, vjust=1, hjust=0, label = "Zebra finch (bTaeGut7)", color = "#000000", size = 1.94, fontface="bold") +
  geom_line(data=df[df$Dist<0,], mapping=aes(x=Dist, y=Mean, group=Assembly), linewidth=0.1, color="#FF02EB") +
  geom_line(data=df[df$Dist>=0,], mapping=aes(x=Dist, y=Mean, group=Assembly), linewidth=0.1, color="black") +
  #annotate("text", y=df[df$Dist==0,]$Mean + 5, x=0, vjust=0, hjust=.5, label = "▼", color = "black", size = 2.25) +
  annotate("text", y=df[df$Dist==-2000,]$Mean+10, x=-2000+100, vjust=0, hjust=0, label = "5’", color = "#FF02EB", size = 1.94) +
  annotate("text", y=df[df$Dist==2000,]$Mean+10, x=2000-100, vjust=0, hjust=1, label = "3’", color = "black", size = 1.94) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,25,50,75,100), name="\nMean MP (%)", expand=c(0,0)) +
  scale_x_continuous(limits=c(-2000, 2000), breaks=c(-2000,0,2000), labels=label_comma(), name="Distance from TSS (bp)", expand=c(0,0)) +
  theme(plot.background=element_blank(), panel.background = element_blank(), panel.grid=element_blank(),
        axis.line=element_line(color="black", linewidth=0.1), axis.ticks = element_blank(),
        axis.title=element_text(color="black", size=5), axis.text=element_text(color="black", size=5))->zf_basic

# Addt'l alt

histones=svgparser::read_svg("svgs/histones.svg")
red_arrow=svgparser::read_svg("svgs/red_arrow.svg")
blue_arrow=svgparser::read_svg("svgs/blue_arrow.svg")
#dna_methylation_tip=svgparser::read_svg("svgs/dna_methylation_tip.svg")
ggplot(data.frame(c(1), c(1)))+ 
  coord_cartesian(clip = "off") +
  geom_rect(xmin=0, xmax=2.5, ymin=0, ymax=10, fill="#E8E8E8") +
  geom_rect(xmin=7.5, xmax=10, ymin=0, ymax=10, fill="#E8E8E8") +
  geom_rect(xmin=3.4375, xmax=6.25, ymin=0, ymax=7.5, fill="#C7E8FB") +
  geom_rect(xmin=3.75, xmax=5, ymin=0, ymax=5, fill="#FFD1DC") +
  geom_segment(yend=7, y=5.5, xend=0.25, x=0.25, linewidth=0.1, color="darkgrey")+
  geom_segment(yend=5.5, y=5.5, xend=0.25, x=0.375, linewidth=0.1, color="darkgrey")+
  geom_segment(yend=7, y=4.5, xend=0.1, x=0.1, linewidth=0.1)+
  geom_segment(yend=4.5, y=4.5, xend=0.1, x=0.375, linewidth=0.1)+
  annotate("text", y=5.5, vjust=.5, x=0.375, hjust=0, label="Histone", size=1.94, color="darkgrey") +
  annotate("text", y=4.5, vjust=.5, x=0.375, hjust=0, label="Genomic DNA", size=1.94) +
  annotate("text", y=3.5, vjust=.5, x=0.375, hjust=0, label="DNA methylation", size=1.94, color="navy") +
  annotation_custom(histones, xmin=0, xmax=10, ymin=0, ymax=10) +
  #annotation_custom(blue_arrow, xmin=3.4375, xmax=6.25, ymin=7, ymax=7.5) +
  #annotation_custom(red_arrow, xmin=3.75, xmax=5, ymin=4.5, ymax=5) +
  #annotation_custom(dna_methylation_tip, xmin=.1375, xmax=.3375, ymin=3, ymax=4) +
  annotate("text", y=9.75, x=1.25, vjust=1, hjust=0.5, label="Heterochromatin\n(inactive)", color="darkgrey", size=1.94) +
  annotate("text", y=9.75, x=5, vjust=1, hjust=0.5, label="Euchromatin\n(active)", color="black", size=1.94) +
  annotate("text", y=9.75, x=8.75, vjust=1, hjust=0.5, label="Heterochromatin\n(inactive)", color="darkgrey", size=1.94) +
  annotate("text", y=0.25, vjust=0, x=0.125, hjust=0, label="5’", size=2.46, color="#FF02EB") +
  annotate("text", y=0.25, vjust=0, x=9.875, hjust=1, label="3’", size=2.46, color="black") +
  geom_segment(x=0, xend=10, y=-0.625, yend=-0.625, linewidth=1, color="#FF02EB") +
  geom_rect(ymin=-1,ymax=-.25,xmin=5,xmax=10,fill="black")+
  annotate("text", y=-0.625, vjust=.5, x=7.5, hjust=.5, label="Gene body", size=1.94, color="white")+
  geom_segment(x=3.4375, xend=3.4375, y=0, yend=7.5, linewidth=0.2, linetype="dashed", alpha=1, color="#009FE9")+
  geom_segment(x=6.25, xend=6.25, y=0, yend=7.5, linewidth=0.2, linetype="dashed", alpha=1, color="#009FE9")+
  annotate("text", y=7.5, x=3.4375+(6.25-3.4375)/2, vjust=1.75, hjust=0.5, label="Promoter", color="#009FE9", size=1.94) +
  geom_segment(x=3.75, xend=3.75, y=0, yend=5, linewidth=0.2, linetype="dashed", alpha=1, color="#DC143C")+
  geom_segment(x=5, xend=5, y=0, yend=5, linewidth=0.2, linetype="dashed", alpha=1, color="#DC143C")+
  annotate("text", y=5, x=3.75+(5-3.75)/2, vjust=-.5, hjust=0.5, label="Core-promoter", fontface="bold", color="#DC143C", size=1.94) +
  #annotate("text", x = 5, y = -1, label = "▲", vjust=1, hjust = .5, fontface="bold", size=2.25) +
  annotate("text", x = 5, y = -1, label = "Transcription start site (TSS)   ", vjust=1, hjust = 1, size=1.94, fontface="bold") +
  scale_y_continuous(limits=c(-1,10), breaks=c(), expand=c(0,0))+
  scale_x_continuous(limits=c(0,10), breaks=c(5), labels=c("Transcription start site (TSS)  "), expand=c(0,0))+
  theme(plot.background=element_blank(), panel.background = element_blank(), panel.grid = element_blank(),
        axis.title=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank() #, axis.ticks.x=element_arrow()
        )->example

example=plot_grid(example, labels=c("c"), label_size=8, hjust=1);example

# Fnct'ns

setwd("D:/R/WD/MP_TSS/fig7")

# bin_size=20; bound=2000
bin_size=10; bound=2000
blank = ggplot() + theme_void()

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
plot_pvalue_heatmap <- function(results_df, bin_col = "Bin", pvalue_col = "p_value", larger_dir_col = "Larger_Dir") {
  #results_df=results_df[results_df$Bin<=bound,]
  df <- results_df %>%
    mutate(
      fill_color = case_when(
        results_df$Larger_Dir == "5'" ~ 1 - .data[[pvalue_col]],   # invert p-value for red gradient
        results_df$Larger_Dir == "3'" ~ -(1 - .data[[pvalue_col]]), # invert p-value for blue gradient as negative
        TRUE ~ NA_real_
      )
    )
  ggplot(df, aes(x = Bin, y = 1, fill = fill_color)) +
    geom_tile(color="transparent") +
    scale_fill_gradientn(
      colors = c("#FF02EB", "white", "black"),
      values = scales::rescale(c(-1, 0, 1)),
      limits = c(-1, 1),
      breaks = c(-1, 0, 1),
      labels = c("0 (3'>5')", "1", "0 (3'>5')"),
      name = "p-value"
    )+
    theme_minimal() +
    #coord_fixed(ratio=10) +
    theme(
      plot.background = element_blank(), panel.background = element_blank(), panel.grid = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.x=element_blank(),
      axis.title.x=element_blank(),
      legend.position="none"
    ) +
    labs(x = "Distance from TSS (bp)") +
    scale_x_continuous(expand=c(0,0))
}
plot_VSline <- function(both, species, cpsize, n, m, SVG){ # m is no. of TSSs
  #left=df[df$Dir=="5'"]; right=df[df$Dir=="3'"]
  #left=data.frame(summarise(group_by(left, Bin),"Mean"=mean(MP),"CI_lower"=ci(MP)[2],"CI_upper"=ci(MP)[3])); right=data.frame(summarise(group_by(right, Bin),"Mean"=mean(MP),"CI_lower"=ci(MP)[2],"CI_upper"=ci(MP)[3]))
  #both=data.frame(); both=rbind(both, cbind(left, "Dir"="5'")); both=rbind(both, cbind(right, "Dir"="3'"))
  #both=both[both$Bin<=bound,]
  ggplot(both)+
    #geom_rect(aes(xmin=0,xmax=as.numeric(cpsize),ymin=trunc(min(both$CI_lower)),ymax=trunc(max(both$CI_upper))+1 ), fill="#FFD1DC")+
    #geom_segment(x=as.numeric(cpsize), xend=as.numeric(cpsize), y=trunc(min(both$CI_lower)), yend=trunc(max(both$CI_upper))+1 , color="#DC143C", linewidth=0.2, linetype="dashed") +
    #annotation_raster(image, xmin=bound-190, xmax=bound, ymin=min(both$CI_lower), ymax=(max(both$CI_upper)-min(both$CI_lower))/2 + min(both$CI_lower)) + 
    geom_ribbon(data=both[both$Dir=="5'",], mapping=aes(y=Mean,x=Bin,ymin=CI_lower,ymax=CI_upper),color="transparent",fill="#FF02EB",alpha=0.25) +
    geom_ribbon(data=both[both$Dir=="3'",], mapping=aes(y=Mean,x=Bin,ymin=CI_lower,ymax=CI_upper),color="transparent",fill="black",alpha=0.25) +
    geom_line(mapping=aes(x=Bin, y=Mean, group=Dir, color=Dir),linewidth=.25)+
    annotation_custom(SVG, xmin = 1300, xmax = 2000, ymin = trunc(min(both$CI_lower)), ymax = trunc(min(both$CI_lower)) + (trunc(max(both$CI_upper))+1 - trunc(min(both$CI_lower))) * 0.66) +
    #annotate("text", y=trunc(max(both$CI_upper))+1, x=0, hjust=0, vjust=1, label = species, color = "black", size = 1.94, fontface="bold") +
    #annotate("text", y=trunc(max(both$CI_upper))+1, x=0, hjust=0, vjust=1, label = paste("\nnTSS =", m), color = "black", size = 1.94) +
    scale_x_continuous(expand=c(0,0))+
    scale_color_manual(name="Location respect to TSS", values=c("5'"="#FF02EB","3'"="black"),na.value="transparent", expand=c(0,0)) +
    scale_y_continuous(limits = c(trunc(min(both$CI_lower)),  trunc(max(both$CI_upper))+1), breaks = c(trunc(min(both$CI_lower)), trunc(max(both$CI_upper)) + 1), name = "\nMean MP (%)", expand=c(0,0))+
    theme(plot.background=element_blank(), panel.background = element_blank(), panel.grid = element_blank(), axis.ticks.x=element_blank(),
          axis.line.y.left = element_line(color="black", linewidth=0.1),    # Hide left y axis line
          axis.text.y.left = element_text(color="black", size=5),    # Hide left y axis texts
          axis.title.y.left = element_text(color="black", size=5),   # Hide left y axis title
          axis.ticks.y=element_blank(),
          #axis.line.y.right = element_line(linewidth = 0.1, color = "black"), # Show right y axis line
          #axis.ticks.y.right = element_blank(),   # Show right y axis ticks
          #axis.text.y.right = element_text(size = 5, color = "black"), # Show right y axis labels
          #axis.title.y.right = element_text(size = 5, color = "black"), # Show right y axis title
          axis.title.x=element_blank(), axis.text.x=element_blank(),
          legend.position="none")
}
plot_fakeVSline <- function(both){ # M is no. of CpGs, m is no. of TSSs
  #left=both[both$Dir=="5'"]; right=both[both$Dir=="3'"]
  #left=data.frame(summarise(group_by(left, Bin),"Mean"=mean(MP),"CI_lower"=ci(MP)[2],"CI_upper"=ci(MP)[3])); right=data.frame(summarise(group_by(right, Bin),"Mean"=mean(MP),"CI_lower"=ci(MP)[2],"CI_upper"=ci(MP)[3]))
  #both=data.frame(); both=rbind(both, cbind(left, "Dir"="5'")); both=rbind(both, cbind(right, "Dir"="3'"))
  #both=both[both$Bin<=bound,]
  ggplot(both)+
    annotate("text", y=both[both$Bin==bound&both$Dir=="5'",]$Mean, x=bound, hjust=.5, vjust=1, label="5’", color="#FF02EB", size=1.94)+
    annotate("text", y=both[both$Bin==bound&both$Dir=="3'",]$Mean, x=bound, hjust=.5, vjust=0, label="3’", color="black", size=1.94)+
    coord_cartesian(clip = "off") +
    scale_x_continuous(limits=c(bound,bound+80), expand=c(0,0))+
    scale_y_continuous(limits = c(trunc(min(both$CI_lower)), trunc(max(both$CI_upper)) + 1), breaks = c(trunc(min(both$CI_lower)), trunc(max(both$CI_upper)) + 1), name = "\nMean MP (%)", expand=c(0,0))+
    theme(plot.background=element_blank(), panel.background = element_blank(), panel.grid = element_blank(), axis.ticks=element_blank(),
          axis.title=element_blank(), axis.text=element_blank(),
          legend.position="none")
}
plot_fakePvalue_heatmap <- function(results_df, bin_col = "Bin", pvalue_col = "p_value", larger_dir_col = "Larger_Dir") {
  #results_df=results_df[results_df$Bin<=bound,]
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
      colors = c("red", "white", "grey", "white", "blue"),
      values = scales::rescale(c(-1, -0.95, 0, 0.95, 1)),
      limits = c(-1, 1),
      breaks = c(-1, -0.95, 0, 0.95, 1),
      labels = c("", "0.05 (3'>5')", "0", "0.05 5'>3'", ""),
      name = "p-value"
    )+
    theme_minimal() +
    #coord_fixed(ratio=10) +
    theme(
      plot.background = element_blank(), panel.background = element_blank(), panel.grid = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.x=element_blank(),
      axis.title.x=element_blank(),
      legend.position="none"
    ) +
    labs(x = "Distance from GSS (bp)") +
    scale_x_continuous(limits=c(-80,0), expand=c(0,0))
}
plot_diff <- function(both, col, cpsize){
  both$Cont = ifelse(both$Dir=="3'", both$Mean * -1, both$Mean)
  both = data.frame(summarise(group_by(both, Bin),"Diff"=sum(Cont)))
  ggplot(both)+
    #annotate("text", y=trunc(min(both$Diff))+(trunc(max(both$Diff))+1 - trunc(min(both$Diff)))/8.3*10, x=0, hjust=0, vjust=1, label = paste("Estimated core promoter size =", cpsize, "bp"), color = "#DC143C", size = 1.94) +
    geom_rect(xmin=0,xmax=as.numeric(cpsize),ymin=trunc(min(both$Diff))-1,ymax=trunc(max(both$Diff))+1, fill="#FFD1DC")+
    geom_segment(x=0,xend=bound,y=0,yend=0, color="grey", linewidth=0.2, linetype="dotted") +
    geom_segment(x=as.numeric(cpsize), xend=as.numeric(cpsize), y=trunc(min(both$Diff))-1, yend=trunc(max(both$Diff))+1, color="#DC143C", linewidth=0.2, linetype="dashed") +
    geom_line(mapping=aes(x=Bin, y=Diff), color=col, linewidth=.25)+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(limits = c(trunc(min(both$Diff))-1, trunc(max(both$Diff))+1 ), breaks = c(trunc(min(both$Diff))-1, 0, trunc(max(both$Diff))+1), name = "\n5’ vs. 3’ asymmetry", expand=c(0,0))+
    theme(plot.background=element_blank(), panel.background = element_blank(), panel.grid = element_blank(), axis.ticks.x=element_blank(),
          axis.line.y.left = element_line(linewidth=0.1, color="black"),    # Hide left y axis line
          axis.text.y.left = element_text(color="black", size=5),    # Hide left y axis texts
          axis.title.y.left = element_text(color="black", size=5),   # Hide left y axis title
          axis.ticks.y=element_blank(),
          #axis.line.y.right = element_line(linewidth = 0.1, color = "black"), # Show right y axis line
          #axis.ticks.y.right = element_blank(),   # Show right y axis ticks
          #axis.text.y.right = element_text(size = 5, color = "black"), # Show right y axis labels
          #axis.title.y.right = element_text(size = 5, color = "black"), # Show right y axis title
          axis.title.x=element_blank(), axis.text.x=element_blank(),
          legend.position="none")
}
plotMedian = function(filename, a, b, color, cpsize){
  df=fread(filename); colnames(df)=c("Dist", "Median"); df=df[df$Dist>=bound*-1&Dist<=bound,]
  realYmax=min(df$Median)+((max(df$Median)-min(df$Median))/8.5*10)
  labelYmax=max(df$Median)
  ggplot(df)+
    coord_cartesian(clip = "off") +
    geom_rect(aes(xmin=a, xmax=b, ymin=min(df$Median), ymax=max(df$Median)), fill="#C7E8FB", show.legend=FALSE)+
    geom_rect(aes(xmin=0-as.numeric(cpsize), xmax=0, ymin=min(df$Median), ymax=min(df$Median) + (max(df$Median)-min(df$Median))/2 ), fill="#FFD1DC", show.legend=FALSE)+
    geom_line(data=df[df$Dist<0,],mapping=aes(x=Dist, y=Median), color="#FF02EB", linewidth=0.1)+
    geom_line(data=df[df$Dist>=0,],mapping=aes(x=Dist, y=Median), color="black", linewidth=0.1)+
    annotate("text", y=max(df[df$Dist<0,]$Median) +(realYmax-labelYmax), vjust=1, x=-2000+100, label="5’", color="#FF02EB", size=1.94)+
    annotate("text", y=max(df[df$Dist>=0,]$Median) +(realYmax-labelYmax), vjust=1, x=2000-100, label="3’", color="black", size=1.94)+
    geom_segment(x=0-as.numeric(cpsize), xend=0-as.numeric(cpsize), y=min(df$Median), yend=min(df$Median) + (max(df$Median)-min(df$Median))/2, color="#DC143C", linewidth=0.2, linetype="dashed") +
    geom_segment(x=0, xend=0, y=min(df$Median), yend=min(df$Median) + (max(df$Median)-min(df$Median))/2, color="#DC143C", linewidth=0.2, linetype="dashed") +
    geom_segment(x=a, xend=a, y=min(df$Median), yend=max(df$Median), color="#009FE9", linewidth=0.2, linetype="dashed") +
    geom_segment(x=b, xend=b, y=min(df$Median), yend=max(df$Median), color="#009FE9", linewidth=0.2, linetype="dashed") +
    annotate("text", y=min(df$Median)+((max(df$Median)-min(df$Median))/8.5*10), x=a+(b-a)/2, hjust=0.5, vjust=1, label = paste(format(b-a, big.mark = ",", scientific = FALSE), "bp"), color = "#009FE9", size = 1.94) +
    annotate("text", y=min(df$Median) + (max(df$Median)-min(df$Median))/2 + (realYmax-labelYmax), x=0-(as.numeric(cpsize)/2), hjust=0.5, vjust=1, label = paste(format(cpsize, big.mark = ",", scientific = FALSE), "bp"), color = "#DC143C", size = 1.94, fontface="bold") +
    scale_y_continuous(limits=c(min(df$Median), realYmax), breaks=c(min(df$Median),labelYmax), expand=c(0,0), name="\nMedian MP (%)")+
    scale_x_continuous(limits=c(-2000,2000), expand=c(0,0))+
    theme(plot.background=element_blank(), panel.background=element_blank(), panel.grid=element_blank(), axis.ticks.x=element_blank(),
          axis.line.y.right=element_blank(), axis.text.y.right=element_blank(), axis.title.y.right=element_blank(),
          axis.line.y.left=element_line(linewidth=0.1, color="black"), axis.text.y.left=element_text(size=5, color="black"), axis.title.y.left=element_text(size=5, color="black"), axis.ticks.y=element_blank(),
          axis.line.x=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank(),
          plot.margin=margin(2.5,0,0,0)
    )
}

# Human

human_svg=svgparser::read_svg("https://upload.wikimedia.org/wikipedia/commons/1/1d/No_image.svg")
# cp="140"; nGenome="1"; nTSS="77,453"; pstart=-874; pend=948; col="#000000"
cp="150"; nGenome="1"; nTSS="77,453"; pstart=-874; pend=948; col="#000000"
# load("T2Thuman_tvd_data_transcript.Rdata"); load("T2Thuman_pvalues_transcript.RData")
load("T2Thuman_10_tvd_data_transcript.Rdata"); load("T2Thuman_10_pvalues_transcript.RData")
T2Thuman_VSline=plot_VSline(T2Thuman_tvd_data, "Human (T2T-CHM13v2.0)", cp, nGenome, nTSS, human_svg)
T2Thuman_tile=plot_pvalue_heatmap(T2Thuman_pvalues)
T2Thuman_fakeVSline=plot_fakeVSline(T2Thuman_tvd_data)
T2Thuman_fakeTile=plot_fakePvalue_heatmap(T2Thuman_pvalues)
T2Thuman_diffLine=plot_diff(T2Thuman_tvd_data, col, cp)
T2Thuman_median=plotMedian("medianMP_transcriptTSS_T2Thuman.tsv", pstart, pend, col, cp)

T2Thuman_fakeNrealVSline = align_plots(T2Thuman_VSline + theme(plot.margin = margin(t=10, l=6, r=0, b=0)), T2Thuman_fakeVSline + theme(plot.margin = margin(t=10, l=0, r=0, b=0)), align="h", axis="tb")
T2Thuman_diffNtile = align_plots(T2Thuman_diffLine + theme(plot.margin = margin(t=10, b=0, l=0, r=0)), T2Thuman_tile + theme(plot.margin = margin(t=0, l=0, r=0, b=0)), align="v", axis="lr")
T2Thuman_soloMedian=align_plots(T2Thuman_median + theme(plot.margin = margin(t=10, l=6, r=7, b=0)))
col1=plot_grid(T2Thuman_fakeNrealVSline[[1]], T2Thuman_fakeNrealVSline[[2]], ncol=2, rel_widths=c(44.5, 2))
col2=plot_grid(T2Thuman_diffNtile[[1]], T2Thuman_diffNtile[[2]], ncol=1, rel_heights=c(90,10))
T2Thuman=plot_grid(col1, col2, T2Thuman_soloMedian[[1]], ncol=3, rel_widths=c(46.5, 44.5, 89), labels=c("d", "", "f"), label_size=8)
ggsave("T2Thuman.png", T2Thuman, height=17.85/105*99 * (111.8/90), width=180, units="mm") # height 17.85 mm corresponds to sum rel_heights of 105

# ZF

zf_svg=svgparser::read_svg("https://upload.wikimedia.org/wikipedia/commons/1/1d/No_image.svg")
# cp="180"; nGenome="1"; nTSS="41,039"; pstart=-1684; pend=1488; col="#000000"
cp="200"; nGenome="1"; nTSS="41,039"; pstart=-1684; pend=1488; col="#000000"
# load("bTaeGut7_tvd_data_transcript.Rdata"); load("bTaeGut7_pvalues_transcript.RData")
load("bTaeGut7_10_tvd_data_transcript.Rdata"); load("bTaeGut7_10_pvalues_transcript.RData")
bTaeGut7_VSline=plot_VSline(bTaeGut7_tvd_data, "Zebra finch (bTaeGut7)", cp, nGenome, nTSS, zf_svg)
bTaeGut7_tile=plot_pvalue_heatmap(bTaeGut7_pvalues)
bTaeGut7_fakeVSline=plot_fakeVSline(bTaeGut7_tvd_data)
bTaeGut7_fakeTile=plot_fakePvalue_heatmap(bTaeGut7_pvalues)
bTaeGut7_diffLine=plot_diff(bTaeGut7_tvd_data, col, cp)
bTaeGut7_median=plotMedian("medianMP_transcriptTSS_bTaeGut7.tsv", pstart, pend, col, cp)

bTaeGut7_fakeNrealVSline = align_plots(bTaeGut7_VSline + theme(plot.margin = margin(t=10, l=6, r=0, b=0)), bTaeGut7_fakeVSline + theme(plot.margin = margin(t=10, l=0, r=0, b=0)), align="h", axis="tb")
bTaeGut7_diffNtile = align_plots(bTaeGut7_diffLine + theme(plot.margin = margin(t=10, b=0, l=0, r=0)), bTaeGut7_tile + theme(plot.margin = margin(t=0, l=0, r=0, b=0)), align="v", axis="lr")
bTaeGut7_soloMedian=align_plots(bTaeGut7_median + theme(plot.margin = margin(t=10, l=6, r=7, b=0)))
col1=plot_grid(bTaeGut7_fakeNrealVSline[[1]], bTaeGut7_fakeNrealVSline[[2]], ncol=2, rel_widths=c(44.5, 2))
col2=plot_grid(bTaeGut7_diffNtile[[1]], bTaeGut7_diffNtile[[2]], ncol=1, rel_heights=c(90,10))
bTaeGut7=plot_grid(col1, col2, bTaeGut7_soloMedian[[1]], ncol=3, rel_widths=c(46.5, 44.5, 89), labels=c("e", "", "g"), label_size=8)
ggsave("bTaeGut7.png", bTaeGut7, height=17.85/105*99 * (111.8/90), width=180, units="mm") # height 17.85 mm corresponds to sum rel_heights of 105



### --- ###

## Legend provision

# For VSlines
emptyPlotMedian=function(){
  ggplot()+
    scale_x_continuous(limits=c(0,bound), breaks=c(0,bound), expand=c(0,0), name="Distance from TSS (bp)", labels = scales::label_comma())+
    theme_minimal() +
    theme(
      plot.background = element_blank(), 
      panel.background = element_blank(), 
      panel.grid = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      # Show x axis text and title
      axis.line.x=element_line(color="black", linewidth=0.1),
      axis.text.x = element_text(color = "black", size = 5),
      axis.title.x = element_text(color = "black", size = 5),
      # Show x axis ticks (needed for marking the positions of texts)
      axis.ticks.x = element_blank(),
      legend.position = "none",
      plot.margin=margin(0,0,0,0)
    )
}
legend=emptyPlotMedian()
legend_VSline = align_plots(T2Thuman_VSline + theme(plot.margin = margin(t=0, l=6, r=0, b=0)), legend + theme(plot.margin = margin(t=0, l=6, r=0, b=0)), align="v", axis="lr")
blank=ggplot() + theme_void()
last_row_1 = plot_grid(legend_VSline[[2]], blank, ncol=2, rel_widths=c(44.5, 2))

# For difflines
plot_pvalue_heatmap_minimal <- function(results_df, bin_col = "Bin", pvalue_col = "p_value", larger_dir_col = "Larger_Dir") {
  df <- results_df %>%
    mutate(
      fill_color = case_when(
        .data[[larger_dir_col]] == "5'" ~ 1 - .data[[pvalue_col]],
        .data[[larger_dir_col]] == "3'" ~ -(1 - .data[[pvalue_col]]),
        TRUE ~ NA_real_
      )
    )
  ggplot(df, aes(x = .data[[bin_col]], y = 1, fill = fill_color)) +
    # Remove geom_tile to hide the heatmap cells
    # geom_tile(color = "black") +
    
    # Remove color scale as coloring is irrelevant without tiles
    # scale_fill_gradientn(...) +
    
    theme_minimal() +
    theme(plot.background = element_blank(), panel.background = element_blank(), panel.grid = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_text(color = "black", size = 5), axis.title.x = element_text(color = "black", size = 5), axis.ticks.x = element_blank(),
          legend.position = "none"
    ) +
    labs(x = "Distance from TSS (bp)") +
    scale_x_continuous(labels = scales::label_comma(), limits=c(0,bound), breaks=c(0,bound), expand=c(0,0)) +
    scale_y_continuous(expand = c(0, 0))  # Optional: remove extra space on y-axis
}
legend=plot_pvalue_heatmap_minimal(T2Thuman_pvalues)
legend_tile = align_plots(T2Thuman_diffLine + theme(plot.margin = margin(b=0, l=0, r=0)), T2Thuman_tile + theme(plot.margin = margin(t=0, l=0, r=0, b=0)), legend + theme(plot.margin = margin(t=0, l=0, r=0, b=0)), align="v", axis="lr")
last_row_2=plot_grid(legend_tile[[3]])

# For medlines
emptyPlotMedian=function(){
  ggplot()+
    scale_x_continuous(limits=c(-2000,2000), breaks=c(-2000,0,2000), expand=c(0,0), name="Distance from TSS (bp)", labels = scales::label_comma())+
    theme_minimal() +
    theme(
      plot.background = element_blank(), 
      panel.background = element_blank(), 
      panel.grid = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      # Show x axis text and title
      axis.line.x=element_line(color="black", linewidth=0.1),
      axis.text.x = element_text(color = "black", size = 5),
      axis.title.x = element_text(color = "black", size = 5),
      # Show x axis ticks (needed for marking the positions of texts)
      axis.ticks.x = element_blank(),
      legend.position = "none",
      plot.margin=margin(0,0,0,0)
    )
}
emptyPlotMedian()->empty_med
meds=align_plots(T2Thuman_median + theme(plot.margin = margin(0,6,0,6)), empty_med + theme(plot.margin = margin(0,6,0,7)), align="v", axis="lr")

legends = plot_grid(last_row_1, last_row_2, meds[[2]], ncol=3, rel_widths=c(46.5, 44.5, 89))

# Colory legend provision

legend_data=data.frame("x"=c(-1,0,1),"y"=c(-1,0,1),"Location"=c("5’", "3’", "3’")); legend_data$"Location"=factor(legend_data$"Location", levels=c("5’", "3’"))
color_legend=as_ggplot(get_legend(ggplot(legend_data) + geom_line(mapping=aes(x=x,y=y,color=Location,group=Location)) +
                                   scale_color_manual(values=c("5’"="#FF02EB", "3’"="black"), name="Position respect to TSS:", guide = guide_legend(order = 1)) +
                                   new_scale_fill() +
                                   geom_tile(mapping=aes(x=x,y=y,fill=x)) +
                                   scale_fill_gradientn(colors=c("#FF02EB", "white", "black"), values = scales::rescale(c(-1, 0, 1)), breaks=c(-1, 0, 1), labels=c("0 (5’ < 3’)", "1", "0 (3’ < 5’)"), name="P-value\n\n\n", guide  = guide_colorbar(ticks = FALSE)) +
                                   theme(axis.ticks.x = element_line(color="black", linewidth=0.1), legend.position="top", legend.title=element_text(size=5, color="black"), legend.text=element_text(size=5, color="black"), legend.key=element_rect(fill="transparent", color=NA), legend.key.width = unit(6, "mm"), legend.key.height=unit(2,"mm"), legend.background=element_rect(fill="transparent"), legend.margin=margin(b=-10), legend.ticks = element_blank())
                                  ))

# Promoter legend provision

colors=c("Core promoter"="#DC143C","Promoter"="#009FE9")
fill_colors=c("Core promoter"="#FFD1DC","Promoter"="#C7E8FB")
as_ggplot(get_legend(ggplot(data.frame(x=c(1,2),col=c("Core promoter","Promoter")), aes(x = x, color = col, fill=col, linetype=col, linewidth=col)) +
  geom_density(alpha = 0.5) +
  scale_color_manual(values=colors,na.value="black")+
  scale_fill_manual(values=fill_colors,na.value="black")+
  scale_linetype_manual(values = c("Core promoter"="dashed", "Promoter"="dashed"))+
  scale_linewidth_manual(values = c("Core promoter"=.1, "Promoter"=.1))+
  theme(legend.position="top",
        legend.title=element_blank(),
        legend.text=element_text(size=5, color="black"),
        legend.key.size=unit(.5,"lines"),
        legend.key.spacing.x=unit(1,"lines"),
        legend.background=element_blank())
))->promoter_legend;promoter_legend

comp_legend = plot_grid(color_legend, promoter_legend, ncol=2, rel_widths=c(91, 89))


# Final arrangement

human_basics = align_plots(human_basic+theme(plot.margin=margin(10,10,0,0)), T2Thuman_VSline+theme(0,0,0,6), align="v", axis="l")
zf_basics=align_plots(zf_basic+theme(plot.margin=margin(10,10,0,0)), T2Thuman_diffLine+theme(0,0,0,0), align="v", axis="l")

rel_heights=c(180)
row=plot_grid(human_basics[[1]], zf_basics[[1]], blank, example, ncol=4, labels=c("a", "b", "", ""), rel_widths=c(44.5, 44.5, 2, 89), label_size=8)
ggsave("row.png", row, height=108.8/640*sum(rel_heights) * (111.8/90), width=180, units="mm")

rel_heights=c(180,99,99,20,70)
comp=plot_grid(row, T2Thuman, bTaeGut7, legends, comp_legend, ncol=1, rel_heights=rel_heights)
ggsave("comp.png", comp, height=108.8/640*sum(rel_heights) * (111.8/90), width=180, units="mm")
ggsave("comp.pdf", comp, height=108.8/640*sum(rel_heights) * (111.8/90), width=180, units="mm")







## Addt'l to fig. 5

df=fread("D:/Downloads/Supplementary Table 1_ Species data - Sheet1.tsv", dec = ".", fill = TRUE)
for(i in 1:length(df$Class)){
  if(df$Class[i]==""){
    df$Class[i]=df$Class[i-1]
  }
}
df$`Genome size (bp)`=as.numeric(gsub(",", "", df$`Genome size (bp)`))
df$`Promoter size (bp)`=as.numeric(gsub(",", "", df$`Promoter size (bp)`))
df$`Core promoter size (bp)`=as.numeric(df$`Core promoter size (bp)`)
df=df[df$`Promoter size (bp)`<10000,]
df=df[!is.na(df$`Core promoter size (bp)`),]
df=df[df$Outlier!="Yes",]

science=c("Mammalia"="#E69F00", "Aves"="#00796B", "Reptilia"="#A6D609", "Amphibia"="#984EA3", "Sarcopterygii"="#A6761D", "Actinopterygii"="#56B4E9", "Chondrichthyes"="#0072B2")
casual=c("Mammal"="#E69F00", "Bird"="#00796B", "Reptile"="#A6D609", "Amphibian"="#984EA3", "Lobe-finned fish"="#A6761D", "Ray-finned fish"="#56B4E9", "Shark"="#0072B2")
legend_data=data.frame("Class"=c("Mammal", "Bird", "Reptile", "Amphibian", "Lobe-finned fish", "Ray-finned fish", "Shark"), "X"=seq(1,7))
legend_data$Class = factor(legend_data$Class, levels=c("Mammal", "Bird", "Reptile", "Amphibian", "Lobe-finned fish", "Ray-finned fish", "Shark"))

ALPHA=0.14

ggplot(df) +
  coord_cartesian(clip = "off") +
  #stat_density_2d(data=df[df$Class=="Sarcopterygii",], mapping=aes(x=`Genome size (bp)`,y=`Promoter size (bp)`,fill=..level..), geom="polygon", alpha=ALPHA) + scale_fill_gradient(low="white", high="#A6761D")+ new_scale_fill() +
  stat_density_2d(data=df[df$Class=="Chondrichthyes",], mapping=aes(x=`Genome size (bp)`,y=`Promoter size (bp)`,fill=..level..), geom="polygon", alpha=ALPHA) + scale_fill_gradient(low="white", high="#0072B2")+ new_scale_fill() +
  stat_density_2d(data=df[df$Class=="Amphibia",], mapping=aes(x=`Genome size (bp)`,y=`Promoter size (bp)`,fill=..level..), geom="polygon", alpha=ALPHA) + scale_fill_gradient(low="white", high="#984EA3")+ new_scale_fill() +
  stat_density_2d(data=df[df$Class=="Reptilia",], mapping=aes(x=`Genome size (bp)`,y=`Promoter size (bp)`,fill=..level..), geom="polygon", alpha=ALPHA) + scale_fill_gradient(low="white", high="#A6D609")+ new_scale_fill() +
  stat_density_2d(data=df[df$Class=="Actinopterygii",], mapping=aes(x=`Genome size (bp)`,y=`Promoter size (bp)`,fill=..level..), geom="polygon", alpha=ALPHA) + scale_fill_gradient(low="white", high="#56B4E9")+ new_scale_fill() +
  stat_density_2d(data=df[df$Class=="Aves",], mapping=aes(x=`Genome size (bp)`,y=`Promoter size (bp)`,fill=..level..), geom="polygon", alpha=ALPHA) + scale_fill_gradient(low="white", high="#00796B")+ new_scale_fill() +
  stat_density_2d(data=df[df$Class=="Mammalia",], mapping=aes(x=`Genome size (bp)`,y=`Promoter size (bp)`,fill=..level..), geom="polygon", alpha=ALPHA) + scale_fill_gradient(low="white", high="#E69F00")+ new_scale_fill() +
  geom_point(mapping=aes(x=`Genome size (bp)`,y=`Promoter size (bp)`,color=`Class`),size=.3)+
  scale_y_continuous(name="Promoter size (kb)", labels=scales::label_number(scale=1e-3), expand=c(0,0), limits=c(0, max(df$`Promoter size (bp)`)))+
  scale_x_continuous(name="Genome size (Gb)", labels=scales::label_number(scale=1e-9), expand=c(0,0), limits=c(0, max(df$`Genome size (bp)`)))+
  scale_color_manual(values=science, na.value="pink")+
  theme(plot.background=element_blank(),
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.ticks=element_blank(),
        axis.line.y.left=element_line(linewidth=0.1, color="black"),
        axis.line.x.bottom = element_line(linewidth=0.1, color="black"),
        axis.title.y=element_text(size=5, color="black"),
        axis.text.y=element_text(size=5, color="black"),
        axis.title.x=element_text(size=5, color="black"),
        axis.text.x=element_text(size=5, color="black"),
        legend.position="none")->pg;pg

ggplot(df) +
  coord_cartesian(clip = "off") +
  #stat_density_2d(data=df[df$Class=="Sarcopterygii",], mapping=aes(x=`Genome size (bp)`,y=`Core promoter size (bp)`,fill=..level..), geom="polygon", alpha=ALPHA) + scale_fill_gradient(low="white", high="#A6761D")+ new_scale_fill() +
  stat_density_2d(data=df[df$Class=="Chondrichthyes",], mapping=aes(x=`Genome size (bp)`,y=`Core promoter size (bp)`,fill=..level..), geom="polygon", alpha=ALPHA) + scale_fill_gradient(low="white", high="#0072B2")+ new_scale_fill() +
  stat_density_2d(data=df[df$Class=="Amphibia",], mapping=aes(x=`Genome size (bp)`,y=`Core promoter size (bp)`,fill=..level..), geom="polygon", alpha=ALPHA) + scale_fill_gradient(low="white", high="#984EA3")+ new_scale_fill() +
  stat_density_2d(data=df[df$Class=="Reptilia",], mapping=aes(x=`Genome size (bp)`,y=`Core promoter size (bp)`,fill=..level..), geom="polygon", alpha=ALPHA) + scale_fill_gradient(low="white", high="#A6D609")+ new_scale_fill() +
  stat_density_2d(data=df[df$Class=="Actinopterygii",], mapping=aes(x=`Genome size (bp)`,y=`Core promoter size (bp)`,fill=..level..), geom="polygon", alpha=ALPHA) + scale_fill_gradient(low="white", high="#56B4E9")+ new_scale_fill() +
  stat_density_2d(data=df[df$Class=="Aves",], mapping=aes(x=`Genome size (bp)`,y=`Core promoter size (bp)`,fill=..level..), geom="polygon", alpha=ALPHA) + scale_fill_gradient(low="white", high="#00796B")+ new_scale_fill() +
  stat_density_2d(data=df[df$Class=="Mammalia",], mapping=aes(x=`Genome size (bp)`,y=`Core promoter size (bp)`,fill=..level..), geom="polygon", alpha=ALPHA) + scale_fill_gradient(low="white", high="#E69F00")+ new_scale_fill() +
  geom_point(mapping=aes(x=`Genome size (bp)`,y=`Core promoter size (bp)`,color=`Class`),size=.3)+
  scale_y_continuous(name="Core promoter size (bp)", expand=c(0,0), limits=c(0, max(df$`Core promoter size (bp)`)))+
  scale_x_continuous(name="Genome size (Gb)", labels=scales::label_number(scale=1e-9), expand=c(0,0), limits=c(0, max(df$`Genome size (bp)`)))+
  scale_color_manual(values=science, na.value="pink")+
  theme(plot.background=element_blank(),
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.ticks=element_blank(),
        axis.line.y.left=element_line(linewidth=0.1, color="black"),
        axis.line.x.bottom = element_line(linewidth=0.1, color="black"),
        axis.title.y=element_text(size=5, color="black"),
        axis.text.y=element_text(size=5, color="black"),
        axis.title.x=element_text(size=5, color="black"),
        axis.text.x=element_text(size=5, color="black"),
        legend.position="none")->cg;cg

ggplot(df) +
  coord_cartesian(clip = "off") +
  #stat_density_2d(data=df[df$Class=="Sarcopterygii",], mapping=aes(x=`Promoter size (bp)`,y=`Core promoter size (bp)`,fill=..level..), geom="polygon", alpha=ALPHA) + scale_fill_gradient(low="white", high="#A6761D")+ new_scale_fill() +
  stat_density_2d(data=df[df$Class=="Chondrichthyes",], mapping=aes(x=`Promoter size (bp)`,y=`Core promoter size (bp)`,fill=..level..), geom="polygon", alpha=ALPHA) + scale_fill_gradient(low="white", high="#0072B2")+ new_scale_fill() +
  stat_density_2d(data=df[df$Class=="Amphibia",], mapping=aes(x=`Promoter size (bp)`,y=`Core promoter size (bp)`,fill=..level..), geom="polygon", alpha=ALPHA) + scale_fill_gradient(low="white", high="#984EA3")+ new_scale_fill() +
  stat_density_2d(data=df[df$Class=="Reptilia",], mapping=aes(x=`Promoter size (bp)`,y=`Core promoter size (bp)`,fill=..level..), geom="polygon", alpha=ALPHA) + scale_fill_gradient(low="white", high="#A6D609")+ new_scale_fill() +
  stat_density_2d(data=df[df$Class=="Actinopterygii",], mapping=aes(x=`Promoter size (bp)`,y=`Core promoter size (bp)`,fill=..level..), geom="polygon", alpha=ALPHA) + scale_fill_gradient(low="white", high="#56B4E9")+ new_scale_fill() +
  stat_density_2d(data=df[df$Class=="Aves",], mapping=aes(x=`Promoter size (bp)`,y=`Core promoter size (bp)`,fill=..level..), geom="polygon", alpha=ALPHA) + scale_fill_gradient(low="white", high="#00796B")+ new_scale_fill() +
  stat_density_2d(data=df[df$Class=="Mammalia",], mapping=aes(x=`Promoter size (bp)`,y=`Core promoter size (bp)`,fill=..level..), geom="polygon", alpha=ALPHA) + scale_fill_gradient(low="white", high="#E69F00")+ new_scale_fill() +
  geom_point(mapping=aes(x=`Promoter size (bp)`,y=`Core promoter size (bp)`,color=`Class`),size=.3)+
  scale_y_continuous(name="Core promoter size (bp)", expand=c(0,0), limits=c(0, max(df$`Core promoter size (bp)`)))+
  scale_x_continuous(name="Promoter size (kb)", labels=scales::label_number(scale=1e-3), expand=c(0,0), limits=c(0, max(df$`Promoter size (bp)`)))+
  scale_color_manual(values=science, na.value="pink")+
  theme(plot.background=element_blank(),
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.ticks=element_blank(),
        axis.line.y.left=element_line(linewidth=0.1, color="black"),
        axis.line.x.bottom = element_line(linewidth=0.1, color="black"),
        axis.title.y=element_text(size=5, color="black"),
        axis.text.y=element_text(size=5, color="black"),
        axis.title.x=element_text(size=5, color="black"),
        axis.text.x=element_text(size=5, color="black"),
        legend.position="none")->cp;cp

genome_size_aligned=align_plots(pg + theme(plot.margin = margin(5,0,0,0)), cg + theme(plot.margin = margin(5,0,0,0)), align="v", axis="lr")
core_aligned=align_plots(cp + theme(plot.margin = margin(5,0,0,0)), cg + theme(plot.margin = margin(5,0,0,0)), align="h", axis="tb")

legend<-as_ggplot(get_legend(ggplot(legend_data)+geom_point(mapping=aes(color=Class,x=X,y=X))+
                               scale_color_manual(name="Class:", values=casual)+
                               guides(color = guide_legend(ncol = 1, override.aes = list(size = 2))) +
                               theme(plot.background = element_blank(),
                                     panel.background=element_blank(),
                                     legend.background=element_blank(),
                                     legend.text=element_text(color="black", size=5.5),
                                     legend.title=element_text(color="black", size=6, face="bold"),
                                     legend.key.size = unit(.75, "lines")
                               )))

plot_grid(genome_size_aligned[[1]], core_aligned[[2]], core_aligned[[1]], legend, ncol=4,
          labels=c("h", "i", "j"), label_size=8)->lengths


rel_heights=c(180,99,99,20,70,200)
compNlengths=plot_grid(comp, lengths, ncol=1, rel_heights=c(468,200))

ggsave("compNlengths.png", compNlengths, height=108.8/640*sum(rel_heights) * (111.8/90), width=180, units="mm")
ggsave("compNlengths.pdf", compNlengths, height=108.8/640*sum(rel_heights) * (111.8/90), width=180, units="mm")



















### --- ###

# Next picture which corresponds to the picture that has row pertaining to ea. clade.

plot_VSline <- function(both, species, cpsize, n, m, SVG, vjuster, hjuster){ # m is no. of TSSs
  #left=df[df$Dir=="5'"]; right=df[df$Dir=="3'"]
  #left=data.frame(summarise(group_by(left, Bin),"Mean"=mean(MP),"CI_lower"=ci(MP)[2],"CI_upper"=ci(MP)[3])); right=data.frame(summarise(group_by(right, Bin),"Mean"=mean(MP),"CI_lower"=ci(MP)[2],"CI_upper"=ci(MP)[3]))
  #both=data.frame(); both=rbind(both, cbind(left, "Dir"="5'")); both=rbind(both, cbind(right, "Dir"="3'"))
  #both=both[both$Bin<=bound,]
  ggplot(both)+
    #geom_rect(aes(xmin=0,xmax=as.numeric(cpsize),ymin=trunc(min(both$CI_lower)),ymax=trunc(max(both$CI_upper))+1 ), fill="#FFD1DC")+
    #geom_segment(x=as.numeric(cpsize), xend=as.numeric(cpsize), y=trunc(min(both$CI_lower)), yend=trunc(max(both$CI_upper))+1 , color="#DC143C", linewidth=0.2, linetype="dashed") +
    #annotation_raster(image, xmin=bound-190, xmax=bound, ymin=min(both$CI_lower), ymax=(max(both$CI_upper)-min(both$CI_lower))/2 + min(both$CI_lower)) + 
    geom_ribbon(data=both[both$Dir=="5'",], mapping=aes(y=Mean,x=Bin,ymin=CI_lower,ymax=CI_upper),color="transparent",fill="#FF02EB",alpha=0.25) +
    geom_ribbon(data=both[both$Dir=="3'",], mapping=aes(y=Mean,x=Bin,ymin=CI_lower,ymax=CI_upper),color="transparent",fill="black",alpha=0.25) +
    geom_line(mapping=aes(x=Bin, y=Mean, group=Dir, color=Dir),linewidth=.25)+
    annotation_custom(SVG, xmin = 1300, xmax = 2000, ymin = trunc(min(both$CI_lower)), ymax = trunc(min(both$CI_lower)) + (trunc(max(both$CI_upper))+1 - trunc(min(both$CI_lower))) * 0.66) +
    annotate("text", y=trunc(min(both$CI_lower)) + (trunc(max(both$CI_upper))+1 - trunc(min(both$CI_lower))) * 0.33, x=1650, vjust=vjuster, hjust=hjuster, label = n, color = "white", size = 1.94) +
    #annotate("text", y=trunc(max(both$CI_upper))+1, x=0, hjust=0, vjust=1, label = species, color = "black", size = 1.94, fontface="bold") +
    #annotate("text", y=trunc(max(both$CI_upper))+1, x=0, hjust=0, vjust=1, label = paste("\nnTSS =", m), color = "black", size = 1.94) +
    scale_x_continuous(expand=c(0,0))+
    scale_color_manual(name="Location respect to TSS", values=c("5'"="#FF02EB","3'"="black"),na.value="transparent", expand=c(0,0)) +
    scale_y_continuous(limits = c(trunc(min(both$CI_lower)),  trunc(max(both$CI_upper))+1), breaks = c(trunc(min(both$CI_lower)), trunc(max(both$CI_upper)) + 1), name = "\nMean MP (%)", expand=c(0,0))+
    theme(plot.background=element_blank(), panel.background = element_blank(), panel.grid = element_blank(), axis.ticks.x=element_blank(),
          axis.line.y.left = element_line(color="black", linewidth=0.1),    # Hide left y axis line
          axis.text.y.left = element_text(color="black", size=5),    # Hide left y axis texts
          axis.title.y.left = element_text(color="black", size=5),   # Hide left y axis title
          axis.ticks.y=element_line(color="black", linewidth=0.1),
          #axis.line.y.right = element_line(linewidth = 0.1, color = "black"), # Show right y axis line
          #axis.ticks.y.right = element_blank(),   # Show right y axis ticks
          #axis.text.y.right = element_text(size = 5, color = "black"), # Show right y axis labels
          #axis.title.y.right = element_text(size = 5, color = "black"), # Show right y axis title
          axis.title.x=element_blank(), axis.text.x=element_blank(),
          legend.position="none")
}

# Mammal

mammal_svg=svgparser::read_svg("svgs/mammal.svg")
plot_fakeVSline <- function(both){ # M is no. of CpGs, m is no. of TSSs
  #left=both[both$Dir=="5'"]; right=both[both$Dir=="3'"]
  #left=data.frame(summarise(group_by(left, Bin),"Mean"=mean(MP),"CI_lower"=ci(MP)[2],"CI_upper"=ci(MP)[3])); right=data.frame(summarise(group_by(right, Bin),"Mean"=mean(MP),"CI_lower"=ci(MP)[2],"CI_upper"=ci(MP)[3]))
  #both=data.frame(); both=rbind(both, cbind(left, "Dir"="5'")); both=rbind(both, cbind(right, "Dir"="3'"))
  #both=both[both$Bin<=bound,]
  ggplot(both)+
    annotate("text", y=both[both$Bin==bound&both$Dir=="5'",]$Mean, x=bound, hjust=.5, vjust=1, label="5’", color="#FF02EB", size=1.94)+
    annotate("text", y=both[both$Bin==bound&both$Dir=="3'",]$Mean, x=bound, hjust=.5, vjust=0, label="3’", color="black", size=1.94)+
    coord_cartesian(clip = "off") +
    scale_x_continuous(limits=c(bound,bound+80), expand=c(0,0))+
    scale_y_continuous(limits = c(trunc(min(both$CI_lower)), trunc(max(both$CI_upper)) + 1), breaks = c(trunc(min(both$CI_lower)), trunc(max(both$CI_upper)) + 1), name = "\nMean MP (%)", expand=c(0,0))+
    theme(plot.background=element_blank(), panel.background = element_blank(), panel.grid = element_blank(), axis.ticks=element_blank(),
          axis.title=element_blank(), axis.text=element_blank(),
          legend.position="none")
} # Same as fig. 7
plotMedian = function(filename, a, b, color, cpsize){
  df=fread(filename); colnames(df)=c("Dist", "Median"); df=df[df$Dist>=bound*-1&Dist<=bound,]
  realYmax=min(df$Median)+((max(df$Median)-min(df$Median))/8.5*10)
  labelYmax=max(df$Median)
  ggplot(df)+
    coord_cartesian(clip = "off") +
    geom_rect(aes(xmin=a, xmax=b, ymin=min(df$Median), ymax=max(df$Median)), fill="#C7E8FB", show.legend=FALSE)+
    geom_rect(aes(xmin=0-as.numeric(cpsize), xmax=0, ymin=min(df$Median), ymax=min(df$Median) + (max(df$Median)-min(df$Median))/2 ), fill="#FFD1DC", show.legend=FALSE)+
    geom_line(data=df,mapping=aes(x=Dist, y=Median), color=color, linewidth=0.1)+
    annotate("text", y=max(df[df$Dist<0,]$Median) +(realYmax-labelYmax), vjust=1, x=-2000+100, label="5’", color="black", size=1.94)+
    annotate("text", y=max(df[df$Dist>=0,]$Median) +(realYmax-labelYmax), vjust=1, x=2000-100, label="3’", color="black", size=1.94)+
    geom_segment(x=0-as.numeric(cpsize), xend=0-as.numeric(cpsize), y=min(df$Median), yend=min(df$Median) + (max(df$Median)-min(df$Median))/2, color="#DC143C", linewidth=0.2, linetype="dashed") +
    geom_segment(x=0, xend=0, y=min(df$Median), yend=min(df$Median) + (max(df$Median)-min(df$Median))/2, color="#DC143C", linewidth=0.2, linetype="dashed") +
    geom_segment(x=a, xend=a, y=min(df$Median), yend=max(df$Median), color="#009FE9", linewidth=0.2, linetype="dashed") +
    geom_segment(x=b, xend=b, y=min(df$Median), yend=max(df$Median), color="#009FE9", linewidth=0.2, linetype="dashed") +
    annotate("text", y=min(df$Median)+((max(df$Median)-min(df$Median))/8.5*10), x=a+(b-a)/2, hjust=0.5, vjust=1, label = paste(format(b-a, big.mark = ",", scientific = FALSE), "bp"), color = "#009FE9", size = 1.94) +
    annotate("text", y=min(df$Median) + (max(df$Median)-min(df$Median))/2 + (realYmax-labelYmax), x=0-(as.numeric(cpsize)/2), hjust=0.5, vjust=1, label = paste(format(cpsize, big.mark = ",", scientific = FALSE), "bp"), color = "#DC143C", size = 1.94, fontface="bold") +
    scale_y_continuous(limits=c(min(df$Median), realYmax), breaks=c(min(df$Median),labelYmax), expand=c(0,0), name="\nMedian MP (%)")+
    scale_x_continuous(limits=c(-2000,2000), expand=c(0,0))+
    theme(plot.background=element_blank(), panel.background=element_blank(), panel.grid=element_blank(), axis.ticks.x=element_blank(),
          axis.line.y.right=element_blank(), axis.text.y.right=element_blank(), axis.title.y.right=element_blank(),
          axis.line.y.left=element_line(linewidth=0.1, color="black"), axis.text.y.left=element_text(size=5, color="black"), axis.title.y.left=element_text(size=5, color="black"), axis.ticks.y=element_line(color="black", linewidth=0.1),
          axis.line.x=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank(),
          plot.margin=margin(2.5,0,0,0)
    )
} # 색깔 생김

name="Mammal"; cp="170"; nGenome="34"; nTSS="1,695,853"; pstart=-988; pend=1158; col="#E69F00"
# load("Mammal_tvd_data_transcript.Rdata"); load("Mammal_pvalues_transcript.RData"); medians="medianMP_transcriptTSS_Mammal.tsv"
load("Mammal_10_tvd_data_transcript.Rdata"); load("Mammal_10_pvalues_transcript.RData"); medians="medianMP_transcriptTSS_Mammal.tsv"
tvd_data = Mammal_tvd_data; pvalues = Mammal_pvalues
Mammal_VSline=plot_VSline(tvd_data, "", cp, nGenome, nTSS, mammal_svg, -.125, 0)
Mammal_tile=plot_pvalue_heatmap(pvalues)
Mammal_fakeVSline=plot_fakeVSline(tvd_data)
Mammal_fakeTile=plot_fakePvalue_heatmap(pvalues)
Mammal_diffLine=plot_diff(tvd_data, col, cp)
Mammal_median=plotMedian(medians, pstart, pend, col, cp)

# Bird

bird_svg=svgparser::read_svg("svgs/bird.svg")
plot_fakeVSline <- function(both){ # M is no. of CpGs, m is no. of TSSs
  #left=both[both$Dir=="5'"]; right=both[both$Dir=="3'"]
  #left=data.frame(summarise(group_by(left, Bin),"Mean"=mean(MP),"CI_lower"=ci(MP)[2],"CI_upper"=ci(MP)[3])); right=data.frame(summarise(group_by(right, Bin),"Mean"=mean(MP),"CI_lower"=ci(MP)[2],"CI_upper"=ci(MP)[3]))
  #both=data.frame(); both=rbind(both, cbind(left, "Dir"="5'")); both=rbind(both, cbind(right, "Dir"="3'"))
  #both=both[both$Bin<=bound,]
  ggplot(both)+
    annotate("text", y=both[both$Bin==bound&both$Dir=="5'",]$Mean, x=bound, hjust=.5, vjust=.5, label="5’", color="#FF02EB", size=1.94)+
    annotate("text", y=both[both$Bin==bound&both$Dir=="3'",]$Mean, x=bound, hjust=.5, vjust=.5, label="3’", color="black", size=1.94)+
    coord_cartesian(clip = "off") +
    scale_x_continuous(limits=c(bound,bound+80), expand=c(0,0))+
    scale_y_continuous(limits = c(trunc(min(both$CI_lower)), trunc(max(both$CI_upper)) + 1), breaks = c(trunc(min(both$CI_lower)), trunc(max(both$CI_upper)) + 1), name = "\nMean MP (%)", expand=c(0,0))+
    theme(plot.background=element_blank(), panel.background = element_blank(), panel.grid = element_blank(), axis.ticks=element_blank(),
          axis.title=element_blank(), axis.text=element_blank(),
          legend.position="none")
} # vjusted to 0.5 for both sides
plotMedian = function(filename, a, b, color, cpsize){
  df=fread(filename); colnames(df)=c("Dist", "Median"); df=df[df$Dist>=bound*-1&Dist<=bound,]
  realYmax=min(df$Median)+((max(df$Median)-min(df$Median))/8.5*10)
  labelYmax=max(df$Median)
  ggplot(df)+
    coord_cartesian(clip = "off") +
    geom_rect(aes(xmin=a, xmax=b, ymin=min(df$Median), ymax=max(df$Median)), fill="#C7E8FB", show.legend=FALSE)+
    geom_rect(aes(xmin=0-as.numeric(cpsize), xmax=0, ymin=min(df$Median), ymax=min(df$Median) + (max(df$Median)-min(df$Median))/2 ), fill="#FFD1DC", show.legend=FALSE)+
    geom_line(data=df,mapping=aes(x=Dist, y=Median), color=color, linewidth=0.1)+
    geom_segment(x=0-as.numeric(cpsize), xend=0-as.numeric(cpsize), y=min(df$Median), yend=min(df$Median) + (max(df$Median)-min(df$Median))/2, color="#DC143C", linewidth=0.2, linetype="dashed") +
    geom_segment(x=0, xend=0, y=min(df$Median), yend=min(df$Median) + (max(df$Median)-min(df$Median))/2, color="#DC143C", linewidth=0.2, linetype="dashed") +
    geom_segment(x=a, xend=a, y=min(df$Median), yend=max(df$Median), color="#009FE9", linewidth=0.2, linetype="dashed") +
    geom_segment(x=b, xend=b, y=min(df$Median), yend=max(df$Median), color="#009FE9", linewidth=0.2, linetype="dashed") +
    annotate("text", y=min(df$Median)+((max(df$Median)-min(df$Median))/8.5*10), x=a+(b-a)/2, hjust=0.5, vjust=1, label = paste(format(b-a, big.mark = ",", scientific = FALSE), "bp"), color = "#009FE9", size = 1.94) +
    annotate("text", y=min(df$Median) + (max(df$Median)-min(df$Median))/2 + (realYmax-labelYmax), x=0-(as.numeric(cpsize)/2), hjust=0.5, vjust=1, label = paste(format(cpsize, big.mark = ",", scientific = FALSE), "bp"), color = "#DC143C", size = 1.94, fontface="bold") +
    scale_y_continuous(limits=c(min(df$Median), realYmax), breaks=c(min(df$Median),labelYmax), expand=c(0,0), name="\nMedian MP (%)")+
    scale_x_continuous(limits=c(-2000,2000), expand=c(0,0))+
    theme(plot.background=element_blank(), panel.background=element_blank(), panel.grid=element_blank(), axis.ticks.x=element_blank(),
          axis.line.y.right=element_blank(), axis.text.y.right=element_blank(), axis.title.y.right=element_blank(),
          axis.line.y.left=element_line(linewidth=0.1, color="black"), axis.text.y.left=element_text(size=5, color="black"), axis.title.y.left=element_text(size=5, color="black"), axis.ticks.y=element_line(color="black", linewidth=0.1),
          axis.line.x=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank(),
          plot.margin=margin(2.5,0,0,0)
    )
} # 5' & 3' 표시 없어짐

name="Bird"; cp="310"; nGenome="18"; nTSS="561,713"; pstart=-1200; pend=1075; col="#00796B"
# load("Bird_tvd_data_transcript.Rdata"); load("Bird_pvalues_transcript.RData"); medians="medianMP_transcriptTSS_Bird.tsv"
load("Bird_10_tvd_data_transcript.Rdata"); load("Bird_10_pvalues_transcript.RData"); medians="medianMP_transcriptTSS_Bird.tsv"
# tvd_data = Mammal_tvd_data; pvalues = Mammal_pvalues # Was saved as "Mammal" srry.
tvd_data = Bird_tvd_data; pvalues = Bird_pvalues
Bird_VSline=plot_VSline(tvd_data, "", cp, nGenome, nTSS, bird_svg, .25, 0)
Bird_tile=plot_pvalue_heatmap(pvalues)
Bird_fakeVSline=plot_fakeVSline(tvd_data)
Bird_fakeTile=plot_fakePvalue_heatmap(pvalues)
Bird_diffLine=plot_diff(tvd_data, col, cp)
Bird_median=plotMedian(medians, pstart, pend, col, cp)

# Reptile

reptile_svg=svgparser::read_svg("svgs/reptile.svg")
plot_fakeVSline <- function(both){ # M is no. of CpGs, m is no. of TSSs
  #left=both[both$Dir=="5'"]; right=both[both$Dir=="3'"]
  #left=data.frame(summarise(group_by(left, Bin),"Mean"=mean(MP),"CI_lower"=ci(MP)[2],"CI_upper"=ci(MP)[3])); right=data.frame(summarise(group_by(right, Bin),"Mean"=mean(MP),"CI_lower"=ci(MP)[2],"CI_upper"=ci(MP)[3]))
  #both=data.frame(); both=rbind(both, cbind(left, "Dir"="5'")); both=rbind(both, cbind(right, "Dir"="3'"))
  #both=both[both$Bin<=bound,]
  ggplot(both)+
    annotate("text", y=both[both$Bin==bound&both$Dir=="5'",]$Mean, x=bound, hjust=.5, vjust=.5, label="5’", color="#FF02EB", size=1.94)+
    annotate("text", y=both[both$Bin==bound&both$Dir=="3'",]$Mean, x=bound, hjust=.5, vjust=0, label="3’", color="black", size=1.94)+
    coord_cartesian(clip = "off") +
    scale_x_continuous(limits=c(bound,bound+80), expand=c(0,0))+
    scale_y_continuous(limits = c(trunc(min(both$CI_lower)), trunc(max(both$CI_upper)) + 1), breaks = c(trunc(min(both$CI_lower)), trunc(max(both$CI_upper)) + 1), name = "\nMean MP (%)", expand=c(0,0))+
    theme(plot.background=element_blank(), panel.background = element_blank(), panel.grid = element_blank(), axis.ticks=element_blank(),
          axis.title=element_blank(), axis.text=element_blank(),
          legend.position="none")
}

name="Reptile"; cp="180"; nGenome="9"; nTSS="317,121"; pstart=-679; pend=816; col="#AEE009"
# load("Reptile_tvd_data_transcript.Rdata"); load("Reptile_pvalues_transcript.RData"); medians="medianMP_transcriptTSS_Reptile.tsv"
load("Reptile_10_tvd_data_transcript.Rdata"); load("Reptile_10_pvalues_transcript.RData"); medians="medianMP_transcriptTSS_Reptile.tsv"
# tvd_data = Mammal_tvd_data; pvalues = Mammal_pvalues # Was saved as "Mammal" srry.
tvd_data = Reptile_tvd_data; pvalues = Reptile_pvalues
Reptile_VSline=plot_VSline(tvd_data, "", cp, nGenome, nTSS, reptile_svg, .125, .5)
Reptile_tile=plot_pvalue_heatmap(pvalues)
Reptile_fakeVSline=plot_fakeVSline(tvd_data)
Reptile_fakeTile=plot_fakePvalue_heatmap(pvalues)
Reptile_diffLine=plot_diff(tvd_data, col, cp)
Reptile_median=plotMedian(medians, pstart, pend, col, cp)

# Amphibian

amphibian_svg=svgparser::read_svg("svgs/amphibian.svg")
plot_fakeVSline <- function(both){ # M is no. of CpGs, m is no. of TSSs
  #left=both[both$Dir=="5'"]; right=both[both$Dir=="3'"]
  #left=data.frame(summarise(group_by(left, Bin),"Mean"=mean(MP),"CI_lower"=ci(MP)[2],"CI_upper"=ci(MP)[3])); right=data.frame(summarise(group_by(right, Bin),"Mean"=mean(MP),"CI_lower"=ci(MP)[2],"CI_upper"=ci(MP)[3]))
  #both=data.frame(); both=rbind(both, cbind(left, "Dir"="5'")); both=rbind(both, cbind(right, "Dir"="3'"))
  #both=both[both$Bin<=bound,]
  ggplot(both)+
    annotate("text", y=both[both$Bin==bound&both$Dir=="5'",]$Mean, x=bound, hjust=.5, vjust=0, label="5’", color="#FF02EB", size=1.94)+
    annotate("text", y=both[both$Bin==bound&both$Dir=="3'",]$Mean, x=bound, hjust=.5, vjust=1, label="3’", color="black", size=1.94)+
    coord_cartesian(clip = "off") +
    scale_x_continuous(limits=c(bound,bound+80), expand=c(0,0))+
    scale_y_continuous(limits = c(trunc(min(both$CI_lower)), trunc(max(both$CI_upper)) + 1), breaks = c(trunc(min(both$CI_lower)), trunc(max(both$CI_upper)) + 1), name = "\nMean MP (%)", expand=c(0,0))+
    theme(plot.background=element_blank(), panel.background = element_blank(), panel.grid = element_blank(), axis.ticks=element_blank(),
          axis.title=element_blank(), axis.text=element_blank(),
          legend.position="none")
}
plotMedian = function(filename, a, b, color, cpsize){
  df=fread(filename); colnames(df)=c("Dist", "Median"); df=df[df$Dist>=bound*-1&Dist<=bound,]
  realYmax=min(df$Median)+((max(df$Median)-min(df$Median))/8.5*10)
  labelYmax=max(df$Median)
  ggplot(df)+
    coord_cartesian(clip = "off") +
    geom_rect(aes(xmin=a, xmax=b, ymin=min(df$Median), ymax=max(df$Median)), fill="#C7E8FB", show.legend=FALSE)+
    geom_rect(aes(xmin=0-as.numeric(cpsize), xmax=0, ymin=min(df$Median), ymax=min(df$Median) + (max(df$Median)-min(df$Median))/2 ), fill="#FFD1DC", show.legend=FALSE)+
    geom_line(data=df,mapping=aes(x=Dist, y=Median), color=color, linewidth=0.1)+
    geom_segment(x=0-as.numeric(cpsize), xend=0-as.numeric(cpsize), y=min(df$Median), yend=min(df$Median) + (max(df$Median)-min(df$Median))/2, color="#DC143C", linewidth=0.2, linetype="dashed") +
    geom_segment(x=0, xend=0, y=min(df$Median), yend=min(df$Median) + (max(df$Median)-min(df$Median))/2, color="#DC143C", linewidth=0.2, linetype="dashed") +
    geom_segment(x=a, xend=a, y=min(df$Median), yend=max(df$Median), color="#009FE9", linewidth=0.2, linetype="dashed") +
    geom_segment(x=b, xend=b, y=min(df$Median), yend=max(df$Median), color="#009FE9", linewidth=0.2, linetype="dashed") +
    annotate("text", y=min(df$Median)+((max(df$Median)-min(df$Median))/8.5*10), x=a+(b-a)/2, hjust=0.5, vjust=1, label = paste(format(b-a, big.mark = ",", scientific = FALSE), "bp"), color = "#009FE9", size = 1.94) +
    annotate("text", y=min(df$Median) + ((min(df$Median) + (max(df$Median)-min(df$Median))/2)-(min(df$Median)))/2, x=a, hjust=1.25, vjust=0.5, label = paste(format(cpsize, big.mark = ",", scientific = FALSE), "bp"), color = "#DC143C", size = 1.94, fontface="bold") +
    scale_y_continuous(limits=c(min(df$Median), realYmax), breaks=c(min(df$Median),labelYmax), expand=c(0,0), name="\nMedian MP (%)")+
    scale_x_continuous(limits=c(-2000,2000), expand=c(0,0))+
    theme(plot.background=element_blank(), panel.background=element_blank(), panel.grid=element_blank(), axis.ticks.x=element_blank(),
          axis.line.y.right=element_blank(), axis.text.y.right=element_blank(), axis.title.y.right=element_blank(),
          axis.line.y.left=element_line(linewidth=0.1, color="black"), axis.text.y.left=element_text(size=5, color="black"), axis.title.y.left=element_text(size=5, color="black"), axis.ticks.y=element_line(color="black", linewidth=0.1),
          axis.line.x=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank(),
          plot.margin=margin(2.5,0,0,0)
    )
} # 잠깐 core promoter 수치 왼쪽으로 옮김

name="Amphibian"; cp="80"; nGenome="6"; nTSS="387,075"; pstart=-200; pend=197; col="#984EA3"
# load("Amphibian_tvd_data_transcript.Rdata"); load("Amphibian_pvalues_transcript.RData"); medians="medianMP_transcriptTSS_Amphibian.tsv"
load("Amphibian_10_tvd_data_transcript.Rdata"); load("Amphibian_10_pvalues_transcript.RData"); medians="medianMP_transcriptTSS_Amphibian.tsv"
# tvd_data = Amphibian_tvd_data; pvalues = Amphibian_pvalues # 乃
tvd_data = Amphibian_tvd_data; pvalues = Amphibian_pvalues
Amphibian_VSline=plot_VSline(tvd_data, "", cp, nGenome, nTSS, amphibian_svg, -.25, .5)
Amphibian_tile=plot_pvalue_heatmap(pvalues)
Amphibian_fakeVSline=plot_fakeVSline(tvd_data)
Amphibian_fakeTile=plot_fakePvalue_heatmap(pvalues)
Amphibian_diffLine=plot_diff(tvd_data, col, cp)
Amphibian_median=plotMedian(medians, pstart, pend, col, cp)

# Coel

coel_svg=svgparser::read_svg("svgs/coel.svg")
plot_fakeVSline <- function(both){ # M is no. of CpGs, m is no. of TSSs
  #left=both[both$Dir=="5'"]; right=both[both$Dir=="3'"]
  #left=data.frame(summarise(group_by(left, Bin),"Mean"=mean(MP),"CI_lower"=ci(MP)[2],"CI_upper"=ci(MP)[3])); right=data.frame(summarise(group_by(right, Bin),"Mean"=mean(MP),"CI_lower"=ci(MP)[2],"CI_upper"=ci(MP)[3]))
  #both=data.frame(); both=rbind(both, cbind(left, "Dir"="5'")); both=rbind(both, cbind(right, "Dir"="3'"))
  #both=both[both$Bin<=bound,]
  ggplot(both)+
    annotate("text", y=both[both$Bin==bound&both$Dir=="5'",]$Mean, x=bound, hjust=.5, vjust=-.25, label="5’", color="#FF02EB", size=1.94)+
    annotate("text", y=both[both$Bin==bound&both$Dir=="3'",]$Mean, x=bound, hjust=.5, vjust=1, label="3’", color="black", size=1.94)+
    coord_cartesian(clip = "off") +
    scale_x_continuous(limits=c(bound,bound+80), expand=c(0,0))+
    scale_y_continuous(limits = c(trunc(min(both$CI_lower)), trunc(max(both$CI_upper)) + 1), breaks = c(trunc(min(both$CI_lower)), trunc(max(both$CI_upper)) + 1), name = "\nMean MP (%)", expand=c(0,0))+
    theme(plot.background=element_blank(), panel.background = element_blank(), panel.grid = element_blank(), axis.ticks=element_blank(),
          axis.title=element_blank(), axis.text=element_blank(),
          legend.position="none")
}
plotMedian = function(filename, a, b, color, cpsize){
  df=fread(filename); colnames(df)=c("Dist", "Median"); df=df[df$Dist>=bound*-1&Dist<=bound,]
  realYmax=min(df$Median)+((max(df$Median)-min(df$Median))/8.5*10)
  labelYmax=max(df$Median)
  ggplot(df)+
    coord_cartesian(clip = "off") +
    geom_rect(aes(xmin=a, xmax=b, ymin=min(df$Median), ymax=max(df$Median)), fill="#C7E8FB", show.legend=FALSE)+
    geom_rect(aes(xmin=0-as.numeric(cpsize), xmax=0, ymin=min(df$Median), ymax=min(df$Median) + (max(df$Median)-min(df$Median))/2 ), fill="#FFD1DC", show.legend=FALSE)+
    geom_line(data=df,mapping=aes(x=Dist, y=Median), color=color, linewidth=0.1)+
    geom_segment(x=0-as.numeric(cpsize), xend=0-as.numeric(cpsize), y=min(df$Median), yend=min(df$Median) + (max(df$Median)-min(df$Median))/2, color="#DC143C", linewidth=0.2, linetype="dashed") +
    geom_segment(x=0, xend=0, y=min(df$Median), yend=min(df$Median) + (max(df$Median)-min(df$Median))/2, color="#DC143C", linewidth=0.2, linetype="dashed") +
    geom_segment(x=a, xend=a, y=min(df$Median), yend=max(df$Median), color="#009FE9", linewidth=0.2, linetype="dashed") +
    geom_segment(x=b, xend=b, y=min(df$Median), yend=max(df$Median), color="#009FE9", linewidth=0.2, linetype="dashed") +
    annotate("text", y=min(df$Median)+((max(df$Median)-min(df$Median))/8.5*10), x=a+(b-a)/2, hjust=0.5, vjust=1, label = paste(format(b-a, big.mark = ",", scientific = FALSE), "bp"), color = "#009FE9", size = 1.94) +
    annotate("text", y=min(df$Median) + (max(df$Median)-min(df$Median))/2 + (realYmax-labelYmax), x=0-as.numeric(cpsize), hjust=0, vjust=1, label = paste(format(cpsize, big.mark = ",", scientific = FALSE), "bp"), color = "#DC143C", size = 1.94, fontface="bold") +
    scale_y_continuous(limits=c(min(df$Median), realYmax), breaks=c(min(df$Median),labelYmax), expand=c(0,0), name="\nMedian MP (%)")+
    scale_x_continuous(limits=c(-2000,2000), expand=c(0,0))+
    theme(plot.background=element_blank(), panel.background=element_blank(), panel.grid=element_blank(), axis.ticks.x=element_blank(),
          axis.line.y.right=element_blank(), axis.text.y.right=element_blank(), axis.title.y.right=element_blank(),
          axis.line.y.left=element_line(linewidth=0.1, color="black"), axis.text.y.left=element_text(size=5, color="black"), axis.title.y.left=element_text(size=5, color="black"), axis.ticks.y=element_line(color="black", linewidth=0.1),
          axis.line.x=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank(),
          plot.margin=margin(2.5,0,0,0)
    )
} # 잠깐 core promoter 수치 왼쪽 align 함

name="Lobe-finned fish"; cp="80"; nGenome="1"; nTSS="33,601"; pstart=-188; pend=469; col="#A6761D"
# load("Coel_tvd_data_transcript.Rdata"); load("Coel_pvalues_transcript.RData"); medians="medianMP_transcriptTSS_Coel.tsv"
load("Coel_10_tvd_data_transcript.Rdata"); load("Coel_10_pvalues_transcript.RData"); medians="medianMP_transcriptTSS_Coel.tsv"
# tvd_data = Coel_tvd_data; pvalues = Coel_pvalues # 奈
tvd_data = Coel_tvd_data; pvalues = Coel_pvalues
Coel_VSline=plot_VSline(tvd_data, "", cp, nGenome, nTSS, coel_svg, .5, -.5)
Coel_tile=plot_pvalue_heatmap(pvalues)
Coel_fakeVSline=plot_fakeVSline(tvd_data)
Coel_fakeTile=plot_fakePvalue_heatmap(pvalues)
Coel_diffLine=plot_diff(tvd_data, col, cp)
Coel_median=plotMedian(medians, pstart, pend, col, cp)

# Fish

fish_svg=svgparser::read_svg("svgs/fish.svg")
plot_fakeVSline <- function(both){ # M is no. of CpGs, m is no. of TSSs
  #left=both[both$Dir=="5'"]; right=both[both$Dir=="3'"]
  #left=data.frame(summarise(group_by(left, Bin),"Mean"=mean(MP),"CI_lower"=ci(MP)[2],"CI_upper"=ci(MP)[3])); right=data.frame(summarise(group_by(right, Bin),"Mean"=mean(MP),"CI_lower"=ci(MP)[2],"CI_upper"=ci(MP)[3]))
  #both=data.frame(); both=rbind(both, cbind(left, "Dir"="5'")); both=rbind(both, cbind(right, "Dir"="3'"))
  #both=both[both$Bin<=bound,]
  ggplot(both)+
    annotate("text", y=both[both$Bin==bound&both$Dir=="5'",]$Mean, x=bound, hjust=.5, vjust=.5, label="5’", color="#FF02EB", size=1.94)+
    annotate("text", y=both[both$Bin==bound&both$Dir=="3'",]$Mean, x=bound, hjust=.5, vjust=.5, label="3’", color="black", size=1.94)+
    coord_cartesian(clip = "off") +
    scale_x_continuous(limits=c(bound,bound+80), expand=c(0,0))+
    scale_y_continuous(limits = c(trunc(min(both$CI_lower)), trunc(max(both$CI_upper)) + 1), breaks = c(trunc(min(both$CI_lower)), trunc(max(both$CI_upper)) + 1), name = "\nMean MP (%)", expand=c(0,0))+
    theme(plot.background=element_blank(), panel.background = element_blank(), panel.grid = element_blank(), axis.ticks=element_blank(),
          axis.title=element_blank(), axis.text=element_blank(),
          legend.position="none")
} # vjusted to 0.5 for both sides
plotMedian = function(filename, a, b, color, cpsize){
  df=fread(filename); colnames(df)=c("Dist", "Median"); df=df[df$Dist>=bound*-1&Dist<=bound,]
  realYmax=min(df$Median)+((max(df$Median)-min(df$Median))/8.5*10)
  labelYmax=max(df$Median)
  ggplot(df)+
    coord_cartesian(clip = "off") +
    geom_rect(aes(xmin=a, xmax=b, ymin=min(df$Median), ymax=max(df$Median)), fill="#C7E8FB", show.legend=FALSE)+
    geom_rect(aes(xmin=0-as.numeric(cpsize), xmax=0, ymin=min(df$Median), ymax=min(df$Median) + (max(df$Median)-min(df$Median))/2 ), fill="#FFD1DC", show.legend=FALSE)+
    geom_line(data=df,mapping=aes(x=Dist, y=Median), color=color, linewidth=0.1)+
    geom_segment(x=0-as.numeric(cpsize), xend=0-as.numeric(cpsize), y=min(df$Median), yend=min(df$Median) + (max(df$Median)-min(df$Median))/2, color="#DC143C", linewidth=0.2, linetype="dashed") +
    geom_segment(x=0, xend=0, y=min(df$Median), yend=min(df$Median) + (max(df$Median)-min(df$Median))/2, color="#DC143C", linewidth=0.2, linetype="dashed") +
    geom_segment(x=a, xend=a, y=min(df$Median), yend=max(df$Median), color="#009FE9", linewidth=0.2, linetype="dashed") +
    geom_segment(x=b, xend=b, y=min(df$Median), yend=max(df$Median), color="#009FE9", linewidth=0.2, linetype="dashed") +
    annotate("text", y=min(df$Median)+((max(df$Median)-min(df$Median))/8.5*10), x=a+(b-a)/2, hjust=0.5, vjust=1, label = paste(format(b-a, big.mark = ",", scientific = FALSE), "bp"), color = "#009FE9", size = 1.94) +
    annotate("text", y=min(df$Median) + (max(df$Median)-min(df$Median))/2 + (realYmax-labelYmax), x=0-(as.numeric(cpsize)/2), hjust=0.5, vjust=1, label = paste(format(cpsize, big.mark = ",", scientific = FALSE), "bp"), color = "#DC143C", size = 1.94, fontface="bold") +
    scale_y_continuous(limits=c(min(df$Median), realYmax), breaks=c(min(df$Median),labelYmax), expand=c(0,0), name="\nMedian MP (%)")+
    scale_x_continuous(limits=c(-2000,2000), expand=c(0,0))+
    theme(plot.background=element_blank(), panel.background=element_blank(), panel.grid=element_blank(), axis.ticks.x=element_blank(),
          axis.line.y.right=element_blank(), axis.text.y.right=element_blank(), axis.title.y.right=element_blank(),
          axis.line.y.left=element_line(linewidth=0.1, color="black"), axis.text.y.left=element_text(size=5, color="black"), axis.title.y.left=element_text(size=5, color="black"), axis.ticks.y=element_line(color="black", linewidth=0.1),
          axis.line.x=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank(),
          plot.margin=margin(2.5,0,0,0)
    )
} # 원래 꺼 (5' & 3' 표시 없어짐까지)로

name="Ray-finned fish"; cp="80"; nGenome="16"; nTSS="726,731"; pstart=-256; pend=679; col="#56B4E9"
# load("Fish_tvd_data_transcript.Rdata"); load("Fish_pvalues_transcript.RData"); medians="medianMP_transcriptTSS_Fish.tsv"
load("Fish_10_tvd_data_transcript.Rdata"); load("Fish_10_pvalues_transcript.RData"); medians="medianMP_transcriptTSS_Fish.tsv"
# tvd_data = Mammal_tvd_data; pvalues = Mammal_pvalues # Was saved as Mammal. srry.
tvd_data = Fish_tvd_data; pvalues = Fish_pvalues
Fish_VSline=plot_VSline(tvd_data, "", cp, nGenome, nTSS, fish_svg, .5, 0)
Fish_tile=plot_pvalue_heatmap(pvalues)
Fish_fakeVSline=plot_fakeVSline(tvd_data)
Fish_fakeTile=plot_fakePvalue_heatmap(pvalues)
Fish_diffLine=plot_diff(tvd_data, col, cp)
Fish_median=plotMedian(medians, pstart, pend, col, cp)

# Shark

shark_svg=svgparser::read_svg("svgs/shark.svg")
plot_fakeVSline <- function(both){ # M is no. of CpGs, m is no. of TSSs
  #left=both[both$Dir=="5'"]; right=both[both$Dir=="3'"]
  #left=data.frame(summarise(group_by(left, Bin),"Mean"=mean(MP),"CI_lower"=ci(MP)[2],"CI_upper"=ci(MP)[3])); right=data.frame(summarise(group_by(right, Bin),"Mean"=mean(MP),"CI_lower"=ci(MP)[2],"CI_upper"=ci(MP)[3]))
  #both=data.frame(); both=rbind(both, cbind(left, "Dir"="5'")); both=rbind(both, cbind(right, "Dir"="3'"))
  #both=both[both$Bin<=bound,]
  ggplot(both)+
    annotate("text", y=both[both$Bin==bound&both$Dir=="5'",]$Mean, x=bound, hjust=.5, vjust=0, label="5’", color="#FF02EB", size=1.94)+
    annotate("text", y=both[both$Bin==bound&both$Dir=="3'",]$Mean, x=bound, hjust=.5, vjust=.5, label="3’", color="black", size=1.94)+
    coord_cartesian(clip = "off") +
    scale_x_continuous(limits=c(bound,bound+80), expand=c(0,0))+
    scale_y_continuous(limits = c(trunc(min(both$CI_lower)), trunc(max(both$CI_upper)) + 1), breaks = c(trunc(min(both$CI_lower)), trunc(max(both$CI_upper)) + 1), name = "\nMean MP (%)", expand=c(0,0))+
    theme(plot.background=element_blank(), panel.background = element_blank(), panel.grid = element_blank(), axis.ticks=element_blank(),
          axis.title=element_blank(), axis.text=element_blank(),
          legend.position="none")
}

name="Cartilaginous fish"; cp="80"; nGenome="4"; nTSS="181,281"; pstart=-440; pend=705; col="#0072B2"
# load("Shark_tvd_data_transcript.Rdata"); load("Shark_pvalues_transcript.RData"); medians="medianMP_transcriptTSS_Shark.tsv"
load("Shark_10_tvd_data_transcript.Rdata"); load("Shark_10_pvalues_transcript.RData"); medians="medianMP_transcriptTSS_Shark.tsv"
# tvd_data = Shark_tvd_data; pvalues = Shark_pvalues # 果
tvd_data = Shark_tvd_data; pvalues = Shark_pvalues
Shark_VSline=plot_VSline(tvd_data, "", cp, nGenome, nTSS, shark_svg, .5, -1)
Shark_tile=plot_pvalue_heatmap(pvalues)
Shark_fakeVSline=plot_fakeVSline(tvd_data)
Shark_fakeTile=plot_fakePvalue_heatmap(pvalues)
Shark_diffLine=plot_diff(tvd_data, col, cp)
Shark_median=plotMedian(medians, pstart, pend, col, cp)

# Legend provision (for alignment only)

emptyPlotMedian=function(){
  ggplot()+
    scale_x_continuous(limits=c(0,bound), breaks=c(0,bound), expand=c(0,0), name="Distance from TSS (bp)", labels = scales::label_comma())+
    theme_minimal() +
    theme(
      plot.background = element_blank(), 
      panel.background = element_blank(), 
      panel.grid = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      # Show x axis text and title
      axis.line.x=element_line(color="black", linewidth=0.1),
      axis.text.x = element_text(color = "black", size = 5),
      axis.title.x = element_text(color = "black", size = 5),
      # Show x axis ticks (needed for marking the positions of texts)
      axis.ticks.x = element_line(color="black", linewidth=0.1),
      legend.position = "none",
      plot.margin=margin(0,0,0,0)
    )
}
legend_fakeNrealVSline=emptyPlotMedian()

plot_pvalue_heatmap_minimal <- function(results_df, bin_col = "Bin", pvalue_col = "p_value", larger_dir_col = "Larger_Dir") {
  df <- results_df %>%
    mutate(
      fill_color = case_when(
        .data[[larger_dir_col]] == "5'" ~ 1 - .data[[pvalue_col]],
        .data[[larger_dir_col]] == "3'" ~ -(1 - .data[[pvalue_col]]),
        TRUE ~ NA_real_
      )
    )
  ggplot(df, aes(x = .data[[bin_col]], y = 1, fill = fill_color)) +
    # Remove geom_tile to hide the heatmap cells
    # geom_tile(color = "black") +
    
    # Remove color scale as coloring is irrelevant without tiles
    # scale_fill_gradientn(...) +
    
    theme_minimal() +
    theme(plot.background = element_blank(), panel.background = element_blank(), panel.grid = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_text(color = "black", size = 5), axis.title.x = element_text(color = "black", size = 5), axis.ticks.x = element_line(color="black", linewidth=0.1),
          legend.position = "none"
    ) +
    labs(x = "Distance from TSS (bp)") +
    scale_x_continuous(labels = scales::label_comma(), limits=c(0,bound), breaks=c(0,bound), expand=c(0,0)) +
    scale_y_continuous(expand = c(0, 0))  # Optional: remove extra space on y-axis
}
legend_diffNtile=plot_pvalue_heatmap_minimal(T2Thuman_pvalues)

emptyPlotMedian=function(){
  ggplot()+
    scale_x_continuous(limits=c(-2000,2000), breaks=c(-2000,0,2000), expand=c(0,0), name="Distance from TSS (bp)", labels = scales::label_comma())+
    theme_minimal() +
    theme(
      plot.background = element_blank(), 
      panel.background = element_blank(), 
      panel.grid = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      # Show x axis text and title
      axis.line.x=element_line(color="black", linewidth=0.1),
      axis.text.x = element_text(color = "black", size = 5),
      axis.title.x = element_text(color = "black", size = 5),
      # Show x axis ticks (needed for marking the positions of texts)
      axis.ticks.x = element_line(color="black", linewidth=0.1),
      legend.position = "none",
      plot.margin=margin(0,0,0,0)
    )
}
emptyPlotMedian()->legend_median

### --- ###

blank = ggplot()+theme_void()

all_VSlines=align_plots(Mammal_VSline + theme(plot.margin = margin(t=10, l=6, r=0, b=0)),
                                Bird_VSline + theme(plot.margin = margin(t=10, l=6, r=0, b=0)),
                                Reptile_VSline + theme(plot.margin = margin(t=10, l=6, r=0, b=0)),
                                Amphibian_VSline + theme(plot.margin = margin(t=10, l=6, r=0, b=0)),
                                Coel_VSline + theme(plot.margin = margin(t=10, l=6, r=0, b=0)),
                                Fish_VSline + theme(plot.margin = margin(t=10, l=6, r=0, b=0)),
                                Shark_VSline + theme(plot.margin = margin(t=10, l=6, r=0, b=0)),
                                legend_fakeNrealVSline + theme(plot.margin = margin(t=0, l=6, r=0, b=0)),
                                align="v", axis="lr")
all_fakeVSlines=align_plots(Mammal_fakeVSline + theme(plot.margin = margin(t=10, l=0, r=0, b=0)),
                            Bird_fakeVSline + theme(plot.margin = margin(t=10, l=0, r=0, b=0)),
                            Reptile_fakeVSline + theme(plot.margin = margin(t=10, l=0, r=0, b=0)),
                            Amphibian_fakeVSline + theme(plot.margin = margin(t=10, l=0, r=0, b=0)),
                            Coel_fakeVSline + theme(plot.margin = margin(t=10, l=0, r=0, b=0)),
                            Fish_fakeVSline + theme(plot.margin = margin(t=10, l=0, r=0, b=0)),
                            Shark_fakeVSline + theme(plot.margin = margin(t=10, l=0, r=0, b=0)),
                            blank + theme(plot.margin = margin(t=10, l=0, r=0, b=0)),
                            align="v", axis="lr")
all_diffNtiles=align_plots(Mammal_diffLine + theme(plot.margin = margin(t=10, b=0, l=0, r=0)), Mammal_tile + theme(plot.margin = margin(t=0, l=0, r=0, b=0)),
                     Bird_diffLine + theme(plot.margin = margin(t=10, b=0, l=0, r=0)), Bird_tile + theme(plot.margin = margin(t=0, l=0, r=0, b=0)),
                     Reptile_diffLine + theme(plot.margin = margin(t=10, b=0, l=0, r=0)), Reptile_tile + theme(plot.margin = margin(t=0, l=0, r=0, b=0)),
                     Amphibian_diffLine + theme(plot.margin = margin(t=10, b=0, l=0, r=0)), Amphibian_tile + theme(plot.margin = margin(t=0, l=0, r=0, b=0)),
                     Coel_diffLine + theme(plot.margin = margin(t=10, b=0, l=0, r=0)), Coel_tile + theme(plot.margin = margin(t=0, l=0, r=0, b=0)),
                     Fish_diffLine + theme(plot.margin = margin(t=10, b=0, l=0, r=0)), Fish_tile + theme(plot.margin = margin(t=0, l=0, r=0, b=0)),
                     Shark_diffLine + theme(plot.margin = margin(t=10, b=0, l=0, r=0)), Shark_tile + theme(plot.margin = margin(t=0, l=0, r=0, b=0)),
                     legend_diffNtile + theme(plot.margin = margin(t=0, l=0, r=0, b=0)),
                     align="v", axis="lr")
all_soloMedians=align_plots(Mammal_median + theme(plot.margin = margin(t=10, l=6, r=7, b=0)),
                            Bird_median + theme(plot.margin = margin(t=10, l=6, r=7, b=0)),
                            Reptile_median + theme(plot.margin = margin(t=10, l=6, r=7, b=0)),
                            Amphibian_median + theme(plot.margin = margin(t=10, l=6, r=7, b=0)),
                            Coel_median + theme(plot.margin = margin(t=10, l=6, r=7, b=0)),
                            Fish_median + theme(plot.margin = margin(t=10, l=6, r=7, b=0)),
                            Shark_median + theme(plot.margin = margin(t=10, l=6, r=7, b=0)),
                            legend_median + theme(plot.margin = margin(0,6,0,7)),
                            align="v", axis="lr")

col1=plot_grid(all_VSlines[[1]], all_fakeVSlines[[1]], ncol=2, rel_widths=c(44.5, 2))
col2=plot_grid(all_diffNtiles[[1]], all_diffNtiles[[2]], ncol=1, rel_heights=c(90,10))
Mammal=plot_grid(col1, col2, all_soloMedians[[1]], ncol=3, rel_widths=c(46.5, 44.5, 89))
ggsave("Mammal.png", Mammal, height=17.85/105*99 * (111.8/90), width=180, units="mm") # height 17.85 mm corresponds to sum rel_heights of 105

col1=plot_grid(all_VSlines[[2]], all_fakeVSlines[[2]], ncol=2, rel_widths=c(44.5, 2))
col2=plot_grid(all_diffNtiles[[3]], all_diffNtiles[[4]], ncol=1, rel_heights=c(90,10))
Bird=plot_grid(col1, col2, all_soloMedians[[2]], ncol=3, rel_widths=c(46.5, 44.5, 89))
ggsave("Bird.png", Bird, height=17.85/105*99 * (111.8/90), width=180, units="mm") # height 17.85 mm corresponds to sum rel_heights of 105

col1=plot_grid(all_VSlines[[3]], all_fakeVSlines[[3]], ncol=2, rel_widths=c(44.5, 2))
col2=plot_grid(all_diffNtiles[[5]], all_diffNtiles[[6]], ncol=1, rel_heights=c(90,10))
Reptile=plot_grid(col1, col2, all_soloMedians[[3]], ncol=3, rel_widths=c(46.5, 44.5, 89))
ggsave("Reptile.png", Reptile, height=17.85/105*99 * (111.8/90), width=180, units="mm") # height 17.85 mm corresponds to sum rel_heights of 105

col1=plot_grid(all_VSlines[[4]], all_fakeVSlines[[4]], ncol=2, rel_widths=c(44.5, 2))
col2=plot_grid(all_diffNtiles[[7]], all_diffNtiles[[8]], ncol=1, rel_heights=c(90,10))
Amphibian=plot_grid(col1, col2, all_soloMedians[[4]], ncol=3, rel_widths=c(46.5, 44.5, 89))
ggsave("Amphibian.png", Amphibian, height=17.85/105*99 * (111.8/90), width=180, units="mm") # height 17.85 mm corresponds to sum rel_heights of 105

col1=plot_grid(all_VSlines[[5]], all_fakeVSlines[[5]], ncol=2, rel_widths=c(44.5, 2))
col2=plot_grid(all_diffNtiles[[9]], all_diffNtiles[[10]], ncol=1, rel_heights=c(90,10))
Coel=plot_grid(col1, col2, all_soloMedians[[5]], ncol=3, rel_widths=c(46.5, 44.5, 89))
ggsave("Coel.png", Coel, height=17.85/105*99 * (111.8/90), width=180, units="mm") # height 17.85 mm corresponds to sum rel_heights of 105

col1=plot_grid(all_VSlines[[6]], all_fakeVSlines[[6]], ncol=2, rel_widths=c(44.5, 2))
col2=plot_grid(all_diffNtiles[[11]], all_diffNtiles[[12]], ncol=1, rel_heights=c(90,10))
Fish=plot_grid(col1, col2, all_soloMedians[[6]], ncol=3, rel_widths=c(46.5, 44.5, 89))
ggsave("Fish.png", Fish, height=17.85/105*99 * (111.8/90), width=180, units="mm") # height 17.85 mm corresponds to sum rel_heights of 105

col1=plot_grid(all_VSlines[[7]], all_fakeVSlines[[7]], ncol=2, rel_widths=c(44.5, 2))
col2=plot_grid(all_diffNtiles[[13]], all_diffNtiles[[14]], ncol=1, rel_heights=c(90,10))
Shark=plot_grid(col1, col2, all_soloMedians[[7]], ncol=3, rel_widths=c(46.5, 44.5, 89))
ggsave("Shark.png", Shark, height=17.85/105*99 * (111.8/90), width=180, units="mm") # height 17.85 mm corresponds to sum rel_heights of 105

# Legend (Plot)

legends = plot_grid(all_VSlines[[8]], all_fakeVSlines[[8]], all_diffNtiles[[15]], all_soloMedians[[8]], ncol=4, rel_widths=c(44.5, 2, 44.5, 89))

# Colory legend provision

legend_data=data.frame("x"=c(-1,0,1),"y"=c(-1,0,1),"Location"=c("5’", "3’", "3’")); legend_data$"Location"=factor(legend_data$"Location", levels=c("5’", "3’"))
color_legend=as_ggplot(get_legend(ggplot(legend_data) + geom_line(mapping=aes(x=x,y=y,color=Location,group=Location)) +
                                    scale_color_manual(values=c("5’"="#FF02EB", "3’"="black"), name="Position respect to TSS:", guide = guide_legend(order = 1)) +
                                    new_scale_fill() +
                                    geom_tile(mapping=aes(x=x,y=y,fill=x)) +
                                    scale_fill_gradientn(colors=c("#FF02EB", "white", "black"), values = scales::rescale(c(-1, 0, 1)), breaks=c(-1, 0, 1), labels=c("0 (5’ < 3’)", "1", "0 (3’ < 5’)"), name="P-value\n\n\n") +
                                    theme(axis.ticks.x = element_line(color="black", linewidth=0.1), legend.position="top", legend.title=element_text(size=5, color="black"), legend.text=element_text(size=5, color="black"), legend.key=element_rect(fill="transparent", color=NA), legend.key.width = unit(6, "mm"), legend.key.height=unit(2,"mm"), legend.background=element_rect(fill="transparent"), legend.margin=margin(b=-10))
))





# Final arrangement (Clades)

rel_heights=c(99,99,99,99,99,99,99,20,70)
clades=plot_grid(Mammal, Bird, Reptile, Amphibian, Coel, Fish, Shark, legends, color_legend, ncol=1, rel_heights=rel_heights, labels=c("a", "b", "c", "d", "e", "f", "g"), label_size=8)
ggsave("clades.png", clades, height=108.8/640*sum(rel_heights) * (111.8/90), width=164.5, units="mm")

# Rotated name addition

Mammal_text = ggdraw() + draw_grob(textGrob("Mammals", rot=270, gp=gpar(fontsize=6, fontface="bold", col="#E69F00")))
Bird_text = ggdraw() + draw_grob(textGrob("Birds", rot=270, gp=gpar(fontsize=6, fontface="bold", col="#00796B")))
Reptile_text = ggdraw() + draw_grob(textGrob("Reptiles", rot=270, gp=gpar(fontsize=6, fontface="bold", col="#AEE009")))
Amphibian_text = ggdraw() + draw_grob(textGrob("Amphibians", rot=270, gp=gpar(fontsize=6, fontface="bold", col="#984EA3")))
Coel_text = ggdraw() + draw_grob(textGrob("Lobe-finned fishes", rot=270, gp=gpar(fontsize=6, fontface="bold", col="#A6761D")))
Fish_text = ggdraw() + draw_grob(textGrob("Ray-finned fishes", rot=270, gp=gpar(fontsize=6, fontface="bold", col="#56B4E9")))
Shark_text = ggdraw() + draw_grob(textGrob("Cartilaginous fishes", rot=270, gp=gpar(fontsize=6, fontface="bold", col="#0072B2")))

Mammal_texted = plot_grid(Mammal, Mammal_text, ncol=2, rel_widths=c(164.5, 5.5))
Bird_texted = plot_grid(Bird, Bird_text, ncol=2, rel_widths=c(164.5, 5.5))
Reptile_texted = plot_grid(Reptile, Reptile_text, ncol=2, rel_widths=c(164.5, 5.5))
Amphibian_texted = plot_grid(Amphibian, Amphibian_text, ncol=2, rel_widths=c(164.5, 5.5))
Coel_texted = plot_grid(Coel, Coel_text, ncol=2, rel_widths=c(164.5, 5.5))
Fish_texted = plot_grid(Fish, Fish_text, ncol=2, rel_widths=c(164.5, 5.5))
Shark_texted = plot_grid(Shark, Shark_text, ncol=2, rel_widths=c(164.5, 5.5))
legends_texted = plot_grid(legends, blank, ncol=2, rel_widths=c(164.5, 5.5))
clades_texted = plot_grid(Mammal_texted, Bird_texted, Reptile_texted, Amphibian_texted, Coel_texted, Fish_texted, Shark_texted, legends_texted, color_legend, ncol=1, rel_heights=rel_heights, labels=c("a", "b", "c", "d", "e", "f", "g"), label_size=8)

# Dendogram addition

load("plannarized_dendogram.RData")
more_plannarized_dendogram = plot_grid(plannarized_dendogram, blank, ncol=1, rel_heights=c(7*99, 20+70))
clades_dendogram = plot_grid(clades_texted, more_plannarized_dendogram, ncol=2, rel_widths=c(170, 10))
ggsave("clades_dendogram.pdf", clades_dendogram, height=108.8/640*sum(rel_heights) * (111.8/90), width=180, units="mm")
ggsave("clades_dendogram.png", clades_dendogram, height=108.8/640*sum(rel_heights) * (111.8/90), width=180, units="mm")

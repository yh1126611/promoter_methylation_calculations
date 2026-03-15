class_color <- c("Human"="#000000", "Mammalia"="#E69F00","Aves"="#00796B","Reptilia"="#A6D609",
                 "Amphibia"="#984EA3","Sarcopterygii"="#A6761D","Actinopterygii"="#56B4E9","Chondrichthyes"="#0072B2")
class_order <- c("Human", "Mammalia", "Aves", "Reptilia", "Amphibia", "Sarcopterygii", "Actinopterygii", "Chondrichthyes")

plot_silhouettes <- function(df, lengths_df, START,END){
  maximilian = max(df$sil_width)
  df=df[START<=df$Xposition & df$Xposition<=END,]
  lengths_df=lengths_df[START<lengths_df$start & lengths_df$end<END,]
  lengths_df$start=lengths_df$start-START+1
  lengths_df$end=lengths_df$end-START+1
  ggplot(df)+
    geom_bar(stat="identity", aes(x=reorder(Number, Xposition), y=sil_width, fill=class_fct)) +
    annotate("segment", linewidth=.25, linetype="dashed", color="black", x=lengths_df[-1,]$start-1, xend=lengths_df[-1,]$start-1, y=maximilian, yend=0) +
    annotate("text", size=1.94, color="black", x=lengths_df$start+(lengths_df$end-lengths_df$start)/2, hjust=0.5, y=maximilian*1.1, vjust=.5, label=lengths_df$cluster)+
    scale_y_continuous(name="Silhouette width (Si)", labels = scales::number_format(accuracy = 0.01), breaks=c(0, maximilian), limits=c(0, maximilian*1.25), expand=c(0,0))+
    scale_fill_manual(values=class_color, na.value="lightgrey")+
    theme(legend.position="none",
          plot.background=element_blank(), panel.background=element_blank(), panel.grid=element_blank(),
          axis.ticks=element_blank(), axis.line = element_line(linewidth = 0.25, lineend="square"),
          axis.title.y=element_text(color="black", size=5), axis.text.y=element_text(color="black", size=5), 
          axis.title.x=element_blank(), axis.text.x=element_text(color="black", size=5),
          plot.margin = unit(c(2.5,0,0,0), "pt")
    )
}


plot_UMAP<-function(umap_df, sil_df, TITLE, COLORBY, COLORBY_CHART){
  umap_df$cluster = sil_df$cluster[match(umap_df$Assembly, sil_df$Assembly)]
  ggplot(umap_df)+
    coord_fixed(ratio=(max(umap_df$UMAP1)-min(umap_df$UMAP1)) / (max(umap_df$UMAP2)-min(umap_df$UMAP2)))+
    geom_point(mapping=aes(x=UMAP1, y=UMAP2, color=COLORBY), size=0.5)+
    geom_point(data=(umap_df[umap_df$Class=="Human",]), mapping=aes(x=UMAP1, y=UMAP2, color=Class), size=0.5)+
    scale_color_manual(values=COLORBY_CHART) +
    scale_x_continuous(limits=c(min(umap_df$UMAP1), max(umap_df$UMAP1)))+
    scale_y_continuous(limits=c(min(umap_df$UMAP2), max(umap_df$UMAP2)))+
    labs(x = "UMAP1", y = "UMAP2", title = TITLE) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 5.5),
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.line.y.left = element_line(),
      axis.line.x.bottom = element_line(),
      panel.grid = element_blank(),
      axis.line = element_line(linewidth = 0.25, lineend="square"),
      axis.ticks = element_blank(),
      axis.title = element_text(size = 5.5),
      axis.text = element_text(size = 5, color = "black"),
      plot.margin = unit(c(5.5,5.5,0,0), "pt")
    )
}

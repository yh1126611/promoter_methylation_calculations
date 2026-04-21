library(ggplot2)

# A vector for color coding scheme must be defined, e.g.
class_color <- c(
  "Human"="#000000", "Mammalia"="#E69F00","Aves"="#00796B","Reptilia"="#A6D609",
  "Amphibia"="#984EA3","Actinopterygii"="#56B4E9","Sarcopterygii"="#A6761D","Chondrichthyes"="#0072B2"
)

plot_UMAP<-function(umap_df, TITLE, COLORBY, COLORBY_CHART, FACE="bold"){
ggplot(umap_df)+
    coord_fixed(ratio= (max(umap_df$UMAP1)-min(umap_df$UMAP1)) / (max(umap_df$UMAP2)-min(umap_df$UMAP2)) )+
  geom_point(mapping=aes(x=UMAP1, y=UMAP2, color=COLORBY), size=0.5)+
  geom_point(data=(umap_df[umap_df$Class=="Human",]), mapping=aes(x=UMAP1, y=UMAP2, color=Class), size=0.5)+
  scale_color_manual(values=COLORBY_CHART, na.value="#CCCCCC") +
    scale_x_continuous(limits=c(min(umap_df$UMAP1), max(umap_df$UMAP1)))+
    scale_y_continuous(limits=c(min(umap_df$UMAP2), max(umap_df$UMAP2)), labels=label_number(accuracy = 0.1))+
  labs(x = "UMAP1", y = "UMAP2", title = TITLE) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 5.5, face=FACE),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y.left = element_line(),
    axis.line.x.bottom = element_line(),
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.25, lineend="square"),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 5.5),
    axis.text = element_text(size = 5, color = "black"),
    plot.margin = unit(c(5.5,5.5,5.5,5.5), "pt")
  )
}

ggplot(umap_df)+
    coord_fixed(ratio= (max(umap_df$UMAP1)-min(umap_df$UMAP1)) / (max(umap_df$UMAP2)-min(umap_df$UMAP2)) )+
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
    plot.margin = unit(c(0,0,0,0), "pt")
  )

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

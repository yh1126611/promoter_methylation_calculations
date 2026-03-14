# Enter median MP information of species/class as a data frame with a, b (broader promoter delination point), custom color (hex code) and cpsize (core promoter size).

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

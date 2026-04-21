library(data.table)
library(ggplot2)
library(ggdensity)
library(cowplot)
library(scales)

class_levels <- c("Mammalia", "Aves", "Reptilia", "Amphibia", "Actinopterygii", "Chondrichthyes")
class_colors <- c(Mammalia = "#E69F00", Aves = "#00796B", Reptilia = "#A6D609", 
                  Amphibia = "#984EA3", Actinopterygii = "#56B4E9", Chondrichthyes = "#0072B2")

ALPHA=0.14
ggplot(df) +
  coord_cartesian(clip = "off") +
geom_point(mapping=aes(x=`Genome.size`,y=`Promoter.size`,color=`Class`),size=.3)+
  scale_y_continuous(name="Promoter size (kb)", labels=scales::label_number(scale=1e-3), expand=c(0,0), limits=c(0, max(df$`Promoter.size`)))+
  scale_x_continuous(name="Genome size (Gb)", labels=scales::label_number(scale=1e-9), expand=c(0,0), limits=c(0, max(df$`Genome.size`)))+
  scale_color_manual(values=class_colors, na.value="pink")+
  theme(plot.background=element_blank(),
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.ticks=element_blank(),
        axis.line.y.left=element_line(linewidth=0.25, lineend="square", color="black"),
        axis.line.x.bottom = element_line(linewidth=0.25, lineend="square", color="black"),
        axis.title.y=element_text(size=5, color="black"),
        axis.text.y=element_text(size=5, color="black"),
        axis.title.x=element_text(size=5, color="black"),
        axis.text.x=element_text(size=5, color="black"),
        legend.position="none")

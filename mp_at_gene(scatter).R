library(ggplot2)
library(data.table)
library(dplyr)
library(stringr)
library(cowplot)
library(scales)
library(png)
library(grid)
library(ggpubr)
setwd("D:/R/WD/MP_TSS/mp_at_gene")

blank = ggplot() + theme_void()

files_GAPDH <- list.files(pattern = "^GAPDH_.*\\.tsv$")
files_ACTB <- list.files(pattern = "^ACTB_.*\\.tsv$")

preprocessor <- function(file){
  species_name = str_match(file, "_(.*)\\.tsv$")[,2]
  dt = bind_cols(fread(file), species_name)
  colnames(dt) = c("Distance", "MP", "Strand", "Assembly")
  dt[Strand == "-", Distance := Distance * -1]
}

data_matrix_GAPDH = bind_rows(lapply(files_GAPDH, preprocessor))
data_matrix_ACTB = bind_rows(lapply(files_ACTB, preprocessor))

comprehensive <- fread("comprehensive.tsv")
data_matrix_GAPDH_merged = merge(data_matrix_GAPDH, comprehensive, by="Assembly", all.x=TRUE)
data_matrix_ACTB_merged = merge(data_matrix_ACTB, comprehensive, by="Assembly", all.x=TRUE)

#ggplot(data_matrix_GAPDH_merged, aes(x=Distance, y=MP)) + geom_point(aes(color=Color), alpha=0.1, stroke=0) +  scale_color_identity()

plot_gene_mp <- function(MATRIX, COLUMN, TYPE, TEXT, ALPHA) {
  class_colors=c("Amphibia"="#984EA3", "Aves"="#00796B", "Reptilia"="#A6D609", "Actinopterygii"="#56B4E9", "Sarcopterygii"="#A6761D", "Chondrichthyes"="#0072B2", "Mammalia"="#E69F00")
  ggplot(MATRIX[MATRIX[[COLUMN]] %in% TYPE, ],aes(x = Distance, y = MP)) +
    #annotation_raster(PIC, xmin = -5000, xmax = Inf, ymin = 25, ymax = Inf) +
    geom_point(aes(color=`Class (강)`), alpha = ALPHA, stroke = NA, shape = 19, size=0.25) +
    scale_color_manual(values=class_colors, na.value="black")+
    theme(
      plot.background = element_blank(),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.line.y.left = element_line(linewidth = 0.1),
      axis.line.x.bottom = element_line(linewidth = 0.1),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 5.5, color="black"),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 5.5, color="black"),
      legend.position = "none"
    ) +
    scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
    scale_x_continuous(labels = c("-10", "0", "10"), limits = c(-10000, 10000), breaks = c(-10000, 0, 10000)) +
    # annotate("text", y=25, x=10000, vjust=1, hjust="inward", label = TEXT, color = "black", size = 1.94) +
    coord_cartesian(clip = "off")
}

blank <- ggplot() + theme_void()
none=readPNG("none.png"); mammal=readPNG("mammal.png"); bird=readPNG("bird.png"); reptile=readPNG("reptile.png"); amphibian=readPNG("amphibian.png"); fish=readPNG("fish.png")

# Legend

class_order <- c("Mammal", "Bird", "Reptile", "Amphibian")
df <- data.frame(X = 1,Y = 1,Class = c("Mammal", "Bird", "Reptile", "Amphibian"))
df$Class <- factor(df$Class, levels = class_order)
fake <- ggplot(df, aes(x = X, y = Y, color = Class)) +
  geom_point(size = 3) +
  scale_color_manual(name="Species:", values = c("Mammal"="#E69F00", "Bird"="#00796B", "Reptile"="#A6D609", "Amphibian"="#984EA3")) +
  guides(color = guide_legend(ncol = 1)) +
  theme(
    #legend.position = "top",
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill="transparent", color=NA),
    legend.title = element_text(size=6, face="bold"),
    legend.text = element_text(size=6),
    legend.key.size=unit(.75, "lines")
  )
legend_actb = as_ggplot(get_legend(fake))

class_order <- c("Mammal (31, 32)", "Bird (17, 17)", "Reptile (8, 9)", "Amphibian (3, 4)", "Lobe-finned fish (0, 1)", "Ray-finned fish (0, 16)")
df <- data.frame(X = 1,Y = 1,Class = c("Mammal (31, 32)", "Bird (17, 17)", "Reptile (8, 9)", "Amphibian (3, 4)", "Lobe-finned fish (0, 1)", "Ray-finned fish (0, 16)"))
df$Class <- factor(df$Class, levels = class_order)
fake <- ggplot(df, aes(x = X, y = Y, color = Class)) +
  geom_point(size = 2) +
  scale_color_manual(name="Genomes (nACTB, nGAPDH):", values = c("Mammal (31, 32)"="#E69F00", "Bird (17, 17)"="#00796B", "Reptile (8, 9)"="#A6D609", "Amphibian (3, 4)"="#984EA3", "Lobe-finned fish (0, 1)"="#A6761D", "Ray-finned fish (0, 16)"="#56B4E9")) +
  guides(color = guide_legend(ncol = 1)) +
  theme(
    #legend.position = "top",
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill="transparent", color=NA),
    legend.title = element_text(size=5.5, face="bold"),
    legend.text = element_text(size=5.5),
    legend.key.size=unit(.75, "lines")
  )
legend_gapdh = as_ggplot(get_legend(fake))

# 본 Plot

plot_gene_mp(data_matrix_ACTB_merged, "Class (강)", c("Mammalia", "Aves", "Reptilia", "Amphibia", "Actinopterygii", "Sarcopterygii", "Chondrichthyes"), "nGenome = 59\nnCpG = 34,501", .25)->all_gene
plot_gene_mp(data_matrix_ACTB_merged, "Class (강)", c("Mammalia"), "nGenome = 31\nnCpG = 19,432", 1)->mammal_gene
plot_gene_mp(data_matrix_ACTB_merged, "Class (강)", c("Aves"), "nGenome = 17\nnCpG = 10,165", 1)->bird_gene
plot_gene_mp(data_matrix_ACTB_merged, "Class (강)", c("Reptilia"), "nGenome = 8\nnCpG = 3,845", 1)->reptile_gene
plot_gene_mp(data_matrix_ACTB_merged, "Class (강)", c("Amphibia"), "nGenome = 3\nnCpG = 1,059", 1)->amphibian_gene
#plot_gene_mp(data_matrix_ACTB_merged, "Class (강)", c("Sarcopterygii"), "nGenome = 0\nCpG = 0", 1)->coel_gene
#plot_gene_mp(data_matrix_ACTB_merged, "Class (강)", c("Actinopterygii"), "nGenome = 0\nCpG = 0", 1)->fish_gene
#plot_gene_mp(data_matrix_ACTB_merged, "Class (강)", c("Chondrichthyes"), "nGenome = 0\nCpG = 0", 1)->shark_gene
ACTBleft=plot_grid(mammal_gene,bird_gene,reptile_gene,amphibian_gene,ncol=2)
#ACTB = plot_grid(ACTBleft, legend_gapdh, ncol=2, rel_widths=c(2,1))

plot_gene_mp(data_matrix_GAPDH_merged, "Class (강)", c("Mammalia", "Aves", "Reptilia", "Amphibia", "Actinopterygii", "Sarcopterygii", "Chondrichthyes"), "nGenome = 79\nnCpG = 38,575", .25)->all_gene
plot_gene_mp(data_matrix_GAPDH_merged, "Class (강)", c("Mammalia"), "nGenome = 32\nnCpG = 19,377", 1)->mammal_gene
plot_gene_mp(data_matrix_GAPDH_merged, "Class (강)", c("Aves"), "nGenome = 17\nnCpG = 8,174", 1)->bird_gene
plot_gene_mp(data_matrix_GAPDH_merged, "Class (강)", c("Reptilia"), "nGenome = 9\nnCpG = 3,084", 1)->reptile_gene
plot_gene_mp(data_matrix_GAPDH_merged, "Class (강)", c("Amphibia"), "nGenome = 4\nnCpG = 1,714", 1)->amphibian_gene
plot_gene_mp(data_matrix_GAPDH_merged, "Class (강)", c("Sarcopterygii"), "nGenome = 1\nCpG = 314", 1)->coel_gene
plot_gene_mp(data_matrix_GAPDH_merged, "Class (강)", c("Actinopterygii"), "nGenome = 16\nCpG = 5,912", 1)->fish_gene
#plot_gene_mp(data_matrix_GAPDH_merged, "Class (강)", c("Chondrichthyes"), "nGenome = 0\nCpG = 0", 1)->shark_gene
#GAPDHtop=plot_grid(all_gene,mammal_gene,bird_gene,reptile_gene,amphibian_gene,coel_gene, ncol=3)
#GAPDHbot=plot_grid(fish_gene, legend_gapdh, ncol=2, rel_widths=c(1,2))
#GAPDH=plot_grid(GAPDHtop, GAPDHbot, ncol=1, rel_heights=c(20,10));GAPDH
GAPDH = plot_grid(mammal_gene, bird_gene, reptile_gene, amphibian_gene, coel_gene, fish_gene, ncol=3)

# overall labels addition

rotated_label <- ggdraw() + draw_grob(textGrob("MP (%)", rot = 90, y = 0.5, just = "center", gp = gpar(fontsize = 5.5, fontface = "plain")))

ACTB_quarterTitled=plot_grid(rotated_label, ACTBleft, ncol=3, rel_widths=c(1,20))
ACTB_halfTitled=plot_grid(ACTB_quarterTitled, blank, ncol=1, rel_heights=c(20,1),
                      labels=c("", "Distance from ACTB TSS (kbp)"), label_size=5.5, label_fontface="plain", label_y=0.5, vjust=0.5, label_x=0.5, hjust=0.5)
ACTB_titled=plot_grid(legend_gapdh, ACTB_halfTitled, ncol=2, rel_widths=c(10,21),
                      labels=c("", "c"), label_size=8)

GAPDH_halfTitled=plot_grid(rotated_label, GAPDH, ncol=2, rel_widths=c(1,30))
GAPDH_titled=plot_grid(GAPDH_halfTitled, blank, ncol=1, rel_heights=c(30,1),
                      labels=c("", "Distance from GAPDH TSS (kbp)"), label_size=5.5, label_fontface="plain", label_y=0.5, vjust=0.5, label_x=0.5, hjust=0.5)
ACTB_titled
GAPDH_titled










#Human only

human = readPNG("humanistic.png")

actb = preprocessor("ACTB_T2T-CHM13v2.0.tsv")
actb_exons = fread("actb_exons.txt"); colnames(actb_exons)=c("Xmin","Xmax")
actb_genes = fread("actb_genes.txt"); colnames(actb_genes)=c("Xmin", "Xmax", "gene_name"); actb_genes$Xmin[actb_genes$Xmin< -10000]=-10000; actb_genes$Xmax[actb_genes$Xmax>10000]=10000
actb_transcripts = fread("actb_transcripts.txt"); colnames(actb_transcripts)=c("Xmin", "Xmax"); actb_transcripts$Xmin[actb_transcripts$Xmin< -10000]=-10000; actb_transcripts$Xmax[actb_transcripts$Xmax>10000]=10000
actb_genes$time = (round((actb_genes$Xmax-actb_genes$Xmin)/1000)-1); actb_genes$space=(actb_genes$Xmax-actb_genes$Xmin)/round((actb_genes$Xmax-actb_genes$Xmin)/1000)
dt <- as.data.table(actb_genes)
dt <- dt[, {
  positions <- seq(Xmin+space, by = space, length.out = time)
  .(gene_name = rep(gene_name, time), position = positions)
}, by = 1:nrow(dt)]

ggplot()+
  #annotation_raster(human, xmin=-10000, xmax=10000, ymin=0, ymax=100)+
  geom_point(data=actb, mapping=aes(x=Distance, y=MP), size=.01)+
  scale_x_continuous(limits=c(-10000, 10000), breaks=c(-10000,10000), labels=c("5,540,601", "5,520,601"), name="Chr. 7", position="top")+
  scale_y_continuous(limits=c(0,100), breaks=c(0,25,50,75,100), name="MP (%)")+
  theme(plot.background = element_blank(), panel.background = element_blank(), panel.grid = element_blank(),
        axis.line.y.left = element_line(linewidth=0.1), axis.line.x.top=element_line(linewidth=0.1),
        axis.text=element_text(size=5.5, color="black"), axis.title=element_text(size=5.5, color="black"), axis.ticks=element_blank())->actb_plot

ggplot()+
  geom_rect(data=actb_exons, mapping=aes(xmax=Xmax, xmin=Xmin, ymax=-10, ymin=-20), fill="#384860")+
  geom_rect(data=actb_genes, mapping=aes(xmax=Xmax, xmin=Xmin, ymax=-11.25, ymin=-18.75), fill="#384860")+
  geom_rect(data=actb_transcripts, mapping=aes(xmax=Xmax, xmin=Xmin, ymax=-2.5, ymin=-7.5), fill="#97a6c4")+
  geom_text(data=actb_genes, mapping=aes(x=Xmin+(Xmax-Xmin)/2, y=-15, label=gene_name), color="white", size=2)+
  geom_text(data=dt[dt$gene_name=="ACTB",], mapping=aes(x=position), y=-5, vjust=0.5, hjust=0.5, label=">", color="white", size=2)+
  geom_text(data=dt[dt$gene_name=="LOC221946",], mapping=aes(x=position), y=-5, vjust=0.5, hjust=0.5, label="<", color="white", size=2)+
  scale_x_continuous(labels = label_comma(), limits=c(-10000, 10000), breaks=c(-10000,0,10000), name="Distance from ACTB TSS (bp)")+
  scale_y_continuous(limits=c(-20,0), breaks=c(-5, -15), labels=c("Transcript", "Gene"))+
  theme(plot.background = element_blank(), panel.background = element_blank(), panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.text=element_text(size=5.5, color="black"), axis.title.x=element_text(size=5.5, color="black"), axis.title.y=element_blank(), axis.ticks=element_blank(),
        axis.text.y=element_text(size=5.5, color="black"),
        axis.line.x.bottom=element_line(linewidth=0.1))->actb_plot_genes;actb_plot_genes

align_plots(actb_plot + theme(plot.margin = margin(5.5,5.5,0,5.5)), actb_plot_genes + theme(plot.margin = margin(0,5.5,5.5,5.5)), align="v", axis="lr")->actb_aligned





gapdh = preprocessor("GAPDH_T2T-CHM13v2.0.tsv")
gapdh_exons = fread("gapdh_exons.txt"); colnames(gapdh_exons)=c("Xmin","Xmax", "gene_name")
gapdh_genes = fread("gapdh_genes.txt"); colnames(gapdh_genes)=c("Xmin", "Xmax", "gene_name"); gapdh_genes$Xmin[gapdh_genes$Xmin< -10000]=-10000; gapdh_genes$Xmax[gapdh_genes$Xmax>10000]=10000
gapdh_transcripts = fread("gapdh_transcripts.txt"); colnames(gapdh_transcripts)=c("Xmin", "Xmax", "gene_name"); gapdh_transcripts$Xmin[gapdh_transcripts$Xmin< -10000]=-10000; gapdh_transcripts$Xmax[gapdh_transcripts$Xmax>10000]=10000
gapdh_transcripts$time = (round((gapdh_transcripts$Xmax-gapdh_transcripts$Xmin)/1000)-1); gapdh_transcripts$space=(gapdh_transcripts$Xmax-gapdh_transcripts$Xmin)/round((gapdh_transcripts$Xmax-gapdh_transcripts$Xmin)/1000)

dt <- as.data.table(gapdh_transcripts[gapdh_transcripts$gene_name=="GAPDH",][1,])
dt_gapdh1 <- dt[, {
  positions <- seq(Xmin+space, by = space, length.out = time)
  .(gene_name = rep(gene_name, time), position = positions)
}, by = 1:nrow(dt)]

dt <- as.data.table(gapdh_transcripts[gapdh_transcripts$gene_name=="GAPDH",][2,])
dt_gapdh2 <- dt[, {
  positions <- seq(Xmin+space, by = space, length.out = time)
  .(gene_name = rep(gene_name, time), position = positions)
}, by = 1:nrow(dt)]

dt <- as.data.table(gapdh_transcripts[gapdh_transcripts$gene_name=="IFFO1",][1,])
dt_iffo11 <- dt[, {
  positions <- seq(Xmin+space, by = space, length.out = time)
  .(gene_name = rep(gene_name, time), position = positions)
}, by = 1:nrow(dt)]

dt <- as.data.table(gapdh_transcripts[gapdh_transcripts$gene_name=="IFFO1",][2,])
dt_iffo12 <- dt[, {
  positions <- seq(Xmin+space, by = space, length.out = time)
  .(gene_name = rep(gene_name, time), position = positions)
}, by = 1:nrow(dt)]

dt <- as.data.table(gapdh_transcripts[gapdh_transcripts$gene_name=="IFFO1",][3,])
dt_iffo13 <- dt[, {
  positions <- seq(Xmin+space, by = space, length.out = time)
  .(gene_name = rep(gene_name, time), position = positions)
}, by = 1:nrow(dt)]

dt <- as.data.table(gapdh_transcripts[gapdh_transcripts$gene_name=="IFFO1",][4,])
dt_iffo14 <- dt[, {
  positions <- seq(Xmin+space, by = space, length.out = time)
  .(gene_name = rep(gene_name, time), position = positions)
}, by = 1:nrow(dt)]

dt <- as.data.table(gapdh_transcripts[gapdh_transcripts$gene_name=="NCAPD2",])
dt_ncapd2 <- dt[, {
  positions <- seq(Xmin+space, by = space, length.out = time)
  .(gene_name = rep(gene_name, time), position = positions)
}, by = 1:nrow(dt)]

ggplot()+
  #annotation_raster(human, xmin=-10000, xmax=10000, ymin=0, ymax=100)+
  geom_point(data=gapdh, mapping=aes(x=Distance, y=MP), size=.01)+
  scale_x_continuous(limits=c(-10000, 10000), breaks=c(-10000,10000), labels=c("6,524,517", "6,544,517"), name="Chr. 12", position="top")+
  scale_y_continuous(limits=c(0,100), breaks=c(0,25,50,75,100), name="MP (%)")+
  theme(plot.background = element_blank(), panel.background = element_blank(), panel.grid = element_blank(),
        axis.line.y.left = element_line(linewidth=0.1), axis.line.x.top=element_line(linewidth=0.1),
        axis.text=element_text(size=5.5, color="black"), axis.title=element_text(size=5.5, color="black"), axis.ticks=element_blank())->gapdh_plot

# ggplot()+
#   geom_rect(data=gapdh_exons[gapdh_exons$gene_name!="GAPDH-DT",], mapping=aes(xmax=Xmax, xmin=Xmin, ymax=-10, ymin=-14), fill="#384860")+
#   geom_rect(data=gapdh_exons[gapdh_exons$gene_name=="GAPDH-DT",], mapping=aes(xmax=Xmax, xmin=Xmin, ymax=-16, ymin=-20), fill="#384860")+
#   geom_rect(data=gapdh_genes[gapdh_genes$gene_name!="GAPDH-DT",], mapping=aes(xmax=Xmax, xmin=Xmin, ymax=-10.5, ymin=-13.5), fill="#384860")+
#   geom_rect(data=gapdh_genes[gapdh_genes$gene_name=="GAPDH-DT",], mapping=aes(xmax=Xmax, xmin=Xmin, ymax=-16.5, ymin=-19.5), fill="#384860")+
#   geom_rect(data=gapdh_transcripts[gapdh_transcripts$gene_name!="GAPDH-DT",], mapping=aes(xmax=Xmax, xmin=Xmin, ymax=-1, ymin=-4), fill="#97a6c4")+
#   geom_rect(data=gapdh_transcripts[gapdh_transcripts$gene_name=="GAPDH-DT",], mapping=aes(xmax=Xmax, xmin=Xmin, ymax=-6, ymin=-9), fill="#97a6c4")+
#   geom_text(data=dt[dt$gene_name!="GAPDH-DT",], mapping=aes(x=position), y=-2.5, vjust=0.5, hjust=0.5, label=">", color="white", size=1.4)+
#   geom_text(data=dt[dt$gene_name=="GAPDH-DT",], mapping=aes(x=position), y=-7.5, vjust=0.5, hjust=0.5, label=">", color="white", size=1.4)+
#   geom_text(data=gapdh_genes[gapdh_genes$gene_name!="GAPDH-DT",], mapping=aes(x=Xmin+(Xmax-Xmin)/2, y=-7, label=gene_name), color="black", size=2)+
#   geom_text(data=gapdh_genes[gapdh_genes$gene_name=="GAPDH-DT",], mapping=aes(x=Xmin+(Xmax-Xmin)/2, y=-14, label=gene_name), color="black", size=1.4)+
#   scale_x_continuous(labels = label_comma(), limits=c(-10000, 10000), breaks=c(-10000,0,10000), name="Distance from GAPDH TSS (bp)")+
#   scale_y_continuous(limits=c(-20,0), breaks=c(-5, -15), labels=c("Transcript", "Gene"))+
#   theme(plot.background = element_blank(), panel.background = element_blank(), panel.grid = element_blank(),
#         axis.line = element_blank(),
#         axis.text=element_text(size=5.5, color="black"), axis.title.x=element_text(size=5.5, color="black"), axis.title.y=element_blank(), axis.ticks=element_blank(),
#         axis.text.y=element_text(size=5.5, color="black"),
#         axis.line.x.bottom=element_line(linewidth=0.1))->gapdh_plot_genes

ggplot()+
  #NCAPD2
  geom_rect(data=gapdh_exons[gene_name=="NCAPD2",], mapping=aes(xmax=Xmax, xmin=Xmin, ymax=-15, ymin=-20), fill="#384860")+
  geom_rect(data=gapdh_genes[gene_name=="NCAPD2",], mapping=aes(xmax=Xmax, xmin=Xmin, ymax=-15.625, ymin=-19.375), fill="#384860")+
  geom_rect(data=gapdh_transcripts[gene_name=="NCAPD2",], mapping=aes(xmax=Xmax, xmin=Xmin, ymax=-5, ymin=-8.75), fill="#97a6c4")+
  geom_text(data=gapdh_genes[gene_name=="NCAPD2",], mapping=aes(x=Xmin+(Xmax-Xmin)/2, y=-17.5, label=gene_name), color="white", size=1.5)+
  geom_text(data=dt_ncapd2, mapping=aes(x=position), y=-6.875, vjust=0.5, hjust=0.5, label=">", color="white", size=1.5)+
  #IFFO1
  geom_rect(mapping=aes(xmin=7134, xmax=10000, ymax=0, ymin=-2.5), fill="#97a6c4")+
  geom_rect(mapping=aes(xmin=5022, xmax=10000, ymax=-3.75, ymin=-6.25), fill="#97a6c4")+
  geom_rect(mapping=aes(xmin=4444, xmax=10000, ymax=-7.5, ymin=-10), fill="#97a6c4")+
  geom_rect(mapping=aes(xmin=4444, xmax=10000, ymax=-11.25, ymin=-13.75), fill="#97a6c4")+
  geom_rect(data=gapdh_exons[gene_name=="IFFO1",], mapping=aes(xmax=Xmax, xmin=Xmin, ymax=-15, ymin=-20), fill="#384860")+
  geom_rect(data=gapdh_genes[gene_name=="IFFO1",], mapping=aes(xmax=Xmax, xmin=Xmin, ymax=-15.625, ymin=-19.375), fill="#384860")+
  geom_text(data=gapdh_genes[gene_name=="IFFO1",], mapping=aes(x=Xmin+(Xmax-Xmin)/2, y=-17.5, label=gene_name), color="white", size=1.5)+
  geom_text(data=dt_iffo14, mapping=aes(x=position), y=-1.25, vjust=0.5, hjust=0.5, label="<", color="white", size=1.5)+
  geom_text(data=dt_iffo13, mapping=aes(x=position), y=-5, vjust=0.5, hjust=0.5, label="<", color="white", size=1.5)+
  geom_text(data=dt_iffo12, mapping=aes(x=position), y=-8.75, vjust=0.5, hjust=0.5, label="<", color="white", size=1.5)+
  geom_text(data=dt_iffo11, mapping=aes(x=position), y=-12.5, vjust=0.5, hjust=0.5, label="<", color="white", size=1.5)+
  #GAPDH
  geom_rect(mapping=aes(xmin=747, xmax=3854, ymax=-1.875, ymin=-5.625), fill="#97a6c4")+
  geom_rect(mapping=aes(xmin=0, xmax=3854, ymax=-8.125, ymin=-11.875), fill="#97a6c4")+
  geom_text(data=dt_gapdh2, mapping=aes(x=position), y=-3.75, vjust=0.5, hjust=0.5, label=">", color="white", size=1.5)+
  geom_text(data=dt_gapdh1, mapping=aes(x=position), y=-10, vjust=0.5, hjust=0.5, label=">", color="white", size=1.5)+
  geom_rect(data=gapdh_exons[gene_name=="GAPDH",], mapping=aes(xmax=Xmax, xmin=Xmin, ymax=-15, ymin=-20), fill="#384860")+
  geom_rect(data=gapdh_genes[gene_name=="GAPDH",], mapping=aes(xmax=Xmax, xmin=Xmin, ymax=-15.625, ymin=-19.375), fill="#384860")+
  geom_text(data=gapdh_genes[gene_name=="GAPDH",], mapping=aes(x=Xmin+(Xmax-Xmin)/2, y=-17.5, label=gene_name), color="white", size=1.5)+
  
  scale_x_continuous(labels = label_comma(), limits=c(-10000, 10000), breaks=c(-10000,0,10000), name="Distance from GAPDH TSS (bp)")+
  scale_y_continuous(limits=c(-20,0), breaks=c(-6.875, -17.5), labels=c("Transcript", "Gene"))+
  theme(plot.background = element_blank(), panel.background = element_blank(), panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.text=element_text(size=5.5, color="black"), axis.title.x=element_text(size=5.5, color="black"), axis.title.y=element_blank(), axis.ticks=element_blank(),
        axis.text.y=element_text(size=5.5, color="black"),
        axis.line.x.bottom=element_line(linewidth=0.1))->gapdh_plot_genes;gapdh_plot_genes

align_plots(gapdh_plot + theme(plot.margin = margin(5.5,5.5,0,5.5)), gapdh_plot_genes + theme(plot.margin = margin(0,5.5,5.5,5.5)), align="v", axis="lr")->gapdh_aligned







plot_grid(blank, actb_aligned[[1]], gapdh_aligned[[1]], blank, actb_aligned[[2]], gapdh_aligned[[2]], blank, ACTB_titled, GAPDH_titled,
          ncol=3, rel_heights=c(11, 7, 21), rel_widths=c(.5,10,10),
          labels=c("a", "", "b", "", "", "", "", "", "d"), label_size=8) -> complete

comp_height= 11 + 7 + 21
comp_width=1+30+1+30+1
ggsave("complete.png", complete, height=180/comp_width*comp_height, width=180, units="mm")
#ggsave("complete.svg", complete, height=180/comp_width*comp_height, width=180, units="mm")
ggsave("complete.pdf", complete, height=180/comp_width*comp_height, width=180, units="mm")

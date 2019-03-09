library("arghqtl")
library("qtl")

mloc <- read.table('data_raw/genetic_map_agren_etal_2013.csv', sep = '\t', header=T)

allqtl_clusters <- readRDS("output/clusters_all_models.rds")
qtlmodels <- readRDS("output/qtl_stepwise_models.rds")
qtlfits   <- readRDS("output/qtl_model_fits.rds")


# set up plot object
pos_plot <- plotQTL(1:5, mloc, nlanes = 10, left_gap = 2, right_gap = 0)


colvec <- rep(c("red3", "blue4"), each=10)
# Details of pleiotropy boxes
boxes <- allqtl_clusters$boxes[allqtl_clusters$boxes$chr, ]
ix    <- match(allqtl_clusters$summary$chr[allqtl_clusters$summary$n.qtl > 1], colnames(pos_plot$lane_margins))
pos   <- 0.5 * (allqtl_clusters$summary$pos_min[allqtl_clusters$summary$n.qtl > 1] + 
                allqtl_clusters$summary$pos_max[allqtl_clusters$summary$n.qtl > 1])

ts <- 0.2
par(mfrow=c(2,1), mar=c(0,1,0,0)+0.1, oma=c(0,1,0,0), cex=0.8)
# Empty plot
plotQTL_base(pos_plot, by=20, maxy = 20)
text(-6.5, -40, "Marker distance (cM)", srt=90, xpd=NA)
# indicate extent of overlapping QTL with grey boxes.
plotQTL_boxes(
  pos_plot,
  allqtl_clusters$boxes,
  border = F,
  col = 'gray75',
  lwd = 2
)

pos_counter <- 1
for(t in c("frut", "seed", "tofu", "mass", "surv")){
  # plot QTL positions
  for(i in c("it2010", "it2011")){
    plotQTL_qtl(pos_plot, qtlmodels[[t]][[i]], qtlfits[[t]][[i]],
                fill_cutoff = 15.2, lane = pos_counter, col='red3')
    pos_counter <- pos_counter + 1
  }
}
# Lane dividers
segments(pos_plot$lane_margins, 0,
         pos_plot$lane_margins, rep(-pos_plot$track_lengths, each=11),
         lty=2, col='gray50')
segments(pos_plot$lane_margins[seq(1,11,2),], 0,
         pos_plot$lane_margins[seq(1,11,2),], rep(-pos_plot$track_lengths, each=6),
         lty=1, col=1)
# Lane labels
text(pos_plot$lane_margins[seq(2,11,2),], 15,
     c("Fr/RP", "Sd/Fr", "Sd/RP", "Sd mass", "Surv"), col=1, cex=0.8, srt=90, adj=c(0.5,0.5))
text(pos_plot$lane_centres, 1, c("10", "11"), cex=0.6, adj=c(0,0.5), srt=90, col='gray30')
# Label pleiotropy boxes
text(pos_plot$lane_margins[1,ix]- 0.2,
     -pos,
     paste("Q",1:13, sep=""),
     adj=c(0.5,0),
     cex=0.7,
     srt=90)



plotQTL_base(pos_plot, by=20, maxy = 20)
text(-6.5, -40, "Marker distance (cM)", srt=90, xpd=NA)
plotQTL_boxes(
  pos_plot,
  allqtl_clusters$boxes,
  border = F,
  col = 'gray75',
  lwd = 2
)
pos_counter <- 1
for(t in c("frut", "seed", "tofu", "mass", "surv")){
  # plot QTL positions
  for(i in c("sw2010", "sw2011")){
    plotQTL_qtl(pos_plot, qtlmodels[[t]][[i]], qtlfits[[t]][[i]],
                fill_cutoff = 15.2, lane = pos_counter, col='blue4')
    pos_counter <- pos_counter + 1
  }
}
# Lanes dividers
segments(pos_plot$lane_margins, 0,
         pos_plot$lane_margins, rep(-pos_plot$track_lengths, each=11),
         lty=2, col='gray50')
segments(pos_plot$lane_margins[seq(1,11,2),], 0,
         pos_plot$lane_margins[seq(1,11,2),], rep(-pos_plot$track_lengths, each=6),
         lty=1, col=1)
# Lane labels
text(pos_plot$lane_margins[seq(2,11,2),], 15,
     c("Fr/RP", "Sd/Fr", "Sd/RP", "Sd mass", "Surv"), col=1, cex=0.8, srt=90, adj=c(0.5,0.5))
text(pos_plot$lane_centres, 1, c("10", "11"), cex=0.6, adj=c(0,0.5), srt=90, col='gray30')
# Label pleiotropy boxes
text(pos_plot$lane_margins[1,ix]-0.2,
     -pos,
     paste("Q",1:13, sep=""),
     adj=c(0.5,0),
     cex=0.7,
     srt=90)


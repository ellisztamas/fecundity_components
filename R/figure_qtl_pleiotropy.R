mloc <- read.table('data_files/genetic_map_agren_etal_2013.csv', sep = '\t', header=T)

allqtl_clusters <- readRDS("output/clusters_all_models.rds")

colvec=c("red3", "red3","blue4","blue4")
#setEPS()
#postscript('../figures/QTL_pleiotropy.eps', width=25/2.56, height = 16/2.54)
pdf('figures/04_QTL_pleiotropy.pdf', width=25/2.56, height = 16.9/2.54)

ts <- 0.2
par(mfrow=c(2,1), mar=c(1,3,0,0)+0.1, oma=c(0,0,0,0), cex=0.8)

chr123 <- plotQTL(1:3, mloc, nlanes = 20)
plotQTL_base(chr123, by=20, maxy = 20)

pos_counter <- 1
# indicate extent of overlapping QTL with grey boxes.
boxes <- allqtl_clusters$boxes[allqtl_clusters$boxes$chr %in% 1:3, ]
plotQTL_boxes(
  chr123,
  allqtl_clusters$boxes,
  border = F,
  col = 'gray75',
  lwd = 2
)
for(t in c("frut", "seed", "tofu", "mass", "surv")){
  # plot QTL positions
  for(i in 1:4){
    plotQTL_qtl(chr123, qtlmodels[[t]][[i]], qtlfits[[t]][[i]],
                fill_cutoff = 15.2, lane = i+pos_counter-1, col=colvec[i])
  }
  pos_counter <- pos_counter + 4
}
# Lanes dividers
segments(chr123$lane_margins, 0,
         chr123$lane_margins, rep(-chr123$track_lengths, each=21),
         lty=2, col='gray50')
segments(chr123$lane_margins[seq(1,22,4),], 00,
         chr123$lane_margins[seq(1,22,4),], rep(-chr123$track_lengths, each=6),
         lty=1, col=1)
# plotQTL_lanes(chr123, lty=0)
# Lane labels
text(chr123$lane_margins[seq(3,21,4),], 15,
     c("Fruits\nper RP", "Seeds\nper fruit", "Seeds\nper RP", "Seed\nmass", "Survival"), col=1, cex=0.7, srt=-90, adj=c(0.5,0.5))
text(chr123$lane_centres, 2, c("10", "11"), cex=0.6, adj=c(0.5,0.5), srt=0, col='gray30')
text(chr123$lane_margins[seq(2,21,2),], 6, c("IT", "SW"), cex=0.6, adj=c(0.5,0.5), srt=0, col=c('red3','blue4'))
os <- 0.2
segments(chr123$lane_margins[seq(1,19,2),] + os, 4, 
         chr123$lane_margins[seq(3,21,2),] - os, 4, col=c('red3','blue4'), lwd=1)
# label pleiotropic regions
# text(chr123$lane_margins[21, c(1,1,1,1,2,3)] + ts,
#      -c(11, 23, 60, 80, 53, 63),
#      c("fm", "m","fs","fs","ms","f"), adj = 0)



chr45 <- plotQTL(4:5, mloc, nlanes = 20)
plotQTL_base(chr45, by=20, maxy = 20)

# indicate extent of overlapping QTL with grey boxes.
boxes <- allqtl_clusters$boxes[allqtl_clusters$boxes$chr %in% 4:5, ]
plotQTL_boxes(
  chr45,
  allqtl_clusters$boxes,
  border = F,
  col = 'gray75',
  lwd = 2
)

pos_counter <- 1
for(t in c("frut", "seed", "tofu", "mass", "surv")){
  # plot QTL positions
  for(i in 1:4){
    plotQTL_qtl(chr45, qtlmodels[[t]][[i]], qtlfits[[t]][[i]],
                fill_cutoff = 15.2, lane = i+pos_counter-1, col=colvec[i])
  }
  pos_counter <- pos_counter + 4
}

segments(chr45$lane_margins, 0,
         chr45$lane_margins, rep(-chr45$track_lengths, each=21),
         lty=2, col='gray50')
segments(chr45$lane_margins[seq(1,22,4),], 00,
         chr45$lane_margins[seq(1,22,4),], rep(-chr45$track_lengths, each=6),
         lty=1, col=1)
# plotQTL_lanes(chr45, lty=0)
# Lane labels
text(chr45$lane_margins[seq(3,21,4),], 15,
     c("Fruits\nper RP", "Seeds\nper fruit", "Seeds\nper RP", "Seed\nmass", "Survival"), col=1, cex=0.7, srt=-90, adj=c(0.5,0.5))
text(chr45$lane_centres, 2, c("10", "11"), cex=0.6, adj=c(0.5,0.5), srt=0, col='gray30')
text(chr45$lane_margins[seq(2,21,2),], 6, c("IT", "SW"), cex=0.6, adj=c(0.5,0.5), srt=0, col=c('red3','blue4'))
os <- 0.2
segments(chr45$lane_margins[seq(1,19,2),] + os, 4, 
         chr45$lane_margins[seq(3,21,2),] - os, 4, col=c('red3','blue4'), lwd=1)

# label pleiotropic regions
# text(chr45$lane_margins[21, c("4","5","5","5","5","5")] + ts,
#      -c(52, 1, 9, 16,55,74),
#      c("f","ms", "fs","m","fs", "fs"), adj = 0)
dev.off()

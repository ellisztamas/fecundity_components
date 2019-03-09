library("arghqtl")
library("qtl")

mloc <- read.table('data_raw/genetic_map_agren_etal_2013.csv', sep = '\t', header=T)
# Import data on fitness boxes from Agren 2013.
boxes <- read.csv("data_raw/fitness_boxes.csv")
boxes$lower <- -boxes$lower
boxes$upper <- -boxes$upper
# Set up colours for each lane
colvec=c("red3", "red3","blue4","blue4")

par(mfrow=c(1,1), mar=c(1,3,0,0)+0.1, oma=c(0,0,0,0))

pos_plot <- plotQTL(1:5, mloc, nlanes = 4)
plotQTL_base(pos_plot, by=20, maxy = 20)
text(35, -95, "Chromosome")

# Plot QTL fitness boxes from Agren 2013.
plotQTL_boxes(pos_plot, boxes[,c("chr", "lower", "upper")], col='gray80', border=F)
# Plot QTL in each lane
pos_counter <- 1
for(i in c("it2010", "it2011", "sw2010", "sw2011")){
  plotQTL_qtl(pos_plot, qtlmodels$tfit[[i]], qtlfits$tfit[[i]],
              fill_cutoff = 15.2, lane = pos_counter, col=colvec[pos_counter])
  pos_counter <- pos_counter + 1
}
# Lanes dividers
segments(pos_plot$lane_margins, 0,
         pos_plot$lane_margins, rep(-pos_plot$track_lengths, each=pos_plot$nlanes+1),
         lty=c(1,2,1,2,1), col=c(1, 'gray50',1,'gray50',1))
# Lane labels
text(pos_plot$lane_centres, 2, c("10", "11"), cex=0.7, adj=c(0.5,0.5), srt=0, col='gray30')
text(pos_plot$lane_margins[seq(2,5,2),], 6, c("It", "Sw"), cex=1, adj=c(0.5,0.5), col=c('red3','blue4'))

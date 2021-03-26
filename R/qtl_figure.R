qtl_figure <- function(stepwise, fits, clusters, marker_locations, expt, title, chr_labels, track_labels, col,bg){
  
  if(length(stepwise) != length(fits)) stop("Lists of stepwise-qtl and fit objects are of different lengths.")
  if(any(! names(stepwise) %in% names(fits))) stop("Entries in lists of stepwise-qtl and fit objects do not match.")
  
  # Set up the plot
  # Create a plotQTL object with all the coordinates of the plot.
  pos_plot <- plotQTL(
    chr = 1:5,
    marker_locations = marker_locations, 
    nlanes = length(stepwise),
    left_gap = 2,
    right_gap = 2
  )
  # Create the empty plot
  plotQTL_base(
    pos_plot,
    by = 20,
    maxy = 0,
    track_labels = chr_labels,
    ylab = ""
  )
  # Y-axis label
  # text(-6.5, -40, "Marker distance (cM)", srt=90, xpd=NA, cex=0.9)
  
  # Annotate the plot.
  # Subplot titles
  # text(-5,15, title, xpd=NA, adj=c(0, 1), cex=1)#, col="red3")
  mtext(title, side=4, line = -2)
  
  # Lane labels
  text(pos_plot$lane_centres, 3,
       track_labels,
       col=1, cex=0.7, srt=90, adj=c(0,0.5))
  
  # indicate extent of overlapping QTL with grey boxes.
  # Details of pleiotropy boxes
  boxes <- clusters$boxes[clusters$summary$n.qtl > 1,][-6,]
  ix    <- match(boxes$chr, colnames(pos_plot$lane_margins))
  pos   <- (boxes$lower + boxes$upper) /2
  # Plot boxes
  plotQTL_boxes(
    pos_plot,
    boxes,
    border = 'gray85',
    col = 'gray85',
    lwd = 2
  )
  # Label pleiotropy boxes
  text(pos_plot$lane_margins[1,ix]- 0.7,
       -pos,
       paste("Q",1:nrow(boxes), sep=""),
       adj=c(0.5,0),
       cex=0.7,
       srt=90)
  plotQTL_ladder(pos_plot)
  
  # Lane dividers
  segments(pos_plot$lane_margins,
           0,
           pos_plot$lane_margins,
           rep(-pos_plot$track_lengths, each=8),
           lty=1, col='gray90', lwd=0.5)
  
  # plot QTL positions for each traits
  pos_counter <- 1
  for(t in c("frut", "seed", "tofu", "mass", "surv", "ffit", "tfit")){
    plotQTL_qtl(
      plotQTL = pos_plot,
      qtl_object = qtl_stepwise[[t]][[expt]],
      model_fit = qtl_fits[[t]][[expt]],
      # fill_cutoff = 15.2,
      lane = pos_counter,
      cex = 0.6,
      col = col[[t]],
      bg=bg[[t]]
    )
    pos_counter <- pos_counter + 1
  }
  
  # # Labels for QTL found for seeds/sdl, but not fruits/sdl
  # os <- 1
  # text(pos_plot$lane_margins[15,3] + os, -19.705, "+")
  # # text(pos_plot$lane_margins[15,4] + os, -20.7, "+")
  # text(pos_plot$lane_margins[15,2] + os, -31.7, "+")
  # text(pos_plot$lane_margins[15,5] + os, -18.5, "+")
  # # Labels for QTL found for fruits/sdl, but not seeds/sdl
  # # text(pos_plot$lane_margins[15,4] + os, -50.0, "-")
  # text(pos_plot$lane_margins[15,3] + os, -1.7, "-")
  
}
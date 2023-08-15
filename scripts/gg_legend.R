# function to extract a legend from a ggplot for plotting separately
# e.g., combining with other panels using cowplot
# get the legend to plot separately
gg_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

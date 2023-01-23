TileDotplot <- function(object, features, col.max = 2, scale.by = "size", scale = T, ...) {
  
  DotPlot(object, 
          features = features, 
          col.max = col.max,
          scale.by = scale.by, 
          scale = scale,
          ...
  ) +
    geom_tile(aes_string(fill = "avg.exp.scaled")) +
    geom_point(aes_string(size = "pct.exp"), shape = 21, colour="white", stroke=1) +
    # scale_colour_viridis(option="magma") +
    guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + 
    RotatedAxis() + 
    coord_flip() +
    coord_fixed(ratio = 1) +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
          axis.title = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
    )+
    scale_color_distiller(palette = "RdBu") +
    scale_fill_distiller(palette = "RdBu") +
    #scale_size(range = c(1,6)) +
    #NoLegend() +
    # theme(legend.position = "none",
    #      
    #     axis.text.x = element_blank()) +
    NULL
}
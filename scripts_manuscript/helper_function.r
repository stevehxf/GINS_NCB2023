## get_gene_cutoff ----
## get the gene that express over a cutoff percentage of specific idents
get_gene_cutoff <- function(object, idents, pct.cutoff){
  # the cell types should be included for percentage calculation
  data <- object@assays$RNA@data %>% as.matrix() %>% as.data.frame()
  meta <- object@meta.data %>% dplyr::select(orig.ident, cell.type)
  pct <- pct.cutoff
  meta.list <- meta %>% 
    split(meta$cell.type)
  meta.list <- meta.list[idents]
  cell.name.list <- meta.list %>% 
    lapply(rownames)
  data.list <- meta.list
  for(i in seq_along(data.list)){
    data.list[[i]] <- data[cell.name.list[[i]]]
  }
  gene.cutoff <- data.list %>% 
    lapply(function(x){
      idx <- rowSums(x>0)/ncol(x) >  pct
      rownames(x)[idx] %>% return()
    }) %>% 
    unlist() %>% 
    unique()
}

## TileDotplot ----
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

## gene.pct --------
## calculate percentage expression in each cell.type
gene.pct <- function(object,pct,group_by){
  list <- SplitObject(object, split.by = group_by) # split the object by group_by into a list
  gene.list.pct <- function(x,pct){
    df <- as.data.frame(x@assays$RNA@data) # fetch the data and convert into a df
    nonzero.cellCount <- rowSums(df != 0) # calculate the number of cells that express each genes.
    total.cellCount <- nrow(x@meta.data) # calculate the total cell number
    cutoff.cellCount <- total.cellCount*pct # calculate the cutoff cell number
    cutoff.gene <- nonzero.cellCount[nonzero.cellCount > cutoff.cellCount] %>% names() # get the genes that pass the cutoff
    return(cutoff.gene)
  }
  pct.gene.list <- lapply(list, gene.list.pct, pct) # get the genes for each group of cells
  unique.genes <- pct.gene.list %>% unlist() %>% unique() # merge 
  return(unique.genes)
}

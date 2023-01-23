## read in the objects for analysis ----
# integration of GINS organoids form timecourse dataset and islets.
cgins.vitro.time.muRM <- readRDS("rds_manuscript/cgins.vitro.time.muRM.rds")
integrated_islet_sub <- readRDS("rds_manuscript/integrated_islet.x_sub.rds")
integrated_islet_sub$source <- integrated_islet_sub$donor
Idents(integrated_islet_sub) <- "cell.type"
integrated_islet_sub.abd <- subset(integrated_islet_sub, idents = c("pp","epsilon"), invert = T)
## combine list containing with cgins.vitro and islets from individual donor
islet.list <- SplitObject(integrated_islet_sub.abd, split.by = "donor")

list_cginsTime.vitro_islet <- c(list(cgins.vitro = cgins.vitro.time.muRM), islet.list[1])

## SCTs normalization
list_cginsTime.vitro_islet <- list_cginsTime.vitro_islet %>% 
  lapply(FUN = SCTransform, method = "glmGamPoi")

## islets+cgins.vitro prepare integration
time.islet.features <- list_cginsTime.vitro_islet %>% SelectIntegrationFeatures()
saveRDS(time.islet.features,"rds_manuscript/time.islet.features.rds")
features <- time.islet.features
list_cginsTime.vitro_islet <- list_cginsTime.vitro_islet %>% 
  PrepSCTIntegration(anchor.features = features)

## RPCA ----
list_cginsTime.vitro_islet <- list_cginsTime.vitro_islet %>% 
  lapply(FUN = RunPCA, features = features)
rpca.anchors <- list_cginsTime.vitro_islet %>%
  FindIntegrationAnchors(
    normalization.method = "SCT",
    anchor.features = features, 
    dims = 1:30, 
    reduction = "rpca", k.anchor = 8)
rpca.integrated_cginsTime.vitro_islet <- 
  IntegrateData(anchorset = rpca.anchors, 
                normalization.method = "SCT", 
                dims = 1:30)

dimensions <- 1:28
rpca.integrated_cginsTime.vitro_islet <- RunPCA(rpca.integrated_cginsTime.vitro_islet, verbose = FALSE)
rpca.integrated_cginsTime.vitro_islet <- RunUMAP(rpca.integrated_cginsTime.vitro_islet, reduction = "pca", dims = dimensions)

p3 <- DimPlot(rpca.integrated_cginsTime.vitro_islet, reduction = "umap", label = F, cols = (my36colors), group.by = "cell.type", split.by = "source", ncol = 4) + NoLegend() + coord_equal()
p4 <- DimPlot(rpca.integrated_cginsTime.vitro_islet, reduction = "umap", label = F, cols = (my36colors), group.by = "cell.type") + coord_equal()
p5 <- DimPlot(rpca.integrated_cginsTime.vitro_islet, reduction = "umap", label = T, cols = (my36colors), group.by = "cell.type", split.by = "source", ncol = 4) + NoLegend() + coord_equal()
p6 <- DimPlot(rpca.integrated_cginsTime.vitro_islet, reduction = "umap", label = T, cols = (my36colors), group.by = "cell.type") + coord_equal()
cowplot::plot_grid(p3,p4,p5,p6, nrow = 2, rel_widths = c(2,2), align = "h")

rpca.integrated_cginsTime.vitro_islet <- rpca.integrated_cginsTime.vitro_islet %>% 
  FindNeighbors(dims = dimensions)

rpca.integrated_cginsTime.vitro_islet <- rpca.integrated_cginsTime.vitro_islet %>% 
  FindClusters(resolution = 0.3)
p1 <- DimPlot(rpca.integrated_cginsTime.vitro_islet, reduction = "umap", label = T, cols = my36colors) + NoLegend()  + coord_equal()
p2 <- DimPlot(rpca.integrated_cginsTime.vitro_islet, reduction = "umap", label = T, cols = (my36colors), group.by = "cell.type") + coord_equal()  + NoLegend() 
plot_grid(p1,p2, align = "v")
## rds ----
saveRDS(rpca.integrated_cginsTime.vitro_islet, "rds_manuscript/rpca.integrated_cginsTime.vitro_islet.rds")
rpca.integrated_cginsTime.vitro_islet <- readRDS("rds_manuscript/rpca.integrated_cginsTime.vitro_islet.rds")

## visualization for publications ----
# umap
# time course umap colored by cell types
p1 <- DimPlot(rpca.integrated_cginsTime.vitro_islet, reduction = "umap", 
              label = F, cols = orig.color, group.by = "orig.ident") + 
  NoLegend() + 
  coord_fixed() + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), 
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
  )

ggsave(plot = p1, "figures_manuscript/fig5a.timeCourse.muRM.islet.celltype.umap.orig.jpeg", height = 4, width = 4, dpi = 500)

order <- c("Stem","Endocrine_1","Endocrine_2","GINS_precursors","delta","Delta_like","Epsilon_like","alpha","Alpha_like", "Beta_like", "beta")
p2 <- rpca.integrated_cginsTime.vitro_islet %>% DimPlot(reduction = "umap", cols = islet.timecourse.cell.type.color, group.by = "cell.type", order = rev(order)) + 
  NoLegend() + 
  coord_fixed() + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), 
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
  )
p1+p2
ggsave(plot = p2, "figures_manuscript/fig5b.timeCourse.muRM.islet.celltype.umap.jpeg", height = 4, width = 4, dpi = 500)

# INS and GAL on UMAP
DefaultAssay(rpca.integrated_cginsTime.vitro_islet) <- "RNA"
markers <- c("GAL", "INS")
marker.parameter.list <- list(c(0.4, 5.5), c(2.5, NA))
names(marker.parameter.list) <- markers
for (marker in names(marker.parameter.list)) {
  rpca.integrated_cginsTime.vitro_islet %>% 
    FeaturePlot(features = marker, 
                min.cutoff = marker.parameter.list[[marker]][1], 
                max.cutoff = marker.parameter.list[[marker]][2]) +
    NoLegend() + 
    coord_fixed() + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), 
          plot.title = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()
          )
  ggsave(paste0("figures_manuscript/Fig6a.Time.mucusRM.islet.",marker,".featureplot.jpeg"), dpi = 500, height = 6, width = 6)
}

## cell type markers in umap
Idents(rpca.integrated_cginsTime.vitro_islet) <- "cell.type"
DefaultAssay(rpca.integrated_cginsTime.vitro_islet) <- "RNA"
markers <- c("LGR5",
             "SOX9","MKI67","TFF2","TOP2A", # STEM
             "SOX4","HES6","ASCL1","CHGA","SLC38A11", # ENDOCRINE
             "SSTR2","GAL","EPHA4", # GINSpre
             "SST", "HHEX", # Delta,
             "INS")
marker.parameter.list <- list(c(0, 0.8), 
                              c(0, 2.8), c(0, 2.5), c(1.5, NA), c(NA, NA), 
                              c(1, 5), c(0.5, 5), c(0, 4), c(1.5, 6.5), c(0, 3),
                              c(0.5, 3.5), c(0.5, 6), c(0, 3),
                              c(1.5,NA), c(NA,NA),
                              c(2.9, NA)
)
names(marker.parameter.list) <- markers
for (marker in names(marker.parameter.list)) {
  FeaturePlot(rpca.integrated_cginsTime.vitro_islet, 
              features = marker, 
              min.cutoff = marker.parameter.list[[marker]][1], 
              max.cutoff = marker.parameter.list[[marker]][2]) +
    NoLegend() + 
    coord_fixed() + ``
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), 
          plot.title = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()
    )
  ggsave(paste0("figures_manuscript/ExtendedDataFig8.Time.mucusRM.wBeta.",marker,".featureplot.png"), dpi = 400, height = 2.4, width = 2.4, bg = "white")
}

Idents(rpca.integrated_cginsTime.vitro_islet) <- "cell.type"
## remove alpha and delta cells
rpca.integrated_cginsTime.vitro_isletbeta <- subset(rpca.integrated_cginsTime.vitro_islet, idents=c("alpha","delta"), invert=T)
## processing ----
DefaultAssay(rpca.integrated_cginsTime.vitro_isletbeta) <- "RNA"
rpca.integrated_cginsTime.vitro_isletbeta <- rpca.integrated_cginsTime.vitro_isletbeta %>% 
  NormalizeData(normalization.method = "LogNormalize") %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(vars.to.regress = c("percent.mt", "S.Score", "G2M.Score") ) %>% 
  RunPCA()

DimPlot(rpca.integrated_cginsTime.vitro_isletbeta, reduction = "pca" ,group.by = "cell.type", cols = my36colors) + coord_fixed()


## rds ----
rpca.integrated_cginsTime.vitro_isletbeta %>% saveRDS("rds_manuscript/rpca.integrated_cginsTime.vitro_isletbeta.rds")


## monocle trajectory analysis ---------------------------------------------------------------
#Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(rpca.integrated_cginsTime.vitro_isletbeta@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = rpca.integrated_cginsTime.vitro_isletbeta@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
# construct monocle cds
rpca.integrated_cginsTime.vitro_isletbeta.cds <- newCellDataSet(data,
                                                                phenoData = pd,
                                                                featureData = fd,
                                                                lowerDetectionLimit = 0.5,
                                                                expressionFamily = negbinomial.size())
# estimate size factors and dispersion
rpca.integrated_cginsTime.vitro_isletbeta.cds <- 
  rpca.integrated_cginsTime.vitro_isletbeta.cds %>% 
  estimateSizeFactors() %>% 
  estimateDispersions()
# cell/gene statistics and filtering
rpca.integrated_cginsTime.vitro_isletbeta.cds <- 
  rpca.integrated_cginsTime.vitro_isletbeta.cds %>% 
  detectGenes(min_expr = 0.1)

# ordering based on genes that differ between cell.type
Idents(rpca.integrated_cginsTime.vitro_isletbeta) <- "cell.type"
rpca.integrated_cginsTime.vitro_isletbeta.markers <- FindAllMarkers(rpca.integrated_cginsTime.vitro_isletbeta, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3)
table(rpca.integrated_cginsTime.vitro_isletbeta.markers$cluster)
top.rpca.integrated_cginsTime.vitro_isletbeta.markers <- rpca.integrated_cginsTime.vitro_isletbeta.markers %>% 
  group_by(cluster) %>% top_n(n = -100, wt = p_val_adj) %>% 
  top_n(n = 100, wt = avg_log2FC)
ordering_genes <- top.rpca.integrated_cginsTime.vitro_isletbeta.markers$gene
rpca.integrated_cginsTime.vitro_isletbeta.cds <- 
  rpca.integrated_cginsTime.vitro_isletbeta.cds %>% 
  setOrderingFilter(ordering_genes)
plot_ordering_genes(rpca.integrated_cginsTime.vitro_isletbeta.cds)

# reduce data dimensionality
rpca.integrated_cginsTime.vitro_isletbeta.cds <- 
  rpca.integrated_cginsTime.vitro_isletbeta.cds %>% 
  reduceDimension(max_components = 2,
                  method = 'DDRTree')

# order cells along the trajectory
rpca.integrated_cginsTime.vitro_isletbeta.cds <- 
  rpca.integrated_cginsTime.vitro_isletbeta.cds %>% 
  orderCells(reverse = TRUE)

# export to RDS
saveRDS(rpca.integrated_cginsTime.vitro_isletbeta.cds, "rds_manuscript/rpca.integrated_cginsTime.vitro_isletbeta.cds.rds")
#rpca.integrated_cginsTime.vitro_isletbeta.cds <- readRDS("rds_manuscript/rpca.integrated_cginsTime.vitro_isletbeta.cds.rds")

# visual
p1 <- plot_cell_trajectory(rpca.integrated_cginsTime.vitro_isletbeta.cds, color_by = "source") + 
  scale_color_manual(values = my36colors) +
  theme(legend.position = "right")
p2 <- plot_cell_trajectory(rpca.integrated_cginsTime.vitro_isletbeta.cds, color_by = "State") + 
  scale_color_manual(values = my36colors) +
  theme(legend.position = "right")
p3 <- plot_cell_trajectory(rpca.integrated_cginsTime.vitro_isletbeta.cds, color_by = "cell.type") + 
  scale_color_manual(values = my36colors)+
  theme(legend.position = "right")
p4 <- plot_cell_trajectory(rpca.integrated_cginsTime.vitro_isletbeta.cds, color_by = "Pseudotime") + 
  theme(legend.position = "right")
cowplot::plot_grid(p1,p2,p3,p4, align = "vh")

#### differential expression along the pseudotime ------------------------------
# filter genes by number of cells expressed
expressed_genes <- 
  row.names(subset(fData(rpca.integrated_cginsTime.vitro_isletbeta.cds),
                   num_cells_expressed >= nrow(pData(rpca.integrated_cginsTime.vitro_isletbeta.cds))*0.02))
# get DEG of pseudotime
diff_test_res.pseudotime.with.beta <- differentialGeneTest(rpca.integrated_cginsTime.vitro_isletbeta.cds[expressed_genes,], 
                                                           fullModelFormulaStr = "~sm.ns(Pseudotime)")
write_csv(diff_test_res.pseudotime.with.beta, "markers/diff_test_res.pseudotime.with.beta.csv")
#diff_test_res.pseudotime.with.beta <- read_csv("markers/diff_test_res.pseudotime.with.beta.csv")
sig.DEG.pseudotime <- diff_test_res.pseudotime.with.beta %>% filter(qval < 0.001)


## pseudotime deg heatmap ----
deg.pseudotime.heatmap <- plot_pseudotime_heatmap(rpca.integrated_cginsTime.vitro_isletbeta.cds[sig.DEG.pseudotime$gene_short_name,],
                                                  cluster_rows = TRUE,
                                                  num_clusters = 7,
                                                  show_rownames = T,
                                                  hmcols = rev(hmcols),
                                                  return_heatmap = TRUE)

plot_grid(deg.pseudotime.heatmap$gtable)
ggsave("figures_manuscript/deg.pseudotime.heatmap.pdf", width = 8, height = 10)
saveRDS(deg.pseudotime.heatmap, "rds_manuscript/deg.pseudotime.heatmap.rds")
deg.pseudotime.heatmap <- readRDS("rds_manuscript/deg.pseudotime.heatmap.rds")

## re-order the heatmap ##
## get the tree
pseudotime.tree <- deg.pseudotime.heatmap$tree_row
## cut the tree into 7 clusters
lbl <- cutree(pseudotime.tree, 7)
## get the cluster and order information
pseudotime.order <- tibble(gene=names(lbl),cluster=lbl, order=pseudotime.tree$order) %>% arrange(cluster,order)
## re-order the cluster
pseudotime.order$cluster <- factor(pseudotime.order$cluster, levels = c(4,3,7,6,1,2,5))
pseudotime.order <- pseudotime.order %>% arrange(cluster) 
pseudotime.order$cluster %>% table()
pseudotime.order.list <- pseudotime.order %>% group_by(cluster) %>% group_split()
names(pseudotime.order.list) <- paste0('cluster',unique(pseudotime.order$cluster))

# create a vector of genes for heatmap 
## the order of some gene cluster need to be flipped vertically
sapply(pseudotime.order.list, nrow)
gene.order <- c(
  
  rev(pseudotime.order.list$cluster4$gene),
  rev(pseudotime.order.list$cluster3$gene),
  rev(pseudotime.order.list$cluster7$gene),
  rev(pseudotime.order.list$cluster6$gene),
  (pseudotime.order.list$cluster1$gene),
  (pseudotime.order.list$cluster2$gene),
  rev(pseudotime.order.list$cluster5$gene)
)

# reorder the cds and plot pseudotime heatmap
ordered.pseudotime.heatmap <- plot_pseudotime_heatmap(rpca.integrated_cginsTime.vitro_isletbeta.cds[gene.order,],
                                                      cluster_rows = FALSE,
                                                      cores = 1,
                                                      show_rownames = F,
                                                      hmcols = rev(hmcols),
                                                      return_heatmap = TRUE)

plot_grid(ordered.pseudotime.heatmap$gtable) 
ggsave(filename = "figures_manuscript/Fig5f.ordered.pseudotime.heatmap.withBeta.pdf", width = 1.8, height = 2)
ggsave(filename = "figures_manuscript/Fig5f.ordered.pseudotime.heatmap.withBeta.png", dpi = 500, width = 1.8, height = 2)
# save RDS of the heatmap file
saveRDS(ordered.pseudotime.heatmap, "rds_manuscript/ordered.pseudotime.heatmap.rds")
#ordered.pseudotime.heatmap <- readRDS("rds_manuscript/ordered.pseudotime.heatmap.rds")

## functional annotation of each cluster
# check the list
pseudotime.order.list %>% sapply(nrow)
enrichment.list <- pseudotime.order.list %>% 
  lapply(function(x){
    x <- x$gene
  })
for (i in seq_along(enrichment.list)) {
  write.csv(enrichment.list[[i]], paste0("markers/enrichment.list_",names(enrichment.list)[i],".csv"))
}
enrichment.list %>% sapply(length)
names(enrichment.list) <- paste0("cluster",1:length(enrichment.list))
## GO term analysis for each cluster
GO.list <- enrichment.list %>% 
  lapply(enrichGO, 
         keyType = "SYMBOL",
         OrgDb = org.Hs.eg.db, 
         ont = "BP", 
         pAdjustMethod = "BH")

## reduce the redundancy of GO terms
GO.list.simplified <- GO.list %>% lapply(simplify, cutoff=0.7)
GO.list.fitlered <- lapply(GO.list.simplified, 
                           function(x){
                             if(min(x@result$p.adjust) > 0.1){x <- NULL}
                             else{x <- x}
                           })
idx <- GO.list.fitlered %>% sapply(is.null)
GO.list.fitlered <- GO.list.fitlered[!idx]
### extract results
dir.create("functional_enrichment")
output_path <- "functional_enrichment/go.results.withBeta.txt"
GO.result <- 
  GO.list.fitlered %>% 
  lapply(function(x){
    x <- x@result %>% 
      arrange(p.adjust) 
  }) %>% 
  bind_rows(.id = "cluster") %>% 
  mutate(x=1) 
GO.result$order <- 1:nrow(GO.result)
idx <- GO.result$Description %>% duplicated()
GO.result$order[idx] <- NA
GO.result$GeneRatio <- GO.result$GeneRatio %>% parse_ratio()
GO.result %>% write.table(output_path,
                          sep = '\t',
                          quote = FALSE, row.names = TRUE)

#TFs
TFs <- read.csv("gene_set/TF_names.txt", header = F)
TFs <- TFs$V1
TFs.order <- gene.order[gene.order %in% TFs]
write.table(TFs.order, "markers/TFs.order.txt")
TFs.heatmap <- plot_pseudotime_heatmap(rpca.integrated_cginsTime.vitro_isletbeta.cds[TFs.order,],
                                       cluster_rows = T,
                                       cores = 1,
                                       num_clusters = 7,
                                       show_rownames = T,
                                       hmcols = rev(hmcols),
                                       return_heatmap = TRUE)
plot_grid(TFs.heatmap$gtable) 
ggsave(filename = "figures_manuscript/ExtendedDataFig8c.ordered.TFs.pseudotime.heatmap.withBeta2.png", dpi = 400, bg="white", width = 1.8, height = 2)
ggsave(filename = "figures_manuscript/ExtendedDataFig8c.ordered.TFs.pseudotime.heatmap.withBeta.pdf", width = 5, height = 45)


## re-order the heatmap ##
## get the tree
TFs.tree <- TFs.heatmap$tree_row
## cut the tree into 7 clusters
lbl <- cutree(TFs.tree, 7)
## get the cluster and order information
TFs.order <- tibble(gene=names(lbl),cluster=lbl, order=TFs.tree$order) %>% arrange(cluster,order)
## re-order the cluster
TFs.order$cluster <- factor(TFs.order$cluster, levels = c(1,2,4,5,3,6,7))
TFs.order <- TFs.order %>% arrange(cluster) 
TFs.order$cluster %>% table()
TFs.order.list <- TFs.order %>% group_by(cluster) %>% group_split()
names(TFs.order.list) <- paste0('cluster',unique(TFs.order$cluster))

# create a vector of genes for heatmap 
## the order of some gene cluster need to be flipped vertically
sapply(TFs.order.list, nrow)
gene.order <- c(
  
  rev(TFs.order.list$cluster1$gene),
  (TFs.order.list$cluster2$gene),
  (TFs.order.list$cluster4$gene),
  (TFs.order.list$cluster5$gene),
  rev(TFs.order.list$cluster3$gene),
  (TFs.order.list$cluster6$gene),
  (TFs.order.list$cluster7$gene)
)

# reorder the cds and plot pseudotime heatmap
ordered.TFs.heatmap <- plot_pseudotime_heatmap(rpca.integrated_cginsTime.vitro_isletbeta.cds[gene.order,],
                                               cluster_rows = FALSE,
                                               cores = 1,
                                               show_rownames = T,
                                               hmcols = rev(hmcols),
                                               return_heatmap = TRUE)
saveRDS(ordered.TFs.heatmap, "rds_manuscript/ordered.TFs.heatmap.rds")

plot_grid(ordered.TFs.heatmap$gtable) 
ggsave("figures_manuscript/ordered.tfs.heatmap.pdf", height = 40, width = 5)
ggsave("figures_manuscript/ordered.tfs.heatmap.jpeg", height = 8, width = 5, bg="white", dpi = 300)


## cluster marker visualization ----
## assign identities
Idents(rpca.integrated_cginsTime.vitro_isletbeta) <- "cell.type"
## find markers
rpca.integrated_cginsTime.vitro_isletbeta.markers <- FindAllMarkers(rpca.integrated_cginsTime.vitro_isletbeta, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3)

## get top markers
top.rpca.integrated_cginsTime.vitro_isletbeta.markers <- rpca.integrated_cginsTime.vitro_isletbeta.markers %>% group_by(cluster) %>% top_n(n = 12, wt = avg_log2FC)
top.rpca.integrated_cginsTime.vitro_isletbeta.markers$cluster <- top.rpca.integrated_cginsTime.vitro_isletbeta.markers$cluster %>% 
  factor(levels = c("Stem","Endocrine_1","Endocrine_2","GINS_precursors","Delta_like","Epsilon_like","Alpha_like","Beta_like", "beta"))
top.rpca.integrated_cginsTime.vitro_isletbeta.markers <- top.rpca.integrated_cginsTime.vitro_isletbeta.markers %>% arrange(cluster)

## tiled dotplot for markers
TileDotplot(rpca.integrated_cginsTime.vitro_isletbeta, 
            features = rev(unique(top.rpca.integrated_cginsTime.vitro_isletbeta.markers$gene)), 
            col.max = 1,
            scale.by = "size", 
            scale = T) +
  NoLegend()

## tiled dotplot for specific markers
order <- c("Stem","Endocrine_1","Endocrine_2","GINS_precursors","Delta_like","Epsilon_like","Alpha_like","Beta_like", "beta")
Idents(rpca.integrated_cginsTime.vitro_isletbeta) <- "cell.type"
rpca.integrated_cginsTime.vitro_isletbeta$cell.type <- rpca.integrated_cginsTime.vitro_isletbeta$cell.type %>% 
  factor(levels = rev(order))

cherry <- c(
  "INS","SCG5",'CPE',
  "GCG","TTR",
  'GHRL','HHEX','SST',
  'SSTR2','GAL',
  'CHGA','HES6','SOX4','LCN2','TFF3',
  'TFF2','TFF1','SPINK4','SOX9','MKI67'
)

TileDotplot(rpca.integrated_cginsTime.vitro_isletbeta, 
            features = rev(cherry), 
            col.max = 1,
            col.min = -1,
            scale = T, 
            idents = rev(order))  +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.title = element_blank(),
        legend.title = element_blank())

ggsave("figures_manuscript/timecourse.mucusRM.beta.markers.tiledotplot.pdf", width = 8, height = 4)

## cluster DEGs pseudotime visualization ----
DEGs <- FindAllMarkers(rpca.integrated_cginsTime.vitro_isletbeta, only.pos = FALSE, min.pct = 0.3, logfc.threshold = 0.3)
DEGs.filtered <- DEGs %>% filter(p_val_adj <0.001)
DEGs.unique <- DEGs.filtered$gene %>% unique()
cluster.deg.pseudotime.heatmap <- 
  plot_pseudotime_heatmap(rpca.integrated_cginsTime.vitro_isletbeta.cds[DEGs.unique,],
                          #hclust_method = "complete",
                          cluster_rows = TRUE,
                          num_clusters = 8,
                          show_rownames = T,
                          hmcols = rev(hmcols),
                          return_heatmap = TRUE)

plot_grid(cluster.deg.pseudotime.heatmap$gtable)
ggsave("figures_manuscript/cluster.deg.pseudotime.heatmap.pdf", width = 8, height = 10)
saveRDS(cluster.deg.pseudotime.heatmap, "rds_manuscript/cluster.deg.pseudotime.heatmap.rds")

## get the tree
pseudotime.tree <- cluster.deg.pseudotime.heatmap$tree_row
## cut the tree into clusters
lbl <- cutree(pseudotime.tree, 8)
## get the cluster and order information
pseudotime.order <- tibble(gene=names(lbl),cluster=lbl, order=pseudotime.tree$order) %>% arrange(cluster,order)
pseudotime.order.list <- pseudotime.order %>% group_by(cluster) %>% group_split()
names(pseudotime.order.list) <- paste0('cluster',unique(pseudotime.order$cluster))
for (i in seq_along(pseudotime.order.list)) {
  write.csv(pseudotime.order.list[[i]], paste0("markers/clusterDEG_",names(enrichment.list)[i],".csv"))
}

## GO term analysis for each cluster
enrichment.list <- pseudotime.order.list %>% 
  lapply(function(x){
    x <- x$gene
  })
GO.list <- enrichment.list %>% 
  lapply(enrichGO, 
         keyType = "SYMBOL",
         OrgDb = org.Hs.eg.db, 
         ont = "BP", 
         pAdjustMethod = "BH")

## reduce the redundancy of GO terms
GO.list.simplified <- GO.list %>% lapply(simplify, cutoff=0.7)
GO.list.fitlered <- lapply(GO.list.simplified, 
                           function(x){
                             if(min(x@result$p.adjust) > 0.1){x <- NULL}
                             else{x <- x}
                           })
idx <- GO.list.fitlered %>% sapply(is.null)
GO.list.fitlered <- GO.list.fitlered[!idx]


## pseudotime trajectory visualization ---------------------------------------
## get meta
cds.meta <- pData(rpca.integrated_cginsTime.vitro_isletbeta.cds) # orignal meta
cds.meta <- cds.meta %>% arrange(Pseudotime) # sort by pseudotime
cds.meta$Pseudotime.order <- 1:nrow(cds.meta) # store the order
number <- length(cds.meta$Pseudotime.order)/100 # bin the pseudotime
bin <- as.integer(cds.meta$Pseudotime.order/number) # bin the pseudotime
cds.meta$Pseudotime.bin <- bin # store the bin
cds.meta$Pseudotime.bin %>% table()
order <- c("Stem","Endocrine_1","Endocrine_2","GINS_precursors","Delta_like","Epsilon_like","Alpha_like","Beta_like", "beta")

cds.meta$cell.type <- cds.meta$cell.type %>% 
  factor(levels = (order))


## plot pseudotime heatmap of the markers  -----------------
# plot a density plot above the heatmap to show cell type distibution
dens.plot <- cds.meta %>% 
  ggplot(aes(x=Pseudotime.bin, fill=cell.type)) +
  geom_density(size=0.25, adjust=1.5,alpha =0.8)+
  scale_fill_manual(values = c(timecourse.cell.type.color, beta = "grey")) + # add transparency value
  #scale_color_manual(values = timecourse.cell.type.color) +
  theme_void()+
  theme(legend.position = "top",
        legend.title = element_blank(),
  )+ 
  NULL
dens.plot
ggsave("figures_manuscript/pseudotime.dens.plot.with.beta.pdf", height = 1, width = 4)

## add meta to seurat
rpca.integrated_cginsTime.vitro_isletbeta <- 
  AddMetaData(rpca.integrated_cginsTime.vitro_isletbeta,
              cds.meta[,c("Pseudotime","Pseudotime.bin")])

## a vlnplot to show combined transgene pdx1/mafa
rpca.integrated_cginsTime.vitro_isletbeta.temp <- rpca.integrated_cginsTime.vitro_isletbeta
rpca.integrated_cginsTime.vitro_isletbeta.temp[['Transgenes']] <- rpca.integrated_cginsTime.vitro_isletbeta.temp@assays$RNA@data["PDX1MUS",] + rpca.integrated_cginsTime.vitro_isletbeta.temp@assays$RNA@data["MAFAMUS",]
rpca.integrated_cginsTime.vitro_isletbeta.temp$cell.type <- 
  rpca.integrated_cginsTime.vitro_isletbeta.temp$cell.type %>% 
  factor(levels = (order))
Idents(rpca.integrated_cginsTime.vitro_isletbeta.temp) <- "cell.type"
p1 <- VlnPlot(rpca.integrated_cginsTime.vitro_isletbeta.temp, "Transgenes", cols = islet.timecourse.cell.type.color, pt.size = 0) +
  theme_classic() + 
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        title = element_blank(),
        axis.ticks.x = element_blank()) +
  NoLegend() +
  NULL
p1
## scatter plot of pseudotime over celltype
p2 <- cds.meta %>% 
  ggplot() + 
  geom_jitter(aes(x=cell.type, y=Pseudotime, color = cell.type), size=0.05, alpha=0.5) + 
  geom_smooth(aes(x=as.numeric(cell.type), y=Pseudotime), 
              method = "loess", se = FALSE, color = "black", size = 0.5) +
  scale_color_manual(values = islet.timecourse.cell.type.color) + # add transparency valuescale_color_manual(values = timecourse.cell.type.color) + # add transparency value
  theme_classic() + 
  NoLegend() +
  theme(axis.title = element_blank(),
        axis.text = element_blank()) +
  NULL
p2
ggsave("figures_manuscript/Fig5d.scatter.pseudotime.png", dpi = 400, width = 4, height = 3)
## combined the two plots above
p3 <- cowplot::plot_grid(p2,p1, ncol = 1, 
                         align = "v",
                         axis = "tblr",
                         rel_heights = c(4,1))
p3
ggsave("figures_manuscript/Fig5d.gins.islet_peudotime_transgene.pdf", width = 6, height = 8)
ggsave("figures_manuscript/Fig5d.gins.islet_peudotime_transgene.jpeg", dpi = 500, width = 6, height = 8)


## gene expression along pseudotime visualization -----------------
my_genes <-  c("MKI67",
               "SOX4", "HES6", "CHGA",
               "GAL",
               "G6PC2","ERO1B","CPE")

for (gene in my_genes) {
  cds_subset <- rpca.integrated_cginsTime.vitro_isletbeta.cds[gene,]
  plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime") + 
    scale_size_manual( values = 3 ) +
    scale_color_gradientn(colors = grad)  +
    theme(#axis.title = element_blank(),
      #axis.text.x = element_blank(),
      #strip.text = element_text(size = 0),
      strip.background = element_blank(),
      legend.position = "none")
  ggsave(paste0("figures_manuscript/ExtendedDataFig8b_",gene,"_in_pseudotime_wBeta.jpeg"), width = 2.4, height = 1.2, dpi = 400)
}

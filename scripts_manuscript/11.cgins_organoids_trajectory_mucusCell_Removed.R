# This script pre-processed the corpus#6 gins organoids differntiation trajectory data without mucus-secreting cells.

## loading data and remove mucus cells #### 
cgins.vitro.time <- readRDS("rds_manuscript/cgins.vitro.time.rds")
Idents(cgins.vitro.time) <- "cell.type"
table(cgins.vitro.time$cell.type)
cgins.vitro.time.muRM <- cgins.vitro.time %>% subset(idents = "Mucoid", invert = TRUE)

## processing ----
cgins.vitro.time.muRM <- cgins.vitro.time.muRM %>% 
  NormalizeData(normalization.method = "LogNormalize") %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(vars.to.regress = c("percent.mt", "S.Score", "G2M.Score") ) %>% 
  RunPCA()

DimHeatmap(cgins.vitro.time.muRM, dims = 1:30, cells = 400, balanced = TRUE, ncol = 3)
ElbowPlot(cgins.vitro.time.muRM, ndims = 50)
# umap and tsne
dimensions <- 1:16
cgins.vitro.time.muRM <- 
  cgins.vitro.time.muRM %>% 
  RunUMAP(dims = dimensions)
p1 <- DimPlot(cgins.vitro.time.muRM, group.by = "Phase", reduction = "umap") + coord_fixed()
p2 <- DimPlot(cgins.vitro.time.muRM, group.by = "orig.ident", reduction = "umap") + coord_fixed()
p3 <- DimPlot(cgins.vitro.time.muRM, group.by = "cell.type", reduction = "umap", cols = my9colors) + coord_fixed()
p1+p2+p3

## rds ----
# cgins.vitro.time.muRM %>% saveRDS("rds_manuscript/cgins.vitro.time.muRM.rds")
cgins.vitro.time.muRM <- readRDS("rds_manuscript/cgins.vitro.time.muRM.rds")

## monocle3 trajectory ---------------------------------------------------------------------------
library(monocle3)
# convert to monocle3 object
cgins.vitro.time.muRM.mono3cds <- as.cell_data_set(cgins.vitro.time.muRM)
# cluster partitions
cgins.vitro.time.muRM.mono3cds <- monocle3::cluster_cells(cgins.vitro.time.muRM.mono3cds, reduction_method = "UMAP", k=20)
cgins.vitro.time.muRM.mono3cds <- monocle3::learn_graph(cgins.vitro.time.muRM.mono3cds, use_partition = TRUE)
cgins.vitro.time.muRM.mono3cds <- monocle3::order_cells(cgins.vitro.time.muRM.mono3cds, reduction_method = "UMAP")
p1 <- plot_cells(cgins.vitro.time.muRM.mono3cds, show_trajectory_graph = TRUE)
p2 <- plot_cells(cgins.vitro.time.muRM.mono3cds, color_cells_by = "partition", show_trajectory_graph = FALSE, label_cell_groups = TRUE)
wrap_plots(p1, p2)
monocle3::plot_cells(cgins.vitro.time.muRM.mono3cds, color_cells_by = "pseudotime")
cgins.vitro.time.muRM.mono3cds.p2 <- 
  cgins.vitro.time.muRM.mono3cds %>% 
  as.Seurat() %>% 
  subset(monocle3_partitions == 2) %>% 
  as.cell_data_set()

cgins.vitro.time.muRM.mono3cds.p2 <- monocle3::cluster_cells(cgins.vitro.time.muRM.mono3cds.p2, reduction_method = "UMAP", k=20)
cgins.vitro.time.muRM.mono3cds.p2 <- monocle3::learn_graph(cgins.vitro.time.muRM.mono3cds.p2, use_partition = TRUE)
cgins.vitro.time.muRM.mono3cds.p2 <- monocle3::order_cells(cgins.vitro.time.muRM.mono3cds.p2, reduction_method = "UMAP")
monocle3::plot_cells(cgins.vitro.time.muRM.mono3cds.p2, 
                     color_cells_by = "pseudotime", 
                     label_cell_groups = FALSE, 
                     label_leaves = FALSE, 
                     label_branch_points = FALSE,
                     label_roots = FALSE) + coord_fixed()
ggsave("figures_manuscript/ExtendedDataFig9c.mono3.trajectory.jpeg", dpi = 400, width = 4, height = 3)
# get pseudotime data
pseudotime.data <- pseudotime(cgins.vitro.time.muRM.mono3cds.p2, reduction_method = "UMAP")
# convert to monocle3 object
cgins.vitro.time.muRM.p2 <- as.Seurat(cgins.vitro.time.muRM.mono3cds.p2)
cgins.vitro.time.muRM.p2 <- AddMetaData(cgins.vitro.time.muRM.p2, metadata=pseudotime.data, col.name = "monocle3.pseudotime")
# get count data of genes
data <- cgins.vitro.time.muRM.p2@assays$RNA@data[c("INS","GAL"),] %>% as.matrix()  %>% t() %>% as_data_frame() %>% as_tibble()
# add pseudotime
data$pseudotime <- cgins.vitro.time.muRM.p2$monocle3.pseudotime
data %>% ggplot(aes(x = log10(GAL+1), y = log10(INS+1), color = pseudotime)) +
  geom_jitter()
## umap visualization --------------------------------------------------
cgins.vitro.time.muRM$cell.type <- cgins.vitro.time.muRM$cell.type %>% 
  factor(levels = c("Stem","Endocrine_1","Endocrine_2",
                    "GINS_precursors","Delta_like","Epsilon_like",
                    "Alpha_like","Beta_like"))
# time course umap colored by cell types
cgins.vitro.time.muRM %>% DimPlot(reduction = "umap", cols = timecourse.cell.type.color, group.by = "cell.type") + 
  coord_fixed() + theme_void() + NoLegend()
ggsave("figures_manuscript/ExtendedDataFig9b.timeCourse.muRM.celltype.umap.jpeg", height = 4, width = 4, dpi = 400)

## cluster marker visualization ----
## assign identities
Idents(cgins.vitro.time.muRM) <- "cell.type"
## find markers
cgins.vitro.time.muRM.markers <- FindAllMarkers(cgins.vitro.time.muRM, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3)

## get top markers
top.cgins.vitro.time.muRM.markers <- cgins.vitro.time.muRM.markers %>% group_by(cluster) %>% top_n(n = 12, wt = avg_log2FC)
DoHeatmap(cgins.vitro.time.muRM,features = top.cgins.vitro.time.muRM.markers$gene)
top.cgins.vitro.time.muRM.markers$cluster <- top.cgins.vitro.time.muRM.markers$cluster %>% 
  factor(levels = c("Stem","Endocrine_1","Endocrine_2","GINS_precursors","Delta_like","Epsilon_like","Alpha_like","Beta_like"))
top.cgins.vitro.time.muRM.markers <- top.cgins.vitro.time.muRM.markers %>% arrange(cluster)

## tiled dotplot for markers
## tiled dotplot for specific markers
order <- c("Stem","Endocrine_1","Endocrine_2","GINS_precursors","Delta_like","Epsilon_like","Alpha_like","Beta_like")
Idents(cgins.vitro.time.muRM) <- "cell.type"
cgins.vitro.time.muRM$cell.type <- cgins.vitro.time.muRM$cell.type %>% 
  factor(levels = rev(order))

cherry <- c(
  "INS","SCG5",'CPE',
  "GCG","TTR",
  'GHRL','HHEX','SST',
  'SSTR2','GAL',
  'CHGA','HES6','SOX4','LCN2','TFF3',
  'TFF2','TFF1','SPINK4','SOX9','MKI67'
)

TileDotplot(cgins.vitro.time.muRM, 
            features = rev(cherry), 
            col.max = 1,
            col.min = -1,
            scale = T, 
            idents = rev(order))  +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.title = element_blank(),
        legend.title = element_blank())

ggsave("figures_manuscript/Fig5c.timecourse.mucusRM.markers.tiledotplot.pdf", width = 8, height = 4)

## plot markers on umap ----
## umap of gins precursor markers

## regulons ----
## Please run Regulon analysis first!
DefaultAssay(cgins.vitro.time.muRM) <- "AUC"
all.auc.feature <- rownames(cgins.vitro.time.muRM)
cgins.vitro.time.muRM <- cgins.vitro.time.muRM %>% ScaleData(feature = all.auc.feature, verbose = FALSE) 
cgins.vitro.time.muRM <- FindVariableFeatures(cgins.vitro.time.muRM)
var.auc.feature <- cgins.vitro.time.muRM@assays$AUC@var.features
cgins.vitro.time.muRM <- cgins.vitro.time.muRM %>% RunUMAP(reduction.name="UMAP_auc", assay = "AUC", features=var.auc.feature)
DimPlot(cgins.vitro.time.muRM, reduction = "UMAP_auc", cols = timecourse.cell.type.color) + 
  coord_fixed() + 
  theme_void() +
  NoLegend()
ggsave("figures_manuscript/muRM.auc.umap.jpeg", dpi = 400, width = 4, height = 3)
Idents(cgins.vitro.time.muRM) <- "cell.type"
cgins.vitro.time.muRM$cell.type <- cgins.vitro.time.muRM$cell.type %>% 
  factor(levels = c("Stem","Endocrine_1","Endocrine_2",
                    "GINS_precursors","Delta_like","Epsilon_like",
                    "Alpha_like","Beta_like"))

## differential regulons
AUC.marker <- FindAllMarkers(cgins.vitro.time.muRM, assay = "AUC", logfc.threshold = 0.06, min.pct = 0.2)
## reserve extended only to remove redundancy
regulon <- AUC.marker$gene %>% unique()
sce.subset <- cgins.vitro.time.muRM %>% subset(downsample=80)
cells.subset <- sce.subset@meta.data %>% arrange(cell.type, Pseudotime) %>% rownames()
cellInfo <- sce.subset@meta.data[,c("cell.type","Pseudotime")]

## regulon auc heatmap and umap 
DefaultAssay(cgins.vitro.time.muRM) <- "AUC"

## get matrix to plot heatmap
auc <- cgins.vitro.time.muRM[AUC.marker$gene, cells.subset]@assays$AUC@scale.data
auc <- auc %>% as.data.frame()
p1 <- auc[,cells.subset] %>% 
  pheatmap(breaks = seq(-1, 1, length.out = 100),
           annotation_col= cellInfo[,c(1,2)],
           annotation_colors = ann_colors,
           clustering_method = "ward.D2",
           show_colnames = F,
           cluster_cols = F,
           fontsize = 4,
           legend = F,
           annotation_legend = FALSE
  )
plot_grid(p1$gtable)
ggsave("figures_manuscript/ExtendedDataFig9a.differential.AUC.heatmap.jpeg", dpi = 400, width = 8, height = 6, bg="white")

# final rds ----
cgins.vitro.time.muRM %>% saveRDS("rds_manuscript/cgins.vitro.time.muRM.rds")

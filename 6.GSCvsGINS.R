# This script processes the vivo cgins data into normalized seurat cgins.vitro.time.gsc.gins
## loading data and remove mucus cells #### 
cgins.vitro.time <- readRDS("rds_manuscript/cgins.vitro.time.rds")
Idents(cgins.vitro.time) <- "orig.ident"
cgins.vitro.time.gsc.gins <- cgins.vitro.time %>% subset(idents = c("GSC","Late_GINS"))
Idents(cgins.vitro.time.gsc.gins) <- "cell.type"
cgins.vitro.time.gsc.gins <- 
  cgins.vitro.time.gsc.gins %>% 
  subset(idents = c("Stem","Mucoid",
                    "Delta_like","Epsilon_like",
                    "Alpha_like","Beta_like"))
## processing ----
cgins.vitro.time.gsc.gins <- cgins.vitro.time.gsc.gins %>% 
  NormalizeData(normalization.method = "LogNormalize") %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
  ScaleData(vars.to.regress = c("percent.mt", "S.Score", "G2M.Score") ) %>% 
  RunPCA()
DimHeatmap(cgins.vitro.time.gsc.gins, dims = 1:18, cells = 400, balanced = TRUE, ncol = 3)
ElbowPlot(cgins.vitro.time.gsc.gins, ndims = 50)
# umap and tsne
dimensions <- 1:15
cgins.vitro.time.gsc.gins <- 
  cgins.vitro.time.gsc.gins %>% 
  RunUMAP(dims = dimensions)
p1 <- DimPlot(cgins.vitro.time.gsc.gins, group.by = "Phase", reduction = "umap") + coord_fixed()
p2 <- DimPlot(cgins.vitro.time.gsc.gins, group.by = "orig.ident", reduction = "umap") + coord_fixed()
p3 <- DimPlot(cgins.vitro.time.gsc.gins, group.by = "cell.type", reduction = "umap", cols = my9colors) + coord_fixed()
p1+p2+p3

## clustering analysis ---------------------------------------------------------
cgins.vitro.time.gsc.gins <- cgins.vitro.time.gsc.gins %>%   
  FindNeighbors(dims = dimensions)
for (resolution in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5)) {
  cgins.vitro.time.gsc.gins <- cgins.vitro.time.gsc.gins %>%  
    FindClusters(resolution = resolution)}
p1 <- DimPlot(cgins.vitro.time.gsc.gins, reduction = "umap", label = T, cols = my36colors, group.by = "RNA_snn_res.0.1") + coord_equal() + NoLegend()
p2 <- DimPlot(cgins.vitro.time.gsc.gins, reduction = "umap", label = T, cols = my36colors, group.by = "RNA_snn_res.0.2") + coord_equal() + NoLegend()
p3 <- DimPlot(cgins.vitro.time.gsc.gins, reduction = "umap", label = T, cols = my36colors, group.by = "RNA_snn_res.0.3") + coord_equal() + NoLegend()
p4 <- DimPlot(cgins.vitro.time.gsc.gins, reduction = "umap", label = T, cols = my36colors, group.by = "RNA_snn_res.0.4") + coord_equal() + NoLegend()
p5 <- DimPlot(cgins.vitro.time.gsc.gins, reduction = "umap", label = T, cols = my36colors, group.by = "RNA_snn_res.0.5") + coord_equal() + NoLegend()
p6 <- DimPlot(cgins.vitro.time.gsc.gins, reduction = "umap", label = T, cols = my36colors, group.by = "RNA_snn_res.0.7") + coord_equal() + NoLegend()
p7 <- DimPlot(cgins.vitro.time.gsc.gins, reduction = "umap", label = T, cols = my36colors, group.by = "RNA_snn_res.1") + coord_equal() + NoLegend()
p8 <- DimPlot(cgins.vitro.time.gsc.gins, reduction = "umap", label = T, cols = my36colors, group.by = "RNA_snn_res.1.5") + coord_equal() + NoLegend()
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,nrow = 2)
## rds ----
cgins.vitro.time.gsc.gins %>% saveRDS("rds_manuscript/cgins.vitro.time.gsc.gins.rds")
cgins.vitro.time.gsc.gins <- readRDS("rds_manuscript/cgins.vitro.time.gsc.gins.rds")

## markers heatmap -----
Idents(cgins.vitro.time.gsc.gins) <- "cell.type"
cgins.vitro.time.gsc.gins.markers <- FindAllMarkers(cgins.vitro.time.gsc.gins, only.pos = TRUE)
top.cgins.vitro.time.gsc.gins.markers <- 
  cgins.vitro.time.gsc.gins.markers %>% 
  filter(p_val_adj < 0.01) %>% 
  group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
DoHeatmap(cgins.vitro.time.gsc.gins, 
          features = c(top.cgins.vitro.time.gsc.gins.markers$gene),
          #disp.min = -1.2, 
          #disp.max = 2
)# +  scale_fill_viridis(option = "B")

## UMAP ----
p1 <- DimPlot(cgins.vitro.time.gsc.gins, group.by = "Phase", reduction = "umap") + coord_fixed()
p2 <- DimPlot(cgins.vitro.time.gsc.gins, group.by = "orig.ident", reduction = "umap", cols = my5colors[c(4,1)]) + coord_fixed()
p3 <- DimPlot(cgins.vitro.time.gsc.gins, group.by = "cell.type", reduction = "umap", cols = timecourseGSC.GINS.cell.type.color) + coord_fixed()
ggsave(plot = p2 + NoLegend(), "figures_manuscript/ExtendedDataFig4a.top.jpeg", dpi = 400, height = 4, width = 4)
ggsave(plot = p3 + NoLegend(), "figures_manuscript/ExtendedDataFig4a.middle.jpeg", dpi = 400, height = 4, width = 4)

## proliferation vlnplot ----
cgins.vitro.time.gsc.gins$cell.type <- 
  cgins.vitro.time.gsc.gins$cell.type %>% 
  factor(levels = c(
    "Stem",
    "Mucoid",
    "Delta_like",
    "Epsilon_like",
    "Alpha_like",
    "Beta_like"
  ))
Idents(cgins.vitro.time.gsc.gins) <- "cell.type"
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
idx <- c(s.genes, g2m.genes) %>% str_starts("MCM")
genes.to.plot <- c("MKI67","PCNA","MCM4","MCM5","MCM6","MCM7")

VlnPlot(cgins.vitro.time.gsc.gins, 
        genes.to.plot, 
        stack = T, flip = T, 
        cols = timecourseGSC.GINS.cell.type.color, 
        fill.by = "ident" ) + 
  theme(axis.title = element_blank()) +
  NoLegend()
ggsave("figures_manuscript/ExtendedDataFig4e.proliferationMarkers.pdf", width = 2, height = 3)

## cell type markers tileDotplot -----
source("helper_functions/TileDotplot.R")
cgins.vitro.time.gsc.gins$cell.type <- 
  cgins.vitro.time.gsc.gins$cell.type %>% 
  factor(levels = rev(c("Beta_like",
                        "Alpha_like",
                        "Epsilon_like",
                        "Delta_like",
                        "Mucoid",
                        "Stem")))
Idents(cgins.vitro.time.gsc.gins) <- "cell.type"

markers.to.plot <- top.cgins.vitro.time.gsc.gins.markers %>% 
  arrange(desc(cluster))
TileDotplot(cgins.vitro.time.gsc.gins, 
            features = unique(markers.to.plot$gene), 
            col.max = 1.8,
            scale.by = "size", 
            scale = T) +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_text(size = 20),
    axis.text.y = element_blank(),
    legend.position = "bottom"
  )

ggsave("figures_manuscript/ExtendedDataFig4b.markers.tiledotplot.pdf", 
       width = length(unique(markers.to.plot$gene))/3.5, height = 5)

## cell type markers in umap ----
markers <- c("LGR5","MKI67","SOX9","AXIN2", # STEM
             "MUC5AC","MUC1", # MUCUS
             "SST","HHEX","GHRL","GCG", # ENDOCRINE
             "TTR","INS","SCG2","ERO1B" # ENDOCRINE
)
marker.parameter.list <- list(c(0,0.8), c(NA,NA), c(NA,NA), c(NA,0.8),
                              c(2.5,NA), c(2.6,NA), 
                              c(2,NA), c(NA,2.5), c(2, NA), c(0.5,3),
                              c(3, NA), c(2.8, NA), c(NA, NA), c(NA, NA))
names(marker.parameter.list) <- markers

for (marker in names(marker.parameter.list)) {
  FeaturePlot(cgins.vitro.time.gsc.gins, marker, 
              min.cutoff = marker.parameter.list[[marker]][1], 
              max.cutoff = marker.parameter.list[[marker]][2]) +
    NoLegend() + 
    coord_fixed() + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), 
          axis.title = element_blank(), 
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          title = element_text(size =  0))

  ggsave(paste0("figures_manuscript/ExtendedDataFig4a.bottom.GSCandGINS.",marker,".featureplot.pdf"), height = 3, width = 3)
}



# This script pre-processed the corpus#6 gins organoids


#### cgins.vitro loading data --------------------------------------------------------------
# read in the 10X matrices from your folder locations (you'll need to change)
cgins.data <- Read10X_h5(filename = 'data/CGINS_3/filtered_feature_bc_matrix.h5')

# create the seurat objects
cgins.vitro <- CreateSeuratObject(cgins.data, min.features=300, min.cells=10, project="C-GINSC")
cgins.data <- NULL
cgins.gene <- rownames(cgins.data)
#### cgins.vitro QC ------------------------------------------------------------------------
# create mito content
cgins.vitro[["percent.mt"]] <- PercentageFeatureSet(cgins.vitro, pattern = "^MT-")

# QC scatterplots 
p1 <- FeatureScatter(cgins.vitro, feature1 = "nCount_RNA", feature2 = "percent.mt")
p2 <- FeatureScatter(cgins.vitro, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p1+p2

# QC plots and make a decision as to cutoffs/filters
VlnPlot(cgins.vitro, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3,pt.size = 0.01, group.by = "orig.ident")

# quantile
quantile(cgins.vitro$nFeature_RNA, 0.03)
quantile(cgins.vitro$nFeature_RNA, 0.97)
quantile(cgins.vitro$nCount_RNA, 0.03)
quantile(cgins.vitro$nCount_RNA, 0.97)

# filter seurat objects - at this stage rerun the previous plts to see how they look now
cgins.vitro <- subset(cgins.vitro, subset = nFeature_RNA > 500 &  nFeature_RNA < 6000 & percent.mt < 15 & nCount_RNA < 24000 & nCount_RNA > 2000)

# pre-process data using standard workflow for QC
cgins.vitro <- 
  cgins.vitro %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA()
DimHeatmap(cgins.vitro, dims = 1:24, cells = 500, balanced = TRUE)
ElbowPlot(cgins.vitro, ndims = 50)
dimensions <- 1:18
cgins.vitro <- 
  cgins.vitro %>% 
  RunUMAP(dims = dimensions) %>% 
  RunTSNE(dims = dimensions) %>% 
  FindNeighbors(dims = dimensions) %>% 
  FindClusters(resolution = 0.5)
DimPlot(cgins.vitro, reduction = "umap", label = T, cols = my36colors) + NoLegend()
VlnPlot(cgins.vitro, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3,pt.size = 0.01)

# find markers for QC 
## cluster markers
cgins.vitro.markers <- FindAllMarkers(cgins.vitro, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3)
top.cgins.vitro.markers <- cgins.vitro.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

## heatmap for markers 
DefaultAssay(cgins.vitro) <- "RNA"
DoHeatmap(cgins.vitro, top.cgins.vitro.markers$gene, angle = 90) 

## is cluster-1 doublet? it has all the markers from different clusters and high UMI. run doublet finder

#### cgins.vitro doublets finder ---------------------------------------------------------------------------------------
# pK Identification (no ground-truth) 
sweep.res.list <- paramSweep_v3(cgins.vitro, PCs = dimensions, sct = F)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

# Homotypic Doublet Proportion Estimate 
DoubletRate <- 0.2                
homotypic.prop <- modelHomotypic(cgins.vitro$RNA_snn_res.0.5) 
nExp_poi <- round(DoubletRate*ncol(cgins.vitro)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder with varying classification stringencies
cgins.vitro <- doubletFinder_v3(cgins.vitro, PCs = dimensions, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = F, sct = F)

# visuzalization
cgins.vitro@meta.data %>% colnames()
p1 <- DimPlot(cgins.vitro, reduction = "tsne", group.by = "DF.classifications_0.25_0.26_1271")
p2 <- DimPlot(cgins.vitro, reduction = "tsne", group.by = "RNA_snn_res.0.5")
p1+p2
VlnPlot(cgins.vitro, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3,pt.size = 0.01, group.by = "DF.classifications_0.25_0.26_1271" )

# remove the doublets
Idents(cgins.vitro) <- "DF.classifications_0.25_0.26_1271"
cgins.vitro <- subset(cgins.vitro, ident = "Singlet")

#### cgins.vitro re-processing after doublet removal ---------------------------------------
DefaultAssay(cgins.vitro) <- "RNA"
cgins.vitro <- cgins.vitro %>% 
  # log-normalization 
  NormalizeData(normalization.method = "LogNormalize") %>% 
  # variable features
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  # scaling
  ScaleData(vars.to.regress = "percent.mt") %>% 
  # PCA
  RunPCA()

DimHeatmap(cgins.vitro, dims = 1:27, cells = 500, balanced = TRUE)
ElbowPlot(cgins.vitro, ndims = 50)
# umap and tsne
DimHeatmap(cgins.vitro, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(cgins.vitro, ndims = 50)
# umap and tsne

#dimensions <- 1:9
dimensions <- 1:13
cgins.vitro <- 
  cgins.vitro %>% 
  RunUMAP(dims = dimensions) %>% 
  RunTSNE(dims = dimensions) %>% 
  FindNeighbors(dims = dimensions)
p1 <- DimPlot(cgins.vitro, group.by = "cell.type", cols = my5colors) + coord_fixed()
p2 <- DimPlot(cgins.vitro, reduction = "tsne", group.by = "cell.type", cols = my5colors) + coord_fixed()
p1+p2

DimPlot(cgins.vitro, reduction = "tsne", group.by = "integrated.ident", cols = my5colors) + coord_fixed()
dimensions <- 1:6
cgins.vitro <- 
  cgins.vitro %>% 
  RunUMAP(dims = dimensions)
FeaturePlot(cgins.vitro, reduction = "tsne", features = c("INS","GCG","GHRL","SST","GAL"))

#### cgins.vitro clustering analysis ---------------------------------------------------------
for (resolution in c(0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.2)) {
  cgins.vitro <- cgins.vitro %>%  
    FindClusters(resolution = resolution)}
DimPlot(cgins.vitro, reduction = "tsne", label = T, cols = my36colors) + coord_equal() + NoLegend()
FeaturePlot(cgins.vitro, c("GCG","INS"), reduction = "tsne")

p1 <- DimPlot(cgins.vitro, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.2") + coord_equal() + NoLegend()
p2 <- DimPlot(cgins.vitro, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.3") + coord_equal() + NoLegend()
p3 <- DimPlot(cgins.vitro, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.4") + coord_equal() + NoLegend()
p4 <- DimPlot(cgins.vitro, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.5") + coord_equal() + NoLegend()
p5 <- DimPlot(cgins.vitro, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.7") + coord_equal() + NoLegend()
p6 <- DimPlot(cgins.vitro, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.1") + coord_equal() + NoLegend()
p7 <- DimPlot(cgins.vitro, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.1.2") + coord_equal() + NoLegend()
plot_grid(p1,p2,p3,p4,p5,p6,p7,ncol = 4)


# select resolution
sel.cluster <- "RNA_snn_res.0.4"
cgins.vitro <- cgins.vitro %>% SetIdent(value = sel.cluster)
DimPlot(cgins.vitro, reduction = "tsne", label = T, cols = my9colors) + coord_equal()
Idents(cgins.vitro) %>% table()

#### cgins.vitro marker identifications -----------------------------------------------------
## cluster markers
cgins.vitro.markers <- FindAllMarkers(cgins.vitro, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3)
top.cgins.vitro.markers <- cgins.vitro.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
## heatmap for markers 
DefaultAssay(cgins.vitro) <- "RNA"
DoHeatmap(cgins.vitro, top.cgins.vitro.markers$gene, angle = 90) 

#### cgins.vitro annotation ---------------------------------------------------------
## select resolution
Idents(cgins.vitro) <- "RNA_snn_res.1.2"
## annotation
cgins.vitro <- 
  cgins.vitro %>%
  RenameIdents("0" = "beta_like",
               "1" = "beta_like",
               "2" = "beta_like",
               "3" = "beta_like",
               "4" = "beta_like",
               "5" = "epsilon_like",
               "6" = "beta_like",
               "7" = "beta_like",
               "8" = "beta_like",
               "9" = "alpha_like",
               "10" = "beta_like",
               "11" = "beta_like",
               "12" = "epsilon_like",
               "13" = "delta_like")

## store the annotation
cgins.vitro$cell.type <- Idents(cgins.vitro)
Idents(cgins.vitro) <- factor(cgins.vitro$cell.type, levels = c("beta_like", "alpha_like", "delta_like", "epsilon_like"))

#### cgins.vitro marker re-do after annotation -----------------------------------------------------
## cluster markers
Idents(cgins.vitro) <- "cell.type"
cgins.vitro.markers <- FindAllMarkers(cgins.vitro, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3)
top.cgins.vitro.markers <- cgins.vitro.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
## visual for markers 
DefaultAssay(cgins.vitro) <- "RNA"
DoHeatmap(cgins.vitro, top.cgins.vitro.markers$gene, angle = 90) 
VlnPlot(cgins.vitro, top.cgins.vitro.markers$gene, stack = T, flip = T) + NoLegend()

#### rds -----------------------------------------------------
# cgins.vitro rds
saveRDS(cgins.vitro, file = "rds_manuscript/cgins.vitro.rds")
cgins.vitro <- readRDS(file = "rds_manuscript/cgins.vitro.rds")



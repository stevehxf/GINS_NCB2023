# This script processes the vivo cgins data into normalized seurat object

#### loading data -------
# read in the 10X matrices from your folder locations (you'll need to change)
cgins.vivo.data <- Read10X_h5(filename = 'data/VIVO_CGINSC_S9_GRCh38_mm10/filtered_feature_bc_matrix.h5')

# create the seurat objects
cgins.vivo <- CreateSeuratObject(cgins.vivo.data, min.features=300, min.cells=10, project="C-GINSC-VIVO")
cgins.vivo.data <- NULL

#### QC and ---------------------
# remove murine cell contaminants
## extract the count
counts <- GetAssayData(cgins.vivo, assay = "RNA")

## extract human and mouse gene counts
hg.idx <- rownames(counts) %>% str_detect("GRCh38-")
hg.counts <- counts[hg.idx,]
hg.counts.sum <- hg.counts %>% colSums()
mm.counts <- counts[!hg.idx,]
mm.counts.sum <- mm.counts %>% colSums()
## make a dataframe storing the feature counts from human or human for each cell
cell.counts.sum <- data.frame(hg=hg.counts.sum,mm=mm.counts.sum)
## extract the cell name whose human gene counts more than mouse gene counts
idx <- cell.counts.sum$hg > cell.counts.sum$mm
keep.cell <- rownames(cell.counts.sum)[idx]

# remove murine cells
cgins.vivo <- subset(cgins.vivo, cells = keep.cell)

# remove mouse genes except Ngn3, Pdx1 and Mafa
hg.gene <- rownames(cgins.vivo) %>% str_subset("GRCh38-")
cgins.vivo <- subset(cgins.vivo, features = c(hg.gene, "mm10---Pdx1", "mm10---Mafa", "mm10---Neurog3"))

# rename the genes
re.name <- row.names(cgins.vivo) %>% str_remove("GRCh38-")
tail(re.name)
length(re.name)
re.name[c(17813:17815)] <- c("PDX1MUS","MAFAMUS","NGN3ERMUS")
cgins.vivo@assays$RNA@counts@Dimnames[[1]] <- re.name
cgins.vivo@assays$RNA@data@Dimnames[[1]] <- re.name
rownames(cgins.vivo@assays$RNA@meta.features) <- re.name

# QC #
## create mito content
cgins.vivo[["percent.mt"]] <- PercentageFeatureSet(cgins.vivo, pattern = "^MT-")

## QC scatterplots
### cgins.vivo
p1 <- FeatureScatter(cgins.vivo, feature1 = "nCount_RNA", feature2 = "percent.mt")
p2 <- FeatureScatter(cgins.vivo, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p1+p2

### QC plots and make a decision as to cutoffs/filters
VlnPlot(cgins.vivo, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3,pt.size = 0.01, group.by = "orig.ident")

### quantile
quantile(cgins.vivo$nFeature_RNA, 0.03)
quantile(cgins.vivo$nFeature_RNA, 0.97)
quantile(cgins.vivo$nCount_RNA, 0.03)
quantile(cgins.vivo$nCount_RNA, 0.97)

## filter seurat objects - at this stage rerun the previous plts to see how they look now
## low feature cutoff 500 accroding to mapping qc.
VlnPlot(cgins.vivo, features = "percent.mt") + ylim(0,20)
VlnPlot(cgins.vivo, features = "nCount_RNA") + ylim(0,10000)
cgins.vivo <- subset(cgins.vivo, subset = nFeature_RNA > 1000 & nFeature_RNA < 5400 & percent.mt < 15 & nCount_RNA < 28000 & nCount_RNA > 3000)
VlnPlot(cgins.vivo, features = c("nCount_RNA","nFeature_RNA","percent.mt"))


# pre-process data using standard workflow for QC
cgins.vivo <- 
  cgins.vivo %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA()
DimHeatmap(cgins.vivo, dims = 1:21, cells = 500, balanced = TRUE)
ElbowPlot(cgins.vivo, ndims = 30)
dimensions <- 1:11
cgins.vivo <- 
  cgins.vivo %>% 
  RunTSNE(dims = dimensions) %>% 
  FindNeighbors(dims = dimensions) %>% 
  FindClusters(resolution = 0.5)
DimPlot(cgins.vivo, reduction = "tsne", label = T, cols = my36colors) + NoLegend()
VlnPlot(cgins.vivo, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3,pt.size = 0.01)

# find markers for QC 
## cluster markers
cgins.vivo.markers <- FindAllMarkers(cgins.vivo, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3)
top.cgins.vivo.markers <- cgins.vivo.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

## heatmap for markers 
DefaultAssay(cgins.vivo) <- "RNA"
DoHeatmap(cgins.vivo, top.cgins.vivo.markers$gene, angle = 90) 

## is cluster-1 doublet? it has all the markers from different clusters and high UMI. run doublet finder

#### cgins.vivo doublets finder ---------------------------------------------------------------------------------------
# pK Identification (no ground-truth) 
sweep.res.list <- paramSweep_v3(cgins.vivo, PCs = dimensions, sct = F)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

# Homotypic Doublet Proportion Estimate 
DoubletRate <- 0.05                
homotypic.prop <- modelHomotypic(cgins.vivo$RNA_snn_res.0.5) 
nExp_poi <- round(DoubletRate*ncol(cgins.vivo)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder with varying classification stringencies
cgins.vivo <- doubletFinder_v3(cgins.vivo, PCs = dimensions, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = F, sct = F)

# visuzalization
p1 <- DimPlot(cgins.vivo, reduction = "tsne", group.by = colnames(cgins.vivo@meta.data)[8])
p2 <- DimPlot(cgins.vivo, reduction = "tsne", group.by = "RNA_snn_res.0.5")
p1+p2
VlnPlot(cgins.vivo, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3,pt.size = 0.01, group.by = "DF.classifications_0.25_0.3_150" )

# remove the doublets
Idents(cgins.vivo) <- colnames(cgins.vivo@meta.data)[8]
cgins.vivo <- subset(cgins.vivo, ident = "Singlet")


#### cgins.vivo re-processing after doublet removal -----
DefaultAssay(cgins.vivo) <- "RNA"
cgins.vivo <- cgins.vivo %>% 
  # log-normalization 
  NormalizeData(normalization.method = "LogNormalize") %>% 
  # variable features
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  # scaling
  ScaleData(vars.to.regress = "percent.mt", features = c(cgins.vivo@assays$RNA@var.features, "INS")) %>% 
  # PCA
  RunPCA()

DimHeatmap(cgins.vivo, dims = 1:30, cells = 500, balanced = TRUE)
ElbowPlot(cgins.vivo, ndims = 30)
# umap and tsne
dimensions <- 1:8
cgins.vivo <- 
  cgins.vivo %>% 
#  RunUMAP(dims = dimensions) %>% 
  RunTSNE(dims = dimensions, seed.use = 5) 

#### cgins.vivo clustering analysis --------
cgins.vivo <- cgins.vivo %>%   
  FindNeighbors(dims = dimensions)

for (resolution in c(0.2, 0.3, 0.5, 0.7, 1.0, 1.3, 1.7, 2.0)) {
  cgins.vivo <- cgins.vivo %>%  
    FindClusters(resolution = resolution)}

#p1 <- DimPlot(cgins.vivo, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.1") + coord_equal() + NoLegend()
p2 <- DimPlot(cgins.vivo, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.2") + coord_equal() + NoLegend()
p3 <- DimPlot(cgins.vivo, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.3") + coord_equal() + NoLegend()
#p4 <- DimPlot(cgins.vivo, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.4") + coord_equal() + NoLegend()
p5 <- DimPlot(cgins.vivo, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.5") + coord_equal() + NoLegend()
p6 <- DimPlot(cgins.vivo, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.7") + coord_equal() + NoLegend()
p7 <- DimPlot(cgins.vivo, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.1") + coord_equal() + NoLegend()
p8 <- DimPlot(cgins.vivo, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.1.3") + coord_equal() + NoLegend()
p9 <- DimPlot(cgins.vivo, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.1.7") + coord_equal() + NoLegend()
p10 <- DimPlot(cgins.vivo, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.2") + coord_equal() + NoLegend()
plot_grid(p2,p3,p5,p6,p7,p8,p9,p10,ncol = 4)

DimPlot(cgins.vivo, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.1.7") + coord_equal()
FeaturePlot(cgins.vivo, c("SST","HHEX","TPH1","INS","ENTPD3"), reduction = "tsne")
# select resolution
sel.cluster <- "RNA_snn_res.1.3"
cgins.vivo <- cgins.vivo %>% SetIdent(value = sel.cluster)
DimPlot(cgins.vivo, reduction = "tsne", label = T, cols = my36colors) + coord_equal()
Idents(cgins.vivo) %>% table()

#### marker identifications ------
# marker identifications
## cluster markers
cgins.vivo.markers <- FindAllMarkers(cgins.vivo, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3)
top.cgins.vivo.markers <- cgins.vivo.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

#### Seurat rds export --------------------------------------------------------------
# export rds
saveRDS(cgins.vivo, file = "rds_manuscript/cgins.vivo.manuscript.rds")
cgins.vivo <- readRDS("rds_manuscript/cgins.vivo.manuscript.rds")

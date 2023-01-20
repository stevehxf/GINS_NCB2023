# This script processes the vivo cgins data into normalized seurat cgins.vitro.time
# run libraries_setup_and_colors.R first

## loading data #### 
# read in the 10X matrices from your folder locations (you'll need to change)
subset.gene <- function(x){
  x <- x$`Gene Expression`
}
cgins.vitro.time1.data <- Read10X_h5(filename = 'data/CGINS_TimeCourse/Zhou_XH_11617_2021_12_22/Sample1/count/sample_feature_bc_matrix.h5') %>% subset.gene()
cgins.vitro.time2.data <- Read10X_h5(filename = 'data/CGINS_TimeCourse/Zhou_XH_11617_2021_12_22/Sample2/count/sample_feature_bc_matrix.h5') %>% subset.gene()
cgins.vitro.time3.data <- Read10X_h5(filename = 'data/CGINS_TimeCourse/Zhou_XH_11617_2021_12_22/Sample3/count/sample_feature_bc_matrix.h5') %>% subset.gene()
cgins.vitro.time4.data <- Read10X_h5(filename = 'data/CGINS_TimeCourse/Zhou_XH_11617_2021_12_22/Sample4/count/sample_feature_bc_matrix.h5') %>% subset.gene()

# create the seurat cgins.vitro.times
cgins.vitro.time1 <- CreateSeuratObject(cgins.vitro.time1.data, min.features=300, min.cells=10, project="GSC")
cgins.vitro.time2 <- CreateSeuratObject(cgins.vitro.time2.data, min.features=300, min.cells=10, project="Endocrine")
cgins.vitro.time3 <- CreateSeuratObject(cgins.vitro.time3.data, min.features=300, min.cells=10, project="Early_GINS")
cgins.vitro.time4 <- CreateSeuratObject(cgins.vitro.time4.data, min.features=300, min.cells=10, project="Late_GINS")

# merge all the cgins.vitro.times
cgins.vitro.time <- merge(cgins.vitro.time1, y = c(cgins.vitro.time2, cgins.vitro.time3, cgins.vitro.time4), project = "GINS.TimeCourse")

## remove murine contaminates  ---------------------

## extract the count
counts <- GetAssayData(cgins.vitro.time, assay = "RNA")

# extract human and mouse gene counts
hg.idx <- rownames(counts) %>% str_detect("GRCh38-")
hg.counts <- counts[hg.idx,]
hg.counts.sum <- hg.counts %>% colSums()
mm.counts <- counts[!hg.idx,]
mm.counts.sum <- mm.counts %>% colSums()
# make a dataframe storing the feature counts from human or human for each cell
cell.counts.sum <- data.frame(hg=hg.counts.sum,mm=mm.counts.sum)
# extract the cell name whose human gene counts more than mouse gene counts
idx <- cell.counts.sum$hg > cell.counts.sum$mm
keep.cell <- rownames(cell.counts.sum)[idx]

# remove murine cells
cgins.vitro.time <- subset(cgins.vitro.time, cells = keep.cell)

# remove mouse genes except Ngn3, Pdx1 and Mafa
hg.gene <- rownames(cgins.vitro.time) %>% str_subset("GRCh38-")
cgins.vitro.time <- subset(cgins.vitro.time, features = c(hg.gene, "mm10---Pdx1", "mm10---Mafa", "mm10---Neurog3"))

# rename the genes
re.name <- row.names(cgins.vitro.time) %>% str_remove("GRCh38-")
tail(re.name)
length(re.name)
re.name[c((length(re.name)-2):length(re.name))] <- c("PDX1MUS","MAFAMUS","NGN3ERMUS")
cgins.vitro.time@assays$RNA@counts@Dimnames[[1]] <- re.name
cgins.vitro.time@assays$RNA@data@Dimnames[[1]] <- re.name
rownames(cgins.vitro.time@assays$RNA@meta.features) <- re.name


## QC ----------------
## create mito content
cgins.vitro.time[["percent.mt"]] <- PercentageFeatureSet(cgins.vitro.time, pattern = "^MT-")

## QC scatterplots
p1 <- FeatureScatter(cgins.vitro.time, feature1 = "nCount_RNA", feature2 = "percent.mt")
p2 <- FeatureScatter(cgins.vitro.time, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p1+p2

## QC plots and make a decision as to cutoffs/filters
cgins.vitro.time$orig.ident <-  factor(cgins.vitro.time$orig.ident, levels = c("GSC", "Endocrine", "Early_GINS", "Late_GINS"))
VlnPlot(cgins.vitro.time, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3,pt.size = 0.01, group.by = "orig.ident")

## split the object into a list
cgins.timeCourse.list <- SplitObject(cgins.vitro.time)

## quantile
nFeature_RNA.quantile <- function(object, n){quantile(object$nFeature_RNA, n)}
nCount_RNA.quantile <- function(object, n){quantile(object$nCount_RNA, n)}
mt.quantile <- function(object, n){quantile(object$percent.mt, n)}
sapply(cgins.timeCourse.list, nFeature_RNA.quantile, 0.05)
sapply(cgins.timeCourse.list, nFeature_RNA.quantile, 0.95)
sapply(cgins.timeCourse.list, nCount_RNA.quantile, 0.05)
sapply(cgins.timeCourse.list, nCount_RNA.quantile, 0.95)
sapply(cgins.timeCourse.list, mt.quantile, 0.05)
sapply(cgins.timeCourse.list, mt.quantile, 0.8)

## subset out low quality cells
cgins.timeCourse.list$GSC <- subset(cgins.timeCourse.list$GSC, subset = nFeature_RNA > 1300 & nFeature_RNA < 6100 & nCount_RNA < 25000 & nCount_RNA > 7000 & percent.mt < 15)
cgins.timeCourse.list$Endocrine <- subset(cgins.timeCourse.list$Endocrine, subset = nFeature_RNA > 640 & nFeature_RNA < 3400 & nCount_RNA < 10000 & nCount_RNA > 1600 & percent.mt < 20)
cgins.timeCourse.list$Early_GINS <- subset(cgins.timeCourse.list$Early_GINS, subset = nFeature_RNA > 400 & nFeature_RNA < 2800 & nCount_RNA < 7000 & nCount_RNA > 1100 & percent.mt < 20)
cgins.timeCourse.list$Late_GINS <- subset(cgins.timeCourse.list$Late_GINS, subset = nFeature_RNA > 360 & nFeature_RNA < 3200 & nCount_RNA < 10000 & nCount_RNA > 1200 & percent.mt < 20)

## merge
cgins.vitro.time <- merge(cgins.timeCourse.list$GSC, y = c(cgins.timeCourse.list$Endocrine, cgins.timeCourse.list$Early_GINS, cgins.timeCourse.list$Late_GINS), project = "GINS.TimeCourse")
VlnPlot(cgins.vitro.time, features = c("nCount_RNA","nFeature_RNA","percent.mt"))
table(cgins.vitro.time$orig.ident)

## pre-process data using standard workflow for QC
cgins.vitro.time <- 
  cgins.vitro.time %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA()
DimHeatmap(cgins.vitro.time, dims = 1:30, cells = 500, balanced = TRUE)
ElbowPlot(cgins.vitro.time, ndims = 50)
dimensions <- 1:15
cgins.vitro.time <- 
  cgins.vitro.time %>% 
  RunTSNE(dims = dimensions) %>% 
  FindNeighbors(dims = dimensions) %>% 
  FindClusters(resolution = 1)
p1 <- DimPlot(cgins.vitro.time, reduction = "tsne", label = T, cols = my36colors) + NoLegend()
p2 <- DimPlot(cgins.vitro.time, reduction = "tsne", label = T, cols = my36colors, group.by = "orig.ident") + NoLegend()
p1+p2
Idents(cgins.vitro.time) <- "RNA_snn_res.0.5"
VlnPlot(cgins.vitro.time, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3,pt.size = 0.01)
VlnPlot(cgins.vitro.time, c("INS","SST","GHRL","GCG"))

## cluster markers
cgins.vitro.time.markers <- FindAllMarkers(cgins.vitro.time, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5)
top.cgins.vitro.time.markers <- cgins.vitro.time.markers %>% group_by(cluster) %>% top_n(n = 8, wt = avg_log2FC)

## heatmap for markers 
DoHeatmap(cgins.vitro.time, top.cgins.vitro.time.markers$gene, angle = 90) 

## cgins.vitro.time doublets finder ---------------------------------------------------------------------------------------
# pK Identification (no ground-truth) 
sweep.res.list <- paramSweep_v3(cgins.vitro.time, PCs = dimensions, sct = F)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

# Homotypic Doublet Proportion Estimate 
DoubletRate <- 0.08                
homotypic.prop <- modelHomotypic(cgins.vitro.time$RNA_snn_res.0.5) 
nExp_poi <- round(DoubletRate*ncol(cgins.vitro.time)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder with varying classification stringencies
cgins.vitro.time <- doubletFinder_v3(cgins.vitro.time, PCs = dimensions, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = F, sct = F)

# visualization for putative doublets
cgins.vitro.time@meta.data %>% colnames()
p1 <- DimPlot(cgins.vitro.time, reduction = "tsne", group.by = "DF.classifications_0.25_0.18_850")
p2 <- DimPlot(cgins.vitro.time, reduction = "tsne", group.by = "RNA_snn_res.0.5")
p1+p2
VlnPlot(cgins.vitro.time, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3,pt.size = 0.01, group.by = "DF.classifications_0.25_0.18_850" )

# remove the doublets
Idents(cgins.vitro.time) <- "DF.classifications_0.25_0.18_850"
cgins.vitro.time <- subset(cgins.vitro.time, ident = "Singlet")

## cgins.vitro.time re-processing after doublet removal ---------------------------------------
cgins.vitro.time$orig.ident <-  factor(cgins.vitro.time$orig.ident, levels = c("GSC", "Endocrine", "Early_GINS", "Late_GINS"))
# cell cylce #
DefaultAssay(cgins.vitro.time) <- "RNA"
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
cgins.vitro.time <- CellCycleScoring(cgins.vitro.time, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
cgins.vitro.time <- cgins.vitro.time %>% 
  # log-normalization 
  NormalizeData(normalization.method = "LogNormalize") %>% 
  # variable features
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  # scaling
  #ScaleData(vars.to.regress = c("percent.mt")) %>% 
  ScaleData(vars.to.regress = c("percent.mt", "S.Score", "G2M.Score") ) %>% 
  # PCA
  RunPCA()

DimHeatmap(cgins.vitro.time, dims = 1:30, cells = 400, balanced = TRUE, ncol = 3)
ElbowPlot(cgins.vitro.time, ndims = 50)
# umap and tsne
dimensions <- 1:15
cgins.vitro.time <- 
  cgins.vitro.time %>% 
  RunUMAP(dims = dimensions) %>% 
  RunTSNE(dims = dimensions, seed.use = 5) 

p1 <- DimPlot(cgins.vitro.time, group.by = "Phase", reduction = "umap") + coord_fixed()
p2 <- DimPlot(cgins.vitro.time, group.by = "orig.ident", reduction = "umap") + coord_fixed()
p1+p2
p1 <- DimPlot(cgins.vitro.time, group.by = "Phase", reduction = "tsne") + coord_fixed()
p2 <- DimPlot(cgins.vitro.time, group.by = "orig.ident", reduction = "tsne") + coord_fixed()
p1+p2

## cgins.vitro.time clustering analysis ---------------------------------------------------------
cgins.vitro.time <- cgins.vitro.time %>%   
  FindNeighbors(dims = dimensions)
for (resolution in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5)) {
  cgins.vitro.time <- cgins.vitro.time %>%  
    FindClusters(resolution = resolution)}
p1 <- DimPlot(cgins.vitro.time, reduction = "umap", label = T, cols = my36colors, group.by = "RNA_snn_res.0.1") + coord_equal() + NoLegend()
p2 <- DimPlot(cgins.vitro.time, reduction = "umap", label = T, cols = my36colors, group.by = "RNA_snn_res.0.2") + coord_equal() + NoLegend()
p3 <- DimPlot(cgins.vitro.time, reduction = "umap", label = T, cols = my36colors, group.by = "RNA_snn_res.0.3") + coord_equal() + NoLegend()
p4 <- DimPlot(cgins.vitro.time, reduction = "umap", label = T, cols = my36colors, group.by = "RNA_snn_res.0.4") + coord_equal() + NoLegend()
p5 <- DimPlot(cgins.vitro.time, reduction = "umap", label = T, cols = my36colors, group.by = "RNA_snn_res.0.5") + coord_equal() + NoLegend()
p6 <- DimPlot(cgins.vitro.time, reduction = "umap", label = T, cols = my36colors, group.by = "RNA_snn_res.0.7") + coord_equal() + NoLegend()
p7 <- DimPlot(cgins.vitro.time, reduction = "umap", label = T, cols = my36colors, group.by = "RNA_snn_res.1") + coord_equal() + NoLegend()
p8 <- DimPlot(cgins.vitro.time, reduction = "umap", label = T, cols = my36colors, group.by = "RNA_snn_res.1.5") + coord_equal() + NoLegend()
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,nrow = 2)
DimPlot(cgins.vitro.time, reduction = "umap", label = T, cols = my36colors, group.by = "RNA_snn_res.0.5", split.by = "orig.ident") + coord_equal() + NoLegend()
gene.of.interest <- c("INS","GCG","SST","GHRL",
                      "ABCC9","ENTPD3", "UCN3",
                      "PCSK1","PCSK2","G6PC2","SLC30A8",
                      "MAFB","MAFA","HOPX","IAPP")
FeaturePlot(cgins.vitro.time, c("SST","HHEX","SOX4","INS","GCG","G6PC2","SOX9"), reduction = "umap")
FeaturePlot(cgins.vitro.time, gene.of.interest, reduction = "umap")
FeaturePlot(cgins.vitro.time, c("INS"), reduction = "umap", min.cutoff = 3, max.cutoff = 5)
VlnPlot(cgins.vitro.time,c("INS"), group.by = "orig.ident")

# select resolution
sel.cluster <- "RNA_snn_res.0.5"
cgins.vitro.time <- cgins.vitro.time %>% SetIdent(value = sel.cluster)
cgins.vitro.time <- BuildClusterTree(cgins.vitro.time, reorder.numeric = T, reorder = T, dims = dimensions)
PlotClusterTree(cgins.vitro.time)
cgins.vitro.time[['tree.ident']] <- Idents(cgins.vitro.time)
p1 <- DimPlot(cgins.vitro.time, reduction = "umap", label = T, cols = my36colors, group.by = "RNA_snn_res.0.5") + coord_equal()
p2 <- DimPlot(cgins.vitro.time, reduction = "umap", label = T, cols = my36colors, group.by = "tree.ident") + coord_equal()
p3 <- DimPlot(cgins.vitro.time, reduction = "umap", label = T, cols = my36colors, group.by = "Phase") + coord_equal()
p4 <- DimPlot(cgins.vitro.time, reduction = "umap", label = T, cols = my36colors, group.by = "orig.ident") + coord_equal()
p1+p2+p3+p4 & NoLegend()
DimPlot(cgins.vitro.time, reduction = "umap", label = T, cols = my36colors, group.by = "tree.ident", split.by = "orig.ident") + coord_equal()
# find markers
Idents(cgins.vitro.time) <- "tree.ident"
cgins.vitro.time.markers <- FindAllMarkers(cgins.vitro.time, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
top.cgins.vitro.time.markers <- cgins.vitro.time.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(cgins.vitro.time, top.cgins.vitro.time.markers$gene, angle = 90) + NoLegend()

## annotation ------------
Idents(cgins.vitro.time) <- "tree.ident"
cgins.vitro.time <- 
  cgins.vitro.time %>%
  RenameIdents("1" = "Mucoid",
               "4" = "Stem",
               
               "2" = "Endocrine_1",
               "3" = "Endocrine_2",
               
               "5" = "Beta_like",
               "6" = "Epsilon_like",
               "7" = "GINS_precursors",
               "8" = "Delta_like"
  )
## store the annotation
cgins.vitro.time$cell.type <- Idents(cgins.vitro.time)

## split to annotate seperately ---------------------
# split the objects into a list
cgins.timeCourse.list <- SplitObject(cgins.vitro.time, split.by = "orig.ident")
# normliazation/scaling/pca
cgins.timeCourse.list <- cgins.timeCourse.list %>% 
  lapply(FUN = function(x) {
    x <- x %>% 
      # log-normalization 
      NormalizeData(normalization.method = "LogNormalize") %>% 
      # variable features
      FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
      # scaling
      #ScaleData(vars.to.regress = c("percent.mt")) %>% 
      ScaleData(vars.to.regress = c("percent.mt", "S.Score", "G2M.Score")) %>% 
      # PCA
      RunPCA()
  })
cgins.timeCourse.list %>% lapply(DimHeatmap, dims = 1:24, cells = 500, balanced = TRUE, ncol =3)
cgins.timeCourse.list %>% lapply(ElbowPlot, ndims = 30)

# tsne
cgins.timeCourse.list$GSC <- cgins.timeCourse.list$GSC %>% RunTSNE(dims = 1:10, seed.use = 5)
DimPlot(cgins.timeCourse.list$GSC, reduction = "tsne")
DimPlot(cgins.timeCourse.list$GSC, reduction = "tsne", group.by = "Phase")
cgins.timeCourse.list$Endocrine <- cgins.timeCourse.list$Endocrine %>% RunTSNE(dims = 1:10, seed.use = 5)
DimPlot(cgins.timeCourse.list$GSC, reduction = "tsne")
cgins.timeCourse.list$Early_GINS <- cgins.timeCourse.list$Early_GINS %>% RunTSNE(dims = 1:9, seed.use = 5)
DimPlot(cgins.timeCourse.list$Early_GINS, reduction = "tsne")
cgins.timeCourse.list$Late_GINS <- cgins.timeCourse.list$Late_GINS %>% RunTSNE(dims = 1:8, seed.use = 5)
DimPlot(cgins.timeCourse.list$Late_GINS, reduction = "tsne")
# FindNeighbors
cgins.timeCourse.list <- cgins.timeCourse.list %>% 
  lapply(FindNeighbors, dims = dimensions)

# Clustering
for (resolution in c(0.1, 0.2, 0.3, 0.4, 0.5)) {
  cgins.timeCourse.list$GSC <- cgins.timeCourse.list$GSC %>%  
    FindClusters(resolution = resolution)}
for (resolution in c(0.1, 0.2, 0.3, 0.4, 0.5)) {
  cgins.timeCourse.list$Endocrine <- cgins.timeCourse.list$Endocrine %>%  
    FindClusters(resolution = resolution)}
for (resolution in c(0.1, 0.2, 0.3, 0.4, 0.5)) {
  cgins.timeCourse.list$Early_GINS <- cgins.timeCourse.list$Early_GINS %>%  
    FindClusters(resolution = resolution)}
for (resolution in c(0.1, 0.2, 0.3, 0.4, 0.5,0.7,1)) {
  cgins.timeCourse.list$Late_GINS <- cgins.timeCourse.list$Late_GINS %>%  
    FindClusters(resolution = resolution)}
saveRDS(cgins.timeCourse.list, file = "rds_manuscript/cgins.timeCourse.list.rds")
cgins.timeCourse.list <- readRDS("rds_manuscript/cgins.timeCourse.list.rds")

# GSC dimplot
p1 <- DimPlot(cgins.timeCourse.list$GSC, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.1") + coord_equal() + NoLegend()
p2 <- DimPlot(cgins.timeCourse.list$GSC, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.2") + coord_equal() + NoLegend()
p3 <- DimPlot(cgins.timeCourse.list$GSC, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.3") + coord_equal() + NoLegend()
p4 <- DimPlot(cgins.timeCourse.list$GSC, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.4") + coord_equal() + NoLegend()
p5 <- DimPlot(cgins.timeCourse.list$GSC, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.5") + coord_equal() + NoLegend()
plot_grid(p1,p2,p3,p4,p5, nrow = 1)
DimPlot(cgins.timeCourse.list$GSC, reduction = "tsne", group.by = "Phase")
# Endocrine dimplot
p1 <- DimPlot(cgins.timeCourse.list$Endocrine, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.1") + coord_equal() + NoLegend()
p2 <- DimPlot(cgins.timeCourse.list$Endocrine, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.2") + coord_equal() + NoLegend()
p3 <- DimPlot(cgins.timeCourse.list$Endocrine, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.3") + coord_equal() + NoLegend()
p4 <- DimPlot(cgins.timeCourse.list$Endocrine, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.4") + coord_equal() + NoLegend()
p5 <- DimPlot(cgins.timeCourse.list$Endocrine, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.5") + coord_equal() + NoLegend()
plot_grid(p1,p2,p3,p4,p5, nrow = 1)

# Early_GINS dimplot
p1 <- DimPlot(cgins.timeCourse.list$Early_GINS, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.1") + coord_equal() + NoLegend()
p2 <- DimPlot(cgins.timeCourse.list$Early_GINS, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.2") + coord_equal() + NoLegend()
p3 <- DimPlot(cgins.timeCourse.list$Early_GINS, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.3") + coord_equal() + NoLegend()
p4 <- DimPlot(cgins.timeCourse.list$Early_GINS, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.4") + coord_equal() + NoLegend()
p5 <- DimPlot(cgins.timeCourse.list$Early_GINS, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.5") + coord_equal() + NoLegend()
plot_grid(p1,p2,p3,p4,p5, nrow = 1)

# Late_GINS dimplot
p1 <- DimPlot(cgins.timeCourse.list$Late_GINS, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.1") + coord_equal() + NoLegend()
p2 <- DimPlot(cgins.timeCourse.list$Late_GINS, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.2") + coord_equal() + NoLegend()
p3 <- DimPlot(cgins.timeCourse.list$Late_GINS, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.3") + coord_equal() + NoLegend()
p4 <- DimPlot(cgins.timeCourse.list$Late_GINS, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.4") + coord_equal() + NoLegend()
p5 <- DimPlot(cgins.timeCourse.list$Late_GINS, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.5") + coord_equal() + NoLegend()
p6 <- DimPlot(cgins.timeCourse.list$Late_GINS, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.7") + coord_equal() + NoLegend()
p7 <- DimPlot(cgins.timeCourse.list$Late_GINS, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.1") + coord_equal() + NoLegend()
plot_grid(p1,p2,p3,p4,p5,p6,p7, nrow = 2)
p1 <- DimPlot(cgins.timeCourse.list$Late_GINS, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.5") + coord_equal() + NoLegend()
p2 <- DimPlot(cgins.timeCourse.list$Late_GINS, reduction = "tsne", label = T, cols = my36colors, group.by = "cell.type") + coord_equal() + NoLegend()
p1+p2
FeaturePlot(cgins.timeCourse.list$Late_GINS, c("GCG","GC"), reduction = "tsne")
# Findmarker
marker.list <- cgins.timeCourse.list
marker.list <- cgins.timeCourse.list %>% 
  lapply(FindAllMarkers, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3)
for (i in seq_along(marker.list)) {
  marker.list[[i]] %>% write_csv(file = paste0("markers/timeCourse.",names(marker.list)[i],".csv"))
  
}
top.marker.list <- marker.list %>% 
  lapply(FUN = function(x) {
    x <- x %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  })
DoHeatmap(cgins.timeCourse.list$GSC, top.marker.list$GSC$gene)
DoHeatmap(cgins.timeCourse.list$Endocrine, top.marker.list$Endocrine$gene)
DoHeatmap(cgins.timeCourse.list$Early_GINS, top.marker.list$Early_GINS$gene)
DoHeatmap(cgins.timeCourse.list$Late_GINS, top.marker.list$Late_GINS$gene)
## Re-annotation -------------------------------------------------------
Idents(cgins.timeCourse.list$Late_GINS) <- "RNA_snn_res.1"
table(cgins.timeCourse.list$Late_GINS$cell.type, cgins.timeCourse.list$Late_GINS$RNA_snn_res.1)
## annotate cluster 8 as Alpha_like
idx <- cgins.timeCourse.list$Late_GINS$RNA_snn_res.1 == "8"
cgins.timeCourse.list$Late_GINS$cell.type <- as.character(cgins.timeCourse.list$Late_GINS$cell.type)
cgins.timeCourse.list$Late_GINS$cell.type[idx] <- "Alpha_like"
cgins.timeCourse.list$Late_GINS$cell.type %>% table()
## reannoate the merged object accordingly
cell.idx <- names(cgins.timeCourse.list$Late_GINS$cell.type[cgins.timeCourse.list$Late_GINS$cell.type == "Alpha_like"])
cgins.vitro.time$cell.type <- as.character(cgins.vitro.time$cell.type)
cgins.vitro.time@meta.data[cell.idx,]$cell.type <- "Alpha_like"
cgins.vitro.time$cell.type %>% table()
## saveRDS -----------------------------------------------------
saveRDS(cgins.vitro.time, "rds_manuscript/cgins.vitro.time.rds")
cgins.vitro.time <- readRDS("rds_manuscript/cgins.vitro.time.rds")
## cell.type markers identification and visualization ----
# scale all genes
Idents(cgins.vitro.time) <- "cell.type"
deg <- cgins.vitro.time %>% FindAllMarkers(only.pos = TRUE, logfc.threshold = 0.3, min.pct = 0.3)
deg.sig <- deg %>% filter(p_val_adj < 0.01)
top.markers <- deg.sig %>% 
  group_by(cluster) %>% 
  top_n(15, wt = avg_log2FC)
var.genes <- cgins.vitro.time@assays$RNA@var.features
cgins.vitro.time <- cgins.vitro.time %>% 
  ScaleData(features = unique(c(var.genes, top.markers$gene)),
            vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
DoHeatmap(cgins.vitro.time, features = top.markers$gene)

## Seurat object visualization --------------------------------------------------
#cgins.vitro.time$cell.type <- Idents(cgins.vitro.time)
cgins.vitro.time$cell.type <- cgins.vitro.time$cell.type %>% 
  factor(levels = c("Stem","Mucoid","Endocrine_1","Endocrine_2","GINS_precursors","Delta_like","Epsilon_like","Alpha_like","Beta_like"))

# time course umap
p1 <- cgins.vitro.time %>% DimPlot(reduction = "umap", cols = my36colors, group.by = "cell.type") + coord_fixed() + ggtitle("Cell types") + theme_void() + NoLegend()
## cell cycle
p2 <- cgins.vitro.time %>% DimPlot(reduction = "umap", cols = my36colors, group.by = "Phase") + coord_fixed() + ggtitle("Cell cycles") + theme_void()
## the original multiplexed sample
p3 <- cgins.vitro.time %>% DimPlot(reduction = "umap", cols = my36colors, group.by = "orig.ident") + coord_fixed() + ggtitle("Original samples") + theme_void() + NoLegend()
## put together
plot_grid(p1, p2, p3, align = "vh", nrow = 1)

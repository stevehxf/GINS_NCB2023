# process islets scRNAseq from GEO datasets GSE114297

# loading data and removing ambient RNA by SoupX ----

## create a directory for the dataset
dir <- "data/islets_020422/"
if(!exists(dir)){
  dir.create(dir)
} 

## organize the data for reading into seurat (donor 1 2 3 and 9)
dir.list <- list.files(dir)
dir.list <- dir.list[c(1,5,6,12)]
dir.filt.list <- paste0(dir, dir.list, "/outs/filtered_feature_bc_matrix.h5") %>% as.list()
dir.raw.list <- paste0(dir, dir.list, "/outs/raw_feature_bc_matrix.h5") %>% as.list()

## create seurat object
islet.filter.data.list <- dir.filt.list %>% lapply(Read10X_h5)
islet.raw.data.list <- dir.raw.list %>% lapply(Read10X_h5)

names(islet.filter.data.list) <- dir.list %>% str_remove("_GRCh38")
names(islet.raw.data.list) <- dir.list %>% str_remove("_GRCh38")

## create soup object
soup.list <- vector(mode = "list", length = length(islet.filter.data.list))
names(soup.list) <- names(islet.filter.data.list)
for (i in seq_along(soup.list)) {
  soup.list[[i]] <-  SoupChannel(islet.raw.data.list[[i]],islet.filter.data.list[[i]])
}

seurat.list <- vector(mode = "list", length = length(islet.filter.data.list))
seurat.list <- islet.filter.data.list %>% lapply(CreateSeuratObject)
## clustering for soup
seurat.list <- 
  seurat.list %>% lapply(function(x){
    x %>% 
      SCTransform(verbose = F) %>% 
      RunPCA(verbose = F) %>% 
      RunUMAP(dims = 1:30, verbose = F) %>% 
      FindNeighbors(dims = 1:30, verbose = F) %>% 
      FindClusters(verbose = T)
  })

## add cluster info to the channel using setClusters.calculating ambient RNA profile.
for(i in seq_along(soup.list)) {
  meta    <- seurat.list[[i]]@meta.data
  umap    <- seurat.list[[i]]@reductions$umap@cell.embeddings
  soup.list[[i]]  <- setClusters(soup.list[[i]], setNames(meta$seurat_clusters, rownames(meta)))
  soup.list[[i]]  <- setDR(soup.list[[i]], umap)
  head(meta)
  # calculating ambient RNA profile.
  soup.list[[i]]  <- autoEstCont(soup.list[[i]])
}

head(soup.list[[1]]$soupProfile[order(soup.list[[1]]$soupProfile$est, decreasing = T), ], n = 30)
head(soup.list[[2]]$soupProfile[order(soup.list[[2]]$soupProfile$est, decreasing = T), ], n = 30)

## output integer matrix
adj.matrix.list <- soup.list
for(i in seq_along(soup.list)) {
  adj.matrix.list[[i]]  <- adjustCounts(soup.list[[i]], roundToInt = T)
  #DropletUtils:::write10xCounts(paste0(dir, dir.list,"_soupX_pbmc10k_filt", adj.matrix[[i]]))
}
adj.matrix.list %>% saveRDS("rds_manuscript/adj.matrix.list.rds")

## create seurat object from SoupX cleaned matrix

## create the seurat objects
islet.list.x <- adj.matrix.list %>% lapply(CreateSeuratObject, min.features=300, min.cells=10)
for (i in seq_along(islet.list.x)) {
  islet.list.x[[i]][["donor"]] <- names(islet.list.x)[i]
}
## merge excluding donor6, 10-12
islets.combined.x <- merge(x=islet.list.x$Donor_1, 
                           y=c(islet.list.x$Donor_2, 
                               islet.list.x$Donor_3,
                               islet.list.x$Donor_9), 
                           project = "islet")

# QC and  cell cycle scoring #### 
## create mito content
islets.combined.x[["percent.mt"]] <- PercentageFeatureSet(islets.combined.x, pattern = "^MT-")
## QC scatterplots
p1 <- FeatureScatter(islets.combined.x, feature1 = "nCount_RNA", feature2 = "percent.mt")
p2 <- FeatureScatter(islets.combined.x, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p1+p2
## QC plots and make a decision as to cutoffs/filters
VlnPlot(islets.combined.x, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3, pt.size = 0.01, group.by = "orig.ident")

## quantile
quantile(islets.combined.x$nFeature_RNA, 0.05)
quantile(islets.combined.x$nFeature_RNA, 0.95)
quantile(islets.combined.x$nCount_RNA, 0.05)
quantile(islets.combined.x$nCount_RNA, 0.95)

## filter seurat objects - at this stage rerun the previous plts to see how they look now
islets.combined.x <- subset(islets.combined.x, subset = nFeature_RNA > 450 &  nFeature_RNA < 3000 & percent.mt < 15 & nCount_RNA < 20000 & nCount_RNA > 1300)
# cell cylce #
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
islets.combined.x <- CellCycleScoring(islets.combined.x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
islet.list.x <- SplitObject(islets.combined.x, split.by = "donor")
#### Normalization and dimension deduction #### 
islet.list.x <- 
  islet.list.x %>% 
  lapply(NormalizeData) %>% 
  lapply(FindVariableFeatures) %>% 
  lapply(ScaleData) %>% 
  lapply(RunPCA)

dimensions <- 1:15
islet.list.x <- 
  islet.list.x %>% 
  lapply(RunUMAP, dim = dimensions) %>% 
  lapply(FindNeighbors, dim = dimensions) %>% 
  lapply(FindClusters, resolution = 0.1)
plot.list <- islet.list.x %>% lapply(DimPlot, label = T, cols = my36colors)
cowplot::plot_grid(plotlist = plot.list, ncol = 4)

# doublets finder ---------------------------------------------------------------------------------------
# pK Identification (no ground-truth) 
bcmvn <- islet.list.x %>% 
  lapply(paramSweep_v3, PCs = dimensions, sct = F) %>%
  lapply(summarizeSweep, GT = FALSE) %>% 
  lapply(find.pK)
pK_bcmvn <- bcmvn
for(i in seq_along(bcmvn)) {
  pK_bcmvn[[i]] <- bcmvn[[i]]$pK[which.max(bcmvn[[i]]$BCmetric)] %>% as.character() %>% as.numeric()
}

# Homotypic Doublet Proportion Estimate
DoubletRate <- 0.05
homotypic.prop <- islet.list.x
nExp_poi <- islet.list.x
nExp_poi.adj <- islet.list.x
for (i in seq_along(islet.list.x)) {
  homotypic.prop[[i]] <- modelHomotypic(islet.list.x[[i]]$RNA_snn_res.0.1) 
  nExp_poi[[i]] <- round(DoubletRate*ncol(islet.list.x[[i]])) 
  nExp_poi.adj[[i]] <- round(nExp_poi[[i]]*(1-homotypic.prop[[i]]))
  # Run DoubletFinder with varying classification stringencies
  islet.list.x[[i]] <- doubletFinder_v3(islet.list.x[[i]], 
                                        PCs = dimensions, 
                                        pN = 0.25, 
                                        pK = pK_bcmvn[[i]], 
                                        nExp = nExp_poi.adj[[i]], 
                                        reuse.pANN = F, 
                                        sct = F)
}

# visuzalization
plot.list <- islet.list.x %>% lapply(function(x){
  DF.group <- x@meta.data %>% names() %>% str_subset("DF.classifications")
  x %>% DimPlot(group.by = DF.group, label = T) + NoLegend()
})
cowplot::plot_grid(plotlist = plot.list)

# remove the doublets
islet.list.x <- islet.list.x %>% 
  lapply(function(x){
    DF.group <- x@meta.data %>% names() %>% str_subset("DF.classifications")
    Idents(x) <- DF.group
    x <- subset(x, ident = "Singlet")
  })

# Re-processing -----------
dimensions <- 1:15
resolution <- 0.1
islet.list.x <- 
  islet.list.x %>% 
  lapply(function(x){
    x <- x %>% 
      NormalizeData() %>% 
      FindVariableFeatures() %>% 
      ScaleData(vars.to.regress = "percent.mt") %>% 
      RunPCA() %>% 
      RunUMAP(dims = dimensions) %>% 
      FindNeighbors(dims = dimensions) %>% 
      FindClusters(resolution = resolution)
  })

# singleR annotation ####
library(scRNAseq)
sceM <- MuraroPancreasData()
# One should normally do cell-based quality control at this point, but for
# brevity's sake, we will just remove the unlabelled libraries here.
sceM <- sceM[,!is.na(sceM$label)]
# SingleR() expects reference datasets to be normalized and log-transformed.
library(scuttle)
sceM <- logNormCounts(sceM)
# annotate islets.combined.x with singleR
## data conversion and annotation
islets.sce.list <- islet.list.x %>% lapply(as.SingleCellExperiment)
rownames(sceM) <- rownames(sceM) %>% str_remove("__chr[1-9,X,Y]")
islets.sce.list <- islets.sce.list %>% lapply(SingleR, ref=sceM, labels=sceM$label, de.method="wilcox")
for (i in seq_along(islets.sce.list)){
  islet.list.x[[i]][["SingleR.labels"]] <- islets.sce.list[[i]]$labels
}
plotScoreHeatmap(islets.sce.list$Donor_1)

# visuzalization
plot.list <- islet.list.x %>% lapply(function(x){
  x %>% DimPlot(group.by = "SingleR.labels", cols = my36colors, label = T) + 
    theme(text = element_blank()) +
    NoLegend() +
    NULL
})
cowplot::plot_grid(plotlist = plot.list, ncol = 4)

## remove non-endocrine cells
for (i in seq_along(islet.list.x)) {
  Idents(islet.list.x[[i]]) <- "SingleR.labels"
  islet.list.x[[i]] <- subset(islet.list.x[[i]], idents = c("alpha","beta","delta","epsilon","pp"))
}


#### Perform integration -----------------------------------------------------
features <- SelectIntegrationFeatures(object.list = islet.list.x)
anchors <- FindIntegrationAnchors(object.list = islet.list.x,
                                  anchor.features = features,
                                  k.anchor = 5)

integrated_islet.x <- IntegrateData(anchorset = anchors)
DefaultAssay(integrated_islet.x) <- "integrated"

integrated_islet.x <- integrated_islet.x %>% 
  ScaleData(vars.to.regress = c("percent.mt", "S.Score", "G2M.Score")) %>% 
  RunPCA(verbose = FALSE)

dimensions <- 1:13
integrated_islet.x <- integrated_islet.x %>% 
  RunUMAP(dims = dimensions, seed.use = 5) 
DimPlot(integrated_islet.x, cols = my9colors)

## PCA heatmap
DimHeatmap(integrated_islet.x, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(integrated_islet.x, dims = 16:30, cells = 500, balanced = TRUE)
ElbowPlot(integrated_islet.x, ndims = 30)

## clustering analysis
integrated_islet.x <- integrated_islet.x %>%   
  FindNeighbors(dims = dimensions, k.param = 15)
for (resolution in c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1, 1.5)) {
  integrated_islet.x <- integrated_islet.x %>%  
    FindClusters(resolution = resolution)}

p1 <- DimPlot(integrated_islet.x, label = T, cols = my36colors, group.by = "integrated_snn_res.0.05") + coord_equal()
p2 <- DimPlot(integrated_islet.x, label = T, cols = my36colors, group.by = "integrated_snn_res.0.1") + coord_equal()
p3 <- DimPlot(integrated_islet.x, label = T, cols = my36colors, group.by = "integrated_snn_res.0.2") + coord_equal()
p4 <- DimPlot(integrated_islet.x, label = T, cols = my36colors, group.by = "integrated_snn_res.0.3") + coord_equal()
p5 <- DimPlot(integrated_islet.x, label = T, cols = my36colors, group.by = "integrated_snn_res.0.5") + coord_equal()
p6 <- DimPlot(integrated_islet.x, label = T, cols = my36colors, group.by = "integrated_snn_res.0.8") + coord_equal()
p7 <- DimPlot(integrated_islet.x, label = T, cols = my36colors, group.by = "integrated_snn_res.1") + coord_equal()
p8 <- DimPlot(integrated_islet.x, label = T, cols = my36colors, group.by = "integrated_snn_res.1.5") + coord_equal()
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8, ncol = 4)

p1 <- DimPlot(integrated_islet.x, label = F, cols = my36colors, group.by = "SingleR.labels") + coord_equal()
p2 <- DimPlot(integrated_islet.x, label = T, cols = my36colors, group.by = "integrated_snn_res.0.3") + coord_equal()
p1+p2

## QC plots and make a decision as to cutoffs/filters
VlnPlot(integrated_islet.x, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3, group.by = "integrated_snn_res.0.3",pt.size = 0.01)

## islet annotation -----
Idents(integrated_islet.x) <- "integrated_snn_res.0.3"
## Annotation
integrated_islet.x <- 
  integrated_islet.x %>%
  RenameIdents("0" = "alpha",
               "1" = "beta",
               "2" = "beta",
               "3" = "delta",
               "4" = "unclear",
               "5" = "pp",
               "6" = "unclear")
integrated_islet.x$cell.type <- Idents(integrated_islet.x) %>% as.character()
DefaultAssay(integrated_islet.x) <- 'RNA'
VlnPlot(integrated_islet.x, "GHRL")
sub <- subset(integrated_islet.x,subset = GHRL > 1)
epi.cell <- colnames(sub)
integrated_islet.x$cell.type[epi.cell] <- "epsilon"
DimPlot(integrated_islet.x, group.by = "cell.type")
FeaturePlot(integrated_islet.x, c("INS","GCG","SST","GHRL","PPY"), split.by = "donor")
VlnPlot(integrated_islet.x, c("INS","GCG","SST","GHRL","PPY"), group.by = "cell.type")
## remove unclear cells
integrated_islet.x_sub <- integrated_islet.x %>% subset(idents="unclear", invert=T)
dimensions <- 1:15
DefaultAssay(integrated_islet.x_sub) <- "integrated"
integrated_islet.x_sub <- integrated_islet.x_sub %>% 
  RunUMAP(dims = dimensions, seed.use = 5) 
## re-clustering
integrated_islet.x_sub <- integrated_islet.x_sub %>%   
  FindNeighbors(dims = dimensions, k.param = 15)
for (resolution in c( 0.1, 0.2, 0.3, 0.5)) {
  integrated_islet.x_sub <- integrated_islet.x_sub %>%  
    FindClusters(resolution = resolution)}
DimPlot(integrated_islet.x_sub, cols = my36colors)

DimPlot(integrated_islet.x_sub, group.by = "cell.type")
## find markers
DefaultAssay(integrated_islet.x_sub) <- 'RNA'
markers <- FindAllMarkers(integrated_islet.x_sub, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3)
FeaturePlot(integrated_islet.x_sub, c("GCG","INS","PPY","SST", "GHRL"))
DimPlot(integrated_islet.x_sub, cols = my9colors, group.by = "cell.type")
VlnPlot(integrated_islet.x_sub, c("GCG","INS","PPY","SST"), group.by = "cell.type")
# saveRDS ----
saveRDS(integrated_islet.x, "rds_manuscript/integrated_islet.x.rds")
saveRDS(integrated_islet.x_sub, "rds_manuscript/integrated_islet.x_sub.rds")

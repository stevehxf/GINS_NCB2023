#### the script to pre-process the agins in vitro data and integrate with agins in vitro and islets
#### author: Xiaofeng Huang
#### date: 11-23-2021
#### agins.vitro loading data --------------------------------------------------------------
# read in the 10X matrices from your folder locations (you'll need to change)
agins.data <- Read10X_h5(filename = 'data/AGINS_4/filtered_feature_bc_matrix.h5')

# create the seurat objects
agins.vitro <- CreateSeuratObject(agins.data, min.features=300, min.cells=10, project="agins.vitro")
agins.data <- NULL

#### agins.vitro QC ------------------------------------------------------------------------
# create mito content
agins.vitro[["percent.mt"]] <- PercentageFeatureSet(agins.vitro, pattern = "^MT-")

# QC scatterplots
p1 <- FeatureScatter(agins.vitro, feature1 = "nCount_RNA", feature2 = "percent.mt")
p2 <- FeatureScatter(agins.vitro, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p1+p2

# QC plots and make a decision as to cutoffs/filters
VlnPlot(agins.vitro, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3,pt.size = 0.01, group.by = "orig.ident")
VlnPlot(agins.vitro, features = c("percent.mt"), ncol = 3,pt.size = 0.01, group.by = "orig.ident") + ylim(c(0,15))
VlnPlot(agins.vitro, features = c("nCount_RNA"), ncol = 3,pt.size = 0.01, group.by = "orig.ident") + ylim(c(0,40000))

# quantile
quantile(agins.vitro$nFeature_RNA, 0.03)
quantile(agins.vitro$nFeature_RNA, 0.97)
quantile(agins.vitro$nCount_RNA, 0.03)
quantile(agins.vitro$nCount_RNA, 0.95)

# filter seurat objects - at this stage rerun the previous plts to see how they look now
agins.vitro <- subset(agins.vitro, subset = nFeature_RNA > 1000 &  nFeature_RNA < 6500 & percent.mt < 13 & nCount_RNA < 35000 & nCount_RNA > 2400)

# pre-process data using standard workflow for QC
agins.vitro <- 
  agins.vitro %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA()
DimHeatmap(agins.vitro, dims = 1:24, cells = 500, balanced = TRUE)
ElbowPlot(agins.vitro, ndims = 50)
dimensions <- 1:18
agins.vitro <- 
  agins.vitro %>% 
  RunUMAP(dims = dimensions) %>% 
  RunTSNE(dims = dimensions) %>% 
  FindNeighbors(dims = dimensions) %>% 
  FindClusters(resolution = 0.5)
DimPlot(agins.vitro, reduction = "tsne", label = T, cols = my36colors) + NoLegend()
VlnPlot(agins.vitro, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3,pt.size = 0.01)

# find markers for QC 
## cluster markers
agins.vitro.markers <- FindAllMarkers(agins.vitro, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3)
top.agins.vitro.markers <- agins.vitro.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

## heatmap for markers 
DefaultAssay(agins.vitro) <- "RNA"
DoHeatmap(agins.vitro, top.agins.vitro.markers$gene, angle = 90) 

## no clear issue

#### agins.vitro doublets finder ---------------------------------------------------------------------------------------
# pK Identification (no ground-truth) 
sweep.res.list <- paramSweep_v3(agins.vitro, PCs = dimensions, sct = F)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

# Homotypic Doublet Proportion Estimate 
DoubletRate <- 0.05               
homotypic.prop <- modelHomotypic(agins.vitro$RNA_snn_res.0.5) 
nExp_poi <- round(DoubletRate*ncol(agins.vitro)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder with varying classification stringencies
agins.vitro <- doubletFinder_v3(agins.vitro, PCs = dimensions, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = F, sct = F)

# visuzalization
agins.vitro@meta.data %>% colnames()
p1 <- DimPlot(agins.vitro, reduction = "tsne", group.by = colnames(agins.vitro@meta.data)[8])
p2 <- DimPlot(agins.vitro, reduction = "tsne", group.by = "RNA_snn_res.0.5")
p1+p2
VlnPlot(agins.vitro, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3,pt.size = 0.01, group.by = "DF.classifications_0.25_0.16_223" )

# remove the doublets
Idents(agins.vitro) <-  colnames(agins.vitro@meta.data)[8]
agins.vitro <- subset(agins.vitro, ident = "Singlet")

#### agins.vitro re-processing after doublet removal ---------------------------------------
DefaultAssay(agins.vitro) <- "RNA"
agins.vitro <- agins.vitro %>% 
  # log-normalization 
  NormalizeData(normalization.method = "LogNormalize") %>% 
  # variable features
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) 
var.gene <- VariableFeatures(agins.vitro)
agins.vitro <- agins.vitro %>% 
  # scaling
  ScaleData(vars.to.regress = "percent.mt", features = c(var.gene, "INS")) %>% 
  # PCA
  RunPCA()

ElbowPlot(agins.vitro, ndims = 30)
# umap and tsne
dimensions <- 1:15
agins.vitro <- 
  agins.vitro %>% 
  RunUMAP(dims = dimensions) %>% 
  RunTSNE(dims = dimensions) %>% 
  FindNeighbors(dims = dimensions)

#### agins.vitro clustering analysis ---------------------------------------------------------
agins.vitro <- agins.vitro %>%   
  FindNeighbors(dims = dimensions)

for (resolution in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.2)) {
  agins.vitro <- agins.vitro %>%  
    FindClusters(resolution = resolution)}

p1 <- DimPlot(agins.vitro, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.1") + coord_equal()
p2 <- DimPlot(agins.vitro, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.2") + coord_equal()
p3 <- DimPlot(agins.vitro, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.3") + coord_equal()
p4 <- DimPlot(agins.vitro, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.4") + coord_equal()
p5 <- DimPlot(agins.vitro, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.5") + coord_equal()
p6 <- DimPlot(agins.vitro, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.7") + coord_equal()
p7 <- DimPlot(agins.vitro, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.1") + coord_equal()
p8 <- DimPlot(agins.vitro, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.1.2") + coord_equal()
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,ncol = 4)

# select resolution
sel.cluster <- "RNA_snn_res.0.7"
agins.vitro <- agins.vitro %>% SetIdent(value = sel.cluster)
DimPlot(agins.vitro, reduction = "tsne", label = T, cols = my36colors) + coord_equal()
Idents(agins.vitro) %>% table()

#### agins.vitro marker identifications -----------------------------------------------------
## cluster markers
agins.vitro.markers <- FindAllMarkers(agins.vitro, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3)
top.agins.vitro.markers <- agins.vitro.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
## heatmap for markers 
DefaultAssay(agins.vitro) <- "RNA"
DoHeatmap(agins.vitro, top.agins.vitro.markers$gene, angle = 90) 

#### annotate agins and integrated agins according to the integration result ------------------------------
## annotation of agins
Idents(agins.vitro) <- "RNA_snn_res.0.7"
agins.vitro <- 
  agins.vitro %>%
  RenameIdents("0" = "beta_like",
               "1" = "beta_like",
               "2" = "beta_like",
               "3" = "beta_like",
               "4" = "beta_like",
               "5" = "beta_like",
               "6" = "beta_like",
               "7" = "G_like",
               "8" = "beta_like",
               "9" = "delta_like")
DimPlot(agins.vitro, reduction = "tsne", cols = my36colors)

VlnPlot(agins.vitro, c("INS","GCG","SST","GHRL","GAL","TAC3","GAST","GC"), stack = T, flip = T) + NoLegend()
FeaturePlot(agins.vitro, c("INS","GCG","SST","GHRL","GAL","TAC3","GAST","GC"), reduction = "tsne")
### store the annotation
agins.vitro$cell.type <- Idents(agins.vitro)
#### agins.vitro marker identifications after annotation -----------------------------------------------------
## cluster markers
agins.vitro.markers <- FindAllMarkers(agins.vitro, only.pos = TRUE, min.pct = 0.4, logfc.threshold = 0.4)
top.agins.vitro.markers <- agins.vitro.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
## heatmap for markers 
DefaultAssay(agins.vitro) <- "RNA"
theme_set(theme_classic())
DoHeatmap(agins.vitro, top.agins.vitro.markers$gene, angle = 90) 

#### rds  -----------------------------------------------------
saveRDS(agins.vitro, file = "rds_manuscript/agins.vitro.rds")
#### read in cgins.vitro and islet endocrine cells for integration analysis -----
agins.vitro <- readRDS(file = "rds_manuscript/agins.vitro.rds")
cgins.vitro <- readRDS(file = "rds_manuscript/cgins.vitro.rds")
agins.vitro$tissue <- "Antrum"
cgins.vitro$tissue <- "Corpus"
#### create a list containing cgins.vitro and islets and log-normalize -----------------------------------------------------
## combine list containing with cgins.vitro, cgins.vivo and islets from individual donor
list_cgins.vitro_agins.vitro <- list(cgins.organoid = cgins.vitro, agins.organoid = agins.vitro)

## log normalization
list_cgins.vitro_agins.vitro <- list_cgins.vitro_agins.vitro %>% 
  lapply(FUN = function(x) {
    DefaultAssay(x) <- "RNA"
    x <- x %>% 
      NormalizeData() %>% 
      FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
  })

#### cgins.vitro + agins.vitro Perform integration -----------------------------------------------------
features <- SelectIntegrationFeatures(object.list = list_cgins.vitro_agins.vitro)
anchors <- FindIntegrationAnchors(object.list = list_cgins.vitro_agins.vitro,
                                  anchor.features = features,
                                  k.anchor = 1)

integrated_cgins.agins.vitro <- IntegrateData(anchorset = anchors)

DefaultAssay(integrated_cgins.agins.vitro) <- "integrated"

dimensions <- 1:20
integrated_cgins.agins.vitro <- integrated_cgins.agins.vitro %>% 
  ScaleData(vars.to.regress = "percent.mt", verbose = FALSE) %>% 
  RunPCA(verbose = FALSE) %>% 
  #  RunUMAP(dims = dimensions) %>% 
  RunTSNE(dims = dimensions)

integrated_cgins.agins.vitro[["inte.idents"]] <- Idents(integrated_cgins.agins.vitro)
DimPlot(integrated_cgins.agins.vitro, reduction = "tsne", label = T, cols = my36colors, group.by = "inte.idents", split.by = "tissue")
DimPlot(integrated_cgins.agins.vitro, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.7",split.by = "tissue", ncol = 4)
DimPlot(integrated_cgins.agins.vitro, reduction = "tsne", cols = my36colors, group.by = "inte.idents")
DimPlot(integrated_cgins.agins.vitro, reduction = "tsne", label = T, cols = my36colors, group.by = "tissue")
DimPlot(integrated_cgins.agins.vitro, reduction = "tsne", label = T, cols = my36colors, group.by = "inte.idents", split.by = "tissue")

## PCA heatmap
DimHeatmap(integrated_cgins.agins.vitro, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(integrated_cgins.agins.vitro, dims = 16:30, cells = 500, balanced = TRUE)
ElbowPlot(integrated_cgins.agins.vitro, ndims = 50)


## clustering analysis
integrated_cgins.agins.vitro <- integrated_cgins.agins.vitro %>%   
  FindNeighbors(dims = dimensions)

for (resolution in c(0.3, 0.5, 0.7, 1.0)) {
  integrated_cgins.agins.vitro <- integrated_cgins.agins.vitro %>%  
    FindClusters(resolution = resolution)}

## cluster visualization
#clustree(integrated_cgins.agins.vitro@meta.data, prefix = "integrated_snn_res.")
#ggsave("figures/inte.clustree.pdf", width = 8, height = 9)
p1 <- DimPlot(integrated_cgins.agins.vitro, reduction = "tsne", label = T, cols = my36colors, group.by = "integrated_snn_res.0.3") + coord_equal()
p2 <- DimPlot(integrated_cgins.agins.vitro, reduction = "tsne", label = T, cols = my36colors, group.by = "integrated_snn_res.0.5") + coord_equal()
p3 <- DimPlot(integrated_cgins.agins.vitro, reduction = "tsne", label = T, cols = my36colors, group.by = "integrated_snn_res.0.7") + coord_equal()
p4 <- DimPlot(integrated_cgins.agins.vitro, reduction = "tsne", label = T, cols = my36colors, group.by = "integrated_snn_res.1") + coord_equal()
plot_grid(p1,p2,p3,p4)

#### visualization for publications ----------------------------------------------------------------------------
DefaultAssay(integrated_cgins.agins.vitro) <- "RNA"
Idents(integrated_cgins.agins.vitro) <- "cell.type"

integrated_cgins.agins.vitro <- integrated_cgins.agins.vitro %>% 
  NormalizeData() %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(vars.to.regress = "percent.mt", verbose = FALSE)
inte.split <- SplitObject(integrated_cgins.agins.vitro, split.by = "tissue")

# dimplot
set.theme <- theme_set(
  theme_classic()+
    theme(plot.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
    )
)

p1 <- DimPlot(integrated_cgins.agins.vitro, 
              reduction = "tsne", 
              cols = gins.color, 
              label = F) + 
  coord_fixed() +
  xlim(-40,40) +
  ylim(-40,40) +
  set.theme +
  NoLegend() +
  NULL

p2 <- DimPlot(inte.split$Corpus,  
              reduction = "tsne", 
              cols = gins.color, 
              label = F) + 
  coord_fixed() +
  xlim(-40,40) +
  ylim(-40,40) +
  set.theme +
  NoLegend() +
  NULL

p3 <- DimPlot(inte.split$Antrum,  
              reduction = "tsne", 
              cols = gins.color, 
              label = F) + 
  coord_fixed() +
  xlim(-40,40) +
  ylim(-40,40) +
  set.theme +
  NoLegend() +
  NULL
cowplot::plot_grid(p1,p2,p3, nrow = 1)
ggsave("figures_manuscript/ExtendedDataFig5d.tsne.jpeg", dpi = 400, width = 9, height = 3)

# paste tissue identity with cell identity
integrated_cgins.agins.vitro$tissue_cell <- paste(integrated_cgins.agins.vitro$tissue, integrated_cgins.agins.vitro$cell.type, sep = "_")
Idents(integrated_cgins.agins.vitro) <- "tissue_cell"
Idents(integrated_cgins.agins.vitro) %>% table()
Idents(integrated_cgins.agins.vitro) <- 
  factor(integrated_cgins.agins.vitro$tissue_cell, 
         levels = c("Corpus_beta_like", "Antrum_beta_like",
                    "Corpus_delta_like", "Antrum_delta_like", 
                    "Corpus_alpha_like", 
                    "Corpus_epsilon_like",
                    "Antrum_G_like"
         ))

## plot pie chart ----
### create a list for pie chart
inte.split <- SplitObject(integrated_cgins.agins.vitro, split.by = "tissue")
pie.tb.list <- inte.split
for (i in seq(inte.split)) {
  pie.tb.list[[i]] <- inte.split[[i]]@meta.data %>% 
    dplyr::select(cell.type) %>% 
    as_tibble() %>% 
    group_by(cell.type) %>% 
    summarise(n=n()) %>% 
    mutate(prop = round(n/sum(n)*100, digits = 4))
  pie.tb.list[[i]]$cell.type <- pie.tb.list[[i]]$cell.type %>% as.character()
  
}
pie.tb.list$Corpus$n %>% sum()
pie.tb.list$Antrum$n %>% sum()
### pie chart parameters
blank_theme <- theme_void()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )
### plot pie
p1 <- pie.tb.list$Corpus %>% 
  ggplot(aes(x="", y=n, fill=cell.type)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = gins.color)+
  blank_theme+
  theme(legend.position = "none") +
  NULL
p2 <- pie.tb.list$Antrum %>% 
  ggplot(aes(x="", y=n, fill=cell.type)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values= gins.color[c(1,3,5)])+
  blank_theme+
  theme(legend.position = "none") +
  NULL
ggsave(plot=p1,"figures_manuscript/ExtendedDataFig5d.corpus.pie.ident.pdf", height = 2, width = 2)
ggsave(plot=p2,"figures_manuscript/ExtendedDataFig5d.antral.pie.ident.pdf", height = 2, width = 2)


DefaultAssay(integrated_cgins.agins.vitro) <- "RNA"
integrated_cgins.agins.vitro$tissue_cell <- 
  factor(integrated_cgins.agins.vitro$tissue_cell, 
         levels = rev(c("Corpus_beta_like", "Antrum_beta_like",
                        "Corpus_delta_like", "Antrum_delta_like", 
                        "Corpus_alpha_like", 
                        "Corpus_epsilon_like",
                        "Antrum_G_like"
         )))
Idents(integrated_cgins.agins.vitro) <- "tissue_cell"



## subset beta like cells for beta cell genes comparison --------
Idents(integrated_cgins.agins.vitro) <- "cell.type"
cgins.agins.beta <- subset(integrated_cgins.agins.vitro, ident="beta_like")
DefaultAssay(cgins.agins.beta) <- "RNA"
Idents(cgins.agins.beta) <- "orig.ident"

# VlnPlot
identity <- c("INS", "GCK","PAX4","PAX6", "PDX1", "RFX6","MNX1",
              "NEUROD1", "NKX2-2","NKX6-1", "PCSK1",
              "GLP1R","GIPR","GJD2",
              "ABCC8","KCNJ11","SLC30A8","TCF7L2") 
metabolism <- c("PCSK1","CPE", "IAPP", "MIF","CHGA","CHGB",
                "SCG2","SCG3","SCG5","SLC30A8","PTPRN","PTPRN2","PAM",
                "CLCN3","KCNJ11","ABCC8", "RAB27A","RAP1GDS1", 
                "RAB1A", "RAB2A", "RAB3A", "NUCB2","PDIA3","YWHAZ",
                "ESRRG")
exocytosis <- c("INS","STX1A","STX1B", "VAMP2", "SNAP25","STX4","STXBP1","CDC42")
insSecretion <- c("CACNA1D","PCSK1N", "CHGB", 
                  "CAMK2N1","KCNJ3","PTPRN","HADH",
                  "CASR", "KCNK3","GRN","TTR")
mody <- c("ABCC8", "APPL1", "BLK", "CEL", "GCK", 
          "HNF1A", "HNF1B", "HNF4A", "INS", 
          "KCNJ11", "KLF11", "NEUROD1", "PAX4", "PDX1") 

Vln.list <- list(
  "metabolism" = metabolism,
  "identity" = identity,
  "mody" = mody,
  "ins_secreting" = insSecretion,
  "exocytosis" = exocytosis
)

DefaultAssay(cgins.agins.beta) <- 'RNA'
for (i in seq_along(Vln.list)) {
  VlnPlot(cgins.agins.beta, 
          Vln.list[[i]], 
          stack = T, flip = T, 
          cols = cVSa, 
          fill.by = "ident" )  +
    theme(axis.title = element_blank(),
          axis.text.x = element_blank(),
          strip.text = element_text(size = 0),
          strip.background = element_blank()
    ) +
    NoLegend()
  file_name <- paste0("figures_manuscript/ExtendedDataFig5e_cVSa_", names(Vln.list)[i], "_vln.pdf")
  ggsave(file_name, width = 1, height = length(Vln.list[[i]])/4.5)
}


#### rds: cgins.vitro and agins.vitro integration -----------------------------------------------------
saveRDS(integrated_cgins.agins.vitro, file = "rds_manuscript/integrated_cgins.agins.vitro.rds")
integrated_cgins.agins.vitro <- readRDS("rds_manuscript/integrated_cgins.agins.vitro.rds")

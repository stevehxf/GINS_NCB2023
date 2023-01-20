
#### loading seurat object rds -----------------------------------------------------
integrated_islet_sub <- readRDS(file = "rds_manuscript/integrated_islet.x_sub.rds")
cgins.vitro <- readRDS(file = "rds_manuscript/cgins.vitro.rds")

#### islets+cgins.vitro: create a list containing cgins.vitro and islets and log-normalize -----------------------------------------------------
integrated_islet_sub$orig.ident <- paste0(integrated_islet_sub$donor, "_islets")
integrated_islet_sub$source <- "islets"
cgins.vitro$source <- "gsins"
## combine list containing with cgins.vitro and islets from individual donor
islet.list <- SplitObject(integrated_islet_sub, split.by = "orig.ident")
list_cgins.vitro_islet <- c(list(cgins.vitro = cgins.vitro), islet.list)

## log normalization
list_cgins.vitro_islet <- list_cgins.vitro_islet %>% 
  lapply(FUN = function(x) {
    DefaultAssay(x) <- "RNA"
    x <- x %>% 
      NormalizeData() %>% 
      FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
  })



#### islets+cgins.vitro: Perform integration -----------------------------------------------------
features <- SelectIntegrationFeatures(object.list = list_cgins.vitro_islet)
anchors <- FindIntegrationAnchors(object.list = list_cgins.vitro_islet,
                                  anchor.features = features,
                                  k.anchor = 2)

integrated_cgins.vitro_islet <- IntegrateData(anchorset = anchors)

DefaultAssay(integrated_cgins.vitro_islet) <- "integrated"


integrated_cgins.vitro_islet <- integrated_cgins.vitro_islet %>% 
  ScaleData(verbose = FALSE) %>% 
  #ScaleData(vars.to.regress = "percent.mt", verbose = FALSE) %>% 
  RunPCA(verbose = FALSE, npcs = 30)

## PCA heatmap
#DimHeatmap(integrated_cgins.vitro_islet, dims = 1:15, cells = 500, balanced = TRUE)
#DimHeatmap(integrated_cgins.vitro_islet, dims = 16:30, cells = 500, balanced = TRUE)
#ElbowPlot(integrated_cgins.vitro_islet, ndims = 50)

#dimensions <- 1:7

dimensions <- 1:20
integrated_cgins.vitro_islet <- integrated_cgins.vitro_islet %>% 
  RunUMAP(dims = dimensions, seed.use = 1)

#integrated_cgins.vitro_islet <- integrated_cgins.vitro_islet %>% RunTSNE(dims = dimensions, seed.use = 5)
p1 <- DimPlot(integrated_cgins.vitro_islet, reduction = "umap", label = T, cols = my9colors, group.by = "cell.type", split.by = "source", ncol = 4) + NoLegend()
p2 <- DimPlot(integrated_cgins.vitro_islet, reduction = "umap", label = T, cols = my9colors, group.by = "cell.type") + NoLegend()
p1+p2

##  clustering analysis ---------------------------------------------------------
integrated_cgins.vitro_islet <- integrated_cgins.vitro_islet %>%   
  FindNeighbors(dims = dimensions)
for (resolution in c(0.05, 0.2)) {
  integrated_cgins.vitro_islet <- integrated_cgins.vitro_islet %>%  
    FindClusters(resolution = resolution)}
p1 <- DimPlot(integrated_cgins.vitro_islet, reduction = "umap", label = T, cols = my36colors, group.by = "integrated_snn_res.0.05", split.by = "source") + coord_equal() + NoLegend()
p2 <- DimPlot(integrated_cgins.vitro_islet, reduction = "umap", label = T, cols = my36colors, group.by = "integrated_snn_res.0.2", split.by = "source") + coord_equal() + NoLegend()
p1+p2

#### rds ----------
saveRDS(integrated_cgins.vitro_islet, "rds_manuscript/integrated_cgins.vitro_islet.rds")
#integrated_cgins.vitro_islet <- readRDS("rds_manuscript/integrated_cgins.vitro_islet.rds")

#### integration visualization for publication-------------------------------------------
# dimplot
Idents(integrated_cgins.vitro_islet) <- "cell.type"
integrated_cgins.vitro_islet@meta.data$cell.type %>% table()

inte.gins.order <- c("alpha", "alpha_like", "beta", "beta_like",  "delta", "delta_like", "epsilon", "epsilon_like","pp")
inte.islet.order <- c("alpha_like", "alpha","beta_like","beta",  "delta_like", "delta", "epsilon", "epsilon_like","pp")
p1 <- DimPlot(integrated_cgins.vitro_islet, 
              reduction = "umap", 
              group.by = "cell.type",
              cols = inte.color, 
              order = inte.islet.order,
              label = F) + 
  theme_void() +
  theme(plot.title = element_blank()) +
  coord_fixed(ratio = 0.5) +
  NoLegend() +
  NULL

p2 <- DimPlot(integrated_cgins.vitro_islet, 
              reduction = "umap", 
              group.by = "cell.type",
              cols = c(
                alpha = "#EBEBEB", # grey
                alpha_like = "#ffca3a", # yellow
                beta = "#EBEBEB", # grey 
                beta_like = "#ff595e", # red
                delta = "#EBEBEB", # grey 
                delta_like = "#8ac926", # green
                epsilon = "#EBEBEB", # grey
                epsilon_like = "#1982c4", # blue
                pp = "#EBEBEB" # grey
              ), 
              order = inte.islet.order,
              label = F) + 
  theme_void() +
  theme(plot.title = element_blank()) +
  coord_fixed(ratio = 0.5) +
  NoLegend() +
  NULL

p3 <- DimPlot(integrated_cgins.vitro_islet, 
              reduction = "umap", 
              group.by = "cell.type",
              cols = c(
                alpha = "#ffca3a", # yellow
                alpha_like = "#EBEBEB", # grey
                beta = "#ff595e", # red
                beta_like = "#EBEBEB", # grey
                delta = "#8ac926", # green
                delta_like = "#EBEBEB", # grey
                epsilon = "#1982c4", # blue
                epsilon_like = "#EBEBEB", # grey
                pp = "#6a4c93"
              ),
              order = inte.gins.order,
              label = F) + 
  theme_void() +
  theme(plot.title = element_blank()) +
  coord_fixed(ratio = 0.5) +
  NoLegend() +
  NULL

ggsave(plot = p1, "figures_manuscript/Fig3a.combined.jpeg", dpi = 400, width = 3, height = 3, bg = "white")
ggsave(plot = p2, "figures_manuscript/Fig3a.gins.jpeg", dpi = 400, width = 3, height = 3, bg = "white")
ggsave(plot = p3, "figures_manuscript/Fig3a.islet.jpeg", dpi = 400, width = 3, height = 3, bg = "white")


#### cgins marker visualization -----------------------------------------------------
## cluster markers
Idents(cgins.vitro) <- "cell.type"
cgins.vitro.markers <- FindAllMarkers(cgins.vitro, only.pos = TRUE)
top.cgins.vitro.markers <- 
  cgins.vitro.markers %>% 
  filter(p_val_adj < 0.01) %>% 
  group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
DoHeatmap(cgins.vitro, 
          features = c(top.cgins.vitro.markers$gene),
          #disp.min = -1.2, 
          #disp.max = 2
)# +  scale_fill_viridis(option = "B")


#### cgins.vitro visualization for publication-------------------------------------
# genes dotplot visualization
cgins.vitro$cell.type <- cgins.vitro$cell.type %>% factor(levels = rev(c("beta_like", "alpha_like", "delta_like", "epsilon_like")))
Idents(cgins.vitro) <- "cell.type"
feature.to.plot <- c("INS", "SCG2", "RGS16",
                     "GCG","TTR","GC","ARX",
                     "SST", "HHEX",
                     "GHRL"
)

## heatmap with dotplot
TileDotplot(cgins.vitro, 
            features = feature.to.plot, 
            scale.by = "size",
            col.max = 1,
            col.min = -1
) + 
  coord_flip() +
  coord_fixed()

ggsave("figures_manuscript/Fig3b.tiledotplot.pdf", width = 6, height = 6)

#### gene visualization for publication --------------------------------------------
# store orig.ident + cell.type as orig.cell.type
integrated_cgins.vitro_islet$orig.cell.type <- paste(integrated_cgins.vitro_islet$orig.ident, 
                                                     integrated_cgins.vitro_islet$cell.type,
                                                     sep="_")
Idents(integrated_cgins.vitro_islet) <- "cell.type"

# re-order the cell.type factor
order <- c("beta","beta_like",
           "alpha","alpha_like",
           "delta","delta_like",
           "epsilon","epsilon_like",
           "pp")
integrated_cgins.vitro_islet$cell.type <- factor(integrated_cgins.vitro_islet$cell.type, levels = order)

#### comparison between beta and beta_like in vitro ------------------------------------------
## 1. beta cell markers vs other islet endocrines
Idents(integrated_cgins.vitro_islet) <- "cell.type"
DefaultAssay(integrated_cgins.vitro_islet) <- "RNA"
DEG.beta.VSotherEndo <- FindMarkers(integrated_cgins.vitro_islet, ident.1 = "beta", ident.2 = c("alpha","delta","epsilon","pp"))
beta.markers <- DEG.beta.VSotherEndo %>% 
  filter(p_val_adj < 0.05 & avg_log2FC > 0) %>% 
  rownames_to_column("gene") %>% 
  as_tibble()

## 2. DEG analysis between beta and beta_like
Idents(integrated_cgins.vitro_islet) <- "cell.type"
DEG.betaVSbetaLike <- FindMarkers(integrated_cgins.vitro_islet, 
                                  ident.1 = "beta", ident.2 = "beta_like",
                                  min.pct = 0.3,
                                  max.cells.per.ident = 1000, # down sample to 1000 cells per ident
                                  logfc.threshold = 0.1)


DEG.betaVSbetaLike %>% head()
DEG.betaVSbetaLike <- DEG.betaVSbetaLike %>%  
  rownames_to_column("gene") %>% 
  as_tibble()

## add a colume to label beta_gene, beta_like_gene, or other
padj.cutoff <- 0.01
log2fc.cutoff <- 1
DEG.betaVSbetaLike <- DEG.betaVSbetaLike %>% 
  mutate(threshold = case_when(avg_log2FC > log2fc.cutoff & p_val_adj < padj.cutoff ~ "beta",
                               avg_log2FC < -log2fc.cutoff & p_val_adj < padj.cutoff ~ "beta_like",
                               TRUE ~ "other"))

table(DEG.betaVSbetaLike$threshold)

## 3. overlapping between DEG and beta.makers
idx <- DEG.betaVSbetaLike$gene[DEG.betaVSbetaLike$threshold == "beta"] %in% beta.markers$gene
gins.deficient.betaGene <- DEG.betaVSbetaLike$gene[DEG.betaVSbetaLike$threshold == "beta"] [idx]
gins.deficient.betaGene %>% 
  write.table(file = "output/gins.deficient.betaGene.txt", sep = "\t")
DEG.betaVSbetaLike$gene[DEG.betaVSbetaLike$threshold == "beta_like"] %>% 
  write.table(file = "output/beta.like.txt", sep = "\t")

# beta like vs beta cell volcano plot
## Volcano plot function
plotVol <- function(DEG.table,ceiling = 300) {
  #Move the points out of the range to the ceiling.
  table_tb <- DEG.table
  idx <- which(-log10(table_tb$p_val_adj) > ceiling)
  table_tb[idx,]$p_val_adj <- 10^(-ceiling)
  
  ggplot(table_tb) +
    geom_point(aes(x = avg_log2FC, 
                   y = -log10(p_val_adj), 
                   colour = threshold),
               size=1) +
    xlab("log2 fold change") + 
    ylab("-log10 adjusted p-value") +
    #scale_y_continuous(limits=c(0,40)) +
    #scale_x_continuous(limits=c(-5,5)) +
    theme_light() + 
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.border = element_rect(linetype = "solid", 
                                  fill = NA, 
                                  size = 0.5, 
                                  colour = "black"),
      legend.position = "none",
      plot.title = element_text(size = rel(1.5), hjust = 0.5),
      axis.title = element_text(size = rel(1.25))
    )
}

## volcano plot ----
idx <- DEG.betaVSbetaLike$gene[DEG.betaVSbetaLike$threshold == "beta"] %in% beta.markers$gene
background_beta_enriched <- DEG.betaVSbetaLike$gene[DEG.betaVSbetaLike$threshold == "beta"] [!idx]

DEGtableForVol <- DEG.betaVSbetaLike %>% filter(!gene %in% background_beta_enriched)
table(DEGtableForVol$threshold)
p1 <- plotVol(DEG.betaVSbetaLike, ceiling = 300) + 
  scale_color_manual(values=c("#ff595e", "#1982c4", "black")) +
  geom_hline(yintercept=-log10(padj.cutoff), linetype="dashed", 
             color = "black", size=0.2) +
  geom_vline(xintercept=c(-log2fc.cutoff,log2fc.cutoff), linetype="dashed", 
             color = "black", size=0.2) +
  ylim(0,300) +
  NULL
p1

ggsave(p1, filename = "figures_manuscript/ExtendedDataFig6a.jpeg", height = 5, width = 5, dpi = 400)

DEG.betaVSbetaLike.df <- DEG.betaVSbetaLike %>% column_to_rownames("gene")

genes.up <- c("RPS17","RPL17","RPL36A", "NGRN")
genes.down <- 
  c("RAMP1","CD81","PRKCI","ITGB1","ANXA2","PKP2","AGR2","CLTC", # protein localization to plasma membrane
    "TXNIP","CALM1","HSPA5", # response to calcium ion
    "DYNC1H1","CTSD","CPNE3","CST3", # neutrophil degranulation
    "ATP2B4","ATP1A1","FXYD3","DSP" # multicellular organismal signaling
  )

genes.to.label <- c(genes.up, genes.down)
EnhancedVolcano(DEG.betaVSbetaLike.df,
                lab = rownames(DEG.betaVSbetaLike.df),
                selectLab = genes.to.label,
                pCutoff = padj.cutoff,
                FCcutoff = log2fc.cutoff,
                labSize = 6,
                col = c('black', 'pink', 'purple4', 'red3'),
                drawConnectors = TRUE,
                boxedLabels = F,
                x = 'avg_log2FC',
                y = 'p_val_adj') +
  theme_classic() +
  theme(legend.position = "none",
        title = element_blank())

#### functional annotation ----
## GO term analysis beta vs beta_like ----
gins.deficient.beta.GO <- enrichGO(gene          = gins.deficient.betaGene,
                                   OrgDb         = org.Hs.eg.db,
                                   ont           = "BP",
                                   keyType       = "SYMBOL",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.1,
                                   qvalueCutoff  = 0.1)
gins.deficient.beta.GO@result %>% write_csv("functional_enrichment/ExtendedDataFig6b.gins.deficient.beta.GO.bp.csv")

bete.like.GO <- enrichGO(gene          = DEG.betaVSbetaLike$gene[DEG.betaVSbetaLike$threshold == "beta_like"],
                         OrgDb         = org.Hs.eg.db,
                         ont           = "BP",
                         keyType       = "SYMBOL",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.1,
                         qvalueCutoff  = 0.1)
bete.like.GO@result %>% write_csv("functional_enrichment/ExtendedDataFig6b.bete.like.GO.bp.csv")


# gene expression visualization for publication ----
## beta vs beta_like
DefaultAssay(integrated_cgins.vitro_islet) <- "RNA"
Idents(integrated_cgins.vitro_islet) <- "cell.type"
integrated_cgins.vitro_islet_beta <- subset(integrated_cgins.vitro_islet, 
                                            idents = c("beta","beta_like"))
# Vlnplot 
genes.to.plot <- c("INS", 
                   "GCK", "ABCC8", # exocytosis
                   "PAX6", "NKX2-2",
                   "PCSK1","G6PC2")
VlnPlot(integrated_cgins.vitro_islet_beta, 
        features = genes.to.plot, 
        stack = T, flip = T, 
        group.by = "source", 
        split.by = "source") + 
  NoLegend() +
  scale_fill_manual(values = c("#1982c4", # blue,
                               "#ff595e" # red
  )) +
  theme(axis.title.x = element_blank(),
        text = element_text(size = 20))
ggsave("figures_manuscript/Fig3c.pdf", width = 3.5, height = 6)


# expression of beta genes: beta VS betalike--------
DefaultAssay(integrated_cgins.vitro_islet_beta) <- "RNA"
source("helper_functions/TileDotplot.R")
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
exocytosis <- c("INS","STX1A","STX1B", "VAMP2", "SNAP25","SLC17A6","STX4","STXBP1","CDC42")
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

DefaultAssay(integrated_cgins.vitro_islet_beta) <- 'RNA'
for (i in seq_along(Vln.list)) {
  VlnPlot(integrated_cgins.vitro_islet_beta, 
          Vln.list[[i]], 
          stack = T, flip = T, 
          cols = ginsVSbeta, 
          fill.by = "ident" )  +
    theme(axis.title = element_blank(),
          axis.text.x = element_blank(),
          strip.text = element_text(size = 0),
          strip.background = element_blank()
    ) +
    NoLegend()
  file_name <- paste0("figures_manuscript/ExtendedDataFig4c.", names(Vln.list)[i], "_vln.pdf")
  ggsave(file_name, width = 1, height = length(Vln.list[[i]])/4.5)
}


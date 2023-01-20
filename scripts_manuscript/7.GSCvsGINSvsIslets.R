# Question: Is there a transcriptome-level adoption of beta cell fate.
## dataset loading and merging ----
cgins.vitro.time <- readRDS("rds_manuscript/cgins.vitro.time.rds")
cgins.vitro.time$source <- "timecourse"
integrated_islet_sub <- readRDS("rds_manuscript/integrated_islet.x_sub.rds")
cgins.vitro.time.list <- SplitObject(cgins.vitro.time, split.by = "orig.ident")
gsc.gins.beta <- merge(x = integrated_islet_sub, y = c(cgins.vitro.time.list$GSC, cgins.vitro.time.list$Late_GINS))

# subset dataset based on cell type
## gsc.gins.beta
Idents(gsc.gins.beta) <- "cell.type"
gsc.gins.beta <- gsc.gins.beta %>% subset(idents = c("beta","Beta_like", "Stem", "Mucoid"))
gsc.gins.beta$cell.type %>% table()
gsc.gins.beta <- 
  gsc.gins.beta %>%
  RenameIdents("beta" = "beta",
               "Mucoid" = "mucus",
               "Beta_like" = "beta_like",
               "Stem" = "stem")
gsc.gins.beta$cell.type <- Idents(gsc.gins.beta)

## processing and dimension reduction ----
DefaultAssay(gsc.gins.beta) <- "RNA"
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
gsc.gins.beta <- CellCycleScoring(gsc.gins.beta, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
gsc.gins.beta <- gsc.gins.beta %>% 
  # log-normalization 
  NormalizeData(normalization.method = "LogNormalize") %>% 
  # variable features
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  # scaling
  ScaleData(vars.to.regress = c("percent.mt", "S.Score", "G2M.Score") ) %>% 
  # PCA
  RunPCA()
#DimHeatmap(gsc.gins.beta, dims = 1:30, cells = 400, balanced = TRUE, ncol = 3)
#ElbowPlot(gsc.gins.beta, ndims = 30)
# umap and tsne
dimensions <- 1:10
gsc.gins.beta <- 
  gsc.gins.beta  %>% 
  # RunTSNE(dims = dimensions, seed.use = 5) %>% 
  RunUMAP(dims = dimensions)

gsc.gins.beta.color <- 
  c("beta" = "#ff595e", # red
    "beta_like"= "#1982c4", # blue
    "stem" = "#8ac926", # green
    "mucus" = "#ffca3a" # yellow
  )
gsc.gins.beta.order <- 
  c("mucus","stem","beta","beta_like") %>% rev()
p1 <- DimPlot(gsc.gins.beta, group.by = "cell.type", 
              reduction = "umap", 
              cols = gsc.gins.beta.color, 
              order = gsc.gins.beta.order) + 
  coord_equal() +
  NULL
p2 <- DimPlot(gsc.gins.beta, group.by = "orig.ident", 
              reduction = "umap", 
              cols = my9colors, 
              order = gsc.gins.beta.order) + 
  coord_equal() +
  NULL
p1+p2

## rds ----
saveRDS(gsc.gins.beta,"rds_manuscript/gsc.gins.beta.rds")
#gsc.gins.beta <- readRDS("rds_manuscript/gsc.gins.beta.rds")

## signature scoring using gene sets from public database ----
# get cell signature from msigdbr database
m_c8 <- msigdbr(species = "Homo sapiens", category = "C8") %>% 
  dplyr::select(gs_name, gene_symbol)
# extract gastric signature
idx <- m_c8$gs_name %>% str_which("GASTRIC")
m_c8_gastric <- m_c8[idx,]
# split gene symbols of different cell types into a list
gastric.list <- 
  m_c8_gastric %>% 
  split(f=m_c8_gastric$gs_name) %>% 
  lapply(function(x){x <- x$gene_symbol})
# remove unwanted list
gastric.signature <- gastric.list[c(2,5,7,9,11,13)] %>% unlist() %>% unique()
write.csv(gastric.signature,"gene_set/gastric.signature.csv")
gastric.signature <- read.csv("gene_set/gastric.signature.csv", row.names = 1)
gastric.signature <- gastric.signature$x

# extract beta cell signature
idx <- m_c8$gs_name %>% str_which("BETA")
m_c8_beta <- m_c8[idx,]
beta.list <- 
  m_c8_beta %>% 
  split(f=m_c8_beta$gs_name) %>% 
  lapply(function(x){x <- x$gene_symbol})
beta.signature <- beta.list[[1]]
write.csv(beta.signature,"gene_set/beta.signature.csv")
beta.signature <- read.csv("gene_set/beta.signature.csv", row.names = 1)
beta.signature <- beta.signature$x

msig.list <- list(gastric.signature=gastric.signature, beta.signature=beta.signature)
all.idents <- gsc.gins.beta$cell.type %>% unique() %>% as.character()
gene.expr <- get_gene_cutoff(gsc.gins.beta, idents = all.idents, pct.cutoff = 0.05)
msig.list.filtered <- msig.list %>% 
  lapply(
    function(x){
      idx <- x %in% gene.expr
      x <- x[idx]
    }
  )

gsc.gins.beta <- gsc.gins.beta %>% AddModuleScore(msig.list.filtered)
idx <- names(gsc.gins.beta@meta.data) %>% str_which("Cluster")
names(gsc.gins.beta@meta.data)[idx] <- names(msig.list.filtered)
Idents(gsc.gins.beta) <- "cell.type"
Idents(gsc.gins.beta) <- Idents(gsc.gins.beta) %>% factor(levels = c("mucus","stem","beta_like","beta"))
VlnPlot(gsc.gins.beta, 
        features = names(msig.list.filtered),
        cols = c("#BEBEBE","#cab2d6",'#fb9a99','#ff595e'),
        pt.size = 0,
        stack = T,
        flip = T,
        fill.by = "ident") + 
  theme_bw() +
  NoLegend() +
  theme(axis.title = element_blank(),
        strip.text = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
        ) +
  NULL
ggsave("figures_manuscript/Fig3d.Scorecards.pdf", width = 2, height = 3)

## disallowed genes visualization ----
source("helper_functions/TileDotplot.R")
disallowed.genes <- c("HK1","LDHA","SLC16A1")
TileDotplot(gsc.gins.beta, 
            features = (disallowed.genes),
            scale.by = "radius") +
  coord_fixed()
ggsave("figures_manuscript/ExtendedDataFig4d.disallowed.pdf", height = 4, width = 4)

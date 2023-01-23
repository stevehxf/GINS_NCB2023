## cgins.vitro and cgins.vivo integration ##
## comparing them with human islets ##

# loading seurat object rds ####
cgins.vitro <- readRDS(file = "rds_manuscript/cgins.vitro.rds")
cgins.vivo <- readRDS(file = "rds_manuscript/cgins.vivo.manuscript.rds")

#### cgins vitro and vivo integration ####
# combine list containing with cgins.vitro, cgins.vivo and islets from individual donor
list_cgins <- list(cgins.vitro = cgins.vitro, cgins.vivo = cgins.vivo)

# cgins log normalization
list_cgins <- list_cgins %>% 
  lapply(FUN = function(x) {
    DefaultAssay(x) <- "RNA"
    x <- x %>% 
      NormalizeData() %>% 
      FindVariableFeatures(selection.method = "mvp", nfeatures = 2000)
  })

#### integrated_cgins.vitro_vivo perform integration ####
cgins.features <- SelectIntegrationFeatures(object.list = list_cgins)
cgins.anchors <- FindIntegrationAnchors(object.list = list_cgins,
                                        anchor.features = cgins.features,
                                        k.anchor = 5)

integrated_cgins.vitro_vivo <- IntegrateData(anchorset = cgins.anchors)

#### integrated_cgins.vitro_vivo: the standard workflow for visualization and clustering  ####
DefaultAssay(integrated_cgins.vitro_vivo) <- "integrated"

integrated_cgins.vitro_vivo <- integrated_cgins.vitro_vivo %>% 
  ScaleData(vars.to.regress = "percent.mt", verbose = FALSE) %>% 
  RunPCA(verbose = FALSE)

DimHeatmap(integrated_cgins.vitro_vivo, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(integrated_cgins.vitro_vivo, dims = 16:30, cells = 500, balanced = TRUE)
ElbowPlot(integrated_cgins.vitro_vivo, ndims = 30)

dimensions <- 1:13
integrated_cgins.vitro_vivo <- 
  integrated_cgins.vitro_vivo %>% 
  RunTSNE(dims = dimensions, seed.use = 1)

DimPlot(integrated_cgins.vitro_vivo, cols = my36colors, group.by = "cell.type",split.by = "orig.ident")
DimPlot(integrated_cgins.vitro_vivo, reduction = "tsne", cols = my36colors, group.by = "RNA_snn_res.0.7",split.by = "orig.ident")

## clustering analysis
integrated_cgins.vitro_vivo <- integrated_cgins.vitro_vivo %>%   
  FindNeighbors(dims = dimensions)

for (resolution in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.8, 1.0, 2.0)) {
  integrated_cgins.vitro_vivo <- integrated_cgins.vitro_vivo %>%  
    FindClusters(resolution = resolution)}
DimPlot(integrated_cgins.vitro_vivo, reduction = "tsne", label = T, cols = my36colors, group.by = "RNA_snn_res.0.5", split.by = "orig.ident") + coord_equal()+ NoLegend()
DimPlot(integrated_cgins.vitro_vivo, reduction = "tsne", label = T, cols = my36colors, group.by = "integrated_snn_res.2", split.by = "orig.ident") + coord_equal()+ NoLegend()

## cluster visualization
p1 <- DimPlot(integrated_cgins.vitro_vivo, reduction = "tsne", label = T, cols = my36colors, group.by = "integrated_snn_res.0.1") + coord_equal()+ NoLegend()
p2 <- DimPlot(integrated_cgins.vitro_vivo, reduction = "tsne", label = T, cols = my36colors, group.by = "integrated_snn_res.0.2") + coord_equal()+ NoLegend()
p3 <- DimPlot(integrated_cgins.vitro_vivo, reduction = "tsne", label = T, cols = my36colors, group.by = "integrated_snn_res.0.3") + coord_equal()+ NoLegend()
p4 <- DimPlot(integrated_cgins.vitro_vivo, reduction = "tsne", label = T, cols = my36colors, group.by = "integrated_snn_res.0.4") + coord_equal()+ NoLegend()
p5 <- DimPlot(integrated_cgins.vitro_vivo, reduction = "tsne", label = T, cols = my36colors, group.by = "integrated_snn_res.0.5") + coord_equal()+ NoLegend()
p6 <- DimPlot(integrated_cgins.vitro_vivo, reduction = "tsne", label = T, cols = my36colors, group.by = "integrated_snn_res.0.8") + coord_equal()+ NoLegend()
p7 <- DimPlot(integrated_cgins.vitro_vivo, reduction = "tsne", label = T, cols = my36colors, group.by = "integrated_snn_res.1") + coord_equal()+ NoLegend()
p8 <- DimPlot(integrated_cgins.vitro_vivo, reduction = "tsne", label = T, cols = my36colors, group.by = "integrated_snn_res.2") + coord_equal()+ NoLegend()

plot_grid(p1,p2,p3,p4,p5,p6,p7,p8, nrow = 2)

DimPlot(integrated_cgins.vitro_vivo, reduction = "tsne", cols = my36colors, group.by = "cell.type", split.by = "orig.ident")
DimPlot(integrated_cgins.vitro_vivo, reduction = "tsne", cols = my36colors, group.by = "cell.type", split.by = "orig.ident")
feature.to.plot <- c("INS","GCG","SST","GHRL")
FeaturePlot(integrated_cgins.vitro_vivo, features = feature.to.plot)

#### integrated_cgins.vitro_vivo annotation ####
# annotation
Idents(integrated_cgins.vitro_vivo) <- "integrated_snn_res.2"
p1 <- DimPlot(integrated_cgins.vitro_vivo, reduction = "tsne", label = T, cols = my36colors, group.by = "cell.type", split.by = "orig.ident") + coord_equal()+ NoLegend()
p2 <- DimPlot(integrated_cgins.vitro_vivo, reduction = "tsne", label = T, cols = my36colors, split.by = "orig.ident") + coord_equal()+ NoLegend()
p1+p2

integrated_cgins.vitro_vivo <- 
  integrated_cgins.vitro_vivo %>% 
  RenameIdents(
    "9" = "delta_like",
    
    "12" = "alpha_like",
    
    "5" = "epsilon_like",
    "19" = "epsilon_like",
    
    "0" = "beta_like",
    "1" = "beta_like",
    "2" = "beta_like",
    "3" = "beta_like",
    "4" = "beta_like",
    "5" = "beta_like",
    "6" = "beta_like",
    "7" = "beta_like",
    "8" = "beta_like",
    "10" = "beta_like",
    "11" = "beta_like",
    "13" = "beta_like", 
    "14" = "beta_like",
    "15" = "beta_like",
    "16" = "beta_like",
    "17" = "beta_like",
    "18" = "beta_like",
    "20" = "beta_like",
    "21" = "beta_like")


table(Idents(integrated_cgins.vitro_vivo))

# store the annotation
integrated_cgins.vitro_vivo$integrated.cell.type <- Idents(integrated_cgins.vitro_vivo)

## RDS ---- 
saveRDS(integrated_cgins.vitro_vivo, "rds_manuscript/integrated_cgins.vitro_vivo.rds" )
#integrated_cgins.vitro_vivo <- readRDS("rds_manuscript/integrated_cgins.vitro_vivo.rds" )

# integration visualization for publications ----
## tsne
DimPlot(integrated_cgins.vitro_vivo, 
        reduction = "tsne", 
        group.by = "integrated.cell.type",
        cols = vivo.vitro.cell.color, 
        #order = inte.islet.order,
        label = F,
        split.by = "orig.ident") + 
  theme_void() +
  theme(plot.title = element_blank(),
        aspect.ratio = 1) +
  NoLegend() +
  NULL
ggsave("figures_manuscript/ExtendedDataFig7f_integrated.vitro.vivo.tsne.jpeg", dpi = 400, width = 3.5,height = 1.7)

## proportion bar graph
Idents(integrated_cgins.vitro_vivo) <- "integrated.cell.type"
cell.type.stats <- table(integrated_cgins.vitro_vivo$orig.ident, 
                         integrated_cgins.vitro_vivo$integrated.cell.type) %>% 
  as.data.frame()
names(cell.type.stats) <- c("orig.ident","integrated.cell.type","freq")
cell.type.stats %>% ggplot(aes(x=orig.ident, y=freq, fill=integrated.cell.type)) +
  geom_col(position = "fill", width = 0.1) +
  scale_fill_manual(values = vivo.vitro.cell.color) +
  theme_void() +
  theme(legend.position = "none")
ggsave("figures_manuscript/ExtendedDataFig7f_integrated.cell.type.prop.pdf", width = 2, height = 2)

## feature plot
DefaultAssay(integrated_cgins.vitro_vivo) <- 'RNA'
FeaturePlot(integrated_cgins.vitro_vivo, "INS", 
            split.by = "orig.ident", 
            coord.fixed = TRUE,
            min.cutoff = 3) &
  scale_color_gradientn(colours = viridis(n=100, option = "A")) &
  theme_void() &
  theme(text = element_blank(),
        aspect.ratio = 1) +
  NoLegend()
ggsave("figures_manuscript/INS.vivo.vitro.jpeg", dpi = 400, width = 4, height = 2)

## vlnplot of beta-like
gene.of.interest <- c("INS","GCG","SST","GHRL","GAL",
                      "ENTPD3", "UCN3", "PAX6","SLC30A8","NKX2.2",
                      "CHGA","GCK")
cell.type.vitroVSvivo.DEG.df <- cell.type.vitroVSvivo.DEG.df %>% 
  mutate(label=gene %in% gene.of.interest) %>% 
  arrange(label) %>% 
  mutate(genelabels = "")
cell.type.vitroVSvivo.DEG.df[cell.type.vitroVSvivo.DEG.df$label,]$genelabels <- 
  cell.type.vitroVSvivo.DEG.df[cell.type.vitroVSvivo.DEG.df$label,]$gene
beta.vitroVSvivo.DEG.df <- cell.type.vitroVSvivo.DEG.df %>% filter(cluster == "beta_like")

EnhancedVolcano(beta.vitroVSvivo.DEG.df,
                lab = beta.vitroVSvivo.DEG.df$gene,
                selectLab = gene.of.interest,
                pCutoff = 10e-6,
                labSize = 6,
                col = c('black', 'pink', 'purple4', 'red3'),
                FCcutoff = 0.3,
                drawConnectors = TRUE,
                boxedLabels = F,
                x = 'avg_log2FC',
                y = 'p_val_adj') +
  theme_classic() +
  theme(legend.position = "none",
        title = element_blank())

ggsave("figures_manuscript/vivo.vitro.beta.vol.pdf", width = 4, height = 3)

# integrate cgins.vitro, cgins.vivo and islets dataset

# loading seurat object rds and create combined object ####
integrated_cgins.vitro_vivo <- readRDS("rds_manuscript/integrated_cgins.vitro_vivo.rds" )
integrated_islet_sub <- readRDS("rds_manuscript/integrated_islet.x_sub.rds")
# create a list containing cgins.vitro and cgins.vivo and islets and normalize
integrated_islet_sub$orig.ident <- paste0(integrated_islet_sub$donor, "_islets")
integrated_islet_sub$source <- "islets"
integrated_islet_sub$integrated.cell.type <- integrated_islet_sub$cell.type
cgins.list <- SplitObject(integrated_cgins.vitro_vivo, split.by = "orig.ident")
cgins.list$`C-GINSC-VIVO`$source <- "graft"
islet.list <- SplitObject(integrated_islet_sub, split.by = "orig.ident")
list_gins.vitro.vivo.islets <- c(cgins.list, islet.list)
names(list_gins.vitro.vivo.islets)[c(1,2)] <- c("gins_organoids","gins_grafts") 
gins.vitro.vivo.islets <- merge(x=list_gins.vitro.vivo.islets$gins_organoids,
                                y=c(list_gins.vitro.vivo.islets$gins_grafts,
                                    list_gins.vitro.vivo.islets$Donor_1_islets,
                                    list_gins.vitro.vivo.islets$Donor_2_islets,
                                    list_gins.vitro.vivo.islets$Donor_3_islets,
                                    list_gins.vitro.vivo.islets$Donor_9_islets
                                    )
)

## subset beta and beta-like cells for trajectory analysis ----
# subset
Idents(gins.vitro.vivo.islets) <- "cell.type"
beta.betaLike <- subset(gins.vitro.vivo.islets, idents = c("beta", "beta_like"))
Idents(beta.betaLike) <- "source"
DefaultAssay(beta.betaLike) <- 'RNA'

beta.betaLike <- beta.betaLike %>% 
  NormalizeData() %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000)
beta.betaLike <- beta.betaLike %>% 
  ScaleData(#features = var.features.combined, 
    vars.to.regress = "percent.mt", verbose = FALSE) 

#saveRDS(beta.betaLike,"rds_manuscript/beta.betaLike.rds")
#beta.betaLike <- readRDS("rds_manuscript/beta.betaLike.rds")

# correlation coefficient of downsampled beta/beta-like cells  ----

Idents(gins.vitro.vivo.islets) <- "cell.type"
gins.vitro.vivo.islets.beta <- subset(gins.vitro.vivo.islets, idents=c("beta","beta_like"))
gins.vitro.vivo.islets.beta$cell.type %>% table()

## get top 2000 variable genes data
DefaultAssay(gins.vitro.vivo.islets.beta) <- "RNA"
gins.vitro.vivo.islets.beta <- 
  gins.vitro.vivo.islets.beta %>% 
  NormalizeData() 

gins.vitro.vivo.islets.beta <- 
  gins.vitro.vivo.islets.beta %>% 
  FindVariableFeatures(selection.method = "mvp", nfeatures = 2000)
gins.vitro.vivo.islets.beta <- 
  gins.vitro.vivo.islets.beta %>% 
  ScaleData(vars.to.regress = "percent.mt", verbose = FALSE)
variable.features <- gins.vitro.vivo.islets.beta@assays$RNA@var.features

## downsample group by source
Idents(gins.vitro.vivo.islets.beta) <- "source"
gins.vitro.vivo.islets.beta.bySource <- gins.vitro.vivo.islets.beta %>% SplitObject(split.by = "ident")
data.list <- vector(mode = "list", length = length(gins.vitro.vivo.islets.beta.bySource))
data.cor.list <- data.list
meta.list <- data.list
for (i in seq_along(gins.vitro.vivo.islets.beta.bySource)) {
  gins.vitro.vivo.islets.beta.bySource[[i]] <- gins.vitro.vivo.islets.beta.bySource[[i]] %>% subset(downsample=300)
  data.list[[i]] <- gins.vitro.vivo.islets.beta.bySource[[i]]@assays$RNA@data %>% as.data.frame()
  data.list[[i]] <- data.list[[i]][variable.features,]
  data.cor.list[[i]] <- data.list[[i]] %>% cor(method = "pearson")
}

legend_breaks <-  c(0,1)
breaks <- seq(0, 1, length.out = 100)
color <- rev(colorRampPalette(RColorBrewer::brewer.pal(9,"RdYlBu"))(100))

clustering_method  <-  "ward.D"
cellwidth <- 1
cellheight <- 1
border_color <- NA

p <- pheatmap(data.cor.list[[1]],
               color = color,
               clustering_method = clustering_method,
               scale = "none",
               legend_breaks = legend_breaks,
               breaks = breaks,
               legend = FALSE,
               treeheight_col = 0,
               treeheight_row = 0,
               cellwidth = cellwidth,
               cellheight = cellheight,
               border_color = border_color,
               show_rownames = F,
               show_colnames = F)
plot_grid(p$gtable)
ggsave("figures_manuscript/Fig4g.cor.gins.organoid.tiff", width = 3, height = 3, dpi = 400,  bg = "transparent")

p <- pheatmap(data.cor.list[[2]],
              color = color,
              clustering_method = clustering_method,
              scale = "none",
              legend_breaks = legend_breaks,
              breaks = breaks,
              legend = FALSE,
              treeheight_col = 0,
              treeheight_row = 0,
              cellwidth = cellwidth,
              cellheight = cellheight,
              border_color = border_color,
              show_rownames = F,
              show_colnames = F)
plot_grid(p$gtable)
ggsave("figures_manuscript/Fig4g.cor.gins.graft.tiff", width = 3, height = 3, dpi = 400,  bg = "transparent")

p <- pheatmap(data.cor.list[[3]],
              color = color,
              clustering_method = clustering_method,
              scale = "none",
              legend_breaks = legend_breaks,
              breaks = breaks,
              legend = FALSE,
              treeheight_col = 0,
              treeheight_row = 0,
              cellwidth = cellwidth,
              cellheight = cellheight,
              border_color = border_color,
              show_rownames = F,
              show_colnames = F)
plot_grid(p$gtable)
ggsave("figures_manuscript/Fig4g.cor.islet.tiff", width = 3, height = 3, dpi = 400,  bg = "transparent")

## genes for vlnplot
gene.of.interest <- c("TFF2", "GHRL",
                      "NEUROD1",'NKX2-2','PAX6', 'UCN3','ENTPD3', "INS")
VlnPlot(beta.betaLike, features = rev(gene.of.interest), pt.size = 0, 
        cols = c("#1982c4", "#ff595e", "grey"),
        fill.by = "ident",
        stack = T, flip = T) + 
  NoLegend()

ggsave("figures_manuscript/Fig4f.vlnplot.pdf", width = 3, height = 5)

## GAL vlnplot
Idents(beta.betaLike) <- "source"
beta.betaLike.sub <- subset(beta.betaLike,downsample=2007)
VlnPlot(beta.betaLike.sub, features = "GAL", pt.size = 0.1, 
        cols = c("#1982c4", "#ff595e", "grey"),
        fill.by = "ident") + 
  NoLegend()
ggsave("figures_manuscript/gal.gins.graft.islet1.png", width = 5, height = 5, dpi = 400,  bg = "transparent")

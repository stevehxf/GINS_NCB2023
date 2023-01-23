## dataset loading and merging ----
cgins.vitro.time <- readRDS("rds_manuscript/cgins.vitro.time.rds")
gsc.vel <- ReadVelocity("data/cgins_timeCourse_loom/Sample1.loom") %>% as.Seurat()
endo.vel <- ReadVelocity("data/cgins_timeCourse_loom/Sample2.loom") %>% as.Seurat()
ginsPre.vel <- ReadVelocity("data/cgins_timeCourse_loom/Sample3.loom") %>% as.Seurat()
gins.vel <- ReadVelocity("data/cgins_timeCourse_loom/Sample4.loom") %>% as.Seurat()
cgins.vitro.time.vel <- merge(gsc.vel, c(endo.vel, ginsPre.vel, gins.vel))
cgins.vitro.time.muRM <- readRDS("rds_manuscript/cgins.vitro.time.muRM.rds")

## subset the cells of cgins.vitro.time.vel according to cgins.vitro.time ----
# rename
name_replaced <- 
  colnames(cgins.vitro.time.vel) %>% 
  str_replace("x","-1") %>% 
  str_remove("Sample1:") %>% 
  str_remove("Sample2:") %>% 
  str_remove("Sample3:") %>% 
  str_remove("Sample4:") 
cgins.vitro.time.vel <- cgins.vitro.time.vel %>% 
  RenameCells(new.names = name_replaced)
# subset
name_subset <- colnames(cgins.vitro.time)
cgins.vitro.time.vel <- subset(cgins.vitro.time.vel, cells = name_subset)

# remove mouse genes except Ngn3, Pdx1 and Mafa
hg.gene <- rownames(cgins.vitro.time.vel) %>% str_subset("GRCh38-")
cgins.vitro.time.vel <- subset(cgins.vitro.time.vel, features = c(hg.gene, "mm10---Pdx1", "mm10---Mafa", "mm10---Neurog3"))

# rename the genes
re.name <- row.names(cgins.vitro.time.vel) %>% str_remove("GRCh38-")
tail(re.name)
length(re.name)
re.name[c((length(re.name)-2):length(re.name))] <- c("PDX1MUS","MAFAMUS","NGN3ERMUS")

cgins.vitro.time.vel@assays$spliced@counts@Dimnames[[1]] <- re.name
cgins.vitro.time.vel@assays$spliced@data@Dimnames[[1]] <- re.name
rownames(cgins.vitro.time.vel@assays$spliced@meta.features) <- re.name

cgins.vitro.time.vel@assays$unspliced@counts@Dimnames[[1]] <- re.name
cgins.vitro.time.vel@assays$unspliced@data@Dimnames[[1]] <- re.name
rownames(cgins.vitro.time.vel@assays$unspliced@meta.features) <- re.name

cgins.vitro.time.vel@assays$ambiguous@counts@Dimnames[[1]] <- re.name
cgins.vitro.time.vel@assays$ambiguous@data@Dimnames[[1]] <- re.name
rownames(cgins.vitro.time.vel@assays$ambiguous@meta.features) <- re.name

## Add metadata from cgins.vitro.time to cgins.vitro.time.vel
cgins.vitro.time.vel <- cgins.vitro.time.vel %>% 
  AddMetaData(metadata = cgins.vitro.time@meta.data)

## add spliced/unspliced assay back to rmMu seurat object
cgins.vitro.time.muRM <- readRDS("rds_manuscript/cgins.vitro.time.muRM.rds")
Idents(cgins.vitro.time.vel) <- "cell.type"
cgins.vitro.time.muRM.vel <- cgins.vitro.time.vel %>% subset(idents = "Mucoid", invert = TRUE)
all(colnames(cgins.vitro.time.muRM.vel) == colnames(cgins.vitro.time.muRM))
cgins.vitro.time.muRM[["spliced"]] <- cgins.vitro.time.muRM.vel[["spliced"]]
cgins.vitro.time.muRM[["unspliced"]] <- cgins.vitro.time.muRM.vel[["unspliced"]]
cgins.vitro.time.muRM[["ambiguous"]] <- cgins.vitro.time.muRM.vel[["ambiguous"]]

DefaultAssay(cgins.vitro.time.muRM) <- "spliced"
cgins.vitro.time.muRM <- cgins.vitro.time.muRM %>% RunVelocity(deltaT = 1, kCells = 25, fit.quantile = 0.02)

ident.colors <- timecourse.cell.type.color
cell.colors <- ident.colors[as.character(cgins.vitro.time.muRM$cell.type)]
names(x = cell.colors) <- colnames(x = cgins.vitro.time.muRM)

pdf("figures_manuscript/velocity.pdf", width = 8, height = 8)
show.velocity.on.embedding.cor(emb = Embeddings(object = cgins.vitro.time.muRM, reduction = "umap"), 
                               vel = Tool(object = cgins.vitro.time.muRM, 
                                          slot = "RunVelocity"), 
                               n = 200, 
                               scale = "sqrt", 
                               cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, 
                               arrow.scale = 10, 
                               show.grid.flow = TRUE,
                               min.grid.cell.mass = 0.5, 
                               grid.n = 70, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)
dev.off()
### rds
cgins.vitro.time.muRM %>% saveRDS("rds_manuscript/cgins.vitro.time.muRM.rds")

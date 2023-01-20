## server setup ####
library(tidyverse)
library(Seurat)
library(SingleR)
library(SoupX)
library(DoubletFinder)
#library(monocle)
library(monocle3)
library(EnhancedVolcano)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(annotables)
library(enrichplot)
library(DOSE)
library(RColorBrewer)
#library(velocyto.R)
library(ggrepel)
library(Seurat)
library(SeuratWrappers)
library(cowplot)
library(patchwork)
library(viridis)
library(pheatmap)
library(RcisTarget)
library(AUCell)
library(GENIE3)
library(SCENIC)
devtools::load_all("C:/Users/zhoulab/AppData/Local/R/win-library/4.2/monocle/")

## in-house functions
source("scripts_manuscript/helper_function.r")

## colors of collection ------------------------------------------------------
my36colors <- c("#E5D2DD", "#53A85F", "#F1BB72", "#F3B1A0", "#D6E7A3", "#57C3F3", "#476D87",
                "#E95C59", "#E59CC4", "#AB3282", "#23452F", "#BD956A", "#8C549C", "#585658",
                "#9FA3A8", "#E0D4CA", "#5F3D69", "#C5DEBA", "#58A4C3", "#E4C755", "#F7F398",
                "#AA9A59", "#E63863", "#E39A35", "#C1E6F3", "#6778AE", "#91D0BE", "#B53E2B",
                "#712820", "#DCC1DD", "#CCE0F5", "#CCC9E6", "#625D9E", "#68A180", "#3A6963",
                "#968175")
my9colors <- c("#fdbf6f", # yellow
               "#a6cee3", # soft blue
               "#0096c7", # cyan
               "#b2df8a", # lightgreen
               "#33a02c", # darkgreen
               "#cab2d6", # grayish violet
               "#fb9a99", # pink
               "#e63946", # brightred
               "#6a3d9a"  # moderate violet
)
my5colors <- c("#ff595e", # red
               "#ffca3a", # yellow
               "#8ac926", # green
               "#1982c4", # blue
               "#6a4c93" # purple
)

#
grad <- colorRampPalette(RColorBrewer::brewer.pal(5,"RdPu"))(50) 
hmcols <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100)
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(100)

## gins cell type color
gins.color <- c("beta_like" = "#ff595e", # red
                "alpha_like" = "#ffca3a", # yellow
                "delta_like" = "#8ac926", # green
                "epsilon_like" = "#1982c4", # blue
                "G_like" = "#6a4c93" # purple"
)

orig.color <- c(
  "Late_GINS" = brewer.pal(12, "Paired")[6], # brightred
  "Early_GINS" = brewer.pal(12, "Paired")[2], #  blue,
  "Endocrine" = brewer.pal(12, "Paired")[4], # darkgreen
  "GSC" = brewer.pal(12, "Paired")[7], # light orange
  "SeuratProject" = "grey" 
)

timecourse.cell.type.color <- c(
  "Mucoid"="#cab2d6", # grayish violet
  "Stem"="#fdbf6f", # light orange
  "Endocrine_1"="#b2df8a", # lightgreen"
  "Endocrine_2"="#33a02c", # darkgreen
  "Delta_like"="#a6cee3", # soft blue
  "GINS_precursors"="#0096c7", # cyan
  "Beta_like"="#e63946", # brightred
  "Epsilon_like"="#fb9a99", # pink
  "Alpha_like"="#cab2d6"  # moderate violet
)

islet.timecourse.cell.type.color <- c(
  
  "Stem" = brewer.pal(12, "Paired")[7], # light orange
  "Endocrine_1" = brewer.pal(12, "Paired")[3], # lightgreen"
  "Endocrine_2" = brewer.pal(12, "Paired")[4], # darkgreen

  "GINS_precursors" = brewer.pal(12, "Paired")[2], #  blue,
  
  "beta" = brewer.pal(12, "Paired")[5], # pink
  "Beta_like" = brewer.pal(12, "Paired")[6], # brightred
  
  
  "Alpha_like" = brewer.pal(12, "Paired")[10],  # moderate violet
  "alpha" = brewer.pal(12, "Paired")[9], #  violet
  
  "Delta_like" = brewer.pal(12, "Paired")[1], # soft blue
  "delta" = brewer.pal(12, "Paired")[12], # brown
  
  "Epsilon_like" = "grey" # grayish violet
)


timecourseGSC.GINS.cell.type.color <- 
  c(
    "Mucoid"="#BEBEBE", # grey
    "Stem"="#6a4c93", # light orange
    "Beta_like" = "#ff595e", # red
    "Alpha_like" = "#ffca3a", # yellow
    "Delta_like" = "#8ac926", # green
    "Epsilon_like" = "#1982c4"
    
  )

ginsVSbeta <- c("beta_like" = "#1982c4", # red
                "beta" = "#ff595e" # blue
)

cVSa <- c("agins.vitro" = "#1982c4", # red
          "C-GINSC" = "#ff595e" # blue
)

vitroVSvivo <- c("C-GINSC-VIVO" = "#1982c4", # red
                 "C-GINSC" = "#ff595e" # blue
)
  
vivo.vitro.cell.color <- c("beta_like" = "#EC7E7A", # subtle red
                           "alpha_like" = "#ffca3a", # yellow
                           "delta_like" = "#8ac926", # green
                           "epsilon_like" = "#1982c4" # blue
)


inte.color <- c(
  alpha = "#ffca3a", # yellow
  alpha_like = "#ffca3a", # yellow
  beta = "#ff595e", # red
  beta_like = "#ff595e", # red
  delta = "#8ac926", # green
  delta_like = "#8ac926", # green
  epsilon = "#1982c4", # blue
  epsilon_like = "#1982c4", # blue
  pp = "#6a4c93"
)

# heatmap annotation color for regulon
# color
ann_colors <-  list(
  cell.type = c(  "Stem"="#fdbf6f", # light orange
                  "Endocrine_1"="#b2df8a", # lightgreen"
                  "Endocrine_2"="#33a02c", # darkgreen
                  "Delta_like"="#a6cee3", # soft blue
                  "GINS_precursors"="#0096c7", # cyan
                  "Beta_like"="#e63946", # brightred
                  "Epsilon_like"="#fb9a99", # pink
                  "Alpha_like"="#cab2d6"  # moderate violet
  ),
  Pseudotime = colorRampPalette(RColorBrewer::brewer.pal(6,"OrRd"))(25)
)

# directories for output ----
## create a directory for the dataset
dir <- "figures_manuscript"
if(!exists(dir)){
  dir.create(dir)
}  


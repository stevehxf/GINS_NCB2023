# GINS_NCB

###### Author: Xiaofeng Huang<sup>1,2,4</sup>, Qiao Zhou<sup>1,3,4</sup>
###### <sup>1</sup> Division of Regenerative Medicine & Ansary Stem Cell Institute, Department of Medicine, Weill Cornell Medicine, 1300 York Avenue, New York, NY 10065, USA 
###### <sup>2</sup> E-mail: xih4001@med.cornell.edu
###### <sup>3</sup> E-mail: jqz4001@med.cornell.edu 
###### <sup>4</sup> Corresponding

## Introduction
The repo contains all the data and scripts you need to reproduce the figures in the paper.

### Publication associated with this repo
Stomach-derived human insulin-secreting organoids restore glucose homeostasis, Nature Cell Biology, 2023, DOI:

## Scripts
#### NOTE: Run Script 1 to load all the required libraries and helper scripts.
* Script 1: libraries loading; scripts loading; directories; colors
* Script 2: data processing for GINS organoids
* Script 3: data processing for islet
* Script 4: integrating GINS organoids and human islets.
* Script 5: data processing for GINS differentiation time course
* Script 6: comparing hGSC and GINS organoids
* Script 7: comparing hGSC, GINS organoids, and islets
* Script 8: comparing antral and corpus GINS organoids
* Script 9: data processing for GINS graft
* Script 10: comparison GINS organoids, GINS grafts and human islets
* Script 11: visualization of GINS differentiation trajectory
* Script 12: visualization of GINS differentiation trajectory with islet datasets
* Script 13: Regulon analysis

## Figures output
* Script 4: Fig. 3a-c; Extended Data Fig. 4c, 6a
* Script 6: Extended Data Fig. 4a, 4b, 4e 
* Script 7: Fig. 3d; Extended Data Fig. 4d
* Script 8: Extended Data Fig 5d, 5e
* Script 10: Fig. 4f-g, 6g; Extended Data Fig 7f
* Script 11: Fig. 5c, 6a; Extended Data Fig 9a-c
* Script 12: Fig. 5a-b,d-f; Extended Data Fig 8a-c

## Data output (no graph plotting)
* Script 4: Extended Data Fig. 6b (pathway enrichment)
* Script 7: Fig. 3d (Signature gene sets)
* Script 12: Fig. 5f (pathway enrichment)

## SessionInfo 
R version 4.2.2 (2022-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19042)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
[3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.utf8    

attached base packages:
[1] splines   stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] monocle_2.24.0              testthat_3.1.4              DDRTree_0.1.5              
 [4] irlba_2.3.5                 VGAM_1.1-7                  Matrix_1.5-1               
 [7] SCENIC_1.3.1                GENIE3_1.18.0               AUCell_1.18.1              
[10] RcisTarget_1.16.0           pheatmap_1.0.12             viridis_0.6.2              
[13] viridisLite_0.4.1           patchwork_1.1.2             cowplot_1.1.1              
[16] SeuratWrappers_0.3.0        RColorBrewer_1.1-3          DOSE_3.22.0                
[19] enrichplot_1.16.1           annotables_0.1.91           org.Hs.eg.db_3.15.0        
[22] AnnotationDbi_1.58.0        clusterProfiler_4.4.4       msigdbr_7.5.1              
[25] EnhancedVolcano_1.14.0      ggrepel_0.9.1               monocle3_1.2.9             
[28] SingleCellExperiment_1.18.0 DoubletFinder_2.0.3         SoupX_1.6.1                
[31] SingleR_1.10.0              SummarizedExperiment_1.26.1 Biobase_2.56.0             
[34] GenomicRanges_1.48.0        GenomeInfoDb_1.32.3         IRanges_2.30.1             
[37] S4Vectors_0.34.0            BiocGenerics_0.42.0         MatrixGenerics_1.8.1       
[40] matrixStats_0.62.0          sp_1.5-0                    SeuratObject_4.1.0         
[43] Seurat_4.1.1                forcats_0.5.2               stringr_1.4.1              
[46] dplyr_1.0.9                 purrr_0.3.4                 readr_2.1.2                
[49] tidyr_1.2.0                 tibble_3.1.8                ggplot2_3.3.6              
[52] tidyverse_1.3.2            

loaded via a namespace (and not attached):
  [1] rsvd_1.0.5                ica_1.0-3                 ps_1.7.1                 
  [4] lmtest_0.9-40             rprojroot_2.0.3           crayon_1.5.1             
  [7] spatstat.core_2.4-4       MASS_7.3-58.1             nlme_3.1-160             
 [10] backports_1.4.1           qlcMatrix_0.9.7           reprex_2.0.2             
 [13] GOSemSim_2.22.0           rlang_1.0.4               XVector_0.36.0           
 [16] ROCR_1.0-11               readxl_1.4.1              nloptr_2.0.3             
 [19] callr_3.7.2               limma_3.52.2              BiocParallel_1.30.3      
 [22] bit64_4.0.5               glue_1.6.2                sctransform_0.3.4        
 [25] vipor_0.4.5               parallel_4.2.2            processx_3.7.0           
 [28] spatstat.sparse_2.1-1     spatstat.geom_2.4-0       haven_2.5.0              
 [31] tidyselect_1.1.2          usethis_2.1.6             fitdistrplus_1.1-8       
 [34] XML_3.99-0.10             zoo_1.8-10                xtable_1.8-4             
 [37] magrittr_2.0.3            cli_3.3.0                 zlibbioc_1.42.0          
 [40] rstudioapi_0.13           miniUI_0.1.1.1            rpart_4.1.19             
 [43] fastmatch_1.1-3           treeio_1.20.2             shiny_1.7.2              
 [46] BiocSingular_1.12.0       leidenbase_0.1.11         pkgbuild_1.3.1           
 [49] cluster_2.1.4             tidygraph_1.2.1           KEGGREST_1.36.3          
 [52] ape_5.6-2                 listenv_0.8.0             Biostrings_2.64.1        
 [55] png_0.1-7                 future_1.27.0             withr_2.5.0              
 [58] bitops_1.0-7              slam_0.1-50               ggforce_0.3.4            
 [61] RBGL_1.72.0               plyr_1.8.7                cellranger_1.1.0         
 [64] GSEABase_1.58.0           sparsesvd_0.2             pillar_1.8.1             
 [67] biocViews_1.64.1          cachem_1.0.6              fs_1.5.2                 
 [70] RUnit_0.4.32              DelayedMatrixStats_1.18.0 vctrs_0.4.1              
 [73] ellipsis_0.3.2            generics_0.1.3            devtools_2.4.4           
 [76] tools_4.2.2               beeswarm_0.4.0            munsell_0.5.0            
 [79] tweenr_2.0.0              fgsea_1.22.0              proxy_0.4-27             
 [82] DelayedArray_0.22.0       fastmap_1.1.0             compiler_4.2.2           
 [85] HSMMSingleCell_1.16.0     pkgload_1.3.0             abind_1.4-5              
 [88] httpuv_1.6.5              sessioninfo_1.2.2         plotly_4.10.0            
 [91] rgeos_0.5-9               GenomeInfoDbData_1.2.8    gridExtra_2.3            
 [94] lattice_0.20-45           deldir_1.0-6              utf8_1.2.2               
 [97] later_1.3.0               jsonlite_1.8.0            scales_1.2.1             
[100] docopt_0.7.1              graph_1.74.0              ScaledMatrix_1.4.0       
[103] tidytree_0.4.0            pbapply_1.5-0             sparseMatrixStats_1.8.0  
[106] lazyeval_0.2.2            promises_1.2.0.1          R.utils_2.12.0           
[109] goftest_1.2-3             spatstat.utils_2.3-1      reticulate_1.25          
[112] Rtsne_0.16                downloader_0.4            uwot_0.1.13              
[115] igraph_1.3.4              survival_3.4-0            htmltools_0.5.3          
[118] memoise_2.0.1             profvis_0.3.7             graphlayouts_0.8.1       
[121] arrow_9.0.0               digest_0.6.29             assertthat_0.2.1         
[124] mime_0.12                 RSQLite_2.2.16            yulab.utils_0.0.5        
[127] future.apply_1.9.0        remotes_2.4.2             data.table_1.14.2        
[130] urlchecker_1.0.1          blob_1.2.3                R.oo_1.25.0              
[133] labeling_0.4.2            fastICA_1.2-3             googledrive_2.0.0        
[136] RCurl_1.98-1.8            broom_1.0.0               hms_1.1.2                
[139] modelr_0.1.9              colorspace_2.0-3          BiocManager_1.30.18      
[142] ggbeeswarm_0.6.0          aplot_0.1.6               ggrastr_1.0.1            
[145] Rcpp_1.0.9                RANN_2.6.1                fansi_1.0.3              
[148] tzdb_0.3.0                brio_1.1.3                parallelly_1.32.1        
[151] R6_2.5.1                  grid_4.2.2                ggridges_0.5.3           
[154] lifecycle_1.0.1           googlesheets4_1.0.1       minqa_1.2.4              
[157] leiden_0.4.2              DO.db_2.9                 qvalue_2.28.0            
[160] desc_1.4.1                RcppAnnoy_0.0.19          htmlwidgets_1.5.4        
[163] beachmat_2.12.0           polyclip_1.10-0           shadowtext_0.1.2         
[166] gridGraphics_0.5-1        terra_1.6-7               rvest_1.0.3              
[169] mgcv_1.8-41               globals_0.16.0            spatstat.random_2.2-0    
[172] progressr_0.10.1          codetools_0.2-18          lubridate_1.8.0          
[175] GO.db_3.15.0              prettyunits_1.1.1         dbplyr_2.2.1             
[178] R.methodsS3_1.8.2         gtable_0.3.0              DBI_1.1.3                
[181] ggfun_0.0.6               tensor_1.5                httr_1.4.4               
[184] KernSmooth_2.23-20        stringi_1.7.8             reshape2_1.4.4           
[187] farver_2.1.1              annotate_1.74.0           ggtree_3.4.2             
[190] xml2_1.3.3                combinat_0.0-8            boot_1.3-28              
[193] BiocNeighbors_1.14.0      lme4_1.1-30               ggplotify_0.1.0          
[196] scattermore_0.8           bit_4.0.4                 scatterpie_0.1.7         
[199] spatstat.data_2.2-0       ggraph_2.0.6              pkgconfig_2.0.3          
[202] babelgene_22.3            gargle_1.2.0      

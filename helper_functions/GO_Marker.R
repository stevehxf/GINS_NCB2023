# DEG.table: input is the raw DEG.table output from FindMarkers function
# padj.cutoff: adjusted p value cutoff for significant DEGs
# log2fc.cutoff: log2 fold change cutoff for significant DEGs
# simplify.cutoff: cutoff to remove redundancy of GO terms
# species: either mouse or human
# ont: ontology from GO terms

GO_Marker <- function(DEG.table, padj.cutoff, log2fc.cutoff, simplify.cutoff, species, ont){
        
        require(clusterProfiler)
        require(tidyverse)
        require(org.Hs.eg.db)
        require(org.Mm.eg.db)
        require(annotables)
        require(enrichplot)
        
        ## begin
        cat('Processing: ')
        cat(paste(unique(DEG.table$cluster), "\n"))
        ## extract signifant up-regulated genes
        up_sig <- DEG.table %>% 
                filter(p_val_adj < padj.cutoff) %>% 
                filter(avg_log2FC > log2fc.cutoff)
        up_sig <- rownames(up_sig)
        cat('Up-regulated: ')
        print(up_sig)
        
        ## extract signifant down-regulated genes
        down_sig <- DEG.table %>% 
                filter(p_val_adj < padj.cutoff) %>% 
                filter(avg_log2FC < -log2fc.cutoff)
        down_sig <- rownames(down_sig)
        cat('Down-regulated: ')
        print(down_sig)
        
        ## species
        if(species == "mouse"){
                OrgDb <- "org.Mm.eg.db"
        }
        if(species == "human"){
                OrgDb <- "org.Hs.eg.db"
        }
        cat(paste('The species is', species, "\n"))
        ## GO ----
        if(length(up_sig) > 0){
                up_ego <- enrichGO(gene          = up_sig,
                                   OrgDb         = OrgDb,
                                   ont           = ont,
                                   keyType       = "SYMBOL",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.1,
                                   qvalueCutoff  = 0.1)
        }
        
        if(length(down_ego) < 0){
                down_ego <- enrichGO(gene          = down_sig,
                                     OrgDb         = OrgDb,
                                     ont           = ont,
                                     keyType       = "SYMBOL",
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.1,
                                     qvalueCutoff  = 0.1)
        }
        
        
        if(is.null(up_ego) & is.null(down_ego)){
                break()     
        }else if(is.null(up_ego)){
                ego.list <- list(down=down_ego)   
        }else if(is.null(down_ego)){
                ego.list <- list(down=up_ego)   
        }else{
                ego.list <- list(up=up_ego, down=down_ego)
        }
        
        ego.list.sim <- ego.list %>% lapply(simplify, cutoff=simplify.cutoff)
        
        return(ego.list.sim)}

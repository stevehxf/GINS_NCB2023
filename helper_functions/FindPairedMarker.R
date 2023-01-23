FindPairedMarker <- function(seurat.object, 
                             split.by, group.by, 
                             ident.1, ident.2, 
                             min.pct = 0.3, 
                             logfc.threshold = 0.1, ...) {
        ## split by cell type and group by treatment would compare treatment vs control across all the cell types
        ## library
        require(tidyverse)
        require(Seurat)
        ## split idents (e.g. split by cell.type)
        Idents(seurat.object) <- group.by
        
        cluster.list <-  SplitObject(seurat.object, split.by = split.by)
        ## subset the cluster.list to make sure each group contains at least 10 cells, other wise return NULL
        meta.list <- cluster.list
        for (i in seq_along(cluster.list)) {
                meta.list[[i]] <- cluster.list[[i]]@meta.data
        }
        
        list.idx <- meta.list
        for(i in seq_along(meta.list)){
                list.idx[[i]] <- all(table(meta.list[[i]][[group.by]]) > 10)
        }
        idx <- unlist(list.idx)
        cluster.list <- cluster.list[idx]
        
        ## Differential expression
        DEG.list <- cluster.list %>% lapply(FindMarkers, ident.1 = ident.1, ident.2 = ident.2, min.pct = min.pct, logfc.threshold = logfc.threshold, ...)
        
        ## Add a variable for cluster
        for(i in seq_along(DEG.list)){
                DEG.list[[i]]$cluster <- names(DEG.list)[i]
        }
        return(DEG.list)
}

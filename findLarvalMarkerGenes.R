
library(Seurat)
library(HandyPack)

rm(list=ls())
## ####################################################
## ####################################################
cleanupCellTypes = function(types)
{
    ugh = c(' ',':','/','-')
    for(u in ugh)
        types = str_replace_all(types,u,'_')

    return(types)
}

## ####################################################
filesIn = c(merged='intermediate/larval_CellType.rds',
            integrated='intermediate/larval_integrated_CellType.rds')

for(n in names(filesIn))
{
    dir = nameAndMakeDir(paste0('larvalMarkerGenes_',n))

    f = readRDS(filesIn[n])
    f = f[,f$CellType != '*']
    f$CellType = cleanupCellTypes(f$CellType)

    clusters = unique(f$CellType)
    
    for(cluster in clusters)
    {
        Tic(paste(n,cluster))
        
        groupMarkerDF = FindMarkers(f,
                                    assay=f@active.assay,
                                    only.pos=TRUE,
                                    group.by=f$CellType,
                                    ident.1=cluster)
        groupMarkerDF = cbind(data.frame(gene=rownames(groupMarkerDF),
                                         stringsAsFactors=FALSE),
                              groupMarkerDF)
        
        cutoff = 0.05
        idx = groupMarkerDF$p_val_adj <= cutoff
        groupMarkerDF = groupMarkerDF[idx,]
        
        if(nrow(groupMarkerDF) == 0)
            next
        
        fileName = paste0(dir,
                          '/markerGenes_',
                          cluster,
                          '.txt')
        
        Write.Table(groupMarkerDF,
                    fileName)
        toc()
    }
}
 

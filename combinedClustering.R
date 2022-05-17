
library(Seurat)
library(ggplot2)
library(HandyPack)
library(tictoc)
library(stringr)
library(plotly)

rm(list=ls())
graphics.off()

## ####################################################
## ####################################################
combined = readRDS('intermediate/combined.rds')

## Get the orig.ident figure:
saveDir = nameAndMakeDir('clusterCombined')
df = FetchData(combined,c('UMAP_1','UMAP_2','orig.ident','CellType'))
g = ggplot(df,aes(x=UMAP_1,y=UMAP_2,color=orig.ident)) +
    geom_point() +
    ggtitle('Combined, data source on UMAP')
print(g)
GGSave(plot=g,
       filename='clusterCombined/orig.identOnUMAP.jpg',
       height=10,width=10)

N = length(unique(df$CellType))
colors = polychrome(N)

h = ggplot(df,aes(x=UMAP_1,y=UMAP_2,color=CellType)) +
    geom_point() +
    scale_color_manual(values=colors) +
    ggtitle('Combined, CellType on UMAP')
dev.new()
print(h)
fileName = 'clusterCombined/CellTypeOnUMAP.jpg'
GGSave(plot=h,
       filename=fileName,
       height=12,width=16)
p = ggplotly(h)
print(p)
fileName = str_replace(fileName,'jpg','html')
htmlwidgets::saveWidget(as_widget(p), fileName)

## ####################################################
resolutions = c(.5,1,1.5)

for(res in resolutions)
{
    f = FindClusters(combined,resolution=res)
    saveDir = c('clusterCombined',
                paste0('combined_',res))
    saveDir = nameAndMakeDir(saveDir)
    saveRDS(f,
            paste0(saveDir,'/combined_',res,'.rds'))

    df = FetchData(f,c('UMAP_1','UMAP_2','seurat_clusters','orig.ident','CellType'))
    N = length(unique(df$seurat_clusters))
    colors = polychrome(N)

    g = ggplot(df,aes(x=UMAP_1,y=UMAP_2,color=seurat_clusters,label=CellType)) +
        geom_point() +
        scale_color_manual(values=colors) +
        ggtitle(paste('Combined, seurat_clusters on UMAP, res',res))
    dev.new()
    print(g)
    fileName = paste0(saveDir,'/combined_seurat_clustersOnUMAP_res',res,'.jpg')
    GGSave(plot=g,
           filename=fileName,
           heigh=12,width=16)
    q = ggplotly(g,tooltips=c('seurat_clusters','label'))
    print(q)

    fileName = str_replace(fileName,'jpg','html')
    htmlwidgets::saveWidget(as_widget(q), fileName)

    saveDir = c('clusterCombined',
                paste0('combined_',res),
                'markerGenes')
    saveDir = nameAndMakeDir(saveDir)

    
    clusters = unique(f$seurat_clusters)
    
    for(cluster in clusters)
    {
        Tic(paste(res,cluster))
        
        groupMarkerDF = FindMarkers(f,
                                    assay=f@active.assay,
                                    only.pos=TRUE,
                                    group.by=f$seurat_clusters,
                                    ident.1=cluster)
        groupMarkerDF = cbind(data.frame(gene=rownames(groupMarkerDF),
                                         stringsAsFactors=FALSE),
                              groupMarkerDF)
        
        cutoff = 0.05
        idx = groupMarkerDF$p_val_adj <= cutoff
        groupMarkerDF = groupMarkerDF[idx,]
        
        if(nrow(groupMarkerDF) == 0)
            next
        
        fileName = paste0(saveDir,
                          '/markerGenes_',
                          cluster,
                          '.txt')
        
        Write.Table(groupMarkerDF,
                    fileName)
        toc()
    }
}

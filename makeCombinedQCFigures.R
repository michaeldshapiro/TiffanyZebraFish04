
library(Seurat)
library(HandyPack)
library(ggplot2)
library(plotly)
library(stringr)
library(Polychrome)

rm(list=ls())
graphics.off()

source('CopyWare.R')
source('IntegrationWare.R')
## ####################################################
## ####################################################

saveDir = nameAndMakeDir('combinedQCFigures')

filesIn = Sys.glob('intermediate/combined*rds')

for(file in filesIn)
{
    combined = readRDS(file)

    
    df = FetchData(combined,c('UMAP_1','UMAP_2','CellType','orig.ident'))
    types = unique(df$CellType)
    types = types[order(types)]

    ## Make orig.ident figures:
    title = file
    title = str_replace(title,
                        'intermediate/',
                        'orig.ident ')
    
    g = ggplot(df,aes(x=UMAP_1,y=UMAP_2,color=orig.ident)) +
        geom_point() +
        ggtitle(title)
    dev.new()
    print(g)
    fileName = str_replace(file,
                           'intermediate/',
                           'combinedQCFigures/orig.identOnUMAP_')
    fileName = str_replace(fileName,'rds','jpg')
    
    GGSave(plot=g,
           filename=fileName)

    title = str_replace(title,'orig.ident','CellType')
    h = ggplot(df,aes(x=UMAP_1,y=UMAP_2,color=CellType)) +
        geom_point() +
        ggtitle(title)
    dev.new()
    print(h)
    fileName = str_replace(file,
                           'intermediate/',
                           'combinedQCFigures/CellTypeOnUMAP_')
    fileName = str_replace(fileName,'rds','jpg')
    
    GGSave(plot=h,
           filename=fileName,
           height=12,width=16)

    p = ggplotly(h)
    print(p)
    fileName = str_replace(fileName,'jpg','html')
    htmlwidgets::saveWidget(as_widget(p), fileName)
}



library(Seurat)
library(HandyPack)
library(ggplot2)
library(plotly)
library(stringr)

rm(list=ls())
graphics.off()

source('CopyWare.R')
source('IntegrationWare.R')
## ####################################################
## ####################################################

saveDir = nameAndMakeDir('larvalQCFigures')

larval = readRDS('intermediate/larval_CellType.rds')
oldLarval = getOld('larval')
oldLarval = runBasicAnalyses(oldLarval)
types = unique(oldLarval$CellType)
types = types[order(types)]
types = c('*',types)

larval$CellType = factor(larval$CellType,levels=types)
oldLarval$CellType = factor(oldLarval$CellType,levels=types)

df = FetchData(larval,c('UMAP_1','UMAP_2','CellType','orig.ident'))
oldDF = FetchData(oldLarval,c('UMAP_1','UMAP_2','CellType','orig.ident'))


## New orig.ident:
N = length(unique(df$orig.ident))
colors = polychrome(N)
gOrig = ggplot(df,aes(x=UMAP_1,y=UMAP_2,color=orig.ident)) +
    geom_point() +
    scale_color_manual(values=colors) +
    ggtitle('Larval data showing data set')

print(gOrig)
pOrig = ggplotly(gOrig)
print(pOrig)

fileName = paste0(saveDir,
                  '/larval_orig.ident.jpg')
ggsave(plot=gOrig,
       filename=fileName)
fileName = str_replace(fileName,'jpg','html')
htmlwidgets::saveWidget(as_widget(pOrig), fileName)

## New CellType:
NType = length(unique(df$CellType))
typeColors = polychrome(NType)

gType = ggplot(df,aes(x=UMAP_1,y=UMAP_2,color=CellType,label=orig.ident)) +
    geom_point() +
    scale_color_manual(values=typeColors) +
    ggtitle('Larval data showing CellType')

dev.new()
print(gType)
pType = ggplotly(gType,tooltips=c('CellType','orig.ident'))
print(pType)


fileName = paste0(saveDir,
                  '/larval_CellType.jpg')
ggsave(plot=gType,
       filename=fileName,
       width=20,height=12,units='in')
fileName = str_replace(fileName,'jpg','html')
htmlwidgets::saveWidget(as_widget(pType), fileName)

## ####################################################

## Old orig.ident:
N = length(unique(oldDF$orig.ident))
colors = polychrome(N)
hOrig = ggplot(oldDF,aes(x=UMAP_1,y=UMAP_2,color=orig.ident)) +
    geom_point() +
    scale_color_manual(values=colors) +
    ggtitle('Old larval data showing data set')

print(hOrig)
qOrig = ggplotly(hOrig)
dev.new()
print(qOrig)

fileName = paste0(saveDir,
                  '/oldLarval_orig.ident.jpg')
ggsave(plot=hOrig,
       filename=fileName)
fileName = str_replace(fileName,'jpg','html')
htmlwidgets::saveWidget(as_widget(qOrig), fileName)

## Old CellType:
hType = ggplot(oldDF,aes(x=UMAP_1,y=UMAP_2,color=CellType,label=orig.ident)) +
    geom_point() +
    scale_color_manual(values=typeColors) +
    ggtitle('Old larval data showing CellType')

dev.new()
print(hType)
qType = ggplotly(hType,tooltips=c('CellType','orig.ident'))
print(qType)


fileName = paste0(saveDir,
                  '/oldLarval_CellType.jpg')
ggsave(plot=hType,
       filename=fileName,
       width=20,height=12,units='in')
fileName = str_replace(fileName,'jpg','html')
htmlwidgets::saveWidget(as_widget(qType), fileName)

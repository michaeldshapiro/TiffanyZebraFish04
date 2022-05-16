
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

saveDir = nameAndMakeDir('adultQCFigures')

adult = readRDS('intermediate/adult_CellType.rds')
types = unique(adult$CellType)
types = types[order(types)]


adult$CellType = factor(adult$CellType,levels=types)

oldAdult = getOld('adult')
oldAdult$CellType = oldAdult$shortName

df = FetchData(adult,c('UMAP_1','UMAP_2','CellType','orig.ident'))
oldDF = FetchData(oldAdult,c('UMAP_1','UMAP_2','CellType','orig.ident'))

## New orig.ident:
N = length(unique(df$orig.ident))
colors = polychrome(N)
gOrig = ggplot(df,aes(x=UMAP_1,y=UMAP_2,color=orig.ident)) +
    geom_point() +
    scale_color_manual(values=colors) +
    ggtitle('Adult data showing data set')

print(gOrig)
pOrig = ggplotly(gOrig)
print(pOrig)

fileName = paste0(saveDir,
                  '/adult_orig.ident.jpg')
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
    ggtitle('Adult data showing CellType')

dev.new()
print(gType)
pType = ggplotly(gType,tooltips=c('CellType','orig.ident'))
print(pType)


fileName = paste0(saveDir,
                  '/adult_CellType.jpg')
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
    ggtitle('Old adult data showing data set')

print(hOrig)
qOrig = ggplotly(hOrig)
dev.new()
print(qOrig)

fileName = paste0(saveDir,
                  '/oldAdult_orig.ident.jpg')
ggsave(plot=hOrig,
       filename=fileName)
fileName = str_replace(fileName,'jpg','html')
htmlwidgets::saveWidget(as_widget(qOrig), fileName)

## Old CellType:
typeColors = colors[2:length(typeColors)]
hType = ggplot(oldDF,aes(x=UMAP_1,y=UMAP_2,color=CellType,label=orig.ident)) +
    geom_point() +
    scale_color_manual(values=typeColors) +
    ggtitle('Old adult data showing CellType')

dev.new()
print(hType)
qType = ggplotly(hType,tooltips=c('CellType','orig.ident'))
print(qType)


fileName = paste0(saveDir,
                  '/oldAdult_CellType.jpg')
ggsave(plot=hType,
       filename=fileName,
       width=20,height=12,units='in')
fileName = str_replace(fileName,'jpg','html')
htmlwidgets::saveWidget(as_widget(qType), fileName)

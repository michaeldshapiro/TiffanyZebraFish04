
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

saveDir = nameAndMakeDir('clusterCombined')

df = FetchData(combined,c('PC_1','PC_2','PC_3','orig.ident','CellType'))

N = length(unique(df$CellType))
colors = polychrome(N)

g = ggplot(df,aes(x=PC_1,y=PC_2,color=CellType)) +
    geom_point() +
    scale_color_manual(values=colors) +
    ggtitle('Combined, CellType on PCA')
print(g)
p = ggplotly(g)
print(p)
fileName = 'clusterCombined/CellTypeOnPCs1_2.jpg'
GGSave(plot=g,
       filename=fileName,
       height=12,width=16)
fileName = str_replace(fileName,'jpg','html')
htmlwidgets::saveWidget(as_widget(p), fileName)


h = ggplot(df,aes(x=PC_1,y=PC_3,color=CellType)) +
    geom_point() +
    scale_color_manual(values=colors) +    
    ggtitle('Combined, CellType on PCA')
dev.new()
print(h)
q = ggplotly(h)
print(q)
fileName = 'clusterCombined/CellTypeOnPCs1_3.jpg'

GGSave(plot=h,
       filename=fileName,
       height=12,width=16)
fileName = str_replace(fileName,'jpg','html')
htmlwidgets::saveWidget(as_widget(q), fileName)

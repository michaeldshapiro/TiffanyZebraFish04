
library(HandyPack)
library(Seurat)

source('CopyWare.R')
## ####################################################
## ####################################################

larval = readRDS('intermediate/larval_CellType.rds')
oldLarval = getOld('larval')

df = data.frame(table(larval$CellType))
names(df) = c('CellType','count')
Write.Table(df,
            'tables/larvalCellTypes_new.txt')

oldDF = data.frame(table(oldLarval$CellType))
names(oldDF) = c('CellType','count')
Write.Table(oldDF,
            'tables/larvalCellTypes_old.txt')


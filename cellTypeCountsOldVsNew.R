
library(Seurat)
library(ggplot2)

rm(list=ls())
graphics.off()

source('CopyWare.R')
## ####################################################
## ####################################################
oldCombined = getOld('combined')
newCombined = readRDS('intermediate/combined_integratedLarval_keepCells2_adult20Only.rds')

